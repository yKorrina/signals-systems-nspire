-- ===================== OUTPUT COLLECTOR =====================
local Out = {}
Out.__index = Out
local function newOut() return setmetatable({lines = {}, plots = {}}, Out) end
function Out:h(s) self.lines[#self.lines + 1] = {k = "h", s = s} end
function Out:t(s) self.lines[#self.lines + 1] = {k = "t", s = s} end
function Out:r(s) self.lines[#self.lines + 1] = {k = "r", s = s} end
function Out:e(s) self.lines[#self.lines + 1] = {k = "e", s = s} end
function Out:plot(p) self.plots[#self.plots + 1] = p end

-- ===================== FORWARD TRANSFORMS (CT) =====================
-- F term: c * e^{-d v} * (v - p)^(-m)   (m may be negative => (v)^|m| with p = 0)
-- freq impulse: {c, w0} -> c δ(ω - w0)
local function addFT(list, c, d, p, m)
  c = tocx(c)
  if cabs(c) < 1e-13 then return end
  for _, t in ipairs(list) do
    if abs(t.d - d) < 1e-10 and t.m == m and near(t.p, p, 1e-9) then t.c = t.c + c; return end
  end
  list[#list + 1] = {c = c, d = d, p = tocx(p), m = m}
end
local function cleanFT(list)
  local r = {}
  for _, t in ipairs(list) do if cabs(t.c) > 1e-11 then t.c = ccln(t.c); t.d = cln(t.d); r[#r + 1] = t end end
  table.sort(r, function(a, b) if abs(a.d - b.d) > 1e-12 then return a.d < b.d end; if a.m ~= b.m then return a.m < b.m end; return a.p.re < b.p.re end)
  return r
end

-- describe a monomial piece c t^k e^{st} on (lo,hi)
local function pieceStr(c, k, s, lo, hi, var, dt)
  local S = {M(c, k, s)}
  local e = segStr(S, var, dt)
  return e .. ",  " .. rangeStr(lo, hi, var, dt)
end

-- CT forward transform. mode "F" (Fourier) or "L" (Laplace)
-- returns terms, impulses, roc{lo,hi}, steps(list of strings)
local function ctTransform(sig, mode)
  local terms, imps, steps = {}, {}, {}
  local rlo, rhi = -INF, INF
  local vs = (mode == "F") and "jω" or "s"
  for _, im in ipairs(sig.imps) do
    addFT(terms, im.c, im.t0, ZERO, -im.m)
    steps[#steps + 1] = impStr(im, "t") .. "  →  " .. cmul(im.c, (im.m > 0 and ("(" .. vs .. ")" .. (im.m > 1 and ("^" .. im.m) or "")) or "") .. (abs(im.t0) > 1e-12 and ("e^(" .. expoArg(-im.t0, (mode == "F") and "jω" or "s") .. ")") or "")) .. "   (sifting)"
  end
  for _, pc in ipairs(sig:monoWindows()) do
    for _, m in ipairs(pc.S) do
      local lo, hi = pc.lo, pc.hi
      local s, k, c = m.s, m.k, m.c
      local pl = {}
      local note = ""
      local function anti(at, sign)
        -- adds sign * [antiderivative of c t^k e^{(s-v)t}] at t=at
        local fk = fact(k)
        for i = 0, k do
          local coef = c * cexp(s * at) * cx(at ^ (k - i)) * (((i % 2 == 1) and -1 or 1) * fk / fact(k - i)) * (((i + 1) % 2 == 1) and -1 or 1)
          addFT(pl, coef * sign, at, s, i + 1)
        end
      end
      if lo > -INF and hi < INF then
        anti(hi, 1); anti(lo, -1)
      elseif lo > -INF then
        if mode == "F" then
          if s.re < -1e-12 then anti(lo, -1)
          elseif abs(s.re) <= 1e-12 and k == 0 then
            anti(lo, -1); imps[#imps + 1] = {c = c * pi, w0 = s.im}
            note = " (+ πδ term)"
          else error("Fourier transform does not exist (signal grows); try Laplace", 0) end
        else
          anti(lo, -1); rlo = math.max(rlo, s.re)
        end
      elseif hi < INF then
        if mode == "F" then
          if s.re > 1e-12 then anti(hi, 1)
          elseif abs(s.re) <= 1e-12 and k == 0 then
            anti(hi, 1); imps[#imps + 1] = {c = c * pi, w0 = s.im}
            note = " (+ πδ term)"
          else error("Fourier transform does not exist (signal grows); try Laplace", 0) end
        else
          anti(hi, 1); rhi = math.min(rhi, s.re)
        end
      else
        if mode == "F" and abs(s.re) <= 1e-12 and k == 0 then
          imps[#imps + 1] = {c = c * 2 * pi, w0 = s.im}
          note = " = 2π·c·δ(ω - ω0)"
        else
          error("transform of a two-sided " .. ((mode == "F") and "non-periodic growing signal does not exist" or "signal has no ROC"), 0)
        end
      end
      local lim1 = (lo == -INF) and "-∞" or fmtR(lo)
      local lim2 = (hi == INF) and "∞" or fmtR(hi)
      local integrand = segStr({M(c, k, s - ((mode == "F") and J or ONE) * 0)}, "t", false)
      local part = pl
      for _, tt in ipairs(pl) do addFT(terms, tt.c, tt.d, tt.p, tt.m) end
      steps[#steps + 1] = "∫[" .. lim1 .. "→" .. lim2 .. "] " .. attach(segStr({M(c, k, s)}, "t", false), "") .. "·e^(-" .. vs .. "t) dt" .. note
      local lam = vfac(s, mode)
      local cs = segStr({M(c, k, 0)}, "t", false)
      local cpre = (cs == "1") and "" or ((topSum(cs) and ("(" .. cs .. ")") or cs) .. "·")
      if cs == "-1" then cpre = "-" end
      if k > 0 then cpre = segStr({M(c, 0, 0)}, "t", false); cpre = (cpre == "1") and "" or (cpre == "-1" and "-" or cpre .. "·") end
      steps[#steps + 1] = "   = " .. cpre .. "∫[" .. lim1 .. "→" .. lim2 .. "] " .. (k > 0 and (powStr("t", k) .. "·") or "") .. "e^(-" .. lam .. "t) dt   (combine exponents)"
      if k == 0 then
        steps[#steps + 1] = "   = " .. cpre .. "[ e^(-" .. lam .. "t) / (-" .. lam .. ") ] from " .. lim1 .. " to " .. lim2
      else
        steps[#steps + 1] = "   use ∫t^k e^(λt)dt = e^(λt)·Σ(-1)^i k!/(k-i)!·t^(k-i)/λ^(i+1),  λ = -" .. lam .. "  (by parts)"
      end
      local sr = (mode == "F") and "Re" or "Re"
      if hi == INF then steps[#steps + 1] = "   at t→∞ the term → 0 because " .. ((mode == "F") and ("Re(" .. fmtC(s) .. ") = " .. fmtR(s.re) .. " < 0") or ("Re{s} > " .. fmtR(s.re) .. " (ROC)")) end
      if lo == -INF then steps[#steps + 1] = "   at t→-∞ the term → 0 because " .. ((mode == "F") and ("Re(" .. fmtC(s) .. ") = " .. fmtR(s.re) .. " > 0") or ("Re{s} < " .. fmtR(s.re) .. " (ROC)")) end
      if lo == 0 and hi == INF then
        steps[#steps + 1] = "   (table pair: t^k e^(-at)u(t) <-> k!/(a+" .. ((mode == "F") and "jω" or "s") .. ")^(k+1))"
      end
      local pstrs = {}
      for _, tt in ipairs(cleanFT(part)) do pstrs[#pstrs + 1] = ftermStr(tt, mode) end
      if #pstrs > 0 then steps[#steps + 1] = "   = " .. joinsum(pstrs) end
    end
  end
  return cleanFT(terms), imps, {lo = rlo, hi = rhi}, steps
end

-- print (v - p) factor
local function vfac(p, mode)
  p = tocx(p)
  local v = (mode == "F") and "jω" or "s"
  if isz(p, 1e-12) then return v, true end
  if mode == "F" then
    local a = -p.re
    local w0 = p.im
    if abs(w0) < 1e-12 then
      if a > 0 then return "(" .. fmtR(a) .. " + jω)" end
      return "(jω - " .. fmtR(-a) .. ")"
    end
    local inner = "j(ω " .. (w0 > 0 and "- " or "+ ") .. fmtR(abs(w0)) .. ")"
    if abs(a) < 1e-12 then return "(" .. inner .. ")" end
    return "(" .. fmtR(a) .. " + " .. inner .. ")"
  end
  local out = "s"
  local a, b = -p.re, -p.im
  if abs(a) > 1e-12 then out = out .. (a > 0 and " + " or " - ") .. fmtR(abs(a)) end
  if abs(b) > 1e-12 then out = out .. (b > 0 and " + " or " - ") .. fmtImag(abs(b)) end
  return "(" .. out .. ")"
end
local function delayStr(d, mode)
  if abs(d) < 1e-12 then return "" end
  if mode == "F" then return "e^(" .. ((d > 0) and "-" or "") .. "j" .. ((abs(d) == 1) and "" or fmtR(abs(d))) .. "ω)" end
  if mode == "L" then return "e^(" .. expoArg(-d, "s") .. ")" end
end
local function ftermStr(t, mode)
  local num = delayStr(t.d, mode)
  local den = ""
  if t.m > 0 then
    local f, bare = vfac(t.p, mode)
    if t.m > 1 then
      if bare and mode == "F" then f = "(" .. f .. ")" end
      f = f .. "^" .. t.m
    elseif bare and mode == "F" then f = "(" .. f .. ")" end
    den = f
  elseif t.m < 0 then
    local v = (mode == "F") and "jω" or "s"
    num = joinf(((-t.m) > 1) and ("(" .. v .. ")^" .. (-t.m)) or v, num)
  end
  if den == "" then return cmul(t.c, num == "" and "1" or num) end
  local c = t.c
  local neg = false
  if isreal(c) and c.re < 0 then neg = true; c = -c end
  local cs = fmtC(c)
  local top
  if num == "" then top = cs
  else
    if cs == "1" then top = num else
      if needpar(cs) or not isreal(c) then cs = "(" .. cs .. ")" end
      top = cs .. num
    end
  end
  local p, q = cs:match("^(%d+)/(%d+)$")
  if p and num == "" then return (neg and "-" or "") .. p .. "/(" .. q .. den .. ")" end
  return (neg and "-" or "") .. top .. "/" .. den
end
local function fimpStr(im, isDT)
  local w = isDT and "Ω" or "ω"
  local a = (abs(im.w0) < 1e-12) and w or (w .. (im.w0 > 0 and " - " or " + ") .. fmtR(abs(im.w0)))
  local s = cmul(im.c, "δ(" .. a .. ")")
  return s
end
-- combine terms with same delay into single rational N/D (in v)
local function combineFT(terms)
  local groups = {}
  for _, t in ipairs(terms) do
    local g
    for _, x in ipairs(groups) do if abs(x.d - t.d) < 1e-10 then g = x end end
    if not g then g = {d = t.d, N = pconst(0), D = pconst(1)}; groups[#groups + 1] = g end
    local tn, td
    if t.m >= 0 then
      tn = pconst(t.c); td = pconst(1)
      for _ = 1, t.m do td = pmul(td, pnew({-t.p, 1})) end
    else
      tn = pconst(t.c); for _ = 1, -t.m do tn = pmul(tn, pnew({0, 1})) end; td = pconst(1)
    end
    -- g.N/g.D + tn/td  (avoid repeated factors: use lcm-ish by checking divisibility)
    local q, r = pdivmod(g.D, td)
    if piszero(r) and g.D.n >= td.n then
      g.N = padd(g.N, pmul(tn, q))
    else
      local q2, r2 = pdivmod(td, g.D)
      if piszero(r2) then g.N = padd(pmul(g.N, q2), tn); g.D = td
      else g.N = padd(pmul(g.N, td), pmul(tn, g.D)); g.D = pmul(g.D, td) end
    end
  end
  return groups
end
local function ratStr(N, D, mode, offN)
  local v = (mode == "F") and "jω" or ((mode == "Z") and "z^-1" or ((mode == "D") and "e^(-jΩ)" or "s"))
  local wrap = (mode == "F") or (mode == "Z") or (mode == "D")
  local lead = D[pdeg(D)]
  -- normalize so that denominator lowest/ leading coefficient is 1
  local norm = (mode == "Z" or mode == "D") and D[0] or lead
  if cabs(norm) < 1e-12 then norm = lead end
  N = pscale(N, ONE / norm); D = pscale(D, ONE / norm)
  local ns, ds
  if mode == "Z" or mode == "D" then ns, ds = pstrQ(N, mode, offN or 0), pstrQ(D, mode)
  else ns, ds = pstr(N, v, wrap), pstr(D, v, wrap) end
  if pdeg(D) == 0 and near(D[0], ONE) then return ns end
  if ns:find(" ") then ns = "(" .. ns .. ")" end
  return ns .. "/(" .. ds .. ")"
end

local function ftSumStr(terms, imps, mode, isDT)
  local parts = {}
  for _, t in ipairs(terms) do parts[#parts + 1] = (isDT and dtermStr or ftermStr)(t, mode) end
  for _, im in ipairs(imps or {}) do parts[#parts + 1] = fimpStr(im, isDT) end
  return joinsum(parts)
end
-- evaluate CT transform at v (complex)
local function ftEval(terms, v)
  local s = ZERO
  for _, t in ipairs(terms) do
    local f = t.c * cexp(-v * t.d)
    if t.m > 0 then f = f / cpow(v - t.p, t.m) elseif t.m < 0 then f = f * cpow(v, -t.m) end
    s = s + f
  end
  return s
end

-- ===================== FORWARD TRANSFORMS (DT) =====================
-- D term: c * q^d * (1 - p q^dir)^(-m),  q = e^{-jΩ} (or z^-1)
local function addDT(list, c, d, p, m, dir)
  c = tocx(c); d = d + 0
  if cabs(c) < 1e-13 then return end
  if m == 0 then p = ZERO; dir = 1 end
  for _, t in ipairs(list) do
    if t.d == d and t.m == m and t.dir == dir and near(t.p, p, 1e-9) then t.c = t.c + c; return end
  end
  list[#list + 1] = {c = c, d = d, p = tocx(p), m = m, dir = dir}
end
-- coefficients beta_m with (x+a)^k = sum_{m=1}^{k+1} beta_m C(x+m-1, m-1)
local function binomBasis(k, a)
  local A, b = {}, {}
  for x = 0, k do
    A[x + 1] = {}
    for m = 1, k + 1 do A[x + 1][m] = cx(binom(x + m - 1, m - 1)) end
    b[x + 1] = cx((x + a) ^ k)
  end
  return solveLin(A, b)
end
local function dtTransform(sig, mode)
  local terms, imps, steps = {}, {}, {}
  local rlo, rhi = 0, INF
  local flo, fhi = sig:dtSupport()
  if flo and fhi - flo <= 40 then
    -- finite-length: write the sum out term by term
    local qs = function(n)
      if n == 0 then return "" end
      if mode == "Z" then return (n == -1) and "z" or ("z^" .. ((n > 0) and "-" or "") .. abs(n)) end
      return "e^(" .. ((n > 0) and "-" or "") .. "j" .. ((abs(n) == 1) and "" or abs(n)) .. "Ω)"
    end
    steps[#steps + 1] = "x is nonzero only for " .. fmtR(flo) .. " ≤ n ≤ " .. fmtR(fhi) .. ", so the sum has a finite number of terms:"
    local raw, vals = {}, {}
    for n = flo, fhi do
      local v = sig:eval(n)
      raw[#raw + 1] = "x[" .. n .. "]" .. ((n ~= 0) and ("·" .. qs(n)) or "")
      vals[#vals + 1] = "x[" .. n .. "]=" .. fmtC(v)
      addDT(terms, v, n, ZERO, 0, 1)
    end
    steps[#steps + 1] = "   " .. table.concat(raw, " + ")
    steps[#steps + 1] = "   values: " .. table.concat(vals, ", ")
    local parts = {}
    for n = flo, fhi do
      local v = sig:eval(n)
      if cabs(v) > 1e-12 then parts[#parts + 1] = cmul(v, (n == 0) and "1" or qs(n)) end
    end
    steps[#steps + 1] = "   = " .. joinsum(parts) .. "   (each x[k]δ[n-k] gives x[k]·" .. ((mode == "Z") and "z^-k" or "e^(-jΩk)") .. ")"
    local r = {}
    for _, t in ipairs(terms) do if cabs(t.c) > 1e-11 then t.c = ccln(t.c); r[#r + 1] = t end end
    table.sort(r, function(a, b) return a.d < b.d end)
    return r, imps, {lo = 0, hi = INF}, steps
  end
  for _, pc in ipairs(sig:monoWindows()) do
    for _, m in ipairs(pc.S) do
      local lo, hi = pc.lo, pc.hi
      local r, k, c = m.s, m.k, m.c
      local desc = pieceStr(c, k, r, lo, hi, "n", true)
      local pl = {}
      local q = (mode == "Z") and "z^-1" or "e^(-jΩ)"
      local qi = (mode == "Z") and "z" or "e^(jΩ)"
      local rs = fmtC(r); if not rs:find("^[%d%.]+$") then rs = "(" .. rs .. ")" end
      local w = (near(r, ONE) and q or (rs .. q))
      if lo > -INF and hi < INF then
        for n = lo, hi - 1 do addDT(pl, c * cx(n ^ k) * cpow(r, n), n, ZERO, 0, 1) end
        steps[#steps + 1] = "Σ[n=" .. fmtR(lo) .. "→" .. fmtR(hi - 1) .. "] " .. desc:match("^(.-),") .. "·" .. ((mode == "Z") and "z^-n" or "e^(-jΩn)") .. "   (finite: write out each term)"
        local lst = {}
        for n = lo, math.min(hi - 1, lo + 5) do lst[#lst + 1] = "x[" .. n .. "]=" .. fmtC(c * cx(n ^ k) * cpow(r, n)) end
        steps[#steps + 1] = "   " .. table.concat(lst, ", ") .. ((hi - lo > 6) and ", ..." or "")
        if k == 0 and hi - lo > 1 then
          steps[#steps + 1] = "   closed form (geometric): " .. cmul(c, "[" .. w .. "]^" .. fmtR(lo)) .. "·(1 - [" .. w .. "]^" .. fmtR(hi - lo) .. ")/(1 - " .. w .. ")"
        end
      elseif lo > -INF then
        local a = lo
        local ar = cabs(r)
        if mode == "D" and ar > 1 + 1e-12 then error("DTFT does not exist (|r| > 1); use Z-transform", 0) end
        if mode == "D" and abs(ar - 1) <= 1e-12 then
          if k > 0 then error("DTFT of n^k·r^n·u[n] with |r|=1 not supported", 0) end
          imps[#imps + 1] = {c = c * pi, w0 = carg(r), periodic = true}
        end
        local beta = binomBasis(k, a)
        for mm = 1, k + 1 do addDT(pl, c * cpow(r, a) * beta[mm], a, r, mm, 1) end
        rlo = math.max(rlo, ar)
        steps[#steps + 1] = "Σ[n=" .. fmtR(lo) .. "→∞] " .. desc:match("^(.-),") .. "·" .. ((mode == "Z") and "z^-n" or "e^(-jΩn)")
        if k == 0 then
          steps[#steps + 1] = "   = " .. cmul(c, "Σ[n=" .. fmtR(lo) .. "→∞] (" .. w .. ")^n") .. "   (ratio w = " .. w .. ")"
          steps[#steps + 1] = "   = " .. cmul(c, "w^" .. fmtR(lo) .. "/(1 - w)") .. "   (geometric series Σ w^n = w^a/(1-w), needs |w|<1 ⇔ " .. ((mode == "Z") and ("|z| > " .. fmtR(ar)) or ("|r| = " .. fmtR(ar) .. " < 1")) .. ")"
        else
          steps[#steps + 1] = "   n^k factor: use Σ C(n+m-1, m-1)w^n = 1/(1-w)^m (derivative of geometric series), w = " .. w
        end
      elseif hi < INF then
        local ap = 1 - hi
        local ar = cabs(r)
        if mode == "D" and ar < 1 - 1e-12 then error("DTFT does not exist (left-sided with |r| < 1)", 0) end
        if mode == "D" and abs(ar - 1) <= 1e-12 then
          if k > 0 then error("not supported", 0) end
          -- c r^n u[hi-1-n] = c r^n - c r^n u[n-hi]
          imps[#imps + 1] = {c = c * pi, w0 = carg(r), periodic = true}
          addDT(pl, -c * cpow(r, hi), hi, r, 1, 1)
        else
          local rinv = ONE / r
          local beta = binomBasis(k, ap)
          for mm = 1, k + 1 do
            addDT(pl, c * cx(((k % 2 == 1) and -1 or 1)) * beta[mm] * cpow(rinv, ap), -ap, rinv, mm, -1)
          end
        end
        rhi = math.min(rhi, ar)
        steps[#steps + 1] = "Σ[n=-∞→" .. fmtR(hi - 1) .. "] " .. desc:match("^(.-),") .. "·" .. ((mode == "Z") and "z^-n" or "e^(-jΩn)")
        steps[#steps + 1] = "   substitute m = -n: Σ[m=" .. fmtR(1 - hi) .. "→∞] (" .. fmtC(ONE / r) .. "·" .. qi .. ")^m" .. (k > 0 and " · (-m)^" .. k or "")
        steps[#steps + 1] = "   geometric series needs |" .. fmtC(ONE / r) .. qi .. "| < 1 ⇔ " .. ((mode == "Z") and ("|z| < " .. fmtR(ar)) or ("|r| = " .. fmtR(ar) .. " > 1"))
      else
        if mode == "D" and abs(cabs(r) - 1) <= 1e-12 and k == 0 then
          imps[#imps + 1] = {c = c * 2 * pi, w0 = carg(r), periodic = true}
          steps[#steps + 1] = desc .. "  (periodic complex exponential, table pair e^(jΩ0 n) <-> 2πΣδ(Ω-Ω0-2πl))"
          steps[#steps + 1] = "   → " .. cmul(c * 2 * pi, "δ(Ω " .. (carg(r) >= 0 and "- " or "+ ") .. fmtR(abs(carg(r))) .. ")") .. "  (+ copies every 2π)"
        else
          error("transform does not exist for two-sided " .. desc, 0)
        end
      end
      local ps = {}
      for _, t in ipairs(pl) do ps[#ps + 1] = dtermStr(t, mode); addDT(terms, t.c, t.d, t.p, t.m, t.dir) end
      if #ps > 0 then steps[#steps + 1] = "   = " .. joinsum(ps) end
    end
  end
  local r = {}
  for _, t in ipairs(terms) do if cabs(t.c) > 1e-11 then t.c = ccln(t.c); r[#r + 1] = t end end
  table.sort(r, function(a, b) if a.m ~= b.m then return a.m < b.m end; return a.d < b.d end)
  return r, imps, {lo = rlo, hi = rhi}, steps
end
function dtermStr(t, mode)
  local function qp(d)
    if d == 0 then return "" end
    if mode == "Z" then
      if d == -1 then return "z" end
      return "z^" .. ((d > 0) and "-" or "") .. tostring(abs(d))
    end
    return "e^(" .. ((d > 0) and "-" or "") .. "j" .. (abs(d) == 1 and "" or tostring(abs(d))) .. "Ω)"
  end
  local num = qp(t.d)
  if t.m == 0 then return cmul(t.c, num == "" and "1" or num) end
  local inner = "1 - " .. cmul(t.p, qp(t.dir))
  inner = inner:gsub("%- %-", "+ ")
  local den = "(" .. inner .. ")" .. ((t.m > 1) and ("^" .. t.m) or "")
  local c = t.c
  local neg = isreal(c) and c.re < 0
  if neg then c = -c end
  local cs = fmtC(c)
  local top
  if num == "" then top = cs:find("/") and ("(" .. cs .. ")") or cs else
    if cs == "1" then top = num else top = ((needpar(cs) or not isreal(c)) and ("(" .. cs .. ")") or cs) .. num end
  end
  return (neg and "-" or "") .. top .. "/" .. den
end
local function dtEval(terms, Om, zmode)
  local q = zmode and (ONE / Om) or cexp(cx(0, -Om))
  local s = ZERO
  for _, t in ipairs(terms) do
    local f = t.c * cpow(q, t.d)
    if t.m > 0 then f = f / cpow(ONE - t.p * cpow(q, t.dir), t.m) end
    s = s + f
  end
  return s
end
-- combine DT terms into single rational in q
local function combineDT(terms)
  -- multiply everything to polynomials in q with offset
  local N, D = pconst(0), pconst(1)
  local minpow = 0
  for _, t in ipairs(terms) do if t.d < minpow then minpow = t.d end; if t.dir == -1 then minpow = math.min(minpow, t.d - t.m) end end
  for _, t in ipairs(terms) do
    local tn = {n = 0, [0] = ZERO}
    local e = t.d - minpow
    tn = {n = e}; for i = 0, e do tn[i] = ZERO end; tn[e] = t.c
    local td = pconst(1)
    if t.m > 0 then
      local f
      if t.dir == 1 then f = pnew({1, -t.p}) else f = pnew({-t.p, 1}) end -- (1 - p q^-1) = (q - p)/q
      for _ = 1, t.m do td = pmul(td, f) end
      if t.dir == -1 then local sh = {n = t.m}; for i = 0, t.m do sh[i] = ZERO end; sh[t.m] = ONE; tn = pmul(tn, sh) end
    end
    local q, r = pdivmod(D, td)
    if piszero(r) and D.n >= td.n then N = padd(N, pmul(tn, q))
    else
      local q2, r2 = pdivmod(td, D)
      if piszero(r2) then N = padd(pmul(N, q2), tn); D = td
      else N = padd(pmul(N, td), pmul(tn, D)); D = pmul(D, td) end
    end
  end
  return N, D, minpow
end
