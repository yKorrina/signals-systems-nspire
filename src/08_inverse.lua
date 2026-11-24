-- ===================== FREQUENCY-DOMAIN EXPRESSIONS & INVERSE =====================
-- element: {c, d, N, D, P, win={lo,hi}, imp=w0}
-- mode: "F" (ω, v=jω), "L" (s), "D" (Ω, q=e^{-jΩ}), "Z" (z, q=z^-1)
local function fe(c) return {c = tocx(c), d = 0, N = pconst(1), D = pconst(1), P = pconst(1)} end
local function feCopy(e)
  return {c = e.c, d = e.d, N = e.N, D = e.D, P = e.P, win = e.win and {lo = e.win.lo, hi = e.win.hi}, imp = e.imp}
end
local FVAR = {F = "w", L = "s", D = "W", Z = "z"}

-- evaluate element at frequency w (real for F/D) ignoring impulse/window
local function feEval(e, mode, w)
  local v, q
  if mode == "F" then v = cx(0, w) elseif mode == "L" then v = tocx(w) end
  local val
  if mode == "F" or mode == "L" then
    val = e.c * cexp(-v * e.d) * peval(e.N, v) / peval(e.D, v)
  else
    q = (mode == "D") and cexp(cx(0, -w)) or (ONE / tocx(w))
    val = e.c * cpow(q, e.d) * peval(e.N, q) / peval(e.D, q) * peval(e.P, cx(w))
  end
  return val
end

local function feMul(a, b, mode)
  local r = {c = a.c * b.c, d = a.d + b.d, N = pmul(a.N, b.N), D = pmul(a.D, b.D), P = pmul(a.P, b.P)}
  if a.win or b.win then
    local lo, hi = -INF, INF
    if a.win then lo = math.max(lo, a.win.lo); hi = math.min(hi, a.win.hi) end
    if b.win then lo = math.max(lo, b.win.lo); hi = math.min(hi, b.win.hi) end
    if hi <= lo then return nil end
    r.win = {lo = lo, hi = hi}
  end
  if a.imp and b.imp then error("product of impulses", 0) end
  local im = a.imp or b.imp
  if im then
    local other = a.imp and b or a
    local w0 = im
    local o2 = feCopy(other); o2.imp = nil
    local v = feEval(o2, mode, w0)
    if other.win and (w0 < other.win.lo or w0 > other.win.hi) then v = ZERO end
    local base = a.imp and a or b
    r = {c = base.c * v, d = 0, N = pconst(1), D = pconst(1), P = pconst(1), imp = w0}
  end
  return r
end
local function listMul(A, B, mode)
  local out = {}
  for _, a in ipairs(A) do for _, b in ipairs(B) do
    local r = feMul(a, b, mode)
    if r and cabs(r.c) > 1e-14 then out[#out + 1] = r end
  end end
  return out
end
local function listScale(A, c) local out = {}; for _, a in ipairs(A) do local e = feCopy(a); e.c = e.c * c; out[#out + 1] = e end; return out end
-- shift integer delay into polynomial (DT)
local function absorbQ(e)
  if (e.d ~= 0) and isint(e.d) then
    local k = round(e.d)
    if k > 0 then
      local sh = {n = k}; for i = 0, k do sh[i] = ZERO end; sh[k] = ONE
      e.N = pmul(e.N, sh)
    else
      local sh = {n = -k}; for i = 0, -k do sh[i] = ZERO end; sh[-k] = ONE
      e.D = pmul(e.D, sh)
    end
    e.d = 0
  end
  return e
end
-- combine a list into a single element (for denominators)
local function listToElem(A, mode)
  if #A == 0 then return fe(0) end
  local disc = (mode == "D" or mode == "Z")
  local res = nil
  for _, a in ipairs(A) do
    if a.win or a.imp then error("cannot divide by windowed/impulsive expression", 0) end
    local e = feCopy(a)
    if disc then absorbQ(e) end
    e.N = pscale(e.N, e.c); e.c = ONE
    if not res then res = e
    else
      if abs(res.d - e.d) > 1e-12 then error("cannot combine different delays in a denominator", 0) end
      if pdeg(res.P) > 0 or pdeg(e.P) > 0 then error("polynomial in Ω in denominator not supported", 0) end
      local q, r = pdivmod(res.D, e.D)
      if piszero(r) and pdeg(res.D) >= pdeg(e.D) then res.N = padd(res.N, pmul(e.N, q))
      else
        local q2, r2 = pdivmod(e.D, res.D)
        if piszero(r2) then res.N = padd(pmul(res.N, q2), e.N); res.D = e.D
        else
          res.N = padd(pmul(res.N, e.D), pmul(e.N, res.D))
          res.D = pmul(res.D, e.D)
        end
      end
    end
  end
  return res
end

local function feList(node, mode, env)
  local fv = FVAR[mode]
  local k = node.k
  local function rec(n) return feList(n, mode, env) end
  if k == "num" then return {fe(node.v)} end
  if k == "var" then
    local nm = node.name
    if nm == fv or (mode == "F" and nm == "omega") or (mode == "D" and nm == "Omega") then
      local e = fe(1)
      if mode == "F" then e.N = pnew({0, cx(0, -1)})
      elseif mode == "L" then e.N = pnew({0, 1})
      elseif mode == "D" then e.P = pnew({0, 1})
      else e.d = -1 end
      return {e}
    end
    return {fe(evalN(node, env))}
  end
  if k == "add" then local a = rec(node.a); for _, x in ipairs(rec(node.b)) do a[#a + 1] = x end; return a end
  if k == "sub" then local a = rec(node.a); for _, x in ipairs(listScale(rec(node.b), -1)) do a[#a + 1] = x end; return a end
  if k == "neg" then return listScale(rec(node.a), -1) end
  if k == "mul" then return listMul(rec(node.a), rec(node.b), mode) end
  if k == "div" then
    local B = listToElem(rec(node.b), mode)
    local inv = {c = ONE / B.c, d = -B.d, N = B.D, D = B.N, P = pconst(1)}
    if piszero(B.N) then error("division by zero", 0) end
    return listMul(rec(node.a), {inv}, mode)
  end
  if k == "pow" then
    local ce = isconst(node.b, fv, env)
    local cb = isconst(node.a, fv, env)
    if ce and cb then return {fe(cpow(cb, ce))} end
    if ce then
      local n = realof(ce, "exponent")
      if not isint(n) then
        if mode == "Z" and node.a.k == "var" then local e = fe(1); e.d = -n; return {e} end
        error("non-integer power", 0)
      end
      n = round(n)
      local B = rec(node.a)
      if n < 0 then
        local E = listToElem(B, mode)
        B = {{c = ONE / E.c, d = -E.d, N = E.D, D = E.N, P = pconst(1)}}
        n = -n
      end
      local r = {fe(1)}
      for _ = 1, n do r = listMul(r, B, mode) end
      return r
    end
    if cb then
      -- c^(linear) = e^(ln c * linear)
      return rec(N("call", {name = "exp", args = {N("mul", {a = NUM(clog(cb)), b = node.b})}}))
    end
    error("unsupported power", 0)
  end
  if k == "call" then
    local nm, A = node.name, node.args
    if nm == "exp" or nm == "cos" or nm == "sin" then
      local c0 = isconst(A[1], fv, env)
      if c0 then return {fe(evalN(node, env))} end
      local a, b = linearOf(A[1], fv, env)
      if not a then error(nm .. "() argument must be linear in " .. (DISP[fv] or fv), 0) end
      local function E(mult)
        local A2, B2 = a * mult, b * mult
        local e = fe(cexp(B2))
        local d
        if mode == "L" then d = -A2 elseif mode == "Z" then error("use z^k for delays", 0) else d = J * A2 end
        if abs(d.im) > 1e-9 then error("exponent must be of form jωd (real delay)", 0) end
        e.d = d.re
        if (mode == "D") then absorbQ(e) end
        return e
      end
      if nm == "exp" then return {E(ONE)} end
      if nm == "cos" then local x, y = E(J), E(-J); x.c = x.c / 2; y.c = y.c / 2; return {x, y} end
      local x, y = E(J), E(-J); x.c = x.c / cx(0, 2); y.c = -y.c / cx(0, 2); return {x, y}
    end
    if nm == "delta" then
      local a, b = linearOf(A[1], fv, env)
      if not a then error("δ() argument must be linear", 0) end
      a, b = realof(a), realof(b)
      local e = fe(1 / abs(a)); e.imp = -b / a
      return {e}
    end
    if nm == "rect" then
      local a, b = linearOf(A[1], fv, env)
      a, b = realof(a), realof(b)
      local x1, x2 = (-0.5 - b) / a, (0.5 - b) / a
      local e = fe(1); e.win = {lo = math.min(x1, x2), hi = math.max(x1, x2)}
      return {e}
    end
    local c0 = isconst(node, fv, env)
    if c0 then return {fe(c0)} end
    error(nm .. "() not supported here", 0)
  end
  if k == "cmp" then
    local cv = mkConv(fv, false, env, {})
    local S = cv:sig(node)
    local out = {}
    for _, pc in ipairs(S:pieces()) do
      local mid = (pc.lo == -INF) and (pc.hi - 1) or ((pc.hi == INF) and (pc.lo + 1) or (pc.lo + pc.hi) / 2)
      if pc.lo == -INF and pc.hi == INF then mid = 0 end
      local e = fe(segEval(pc.S, mid, false)); e.win = {lo = pc.lo, hi = pc.hi}
      out[#out + 1] = e
    end
    return out
  end
  error("unsupported expression", 0)
end

-- textual substitution of a variable letter in a formatted expression
function substVar(s, v, rep)
  local out = s:gsub("%f[%a]" .. v .. "%f[%A]", (rep:gsub("%%", "%%%%")))
  out = out:gsub("%((%b())%)", "%1")
  return out
end
-- C(n+m-1, m-1) as polynomial coefficients in n
local function binomPolyN(m)
  local p = pconst(1)
  for i = 1, m - 1 do p = pmul(p, pnew({i, 1})) end
  return pscale(p, cx(1 / fact(m - 1)))
end

-- windowed synthesis integral: (1/2π)∫_{lo}^{hi} P(w) e^{j w (t - d)} dw  ; returns list {c,W,pw,d}
local function synthTerms(Pw, lo, hi, d, coef)
  local res = {}
  local function add(c, W, pw)
    for _, r in ipairs(res) do if abs(r.W - W) < 1e-12 and r.pw == pw and abs(r.d - d) < 1e-12 then r.c = r.c + c; return end end
    res[#res + 1] = {c = c, W = W, pw = pw, d = d}
  end
  Pw = ptrim(Pw)
  for kk = 0, Pw.n do
    local ck = Pw[kk] * coef / (2 * pi)
    if cabs(ck) > 1e-14 then
      local fk = fact(kk)
      for _, lim in ipairs({{hi, 1}, {lo, -1}}) do
        local W, sg = lim[1], lim[2]
        for i = 0, kk do
          local c = ck * sg * cx(W ^ (kk - i)) * (((i % 2 == 1) and -1 or 1) * fk / fact(kk - i)) / cpow(J, i + 1)
          add(c, W, i + 1)
        end
      end
    end
  end
  return res
end
local function synthStr(terms, var)
  local used, parts = {}, {}
  for i, t in ipairs(terms) do
    if not used[i] and cabs(t.c) > 1e-12 then
      used[i] = true
      local tau = (abs(t.d) < 1e-12) and var or ("(" .. shiftStr(var, t.d) .. ")")
      local den = (t.pw == 1) and tau or (tau .. "^" .. t.pw)
      local A, B = t.c, ZERO
      local W = abs(t.W)
      if abs(t.W) > 1e-12 then
        for j = i + 1, #terms do
          local o = terms[j]
          if not used[j] and o.pw == t.pw and abs(o.W + t.W) < 1e-12 and abs(o.d - t.d) < 1e-12 then
            used[j] = true
            local cp, cm = t.c, o.c
            if t.W < 0 then cp, cm = o.c, t.c end
            A = cp + cm; B = (cp - cm) * J
            break
          end
        end
        if near(B, ZERO, 1e-12) and not near(A, ZERO, 1e-12) and A == t.c then
          -- unpaired
          parts[#parts + 1] = cmul(t.c, "e^(j" .. expoArg(t.W, tau) .. ")/" .. den)
          A = ZERO
        end
        local arg = expoArg(W, tau)
        if cabs(A) > 1e-12 then parts[#parts + 1] = cmul(A, "cos(" .. arg .. ")/" .. den) end
        if cabs(B) > 1e-12 then parts[#parts + 1] = cmul(B, "sin(" .. arg .. ")/" .. den) end
      else
        parts[#parts + 1] = cmul(t.c, "1/" .. den)
      end
    end
  end
  return joinsum(parts)
end

-- inverse transform. returns Sig (may be nil), synth string (or nil), steps (Out)
local function inverseXform(src, mode, env, roc, out)
  local fv = FVAR[mode]
  local dt = (mode == "D" or mode == "Z")
  local tv = dt and "n" or "t"
  local vsym = ({F = "jω", L = "s", D = "e^(-jΩ)", Z = "z^-1"})[mode]
  local node = parse(src)
  local list = feList(node, mode, env)
  local result = newSig(dt)
  local synth = {}
  local plain, wins = {}, {}
  for _, e in ipairs(list) do
    if e.imp then
      if mode == "L" or mode == "Z" then error("impulses not allowed here", 0) end
      local c = e.c / (2 * pi)
      local s = dt and cexp(cx(0, e.imp)) or cx(0, e.imp)
      local w = dt and "Ω" or "ω"
      local ds = (abs(e.imp) < 1e-12) and w or (w .. (e.imp >= 0 and " - " or " + ") .. fmtR(abs(e.imp)))
      local es = (abs(e.imp) < 1e-12) and "1" or ("e^(" .. (e.imp < 0 and "-" or "") .. "j" .. expoArg(abs(e.imp), tv) .. ")")
      out:t(cmul(e.c, "δ(" .. ds .. ")") .. "  <->  " .. cmul(c, es))
      result = sadd(result, sexp(s, dt):scale(c))
    elseif e.win or (dt and (pdeg(e.P) > 0 or not isint(e.d))) or (mode == "F" and pdeg(e.D) == 0 and false) then
      wins[#wins + 1] = e
    else
      plain[#plain + 1] = e
    end
  end
  -- group plain by delay (CT) / all (DT)
  local groups = {}
  for _, e in ipairs(plain) do
    local x = feCopy(e)
    if dt then absorbQ(x) end
    local g
    for _, gg in ipairs(groups) do if abs(gg.d - x.d) < 1e-12 then g = gg end end
    if not g then g = {d = x.d, items = {}}; groups[#groups + 1] = g end
    g.items[#g.items + 1] = x
  end
  for _, g in ipairs(groups) do
    local E = listToElem(g.items, mode)
    local Nn, Dd = pscale(E.N, E.c), E.D
    local off = 0
    if dt then
      -- remove q^k factors from denominator
      Dd = ptrim(Dd)
      while Dd.n > 0 and cabs(Dd[0]) < 1e-14 do
        local r = {n = Dd.n - 1}; for i = 1, Dd.n do r[i - 1] = Dd[i] end; Dd = r; off = off - 1
      end
      Nn = ptrim(Nn)
      while Nn.n > 0 and cabs(Nn[0]) < 1e-14 do
        local r = {n = Nn.n - 1}; for i = 1, Nn.n do r[i - 1] = Nn[i] end; Nn = r; off = off + 1
      end
    end
    local delay = dt and off or g.d
    local hdr = ratStr(Nn, Dd, mode)
    if abs(delay) > 1e-12 then
      out:t("Term with delay factor " .. (dt and (mode == "Z" and ("z^" .. -delay) or ("e^(-j" .. delay .. "Ω)")) or delayStr(delay, mode)) .. " → shift by " .. fmtR(delay) .. " at the end")
    end
    local Q, R = pdivmod(Nn, Dd)
    local part = newSig(dt)
    if pdeg(Dd) > 0 then out:t("Rational part: " .. hdr) end
    if not piszero(Q) then
      if pdeg(Dd) > 0 then out:t("Long division: quotient " .. (dt and pstrQ(Q, mode) or pstr(Q, vsym, true)) .. ", remainder " .. (dt and pstrQ(R, mode) or pstr(R, vsym, true))) end
      for i = 0, Q.n do
        if cabs(Q[i]) > 1e-12 then
          if dt then part = sadd(part, swindow(i, i + 1, true):scale(Q[i]))
          else
            if mode == "L" or mode == "F" then
              local r = newSig(false); r.imps = {{c = Q[i], t0 = 0, m = i}}; part = sadd(part, r)
            end
          end
        end
      end
      out:t("  " .. (dt and pstrQ(Q, mode) or pstr(Q, vsym, true)) .. "  <->  " .. part:compact(tv))
    end
    if not piszero(R) then
      local terms = pfe(R, Dd)
      -- factored denominator & residue formulas
      local mult = {}
      for _, t in ipairs(terms) do
        local key = fmtC(t.p)
        mult[key] = math.max(mult[key] or 0, t.m)
      end
      local groupsR = rootGroups(proots(Dd))
      local fac = {}
      local vname = dt and "q" or (mode == "F" and "jω" or "s")
      for _, g2 in ipairs(groupsR) do
        local f
        if dt then f = "(q - " .. fmtC(g2.p) .. ")" else f = vfac(g2.p, mode) end
        fac[#fac + 1] = f .. ((g2.m > 1) and ("^" .. g2.m) or "")
      end
      local lead = Dd[pdeg(Dd)]
      if not dt then
      out:t("Factor the denominator" .. (dt and " (q = " .. ((mode == "Z") and "z^-1" or "e^(-jΩ)") .. ")" or "") .. ": " .. (near(lead, ONE) and "" or (fmtC(lead) .. "·")) .. table.concat(fac, ""))
      for _, g2 in ipairs(groupsR) do
        local Mx = g2.m
        local function Xred(v)
          -- (v - p)^Mx * R(v)/D(v)
          local num = peval(R, v)
          local den = lead
          for _, g3 in ipairs(groupsR) do if g3 ~= g2 then den = den * cpow(v - g3.p, g3.m) end end
          return num / den
        end
        if Mx == 1 then
          out:t("  residue at " .. vname .. " = " .. fmtC(g2.p) .. ": A = " .. vfac(g2.p, mode) .. "·X|at pole = " .. fmtC(ccln(Xred(g2.p))) .. "   (cover-up)")
        else
          out:t("  pole " .. fmtC(g2.p) .. " has order " .. Mx .. ": A_m = 1/(" .. Mx .. "-m)! · d^(" .. Mx .. "-m)/d" .. vname .. "^(" .. Mx .. "-m) [" .. vfac(g2.p, mode) .. "^" .. Mx .. "·X] at pole")
          out:t("    A_" .. Mx .. " = " .. fmtC(ccln(Xred(g2.p))) .. "  (cover-up for highest power)")
          for _, t in ipairs(terms) do
            if near(t.p, g2.p, 1e-6) and t.m < Mx then out:t("    A_" .. t.m .. " = " .. fmtC(t.A) .. "  (derivative formula)") end
          end
        end
      end
      end
      -- choose sides
      local pfs = {}
      if dt then
        -- convert A/(q - qi)^j -> B/(1 - p q)^j
        for _, t in ipairs(terms) do
          if cabs(t.p) < 1e-12 then error("pole at q=0", 0) end
          local p = ONE / t.p
          local B = t.A * cpow(-p, t.m)
          pfs[#pfs + 1] = {p = p, m = t.m, A = ccln(B)}
        end
        local fstr = {}
        for _, t in ipairs(pfs) do fstr[#fstr + 1] = dtermStr({c = t.A, d = 0, p = t.p, m = t.m, dir = 1}, mode) end
        local polesS = {}
        for _, t in ipairs(pfs) do if t.m == 1 or true then polesS[#polesS + 1] = fmtC(t.p) .. (t.m > 1 and ("(×" .. t.m .. ")") or "") end end
        local fz = {}
        for _, t in ipairs(pfs) do fz[#fz + 1] = "(1 - " .. cmul(t.p, (mode == "Z") and "z^-1" or "e^(-jΩ)") .. ")" .. (t.m > 1 and ("^" .. t.m) or "") end
        out:t("Write denominator as product of (1 - p·" .. ((mode == "Z") and "z^-1" or "e^(-jΩ)") .. ") factors")
        for _, t in ipairs(pfs) do
          out:t("  A for p = " .. fmtC(t.p) .. (t.m > 1 and (" (power " .. t.m .. ")") or "") .. ": multiply by (1 - p·q)" .. (t.m > 1 and ("^" .. t.m) or "") .. ", set q = 1/p → " .. fmtC(t.A))
        end
        out:t("Poles (in z): " .. table.concat(polesS, ", "))
        out:t("Partial fractions: " .. joinsum(fstr))
      else
        pfs = terms
        local fstr, polesS = {}, {}
        for _, t in ipairs(terms) do
          fstr[#fstr + 1] = ftermStr({c = t.A, d = 0, p = t.p, m = t.m}, mode)
          local key = fmtC(t.p)
          local mm = mult[key] or 1
          local dup = false
          for _, q in ipairs(polesS) do if q:sub(1, #key) == key and (q == key or q:sub(#key + 1, #key + 1) == " ") then dup = true end end
          if not dup then polesS[#polesS + 1] = key .. (mm > 1 and (" (order " .. mm .. ")") or "") end
        end
        out:t("Poles: " .. (mode == "F" and "jω = " or "s = ") .. table.concat(polesS, ", "))
        out:t("Partial fractions: " .. joinsum(fstr))
      end
      for _, t in ipairs(pfs) do
        local side
        local mag = dt and cabs(t.p) or t.p.re
        if mode == "F" then side = (mag < -1e-12) and "R" or ((mag > 1e-12) and "L" or "B")
        elseif mode == "D" then side = (mag < 1 - 1e-12) and "R" or ((mag > 1 + 1e-12) and "L" or "B")
        else
          local rl, rh = roc and roc.lo, roc and roc.hi
          if not roc then side = "R"
          elseif rl and mag <= rl + 1e-9 then side = "R"
          elseif rh and mag >= rh - 1e-9 then side = "L"
          else side = "R" end
        end
        local sg
        if dt then
          local bp = binomPolyN(t.m)
          local S = {}
          for i = 0, bp.n do if cabs(bp[i]) > 1e-14 then S[#S + 1] = M(bp[i] * t.A, i, t.p) end end
          if side == "R" then sg = fromPieces(true, {{lo = 0, hi = INF, S = S}})
          elseif side == "L" then
            for _, m in ipairs(S) do m.c = -m.c end
            sg = fromPieces(true, {{lo = -INF, hi = 0, S = S}})
          else error("pole on the unit circle: DTFT needs impulses (not handled)", 0) end
        else
          local S = {M(t.A / fact(t.m - 1), t.m - 1, t.p)}
          if side == "R" then sg = fromPieces(false, {{lo = 0, hi = INF, S = S}})
          elseif side == "L" then S[1].c = -S[1].c; sg = fromPieces(false, {{lo = -INF, hi = 0, S = S}})
          else
            if t.m > 1 then error("repeated pole on jω axis not supported", 0) end
            sg = fromPieces(false, {{lo = 0, hi = INF, S = {M(t.A / 2, 0, t.p)}}, {lo = -INF, hi = 0, S = {M(-t.A / 2, 0, t.p)}}})
          end
        end
        local lhs = dt and dtermStr({c = t.A, d = 0, p = t.p, m = t.m, dir = 1}, mode) or ftermStr({c = t.A, d = 0, p = t.p, m = t.m}, mode)
        local why = ""
        if mode == "F" then why = (side == "R") and "  (Re pole<0 → causal)" or ((side == "L") and "  (Re pole>0 → anticausal)" or "  (pole on axis)") end
        if mode == "D" then why = (side == "R") and "  (|p|<1)" or "  (|p|>1 → left-sided)" end
        out:t("  " .. lhs .. "  <->  " .. sg:compact(tv) .. why)
        part = sadd(part, sg)
      end
    end
    if abs(delay) > 1e-12 then
      local unsh = part:compact(tv)
      part = part:compose(1, -delay)
      out:t("Apply shift " .. tv .. " → " .. shiftStr(tv, delay) .. ": " .. substVar(unsh, tv, "(" .. shiftStr(tv, delay) .. ")"))
    end
    result = sadd(result, part)
  end
  -- windowed parts
  for _, e in ipairs(wins) do
    local lo, hi = -pi, pi
    if e.win then lo, hi = e.win.lo, e.win.hi end
    if not dt then
      if pdeg(e.D) > 0 then error("rational × window not supported", 0) end
      -- N(v) with v = jω -> poly in ω
      local Pw = {n = e.N.n}
      for i = 0, e.N.n do Pw[i] = e.N[i] * cpow(J, i) / e.D[0] end
      out:t("Windowed: x(t) = (1/2π)∫[" .. fmtR(lo) .. "→" .. fmtR(hi) .. "] (" .. pstr(Pw, "ω") .. ")" .. ((abs(e.d) > 1e-12) and ("e^(-j" .. fmtR(e.d) .. "ω)") or "") .. "·e^(jωt) dω")
      local ts = synthTerms(Pw, lo, hi, e.d, e.c)
      for _, t in ipairs(ts) do synth[#synth + 1] = t end
    else
      local x = feCopy(e)
      if pdeg(x.D) > 0 then error("rational × window not supported", 0) end
      for i = 0, x.N.n do
        local cN = x.N[i] / x.D[0]
        if cabs(cN) > 1e-14 then
          local ts = synthTerms(x.P, lo, hi, x.d + i, x.c * cN)
          for _, t in ipairs(ts) do synth[#synth + 1] = t end
        end
      end
      out:t("Windowed: x[n] = (1/2π)∫[" .. fmtR(lo) .. "→" .. fmtR(hi) .. "] X(e^(jΩ))e^(jΩn) dΩ")
    end
  end
  local sstr = nil
  if #synth > 0 then sstr = synthStr(synth, tv) end
  return result:normalize(), sstr, synth
end
local function synthEval(terms, t)
  local v = ZERO
  for _, s in ipairs(terms) do
    local tau = t - s.d
    if abs(tau) < 1e-9 then tau = 1e-9 end
    v = v + s.c * cexp(cx(0, s.W * tau)) / cx(tau ^ s.pw)
  end
  return v
end
