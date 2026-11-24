-- ===================== TOOLS (part A) =====================
local function sampleIn(lo, hi, dt)
  if lo == -INF and hi == INF then return 0 end
  if lo == -INF then return hi - 1 end
  if hi == INF then return dt and lo or lo + 1 end
  return dt and lo or (lo + hi) / 2
end

local TOOLS = {}
local function addTool(t) TOOLS[#TOOLS + 1] = t end
local function V(dt) return dt and "n" or "t" end
local function sigLines(out, sig, var, name)
  name = name or "x"
  local o, c = sig.dt and "[" or "(", sig.dt and "]" or ")"
  out:r(name .. o .. var .. c .. " = " .. sig:compact(var))
  local pw = sig:piecewise(var)
  if #pw > 1 then
    out:t("Piecewise:")
    for _, l in ipairs(pw) do out:t("   " .. l.expr .. ",   " .. l.range) end
  end
  if #sig.imps > 0 and #pw > 1 then out:t("   plus impulses shown above") end
  local lst = sig:dtList()
  if lst and #lst <= 16 and #lst > 1 then
    local parts = {}
    for _, e in ipairs(lst) do parts[#parts + 1] = (e.n == 0) and ("[" .. fmtC(e.v) .. "]") or fmtC(e.v) end
    out:t("Values n=" .. lst[1].n .. ".." .. lst[#lst].n .. " ([ ] = n=0): {" .. table.concat(parts, ", ") .. "}")
  end
end

-- ---------- periodicity ----------
local function periodInfo(sig, out)
  local dt = sig.dt
  local var = V(dt)
  if #sig.bp > 0 or #sig.imps > 0 then
    out:t("Not periodic: signal is not the same for all " .. var .. " (it has a start/stop or pieces).")
    return nil
  end
  local S = sig.seg[1]
  if #S == 0 then out:t("x = 0 (trivially periodic)"); return nil end
  local ws = {}
  for _, m in ipairs(S) do
    if m.k > 0 then out:t("Not periodic: contains " .. var .. "^" .. m.k .. " (grows)."); return nil end
    if dt then
      if abs(cabs(m.s) - 1) > 1e-9 then out:t("Not periodic: |r| = " .. fmtR(cabs(m.s)) .. " ≠ 1 (grows/decays)."); return nil end
      local w = carg(m.s)
      if abs(w) > 1e-12 then ws[#ws + 1] = w end
    else
      if abs(m.s.re) > 1e-12 then out:t("Not periodic: has e^(" .. fmtR(m.s.re) .. var .. ") growth/decay."); return nil end
      if abs(m.s.im) > 1e-12 then ws[#ws + 1] = m.s.im end
    end
  end
  if #ws == 0 then out:t("Constant signal: periodic with any period (no fundamental period)."); return nil end
  out:t("Rewrite as sum of sinusoids / exponentials (Euler, trig identities): x = " .. sig:compact(var))
  out:t("Each e^(jω" .. var .. ") (or cos/sin(ω" .. var .. ")) term has its own period; x is periodic only if all periods share a common multiple:")
  local uniq = {}
  for _, w in ipairs(ws) do
    local aw = abs(w)
    local f = false
    for _, u in ipairs(uniq) do if abs(u - aw) < 1e-9 then f = true end end
    if not f then uniq[#uniq + 1] = aw end
  end
  if dt then
    local N = 1
    for _, w in ipairs(uniq) do
      local p, q = ratapprox(w / (2 * pi), 10000, 1e-9)
      if not p then
        out:t("  Ω = " .. fmtR(w) .. ": Ω/2π = " .. fmtR(w / (2 * pi)) .. " is irrational → NOT periodic")
        out:r("Not periodic")
        return nil
      end
      local Ni = q
      out:t("  Ω = " .. fmtR(w) .. ": 2π/Ω = " .. fmtR(2 * pi / w) .. " = N/k with N = " .. Ni .. " (reduced) → period " .. Ni)
      N = lcmi(N, Ni)
    end
    out:r("Periodic, fundamental period N = " .. N .. (#uniq > 1 and "  (LCM of component periods)" or ""))
    return N
  else
    local base = uniq[1]
    local qs = {}
    for _, w in ipairs(uniq) do
      out:t("  ω = " .. fmtR(w) .. " → period 2π/ω = " .. fmtR(2 * pi / w))
      local p, q = ratapprox(w / base, 1000, 1e-9)
      if not p then
        out:t("  ratio " .. fmtR(w) .. "/" .. fmtR(base) .. " is irrational → NOT periodic")
        out:r("Not periodic")
        return nil
      end
      qs[#qs + 1] = {p, q}
    end
    -- w_i = base * p_i/q_i ; w0 = base * gcd(p_i/q_i) = base * gcd(p)/lcm(q)
    local L = 1
    for _, r in ipairs(qs) do L = lcmi(L, r[2]) end
    local G = 0
    for _, r in ipairs(qs) do G = gcdi(G, r[1] * (L / r[2])) end
    local w0 = base * G / L
    local T = 2 * pi / w0
    out:t("Fundamental frequency ω0 = GCD of frequencies = " .. fmtR(w0))
    out:r("Periodic, fundamental period T = 2π/ω0 = " .. fmtR(T))
    return T
  end
end

-- ---------- even/odd ----------
local function evenOdd(sig, out)
  local var = V(sig.dt)
  local r = sig:compose(-1, 0)
  local o, c = sig.dt and "[" or "(", sig.dt and "]" or ")"
  out:t("x" .. o .. "-" .. var .. c .. " = " .. r:compact(var))
  local even = ssub(sig, r):isZero()
  local odd = sadd(sig, r):isZero()
  if even and odd then out:r("Both even and odd (x = 0)")
  elseif even then out:r("EVEN: x" .. o .. "-" .. var .. c .. " = x" .. o .. var .. c)
  elseif odd then out:r("ODD: x" .. o .. "-" .. var .. c .. " = -x" .. o .. var .. c)
  else
    out:r("Neither even nor odd")
    local xe, xo = sadd(sig, r):scale(0.5), ssub(sig, r):scale(0.5)
    out:t("Even part  xe = [x" .. o .. var .. c .. " + x" .. o .. "-" .. var .. c .. "]/2 = " .. xe:compact(var))
    for _, l in ipairs(xe:piecewise(var)) do if #xe.bp > 0 then out:t("     " .. l.expr .. ",  " .. l.range) end end
    out:t("Odd part   xo = [x" .. o .. var .. c .. " - x" .. o .. "-" .. var .. c .. "]/2 = " .. xo:compact(var))
    for _, l in ipairs(xo:piecewise(var)) do if #xo.bp > 0 then out:t("     " .. l.expr .. ",  " .. l.range) end end
    out:t("Check: xe + xo = x,  xe(-" .. var .. ") = xe,  xo(-" .. var .. ") = -xo")
    out:plot({kind = "sig", sig = xe, label = "even part"})
    out:plot({kind = "sig", sig = xo, label = "odd part"})
  end
  return r
end

local function conjSig(sig)
  local Cj = sig:copy()
  for _, S in ipairs(Cj.seg) do for _, m in ipairs(S) do m.c = cconj(m.c); m.s = cconj(m.s) end end
  for _, im in ipairs(Cj.imps) do im.c = cconj(im.c) end
  return Cj:normalize()
end

-- ===== Tool: Signal analysis =====
addTool({
  name = "Signal analyzer",
  help = "Exact form, sketch, even/odd, periodic, energy. Use defs: x(2t+1), x[3n], x(t)u(t)...",
  fields = {{"x(t) or x[n] =", "e^(-2t)u(t)"}},
  run = function(v, ctx, out)
    local node = P(ctx, v[1])
    local dt = isDTexpr(node, ctx)
    local var = V(dt)
    local ok, sig = pcall(toSig, ctx, node, dt)
    if not ok then
      out:e("Exact form unavailable: " .. tostring(sig))
      out:plot({kind = "fn", f = numFn(ctx, node, var), dt = dt, label = v[1]})
      return
    end
    out:h("Signal")
    local nm = v[1]:match("^%s*([%a][%w]*)%s*[%(%[]") 
    sigLines(out, sig, var, (nm and ctx.defs[nm] and #v[1] < 18) and v[1]:gsub("[%(%[].*$", "") or "x")
    out:plot({kind = "sig", sig = sig, label = v[1]})
    -- time transformation explanation
    if node.k == "call" and ctx.defs[node.name] then
      local a, b = linearOf(node.args[1], var, ctx.env)
      if a then
        a, b = a.re, b.re
        local d = ctx.defs[node.name]
        local base = toSig(ctx, P(ctx, node.name .. (dt and "[" or "(") .. var .. (dt and "]" or ")")), dt)
        out:h("Time transformation " .. v[1])
        out:t("Original: " .. node.name .. (dt and "[" or "(") .. var .. (dt and "]" or ")") .. " = " .. base:compact(var))
        if abs(a - 1) < 1e-12 then
          out:t(abs(b) < 1e-12 and "No change." or ("Pure shift: " .. ((b > 0) and ("LEFT by " .. fmtR(b)) or ("RIGHT by " .. fmtR(-b))) .. " (x(t+b) moves the graph left when b>0)"))
        else
          out:t("Method 1 (shift then scale): first s(" .. var .. ") = x(" .. shiftStr(var, -b) .. ") [" .. ((b > 0) and ("left by " .. fmtR(b)) or ((b < 0) and ("right by " .. fmtR(-b)) or "no shift")) .. "], then replace " .. var .. " by " .. fmtR(a) .. var)
          out:t("  scaling by " .. fmtR(a) .. ": " .. ((abs(a) > 1) and "compress (squeeze) by factor " .. fmtR(abs(a)) or "stretch by factor " .. fmtR(1 / abs(a))) .. ((a < 0) and ", and time-reverse (flip)" or ""))
          out:t("Method 2: x(" .. fmtR(a) .. "(" .. shiftStr(var, -b / a) .. ")) → scale first, then shift by " .. fmtR(-b / a))
        end
        if #base.bp > 0 or #base.imps > 0 then
          local mp = {}
          for _, p in ipairs(base.bp) do mp[#mp + 1] = fmtR(p) .. "→" .. fmtR((p - b) / a) end
          for _, im in ipairs(base.imps) do mp[#mp + 1] = "δ@" .. fmtR(im.t0) .. "→" .. fmtR((im.t0 - b) / a) end
          out:t("Key points move by " .. var .. "_new = (" .. var .. "_old - " .. fmtR(b) .. ")/" .. fmtR(a) .. ":  " .. table.concat(mp, ", "))
          if not dt and abs(a) ~= 1 and #base.imps > 0 then out:t("Impulse areas get multiplied by 1/|a| = " .. fmtR(1 / abs(a)) .. "  (δ(at) = δ(t)/|a|)") end
          if dt and abs(a) > 1 then out:t("DT compression x[" .. fmtR(a) .. "n] keeps only every " .. fmtR(abs(a)) .. "-th sample (others are lost)") end
        end
      end
    end
    out:h("Even / odd")
    out:t("Substitute " .. var .. " → -" .. var .. ": x" .. (dt and "[-n]" or "(-t)") .. " = " .. tostr(subst(node, var, N("neg", {a = N("var", {name = var})}))))
    evenOdd(sig, out)
    out:h("Periodicity")
    local T = periodInfo(sig, out)
    out:h("Causality of the signal")
    local first = sig.bp[1]
    local leftZero = (#sig.seg[1] == 0)
    local negZero = true
    for i, S in ipairs(sig.seg) do
      local lo = sig.bp[i - 1] or -INF
      if #S > 0 and lo < 0 then negZero = false end
    end
    for _, im in ipairs(sig.imps) do if im.t0 < 0 then negZero = false end end
    out:t(negZero and ("x = 0 for " .. var .. " < 0 → causal (right-sided) signal") or ("x ≠ 0 for some " .. var .. " < 0 → non-causal signal"))
    out:h("Energy / power")
    if #sig.imps > 0 then out:t("Contains impulses → infinite energy")
    else
      local ok2, E = pcall(function() return smul(sig, conjSig(sig)):integrate(-INF, INF) end)
      if ok2 then
        out:t("E = " .. (dt and "Σ|x[n]|²" or "∫|x(t)|²dt") .. " = " .. fmtC(E))
        out:r("Energy signal (finite energy), P = 0")
      else
        out:t("E = ∞")
        if T then
          local per = dt and T or T
          local ok3, Pv = pcall(function()
            local x2 = smul(sig, conjSig(sig))
            if dt then return x2:integrate(0, per - 1) / per end
            return x2:integrate(0, per) / per
          end)
          if ok3 then out:r("Power signal: P = (1/T)∫|x|² over a period = " .. fmtC(Pv)) end
        else
          local ok4, Pv = pcall(function()
            local x2 = smul(sig, conjSig(sig))
            local Lv = 2000
            if dt then return x2:integrate(-Lv, Lv) / (2 * Lv + 1) end
            return x2:integrate(-Lv, Lv) / (2 * Lv)
          end)
          if ok4 and cabs(Pv) < 1e6 then out:t("P ≈ " .. fmtC(Pv) .. " (average over a long window)") end
        end
      end
    end
  end,
})

-- ===== Tool: running integral / sum =====
addTool({
  name = "Running integral / sum",
  help = "y(t)=∫[-∞,t] x(τ)dτ  or  y[n]=Σ[k=-∞,n] x[k]",
  fields = {{"x(t) or x[n] =", "sin(2t)u(t)"}},
  run = function(v, ctx, out)
    local node = P(ctx, v[1])
    local dt = isDTexpr(node, ctx)
    local var = V(dt)
    local sig = toSig(ctx, node, dt)
    local dv = dt and "k" or "τ"
    out:h("Input")
    sigLines(out, sig, var)
    out:h("Steps")
    out:t(dt and "y[n] = Σ[k=-∞→n] x[k]" or "y(t) = ∫[-∞→t] x(τ) dτ")
    if #sig.imps > 0 then
      for _, im in ipairs(sig.imps) do
        if im.m == 0 then out:t("∫ of " .. impStr(im, "τ") .. " → " .. cmul(im.c, uStr("t", im.t0, false))) end
      end
    end
    local Y = sig:running(false)
    -- region by region
    local onlyImp = true
    for _, S in ipairs(sig.seg) do if #S > 0 then onlyImp = false end end
    for i, S in ipairs(onlyImp and {} or sig.seg) do
      local lo, hi = sig.bp[i - 1] or -INF, sig.bp[i] or INF
      local yexpr = segStr(Y.seg[Y:segIndex(sampleIn(lo, hi, dt))], var, dt)
      local xs = segStr(S, dv, dt)
      if #S == 0 and i == 1 then
        out:t("For " .. rangeStr(lo, hi, var, dt) .. ": x = 0 so y = 0")
      else
        local start = (lo == -INF) and "-∞" or fmtR(lo)
        out:t("For " .. rangeStr(lo, hi, var, dt) .. ":")
        if lo == -INF then
          out:t(dt and ("  y[n] = Σ[k=-∞→n] " .. xs) or ("  y(t) = ∫[-∞→t] " .. xs .. " dτ"))
        else
          out:t(dt and ("  y[n] = y[" .. fmtR(lo - 1) .. "] + Σ[k=" .. start .. "→n] " .. xs) or ("  y(t) = y(" .. start .. ") + ∫[" .. start .. "→t] " .. (xs == "0" and "0" or xs) .. " dτ"))
        end
        if #S == 0 and lo > -INF then
          out:t("    x = 0 here, so y keeps its value " .. fmtC(dt and Y:eval(lo - 1) or select(1, Y:limits(lo))) .. ((not dt and #sig.imps > 0) and " (plus impulse areas)" or ""))
        end
        if #S > 0 then
          local G = dt and antiDiffS(S) or antiS(S)
          local Gs = segStr(G, dv, dt)
          if dt then
            out:t("    antidifference: G[k] = " .. Gs .. "  (G[k] - G[k-1] = x[k])")
            out:t("    Σ[k=" .. start .. "→n] x[k] = G[n] - G[" .. ((lo == -INF) and "-∞" or fmtR(lo - 1)) .. "]")
          else
            out:t("    antiderivative: ∫ " .. xs .. " dτ = " .. Gs)
            out:t("    ∫[" .. start .. "→t] = [" .. Gs .. "] evaluated " .. start .. " → t" .. ((lo == -INF) and "  (lower limit → 0)" or (" = G(t) - G(" .. start .. ")")))
          end
          if lo > -INF then
            local yprev = dt and Y:eval(lo - 1) or select(1, Y:limits(lo))
            out:t("    y" .. (dt and ("[" .. fmtR(lo - 1) .. "]") or ("(" .. start .. ")")) .. " = " .. fmtC(yprev) .. ((not dt) and ("   (value carried from previous interval" .. (#sig.imps > 0 and ", plus impulse areas at " .. start or "") .. ")") or ""))
          end
        end
        out:t("       = " .. yexpr)
        if dt and lo > -INF and hi < INF and hi - lo <= 8 then
          local vals = {}
          for n = lo, hi - 1 do vals[#vals + 1] = "y[" .. n .. "]=" .. fmtC(Y:eval(n)) end
          out:t("       " .. table.concat(vals, ", "))
        end
      end
    end
    out:h("Result")
    sigLines(out, Y, var, "y")
    out:plot({kind = "sig", sig = Y, label = "y"})
    out:plot({kind = "sig", sig = sig, label = "x"})
  end,
})

-- ===== Tool: definite integral / sum =====
addTool({
  name = "Integral / sum evaluator",
  help = "∫[a,b] f(t)dt with δ sifting & u() limits. Use inf for ∞. For DT f[n] sums n=a..b",
  fields = {{"f(t) or f[n] =", "u(1-t)cos(2pi t)"}, {"from", "0"}, {"to", "inf"}},
  run = function(v, ctx, out)
    local node = P(ctx, v[1])
    local dt = isDTexpr(node, ctx)
    local var = V(dt)
    local a = realof(evalN(P(ctx, v[2]), ctx.env), "limit")
    local b = realof(evalN(P(ctx, v[3]), ctx.env), "limit")
    local ok, sig = pcall(toSig, ctx, node, dt)
    if not ok then
      out:e("Not closed-form (" .. tostring(sig) .. "). Numeric:")
      local f = numFn(ctx, node, var)
      local lo, hi = math.max(a, -200), math.min(b, 200)
      local s
      if dt then s = ZERO; for n = math.ceil(lo), floor(hi) do s = s + (f(n) or ZERO) end
      else s = simpson(function(x) return f(x) or ZERO end, lo, hi, 4000) end
      out:r("≈ " .. fmtC(s))
      return
    end
    out:h("Integrand")
    sigLines(out, sig, var, "f")
    out:h("Steps")
    local lim = (dt and "Σ[n=" or "∫[") .. fmtR(a) .. "→" .. fmtR(b) .. "]"
    for _, im in ipairs(sig.imps) do
      local inside = im.t0 >= math.min(a, b) and im.t0 <= math.max(a, b)
      if im.m == 0 then
        out:t("Sifting: " .. impStr(im, "t") .. " picks value at t = " .. fmtR(im.t0) .. (inside and " (inside limits) → " .. fmtC(im.c) or " (outside limits) → 0"))
      else
        out:t("Doublet " .. impStr(im, "t") .. ": ∫f(t)δ'(t-t0)dt = -f'(t0)" .. (inside and "" or " (outside → 0)"))
      end
    end
    local total = ZERO
    local partsS = {}
    for _, im in ipairs(sig.imps) do
      if im.m == 0 and im.t0 >= math.min(a, b) and im.t0 <= math.max(a, b) then partsS[#partsS + 1] = fmtC(im.c) end
    end
    for i, S in ipairs(sig.seg) do
      local lo, hi = sig.bp[i - 1] or -INF, sig.bp[i] or INF
      local L, H = math.max(lo, math.min(a, b)), math.min(hi, math.max(a, b))
      if dt then H = math.min(hi - 1, math.max(a, b)) end
      if #S > 0 and H >= L and (H > L or dt) then
        local piece = fromPieces(dt, {{lo = L, hi = dt and (H + 1) or H, S = S}})
        local val = piece:integrate(L, H)
        out:t(lim:gsub("%[.-%]", "[" .. fmtR(L) .. "→" .. fmtR(H) .. "]") .. " " .. segStr(S, var, dt) .. (dt and "" or " dt"))
        local G = dt and antiDiffS(S) or antiS(S)
        local function Gv(x)
          if x == INF then return limitS(G, 1, dt) end
          if x == -INF then return limitS(G, -1, dt) end
          return segEval(G, x, dt)
        end
        if dt then
          out:t("   antidifference G[n] = " .. segStr(G, var, dt) .. ";  sum = G[" .. fmtR(H) .. "] - G[" .. fmtR(L - 1) .. "]")
          local a1, a2 = Gv(H), (L == -INF) and Gv(-INF) or segEval(G, L - 1, true)
          if a1 and a2 then out:t("   = " .. fmtC(a1) .. " - (" .. fmtC(a2) .. ") = " .. fmtC(val)) else out:t("   = " .. fmtC(val)) end
        else
          out:t("   antiderivative G(t) = " .. segStr(G, var, dt))
          local a1, a2 = Gv(H), Gv(L)
          local ls = (L == -INF) and "G(-∞)" or ("G(" .. fmtR(L) .. ")")
          local hs = (H == INF) and "G(∞)" or ("G(" .. fmtR(H) .. ")")
          if a1 and a2 then out:t("   = " .. hs .. " - " .. ls .. " = " .. fmtC(a1) .. " - (" .. fmtC(a2) .. ") = " .. fmtC(val)) else out:t("   = " .. fmtC(val)) end
        end
        total = total + val
        partsS[#partsS + 1] = fmtC(val)
      end
    end
    local res = sig:integrate(a, b)
    out:h("Result")
    if #partsS > 1 then out:t("Add the parts: " .. joinsum(partsS)) end
    out:r(lim .. " f = " .. ((#partsS > 1 and not fmtC(res):find("[π√e/]")) and joinsum(partsS) .. " = " or "") .. fmtC(res))
    out:plot({kind = "sig", sig = sig, label = "f"})
  end,
})

-- ===== Tool: convolution =====
-- h(var - dv) expanded, per monomial list
local function shiftedStr(S, var, dv, dt)
  local terms = {}
  for _, m in ipairs(S) do
    for i = 0, m.k do
      local cf = m.c * binom(m.k, i) * (((i % 2) == 1) and -1 or 1)
      local tp = joinf(powStr(var, m.k - i), dt and baseStr(m.s, var) or expStr(m.s, var))
      local dp = joinf(powStr(dv, i), dt and baseStr(ONE / m.s, dv) or expStr(-m.s, dv))
      local f = joinf(tp, dp)
      terms[#terms + 1] = cmul(cf, f == "" and "1" or f)
    end
  end
  return joinsum(terms)
end
local function convSteps(X, H, out)
  local dt = X.dt
  local var = V(dt)
  local dv = dt and "k" or "τ"
  local o, c = dt and "[" or "(", dt and "]" or ")"
  out:h("Step 1: Convolution formula")
  out:t(dt and "y[n] = x[n]*h[n] = Σ[k=-∞→∞] x[k]·h[n-k]" or "y(t) = x(t)*h(t) = ∫[-∞→∞] x(τ)·h(t-τ) dτ")
  out:t("Rename the variable in x: x" .. o .. dv .. c .. " = " .. X:compact(dv))
  out:t("h" .. o .. var .. c .. " = " .. H:compact(var))
  out:h("Step 2: Flip h, then shift by " .. var)
  local Hf = H:compose(-1, 0)
  out:t("Flip: h" .. o .. "-" .. dv .. c .. " = " .. Hf:compact(dv))
  local hs = H:compact(var)
  out:t("Shift: h" .. o .. var .. "-" .. dv .. c .. " = " .. substVar(hs, var, "(" .. var .. " - " .. dv .. ")"))
  for _, pb in ipairs(H:pieces()) do
    local lo = (pb.hi == INF) and "-∞" or shiftStr(var, dt and (pb.hi - 1) or pb.hi)
    local hi = (pb.lo == -INF) and "∞" or shiftStr(var, pb.lo)
    out:t("   = " .. shiftedStr(pb.S, var, dv, dt) .. "   for " .. lo .. (dt and " ≤ " or " < ") .. dv .. (dt and " ≤ " or " < ") .. hi)
  end
  local ed = {}
  for _, b in ipairs(H.bp) do ed[#ed + 1] = dv .. " = " .. shiftStr(var, b) end
  if #ed > 0 and not dt then out:t("Edges of h" .. o .. var .. "-" .. dv .. c .. " sit at " .. table.concat(ed, ", ") .. " (they move right as " .. var .. " grows)") end
  local ex = {}
  for _, b in ipairs(X.bp) do ex[#ex + 1] = fmtR(b) end
  if #ex > 0 then out:t("Edges of x" .. o .. dv .. c .. " are fixed at " .. dv .. " = " .. table.concat(ex, ", ")) end
  out:t("Breakpoints in " .. var .. " happen when an edge of h meets an edge of x.")
  for _, im in ipairs(X.imps) do out:t("Impulse rule: " .. impStr(im, var) .. " * h" .. o .. var .. c .. " = " .. cmul(im.c, "h" .. o .. shiftStr(var, im.t0) .. c) .. "  (x*δ(t-t0) = x(t-t0))") end
  for _, im in ipairs(H.imps) do out:t("Impulse rule: x" .. o .. var .. c .. " * " .. impStr(im, var) .. " = " .. cmul(im.c, "x" .. o .. shiftStr(var, im.t0) .. c)) end
  local Y = sconv(X, H)
  out:h("Step 3: Slide through each interval")
  local PX, PH = X:pieces(), H:pieces()
  local sm = {}
  for i = 1, #Y.seg do
    local lo, hi = Y.bp[i - 1] or -INF, Y.bp[i] or INF
    local tm = sampleIn(lo, hi, dt)
    if dt and lo == -INF and hi ~= INF then tm = hi - 1 end
    out:t("Case " .. rangeStr(lo, hi, var, dt) .. ":")
    local any = false
    for _, pa in ipairs(PX) do
      for _, pb in ipairs(PH) do
        local L, U, Ls, Us
        if dt then
          L = math.max(pa.lo, tm - pb.hi + 1); U = math.min(pa.hi - 1, tm - pb.lo)
          Ls = (pa.lo >= tm - pb.hi + 1) and ((pa.lo == -INF) and "-∞" or fmtR(pa.lo)) or shiftStr("n", pb.hi - 1)
          Us = (pa.hi - 1 <= tm - pb.lo) and ((pa.hi == INF) and "∞" or fmtR(pa.hi - 1)) or shiftStr("n", pb.lo)
          if pb.hi == INF then Ls = (pa.lo == -INF) and "-∞" or fmtR(pa.lo) end
          if pb.lo == -INF then Us = (pa.hi == INF) and "∞" or fmtR(pa.hi - 1) end
        else
          L = math.max(pa.lo, tm - pb.hi); U = math.min(pa.hi, tm - pb.lo)
          Ls = (pa.lo >= tm - pb.hi) and ((pa.lo == -INF) and "-∞" or fmtR(pa.lo)) or shiftStr("t", pb.hi)
          Us = (pa.hi <= tm - pb.lo) and ((pa.hi == INF) and "∞" or fmtR(pa.hi)) or shiftStr("t", pb.lo)
        end
        if U > L or (dt and U >= L) then
          any = true
          local fx = segStr(pa.S, dv, dt)
          local fh = shiftedStr(pb.S, var, dv, dt)
          if topSum(fh) or fh:sub(1, 1) == "-" then fh = "(" .. fh .. ")" end
          local prod
          if fx == "1" then prod = fh elseif fx == "-1" then prod = "-" .. fh
          else prod = (topSum(fx) and ("(" .. fx .. ")") or fx) .. "·" .. fh end
          local SIG = dt and ("Σ[k=" .. Ls .. "→" .. Us .. "] ") or ("∫[" .. Ls .. "→" .. Us .. "] ")
          out:t("  Overlap: " .. Ls .. (dt and " ≤ k ≤ " or " < τ < ") .. Us)
          out:t("  " .. SIG .. prod .. (dt and "" or " dτ"))
          local nterms = #pa.S * #pb.S
          for _, ma in ipairs(pa.S) do
            for _, mb in ipairs(pb.S) do
              local D = dt and (ma.s / mb.s) or (ma.s - mb.s)
              local pre = M(ma.c * mb.c, 0, mb.s)
              local preS = segStr({pre}, var, dt)
              for ii = 0, mb.k do
                local cf = binom(mb.k, ii) * (((ii % 2) == 1) and -1 or 1)
                local p = ma.k + ii
                local outside = segStr({M(ma.c * mb.c * cf, mb.k - ii, mb.s)}, var, dt)
                local inside = segStr({M(1, p, D)}, dv, dt)
                local o1 = (outside == "1") and "" or ((topSum(outside) and ("(" .. outside .. ")") or outside) .. "·")
                if o1 == "-1·" then o1 = "-" end
                out:t("    " .. ((nterms > 1 or mb.k > 0) and "term: " or "= ") .. o1 .. SIG .. inside .. (dt and "" or " dτ"))
                local Lk, Uk
                if dt then
                  Lk = (pa.lo >= tm - pb.hi + 1 or pb.hi == INF) and {kind = "c", v = pa.lo} or {kind = "t", c0 = pb.hi - 1}
                  Uk = (pa.hi - 1 <= tm - pb.lo or pb.lo == -INF) and {kind = "c", v = pa.hi} or {kind = "t", c0 = pb.lo - 1}
                else
                  Lk = (pa.lo >= tm - pb.hi or pb.hi == INF) and {kind = "c", v = pa.lo} or {kind = "t", c0 = pb.hi}
                  Uk = (pa.hi <= tm - pb.lo or pb.lo == -INF) and {kind = "c", v = pa.hi} or {kind = "t", c0 = pb.lo}
                end
                local tval = nil
                local okv = pcall(function()
                  local up = antiAt(p, D, Uk, dt)
                  local dn = antiAt(p, D, Lk, dt)
                  for _, m in ipairs(dn) do m.c = -m.c; up[#up + 1] = m end
                  tval = segStr(mergeS(mulS(up, {M(ma.c * mb.c * cf, mb.k - ii, mb.s)}, dt), dt), var, dt)
                end)
                if dt then
                  if p == 0 and not near(D, ONE, 1e-12) then
                    local Ds = fmtC(D); if not Ds:find("^[%d%.]+$") then Ds = "(" .. Ds .. ")" end
                    if Us == "∞" then out:t("      = " .. o1 .. Ds .. "^(" .. Ls .. ")/(1 - " .. Ds .. ")   (geometric, |ratio|<1)")
                    elseif Ls == "-∞" then out:t("      = " .. o1 .. Ds .. "^(" .. Us .. "+1)/(" .. Ds .. " - 1)   (geometric, |ratio|>1)")
                    else out:t("      = " .. o1 .. "[" .. Ds .. "^(" .. Ls .. ") - " .. Ds .. "^(" .. Us .. "+1)]/(1 - " .. Ds .. ")   (finite geometric sum)") end
                  elseif p == 0 then
                    out:t("      = " .. o1 .. "(number of terms) = " .. o1 .. "((" .. Us .. ") - (" .. Ls .. ") + 1)")
                  else
                    local q, d = antidiffQ(p, D)
                    local Hs = {}
                    for j = 0, d do Hs[#Hs + 1] = M(q[j], j, D) end
                    out:t("      = " .. o1 .. "[" .. segStr(mergeS(Hs, true), "k", true) .. "] from k=" .. Ls .. " to k=" .. Us .. "+1  (antidifference)")
                  end
                else
                  local G = antiS({M(1, p, D)})
                  out:t("      = " .. o1 .. "[" .. segStr(G, "τ", false) .. "] evaluated " .. Ls .. " → " .. Us)
                end
                if okv and tval then out:t("      = " .. tval) end
              end
            end
          end
        end
      end
    end
    if not any then out:t("  No overlap between x and the shifted h  ⇒  y = 0") end
    for _, im in ipairs(X.imps) do out:t("  + impulse term " .. cmul(im.c, "h" .. o .. shiftStr(var, im.t0) .. c)) end
    for _, im in ipairs(H.imps) do out:t("  + impulse term " .. cmul(im.c, "x" .. o .. shiftStr(var, im.t0) .. c)) end
    local ys = segStr(Y.seg[i], var, dt)
    out:r("  y = " .. ys .. "   for " .. rangeStr(lo, hi, var, dt))
    if lo ~= -INF and hi ~= INF and not dt then
      out:t("  check: y(" .. fmtR(lo) .. "+) = " .. fmtC(segEval(Y.seg[i], lo, false)) .. ",  y(" .. fmtR(hi) .. "-) = " .. fmtC(segEval(Y.seg[i], hi, false)))
    end
  end
  out:h("Step 4: Combine the cases")
  sigLines(out, Y, var, "y")
  if not dt then
    local ok, area = pcall(function() return Y:integrate(-INF, INF) end)
    local okx, ax = pcall(function() return X:integrate(-INF, INF) end)
    local okh, ah = pcall(function() return H:integrate(-INF, INF) end)
    if ok and okx and okh then out:t("Check (area rule): area(y) = " .. fmtC(area) .. " = area(x)·area(h) = " .. fmtC(ax) .. "·" .. fmtC(ah)) end
    local function span(S) return (S.bp[1] or -INF), (S.bp[#S.bp] or INF) end
  end
  out:plot({kind = "sig", sig = Y, label = "y"})
  out:plot({kind = "sig", sig = X, label = "x"})
  out:plot({kind = "sig", sig = H, label = "h"})
  return Y
end
addTool({
  name = "Convolution",
  help = "y = x * h with flip-and-slide steps. Works for CT (t) and DT (n).",
  fields = {{"x =", "e^(-t)u(t)"}, {"h =", "e^(-2t)u(t)"}},
  run = function(v, ctx, out)
    local nx, nh = P(ctx, v[1]), P(ctx, v[2])
    local dt = isDTexpr(nx, ctx) or isDTexpr(nh, ctx)
    local X, H = toSig(ctx, nx, dt), toSig(ctx, nh, dt)
    convSteps(X, H, out)
  end,
})
