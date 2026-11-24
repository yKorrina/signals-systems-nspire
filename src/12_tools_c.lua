-- ===================== TOOLS (part C) =====================
-- linear coefficient extraction from difference/differential equation side
-- returns {Y = poly, X = poly} in v (CT: jω power, DT: q = e^{-jΩ} power)
local function deSide(node, dt, env)
  local k = node.k
  local function zero() return {Y = pconst(0), X = pconst(0)} end
  local function addp(a, b) return {Y = padd(a.Y, b.Y), X = padd(a.X, b.X)} end
  local function scl(a, c) return {Y = pscale(a.Y, c), X = pscale(a.X, c)} end
  local function mono(which, pw)
    local p = {n = pw}; for i = 0, pw do p[i] = ZERO end; p[pw] = ONE
    local r = zero(); r[which] = p; return r
  end
  if k == "add" then return addp(deSide(node.a, dt, env), deSide(node.b, dt, env)) end
  if k == "sub" then return addp(deSide(node.a, dt, env), scl(deSide(node.b, dt, env), -1)) end
  if k == "neg" then return scl(deSide(node.a, dt, env), -1) end
  if k == "mul" then
    local ca = isconst(node.a, dt and "n" or "t", env)
    if ca and not contains(node.a, function(n) return n.k == "call" or (n.k == "var" and (n.name:match("^[xy]'*$"))) end) then return scl(deSide(node.b, dt, env), ca) end
    local cb = isconst(node.b, dt and "n" or "t", env)
    if cb and not contains(node.b, function(n) return n.k == "call" or (n.k == "var" and (n.name:match("^[xy]'*$"))) end) then return scl(deSide(node.a, dt, env), cb) end
    error("equation must be linear with constant coefficients", 0)
  end
  if k == "div" then
    local cb = isconst(node.b, dt and "n" or "t", env)
    if cb then return scl(deSide(node.a, dt, env), ONE / cb) end
    error("bad division", 0)
  end
  if k == "var" then
    local base = node.name:gsub("'+$", "")
    local pr = #node.name - #base
    if base == "y" then return mono("Y", pr) end
    if base == "x" then return mono("X", pr) end
    local v = evalN(node, env)
    error("unexpected term " .. node.name, 0)
  end
  if k == "call" then
    if node.name == "der" or node.name == "der2" then
      local inner = deSide(node.args[1], dt, env)
      local ord = (node.name == "der") and 1 or 2
      local sh = {n = ord}; for i = 0, ord do sh[i] = ZERO end; sh[ord] = ONE
      return {Y = pmul(inner.Y, sh), X = pmul(inner.X, sh)}
    end
    if node.name == "y" or node.name == "x" then
      local which = (node.name == "y") and "Y" or "X"
      if dt then
        local a, b = linearOf(node.args[1], "n", env)
        if not a or abs(a.re - 1) > 1e-9 then error("use y[n-k] / x[n-k]", 0) end
        local d = -round(b.re)
        if d < 0 then
          -- advance: q^-1 ; handle by returning negative power marker
          error("use delays only (y[n-k]); rewrite advances", 0)
        end
        return mono(which, d)
      else
        return mono(which, 0)
      end
    end
  end
  if k == "num" then
    if isz(node.v) then return zero() end
  end
  error("cannot read equation term: " .. tostr(node), 0)
end

-- real/imag parts of polynomial evaluated on the unit circle / jω axis, as strings in w
local function reImPoly(p, dt)
  local re, im = {}, {}
  p = ptrim(p)
  for k = 0, p.n do
    local c = p[k]
    if cabs(c) > 1e-12 then
      if dt then
        -- c e^{-jkΩ} = c cos(kΩ) - j c sin(kΩ)   (c real assumed)
        if k == 0 then re[#re + 1] = fmtC(c)
        else
          local arg = (k == 1) and "Ω" or (k .. "Ω")
          re[#re + 1] = cmul(c, "cos(" .. arg .. ")")
          im[#im + 1] = cmul(-c, "sin(" .. arg .. ")")
        end
      else
        -- c (jω)^k = c j^k ω^k
        local jk = cpow(J, k) * c
        local wk = (k == 0) and "" or ((k == 1) and "ω" or ("ω^" .. k))
        if abs(jk.re) > 1e-12 then re[#re + 1] = cmul(cx(jk.re), wk == "" and "1" or wk) end
        if abs(jk.im) > 1e-12 then im[#im + 1] = cmul(cx(jk.im), wk == "" and "1" or wk) end
      end
    end
  end
  return joinsum(re), joinsum(im)
end
local function magPhase(out, Hterms, dt)
  out:h("Magnitude and phase")
  local f = function(w) if dt then return dtEval(Hterms, w) end return ftEval(Hterms, cx(0, w)) end
  if not dt then
    local g = combineFT(Hterms)
    if #g == 1 then
      local Nn, Dd = g[1].N, g[1].D
      local nr, ni = reImPoly(Nn, false)
      local dr, di = reImPoly(Dd, false)
      out:t("Put jω in: numerator = (" .. nr .. ") + j(" .. ni .. ")")
      out:t("           denominator = (" .. dr .. ") + j(" .. di .. ")")
      out:t("|H(jω)| = √[(" .. nr .. ")² + (" .. ni .. ")²] / √[(" .. dr .. ")² + (" .. di .. ")²]")
      out:t("∠H(jω) = atan2(" .. ni .. ", " .. nr .. ") - atan2(" .. di .. ", " .. dr .. ")" .. ((abs(g[1].d) > 1e-12) and (" - " .. fmtR(g[1].d) .. "ω") or ""))
    else
      local pure = true
      for _, t in ipairs(Hterms) do if t.m ~= 0 then pure = false end end
      if pure then
        local re, im = {}, {}
        for _, t in ipairs(Hterms) do
          if abs(t.d) < 1e-12 then re[#re + 1] = fmtC(t.c)
          else
            local a = expoArg(t.d, "ω")
            re[#re + 1] = cmul(t.c, "cos(" .. a .. ")"); im[#im + 1] = cmul(-t.c, "sin(" .. a .. ")")
          end
        end
        out:t("Use e^(-jωt0) = cos(ωt0) - j sin(ωt0):")
        out:t("  H = (" .. joinsum(re) .. ") + j(" .. joinsum(im) .. ")")
        out:t("|H(jω)| = √[(" .. joinsum(re) .. ")² + (" .. joinsum(im) .. ")²]")
        out:t("∠H(jω) = atan2(" .. joinsum(im) .. ", " .. joinsum(re) .. ")")
        if #Hterms == 2 and near(Hterms[1].c, Hterms[2].c) then
          local d1, d2 = Hterms[1].d, Hterms[2].d
          local mid = (d1 + d2) / 2
          out:t("Symmetric trick: factor e^(-j" .. expoArg(mid, "ω") .. "): H = " .. cmul(Hterms[1].c * 2, "e^(-j" .. expoArg(mid, "ω") .. ")·cos(" .. expoArg(abs(d2 - d1) / 2, "ω") .. ")"))
          out:t("  ⇒ |H| = " .. cmul(cabs(Hterms[1].c) * 2, "|cos(" .. expoArg(abs(d2 - d1) / 2, "ω") .. ")|") .. ",  ∠H = -" .. expoArg(mid, "ω") .. " (+π where cos<0)")
        end
      else
        out:t("H has several delayed parts: |H(jω)| = |Σ terms|; each e^(-jωt0) adds phase -ω·t0")
        out:t("Tip: factor out e^(-jω·t_mid) to make the rest real (e.g. (1-e^(-j2ω))/(jω) = e^(-jω)·2sin(ω)/ω)")
      end
    end
  else
    local ok, Nn, Dd, mp = pcall(combineDT, Hterms)
    if ok then
      local nr, ni = reImPoly(Nn, true)
      local dr, di = reImPoly(Dd, true)
      out:t("Use e^(-jkΩ) = cos(kΩ) - j sin(kΩ):")
      out:t("  numerator = (" .. nr .. ") + j(" .. (ni == "0" and "0" or ni) .. ")")
      out:t("  denominator = (" .. dr .. ") + j(" .. (di == "0" and "0" or di) .. ")")
      out:t("|H| = √(Re²+Im²) of numerator / √(Re²+Im²) of denominator")
      out:t("∠H = ∠numerator - ∠denominator" .. ((mp ~= 0) and (" + " .. (-mp) .. "Ω") or ""))
    end
  end
  local ws = dt and {0, pi / 4, pi / 2, 3 * pi / 4, pi} or {0, 0.5, 1, 2, 5, 10, 100}
  for _, w in ipairs(ws) do
    local okv, v = pcall(f, w)
    if okv and v and v.re == v.re then
      out:t((dt and "  Ω = " or "  ω = ") .. fmtR(w) .. ":  |H| = " .. fmtR(cabs(v)) .. ",  ∠H = " .. fmtR(carg(v)) .. " rad")
    else
      out:t((dt and "  Ω = " or "  ω = ") .. fmtR(w) .. ":  undefined (pole / limit)")
    end
  end
  if not dt and #combineFT(Hterms) == 1 then out:t("  As ω→∞: |H| → " .. fmtR(cln(cabs(f(1e7)) * 1e4) / 1e4) .. " (approx)") end
end
addTool({
  name = "Frequency response H",
  help = "Diff. eq. like y''+3y'+2y=x'+x  /  y[n]-0.5y[n-1]=x[n]  or system y(t)=x(t)+x(t-1). Optional input x.",
  fields = {{"equation", "y'+2y=x'+x"}, {"input x (optional)", "e^(-t)u(t)"}},
  run = function(v, ctx, out)
    local src = v[1]
    local lhs, rhs = splitEq(src)
    if not lhs then error("need an equation with =", 0) end
    local lt = preprocess(lhs):gsub("%s", "")
    local dt = (src:find("%[") ~= nil)
    local var = V(dt)
    local mode = dt and "D" or "F"
    local Hstr, hSig, Hterms
    if lt == "y(t)" or lt == "y[n]" then
      -- explicit system: find h by input δ
      out:h("Impulse response (input δ)")
      local defs = {}
      for k2, d in pairs(ctx.defs) do defs[k2] = d end
      defs.x = {param = var, body = parse(dt and "delta[n]" or "delta(t)"), dt = dt}
      local node = P(ctx, rhs, {x = true})
      local cv = mkConv(var, dt, ctx.env, defs)
      hSig = cv:sig(node)
      out:r("h" .. (dt and "[n]" or "(t)") .. " = " .. hSig:compact(var))
      out:h("Frequency response")
      out:t(dt and "x[n] = e^(jΩn) => y[n] = H(e^(jΩ))e^(jΩn),  H = Σ h[n]e^(-jΩn)" or "x(t) = e^(jωt) => y(t) = H(jω)e^(jωt),  H = ∫ h(t)e^(-jωt)dt")
      local T, I
      if dt then T, I = dtTransform(hSig, "D") else T, I = ctTransform(hSig, "F") end
      Hterms = T
      Hstr = ftSumStr(T, I, mode, dt)
      out:r((dt and "H(e^(jΩ)) = " or "H(jω) = ") .. Hstr)
      if not dt and #T > 1 and #combineFT(T) == 1 then
        local parts = {}
        for _, g in ipairs(combineFT(T)) do
          local ds = delayStr(g.d, "F")
          parts[#parts + 1] = (ds ~= "" and (ds .. "·[") or "") .. ratStr(g.N, g.D, "F") .. (ds ~= "" and "]" or "")
        end
        out:r("       = " .. joinsum(parts))
      end
    else
      local fn = {y = true, x = true}
      local L = deSide(parse(lhs, {fnames = fn}), dt, ctx.env)
      local R = deSide(parse(rhs, {fnames = fn}), dt, ctx.env)
      local A = padd(L.Y, pscale(R.Y, -1))
      local B = padd(R.X, pscale(L.X, -1))
      out:h("Transform both sides")
      if dt then
        out:t("Use x[n-k] <-> e^(-jΩk)X  (delay property)")
        out:t("(" .. pstrQ(A, "D") .. ")Y = (" .. pstrQ(B, "D") .. ")X")
        Hstr = "(" .. pstrQ(B, "D") .. ")/(" .. pstrQ(A, "D") .. ")"
        out:r("H(e^(jΩ)) = Y/X = " .. ratStr(B, A, "D"))
      else
        out:t("Use d^k/dt^k <-> (jω)^k  (differentiation property)")
        out:t("(" .. pstr(A, "jω", true) .. ")Y(jω) = (" .. pstr(B, "jω", true) .. ")X(jω)")
        Hstr = "(" .. pstr(B, "jω", true) .. ")/(" .. pstr(A, "jω", true) .. ")"
        out:r("H(jω) = Y(jω)/X(jω) = " .. ratStr(B, A, "F"))
      end
      out:h("Impulse response h = inverse of H")
      hSig = inverseXform(Hstr, mode, ctx.env, nil, out)
      out:r("h" .. (dt and "[n]" or "(t)") .. " = " .. hSig:compact(var))
      if dt then Hterms = dtTransform(hSig, "D") else Hterms = ctTransform(hSig, "F") end
    end
    pcall(magPhase, out, Hterms, dt)
    out:plot({kind = "freq", f = function(w) if dt then return dtEval(Hterms, w) end return ftEval(Hterms, cx(0, w)) end, label = "|H|", mode = mode, phase = true})
    out:plot({kind = "sig", sig = hSig, label = "h"})
    if v[2] and v[2] ~= "" then
      out:h("Output for input x")
      local X = toSig(ctx, P(ctx, v[2]), dt)
      local T, I
      if dt then T, I = dtTransform(X, "D") else T, I = ctTransform(X, "F") end
      local Xs = ftSumStr(T, I, mode, dt)
      out:t((dt and "X(e^(jΩ)) = " or "X(jω) = ") .. Xs)
      out:t("Y = X·H  (convolution property)")
      local ok, Y = pcall(function()
        local o2 = newOut()
        local sig = inverseXform("(" .. Xs .. ")*(" .. Hstr .. ")", mode, ctx.env, nil, o2)
        for _, l in ipairs(o2.lines) do out.lines[#out.lines + 1] = l end
        return sig
      end)
      if not ok then
        out:t("(using convolution instead: " .. tostring(Y) .. ")")
        Y = sconv(X, hSig)
      end
      out:r("y" .. (dt and "[n]" or "(t)") .. " = " .. Y:compact(var))
      out:plot({kind = "sig", sig = Y, label = "y"})
    end
  end,
})

-- ===== Fourier series =====
addTool({
  name = "Fourier series (CT/DT)",
  help = "Give one period of x (or a sum of sinusoids) and period T/N. Gives a_k with steps.",
  fields = {{"x over one period =", "(abs(t)<1)"}, {"period T or N =", "4"}, {"show k up to", "5"}},
  run = function(v, ctx, out)
    local node = P(ctx, v[1])
    local dt = isDTexpr(node, ctx)
    local var = V(dt)
    local T = realof(evalN(P(ctx, v[2]), ctx.env), "period")
    local K = round(realof(evalN(P(ctx, v[3] ~= "" and v[3] or "5"), ctx.env)))
    local X = toSig(ctx, node, dt)
    local w0 = 2 * pi / T
    out:h("Setup")
    out:t((dt and "N = " or "T = ") .. fmtR(T) .. ",  " .. (dt and "Ω0 = 2π/N = " or "ω0 = 2π/T = ") .. fmtR(w0))
    out:t(dt and "a_k = (1/N) Σ[n=<N>] x[n] e^(-jkΩ0 n)" or "a_k = (1/T) ∫[T] x(t) e^(-jkω0 t) dt")
    local periodicInput = (#X.bp == 0 and #X.imps == 0)
    local coef = {}
    if periodicInput then
      out:t("x is a sum of complex exponentials → read a_k directly (Euler):")
      out:t("x = " .. X:compact(var))
      for _, m in ipairs(X.seg[1]) do
        local w = dt and carg(m.s) or m.s.im
        local kk = w / w0
        if m.k > 0 or (not dt and abs(m.s.re) > 1e-12) or (dt and abs(cabs(m.s) - 1) > 1e-9) or not isint(kk) then
          error("not periodic with this period", 0)
        end
        kk = round(kk)
        if dt then kk = kk % round(T) end
        coef[kk] = (coef[kk] or ZERO) + m.c
      end
      local keys = {}
      for kk in pairs(coef) do keys[#keys + 1] = kk end
      table.sort(keys)
      out:h("Coefficients")
      for _, kk in ipairs(keys) do out:r("a_" .. kk .. " = " .. fmtC(coef[kk]) .. ((not isreal(coef[kk])) and ("   (= " .. fmtPolar(coef[kk]) .. ")") or "")) end
      out:t("all other a_k = 0")
      local pw = 0
      for _, kk in ipairs(keys) do pw = pw + cabs(coef[kk]) ^ 2 end
      out:t("Parseval: average power P = Σ|a_k|² = " .. fmtR(pw))
      return
    end
    -- one period: restrict to [0,T) or the given support
    local lo, hi
    if dt then lo, hi = 0, round(T) - 1 else lo, hi = 0, T end
    -- if signal support is centered (e.g. |t|<1 within (-T/2,T/2)) use it
    local sLo, sHi = X.bp[1] or -INF, X.bp[#X.bp] or INF
    if #X.seg[1] == 0 and #X.seg[#X.seg] == 0 and sLo >= -T / 2 - 1e-9 and sHi <= T / 2 + 1e-9 and sLo < 0 then lo, hi = -T / 2, T / 2 end
    if dt and #X.seg[1] == 0 and #X.seg[#X.seg] == 0 and sLo < 0 then lo = round(sLo); hi = lo + round(T) - 1 end
    out:t("Integrate over one period: " .. fmtR(lo) .. " → " .. fmtR(hi))
    local Xp = smul(X, dt and swindow(lo, hi + 1, true) or swindow(lo, hi, false))
    if dt then
      local parts = {}
      for n = lo, hi do parts[#parts + 1] = fmtC(Xp:eval(n)) end
      out:t("x[n] over period: {" .. table.concat(parts, ", ") .. "}")
    else
      out:t("x(t) on period = " .. Xp:compact("t"))
      -- general formula when pieces are exponentials (k=0)
      local simple = true
      for _, pc in ipairs(Xp:pieces()) do for _, m in ipairs(pc.S) do if m.k > 0 then simple = false end end end
      if simple then
        local terms = {}
        for _, pc in ipairs(Xp:pieces()) do
          for _, m in ipairs(pc.S) do
            local lam = "(" .. (isz(m.s) and "" or (fmtC(m.s) .. " ")) .. "- jk" .. fmtR(w0) .. ")"
            if isz(m.s) then lam = "(-jk" .. fmtR(w0) .. ")" end
            terms[#terms + 1] = cmul(m.c, "[e^(" .. lam .. "·" .. fmtR(pc.hi) .. ") - e^(" .. lam .. "·" .. fmtR(pc.lo) .. ")]/" .. lam)
          end
        end
        out:t("General (k≠0): a_k = (1/" .. fmtR(T) .. ")·{" .. joinsum(terms) .. "}")
      end
    end
    out:h("Coefficients")
    local acoef = {}
    local k1, k2 = -K, K
    if dt then k1, k2 = 0, round(T) - 1 end
    for kk = k1, k2 do
      do
        local val
        if dt then
          val = ZERO
          for n = lo, hi do val = val + Xp:eval(n) * cexp(cx(0, -kk * w0 * n)) end
          val = val / T
        else
          local e = sexp(cx(0, -kk * w0), false)
          val = smul(Xp, e):integrate(lo, hi) / T
          for _, im in ipairs(Xp.imps) do end
        end
        val = ccln(val)
        local ps = fmtPolar(val)
        out:r("a_" .. kk .. " = " .. fmtC(val) .. ((not isreal(val) and not isreal(val * J)) and ("  = " .. ps) or ""))
        acoef[#acoef + 1] = {k = kk, v = val}
      end
    end
    if dt then out:t("(a_k is periodic with period N)") end
    -- Parseval & synthesis
    local okP, Pv = pcall(function()
      local x2 = smul(Xp, conjSig(Xp))
      if dt then return x2:integrate(lo, hi) / T end
      return x2:integrate(lo, hi) / T
    end)
    if okP then out:t("Parseval: P = (1/" .. (dt and "N" or "T") .. ")" .. (dt and "Σ|x|²" or "∫|x|²") .. " over a period = " .. fmtC(Pv) .. " = Σ over all k of |a_k|²") end
    local sp = 0
    for _, e in ipairs(acoef) do sp = sp + cabs(e.v) ^ 2 end
    out:t("Sum of |a_k|² for the listed k = " .. fmtR(sp) .. (dt and "" or " (approaches P as more k are added)"))
    out:t("Synthesis: x = Σ a_k e^(jk" .. (dt and "Ω0n" or "ω0t") .. ")  (partial sum plotted)")
    out:plot({kind = "fn", dt = dt, label = "partial FS sum", f = function(tt)
      local s = ZERO
      for _, e in ipairs(acoef) do s = s + e.v * cexp(cx(0, e.k * w0 * tt)) end
      return s
    end})
  end,
})

-- ===== Complex numbers =====
addTool({
  name = "Complex numbers",
  help = "Cartesian/polar with Euler steps. e.g. (3+j)/(1-j), 2e^(j3pi/4)",
  fields = {{"z =", "(-2+j)/(1+j)"}},
  run = function(v, ctx, out)
    local node = P(ctx, v[1])
    out:h("Steps")
    -- Euler expansions
    contains(node, function(n)
      if n.k == "call" and n.name == "exp" then
        local ok, a = pcall(evalN, n.args[1], ctx.env)
        if ok and abs(a.im) > 1e-12 then
          local th = a.im
          local pre = abs(a.re) > 1e-12 and ("e^(" .. fmtR(a.re) .. ")·") or ""
          out:t(tostr(n) .. " = " .. pre .. "(cos(" .. fmtR(th) .. ") + j sin(" .. fmtR(th) .. ")) = " .. fmtC(cexp(a)))
        end
      end
      return false
    end)
    if node.k == "div" then
      local a, b = evalN(node.a, ctx.env), evalN(node.b, ctx.env)
      out:t("Numerator = " .. fmtC(a) .. ",  denominator = " .. fmtC(b))
      out:t("Multiply top & bottom by conj(den) = " .. fmtC(cconj(b)) .. ":")
      out:t("  = (" .. fmtC(a * cconj(b)) .. ")/" .. fmtR(cabs(b) ^ 2))
      out:t("Polar: |num|/|den| = " .. fmtR(cabs(a)) .. "/" .. fmtR(cabs(b)) .. ",  ∠num - ∠den = " .. fmtR(carg(a)) .. " - (" .. fmtR(carg(b)) .. ")")
    end
    local z = evalN(node, ctx.env)
    out:h("Result")
    out:r("Cartesian: " .. fmtC(z))
    local m, ang = cabs(z), carg(z)
    local function sq(x) local s = fmtR(x); if x < 0 or s:find("[/ ]") then s = "(" .. s .. ")" end; return s .. "²" end
    out:t("|z| = √(" .. sq(z.re) .. " + " .. sq(z.im) .. ") = " .. fmtR(m))
    local q = ""
    if z.re < 0 then q = "  (Re<0: arctan(" .. fmtR(z.im / z.re) .. ") " .. (z.im >= 0 and "+ π" or "- π") .. ")"
    elseif abs(z.re) > 1e-12 then q = "  (= arctan(" .. fmtR(z.im / z.re) .. "))" end
    out:t("∠z = " .. fmtR(ang) .. " rad = " .. fmtR(ang * 180 / pi) .. "°" .. q)
    out:r("Polar: " .. fmtPolar(z))
    out:t("Re = " .. fmtR(z.re) .. ",  Im = " .. fmtR(z.im) .. ",  conj = " .. fmtC(cconj(z)))
  end,
})

-- ===== Graph =====
addTool({
  name = "Graph signals",
  help = "Plot up to 3 signals (CT in t, DT in n). Impulses drawn as arrows.",
  fields = {{"f1 =", "u(t+1)-u(t-1)"}, {"f2 =", ""}, {"f3 =", ""}},
  run = function(v, ctx, out)
    out:h("Plots")
    for i = 1, 3 do
      if v[i] and v[i] ~= "" then
        local node = P(ctx, v[i])
        local dt = isDTexpr(node, ctx)
        local ok, sig = pcall(toSig, ctx, node, dt)
        if ok then
          out:t("f" .. i .. " = " .. sig:compact(V(dt)))
          out:plot({kind = "sig", sig = sig, label = "f" .. i})
        else
          out:t("f" .. i .. ": numeric plot (" .. tostring(sig) .. ")")
          out:plot({kind = "fn", f = numFn(ctx, node, V(dt)), dt = dt, label = "f" .. i})
        end
      end
    end
    out:t("Press G / Menu > Graph to view.")
  end,
})
