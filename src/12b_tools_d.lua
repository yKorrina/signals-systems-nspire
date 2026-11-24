-- ===================== TOOLS (part D) =====================
-- CT signal sampled at t = n Ts -> DT signal
local function sampleSig(X, Ts)
  local pieces = {}
  for _, pc in ipairs(X:pieces()) do
    local lo = (pc.lo == -INF) and -INF or math.ceil(pc.lo / Ts - 1e-9)
    local hi = (pc.hi == INF) and INF or math.ceil(pc.hi / Ts - 1e-9)
    if hi > lo then
      local S = {}
      for _, m in ipairs(pc.S) do S[#S + 1] = M(m.c * cx(Ts ^ m.k), m.k, cexp(m.s * Ts)) end
      pieces[#pieces + 1] = {lo = lo, hi = hi, S = S}
    end
  end
  return fromPieces(true, pieces)
end

addTool({
  name = "Sampling & aliasing",
  help = "x(t) (sum of sinusoids or any signal) sampled every Ts: Nyquist check, x[n] = x(nTs), aliased frequencies, ideal reconstruction.",
  fields = {{"x(t) =", "cos(2pi*3t)+sin(2pi*7t)"}, {"Ts (or fs=... )", "fs=10"}},
  run = function(v, ctx, out)
    local X = toSig(ctx, P(ctx, v[1]), false)
    local s2 = v[2]:gsub("%s", "")
    local Ts
    if s2:match("^fs=") then Ts = 1 / realof(evalN(P(ctx, s2:sub(4)), ctx.env), "fs")
    elseif s2:match("^ws=") or s2:match("^ωs=") then Ts = 2 * pi / realof(evalN(P(ctx, (s2:gsub("^.-=", ""))), ctx.env), "ws")
    else Ts = realof(evalN(P(ctx, (s2:gsub("^Ts=", ""))), ctx.env), "Ts") end
    local ws = 2 * pi / Ts
    out:h("Sampling rate")
    out:t("Ts = " .. fmtR(Ts) .. " s,  fs = 1/Ts = " .. fmtR(1 / Ts) .. " Hz,  ωs = 2π/Ts = " .. fmtR(ws) .. " rad/s")
    out:h("Signal content")
    out:t("x(t) = " .. X:compact("t"))
    local periodicLike = (#X.bp == 0 and #X.imps == 0)
    local wmax = 0
    local comps = {}
    if periodicLike then
      for _, m in ipairs(X.seg[1]) do
        if abs(m.s.re) < 1e-12 and m.k == 0 then
          wmax = math.max(wmax, abs(m.s.im))
          comps[#comps + 1] = m
        end
      end
      local fl = {}
      for _, m in ipairs(comps) do if m.s.im >= 0 then fl[#fl + 1] = "ω = " .. fmtR(m.s.im) .. " (f = " .. fmtR(m.s.im / (2 * pi)) .. " Hz)" end end
      out:t("Frequencies present: " .. table.concat(fl, ", "))
      out:t("Highest frequency ωmax = " .. fmtR(wmax) .. "  (fmax = " .. fmtR(wmax / (2 * pi)) .. " Hz)")
      out:h("Nyquist check")
      out:t("Need ωs > 2ωmax = " .. fmtR(2 * wmax) .. "  (Nyquist rate, fs > " .. fmtR(wmax / pi) .. " Hz)")
      if ws > 2 * wmax + 1e-9 then out:r("ωs = " .. fmtR(ws) .. " > " .. fmtR(2 * wmax) .. " → NO aliasing; x(t) recoverable")
      elseif abs(ws - 2 * wmax) < 1e-9 then out:r("ωs = 2ωmax exactly → borderline (sine terms at ωmax can vanish)")
      else out:r("ωs = " .. fmtR(ws) .. " < " .. fmtR(2 * wmax) .. " → ALIASING occurs") end
    else
      out:t("Signal is not a finite sum of sinusoids → not band-limited in general (check |X(jω)|).")
    end
    out:h("Sampled sequence x[n] = x(n·Ts)")
    local Xd = sampleSig(X, Ts)
    out:t("Replace t by n·" .. fmtR(Ts) .. ":  ω t → (ω·Ts) n = Ω n")
    out:r("x[n] = " .. Xd:compact("n"))
    out:plot({kind = "sig", sig = Xd, label = "x[n]"})
    out:plot({kind = "sig", sig = X, label = "x(t)"})
    if periodicLike and #comps > 0 then
      out:h("Digital frequencies and aliasing")
      local rec = {}
      for _, m in ipairs(comps) do
        local Om = m.s.im * Ts
        local Oa = Om - 2 * pi * round(Om / (2 * pi))
        if abs(Oa + pi) < 1e-9 then Oa = pi end
        if m.s.im >= 0 then
          out:t("ω = " .. fmtR(m.s.im) .. ": Ω = ω·Ts = " .. fmtR(Om) .. (abs(Om - Oa) > 1e-9 and (" ≡ " .. fmtR(Oa) .. " (subtract multiples of 2π)") or "") .. " → appears at ω = Ω/Ts = " .. fmtR(abs(Oa) / Ts) .. (abs(abs(Oa) / Ts - m.s.im) > 1e-9 and "  ← ALIASED" or ""))
        end
        rec[#rec + 1] = M(m.c, 0, cx(0, Oa / Ts))
      end
      for _, m in ipairs(X.seg[1]) do if isz(m.s) then rec[#rec + 1] = M(m.c, 0, 0) end end
      local R = newSig(false); R.seg = {mergeS(rec)}; R:normalize()
      out:h("Ideal reconstruction (low-pass, cutoff ωs/2 = " .. fmtR(ws / 2) .. ")")
      out:r("x_r(t) = " .. R:compact("t"))
      if ssub(R, X):isZero() then out:t("x_r(t) = x(t): perfect reconstruction") else out:t("x_r(t) ≠ x(t) because of aliasing") end
      out:plot({kind = "sig", sig = R, label = "reconstructed"})
    end
  end,
})

addTool({
  name = "Interconnected systems",
  help = "Cascade (series) or parallel LTI systems: h = h1*h2 or h1+h2, H = H1·H2 or H1+H2.",
  fields = {{"h1 =", "e^(-t)u(t)"}, {"h2 =", "e^(-2t)u(t)"}, {"cascade or parallel", "cascade"}},
  run = function(v, ctx, out)
    local n1, n2 = P(ctx, v[1]), P(ctx, v[2])
    local dt = isDTexpr(n1, ctx) or isDTexpr(n2, ctx)
    local var = V(dt)
    local H1, H2 = toSig(ctx, n1, dt), toSig(ctx, n2, dt)
    local par = v[3]:lower():find("par") ~= nil
    out:h(par and "Parallel connection" or "Cascade (series) connection")
    out:t(par and "Outputs add: y = x*h1 + x*h2 = x*(h1 + h2)  ⇒  h = h1 + h2" or "Output of system 1 feeds system 2: y = (x*h1)*h2 = x*(h1*h2)  ⇒  h = h1*h2 (order does not matter)")
    out:t(par and "Frequency domain: H = H1 + H2" or "Frequency domain: H = H1·H2  (convolution property)")
    local Hh
    if par then
      Hh = sadd(H1, H2)
      out:r("h = " .. Hh:compact(var))
      out:plot({kind = "sig", sig = Hh, label = "h"})
    else
      Hh = convSteps(H1, H2, out)
    end
    ltiProps(Hh, out)
    local ok, T = pcall(function() if dt then return dtTransform(Hh, "D") end return ctTransform(Hh, "F") end)
    if ok then
      local mode = dt and "D" or "F"
      out:h("Overall frequency response")
      out:r((dt and "H(e^(jΩ)) = " or "H(jω) = ") .. ftSumStr(T, {}, mode, dt))
    end
  end,
})

addTool({
  name = "Convolution property check",
  help = "Verify y = x*h two ways: direct convolution and inverse transform of X·H.",
  fields = {{"x =", "e^(-t)u(t)"}, {"h =", "e^(-t)u(t)"}},
  run = function(v, ctx, out)
    local nx, nh = P(ctx, v[1]), P(ctx, v[2])
    local dt = isDTexpr(nx, ctx) or isDTexpr(nh, ctx)
    local var = V(dt)
    local mode = dt and "D" or "F"
    local X, H = toSig(ctx, nx, dt), toSig(ctx, nh, dt)
    out:h("Method 1: direct convolution")
    local Y1 = convSteps(X, H, out)
    out:h("Method 2: transforms")
    local TX, IX, TH, IH
    if dt then TX, IX = dtTransform(X, "D"); TH, IH = dtTransform(H, "D")
    else TX, IX = ctTransform(X, "F"); TH, IH = ctTransform(H, "F") end
    local Xs, Hs = ftSumStr(TX, IX, mode, dt), ftSumStr(TH, IH, mode, dt)
    out:t("X = " .. Xs)
    out:t("H = " .. Hs)
    out:t("Y = X·H = (" .. Xs .. ")·(" .. Hs .. ")")
    local o2 = newOut()
    local ok, Y2 = pcall(inverseXform, "(" .. Xs .. ")*(" .. Hs .. ")", mode, ctx.env, nil, o2)
    for _, l in ipairs(o2.lines) do out.lines[#out.lines + 1] = l end
    if ok then
      out:r("Inverse transform: y = " .. Y2:compact(var))
      local same = ssub(Y1, Y2):isZero()
      out:r(same and "Both methods agree (OK)" or "Results differ (check impulses/ROC)")
    else
      out:e("Inverse step failed: " .. tostring(Y2))
    end
  end,
})
