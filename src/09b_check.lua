-- ===================== NUMERIC SELF-CHECKS =====================
-- Each transform answer is compared with a direct numeric computation at a test point.
local function okLine(out, label, a, b)
  local good = near(a, b, 1e-5) or cabs(a - b) < 1e-7
  if good then
    out:t("Check at " .. label .. " computed directly = " .. decstr(cln(a.re)) .. (abs(a.im) > 1e-9 and ((a.im < 0) and " - j" or " + j") .. decstr(abs(cln(a.im))) or "") .. ", matches  (OK)")
  else
    out:e("Check FAILED at " .. label .. " direct = " .. fmtC(a) .. ", formula = " .. fmtC(b) .. "  (please report this input)")
  end
  return good
end

-- numeric ∫ x(t) e^{-s t} dt over the signal pieces, plus impulses
local function numLaplace(sig, s0)
  local total = ZERO
  for _, pc in ipairs(sig:pieces()) do
    local lo, hi = math.max(pc.lo, -60), math.min(pc.hi, 60)
    if hi > lo then
      local n = math.min(4000, math.max(200, math.ceil((hi - lo) * 60)))
      total = total + simpson(function(t) return segEval(pc.S, t, false) * cexp(-s0 * t) end, lo, hi, n)
    end
  end
  for _, im in ipairs(sig.imps) do total = total + im.c * cpow(s0, im.m) * cexp(-s0 * im.t0) end
  return total
end
local function numZ(sig, z0)
  local lo, hi = -400, 400
  local s = ZERO
  for n = lo, hi do
    local v = sig:eval(n)
    if not isz(v, 0) then s = s + v * cpow(z0, -n) end
  end
  return s
end

function checkForward(sig, mode, T, I, roc, out)
  if I and #I > 0 then
    out:t("Check: (impulses in the answer, numeric check skipped)")
    return
  end
  local ok = pcall(function()
    if mode == "Z" or mode == "D" then
      local z0
      if mode == "D" then z0 = cexp(cx(0, 0.7))
      else
        local r0
        if roc.hi == INF then r0 = math.max(roc.lo * 1.5, roc.lo + 0.5, 1.2) else r0 = (roc.lo + roc.hi) / 2 end
        z0 = cexp(cx(0, 0.7)) * r0
      end
      local direct = numZ(sig, z0)
      local formula = dtEval(T, (mode == "Z") and z0 or 0.7, mode == "Z")
      okLine(out, (mode == "Z") and ("z = " .. fmtR(cabs(z0)) .. "e^(j0.7):") or "Ω = 0.7:", direct, formula)
    else
      local s0
      if mode == "F" then s0 = cx(0, 0.7)
      else
        local sig0
        if roc.lo > -INF and roc.hi < INF then sig0 = (roc.lo + roc.hi) / 2
        elseif roc.lo > -INF then sig0 = roc.lo + 1
        elseif roc.hi < INF then sig0 = roc.hi - 1
        else sig0 = 0 end
        s0 = cx(sig0, 0.7)
      end
      local direct = numLaplace(sig, s0)
      local formula = ftEval(T, s0)
      okLine(out, (mode == "F") and "ω = 0.7:" or ("s = " .. fmtC(s0) .. ":"), direct, formula)
    end
  end)
  if not ok then out:t("Check: (numeric check not possible here)") end
end

-- inverse: transform the answer back and compare with the given X at a test point
function checkInverse(src, mode, env, sig, hasSynth, roc, out)
  if hasSynth then return end
  local node = parse(src)
  if contains(node, function(n) return n.k == "call" and (n.name == "delta" or n.name == "rect") or n.k == "cmp" end) then return end
  local ok, err = pcall(function()
    local dt = (mode == "D" or mode == "Z")
    local fv = FVAR[mode]
    if dt then
      local _, _, rc = dtTransform(sig, "Z")
      local r0
      if rc.hi == INF then r0 = math.max(rc.lo * 1.5, rc.lo + 0.5, 1.2) else r0 = (rc.lo + rc.hi) / 2 end
      if mode == "D" then r0 = 1 end
      local z0 = cexp(cx(0, 0.7)) * r0
      env.vars[fv] = (mode == "D") and cx(0.7) or z0
      local given = evalN(node, env)
      local back = numZ(sig, z0)
      okLine(out, ((mode == "D") and "Ω = 0.7" or ("z = " .. fmtR(r0) .. "e^(j0.7)")) .. ": Σ x[n]z^-n of the answer =", back, given)
    else
      local s0
      if mode == "F" then s0 = cx(0, 0.7); env.vars[fv] = cx(0.7)
      else
        local _, _, rc = ctTransform(sig, "L")
        local sg = (rc.lo > -INF and rc.hi < INF) and (rc.lo + rc.hi) / 2 or ((rc.lo > -INF) and rc.lo + 1 or ((rc.hi < INF) and rc.hi - 1 or 0))
        s0 = cx(sg, 0.7); env.vars[fv] = s0
      end
      local given = evalN(node, env)
      local back = numLaplace(sig, s0)
      okLine(out, ((mode == "F") and "ω = 0.7" or ("s = " .. fmtC(s0))) .. ": ∫x(t)e^(-st)dt of the answer =", back, given)
    end
  end)
end
