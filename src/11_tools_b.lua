-- ===================== TOOLS (part B) =====================
-- ---------- numeric system tester ----------
local function sysRunner(ctx, node, dt)
  local var = V(dt)
  -- returns function(inputFn) -> function(t) -> complex|nil
  return function(xf)
    return function(t)
      local env = ctx.env
      env.vars[var] = cx(t)
      env.tvar = var
      env.dh = 1e-4
      env.intCap = 240
      env.sumInf = 28
      local saved = env.fns.x
      env.fns.x = function(a)
        local r = tocx(a).re
        if dt then
          if abs(r - round(r)) > 1e-6 then return ZERO end
          r = round(r)
        end
        return tocx(xf(r))
      end
      local ok, val = pcall(evalN, node, env)
      env.fns.x = saved
      if not ok then return nil, val end
      return val
    end
  end
end
local function rnd(seed) local s = seed; return function() s = (s * 16807) % 2147483647; return s / 2147483647 end end

local function sysTests(ctx, node, dt)
  local run = sysRunner(ctx, node, dt)
  local R = {}
  local pts = dt and {-7, -3, -2, -1, 0, 1, 2, 3, 5, 8} or {-3.7, -1.3, -0.4, 0, 0.6, 1.1, 2.5, 3.9}
  local function x1(t) return dt and ((t >= -12 and t <= 12) and cx(math.sin(1.3 * t + 0.4) + 0.3 * math.cos(0.7 * t) + 0.2) or ZERO) or cx(math.exp(-(t - 0.3) ^ 2 / 4) * (math.cos(1.3 * t) + 0.5)) end
  local function x2(t) return dt and ((t >= -12 and t <= 12) and cx(0.5 * math.cos(2.1 * t - 1) - 0.4 + 0.1 * t) or ZERO) or cx(1 / (1 + t * t) * (math.sin(0.9 * t) + 0.8)) end
  local function eqf(fa, fb, tolf)
    for _, t in ipairs(pts) do
      local a, ea = fa(t); local b, eb = fb(t)
      if a == nil or b == nil then return nil, (ea or eb) end
      if not near(a, b, tolf or 1e-4) then return false, t, a, b end
    end
    return true
  end
  -- linearity
  local A, B = 2.3, -1.7
  local yc = run(function(t) return x1(t) * A + x2(t) * B end)
  local y1, y2 = run(x1), run(x2)
  local lin, lt, la, lb = eqf(yc, function(t) local a = y1(t); local b = y2(t); if not a or not b then return nil end; return a * A + b * B end)
  R.linear = {ok = lin, t = lt, a = la, b = lb}
  if lin == nil then R.err = lt end
  -- time invariance
  local t0 = dt and 2 or 0.7
  local ys = run(function(t) return x1(t - t0) end)
  local ti, tt, ta, tb = eqf(ys, function(t) return y1(t - t0) end)
  R.ti = {ok = ti, t = tt, a = ta, b = tb, t0 = t0}
  -- causality & memory: perturb future / past
  local function bump(t, c, side, width)
    -- smooth bump on side of c (side=1 future, -1 past), zero within width
    local d = (t - c) * side
    if dt then return (d >= 1) and cx(1.7 + 0.3 * d) or ZERO end
    if d <= width then return ZERO end
    local u = d - width
    return cx(u * u * u * math.exp(-u))
  end
  local causal, memless = true, true
  local cinfo, minfo
  for _, c in ipairs(dt and {-3, 0, 2, 5} or {-2, 0, 0.3, 1.7}) do
    local yf = run(function(t) return x1(t) + bump(t, c, 1, 0) end)
    local a, b = y1(c), yf(c)
    if a and b and not near(a, b, 1e-6) then causal = false; cinfo = cinfo or {c = c, a = a, b = b} end
    local yp = run(function(t) return x1(t) + bump(t, c, -1, dt and 0 or 0.01) end)
    local a2, b2 = y1(c), yp(c)
    if a2 and b2 and not near(a2, b2, 1e-6) then memless = false; minfo = minfo or {c = c, side = "past"} end
    if a and b and not near(a, b, 1e-6) then memless = false; minfo = minfo or {c = c, side = "future"} end
  end
  if contains(node, function(n) return n.k == "call" and (n.name == "der" or n.name == "der2") end) then memless = false; minfo = minfo or {deriv = true} end
  R.causal = {ok = causal, info = cinfo}
  R.memless = {ok = memless, info = minfo}
  -- stability
  local bounded = {
    {name = "u" .. (dt and "[n]" or "(t)"), f = function(t) return cx(t >= 0 and 1 or 0) end},
    {name = "1", f = function(t) return ONE end},
    {name = dt and "(-1)^n" or "cos(t)", f = function(t) return dt and cx((round(t) % 2 == 0) and 1 or -1) or cx(math.cos(t)) end},
  }
  local stable, sinfo = true, nil
  for _, bi in ipairs(bounded) do
    local y = run(bi.f)
    local function maxOn(L)
      local m = 0
      local step = dt and math.max(1, math.floor(L / 20)) or (L / 12)
      local t = -L
      while t <= L + 1e-9 do
        local v = y(t)
        if v then m = math.max(m, cabs(v)) end
        t = t + step
      end
      local v0 = y(0); if v0 then m = math.max(m, cabs(v0)) end
      return m
    end
    local m10, m40 = maxOn(10), maxOn(40)
    if m10 > 1e3 or m40 > 1e3 or (m40 > 3 * m10 + 1e-9 and m40 > 15) then
      stable = false
      sinfo = {name = bi.name, m10 = m10, m40 = m40}
      break
    end
  end
  R.stable = {ok = stable, info = sinfo}
  -- invertibility: search for two different inputs with same output
  local cands = {
    {"x = 0", function() return ZERO end, "x = 1", function() return ONE end},
    {"x", x1, "-x", function(t) return -x1(t) end},
    {"x = 0", function() return ZERO end, dt and "x = (-1)^n" or "x = odd signal t·e^(-t²)", function(t) return dt and cx((round(t) % 2 == 0) and 1 or -1) or cx(t * math.exp(-t * t)) end},
    {"x = 0", function() return ZERO end, "x = even signal", function(t) return cx(math.exp(-t * t / 4)) end},
    {"x", x1, "x + (pulse at negative times)", function(t) return x1(t) + bump(t, -1, -1, 0) end},
    {"x", x1, "x + (pulse at positive times)", function(t) return x1(t) + bump(t, 1, 1, 0) end},
  }
  if dt then
    cands[#cands + 1] = {"x", x1, "x + δ[n-1]", function(t) return x1(t) + cx(t == 1 and 1 or 0) end}
    cands[#cands + 1] = {"x", x1, "x + δ[n-3]", function(t) return x1(t) + cx(t == 3 and 1 or 0) end}
  end
  local inv, iinfo = true, nil
  for _, cd in ipairs(cands) do
    local ya, yb = run(cd[2]), run(cd[4])
    local same = true
    local tp = dt and {-9, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 9} or {-3.3, -2.1, -1.5, -0.7, -0.2, 0, 0.4, 0.9, 1.6, 2.2, 3.1}
    for _, t in ipairs(tp) do
      local a, b = ya(t), yb(t)
      if not a or not b or not near(a, b, 1e-6) then same = false; break end
    end
    if same then inv = false; iinfo = cd[1] .. "  and  " .. cd[3]; break end
  end
  R.invertible = {ok = inv, info = iinfo}
  return R
end

local function xArgs(node)
  local args = {}
  local function walk(n, note)
    if n.k == "call" and n.name == "x" then args[#args + 1] = tostr(n.args[1]) .. note; return end
    if n.k == "call" and (n.name == "int" or n.name == "sum") and #n.args == 4 then
      local dv = tostr(n.args[2])
      walk(n.args[1], " (" .. dv .. " from " .. tostr(n.args[3]) .. " to " .. tostr(n.args[4]) .. ")")
      return
    end
    for kk, v in pairs(n) do
      if type(v) == "table" and v.k then walk(v, note)
      elseif kk == "args" then for _, x in ipairs(v) do walk(x, note) end end
    end
  end
  walk(node, "")
  return args
end

-- symbolic bound of |y| given |x| <= B
local function boundOf(n, ctx, var)
  local k = n.k
  local function rec(m) return boundOf(m, ctx, var) end
  if k == "num" then return fmtR(cabs(n.v)) end
  if k == "var" then
    if n.name == var then return nil end
    local ok, v = pcall(evalN, n, ctx.env)
    if ok then return fmtR(cabs(v)) end
    return nil
  end
  if k == "add" or k == "sub" then
    local a, b = rec(n.a), rec(n.b)
    if a and b then return a .. " + " .. b end
    return nil
  end
  if k == "neg" then return rec(n.a) end
  if k == "mul" then
    local a, b = rec(n.a), rec(n.b)
    if not a or not b then return nil end
    if a == "1" then return b elseif b == "1" then return a end
    local function w(s) return topSum(s) and ("(" .. s .. ")") or s end
    return w(a) .. "·" .. w(b)
  end
  if k == "div" then
    local cb = isconst(n.b, var, ctx.env)
    local a = rec(n.a)
    if cb and a and not isz(cb) then return "(" .. a .. ")/" .. fmtR(cabs(cb)) end
    return nil
  end
  if k == "pow" then
    local ce = isconst(n.b, var, ctx.env)
    local a = rec(n.a)
    if ce and a then return "(" .. a .. ")^" .. fmtR(ce.re) end
    return nil
  end
  if k == "call" then
    local nm = n.name
    if nm == "x" then return "B" end
    if nm == "u" or nm == "rect" or nm == "sgn" or nm == "sign" or nm == "step" then return "1" end
    if nm == "cos" or nm == "sin" then
      if not contains(n.args[1], function(m) return m.k == "var" and (m.name == "j" or m.name == "i") end) then return "1" end
      return nil
    end
    if nm == "exp" then
      -- e^{j·real} has magnitude 1
      local a = n.args[1]
      local imagOnly = false
      if a.k == "mul" and ((a.a.k == "var" and (a.a.name == "j" or a.a.name == "i")) or (a.b.k == "var" and (a.b.name == "j" or a.b.name == "i"))) then imagOnly = true end
      if a.k == "neg" and a.a.k == "mul" and ((a.a.a.k == "var" and a.a.a.name == "j")) then imagOnly = true end
      if imagOnly and not contains(a, function(m) return m.k == "call" and m.name == "x" end) then return "1" end
      local b = rec(a)
      if b then return "e^(" .. b .. ")" end
      return nil
    end
    if nm == "abs" then return rec(n.args[1]) end
    if (nm == "int" or nm == "sum") and #n.args == 4 then
      local a1, b1 = linearOf(n.args[3], var, ctx.env)
      local a2, b2 = linearOf(n.args[4], var, ctx.env)
      if a1 and a2 and near(a1, a2) then
        local Lr = realof(b2 - b1)
        local inner = boundOf(n.args[1], ctx, var)
        if inner and Lr >= 0 then
          local cnt = (nm == "sum") and (Lr + 1) or Lr
          return fmtR(cnt) .. "·" .. (topSum(inner) and ("(" .. inner .. ")") or inner) .. ((nm == "sum") and "   (" .. fmtR(cnt) .. " terms)" or "   (interval length " .. fmtR(Lr) .. ")")
        end
      end
      return nil
    end
    local c0 = isconst(n, var, ctx.env)
    if c0 then return fmtR(cabs(c0)) end
  end
  return nil
end
-- exact output for a given input expression (or nil)
local function exactY(ctx, node, dt, inputSrc)
  local var = dt and "n" or "t"
  local ok, r = pcall(function()
    local defs = {}
    for k2, d in pairs(ctx.defs) do defs[k2] = d end
    defs.x = {param = var, body = parse(inputSrc), dt = dt}
    local cv = mkConv(var, dt, ctx.env, defs)
    return cv:sig(node)
  end)
  if ok then return r end
  return nil
end
-- arguments of x with info
local function xArgList(node, var, ctx)
  local list = {}
  local function walk(n, intinfo)
    if n.k == "call" and n.name == "x" then list[#list + 1] = {arg = n.args[1], int = (intinfo and not intinfo.der) and intinfo or nil, der = intinfo and intinfo.der}; return end
    if n.k == "call" and (n.name == "int" or n.name == "sum") and #n.args == 4 then
      walk(n.args[1], {dv = tostr(n.args[2]), lo = n.args[3], hi = n.args[4], name = n.args[2].name})
      return
    end
    if n.k == "call" and (n.name == "der" or n.name == "der2") then
      walk(n.args[1], {der = true})
      return
    end
    for kk, v in pairs(n) do
      if type(v) == "table" and v.k then walk(v, intinfo)
      elseif kk == "args" then for _, x in ipairs(v) do walk(x, intinfo) end end
    end
  end
  walk(node, nil)
  return list
end

addTool({
  name = "System properties",
  help = "Enter y(t)=... using x(...). Checks memoryless, time-invariant, linear, causal, stable, invertible with full reasoning.",
  fields = {{"y =", "x(t-1) + x(t+1)"}},
  run = function(v, ctx, out)
    local src = v[1]
    local lhs, rhs = splitEq(src)
    if lhs then src = rhs end
    local node = P(ctx, src, {x = true})
    local dt = isDTexpr(node, ctx) or (src:find("%[") ~= nil and not hasvar(node, "t"))
    local var = V(dt)
    local o, c = dt and "[" or "(", dt and "]" or ")"
    local yv = "y" .. o .. var .. c
    local xv = "x" .. o .. var .. c
    out:h("System")
    out:r(yv .. " = " .. tostr(node))
    local R = sysTests(ctx, node, dt)
    if R.err then out:e("Evaluation problem: " .. tostring(R.err)) end
    local args = xArgList(node, var, ctx)
    local t0n = N("var", {name = dt and "n0" or "t0"})
    local t0s = dt and "n0" or "t0"
    local function yn(b) return b and "YES" or "NO" end
    out:h("Summary")
    out:r("Memoryless: " .. yn(R.memless.ok) .. "   Time-inv: " .. yn(R.ti.ok))
    out:r("Linear: " .. yn(R.linear.ok) .. "   Causal: " .. yn(R.causal.ok))
    out:r("Stable (BIBO): " .. yn(R.stable.ok) .. "   Invertible: " .. yn(R.invertible.ok))
    local function showY(label, S) if S then out:t("    " .. label .. " = " .. S:compact(var)) end end
    -- which times does y use?
    out:h("Which input values does " .. yv .. " use?")
    local samples = dt and {-2, 0, 3} or {-2, 0, 2}
    local future, other = false, false
    for _, a in ipairs(args) do
      if a.der then
        out:t("  d/d" .. var .. " of x" .. o .. tostr(a.arg) .. c .. ": needs x just around " .. var .. " (limit of a difference quotient)")
        other = true
      elseif a.int then
        out:t("  x" .. o .. tostr(a.arg) .. c .. " with " .. a.int.dv .. " from " .. tostr(a.int.lo) .. " to " .. tostr(a.int.hi))
        local okh, hv = pcall(function() local env = ctx.env; env.vars[var] = cx(samples[3]); return evalN(a.int.hi, env) end)
        if okh and hv.re > samples[3] + 1e-9 then out:t("     upper limit " .. tostr(a.int.hi) .. " > " .. var .. " for " .. var .. " = " .. samples[3] .. "  → future values"); future = true end
        other = true
      else
        local diffs = {}
        local fut = nil
        for _, tt in ipairs(dt and {-3, -1, 0, 1, 3} or {-3, -1, -0.5, 0.5, 1, 3}) do
          ctx.env.vars[var] = cx(tt)
          local okv, av = pcall(evalN, a.arg, ctx.env)
          if okv then
            if abs(av.re - tt) > 1e-9 then other = true end
            if av.re > tt + 1e-9 then fut = fut or {tt, av.re} end
          end
        end
        local same = not fut and (function()
          for _, tt in ipairs({-3, 0.5, 2}) do ctx.env.vars[var] = cx(tt); local okv, av = pcall(evalN, a.arg, ctx.env); if okv and abs(av.re - tt) > 1e-9 then return false end end
          return true
        end)()
        local msg = "  x" .. o .. tostr(a.arg) .. c .. ": "
        if same then msg = msg .. "argument equals " .. var .. " (present value)"
        elseif fut then msg = msg .. "at " .. var .. " = " .. fmtR(fut[1]) .. " it reads x at " .. fmtR(fut[2]) .. " > " .. fmtR(fut[1]) .. " → FUTURE"; future = true
        else msg = msg .. "argument ≤ " .. var .. " always → past/present" end
        out:t(msg)
      end
    end
    -- (i)
    out:h("(i) Memoryless")
    out:t("Definition: y at time " .. var .. " depends only on x at the same time " .. var .. ".")
    if R.memless.ok then out:t("Every x term is evaluated at " .. var .. " itself → MEMORYLESS.")
    else
      if R.memless.info and R.memless.info.deriv then out:t("A derivative needs values of x at nearby times → has MEMORY (not memoryless).")
      else out:t(yv .. " uses x at times other than " .. var .. " (see list above) → NOT memoryless.") end
    end
    -- (ii)
    out:h("(ii) Time invariance")
    out:t("Test: delay the input by " .. t0s .. " and compare with the delayed output.")
    out:t("1) Input r" .. o .. var .. c .. " = x" .. o .. var .. " - " .. t0s .. c .. ". Replace every x(·) by x(· - " .. t0s .. "):")
    local zn = mapcalls(node, "x", function(a, orig) return N("call", {name = "x", args = {N("sub", {a = a[1], b = t0n})}, br = orig.br}) end)
    out:t("   z" .. o .. var .. c .. " = " .. tostr(zn))
    out:t("2) Delay the output: replace every " .. var .. " by (" .. var .. " - " .. t0s .. "):")
    local ysh = subst(node, var, N("sub", {a = N("var", {name = var}), b = t0n}))
    out:t("   y" .. o .. var .. " - " .. t0s .. c .. " = " .. tostr(ysh))
    if contains(node, function(nn) return nn.k == "call" and (nn.name == "int" or nn.name == "sum") end) then
      out:t("   (for integrals/sums substitute the dummy variable, e.g. τ = σ - " .. t0s .. ", then compare limits)")
    end
    if R.ti.ok then out:t("3) z" .. o .. var .. c .. " = y" .. o .. var .. " - " .. t0s .. c .. " for every input and shift → TIME-INVARIANT.")
    else
      out:t("3) They differ → TIME-VARYING.")
      local function hasTexplicit()
        local stripped = mapcalls(node, "x", function() return NUM(1) end)
        return hasvar(stripped, var)
      end
      if hasTexplicit() then out:t("   Reason: " .. var .. " appears outside of x(·) (a time-varying gain/limit) or x has a scaled/reversed argument.")
      else out:t("   Reason: the argument of x is scaled or reversed (not just " .. var .. " + const).") end
      local inp = dt and "u[n]" or "u(t)"
      local y1 = exactY(ctx, node, dt, dt and "u[n]" or "u(t)")
      local z1 = exactY(ctx, node, dt, dt and "u[n-2]" or "u(t-1)")
      if y1 and z1 then
        local d0 = dt and 2 or 1
        out:t("   Example: x = " .. inp .. " gives y = " .. y1:compact(var))
        out:t("   input " .. (dt and "u[n-2]" or "u(t-1)") .. " gives z = " .. z1:compact(var))
        out:t("   but y" .. o .. shiftStr(var, d0) .. c .. " = " .. y1:compose(1, -d0):compact(var) .. "  ≠ z")
      elseif R.ti.t then out:t("   e.g. " .. t0s .. "=" .. R.ti.t0 .. ", at " .. var .. "=" .. fmtR(R.ti.t) .. ": z = " .. fmtC(R.ti.a) .. ", y shifted = " .. fmtC(R.ti.b)) end
    end
    -- (iii)
    out:h("(iii) Linearity")
    out:t("Test superposition: input x = a·x1 + b·x2 must give a·y1 + b·y2.")
    local ln = mapcalls(node, "x", function(a, orig)
      return N("add", {a = N("mul", {a = N("var", {name = "a"}), b = N("call", {name = "x1", args = a, br = orig.br})}), b = N("mul", {a = N("var", {name = "b"}), b = N("call", {name = "x2", args = a, br = orig.br})})})
    end)
    out:t("  Output for a·x1 + b·x2:  " .. tostr(ln))
    local l1 = mapcalls(node, "x", function(a, orig) return N("call", {name = "x1", args = a, br = orig.br}) end)
    local l2 = mapcalls(node, "x", function(a, orig) return N("call", {name = "x2", args = a, br = orig.br}) end)
    out:t("  a·y1 + b·y2 = a[" .. tostr(l1) .. "] + b[" .. tostr(l2) .. "]")
    if R.linear.ok then
      out:t("  Expanding (distributing a, b through the operations) gives the same expression → LINEAR.")
      out:t("  (also: zero input gives zero output)")
    else
      out:t("  These are NOT equal → NONLINEAR.")
      local one = dt and "1" or "1"
      local y1 = exactY(ctx, node, dt, dt and "u[n]" or "u(t)")
      local y2 = exactY(ctx, node, dt, dt and "2u[n]" or "2u(t)")
      local y0 = exactY(ctx, node, dt, "0")
      if y0 and not y0:isZero() then out:t("  Quick test: x = 0 gives y = " .. y0:compact(var) .. " ≠ 0 (a linear system must give 0)") end
      if y1 and y2 then
        out:t("  Counterexample (homogeneity, a = 2):")
        out:t("    x = " .. (dt and "u[n]" or "u(t)") .. "  → y = " .. y1:compact(var))
        out:t("    x = " .. (dt and "2u[n]" or "2u(t)") .. " → y = " .. y2:compact(var))
        out:t("    2·(first) = " .. y1:scale(2):compact(var) .. "  ≠ second")
      elseif R.linear.t then out:t("  e.g. numeric test at " .. var .. "=" .. fmtR(R.linear.t) .. ": " .. fmtC(R.linear.a) .. " vs " .. fmtC(R.linear.b)) end
    end
    -- (iv)
    out:h("(iv) Causality")
    out:t("Definition: y at time " .. var .. " depends only on present and PAST input values.")
    if R.memless.ok and R.causal.ok then out:t("Memoryless ⇒ causal.")
    elseif R.causal.ok then out:t("All x values used are at times ≤ " .. var .. " (see list above) → CAUSAL.")
    else
      out:t("y uses future values of x (see list above) → NON-CAUSAL.")
      if R.causal.info then out:t("  Check: change x only after " .. var .. " = " .. fmtR(R.causal.info.c) .. "; y(" .. fmtR(R.causal.info.c) .. ") changes from " .. fmtC(R.causal.info.a) .. " to " .. fmtC(R.causal.info.b)) end
    end
    -- (v)
    out:h("(v) BIBO stability")
    out:t("Definition: every bounded input |x| ≤ B must give a bounded output.")
    local bnd = boundOf(node, ctx, var)
    if R.stable.ok then
      if bnd then
        out:t("Assume |x" .. o .. var .. c .. "| ≤ B. Using |a+b| ≤ |a|+|b|, |ab| = |a||b|, |e^(jθ)| = 1:")
        out:t("  |" .. yv .. "| ≤ " .. bnd .. "  < ∞")
        out:t("→ output is bounded → STABLE.")
      else
        out:t("Tested bounded inputs (u, 1, " .. (dt and "(-1)^n" or "cos t") .. "): outputs stay bounded → STABLE.")
      end
    else
      local s = R.stable.info
      out:t("Find one bounded input with an unbounded output:")
      local inp = dt and "u[n]" or "u(t)"
      local ys = exactY(ctx, node, dt, inp)
      if ys and s.name:find("^u") then
        out:t("  x = " .. inp .. " (|x| ≤ 1) gives y = " .. ys:compact(var))
        if #ys.imps > 0 then out:t("  y contains an impulse (infinite height) → unbounded")
        else out:t("  this grows without bound as " .. var .. " → ±∞") end
      else
        out:t("  x = " .. s.name .. " (bounded) gives an output that grows:")
        if s.m10 > 1e3 then out:t("  output contains an impulse/very large values (e.g. derivative of a step).")
        else out:t("  max|y| on ±10: " .. fmtR(s.m10) .. ",  on ±40: " .. fmtR(s.m40)) end
      end
      out:t("→ UNSTABLE.")
    end
    -- (vi)
    out:h("(vi) Invertibility")
    out:t("Definition: different inputs always give different outputs (input can be recovered).")
    if R.invertible.ok then
      out:t("No two different test inputs gave the same output → INVERTIBLE.")
      out:t("(Try to solve for x in terms of y to write the inverse system.)")
    else
      out:t("Different inputs " .. R.invertible.info .. " give the SAME output → NOT invertible.")
      local pairs0 = {{"0", "1"}}
      if R.invertible.info:find("^x = 0") and R.invertible.info:find("x = 1") then
        local a, b = exactY(ctx, node, dt, "0"), exactY(ctx, node, dt, "1")
        if a and b then out:t("  x = 0 → y = " .. a:compact(var) .. ";   x = 1 → y = " .. b:compact(var)) end
      end
    end
    out:t("(Each verdict was also checked numerically with random test signals.)")
  end,
})

-- ===== Tool: LTI system from h =====
local function absSig(sig)
  -- exact |h| when each piece has one monomial with k=0; otherwise triangle bound
  local pieces, exact = {}, true
  for _, pc in ipairs(sig:pieces()) do
    local S = {}
    if #pc.S > 1 then exact = false end
    for _, m in ipairs(pc.S) do
      if m.k > 0 then
        -- t^k sign may change at 0: split
        exact = exact and (pc.lo >= 0 or pc.hi <= 0)
      end
      local s2 = sig.dt and cx(cabs(m.s)) or cx(m.s.re)
      S[#S + 1] = M(cabs(m.c), m.k, s2)
      if m.k % 2 == 1 and pc.hi <= 0 then S[#S].c = -S[#S].c end
    end
    pieces[#pieces + 1] = {lo = pc.lo, hi = pc.hi, S = S}
  end
  local r = fromPieces(sig.dt, pieces)
  r.imps = {}
  for _, im in ipairs(sig.imps) do r.imps[#r.imps + 1] = {c = cx(cabs(im.c)), t0 = im.t0, m = 0}; if im.m > 0 then exact = false end end
  return r, exact
end
local function ltiProps(H, out)
  local dt = H.dt
  local var = V(dt)
  local o, c = dt and "[" or "(", dt and "]" or ")"
  out:h("Causality")
  local nonc = nil
  for i, S in ipairs(H.seg) do
    local lo = H.bp[i - 1] or -INF
    if #S > 0 and lo < 0 then nonc = segStr(S, var, dt) .. " for " .. rangeStr(lo, math.min(H.bp[i] or INF, 0), var, dt) end
  end
  for _, im in ipairs(H.imps) do if im.t0 < 0 then nonc = impStr(im, var) end end
  if nonc then out:r("NOT causal: h ≠ 0 for " .. var .. " < 0  (" .. nonc .. ")")
  else out:r("CAUSAL: h" .. o .. var .. c .. " = 0 for " .. var .. " < 0") end
  out:h("Stability")
  local A, exact = absSig(H)
  out:t((dt and "Σ|h[n]|" or "∫|h(t)|dt") .. (exact and " = " or " ≤ ") .. (dt and "Σ " or "∫ ") .. A:compact(var))
  local ok, val = pcall(function() return A:integrate(-INF, INF) end)
  if ok then
    local pcs = {}
    for _, pc in ipairs(A:pieces()) do
      local piece = fromPieces(dt, {{lo = pc.lo, hi = pc.hi, S = pc.S}})
      local pv = piece:integrate(pc.lo, dt and (pc.hi == INF and INF or pc.hi - 1) or pc.hi)
      pcs[#pcs + 1] = fmtC(pv)
      local G = dt and antiDiffS(pc.S) or antiS(pc.S)
      local lo_s = (pc.lo == -INF) and "-∞" or fmtR(pc.lo)
      local hi_s = (pc.hi == INF) and "∞" or fmtR(dt and pc.hi - 1 or pc.hi)
      if dt then
        out:t("  Σ[n=" .. lo_s .. "→" .. hi_s .. "] " .. segStr(pc.S, var, dt) .. ((pc.hi == INF or pc.lo == -INF) and "  (geometric series)" or "") .. " = " .. fmtC(pv))
      else
        out:t("  ∫[" .. lo_s .. "→" .. hi_s .. "] " .. segStr(pc.S, var, dt) .. " dt = [" .. segStr(G, var, dt) .. "] from " .. lo_s .. " to " .. hi_s .. " = " .. fmtC(pv))
      end
    end
    for _, im in ipairs(A.imps) do pcs[#pcs + 1] = fmtC(im.c) end
    out:t("  = " .. table.concat(pcs, " + ") .. (exact and " = " or " → total ") .. fmtC(val) .. " < ∞")
    out:r("STABLE (absolutely " .. (dt and "summable" or "integrable") .. ")")
  else
    out:t("  = ∞ (" .. tostring(val) .. ")")
    out:r("UNSTABLE")
  end
  out:h("Memoryless")
  local isMem = (#H.imps == 1 and abs(H.imps[1].t0) < 1e-12 and H.imps[1].m == 0)
  local smoothZero = true
  for _, S in ipairs(H.seg) do if #S > 0 then smoothZero = false end end
  if dt then
    isMem = #H.bp == 2 and H.bp[1] == 0 and H.bp[2] == 1 and #H.seg[1] == 0 and #H.seg[3] == 0
    out:t(isMem and "h[n] = K·δ[n] → memoryless" or "h[n] ≠ K·δ[n] → has memory")
  else
    isMem = isMem and smoothZero
    out:t(isMem and "h(t) = K·δ(t) → memoryless" or "h(t) ≠ K·δ(t) → has memory")
  end
end
addTool({
  name = "LTI system from h / step resp.",
  help = "Give h(t)/h[n] (or the step response s). Gives causal, stable, memoryless, step/impulse response, output for x.",
  fields = {{"h =", "e^(-4t)u(t-2)"}, {"or step resp s =", ""}, {"input x (optional) =", ""}},
  run = function(v, ctx, out)
    local H, dt
    if v[1] ~= "" then
      local nh = P(ctx, v[1]); dt = isDTexpr(nh, ctx); H = toSig(ctx, nh, dt)
    elseif v[2] ~= "" then
      local ns = P(ctx, v[2]); dt = isDTexpr(ns, ctx)
      local S = toSig(ctx, ns, dt)
      out:h("Impulse response from step response")
      if dt then
        H = ssub(S, S:compose(1, -1))
        out:t("h[n] = s[n] - s[n-1]")
      else
        H = S:deriv()
        out:t("h(t) = ds(t)/dt  (jumps give δ terms; d/dt u(t) = δ(t))")
      end
    else error("enter h or s", 0) end
    local var = V(dt)
    out:h("Impulse response")
    sigLines(out, H, var, "h")
    out:plot({kind = "sig", sig = H, label = "h"})
    ltiProps(H, out)
    if v[1] ~= "" then
      out:h("Step response")
      out:t(dt and "s[n] = Σ[k=-∞→n] h[k]" or "s(t) = ∫[-∞→t] h(τ)dτ")
      local okS, S = pcall(H.running, H, false)
      if okS then out:r("s" .. (dt and "[n]" or "(t)") .. " = " .. S:compact(var))
      else out:t("diverges (h does not vanish as " .. var .. "→-∞)") end
    end
    if v[3] ~= "" then
      local nx = P(ctx, v[3])
      local X = toSig(ctx, nx, dt)
      out:h("Output y = x * h")
      convSteps(X, H, out)
    end
  end,
})

-- ===== Transform tools =====
local function fwdTool(name, mode, dtMode, help, def)
  addTool({
    name = name, help = help,
    fields = {{dtMode and "x[n] =" or "x(t) =", def}},
    run = function(v, ctx, out)
      local node = P(ctx, v[1])
      local var = V(dtMode)
      local okX, X = pcall(toSig, ctx, node, dtMode)
      if not okX then
        out:e("Not in closed-form class (" .. tostring(X) .. ")")
        out:t("Use table pairs, e.g. sin(Wt)/(πt) <-> 1 for |ω|<W;  shift t-t0 => e^(-jωt0);  scale x(at) => (1/|a|)X(jω/a)")
        if mode == "F" or mode == "D" then
          out:h("Numeric X")
          local f = numFn(ctx, node, var)
          local function Xn(w)
            if dtMode then
              local s = ZERO
              for n = -200, 200 do local xv = f(n); if xv then s = s + xv * cexp(cx(0, -w * n)) end end
              return s
            end
            return simpson(function(t) local xv = f(t); if not xv then return ZERO end; return xv * cexp(cx(0, -w * t)) end, -40, 40, 1600)
          end
          for _, w in ipairs(dtMode and {0, pi / 4, pi / 2, pi} or {0, 0.5, 1, 1.5, 2, 3}) do
            out:t((dtMode and "Ω=" or "ω=") .. fmtR(w) .. ":  X ≈ " .. fmtC(Xn(w)) )
          end
          out:plot({kind = "freq", f = Xn, label = "|X| (numeric)", mode = mode, N = 60})
        end
        return
      end
      out:h("Signal")
      sigLines(out, X, var)
      out:plot({kind = "sig", sig = X, label = "x"})
      out:h("Steps")
      local T, I, roc, steps
      if dtMode then
        out:t(mode == "D" and "X(e^(jΩ)) = Σ x[n]e^(-jΩn)" or "X(z) = Σ x[n]z^(-n)")
        T, I, roc, steps = dtTransform(X, mode)
      else
        out:t(mode == "F" and "X(jω) = ∫ x(t)e^(-jωt) dt" or "X(s) = ∫ x(t)e^(-st) dt")
        T, I, roc, steps = ctTransform(X, mode)
      end
      local fl0 = dtMode and X:dtSupport()
      if not (fl0 and select(2, X:dtSupport()) - fl0 <= 40) then out:t("Split x into exponential pieces:") end
      for _, pc in ipairs((fl0 and select(2, X:dtSupport()) - fl0 <= 40) and {} or X:monoWindows()) do
        for _, m in ipairs(pc.S) do out:t("  " .. pieceStr(m.c, m.k, m.s, pc.lo, pc.hi, var, dtMode)) end
      end
      for _, s in ipairs(steps) do out:t("  " .. s) end
      out:h("Result")
      local lhs = ({F = "X(jω)", L = "X(s)", D = "X(e^(jΩ))", Z = "X(z)"})[mode]
      local main = ftSumStr(T, I, mode, dtMode)
      if dtMode and mode == "Z" then main = main:gsub("Ω", "") end
      out:r(lhs .. " = " .. main)
      if dtMode then
        local ok, N, D, mp = pcall(combineDT, T)
        if ok and #T > 1 and #I == 0 and ratStr(N, D, mode, mp) ~= main then out:r("Combined: " .. lhs .. " = " .. ratStr(N, D, mode, mp)) end
      else
        local groups = combineFT(T)
        if #T > 1 and #groups == 1 then
          local parts = {}
          for _, g in ipairs(groups) do
            local ds = delayStr(g.d, mode)
            local rs = ratStr(g.N, g.D, mode)
            parts[#parts + 1] = (ds ~= "" and (ds .. "·") or "") .. (ds ~= "" and ("[" .. rs .. "]") or rs)
          end
          out:r("Combined: " .. joinsum(parts))
        end
      end
      if mode == "L" then
        out:r("ROC: " .. (roc.lo > -INF and (fmtR(roc.lo) .. " < ") or "") .. "Re{s}" .. (roc.hi < INF and (" < " .. fmtR(roc.hi)) or "") .. ((roc.lo == -INF and roc.hi == INF) and " (entire plane)" or ""))
        if roc.lo >= roc.hi then out:e("ROC is empty → Laplace transform does not exist") end
      elseif mode == "Z" then
        if roc.lo == 0 and roc.hi == INF then
          local pos, neg = false, false
          for _, t in ipairs(T) do if t.d > 0 then pos = true elseif t.d < 0 then neg = true end end
          local ex = {}
          if pos then ex[#ex + 1] = "z = 0" end
          if neg then ex[#ex + 1] = "z = ∞" end
          out:r("ROC: all z" .. ((#ex > 0) and (" except " .. table.concat(ex, " and ")) or "") .. "  (finite-length signal)")
        else
          out:r("ROC: " .. ((roc.lo > 0) and (fmtR(roc.lo) .. " < ") or "") .. "|z|" .. (roc.hi < INF and (" < " .. fmtR(roc.hi)) or ""))
        end
      end
      checkForward(X, mode, T, I, roc, out)
      if mode == "F" then
        out:plot({kind = "freq", f = function(w) return ftEval(T, cx(0, w)) end, label = "|X(jω)|", mode = "F", imps = I})
      elseif mode == "D" then
        out:plot({kind = "freq", f = function(w) return dtEval(T, w) end, label = "|X(e^jΩ)|", mode = "D", imps = I})
      end
    end,
  })
end
fwdTool("Fourier transform (CT)", "F", false, "X(jω) from x(t) via analysis integral", "e^(-2t)u(t)")
fwdTool("DTFT", "D", true, "X(e^jΩ) from x[n]", "(1/2)^n u[n]")
fwdTool("Laplace transform", "L", false, "Bilateral X(s) with ROC", "e^(-t)u(t)")
fwdTool("Z-transform", "Z", true, "X(z) with ROC", "(1/2)^n u[n]")

local function invTool(name, mode, help, def, rocLabel)
  local fields = {{({F = "X(jω) =", L = "X(s) =", D = "X(e^jΩ) =", Z = "X(z) ="})[mode], def}}
  if rocLabel then fields[2] = {rocLabel, ""} end
  addTool({
    name = name, help = help, fields = fields,
    run = function(v, ctx, out)
      local roc = nil
      if rocLabel and v[2] ~= "" then
        -- forms: ">a", "<b", "a<x<b"
        local s = preprocess(v[2]):gsub("%s", ""):gsub("Re%(s%)", "x"):gsub("re%(s%)", "x"):gsub("|z|", "x"):gsub("abs%(z%)", "x"):gsub("[sz]", "x")
        local a, b = s:match("^(.-)<x<(.-)$")
        local function ev(e) return realof(evalN(parse(e), ctx.env)) end
        if a then roc = {lo = ev(a), hi = ev(b)}
        else
          local g = s:match("^x>(.+)$") or s:match("^>(.+)$")
          local l = s:match("^x<(.+)$") or s:match("^<(.+)$")
          if g then roc = {lo = ev(g), hi = INF} elseif l then roc = {lo = -INF, hi = ev(l)} end
        end
      end
      out:h("Steps")
      local sig, syn, synT = inverseXform(v[1], mode, ctx.env, roc, out)
      local dt = (mode == "D" or mode == "Z")
      local var = V(dt)
      out:h("Result")
      local s = sig:isZero() and "" or sig:compact(var)
      if syn then s = (s == "") and syn or (s .. " + " .. syn) end
      if s == "" then s = "0" end
      out:r((dt and "x[n]" or "x(t)") .. " = " .. s)
      if not sig:isZero() then
        local pw = sig:piecewise(var)
        if #pw > 1 then for _, l in ipairs(pw) do out:t("   " .. l.expr .. ",   " .. l.range) end end
      end
      if synT and #synT > 0 then
        out:t("(at " .. var .. " = shift point use the limit; e.g. integrate X directly)")
      end
      checkInverse(v[1], mode, ctx.env, sig, synT and #synT > 0, roc, out)
      local f = function(t)
        local val = sig:eval(t)
        if synT and #synT > 0 then val = val + synthEval(synT, t) end
        return val
      end
      if synT and #synT > 0 then out:plot({kind = "fn", f = f, dt = dt, label = "x", sig = sig})
      else out:plot({kind = "sig", sig = sig, label = "x"}) end
    end,
  })
end
invTool("Inverse Fourier (CT)", "F", "Use jw or w. Supports rational, e^(-jw t0), δ(w-w0), windows like (abs(w)<1).", "1/((jw+1)(jw+2))")
invTool("Inverse DTFT", "D", "Use W for Ω: e^(-jW), δ(W-W0), (abs(W)<pi/4)...", "1/(1-1/2e^(-jW))")
invTool("Inverse Laplace", "L", "Rational in s, e^(-as). ROC like >-1, <2, -1<x<2 (default causal).", "1/((s+1)(s+2))", "ROC (e.g. >-1):")
invTool("Inverse Z", "Z", "Rational in z or z^-1. ROC like >0.5 (default causal).", "z/(z-1/2)", "ROC (e.g. >0.5):")
