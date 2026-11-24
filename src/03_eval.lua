-- ===================== NUMERIC EVALUATOR =====================
local evalN
local NUMINF = 30   -- truncation for infinite numeric integrals
local SUMINF = 60

local function realof(z, what)
  z = tocx(z)
  if abs(z.im) > 1e-9 * math.max(1, abs(z.re)) then error((what or "value") .. " must be real", 0) end
  return z.re
end

local function simpson(f, a, b, n)
  n = n or 200
  if n % 2 == 1 then n = n + 1 end
  local h = (b - a) / n
  local s = f(a) + f(b)
  for i = 1, n - 1 do s = s + f(a + i * h) * ((i % 2 == 1) and 4 or 2) end
  return s * (h / 3)
end

local BUILTIN = {}
BUILTIN.sin = function(z) return csin(z) end
BUILTIN.cos = function(z) return ccos(z) end
BUILTIN.tan = function(z) return csin(z) / ccos(z) end
BUILTIN.sec = function(z) return ONE / ccos(z) end
BUILTIN.csc = function(z) return ONE / csin(z) end
BUILTIN.cot = function(z) return ccos(z) / csin(z) end
BUILTIN.asin = function(z) return cx(math.asin(realof(z))) end
BUILTIN.acos = function(z) return cx(math.acos(realof(z))) end
BUILTIN.atan = function(z) return cx(math.atan(realof(z))) end
BUILTIN.arcsin, BUILTIN.arccos, BUILTIN.arctan = BUILTIN.asin, BUILTIN.acos, BUILTIN.atan
BUILTIN.sinh = function(z) return (cexp(z) - cexp(-z)) / 2 end
BUILTIN.cosh = function(z) return (cexp(z) + cexp(-z)) / 2 end
BUILTIN.tanh = function(z) return (cexp(z) - cexp(-z)) / (cexp(z) + cexp(-z)) end
BUILTIN.exp = cexp
BUILTIN.ln = clog
BUILTIN.log = function(z, b) if b then return clog(z) / clog(b) end return clog(z) / log(10) end
BUILTIN.sqrt = csqrt
BUILTIN.abs = function(z) return cx(cabs(z)) end
BUILTIN.u = function(z) return cx(realof(z, "u() argument") >= 0 and 1 or 0) end
BUILTIN.delta = function(z) local x = realof(z); return cx(0) end
BUILTIN.rect = function(z) local x = abs(realof(z)); return cx(x < 0.5 and 1 or (x == 0.5 and 1 or 0)) end
BUILTIN.tri = function(z) local x = abs(realof(z)); return cx(x < 1 and 1 - x or 0) end
BUILTIN.ramp = function(z) local x = realof(z); return cx(x > 0 and x or 0) end
BUILTIN.sinc = function(z) local x = realof(z); if abs(x) < 1e-12 then return cx(1) end return cx(sin(pi * x) / (pi * x)) end
BUILTIN.sgn = function(z) local x = realof(z); return cx(x > 0 and 1 or (x < 0 and -1 or 0)) end
BUILTIN.sign = BUILTIN.sgn
BUILTIN.re = function(z) return cx(tocx(z).re) end
BUILTIN.real = BUILTIN.re
BUILTIN.im = function(z) return cx(tocx(z).im) end
BUILTIN.imag = BUILTIN.im
BUILTIN.conj = cconj
BUILTIN.arg = function(z) return cx(carg(z)) end
BUILTIN.angle = BUILTIN.arg
BUILTIN.floor = function(z) return cx(floor(realof(z))) end
BUILTIN.ceil = function(z) return cx(math.ceil(realof(z))) end
BUILTIN.round = function(z) return cx(round(realof(z))) end
BUILTIN.max = function(a, b) return cx(math.max(realof(a), realof(b))) end
BUILTIN.min = function(a, b) return cx(math.min(realof(a), realof(b))) end
BUILTIN.mod = function(a, b) return cx(realof(a) % realof(b)) end
BUILTIN.fact = function(a) return cx(fact(round(realof(a)))) end

local CONST = { pi = cx(pi), e = cx(math.exp(1)), j = J, i = J }

local function lookupVar(env, name)
  local v = env.vars and env.vars[name]
  if v then return v end
  v = CONST[name]
  if v then return v end
  error("unknown symbol '" .. (DISP[name] or name) .. "' (set a value in Params)", 0)
end

evalN = function(node, env)
  local k = node.k
  if k == "num" then return node.v
  elseif k == "var" then return lookupVar(env, node.name)
  elseif k == "add" then return evalN(node.a, env) + evalN(node.b, env)
  elseif k == "sub" then return evalN(node.a, env) - evalN(node.b, env)
  elseif k == "mul" then
    local a = evalN(node.a, env)
    if isz(a, 0) then
      -- still evaluate b to propagate errors? skip for speed
      return ZERO
    end
    return a * evalN(node.b, env)
  elseif k == "div" then return evalN(node.a, env) / evalN(node.b, env)
  elseif k == "neg" then return -evalN(node.a, env)
  elseif k == "pow" then return cpow(evalN(node.a, env), evalN(node.b, env))
  elseif k == "cmp" then
    local a, b = realof(evalN(node.a, env)), realof(evalN(node.b, env))
    local op = node.op
    local r = (op == "<" and a < b) or (op == ">" and a > b) or (op == "<=" and a <= b) or (op == ">=" and a >= b)
    return cx(r and 1 or 0)
  elseif k == "call" then
    local nm = node.name
    local A = node.args
    if nm == "int" or nm == "integral" then
      if #A ~= 4 then error("int(f, var, lower, upper)", 0) end
      local var = A[2].name
      local lo, hi = realof(evalN(A[3], env)), realof(evalN(A[4], env))
      local sg = 1
      if lo > hi then lo, hi = hi, lo; sg = -1 end
      if lo == -INF then lo = math.min(-NUMINF, hi - NUMINF) end
      if hi == INF then hi = math.max(NUMINF, lo + NUMINF) end
      if hi - lo < 1e-14 then return ZERO end
      local saved = env.vars[var]
      local f = function(x) env.vars[var] = cx(x); return evalN(A[1], env) end
      local n = math.min(env.intCap or 1200, math.max(40, math.ceil((hi - lo) * 30)))
      local r = simpson(f, lo, hi, n)
      env.vars[var] = saved
      return r * sg
    elseif nm == "sum" then
      if #A ~= 4 then error("sum(f, k, lower, upper)", 0) end
      local var = A[2].name
      local lo, hi = realof(evalN(A[3], env)), realof(evalN(A[4], env))
      local SI = env.sumInf or SUMINF
      if lo == -INF then lo = math.min(-SI, hi - SI) end
      if hi == INF then hi = math.max(SI, lo + SI) end
      local saved = env.vars[var]
      local s = ZERO
      for kk = math.ceil(lo - 1e-9), floor(hi + 1e-9) do
        env.vars[var] = cx(kk); s = s + evalN(A[1], env)
      end
      env.vars[var] = saved
      return s
    elseif nm == "der" or nm == "der2" then
      local var = env.tvar or "t"
      local t0 = lookupVar(env, var)
      local h = env.dh or 1e-4
      local f = function(x) env.vars[var] = cx(x); return evalN(A[1], env) end
      local r
      if nm == "der" then r = (f(t0.re + h) - f(t0.re - h)) / (2 * h)
      else r = (f(t0.re + h) - f(t0.re) * 2 + f(t0.re - h)) / (h * h) end
      env.vars[var] = t0
      return r
    elseif nm == "conv" then
      error("conv() only in exact mode", 0)
    end
    local uf = env.fns and env.fns[nm]
    if uf then
      local argv = {}
      for i2, a in ipairs(A) do argv[i2] = evalN(a, env) end
      if type(uf) == "function" then return uf(table.unpack and table.unpack(argv) or unpack(argv)) end
      local saved = env.vars[uf.param]
      env.vars[uf.param] = argv[1]
      local ok, r = pcall(evalN, uf.body, env)
      env.vars[uf.param] = saved
      if not ok then error(r, 0) end
      return r
    end
    local bf = BUILTIN[nm]
    if not bf then error("unknown function '" .. nm .. "'", 0) end
    local a1 = A[1] and evalN(A[1], env)
    local a2 = A[2] and evalN(A[2], env)
    return bf(a1, a2)
  end
  error("cannot evaluate node " .. tostring(k), 0)
end

-- make evaluation env
local function mkenv(params, fns)
  local env = {vars = {}, fns = fns or {}}
  for kk, v in pairs(params or {}) do env.vars[kk] = v end
  return env
end

-- test whether node is linear in var: returns alpha, beta (reals) or nil
local function linearOf(node, var, env)
  if contains(node, function(n) return n.k == "call" and (n.name == "delta" or n.name == "u" or n.name == "abs" or n.name == "int" or n.name == "sum") end) then
    -- still may be linear if function of no var; check var presence
    if hasvar(node, var) then return nil end
  end
  local saved = env.vars[var]
  local ok, res = pcall(function()
    local f = function(x) env.vars[var] = cx(x); return evalN(node, env) end
    local f0, f1, f2, f3 = f(0), f(1), f(2.7), f(-1.9)
    local b = f0; local a = f1 - f0
    if not (near(f2, b + a * 2.7, 1e-8) and near(f3, b - a * 1.9, 1e-8)) then return nil end
    return {a, b}
  end)
  env.vars[var] = saved
  if not ok or not res then return nil end
  return res[1], res[2]
end
local function isconst(node, var, env)
  if not hasvar(node, var) then
    local ok, v = pcall(evalN, node, env)
    if ok then return v end
  end
  return nil
end
