-- ===================== PARSER =====================
local FUNCS = {
  sin = 1, cos = 1, tan = 1, sec = 1, csc = 1, cot = 1, asin = 1, acos = 1, atan = 1, arcsin = 1, arccos = 1, arctan = 1,
  sinh = 1, cosh = 1, tanh = 1, exp = 1, ln = 1, log = 1, sqrt = 1, abs = 1, u = 1, step = 1, delta = 1,
  dirac = 1, imp = 1, rect = 1, tri = 1, sinc = 1, sgn = 1, sign = 1, re = 1, real = 1, im = 1, imag = 1,
  conj = 1, arg = 1, angle = 1, int = 1, integral = 1, sum = 1, der = 1, der2 = 1, conv = 1, floor = 1, ceil = 1,
  round = 1, max = 1, min = 1, mod = 1, ramp = 1, r = 0,
}
local NAMES = { "arcsin", "arccos", "arctan", "integral", "lambda", "alpha", "omega", "Omega", "theta", "sigma",
  "delta", "dirac", "gamma", "beta", "sinh", "cosh", "tanh", "asin", "acos", "atan", "sqrt", "rect", "sinc",
  "sign", "real", "imag", "conj", "step", "ramp", "angle", "floor", "ceil", "round", "der2", "conv",
  "sin", "cos", "tan", "sec", "csc", "cot", "exp", "abs", "sgn", "tri", "imp", "int", "sum", "der", "arg",
  "max", "min", "mod", "tau", "phi", "inf", "mu", "pi", "ln", "log", "re", "im", "oo" }

local UREP = {
  {"−", "-"}, {"⁻", "-"}, {"‐", "-"}, {"·", "*"}, {"×", "*"}, {"∙", "*"}, {"÷", "/"}, {"ℯ", "e"}, {"𝑖", "j"},
  {"²", "^2"}, {"³", "^3"}, {"√", "sqrt"}, {"π", "pi"}, {"∞", "inf"}, {"≤", "<="}, {"≥", ">="},
  {"∫", "int"}, {"Σ", "sum"}, {"∑", "sum"}, {"δ", "delta"}, {"ω", "w"}, {"Ω", "W"}, {"τ", "tau"},
  {"α", "alpha"}, {"β", "beta"}, {"γ", "gamma"}, {"θ", "theta"}, {"φ", "phi"}, {"σ", "sigma"},
  {"λ", "lambda"}, {"μ", "mu"}, {"ᴇ", "E"}, {"–", "-"}, {"\226\128\178", "'"}, {"’", "'"},
}
-- display names for identifiers
local DISP = { pi = "π", w = "ω", W = "Ω", tau = "τ", alpha = "α", beta = "β", gamma = "γ", theta = "θ",
  phi = "φ", sigma = "σ", lambda = "λ", mu = "μ", omega = "ω", Omega = "Ω", inf = "∞", delta = "δ", oo = "∞" }

local function preprocess(s)
  for _, p in ipairs(UREP) do s = s:gsub(p[1], p[2]) end
  s = s:gsub("d%^2(%a)/dt%^2", "%1''"):gsub("d%^2(%a)/dn%^2", "%1''")
  s = s:gsub("d(%a)/dt", "%1'")
  s = s:gsub("d%^2/dt%^2", "der2 "):gsub("d/dt", "der ")
  return s
end

local function lex(src)
  local s = preprocess(src)
  local toks = {}
  local i, n = 1, #s
  local function push(t, v) toks[#toks + 1] = {t = t, v = v} end
  while i <= n do
    local c = s:sub(i, i)
    if c:match("%s") then i = i + 1
    elseif c:match("[%d%.]") then
      local num = s:match("^%d*%.?%d*", i)
      local j = i + #num
      local ex = s:match("^E%-?%d+", j)
      if ex then num = num .. "e" .. ex:sub(2); j = j + #ex end
      local v = tonumber(num)
      if not v then error("bad number '" .. num .. "'", 0) end
      push("num", v); i = j
    elseif c:match("[%a_]") then
      local run = s:match("^[%a_]+", i)
      local j = i + #run
      local dg0 = s:match("^%d+", j)
      if dg0 and FUNCS[run .. dg0] then run = run .. dg0; j = j + #dg0 end
      -- split run into known names / single letters
      local pieces = {}
      local p = 1
      while p <= #run do
        local matched = nil
        for _, nm in ipairs(NAMES) do
          if run:sub(p, p + #nm - 1) == nm then matched = nm; break end
        end
        if matched and not (#run == p + #matched - 1 + 0 and false) then
          pieces[#pieces + 1] = matched; p = p + #matched
        else
          pieces[#pieces + 1] = run:sub(p, p); p = p + 1
        end
      end
      -- digits attach to final single letter (except j,i,e)
      local digs = s:match("^%d+", j)
      if digs and #pieces[#pieces] == 1 and not pieces[#pieces]:match("[jie]") then
        pieces[#pieces] = pieces[#pieces] .. digs; j = j + #digs
      end
      -- primes
      local pr = s:match("^'+", j)
      if pr then pieces[#pieces] = pieces[#pieces] .. pr; j = j + #pr end
      for _, pc in ipairs(pieces) do push("id", pc) end
      i = j
    else
      local two = s:sub(i, i + 1)
      if two == "<=" or two == ">=" or two == "!=" or two == "==" then push("op", two); i = i + 2
      elseif c:match("[%+%-%*/%^%(%)%[%],|=<>{};!]") then push("op", c); i = i + 1
      else error("unexpected character '" .. c .. "'", 0) end
    end
  end
  push("eof", nil)
  return toks
end

-- AST constructors
local function N(k, t) t.k = k; return t end
local function NUM(v) return N("num", {v = tocx(v)}) end

local Parser = {}
Parser.__index = Parser

local function parse(src, opts)
  opts = opts or {}
  local P = setmetatable({toks = lex(src), pos = 1, fn = opts.fnames or {}, absd = 0}, Parser)
  local e = P:expr()
  if P:peek().t ~= "eof" then error("unexpected '" .. tostring(P:peek().v) .. "'", 0) end
  return e
end

function Parser:peek(o) return self.toks[self.pos + (o or 0)] end
function Parser:next() local t = self.toks[self.pos]; self.pos = self.pos + 1; return t end
function Parser:isop(v) local t = self:peek(); return t.t == "op" and t.v == v end
function Parser:expect(v)
  local t = self:next()
  if t.t ~= "op" or t.v ~= v then error("expected '" .. v .. "'", 0) end
end

function Parser:expr() return self:cmp() end
function Parser:cmp()
  local a = self:add()
  local res = nil
  while true do
    local t = self:peek()
    if t.t == "op" and (t.v == "<" or t.v == ">" or t.v == "<=" or t.v == ">=") then
      self:next()
      local b = self:add()
      local c = N("cmp", {op = t.v, a = a, b = b})
      res = res and N("mul", {a = res, b = c}) or c
      a = b
    else break end
  end
  return res or a
end
function Parser:add()
  local a = self:mul()
  while true do
    if self:isop("+") then self:next(); a = N("add", {a = a, b = self:mul()})
    elseif self:isop("-") then self:next(); a = N("sub", {a = a, b = self:mul()})
    else break end
  end
  return a
end
function Parser:startsPrimary()
  local t = self:peek()
  if t.t == "num" or t.t == "id" then return true end
  if t.t == "op" and (t.v == "(" or t.v == "[" or t.v == "{") then return true end
  if t.t == "op" and t.v == "|" and self.absd == 0 then return true end
  return false
end
function Parser:mul()
  local a = self:unary()
  while true do
    if self:isop("*") then self:next(); a = N("mul", {a = a, b = self:unary()})
    elseif self:isop("/") then self:next(); a = N("div", {a = a, b = self:unary()})
    elseif self:startsPrimary() then a = N("mul", {a = a, b = self:pow(), imp = true})
    else break end
  end
  return a
end
function Parser:unary()
  if self:isop("-") then self:next(); return N("neg", {a = self:unary()}) end
  if self:isop("+") then self:next(); return self:unary() end
  return self:pow()
end
function Parser:pow()
  local b = self:postfix()
  if self:isop("^") then
    self:next()
    local e
    if self:isop("-") then self:next(); e = N("neg", {a = self:pow()}) else e = self:pow() end
    return N("pow", {a = b, b = e})
  end
  return b
end
function Parser:postfix()
  local a = self:primary()
  while self:isop("!") do self:next(); a = N("call", {name = "fact", args = {a}}) end
  return a
end
function Parser:args(close)
  local args = {}
  if self:isop(close) then self:next(); return args end
  local saved = self.absd; self.absd = 0
  while true do
    args[#args + 1] = self:expr()
    if self:isop(",") then self:next() else break end
  end
  self.absd = saved
  self:expect(close)
  return args
end
function Parser:primary()
  local t = self:next()
  if t.t == "num" then return NUM(t.v) end
  if t.t == "op" then
    if t.v == "(" or t.v == "{" then
      local saved = self.absd; self.absd = 0
      local e = self:expr()
      self.absd = saved
      self:expect(t.v == "(" and ")" or "}")
      return e
    elseif t.v == "[" then
      local saved = self.absd; self.absd = 0
      local e = self:expr(); self.absd = saved
      self:expect("]"); return e
    elseif t.v == "|" then
      self.absd = self.absd + 1
      local e = self:expr()
      self.absd = self.absd - 1
      self:expect("|")
      return N("call", {name = "abs", args = {e}})
    end
    error("unexpected '" .. t.v .. "'", 0)
  end
  if t.t == "id" then
    local name = t.v
    local base = name:gsub("'+$", "")
    local primes = #name - #base
    if base == "e" and self:isop("^") then
      self:next()
      local e
      if self:isop("-") then self:next(); e = N("neg", {a = self:pow()}) else e = self:pow() end
      return N("call", {name = "exp", args = {e}})
    end
    local isfun = (FUNCS[base] == 1) or self.fn[base] or (base == "d" and (self:isop("(") or self:isop("[")))
    if (base == "der" or base == "der2") and not self:isop("(") then
      -- der x(t)
      local a = self:pow()
      return N("call", {name = base, args = {a}})
    end
    if isfun and (self:isop("(") or self:isop("[")) then
      local open = self:next().v
      local args = self:args(open == "(" and ")" or "]")
      if base == "d" or base == "dirac" or base == "imp" then base = "delta" end
      if base == "step" then base = "u" end
      local node = N("call", {name = base, args = args, br = open})
      if primes > 0 then
        for _ = 1, primes do node = N("call", {name = "der", args = {node}}) end
      end
      return node
    end
    if self.fn[base] and primes > 0 then
      return N("var", {name = name})
    end
    if FUNCS[base] == 1 and base ~= "u" and base ~= "re" and base ~= "im" then
      -- function without parentheses: sin 2t  -> apply to next factor
      if self:startsPrimary() then
        local a = self:pow()
        return N("call", {name = base, args = {a}})
      end
    end
    if base == "inf" or base == "oo" then return N("num", {v = cx(INF, 0), inf = 1}) end
    return N("var", {name = name})
  end
  error("unexpected end of input", 0)
end

-- split "lhs = rhs" at top-level "="
local function splitEq(s)
  local depth = 0
  for i = 1, #s do
    local c = s:sub(i, i)
    if c == "(" or c == "[" then depth = depth + 1 elseif c == ")" or c == "]" then depth = depth - 1
    elseif c == "=" and depth == 0 then
      local pc, nc = s:sub(i - 1, i - 1), s:sub(i + 1, i + 1)
      if pc ~= "<" and pc ~= ">" and pc ~= "!" and nc ~= "=" then return s:sub(1, i - 1), s:sub(i + 1) end
    end
  end
  return nil, s
end

-- ---------- AST printer ----------
local PREC = {cmp = 1, add = 2, sub = 2, mul = 3, div = 3, neg = 4, pow = 5}
local tostr
local function wrap(node, minp, s)
  local p = PREC[node.k] or 9
  if p < minp then return "(" .. s .. ")" end
  return s
end
local function dname(n) return DISP[n] or n end
tostr = function(node)
  local k = node.k
  if k == "num" then
    if node.inf then return node.v.re > 0 and "∞" or "-∞" end
    if node.v.re == INF then return "∞" elseif node.v.re == -INF then return "-∞" end
    return fmtC(node.v)
  elseif k == "var" then return dname(node.name)
  elseif k == "add" then return tostr(node.a) .. " + " .. wrap(node.b, 2, tostr(node.b))
  elseif k == "sub" then return tostr(node.a) .. " - " .. wrap(node.b, 3, tostr(node.b))
  elseif k == "mul" then
    local a, b = wrap(node.a, 3, tostr(node.a)), wrap(node.b, 3, tostr(node.b))
    if node.a.k == "num" and a:find("[/ ]") then a = "(" .. a .. ")" end
    if b:sub(1, 1) == "-" then b = "(" .. b .. ")" end
    if node.a.k == "var" and (node.a.name == "j" or node.a.name == "i") and node.b.k == "num" then return "j" .. b end
    if node.a.k == "num" and (node.b.k ~= "num") then return a .. b end
    if node.b.k == "num" then return a .. "·" .. b end
    if b:sub(1, 1) == "-" then b = "(" .. b .. ")" end
    if a:find("[%a]$") and b:find("^[%a]") and not (node.a.k == "var" and (node.a.name == "j" or node.a.name == "i")) then return a .. "·" .. b end
    return a .. b
  elseif k == "div" then
    return wrap(node.a, 3, tostr(node.a)) .. "/" .. wrap(node.b, 4, tostr(node.b))
  elseif k == "neg" then
    local s = wrap(node.a, 4, tostr(node.a))
    if s:sub(1, 1) == "-" then s = "(" .. s .. ")" end
    return "-" .. s
  elseif k == "pow" then
    local b = tostr(node.b)
    if not (node.b.k == "num" and node.b.v.re >= 0 or node.b.k == "var") then b = "(" .. b .. ")" end
    return wrap(node.a, 6, tostr(node.a)) .. "^" .. b
  elseif k == "cmp" then
    local op = ({["<="] = "≤", [">="] = "≥"})[node.op] or node.op
    return "[" .. tostr(node.a) .. op .. tostr(node.b) .. "]"
  elseif k == "call" then
    local nm = node.name
    if nm == "exp" then
      local a = tostr(node.args[1])
      if node.args[1].k == "var" or (node.args[1].k == "num" and node.args[1].v.re >= 0) then return "e^" .. a end
      return "e^(" .. a .. ")"
    end
    if nm == "abs" then return "|" .. tostr(node.args[1]) .. "|" end
    if nm == "fact" then return wrap(node.args[1], 9, tostr(node.args[1])) .. "!" end
    if nm == "int" or nm == "integral" then
      local a = node.args
      if #a == 4 then return "∫[" .. tostr(a[3]) .. "→" .. tostr(a[4]) .. "] " .. wrap(a[1], 3, tostr(a[1])) .. " d" .. tostr(a[2]) end
    end
    if nm == "sum" then
      local a = node.args
      if #a == 4 then return "Σ[" .. tostr(a[2]) .. "=" .. tostr(a[3]) .. "→" .. tostr(a[4]) .. "] " .. wrap(a[1], 3, tostr(a[1])) end
    end
    if nm == "der" then return "d/dt{" .. tostr(node.args[1]) .. "}" end
    if nm == "der2" then return "d²/dt²{" .. tostr(node.args[1]) .. "}" end
    local parts = {}
    for i, a in ipairs(node.args) do parts[i] = tostr(a) end
    local o, c = "(", ")"
    if node.br == "[" then o, c = "[", "]" end
    return dname(nm) .. o .. table.concat(parts, ", ") .. c
  end
  return "?"
end

-- substitute variable name -> node (deep copy)
local function subst(node, name, rep)
  if node.k == "var" then
    if node.name == name then return rep end
    return node
  end
  local r = {}
  for kk, v in pairs(node) do
    if type(v) == "table" and v.k then r[kk] = subst(v, name, rep)
    elseif kk == "args" then
      local a = {}
      for i, x in ipairs(v) do a[i] = subst(x, name, rep) end
      r[kk] = a
    else r[kk] = v end
  end
  return r
end
-- replace calls name(arg) with f(argnode) -> node
local function mapcalls(node, name, f)
  if node.k == "call" and node.name == name then
    local a = {}
    for i, x in ipairs(node.args) do a[i] = mapcalls(x, name, f) end
    return f(a, node)
  end
  if node.k == "num" or node.k == "var" then return node end
  local r = {}
  for kk, v in pairs(node) do
    if type(v) == "table" and v.k then r[kk] = mapcalls(v, name, f)
    elseif kk == "args" then
      local a = {}
      for i, x in ipairs(v) do a[i] = mapcalls(x, name, f) end
      r[kk] = a
    else r[kk] = v end
  end
  return r
end
local function contains(node, pred)
  if pred(node) then return true end
  for kk, v in pairs(node) do
    if type(v) == "table" and v.k then if contains(v, pred) then return true end
    elseif kk == "args" then for _, x in ipairs(v) do if contains(x, pred) then return true end end end
  end
  return false
end
local function hasvar(node, name) return contains(node, function(n) return n.k == "var" and n.name == name end) end
