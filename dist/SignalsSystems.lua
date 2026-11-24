platform.apilevel = "2.0"
-- ===================== CORE: complex numbers & formatting =====================
abs, floor, sqrt, exp, log, sin, cos, atan2, pi = math.abs, math.floor, math.sqrt, math.exp, math.log, math.sin, math.cos, math.atan2, math.pi
INF = math.huge
EPS = 1e-9

C = {}
C.__index = C
function cx(re, im) return setmetatable({re = re, im = im or 0}, C) end
function tocx(a) if type(a) == "number" then return cx(a, 0) end return a end
C.__add = function(a, b) a, b = tocx(a), tocx(b); return cx(a.re + b.re, a.im + b.im) end
C.__sub = function(a, b) a, b = tocx(a), tocx(b); return cx(a.re - b.re, a.im - b.im) end
C.__mul = function(a, b) a, b = tocx(a), tocx(b); return cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re) end
C.__div = function(a, b)
  a, b = tocx(a), tocx(b)
  local d = b.re * b.re + b.im * b.im
  if d == 0 then error("division by zero", 0) end
  return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d)
end
C.__unm = function(a) return cx(-a.re, -a.im) end
C.__eq = function(a, b) return a.re == b.re and a.im == b.im end
function cabs(a) a = tocx(a); return sqrt(a.re * a.re + a.im * a.im) end
function carg(a) a = tocx(a); return atan2(a.im, a.re) end
function cexp(a) a = tocx(a); local m = exp(a.re); return cx(m * cos(a.im), m * sin(a.im)) end
function clog(a) a = tocx(a); return cx(log(cabs(a)), carg(a)) end
function cconj(a) a = tocx(a); return cx(a.re, -a.im) end
function cpow(a, b)
  a, b = tocx(a), tocx(b)
  if b.im == 0 and b.re == floor(b.re) and abs(b.re) < 64 then
    local n = b.re; local r = cx(1, 0); local base = a
    if n < 0 then base = cx(1, 0) / a; n = -n end
    for _ = 1, n do r = r * base end
    return r
  end
  if a.re == 0 and a.im == 0 then
    if b.re > 0 then return cx(0, 0) end
    error("0 to a non-positive power", 0)
  end
  return cexp(b * clog(a))
end
function csqrt(a) a = tocx(a); local r = cabs(a); local t = carg(a) / 2; local m = sqrt(r); return cx(m * cos(t), m * sin(t)) end
function csin(a) a = tocx(a); return (cexp(cx(-a.im, a.re)) - cexp(cx(a.im, -a.re))) / cx(0, 2) end
function ccos(a) a = tocx(a); return (cexp(cx(-a.im, a.re)) + cexp(cx(a.im, -a.re))) / 2 end
function isz(a, tol) tol = tol or EPS; a = tocx(a); return abs(a.re) <= tol and abs(a.im) <= tol end
function near(a, b, tol)
  a, b = tocx(a), tocx(b); tol = tol or 1e-7
  local d = cabs(a - b); local m = math.max(1, cabs(a), cabs(b))
  return d <= tol * m
end
function isreal(a, tol) a = tocx(a); return abs(a.im) <= (tol or 1e-9) * math.max(1, abs(a.re)) end
J = cx(0, 1)
ZERO, ONE = cx(0, 0), cx(1, 0)

-- clean tiny float noise
function cln(x)
  if abs(x) < 1e-11 then return 0 end
  local r = floor(x + 0.5)
  if abs(x - r) < 1e-9 * math.max(1, abs(x)) then return r end
  return x
end
function ccln(z) z = tocx(z); return cx(cln(z.re), cln(z.im)) end

-- ---------- number formatting ----------
function ratapprox(x, maxq, tol)
  -- continued fraction rational approx; returns p, q or nil
  tol = tol or 1e-9
  local sgn = x < 0 and -1 or 1; x = abs(x)
  local h0, h1, k0, k1 = 0, 1, 1, 0
  local y = x
  for _ = 1, 40 do
    local a = floor(y)
    local h2, k2 = a * h1 + h0, a * k1 + k0
    if k2 > maxq then return nil end
    h0, h1, k0, k1 = h1, h2, k1, k2
    if abs(x - h1 / k1) <= tol * math.max(1, x) then return sgn * h1, k1 end
    local fr = y - a
    if fr < 1e-15 then break end
    y = 1 / fr
  end
  if abs(x - h1 / k1) <= tol * math.max(1, x) and k1 <= maxq then return sgn * h1, k1 end
  return nil
end

function trimdec(s)
  if s:find("e") then return s end
  if s:find("%.") then s = s:gsub("0+$", ""):gsub("%.$", "") end
  return s
end
function decstr(x)
  local s
  if abs(x) >= 1e5 or (abs(x) < 1e-3 and x ~= 0) then s = string.format("%.4e", x)
  else s = trimdec(string.format("%.5g", x)) end
  return s
end

fracmode = true
function fracstr(p, q)
  if q == 1 then return tostring(p) end
  return p .. "/" .. q
end
SQ = {2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21}
-- format real number; returns string, and "kind"
function fmtR(x)
  if x ~= x then return "undef" end
  if x == INF then return "∞" elseif x == -INF then return "-∞" end
  x = cln(x)
  if x == floor(x) and abs(x) < 1e15 then return string.format("%d", x) end
  if not fracmode then return decstr(x) end
  local p, q = ratapprox(x, 64)
  if p and q > 16 and abs(x * 1000 - floor(x * 1000 + 0.5)) < 1e-6 then return trimdec(string.format("%.3f", x)) end
  if p then return fracstr(p, q) end
  p, q = ratapprox(x / pi, 48)
  if p then
    local ps = (p == 1) and "" or (p == -1 and "-" or tostring(p))
    if q == 1 then return ps .. "π" end
    return ps .. "π/" .. q
  end
  p, q = ratapprox(x / (pi * pi), 24)
  if p then return ((p == 1) and "" or (p == -1 and "-" or tostring(p))) .. "π²" .. ((q == 1) and "" or ("/" .. q)) end
  p, q = ratapprox(x * pi, 48)
  if p then
    local ps = (p == 1 and "1" or p == -1 and "-1" or tostring(p))
    if q == 1 then return ps .. "/π" end
    return ps .. "/(" .. q .. "π)"
  end
  for _, n in ipairs(SQ) do
    p, q = ratapprox(x / sqrt(n), 48)
    if p then
      local ps = (p == 1) and "" or (p == -1 and "-" or tostring(p))
      if q == 1 then return ps .. "√" .. n end
      return ps .. "√" .. n .. "/" .. q
    end
  end
  for k = 1, 12 do
    for _, sg in ipairs({-1, 1}) do
      local e = sg * k
      p, q = ratapprox(x / exp(e), 16, 1e-10)
      if p and abs(p) < 1000 then
        local es = "e^(" .. e .. ")"
        if e == 1 then es = "e" end
        local ps = (p == 1) and "" or (p == -1 and "-" or (p .. "·"))
        if q == 1 then return ps .. es end
        return ps .. es .. "/" .. q
      end
    end
  end
  local x3 = x * 1000
  if abs(x3 - floor(x3 + 0.5)) < 1e-6 then return trimdec(string.format("%.3f", x)) end
  p, q = ratapprox(x, 2000)
  if p and abs(p) < 1e6 then return fracstr(p, q) end
  return decstr(x)
end

-- does a formatted real need parentheses when multiplied / used as coefficient
function needpar(s)
  return s:find("[%+/]") ~= nil or s:find(".%-") ~= nil
end

function fmtImag(y)
  -- y real, returns string like "j", "-j", "j2", "-j(1/2)"
  local s = fmtR(abs(y))
  local sg = y < 0 and "-" or ""
  if s == "1" then return sg .. "j" end
  if s:find("^[%d%.]+$") or s:find("^[%d]*π$") or s:find("^√%d+$") then return sg .. "j" .. s end
  return sg .. "j(" .. s .. ")"
end

polarmode = false
function fmtC(z)
  z = tocx(z)
  local re, im = cln(z.re), cln(z.im)
  if abs(im) <= 1e-10 * math.max(1, abs(re)) then return fmtR(re) end
  if abs(re) <= 1e-10 * math.max(1, abs(im)) then return fmtImag(im) end
  local s = fmtR(re)
  if im < 0 then return s .. " - " .. fmtImag(-im) end
  return s .. " + " .. fmtImag(im)
end
function fmtPolar(z)
  z = tocx(z)
  local m = cabs(z)
  if m < 1e-12 then return "0" end
  local a = carg(z)
  local ms = fmtR(m)
  if abs(a) < 1e-12 then return ms end
  local sg = a < 0 and "-" or ""
  local as = fmtR(abs(a))
  if needpar(ms) then ms = "(" .. ms .. ")" end
  if ms == "1" then ms = "" end
  return ms .. "e^(" .. sg .. "j" .. (as:find("[%+%-]") and ("(" .. as .. ")") or as) .. ")"
end
-- coefficient * expression string, handles 1, -1, parentheses
function cmul(c, expr)
  c = tocx(c)
  if expr == nil or expr == "" or expr == "1" then return fmtC(c) end
  if near(c, ONE, 1e-10) then return expr end
  if near(c, -ONE, 1e-10) then return "-" .. expr end
  local s = fmtC(c)
  if isreal(c) then
    if s:find("^%-?%d+/%d+$") then
      -- p/q * expr -> (p/q)expr  shown as p·expr/q is nicer for single factor
      local p, q = s:match("^(%-?%d+)/(%d+)$")
      if p == "1" then return "(1/" .. q .. ")" .. expr end
      if p == "-1" then return "-(1/" .. q .. ")" .. expr end
      if p:sub(1, 1) == "-" then return "-(" .. p:sub(2) .. "/" .. q .. ")" .. expr end
      return "(" .. s .. ")" .. expr
    end
    if needpar(s) then return "(" .. s .. ")" .. expr end
    if expr:find("^%([%d%-]") or expr:find("^%d") then return s .. "·" .. expr end
    if s:find("[%a]$") or s:find("[π√]") then return s .. "·" .. expr end
    return s .. expr
  end
  local ps = fmtPolar(c)
  if not ps:find("%d%.%d") and ulen(ps) < ulen(s) + 2 then return ps .. "·" .. expr end
  return "(" .. s .. ")" .. expr
end
-- does string have a top-level + or - (outside parentheses)
function topSum(s)
  local d = 0
  for i = 1, #s do
    local c = s:sub(i, i)
    if c == "(" or c == "[" then d = d + 1 elseif c == ")" or c == "]" then d = d - 1
    elseif d == 0 and i > 1 and (c == "+" or c == "-") and s:sub(i - 1, i - 1) == " " then return true end
  end
  return false
end
-- join term strings into sum with proper signs
function joinsum(terms)
  if #terms == 0 then return "0" end
  local out = terms[1]
  for i = 2, #terms do
    local t = terms[i]
    if t:sub(1, 1) == "-" then out = out .. " - " .. t:sub(2) else out = out .. " + " .. t end
  end
  return out
end

-- ---------- misc helpers ----------
function gcdi(a, b) a, b = abs(a), abs(b); while b ~= 0 do a, b = b, a % b end; return a end
function lcmi(a, b) return a / gcdi(a, b) * b end
function fact(n) local r = 1; for i = 2, n do r = r * i end; return r end
function binom(n, k) if k < 0 or k > n then return 0 end; return fact(n) / (fact(k) * fact(n - k)) end
function copy(t) local r = {}; for k, v in pairs(t) do r[k] = v end; return r end
function isint(x) return abs(x - floor(x + 0.5)) < 1e-9 end
function round(x) return floor(x + 0.5) end

-- UTF-8 helpers
function uchars(s)
  local t = {}
  for ch in s:gmatch("[%z\1-\127\194-\244][\128-\191]*") do t[#t + 1] = ch end
  return t
end
function ulen(s) local _, n = s:gsub("[^\128-\191]", ""); return n end
-- ===================== PARSER =====================
FUNCS = {
  sin = 1, cos = 1, tan = 1, sec = 1, csc = 1, cot = 1, asin = 1, acos = 1, atan = 1, arcsin = 1, arccos = 1, arctan = 1,
  sinh = 1, cosh = 1, tanh = 1, exp = 1, ln = 1, log = 1, sqrt = 1, abs = 1, u = 1, step = 1, delta = 1,
  dirac = 1, imp = 1, rect = 1, tri = 1, sinc = 1, sgn = 1, sign = 1, re = 1, real = 1, im = 1, imag = 1,
  conj = 1, arg = 1, angle = 1, int = 1, integral = 1, sum = 1, der = 1, der2 = 1, conv = 1, floor = 1, ceil = 1,
  round = 1, max = 1, min = 1, mod = 1, ramp = 1, r = 0,
}
NAMES = { "arcsin", "arccos", "arctan", "integral", "lambda", "alpha", "omega", "Omega", "theta", "sigma",
  "delta", "dirac", "gamma", "beta", "sinh", "cosh", "tanh", "asin", "acos", "atan", "sqrt", "rect", "sinc",
  "sign", "real", "imag", "conj", "step", "ramp", "angle", "floor", "ceil", "round", "der2", "conv",
  "sin", "cos", "tan", "sec", "csc", "cot", "exp", "abs", "sgn", "tri", "imp", "int", "sum", "der", "arg",
  "max", "min", "mod", "tau", "phi", "inf", "mu", "pi", "ln", "log", "re", "im", "oo" }

UREP = {
  {"−", "-"}, {"⁻", "-"}, {"‐", "-"}, {"·", "*"}, {"×", "*"}, {"∙", "*"}, {"÷", "/"}, {"ℯ", "e"}, {"𝑖", "j"},
  {"²", "^2"}, {"³", "^3"}, {"√", "sqrt"}, {"π", "pi"}, {"∞", "inf"}, {"≤", "<="}, {"≥", ">="},
  {"∫", "int"}, {"Σ", "sum"}, {"∑", "sum"}, {"δ", "delta"}, {"ω", "w"}, {"Ω", "W"}, {"τ", "tau"},
  {"α", "alpha"}, {"β", "beta"}, {"γ", "gamma"}, {"θ", "theta"}, {"φ", "phi"}, {"σ", "sigma"},
  {"λ", "lambda"}, {"μ", "mu"}, {"ᴇ", "E"}, {"–", "-"}, {"\226\128\178", "'"}, {"’", "'"},
}
-- display names for identifiers
DISP = { pi = "π", w = "ω", W = "Ω", tau = "τ", alpha = "α", beta = "β", gamma = "γ", theta = "θ",
  phi = "φ", sigma = "σ", lambda = "λ", mu = "μ", omega = "ω", Omega = "Ω", inf = "∞", delta = "δ", oo = "∞" }

function preprocess(s)
  for _, p in ipairs(UREP) do s = s:gsub(p[1], p[2]) end
  s = s:gsub("d%^2(%a)/dt%^2", "%1''"):gsub("d%^2(%a)/dn%^2", "%1''")
  s = s:gsub("d(%a)/dt", "%1'")
  s = s:gsub("d%^2/dt%^2", "der2 "):gsub("d/dt", "der ")
  return s
end

function lex(src)
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
function N(k, t) t.k = k; return t end
function NUM(v) return N("num", {v = tocx(v)}) end

Parser = {}
Parser.__index = Parser

function parse(src, opts)
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
function splitEq(s)
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
PREC = {cmp = 1, add = 2, sub = 2, mul = 3, div = 3, neg = 4, pow = 5}
function wrap(node, minp, s)
  local p = PREC[node.k] or 9
  if p < minp then return "(" .. s .. ")" end
  return s
end
function dname(n) return DISP[n] or n end
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
function subst(node, name, rep)
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
function mapcalls(node, name, f)
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
function contains(node, pred)
  if pred(node) then return true end
  for kk, v in pairs(node) do
    if type(v) == "table" and v.k then if contains(v, pred) then return true end
    elseif kk == "args" then for _, x in ipairs(v) do if contains(x, pred) then return true end end end
  end
  return false
end
function hasvar(node, name) return contains(node, function(n) return n.k == "var" and n.name == name end) end
-- ===================== NUMERIC EVALUATOR =====================
NUMINF = 30   -- truncation for infinite numeric integrals
SUMINF = 60

function realof(z, what)
  z = tocx(z)
  if abs(z.im) > 1e-9 * math.max(1, abs(z.re)) then error((what or "value") .. " must be real", 0) end
  return z.re
end

function simpson(f, a, b, n)
  n = n or 200
  if n % 2 == 1 then n = n + 1 end
  local h = (b - a) / n
  local s = f(a) + f(b)
  for i = 1, n - 1 do s = s + f(a + i * h) * ((i % 2 == 1) and 4 or 2) end
  return s * (h / 3)
end

BUILTIN = {}
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

CONST = { pi = cx(pi), e = cx(math.exp(1)), j = J, i = J }

function lookupVar(env, name)
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
function mkenv(params, fns)
  local env = {vars = {}, fns = fns or {}}
  for kk, v in pairs(params or {}) do env.vars[kk] = v end
  return env
end

-- test whether node is linear in var: returns alpha, beta (reals) or nil
function linearOf(node, var, env)
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
function isconst(node, var, env)
  if not hasvar(node, var) then
    local ok, v = pcall(evalN, node, env)
    if ok then return v end
  end
  return nil
end
-- ===================== SIGNAL ENGINE =====================
-- Sig: piecewise exponential-polynomial signal
--  CT: segment i = interval (bp[i-1], bp[i]); monomial {c,k,s} = c t^k e^{s t}
--  DT: segment i = integers bp[i-1] <= n < bp[i]; monomial {c,k,s} = c n^k s^n  (s is the base r)
--  imps (CT only): {c, t0, m} = c * delta^(m)(t - t0)
Sig = {}
Sig.__index = Sig
TOL = 1e-9

function M(c, k, s) return {c = tocx(c), k = k or 0, s = tocx(s)} end
function newSig(dt) return setmetatable({dt = dt and true or false, bp = {}, seg = {{}}, imps = {}}, Sig) end

function monoEval(m, t, dt)
  local v = m.c
  if m.k > 0 then v = v * cx(t ^ m.k) end
  if dt then
    if not near(m.s, ONE, 1e-12) then v = v * cpow(m.s, t) end
  else
    if not isz(m.s, 0) then v = v * cexp(m.s * t) end
  end
  return v
end
function segEval(S, t, dt)
  local v = ZERO
  for _, m in ipairs(S) do v = v + monoEval(m, t, dt) end
  return v
end
function sameS(a, b, dt)
  local idA = dt and ONE or ZERO
  return near(a, b, 1e-9)
end
-- merge like monomials, drop zeros
function mergeS(S, dt)
  local out = {}
  for _, m in ipairs(S) do
    local found = false
    for _, o in ipairs(out) do
      if o.k == m.k and near(o.s, m.s, 1e-9) then o.c = o.c + m.c; found = true; break end
    end
    if not found then out[#out + 1] = M(m.c, m.k, m.s) end
  end
  local res = {}
  for _, m in ipairs(out) do
    if cabs(m.c) > 1e-10 then m.c = ccln(m.c); m.s = ccln(m.s); res[#res + 1] = m end
  end
  table.sort(res, function(a, b)
    if a.k ~= b.k then return a.k > b.k end
    if abs(a.s.re - b.s.re) > 1e-12 then return a.s.re > b.s.re end
    return a.s.im > b.s.im
  end)
  return res
end
function eqS(A, B)
  if #A ~= #B then return false end
  for i = 1, #A do
    local a, b = A[i], B[i]
    if a.k ~= b.k or not near(a.s, b.s, 1e-9) or not near(a.c, b.c, 1e-8) then return false end
  end
  return true
end

function Sig:copy()
  local r = newSig(self.dt)
  for i, b in ipairs(self.bp) do r.bp[i] = b end
  r.seg = {}
  for i, S in ipairs(self.seg) do
    local T = {}
    for j, m in ipairs(S) do T[j] = M(m.c, m.k, m.s) end
    r.seg[i] = T
  end
  for i, im in ipairs(self.imps) do r.imps[i] = {c = im.c, t0 = im.t0, m = im.m} end
  return r
end

function Sig:normalize()
  for i, S in ipairs(self.seg) do self.seg[i] = mergeS(S, self.dt) end
  -- merge equal adjacent segments
  local bp, seg = {}, {self.seg[1]}
  for i, b in ipairs(self.bp) do
    local S = self.seg[i + 1]
    if eqS(seg[#seg], S) then
      -- drop breakpoint
    else
      bp[#bp + 1] = b; seg[#seg + 1] = S
    end
  end
  self.bp, self.seg = bp, seg
  -- impulses
  local imps = {}
  for _, im in ipairs(self.imps) do
    local f = false
    for _, o in ipairs(imps) do
      if abs(o.t0 - im.t0) < 1e-9 and o.m == im.m then o.c = o.c + im.c; f = true; break end
    end
    if not f then imps[#imps + 1] = {c = im.c, t0 = im.t0, m = im.m} end
  end
  self.imps = {}
  for _, im in ipairs(imps) do
    if cabs(im.c) > 1e-10 then im.c = ccln(im.c); im.t0 = cln(im.t0); self.imps[#self.imps + 1] = im end
  end
  table.sort(self.imps, function(a, b) return a.t0 < b.t0 end)
  return self
end

-- index of segment containing t (CT: right-continuous at breakpoints)
function Sig:segIndex(t)
  local i = 1
  for j, b in ipairs(self.bp) do if t >= b - 1e-12 then i = j + 1 else break end end
  return i
end
function Sig:eval(t)
  return segEval(self.seg[self:segIndex(t)], t, self.dt)
end
-- left/right limits at t (CT)
function Sig:limits(t)
  local i = self:segIndex(t)
  local right = segEval(self.seg[i], t, self.dt)
  local li = i
  for j, b in ipairs(self.bp) do if abs(b - t) < 1e-12 then li = j end end
  local left = segEval(self.seg[li], t, self.dt)
  return left, right
end

-- split to include breakpoints in list
function Sig:refine(pts)
  local all = {}
  for _, b in ipairs(self.bp) do all[#all + 1] = b end
  for _, p in ipairs(pts) do
    if p > -INF and p < INF then
      local dup = false
      for _, b in ipairs(all) do if abs(b - p) < 1e-12 then dup = true; break end end
      if not dup then all[#all + 1] = p end
    end
  end
  table.sort(all)
  local seg = {}
  for i = 1, #all + 1 do
    local lo = all[i - 1] or -INF
    local hi = all[i] or INF
    local mid
    if lo == -INF and hi == INF then mid = 0 elseif lo == -INF then mid = hi - 1 elseif hi == INF then mid = lo + 1 else mid = (lo + hi) / 2 end
    if self.dt then
      if lo == -INF then mid = hi - 1 else mid = lo end
    end
    local S = self.seg[self:segIndex(mid)]
    local T = {}
    for j, m in ipairs(S) do T[j] = M(m.c, m.k, m.s) end
    seg[i] = T
  end
  local r = newSig(self.dt)
  r.bp, r.seg = all, seg
  for i, im in ipairs(self.imps) do r.imps[i] = {c = im.c, t0 = im.t0, m = im.m} end
  return r
end

-- build signal from pieces {lo, hi, S} (summed where overlapping)
function fromPieces(dt, pieces, imps)
  local pts = {}
  for _, p in ipairs(pieces) do pts[#pts + 1] = p.lo; pts[#pts + 1] = p.hi end
  local r = newSig(dt):refine(pts)
  for i = 1, #r.seg do
    local lo = r.bp[i - 1] or -INF
    local hi = r.bp[i] or INF
    local S = {}
    for _, p in ipairs(pieces) do
      if p.lo <= lo + 1e-12 and p.hi >= hi - 1e-12 and p.hi > p.lo then
        for _, m in ipairs(p.S) do S[#S + 1] = M(m.c, m.k, m.s) end
      end
    end
    r.seg[i] = S
  end
  r.imps = {}
  for i, im in ipairs(imps or {}) do r.imps[i] = {c = im.c, t0 = im.t0, m = im.m} end
  return r:normalize()
end
function Sig:pieces()
  local P = {}
  for i, S in ipairs(self.seg) do
    if #S > 0 then P[#P + 1] = {lo = self.bp[i - 1] or -INF, hi = self.bp[i] or INF, S = S} end
  end
  return P
end

-- maximal windows for each monomial: list {lo, hi, S={m}}
function Sig:monoWindows()
  local open, res = {}, {}
  for i, S in ipairs(self.seg) do
    local lo, hi = self.bp[i - 1] or -INF, self.bp[i] or INF
    local nopen = {}
    for _, m in ipairs(S) do
      local found = nil
      for j, w in ipairs(open) do
        if w.m.k == m.k and near(w.m.s, m.s, 1e-9) and near(w.m.c, m.c, 1e-9) and abs(w.hi - lo) < 1e-12 then found = j; break end
      end
      if found then local w = table.remove(open, found); w.hi = hi; nopen[#nopen + 1] = w
      else nopen[#nopen + 1] = {lo = lo, hi = hi, m = M(m.c, m.k, m.s)} end
    end
    for _, w in ipairs(open) do res[#res + 1] = w end
    open = nopen
  end
  for _, w in ipairs(open) do res[#res + 1] = w end
  table.sort(res, function(a, b) if a.lo ~= b.lo then return a.lo < b.lo end; return a.hi < b.hi end)
  local out = {}
  for _, w in ipairs(res) do out[#out + 1] = {lo = w.lo, hi = w.hi, S = {w.m}} end
  return out
end
function sconst(c, dt) local r = newSig(dt); r.seg = {{M(c, 0, dt and 1 or 0)}}; return r:normalize() end
function svar(dt) local r = newSig(dt); r.seg = {{M(1, 1, dt and 1 or 0)}}; return r end
function swindow(lo, hi, dt, S) return fromPieces(dt, {{lo = lo, hi = hi, S = S or {M(1, 0, dt and 1 or 0)}}}) end
function sexp(s, dt) local r = newSig(dt); r.seg = {{M(1, 0, s)}}; return r end

function Sig:isZero()
  for _, S in ipairs(self.seg) do if #S > 0 then return false end end
  return #self.imps == 0
end
-- constant value if signal is constant everywhere (no impulses)
function Sig:constValue()
  if #self.bp > 0 or #self.imps > 0 then return nil end
  local S = self.seg[1]
  if #S == 0 then return ZERO end
  if #S == 1 and S[1].k == 0 and (self.dt and near(S[1].s, ONE) or (not self.dt and isz(S[1].s))) then return S[1].c end
  return nil
end

function Sig:scale(c)
  local r = self:copy()
  c = tocx(c)
  for _, S in ipairs(r.seg) do for _, m in ipairs(S) do m.c = m.c * c end end
  for _, im in ipairs(r.imps) do im.c = im.c * c end
  return r:normalize()
end

function sadd(A, B)
  local pts = {}
  for _, b in ipairs(B.bp) do pts[#pts + 1] = b end
  local Ar = A:refine(pts)
  local Br = B:refine(Ar.bp)
  for i, S in ipairs(Ar.seg) do
    for _, m in ipairs(Br.seg[i]) do S[#S + 1] = M(m.c, m.k, m.s) end
  end
  for _, im in ipairs(B.imps) do Ar.imps[#Ar.imps + 1] = {c = im.c, t0 = im.t0, m = im.m} end
  return Ar:normalize()
end
function ssub(A, B) return sadd(A, B:scale(-1)) end

function mulS(A, B, dt)
  local out = {}
  for _, a in ipairs(A) do
    for _, b in ipairs(B) do
      out[#out + 1] = M(a.c * b.c, a.k + b.k, dt and (a.s * b.s) or (a.s + b.s))
    end
  end
  return out
end
-- derivative of monomial list (CT)
function derS(S)
  local out = {}
  for _, m in ipairs(S) do
    if m.k > 0 then out[#out + 1] = M(m.c * m.k, m.k - 1, m.s) end
    if not isz(m.s, 0) then out[#out + 1] = M(m.c * m.s, m.k, m.s) end
  end
  return mergeS(out)
end
-- value of n-th derivative of signal at t (smooth part, right side)
function Sig:derivAt(t, n)
  local S = self.seg[self:segIndex(t)]
  for _ = 1, n do S = derS(S) end
  return segEval(S, t, false)
end

function smul(A, B)
  local pts = {}
  for _, b in ipairs(B.bp) do pts[#pts + 1] = b end
  local Ar = A:refine(pts)
  local Br = B:refine(Ar.bp)
  local r = newSig(A.dt)
  r.bp = Ar.bp
  for i, S in ipairs(Ar.seg) do r.seg[i] = mulS(S, Br.seg[i], A.dt) end
  local function impTimes(imps, other)
    for _, im in ipairs(imps) do
      -- f(t) delta^(m)(t-t0) = sum_i (-1)^i C(m,i) f^(i)(t0) delta^(m-i)(t-t0)
      for i = 0, im.m do
        local fv = (i == 0) and other:eval(im.t0) or other:derivAt(im.t0, i)
        if i == 0 then
          -- use average of limits if discontinuity at t0 inside smooth part
          local l, rr = other:limits(im.t0)
          if not near(l, rr) then fv = rr end
        end
        local c = im.c * fv * (((i % 2) == 1) and -1 or 1) * binom(im.m, i)
        r.imps[#r.imps + 1] = {c = c, t0 = im.t0, m = im.m - i}
      end
      if #other.imps > 0 then
        for _, o in ipairs(other.imps) do
          if abs(o.t0 - im.t0) < 1e-12 then error("product of two impulses is undefined", 0) end
        end
      end
    end
  end
  impTimes(A.imps, B)
  impTimes(B.imps, A)
  return r:normalize()
end

-- x(alpha*t + beta)
function Sig:compose(alpha, beta)
  if abs(alpha) < 1e-14 then return sconst(self:eval(beta), self.dt) end
  local dt = self.dt
  local pieces = {}
  local function mapS(S)
    local out = {}
    for _, m in ipairs(S) do
      local base
      if dt then base = m.c * cpow(m.s, beta) else base = m.c * cexp(m.s * beta) end
      local ns = dt and cpow(m.s, alpha) or (m.s * alpha)
      for i = 0, m.k do
        local cf = binom(m.k, i) * (alpha ^ i) * (beta ^ (m.k - i))
        if cf ~= 0 then out[#out + 1] = M(base * cf, i, ns) end
      end
    end
    return out
  end
  if dt then
    if not isint(alpha) then error("x[a·n+b] needs integer a", 0) end
    if not isint(beta) then error("x[a·n+b] needs integer b", 0) end
    alpha, beta = round(alpha), round(beta)
  end
  for i, S in ipairs(self.seg) do
    if #S > 0 then
      local p = self.bp[i - 1] or -INF
      local q = self.bp[i] or INF
      local lo, hi
      if not dt then
        local a1, a2 = (p - beta) / alpha, (q - beta) / alpha
        if p == -INF then a1 = (alpha > 0) and -INF or INF end
        if q == INF then a2 = (alpha > 0) and INF or -INF end
        lo, hi = math.min(a1, a2), math.max(a1, a2)
      else
        -- integers n with p <= alpha n + beta < q
        local function cdiv(x) if x == INF or x == -INF then return x end return math.ceil(x - 1e-9) end
        local function fdiv(x) if x == INF or x == -INF then return x end return floor(x + 1e-9) end
        if alpha > 0 then
          lo = (p == -INF) and -INF or cdiv((p - beta) / alpha)
          hi = (q == INF) and INF or cdiv((q - beta) / alpha)
        else
          lo = (q == INF) and -INF or (fdiv((q - beta) / alpha) + 1)
          hi = (p == -INF) and INF or (fdiv((p - beta) / alpha) + 1)
        end
      end
      if hi > lo then pieces[#pieces + 1] = {lo = lo, hi = hi, S = mapS(S)} end
    end
  end
  local imps = {}
  for _, im in ipairs(self.imps) do
    imps[#imps + 1] = {c = im.c / (abs(alpha) * alpha ^ im.m), t0 = (im.t0 - beta) / alpha, m = im.m}
  end
  return fromPieces(dt, pieces, imps)
end

-- derivative (CT): smooth part + impulses at jumps
function Sig:deriv()
  if self.dt then error("use first difference for DT", 0) end
  local r = newSig(false)
  r.bp = {}
  for i, b in ipairs(self.bp) do r.bp[i] = b end
  for i, S in ipairs(self.seg) do r.seg[i] = derS(S) end
  for _, b in ipairs(self.bp) do
    local l, rr = self:limits(b)
    local jmp = rr - l
    if cabs(jmp) > 1e-10 then r.imps[#r.imps + 1] = {c = jmp, t0 = b, m = 0} end
  end
  for _, im in ipairs(self.imps) do r.imps[#r.imps + 1] = {c = im.c, t0 = im.t0, m = im.m + 1} end
  return r:normalize()
end

-- antiderivative of monomial list (CT): returns S'
function antiS(S)
  local out = {}
  for _, m in ipairs(S) do
    if isz(m.s, 1e-12) then
      out[#out + 1] = M(m.c / (m.k + 1), m.k + 1, ZERO)
    else
      local fk = fact(m.k)
      for i = 0, m.k do
        local cf = ((i % 2 == 1) and -1 or 1) * fk / fact(m.k - i)
        out[#out + 1] = M(m.c * cf / cpow(m.s, i + 1), m.k - i, m.s)
      end
    end
  end
  return mergeS(out)
end
-- polynomial Q for antidifference: r Q(n+1) - Q(n) = n^k  (r~=1, deg k) or Q(n+1)-Q(n)=n^k (r==1, deg k+1, Q(0)=0)
function solveLin(A, b)
  local n = #b
  for i = 1, n do
    local piv = i
    for r = i + 1, n do if cabs(A[r][i]) > cabs(A[piv][i]) then piv = r end end
    A[i], A[piv] = A[piv], A[i]; b[i], b[piv] = b[piv], b[i]
    if cabs(A[i][i]) < 1e-14 then error("singular system", 0) end
    for r = i + 1, n do
      local f = A[r][i] / A[i][i]
      for c = i, n do A[r][c] = A[r][c] - f * A[i][c] end
      b[r] = b[r] - f * b[i]
    end
  end
  local x = {}
  for i = n, 1, -1 do
    local s = b[i]
    for c = i + 1, n do s = s - A[i][c] * x[c] end
    x[i] = s / A[i][i]
  end
  return x
end
-- returns coefficient array q[0..d] of Q
function antidiffQ(k, r)
  r = tocx(r)
  local one = near(r, ONE, 1e-12)
  local d = one and (k + 1) or k
  -- unknowns q_0..q_d ; equations: coefficient of n^p for p=0..d (plus Q(0)=0 if one)
  local nunk = d + 1
  local A, b = {}, {}
  local rows = {}
  -- expression: r*Q(n+1) - Q(n) = sum_j q_j (r (n+1)^j - n^j)
  for p = 0, d do
    local row = {}
    for j = 0, d do
      local v = ZERO
      if j >= p then v = v + r * binom(j, p) end
      if j == p then v = v - ONE end
      row[j + 1] = tocx(v)
    end
    rows[#rows + 1] = {row, cx((p == k) and 1 or 0)}
  end
  if one then
    -- drop equation for highest power p=d (automatically consistent) and add Q(0)=0
    rows[#rows] = nil
    local row = {}
    for j = 0, d do row[j + 1] = cx(j == 0 and 1 or 0) end
    rows[#rows + 1] = {row, ZERO}
  end
  for i, rw in ipairs(rows) do A[i] = rw[1]; b[i] = rw[2] end
  local x = solveLin(A, b)
  local q = {}
  for j = 0, d do q[j] = x[j + 1] end
  return q, d
end
-- DT running-sum antidifference of monomial list: G with G[n]-G[n-1] = S[n]
function antiDiffS(S)
  local out = {}
  for _, m in ipairs(S) do
    local q, d = antidiffQ(m.k, m.s)
    -- H[n] = r^n Q(n);  G[n] = H[n+1] = r * r^n Q(n+1)
    for j = 0, d do
      for i = 0, j do
        local cf = q[j] * binom(j, i)   -- (n+1)^j = sum C(j,i) n^i
        if not isz(cf, 1e-15) then out[#out + 1] = M(m.c * m.s * cf, i, m.s) end
      end
    end
  end
  return mergeS(out)
end
-- limit of monomial list at +inf (dir=1) or -inf (dir=-1); returns value or nil if diverges
function limitS(S, dir, dt)
  local v = ZERO
  for _, m in ipairs(S) do
    local isconstM = (m.k == 0) and (dt and near(m.s, ONE, 1e-12) or (not dt and isz(m.s, 1e-12)))
    if isconstM then v = v + m.c
    else
      local decays
      if dt then
        local a = cabs(m.s)
        decays = (dir > 0 and a < 1 - 1e-12) or (dir < 0 and a > 1 + 1e-12)
      else
        decays = (dir > 0 and m.s.re < -1e-12) or (dir < 0 and m.s.re > 1e-12)
      end
      if not decays then return nil end
    end
  end
  return v
end

-- running integral (CT) / running sum (DT) from -inf.  If lower limit diverges and allowFree, uses arbitrary constant
function Sig:running(allowFree)
  local dt = self.dt
  local src = self
  if not dt then
    local pts = {}
    for _, im in ipairs(self.imps) do pts[#pts + 1] = im.t0 end
    src = self:refine(pts)
  end
  local r = newSig(dt)
  r.bp = {}
  for i, b in ipairs(src.bp) do r.bp[i] = b end
  local diverged = false
  local prevS, prevK = nil, ZERO
  for i, S in ipairs(src.seg) do
    local G = dt and antiDiffS(S) or antiS(S)
    local K
    if i == 1 then
      local lim = limitS(G, -1, dt)
      if lim == nil then
        if #G == 0 then lim = ZERO else diverged = true; lim = ZERO end
      end
      K = -lim
    else
      local b = src.bp[i - 1]
      if dt then
        local prevVal = segEval(prevS, b - 1, true) + prevK
        K = prevVal - segEval(G, b - 1, true)
      else
        local prevVal = segEval(prevS, b, false) + prevK
        for _, im in ipairs(src.imps) do
          if abs(im.t0 - b) < 1e-12 and im.m == 0 then prevVal = prevVal + im.c end
        end
        K = prevVal - segEval(G, b, false)
      end
    end
    local T = {}
    for j, m in ipairs(G) do T[j] = m end
    if not isz(K, 1e-13) then T[#T + 1] = M(K, 0, dt and 1 or 0) end
    r.seg[i] = T
    prevS, prevK = G, K
  end
  for _, im in ipairs(src.imps) do
    if im.m > 0 then r.imps[#r.imps + 1] = {c = im.c, t0 = im.t0, m = im.m - 1} end
  end
  if diverged and not allowFree then error("running integral diverges (signal does not vanish as t→-∞)", 0) end
  return r:normalize(), diverged
end

-- definite integral / sum over [a,b] (CT: a<b reals or ±inf; DT inclusive integers)
function Sig:integrate(a, b)
  local sg = 1
  if a > b then a, b = b, a; sg = -1 end
  local F, div = self:running(true)
  local function at(x, side)
    if x == INF then
      local v = limitS(F.seg[#F.seg], 1, self.dt)
      if v == nil then error("integral diverges at +∞", 0) end
      return v
    elseif x == -INF then
      if div then error("integral diverges at -∞", 0) end
      return ZERO
    end
    if self.dt then return F:eval(x) end
    -- CT: include impulses at endpoints fully
    local l, r = F:limits(x)
    if side == "lo" then return l end
    return r
  end
  local va
  if self.dt then va = (a == -INF) and at(a) or F:eval(a - 1) else va = at(a, "lo") end
  if self.dt and a == -INF then va = at(-INF) end
  local vb = at(b, "hi")
  return (vb - va) * sg
end

-- ---------- convolution ----------
-- antiderivative of tau^p e^{D tau} as list {coef, pow, hasExp}
function antiTerms(p, D, dt)
  local res = {}
  if dt then
    local q, d = antidiffQ(p, D)
    for j = 0, d do res[#res + 1] = {c = q[j], p = j} end
    return res, true
  end
  if isz(D, 1e-12) then res[1] = {c = cx(1 / (p + 1)), p = p + 1}; return res, false end
  local fp = fact(p)
  for i = 0, p do
    res[#res + 1] = {c = cx(((i % 2 == 1) and -1 or 1) * fp / fact(p - i)) / cpow(D, i + 1), p = p - i}
  end
  return res, true
end
-- evaluate antiderivative at a limit: lim = {kind="c", v} or {kind="t", c0} (tau = t - c0)
-- returns monomial list in t (factor e^{D tau} included)
function antiAt(p, D, lim, dt)
  local terms, hasE = antiTerms(p, D, dt)
  if dt then hasE = not near(D, ONE, 1e-12) end
  local out = {}
  if lim.kind == "c" then
    local x = lim.v
    if x == INF or x == -INF then
      -- limit must be zero (decaying)
      local dir = (x == INF) and 1 or -1
      local ok
      if dt then
        local a = cabs(D)
        ok = hasE and ((dir > 0 and a < 1 - 1e-12) or (dir < 0 and a > 1 + 1e-12))
      else
        ok = hasE and ((dir > 0 and D.re < -1e-12) or (dir < 0 and D.re > 1e-12))
      end
      if not ok then error("convolution integral diverges", 0) end
      return {}
    end
    local v = ZERO
    for _, tm in ipairs(terms) do v = v + tm.c * cx(x ^ tm.p) end
    if hasE then v = v * (dt and cpow(D, x) or cexp(D * x)) end
    out[1] = M(v, 0, dt and 1 or 0)
    return out
  end
  -- tau = t - c0
  local c0 = lim.c0
  local ef = hasE and (dt and cpow(D, -c0) or cexp(D * (-c0))) or ONE
  local s = hasE and D or (dt and ONE or ZERO)
  for _, tm in ipairs(terms) do
    for i = 0, tm.p do
      local cf = binom(tm.p, i) * ((-c0) ^ (tm.p - i))
      if cf ~= 0 then out[#out + 1] = M(tm.c * ef * cf, i, s) end
    end
  end
  return out
end

-- convolve windowed monomials. A: {lo,hi,m}, B: {lo,hi,m}
function convPieces(pa, ma, pb, mb, dt, out)
  local a1, b1, a2, b2 = pa.lo, pa.hi, pb.lo, pb.hi
  local D = dt and (ma.s / mb.s) or (ma.s - mb.s)
  local pts
  if dt then pts = {a1 + a2, a1 + b2, b1 + a2, b1 + b2 - 1}
  else pts = {a1 + a2, a1 + b2, b1 + a2, b1 + b2} end
  local fin = {}
  for _, p in ipairs(pts) do if p == p and p > -INF and p < INF then fin[#fin + 1] = p end end
  table.sort(fin)
  local uniq = {}
  for _, p in ipairs(fin) do if #uniq == 0 or abs(uniq[#uniq] - p) > 1e-12 then uniq[#uniq + 1] = p end end
  local lo0 = a1 + a2; if lo0 ~= lo0 then lo0 = -INF end
  local hi0 = dt and (b1 + b2 - 1) or (b1 + b2); if hi0 ~= hi0 then hi0 = INF end
  local edges = {lo0}
  for _, p in ipairs(uniq) do if p > lo0 + 1e-12 and p < hi0 - 1e-12 then edges[#edges + 1] = p end end
  edges[#edges + 1] = hi0
  for ii = 1, #edges - 1 do
    local lo, hi = edges[ii], edges[ii + 1]
    if hi > lo + 1e-12 then
      local tm
      if lo == -INF and hi == INF then tm = 0 elseif lo == -INF then tm = hi - 7.3 elseif hi == INF then tm = lo + 7.3 else tm = (lo + hi) / 2 end
      if dt then if lo == -INF and hi == INF then tm = 0 else tm = (lo == -INF) and (hi - 1) or lo end end
      local L, U, Lk, Uk
      -- lower limit
      local cand = dt and (tm - b2 + 1) or (tm - b2)
      if a1 >= cand then L = {kind = "c", v = a1} else L = {kind = "t", c0 = dt and (b2 - 1) or b2} end
      if b2 == INF then L = {kind = "c", v = a1} end
      local candU = tm - a2
      local ub = dt and (b1 - 1) or b1
      if ub <= candU then U = {kind = "c", v = ub} else U = {kind = "t", c0 = a2} end
      if a2 == -INF then U = {kind = "c", v = ub} end
      -- for DT, sum_{k=L}^{U} = H[U+1]-H[L]
      if dt then
        if U.kind == "c" then U = {kind = "c", v = (U.v == INF) and INF or U.v + 1} else U = {kind = "t", c0 = U.c0 - 1} end
      end
      local S = {}
      for i = 0, mb.k do
        local cf = binom(mb.k, i) * (((i % 2) == 1) and -1 or 1)
        local up = antiAt(ma.k + i, D, U, dt)
        local dn = antiAt(ma.k + i, D, L, dt)
        for _, m in ipairs(dn) do m.c = -m.c; up[#up + 1] = m end
        -- multiply by c1 c2 t^{kb-i} e^{sb t}
        for _, m in ipairs(up) do
          S[#S + 1] = M(m.c * ma.c * mb.c * cf, m.k + mb.k - i, dt and (m.s * mb.s) or (m.s + mb.s))
        end
      end
      out[#out + 1] = {lo = lo, hi = dt and ((hi == INF) and INF or hi) or hi, S = mergeS(S, dt)}
    end
  end
end

function sconv(A, B)
  if A.dt ~= B.dt then error("cannot convolve CT with DT", 0) end
  local dt = A.dt
  local out, imps = {}, {}
  for _, pa in ipairs(A:pieces()) do
    for _, ma in ipairs(pa.S) do
      for _, pb in ipairs(B:pieces()) do
        for _, mb in ipairs(pb.S) do
          local piece = {}
          if dt then
            convPieces({lo = pa.lo, hi = pa.hi}, ma, {lo = pb.lo, hi = pb.hi}, mb, dt, piece)
            -- fix DT hi: edges computed with last valid n = b1+b2-2 so half-open bound b1+b2-1
          else
            convPieces({lo = pa.lo, hi = pa.hi}, ma, {lo = pb.lo, hi = pb.hi}, mb, dt, piece)
          end
          for _, p in ipairs(piece) do out[#out + 1] = p end
        end
      end
    end
  end
  local res = fromPieces(dt, out, {})
  -- impulse parts
  local function impConv(imps1, other)
    local acc = nil
    for _, im in ipairs(imps1) do
      local d = other
      for _ = 1, im.m do d = d:deriv() end
      local sh = d:compose(1, -im.t0):scale(im.c)
      acc = acc and sadd(acc, sh) or sh
    end
    return acc
  end
  if not dt then
    local smoothA = A:copy(); smoothA.imps = {}
    local smoothB = B:copy(); smoothB.imps = {}
    local x1 = impConv(A.imps, smoothB)
    local x2 = impConv(B.imps, smoothA)
    if x1 then res = sadd(res, x1) end
    if x2 then res = sadd(res, x2) end
    for _, ia in ipairs(A.imps) do
      for _, ib in ipairs(B.imps) do
        res.imps[#res.imps + 1] = {c = ia.c * ib.c, t0 = ia.t0 + ib.t0, m = ia.m + ib.m}
      end
    end
    res:normalize()
  end
  return res
end
-- ===================== AST -> Sig conversion =====================
function baseSig(name, dt)
  local one = dt and 1 or 0
  if name == "u" then return swindow(0, INF, dt) end
  if name == "rect" then
    if dt then return swindow(0, 1, dt) end   -- rarely used in DT
    return swindow(-0.5, 0.5, dt)
  end
  if name == "tri" then
    return fromPieces(dt, {{lo = -1, hi = 0, S = {M(1, 0, one), M(1, 1, one)}}, {lo = 0, hi = 1, S = {M(1, 0, one), M(-1, 1, one)}}})
  end
  if name == "ramp" then return fromPieces(dt, {{lo = 0, hi = INF, S = {M(1, 1, one)}}}) end
  if name == "abs" then return fromPieces(dt, {{lo = -INF, hi = 0, S = {M(-1, 1, one)}}, {lo = 0, hi = INF, S = {M(1, 1, one)}}}) end
  if name == "sgn" or name == "sign" then
    return fromPieces(dt, {{lo = -INF, hi = 0, S = {M(-1, 0, one)}}, {lo = dt and 1 or 0, hi = INF, S = {M(1, 0, one)}}})
  end
  if name == "delta" then
    if dt then return swindow(0, 1, dt) end
    local r = newSig(false); r.imps = {{c = ONE, t0 = 0, m = 0}}; return r
  end
  return nil
end

Conv = {}
Conv.__index = Conv
-- ctx: {var=, dt=, env=, defs = {name={param=,body=,dt=}}, cache={}}
function mkConv(var, dt, env, defs)
  return setmetatable({var = var, dt = dt, env = env, defs = defs or {}, cache = {}}, Conv)
end

function Conv:lin(node, var)
  local a, b = linearOf(node, var or self.var, self.env)
  return a, b
end
function cre(z, what) z = tocx(z); if abs(z.im) > 1e-9 then error(what .. " must be real", 0) end; return z.re end

function Conv:userSig(name)
  if self.cache[name] then return self.cache[name] end
  local d = self.defs[name]
  if not d then error("'" .. name .. "' is not defined", 0) end
  if self.busy and self.busy[name] then error("recursive definition of " .. name, 0) end
  self.busy = self.busy or {}
  self.busy[name] = true
  local sub = mkConv(d.param, d.dt, self.env, self.defs)
  sub.cache = self.cache; sub.busy = self.busy
  local ok, s = pcall(sub.sig, sub, d.body)
  self.busy[name] = nil
  if not ok then error(s, 0) end
  self.cache[name] = s
  return s
end

function Conv:sig(node)
  local k = node.k
  local dt, var = self.dt, self.var
  if k == "num" then
    if node.v.re == INF or node.v.re == -INF then error("∞ not allowed here", 0) end
    return sconst(node.v, dt)
  elseif k == "var" then
    if node.name == var then return svar(dt) end
    if self.defs[node.name] then error("use " .. node.name .. "(" .. var .. ")", 0) end
    return sconst(evalN(node, self.env), dt)
  elseif k == "add" then return sadd(self:sig(node.a), self:sig(node.b))
  elseif k == "sub" then return ssub(self:sig(node.a), self:sig(node.b))
  elseif k == "neg" then return self:sig(node.a):scale(-1)
  elseif k == "mul" then
    local A = self:sig(node.a)
    local ca = A:constValue()
    if ca and isz(ca) then return newSig(dt) end
    return smul(A, self:sig(node.b))
  elseif k == "div" then
    local B = self:sig(node.b)
    local cb = B:constValue()
    if not cb or isz(cb) then error("division by a non-constant signal is not closed-form", 0) end
    return self:sig(node.a):scale(ONE / cb)
  elseif k == "pow" then
    local ce = isconst(node.b, var, self.env)
    local cb = isconst(node.a, var, self.env)
    if ce then
      if cb then return sconst(cpow(cb, ce), dt) end
      local n = cre(ce, "exponent")
      if n >= 0 and isint(n) then
        local B = self:sig(node.a)
        local r = sconst(1, dt)
        for _ = 1, round(n) do r = smul(r, B) end
        return r
      end
      error("non-integer power of a signal is not closed-form", 0)
    end
    if cb then
      local a, b = self:lin(node.b)
      if not a then error("exponent must be linear in " .. var, 0) end
      local coef = cpow(cb, b)
      if dt then
        return sexp(cpow(cb, a), true):scale(coef)
      else
        return sexp(clog(cb) * a, false):scale(coef)
      end
    end
    error("unsupported power", 0)
  elseif k == "cmp" then
    local D = ssub(self:sig(node.a), self:sig(node.b))
    return self:indicator(D, node.op)
  elseif k == "call" then
    local nm, A = node.name, node.args
    if self.defs[nm] then
      local a, b = self:lin(A[1])
      if not a then error("argument of " .. nm .. "() must be linear in " .. var, 0) end
      return self:userSig(nm):compose(cre(a, "scale"), cre(b, "shift"))
    end
    if nm == "exp" or nm == "sin" or nm == "cos" or nm == "sinh" or nm == "cosh" then
      local c0 = isconst(A[1], var, self.env)
      if c0 then return sconst(evalN(node, self.env), dt) end
      local a, b = self:lin(A[1])
      if not a and nm == "exp" then
        -- piecewise-linear exponent
        local X = self:sig(A[1])
        local pieces = {}
        for _, pc in ipairs(X:pieces()) do
          local c0, c1 = ZERO, ZERO
          for _, m in ipairs(pc.S) do
            local polyOk = dt and near(m.s, ONE, 1e-12) or (not dt and isz(m.s, 1e-12))
            if not polyOk or m.k > 1 then error("exp() argument must be (piecewise) linear", 0) end
            if m.k == 0 then c0 = c0 + m.c else c1 = c1 + m.c end
          end
          pieces[#pieces + 1] = {lo = pc.lo, hi = pc.hi, S = {M(cexp(c0), 0, dt and cexp(c1) or c1)}}
        end
        -- zero-exponent gaps give e^0 = 1
        local covered = fromPieces(dt, pieces)
        local one = dt and 1 or 0
        local gaps = {}
        for i, S in ipairs(X.seg) do if #S == 0 then gaps[#gaps + 1] = {lo = X.bp[i - 1] or -INF, hi = X.bp[i] or INF, S = {M(1, 0, one)}} end end
        if #gaps > 0 then covered = sadd(covered, fromPieces(dt, gaps)) end
        return covered
      end
      if not a then error(nm .. "() argument must be linear in " .. var, 0) end
      local function E(mult) -- e^{mult*(a v + b)}
        local s = mult * a
        return sexp(dt and cexp(s) or s, dt):scale(cexp(mult * b))
      end
      if nm == "exp" then return E(ONE) end
      if nm == "cos" then return sadd(E(J), E(-J)):scale(0.5) end
      if nm == "sin" then return ssub(E(J), E(-J)):scale(ONE / cx(0, 2)) end
      if nm == "cosh" then return sadd(E(ONE), E(-ONE)):scale(0.5) end
      if nm == "sinh" then return ssub(E(ONE), E(-ONE)):scale(0.5) end
    end
    local base = baseSig(nm, dt)
    if base then
      local c0 = isconst(A[1], var, self.env)
      if c0 then
        if nm == "delta" then
          if dt then return sconst(isz(c0) and 1 or 0, dt) end
          error("δ(constant) is not a signal", 0)
        end
        return sconst(evalN(node, self.env), dt)
      end
      local a, b = self:lin(A[1])
      if a then return base:compose(cre(a, "scale"), cre(b, "shift")) end
      if nm == "abs" or nm == "sgn" or nm == "sign" then
        local X = self:sig(A[1])
        local P = self:indicator(X, ">")
        local Nn = self:indicator(X, "<")
        local sgn = ssub(P, Nn)
        if nm == "abs" then return smul(X, sgn) end
        return sgn
      end
      if nm == "u" then return self:indicator(self:sig(A[1]), ">=") end
      error(nm .. "() argument must be linear in " .. var, 0)
    end
    if nm == "int" or nm == "integral" or nm == "sum" then return self:integral(node) end
    if nm == "der" or nm == "der2" then
      if dt then error("derivative not defined for DT; use x[n]-x[n-1]", 0) end
      local X = self:sig(A[1]):deriv()
      if nm == "der2" then X = X:deriv() end
      return X
    end
    if nm == "conv" then return sconv(self:sig(A[1]), self:sig(A[2])) end
    if nm == "conj" or nm == "re" or nm == "real" or nm == "im" or nm == "imag" then
      local X = self:sig(A[1])
      local Cj = X:copy()
      for _, S in ipairs(Cj.seg) do for _, m in ipairs(S) do m.c = cconj(m.c); m.s = cconj(m.s) end end
      for _, im in ipairs(Cj.imps) do im.c = cconj(im.c) end
      Cj:normalize()
      if nm == "conj" then return Cj end
      if nm == "re" or nm == "real" then return sadd(X, Cj):scale(0.5) end
      return ssub(X, Cj):scale(ONE / cx(0, 2))
    end
    local c0 = isconst(node, var, self.env)
    if c0 then return sconst(c0, dt) end
    error(nm .. "() is not supported in exact mode", 0)
  end
  error("unsupported expression", 0)
end

-- indicator {D op 0} for piecewise-linear D
function Conv:indicator(D, op)
  local dt = self.dt
  local pieces = {}
  local strict = (op == ">" or op == "<")
  local wantPos = (op == ">" or op == ">=")
  local one = dt and 1 or 0
  for i, S in ipairs(D.seg) do
    local lo, hi = D.bp[i - 1] or -INF, D.bp[i] or INF
    -- S must be degree<=1 polynomial
    local c0, c1 = ZERO, ZERO
    for _, m in ipairs(S) do
      local polyOk = dt and near(m.s, ONE, 1e-12) or (not dt and isz(m.s, 1e-12))
      if not polyOk or m.k > 1 then error("condition must be (piecewise) linear", 0) end
      if m.k == 0 then c0 = c0 + m.c else c1 = c1 + m.c end
    end
    local a, b = cre(c1, "condition"), cre(c0, "condition")
    local function test(v)
      if wantPos then return strict and v > 1e-12 or v >= -1e-12 end
      return strict and v < -1e-12 or v <= 1e-12
    end
    if abs(a) < 1e-14 then
      if test(b) then pieces[#pieces + 1] = {lo = lo, hi = hi, S = {M(1, 0, one)}} end
    else
      local r = -b / a
      local L, H
      local posRight = a > 0
      local wantRight = (posRight == wantPos)
      if dt then
        -- integers n in [lo,hi) satisfying
        local rr
        if wantRight then
          -- n >= r (non-strict) or n > r (strict) if sign matches
          rr = (strict and (isint(r) and round(r) + 1 or math.ceil(r))) or math.ceil(r - 1e-9)
          L, H = math.max(lo, rr), hi
        else
          rr = (strict and (isint(r) and round(r) - 1 or floor(r))) or floor(r + 1e-9)
          L, H = lo, math.min(hi, rr + 1)
        end
      else
        if wantRight then L, H = math.max(lo, r), hi else L, H = lo, math.min(hi, r) end
      end
      if H > L then pieces[#pieces + 1] = {lo = L, hi = H, S = {M(1, 0, one)}} end
    end
  end
  return fromPieces(dt, pieces)
end

function Conv:integral(node)
  local A = node.args
  if #A ~= 4 then error(node.name .. "(f, var, lower, upper)", 0) end
  local iv = A[2].name
  if not iv then error("2nd argument must be the integration variable", 0) end
  local dt = self.dt
  if hasvar(A[1], self.var) and iv ~= self.var then error("integrand depends on " .. self.var .. " (not closed-form here)", 0) end
  local sub = mkConv(iv, dt, self.env, self.defs)
  sub.cache = {}
  local f = sub:sig(A[1])
  local F, div = f:running(true)
  local function atLimit(lim, isLower)
    local lc = isconst(lim, self.var, self.env)
    if lc and (lc.re == INF or lc.re == -INF) then
      local dir = lc.re > 0 and 1 or -1
      if dir < 0 then
        if div then error("integral diverges at -∞", 0) end
        return sconst(0, dt)
      end
      local v = limitS(F.seg[#F.seg], 1, dt)
      if v == nil then error("integral diverges at +∞", 0) end
      return sconst(v, dt)
    end
    local a, b = self:lin(lim)
    if not a then error("integration limits must be linear in " .. self.var, 0) end
    a, b = cre(a, "limit"), cre(b, "limit")
    if dt and isLower then b = b - 1 end
    if abs(a) < 1e-14 then return sconst(F:eval(b), dt) end
    return F:compose(a, b)
  end
  local hi = atLimit(A[4], false)
  local lo = atLimit(A[3], true)
  local R = ssub(hi, lo)
  R.imps = {}
  return R:normalize()
end

-- ===================== Sig printing =====================
function expoArg(sv, var)
  -- string for (sv)*var, like "2t", "-t", "t/2", "πt/3"
  sv = cln(sv)
  if sv == 1 then return var end
  if sv == -1 then return "-" .. var end
  local s = fmtR(sv)
  local p, q = s:match("^(%-?%d*)/(%d+)$")
  if not p then p, q = s:match("^(%-?%d*π)/(%d+)$") end
  if p then
    if p == "" then p = "1" end
    if p == "1" then return var .. "/" .. q end
    if p == "-1" then return "-" .. var .. "/" .. q end
    if p == "π" then return "π" .. var .. "/" .. q end
    if p == "-π" then return "-π" .. var .. "/" .. q end
    return p .. var .. "/" .. q
  end
  if s:find("[%+ ]") or s:find(".%-") then return "(" .. s .. ")" .. var end
  if s:find("[√e]") then return s .. "·" .. var end
  return s .. var
end
function expStr(s, var)
  -- e^{s var}, s complex
  s = tocx(s)
  if isz(s, 1e-12) then return "" end
  if abs(s.im) < 1e-12 then
    local a = expoArg(s.re, var)
    if a == var then return "e^" .. var end
    return "e^(" .. a .. ")"
  end
  if abs(s.re) < 1e-12 then
    local a = expoArg(abs(s.im), var)
    return "e^(" .. (s.im < 0 and "-" or "") .. "j" .. a .. ")"
  end
  return "e^((" .. fmtC(s) .. ")" .. var .. ")"
end
function powStr(var, k)
  if k == 0 then return "" elseif k == 1 then return var end
  return var .. "^" .. k
end
function baseStr(r, var)
  r = tocx(r)
  if near(r, ONE, 1e-12) then return "" end
  local s = fmtC(r)
  if isreal(r) and r.re > 0 and not s:find("[/%-%s√π·]") then return s .. "^" .. var end
  return "(" .. s .. ")^" .. var
end
function joinf(...)
  local t = {}
  for _, s in ipairs({...}) do if s ~= "" then t[#t + 1] = s end end
  local out = ""
  for i, s in ipairs(t) do
    if i > 1 and out:find("[%a]$") and s:find("^[%a]") then out = out .. "·" .. s else out = out .. s end
  end
  return out
end

-- format monomial list as a string (combining conjugate pairs into cos/sin)
function segStr(S, var, dt)
  local used = {}
  local terms = {}
  for i, m in ipairs(S) do
    if not used[i] then
      used[i] = true
      local paired = nil
      local hasIm = dt and (abs(carg(m.s)) > 1e-10) or (not dt and abs(m.s.im) > 1e-12)
      if hasIm then
        for j = i + 1, #S do
          local o = S[j]
          local sym = dt and near(o.s, cconj(m.s), 1e-9) and near(cabs(m.s), 1, 1e-9) or near(o.s, cconj(m.s), 1e-9)
          if not used[j] and o.k == m.k and sym and (near(o.c, cconj(m.c), 1e-8) or (isreal(m.c * cx(0,1)) or isreal(m.c)) and (isreal(o.c * cx(0,1)) or isreal(o.c))) then
            paired = j; break
          end
        end
      end
      if paired then
        used[paired] = true
        local mm = m
        local w
        if dt then w = carg(m.s) else w = m.s.im end
        if w < 0 then mm = S[paired]; w = -w end
        -- c1 e^{jw} + c2 e^{-jw} = (c1+c2) cos + j(c1-c2) sin
        local other = (mm == m) and S[paired] or m
        local Ac, Bc = mm.c + other.c, (mm.c - other.c) * J
        local A, B = Ac.re, Bc.re
        local cplx = not (isreal(Ac) and isreal(Bc))
        local env, wstr
        if dt then
          local rho = cabs(mm.s)
          env = near(rho, 1, 1e-12) and "" or baseStr(rho, var)
        else
          env = expStr(cx(mm.s.re, 0), var)
        end
        wstr = expoArg(w, var)
        local tp = powStr(var, mm.k)
        local inner = {}
        -- try amplitude-phase if both nonzero? keep A cos + B sin
        if cplx then A, B = Ac, Bc else A, B = cx(A), cx(B) end
        if cabs(A) > 1e-10 then inner[#inner + 1] = cmul(A, "cos(" .. wstr .. ")") end
        if cabs(B) > 1e-10 then inner[#inner + 1] = cmul(B, "sin(" .. wstr .. ")") end
        local pre = joinf(tp, env)
        if #inner == 1 then
          local s = inner[1]
          local neg = s:sub(1, 1) == "-"
          if neg then s = s:sub(2) end
          if pre ~= "" then
            -- coefficient first
            local cpart, fpart = s:match("^(.-)([cs][oi][sn]%(.*)$")
            s = joinf(cpart or "", pre, fpart or s)
            if cpart == nil then s = joinf(pre, s) end
          end
          terms[#terms + 1] = (neg and "-" or "") .. s
        else
          if pre ~= "" then terms[#terms + 1] = pre .. "[" .. joinsum(inner) .. "]"
          else for _, s in ipairs(inner) do terms[#terms + 1] = s end end
        end
      else
        local f
        if dt then f = joinf(powStr(var, m.k), baseStr(m.s, var))
        else f = joinf(powStr(var, m.k), expStr(m.s, var)) end
        terms[#terms + 1] = cmul(m.c, f)
      end
    end
  end
  return joinsum(terms)
end

function impStr(im, var)
  local d = "δ" .. string.rep("'", im.m)
  local a = (abs(im.t0) < 1e-12) and var or (var .. (im.t0 > 0 and " - " or " + ") .. fmtR(abs(im.t0)))
  return cmul(im.c, d .. "(" .. a .. ")")
end
function shiftStr(var, a)
  -- (var - a)
  if abs(a) < 1e-12 then return var end
  if a > 0 then return var .. " - " .. fmtR(a) end
  return var .. " + " .. fmtR(-a)
end
function uStr(var, a, dt, rev)
  local o, c = dt and "[" or "(", dt and "]" or ")"
  if rev then
    -- u(a - var)
    if abs(a) < 1e-12 then return "u" .. o .. "-" .. var .. c end
    return "u" .. o .. fmtR(a) .. " - " .. var .. c
  end
  return "u" .. o .. shiftStr(var, a) .. c
end
function windowStr(lo, hi, var, dt)
  local o, c = dt and "[" or "(", dt and "]" or ")"
  if dt and hi ~= INF and lo ~= -INF and hi - lo == 1 then return "δ" .. o .. shiftStr(var, lo) .. c end
  if lo == -INF and hi == INF then return "" end
  if hi == INF then return uStr(var, lo, dt) end
  if lo == -INF then
    if dt then return uStr(var, hi - 1, dt, true) end
    return uStr(var, hi, dt, true)
  end
  return "[" .. uStr(var, lo, dt) .. " - " .. uStr(var, hi, dt) .. "]"
end
function attach(expr, w)
  if w == "" then return expr end
  if expr == "1" then return w end
  if expr == "-1" then return "-" .. w end
  local single = not topSum(expr)
  if single and expr:find("^%-?[%d]+/[%d]+$") then return "(" .. expr .. ")" .. w end
  if single then return joinf(expr, w) end
  return "(" .. expr .. ")" .. w
end

-- interval description
function rangeStr(lo, hi, var, dt)
  if dt then
    if lo == -INF and hi == INF then return "all " .. var end
    if lo == -INF then return var .. " ≤ " .. fmtR(hi - 1) end
    if hi == INF then return var .. " ≥ " .. fmtR(lo) end
    if hi - lo == 1 then return var .. " = " .. fmtR(lo) end
    return fmtR(lo) .. " ≤ " .. var .. " ≤ " .. fmtR(hi - 1)
  end
  if lo == -INF and hi == INF then return "all " .. var end
  if lo == -INF then return var .. " < " .. fmtR(hi) end
  if hi == INF then return var .. " > " .. fmtR(lo) end
  return fmtR(lo) .. " < " .. var .. " < " .. fmtR(hi)
end

-- piecewise lines
function Sig:piecewise(var)
  local lines = {}
  for i, S in ipairs(self.seg) do
    local lo, hi = self.bp[i - 1] or -INF, self.bp[i] or INF
    lines[#lines + 1] = {expr = segStr(S, var, self.dt), range = rangeStr(lo, hi, var, self.dt), lo = lo, hi = hi}
  end
  return lines
end

-- compact u() form
function Sig:compact(var)
  var = var or (self.dt and "n" or "t")
  local dt = self.dt
  if dt then
    local lo, hi = self:dtSupport()
    if lo and hi - lo <= 12 then
      local parts = {}
      for n = lo, hi do
        local v = self:eval(n)
        if cabs(v) > 1e-12 then parts[#parts + 1] = cmul(v, "δ[" .. shiftStr(var, n) .. "]") end
      end
      if #parts > 0 then return joinsum(parts) end
    end
  end
  local parts = {}
  -- window form
  for i, S in ipairs(self.seg) do
    if #S > 0 then
      local lo, hi = self.bp[i - 1] or -INF, self.bp[i] or INF
      local w = windowStr(lo, hi, var, dt)
      if dt and lo ~= -INF and hi ~= INF and hi - lo == 1 then
        parts[#parts + 1] = cmul(segEval(S, lo, true), w)
      else
        parts[#parts + 1] = attach(segStr(S, var, dt), w)
      end
    end
  end
  local winForm = joinsum(parts)
  -- step form: S1 + sum (S_i - S_{i-1}) u(t - b_i)
  local sp = {}
  local ok = true
  if #self.seg[1] > 0 then ok = false end
  if ok then
    for i = 2, #self.seg do
      local diff = {}
      for _, m in ipairs(self.seg[i]) do diff[#diff + 1] = M(m.c, m.k, m.s) end
      for _, m in ipairs(self.seg[i - 1]) do diff[#diff + 1] = M(-m.c, m.k, m.s) end
      diff = mergeS(diff, dt)
      if #diff > 0 then
        local b = self.bp[i - 1]
        sp[#sp + 1] = attach(segStr(diff, var, dt), uStr(var, b, dt))
      end
    end
  end
  local stepForm = ok and joinsum(sp) or nil
  local best = winForm
  if stepForm and #self.seg > 2 and ulen(stepForm) < ulen(winForm) then best = stepForm end
  if #self.seg == 1 and #self.seg[1] > 0 then best = segStr(self.seg[1], var, dt) end
  if #parts == 0 then best = "" end
  local ip = {}
  for _, im in ipairs(self.imps) do ip[#ip + 1] = impStr(im, var) end
  if #ip > 0 then
    if best == "" then best = joinsum(ip) else best = best .. " + " .. joinsum(ip):gsub("^%-", "- ", 1) end
    best = best:gsub("%+ %- ", "- ")
  end
  if best == "" then best = "0" end
  return best
end

-- finite DT support listing
function Sig:dtSupport()
  if not self.dt then return nil end
  if #self.bp == 0 or #self.seg[1] > 0 or #self.seg[#self.seg] > 0 then return nil end
  local lo, hi = self.bp[1], self.bp[#self.bp] - 1
  if not lo or hi - lo > 40 then return nil end
  return lo, hi
end
function Sig:dtList()
  local lo, hi = self:dtSupport()
  if not lo then return nil end
  local vals = {}
  for n = lo, hi do vals[#vals + 1] = {n = n, v = self:eval(n)} end
  return vals
end
-- ===================== POLYNOMIALS (complex coeffs, index 0..deg) =====================
Poly = {}
function pnew(t) local p = {}; for i = 0, #t - 1 do p[i] = tocx(t[i + 1]) end; p.n = #t - 1; return p end
function pdeg(p)
  local d = p.n
  while d > 0 and cabs(p[d]) < 1e-12 do d = d - 1 end
  return d
end
function ptrim(p) local d = pdeg(p); local r = {n = d}; for i = 0, d do r[i] = p[i] end; return r end
function pconst(c) return {[0] = tocx(c), n = 0} end
function padd(a, b)
  local n = math.max(a.n, b.n); local r = {n = n}
  for i = 0, n do r[i] = (a[i] or ZERO) + (b[i] or ZERO) end
  return ptrim(r)
end
function pscale(a, c) local r = {n = a.n}; for i = 0, a.n do r[i] = a[i] * c end; return r end
function pmul(a, b)
  local r = {n = a.n + b.n}
  for i = 0, r.n do r[i] = ZERO end
  for i = 0, a.n do for j = 0, b.n do r[i + j] = r[i + j] + a[i] * b[j] end end
  return ptrim(r)
end
function peval(p, x)
  local v = ZERO
  for i = p.n, 0, -1 do v = v * x + p[i] end
  return v
end
function pdivmod(a, b)
  a = ptrim(a); b = ptrim(b)
  local db = b.n
  if db == 0 then return pscale(a, ONE / b[0]), pconst(0) end
  if a.n < db then return pconst(0), a end
  local r = {n = a.n}; for i = 0, a.n do r[i] = a[i] end
  local q = {n = math.max(0, a.n - db)}; for i = 0, q.n do q[i] = ZERO end
  for i = a.n - db, 0, -1 do
    local c = r[i + db] / b[db]
    q[i] = c
    for j = 0, db do r[i + j] = r[i + j] - c * b[j] end
  end
  r.n = math.max(0, db - 1)
  return ptrim(q), ptrim(r)
end
function piszero(p) p = ptrim(p); return p.n == 0 and cabs(p[0]) < 1e-12 end
-- roots via Durand-Kerner + Newton polish
function proots(p)
  p = ptrim(p)
  local n = p.n
  if n == 0 then return {} end
  local lead = p[n]
  local a = {}; for i = 0, n do a[i] = p[i] / lead end
  local mon = {n = n}; for i = 0, n do mon[i] = a[i] end
  if n == 1 then return {-a[0]} end
  if n == 2 then
    local b, c = a[1], a[0]
    local d = csqrt(b * b - c * 4)
    return {(-b + d) / 2, (-b - d) / 2}
  end
  local R = 1; for i = 0, n - 1 do R = math.max(R, 1 + cabs(a[i])) end
  local z = {}
  local seed = cx(0.4, 0.9)
  for i = 1, n do z[i] = cpow(seed, i - 1) * (R * 0.5) end
  for _ = 1, 500 do
    local maxd = 0
    for i = 1, n do
      local num = peval(mon, z[i])
      local den = ONE
      for j = 1, n do if j ~= i then den = den * (z[i] - z[j]) end end
      if cabs(den) < 1e-300 then den = cx(1e-12) end
      local dz = num / den
      z[i] = z[i] - dz
      maxd = math.max(maxd, cabs(dz))
    end
    if maxd < 1e-14 then break end
  end
  for i = 1, n do z[i] = ccln(z[i]); if abs(z[i].im) < 1e-9 then z[i] = cx(z[i].re) end end
  return z
end
-- group roots with multiplicities
function rootGroups(roots)
  local g = {}
  for _, r in ipairs(roots) do
    local f = false
    for _, x in ipairs(g) do
      if cabs(x.p - r) < 1e-5 * math.max(1, cabs(r)) then
        x.sum = x.sum + r; x.m = x.m + 1; x.p = x.sum / x.m; f = true; break
      end
    end
    if not f then g[#g + 1] = {p = r, m = 1, sum = r} end
  end
  for _, x in ipairs(g) do
    x.p = ccln(x.p)
    local rr = ratapprox(x.p.re, 1000, 1e-7); if rr then local p, q = ratapprox(x.p.re, 1000, 1e-7); x.p = cx(p / q, x.p.im) end
    local ii = ratapprox(x.p.im, 1000, 1e-7); if ii then local p, q = ratapprox(x.p.im, 1000, 1e-7); x.p = cx(x.p.re, p / q) end
  end
  return g
end
-- partial fractions of R/D (deg R < deg D). Returns list {p, m, A}: A/(v-p)^m
function pfe(R, D)
  D = ptrim(D)
  local groups = rootGroups(proots(D))
  local lead = D[D.n]
  local basis = {}
  for _, g in ipairs(groups) do
    for j = 1, g.m do
      -- polynomial D/(v-p)^j
      local q = pscale(D, ONE / lead)
      for _ = 1, j do q = pdivmod(q, pnew({-g.p, 1})) end
      basis[#basis + 1] = {p = g.p, m = j, poly = q}
    end
  end
  local n = D.n
  local A, b = {}, {}
  for row = 0, n - 1 do
    A[row + 1] = {}
    for c, bs in ipairs(basis) do A[row + 1][c] = bs.poly[row] or ZERO end
    b[row + 1] = (R[row] or ZERO) / lead
  end
  local x = solveLin(A, b)
  local out = {}
  for c, bs in ipairs(basis) do
    local v = ccln(x[c])
    if cabs(v) > 1e-11 then out[#out + 1] = {p = bs.p, m = bs.m, A = v} end
  end
  return out
end
-- polynomial in q=e^{-jΩ} or z^-1, ascending
function pstrQ(p, mode, off)
  p = ptrim(p); off = off or 0
  local terms = {}
  for i = 0, p.n do
    local c = p[i]
    if cabs(c) > 1e-12 then
      local d = i + off
      local vs
      if d == 0 then vs = ""
      elseif mode == "Z" then vs = "z^" .. (d > 0 and "-" or "") .. abs(d)
      else vs = "e^(" .. (d > 0 and "-" or "") .. "j" .. (abs(d) == 1 and "" or abs(d)) .. "Ω)" end
      terms[#terms + 1] = cmul(c, vs)
    end
  end
  return joinsum(terms)
end
-- polynomial printer in variable string (e.g. "jω", "s", "z^-1")
function pstr(p, v, wrapv)
  p = ptrim(p)
  local terms = {}
  for i = p.n, 0, -1 do
    local c = p[i]
    if cabs(c) > 1e-12 then
      local vs
      if i == 0 then vs = "" elseif i == 1 then vs = v else vs = (wrapv and ("(" .. v .. ")") or v) .. "^" .. i end
      terms[#terms + 1] = cmul(c, vs)
    end
  end
  return joinsum(terms)
end
-- ===================== OUTPUT COLLECTOR =====================
Out = {}
Out.__index = Out
function newOut() return setmetatable({lines = {}, plots = {}}, Out) end
function Out:h(s) self.lines[#self.lines + 1] = {k = "h", s = s} end
function Out:t(s) self.lines[#self.lines + 1] = {k = "t", s = s} end
function Out:r(s) self.lines[#self.lines + 1] = {k = "r", s = s} end
function Out:e(s) self.lines[#self.lines + 1] = {k = "e", s = s} end
function Out:plot(p) self.plots[#self.plots + 1] = p end

-- ===================== FORWARD TRANSFORMS (CT) =====================
-- F term: c * e^{-d v} * (v - p)^(-m)   (m may be negative => (v)^|m| with p = 0)
-- freq impulse: {c, w0} -> c δ(ω - w0)
function addFT(list, c, d, p, m)
  c = tocx(c)
  if cabs(c) < 1e-13 then return end
  for _, t in ipairs(list) do
    if abs(t.d - d) < 1e-10 and t.m == m and near(t.p, p, 1e-9) then t.c = t.c + c; return end
  end
  list[#list + 1] = {c = c, d = d, p = tocx(p), m = m}
end
function cleanFT(list)
  local r = {}
  for _, t in ipairs(list) do if cabs(t.c) > 1e-11 then t.c = ccln(t.c); t.d = cln(t.d); r[#r + 1] = t end end
  table.sort(r, function(a, b) if abs(a.d - b.d) > 1e-12 then return a.d < b.d end; if a.m ~= b.m then return a.m < b.m end; return a.p.re < b.p.re end)
  return r
end

-- describe a monomial piece c t^k e^{st} on (lo,hi)
function pieceStr(c, k, s, lo, hi, var, dt)
  local S = {M(c, k, s)}
  local e = segStr(S, var, dt)
  return e .. ",  " .. rangeStr(lo, hi, var, dt)
end

-- CT forward transform. mode "F" (Fourier) or "L" (Laplace)
-- returns terms, impulses, roc{lo,hi}, steps(list of strings)
function ctTransform(sig, mode)
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
function vfac(p, mode)
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
function delayStr(d, mode)
  if abs(d) < 1e-12 then return "" end
  if mode == "F" then return "e^(" .. ((d > 0) and "-" or "") .. "j" .. ((abs(d) == 1) and "" or fmtR(abs(d))) .. "ω)" end
  if mode == "L" then return "e^(" .. expoArg(-d, "s") .. ")" end
end
function ftermStr(t, mode)
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
function fimpStr(im, isDT)
  local w = isDT and "Ω" or "ω"
  local a = (abs(im.w0) < 1e-12) and w or (w .. (im.w0 > 0 and " - " or " + ") .. fmtR(abs(im.w0)))
  local s = cmul(im.c, "δ(" .. a .. ")")
  return s
end
-- combine terms with same delay into single rational N/D (in v)
function combineFT(terms)
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
function ratStr(N, D, mode, offN)
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

function ftSumStr(terms, imps, mode, isDT)
  local parts = {}
  for _, t in ipairs(terms) do parts[#parts + 1] = (isDT and dtermStr or ftermStr)(t, mode) end
  for _, im in ipairs(imps or {}) do parts[#parts + 1] = fimpStr(im, isDT) end
  return joinsum(parts)
end
-- evaluate CT transform at v (complex)
function ftEval(terms, v)
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
function addDT(list, c, d, p, m, dir)
  c = tocx(c); d = d + 0
  if cabs(c) < 1e-13 then return end
  if m == 0 then p = ZERO; dir = 1 end
  for _, t in ipairs(list) do
    if t.d == d and t.m == m and t.dir == dir and near(t.p, p, 1e-9) then t.c = t.c + c; return end
  end
  list[#list + 1] = {c = c, d = d, p = tocx(p), m = m, dir = dir}
end
-- coefficients beta_m with (x+a)^k = sum_{m=1}^{k+1} beta_m C(x+m-1, m-1)
function binomBasis(k, a)
  local A, b = {}, {}
  for x = 0, k do
    A[x + 1] = {}
    for m = 1, k + 1 do A[x + 1][m] = cx(binom(x + m - 1, m - 1)) end
    b[x + 1] = cx((x + a) ^ k)
  end
  return solveLin(A, b)
end
function dtTransform(sig, mode)
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
function dtEval(terms, Om, zmode)
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
function combineDT(terms)
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
-- ===================== FREQUENCY-DOMAIN EXPRESSIONS & INVERSE =====================
-- element: {c, d, N, D, P, win={lo,hi}, imp=w0}
-- mode: "F" (ω, v=jω), "L" (s), "D" (Ω, q=e^{-jΩ}), "Z" (z, q=z^-1)
function fe(c) return {c = tocx(c), d = 0, N = pconst(1), D = pconst(1), P = pconst(1)} end
function feCopy(e)
  return {c = e.c, d = e.d, N = e.N, D = e.D, P = e.P, win = e.win and {lo = e.win.lo, hi = e.win.hi}, imp = e.imp}
end
FVAR = {F = "w", L = "s", D = "W", Z = "z"}

-- evaluate element at frequency w (real for F/D) ignoring impulse/window
function feEval(e, mode, w)
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

function feMul(a, b, mode)
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
function listMul(A, B, mode)
  local out = {}
  for _, a in ipairs(A) do for _, b in ipairs(B) do
    local r = feMul(a, b, mode)
    if r and cabs(r.c) > 1e-14 then out[#out + 1] = r end
  end end
  return out
end
function listScale(A, c) local out = {}; for _, a in ipairs(A) do local e = feCopy(a); e.c = e.c * c; out[#out + 1] = e end; return out end
-- shift integer delay into polynomial (DT)
function absorbQ(e)
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
function listToElem(A, mode)
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

function feList(node, mode, env)
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
function binomPolyN(m)
  local p = pconst(1)
  for i = 1, m - 1 do p = pmul(p, pnew({i, 1})) end
  return pscale(p, cx(1 / fact(m - 1)))
end

-- windowed synthesis integral: (1/2π)∫_{lo}^{hi} P(w) e^{j w (t - d)} dw  ; returns list {c,W,pw,d}
function synthTerms(Pw, lo, hi, d, coef)
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
function synthStr(terms, var)
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
function inverseXform(src, mode, env, roc, out)
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
function synthEval(terms, t)
  local v = ZERO
  for _, s in ipairs(terms) do
    local tau = t - s.d
    if abs(tau) < 1e-9 then tau = 1e-9 end
    v = v + s.c * cexp(cx(0, s.W * tau)) / cx(tau ^ s.pw)
  end
  return v
end
-- ===================== CONTEXT (definitions & params) =====================
function buildCtx(defsText)
  local ctx = {defs = {}, params = {}, fnames = {}, errs = {}}
  -- first pass: names
  local lines = {}
  for line in (defsText .. "\n"):gmatch("([^\n;]*)[\n;]") do
    line = line:gsub("^%s+", ""):gsub("%s+$", "")
    if line ~= "" then lines[#lines + 1] = line end
  end
  for _, line in ipairs(lines) do
    local lhs = splitEq(line)
    if lhs then
      local nm = lhs:match("^%s*([%a][%w]*)%s*[%(%[]")
      if nm then ctx.fnames[nm] = true end
    end
  end
  ctx.env = mkenv({})
  for _, line in ipairs(lines) do
    local lhs, rhs = splitEq(line)
    if lhs then
      local nm, br, p = lhs:match("^%s*([%a][%w]*)%s*([%(%[])%s*([%a][%w]*)%s*[%)%]]%s*$")
      local ok, err = pcall(function()
        if nm then
          local pre = preprocess(p)
          local body = parse(rhs, {fnames = ctx.fnames})
          ctx.defs[nm] = {param = pre, body = body, dt = (br == "["), src = rhs}
          ctx.env.fns[nm] = {param = pre, body = body}
        else
          local pn = preprocess(lhs):gsub("%s", "")
          local v = evalN(parse(rhs, {fnames = ctx.fnames}), ctx.env)
          ctx.params[pn] = v
          ctx.env.vars[pn] = v
        end
      end)
      if not ok then ctx.errs[#ctx.errs + 1] = line .. ": " .. tostring(err) end
    end
  end
  return ctx
end
function P(ctx, s, extra)
  local fn = {}
  for k in pairs(ctx.fnames) do fn[k] = true end
  for k in pairs(extra or {}) do fn[k] = true end
  return parse(s, {fnames = fn})
end
-- decide DT vs CT for an expression
function isDTexpr(node, ctx)
  if hasvar(node, "n") and not hasvar(node, "t") then return true end
  local dt = false
  contains(node, function(nd)
    if nd.k == "call" and nd.br == "[" then dt = true end
    if nd.k == "call" and ctx.defs[nd.name] and ctx.defs[nd.name].dt then dt = true end
    return false
  end)
  if hasvar(node, "t") then return false end
  return dt
end
function toSig(ctx, node, dt)
  local var = dt and "n" or "t"
  local cv = mkConv(var, dt, ctx.env, ctx.defs)
  return cv:sig(node)
end
-- numeric function of a node (for plotting/fallback)
function numFn(ctx, node, var)
  return function(x)
    ctx.env.vars[var] = cx(x)
    ctx.env.tvar = var
    local ok, v = pcall(evalN, node, ctx.env)
    if ok then return v end
    return nil
  end
end
-- ===================== NUMERIC SELF-CHECKS =====================
-- Each transform answer is compared with a direct numeric computation at a test point.
function okLine(out, label, a, b)
  local good = near(a, b, 1e-5) or cabs(a - b) < 1e-7
  if good then
    out:t("Check at " .. label .. " computed directly = " .. decstr(cln(a.re)) .. (abs(a.im) > 1e-9 and ((a.im < 0) and " - j" or " + j") .. decstr(abs(cln(a.im))) or "") .. ", matches  (OK)")
  else
    out:e("Check FAILED at " .. label .. " direct = " .. fmtC(a) .. ", formula = " .. fmtC(b) .. "  (please report this input)")
  end
  return good
end

-- numeric ∫ x(t) e^{-s t} dt over the signal pieces, plus impulses
function numLaplace(sig, s0)
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
function numZ(sig, z0)
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
-- ===================== TOOLS (part A) =====================
function sampleIn(lo, hi, dt)
  if lo == -INF and hi == INF then return 0 end
  if lo == -INF then return hi - 1 end
  if hi == INF then return dt and lo or lo + 1 end
  return dt and lo or (lo + hi) / 2
end

TOOLS = {}
function addTool(t) TOOLS[#TOOLS + 1] = t end
function V(dt) return dt and "n" or "t" end
function sigLines(out, sig, var, name)
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
function periodInfo(sig, out)
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
function evenOdd(sig, out)
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

function conjSig(sig)
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
function shiftedStr(S, var, dv, dt)
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
function convSteps(X, H, out)
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
-- ===================== TOOLS (part B) =====================
-- ---------- numeric system tester ----------
function sysRunner(ctx, node, dt)
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
function rnd(seed) local s = seed; return function() s = (s * 16807) % 2147483647; return s / 2147483647 end end

function sysTests(ctx, node, dt)
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

function xArgs(node)
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
function boundOf(n, ctx, var)
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
function exactY(ctx, node, dt, inputSrc)
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
function xArgList(node, var, ctx)
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
function absSig(sig)
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
function ltiProps(H, out)
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
function fwdTool(name, mode, dtMode, help, def)
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

function invTool(name, mode, help, def, rocLabel)
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
-- ===================== TOOLS (part C) =====================
-- linear coefficient extraction from difference/differential equation side
-- returns {Y = poly, X = poly} in v (CT: jω power, DT: q = e^{-jΩ} power)
function deSide(node, dt, env)
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
function reImPoly(p, dt)
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
function magPhase(out, Hterms, dt)
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
-- ===================== TOOLS (part D) =====================
-- CT signal sampled at t = n Ts -> DT signal
function sampleSig(X, Ts)
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
-- ===================== EXAMPLE PROBLEMS =====================
-- {cat, title, tool, vals, defs(optional, ';' separated), keys(optional)}
EX = {
  -- Complex numbers
  {"Complex", "2e^(jπ) to Cartesian", "Complex numbers", {"2e^(jpi)"}},
  {"Complex", "e^(j2π) to Cartesian", "Complex numbers", {"e^(j2pi)"}},
  {"Complex", "2e^(j3π/4) to Cartesian", "Complex numbers", {"2e^(j3pi/4)"}},
  {"Complex", "e^(jπ/2)+e^(jπ)", "Complex numbers", {"e^(jpi/2)+e^(jpi)"}},
  {"Complex", "(3+j)/(1-j) divide", "Complex numbers", {"(3+j)/(1-j)"}},
  {"Complex", "-1-j to polar", "Complex numbers", {"-1-j"}},
  {"Complex", "(-2+j)/(1+j) polar", "Complex numbers", {"(-2+j)/(1+j)"}},
  -- Running integral / sum
  {"Running integral", "∫ u(τ) → ramp", "Running integral / sum", {"u(t)"}},
  {"Running integral", "∫ e^(-τ)u(τ)", "Running integral / sum", {"e^(-t)u(t)"}},
  {"Running integral", "∫ sin(2τ)u(τ)", "Running integral / sum", {"sin(2t)u(t)"}},
  {"Running integral", "∫ δ(t+2)-δ(t-2)", "Running integral / sum", {"d(t+2)-d(t-2)"}},
  {"Running integral", "∫ t[δ(t-1)+δ(t-2)+δ(t-3)]", "Running integral / sum", {"t(d(t-1)+d(t-2)+d(t-3))"}},
  {"Running sum", "Σ u[k]", "Running integral / sum", {"u[n]"}},
  {"Running sum", "Σ (-1)^k u[k]", "Running integral / sum", {"(-1)^n u[n]"}},
  {"Running sum", "Σ (u[k]-u[k-5])cos(πk/2)", "Running integral / sum", {"(u[n]-u[n-5])cos(pi/2 n)"}},
  {"Running sum", "Σ (1/2)^k u[k]", "Running integral / sum", {"(1/2)^n u[n]"}},
  {"Running sum", "Σ 2^k (two-sided)", "Running integral / sum", {"2^n"}},
  -- Signal transforms / sketches
  {"Signal sketch", "x(t-1) of piecewise x", "Signal analyzer", {"x(t-1)"}, "x(t)=(t+2)(-2<=t)(t<0)+2(0<=t)(t<1)+(2-t)(1<=t)(t<2)"},
  {"Signal sketch", "x(2t+1) scale & shift", "Signal analyzer", {"x(2t+1)"}, "x(t)=(t+2)(-2<=t)(t<0)+2(0<=t)(t<1)+(2-t)(1<=t)(t<2)"},
  {"Signal sketch", "[x(t)+x(-t)]u(t)", "Signal analyzer", {"(x(t)+x(-t))u(t)"}, "x(t)=(t+2)(-2<=t)(t<0)+2(0<=t)(t<1)+(2-t)(1<=t)(t<2)"},
  {"Signal sketch", "x(t)[δ(t+3/2)-δ(t-3/2)]", "Signal analyzer", {"x(t)(d(t+3/2)-d(t-3/2))"}, "x(t)=(t+2)(-2<=t)(t<0)+2(0<=t)(t<1)+(2-t)(1<=t)(t<2)"},
  {"Signal sketch", "DT x[n-4]", "Signal analyzer", {"x[n-4]"}, "x[n]=-d[n+2]+d[n+1]+d[n]+d[n-1]+1/2d[n-2]+1/2d[n-3]"},
  {"Signal sketch", "DT x[3n] decimation", "Signal analyzer", {"x[3n]"}, "x[n]=-d[n+2]+d[n+1]+d[n]+d[n-1]+1/2d[n-2]+1/2d[n-3]"},
  {"Signal sketch", "DT x[n]u[3-n]", "Signal analyzer", {"x[n]u[3-n]"}, "x[n]=-d[n+2]+d[n+1]+d[n]+d[n-1]+1/2d[n-2]+1/2d[n-3]"},
  {"Signal sketch", "DT x[n-2]δ[n-2]", "Signal analyzer", {"x[n-2]d[n-2]"}, "x[n]=-d[n+2]+d[n+1]+d[n]+d[n-1]+1/2d[n-2]+1/2d[n-3]"},
  {"Signal sketch", "First difference of δ, twice", "Signal analyzer", {"z[n]-z[n-1]"}, "z[n]=d[n]-d[n-1]"},
  {"Signal sketch", "LTI superposition y1(t)-y1(t-2)", "Signal analyzer", {"y1(t)-y1(t-2)"}, "y1(t)=t(0<=t)(t<1)+(2-t)(1<=t)(t<2)"},
  -- Periodicity / even-odd
  {"Periodic", "3cos(5t+π/2)", "Signal analyzer", {"3cos(5t+pi/2)"}},
  {"Periodic", "e^(-jπt)", "Signal analyzer", {"e^(-j pi t)"}},
  {"Periodic", "[cos(2t)]^2", "Signal analyzer", {"cos(2t)^2"}},
  {"Periodic", "e^(1+j2πt)", "Signal analyzer", {"e^(1+j2pi t)"}},
  {"Periodic", "e^((1+j2π)t) not periodic", "Signal analyzer", {"e^((1+j2pi)t)"}},
  {"Periodic", "DT sin(3πn/5)", "Signal analyzer", {"sin(3pi/5 n)"}},
  {"Periodic", "DT sin(13πn/5)", "Signal analyzer", {"sin(13pi/5 n)"}},
  {"Periodic", "DT cos(n) not periodic", "Signal analyzer", {"cos(n)"}},
  {"Periodic", "DT j^n", "Signal analyzer", {"j^n"}},
  {"Periodic", "sum of two sinusoids", "Signal analyzer", {"cos(2t)+sin(3t)"}},
  {"Even/odd", "DT cos(πn/3)", "Signal analyzer", {"cos(pi/3 n)"}},
  {"Even/odd", "DT (n-1)^2", "Signal analyzer", {"(n-1)^2"}},
  {"Even/odd", "DT u[n]-u[-n]", "Signal analyzer", {"u[n]-u[-n]"}},
  {"Even/odd", "even & odd parts of e^(-t)u(t)", "Signal analyzer", {"e^(-t)u(t)"}},
  {"Energy", "energy of e^(-2t)u(t)", "Signal analyzer", {"e^(-2t)u(t)"}},
  -- System properties
  {"System props", "y = e^(jt)x(t)", "System properties", {"y(t)=e^(jt)x(t)"}},
  {"System props", "y = ∫[t-1,t+1] x(τ)^2 dτ", "System properties", {"y(t)=int(x(tau)^2,tau,t-1,t+1)"}},
  {"System props", "y = ∫[-∞,t^2] x(τ)dτ", "System properties", {"y(t)=int(x(tau),tau,-inf,t^2)"}},
  {"System props", "y = dx/dt", "System properties", {"y(t)=der(x(t))"}},
  {"System props", "y[n] = x[-n]", "System properties", {"y[n]=x[-n]"}},
  {"System props", "y[n] = x[n-1]-2x[n]+x[n+1]", "System properties", {"y[n]=x[n-1]-2x[n]+x[n+1]"}},
  {"System props", "y[n] = x[n]u[n]", "System properties", {"y[n]=x[n]u[n]"}},
  {"System props", "y[n] = Σ 2^k x[n-k]^2", "System properties", {"y[n]=sum(2^k x[n-k]^2,k,-inf,inf)"}},
  {"System props", "y = e^(x(t)) invertible", "System properties", {"y(t)=e^(x(t))"}},
  {"System props", "y = (t^2+1)x(t)", "System properties", {"y(t)=(t^2+1)x(t)"}},
  {"System props", "y = x(t)-x(-t)", "System properties", {"y(t)=x(t)-x(-t)"}},
  {"System props", "y = x(2t) (time scaling)", "System properties", {"y(t)=x(2t)"}},
  {"System props", "y = 2x(t)+3 (affine)", "System properties", {"y(t)=2x(t)+3"}},
  {"System props", "y[n] = x[n]x[n-1]", "System properties", {"y[n]=x[n]x[n-1]"}},
  -- Convolution
  {"Convolution", "α^n u[n] * β^n u[n]", "Convolution", {"a^n u[n]", "b^n u[n]"}, "a=1/2;b=1/3"},
  {"Convolution", "α^n u[n] * α^n u[n]", "Convolution", {"a^n u[n]", "a^n u[n]"}, "a=1/2"},
  {"Convolution", "pulse * anticausal 2^n", "Convolution", {"u[n]-u[n-5]", "2^n u[2-n]"}},
  {"Convolution", "e^(-αt)u(t) * e^(-βt)u(t)", "Convolution", {"e^(-a t)u(t)", "e^(-b t)u(t)"}, "a=1;b=2"},
  {"Convolution", "e^(-αt)u(t) * e^(-αt)u(t)", "Convolution", {"e^(-a t)u(t)", "e^(-a t)u(t)"}, "a=1"},
  {"Convolution", "square wave piece * ramp-down", "Convolution", {"x(t)", "h(t)"}, "x(t)=rect(t)-rect(t-1);h(t)=(1-t)(u(t)-u(t-1))"},
  {"Convolution", "u(t)*u(t)*u(t)", "Convolution", {"conv(u(t),u(t))", "u(t)"}},
  {"Convolution", "two rectangles → trapezoid", "Convolution", {"u(t)-u(t-2)", "u(t)-u(t-1)"}},
  {"Convolution", "rect * rect → triangle", "Convolution", {"rect(t)", "rect(t)"}},
  {"Convolution", "impulse train * triangle", "Convolution", {"d(t)+d(t-3)+d(t+3)", "tri(t)"}},
  {"Convolution", "cos(t) * e^(-at)u(t)", "Convolution", {"cos(t)", "e^(-a t)u(t)"}, "a=2"},
  {"Convolution", "e^(-at)u(t) * e^(bt)u(-t)", "Convolution", {"e^(-a t)u(t)", "e^(b t)u(-t)"}, "a=1;b=3"},
  {"Convolution", "DT finite sequences", "Convolution", {"d[n]+2d[n-1]+3d[n-2]", "d[n]-d[n-1]"}},
  -- LTI / step response
  {"LTI", "h from step s = e^(-2t)u(t), input t u(t)", "LTI system from h / step resp.", {"", "e^(-2t)u(t)", "tu(t)"}},
  {"LTI", "h[n]=(1/5)^n u[n]", "LTI system from h / step resp.", {"(1/5)^n u[n]", "", ""}},
  {"LTI", "h[n]=(1/2)^n u[-n]", "LTI system from h / step resp.", {"(1/2)^n u[-n]", "", ""}},
  {"LTI", "h[n]=(-1/2)^n u[n]+1.01^n u[n-1]", "LTI system from h / step resp.", {"(-1/2)^n u[n]+1.01^n u[n-1]", "", ""}},
  {"LTI", "h[n]=(-1/2)^n u[n]+1.01^n u[1-n]", "LTI system from h / step resp.", {"(-1/2)^n u[n]+1.01^n u[1-n]", "", ""}},
  {"LTI", "h(t)=e^(-4t)u(t-2)", "LTI system from h / step resp.", {"e^(-4t)u(t-2)", "", ""}},
  {"LTI", "h(t)=e^(-6t)u(3-t)", "LTI system from h / step resp.", {"e^(-6t)u(3-t)", "", ""}},
  {"LTI", "h(t)=e^(-6|t|)", "LTI system from h / step resp.", {"e^(-6abs(t))", "", ""}},
  {"LTI", "h(t)=t e^(-t)u(t)", "LTI system from h / step resp.", {"t e^(-t)u(t)", "", ""}},
  {"LTI", "cascade of two exponentials", "Interconnected systems", {"e^(-t)u(t)", "e^(-2t)u(t)", "cascade"}},
  {"LTI", "parallel systems", "Interconnected systems", {"d(t)", "-e^(-t)u(t)", "parallel"}},
  -- Integrals with impulses
  {"Integrals", "∫ δ(t)cos(t) dt (sifting)", "Integral / sum evaluator", {"d(t)cos(t)", "-inf", "inf"}},
  {"Integrals", "∫[0,5] sin(2πt)δ(t+3) dt", "Integral / sum evaluator", {"sin(2pi t)d(t+3)", "0", "5"}},
  {"Integrals", "∫[-5,5] u(1-τ)cos(2πτ) dτ", "Integral / sum evaluator", {"u(1-t)cos(2pi t)", "-5", "5"}},
  {"Integrals", "∫ e^(-6|t|) dt", "Integral / sum evaluator", {"e^(-6abs(t))", "-inf", "inf"}},
  {"Integrals", "Σ (1/2)^n u[n]", "Integral / sum evaluator", {"(1/2)^n u[n]", "-inf", "inf"}},
  -- CTFT
  {"Fourier (CT)", "e^(-αt)cos(ω0 t)u(t)", "Fourier transform (CT)", {"e^(-a t)cos(w0 t)u(t)"}, "a=1;w0=3"},
  {"Fourier (CT)", "e^(-3|t|)sin(2t)", "Fourier transform (CT)", {"e^(-3abs(t))sin(2t)"}},
  {"Fourier (CT)", "(1+cos πt) for |t|≤1", "Fourier transform (CT)", {"(1+cos(pi t))(abs(t)<=1)"}},
  {"Fourier (CT)", "impulse train Σ α^k δ(t-kT) (5 terms)", "Fourier transform (CT)", {"d(t)+a d(t-T)+a^2 d(t-2T)+a^3 d(t-3T)+a^4 d(t-4T)"}, "a=1/2;T=1"},
  {"Fourier (CT)", "(t^2+t+1)e^(-at)u(t)", "Fourier transform (CT)", {"(t^2+t+1)e^(-a t)u(t)"}, "a=2"},
  {"Fourier (CT)", "t for |t|<1", "Fourier transform (CT)", {"t(abs(t)<1)"}},
  {"Fourier (CT)", "rect(t/2) → sinc", "Fourier transform (CT)", {"rect(t/2)"}},
  {"Fourier (CT)", "triangle tri(t)", "Fourier transform (CT)", {"tri(t)"}},
  {"Fourier (CT)", "cos(3t)+1 (impulses)", "Fourier transform (CT)", {"1+cos(3t)"}},
  {"Fourier (CT)", "u(t) (with πδ)", "Fourier transform (CT)", {"u(t)"}},
  {"Inverse FT", "cos(2ω-π/4)", "Inverse Fourier (CT)", {"cos(2w-pi/4)"}},
  {"Inverse FT", "3[δ(ω-2)-δ(ω+2)]+5[δ(ω-5)+δ(ω+5)]", "Inverse Fourier (CT)", {"3(d(w-2)-d(w+2))+5(d(w-5)+d(w+5))"}},
  {"Inverse FT", "jω for |ω|<1", "Inverse Fourier (CT)", {"jw(abs(w)<1)"}},
  {"Inverse FT", "2πδ(ω)+πδ(ω∓3)-(3π/j)δ(ω-4)...", "Inverse Fourier (CT)", {"2pi d(w)+pi d(w-3)+pi d(w+3)-3pi/j d(w-4)+3pi/j d(w+4)"}},
  {"Inverse FT", "(1+e^(-jω)-2e^(-2jω))/(jω)", "Inverse Fourier (CT)", {"(1+e^(-jw)-2e^(-2jw))/(jw)"}},
  {"Inverse FT", "(jω)^2/(a+jω)^2", "Inverse Fourier (CT)", {"(jw)^2/(a+jw)^2"}, "a=2"},
  {"Inverse FT", "1/((a+jω)(b-jω)) two-sided", "Inverse Fourier (CT)", {"1/((a+jw)(b-jw))"}, "a=1;b=3"},
  {"Inverse FT", "ideal low-pass |ω|<W", "Inverse Fourier (CT)", {"(abs(w)<W0)"}, "W0=2"},
  {"Conv. property", "e^(-at)u(t)*e^(-at)u(t) both ways", "Convolution property check", {"e^(-a t)u(t)", "e^(-a t)u(t)"}, "a=1"},
  {"Conv. property", "e^(-at)u(t)*e^(bt)u(-t) both ways", "Convolution property check", {"e^(-a t)u(t)", "e^(b t)u(-t)"}, "a=1;b=3"},
  {"Conv. property", "DT a^n u[n]*b^n u[n] both ways", "Convolution property check", {"a^n u[n]", "b^n u[n]"}, "a=1/2;b=1/3"},
  {"Conv. property", "DT δ[n-n0]*δ[n-n1]", "Convolution property check", {"d[n-2]", "d[n-3]"}},
  -- Frequency response
  {"Freq. response", "y = dx/dt", "Frequency response H", {"y(t)=der(x(t))", ""}},
  {"Freq. response", "y = d²x/dt²", "Frequency response H", {"y(t)=der2(x(t))", ""}},
  {"Freq. response", "y = x(t)+x(t-1)", "Frequency response H", {"y(t)=x(t)+x(t-1)", ""}},
  {"Freq. response", "y = ∫[t-2,t] x(τ)dτ", "Frequency response H", {"y(t)=int(x(tau),tau,t-2,t)", ""}},
  {"Freq. response", "y'+2y = x'+x, x=e^(-t)u(t)", "Frequency response H", {"y'+2y=x'+x", "e^(-t)u(t)"}},
  {"Freq. response", "y''+3y'+2y = x", "Frequency response H", {"y''+3y'+2y=x", "u(t)"}},
  {"Freq. response", "DT y[n]-0.5y[n-1]=x[n]", "Frequency response H", {"y[n]-1/2y[n-1]=x[n]", "(1/3)^n u[n]"}},
  {"Freq. response", "DT moving average", "Frequency response H", {"y[n]=(x[n]+x[n-1]+x[n-2])/3", ""}},
  -- DTFT
  {"DTFT", "u[n-2]-u[n-6]", "DTFT", {"u[n-2]-u[n-6]"}},
  {"DTFT", "(1/2)^(-n)u[-n-1]", "DTFT", {"(1/2)^(-n)u[-n-1]"}},
  {"DTFT", "n for -3≤n≤3", "DTFT", {"n(abs(n)<=3)"}},
  {"DTFT", "(n+1)a^n u[n]", "DTFT", {"(n+1)(1/2)^n u[n]"}},
  {"DTFT", "cos(2πn/3) (impulses)", "DTFT", {"cos(2pi/3 n)"}},
  {"Inverse DTFT", "band-pass π/4≤|Ω|≤3π/4", "Inverse DTFT", {"(abs(W)>=pi/4)(abs(W)<=3pi/4)"}},
  {"Inverse DTFT", "1+3e^(-jΩ)+2e^(-j2Ω)-4e^(-j3Ω)+e^(-j10Ω)", "Inverse DTFT", {"1+3e^(-jW)+2e^(-2jW)-4e^(-3jW)+e^(-10jW)"}},
  {"Inverse DTFT", "e^(-jΩ/2) on [-π,π]", "Inverse DTFT", {"e^(-jW/2)"}},
  {"Inverse DTFT", "(e^(-jΩ)-1/5)/(1-e^(-jΩ)/5)", "Inverse DTFT", {"(e^(-jW)-1/5)/(1-1/5e^(-jW))"}},
  {"Inverse DTFT", "partial fractions 2 poles", "Inverse DTFT", {"1/((1-1/2e^(-jW))(1-1/3e^(-jW)))"}},
  {"Inverse DTFT", "inverse system G=1/H", "Inverse DTFT", {"(1-1/2e^(-jW))/(1-1/3e^(-jW))"}},
  -- Laplace & Z
  {"Laplace", "e^(-2t)u(t)+e^t u(-t) (ROC strip)", "Laplace transform", {"e^(-2t)u(t)+e^(t)u(-t)"}},
  {"Laplace", "t e^(-2t)u(t-1)", "Laplace transform", {"t e^(-2t)u(t-1)"}},
  {"Laplace", "sin(3t)u(t)", "Laplace transform", {"sin(3t)u(t)"}},
  {"Inverse Laplace", "repeated pole (s+3)/((s+1)(s+2)^2)", "Inverse Laplace", {"(s+3)/((s+1)(s+2)^2)", ""}},
  {"Inverse Laplace", "complex poles (s+1)/(s^2+2s+5)", "Inverse Laplace", {"(s+1)/(s^2+2s+5)", ""}},
  {"Inverse Laplace", "two-sided with ROC -2<Re<-1", "Inverse Laplace", {"1/((s+1)(s+2))", "-2<x<-1"}},
  {"Inverse Laplace", "delay e^(-2s)/(s+1)", "Inverse Laplace", {"e^(-2s)/(s+1)", ""}},
  {"Z-transform", "(1/2)^n u[n] - 2^n u[-n-1]", "Z-transform", {"(1/2)^n u[n] - 2^n u[-n-1]"}},
  {"Inverse Z", "ROC 1/2<|z|<2", "Inverse Z", {"1/((1-1/2z^-1)(1-2z^-1))", "1/2<x<2"}},
  {"Inverse Z", "z/(z-1/2) causal", "Inverse Z", {"z/(z-1/2)", ""}},
  -- Fourier series
  {"Fourier series", "cos(2πt)+sin(4πt), T=1", "Fourier series (CT/DT)", {"cos(2pi t)+sin(4pi t)", "1", "3"}},
  {"Fourier series", "square pulse |t|<1, T=4", "Fourier series (CT/DT)", {"(abs(t)<1)", "4", "5"}},
  {"Fourier series", "sawtooth t on [0,1), T=1", "Fourier series (CT/DT)", {"t(0<=t)(t<1)", "1", "4"}},
  {"Fourier series", "DT pulse u[n]-u[n-3], N=5", "Fourier series (CT/DT)", {"u[n]-u[n-3]", "5", ""}},
  -- Sampling
  {"Sampling", "aliasing 3 Hz + 7 Hz at fs=10", "Sampling & aliasing", {"cos(2pi*3t)+sin(2pi*7t)", "fs=10"}},
  {"Sampling", "no aliasing cos(2t) at Ts=0.5", "Sampling & aliasing", {"cos(2t)", "0.5"}},
}
EXAMPLES = {}
for _, e in ipairs(EX) do
  EXAMPLES[#EXAMPLES + 1] = {cat = e[1], title = e[2], tool = e[3], vals = e[4], defs = e[5]}
end
-- ===================== REFERENCE TABLES =====================
REFS = {
  {title = "Basic signals & impulse rules", lines = {
    "u(t) = 1 for t>0, 0 for t<0;   u[n] = 1 for n≥0",
    "δ(t): area 1 at t=0;  δ[n] = 1 only at n=0",
    "du(t)/dt = δ(t);   ∫[-∞,t] δ(τ)dτ = u(t);   u[n]-u[n-1] = δ[n]",
    "Sifting: ∫ x(t)δ(t-t0)dt = x(t0);  x(t)δ(t-t0) = x(t0)δ(t-t0)",
    "Scaling: δ(at) = δ(t)/|a|;   x(t)*δ(t-t0) = x(t-t0)",
    "Doublet: ∫ x(t)δ'(t-t0)dt = -x'(t0)",
    "rect(t) = 1 for |t|<1/2;  tri(t) = 1-|t| for |t|<1;  ramp r(t) = t u(t)",
    "x(at+b): shift by b first (left if b>0), then scale by a; or x(a(t+b/a))",
    "DT x[kn]: keep every k-th sample (decimation)",
  }},
  {title = "Even/odd, energy, power, periodicity", lines = {
    "Even: x(-t)=x(t).  Odd: x(-t)=-x(t)",
    "xe = [x(t)+x(-t)]/2,  xo = [x(t)-x(-t)]/2,  x = xe + xo",
    "Energy E = ∫|x|²dt (Σ|x[n]|²);  energy signal: 0<E<∞, P=0",
    "Power P = lim (1/2T)∫[-T,T]|x|²dt; periodic: P = (1/T0)∫[T0]|x|²dt",
    "CT: cos(ω0t+φ), e^(jω0t) periodic with T = 2π/|ω0|",
    "Sum of periodic CT signals periodic iff ratio of periods rational; T0 = LCM",
    "DT e^(jΩ0n) periodic iff Ω0/2π = k/N rational; N (reduced) = fundamental",
    "e^((σ+jω)t) with σ≠0 is NOT periodic",
  }},
  {title = "Complex numbers & identities", lines = {
    "Euler: e^(jθ) = cosθ + j sinθ;  |e^(jθ)| = 1",
    "cosθ = (e^(jθ)+e^(-jθ))/2;  sinθ = (e^(jθ)-e^(-jθ))/(2j)",
    "z = a+jb = r e^(jθ),  r = √(a²+b²),  θ = atan2(b,a)",
    "θ: a>0 → arctan(b/a); a<0 → arctan(b/a) ± π",
    "Division: (a+jb)/(c+jd) = (a+jb)(c-jd)/(c²+d²)",
    "cos²θ = (1+cos2θ)/2;  sin²θ = (1-cos2θ)/2;  sin2θ = 2sinθcosθ",
    "cos(a±b) = cos a cos b ∓ sin a sin b;  sin(a±b) = sin a cos b ± cos a sin b",
    "cos a cos b = [cos(a-b)+cos(a+b)]/2",
  }},
  {title = "Series & integrals", lines = {
    "Σ[n=0,N-1] a^n = (1-a^N)/(1-a)  (a≠1);  = N if a = 1",
    "Σ[n=k,∞] a^n = a^k/(1-a), |a|<1",
    "Σ[n=0,∞] n a^n = a/(1-a)²;  Σ (n+1)a^n = 1/(1-a)²",
    "∫ e^(at)dt = e^(at)/a;   ∫ t e^(at)dt = e^(at)(t/a - 1/a²)",
    "∫ t^k e^(at)dt = e^(at) Σ(-1)^i k!/(k-i)! t^(k-i)/a^(i+1)",
    "∫ e^(-at²)... : ∫[-∞,∞] e^(-t²+zt)dt = √π e^(z²/4)",
    "By parts: ∫u dv = uv - ∫v du",
  }},
  {title = "System properties (definitions & tests)", lines = {
    "Memoryless: y(t) depends only on x(t) at the same t",
    "Time-invariant: input x(t-t0) gives output y(t-t0). Test: z(t) vs y(t-t0)",
    "  red flags: t outside x(·), x(2t), x(-t), limits like t²",
    "Linear: a x1 + b x2 → a y1 + b y2 (homogeneity + additivity); zero in → zero out",
    "  red flags: x², |x|, e^x, + constant, x·x(t-1)",
    "Causal: y(t) uses x(τ) only for τ ≤ t. Future args (t+1, -t for t<0) → noncausal",
    "BIBO stable: |x|≤B ⇒ |y|≤C. Integrator ∫[-∞,t], derivative, t·x(t) are unstable",
    "Invertible: distinct inputs → distinct outputs. Not: x², d/dt, x(t)-x(-t), x[2n]",
    "LTI: causal ⇔ h(t)=0 for t<0;  stable ⇔ ∫|h|dt < ∞ (Σ|h[n]| < ∞)",
    "LTI memoryless ⇔ h = Kδ;  step response s(t) = ∫h, h = ds/dt (h[n]=s[n]-s[n-1])",
  }},
  {title = "Convolution properties", lines = {
    "y = x*h = ∫x(τ)h(t-τ)dτ  (Σ x[k]h[n-k])",
    "Commutative, associative, distributive",
    "x*δ(t-t0) = x(t-t0);   x*u(t) = ∫[-∞,t] x(τ)dτ",
    "Support: [a1,b1]*[a2,b2] → [a1+a2, b1+b2];  length L1+L2-1 for DT",
    "Area: ∫y = (∫x)(∫h);   derivative: (x*h)' = x'*h",
    "e^(-at)u * e^(-bt)u = (e^(-at)-e^(-bt))/(b-a) u(t);  a=b: t e^(-at)u(t)",
    "a^n u[n] * b^n u[n] = (b^(n+1)-a^(n+1))/(b-a) u[n];  a=b: (n+1)a^n u[n]",
    "u(t)*u(t) = t u(t);  k-fold: t^(k-1)/(k-1)! u(t)",
    "Cascade: h1*h2;  parallel: h1+h2",
  }},
  {title = "CT Fourier transform pairs", lines = {
    "X(jω) = ∫x(t)e^(-jωt)dt;  x(t) = (1/2π)∫X(jω)e^(jωt)dω",
    "δ(t) <-> 1;   1 <-> 2πδ(ω);   δ(t-t0) <-> e^(-jωt0)",
    "e^(jω0t) <-> 2πδ(ω-ω0);   cos ω0t <-> π[δ(ω-ω0)+δ(ω+ω0)]",
    "sin ω0t <-> (π/j)[δ(ω-ω0)-δ(ω+ω0)]",
    "u(t) <-> 1/(jω) + πδ(ω);   sgn(t) <-> 2/(jω)",
    "e^(-at)u(t) <-> 1/(a+jω), Re a>0",
    "t e^(-at)u(t) <-> 1/(a+jω)²;  t^(n-1)/(n-1)! e^(-at)u(t) <-> 1/(a+jω)^n",
    "e^(-a|t|) <-> 2a/(a²+ω²)",
    "rect: 1 for |t|<T1 <-> 2 sin(ωT1)/ω",
    "sin(Wt)/(πt) <-> 1 for |ω|<W (ideal low-pass)",
    "e^(-at)cos(ω0t)u(t) <-> (a+jω)/((a+jω)²+ω0²)",
    "e^(-at)sin(ω0t)u(t) <-> ω0/((a+jω)²+ω0²)",
    "Σ δ(t-kT) <-> (2π/T)Σ δ(ω-2πk/T)",
    "e^(-t²/(2σ²)) <-> σ√(2π) e^(-σ²ω²/2)",
  }},
  {title = "CT Fourier transform properties", lines = {
    "Linearity: ax+by <-> aX+bY",
    "Shift: x(t-t0) <-> e^(-jωt0)X(jω);   e^(jω0t)x(t) <-> X(j(ω-ω0))",
    "Scaling: x(at) <-> (1/|a|)X(jω/a);   x(-t) <-> X(-jω)",
    "Conjugate: x*(t) <-> X*(-jω); real x ⇒ X(-jω) = X*(jω), |X| even, ∠X odd",
    "Differentiation: dx/dt <-> jωX;   t x(t) <-> j dX/dω",
    "Integration: ∫[-∞,t]x <-> X/(jω) + πX(0)δ(ω)",
    "Convolution: x*h <-> X·H;   Multiplication: x·p <-> (1/2π) X*P",
    "Modulation: x(t)cos ω0t <-> [X(j(ω-ω0)) + X(j(ω+ω0))]/2",
    "Duality: if x(t) <-> X(jω) then X(jt) <-> 2π x(-ω)",
    "Parseval: ∫|x|²dt = (1/2π)∫|X|²dω",
    "LTI: x=e^(jωt) ⇒ y = H(jω)e^(jωt);  diff. eq: replace d^k/dt^k by (jω)^k",
  }},
  {title = "DTFT pairs", lines = {
    "X(e^(jΩ)) = Σ x[n]e^(-jΩn) (2π-periodic);  x[n] = (1/2π)∫[2π] X e^(jΩn)dΩ",
    "δ[n] <-> 1;   δ[n-n0] <-> e^(-jΩn0)",
    "a^n u[n] <-> 1/(1-a e^(-jΩ)), |a|<1",
    "(n+1)a^n u[n] <-> 1/(1-a e^(-jΩ))²",
    "-a^n u[-n-1] <-> 1/(1-a e^(-jΩ)), |a|>1",
    "1 for |n|≤N1 <-> sin(Ω(N1+1/2))/sin(Ω/2)",
    "sin(Wn)/(πn) <-> 1 for |Ω|<W (periodic)",
    "e^(jΩ0n) <-> 2π Σ δ(Ω-Ω0-2πl);   1 <-> 2π Σ δ(Ω-2πl)",
    "cos Ω0n <-> π Σ[δ(Ω-Ω0-2πl)+δ(Ω+Ω0-2πl)]",
    "u[n] <-> 1/(1-e^(-jΩ)) + π Σ δ(Ω-2πl)",
  }},
  {title = "DTFT properties", lines = {
    "Shift: x[n-n0] <-> e^(-jΩn0)X;   e^(jΩ0n)x[n] <-> X(e^(j(Ω-Ω0)))",
    "Reversal: x[-n] <-> X(e^(-jΩ));  conjugate symmetry for real x",
    "Difference: x[n]-x[n-1] <-> (1-e^(-jΩ))X",
    "Accumulation: Σ[k≤n] x[k] <-> X/(1-e^(-jΩ)) + πX(e^(j0))Σδ(Ω-2πl)",
    "Freq. differentiation: n x[n] <-> j dX/dΩ",
    "Convolution: x*h <-> X·H;  multiplication: x·y <-> (1/2π)∫[2π] X(θ)Y(Ω-θ)dθ",
    "Parseval: Σ|x[n]|² = (1/2π)∫[2π]|X|²dΩ",
    "Difference equation: y[n-k] → e^(-jkΩ)Y",
    "Inverse system: G(e^(jΩ)) = 1/H(e^(jΩ))",
  }},
  {title = "Laplace transform", lines = {
    "X(s) = ∫x(t)e^(-st)dt, ROC = strip in Re{s}",
    "δ(t) <-> 1 (all s);  u(t) <-> 1/s, Re{s}>0;  -u(-t) <-> 1/s, Re{s}<0",
    "e^(-at)u(t) <-> 1/(s+a), Re{s}>-a;  -e^(-at)u(-t) <-> 1/(s+a), Re{s}<-a",
    "t^(n-1)/(n-1)! e^(-at)u(t) <-> 1/(s+a)^n",
    "cos ω0t u(t) <-> s/(s²+ω0²);  sin ω0t u(t) <-> ω0/(s²+ω0²)",
    "e^(-at)cos ω0t u(t) <-> (s+a)/((s+a)²+ω0²)",
    "Shift x(t-t0) <-> e^(-st0)X;  dx/dt <-> sX;  ∫x <-> X/s;  x*h <-> XH",
    "Right-sided ⇒ ROC right of rightmost pole; left-sided ⇒ left of leftmost",
    "Causal & stable LTI ⇔ all poles in left half-plane (ROC includes jω axis)",
    "X(jω) = X(s) at s = jω when ROC contains the jω axis",
  }},
  {title = "Z transform", lines = {
    "X(z) = Σ x[n]z^(-n); ROC is a ring |z| between radii",
    "δ[n] <-> 1;  a^n u[n] <-> 1/(1-a z^-1), |z|>|a|;  -a^n u[-n-1] <-> same, |z|<|a|",
    "n a^n u[n] <-> a z^-1/(1-a z^-1)²",
    "x[n-n0] <-> z^(-n0)X;  x*h <-> XH;  n x[n] <-> -z dX/dz",
    "Causal & stable ⇔ poles inside unit circle;  DTFT = X(z) at z = e^(jΩ)",
  }},
  {title = "Fourier series", lines = {
    "CT: x(t) = Σ a_k e^(jkω0t),  a_k = (1/T)∫[T] x(t)e^(-jkω0t)dt,  ω0 = 2π/T",
    "DT: x[n] = Σ[k=<N>] a_k e^(jk(2π/N)n),  a_k = (1/N)Σ[n=<N>] x[n]e^(-jk(2π/N)n)",
    "a_0 = average value;  real x ⇒ a_-k = a_k*",
    "Square pulse |t|<T1, period T: a_k = sin(kω0T1)/(kπ), a_0 = 2T1/T",
    "Shift: x(t-t0) → a_k e^(-jkω0t0);  derivative → jkω0 a_k",
    "Parseval: (1/T)∫|x|² = Σ|a_k|²",
    "LTI: y = Σ a_k H(jkω0) e^(jkω0t)",
  }},
  {title = "Sampling", lines = {
    "x[n] = x(nTs),  ωs = 2π/Ts,  Ω = ωTs",
    "Nyquist: no aliasing if ωs > 2ωmax (fs > 2fmax)",
    "Aliasing: frequency ω appears at ω - kωs (fold into |ω| < ωs/2)",
    "Ideal reconstruction: low-pass cutoff ωs/2, gain Ts",
    "Spectrum of sampled signal: (1/Ts)Σ X(j(ω - kωs))",
  }},
}
-- ===================== NSPIRE USER INTERFACE =====================
DEFAULT_DEFS = {"x(t)=e^(-t)u(t)", "h(t)=e^(-2t)u(t)", "a=2", "", "", ""}
App = {
  screen = "home", sel = 1, top = 1,
  tool = nil, field = 1, cur = {}, vals = {}, defs = {},
  out = nil, rl = nil, rscroll = 0, rwidth = 0,
  plots = {}, pidx = 1, g = nil, busy = false, msg = nil, phase = false,
}
for i, d in ipairs(DEFAULT_DEFS) do App.defs[i] = d end

HOME_ITEMS = {}
function buildHome()
  HOME_ITEMS = {
    {label = "Search tools / examples / tables", kind = "search"},
    {label = "Example problems (" .. #EXAMPLES .. ")", kind = "examples"},
    {label = "Reference tables (pairs, properties)", kind = "refs"},
    {label = "Definitions  x(t), h(t), params", kind = "defs"},
  }
  for i, t in ipairs(TOOLS) do HOME_ITEMS[#HOME_ITEMS + 1] = {label = i .. ". " .. t.name, kind = "tool", idx = i} end
  HOME_ITEMS[#HOME_ITEMS + 1] = {label = "Syntax help & keys", kind = "help"}
end

HELP = {
  {k = "h", s = "Syntax"},
  {k = "t", s = "u(t) step, d(t) or delta(t) impulse, rect(t), tri(t), ramp(t), sgn(t), abs(t) or |t|"},
  {k = "t", s = "e^(-2t), sin, cos, pi, j (or i) imaginary, inf for ∞"},
  {k = "t", s = "DT: use n and brackets: u[n], d[n-2], (1/2)^n u[n], x[n-1]"},
  {k = "t", s = "Conditions: (abs(t)<1) or (0<=t<=2) are 1/0 windows"},
  {k = "t", s = "Integral: int(f, tau, a, b)   Sum: sum(f, k, a, b)"},
  {k = "t", s = "Derivative: der(x(t)) or x'(t) or d/dt x(t)"},
  {k = "t", s = "Convolution inline: conv(f(t), g(t))"},
  {k = "t", s = "Frequency: w = ω (jw), W = Ω (e^(-jW)), s, z"},
  {k = "h", s = "Definitions"},
  {k = "t", s = "x(t)=..., h(t)=..., x[n]=..., a=2 (params). Use x(2t+1), x[n-4] etc. anywhere."},
  {k = "h", s = "Finding things"},
  {k = "t", s = "Home: start typing letters to SEARCH tools, example problems and reference tables."},
  {k = "t", s = "Example problems: pick one, it fills the tool (and its definitions) and solves it. Esc to edit the inputs."},
  {k = "t", s = "Reference tables: transform pairs & properties, system tests, identities, series."},
  {k = "h", s = "Keys"},
  {k = "t", s = "Home: up/dn select, enter open.  Tool: tab/up/dn field, enter solve, esc back"},
  {k = "t", s = "Results: up/dn scroll, l/r page, G graph, esc edit"},
  {k = "t", s = "Graph: l/r trace, up/dn next plot, + - zoom, 4 6 8 2 pan, A auto, P phase, esc back"},
}

-- ---------- helpers ----------
function W() return platform.window:width() end
function Hh() return platform.window:height() end
function invalidate() platform.window:invalidate() end
function setFont(gc, st, sz) gc:setFont("sansserif", st or "r", sz or 9) end
function col(gc, r, g, b) gc:setColorRGB(r, g, b) end
function sw(gc, s) return gc:getStringWidth(s) end

function wrapText(gc, s, maxw)
  local out = {}
  local line = ""
  local chars = uchars(s)
  local i = 1
  local lastSpace = nil
  local cur = {}
  local function flush(upto)
    out[#out + 1] = table.concat(cur, "", 1, upto)
    local rest = {}
    for k = upto + 1, #cur do rest[#rest + 1] = cur[k] end
    cur = rest
  end
  for _, ch in ipairs(chars) do
    cur[#cur + 1] = ch
    if sw(gc, table.concat(cur)) > maxw then
      -- find break
      local br = nil
      for k = #cur - 1, 2, -1 do if cur[k] == " " then br = k; break end end
      if br and br > #cur / 3 then flush(br) else flush(#cur - 1) end
      -- indent continuation
      table.insert(cur, 1, "  ")
    end
  end
  if #cur > 0 then out[#out + 1] = table.concat(cur) end
  return out
end

function currentTool() return TOOLS[App.tool] end
function fieldsOf()
  if App.screen == "defs" then
    local f = {}
    for i = 1, #App.defs do f[i] = {"def " .. i, ""} end
    return f, App.defs
  end
  local t = currentTool()
  if not App.vals[App.tool] then
    App.vals[App.tool] = {}
    for i, f in ipairs(t.fields) do App.vals[App.tool][i] = f[2] end
  end
  return t.fields, App.vals[App.tool]
end
function curPos(i, s)
  local key = (App.screen == "defs" and "d" or App.tool) .. ":" .. i
  if App.cur[key] == nil then App.cur[key] = #uchars(s) end
  return App.cur[key], key
end

-- ---------- run ----------
function defsText() return table.concat(App.defs, "\n") end
function runTool()
  local t = currentTool()
  local _, vals = fieldsOf()
  local ctx = buildCtx(defsText())
  local out = newOut()
  for _, e in ipairs(ctx.errs) do out:e("Definition error: " .. e) end
  if App.exTitle then out:h("Example — " .. App.exTitle) end
  if App.exNote then out:t(App.exNote) end
  App.exTitle = nil; App.exNote = nil
  local ok, err = pcall(t.run, vals, ctx, out)
  if not ok then out:e("Error: " .. tostring(err)) end
  App.out = out
  App.rl = nil
  App.rscroll = 0
  App.plots = out.plots
  App.pidx = 1
  App.g = nil
  App.screen = "result"
end

-- ---------- painting: header/footer ----------
function header(gc, title, right)
  col(gc, 0, 0, 0); setFont(gc, "b", 10)
  gc:drawString(title, 4, 1, "top")
  if right then col(gc, 110, 110, 110); setFont(gc, "r", 7); gc:drawString(right, W() - sw(gc, right) - 4, 4, "top") end
  col(gc, 160, 160, 160); gc:drawLine(0, 17, W(), 17)
end
function footer(gc, s)
  col(gc, 255, 255, 255); gc:fillRect(0, Hh() - 13, W(), 13)
  col(gc, 160, 160, 160); gc:drawLine(0, Hh() - 13, W(), Hh() - 13)
  col(gc, 90, 90, 90); setFont(gc, "r", 7)
  gc:drawString(s, 3, Hh() - 12, "top")
end

function paintHome(gc)
  header(gc, "Signals & Systems", nil)
  setFont(gc, "r", 9)
  local lh = 15
  local rows = math.floor((Hh() - 32) / lh)
  if App.sel < App.top then App.top = App.sel end
  if App.sel >= App.top + rows then App.top = App.sel - rows + 1 end
  for i = App.top, math.min(#HOME_ITEMS, App.top + rows - 1) do
    local y = 19 + (i - App.top) * lh
    if i == App.sel then col(gc, 215, 215, 215); gc:fillRect(2, y, W() - 4, lh) end
    col(gc, 0, 0, 0)
    gc:drawString(HOME_ITEMS[i].label, 6, y, "top")
  end
  footer(gc, "up/dn select   enter open   menu: tools")
end

function paintFields(gc, title, help)
  header(gc, title, App.screen == "defs" and "enter=save" or "enter=solve")
  local fields, vals = fieldsOf()
  local y = 20
  if help then
    setFont(gc, "r", 7); col(gc, 90, 90, 90)
    for _, l in ipairs(wrapText(gc, help, W() - 8)) do gc:drawString(l, 4, y, "top"); y = y + 10 end
    y = y + 2
  end
  for i, f in ipairs(fields) do
    local active = (i == App.field)
    setFont(gc, "r", 9); col(gc, 40, 40, 40)
    if App.screen ~= "defs" then gc:drawString(f[1], 4, y, "top"); y = y + 12 end
    local bx, bw, bh = 4, W() - 8, 18
    col(gc, 255, 255, 255); gc:fillRect(bx, y, bw, bh)
    col(gc, active and 0 or 170, active and 0 or 170, active and 0 or 170); gc:drawRect(bx, y, bw, bh)
    setFont(gc, "r", 10); col(gc, 0, 0, 0)
    local s = vals[i] or ""
    local ch = uchars(s)
    local cp = curPos(i, s)
    if cp > #ch then cp = #ch end
    -- horizontal scroll so cursor visible
    local start = 1
    local maxw = bw - 8
    while start < cp and sw(gc, table.concat(ch, "", start, cp)) > maxw do start = start + 1 end
    local vis = {}
    for k = start, #ch do
      vis[#vis + 1] = ch[k]
      if sw(gc, table.concat(vis)) > maxw then vis[#vis] = nil; break end
    end
    gc:drawString(table.concat(vis), bx + 4, y + 1, "top")
    if active then
      local cx0 = bx + 4 + sw(gc, table.concat(ch, "", start, cp))
      if cp < start then cx0 = bx + 4 end
      col(gc, 0, 0, 0); gc:fillRect(cx0, y + 3, 1, bh - 6)
    end
    y = y + bh + 4
  end
  if App.busy then
    col(gc, 255, 255, 255); gc:fillRect(W() / 2 - 60, Hh() / 2 - 12, 120, 24)
    col(gc, 0, 0, 0); gc:drawRect(W() / 2 - 60, Hh() / 2 - 12, 120, 24)
    setFont(gc, "r", 10); gc:drawString("Working...", W() / 2 - 30, Hh() / 2 - 8, "top")
  end
  footer(gc, App.screen == "defs" and "up/dn/tab line   enter save   esc back" or "tab/up/dn field   enter solve   esc home")
end

KCOL = {h = {0, 0, 0}, t = {0, 0, 0}, r = {0, 0, 0}, e = {0, 0, 0}}
function paintResult(gc)
  local t = currentTool()
  local isView = (App.screen == "help" or App.screen == "view")
  header(gc, isView and (App.view and App.view.title or "Help") or t.name, (not isView and #App.plots > 0) and "G=graph" or nil)
  local maxw = W() - 10
  if not App.rl or App.rwidth ~= maxw then
    App.rl = {}
    App.rwidth = maxw
    local src = isView and (App.view and App.view.lines or HELP) or App.out.lines
    for _, l in ipairs(src) do
      local st = (l.k == "h" or l.k == "r") and "b" or "r"
      setFont(gc, st, 9)
      local txt = l.s
      if l.k == "e" then txt = "Error: " .. txt:gsub("^Error: ", "") end
      for _, w in ipairs(wrapText(gc, txt, maxw)) do App.rl[#App.rl + 1] = {k = l.k, s = w} end
    end
  end
  local lh = 13
  local rows = math.floor((Hh() - 32) / lh)
  local maxs = math.max(0, #App.rl - rows)
  if App.rscroll > maxs then App.rscroll = maxs end
  if App.rscroll < 0 then App.rscroll = 0 end
  for i = 1, rows do
    local l = App.rl[App.rscroll + i]
    if not l then break end
    local y = 19 + (i - 1) * lh
    if l.k == "h" and i > 1 then col(gc, 200, 200, 200); gc:drawLine(0, y, W(), y) end
    local c = KCOL[l.k] or KCOL.t
    col(gc, c[1], c[2], c[3])
    setFont(gc, (l.k == "h" or l.k == "r") and "b" or "r", 9)
    gc:drawString(l.s, 4, y, "top")
  end
  -- scrollbar
  if #App.rl > rows then
    local th = Hh() - 32
    local bh = math.max(10, th * rows / #App.rl)
    local by = 18 + (th - bh) * App.rscroll / maxs
    col(gc, 150, 150, 150); gc:fillRect(W() - 3, by, 2, bh)
  end
  footer(gc, isView and "up/dn scroll  l/r page  esc back" or "up/dn scroll  l/r page  G graph  esc edit")
end

-- ---------- graph ----------
function niceStep(range, n)
  local raw = range / n
  local mag = 10 ^ math.floor(math.log10(raw))
  local r = raw / mag
  local st
  if r < 1.5 then st = 1 elseif r < 3 then st = 2 elseif r < 7 then st = 5 else st = 10 end
  return st * mag
end
function plotRange(p)
  if p.kind == "freq" then
    if p.mode == "D" then return -pi, pi end
    return -10, 10
  end
  local dt = p.dt or (p.sig and p.sig.dt)
  local pts = {}
  local sig = p.sig
  if sig then
    for _, b in ipairs(sig.bp) do pts[#pts + 1] = b end
    for _, im in ipairs(sig.imps) do pts[#pts + 1] = im.t0 end
  end
  local lo, hi
  if #pts == 0 then
    if dt then lo, hi = -10, 10 else lo, hi = -2 * pi, 2 * pi end
    -- periodic-looking or exp: default
  else
    lo, hi = pts[1], pts[1]
    for _, x in ipairs(pts) do lo = math.min(lo, x); hi = math.max(hi, x) end
    local span = math.max(hi - lo, dt and 6 or 2)
    local padL, padR = span * 0.3, span * 0.3
    if sig and #sig.seg[1] > 0 then padL = math.max(padL, dt and 8 or 4) end
    if sig and #sig.seg[#sig.seg] > 0 then padR = math.max(padR, dt and 8 or 4) end
    lo, hi = lo - padL, hi + padR
  end
  if dt then lo, hi = math.floor(lo), math.ceil(hi) end
  return lo, hi
end
function evalPlot(p, x)
  local g = App.g
  if g and p.kind ~= "sig" then
    g.cache = g.cache or {}
    local key = string.format("%.9g", x)
    local v = g.cache[key]
    if v == nil then v = evalPlotRaw(p, x) or false; g.cache[key] = v end
    return v or nil
  end
  return evalPlotRaw(p, x)
end
evalPlotRaw = function(p, x)
  if p.kind == "sig" then return p.sig:eval(x) end
  if p.kind == "fn" then return p.f(x) end
  if p.kind == "freq" then
    local ok, v = pcall(p.f, x)
    if not ok or not v then return nil end
    return v
  end
end
function graphInit()
  local p = App.plots[App.pidx]
  local lo, hi = plotRange(p)
  App.g = {x0 = lo, x1 = hi, tr = nil, y0 = nil, y1 = nil}
end
function graphAutoY(p, g, N, xs)
  local ymin, ymax = INF, -INF
  local dt = (p.kind ~= "freq") and (p.dt or (p.sig and p.sig.dt))
  for _, x in ipairs(xs) do
    local v = evalPlot(p, x)
    if v then
      local y
      if p.kind == "freq" then y = App.phase and carg(v) or cabs(v) else y = v.re end
      if y == y and abs(y) < 1e9 then ymin = math.min(ymin, y); ymax = math.max(ymax, y)
        if p.kind ~= "freq" and abs(v.im) > 1e-9 then ymin = math.min(ymin, v.im); ymax = math.max(ymax, v.im) end
      end
    end
  end
  if p.kind == "sig" then
    for _, im in ipairs(p.sig.imps) do ymin = math.min(ymin, im.c.re, 0); ymax = math.max(ymax, im.c.re, 0) end
  end
  if p.imps then for _, im in ipairs(p.imps) do ymax = math.max(ymax, cabs(im.c)) end end
  if ymin == INF then ymin, ymax = -1, 1 end
  ymin = math.min(ymin, 0); ymax = math.max(ymax, 0)
  if ymax - ymin < 1e-9 then ymax = ymax + 1; ymin = ymin - 1 end
  local pad = (ymax - ymin) * 0.12
  return ymin - pad, ymax + pad
end

function paintGraph(gc)
  local p = App.plots[App.pidx]
  if not p then App.screen = "result"; return paintResult(gc) end
  if not App.g then graphInit() end
  local g = App.g
  local dt = (p.kind ~= "freq") and (p.dt or (p.sig and p.sig.dt)) or false
  local title = p.label .. "  (" .. App.pidx .. "/" .. #App.plots .. ")"
  if p.kind == "freq" then title = (App.phase and ("phase ∠" .. p.label:gsub("|", "")) or p.label) .. "  (" .. App.pidx .. "/" .. #App.plots .. ")" end
  header(gc, title, "up/dn plot  A auto")
  local L, T, R, B = 30, 20, W() - 6, Hh() - 26
  local xs = {}
  local N = dt and (g.x1 - g.x0) or (p.N or 260)
  if dt then for n = math.ceil(g.x0), math.floor(g.x1) do xs[#xs + 1] = n end
  else for i = 0, N do xs[#xs + 1] = g.x0 + (g.x1 - g.x0) * i / N end end
  if not g.y0 then g.y0, g.y1 = graphAutoY(p, g, N, xs) end
  local function X(x) return L + (x - g.x0) / (g.x1 - g.x0) * (R - L) end
  local function Y(y) return B - (y - g.y0) / (g.y1 - g.y0) * (B - T) end
  -- grid & axes
  setFont(gc, "r", 7)
  local xst = niceStep(g.x1 - g.x0, 6)
  if p.kind == "freq" and p.mode == "D" then xst = pi / 2 end
  local ystp = niceStep(g.y1 - g.y0, 5)
  col(gc, 235, 235, 235)
  local xv = math.ceil(g.x0 / xst) * xst
  while xv <= g.x1 do
    col(gc, 0, 0, 0); gc:drawLine(X(xv), B, X(xv), B - 3)
    col(gc, 0, 0, 0)
    local lbl = (p.kind == "freq" and p.mode == "D") and fmtR(xv) or decstr(cln(xv))
    gc:drawString(lbl, X(xv) - sw(gc, lbl) / 2, B + 1, "top")
    col(gc, 235, 235, 235)
    xv = xv + xst
  end
  local yv = math.ceil(g.y0 / ystp) * ystp
  while yv <= g.y1 do
    col(gc, 0, 0, 0); gc:drawLine(L, Y(yv), L + 3, Y(yv))
    col(gc, 0, 0, 0)
    local lbl = decstr(cln(yv))
    gc:drawString(lbl, L - sw(gc, lbl) - 2, Y(yv) - 5, "top")
    col(gc, 235, 235, 235)
    yv = yv + ystp
  end
  col(gc, 0, 0, 0)
  if g.y0 <= 0 and g.y1 >= 0 then gc:drawLine(L, Y(0), R, Y(0)) end
  if g.x0 <= 0 and g.x1 >= 0 then gc:drawLine(X(0), T, X(0), B) end
  gc:drawRect(L, T, R - L, B - T)
  local function clipY(y) return math.max(T, math.min(B, y)) end
  -- curve
  local function getY(v)
    if not v then return nil end
    if p.kind == "freq" then return App.phase and carg(v) or cabs(v) end
    return v.re
  end
  if dt then
    for _, n in ipairs(xs) do
      local v = evalPlot(p, n)
      local y = getY(v)
      if y and y == y then
        local px, py = X(n), clipY(Y(y))
        col(gc, 0, 0, 160)
        gc:drawLine(px, clipY(Y(0)), px, py)
        gc:fillArc(px - 2, py - 2, 5, 5, 0, 360)
        if v and abs(v.im) > 1e-9 then col(gc, 160, 0, 0); gc:fillArc(px - 2, clipY(Y(v.im)) - 2, 5, 5, 0, 360) end
      end
    end
  else
    local prevx, prevy, prevseg, previm = nil, nil, nil, nil
    for _, x in ipairs(xs) do
      local v = evalPlot(p, x)
      local y = getY(v)
      local seg = (p.kind == "sig") and p.sig:segIndex(x) or nil
      if y and y == y and abs(y) < 1e12 then
        local px, py = X(x), Y(y)
        if prevx then
          if seg ~= prevseg then col(gc, 0, 0, 160) else col(gc, 0, 0, 160) end
          gc:drawLine(prevx, clipY(prevy), px, clipY(py))
          if v and p.kind ~= "freq" and previm and abs(v.im) > 1e-9 then
            col(gc, 160, 0, 0); gc:drawLine(prevx, clipY(previm), px, clipY(Y(v.im)))
          end
        end
        prevx, prevy, prevseg = px, py, seg
        previm = (v and p.kind ~= "freq") and Y(v.im) or nil
      else prevx = nil end
    end
  end
  -- impulses
  local function arrow(x, h, label)
    if x < g.x0 or x > g.x1 then return end
    local px = X(x)
    local y0p, y1p = clipY(Y(0)), clipY(Y(h))
    col(gc, 0, 0, 0)
    gc:drawLine(px, y0p, px, y1p)
    gc:drawLine(px + 1, y0p, px + 1, y1p)
    local d = (h >= 0) and 1 or -1
    gc:fillPolygon({px - 3, y1p + 5 * d, px + 4, y1p + 5 * d, px + 0.5, y1p, px - 3, y1p + 5 * d})
    setFont(gc, "r", 7)
    gc:drawString("(" .. label .. ")", px + 3, y1p - (d > 0 and 2 or 8), "top")
  end
  if p.kind == "sig" then for _, im in ipairs(p.sig.imps) do arrow(im.t0, im.c.re, fmtC(im.c) .. (im.m > 0 and "δ'" or "")) end end
  if p.imps and not App.phase then
    for _, im in ipairs(p.imps) do
      local per = (p.mode == "D") and 2 * pi or nil
      if per then
        for l = -2, 2 do arrow(im.w0 + l * per, cabs(im.c), fmtC(im.c)) end
      else arrow(im.w0, cabs(im.c), fmtC(im.c)) end
    end
  end
  -- trace
  if g.tr then
    local v = evalPlot(p, g.tr)
    col(gc, 120, 120, 120)
    gc:drawLine(X(g.tr), T, X(g.tr), B)
    local s = (dt and "n=" or ((p.kind == "freq") and (p.mode == "D" and "Ω=" or "ω=") or "t=")) .. decstr(cln(g.tr))
    if v then
      if p.kind == "freq" then s = s .. "  |X|=" .. decstr(cabs(v)) .. "  ∠=" .. decstr(carg(v))
      else s = s .. "  x=" .. fmtC(v) end
    end
    footer(gc, s)
  else
    footer(gc, "l/r trace  +/- zoom  4 6 8 2 pan" .. ((p.kind == "freq") and "  P phase" or "") .. "  esc")
  end
end

function paintList(gc)
  local L = App.list
  header(gc, L.title, #L.items .. " items")
  setFont(gc, "r", 9)
  local lh = 14
  local rows = math.floor((Hh() - 32) / lh)
  if L.sel < L.top then L.top = L.sel end
  if L.sel >= L.top + rows then L.top = L.sel - rows + 1 end
  if L.top > 1 and L.items[L.top - 1] and L.items[L.top - 1].kind == "hdr" and L.sel - L.top + 2 <= rows then L.top = L.top - 1 end
  for i = L.top, math.min(#L.items, L.top + rows - 1) do
    local it = L.items[i]
    local y = 19 + (i - L.top) * lh
    if it.kind == "hdr" then
      col(gc, 0, 0, 0); setFont(gc, "b", 9)
    else
      if i == L.sel then col(gc, 215, 215, 215); gc:fillRect(2, y, W() - 4, lh) end
      col(gc, 0, 0, 0); setFont(gc, "r", 9)
    end
    gc:drawString(it.label, 6, y, "top")
  end
  footer(gc, "up/dn select   enter open   esc back")
end
function paintSearch(gc)
  header(gc, "Search", "type words")
  local bx, by, bw, bh = 4, 20, W() - 8, 18
  col(gc, 255, 255, 255); gc:fillRect(bx, by, bw, bh)
  col(gc, 0, 0, 0); gc:drawRect(bx, by, bw, bh)
  setFont(gc, "r", 10); col(gc, 0, 0, 0)
  local q = App.query or ""
  gc:drawString(q, bx + 4, by + 1, "top")
  col(gc, 0, 0, 0); gc:fillRect(bx + 4 + sw(gc, q), by + 3, 1, bh - 6)
  local res = App.sres or {}
  setFont(gc, "r", 9)
  local lh = 14
  local rows = math.floor((Hh() - 56) / lh)
  App.ssel = math.max(1, math.min(App.ssel or 1, #res))
  App.stop = App.stop or 1
  if App.ssel < App.stop then App.stop = App.ssel end
  if App.ssel >= App.stop + rows then App.stop = App.ssel - rows + 1 end
  if #res == 0 then
    col(gc, 110, 110, 110)
    gc:drawString((q == "") and "e.g. convolution, causal, dtft, laplace, aliasing" or "no matches", 6, 42, "top")
  end
  for i = App.stop, math.min(#res, App.stop + rows - 1) do
    local y = 42 + (i - App.stop) * lh
    if i == App.ssel then col(gc, 215, 215, 215); gc:fillRect(2, y, W() - 4, lh) end
    col(gc, 0, 0, 0)
    gc:drawString(res[i].label, 6, y, "top")
  end
  footer(gc, #res .. " results   up/dn select   enter open   esc back")
end

-- ---------- event handlers ----------
function on.paint(gc)
  if App.screen == "home" then paintHome(gc)
  elseif App.screen == "tool" then local t = currentTool(); paintFields(gc, t.name, t.help)
  elseif App.screen == "defs" then paintFields(gc, "Definitions", "Lines like x(t)=..., h[n]=..., a=2. Used by every tool.")
  elseif App.screen == "result" or App.screen == "help" or App.screen == "view" then paintResult(gc)
  elseif App.screen == "list" then paintList(gc)
  elseif App.screen == "search" then paintSearch(gc)
  elseif App.screen == "graph" then paintGraph(gc) end
end

function loadExample(ex)
  local idx
  for i, t in ipairs(TOOLS) do if t.name == ex.tool then idx = i end end
  if not idx then return end
  if ex.defs then
    local d = {}
    for part in (ex.defs .. ";"):gmatch("([^;]*);") do if part ~= "" then d[#d + 1] = part end end
    while #d < 6 do d[#d + 1] = "" end
    App.defs = d
  end
  local vals = {}
  for i = 1, #TOOLS[idx].fields do vals[i] = ex.vals[i] or "" end
  App.vals[idx] = vals
  for k in pairs(App.cur) do if tostring(k):match("^" .. idx .. ":") then App.cur[k] = nil end end
  App.tool = idx; App.screen = "tool"; App.field = 1
  App.exNote = ex.defs and ("Example definitions loaded: " .. ex.defs:gsub(";", ";  ")) or nil
  App.exTitle = ex.cat .. ": " .. ex.title
  App.busy = true
  invalidate()
  timer.start(0.05)
end
function openList(title, items, back)
  App.list = {title = title, items = items, sel = 1, top = 1, back = back or "home"}
  for i, it in ipairs(items) do if it.kind ~= "hdr" then App.list.sel = i; break end end
  App.screen = "list"
  invalidate()
end
function exampleItems()
  local items, last = {}, nil
  for _, ex in ipairs(EXAMPLES) do
    if ex.cat ~= last then items[#items + 1] = {label = ex.cat, kind = "hdr"}; last = ex.cat end
    items[#items + 1] = {label = "  " .. ex.title, kind = "example", ex = ex}
  end
  return items
end
function refItems()
  local items = {}
  for _, r in ipairs(REFS) do items[#items + 1] = {label = r.title, kind = "ref", ref = r} end
  return items
end
function openView(title, lines, back)
  App.view = {title = title, lines = lines, back = back}
  App.screen = "view"; App.rl = nil; App.rscroll = 0
  invalidate()
end
function refLines(r)
  local L = {{k = "h", s = r.title}}
  for _, s in ipairs(r.lines) do L[#L + 1] = {k = "t", s = s} end
  return L
end
function lower(s) return (s:gsub("%u", string.lower)) end
function searchItems(q)
  local words = {}
  for w in lower(preprocess(q)):gmatch("%S+") do words[#words + 1] = w end
  local res = {}
  local function match(text)
    local t = lower(preprocess(text))
    for _, w in ipairs(words) do if not t:find(w, 1, true) then return false end end
    return true
  end
  if #words == 0 then return res end
  for i, t in ipairs(TOOLS) do
    if match(t.name .. " " .. (t.help or "")) then res[#res + 1] = {label = "[tool] " .. t.name, kind = "tool", idx = i} end
  end
  for _, ex in ipairs(EXAMPLES) do
    if match(ex.cat .. " " .. ex.title .. " " .. ex.tool .. " " .. table.concat(ex.vals, " ")) then
      res[#res + 1] = {label = "[ex] " .. ex.cat .. ": " .. ex.title, kind = "example", ex = ex}
    end
  end
  for _, r in ipairs(REFS) do
    local hit = match(r.title)
    for _, l in ipairs(r.lines) do if match(l) then hit = true end end
    if hit then res[#res + 1] = {label = "[table] " .. r.title, kind = "ref", ref = r} end
  end
  return res
end
openItem = function(it)
  if it.kind == "tool" then App.tool = it.idx; App.screen = "tool"; App.field = 1; App.exNote = nil; App.exTitle = nil
  elseif it.kind == "defs" then App.screen = "defs"; App.field = 1
  elseif it.kind == "help" then openView("Syntax help", HELP, "home"); return
  elseif it.kind == "examples" then openList("Example problems", exampleItems(), "home"); return
  elseif it.kind == "refs" then openList("Reference tables", refItems(), "home"); return
  elseif it.kind == "search" then App.screen = "search"; App.query = App.query or ""; App.sres = searchItems(App.query); App.ssel = 1; App.stop = 1
  elseif it.kind == "example" then loadExample(it.ex); return
  elseif it.kind == "ref" then openView(it.ref.title, refLines(it.ref), App.screen); return
  end
  invalidate()
end

function editField(fn)
  local fields, vals = fieldsOf()
  local i = App.field
  local s = vals[i] or ""
  local ch = uchars(s)
  local cp, key = curPos(i, s)
  if cp > #ch then cp = #ch end
  local nch, ncp = fn(ch, cp)
  vals[i] = table.concat(nch)
  App.cur[key] = ncp
  invalidate()
end

function searchUpdate()
  App.sres = searchItems(App.query or "")
  App.ssel, App.stop = 1, 1
  invalidate()
end
function on.charIn(c)
  if App.screen == "search" then
    App.query = (App.query or "") .. c
    searchUpdate()
    return
  end
  if App.screen == "home" and c:match("^%a") then
    App.screen = "search"; App.query = c; searchUpdate(); return
  end
  if App.screen == "tool" or App.screen == "defs" then
    editField(function(ch, cp)
      local ins = uchars(c)
      for k = #ins, 1, -1 do table.insert(ch, cp + 1, ins[k]) end
      return ch, cp + #ins
    end)
  elseif App.screen == "result" then
    if (c == "g" or c == "G") and #App.plots > 0 then App.screen = "graph"; App.pidx = 1; App.g = nil; invalidate() end
    if c == "e" or c == "E" then App.screen = "tool"; invalidate() end
  elseif App.screen == "graph" then
    local g = App.g
    if not g then return end
    local span = g.x1 - g.x0
    local c0 = g.tr or (g.x0 + g.x1) / 2
    if c == "+" or c == "*" then g.x0 = c0 - span / 4; g.x1 = c0 + span / 4; g.y0 = nil
    elseif c == "-" or c == "/" then g.x0 = c0 - span; g.x1 = c0 + span; g.y0 = nil
    elseif c == "4" then g.x0 = g.x0 - span / 5; g.x1 = g.x1 - span / 5
    elseif c == "6" then g.x0 = g.x0 + span / 5; g.x1 = g.x1 + span / 5
    elseif c == "8" then local d = (g.y1 - g.y0) / 5; g.y0 = (g.y0 or 0) + d; g.y1 = g.y1 + d
    elseif c == "2" then local d = (g.y1 - g.y0) / 5; g.y0 = (g.y0 or 0) - d; g.y1 = g.y1 - d
    elseif c == "a" or c == "A" then graphInit()
    elseif c == "p" or c == "P" then App.phase = not App.phase; g.y0 = nil end
    if App.plots[App.pidx] and (App.plots[App.pidx].dt or (App.plots[App.pidx].sig and App.plots[App.pidx].sig.dt)) then
      g.x0, g.x1 = math.floor(g.x0), math.ceil(g.x1)
    end
    invalidate()
  elseif App.screen == "home" then
    local n = tonumber(c)
    if n then for _, it in ipairs(HOME_ITEMS) do if it.kind == "tool" and it.idx == n then openItem(it); return end end end
  end
end
function on.backspaceKey()
  if App.screen == "search" then
    local ch = uchars(App.query or "")
    ch[#ch] = nil
    App.query = table.concat(ch)
    searchUpdate()
    return
  end
  if App.screen == "tool" or App.screen == "defs" then
    editField(function(ch, cp)
      if cp > 0 then table.remove(ch, cp); cp = cp - 1 end
      return ch, cp
    end)
  elseif App.screen ~= "home" then on.escapeKey() end
end
function on.mouseUp() end
function on.deleteKey() on.backspaceKey() end
function on.clearKey()
  if App.screen == "search" then App.query = ""; searchUpdate(); return end
  if App.screen == "tool" or App.screen == "defs" then editField(function() return {}, 0 end) end
end
function moveField(d)
  local fields = fieldsOf()
  App.field = ((App.field - 1 + d) % #fields) + 1
  invalidate()
end
function on.tabKey() if App.screen == "tool" or App.screen == "defs" then moveField(1) elseif App.screen == "graph" then on.arrowKey("down") end end
function on.backtabKey() if App.screen == "tool" or App.screen == "defs" then moveField(-1) end end
function on.arrowKey(k)
  local s = App.screen
  if s == "home" then
    if k == "up" then App.sel = (App.sel - 2) % #HOME_ITEMS + 1 elseif k == "down" then App.sel = App.sel % #HOME_ITEMS + 1 end
  elseif s == "tool" or s == "defs" then
    if k == "up" then moveField(-1) elseif k == "down" then moveField(1)
    elseif k == "left" then editField(function(ch, cp) return ch, math.max(0, cp - 1) end)
    elseif k == "right" then editField(function(ch, cp) return ch, math.min(#ch, cp + 1) end) end
  elseif s == "list" then
    local L = App.list
    local d = (k == "up" or k == "left") and -1 or 1
    local steps = (k == "left" or k == "right") and 8 or 1
    for _ = 1, steps do
      local i = L.sel
      repeat i = i + d until i < 1 or i > #L.items or L.items[i].kind ~= "hdr"
      if i >= 1 and i <= #L.items then L.sel = i end
    end
  elseif s == "search" then
    if k == "up" then App.ssel = math.max(1, (App.ssel or 1) - 1) elseif k == "down" then App.ssel = (App.ssel or 1) + 1 end
  elseif s == "result" or s == "help" or s == "view" then
    local rows = math.floor((Hh() - 32) / 13)
    if k == "up" then App.rscroll = App.rscroll - 1 elseif k == "down" then App.rscroll = App.rscroll + 1
    elseif k == "left" then App.rscroll = App.rscroll - rows + 1 elseif k == "right" then App.rscroll = App.rscroll + rows - 1 end
  elseif s == "graph" then
    local g = App.g
    if k == "up" or k == "down" then
      App.pidx = (App.pidx - 1 + (k == "down" and 1 or -1)) % #App.plots + 1
      App.g = nil
    else
      local p = App.plots[App.pidx]
      local dt = p and (p.kind ~= "freq") and (p.dt or (p.sig and p.sig.dt))
      local step = dt and 1 or (g.x1 - g.x0) / 100
      if not g.tr then g.tr = dt and 0 or (g.x0 + g.x1) / 2 end
      g.tr = g.tr + ((k == "right") and step or -step)
      if g.tr < g.x0 then g.tr = g.x0 elseif g.tr > g.x1 then g.tr = g.x1 end
    end
  end
  invalidate()
end
function on.enterKey()
  local s = App.screen
  if s == "home" then openItem(HOME_ITEMS[App.sel])
  elseif s == "list" then local it = App.list.items[App.list.sel]; if it and it.kind ~= "hdr" then openItem(it) end
  elseif s == "search" then local it = App.sres and App.sres[App.ssel]; if it then openItem(it) end
  elseif s == "tool" then
    App.busy = true
    invalidate()
    timer.start(0.05)
  elseif s == "defs" then App.screen = "home"; invalidate()
  elseif s == "result" then
    if #App.plots > 0 then App.screen = "graph"; App.pidx = 1; App.g = nil end
    invalidate()
  end
end
function on.timer()
  timer.stop()
  if App.busy then
    App.busy = false
    runTool()
    invalidate()
  end
end
function on.escapeKey()
  local s = App.screen
  if s == "list" then App.screen = App.list.back or "home"
  elseif s == "search" then App.screen = "home"
  elseif s == "view" or s == "help" then
    local b = App.view and App.view.back or "home"
    App.screen = (b == "list" and App.list) and "list" or ((b == "search") and "search" or "home")
  elseif s == "tool" or s == "defs" then App.screen = "home"
  elseif s == "result" then App.screen = "tool"
  elseif s == "graph" then App.screen = "result" end
  invalidate()
end
function on.mouseDown(x, y)
  if App.screen == "list" then
    local i = App.list.top + math.floor((y - 19) / 14)
    local it = App.list.items[i]
    if it and it.kind ~= "hdr" then App.list.sel = i; openItem(it) end
    return
  end
  if App.screen == "search" then
    local i = (App.stop or 1) + math.floor((y - 42) / 14)
    if App.sres and App.sres[i] then App.ssel = i; openItem(App.sres[i]) end
    return
  end
  if App.screen == "home" then
    local i = App.top + math.floor((y - 19) / 15)
    if HOME_ITEMS[i] then App.sel = i; openItem(HOME_ITEMS[i]) end
  elseif App.screen == "tool" or App.screen == "defs" then
    -- pick field by y
    local fields = fieldsOf()
    local best, bd = App.field, 1e9
    -- approximate positions
    for i = 1, #fields do
      local fy = 40 + (i - 1) * (App.screen == "defs" and 22 or 34)
      if abs(fy - y) < bd then bd = abs(fy - y); best = i end
    end
    App.field = best; invalidate()
  end
end
function on.save() return {defs = App.defs, vals = App.vals, v = 1} end
function on.restore(st)
  if type(st) == "table" then
    if st.defs then App.defs = st.defs end
    if st.vals then App.vals = st.vals end
  end
end
function on.resize() App.rl = nil; invalidate() end
function on.construction() buildHome() end
function on.paste()
  if (App.screen == "tool" or App.screen == "defs") and clipboard then
    local s = clipboard.getText()
    if s then on.charIn(s) end
  end
end
function on.copy()
  if clipboard and App.screen == "result" and App.out then
    for _, l in ipairs(App.out.lines) do if l.k == "r" then clipboard.addText(l.s); return end end
  end
end

-- menu
function registerMenu()
  if not toolpalette then return end
  local menu = {}
  local groups = {
    {"Signals", {1, 2, 3, 4, 21}},
    {"Systems", {5, 6, 15, 20}},
    {"Transforms", {7, 8, 9, 10, 11, 12, 13, 14, 16}},
    {"Other", {17, 18, 19}},
  }
  for _, gpair in ipairs(groups) do
    local sub = {gpair[1]}
    for _, i in ipairs(gpair[2]) do
      if TOOLS[i] then
        sub[#sub + 1] = {TOOLS[i].name, function() App.tool = i; App.screen = "tool"; App.field = 1; invalidate() end}
      end
    end
    menu[#menu + 1] = sub
  end
  menu[#menu + 1] = {"App", {"Home", function() App.screen = "home"; invalidate() end},
    {"Search", function() App.screen = "search"; App.query = ""; App.sres = {}; invalidate() end},
    {"Examples", function() openList("Example problems", exampleItems(), "home") end},
    {"Reference tables", function() openList("Reference tables", refItems(), "home") end},
    {"Definitions", function() App.screen = "defs"; App.field = 1; invalidate() end},
    {"Graph results", function() if #App.plots > 0 then App.screen = "graph"; App.g = nil; invalidate() end end},
    {"Syntax help", function() openView("Syntax help", HELP, "home") end},
    {"Fractions on/off", function() fracmode = not fracmode; invalidate() end}}
  toolpalette.register(menu)
end
buildHome()
registerMenu()
