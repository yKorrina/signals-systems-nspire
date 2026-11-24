-- ===================== CORE: complex numbers & formatting =====================
local abs, floor, sqrt, exp, log, sin, cos, atan2, pi = math.abs, math.floor, math.sqrt, math.exp, math.log, math.sin, math.cos, math.atan2, math.pi
local INF = math.huge
local EPS = 1e-9

local C = {}
C.__index = C
local function cx(re, im) return setmetatable({re = re, im = im or 0}, C) end
local function tocx(a) if type(a) == "number" then return cx(a, 0) end return a end
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
local function cabs(a) a = tocx(a); return sqrt(a.re * a.re + a.im * a.im) end
local function carg(a) a = tocx(a); return atan2(a.im, a.re) end
local function cexp(a) a = tocx(a); local m = exp(a.re); return cx(m * cos(a.im), m * sin(a.im)) end
local function clog(a) a = tocx(a); return cx(log(cabs(a)), carg(a)) end
local function cconj(a) a = tocx(a); return cx(a.re, -a.im) end
local function cpow(a, b)
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
local function csqrt(a) a = tocx(a); local r = cabs(a); local t = carg(a) / 2; local m = sqrt(r); return cx(m * cos(t), m * sin(t)) end
local function csin(a) a = tocx(a); return (cexp(cx(-a.im, a.re)) - cexp(cx(a.im, -a.re))) / cx(0, 2) end
local function ccos(a) a = tocx(a); return (cexp(cx(-a.im, a.re)) + cexp(cx(a.im, -a.re))) / 2 end
local function isz(a, tol) tol = tol or EPS; a = tocx(a); return abs(a.re) <= tol and abs(a.im) <= tol end
local function near(a, b, tol)
  a, b = tocx(a), tocx(b); tol = tol or 1e-7
  local d = cabs(a - b); local m = math.max(1, cabs(a), cabs(b))
  return d <= tol * m
end
local function isreal(a, tol) a = tocx(a); return abs(a.im) <= (tol or 1e-9) * math.max(1, abs(a.re)) end
local J = cx(0, 1)
local ZERO, ONE = cx(0, 0), cx(1, 0)

-- clean tiny float noise
local function cln(x)
  if abs(x) < 1e-11 then return 0 end
  local r = floor(x + 0.5)
  if abs(x - r) < 1e-9 * math.max(1, abs(x)) then return r end
  return x
end
local function ccln(z) z = tocx(z); return cx(cln(z.re), cln(z.im)) end

-- ---------- number formatting ----------
local function ratapprox(x, maxq, tol)
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

local function trimdec(s)
  if s:find("e") then return s end
  if s:find("%.") then s = s:gsub("0+$", ""):gsub("%.$", "") end
  return s
end
local function decstr(x)
  local s
  if abs(x) >= 1e5 or (abs(x) < 1e-3 and x ~= 0) then s = string.format("%.4e", x)
  else s = trimdec(string.format("%.5g", x)) end
  return s
end

local fracmode = true
local function fracstr(p, q)
  if q == 1 then return tostring(p) end
  return p .. "/" .. q
end
local SQ = {2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21}
-- format real number; returns string, and "kind"
local function fmtR(x)
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
local function needpar(s)
  return s:find("[%+/]") ~= nil or s:find(".%-") ~= nil
end

local function fmtImag(y)
  -- y real, returns string like "j", "-j", "j2", "-j(1/2)"
  local s = fmtR(abs(y))
  local sg = y < 0 and "-" or ""
  if s == "1" then return sg .. "j" end
  if s:find("^[%d%.]+$") or s:find("^[%d]*π$") or s:find("^√%d+$") then return sg .. "j" .. s end
  return sg .. "j(" .. s .. ")"
end

local polarmode = false
local function fmtC(z)
  z = tocx(z)
  local re, im = cln(z.re), cln(z.im)
  if abs(im) <= 1e-10 * math.max(1, abs(re)) then return fmtR(re) end
  if abs(re) <= 1e-10 * math.max(1, abs(im)) then return fmtImag(im) end
  local s = fmtR(re)
  if im < 0 then return s .. " - " .. fmtImag(-im) end
  return s .. " + " .. fmtImag(im)
end
local function fmtPolar(z)
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
local function cmul(c, expr)
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
local function topSum(s)
  local d = 0
  for i = 1, #s do
    local c = s:sub(i, i)
    if c == "(" or c == "[" then d = d + 1 elseif c == ")" or c == "]" then d = d - 1
    elseif d == 0 and i > 1 and (c == "+" or c == "-") and s:sub(i - 1, i - 1) == " " then return true end
  end
  return false
end
-- join term strings into sum with proper signs
local function joinsum(terms)
  if #terms == 0 then return "0" end
  local out = terms[1]
  for i = 2, #terms do
    local t = terms[i]
    if t:sub(1, 1) == "-" then out = out .. " - " .. t:sub(2) else out = out .. " + " .. t end
  end
  return out
end

-- ---------- misc helpers ----------
local function gcdi(a, b) a, b = abs(a), abs(b); while b ~= 0 do a, b = b, a % b end; return a end
local function lcmi(a, b) return a / gcdi(a, b) * b end
local function fact(n) local r = 1; for i = 2, n do r = r * i end; return r end
local function binom(n, k) if k < 0 or k > n then return 0 end; return fact(n) / (fact(k) * fact(n - k)) end
local function copy(t) local r = {}; for k, v in pairs(t) do r[k] = v end; return r end
local function isint(x) return abs(x - floor(x + 0.5)) < 1e-9 end
local function round(x) return floor(x + 0.5) end

-- UTF-8 helpers
local function uchars(s)
  local t = {}
  for ch in s:gmatch("[%z\1-\127\194-\244][\128-\191]*") do t[#t + 1] = ch end
  return t
end
local function ulen(s) local _, n = s:gsub("[^\128-\191]", ""); return n end
