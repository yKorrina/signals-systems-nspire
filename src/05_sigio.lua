-- ===================== AST -> Sig conversion =====================
local function baseSig(name, dt)
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

local Conv = {}
Conv.__index = Conv
-- ctx: {var=, dt=, env=, defs = {name={param=,body=,dt=}}, cache={}}
local function mkConv(var, dt, env, defs)
  return setmetatable({var = var, dt = dt, env = env, defs = defs or {}, cache = {}}, Conv)
end

function Conv:lin(node, var)
  local a, b = linearOf(node, var or self.var, self.env)
  return a, b
end
local function cre(z, what) z = tocx(z); if abs(z.im) > 1e-9 then error(what .. " must be real", 0) end; return z.re end

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
local function expoArg(sv, var)
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
local function expStr(s, var)
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
local function powStr(var, k)
  if k == 0 then return "" elseif k == 1 then return var end
  return var .. "^" .. k
end
local function baseStr(r, var)
  r = tocx(r)
  if near(r, ONE, 1e-12) then return "" end
  local s = fmtC(r)
  if isreal(r) and r.re > 0 and not s:find("[/%-%s√π·]") then return s .. "^" .. var end
  return "(" .. s .. ")^" .. var
end
local function joinf(...)
  local t = {}
  for _, s in ipairs({...}) do if s ~= "" then t[#t + 1] = s end end
  local out = ""
  for i, s in ipairs(t) do
    if i > 1 and out:find("[%a]$") and s:find("^[%a]") then out = out .. "·" .. s else out = out .. s end
  end
  return out
end

-- format monomial list as a string (combining conjugate pairs into cos/sin)
local function segStr(S, var, dt)
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

local function impStr(im, var)
  local d = "δ" .. string.rep("'", im.m)
  local a = (abs(im.t0) < 1e-12) and var or (var .. (im.t0 > 0 and " - " or " + ") .. fmtR(abs(im.t0)))
  return cmul(im.c, d .. "(" .. a .. ")")
end
local function shiftStr(var, a)
  -- (var - a)
  if abs(a) < 1e-12 then return var end
  if a > 0 then return var .. " - " .. fmtR(a) end
  return var .. " + " .. fmtR(-a)
end
local function uStr(var, a, dt, rev)
  local o, c = dt and "[" or "(", dt and "]" or ")"
  if rev then
    -- u(a - var)
    if abs(a) < 1e-12 then return "u" .. o .. "-" .. var .. c end
    return "u" .. o .. fmtR(a) .. " - " .. var .. c
  end
  return "u" .. o .. shiftStr(var, a) .. c
end
local function windowStr(lo, hi, var, dt)
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
local function attach(expr, w)
  if w == "" then return expr end
  if expr == "1" then return w end
  if expr == "-1" then return "-" .. w end
  local single = not topSum(expr)
  if single and expr:find("^%-?[%d]+/[%d]+$") then return "(" .. expr .. ")" .. w end
  if single then return joinf(expr, w) end
  return "(" .. expr .. ")" .. w
end

-- interval description
local function rangeStr(lo, hi, var, dt)
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
