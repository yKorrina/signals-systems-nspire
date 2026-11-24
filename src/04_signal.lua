-- ===================== SIGNAL ENGINE =====================
-- Sig: piecewise exponential-polynomial signal
--  CT: segment i = interval (bp[i-1], bp[i]); monomial {c,k,s} = c t^k e^{s t}
--  DT: segment i = integers bp[i-1] <= n < bp[i]; monomial {c,k,s} = c n^k s^n  (s is the base r)
--  imps (CT only): {c, t0, m} = c * delta^(m)(t - t0)
local Sig = {}
Sig.__index = Sig
local TOL = 1e-9

local function M(c, k, s) return {c = tocx(c), k = k or 0, s = tocx(s)} end
local function newSig(dt) return setmetatable({dt = dt and true or false, bp = {}, seg = {{}}, imps = {}}, Sig) end

local function monoEval(m, t, dt)
  local v = m.c
  if m.k > 0 then v = v * cx(t ^ m.k) end
  if dt then
    if not near(m.s, ONE, 1e-12) then v = v * cpow(m.s, t) end
  else
    if not isz(m.s, 0) then v = v * cexp(m.s * t) end
  end
  return v
end
local function segEval(S, t, dt)
  local v = ZERO
  for _, m in ipairs(S) do v = v + monoEval(m, t, dt) end
  return v
end
local function sameS(a, b, dt)
  local idA = dt and ONE or ZERO
  return near(a, b, 1e-9)
end
-- merge like monomials, drop zeros
local function mergeS(S, dt)
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
local function eqS(A, B)
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
local function fromPieces(dt, pieces, imps)
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
local function sconst(c, dt) local r = newSig(dt); r.seg = {{M(c, 0, dt and 1 or 0)}}; return r:normalize() end
local function svar(dt) local r = newSig(dt); r.seg = {{M(1, 1, dt and 1 or 0)}}; return r end
local function swindow(lo, hi, dt, S) return fromPieces(dt, {{lo = lo, hi = hi, S = S or {M(1, 0, dt and 1 or 0)}}}) end
local function sexp(s, dt) local r = newSig(dt); r.seg = {{M(1, 0, s)}}; return r end

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

local function sadd(A, B)
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
local function ssub(A, B) return sadd(A, B:scale(-1)) end

local function mulS(A, B, dt)
  local out = {}
  for _, a in ipairs(A) do
    for _, b in ipairs(B) do
      out[#out + 1] = M(a.c * b.c, a.k + b.k, dt and (a.s * b.s) or (a.s + b.s))
    end
  end
  return out
end
-- derivative of monomial list (CT)
local function derS(S)
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

local function smul(A, B)
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
local function antiS(S)
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
local function solveLin(A, b)
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
local function antidiffQ(k, r)
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
local function antiDiffS(S)
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
local function limitS(S, dir, dt)
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
local function antiTerms(p, D, dt)
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
local function antiAt(p, D, lim, dt)
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
local function convPieces(pa, ma, pb, mb, dt, out)
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

local function sconv(A, B)
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
