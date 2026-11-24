-- ===================== POLYNOMIALS (complex coeffs, index 0..deg) =====================
local Poly = {}
local function pnew(t) local p = {}; for i = 0, #t - 1 do p[i] = tocx(t[i + 1]) end; p.n = #t - 1; return p end
local function pdeg(p)
  local d = p.n
  while d > 0 and cabs(p[d]) < 1e-12 do d = d - 1 end
  return d
end
local function ptrim(p) local d = pdeg(p); local r = {n = d}; for i = 0, d do r[i] = p[i] end; return r end
local function pconst(c) return {[0] = tocx(c), n = 0} end
local function padd(a, b)
  local n = math.max(a.n, b.n); local r = {n = n}
  for i = 0, n do r[i] = (a[i] or ZERO) + (b[i] or ZERO) end
  return ptrim(r)
end
local function pscale(a, c) local r = {n = a.n}; for i = 0, a.n do r[i] = a[i] * c end; return r end
local function pmul(a, b)
  local r = {n = a.n + b.n}
  for i = 0, r.n do r[i] = ZERO end
  for i = 0, a.n do for j = 0, b.n do r[i + j] = r[i + j] + a[i] * b[j] end end
  return ptrim(r)
end
local function peval(p, x)
  local v = ZERO
  for i = p.n, 0, -1 do v = v * x + p[i] end
  return v
end
local function pdivmod(a, b)
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
local function piszero(p) p = ptrim(p); return p.n == 0 and cabs(p[0]) < 1e-12 end
-- roots via Durand-Kerner + Newton polish
local function proots(p)
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
local function rootGroups(roots)
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
local function pfe(R, D)
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
local function pstrQ(p, mode, off)
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
local function pstr(p, v, wrapv)
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
