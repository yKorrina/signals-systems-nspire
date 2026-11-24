-- ===================== NSPIRE USER INTERFACE =====================
local DEFAULT_DEFS = {"x(t)=e^(-t)u(t)", "h(t)=e^(-2t)u(t)", "a=2", "", "", ""}
local App = {
  screen = "home", sel = 1, top = 1,
  tool = nil, field = 1, cur = {}, vals = {}, defs = {},
  out = nil, rl = nil, rscroll = 0, rwidth = 0,
  plots = {}, pidx = 1, g = nil, busy = false, msg = nil, phase = false,
}
for i, d in ipairs(DEFAULT_DEFS) do App.defs[i] = d end

local HOME_ITEMS = {}
local function buildHome()
  HOME_ITEMS = {
    {label = "Search tools / examples / tables", kind = "search"},
    {label = "Example problems (" .. #EXAMPLES .. ")", kind = "examples"},
    {label = "Reference tables (pairs, properties)", kind = "refs"},
    {label = "Definitions  x(t), h(t), params", kind = "defs"},
  }
  for i, t in ipairs(TOOLS) do HOME_ITEMS[#HOME_ITEMS + 1] = {label = i .. ". " .. t.name, kind = "tool", idx = i} end
  HOME_ITEMS[#HOME_ITEMS + 1] = {label = "Syntax help & keys", kind = "help"}
end

local HELP = {
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
local function W() return platform.window:width() end
local function Hh() return platform.window:height() end
local function invalidate() platform.window:invalidate() end
local function setFont(gc, st, sz) gc:setFont("sansserif", st or "r", sz or 9) end
local function col(gc, r, g, b) gc:setColorRGB(r, g, b) end
local function sw(gc, s) return gc:getStringWidth(s) end

local function wrapText(gc, s, maxw)
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

local function currentTool() return TOOLS[App.tool] end
local function fieldsOf()
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
local function curPos(i, s)
  local key = (App.screen == "defs" and "d" or App.tool) .. ":" .. i
  if App.cur[key] == nil then App.cur[key] = #uchars(s) end
  return App.cur[key], key
end

-- ---------- run ----------
local function defsText() return table.concat(App.defs, "\n") end
local function runTool()
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
local function header(gc, title, right)
  col(gc, 0, 0, 0); setFont(gc, "b", 10)
  gc:drawString(title, 4, 1, "top")
  if right then col(gc, 110, 110, 110); setFont(gc, "r", 7); gc:drawString(right, W() - sw(gc, right) - 4, 4, "top") end
  col(gc, 160, 160, 160); gc:drawLine(0, 17, W(), 17)
end
local function footer(gc, s)
  col(gc, 255, 255, 255); gc:fillRect(0, Hh() - 13, W(), 13)
  col(gc, 160, 160, 160); gc:drawLine(0, Hh() - 13, W(), Hh() - 13)
  col(gc, 90, 90, 90); setFont(gc, "r", 7)
  gc:drawString(s, 3, Hh() - 12, "top")
end

local function paintHome(gc)
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

local function paintFields(gc, title, help)
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

local KCOL = {h = {0, 0, 0}, t = {0, 0, 0}, r = {0, 0, 0}, e = {0, 0, 0}}
local function paintResult(gc)
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
local function niceStep(range, n)
  local raw = range / n
  local mag = 10 ^ math.floor(math.log10(raw))
  local r = raw / mag
  local st
  if r < 1.5 then st = 1 elseif r < 3 then st = 2 elseif r < 7 then st = 5 else st = 10 end
  return st * mag
end
local function plotRange(p)
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
local evalPlotRaw
local function evalPlot(p, x)
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
local function graphInit()
  local p = App.plots[App.pidx]
  local lo, hi = plotRange(p)
  App.g = {x0 = lo, x1 = hi, tr = nil, y0 = nil, y1 = nil}
end
local function graphAutoY(p, g, N, xs)
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

local function paintGraph(gc)
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

local function paintList(gc)
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
local function paintSearch(gc)
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

local openItem
local function loadExample(ex)
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
local function openList(title, items, back)
  App.list = {title = title, items = items, sel = 1, top = 1, back = back or "home"}
  for i, it in ipairs(items) do if it.kind ~= "hdr" then App.list.sel = i; break end end
  App.screen = "list"
  invalidate()
end
local function exampleItems()
  local items, last = {}, nil
  for _, ex in ipairs(EXAMPLES) do
    if ex.cat ~= last then items[#items + 1] = {label = ex.cat, kind = "hdr"}; last = ex.cat end
    items[#items + 1] = {label = "  " .. ex.title, kind = "example", ex = ex}
  end
  return items
end
local function refItems()
  local items = {}
  for _, r in ipairs(REFS) do items[#items + 1] = {label = r.title, kind = "ref", ref = r} end
  return items
end
local function openView(title, lines, back)
  App.view = {title = title, lines = lines, back = back}
  App.screen = "view"; App.rl = nil; App.rscroll = 0
  invalidate()
end
local function refLines(r)
  local L = {{k = "h", s = r.title}}
  for _, s in ipairs(r.lines) do L[#L + 1] = {k = "t", s = s} end
  return L
end
local function lower(s) return (s:gsub("%u", string.lower)) end
local function searchItems(q)
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

local function editField(fn)
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

local function searchUpdate()
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
local function moveField(d)
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
local function registerMenu()
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
