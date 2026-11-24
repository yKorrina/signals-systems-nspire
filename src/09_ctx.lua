-- ===================== CONTEXT (definitions & params) =====================
local function buildCtx(defsText)
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
local function P(ctx, s, extra)
  local fn = {}
  for k in pairs(ctx.fnames) do fn[k] = true end
  for k in pairs(extra or {}) do fn[k] = true end
  return parse(s, {fnames = fn})
end
-- decide DT vs CT for an expression
local function isDTexpr(node, ctx)
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
local function toSig(ctx, node, dt)
  local var = dt and "n" or "t"
  local cv = mkConv(var, dt, ctx.env, ctx.defs)
  return cv:sig(node)
end
-- numeric function of a node (for plotting/fallback)
local function numFn(ctx, node, var)
  return function(x)
    ctx.env.vars[var] = cx(x)
    ctx.env.tvar = var
    local ok, v = pcall(evalN, node, ctx.env)
    if ok then return v end
    return nil
  end
end
