-- Runs every built-in example through its tool and reports failures.
-- Usage (from the repo root, after ./build.sh):
--   cat tests/stubs.lua dist/SignalsSystems.lua tests/run_examples.lua > /tmp/t.lua && lua5.1 /tmp/t.lua
local bad = 0
for i, ex in ipairs(EXAMPLES) do
  local ctx = buildCtx(ex.defs or "")
  local tool
  for _, t in ipairs(TOOLS) do if t.name == ex.tool then tool = t end end
  local out = newOut()
  local vals = {}
  for k = 1, #tool.fields do vals[k] = ex.vals[k] or "" end
  local ok, err = pcall(tool.run, vals, ctx, out)
  local last = ""
  for _, l in ipairs(out.lines) do
    if l.k == "e" then ok = false; err = l.s end
    if l.k == "r" then last = l.s end
  end
  if not ok then bad = bad + 1 end
  print(string.format("%3d %-18s %s", i, ex.cat, ok and last or ("FAILED: " .. tostring(err))))
end
print(bad == 0 and "all examples ok" or (bad .. " failed"))
