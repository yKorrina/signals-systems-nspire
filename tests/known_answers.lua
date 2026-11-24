-- Checks a set of problems against answers worked out by hand.
-- Usage (from the repo root, after ./build.sh):
--   cat tests/stubs.lua dist/SignalsSystems.lua tests/known_answers.lua > /tmp/t.lua && lua5.1 /tmp/t.lua
local cases = {
  {"Z-transform", {"d[n]+2d[n-1]+2d[n-2]+d[n-3]"}, "X(z) = 1 + 2z^-1 + 2z^-2 + z^-3"},
  {"DTFT", {"(1/2)^n u[n]"}, "X(e^(jΩ)) = 1/(1 - (1/2)e^(-jΩ))"},
  {"DTFT", {"u[n-2]-u[n-6]"}, "X(e^(jΩ)) = e^(-j2Ω) + e^(-j3Ω) + e^(-j4Ω) + e^(-j5Ω)"},
  {"Fourier transform (CT)", {"e^(-2t)u(t)"}, "X(jω) = 1/(2 + jω)"},
  {"Laplace transform", {"e^(-t)u(t)"}, "X(s) = 1/(s + 1)"},
  {"Convolution", {"e^(-t)u(t)", "e^(-2t)u(t)"}, "y(t) = (e^(-t) - e^(-2t))u(t)"},
  {"Convolution", {"(1/2)^n u[n]", "(1/3)^n u[n]"}, "y[n] = (3·(1/2)^n - 2·(1/3)^n)u[n]"},
  {"Convolution", {"u(t)", "u(t)"}, "y(t) = t·u(t)"},
  {"Running integral / sum", {"e^(-t)u(t)"}, "y(t) = (1 - e^(-t))u(t)"},
  {"Running integral / sum", {"(1/2)^n u[n]"}, "y[n] = (2 - (1/2)^n)u[n]"},
  {"Inverse Fourier (CT)", {"1/((jw+1)(jw+2))"}, "x(t) = (e^(-t) - e^(-2t))u(t)"},
  {"Inverse Laplace", {"(s+1)/(s^2+2s+5)", ""}, "x(t) = e^(-t)cos(2t)u(t)"},
  {"Inverse Z", {"z/(z-1/2)", ""}, "x[n] = (1/2)^n·u[n]"},
  {"Inverse DTFT", {"1/((1-1/2e^(-jW))(1-1/3e^(-jW)))"}, "x[n] = (3·(1/2)^n - 2·(1/3)^n)u[n]"},
  {"Complex numbers", {"(3+j)/(1-j)"}, "Cartesian: 1 + j2"},
  {"Integral / sum evaluator", {"e^(-6abs(t))", "-inf", "inf"}, "∫[-∞→∞] f = 1/3"},
  {"System properties", {"y(t)=x(t-1)+x(t+1)"}, "Linear: YES   Causal: NO"},
  {"System properties", {"y[n]=x[n]u[n]"}, "Memoryless: YES   Time-inv: NO"},
  {"Frequency response H", {"y'+2y=x'+x", ""}, "H(jω) = Y(jω)/X(jω) = (jω + 1)/(jω + 2)"},
}
local bad = 0
for _, c in ipairs(cases) do
  local tool
  for _, t in ipairs(TOOLS) do if t.name == c[1] then tool = t end end
  local out = newOut()
  local vals = {}
  for k = 1, #tool.fields do vals[k] = c[2][k] or "" end
  pcall(tool.run, vals, buildCtx(""), out)
  local found, failedCheck = false, false
  for _, l in ipairs(out.lines) do
    if l.s == c[3] then found = true end
    if l.s:find("Check FAILED", 1, true) then failedCheck = true end
  end
  if not found or failedCheck then bad = bad + 1; print("FAIL: " .. c[1] .. "  " .. table.concat(c[2], " | ")) end
end
print(bad == 0 and ("all " .. #cases .. " known answers ok") or (bad .. " failed"))
