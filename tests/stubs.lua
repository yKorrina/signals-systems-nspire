-- Minimal stand-ins for the Nspire API so the app loads in desktop Lua 5.1.
platform = {apilevel = "2.0", window = {
  width = function() return 318 end, height = function() return 212 end,
  invalidate = function() end}}
on = {}
timer = {start = function() end, stop = function() end}
toolpalette = {register = function() end}
clipboard = {getText = function() return "" end, addText = function() end}
