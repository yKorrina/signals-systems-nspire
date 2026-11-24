#!/bin/sh
# Builds dist/SignalsSystems.lua and dist/SignalsSystems.tns from src/.
#
# The source is split into files that share top-level locals. Lua 5.1 only
# allows 200 locals per chunk, so when the files are joined the top-level
# "local" keywords are dropped. Comments and indentation are also stripped
# from the copy that goes into the .tns to keep it small.
#
# Needs Luna (https://github.com/ndless-nspire/Luna) for the .tns step.
# Set LUNA=/path/to/luna if it is not on your PATH.

set -e
cd "$(dirname "$0")"
LUNA=${LUNA:-luna}
mkdir -p dist

strip_locals() {
  sed -E '/^local [A-Za-z_][A-Za-z0-9_, ]*$/d' | sed -E 's/^local (function )?/\1/'
}

{
  echo 'platform.apilevel = "2.0"'
  cat src/0*.lua src/10_*.lua src/11_*.lua src/12_*.lua src/12b_*.lua \
      src/14_*.lua src/15_*.lua src/13_*.lua | strip_locals
} > dist/SignalsSystems.lua

sed -E 's/^[[:space:]]+//' dist/SignalsSystems.lua | sed -E '/^--/d' | sed -E '/^$/d' > dist/.min.lua

if command -v "$LUNA" >/dev/null 2>&1; then
  "$LUNA" dist/.min.lua dist/SignalsSystems.tns
  echo "built dist/SignalsSystems.tns"
else
  echo "luna not found, skipped .tns (dist/SignalsSystems.lua was still built)"
fi
rm -f dist/.min.lua
