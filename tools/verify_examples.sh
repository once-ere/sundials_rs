#!/bin/bash
# verify_examples.sh [crate|all] — build + run every [[example]] of the given
# crate(s) in release mode, diff each against the upstream SUNDIALS reference
# .out file, and write one line per example to logs/summary.txt:
#   <crate>/<example>  IDENTICAL | DIFF(n lines) | FAIL(code) | NOREF
# Read logs/summary.txt after running; open logs/<name>.diff only for
# non-IDENTICAL entries.
set -u
cd "$(dirname "$0")/.."
REFROOT="../sundials-7.7.0/examples"
mkdir -p logs
: > logs/summary.txt

refdir_for() {
  case "$1" in
    cvode_rs)  echo "$REFROOT/cvode/serial" ;;
    cvodes_rs) echo "$REFROOT/cvodes/serial" ;;
    kinsol_rs) echo "$REFROOT/kinsol/serial" ;;
    ida_rs)    echo "$REFROOT/ida/serial" ;;
    idas_rs)   echo "$REFROOT/idas/serial" ;;
    arkode_rs) echo "$REFROOT/arkode/C_serial" ;;
    *) echo "" ;;
  esac
}

# Filter machine-dependent noise symmetrically from both sides of the diff.
# Currently: nothing filtered (extend here if an example prints timings).
flt() { cat "$1"; }

verify_crate() {
  local crate="$1"
  local refdir; refdir="$(refdir_for "$crate")"
  local names
  names=$(grep -A1 '^\[\[example\]\]' "crates/$crate/Cargo.toml" 2>/dev/null \
          | grep 'name *= *"' | sed 's/.*"\(.*\)".*/\1/')
  [ -z "$names" ] && { echo "$crate: no examples" >> logs/summary.txt; return; }
  cargo build --release --examples -p "$crate" > "logs/build_$crate.log" 2>&1
  if [ $? -ne 0 ]; then
    echo "$crate: EXAMPLE BUILD FAILED (see logs/build_$crate.log)" >> logs/summary.txt
    return
  fi
  for name in $names; do
    local bin="target/release/examples/$name"
    if [ ! -x "$bin" ]; then
      echo "$crate/$name  FAIL(no binary)" >> logs/summary.txt; continue
    fi
    "$bin" > "logs/$name.out" 2>&1
    local code=$?
    if [ $code -ne 0 ]; then
      echo "$crate/$name  FAIL($code)" >> logs/summary.txt; continue
    fi
    local ref="$refdir/$name.out"
    if [ ! -f "$ref" ]; then
      echo "$crate/$name  NOREF" >> logs/summary.txt; continue
    fi
    if diff <(flt "logs/$name.out") <(flt "$ref") > "logs/$name.diff" 2>&1; then
      echo "$crate/$name  IDENTICAL" >> logs/summary.txt
      rm -f "logs/$name.diff"
    else
      echo "$crate/$name  DIFF($(wc -l < "logs/$name.diff" | tr -d ' ') lines)" >> logs/summary.txt
    fi
  done
}

if [ "${1:-all}" = "all" ]; then
  for c in cvode_rs cvodes_rs kinsol_rs ida_rs idas_rs arkode_rs; do
    grep -q '^\[\[example\]\]' "crates/$c/Cargo.toml" 2>/dev/null && verify_crate "$c"
  done
else
  verify_crate "$1"
fi
echo "---- summary written to logs/summary.txt ----"
