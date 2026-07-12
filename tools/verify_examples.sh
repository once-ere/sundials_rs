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
# Optional second reference tree: outputs of the locally-built C library
# (same layout as $REFROOT). When a run differs from the shipped .out but
# matches the local C build byte-for-byte, report LOCAL-C instead of DIFF.
# The committed localref/ tree holds these outputs (Release,
# -ffp-contract=off, no LAPACK/KLU — see VERIFICATION.md); earlier
# sessions kept them in an ephemeral scratchpad, which did not survive.
LOCALREF="${SUNDIALS_LOCALREF:-localref}"
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
  # Run one (example, args) pair and append its summary line.
  # $1 = crate, $2 = example name, $3 = reference file, $4 = output tag
  # (name or name_argsuffix), $5.. = decoded command-line args.
  verify_one() {
    local crate="$1" name="$2" ref="$3" tag="$4"; shift 4
    local bin="target/release/examples/$name"
    "$bin" "$@" > "logs/$tag.out" 2>&1
    local code=$?
    if [ $code -ne 0 ]; then
      echo "$crate/$tag  FAIL($code)" >> logs/summary.txt; return
    fi
    if diff <(flt "logs/$tag.out") <(flt "$ref") > "logs/$tag.diff" 2>&1; then
      echo "$crate/$tag  IDENTICAL" >> logs/summary.txt
      rm -f "logs/$tag.diff"
    else
      local lref="$LOCALREF/${refdir#$REFROOT/}/$(basename "$ref")"
      if [ -f "$lref" ] && diff <(flt "logs/$tag.out") <(flt "$lref") > /dev/null 2>&1; then
        echo "$crate/$tag  LOCAL-C (shipped-ref diff: $(wc -l < "logs/$tag.diff" | tr -d ' ') lines)" >> logs/summary.txt
      else
        echo "$crate/$tag  DIFF($(wc -l < "logs/$tag.diff" | tr -d ' ') lines)" >> logs/summary.txt
      fi
    fi
  }

  for name in $names; do
    local bin="target/release/examples/$name"
    if [ ! -x "$bin" ]; then
      echo "$crate/$name  FAIL(no binary)" >> logs/summary.txt; continue
    fi

    # Plain (no-argument) reference, when one is shipped.
    local ref="$refdir/$name.out"
    local have_any=0
    if [ -f "$ref" ]; then
      have_any=1
      verify_one "$crate" "$name" "$ref" "$name"
    fi

    # Argument-encoded reference variants: <name>_<argv-with-underscores>.out
    # (e.g. idasRoberts_FSA_dns_-sensi_stg_t.out -> args "-sensi stg t").
    # The underscore join is lossy for keys that contain underscores
    # (idas.init_step, kinsol m_aa, ...): tools/verify_args.map overrides
    # the decoding per reference basename.  Variant scanning is enabled
    # for idas_rs/arkode_rs only; the older crates' variants were
    # verified manually (see VERIFICATION.md) and several use ambiguous
    # encodings.  Skip files that belong to a LONGER example name
    # sharing this prefix.
    if [ "$crate" = "idas_rs" ] || [ "$crate" = "arkode_rs" ] || [ "$crate" = "cvodes_rs" ]; then
      local vref
      for vref in "$refdir/${name}_"*.out; do
        [ -f "$vref" ] || continue
        local base; base="$(basename "$vref" .out)"
        local suffix="${base#${name}_}"
        # skip refs owned by a LONGER example name sharing this prefix
        # (any number of extra underscore-joined tokens)
        local owned_by_longer=0 other
        for other in $names; do
          [ "$other" = "$name" ] && continue
          case "$base" in
            "$other" | "${other}_"*)
              if [ ${#other} -gt ${#name} ]; then owned_by_longer=1; fi ;;
          esac
        done
        # ... including LONGER C example names that are not ported yet
        # (e.g. ark_brusselator1D_imexmri_*.out are NOT arg variants of
        # ark_brusselator1D)
        if [ $owned_by_longer -eq 0 ]; then
          local csrc cbase
          for csrc in "$refdir/${name}_"*.c; do
            [ -f "$csrc" ] || continue
            cbase="$(basename "$csrc" .c)"
            [ "$cbase" = "$name" ] && continue
            case "$base" in
              "$cbase" | "${cbase}_"*)
                if [ ${#cbase} -gt ${#name} ]; then owned_by_longer=1; fi ;;
            esac
          done
        fi
        [ $owned_by_longer -eq 1 ] && continue
        local args
        args="$(awk -v k="$base" '$1 == k { $1=""; print substr($0,2); exit }' \
                tools/verify_args.map 2>/dev/null)"
        [ -z "$args" ] && args="$(printf '%s' "$suffix" | tr '_' ' ')"
        have_any=1
        # shellcheck disable=SC2086
        verify_one "$crate" "$name" "$vref" "$base" $args
      done
    fi

    if [ $have_any -eq 0 ]; then
      echo "$crate/$name  NOREF" >> logs/summary.txt
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
