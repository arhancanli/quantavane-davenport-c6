#!/bin/sh
# One command:  sh check.sh            (the certificate, both of our programs, our red controls)
#               sh check.sh --referee  (also the referee's two programs and twelve controls)
# Needs a C compiler (cc) and Python 3 with numpy (set PYTHON=/path/to/python if numpy is not in python3).
# Program A needs about 3.3 GB of memory and 25 s; the whole run takes about 3 minutes on one core.
set -eu
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
fail() { echo "CHECK FAILED: $*"; exit 1; }
"$PY" -c 'import numpy' 2>/dev/null || fail "Python with numpy is needed (pip install numpy, or PYTHON=/path/to/python sh check.sh)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

echo "1. File hashes"
"$PY" -I tools/sha256check.py certs/SHA256SUMS || fail "certificate hashes"

echo "2. Rebuild the 66 terms from the sequences printed in arXiv 2611.24775 and compare byte for byte"
mkdir -p "$T/b/src" "$T/b/seeds" "$T/b/certs"
cp src/build_c6_13.py "$T/b/src/"; cp seeds/paper_sequences_2611.24775.json "$T/b/seeds/"
"$PY" -I "$T/b/src/build_c6_13.py" > "$T/build.log" 2>&1 || { cat "$T/build.log"; fail "build"; }
cmp -s "$T/b/certs/c6_13_len66.txt" certs/c6_13_len66.txt || fail "rebuilt .txt differs"
cmp -s "$T/b/certs/c6_13_len66.json" certs/c6_13_len66.json || fail "rebuilt .json differs"
echo "   identical"

echo "3. Program A (subset sums over the whole group C_6^13)"
cc -O2 -o "$T/checkA" src/checkA_fulldp.c
"$T/checkA" 13 certs/c6_13_len66.txt 2>&1 | grep -v "added$" | tee "$T/A.txt"
grep -q "ZERO-SUMFREE; distinct nonempty subsums 13060694015 of 13060694016" "$T/A.txt" || fail "Program A"

echo "4. Program B (exact enumeration over the three blocks)"
"$PY" -I src/checkB_blocks.py certs/c6_13_len66.json | tee "$T/B.txt"
grep -q "^checkB: ZERO-SUMFREE" "$T/B.txt" || fail "Program B"

echo "5. Red controls: each must be refused by both programs"
for c in flipz vzero append oneline; do
  "$T/checkA" 13 controls/ctl_$c.txt 2>&1 | grep -q "NOT zero-sumfree" || fail "Program A accepted control $c"
  "$PY" -I src/checkB_blocks.py controls/ctl_$c.json | grep -q "^checkB: NOT zero-sumfree" || fail "Program B accepted control $c"
  echo "   $c: refused by A and B"
done

if [ "${1:-}" = "--referee" ]; then
  echo "6. The referee's programs (written separately, before reading ours)"
  cc -O2 -o "$T/zsf_dp" referee/zsf_dp.c -lpthread
  "$T/zsf_dp" certs/c6_13_len66.txt --threads 2 2>/dev/null | tee "$T/R.txt" | grep "VERDICT"
  grep -q "^VERDICT: ZERO-SUMFREE" "$T/R.txt" || fail "referee zsf_dp"
  "$PY" -I referee/mitm.py certs/c6_13_len66.txt | tee "$T/M.txt" | grep "VERDICT"
  grep -q "^VERDICT (mitm): ZERO-SUMFREE" "$T/M.txt" || fail "referee mitm"
  for f in referee/controls/*.txt; do
    "$T/zsf_dp" "$f" --threads 2 2>/dev/null | grep -q "^VERDICT: NOT zero-sumfree" || fail "zsf_dp accepted $f"
    "$PY" -I referee/mitm.py "$f" | grep -q "^VERDICT (mitm): NOT zero-sumfree" || fail "mitm accepted $f"
  done
  echo "   referee programs: witness ZERO-SUMFREE, all 12 controls refused"
fi
echo "ALL CHECKS PASSED: the 66-term sequence over C_6^13 is zero-sumfree, so D(C_6^13) >= 67."
