#!/usr/bin/env bash
# BRIEF_to_gauge_close -- the one queue for this box (32 GB: ONE model-1024
# MATLAB at a time).  Supersedes closechain2.sh, which was stopped to insert
# vampseq (item 1's fifth leg; see that script for the measurement that made
# it necessary).  Remaining order:
#   vampseq   item 1: v_cal 'amp', the reference frames every bench takes
#   gateseq3  item 3: the tail tuner's winner gate, two verifies
#   item2seq  item 2: unwrapped ladder, lens control, loop, descent
set -u
Z=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../zwfs_dm96/runs" && pwd)
T=$(cd "$(dirname "${BASH_SOURCE[0]}")/." && pwd)
n=0
while ! grep -q '\[vfitseq\] done' "$Z/vfitseq.nohup" 2>/dev/null; do
    n=$((n+1)); [ $n -gt 240 ] && { echo "[closechain3] ABORT: vfitseq never finished (2 h)"; exit 1; }
    sleep 30
done
echo "[closechain3] fit leg done -> the amp leg"
"$Z/vampseq.sh"  > "$Z/vampseq.nohup"   2>&1
echo "[closechain3] item 1 done -> item 3 gate verifies"
"$T/gateseq3.sh" > "$T/gateseq3.nohup"  2>&1
echo "[closechain3] item 3 done -> item 2 runs"
"$T/item2seq.sh" > "$T/item2seq.nohup"  2>&1
echo "[closechain3] all queued items finished"
