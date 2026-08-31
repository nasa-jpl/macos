#!/bin/bash
# Sync Dave's editing copy + baseline from the current build -- GATED.
# Refuses when the edit deck differs from the baseline (unrecovered edits);
# fold them into the md/sidecar first, then re-run.  --force overrides ONLY
# after Dave has confirmed the on-disk edits are disposable.
cd "$(dirname "$0")"
if [ -f ".~lock.deck_keysight_edit.pptx#" ]; then
    echo "REFUSED: deck_keysight_edit.pptx is OPEN in Impress (lock file present)."
    echo "Have Dave close it (or Save As a side file) first."
    exit 1
fi
T=$(python3 pptx_text_diff.py baseline_pass3.pptx deck_keysight_edit.pptx | wc -l)
G=$(python3 pptx_geo_diff.py baseline_pass3.pptx deck_keysight_edit.pptx | grep -cE "moved|font|only in")
if [ "$T" -ne 0 ] || [ "$G" -ne 0 ]; then
    if [ "$1" != "--force" ] && [ "$1" != "--folded" ]; then
        echo "REFUSED: edit deck has unrecovered changes (text $T lines, geo $G deliberate)."
        echo "Fold into md/sidecar + rebuild, then re-run with --folded."
        echo "--force DISCARDS them (only on Dave's say-so)."
        exit 1
    fi
fi
cp deck_keysight.pptx deck_keysight_edit.pptx
cp deck_keysight.pptx baseline_pass3.pptx
echo "synced (text $T, geo $G at check time)"
