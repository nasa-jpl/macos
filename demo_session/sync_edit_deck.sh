#!/bin/bash
# Sync Dave's editing copy + baseline from the current build -- GATED.
#   ./sync_edit_deck.sh <stem> [--folded | --force]      e.g. deck_gauges, deck_keysight
# Files: <stem>.pptx (the build), <stem>_edit.pptx (Dave edits this one),
# <stem>_baseline.pptx (what the edit copy was last synced from; the diff
# reference).  deck_keysight keeps its historical baseline name.
# Refuses when the edit deck differs from the baseline (unrecovered edits):
# fold them into the md / <stem>.geo.json sidecar first, rebuild, then re-run
# with --folded.  --force DISCARDS them -- only after Dave says the on-disk
# edits are disposable.  Never chain cp after a diff by hand: use this.
cd "$(dirname "$0")"
stem="${1:?usage: sync_edit_deck.sh <stem> [--folded|--force]}"; mode="${2:-}"
edit="${stem}_edit.pptx"; build="${stem}.pptx"; base="${stem}_baseline.pptx"
[ "$stem" = "deck_keysight" ] && base="baseline_pass3.pptx"
if [ -f ".~lock.${edit}#" ]; then
    echo "REFUSED: ${edit} is OPEN in Impress (lock file present)."
    echo "Have Dave close it (or Save As a side file) first."
    exit 1
fi
[ -f "$build" ] || { echo "no build ${build}"; exit 1; }
if [ -f "$base" ] && [ -f "$edit" ]; then
    T=$(python3 pptx_text_diff.py "$base" "$edit" | wc -l)
    G=$(python3 pptx_geo_diff.py "$base" "$edit" | grep -cE "moved|font|only in")
    if [ "$T" -ne 0 ] || [ "$G" -ne 0 ]; then
        if [ "$mode" != "--force" ] && [ "$mode" != "--folded" ]; then
            echo "REFUSED: ${edit} has unrecovered changes (text $T lines, geo $G deliberate)."
            echo "Fold into the md/sidecar + rebuild, then re-run with --folded."
            echo "--force DISCARDS them (only on Dave's say-so)."
            exit 1
        fi
    fi
else
    T=0; G=0
fi
cp -p "$build" "$edit"
cp -p "$build" "$base"
echo "synced ${stem}: edit + baseline = build (text $T, geo $G at check time)"
