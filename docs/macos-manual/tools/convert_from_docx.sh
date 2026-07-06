#!/bin/bash
# One-time migration (2026-07-06): styled .docx -> Markdown source tree.
# Kept for provenance; src/ is now the canonical source and this script
# should not normally be re-run (it would overwrite src/ edits).
#
# Usage: bash tools/convert_from_docx.sh <input.docx> <workdir>
set -euo pipefail
IN=${1:?input docx}
WORK=${2:?work dir}
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$WORK" && cd "$WORK"
pandoc -f docx+styles \
       --lua-filter="$HERE/style_map.lua" \
       --extract-media=. \
       --wrap=preserve \
       -t markdown-smart \
       -o full_mapped.md \
       "$IN"
python3 "$HERE/split_sections.py" full_mapped.md "$HERE/../src"
cp -r media "$HERE/../src/"
python3 "$HERE/gen_appendix_c.py"
echo "src/ regenerated from $IN"
