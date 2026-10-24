#!/usr/bin/env python3
"""Generate Part II (the template catalog) of the mmacos manual: one page per template, from its
README (or, for a template without one, the README text the manual keeps in part2/<name>.md), its
spec (tools/part2/specs/<name>.json) and the artefacts its driver tools/transcripts/p2_<name>.m
produced (make transcripts).  Nothing on a page is typed except three short fields per spec:
'caption' (the figure's caption), 'number' (the sentence that states the reproduced number, with
%pin:...% references) and 'prompt' (a request a reader can paste to an assistant, one or two
sentences).  The page order is tools/part2/order.txt (the index's order).

  python3 tools/gen_part2.py        (from docs/mmacos-manual; writes src/200_part2.md, src/2NN_<name>.md)
"""
import json, re, sys
from pathlib import Path

MAN = Path(__file__).resolve().parent.parent
SPECS = MAN / "tools" / "part2" / "specs"
ORDER = MAN / "tools" / "part2" / "order.txt"
KEPT = MAN / "part2"                               # README text for templates without one
TPL = Path.home() / "dev" / "MACOS_resources" / "mmacos" / "templates"
FIRST_FIGURE = 20                                  # Part I ends at FIGURE 19

INTRO = """## SECTION 8: The Template Catalog

Part I followed one telescope through the library. This part is the library itself, one page per
design template, in the order of the template index (`templates/00_INDEX.md`). Each page says what
the template is (from its own README), shows the knobs it exposes, what it prints and draws when it
runs, and the number it reproduces, and ends with a request you can paste to an assistant to adapt
it. Every listing, figure and number on these pages comes from running the template, unchanged, the
way Part I's did (the front matter, *How this manual is made*).
"""


def readme_intro(name, section, start=None):
    """The README's title (after the dash) and its opening paragraphs, up to the first heading.
    Blockquoted notes (status banners, cross-references) are left out; a spec's 'intro_from' (a
    regex) starts the text at the first line it matches, for a README that opens with a table."""
    f = TPL / section / name / "README.md"
    kept = not f.is_file()
    if kept:
        f = KEPT / f"{name}.md"
    text = f.read_text()
    lines = text.splitlines()
    title = re.sub(r"^#\s*", "", lines[0]).strip()
    title = re.split(r"\s+[—-]\s+", title, maxsplit=1)[-1]
    rest = lines[1:]
    if start:
        k = next(i for i, ln in enumerate(rest) if re.search(start, ln))
        rest = rest[k:]
    body = []
    for ln in rest:
        if ln.startswith("#"):
            break
        if not ln.startswith(">"):
            body.append(ln)
    para = re.sub(r"\n{3,}", "\n\n", "\n".join(body)).strip()
    para = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", para)   # links -> their text (paths are not live here)
    if kept:
        para += ("\n\n*This template has no README; the text above is generated from its entry script's "
                 f"header (`docs/mmacos-manual/part2/{name}.md`, `tools/gen_readme_part2.py`).*")
    return title, para


def page(name, fig_no):
    s = json.loads((SPECS / f"{name}.json").read_text())
    title, intro = readme_intro(name, s["section"], s.get("intro_from"))
    out = [f"### `{name}`: {title}", "", intro, ""]
    out += ["**Run it.**", "", "```matlab",
            f"cd('<path-to>/mmacos/templates/{s['section']}/{name}');", s["entry"], "```", ""]
    if "run_note" in s:
        out += [s["run_note"], ""]
    if "knobs" in s:
        out += [f"**The knobs** (`{s['knobs']['file']}`):", "", f"<!-- transcript: p2_{name}_knobs -->", ""]
    if "run_ref" in s:                             # a template Part I already ran: point to it
        out += [f"**What it prints:** {s['run_ref']}", ""]
    else:
        out += ["**What it prints:**", "", f"<!-- transcript: p2_{name}_run -->", ""]
    attr = s.get("fig_attr", "")                   # e.g. "{height=8in}" for a tall figure
    fig = s.get("fig_media", f"p2_{name}_fig.png")  # or a figure a Part I driver generated
    out += [f"![](media/{fig}){attr}", "", f"**FIGURE {fig_no}** {s['caption']}", ""]
    out += [f"**The number it reproduces.** {s['number']}", ""]
    out += [f"**Adapt it.** *\"{s['prompt']}\"*", ""]
    return "\n".join(out)


def main():
    names = [n.strip() for n in ORDER.read_text().splitlines() if n.strip() and not n.startswith("#")]
    (MAN / "src" / "200_part2.md").write_text(INTRO)
    for old in (MAN / "src").glob("2[0-9][0-9]_*.md"):
        if old.name != "200_part2.md":
            old.unlink()
    fig = FIRST_FIGURE
    for i, n in enumerate(names, start=1):
        (MAN / "src" / f"2{i:02d}_{n}.md").write_text(page(n, fig) + "\n")
        fig += 1
    print(f"wrote src/200_part2.md and {len(names)} pages (FIGURE {FIRST_FIGURE}..{fig - 1})")


if __name__ == "__main__":
    sys.exit(main())
