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
        if not ln.startswith(">") and not ln.lstrip().startswith("|"):   # notes and tables stay in the README
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
    cd = f"cd('<path-to>/mmacos/templates/{s['section']}/{name}');"
    if "refusal" in s:                             # the bare call refuses today: that is printed first
        out += ["**Run it, as it stands today.**", "", "```matlab", cd, s["refusal"]["entry"], "```", ""]
        out += ["**It refuses:**", "", f"<!-- transcript: p2_{name}_refusal -->", ""]
    if "entry" in s:
        lead = "**Run it with the README's own value.**" if "refusal" in s else "**Run it.**"
        out += [lead, "", "```matlab", cd, s["entry"], "```", ""]
    if "run_note" in s:
        out += [s["run_note"], ""]
    if "knobs" in s:
        out += [f"**The knobs** (`{s['knobs']['file']}`):", "", f"<!-- transcript: p2_{name}_knobs -->", ""]
    if "run_ref" in s:                             # a template Part I already ran: point to it
        out += [f"**What it prints:** {s['run_ref']}", ""]
    elif "run" in s:
        out += ["**What it prints:**", "", f"<!-- transcript: p2_{name}_run -->", ""]
    if "figure_ref" in s:                          # its figure is already in Part I: point to it
        out += [f"**Its figure:** {s['figure_ref']}", ""]
    else:
        attr = s.get("fig_attr", "")               # e.g. "{height=8in}" for a tall figure
        fig = s.get("fig_media", f"p2_{name}_fig.png")  # or a figure a Part I driver generated
        out += [f"![](media/{fig}){attr}", "", f"**FIGURE {fig_no}** {s['caption']}", ""]
    out += [f"**The number it reproduces.** {s['number']}", ""]
    out += [f"**Adapt it.** *\"{s['prompt']}\"*", ""]
    return "\n".join(out)


TOOLS = MAN / "tools" / "part2" / "tools.txt"
LIBDIRS = [("design/src", "the design layer's helpers"), ("sensitivities", "the sensitivity drivers")]

TOOLS_INTRO = """## SECTION 9: The Tool Templates

The other seventeen templates are tools, not designs: examples of an interface, a method, or a driver that
the design templates and the runners call. Two are the entry points of Section 1's thread and get a page each;
the table lists all seventeen, each described by its own README's title or, where it has none, by its main
script's header. A last table lists the library functions no template calls.
"""


def describe(section, name):
    """A tool's one-line description: its README's title, else its main script's heading (generated)."""
    d = TPL / section / name
    if (d / "README.md").is_file():
        return readme_intro(name, section)[0], "README"
    sys.path.insert(0, str(MAN / "tools"))
    import gen_readme_part2 as gr
    ms = sorted(d.glob("*.m"), key=lambda f: (name not in f.stem, "example" not in f.stem, f.name))  # its own, then its example
    hl = gr.header_lines(ms[0])
    h1 = re.match(r"^([A-Z][A-Z0-9_]+)\s+(\S.*)$", hl[0]) if hl else None
    strip = lambda x: re.sub(r"^\S+\.m\s*(\([^)]*\)\s*)?(--\s*)?", "", x.strip()).strip().rstrip(".")   # "name.m -- what"
    text = h1.group(2) if h1 else next(strip(x) for x in hl if strip(x) and not re.fullmatch(r"[=\-%\s]+", x))
    text = strip(text) if not h1 else text.strip().rstrip(".")
    text = gr.heading_case(text) if re.sub(r"\([^)]*\)", "", text).isupper() else text
    return text, f"`{ms[0].name}` header"


def h1_line(f):
    """A MATLAB function's H1 line (the first comment line, '%NAME  summary'), '' if it has none."""
    for ln in f.read_text(errors="ignore").splitlines():
        t = ln.strip()
        if t.startswith("%"):
            m = re.match(r"^%\s*[A-Za-z0-9_.]+\s+(.*)$", t)
            return (m.group(1) if m else t.lstrip("% ")).strip()
        if t and not t.startswith("function"):
            return ""
    return ""


def uncalled():
    """Library functions no template's .m names (word-boundary match), by folder."""
    tpl = "\n".join(f.read_text(errors="ignore") for f in TPL.rglob("*.m"))
    out = []
    for d, label in LIBDIRS:
        fs = sorted((TPL.parent / d).glob("*.m"))
        un = [f for f in fs if not re.search(r"\b" + re.escape(f.stem) + r"\b", tpl)]
        out.append((d, label, len(fs), un))
    return out


def tools_tables():
    rows = [n.strip() for n in TOOLS.read_text().splitlines() if n.strip() and not n.startswith("#")]
    out = [TOOLS_INTRO, "### The seventeen tools", "", "| Template | Section | What it is | Described by |", "|---|---|---|---|"]
    for r in rows:
        sec, name = r.split("/")
        what, src = describe(sec, name)
        out.append(f"| `{name}` | {sec} | {what.replace('|', '/')} | {src} |")
    out += ["", "`40_benches/vsg_wip` is a planning folder (a VSG2 model plan and review notes), not runnable; it has no page.", "", ""]
    hl = ["### The functions no template calls", ""]
    parts = uncalled()
    hl.append("Counted by a word-boundary search of every template script for each function's name: "
              + "; ".join(f"`{d}`, {len(un)} of {n}" for d, _, n, un in parts) + ". They are documented here, not "
              "deleted (the Part I report's default); a runner or a test may still call them.")
    hl.append("")
    for d, label, n, un in parts:
        hl += [f"**`{d}`** ({label}, {len(un)} of {n}):", "", "| Function | What its help line says |", "|---|---|"]
        hl += [f"| `{f.stem}` | {h1_line(f).replace('|', '/') or '(no help line)'} |" for f in un]
        hl.append("")
    return "\n".join(out), "\n".join(hl)


def main():
    names = [n.strip() for n in ORDER.read_text().splitlines() if n.strip() and not n.startswith("#")]
    (MAN / "src" / "200_part2.md").write_text(INTRO)
    for old in (MAN / "src").glob("2[0-9][0-9]_*.md"):
        if old.name != "200_part2.md":
            old.unlink()
    fig = FIRST_FIGURE
    for i, n in enumerate(names, start=1):
        (MAN / "src" / f"2{i:02d}_{n}.md").write_text(page(n, fig) + "\n")
        if "figure_ref" not in json.loads((SPECS / f"{n}.json").read_text()):
            fig += 1
    tools_md, helpers_md = tools_tables()
    (MAN / "src" / "290_tools.md").write_text(tools_md + "\n")
    tnames = ["design_layer_api", "pupil_id"]
    for i, n in enumerate(tnames, start=1):
        (MAN / "src" / f"29{i}_{n}.md").write_text(page(n, fig).replace("### `", "### `", 1) + "\n")
        if "figure_ref" not in json.loads((SPECS / f"{n}.json").read_text()):
            fig += 1
    (MAN / "src" / "299_helpers.md").write_text(helpers_md + "\n")
    print(f"wrote src/200_part2.md, {len(names)} pages, the tools section with {len(tnames)} pages (FIGURE {FIRST_FIGURE}..{fig - 1})")


if __name__ == "__main__":
    sys.exit(main())
