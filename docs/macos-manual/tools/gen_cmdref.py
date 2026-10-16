#!/usr/bin/env python3
"""Regenerate the auto-skeleton parts of the MACOS Command Reference:

  cmdref/10_engine_commands.md   from macos_f90/macos_help.inc
  cmdref/20_bindings.md          from mmacos +macos/*.m and pymacos
                                 macos.py/toolbox.py (paired by name)
  cmdref/30_higher_level.md      from mmacos +design and +channels

Hand-written prose survives regeneration ONLY in blocks bracketed by
  <!-- BEGIN NOTES <entry-key> --> ... <!-- END NOTES <entry-key> -->
Everything else in these three files is overwritten.  The orientation
and per-surface syntax chapters (00_-04_) are hand-maintained and not
touched.

Run from docs/macos-manual:  python3 tools/gen_cmdref.py
"""
import ast
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANUAL = HERE.parent
CMDREF = MANUAL / "cmdref"
ROOT = MANUAL.parent.parent                      # repo root
HELP_INC = ROOT / "macos_f90" / "macos_help.inc"
API_MOD = ROOT / "macos_f90" / "macos_api_mod.F90"
MMACOS = ROOT.parent / "MACOS_resources" / "mmacos" / "src" / "+macos"
PYMACOS = ROOT.parent / "MACOS_resources" / "pymacos" / "src" / "pymacos"
CLI_CMD = ROOT / "cli_tests" / "cmd"          # the CLI suite's journals (Phase B)

TODO = "*TODO: expand — dialog details, behavior notes, related commands.*"


# ----------------------------------------------------------------- notes
def load_notes(path):
    """Harvest hand-written NOTES blocks from an existing generated file
    so regeneration preserves them."""
    notes = {}
    if not path.exists():
        return notes
    pat = re.compile(
        r"<!-- BEGIN NOTES (\S+) -->\n(.*?)<!-- END NOTES \1 -->",
        re.S)
    for m in pat.finditer(path.read_text()):
        body = m.group(2).rstrip("\n")
        if body.strip():          # empty block = unfilled; keep its TODO
            # keys are matched case-insensitively so a capitalization
            # fix in macos_help.inc can't orphan an entry's notes
            notes[m.group(1).lower()] = body
    return notes


def notes_block(key, notes):
    body = notes.get(key.lower(), "")
    return (f"<!-- BEGIN NOTES {key} -->\n{body}\n"
            f"<!-- END NOTES {key} -->" if body else
            f"<!-- BEGIN NOTES {key} -->\n<!-- END NOTES {key} -->")


# ---------------------------------------------------------- engine part
def parse_help_inc():
    """Return [(category, [(names, desc)])] preserving help order."""
    line_pat = re.compile(r"WRITE\(\*,\*\)'(.*)'\s*$")
    cats, cur, cur_items = [], None, []
    pending = None
    for raw in HELP_INC.read_text().splitlines():
        m = line_pat.search(raw)
        if not m:
            continue
        s = m.group(1).rstrip()
        st = s.strip()
        if not st:
            continue
        if st.isupper() and " - " not in st and len(st.split()) <= 5 \
                and not st.startswith(("MACOS", "TAGS", "CASING")):
            if pending:
                cur_items.append(pending)
                pending = None
            if cur is not None:
                cats.append((cur, cur_items))
            cur, cur_items = st.title(), []
            continue
        if st.startswith(("MACOS command", "Casing:", "Tags:")):
            continue
        cm = re.match(r"^(\S[^-]*?)\s+-\s+(.*)$", st)
        if cm and cur is not None:
            if pending:
                cur_items.append(pending)
            pending = [cm.group(1).strip(), cm.group(2).strip()]
        elif pending:
            pending[1] += " " + st            # continuation line
    if pending:
        cur_items.append(pending)
    if cur is not None:
        cats.append((cur, cur_items))
    return cats


def _entry_key(names):
    return "cmd-" + re.sub(r"[^A-Za-z0-9]+", "-",
                           names.split(",")[0]).strip("-")


def journal_index(cats):
    """Phase B (PLAN_CONSOLIDATION 1.5): map each engine entry to the CLI-suite
    journals (cli_tests/cmd/*.jou) that issue it.  A journal line that is not a
    comment is a command or a prompt answer; its first token is matched the way
    the dispatcher matches -- it must start with the entry's minimum-match
    prefix (the capitalized run of the HELP name) and be a prefix of the full
    name.  Returns {entry key (lower): [journal names]}."""
    entries = []
    for _cat, items in cats:
        for names, _desc in items:
            forms = []
            for nm in names.split(","):
                word = nm.strip().split()[0] if nm.strip() else ""
                word = re.sub(r"[<>\[\]]", "", word)
                if not word:
                    continue
                mm = re.match(r"^[A-Z!]+", word)
                mm = mm.group(0) if mm else word.upper()
                forms.append((mm.upper(), word.upper()))
            entries.append((_entry_key(names).lower(), forms))
    index = {}
    if not CLI_CMD.is_dir():
        return index
    for jou in sorted(CLI_CMD.glob("*.jou")):
        hits = set()
        for line in jou.read_text(errors="replace").splitlines():
            st = line.strip()
            if not st or st.startswith("#"):
                continue
            tok = st.split()[0]
            if not re.match(r"^[A-Za-z!][A-Za-z]*$", tok):
                continue
            tu = tok.upper()
            for key, forms in entries:
                for mm, full in forms:
                    if len(tu) >= len(mm) and full.startswith(tu):
                        hits.add(key)
        for key in hits:
            index.setdefault(key, []).append(jou.name)
    return index


def gen_engine(notes):
    cats = parse_help_inc()
    jindex = journal_index(cats)
    uncovered = []
    out = ["<!-- GENERATED by tools/gen_cmdref.py from"
           " macos_f90/macos_help.inc.",
           "     Edit prose only inside BEGIN/END NOTES blocks;"
           " everything else is overwritten. -->",
           "",
           "## Part I — Engine command catalog",
           "",
           "Commands below are available from the interactive CLI and, "
           "with few exceptions, from a `CALL SMACOS(command, ...)` "
           "invocation (see the SMACOS chapter for the argument "
           "convention; purely interactive conveniences — shell/file "
           "utilities like `!`, PWD, CD, RX/LS/LL, VI/EMAcs, and a few "
           "others noted per entry — exist only in the CLI build).  "
           "The **capitalized prefix** of each name is the "
           "minimum-match abbreviation.  Prerequisite tags: \\[Rx\\] "
           "needs a loaded prescription (OLD/NEW), \\[BLD\\] a built "
           "linear model (BUILD), \\[DIFF\\] a diffraction output "
           "(the command itself propagates the wavefront; there is no "
           "separate PROpagate command).",
           ""]
    n = 0
    for cat, items in cats:
        out.append(f"### {cat}")
        out.append("")
        for names, desc in items:
            n += 1
            key = "cmd-" + re.sub(r"[^A-Za-z0-9]+", "-",
                                  names.split(",")[0]).strip("-")
            tags = " ".join(re.findall(r"\[(?:Rx|BLD|DIFF)\]", desc))
            desc_clean = re.sub(r"\s*\[(?:Rx|BLD|DIFF)\]", "", desc)
            out.append(f"#### {names}")
            out.append("")
            meta = f"*Minimum match:* `{names.split(',')[0].strip()}`"
            if tags:
                meta += f" — *needs:* {tags}"
            out.append(meta)
            out.append("")
            out.append(desc_clean[:1].upper() + desc_clean[1:]
                       + ("" if desc_clean.endswith(".") else "."))
            out.append("")
            jl = jindex.get(key.lower())
            if jl:
                out.append("*CLI suite:* " + ", ".join(f"`{j}`" for j in jl)
                           + " (`cli_tests/cmd/`; each journal is the "
                           "entry's example and its regression test).")
            else:
                out.append("*CLI suite:* none yet.")
                uncovered.append((cat, names))
            out.append("")
            nb = notes_block(key, notes)
            out.append(nb)
            if key.lower() not in notes:
                out.append(TODO)
            out.append("")
    (CMDREF / "10_engine_commands.md").write_text("\n".join(out))
    # the Phase B worklist: engine commands no journal exercises yet
    cov = ["<!-- GENERATED by tools/gen_cmdref.py from cli_tests/cmd/*.jou"
           " and macos_f90/macos_help.inc; do not edit. -->",
           "",
           "## Appendix — CLI-suite coverage of the engine commands",
           "",
           f"Of the {n} engine commands, {n - len(uncovered)} are issued by at "
           f"least one journal of the CLI regression suite "
           f"(`cli_tests/cmd/*.jou`, run by `cli_tests/run_cli_tests.sh cmd` "
           f"on both compilers); each entry of Part I names its journals.  "
           f"The {len(uncovered)} commands below have no journal yet — the "
           f"Phase B worklist.  A journal is one pty transcript with "
           f"`# EXPECT:` / `# EXPECT-NOT:` assertions, so writing one is "
           f"writing the entry's example and its test at once.",
           ""]
    bycat = {}
    for cat, names in uncovered:
        bycat.setdefault(cat, []).append(names)
    for cat in [c for c, _ in cats if c in bycat]:
        cov.append(f"- **{cat}:** " + ", ".join(bycat[cat]))
    cov.append("")
    (CMDREF / "40_suite_coverage.md").write_text("\n".join(cov))
    print(f"cli-suite coverage: {n - len(uncovered)} of {n} engine commands")
    return n


# --------------------------------------------------------- bindings part
def parse_api_mod():
    names = set()
    for m in re.finditer(r"^\s*(?:module\s+)?subroutine\s+([a-z0-9_]+)",
                         API_MOD.read_text(), re.M | re.I):
        names.add(m.group(1).lower())
    return names


def parse_mmacos_dir(d):
    """{name: (signature, [help lines])} for every top-level .m file."""
    out = {}
    for f in sorted(d.glob("*.m")):
        text = f.read_text(errors="replace").splitlines()
        sig, helplines = None, []
        for ln in text:
            s = ln.strip()
            if sig is None:
                if s.startswith("function"):
                    sig = re.sub(r"^function\s+", "", s).strip()
                elif s.startswith("classdef"):
                    sig = s
                continue
            if s.startswith("%"):
                helplines.append(s.lstrip("% ").rstrip())
            elif s:
                break
        out[f.stem] = (sig or f.stem, helplines)
    return out


def parse_pymacos():
    """{name: (signature, first docstring paragraph, alias_of)}"""
    out = {}
    for fname in ("macos.py", "toolbox.py"):
        p = PYMACOS / fname
        if not p.exists():
            continue
        tree = ast.parse(p.read_text())
        for node in tree.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            if node.name.startswith("_"):
                continue
            args = [a.arg for a in node.args.args]
            sig = f"{node.name}({', '.join(args)})"
            doc = ast.get_docstring(node) or ""
            para = doc.split("\n\n")[0].replace("\n", " ").strip()
            alias_of = None
            if not doc and len(node.body) == 1 \
                    and isinstance(node.body[0], ast.Return) \
                    and isinstance(node.body[0].value, ast.Call) \
                    and isinstance(node.body[0].value.func, ast.Name):
                alias_of = node.body[0].value.func.id
            out[node.name] = (sig, para, alias_of)
    return out


GROUPS = [
    ("Session, prescription and model I/O",
     r"^(init|load_rx|save_rx|has_rx|modify|num_elt|sys_units|cbm|"
     r"modified_rx|Session|save|reset)"),
    ("Source and wavelength",
     r"(^|_)(src|source)(_|$)|^sxp"),
    ("Element geometry and properties",
     r"^(get_elt|set_elt|elt_|del_elt|find_.*_elts|prb_grp)"),
    ("Gratings", r"grating"),
    ("Ray trace and outputs",
     r"^(trace|spot|opd|intensity|complex_field|apodize|ffp|pfp|fex|"
     r"xps|window|stop|srs|get_ray|draw_rays|dx_at|"
     r"first_order_properties|pupil_quality|get_xp|set_xp|compose)"),
    ("Perturbation", r"^(perturb|prb_elt)"),
    ("Jacobians (sensitivity channels)", r"^dw|^dwd"),
    ("CALIB (native optimizer)", r"^calib"),
    ("Segment and grid bases",
     r"(grid_basis|write_grid_file|zrn_freeform|grid_add|grid_elts)"),
]


def norm(name):
    return re.sub(r"[^a-z0-9]", "", name.lower())


def group_for(name):
    for title, pat in GROUPS:
        if re.search(pat, name):
            return title
    return "Other utilities"


def gen_bindings(notes):
    api = parse_api_mod()
    mm = parse_mmacos_dir(MMACOS)
    py = parse_pymacos()

    # map normalized-name -> pymacos primary + aliases
    py_primary, py_alias = {}, {}
    for name, (sig, para, alias_of) in py.items():
        if alias_of and alias_of in py:
            py_alias.setdefault(alias_of, []).append(name)
        else:
            py_primary[norm(name)] = (name, sig, para)

    def base_of(name):
        return re.sub(r"^(get_?|set_?)", "", name, flags=re.I)

    # same operation, different name across the bindings:
    # norm(pymacos base) -> norm(mmacos base)
    SYNONYMS = {
        "tracerays": "trace",
        "load": "loadrx",
        "save": "saverx",
        "eltgrprm": "deleltgrp",
        "eltgridfnd": "findgridelts",
        "eltzrnfnd": "findzernelts",
        "rxmodified": "modifiedrx",
    }
    EXCLUDE_PY = {"isNumeric", "isfloat", "isinteger", "isnumeric"}

    # one entry per OPERATION: mmacos get_/set_ pairs and pymacos
    # combined getter/setters all merge onto the stripped base name
    entries = {}          # norm(base) -> dict
    for name, (sig, helplines) in mm.items():
        b = norm(base_of(name))
        e = entries.setdefault(b, dict(name=base_of(name), mm=[], py=None))
        e["mm"].append((name, sig, helplines))
    for k, (name, sig, para) in py_primary.items():
        if name in EXCLUDE_PY:
            continue
        b = norm(base_of(name))
        b = SYNONYMS.get(b, b)
        e = entries.setdefault(b, dict(name=base_of(name), mm=[], py=None))
        if e["py"] is None:
            e["py"] = (name, sig, para)
        if not e["mm"]:
            e["name"] = name
    entries = {e["name"]: e for e in entries.values()}

    grouped = {}
    for name, e in entries.items():
        grouped.setdefault(group_for(name), []).append((name, e))

    out = ["<!-- GENERATED by tools/gen_cmdref.py from the mmacos +macos"
           " package and the pymacos python package.",
           "     Edit prose only inside BEGIN/END NOTES blocks;"
           " everything else is overwritten. -->",
           "",
           "## Part II — mmacos / pymacos function catalog",
           "",
           "Each entry shows the MATLAB (mmacos) and Python (pymacos) "
           "form of the same operation where both exist; both wrap the "
           "same engine routine in `macos_api_mod.F90`, so arguments, "
           "units and behavior match.  An entry with only one language "
           "is a coverage gap in the other (noted).",
           ""]
    n = 0
    titles = [t for t, _ in GROUPS] + ["Other utilities"]
    for title in titles:
        if title not in grouped:
            continue
        out.append(f"### {title}")
        out.append("")
        for name, e in sorted(grouped[title]):
            n += 1
            key = "fn-" + re.sub(r"[^A-Za-z0-9]+", "-", name).strip("-")
            out.append(f"#### {name}")
            out.append("")
            for mname, msig, mhelp in sorted(e["mm"]):
                # qualify the function name, not the output binding:
                # "v = get_x(s)" -> "v = macos.get_x(s)"
                qsig = re.sub(r"(^|=\s*)([A-Za-z]\w*\s*\()",
                              r"\1macos.\2", msig, count=1)
                out.append(f"- **mmacos:** `{qsig}`")
            if e["py"]:
                pname, psig, ppara = e["py"]
                alias = py_alias.get(pname)
                a = f"  (alias: {', '.join(alias)})" if alias else ""
                out.append(f"- **pymacos:** `pymacos.{psig}`{a}")
            if e["mm"] and not e["py"]:
                out.append("- **pymacos:** *not available*")
            if e["py"] and not e["mm"]:
                out.append("- **mmacos:** *not available*")
            api_hits = {norm(n2) for n2, _, _ in e["mm"]} | {norm(name)}
            hit = next((h for h in sorted(api_hits) if h in
                        {norm(a2) for a2 in api}), None)
            if hit:
                real = next(a2 for a2 in sorted(api) if norm(a2) == hit)
                out.append(f"- **engine API:** `macos_api_mod::{real}`")
            out.append("")
            desc = ""
            helps = [h for _, _, h in sorted(e["mm"]) if h]
            if helps:
                desc = " ".join(helps[0][:6])
                desc = re.sub(r"^MACOS\.\S+\s+", "", desc)
            elif e["py"] and e["py"][2]:
                desc = e["py"][2]
            if desc:
                out.append(desc)
                out.append("")
            out.append(notes_block(key, notes))
            if key.lower() not in notes and not desc:
                out.append(TODO)
            out.append("")
    (CMDREF / "20_bindings.md").write_text("\n".join(out))
    return n


# ------------------------------------------------------ higher level part
def gen_higher(notes):
    out = ["<!-- GENERATED by tools/gen_cmdref.py from mmacos +design and"
           " +channels.",
           "     Edit prose only inside BEGIN/END NOTES blocks;"
           " everything else is overwritten. -->",
           "",
           "## Part III — Higher-level MATLAB layers (current state)",
           "",
           "These packages sit on top of the Part II bindings.  They are "
           "under active development; interfaces here reflect the "
           "current tree and may evolve faster than Parts I–II.",
           ""]
    n = 0
    for sub, blurb in [
            ("+design",
             "Parametric telescope design layer: build a System from "
             "first-order parameters, add mirrors, solve conics/Zernikes "
             "against multi-field wavefront targets, and emit a MACOS "
             ".in prescription."),
            ("+channels",
             "Sensitivity-channel classes used by the dw_* Jacobian "
             "drivers: each channel enumerates one family of "
             "perturbations (rigid-body, Zernike, grid, source, ...) "
             "and knows how to apply/revert them.")]:
        d = MMACOS / sub
        if not d.exists():
            continue
        out.append(f"### macos.{sub[1:]} — {sub}")
        out.append("")
        out.append(blurb)
        out.append("")
        for name, (sig, helplines) in parse_mmacos_dir(d).items():
            n += 1
            key = "hl-" + re.sub(r"[^A-Za-z0-9]+", "-", name).strip("-")
            out.append(f"#### {sub[1:]}.{name}")
            out.append("")
            out.append(f"`{sig}`" if sig.startswith("classdef")
                       else f"`macos.{sub[1:]}.{sig}`")
            out.append("")
            if helplines:
                h = " ".join(helplines[:5])
                h = re.sub(r"^[A-Z.]*\S*\s+", "", h) if False else h
                out.append(h)
                out.append("")
            out.append(notes_block(key, notes))
            if key.lower() not in notes and not helplines:
                out.append(TODO)
            out.append("")
    (CMDREF / "30_higher_level.md").write_text("\n".join(out))
    return n


def main():
    CMDREF.mkdir(exist_ok=True)
    n1 = gen_engine(load_notes(CMDREF / "10_engine_commands.md"))
    n2 = gen_bindings(load_notes(CMDREF / "20_bindings.md"))
    n3 = gen_higher(load_notes(CMDREF / "30_higher_level.md"))
    print(f"engine commands: {n1}, binding functions: {n2}, "
          f"higher-level: {n3}")


if __name__ == "__main__":
    main()
