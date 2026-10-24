#!/usr/bin/env python3
"""Generate the README text the manual keeps for a design template that has no README.md of its own:
part2/<name>.md, from the header comment of the template's entry script (the spec's entry, the first
.m it runs).  Nothing is typed: the title is the header's heading line(s), the body its first prose
paragraph(s), up to the first all-caps section label or rule.  The templates themselves are not
touched (whether to promote these into them is Dave's question, OPEN_QUESTIONS).

  python3 tools/gen_readme_part2.py name [name ...]     (from docs/mmacos-manual)
"""
import json, re, sys
from pathlib import Path

MAN = Path(__file__).resolve().parent.parent
SPECS = MAN / "tools" / "part2" / "specs"
KEPT = MAN / "part2"
TPL = Path.home() / "dev" / "MACOS_resources" / "mmacos" / "templates"

# words an all-caps heading keeps in capitals, and the proper names it keeps capitalized
KEEP_UPPER = {"RC", "TMA", "FOV", "M1", "M2", "M3", "M4", "FP", "FP2", "DM", "OAP", "PSI", "IFO", "MET",
              "MACOS", "OTA", "APLC", "VIS", "CTB", "HCIT", "VSG2", "J18", "PROPER"}
PROPER = {"RITCHEY": "Ritchey", "CHRETIEN": "Chretien", "KORSCH": "Korsch", "ZERNIKE": "Zernike",
          "CASSEGRAIN": "Cassegrain", "TWYMAN": "Twyman", "GREEN": "Green", "MACOS": "MACOS", "GROOT": "Groot", "JONES": "Jones"}


def heading_case(text):
    out = []
    for i, w in enumerate(re.split(r"(\s+|-|\+)", text)):
        u = re.sub(r"[^A-Z0-9]", "", w.upper())
        if u in KEEP_UPPER:
            out.append(w.upper())
        elif u in PROPER:
            out.append(w.replace(u, PROPER[u]) if u in w else PROPER[u])
        elif w.isupper():
            out.append(w.lower())
        else:
            out.append(w)
    s = "".join(out)
    return s[:1].upper() + s[1:]


def header_lines(script):
    lines = []
    for ln in script.read_text().splitlines():
        if not ln.startswith("%"):
            if lines:
                break
            continue
        lines.append(re.sub(r"^%\s?", "", ln).rstrip())
    return lines


def entry_script(spec, name):
    m = re.search(r"'([A-Za-z0-9_]+\.m)'", spec["entry"])
    if m:
        return m.group(1)
    return re.search(r"([A-Za-z0-9_]+)\(", spec["entry"]).group(1) + ".m"


def build(name):
    spec = json.loads((SPECS / f"{name}.json").read_text())
    script = TPL / spec["section"] / name / entry_script(spec, name)
    hl = header_lines(script)
    rule = lambda s: re.fullmatch(r"[=\-]{5,}", s.strip()) is not None
    h1 = re.match(r"^([A-Z][A-Z0-9_]+)\s+(\S.*)$", hl[0]) if hl else None
    if h1 and not any(rule(x) for x in hl[:4]):  # a MATLAB help header: %NAME  summary, then text
        title = h1.group(2).strip().rstrip(".")
        body = []
        for ln in hl[1:]:
            t_ = ln.strip()
            if not t_:
                if body:
                    break
                continue
            u_ = re.match(r"^\w+\s*=\s*([A-Z][A-Z0-9_]+)\(([^)]*)\)\s*(.*)$", t_)   # "out = NAME(...) runs ..."
            if u_:
                t_ = f"`{u_.group(1).lower()}({u_.group(2)})` {u_.group(3)}".strip()
            body.append(t_)
        para = re.sub(r"\s+", " ", " ".join(body)).strip()
        return write_(name, title, para, script)
    # the heading: the all-caps line(s) between the first two rules (or right after the file-name line)
    i = 0
    while i < len(hl) and not rule(hl[i]):
        i += 1
    i += 1
    caps = lambda t: re.sub(r"[^A-Za-z]", "", t) != "" and re.sub(r"\([^)]*\)", "", t).upper() == re.sub(r"\([^)]*\)", "", t)
    head = []
    while i < len(hl) and not rule(hl[i]) and hl[i].strip() and caps(hl[i].strip()):
        head.append(hl[i].strip())
        i += 1
    if not head and i < len(hl) and not rule(hl[i]):  # a mixed-case heading line ("SPHERE + ZERNIKE unobscured ...")
        head.append(hl[i].strip())
        i += 1
    while i < len(hl) and (rule(hl[i]) or not hl[i].strip()):
        i += 1
    while i < len(hl) and hl[i].strip().startswith("("):   # a parenthetical aside under the heading
        while i < len(hl) and hl[i].strip() and not rule(hl[i]):
            i += 1
        while i < len(hl) and (rule(hl[i]) or not hl[i].strip()):
            i += 1
    # a single subtitle line closed by a rule belongs to the title ("... -- vs 3-mirror (Korsch TMA)")
    sub = ""
    if i + 1 < len(hl) and hl[i].strip() and rule(hl[i + 1]):
        sub = hl[i].strip()
        i += 2
    title = " ".join(h for h in head if not h.startswith("("))
    title = re.sub(r"^MACOS DESIGN LAYER\s*--\s*", "", title)
    title = heading_case(title.strip('"'))
    if sub:
        title = title.rstrip(':') + ("" if title.endswith("?") else ":") + " " + sub
    # the body: prose lines up to the first all-caps label ("THE ...", "GEOMETRY ...") or rule
    label = r"^[A-Z][A-Z0-9 /+-]{4,}[.:]\s"
    body = []
    paras = 0
    for ln in hl[i:]:
        s = ln.strip()
        if not s:                                       # at most the first two paragraphs
            paras += 1 if body else 0
            if paras >= 2:
                break
            continue
        if not body and re.match(label, s):       # a body that opens with its own label ("ARCHITECTURE.")
            s = re.sub(label, "", s)
        elif rule(s) or re.match(label, s) or re.match(r"^[A-Z][A-Z0-9 /+-]{4,}\(", s) or s.startswith("Run:"):
            break
        body.append(re.sub(r"^\*\s+", "", s))            # a bullet's marker, not its text
    para = ""
    for s_ in body:                                     # rejoin, mending words hyphenated at a line break
        if re.search(r"[A-Za-z]-$", para) and s_[:1].islower():
            para += s_
        else:
            para += (" " if para else "") + s_
    para = re.sub(r"\s+", " ", para).strip()
    return write_(name, title, para, script)


def write_(name, title, para, script):
    assert title and para, f"{name}: no heading or no prose in {script.name}'s header"
    KEPT.mkdir(exist_ok=True)
    (KEPT / f"{name}.md").write_text(
        f"# {name} — {title}\n\n{para}\n\n## Source\n\nGenerated by tools/gen_readme_part2.py from the header of "
        f"`{script.name}`; the template has no README.md.\n")
    return title


if __name__ == "__main__":
    for n in sys.argv[1:]:
        print(f"{n}: {build(n)}")
