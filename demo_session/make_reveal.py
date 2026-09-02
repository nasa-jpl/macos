#!/usr/bin/env python3
"""Live-reveal swap for slide 33 ("The live design, revealed").

Harvests the just-finished oi_demo_step(<W>) run bundle, crops its
figures, splices the live numbers into a COPY of the deck
(deck_keysight_live.md -> deck_keysight_live.pptx), and renders the
replaced slide for a pre-flight check.  The presented deck is never
touched; Dave opens deck_keysight_live.pptx on prompt at the reveal.

Run (from demo_session/):  python3 make_reveal.py <width_deg>
Dry run uses the committed 12 deg rehearsal bundle: make_reveal.py 12
"""
import re, sys, shutil, subprocess, time, os, glob
from PIL import Image

t0 = time.time()
W = int(sys.argv[1])
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
BUNDLE = os.path.expanduser(
    sys.argv[2] if len(sys.argv) > 2 else
    "~/dev/MACOS_res_dev/mmacos/templates/10_telescopes/offset_imager/"
    "demo_adjacent")
SLIDE_TITLE = "## The live design, revealed |"

# ---- 1. harvest the NEWEST matching bundle -----------------------
# live runs are timestamped (oi_demo_<W>deg_<stamp>_*); the curated
# rehearsal bundles are bare (oi_demo_<W>deg_*).  Newest verdict wins.
cands = glob.glob(f"{BUNDLE}/oi_demo_{W}deg*_verdict.txt")
assert cands, f"no oi_demo_{W}deg*_verdict.txt in {BUNDLE}"
vfile = max(cands, key=os.path.getmtime)
STEM = vfile[:-len("_verdict.txt")]
print(f"  bundle: {os.path.basename(STEM)}  "
      f"(of {len(cands)} candidate{'s' if len(cands)>1 else ''})")
vtxt = open(vfile).read()
def grab(pat, what):
    m = re.search(pat, vtxt)
    assert m, f"verdict parse failed: {what}"
    return m.groups()
(status,)        = grab(r"\[(PASS|FAIL)\]", "PASS/FAIL")
(solved,)        = grab(r"DENSE MAP MAX\s*:\s*([\d.]+) nm", "solved")
pred, meas, ratio = grab(r"predicted ~([\d.]+) nm -> measured ([\d.]+) nm"
                         r" \(([\d.]+)x\)", "prediction")
(clear,)         = grab(r"clearance floor\s*:\s*([\d.]+) mm", "clearance")
(wall,)          = grab(r"wall time\s*:\s*([\d.]+) min", "wall time")
assert meas == solved, "verdict internal mismatch (map max vs measured)"
print(f"  verdict: [{status}] solved {solved} nm vs ~{pred} nm predicted "
      f"({ratio}x), clearance {clear} mm, {wall} min")

# ---- 2. crop the live figures ------------------------------------
PAD, THRESH, GUT = 14, 248, 26
def trim(im, pad=PAD, thresh=THRESH):
    mask = im.convert("L").point(lambda v: 255 if v < thresh else 0)
    l, t, r, b = mask.getbbox() or (0, 0, im.width, im.height)
    return im.crop((max(0, l - pad), max(0, t - pad),
                    min(im.width, r + pad), min(im.height, b + pad)))

lay = Image.open(f"{STEM}_layout.png").convert("RGB")
px = lay.load()                       # whiten the spanning suptitle band
for y in range(min(80, lay.height)):
    for x in range(int(lay.width * .3), int(lay.width * .7)):
        px[x, y] = (255, 255, 255)
hw, hh = lay.width // 2, lay.height // 2
q = [trim(lay.crop(b)) for b in [(0, 0, hw, hh), (hw, 0, lay.width, hh),
     (0, hh, hw, lay.height), (hw, hh, lay.width, lay.height)]]
cw = [max(q[0].width, q[2].width), max(q[1].width, q[3].width)]
rh = [max(q[0].height, q[1].height), max(q[2].height, q[3].height)]
out = Image.new("RGB", (sum(cw) + GUT, sum(rh) + GUT), (255, 255, 255))
for i, qi in enumerate(q):
    out.paste(qi, ((i % 2) * (cw[0] + GUT), (i // 2) * (rh[0] + GUT)))
out.save("figs/oi_demo_live_layout.png")
trim(Image.open(f"{STEM}_map.png").convert("RGB"),
     pad=12).save("figs/oi_demo_live_map.png")
print(f"  figs: oi_demo_live_layout.png {out.size} + oi_demo_live_map.png")

# ---- 3. splice slide 33 into a COPY of the deck ------------------
md = open("deck_keysight.md").read()
i0 = md.index(SLIDE_TITLE)                       # keep title line verbatim
i1 = md.index("\n", i0) + 1                      # body starts after title
i2 = md.index("\n## ", i1)                       # next slide
body = f"""::: left
- **The ask from the room:** a {W}×{W}° field box.  **The prediction, stated before the solve:** ~{pred} nm, interpolated from the committed frontier rows.
- **The live run [{status}]:** dense-map max **{solved} nm** — {ratio}× the prediction — at a {clear} mm clearance floor, in {wall} minutes, while we talked.
- The run is deterministic to every printed digit, and every ask is a genuine continuation step — never a re-score of a stored design.
::: right
![The room's {W}° design and its field map: {solved} nm max against a ~{pred} nm prediction.](figs/oi_demo_live_layout.png){{h=3.0}}
![](figs/oi_demo_live_map.png){{h=2.2}}
~ Solved live during this talk; artifacts {os.path.basename(STEM)}_* in the demo_adjacent directory.
"""
open("deck_keysight_live.md", "w").write(md[:i1] + body + md[i2:])

# geo overrides are keyed by FILENAME (img:) / text prefix (txt:) --
# remap the reveal slide's keys to the live blocks so the presented
# deck's slide-33 geometry (map under the text, layout large on the
# right) carries over verbatim.
import json
geo = json.load(open("deck_keysight.geo.json"))
sl = geo["The live design, revealed"]
remap = {"img:crop_oi_demo_12deg_layout.png": "img:oi_demo_live_layout.png",
         "img:oi_demo_12deg_map.png":         "img:oi_demo_live_map.png",
         "txt:On demo day the live run's figures":
                                              "txt:Solved live during this talk"}
for old, new in remap.items():
    assert old in sl, f"geo key drift: {old} missing on the reveal slide"
    sl[new] = sl.pop(old)
json.dump(geo, open("deck_keysight_live.geo.json", "w"), indent=1)

# ---- 4. build + render the check ---------------------------------
subprocess.run(["python3", "make_brief_slides.py", "deck_keysight_live.md"],
               check=True, capture_output=True)
subprocess.run(["soffice", "--headless", "--convert-to", "pdf",
                "deck_keysight_live.pptx", "--outdir", "renders"],
               check=True, capture_output=True)
subprocess.run(["pdftoppm", "-f", "33", "-l", "33", "-r", "90", "-png",
                "renders/deck_keysight_live.pdf", "renders/live33"],
               check=True, capture_output=True)
print(f"  built deck_keysight_live.pptx + renders/live33-33.png")
print(f"READY in {time.time()-t0:.0f} s -- prompt Dave to open "
      f"demo_session/deck_keysight_live.pptx at the reveal")
