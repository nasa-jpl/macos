"""Panel crops for deck_gauges: cut a fractional box out of a lane's own
layout PNG, then trim white margins.  Cropping only, never re-rendering
(Dave's rule: deck figures are the tools' own output).
Run: python3 crop_panels.py

SOURCES ARE RUN OUTPUTS, so a retired run tag silently freezes a deck figure:
the OAP panels were still being cut from runs/oapdraw3, a render from before
the substrates were decided, and the lens/oap entries named bare *_vlayout.png
files no run writes any more.  Both now point at runs/lay96_{lens,oap}, whose
producer is tg96_run's figs stage (regenerate with runs/layoutfigs.sh).  If a
source path here does not exist, the crop is SKIPPED and the deck keeps the
stale PNG -- which is why this file now reports every miss at the end."""
from PIL import Image, ImageChops
R = '/home/dcr/dev/MACOS_resources/mmacos/templates/40_benches/'
JOBS = [  # (source, out name, (l, t, r, b) fractions)
    (R+'tg_psi_dm96_oap/runs/lay96_lens/lay96_lens_vlayout.png', 'crop_lens_vlayout_train.png', (0.0, 0.0, 1.0, 0.33)),
    (R+'tg_psi_dm96_oap/runs/lay96_lens/lay96_lens_vlayout.png', 'crop_lens_vlayout_node.png',  (0.2, 0.33, 0.8, 0.66)),
    (R+'tg_psi_dm96_oap/runs/lay96_lens/lay96_lens_vlayout.png', 'crop_lens_vlayout_tail.png',  (0.2, 0.66, 0.8, 1.0)),
    (R+'tg_psi_dm96_oap/runs/lay96_oap/lay96_oap_vlayout.png',  'crop_oap_vlayout_train.png',  (0.0, 0.0, 1.0, 0.33)),
    (R+'tg_psi_dm96_oap/runs/lay96_oap/lay96_oap_vlayout.png',  'crop_oap_vlayout_tail.png',   (0.2, 0.66, 0.8, 1.0)),
    (R+'pdi_dm96/psri_layout.png',         'crop_psri_layout_train.png',  (0.0, 0.0, 1.0, 0.44)),
    (R+'pdi_dm96/psri_layout.png',         'crop_psri_layout_arm.png',    (0.0, 0.44, 1.0, 0.88)),
    (R+'pdi_dm96/pdi_layout.png',          'crop_pdi_layout_train.png',   (0.0, 0.0, 1.0, 0.37)),
    (R+'pdi_dm96/pdi_layout.png',          'crop_pdi_layout_tail.png',    (0.0, 0.37, 1.0, 1.0)),
    (R+'zwfs_dm96/zwfs_vlayout.png',       'crop_zwfs_vlayout_train.png', (0.0, 0.0, 1.0, 0.5)),
    (R+'zwfs_dm96/zwfs_vlayout.png',       'crop_zwfs_vlayout_tail.png',  (0.0, 0.5, 1.0, 1.0)),
    (R+'zwfs_dm96/bench_bs30.png',         'crop_bench_bs30_train.png',   (0.15, 0.0, 0.85, 0.40)),
    (R+'zwfs_dm96/bench_bs30.png',         'crop_bench_bs30_node.png',    (0.15, 0.40, 0.85, 1.0)),
    (R+'zwfs_dm96/bench_bs22.png',         'crop_bench_bs22_node.png',    (0.15, 0.40, 0.85, 1.0)),
    (R+'zwfs_dm96/bench_bs7.png',          'crop_bench_bs7_node.png',     (0.15, 0.40, 0.85, 1.0)),
    (R+'zwfs_dm96/bench_bs22.png',         'crop_bench_bs22_train.png',   (0.15, 0.0, 0.85, 0.40)),
    (R+'tg_psi_dm96_oap/runs/lay96_oap/lay96_oap_vlayout.png', 'crop_lay96_oap_train.png', (0.15, 0.05, 0.9, 0.33)),
    (R+'tg_psi_dm96_oap/runs/lay96_oap/lay96_oap_vlayout.png', 'crop_lay96_oap_node.png',  (0.3, 0.33, 0.8, 0.66)),
]
MISSING = []


def trim(im):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    diff = ImageChops.difference(im.convert('RGB'), bg.convert('RGB'))
    bbox = diff.getbbox()
    if bbox:
        l, t, r, b = bbox; m = 12
        im = im.crop((max(0, l-m), max(0, t-m), min(im.width, r+m), min(im.height, b+m)))
    return im
import os
for src, out, (l, t, r, b) in JOBS:
    if not os.path.exists(src):
        # A MISSING SOURCE IS NOT A NO-OP: the old crop stays in figs/ and the
        # deck keeps showing it, which is how the OAP panels stayed two days
        # behind the substrate decision.  Say so, loudly, at the end.
        MISSING.append((out, src))
        print('SKIPPED (source missing):', out)
        continue
    im = Image.open(src).convert('RGB')
    W, H = im.size
    box = (int(l*W), int(t*H), int(r*W), int(b*H))
    im2 = trim(im.crop(box))
    im2.save('figs/' + out)
    print(out, im2.size)

# bench_three_nodes.png: the splitter-angle comparison, 7 / 22.5 / 30 left to
# right.  Tiled from the three node crops above rather than rendered, so it
# cannot drift from them.  It had no producer at all until 2026-09-18 -- it was
# assembled by hand and went stale with the rest.
TILE = ['crop_bench_bs7_node.png', 'crop_bench_bs22_node.png', 'crop_bench_bs30_node.png']
if not MISSING:
    ims = [Image.open('figs/' + f) for f in TILE]
    GAP = 52
    W = sum(i.width for i in ims) + GAP*(len(ims)-1)
    H = max(i.height for i in ims)
    comp = Image.new('RGB', (W, H), (255, 255, 255))
    x = 0
    for i in ims:
        comp.paste(i, (x, (H - i.height)//2));  x += i.width + GAP
    comp.save('figs/bench_three_nodes.png')
    print('bench_three_nodes.png', comp.size, '(tiled from the three node crops)')

if MISSING:
    print('\n%d CROP(S) NOT REBUILT -- the deck still shows the OLD figure:' % len(MISSING))
    for out, src in MISSING:
        print('   %-34s wanted %s' % (out, src))
    raise SystemExit(1)
print('\nall %d crops rebuilt' % len(JOBS))
