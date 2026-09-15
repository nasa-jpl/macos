"""Panel crops for deck_gauges: cut a fractional box out of a lane's own
layout PNG, then trim white margins.  Cropping only, never re-rendering
(Dave's rule: deck figures are the tools' own output).
Run: python3 crop_panels.py"""
from PIL import Image, ImageChops
R = '/home/dcr/dev/MACOS_resources/mmacos/templates/40_benches/'
JOBS = [  # (source, out name, (l, t, r, b) fractions)
    (R+'tg_psi_dm96_oap/lens_vlayout.png', 'crop_lens_vlayout_train.png', (0.0, 0.0, 1.0, 0.33)),
    (R+'tg_psi_dm96_oap/lens_vlayout.png', 'crop_lens_vlayout_node.png',  (0.2, 0.33, 0.8, 0.66)),
    (R+'tg_psi_dm96_oap/lens_vlayout.png', 'crop_lens_vlayout_tail.png',  (0.2, 0.66, 0.8, 1.0)),
    (R+'tg_psi_dm96_oap/oap_vlayout.png',  'crop_oap_vlayout_train.png',  (0.0, 0.0, 1.0, 0.33)),
    (R+'tg_psi_dm96_oap/oap_vlayout.png',  'crop_oap_vlayout_tail.png',   (0.2, 0.66, 0.8, 1.0)),
    (R+'pdi_dm96/psri_layout.png',         'crop_psri_layout_train.png',  (0.0, 0.0, 1.0, 0.44)),
    (R+'pdi_dm96/psri_layout.png',         'crop_psri_layout_arm.png',    (0.0, 0.44, 1.0, 0.88)),
    (R+'pdi_dm96/pdi_layout.png',          'crop_pdi_layout_train.png',   (0.0, 0.0, 1.0, 0.37)),
    (R+'pdi_dm96/pdi_layout.png',          'crop_pdi_layout_tail.png',    (0.0, 0.37, 1.0, 1.0)),
    (R+'zwfs_dm96/zwfs_vlayout.png',       'crop_zwfs_vlayout_train.png', (0.0, 0.0, 1.0, 0.5)),
    (R+'zwfs_dm96/bench_bs30.png',         'crop_bench_bs30_train.png',   (0.15, 0.0, 0.85, 0.40)),
    (R+'zwfs_dm96/bench_bs30.png',         'crop_bench_bs30_node.png',    (0.15, 0.40, 0.85, 1.0)),
    (R+'zwfs_dm96/bench_bs22.png',         'crop_bench_bs22_node.png',    (0.15, 0.40, 0.85, 1.0)),
    (R+'zwfs_dm96/bench_bs7.png',          'crop_bench_bs7_node.png',     (0.15, 0.40, 0.85, 1.0)),
]
def trim(im):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    diff = ImageChops.difference(im.convert('RGB'), bg.convert('RGB'))
    bbox = diff.getbbox()
    if bbox:
        l, t, r, b = bbox; m = 12
        im = im.crop((max(0, l-m), max(0, t-m), min(im.width, r+m), min(im.height, b+m)))
    return im
for src, out, (l, t, r, b) in JOBS:
    im = Image.open(src).convert('RGB')
    W, H = im.size
    box = (int(l*W), int(t*H), int(r*W), int(b*H))
    im2 = trim(im.crop(box))
    im2.save('figs/' + out)
    print(out, im2.size)
