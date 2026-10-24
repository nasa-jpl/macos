#!/usr/bin/env python3
"""White-margin autocrop for deck figures (DECK_STYLE: minimize whitespace).

Reads figs/<name>.png (symlinks to committed artifacts), writes
figs/crop_<name>.png.  Cropping only — never regenerates content.
Run: python3 crop_figs.py   (PIL only, system python)
"""
from PIL import Image

TARGETS = [
    "r3t_s1_layout.png",
    "r3t_s4_layout.png",
    "r3t_s4_map.png",
    "t5_walk_k05_layout.png",
    "t5_walk_k03_map.png",
]

PAD = 12
THRESH = 248


def autocrop(name):
    im = Image.open(f"figs/{name}").convert("RGB")
    gray = im.convert("L")
    mask = gray.point(lambda v: 255 if v < THRESH else 0)
    bbox = mask.getbbox()
    if bbox is None:
        print(f"{name}: blank? skipped")
        return
    l, t, r, b = bbox
    l, t = max(0, l - PAD), max(0, t - PAD)
    r, b = min(im.width, r + PAD), min(im.height, b + PAD)
    out = im.crop((l, t, r, b))
    out.save(f"figs/crop_{name}")
    print(f"{name}: {im.size} -> {out.size} "
          f"({100 - 100 * out.width * out.height // (im.width * im.height)}% removed)")


for n in TARGETS:
    autocrop(n)


# ---- 2x2 panel recompose (DECK_STYLE): blank the suptitle, split into
# quadrants, white-trim each, reassemble with narrow gutters. ----------
def recompose(name, title_frac=0.055, pad=10, gutter=24, bbox_thresh=200):
    # bbox_thresh < THRESH: faint gray label-leader lines are ignored when
    # sizing each panel, so a long leader cannot hold a tall empty margin.
    im = Image.open(f"figs/{name}").convert("RGB")
    im = im.crop((0, int(im.height * title_frac), im.width, im.height))
    hw, hh = im.width // 2, im.height // 2
    quads = [im.crop(b) for b in
             [(0, 0, hw, hh), (hw, 0, im.width, hh),
              (0, hh, hw, im.height), (hw, hh, im.width, im.height)]]
    trimmed = []
    for q in quads:
        mask = q.convert("L").point(lambda v: 255 if v < bbox_thresh else 0)
        bb = mask.getbbox() or (0, 0, q.width, q.height)
        l, t, r, b = bb
        trimmed.append(q.crop((max(0, l - pad), max(0, t - pad),
                               min(q.width, r + pad), min(q.height, b + pad))))
    cw = [max(trimmed[0].width, trimmed[2].width),
          max(trimmed[1].width, trimmed[3].width)]
    rh = [max(trimmed[0].height, trimmed[1].height),
          max(trimmed[2].height, trimmed[3].height)]
    out = Image.new("RGB", (cw[0] + cw[1] + gutter, rh[0] + rh[1] + gutter),
                    (255, 255, 255))
    pos = [(0, 0), (cw[0] + gutter, 0),
           (0, rh[0] + gutter), (cw[0] + gutter, rh[0] + gutter)]
    for q, (x, y) in zip(trimmed, pos):
        out.paste(q, (x, y))
    out.save(f"figs/recomp_{name}")
    print(f"recomp_{name}: {Image.open(f'figs/{name}').size} -> {out.size}")


for n in ["r3t_s4_layout.png", "t5_walk_k05_layout.png"]:
    recompose(n)


# ---- bottom-row pair (iso + side views): the clearance/packaging story,
# as a wide 2-panel strip that fills a slide column. ----------------------
def bottom_pair(name, title_frac=0.055, pad=10, gutter=24, bbox_thresh=200):
    im = Image.open(f"figs/{name}").convert("RGB")
    im = im.crop((0, int(im.height * title_frac), im.width, im.height))
    hw, hh = im.width // 2, im.height // 2
    quads = [im.crop((0, hh, hw, im.height)),
             im.crop((hw, hh, im.width, im.height))]
    trimmed = []
    for q in quads:
        mask = q.convert("L").point(lambda v: 255 if v < bbox_thresh else 0)
        bb = mask.getbbox() or (0, 0, q.width, q.height)
        l, t, r, b = bb
        trimmed.append(q.crop((max(0, l - pad), max(0, t - pad),
                               min(q.width, r + pad), min(q.height, b + pad))))
    h = max(t.height for t in trimmed)
    out = Image.new("RGB", (sum(t.width for t in trimmed) + gutter, h),
                    (255, 255, 255))
    out.paste(trimmed[0], (0, (h - trimmed[0].height) // 2))
    out.paste(trimmed[1], (trimmed[0].width + gutter,
                           (h - trimmed[1].height) // 2))
    out.save(f"figs/pair_{name}")
    print(f"pair_{name}: -> {out.size}")


for n in ["r3t_s4_layout.png", "t5_walk_k05_layout.png",
          "r3t_s1_layout.png"]:
    bottom_pair(n)


# ---- 1x4 row strips (view_std single-row): keep two panels ------------
def row_keep(name, idx=(0, 2), title_frac=0.07, pad=10, gutter=24,
             bbox_thresh=200):
    im = Image.open(f"figs/{name}").convert("RGB")
    im = im.crop((0, int(im.height * title_frac), im.width, im.height))
    qw = im.width // 4
    kept = []
    for i in idx:
        q = im.crop((i * qw, 0, (i + 1) * qw, im.height))
        mask = q.convert("L").point(lambda v: 255 if v < bbox_thresh else 0)
        bb = mask.getbbox() or (0, 0, q.width, q.height)
        l, t, r, b = bb
        kept.append(q.crop((max(0, l - pad), max(0, t - pad),
                            min(q.width, r + pad), min(q.height, b + pad))))
    h = max(k.height for k in kept)
    out = Image.new("RGB", (sum(k.width for k in kept) + gutter, h),
                    (255, 255, 255))
    out.paste(kept[0], (0, (h - kept[0].height) // 2))
    out.paste(kept[1], (kept[0].width + gutter, (h - kept[1].height) // 2))
    out.save(f"figs/rowpair_{name}")
    print(f"rowpair_{name}: -> {out.size}")


for n in ["s1_views.png", "s3_views_pie.png"]:
    row_keep(n)
