#!/usr/bin/env python3
"""Where can a DM surface gauge look at a coronagraph's DM from?

Reads a MACOS deck (the CTB bench, the e2e6m space relay), takes every real
element as a body (sphere of its stated radius + mount) and every leg between
consecutive real elements as a science beam (cylinder), then for each DM scans
a gauge beam about the DM normal: polar angle theta (0 = normal incidence),
azimuth phi about the normal.  The gauge's IN and OUT legs sit at +theta and
-theta in the plane of azimuth phi.  Two checks, both numbers:

  beam   -- min clearance of the two gauge legs (radius = the DM beam) from
            every body except the DM itself, over a length L (mm);
  optic  -- the gauge's first optic, a body of radius r_gauge + mount placed
            at distance L_h along each leg, against every body AND every
            science beam (light crossing light is not a collision; hardware
            in a beam is).

Prints, per DM: the normal-incidence verdict, the best azimuth at each
theta, and the smallest theta that clears by >= MARGIN.  No MATLAB.
Run: python3 gauge_in_coro_clearance.py
"""
import re, sys, math
import numpy as np

R = '/home/dcr/dev/MACOS_resources/mmacos/templates/'
MOUNT = 8.0      # mm, mount allowance on every body (dmg_bench_clearance's)
MARGIN = 25.0    # mm, the bench rule
THETAS = [0, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30]
PHIS = np.arange(0, 360, 5)
L_BEAM = 700.0   # mm, how far the gauge legs are checked
L_OPT = [300.0, 500.0]  # mm, candidate stations for the gauge's first optic

CASES = [
    dict(name='CTB bench (ctb_dcr.in, mm)', deck=R+'30_instruments/bench_ctb/ctb_dcr.in', scale=1.0,
         radius={'OAP': 75.0, 'DM': 22.5, 'default': 25.0}, r_beam=22.5, r_sci=10.62, first=1),   # gauge covers the DM's clear radius; the traced science beam is 21.24 mm across (ctb_beam_probe, TO 2026-09-16)
    dict(name='e2e6m space relay (r1_seg_d040_full.in, m -> mm)', deck=R+'80_end_to_end/e2e6m_r2/r1_seg_d040_full.in',
         scale=1000.0, radius={'OAP': 60.0, 'DM': 30.0, 'M2': 400.0, 'M3': 400.0, 'default': 25.0}, r_beam=30.0, r_sci=23.75, first=20),
]

def parabola_hit(P, u, V, a, R):
    """point where the line P + s u meets the parabola with vertex V, axis a (unit), radius R: z = |h|^2/(2R)"""
    best = None
    for s0 in np.linspace(50, 6000, 120):
        s = s0
        for _ in range(60):
            X = P + s * u; d = X - V; z = d @ a; h = d - z * a
            F = z - (h @ h) / (2 * R)
            dF = u @ a - (h @ (u - (u @ a) * a)) / R
            if abs(dF) < 1e-12: break
            s -= F / dF
            if abs(F) < 1e-9: break
        if s > 1 and abs(F) < 1e-6 and (best is None or s < best): best = s
    return P + best * u if best is not None else None

def read_deck(path, scale, first):
    s = open(path).read()
    hdr = lambda k: np.array([float(x.replace('D', 'E')) for x in re.search(r'\b' + k + r'\s*=\s*([^\n]*)', s).group(1).split()])
    blocks = re.split(r'\n(?=\s*iElt\s*=)', s)
    els = []
    for b in blocks[1:]:
        g = lambda k: (re.search(r'\b' + k + r'\s*=\s*([^\n]*)', b) or [None, ''])[1].strip()
        ie = g('iElt')
        if not ie.isdigit() or int(ie) < first: continue
        typ = g('Element')
        if typ not in ('Reflector', 'Refractor', 'FocalPlane', 'Reference', 'Return'): continue
        if typ == 'Return' and g('Surface') != 'Flat': continue
        nm = g('EltName')
        v = np.array([float(x.replace('D', 'E')) for x in g('VptElt').split()]) * scale
        p = np.array([float(x.replace('D', 'E')) for x in g('psiElt').split()])
        if any(np.linalg.norm(v - e['v']) < 1e-3 for e in els): continue
        if nm.endswith('_start') or nm.endswith('_end') or nm.endswith('_Pst'): continue
        kr = float(g('KrElt').replace('D', 'E')) * scale if g('KrElt') else 1e22
        els.append(dict(i=int(ie), name=nm, type=typ, v=v.copy(), vtx=v.copy(), psi=p / np.linalg.norm(p), R=abs(kr),
                        oap=(typ == 'Reflector' and g('Surface') == 'Conic' and nm.upper().startswith('OAP'))))
    # stations: flats and markers are exact; an OAP's station is where the running chief meets its parabola,
    # and the chief reflects there off the parabola's own normal (grad of z - |h|^2/2R = a - h/R)
    def refl(u, n): return u - 2 * (u @ n) * n
    def para_n(X, e):
        d = X - e['vtx']; z = d @ e['psi']; h = d - z * e['psi']; n = e['psi'] - h / e['R']; return n / np.linalg.norm(n)
    k0 = [i for i, e in enumerate(els) if e['name'].upper().startswith('DM')][0]
    if els[0]['oap'] and k0 == 1:      # CTB: the source is in the header
        P = hdr('ChfRayPos') * scale; u = hdr('ChfRayDir'); u /= np.linalg.norm(u); kstart = 0
    else:                              # space relay: OAP1 by walking backward from DM1 (M2/M3 stay at their vertices)
        d_out = els[k0 + 1]['v'] - els[k0]['v']; d_out /= np.linalg.norm(d_out)
        d_in = refl(d_out, els[k0]['psi'])
        e = els[k0 - 1]; assert e['oap']
        e['v'] = parabola_hit(els[k0]['v'], -d_in, e['vtx'], e['psi'], e['R'])
        P = els[k0]['v']; u = d_out; kstart = k0 + 1
    for k in range(kstart, len(els)):
        e = els[k]
        if e['oap']:
            X = parabola_hit(P, u, e['vtx'], e['psi'], e['R'])
            if X is None: print('  ! no parabola hit for', e['name']); X = e['vtx']
            e['v'] = X; u = refl(u, para_n(X, e)); P = X
        elif e['type'] == 'Reflector':
            P = e['v']; u = refl(u, e['psi'])
        else:
            P = e['v']
    for e in els:
        if e['oap']: print(f"  {e['name']}: pole {np.round(e['v'],1)}  (parent vertex {np.round(e['vtx'],1)}, off-axis {np.linalg.norm(e['v']-e['vtx']):.1f} mm)")
    return els

def body_radius(e, radius):
    for k, r in radius.items():
        if k != 'default' and e['name'].upper().startswith(k.upper()): return r
    return radius['default']

def seg_seg_dist(p1, q1, p2, q2):
    """closest distance between segments p1-q1 and p2-q2 (Ericson 5.1.9)"""
    d1 = q1 - p1; d2 = q2 - p2; r = p1 - p2
    a = d1 @ d1; e = d2 @ d2; f = d2 @ r
    if a < 1e-12 and e < 1e-12: return np.linalg.norm(r)
    if a < 1e-12: s = 0.0; t = min(max(f / e, 0.0), 1.0)
    else:
        c = d1 @ r
        if e < 1e-12: t = 0.0; s = min(max(-c / a, 0.0), 1.0)
        else:
            b = d1 @ d2; den = a * e - b * b
            s = min(max((b * f - c * e) / den, 0.0), 1.0) if den > 1e-12 else 0.0
            t = (b * s + f) / e
            if t < 0: t = 0.0; s = min(max(-c / a, 0.0), 1.0)
            elif t > 1: t = 1.0; s = min(max((b - c) / a, 0.0), 1.0)
    return np.linalg.norm((p1 + d1 * s) - (p2 + d2 * t))

def point_seg_dist(x, p, q):
    d = q - p; t = min(max(((x - p) @ d) / (d @ d), 0.0), 1.0)
    return np.linalg.norm(x - (p + t * d))

def scan(case):
    els = read_deck(case['deck'], case['scale'], case['first'])
    rad = case['radius']; rg = case['r_beam']
    legs = [(els[k], els[k + 1]) for k in range(len(els) - 1)]
    dms = [e for e in els if e['name'].upper().startswith('DM')]
    print('=' * 100); print(case['name']); print('elements:', ', '.join(f"{e['name']}" for e in els))
    for a, b in legs:
        if a['name'].startswith('DM') or b['name'].startswith('DM'):
            print(f"  leg {a['name']:>10} -> {b['name']:<10} {np.linalg.norm(b['v']-a['v']):7.1f} mm")
    for dm in dms:
        n = dm['psi']
        # incidence of the science beam on this DM
        k = [i for i, e in enumerate(els) if e is dm][0]
        din = dm['v'] - els[k - 1]['v']; din /= np.linalg.norm(din)
        aoi = math.degrees(math.acos(abs(din @ n)))
        # frame about the normal
        t1 = np.cross(n, [0, 0, 1.0]);
        if np.linalg.norm(t1) < 0.3: t1 = np.cross(n, [0, 1.0, 0])
        t1 /= np.linalg.norm(t1); t2 = np.cross(n, t1)
        print(f"\n-- {dm['name']}: science AOI {aoi:.1f} deg; gauge beam radius {rg} mm; other bodies + {MOUNT} mm mounts; rule >= {MARGIN} mm")
        print(f"   {'theta':>5} | {'best phi':>8} {'beam clr':>9} {'binds':>14} | optic at {L_OPT[0]:.0f}: {'clr':>7} {'binds':>14} | at {L_OPT[1]:.0f}: {'clr':>7} {'binds':>14}")
        first_ok = None
        for th in THETAS:
            best = None
            for ph in (PHIS if th > 0 else [0]):
                u = math.cos(math.radians(ph)) * t1 + math.sin(math.radians(ph)) * t2
                dirs = [math.cos(math.radians(th)) * n + math.sin(math.radians(th)) * u,
                        math.cos(math.radians(th)) * n - math.sin(math.radians(th)) * u]
                # beam vs bodies
                cb, wb = 1e9, ''
                for d in dirs:
                    p, q = dm['v'], dm['v'] + L_BEAM * d
                    for e in els:
                        if e is dm: continue
                        c = point_seg_dist(e['v'], p, q) - body_radius(e, rad) - rg - MOUNT
                        if c < cb: cb, wb = c, e['name']
                # optic vs bodies + science beams
                co = []
                for Lh in L_OPT:
                    cc, ww = 1e9, ''
                    for d in dirs:
                        x = dm['v'] + Lh * d
                        for e in els:
                            if e is dm: continue
                            c = np.linalg.norm(x - e['v']) - body_radius(e, rad) - rg - MOUNT
                            if c < cc: cc, ww = c, e['name']
                        for a, b in legs:
                            c = point_seg_dist(x, a['v'], b['v']) - case['r_sci'] - rg - MOUNT
                            if c < cc: cc, ww = c, f"beam {a['name']}>{b['name']}"
                    co.append((cc, ww))
                score = min(cb, co[0][0])
                if best is None or score > best[0]: best = (score, ph, cb, wb, co)
            score, ph, cb, wb, co = best
            print(f"   {th:5.0f} | {ph:8.0f} {cb:9.1f} {wb:>14} | {co[0][0]:7.1f} {co[0][1]:>22} | {co[1][0]:7.1f} {co[1][1]:>22}")
            if first_ok is None and cb >= MARGIN and co[0][0] >= MARGIN: first_ok = th
        print(f"   smallest theta clearing both by >= {MARGIN} mm with the optic at {L_OPT[0]:.0f} mm: {first_ok}")

if __name__ == '__main__':
    for c in CASES: scan(c)
