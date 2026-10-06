#!/usr/bin/env python3
"""CLI gate for SPOT centering after a LOAD (Luis, OPTIIX FSM 2e-7 rad, 2026-10-06).
Load Rx_Cass_FarField; TEXT (spot to a file); SPOT 3 TOUT; PERTURB elt 2 by 2e-7 rad
about x; SPOT again.  The MEAN of the written spot (RaySpot) must move with the chief
ray.  Pre-fix, elt_mod_init_vars set spcOption = 0 at every load, so the SPOT branch
took its ELSE (chief-ray centring): RaySpot did not move while the chief did.  With
--chfray the same binary is put in that state by SPCENTER chfray (must-fail leg)."""
import os, pty, select, sys, time, re, shutil, glob
chf = '--chfray' in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith('--')]
MACOS = os.path.expanduser(args[0] if args else '~/dev/macos/build_release_gfortran/bin/macos')
DECK = os.path.expanduser('~/dev/MACOS_resources/pymacos/tests/Rx/Rx_Cass_FarField')
work = os.path.join(os.getcwd(), 'spotgate'); shutil.rmtree(work, ignore_errors=True); os.makedirs(work); os.chdir(work)
pid, fd = pty.fork()
if pid == 0: os.execv(MACOS, [MACOS])
log = b''
mark = 0
def drain(timeout=25, until=None):
    global log, mark; t0 = time.time()
    while time.time() - t0 < timeout:
        r, _, _ = select.select([fd], [], [], 0.3)
        if r:
            try: chunk = os.read(fd, 65536)
            except OSError: break
            log += chunk
            if until and until in log[mark:]: mark = len(log); return True
        elif until is None: mark = len(log); return True
    mark = len(log); return until is None
def send(c): os.write(fd, c.encode() + b'\r')
def spot(tag):
    send('SPOT'); drain(10, b'spot diagram is to be computed'); send('6')
    drain(10, b'coordinate option'); send('TOUT')
    global mark
    t0 = time.time(); start = mark
    while time.time() - t0 < 150:
        drain(2)
        new = log[start:]
        if b'Graphics device/type' in new and b'/null' not in new: send('/null'); start = len(log)
        if new.rstrip().endswith(b'MACOS>'): break
    f = sorted(glob.glob(DECK + '.spot6*'), key=os.path.getmtime)   # SPOTOUT writes beside the deck
    assert f, 'no spot file written: ' + log.decode('utf-8','replace')[-1500:]
    shutil.move(f[-1], tag + '.spot'); return tag + '.spot'
drain(20, b'MACOS model size:'); send(''); drain(30, b'MACOS>')
send('OLD'); drain(10); send(DECK); drain(40, b'MACOS>')
send('TEXT'); drain(5, b'MACOS>')
if chf: send('SPCENTER'); drain(5, b'centering option'); send('chfray'); drain(5, b'MACOS>')
f0 = spot('nominal')
send('PERTURB'); drain(5); send('3'); drain(5); send('GLOBAL'); drain(5); send('2e-7 0 0'); drain(5); send('0 0 0'); drain(20, b'MACOS>')
f1 = spot('perturbed')
send('QUIT'); drain(5)
txt = log.decode('utf-8', 'replace')
rows = re.findall(r'Chief ray location: x=\s*([-\d.D+E]+)\s*y=\s*([-\d.D+E]+)', txt)
f = lambda s: float(s.replace('D', 'E'))
dchief = f(rows[-1][1]) - f(rows[0][1])
def mean_y(fn):
    ys = []
    for line in open(fn):
        p = line.replace('D','E').split()
        try: ys.append(float(p[1]))
        except (IndexError, ValueError): pass
    return sum(ys)/len(ys), len(ys)
m0, n0 = mean_y(f0); m1, n1 = mean_y(f1)
print('SPCENTER %s: chief moved %.3e; written spot mean moved %.3e (%d rays)' % ('chfray' if chf else 'elt (default after load)', dchief, m1 - m0, n1))
ok = abs(dchief) > 1e-9 and abs((m1 - m0) - dchief) < 0.05*abs(dchief)
print('spot follows the chief:', ok)
sys.exit(0 if ok else 1)
