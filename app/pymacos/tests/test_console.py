"""MACOS_CONSOLE: the engine's console output inside a host (emsg_mod, read at the api's init).

Every WRITE(*,...) in the engine writes to the process's console handle.  A host that closes
or redirects that handle while the engine is loaded leaves the Fortran runtime with a write
that fails, and the runtime ends the process (forrtl severe 38; seen at random under pytest's
output capture on Windows, L. Marchen 2026-11-20).  MACOS_CONSOLE=off (or a file name) connects
the engine's unit 6 elsewhere, so the engine never writes to the console handle.

Each leg is a fresh Python process, because the setting is read once per process and the legs
take the console away.  Legs:
  negative control   console on, the handle closed under the engine: an ifx engine dies.
                     gfortran's runtime ignores the failed write, so a gfortran engine survives;
                     the leg asserts whichever applies to the library under test.
  off                the same with MACOS_CONSOLE=off: the process survives and traces, whether
                     the handle was closed after the engine's init or before it (closed before,
                     a gfortran engine faulted until unit 6 was closed ahead of its reopening)
  file               MACOS_CONSOLE=<file>: the engine's lines are in the file, not on stdout
  bad file           a file that cannot be opened: a warning, the console as before, no fort.6
  messages           with the console off the message buffer still holds the load's lines
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SRC = os.environ.get('PYMACOS_SRC', str(HERE.parent / 'src'))      # the side-tree .so too
RX = str(HERE / 'Rx' / 'Rx_Cass_FarField.in')
RX_MSG = str(HERE / 'Rx' / 'Rx_HdrUnknown.in')

CODE = (
    "import os, sys; sys.path.insert(0, %r)\n"
    "import pymacos.macos as m\n"
    "if %r == 'before': os.close(1)\n"               # the host has no console when the engine starts
    "m.init(128)\n"
    "if %r in (True, 'after'): os.close(1)\n"        # the console handle goes away under the engine
    "n = m.load(%r)\n"
    "m.trace_rays(-1)\n"
    "kept = m.messages() if %r else []\n"
    "sys.stderr.write('SURVIVED nElt %%d kept %%d\\n' %% (n, len(kept)))\n"
)


def run(console, close_fd, rx=RX, want_messages=False, cwd=None):
    """close_fd: False, True or 'after' (closed after the engine's init), or 'before' (closed before it)."""
    env = dict(os.environ)
    env.pop('MACOS_CONSOLE', None)
    if console is not None:
        env['MACOS_CONSOLE'] = console
    code = CODE % (SRC, close_fd, close_fd, rx, want_messages)
    return subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=300, env=env, cwd=cwd)


def fortran_runtime():
    """'ifx' or 'gfortran': which Fortran runtime the library under test is linked to (None if neither is named)."""
    import glob
    libs = glob.glob(os.path.join(SRC, 'pymacos', 'pymacosf90*'))
    blob = b''.join(open(f, 'rb').read() for f in libs if not f.endswith(('.py', '.pyf')))
    if b'libifcore' in blob:
        return 'ifx'
    return 'gfortran' if b'libgfortran' in blob else None


def console_on_or_skip():
    """The control: console on, nothing closed.  Its stdout must carry the load and the trace.
    On macOS the engine's unit 6 has reached no stream a test can read (test_messages guards its
    capfd leg for that); if that is so here, the legs that read stdout are skipped, not failed."""
    loud = run(None, close_fd=False)
    assert loud.returncode == 0, loud.stderr[-1500:]
    if len(loud.stdout) <= 200:
        if sys.platform == 'darwin':
            pytest.skip("the engine's console lines reach no stream the parent can read on macOS")
        raise AssertionError('the control prints the load and the trace: %r' % loud.stdout[:300])
    return loud


def test_negative_control_a_closed_console_kills_an_ifx_host():
    r = run(None, close_fd=True)
    rt = fortran_runtime()
    if rt == 'ifx':
        assert r.returncode != 0 and 'SURVIVED' not in r.stderr, (r.returncode, r.stderr[-500:])
    elif rt == 'gfortran':                           # its runtime ignores the failed write
        assert r.returncode == 0 and 'SURVIVED nElt' in r.stderr, (r.returncode, r.stderr[-1500:])
    else:
        pytest.skip('the Fortran runtime of the library could not be told')


@pytest.mark.parametrize('when', ['after', 'before'])
def test_with_the_console_off_the_host_survives_a_closed_console(when):
    r = run('off', close_fd=when)
    assert r.returncode == 0 and 'SURVIVED nElt' in r.stderr, (r.returncode, r.stderr[-1500:])


def test_a_file_that_cannot_be_opened_leaves_the_console_as_it_was(tmp_path):
    r = run(str(tmp_path / 'no_such_folder' / 'engine.log'), close_fd=False, cwd=str(tmp_path))
    assert r.returncode == 0 and 'SURVIVED nElt' in r.stderr, (r.returncode, r.stderr[-1500:])
    assert not (tmp_path / 'fort.6').exists(), 'the engine wrote its console lines to a new file fort.6'
    console_on_or_skip()
    assert 'could not be opened; console output is unchanged' in r.stdout, r.stdout[:300]
    assert len(r.stdout) > 200, 'the load and the trace still print'


def test_with_the_console_off_nothing_is_printed():
    quiet = run('off', close_fd=False)
    assert quiet.returncode == 0, quiet.stderr[-1500:]
    assert quiet.stdout.strip() == '', quiet.stdout[:300]
    console_on_or_skip()                             # ... and with the console on there is something to silence


def test_a_file_name_receives_the_engines_lines(tmp_path):
    log = tmp_path / 'engine.log'
    r = run(str(log), close_fd=False)
    assert r.returncode == 0 and r.stdout.strip() == '', (r.returncode, r.stdout[:300])
    text = log.read_text(errors='replace')
    assert 'Tracing' in text or 'rays' in text.lower(), text[:300]
    again = run(str(log), close_fd=False)                       # a second session appends
    assert again.returncode == 0 and len(log.read_text(errors='replace')) > len(text)


def test_the_message_buffer_still_holds_the_lines_with_the_console_off():
    r = run('off', close_fd=True, rx=RX_MSG, want_messages=True)
    assert r.returncode == 0, r.stderr[-1500:]
    kept = int(r.stderr.split('kept')[-1].split()[0])
    assert kept > 0, r.stderr[-300:]


@pytest.mark.parametrize('value', ['on', 'ON'])
def test_on_is_the_console(value):
    loud = console_on_or_skip()
    r = run(value, close_fd=False)
    assert r.returncode == 0 and len(r.stdout) > 200 and abs(len(r.stdout) - len(loud.stdout)) < 100, (r.returncode, len(r.stdout), len(loud.stdout))
