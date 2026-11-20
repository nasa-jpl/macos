"""MACOS_CONSOLE: the engine's console output inside a host (emsg_mod, read at the api's init).

Every WRITE(*,...) in the engine writes to the process's console handle.  A host that closes
or redirects that handle while the engine is loaded leaves the Fortran runtime with a write
that fails, and the runtime ends the process (forrtl severe 38; seen at random under pytest's
output capture on Windows, L. Marchen 2026-11-20).  MACOS_CONSOLE=off (or a file name) connects
the engine's unit 6 elsewhere, so the engine never writes to the console handle.

Each leg is a fresh Python process, because the setting is read once per process and the legs
take the console away.  Legs:
  negative control   console on, the handle closed under the engine: the process dies
  off                the same with MACOS_CONSOLE=off: the process survives and traces
  file               MACOS_CONSOLE=<file>: the engine's lines are in the file, not on stdout
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
    "m.init(128)\n"
    "if %r: os.close(1)\n"                           # the console handle goes away under the engine
    "n = m.load(%r)\n"
    "m.trace_rays(-1)\n"
    "kept = m.messages() if %r else []\n"
    "sys.stderr.write('SURVIVED nElt %%d kept %%d\\n' %% (n, len(kept)))\n"
)


def run(console, close_fd, rx=RX, want_messages=False):
    env = dict(os.environ)
    env.pop('MACOS_CONSOLE', None)
    if console is not None:
        env['MACOS_CONSOLE'] = console
    code = CODE % (SRC, close_fd, rx, want_messages)
    return subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=300, env=env)


def test_negative_control_a_closed_console_kills_the_host():
    r = run(None, close_fd=True)
    assert r.returncode != 0 and 'SURVIVED' not in r.stderr, (r.returncode, r.stderr[-500:])


def test_with_the_console_off_the_host_survives_a_closed_console():
    r = run('off', close_fd=True)
    assert r.returncode == 0 and 'SURVIVED nElt' in r.stderr, (r.returncode, r.stderr[-1500:])


def test_with_the_console_off_nothing_is_printed():
    quiet = run('off', close_fd=False)
    loud = run(None, close_fd=False)
    assert quiet.returncode == 0 and loud.returncode == 0
    assert len(loud.stdout) > 200, 'the control prints the load and the trace'
    assert quiet.stdout.strip() == '', quiet.stdout[:300]


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
    r = run(value, close_fd=False)
    assert r.returncode == 0 and len(r.stdout) > 200
