"""The engine finds its data files with no MACOS_HOME (issue #75; twin of mmacos tDataDir).

The build records its own data directory (BUILD.md, "Where the engine finds its data
files"), so `init` must work from a directory that holds no macos_param.txt -- with
MACOS_HOME unset, and with MACOS_HOME naming a directory that does not exist (a stale
value in a shell start-up file: the engine warns and falls through).

Each leg is a subprocess: a failed lookup STOPS the process, so it cannot run in the
suite's own interpreter.  Before the fix the built-in directory was dead (a macro name
inside quotes) and both legs stopped at init -- with exit status 0, which is why the
legs wait for INIT_OK and do not trust the status alone.
"""
import os
import subprocess
import sys
from pathlib import Path

from context import pymacos as m

SRC = str(Path(m.__file__).resolve().parents[1])     # the src/ this suite itself loaded
LEG = ("import sys; sys.path.insert(0, {src!r}); import pymacos.macos as m; "
       "m.init(128); print('INIT_OK', flush=True)")


def _init_in(tmp_path, macos_home):
    env = dict(os.environ)
    env.pop('MACOS_HOME', None)
    if macos_home is not None:
        env['MACOS_HOME'] = macos_home
    assert not (tmp_path / 'macos_param.txt').exists()      # the working directory cannot supply it
    r = subprocess.run([sys.executable, '-c', LEG.format(src=SRC)], cwd=tmp_path, env=env,
                       capture_output=True, text=True, timeout=300)
    assert 'INIT_OK' in r.stdout and r.returncode == 0, \
        f'init did not complete (status {r.returncode})\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}'


def test_init_with_macos_home_unset(tmp_path):
    _init_in(tmp_path, None)


def test_init_with_a_stale_macos_home(tmp_path):
    _init_in(tmp_path, str(tmp_path / 'no_such_directory'))
