"""Where this repository's files are, found from this file's own location
(app/pymacos/tests), so a test reads the checkout it lives in whatever the clone
is called, wherever it is, and whatever MACOS_HOME holds."""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
ENGINE_DIR = REPO_ROOT / 'src'      # the engine's Fortran source
