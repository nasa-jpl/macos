## pymacos syntax conventions (Python)

pymacos is a Python package wrapping the same engine through an
f2py-compiled extension.  It works with numpy throughout.

### Setup

Build per the pymacos README (compiles the engine sources into the
`pymacosf90` extension), then:

```python
import pymacos

pymacos.init(256)                    # model size, once per session
pymacos.load('Cassegrain.in')
pymacos.trace_rays(-1)               # trace to the last surface
W = pymacos.opd()                    # numpy (N, N) OPD array
```

### Conventions

- **numpy in, numpy out.**  Vectors and matrices are numpy arrays;
  element indices are 1-based ints, and many functions accept an array
  of element IDs to operate on several at once.
- **Negative element indices** count from the end of the train:
  `-1` is the last surface (image), `-2` the one before it — handy for
  "trace to the exit pupil" idioms like `trace_rays(-2)`.
- **Combined getter/setters.**  Where mmacos has `get_x`/`set_x`
  pairs, pymacos usually has one function whose optional argument
  selects the mode: `elt_vpt(srf)` reads, `elt_vpt(srf, vpt)` writes.
  Part II entries list both languages' forms together.
- **Naming**: snake_case primaries with some camelCase aliases kept
  for backward compatibility (`traceWavefront` → `trace_rays`); use
  the snake_case names in new code.
- **Trace state** works as in mmacos: outputs reflect the most recent
  trace; re-trace after changing the system.
- **Errors** raise Python exceptions carrying the engine message.
- Docstrings are complete — `help(pymacos.trace_rays)` gives the full
  argument contract, including shapes and ranges.

### Coverage relative to mmacos

Part II flags per-entry gaps.  As of this writing the multi-channel
Jacobian drivers, CALIB interface, grid-basis generators, and the
design layer exist on the MATLAB side only; pymacos covers the full
session/source/element/trace/output core plus the FreeForm and
Zernike surface interfaces.
