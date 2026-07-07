## mmacos syntax conventions (MATLAB)

mmacos is a MATLAB package: everything lives in the `macos.` namespace
(the `+macos` package directory), backed by one MEX gateway
(`mmacos.mexa64`) linked against the engine.

### Setup

Build the engine (`makeall.sh` or `makegfortran.sh`), then build the
mex per the mmacos README, and put the mmacos `src/` directory on the
MATLAB path.  Then:

```matlab
macos.init(256);                    % model size, once per session
macos.load_rx('Cassegrain.in');
macos.trace();
W = macos.opd();                    % N x N OPD matrix
```

### Conventions

- **Everything numeric is `double`.**  The MEX dispatcher reads all
  numeric inputs as REAL(8); pass `double(...)` values — never
  `int32` or `logical` types.  Element indices are 1-based doubles.
- **Getter/setter pairs.**  Element and source properties come as
  `get_*`/`set_*` pairs (`get_elt_vpt`/`set_elt_vpt`, ...).  Setters
  write engine state directly; for rigid-body moves of real optics
  prefer `macos.perturb`, which keeps linked coordinate frames
  consistent.
- **Trace state.**  Output functions (`opd`, `spot`, `intensity`, ...)
  read the *most recent trace*; call `macos.trace()` (or a driver that
  traces internally) after any change.  Setters mark the state dirty
  so stale results are not silently reused.
- **Units** follow the loaded prescription's BaseUnits/WaveUnits;
  `macos.sys_units` reports them, `macos.cbm` gives the base-unit→
  meter conversion.
- **Errors**: functions raise MATLAB errors with the engine message;
  most also return an `ok` flag when called with output arguments.
- `macos.Session` provides an object wrapper (load/trace/query in one
  handle) if you prefer OO style; the plain functions are equivalent.

### Higher-level drivers

Beyond one-call wrappers, mmacos adds drivers with no CLI equivalent —
documented in Parts II and III:

- `macos.dw_dx`, `dw_dz_zernike`, `dw_dsurf`, `dw_dgrid` (+ `_multi`
  variants): wavefront Jacobians over perturbation channels.
- `macos.calib*`: the native CALIB optimizer interface.
- `macos.segment_grid_basis`, `gs_zernike_segment_basis`,
  `write_grid_file`: segment/Zernike grid-basis generation.
- `macos.design.*`: parametric telescope construction and wavefront
  solving (Part III).
- `macos.channels.*`: the perturbation-channel classes the Jacobian
  drivers consume (Part III).

Help is embedded: `help macos.opd`, `help macos.dw_dx`, etc.
