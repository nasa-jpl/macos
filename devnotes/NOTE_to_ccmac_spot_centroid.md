# NOTE to CCMac — the spot-centroid fix, branch `mr_merge_spot_centroid` (CC, 2026-11-18)

This branch is one commit on top of `mr_merge` at `b5b37c34`, the merged series you are running the
full gate on. It changes engine source, so it needs a Mac gate of its own, after that one. Read
this note from `devnotes/` in your clone.

## What it is

The tDwDx NaN you saw on `mr_merge_param_print` was not a transient from machine load. It was an
engine defect that predates that branch, and this branch fixes it.

`spot_get` (`src/macos_api_mod.F90`) computed the centroid as the sum of ALL `mRay` rows of
`RaySpot`, divided by `N`. SPOT writes rows `1..N` only. The rows after them held:

- **an earlier, larger spot.** Deterministic. On Linux, `e5hex1.in` at source sampling 31 (552
  rays) after the deck's own sampling (10,245 rays) gave centroid y = −1,156.7 mm where the mean
  of the returned points is −62.32 mm.
- **whatever the allocation returned.** The library never zeroed `RaySpot`. A fresh block is
  zero; a recycled block is not, and if it holds anything that is not a finite number the
  centroid is NaN. Load matters only because it changes what the allocator returns.

The second is my reading of your NaN. I did not reproduce a NaN on Linux. Your three failures
fit it exactly: `tDwDx` has five tests that ask for `dcdx`, and the three that failed are the
three that can notice a NaN. One of the other two checks shapes only; the other compares two
results with `verifyEqual`, which counts NaN equal to NaN.

PR #96 (the print fix) is merged with its text saying the NaN was a transient from load. If you
can still edit that description, change the sentence to: a known engine defect, unrelated to that
change, fixed on `mr_merge_spot_centroid`.

## What changed

1. **`src/macos_api_mod.F90`**: `spot_get` sums rows `1:N`. `init` zeroes `PixArray`, `OPDMat`
   and `RaySpot` when it allocates them, as the command-line program does for its own copies.
2. **`app/mmacos/tests/tSpot.m`**: new `test_centroid_is_the_mean_of_the_returned_points`.
3. **`app/pymacos/tests/test_spot_centroid.py`**: new file, two tests.
4. **`app/mmacos/tests/tDwDx.m`**: `test_multi_compute_los` and
   `test_trans_output_default_is_base_x_cbm` now check that `dcdx` is finite.
5. **`src/notes/hygiene.md`**: the record.

## What to do on the Mac

1. **Before you rebuild**, with the `mr_merge` build of your full gate still in place, confirm
   the defect on arm64. In MATLAB:

   ```matlab
   macos.init(128);  macos.load_rx('<repo>/app/pymacos/tests/Rx/e5hex1.in');
   n = macos.num_elt();
   a = macos.spot(n, 'ref', 'tout', 'at', 'elt');
   macos.set_src_sampling(31);
   b = macos.spot(n, 'ref', 'tout', 'at', 'elt');
   [b.centroid; mean(b.pts, 1)]
   ```

   On Linux the two rows are `[2.07e-05 -1156.73]` and `[3.94e-16 -62.324]`. Send me what the
   Mac prints.
2. `git fetch origin && git switch mr_merge_spot_centroid`.
3. Rebuild the engine in the existing tree (`source ./makems.sh release gfortran`), then the mex
   and pymacos.
4. Repeat step 1. The two rows must now agree.
5. The mmacos fast suite and the pymacos suite, with nothing else running, and then the mmacos
   fast suite once more **beside the same jobs that were running when the NaN appeared** (the
   cmd gate, pymacos, an isolated MATLAB). The second run is the one that matters.

The command-line gates are not needed: the `macos` program does not contain the changed file.

## What I expect on the Mac

- Step 1 shows the same defect as Linux, with the same numbers to the digits shown.
- After the fix, `tDwDx` passes with and without load. If a NaN appears again after the fix, the
  cause is not this one: send the test name and, from the same session,
  `s = macos.spot(...)`, `any(isnan(s.centroid))` and `any(isnan(s.pts(:)))`.
- Counts: one more mmacos test and two more pymacos tests than on `mr_merge`. Linux: mmacos fast
  817 pass, pymacos 7,032 pass and 1 skipped.

## Linux results (CC; `MACOS_HOME` unset, a separate checkout)

| Check | Result |
|---|---|
| Engine rebuild, gfortran and ifx; mmacos mex; pymacos | all succeed |
| New mmacos test against the engine before the fix | fails: centroid y −1,156.73 against −62.324 |
| New pymacos tests against the engine before the fix | both fail, same numbers |
| mmacos fast suite | 817 pass, 0 fail, 0 incomplete (816 before, plus the new test) |
| pymacos suite | 7,032 pass, 0 fail, 1 skipped (7,030 before, plus the two new tests) |
| Command reference regenerated from the changed source | no difference |

The binaries were built on the print-fix commit (`81e8cc5f`) plus this change. `mr_merge` adds
twelve files to that (the GMI regression scripts, `MAC_PORT.md`, two notes), none of them compiled
into the engine, the mex or pymacos. A first mmacos run reported 45 incomplete because SegMirMaker
was not built in the new checkout; it was built and the whole suite run again.

Not run: the command-line gates, the GMI regression, the large-model mmacos groups, the manual
builds. No manual driver calls the spot centroid or `compute_los`.

## Reporting

As before. Open the PR from `mr_merge_spot_centroid` with base `mr_merge`, with the text below and
your Mac lines added. No issue number. Nothing is stacked on it, so any merge method is safe.

## PR text

```markdown
# Spot centroid: sum only the rays of this spot

Branch `mr_merge_spot_centroid`, on top of `mr_merge` (`b5b37c34`, the merged reorganization series). No issue number.

## What changes for a user

`macos.spot` in mmacos and `spot` in pymacos return a spot's points and its centroid. The centroid was wrong whenever an earlier spot in the same session had more rays, and it could be NaN on a first spot. It is now the mean of the points returned by the same call, always.

Measured on `e5hex1.in` (a millimetre deck), model size 128, spot at the focal plane about the element vertex:

| Spot | Rays | Centroid y before | Centroid y now | Mean of the returned points |
|---|---|---|---|---|
| the deck's own source sampling | 10,245 | −62.324 | −62.324 | −62.324 |
| then source sampling 31 | 552 | **−1,156.7** | −62.324 | −62.324 |

The points themselves were always right. So was the command-line program, which does not use this routine.

The centroid sensitivity `dcdx` that `dw_dx(..., 'compute_los', true)` and its relatives return is a difference of two such centroids. Equal errors cancelled, so it was wrong only when the perturbation changed how many rays arrived, or when the error was not a finite number.

## Why

`spot_get` in `src/macos_api_mod.F90` summed every row of the engine's spot array and divided by the number of rays in this spot. The SPOT command writes only that many rows. The rest held one of two things:

- **Rows of an earlier, larger spot.** This is the case in the table. It repeats every run, on both compilers.
- **Rows no spot had written.** The library allocated the array without zeroing it; the command-line program zeroes its own copy. A fresh block of memory is zero and the answer came out right. A recycled block is not.

The second case is our reading of a Mac result on 2026-11-18: three `tDwDx` centroid tests returned NaN in a suite run beside other jobs, and passed when the suite was re-run alone. Machine load cannot change a deterministic calculation, but it can change which block of memory an allocation returns. The NaN itself was not reproduced on Linux; the first case was.

Two things kept the defect out of sight. Every spot test ran on a centred telescope, where the extra rows sum to nearly nothing. And of the five `tDwDx` tests that ask for `dcdx`, two could not notice a NaN: one checks array shapes only, and the other compares two results with `verifyEqual`, which counts NaN as equal to NaN.

## What is in the commit

- **`src/macos_api_mod.F90`**: the centroid sums rows 1 to N. The library zeroes `PixArray`, `OPDMat` and `RaySpot` when it allocates them.
- **`app/mmacos/tests/tSpot.m`**: a new test, on `e5hex1.in`, that the centroid of a small spot taken after a large one is the mean of its own points.
- **`app/pymacos/tests/test_spot_centroid.py`**: the same test, and a second one that a spot's centroid does not depend on which deck was traced before it.
- **`app/mmacos/tests/tDwDx.m`**: the two tests that could not notice a NaN now check that `dcdx` is finite.
- **`src/notes/hygiene.md`**: the record.

## Test results (Linux)

All runs with `MACOS_HOME` unset, in a separate checkout. The binaries were built on the print-fix commit (`81e8cc5f`) plus this change. The branch as pushed sits on `mr_merge`, which adds twelve files to that: the GMI regression scripts, `MAC_PORT.md` and two notes. None of them is compiled into the engine, the mex or pymacos.

| Check | Result |
|---|---|
| Engine rebuild, gfortran and ifx; mmacos mex; pymacos | all succeed |
| New mmacos test against the engine before the fix | fails: centroid y −1,156.73 against −62.324 |
| New pymacos tests against the engine before the fix | both fail, same numbers |
| mmacos fast suite | 817 pass, 0 fail, 0 incomplete (816 before, plus the new test) |
| pymacos suite | 7,032 pass, 0 fail, 1 skipped (7,030 before, plus the two new tests) |
| Command reference regenerated from the changed source | no difference |

A first mmacos run reported 772 pass and 45 incomplete: the SegMirMaker binary had not been built in the new checkout. It was built and the whole suite run again; the table shows that run.

Not run for this change: the command-line gates (the `macos` program does not contain the changed file), the GMI regression, the large-model mmacos groups, the manual builds. No manual driver calls the spot centroid or `compute_los`.

## Not covered

- **Mac and Windows.** Not run here.
- **The NaN on the Mac was not reproduced.** The fix removes the only read of unwritten rows in the library, and zeroes the arrays, so it should not recur. If it does, the cause is elsewhere.
- **Other getters.** The three unzeroed arrays have three other direct readers in the library (`OPD`, the spot points, the pixel array). Each reads a bounded block. Other buffers were not audited.
```
