# NOTE to CCMac — the manuals rebuild, branch `mr_merge_manuals` (CC, 2026-11-18)

Nothing to build or run on the Mac for this one: the manuals are built on Linux (pandoc and
LibreOffice), and no code changed. It is the last branch of the stack:

    #75 ... #82 (`mr_merge_remove_scratch`), `mr_merge_gmi_mac`, `mr_merge_manuals`

## What to do

1. Open the PR from `mr_merge_manuals` with base `mr_merge_gmi_mac` and the text below. If the GMI
   branch is still being worked on, the base stays `mr_merge_gmi_mac`; this branch does not touch
   any file that one does.
2. If you want a check: open `docs/mmacosMan.pdf`, `docs/macosMan.pdf` and `docs/macosCmdRef.pdf`
   on the Mac and confirm they open and that the installation text names one repository.

## What the branch changes outside `docs/`

Twenty-five files under `app/mmacos/` (template READMEs, the template index and one comment in
`mmacos_setup.m`): the placeholder `<path-to>/mmacos/...` is now `<path-to>/macos/app/mmacos/...`.
The top-level `README.md` has the same correction in its Windows and pymacos steps. No test reads
any of these lines.

## PR text

```markdown
# Rebuild the three manuals for the reorganized repository

Branch `mr_merge_manuals`, on top of `mr_merge_gmi_mac`. No issue number; this is the rebuild that was deferred through issues #76 to #82 so that it would be done once.

## What changes for a user

`docs/macosMan.pdf`, `docs/macosCmdRef.pdf` and `docs/mmacosMan.pdf` describe the repository as it now is: one repository, the engine in `src/`, its data in `data/`, the interfaces in `app/`, the stand-alone tools in `tools/`, the tests in `tests/`. Until now the three PDFs still described two repositories.

| Manual | Pages | Was |
|---|---|---|
| MACOS Manual | 244 | 244 |
| MACOS Command Reference | 121 | 120 |
| mmacos Manual | 147 | 147 |

## What is in the commit

- **The three PDFs**, built from the sources on this branch.
- **Command reference regenerated** from the engine, the bindings and the test journals. It gains the newer design-layer functions (`design.Lens`, `add_exit_pupil`, the camera functions) and one journal in the coverage lists; these were missing before the reorganization began.
- **MACOS manual listings and dialogs regenerated** from the engine (`make transcripts`, the ifx build). Every one reproduced; two changed, the start-up banner and the MRESET dialog, which now show the data files under `data/` and 205 glasses.
- **Appendix C and the Zernike tables regenerated**; neither changed.
- **mmacos manual:** the Part II pages are regenerated from the template READMEs and specs. One page changed, picking up a README sentence written at the November freeze. One listing was stale and is regenerated from its driver: the singlet page quoted a template line naming the glass table's old location.
- **A path placeholder corrected in 75 files.** The manuals and the template READMEs tell the reader to run `<path-to>/mmacos/...`; that is now `<path-to>/macos/app/mmacos/...`. A few other leftovers from the old layout are corrected in the README's Windows section, the mmacos manual's installation steps and the command reference's orientation page.
- **Word edit copies of the MACOS manual and command reference** (`docs/macos-manual/edits/`) refreshed from the new build. Both were unedited copies of the previous build.

## Checks

| Check | Result |
|---|---|
| `make transcripts`, MACOS manual, ifx | all listings and dialogs reproduce; 2 files change, as described |
| `make zernike-table`, `make appendix-c` | no change |
| `gen_part2.py`, mmacos manual | 37 of 38 pages identical; 1 picks up a README change |
| mmacos manual driver `p2_e5_seg` | passes: 2 listings, 1 figure, 5 pinned numbers. This is the page whose tool was re-laid for `tools/segmirmaker` in #87 and could not be run then. |
| mmacos manual driver `p2_singlet` | first run refused, correctly: one listing line had moved. After regenerating that listing it passes: 2 listings, 1 figure, 4 pinned numbers. |
| All three builds | no pandoc warning; contents pages populated; no old directory name left in the text of any PDF |

## Not covered

- **The other 52 drivers of the mmacos manual were not re-run.** They take many hours in all, and one chapter needs eight hours of stage outputs. A search of every generated listing found one line that names a repository path, and that one is fixed, so the listings are consistent with the new layout; but only two drivers were actually run in it.
- **No figure was regenerated.** The two drivers that ran reproduced their figures in content, but the machine they ran on has no display for MATLAB and draws them at a different size, so the committed figures are kept.
- **The polarization validation report** is not rebuilt.
- **The Mac.** The manuals are built on Linux only.
```
