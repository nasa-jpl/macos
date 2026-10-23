# mmacos manual -- open questions and findings for Dave (TO, from 2026-10-23)

The report on Dave's return gathers this list (BRIEF_to_mmacos_manual.md, Report).

## Defaults taken so the work could start (overturn in one line)
1. **The thread's telescope = e2e** (`templates/80_end_to_end/e2e`, the 4 m f/18 folded Korsch TMA, all
   seven stages).  The brief first described it as "the 3 m f/20 TMA", which is e2e2 (design stages
   only); CC corrected the brief 2026-10-23.  e2e2 is chapter 1's pointer to the improved design flow.
2. **The forty uncalled helpers** are documented in one table (Part II), not deleted.

## Findings met while writing (no library change made; PLAN_CONSOLIDATION 1.2 if Dave wants them fixed)
- **ch. 1: the e2e template's committed stage-1 outputs are stale.**  A fresh run of the unchanged
  `s1_telescope.m` gives conics K = [-0.98972 -1.51153 -0.97231] (committed report: -0.98946 -1.50791
  -0.98712), worst +-1' WFE 0.0138 waves (committed 0.0158), M1 hole r 0.232 m (committed 0.095 m --
  the `hole_min_r_m` floor in `e2e_params.m` postdates the committed report).  The manual prints the
  live run; the template's `s1_report.txt` / `.in` / `.mat` / PNGs are from August.  Regenerating the
  committed template outputs is Dave's call (stages 2-7 consume them).
- **ch. 1: `design_report`'s "EFL" is the local plate scale at the biased field point** (chief-ray
  image displacement per field angle, +-0.1' about the 2' bias): 72.32 m, vs the engine's SYSPROP EFL
  72.14 m on the emitted deck.  Both are right for what they measure; the report's label could say
  "plate-scale EFL at the field point".

## Questions
- Authorship / acknowledgment block in the front matter: copied from the MACOS manual as a placeholder.
