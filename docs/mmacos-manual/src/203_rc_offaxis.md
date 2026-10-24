### `rc_offaxis`: Unobscured off-axis Ritchey-Chretien

Build an on-axis 2-mirror RC, take it OFF-AXIS as an eccentric-pupil section, clear the beam, and re-optimize to diffraction-limited.

*This template has no README; the text above is generated from its entry script's header (`docs/mmacos-manual/part2/rc_offaxis.md`, `tools/gen_readme_part2.py`).*

**Run it, as it stands today.**

```matlab
cd('<path-to>/mmacos/templates/10_telescopes/rc_offaxis');
run(fullfile(pwd,'example_rc_offaxis.m'));
```

**It refuses:**

<!-- transcript: p2_rc_offaxis_refusal -->

**Its figure:** none: the run stops at its clipping check.

**The number it reproduces.** There is no number to reproduce today. The design reaches its final stage and then fails its own proof: the clipping check finds the secondary cutting a foreign beam, clearance %pin:p2_rc_offaxis_clr% m (finding (rrr)). The template has not changed since 2026-08-19. `rc_unobscured`, its parameterized sibling, runs clear (its page above).

**Adapt it.** *"Using templates/10_telescopes/rc_offaxis, find which beam M2 cuts in the final design and by how much at each field point, and tell me the smallest extra decenter that clears it, before changing anything."*

