# api binding gap — macos 0591a8c (api source last changed 5792b71) / resources 920eb58

Measured 2026-10-16 (TO, BRIEF_to_bindings A.1) from SOURCE by
`MACOS_resources/mmacos/src/api_binding_gap.py` (the mex codegen's own parse of
`macos_api_mod.F90`; 157 `subroutine` lines, all parsed — a raw grep and the parse
agree name-for-name; the 2 `translate*ID` helpers + 6 helper functions are excluded
as internal).  Re-run: `python3 mmacos/src/api_binding_gap.py`.

Bound means: mex — a `do_*` helper in `mmacos_gen.F` / `mmacos_mex.F` CALLs the
routine; f2py — `pymacos_f2py.f90` renames it (`NAME_impl => NAME`) and defines the
forwarding wrapper; veneer — `macos.py` references `api.<name>`.

**Against the 2026-10-10 record:** it read "136 of 146" but named ~21 missing, so the
two disagreed.  Today: 155 api routines, 136 f2py-bound, 19 missing — the plan's
list less two that were miscounted there: `translate` is the internal
`translateSurfaceID` / `translateEltID` helper pair (not an entry point), and the
`elt_srf_zrn_*` family IS f2py-bound (`elt_srf_zrn_get` alone is unreached from
the veneer — the documented `getEltSrfZern` misname).  No routine added since
10-10 is missing (`modes_get`, `calib_stats_get`, `grid_idx_ovf_get` are bound).

- **f2py missing (19):** `elt_z ray_status_get ray_hist_set ray_pos_hist_get
  elt_info_get elt_obs_get src_seg_get met_calc met_get met_geom_get polelt_set
  polelt_get jmat_elt_get xps_cmd draw_rays_cmd draw_rays_get draw_rays3d_get
  deadband_notes_get fwd_root_flips_get`
- **veneer missing (f2py-bound, no `api.<name>` in macos.py) (9):** `set_src_info
  get_src_size set_src_size elt_srf_csys_set elt_srf_csys_get elt_srf_zrn_get
  elt_srf_ff_any elt_srf_mon_zrn_max_modes cfield_get` (`cfield_get` is reached
  through `cfield_plane_get` plane 0 — same buffer)
- **mex missing (1):** `cfield_get` — covered: the hand-written `complex_field`
  calls `cfield_plane_get`, whose plane 0 is exactly `cfield_get`'s buffer.

Mex codegen refactor-safety baseline: `gen_mex_wrappers.py` run at this SHA
reproduces the committed `mmacos_gen.F` and `mmacos_gen_cmds.txt` byte-for-byte
(`cmp` clean; 137 generated, 1 skipped = `elt_csys_get`, rank 3).

api public routines (from source): 155  (internal helpers excluded: 8)

| layer | bound | missing |
|---|---|---|
| mmacos mex (do_* helper CALLs it) | 154 | 1 |
| pymacos f2py (`_impl =>` + wrapper) | 136 | 19 |
| pymacos veneer (`api.<name>` in macos.py) | 127 | 28 |

| # | api routine | mex | f2py | veneer |
|---|---|---|---|---|
| 1 | `modified_rx` | yes | yes | yes |
| 2 | `prb_elt` | yes | yes | yes |
| 3 | `prb_elt_grp` | yes | yes | yes |
| 4 | `set_src_csys` | yes | yes | yes |
| 5 | `get_src_csys` | yes | yes | yes |
| 6 | `set_src_sampling` | yes | yes | yes |
| 7 | `get_src_sampling` | yes | yes | yes |
| 8 | `sys_units` | yes | yes | yes |
| 9 | `src_info` | yes | yes | yes |
| 10 | `get_src_fov` | yes | yes | yes |
| 11 | `src_wvl` | yes | yes | yes |
| 12 | `src_flux` | yes | yes | yes |
| 13 | `sysprop` | yes | yes | yes |
| 14 | `set_src_info` | yes | yes | **no** |
| 15 | `set_src_fov` | yes | yes | yes |
| 16 | `src_size` | yes | yes | yes |
| 17 | `get_src_size` | yes | yes | **no** |
| 18 | `set_src_size` | yes | yes | **no** |
| 19 | `src_finite` | yes | yes | yes |
| 20 | `elt_grp_max_all` | yes | yes | yes |
| 21 | `elt_grp_max` | yes | yes | yes |
| 22 | `elt_grp_any` | yes | yes | yes |
| 23 | `elt_grp_fnd` | yes | yes | yes |
| 24 | `elt_grp_del` | yes | yes | yes |
| 25 | `elt_grp_del_all` | yes | yes | yes |
| 26 | `elt_grp_set` | yes | yes | yes |
| 27 | `elt_grp_get` | yes | yes | yes |
| 28 | `elt_vpt` | yes | yes | yes |
| 29 | `elt_psi` | yes | yes | yes |
| 30 | `elt_rpt` | yes | yes | yes |
| 31 | `elt_kc` | yes | yes | yes |
| 32 | `elt_kr` | yes | yes | yes |
| 33 | `elt_z` | yes | **no** | **no** |
| 34 | `elt_csys_set` | yes | yes | yes |
| 35 | `elt_csys_get` | yes | yes | yes |
| 36 | `elt_csys_rm` | yes | yes | yes |
| 37 | `elt_srf_csys` | yes | yes | yes |
| 38 | `elt_srf_csys_pos` | yes | yes | yes |
| 39 | `elt_srf_csys_dir` | yes | yes | yes |
| 40 | `elt_srf_csys_set` | yes | yes | **no** |
| 41 | `elt_srf_csys_get` | yes | yes | **no** |
| 42 | `elt_srf_grid_any` | yes | yes | yes |
| 43 | `elt_srf_grid_fnd` | yes | yes | yes |
| 44 | `elt_srf_grid_fnd_type` | yes | yes | yes |
| 45 | `elt_srf_grid_size_max` | yes | yes | yes |
| 46 | `elt_srf_grid_size` | yes | yes | yes |
| 47 | `elt_srf_grid_spacing` | yes | yes | yes |
| 48 | `elt_srf_grid_data_scale` | yes | yes | yes |
| 49 | `elt_srf_grid_data` | yes | yes | yes |
| 50 | `elt_srf_grid_data_add` | yes | yes | yes |
| 51 | `elt_srf_grating_any` | yes | yes | yes |
| 52 | `elt_srf_grating_fnd` | yes | yes | yes |
| 53 | `elt_srf_grating_rule_width` | yes | yes | yes |
| 54 | `elt_srf_grating_order` | yes | yes | yes |
| 55 | `elt_srf_grating_type` | yes | yes | yes |
| 56 | `elt_srf_grating_rule_dir` | yes | yes | yes |
| 57 | `elt_srf_grating_params` | yes | yes | yes |
| 58 | `elt_srf_zrn_any` | yes | yes | yes |
| 59 | `elt_srf_zrn_fnd` | yes | yes | yes |
| 60 | `elt_srf_zrn_coef` | yes | yes | yes |
| 61 | `elt_srf_zrn_type` | yes | yes | yes |
| 62 | `elt_srf_zrn_norm_radius` | yes | yes | yes |
| 63 | `elt_srf_zrn_set` | yes | yes | yes |
| 64 | `elt_srf_zrn_mode_set` | yes | yes | yes |
| 65 | `elt_srf_zrn_get` | yes | yes | **no** |
| 66 | `elt_srf_ff_any` | yes | yes | **no** |
| 67 | `elt_srf_ff_fnd` | yes | yes | yes |
| 68 | `elt_srf_mon_zrn_max_modes` | yes | yes | **no** |
| 69 | `elt_srf_mon_zrn_coef` | yes | yes | yes |
| 70 | `elt_srf_ff_zrn_coef` | yes | yes | yes |
| 71 | `elt_srf_zrn_FreeForm` | yes | yes | yes |
| 72 | `ray_info_get` | yes | yes | yes |
| 73 | `ray_info_set` | yes | yes | yes |
| 74 | `ray_status_get` | yes | **no** | **no** |
| 75 | `ray_hist_set` | yes | **no** | **no** |
| 76 | `ray_pos_hist_get` | yes | **no** | **no** |
| 77 | `elt_info_get` | yes | **no** | **no** |
| 78 | `elt_obs_get` | yes | **no** | **no** |
| 79 | `src_seg_get` | yes | **no** | **no** |
| 80 | `met_calc` | yes | **no** | **no** |
| 81 | `met_get` | yes | **no** | **no** |
| 82 | `met_geom_get` | yes | **no** | **no** |
| 83 | `trace_rays` | yes | yes | yes |
| 84 | `opd_val` | yes | yes | yes |
| 85 | `opd_mask_get` | yes | yes | yes |
| 86 | `opd_ref_set` | yes | yes | yes |
| 87 | `opd_ref_get` | yes | yes | yes |
| 88 | `obs_set` | yes | yes | yes |
| 89 | `spot_cmd` | yes | yes | yes |
| 90 | `spot_get` | yes | yes | yes |
| 91 | `int_cmd` | yes | yes | yes |
| 92 | `int_get` | yes | yes | yes |
| 93 | `compose_start` | yes | yes | yes |
| 94 | `compose_add` | yes | yes | yes |
| 95 | `compose_get` | yes | yes | yes |
| 96 | `window_set` | yes | yes | yes |
| 97 | `window_off` | yes | yes | yes |
| 98 | `beam_set` | yes | yes | yes |
| 99 | `beam_get` | yes | yes | yes |
| 100 | `pol_set` | yes | yes | yes |
| 101 | `pol_get` | yes | yes | yes |
| 102 | `vecdif_set` | yes | yes | yes |
| 103 | `coat_set` | yes | yes | yes |
| 104 | `coat_get` | yes | yes | yes |
| 105 | `polelt_set` | yes | **no** | **no** |
| 106 | `polelt_get` | yes | **no** | **no** |
| 107 | `jmat_elt_get` | yes | **no** | **no** |
| 108 | `rayfield_get` | yes | yes | yes |
| 109 | `ffp` | yes | yes | yes |
| 110 | `pfp` | yes | yes | yes |
| 111 | `cfield_cmd` | yes | yes | yes |
| 112 | `elt_dx_get` | yes | yes | yes |
| 113 | `base_unit_to_metres` | yes | yes | yes |
| 114 | `cfield_get` | **no** | yes | **no** |
| 115 | `cfield_plane_get` | yes | yes | yes |
| 116 | `cfield_apodize` | yes | yes | yes |
| 117 | `cfield_apodize_complex` | yes | yes | yes |
| 118 | `perturb_elt` | yes | yes | yes |
| 119 | `perturb_src` | yes | yes | yes |
| 120 | `xp_set` | yes | yes | yes |
| 121 | `xp_get` | yes | yes | yes |
| 122 | `xp_fnd` | yes | yes | yes |
| 123 | `sxp_fnd` | yes | yes | yes |
| 124 | `xps_cmd` | yes | **no** | **no** |
| 125 | `draw_rays_cmd` | yes | **no** | **no** |
| 126 | `draw_rays_get` | yes | **no** | **no** |
| 127 | `draw_rays3d_get` | yes | **no** | **no** |
| 128 | `ors_run` | yes | yes | yes |
| 129 | `srs_run` | yes | yes | yes |
| 130 | `calib_run` | yes | yes | yes |
| 131 | `modes_get` | yes | yes | yes |
| 132 | `calib_buffer_dims` | yes | yes | yes |
| 133 | `calib_stats_get` | yes | yes | yes |
| 134 | `deadband_notes_get` | yes | **no** | **no** |
| 135 | `fwd_root_flips_get` | yes | **no** | **no** |
| 136 | `grid_dropped_get` | yes | yes | yes |
| 137 | `grid_idx_ovf_get` | yes | yes | yes |
| 138 | `calib_set_var_elt` | yes | yes | yes |
| 139 | `calib_clear_var_elts` | yes | yes | yes |
| 140 | `calib_set_iter` | yes | yes | yes |
| 141 | `calib_set_tol` | yes | yes | yes |
| 142 | `calib_set_target` | yes | yes | yes |
| 143 | `calib_set_beam` | yes | yes | yes |
| 144 | `calib_set_beam_pos_fov` | yes | yes | yes |
| 145 | `calib_set_beam_wt` | yes | yes | yes |
| 146 | `ffcut_set` | yes | yes | yes |
| 147 | `ffcut_get` | yes | yes | yes |
| 148 | `elt_asph_get` | yes | yes | yes |
| 149 | `stop_info_get` | yes | yes | yes |
| 150 | `stop_info_set` | yes | yes | yes |
| 151 | `stop_obj_set` | yes | yes | yes |
| 152 | `n_elt` | yes | yes | yes |
| 153 | `init` | yes | yes | yes |
| 154 | `load_rx` | yes | yes | yes |
| 155 | `save_rx` | yes | yes | yes |
