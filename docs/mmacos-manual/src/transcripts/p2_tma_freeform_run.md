    ====================================================================
     Off-axis unobscured 3-mirror freeform | D=1.00 m | half-field +-3'
    ====================================================================

    [1] coaxial Korsch TMA seeded: on-axis 0.0103 waves (4 elts)
    [2] unobscured off-axis: decenter 0.89 m (0.89 D); eccentric-pupil WFE 3.8 waves
    [3] 3-conic + rigid optimize: worst-field 0.0213 waves (3-mirror limit)
    [4] + freeform M[2 3] (modes [5 6 7 8 11], CALIB OptZern, radii/conics held):
        worst-field 0.0213 -> 0.0176 waves  (1.2x) -> DIFFRACTION-LIMITED
    [5] clearance: 4/4 optics clear -> UNOBSCURED
    [6] WFE-across-field map: <work>/templates/10_telescopes/tma_freeform/tma_freeform_wfe.png
        per-field before/after: <work>/templates/10_telescopes/tma_freeform/tma_freeform_fields.png
        orthographic layout: <work>/templates/10_telescopes/tma_freeform/tma_freeform_layout.png
    [7] saved: <work>/templates/10_telescopes/tma_freeform/tma_freeform.in
               + <work>/templates/10_telescopes/tma_freeform/tma_freeform.mat

    --------------------------------------------------------------------
     D=1.00 m | off-axis 0.89 D (unobscured) | half-field +-3'
     worst-field RMS WFE: 3-conic limit 0.0213 -> + freeform 0.0176 waves  (1.2x)
     freeform on M[2 3] holds all radii/conics; the 3rd mirror gives the field-
     dependent control a 2-mirror cannot (cf. wf2_freeform: ~1x over a wide field)
    ====================================================================
