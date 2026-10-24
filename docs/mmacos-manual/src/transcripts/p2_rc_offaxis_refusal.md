    [6] clearance check (no body sits in any beam):
    check_clipping  (family=ritchey_chretien, 3 elements)
      name       kind            ap_r    foot_r    margin  clearnce obstruct  status
      M1         Reflector        0.5    0.4999 9.901e-05    0.9012        0  OK
      M2         Reflector        0.3    0.1034    0.1966  -0.01573       11  ** CLIP
      FP         FocalPlane       0.2 1.627e-05       0.2        --        0  OK
      => layout has CONFLICTS: margin<0 = own beam overfills aperture; clearance<0 = body cuts a foreign beam

    REFUSED: unexpected clipping in the final design
