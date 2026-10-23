    bench_layout: 13 elements, chief path 2286.320 mm
       1  Baffle       Obscuration s=  250.000  vpt=[  250.000     0.000     0.000]
       2  L1flat       Refractor   s=  496.193  vpt=[  496.193     0.000     0.000]
       3  L1pow        Refractor   s=  500.000  vpt=[  500.000     0.000     0.000]
       4  BSrefl       Reflector   s=  650.000  vpt=[  650.000     0.000     0.000]
       5  DM           Reflector   s=  900.000  vpt=[  650.000  -250.000     0.000]
       6  BStxf        Refractor   s= 1150.000  vpt=[  650.000     0.000     0.000]
       7  BStxb        Refractor   s= 1159.071  vpt=[  652.633     8.681     0.000]
       8  Fold1        Reflector   s= 1459.071  vpt=[  652.633   308.681     0.000]
       9  Fold2        Reflector   s= 1759.071  vpt=[  952.633   308.681     0.000]
      10  L2pow        Refractor   s= 1959.071  vpt=[  952.633   108.681     0.000]
      11  L2flat       Refractor   s= 1964.609  vpt=[  952.633   103.143     0.000]
      12  FocalMask    Reference   s= 2209.071  vpt=[  952.633  -141.319     0.000]
      13  Detector     FocalPlane  s= 2286.320  vpt=[  952.633  -218.569     0.000]
    seed prescription: <work>/bench_layout_seed.in
    parameter schematic: <work>/bench_layout_params.png
    rays: 3002 at baffle, 3002 at detector
    builder-vs-engine chief agreement: max < 1e-09 mm

    === Stage A: collimate at DM (elt 5), solve L1 ===
      L1: Kr 250 -> 236.866, Kc 0 -> -0.5829;  ray spread 1.71 -> < 1e-06 mrad rms

    === Stage B: focus at mask (elt 12), solve L2 ===
      L2: Kr -125 -> -124.076, Kc -2.25 -> -0.5826;  focus spot 403 -> 0.0177 um rms

    === Stage C: detector at the DM pupil image ===
      thin-lens seed off by -0.824 mm; tilt-induced beam shift 1.07 -> < 1e-05 um
      pupil: 22.47 mm radius at DM -> 6.96 mm at detector (mag 0.310, thin-lens 0.309)
    rendered: <work>/bench_layout_view_rx.png
              <work>/bench_layout_view_std.png

    DONE.  Optimized prescription: <work>/bench_layout_opt.in
