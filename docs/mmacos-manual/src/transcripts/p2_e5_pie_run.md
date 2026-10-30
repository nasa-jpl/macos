    step 1: pie system, 7 segments, width 2666.7 mm, gap 50 mm
            baseline: 12128 src rays, 12128 pass, rmsWFE 3.91181e-05 mm
    step 2: footprint map + boundary overlay -> e5pie_step2_footprints.png
    step 3: apertures emitted (pad 0 mm): center hexagon with 6 vertices, wedges 15 each -> e5pie_polyap.in
    step 4: with apertures: 12008 pass (baseline 12128), rmsWFE 3.88557e-05 mm (-0.67%)
            120 rays clipped; first-fail element(s): [1 2 3 4 5 6 7]
    step 5: 100 mm wedge decenter: 10265 rays pass with apertures, 10351 without --
            the aperture-less trace keeps 86 rays that have no glass under them
    step 6: seg_boundary source = rxpoly (auto-detected from the Rx)
            48 MET beams; launchers sampled from the aperture polygons -> e5pie_step6_met.png
    findings.txt + 5 figures written to <work>/templates/20_segmentation/e5_pie
