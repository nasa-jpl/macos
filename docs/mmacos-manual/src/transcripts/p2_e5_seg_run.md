    [1] segmenting e5mono (rings=1, Hex, gap=50)
    [2] MET truss: 6 launchers/seg + 6 around elt 11 -> 6 fiducials on elt 8
    [3] segmented+met system: 13 elts, 10246 rays, rmsWFE 3.62e-05
    [3b] dwdx FD over 54 DOFs x 10244 OPD samples

    [4] post-control wavefront residual (rms, control = -pinv(dwdx)*w_hat):
       prior (no sensing)    12860.131 nm rms
       edge sensors only        55.322 nm rms
       MET only                472.908 nm rms
       edge + MET                4.704 nm rms
       sensing removes 100.0% of the prior wavefront variance (edge+MET)
       Monte-Carlo (200 draws): 4.759 nm rms  (analytic 4.704) [1.2%]

    [5] artifacts saved beside the script (e5_seg_met.in, .mat, metric + layout png)
