"""The spot centroid is the mean of the points the SAME call returns -- twin of
mmacos tSpot.test_centroid_is_the_mean_of_the_returned_points.

Until 2026-11-18 the api's spot_get summed EVERY row of the engine's spot array
and divided by the number of rays of this spot.  The rows after that number
still held an earlier spot with more rays, so the centroid of a small spot taken
after a large one carried the large one's tail: on e5hex1 at source sampling 31
(552 rays) after the deck's own sampling (10,245 rays) it read -1,156.7 mm in y
where the mean of the returned points is -62.32 mm.  The same rows were never
zeroed by the library, so a first spot could also read whatever the allocation
returned (the dw_dx centroid NaN seen on the Mac under memory pressure).

The deck matters: e5hex1's focal-plane spot sits 62 mm off the element vertex,
so a stale tail is far from zero.  On a centred deck (the Cassegrain fixture of
the other spot tests) the tail sums to nearly nothing and the defect is hidden.
"""
import numpy as np
from context import pymacos as m
from context import rx_path

MODEL = 128
SMALL = 31          # source sampling of the second, smaller spot


def _spot_at_image():
    # element-centred, in the TOUT frame: the call dw_dx makes for its centroid
    pts, cen, shift, csys = m.spot(-1, vpt_center=True, beam_csys=2)
    return pts, cen


def _tol(pts):
    return 1e-10 * np.abs(pts).max()


def test_centroid_is_the_mean_of_the_returned_points():
    m.init(MODEL)
    m.load(str(rx_path('e5hex1.in')))
    big, cen_big = _spot_at_image()
    # the test bites only if the first spot is larger and well off the vertex
    assert abs(big.mean(axis=0)[1]) > 10.0
    np.testing.assert_allclose(cen_big, big.mean(axis=0), rtol=0, atol=_tol(big))

    m.src_sampling(SMALL)
    small, cen_small = _spot_at_image()
    assert small.shape[0] < big.shape[0] // 4
    assert np.isfinite(cen_small).all()
    np.testing.assert_allclose(cen_small, small.mean(axis=0), rtol=0, atol=_tol(small))


def test_centroid_does_not_depend_on_the_spot_before_it():
    # the same small spot, taken after two DIFFERENT large spots: the centred
    # Cassegrain's (11,484 rays) and e5hex1's own (10,245 rays, 62 mm off the vertex)
    m.init(MODEL)
    m.load(str(rx_path('Rx_Cass_FarField.in')))
    _spot_at_image()
    m.load(str(rx_path('e5hex1.in')))
    m.src_sampling(SMALL)
    after_cass, cen_after_cass = _spot_at_image()

    m.load(str(rx_path('e5hex1.in')))
    _spot_at_image()
    m.src_sampling(SMALL)
    after_e5, cen_after_e5 = _spot_at_image()

    assert after_e5.shape == after_cass.shape
    np.testing.assert_allclose(cen_after_e5, cen_after_cass, rtol=0, atol=_tol(after_e5))
