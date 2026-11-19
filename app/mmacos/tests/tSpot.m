classdef tSpot < matlab.unittest.TestCase
%TSPOT  Ray-surface intersection (spot diagram) veneer + engine coverage.
%   Exercises macos.spot across all three reference frames -- 'telt',
%   'tout' and 'beam' -- on the centrally-obscured Cassegrain fixture.
%
%   Regression guard for the 'beam' reference frame (engine fix in
%   tracesub.F LocalCoord): the beam coordinate frame is derived from a
%   chief-ray trace, and the on-axis chief ray of a centrally-obscured
%   telescope is VIGNETTED.  LocalCoord used to treat that obscuration as
%   "chief ray becomes undefined" and abort, so macos.spot(...,'ref',
%   'beam') failed at every element.  Obscuration sets only the flux flag
%   (not the geometric intersection), so the beam frame is well-defined;
%   LocalCoord now tolerates a vignetted-but-geometrically-valid chief
%   ray and fails only on a genuine geometric miss.

    properties (Constant)
        ModelSize = 128
        RxName    = 'Rx_Cass_FarField.in'
    end

    properties
        rx_path
    end

    methods (TestClassSetup)
        function setupClass(testCase)
            testCase.rx_path = rx_fixture_path(testCase.RxName);
            macos.init(testCase.ModelSize);
        end
    end

    methods (TestMethodSetup)
        function setupMethod(testCase)
            macos.load_rx(testCase.rx_path);
            macos.trace();
        end
    end

    methods (Test)
        function test_beam_ref_succeeds_all_elts(testCase)
            % The core regression: 'beam' must succeed at every element of
            % an obscured telescope, not abort on the vignetted chief ray.
            n = macos.num_elt();
            for e = 2:n
                s = macos.spot(e, 'ref', 'beam');
                testCase.verifyGreaterThan(s.nSpot, 0, ...
                    sprintf('beam spot returned no rays at elt %d', e));
            end
        end

        function test_beam_matches_tout_ray_count(testCase)
            % Same obscuration option -> same rays contribute regardless of
            % the reference frame the spot is expressed in.
            n = macos.num_elt();
            sb = macos.spot(n, 'ref', 'beam');
            st = macos.spot(n, 'ref', 'tout');
            testCase.verifyEqual(sb.nSpot, st.nSpot);
        end

        function test_tout_telt_agree(testCase)
            % TOUT and TELT are both stored element frames; on this on-axis
            % Rx their centroids coincide.
            n = macos.num_elt();
            so = macos.spot(n, 'ref', 'tout');
            se = macos.spot(n, 'ref', 'telt');
            testCase.verifyEqual(so.centroid, se.centroid, 'AbsTol', 1e-12);
        end

        function test_centroid_is_the_mean_of_the_returned_points(testCase)
            % The centroid is the mean of the points the SAME call returns.
            % Until 2026-11-18 the api's spot_get summed EVERY row of the
            % engine's spot array and divided by this spot's ray count, so
            % the rows of an earlier, larger spot stayed in the sum: on
            % e5hex1 at source sampling 31 (552 rays) after the deck's own
            % sampling (10,245 rays) the centroid read -1,156.7 mm in y
            % where the mean of the returned points is -62.32 mm.  The
            % library never zeroed those rows either, so a first spot could
            % read whatever the allocation returned (the dw_dx centroid NaN
            % seen on the Mac under memory pressure).
            % The deck matters: e5hex1's focal-plane spot sits 62 mm off the
            % element vertex, so a stale tail is far from zero.  On this
            % class's centred Cassegrain the tail sums to nearly nothing,
            % which is why no test here saw the defect.
            % Twin: pymacos tests/test_spot_centroid.py.
            macos.load_rx(rx_fixture_path('e5hex1.in'));
            n = macos.num_elt();
            big = macos.spot(n, 'ref', 'tout', 'at', 'elt');
            mb  = mean(big.pts, 1);
            testCase.assertGreaterThan(abs(mb(2)), 10, ...
                'the first spot must sit well off the vertex, or the test is blind');
            testCase.verifyEqual(big.centroid, mb, 'AbsTol', 1e-10*max(abs(big.pts(:))));

            macos.set_src_sampling(31);
            small = macos.spot(n, 'ref', 'tout', 'at', 'elt');
            testCase.assertLessThan(small.nSpot, big.nSpot/4, ...
                'the second spot must be much smaller than the first');
            testCase.verifyTrue(all(isfinite(small.centroid)));
            ms = mean(small.pts, 1);
            testCase.verifyEqual(small.centroid, ms, 'AbsTol', 1e-10*max(abs(small.pts(:))), ...
                sprintf('centroid %s vs the mean of its own points %s', ...
                        mat2str(small.centroid, 6), mat2str(ms, 6)));
        end

        function test_struct_shape(testCase)
            n = macos.num_elt();
            s = macos.spot(n, 'ref', 'beam');
            testCase.verifySize(s.pts, [s.nSpot, 2]);
            testCase.verifySize(s.centroid, [1, 2]);
            testCase.verifySize(s.shift, [1, 4]);
            testCase.verifySize(s.csys, [3, 3]);
        end

        function test_at_elt_option(testCase)
            % 'at','elt' (element vertex) vs default 'chief' both trace.
            n = macos.num_elt();
            s = macos.spot(n, 'ref', 'beam', 'at', 'elt');
            testCase.verifyGreaterThan(s.nSpot, 0);
        end

        function test_stop_info_positive(testCase)
            % get_stop_info returns the element + vertex offset once an
            % element-anchored stop is set (negative/raise path is in
            % tSurfInspect).
            macos.stop(2);
            s = macos.get_stop_info();
            testCase.verifyEqual(s.elt, 2);
            testCase.verifySize(s.offset, [1, 2]);
        end

        function test_session_delegates(testCase)
            % macos.Session.spot must delegate to the package function.
            m = macos.Session(testCase.ModelSize);
            m.load_rx(testCase.rx_path);
            m.trace();
            n = m.num_elt();
            s = m.spot(n, 'ref', 'beam');
            testCase.verifyGreaterThan(s.nSpot, 0);
        end
    end
end
