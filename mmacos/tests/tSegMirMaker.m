classdef tSegMirMaker < matlab.unittest.TestCase
    % SegMirMaker batch-mode regression (design-layer Sprint 2D, S0).
    %
    % The macos.design.segmirmaker_run driver scripts the SegMirMaker
    % interactive dialog over stdin.  These tests pin (a) the answer
    % ordering the driver assembles and (b) the generator's output
    % against the committed references in segmirmaker/test_in:
    % a batch run must reproduce e5pie / e5hex2 -- the same text, and
    % every number to 1e-12 (TEXT_NUMBERS_MATCH).  The comparison was
    % byte for byte until 2026-11-17; the references are x86_64 output
    % and the arm64 build prints a sixteenth digit that differs by one
    % (CCMac's Mac run of issue #76), so it compares numbers now, and
    % test_comparator_refuses_real_changes shows what it still refuses.
    %
    % Runs no engine session (the binary is standalone); skips cleanly
    % if the SegMirMaker binary has not been built on this machine.

    properties
        tin     % segmirmaker/test_in
        bin     % SegMirMaker executable
    end

    methods (TestClassSetup)
        function locate(tc)
            here = fileparts(mfilename('fullpath'));            % mmacos/tests
            res_root = fileparts(fileparts(here));              % the repository root
            tc.tin = fullfile(res_root, 'segmirmaker', 'test_in');
            % Locate whichever build tree exists (ifx on Linux, gfortran on
            % macOS / gfortran builds); ifx preferred when both are present.
            smmdir = fullfile(res_root, 'segmirmaker');
            tc.bin = '';
            for tag = ["build_release_ifx", "build_release_gfortran", ...
                       "build_debug_ifx", "build_debug_gfortran"]
                cand = fullfile(smmdir, tag, 'SegMirMaker');
                if isfile(cand), tc.bin = cand; break; end
            end
            tc.assumeTrue(isfolder(tc.tin), ...
                'segmirmaker/test_in not found');
            tc.assumeTrue(~isempty(tc.bin) && isfile(tc.bin), ...
                'SegMirMaker binary not built (source ./makesegmirmaker.sh)');
        end
    end

    methods (Test)
        function test_pie_reproduces_committed(tc)
            out = macos.design.segmirmaker_run( ...
                fullfile(tc.tin, 'e5mono.in'), ...
                'name', 'e5pie', 'grid', 'Pie', 'rings', 1, ...
                'gap', 50, 'dofs', 6, 'meas_config', 1);
            tc.verify_matches_(out.presc, 'e5pie.presc');
            tc.verify_matches_(out.hx, 'e5pieHx.m');
        end

        function test_hex2_reproduces_committed(tc)
            % 2 rings, all-adjacency edge sensors (config 2).  Note the
            % grid-type string is echoed verbatim into the .presc
            % header, so the committed lowercase 'hex' must be matched.
            out = macos.design.segmirmaker_run( ...
                fullfile(tc.tin, 'e5mono.in'), ...
                'name', 'e5hex2', 'grid', 'hex', 'rings', 2, ...
                'gap', 50, 'dofs', 6, 'meas_config', 2);
            tc.verify_matches_(out.presc, 'e5hex2.presc');
            tc.verify_matches_(out.hx, 'e5hex2Hx.m');
            % Hx dimension sanity from the committed run: 19 segments,
            % 42 shared edges x 2 locations x 3 axes (2026-07-19 model,
            % no anchor row) -> nMeas=252, nState=114.
            hxtxt = fileread(out.hx);
            tc.verifySubstring(hxtxt, 'nMeas=   252;');
            tc.verifySubstring(hxtxt, 'nState=   114;');
        end

        function test_comparator_refuses_real_changes(tc)
            % The must-fail legs of the comparison the two tests above rely on, on the
            % committed e5pie reference itself: round-off passes, a real change does not.
            ref = fileread(fullfile(tc.tin, 'e5pie.presc'));
            tc.verifyTrue(text_numbers_match(ref, ref), 'a file matches itself');
            v = regexp(ref, '\d\.\d{15}E[-+]\d\d', 'match', 'once');       % a 16-digit value
            tc.assertNotEmpty(v, 'the reference prints 16-digit values');
            last = v(17);  bump = char('0' + mod(last - '0' + 1, 10));
            ulp = [v(1:16) bump v(18:end)];                                    % one unit in the last digit
            tc.verifyTrue(text_numbers_match(strrep(ref, v, ulp), ref), 'one unit in the 16th digit is round-off');
            d7 = v;  d7(8) = char('0' + mod(v(8) - '0' + 1, 10));              % the 7th significant digit
            [ok, msg] = text_numbers_match(strrep(ref, v, d7), ref);
            tc.verifyFalse(ok, 'a change in the 7th digit must be refused');
            tc.verifySubstring(msg, 'against', 'and the message names the pair');
            tc.verifyFalse(text_numbers_match(regexprep(ref, 'nSeg=\s*7', 'nSeg=  8', 'once'), ref), 'an integer must be equal');
            tc.verifyFalse(text_numbers_match(strrep(ref, 'Segment', 'Segmint'), ref), 'the text must be equal');
            tc.verifyTrue(text_numbers_match(' -0.0E+00  1.7E-16', '  0.0E+00 -3.1E-17'), 'numerical zeros are equal, sign and digits aside');
            tc.verifyFalse(text_numbers_match('5.01E+01', '5.0E+01'), 'a short print gets no last-digit allowance');
        end

        function test_size_args_mutually_exclusive(tc)
            tc.verifyError(@() macos.design.segmirmaker_run( ...
                fullfile(tc.tin, 'e5mono.in'), ...
                'seg_size', 100, 'ap_diam', 8000), ...
                'macos:design:segmirmaker_run:sizeargs');
        end
    end

    methods (Access = private)
        function verify_matches_(tc, got_file, ref_name)
            [ok, msg] = text_numbers_match(fileread(got_file), fileread(fullfile(tc.tin, ref_name)));
            tc.verifyTrue(ok, sprintf('%s does not match the committed reference: %s', ref_name, msg));
        end

    end
end
