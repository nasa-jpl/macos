function result = test07_tolerance_refuses(opts)
% The reference tolerance must still refuse a real change (non-vacuity).
% No engine call: for each prescription, a map that differs from another
% by HALF the tolerance must pass and one that differs by TWICE the
% tolerance must fail.  If someone loosens the rule until it passes
% everything, this test goes red.

    name = 'test07_tolerance_refuses';
    result = make_result(name);
    rxs = {'e2e_pie', 'e5hex1'};
    for k = 1:numel(rxs)
        [tol, L] = ref_tolerance(rxs{k}, opts);
        a = zeros(4);
        b = a;  b(2, 3) = 0.5 * tol;
        c = a;  c(2, 3) = 2.0 * tol;
        ok_small = compare_within(a, b, tol, 'half the tolerance');
        ok_large = compare_within(a, c, tol, 'twice the tolerance');
        if ~ok_small || ok_large || ~(tol > 0) || tol > 1e-9 * L
            result.pass = false;
            result.msg = sprintf(['[%s] FAIL  %s: tol %.3e (L %.6g): half passes = %d, ' ...
                                  'twice passes = %d'], name, rxs{k}, tol, L, ok_small, ok_large);
            return
        end
    end
    result.pass = true;
    result.msg = sprintf('[%s] PASS  half the tolerance passes and twice is refused, both prescriptions', name);
end

function r = make_result(name)
    r.name = name;
    r.pass = false;
    r.msg = sprintf('[%s] not run', name);
end
