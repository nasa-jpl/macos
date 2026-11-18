function [ok, msg] = text_numbers_match(got, ref, rtol, ztol)
%TEXT_NUMBERS_MATCH  Two generated text files agree: the same text, and the same numbers to a tolerance.
%   [OK, MSG] = TEXT_NUMBERS_MATCH(GOT, REF) compares two texts that a program wrote (a prescription,
%   a generated .m file).  The two must have the same lines; on each line everything that is not a
%   number must be identical, blanks aside, and the numbers are compared one by one:
%     - an integer must be equal;
%     - two values both below ZTOL in magnitude are numerical zeros and are equal (their digits and
%       their sign are round-off);
%     - any other value must agree to RTOL times the SCALE OF ITS LINE -- the largest magnitude
%       among the reference's values on that line, or the value itself if it is larger -- or to one
%       unit in its last printed digit when it is printed to 8 significant digits or more (a
%       full-precision print can differ by that much from one unit of round-off).
%   Defaults RTOL = 1e-12, ZTOL = 1e-9.
%
%   Why the scale of the line: a line holds one vector or one matrix row, and a small entry beside
%   large ones is a difference of large terms, good to their round-off and no better.  e5hex2Hx.m
%   row 62 carries -0.364217987D-08 beside 0.399952076D+03; arm64 prints -0.364217951D-08, which is
%   3.6e-16 away -- one unit of round-off of the terms it came from, but 1e-7 of the entry itself
%   (CCMac's Mac run of issue #77).  Against the line's scale it is 9e-19; a change of that entry
%   to -0.364D-06 is 9e-10 and is refused.
%
%   Why not byte for byte: the same source on arm64 and on x86_64 prints a sixteenth significant
%   digit that differs by one (CCMac's Mac run of issue #76, 2026-11-17), and Dave's ruling of
%   2026-10-19 is that a reference is never pinned to digits the two architectures do not share.
%   MSG names the first difference.
if nargin < 3 || isempty(rtol), rtol = 1e-12; end
if nargin < 4 || isempty(ztol), ztol = 1e-9;  end
ok = false;
pat = '(?<![A-Za-z_\d.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?(?![A-Za-z_\d.])';
la = regexp(regexprep(char(got), '\r', ''), '\n', 'split');
lb = regexp(regexprep(char(ref), '\r', ''), '\n', 'split');
if numel(la) ~= numel(lb)
    msg = sprintf('%d lines against %d in the reference', numel(la), numel(lb));  return
end
for L = 1:numel(lb)
    [na, ta] = regexp(la{L}, pat, 'match', 'split');
    [nb, tb] = regexp(lb{L}, pat, 'match', 'split');
    ta = regexprep(ta, '[ \t]+', '');  tb = regexprep(tb, '[ \t]+', '');
    if numel(na) ~= numel(nb)
        msg = sprintf('line %d: %d numbers against %d in the reference', L, numel(na), numel(nb));  return
    end
    k = find(~strcmp(ta, tb), 1);
    if ~isempty(k)
        msg = sprintf('line %d: text differs: "%s" against "%s"', L, clip_(ta{k}), clip_(tb{k}));  return
    end
    n = numel(nb);  a = zeros(1, n);  b = a;  ia = false(1, n);  ib = ia;  ua = a;  ub = a;  da = a;  db = a;
    for k = 1:n
        [a(k), ia(k), ua(k), da(k)] = parse_(na{k});  [b(k), ib(k), ub(k), db(k)] = parse_(nb{k});
    end
    scale = max([0, abs(b(~ib))]);                  % the largest real value on the reference's line
    for k = 1:n
        if ia(k) && ib(k)
            same = a(k) == b(k);
        elseif abs(a(k)) < ztol && abs(b(k)) < ztol
            same = true;
        else
            tol = rtol * max([abs(a(k)), abs(b(k)), scale]);
            if da(k) >= 8 && db(k) >= 8, tol = max(tol, 1.0000001 * max(ua(k), ub(k))); end
            same = abs(a(k) - b(k)) <= tol;
        end
        if ~same
            msg = sprintf('line %d, number %d: %s against %s in the reference', L, k, na{k}, nb{k});  return
        end
    end
end
ok = true;  msg = '';
end

function [v, isint, unit, ndig] = parse_(tok)
% value, whether it is an integer token, one unit in its last printed digit, its significant digits
t = regexprep(tok, '[Dd]', 'E');
v = str2double(t);
m = regexp(t, '^[-+]?(\d*)\.?(\d*)', 'tokens', 'once');           % the mantissa's two digit runs
e = regexp(t, '[Ee]([-+]?\d+)$', 'tokens', 'once');               % the exponent, if it has one
isint = isempty(regexp(t, '[.Ee]', 'once'));
ex = 0;  if ~isempty(e), ex = str2double(e{1}); end
unit = 10 ^ (ex - numel(m{2}));
ndig = numel(regexprep([m{1} m{2}], '^0+', ''));
end

function s = clip_(s)
s = strrep(s, newline, '\n');
if numel(s) > 60, s = [s(1:57) '...']; end
end
