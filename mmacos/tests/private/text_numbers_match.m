function [ok, msg] = text_numbers_match(got, ref, rtol, ztol)
%TEXT_NUMBERS_MATCH  Two generated text files agree: the same text, and the same numbers to a tolerance.
%   [OK, MSG] = TEXT_NUMBERS_MATCH(GOT, REF) compares two texts that a program wrote (a prescription,
%   a generated .m file).  Everything that is not a number must be identical, blanks aside.  The
%   numbers are compared one by one:
%     - an integer must be equal;
%     - two values both below ZTOL in magnitude are numerical zeros and are equal (their digits and
%       their sign are round-off);
%     - any other value must agree to RTOL, or to one unit in its last printed digit when it is
%       printed to 8 significant digits or more (a full-precision print can differ by that much from
%       one unit of round-off).
%   Defaults RTOL = 1e-12, ZTOL = 1e-9.
%
%   Why not byte for byte: the same source on arm64 and on x86_64 prints a sixteenth significant
%   digit that differs by one (CCMac's Mac run of issue #76, 2026-11-17), and Dave's ruling of
%   2026-10-19 is that a reference is never pinned to digits the two architectures do not share.
%   MSG names the first difference.
if nargin < 3 || isempty(rtol), rtol = 1e-12; end
if nargin < 4 || isempty(ztol), ztol = 1e-9;  end
ok = false;
pat = '(?<![A-Za-z_\d.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?(?![A-Za-z_\d.])';
[na, ta] = regexp(char(got), pat, 'match', 'split');
[nb, tb] = regexp(char(ref), pat, 'match', 'split');
ta = regexprep(ta, '[ \t\r]+', '');  tb = regexprep(tb, '[ \t\r]+', '');
if numel(na) ~= numel(nb)
    msg = sprintf('%d numbers against %d in the reference', numel(na), numel(nb));  return
end
k = find(~strcmp(ta, tb), 1);
if ~isempty(k)
    msg = sprintf('text differs after number %d: "%s" against "%s"', k - 1, clip_(ta{k}), clip_(tb{k}));  return
end
for k = 1:numel(na)
    [a, ia, ua, da] = parse_(na{k});  [b, ib, ub, db] = parse_(nb{k});
    if ia && ib
        same = a == b;
    elseif abs(a) < ztol && abs(b) < ztol
        same = true;
    else
        tol = rtol * max(abs(a), abs(b));
        if da >= 8 && db >= 8, tol = max(tol, 1.0000001 * max(ua, ub)); end
        same = abs(a - b) <= tol;
    end
    if ~same
        msg = sprintf('number %d: %s against %s in the reference', k, na{k}, nb{k});  return
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
