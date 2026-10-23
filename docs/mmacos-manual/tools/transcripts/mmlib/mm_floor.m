function txt = mm_floor(name, txt, pattern, bound)
%MM_FLOOR  A number set by a numerical zero (a condition number over a noise-floor singular value)
%   prints as the bound it is known to exceed, and is pinned by it (the cross-architecture rule).
%   PATTERN is a regexp whose one token is the number; its first match must exceed BOUND (pin NAME,
%   rule 'min:<bound>'), and every match's number is printed as '> BOUND'.  The mirror of mm_zero.
txt = char(txt);
[tok, ext] = regexp(txt, pattern, 'tokens', 'tokenExtents');
if isempty(tok), error('mm_floor:nomatch', '%s: no match for %s', name, pattern); end
shown = sprintf('> %g', bound);
mm_pin(name, str2double(tok{1}{1}), sprintf('min:%g', bound), shown);
for k = numel(ext):-1:1                      % last first, so earlier extents stay valid
    e = ext{k}(1, :);
    txt = [txt(1:e(1)-1) shown txt(e(2)+1:end)];
end
end
