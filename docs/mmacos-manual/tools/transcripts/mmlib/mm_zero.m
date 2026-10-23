function txt = mm_zero(name, txt, pattern, bound)
%MM_ZERO  A numerical zero prints as its bound and is pinned by it (the cross-architecture rule).
%   PATTERN is a regexp whose one token is the number; its first match must be below BOUND in
%   magnitude (pin NAME, rule 'abs:<bound>'), and every match's number is printed as '< BOUND'.
txt = char(txt);
[tok, ext] = regexp(txt, pattern, 'tokens', 'tokenExtents');
if isempty(tok), error('mm_zero:nomatch', '%s: no match for %s', name, pattern); end
shown = sprintf('< %g', bound);
mm_pin(name, str2double(tok{1}{1}), sprintf('abs:%g', bound), shown);
for k = numel(ext):-1:1                      % last first, so earlier extents stay valid
    e = ext{k}(1, :);
    txt = [txt(1:e(1)-1) shown txt(e(2)+1:end)];
end
end
