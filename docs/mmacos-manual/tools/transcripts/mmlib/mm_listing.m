function mm_listing(name, file, rows)
%MM_LISTING  Emit a source or output FILE (verbatim; ROWS = [first last] or all) as transcript NAME.
t = splitlines(string(fileread(file)));
if nargin > 2 && ~isempty(rows), t = t(rows(1):min(rows(2), numel(t))); end
mm_transcript(name, strjoin(t, newline));
end
