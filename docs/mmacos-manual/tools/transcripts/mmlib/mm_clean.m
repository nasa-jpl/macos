function txt = mm_clean(txt)
%MM_CLEAN  Remove what differs by machine and not by the design from captured output: MATLAB's
%   bracketed warning blocks (graphics hardware, the warning's stack, links) and compute-time lines.
bs = char(8);                    % -batch brackets a warning as '[' BS ... ']' BS
txt = regexprep(char(txt), ['\[' bs '?Warning:.*?\]' bs '?[ \t]*\n'], '');
txt = regexprep(txt, ['\[' bs '?> In .*?\]' bs '?[ \t]*\n'], '');
txt = regexprep(txt, '[^\n]*(Compute time|wall time|elapsed)[^\n]*\n', '');
end
