function mm_pin(name, value, rule, shown)
%MM_PIN  A number the chapter quotes, checked against tools/pins/<chapter>.tsv on every regen.
%   RULE 'sig5' (a signal: five significant figures, the cross-architecture rule), 'abs:<bound>'
%   (a numerical zero: |value| below the bound) or 'printed' (a number the template prints at its own
%   precision: SHOWN is compared by the listing rule, one unit in its last printed digit).  SHOWN is the text the prose prints (%pin:NAME% in
%   a chapter becomes SHOWN at build time -- tools/pins.lua); default: VALUE at five figures.
%   Re-pins (Dave 2026-10-28, PLAN_CONSOLIDATION 1.2 (zzz)): a move past five figures with a KNOWN cause
%   may be re-pinned by the lane that found it, and is reported (cause, old -> new, relative move);
%   an unexplained move is a finding, not a re-pin.
if nargin < 4 || isempty(shown), shown = sprintf('%.5g', value); end
fprintf('@@MM|P|%s|%.17g|%s|%s\n', name, value, rule, shown);
end
