function mm_pin(name, value, rule, shown)
%MM_PIN  A number the chapter quotes, checked against tools/pins/<chapter>.tsv on every regen.
%   RULE 'sig5' (a signal: five significant figures, the cross-architecture rule) or 'abs:<bound>'
%   (a numerical zero: |value| below the bound).  SHOWN is the text the prose prints (%pin:NAME% in
%   a chapter becomes SHOWN at build time -- tools/pins.lua); default: VALUE at five figures.
if nargin < 4 || isempty(shown), shown = sprintf('%.5g', value); end
fprintf('@@MM|P|%s|%.17g|%s|%s\n', name, value, rule, shown);
end
