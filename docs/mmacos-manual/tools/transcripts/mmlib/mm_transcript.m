function mm_transcript(name, text)
%MM_TRANSCRIPT  Emit TEXT (char/string, lines split on newline) as the generated listing
%   src/transcripts/NAME.md.  The lines go to stdout tagged '@@MM|T|NAME|'; the regenerator
%   (tools/regen_mmacos_transcripts.py) collects them and writes nothing unless the whole driver
%   ran and every pin held.  A chapter pulls the file in with <!-- transcript: NAME -->.
lines = splitlines(string(text));
if ~isempty(lines) && lines(end) == "", lines(end) = []; end
for k = 1:numel(lines)
    fprintf('@@MM|T|%s|%s\n', name, lines(k));
end
fprintf('@@MM|E|%s\n', name);
end
