function mm_figure(name, file)
%MM_FIGURE  Hand the PNG FILE to the regenerator as src/media/NAME.png (copied only on success).
if ~isfile(file), error('mm_figure:missing', 'figure %s: %s was not written', name, file); end
fprintf('@@MM|F|%s|%s\n', name, char(java.io.File(file).getAbsolutePath()));
end
