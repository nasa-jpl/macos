function [here, check] = mm_template(section, name)
%MM_TEMPLATE  Copy a template into the scratch directory and run there (the manual's rule: a
%   committed template's own outputs are never touched).  The copy keeps the library layout
%   (<work>/templates/SECTION/NAME, with src, design, sensitivities and mmacos_setup.m LINKED beside
%   it), so the template's '..'-relative paths resolve.  Returns HERE, the copy (cwd, on the path),
%   and CHECK, a function the driver calls after its runs: it errors if `git status` of
%   MACOS_resources (untracked files included) moved, i.e. if a run wrote into the library.
root = getenv('MMACOS_ROOT');  work = getenv('MMMAN_WORK');
res  = fileparts(root);
gs   = @() evalc(sprintf('!git -C %s status --porcelain --untracked-files=all', res));
before = gs();
dst = fullfile(work, 'templates', section);  mkdir(dst);
copyfile(fullfile(root, 'templates', section, name), fullfile(dst, name));
for top = {'src', 'design', 'sensitivities', 'mmacos_setup.m'}
    system(sprintf('ln -s %s %s', fullfile(root, top{1}), fullfile(work, top{1})));
end
addpath(fullfile(root, 'src'));
here = fullfile(dst, name);  addpath(here);  cd(here);
check = @() assert(strcmp(gs(), before), 'mm_template:library', ...
    '%s: a run changed MACOS_resources (git status moved)', name);
end
