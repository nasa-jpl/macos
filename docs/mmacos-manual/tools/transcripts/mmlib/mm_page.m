function mm_page(name)
%MM_PAGE  Run one template for its Part II page, from its spec tools/part2/specs/<name>.json.
%   The template's folder is copied into the scratch directory (MMMAN_WORK) under its own relative
%   layout (templates/<section>/<name>), with any sibling templates it reads (spec.siblings), and the
%   library's top level (src, design, sensitivities, mmacos_setup.m) LINKED beside it so the
%   template's '..'-relative paths resolve.  A template that climbs past mmacos to MACOS_resources
%   (spec.res_links, e.g. ["segmirmaker"]) is copied under <work>/mmacos instead, with those
%   resources linked at <work>/<name>; its printed paths still read <work>/templates/...  The entry (spec.entry, a MATLAB statement) runs in the
%   copy, unchanged.  The library is checked before and after: if `git status` of MACOS_resources
%   (untracked files included) changed, the page is refused -- a template that writes outside its
%   own folder would have written into the library.
%
%   Emits, through the gated helpers: p2_<name>_knobs (a listing of spec.knobs.file, from the line
%   matching spec.knobs.from through the first line matching spec.knobs.to), p2_<name>_run (the run's
%   own printout, cleaned, from spec.run.from through spec.run.to), p2_<name>_fig (spec.figure),
%   and one pin per spec.pins entry {name, pattern (one token), rule, format}.
spec = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'part2', 'specs', ...
                                    [name '.json'])));
root = getenv('MMACOS_ROOT');  work = getenv('MMMAN_WORK');
res  = fileparts(root);                                         % MACOS_resources
gs   = @() evalc(sprintf('!git -C %s status --porcelain --untracked-files=all', res));
before = gs();

lib = work;                                                     % the copy's mmacos root
lnk = {};  if isfield(spec, 'res_links'), lnk = cellstr(spec.res_links); end
if ~isempty(lnk)          % a template that finds MACOS_resources by '..': copy under <work>/mmacos,
    lib = fullfile(work, 'mmacos');  mkdir(lib);      % link the resources it reads at <work>/<name>
    for k = 1:numel(lnk)
        system(sprintf('ln -s %s %s', fullfile(res, lnk{k}), fullfile(work, lnk{k})));
    end
end
dst = fullfile(lib, 'templates', spec.section);  mkdir(dst);
copyfile(fullfile(root, 'templates', spec.section, name), fullfile(dst, name));
sib = {};  if isfield(spec, 'siblings'), sib = cellstr(spec.siblings); end
for k = 1:numel(sib)            % a sibling in the same section, or '<section>/<name>' in another
    if contains(sib{k}, '/'), sp = sib{k}; else, sp = [spec.section '/' sib{k}]; end
    to = fullfile(lib, 'templates', sp);  mkdir(fileparts(to));
    copyfile(fullfile(root, 'templates', sp), to);
end
for top = {'src', 'design', 'sensitivities', 'mmacos_setup.m'}
    system(sprintf('ln -s %s %s', fullfile(root, top{1}), fullfile(lib, top{1})));
end
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'src'));
addpath(fullfile(root, 'design', 'runners'));
here = fullfile(dst, name);  addpath(here);  cd(here);

t0 = tic;
try
    txt = mm_clean(strrep(run_entry_(spec.entry), lib, work));  % in a workspace of its own; paths read <work>/templates/...
catch ME                       % a failed run is checked against the library too, then reported
    cd(work);
    if ~strcmp(gs(), before)
        error('mm_page:library', '%s: the run failed (%s) AND changed MACOS_resources (git status moved)', name, ME.message);
    end
    rethrow(ME);
end
fprintf('%s wall time %.0f s\n', name, toc(t0));
cd(work);
if ~strcmp(gs(), before)
    error('mm_page:library', '%s: the run changed MACOS_resources (git status moved) -- refused', name);
end
disp(txt)                                                      % the whole run, in the log

if isfield(spec, 'knobs')
    kl = splitlines(string(fileread(fullfile(here, spec.knobs.file))));
    k0 = find(~cellfun(@isempty, regexp(kl, spec.knobs.from, 'once')), 1);
    k1 = k0 - 1 + find(~cellfun(@isempty, regexp(kl(k0:end), spec.knobs.to, 'once')), 1);
    assert(~isempty(k0) && ~isempty(k1), 'mm_page:knobs', '%s: knobs block not found', name);
    mm_listing(['p2_' name '_knobs'], fullfile(here, spec.knobs.file), [k0 k1]);
end
r0 = regexp(txt, spec.run.from, 'once');                       % to '' = through the end of the printout
if isempty(spec.run.to), r1 = numel(txt) - r0 + 1; else, r1 = regexp(txt(r0:end), spec.run.to, 'end', 'once'); end
assert(~isempty(r0) && ~isempty(r1), 'mm_page:run', '%s: the run excerpt was not found', name);
ex = txt(r0:r0+r1-1);
ex = regexprep(ex, '\d{4}-\d{2}-\d{2} \d{2}:\d{2}(:\d{2})?', '<date>');   % a run's wall-clock stamp is not a result
if isfield(spec, 'masks')                                      % other run-time values: {pattern, with}
    for k = 1:numel(spec.masks), ex = regexprep(ex, spec.masks(k).pattern, spec.masks(k).with); end
end
if isfield(spec, 'zeros')                                     % numerical zeros -> their bounds
    for z = 1:numel(spec.zeros)
        q = spec.zeros(z);
        ex = mm_zero(['p2_' name '_' q.name], ex, q.pattern, q.bound);
    end
end
mm_transcript(['p2_' name '_run'], strtrim(ex));
mm_figure(['p2_' name '_fig'], fullfile(here, spec.figure));
for p = 1:numel(spec.pins)
    q = spec.pins(p);
    t = regexp(txt, q.pattern, 'tokens', 'once');
    assert(~isempty(t), 'mm_page:pin', '%s: no match for pin %s', name, q.name);
    mm_pin(['p2_' name '_' q.name], str2double(t{1}), q.rule, sprintf(q.format, str2double(t{1})));
end
end

function txt = run_entry_(entry)
% The template's entry runs here, so its variables cannot clobber mm_page's.
txt = evalc(entry);
end
