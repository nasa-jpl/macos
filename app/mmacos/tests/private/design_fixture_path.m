function p = design_fixture_path(name)
%DESIGN_FIXTURE_PATH  Absolute path to a shared optical-design fixture.
%   The telescope-design reference fixtures live one level up, shared by
%   mmacos and pymacos, at MACOS_resources/optical_design/fixtures/ (see
%   that dir's README + OPTICAL_DESIGN_AGENT_GUIDE.md).  Use this helper
%   so the physical layout can change without touching every test.
%
%   Example:
%       fx = jsondecode(fileread( ...
%               design_fixture_path('telescope_design_fixtures.json')));
    arguments
        name (1,:) char
    end
%   Roots are anchored to THIS clone first (repo-relative), with the
%   canonical ~/dev/MACOS_resources checkout as a fallback.
    roots = { fullfile(repo_root_(), 'optical_design', 'fixtures') };   % this repository only
    for i = 1:numel(roots)
        p = fullfile(roots{i}, name);
        if exist(p, 'file'), return; end
    end
    error('design_fixture_path:notFound', ...
        'optical-design fixture not found: %s (searched %s)', ...
        name, strjoin(roots, ', '));
end
