function r = repo_root_()
%REPO_ROOT_  The root of this repository, found from this file's own location
%   (app/mmacos/tests/private), so a test reads the checkout it lives in whatever
%   the clone is called, wherever it is, and whatever MACOS_HOME holds.
r = fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))));
end
