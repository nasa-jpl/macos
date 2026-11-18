function h = macos_home_()
%MACOS_HOME_  The MACOS_HOME a nested `matlab -batch` is given for macos_param.txt: the parent's, or
%   this checkout's own macos_f90 directory (a system() from MATLAB on macOS drops the parent
%   environment, so the nested engine init found no parameter file; CCMac's Mac gate 2026-10-19).
h = getenv('MACOS_HOME');
if isempty(h)
    h = fullfile(repo_root_(), 'macos_f90');
end
end
