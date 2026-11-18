classdef tDataDir < matlab.unittest.TestCase
%TDATADIR  The engine finds its data files with no MACOS_HOME (issue #75; twin of pymacos
%   tests/test_data_dir.py).
%   The build records its own data directory (BUILD.md, "Where the engine finds its data files"), so
%   macos.init must work from a directory that holds no macos_param.txt -- with MACOS_HOME unset, and
%   with MACOS_HOME naming a directory that does not exist (a stale value in a shell start-up file:
%   the engine warns and falls through).
%
%   THE LEG runs in a matlab -batch SUBPROCESS: a failed lookup STOPS the engine, and inside the mex
%   that takes MATLAB with it.  It clears MACOS_HOME, inits, sets a stale value, and inits again.
%   Before the fix the built-in directory was dead (a macro name inside quotes) and the first init
%   ended the subprocess: no INIT_OK 1.
%   The leg clears the variable ITSELF, inside MATLAB: removing it from the child's environment is
%   not enough, because matlab -batch runs the user's startup.m first, and a startup.m that sets
%   MACOS_HOME (the documented Mac set-up before issue #75) put it back (CCMac's Mac run, 2026-11-17).
%   An empty value is the same as unset to the engine.
    properties (Constant)
        ModelSize = 128
    end
    methods (Test)
        function test_init_with_macos_home_unset_then_stale(tc)
            wd = tempname;  mkdir(wd);  c = onCleanup(@() rmdir(wd, 's'));
            mm = fileparts(fileparts(mfilename('fullpath')));
            scr = fullfile(wd, 'leg.m');  fid = fopen(scr, 'w');
            fprintf(fid, 'run(''%s'');\n', fullfile(mm, 'mmacos_setup.m'));
            fprintf(fid, 'setenv(''MACOS_HOME'', '''');\n');
            fprintf(fid, 'assert(isempty(getenv(''MACOS_HOME'')), ''MACOS_HOME could not be cleared'');\n');
            fprintf(fid, 'assert(~isfile(''macos_param.txt''), ''the working directory must not supply the file'');\n');
            fprintf(fid, 'macos.init(%d);  fprintf(''INIT_OK 1\\n'');\n', tc.ModelSize);
            fprintf(fid, 'setenv(''MACOS_HOME'', fullfile(pwd, ''no_such_directory''));\n');
            fprintf(fid, 'macos.init(%d);  fprintf(''INIT_OK 2\\n'');\n', tc.ModelSize);
            fclose(fid);
            [st, out] = system(sprintf('cd %s && env -u MACOS_HOME "%s" -batch "run(''%s'')" 2>&1', wd, fullfile(matlabroot, 'bin', 'matlab'), scr));
            tc.verifyEqual(st, 0, sprintf('init must not end MATLAB (pre-fix: no parameter file found)\n%s', out));
            tc.verifyTrue(contains(out, 'INIT_OK 1'), sprintf('init with MACOS_HOME unset\n%s', out));
            tc.verifyTrue(contains(out, 'INIT_OK 2'), sprintf('init with a stale MACOS_HOME\n%s', out));
        end
    end
end
