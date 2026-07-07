import os
import sys
import subprocess
from pathlib import Path

# Change into an intermediate directory because gfortran puts the .mod files in the cwd
intermediate = Path('gcc/')
home = Path('../')
path = r'c:\mingw64\bin;c:\python27\scripts'
compileMacos = True

if len(sys.argv) > 1:
    if sys.argv[1] == 'nomacos':
        compileMacos = False

doLink = True
doCompile = True

npsol = {
    'folder': '../npsol/',
    'compiler': ['gcc-4.8.2', '-O2', '-fPIC', '-shared', '-ffixed-line-length-none', '-w', '-c'],
    'files': ['blas/daxpy.f', 'blas/dcopy.f', '...']  # same file list as before
}

macos = {
    'folder': './',
    'compiler': ['gcc-4.8.2', '-O2', '-fPIC', '-shared', '-ffixed-line-length-none', '-w', '-c',
                 '-DCSMACOS', '-DNR_FFT', '-DMSWIN'],
    'files': ['kinds.F90', 'constants.F90', '...']  # same file list as before
}

pymacos = {
    'folder': '',
    'compiler': ['gfortran-4.8.2', '-c', '-O2', '-fPIC', '-shared', '-DCSMACOS',
                 '-xf77-cpp-input', '-ffree-form', '-fimplicit-none',
                 '-ffixed-line-length-none', '-w', '-DNR_FFT', '-DMSWIN'],
    'files': ['pymacos.f90']
}

compileGroups = [npsol]
linkFiles = []

# Use pathlib instead of os.mkdir
intermediate.mkdir(parents=True, exist_ok=True)
os.chdir(intermediate)

# Build an explicit env dict rather than mutating os.environ in place
build_env = os.environ.copy()
build_env['PATH'] = path + ';' + build_env['PATH']

# Compile
for c in compileGroups:
    for f in c['files']:
        src = (home / c['folder'] / f).resolve()
        result = Path(f).stem + '.o'
        linkFiles.append(result)

        # List form — no shell, no injection risk
        cmd = c['compiler'] + ['-o', result, str(src)]

        if doCompile:
            print(' '.join(cmd))
            subprocess.run(cmd, check=True, env=build_env)

if doLink:
    # ar is also a subprocess call — keep it as a list
    cmd = ['ar', 'rcs', './npsol_lib.a'] + linkFiles
    print(' '.join(cmd))
    subprocess.run(cmd, check=True, env=build_env)

os.chdir(home)