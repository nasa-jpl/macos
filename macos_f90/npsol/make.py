import os
import sys

# Change into an intermediate directory because gfortran puts the .mod files in the cwd
intermediate = 'gcc/'
home = '../'
path = 'c:\\mingw64\\bin;c:\\python27\\scripts'
compileMacos = True

if len(sys.argv) > 1:
    if sys.argv[1] == 'nomacos':
        compileMacos = False

doLink = True
doCompile = True
    
npsol = {'folder': '../npsol/',
    'compiler': 'gcc-4.8.2 -O2 -fPIC -shared -ffixed-line-length-none -w -c',
    'files': ['blas/daxpy.f', 'blas/dcopy.f', 'blas/ddot.f', 'blas/dgemm.f', 'blas/dgemv.f', 'blas/dger.f', 'blas/dnrm2.f', 'blas/dscal.f', 'blas/dswap.f', 'blas/dsyr.f', 'blas/dtrmm.f', 'blas/dtrmv.f', 'blas/dtrsv.f', 'blas/idamax.f', 'blas/lsame.f', 'blas/xerbla.f', 'lapack/dbdsqr.f', 'lapack/dgebd2.f', 'lapack/dgebrd.f', 'lapack/dgelq2.f', 'lapack/dgelqf.f', 'lapack/dgeqr2.f', 'lapack/dgeqrf.f', 'lapack/dgesvd.f', 'lapack/dlabrd.f', 'lapack/dlacpy.f', 'lapack/dlamch.f', 'lapack/dlange.f', 'lapack/dlapy2.f', 'lapack/dlarf.f', 'lapack/dlarfb.f', 'lapack/dlarfg.f', 'lapack/dlarft.f', 'lapack/dlartg.f', 'lapack/dlas2.f', 'lapack/dlascl.f', 'lapack/dlaset.f', 'lapack/dlasr.f', 'lapack/dlassq.f', 'lapack/dlasv2.f', 'lapack/dorg2r.f', 'lapack/dorgbr.f', 'lapack/dorgl2.f', 'lapack/dorglq.f', 'lapack/dorgqr.f', 'lapack/dorm2r.f', 'lapack/dormbr.f', 'lapack/dorml2.f', 'lapack/dormlq.f', 'lapack/dormqr.f', 'lapack/ilaenv.f', 'blaso.f', 'chsubs.f', 'cmsubs.f', 'lssubs.f', 'mcsubs.f', 'npsubs.f', 'opsubs.f', 'srsubs.f']
    }

macos = {'folder': './',
    'compiler': 'gcc-4.8.2 -O2 -fPIC -shared -ffixed-line-length-none -w -c -DCSMACOS -DNR_FFT -DMSWIN',
    'files': ['kinds.F90', 'constants.F90', 'mathsub.F', 'export_fits_format.F', 'macos_debug.F', 'param_mod.F', 'elt_mod.F', 'src_mod.F', 'cfiles_mod.F', 'usersub.F', 'lohpars_mod.F', 'nn_util.F', 'lsq.F', 'loh_mod.F', 'traceutil_mod.F', 'dopt_mod.F', 'surfsub.F', 'didesub.F', 'macos_mod.F', 'smacos_vars_mod.F', 'smacosutil.F', 'smacosio.F', 'zern_wf.F90', 'macos_init.F', 'elemsub.F', 'funcsub.F', 'pixsub.F', 'sourcsub.F', 'tracesub.F', 'srtrace.F', 'utilsub.F', 'linsub.F', 'propsub.F', 'pgplotdummy.F', 'sunsub.F', 'nls.F', 'macos_ops.F', 'design_optim.F', 'design_cons_optim.F', 'stop_set.F', 'smacos.F']
    }

pymacos = {'folder': '',
    'compiler': 'gfortran-4.8.2 -c -O2 -fPIC -shared -DCSMACOS -xf77-cpp-input -ffree-form -fimplicit-none -ffixed-line-length-none -w -DNR_FFT -DMSWIN',
    'files': ['pymacos.f90']
    }

compileGroups = [npsol]
linkFiles = []

# link = gcc-4.8.2 -shared -o ../pymacos.dll blaso.o cfiles_mod.o chsubs.o cmsubs.o constants.o daxpy.o dbdsqr.o dcopy.o ddot.o design_cons_optim.o design_optim.o dgebd2.o dgebrd.o dgelq2.o dgelqf.o dgemm.o dgemv.o dgeqr2.o dgeqrf.o dger.o dgesvd.o didesub.o dlabrd.o dlacpy.o dlamch.o dlange.o dlapy2.o dlarf.o dlarfb.o dlarfg.o dlarft.o dlartg.o dlas2.o dlascl.o dlaset.o dlasr.o dlassq.o dlasv2.o dnrm2.o dopt_mod.o dorg2r.o dorgbr.o dorgl2.o dorglq.o dorgqr.o dorm2r.o dormbr.o dorml2.o dormlq.o dormqr.o dscal.o dswap.o dsyr.o dtrmm.o dtrmv.o dtrsv.o elemsub.o elt_mod.o export_fits_format.o funcsub.o idamax.o ilaenv.o jwst2.o kinds.o linsub.o lohpars_mod.o loh_mod.o lsame.o lsq.o lssubs.o macos_debug.o macos_init.o macos_mod.o macos_ops.o mathsub.o mcsubs.o nls.o nn_util.o npsubs.o opsubs.o param_mod.o pgplotdummy.o pixsub.o propsub.o pygate.o smacos.o smacosio.o smacosutil.o smacos_vars_mod.o sourcsub.o src_mod.o srsubs.o srtrace.o stop_set.o sunsub.o surfsub.o tracesub.o traceutil_mod.o usersub.o utilsub.o xerbla.o zern_wf.o -lgfortran -Wl,--out-implib,../pymacos.a

if not os.path.exists(intermediate):
    os.mkdir(intermediate)
os.chdir(intermediate)

oldPath = os.environ['PATH']
os.environ['PATH'] = path + ';' + oldPath;

#compile
for c in compileGroups:
    for f in c['files']:
        name = os.path.abspath( home + c['folder'] + f)
        result = os.path.splitext(os.path.basename(f))[0] + '.o'
        linkFiles.append(result)
        cmd = c['compiler'] + ' -o ' +  result + ' ' +  name        
        if doCompile:
            print(cmd)  
            os.system(cmd)

if doLink:
#link
    cmd = 'ar rcs ./npsol_lib.a '+' '.join(linkFiles)
    print(cmd)
    os.system(cmd)

os.chdir(home)

os.environ['PATH'] = oldPath;
