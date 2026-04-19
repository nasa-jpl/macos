
#--------------------------------------------------------------------------
# Compile GMI
#
# Note:
#
# 1. User must change matlab_version above, or leave it as is if matlab
#    has been setup as /usr/local/bin/matlab
#
#-------------------------------------------------------------------------

cd ~/dev/MACOS_resources/GMI
export macossrc_dir=~/dev/macos/macos_f90

# WAS:
## eval $(ifx -C -traceback -fstack-protector -c  -I$macossrc_dir -I$matlab_version/extern/include -I$matlab_version/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMI.F")

## eval $(ifx -C -traceback -fstack-protector -c  -I$macossrc_dir -I$matlab_version/extern/include -I$matlab_version/simulink/include -nologo -fpic -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer -D__amd64 -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  -DGMI_SVN_REV="''" -DGMI_DATE="'2024-01-18'"  -DMX_COMPAT_32 -O2 -xHOST  "GMIG.F")

## eval $(ifx -C -traceback -fstack-protector -O -shared-intel -shared -Wl,--version-script,$matlab_version/extern/lib/glnxa64/fexport.map  -Wl,--no-undefined -o  "GMI.mexa64"  GMI.o GMIG.o   -Wl,-rpath-link,$matlab_version/bin/glnxa64 -L$matlab_version/bin/glnxa64  -l:libmx.so -l:libmex.so -lmat -L/opt/intel-14.0.0/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64 -lirc -lm -lstdc++  $macossrc_dir/SMACOS_OBJS/Linux-x86_64/smacos_lib.a)

# IS:
#  -shared-intel  -L/opt/intel/oneapi/compiler/latest/lib/compiler/msan -l:libmsan.so\
ifx -check all -traceback -fstack-protector -c  \
  -I$macossrc_dir -I/usr/local/MATLAB/R2025b/extern/include \
  -nologo -fPIC -fpp -132 -gen-interfaces -fp-model strict -fno-omit-frame-pointer \
  -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  \
  -DMX_COMPAT_32 -O2 -xHOST  "GMI.F"

ifx -check all -traceback -fstack-protector -c  \
  -I$macossrc_dir -I/usr/local/MATLAB/R2025b/extern/include \
  -I/usr/local/MATLAB/R2025b/extern/include/fintrf.h \
  -nologo -fPIC -fpp -132 -gen-interfaces -fp-model strict  -fno-omit-frame-pointer \
  -module $macossrc_dir/SMACOS_OBJS/Linux-x86_64  \
  -DMX_COMPAT_32 -O2 -xHOST  "GMIG.F"

ifx -check all -traceback -fstack-protector -O -shared-intel -shared \
  -Wl,--version-script,/usr/local/MATLAB/R2025b/extern/lib/glnxa64/fexport.map   \
  -Wl,--no-undefined -o  "GMI.mexa64"  GMI.o GMIG.o   \
  -Wl,-rpath-link,/usr/local/MATLAB/R2025b/bin/glnxa64 \
  -L/usr/local/MATLAB/R2025b/bin/glnxa64  -l:libmx.so -l:libmex.so -lmat \
  -L/opt/intel/oneapi/compiler/latest/lib \
  -Wl,-rpath,/opt/intel/oneapi/compiler/latest/lib \
  -Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib \
  -Wl,--disable-new-dtags \
  -L/usr/lib/x86_64-linux-gnu/libstdc++.so -lirc -lm -lstdc++  \
  $macossrc_dir/SMACOS_OBJS/Linux-x86_64/smacos_lib.a

cd ~/dev/macos

#--------------------------------------------------------------------------
# Requirements to Compile on a Linux System 
#--------------------------------------------------------------------------
# Luis used following versions of compilers for compiling macos on "marchen":
# 1. bash (run on bash, any version)
# 2. Fortran Intel Compiler oneAPI (2024 or earlier down to 2013) ifort or ifx (ifx early 2024 won't compile GMI)
#    User must install oneAPI HPC Toolbox to compile Fortran and Intel oneAPIBase Toolkit. 
# 3. GCC (4.2 or later, used to compile supporting libraries)
# 4. Matlab 2013a to Latest (2023b, GMI will compile with any of these version with either ifort or ifx 
#    with the exception of the ifx early 2024 Fortran Intel Compiler


