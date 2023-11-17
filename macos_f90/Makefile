#########################################################################
#  	Top-Level Makefile for MACOS
#
#       Contact:     John Z. Lou
#       Address:     JPL 306-472, 4800 Oak Grove Dr., Pasadena CA 91109
#       Email:       John.Z.Lou@jpl.nasa.gov
#       Phone:       (818)354-4870
#########################################################################

PWD := $(shell pwd)

# Our source files are in the current directory
export SRC_DIR := "$(PWD)"

# The base of our build location
export BUILD_LOC := $(SRC_DIR)

# Current Revision:
export SVN_REV := $(shell svnversion)

# Determine our platform and 32-bit vs. 64-bit
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
    PLATFORM := $(shell uname -p)
    ARCH     := $(OS)-$(PLATFORM)
    Make.T   := Make-$(ARCH).T
    FFT_LIB  := NR_FFT
else
    PLATFORM := $(shell uname -i)
    ifeq ($(PLATFORM),x86_64)
	OS      := LNXOS64
	ARCH    := Linux-x86_64
	Make.T  := Make-$(ARCH).T
	FFT_LIB := NR_FFT
	# FFT_LIB = INTEL_FFT_LIB
	# FFT_LIB = FFTW_LIB
    endif

    ifeq ($(PLATFORM),i386)
	OS      := LNXOS32
	ARCH    := Linux-i386
	Make.T  := Make-$(ARCH).T
	FFT_LIB := NR_FFT
	# FFT_LIB = INTEL_FFT_LIB
	# FFT_LIB = FFTW_LIB
    endif

    ifeq ($(findstring SUNW,$(PLATFORM)),SUNW)
	OS     := SUNOS
	ARCH   := SUNW-$(shell uname -p)
	Make.T := Make-Sun.T
	MACOSDIR := /home/jzlou/macos-vY/macos-v3.2b
    endif
endif

export M_OBJS_DIR    := MACOS_OBJS/$(ARCH)
export SM_OBJS_DIR   := SMACOS_OBJS/$(ARCH)
export MV3D_OBJS_DIR := MACOS_V3D_OBJS/$(ARCH)
export smacoslib     := smacos_lib

################################################################################
# Set defaults for compilers and libraries.
# These settings can be overriden as needed in the architecture-specific Makefile.

#### Compilers:
#export FC = f95
export FC = ifort -g -fPIC
#export FC = gfortran -fPIC
export CC = gcc -fPIC
#export CC = icc -fPIC
#export CC = mcc -fPIC

## Linker:
export LD = $(FC)

#### Compiler/Linker Flags:

export FFLAGS += -fs -w -tpp
# For Optimization:
#FFLAGS += -O2
# For Debugging:
FFLAGS += -g -fpen

DEFINES += -D$(OS) -D$(FFT_LIB)
# For Threading support:
DEFINES += -DMT0
# For some extra debugging prints:
#DEFINES += -DDEBUG

#export LDFLAGS += -static -O3 --openmp
export LDFLAGS += -O3 --openmp
# For Debugging
#export LDFLAGS = -g

# Look for source files in the base directory
export VPATH = ../..

# Default names?
export macosexe = macos
export macosexe_v3d = macosv3d
export FLINT = flint
export FMEX = fmex

################################################################################

export MOD_SRCS = kinds.f90 constants.f90 macos_debug.F param_mod.F elt_mod.F \
           src_mod.F cfiles_mod.F lohpars_mod.F nn_util.F \
           traceutil_mod.F macos_mod.F loh_mod.F dopt_mod.F

export MOD_OBJS = kinds.o constants.o macos_debug.o param_mod.o elt_mod.o src_mod.o cfiles_mod.o \
           lohpars_mod.o traceutil_mod.o macos_mod.o nn_util.o loh_mod.o \
           dopt_mod.o

export M_OBJS = macos_vars_mod.o mathsub.o macosio.o macos_init.o \
        elemsub.o surfsub.o funcsub.o \
        pixsub.o sourcsub.o tracesub.o srtrace.o utilsub.o \
        linsub.o propsub.o usersub.o matloh.o didesub.o \
        dftsub.o pgplotsub.o sunsub.o lsq.o nls.o macos_ops.o\
        zern_wf.o design_optim.o design_cons_optim.o stop_set.o psesub.o f90prt.o

export SM_OBJS = smacos_vars_mod.o mathsub.o smacosutil.o smacosio.o \
        macos_init.o elemsub.o surfsub.o funcsub.o \
        pixsub.o sourcsub.o tracesub.o srtrace.o utilsub.o \
        linsub.o propsub.o usersub.o matloh.o didesub.o \
        dftsub.o pgplotdummy.o sunsub.o lsq.o nls.o \
        macos_ops.o zern_wf.o design_optim.o design_cons_optim.o stop_set.o psesub.o

export C_SRCS = utilsub_c.c mhist.c splicubi2.c
export C_OBJS = utilsub_c.o mhist.o splicubi2.o
export C2_SRCS = utilsub_c.c fk.c splicubi2.c
export C2_OBJS = utilsub_c.o fk.o splicubi2.o

################################################################################

default: macos smacos

# Make the object directories:
make-obj-dirs:
	@- mkdir -p "$(M_OBJS_DIR)" "$(MV3D_OBJS_DIR)" "$(SM_OBJS_DIR)" 2>/dev/null
	cp $(Make.T) $(M_OBJS_DIR)/Makefile

# For testing what values are defined.
test-values: make-obj-dirs
	# $(PLATFORM)
	# $(OS)
	# $(Make.T)
	# $(ARCH)
	# $(MACOSDIR)
	# $(M_OBJS_DIR)
	# $(MV3D_OBJS_DIR)
	# $(SM_OBJS_DIR)
	# $(PWD)

# Choice of a MACOS FFT routine
# -DSUNOS ==> use FFT from Sun performance libaray (by default)
# -DUSE_FFTW ==> use FFT from a fftw libaray
# -DNR_FFT ==> use FFT from Numerical Recipes
# When -DSUNOS and -DUSE_FFTW both defined, use fftw

macos: $(M_OBJS_DIR)_Makefile rm_macos make-obj-dirs
	(cd $(M_OBJS_DIR); $(MAKE) DEFINES="$(DEFINES) -DCMACOS" macos)
	cp $(M_OBJS_DIR)/macos "$(PWD)/macos-$(ARCH)"
	/bin/rm -f macos; ln -s "$(PWD)/macos-$(ARCH)" macos
#	cp macos /home/jzlou/Pub/macoslnx

macosv3d: $(MV3D_OBJS_DIR)_Makefile rm_macos rm_macosv3d make-obj-dirs
	touch macos_cmd_loop.inc
	(cd $(MV3D_OBJS_DIR); $(MAKE) DEFINES="$(DEFINES) -DCMACOS -DVIS3D_DEMO" macosv3d)
	cp $(M_OBJS_DIR)/macosv3d "$(PWD)/macosv3d-$(ARCH)"

smacos: $(SM_OBJS_DIR)_Makefile rm_smacos make-obj-dirs
	(cd $(SM_OBJS_DIR); $(MAKE) DEFINES="$(DEFINES) -DCSMACOS" $(smacoslib).o)

smacosdvr: $(SM_OBJS_DIR)_Makefile rm_smacos_dvr
	(cd $(SM_OBJS_DIR); $(MAKE) DEFINES="$(DEFINES) -DCSMACOS" smacos_dvr)

$(M_OBJS_DIR)_Makefile: $(Make.T)
	if [ ! -d $(M_OBJS_DIR) ]; then mkdir -p $(M_OBJS_DIR); fi;
	cp $(Make.T) $(M_OBJS_DIR)/Makefile

$(MV3D_OBJS_DIR)_Makefile: $(Make.T)
	if [ ! -d $(MV3D_OBJS_DIR) ]; then mkdir -p $(MV3D_OBJS_DIR); fi;
	cp $(Make.T) $(MV3D_OBJS_DIR)/Makefile

$(SM_OBJS_DIR)_Makefile: $(Make.T)
	if [ ! -d $(SM_OBJS_DIR) ]; then mkdir -p $(SM_OBJS_DIR); fi;
	cp $(Make.T) $(SM_OBJS_DIR)/Makefile

rm_macos:
rm_macosv3d:
	rm -f macos-$(ARCH) $(M_OBJS_DIR)/macos macosv3d-$(ARCH) $(M_OBJS_DIRS)/macosv3d

rm_smacos:
	rm -f $(SM_OBJS_DIR)/$(smacoslib).o

rm_smacos_dvr:
	rm -f $(SM_OBJS_DIR)/smacos_dvr

clean-smacos:
	-rm -rf $(SM_OBJS_DIR)

clean-macos::
	-rm -f "macos-$(ARCH)"
	-rm -rf $(M_OBJS_DIR)
	-rm -f macos

clean-macosv3d::
	-rm -f macosv3d-$(ARCH)
	-rm -rf $(MV3D_OBJS_DIR)

clean:: clean-smacos clean-macos clean-macosv3d

cleanall:: clean-smacos clean-macos clean-macosv3d
	-rm -rf MACOS_OBJS SMACOS_OBJS MACOS_V3D_OBJS

move::
	-cp -f $(M_OBJS_DIR)/macos ../bin
	-mv -f $(SM_OBJS_DIR)/$(smacoslib).o ../bin

.PHONY: make-obj-dirs test-values clean clean-macos clean-macosv3d move
