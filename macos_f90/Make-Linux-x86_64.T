##############################################################################
# Makefile template (Linux 64-bit version) - for MACOS version 3.34
##############################################################################

#LOCAL_FFLAGS = $(FFLAGS)
# Override FFLAGS with our compiler settings:
LOCAL_FFLAGS = -g
#LOCAL_FFLAGS = -Kpic -fs -w -tpp --wo
# Wth Optimization
#LOCAL_FFLAGS += -O3
#LOCAL_FFLAGS += -O3 --fast -Kfast
# With Debugging and trace
#LOCAL_FFLAGS += -g --chk --trace --info
#LOCAL_FFLAGS += -g -fpe1
#LOCAL_FFLAGS += -O3

# Readline Library Defines
USE_READLINE = 1
ifeq ($(USE_READLINE),1)
  RL_LIBRARY_VERSION = 8.2
  RL_ROOT = $(SRC_DIR)/readline-8.2
  RL_INC = $(RL_ROOT)
  DEFS = -DHAVE_CONFIG_H
  LOCAL_FFLAGS += -DREADLINE_LIBRARY
  LOCAL_CFLAGS += -I$(RL_INC) $(DEFS) -DREADLINE_LIBRARY -DRL_LIBRARY_VERSION='"$(RL_LIBRARY_VERSION)"'
  READLINE_LIB = $(RL_ROOT)/libreadline.a
  #HISTORY_LIB = $(RL_ROOT)/lnx_libs/64/libhistory.a
  RL_LIBS = $(READLINE_LIB) -ltermcap
endif

# Add in Defines:
LOCAL_FFLAGS += $(DEFINES) -DBUILD_LOC="'$(BUILD_LOC)'"
LOCAL_CFLAGS += $(CFLAGS) $(DEFINES) -DBUILD_LOC="\"$(BUILD_LOC)\""

# Add in SVN_REV:
LOCAL_FFLAGS += -DSVN_REV="'$(SVN_REV)'"
LOCAL_CFLAGS += -DSVN_REV="\"$(SVN_REV)\""

################################################################################
# Libraries for MACOS:

MLIBS = /opt/intel-$(intel_version)/lib/intel64/libimf.a \
        /opt/intel-$(intel_version)/lib/intel64/libifport.a \
        /opt/intel-$(intel_version)/lib/intel64/libifcore.a \
        -L/usr/lib64 -lX11 \
        -L$(SRC_DIR)/pgplot -lpgplot \
        -L$(SRC_DIR)/fits_build -lfits \
        $(RL_LIBS)

NPSOL_ROOT=$(SRC_DIR)/npsol
NPSOL_OBJS=$(NPSOL_ROOT)/blas/blas.o $(NPSOL_ROOT)/lapack/lapack.o \
            $(NPSOL_ROOT)/Linux-x86_64/npsol.o
###############################################################################

# MACOS executable
MACOSOBJS = $(MOD_OBJS) $(C_OBJS) $(M_OBJS) macos.o
$(macosexe): $(MACOSOBJS)
	$(FC) $(LDFLAGS) -o $@ $(NPSOL_OBJS) $(MACOSOBJS) $(MLIBS)


SMDRVLIBS = -L/opt/intel/composer_xe_2011_sp1.8.273/compiler/lib/intel64 -limf -lintlc -lifport -lifcore -lc -lm \
            -L/usr/lib64 -lX11

# SMACOS library
SMACOSOBJS = $(MOD_OBJS) $(SM_OBJS) $(C2_OBJS) smacos.o
$(smacoslib).o: $(SMACOSOBJS) $(NPSOL_OBJS)
	ld -r -o $@ $(NPSOL_OBJS) $(SMACOSOBJS)
	ar -r smacos_lib.a $(NPSOL_OBJS) $(SMACOSOBJS)

# SMACOS driver program (example)
SMACOS_DVR_OBJS = $(SMACOSOBJS) smacos_dvr.o
smacos_dvr: $(SMACOS_DVR_OBJS) $(NPSOL_OBJS)
	$(FC) -o $@ $(SMACOS_DVR_OBJS) $(NPSOL_OBJS) $(SMDRVLIBS)

###############################################################################
# Default compilation rules:

%.o:%.f
	$(FC) $(LOCAL_FFLAGS) -c $<

%.o:%.f90
	$(FC) $(LOCAL_FFLAGS) -c $<

%.o:%.F
	$(FC) $(LOCAL_FFLAGS) -c $<

%.o:%.F90
	$(FC) $(LOCAL_FFLAGS) -c $<

%.o:%.c
	$(CC) $(LOCAL_CFLAGS) -c $<

###############################################################################

depend:
	genmake/fmkdep -srcdir "$(SRC_DIR)" $(ALLSRCS)

clean::
	rm -f *.o *~ *.bak

###############################################################################

# begin depends
kinds.o: kinds.f90
constants.o: constants.f90 kinds.f90
dftsub.o: dftsub.F
macos_debug.o: macos_debug.F
param_mod.o: param_mod.F
elt_mod.o: elt_mod.F param_mod.F
src_mod.o: src_mod.F param_mod.F
cfiles_mod.o: cfiles_mod.F param_mod.F
macos_mod.o: macos_mod.F param_mod.F
lohpars_mod.o: lohpars_mod.F
macos_vars_mod.o: macos_vars_mod.F

macos_init.o: macos_init.F
mathsub.o: mathsub.F
elemsub.o: elemsub.F $(MOD_SRCS)
surfsub.o: surfsub.F $(MOD_SRCS) hzdef.inc fzdef.inc
funcsub.o: funcsub.F $(MOD_SRCS) lnk_pert.inc
didesub.o: didesub.F
pgplotsub.o: pgplotsub.F $(MOD_SRCS) fits1.inc mcolors.inc
sunsub.o: sunsub.F $(MOD_SRCS) fits1.inc
pixsub.o: pixsub.F $(MOD_SRCS)
propsub.o: propsub.F $(MOD_SRCS)
sourcsub.o: sourcsub.F $(MOD_SRCS) ssrcray.inc
utilsub.o: utilsub.F lensarr_indexes.inc collim.inc
linsub.o: linsub.F $(MOD_SRCS)
tracesub.o: tracesub.F $(MOD_SRCS) fpp.inc
srtrace.o: srtrace.F $(MOD_SRCS)
usersub.o: usersub.F $(MOD_SRCS) collim.inc
macos.o: macos.F $(MOD_SRCS) fits1.inc glass.inc
macosio.o: macosio.F $(MOD_SRCS) prop_defs.inc msmacosio.inc iosub.inc viewdata_formats.inc
smacosutil.o: smacosutil.F
smacosio.o: smacosio.F $(MOD_SRCS) prop_defs.inc msmacosio.inc iosub.inc viewdata_formats.inc
# matloh.o: matloh.F
pgplotdummy.o: pgplotdummy.F
mhist.o : mhist.c

macos_inc = macos_cmd_loop.inc macos_formats.inc macos_glass.inc \
            fits1.inc glass.inc fits1.inc macos_init.inc efl.inc \
            setbeam.inc
smacos_inc = macos_cmd_loop.inc macos_formats.inc smacos_glass.inc \
            fits1.inc glass.inc fits1.inc macos_init.inc setbeam.inc

stop_set.o: stop_set.F $(MOD_SRCS) design_optim.F

macos.o: macos.F macos_vars_mod.F $(MOD_SRCS) $(macos_inc)
smacos.o: smacos.F smacos_vars_mod.F $(MOD_SRCS) $(smacos_inc)

nn_util.o: nn_util.F
lsq.o: lsq.F
nls.o: nls.F
design_optim.o: design_optim.F smacos_compute.inc
design_cons_optim.o: design_cons_optim.F smacos_compute.inc setbeam_vars.inc setbeam.inc
zern_wf.o: zern_wf.F

utilsub_c.o: utilsub_c.c
splicubi2.o: splicubi2.c

# end depends

###############################################################################

