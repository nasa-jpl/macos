PROJ_DIR := $(shell pwd)/..
NPSOL_DIR := $(PROJ_DIR)/npsol
MACOS_DIR := $(PROJ_DIR)/macos_f90
ARCH := Linux-x86_64
export SMACOS_OBJS := ./gcc
export SMACOS_LIBS := $(SMACOS_OBJS)/smacos_lib.a
PYTHON := python3.3
PYTHON_CMD :=
FFLAGS := 
FC_TAG := gfortran
export FC := gfortran-4.8.2

ifeq ($(FC_TAG),gfortran)
	FCOMPILER := gfortran
	F90EXEC := $(FC)
endif
ifeq ($(FC_TAG),ifort)
	FCOMPILER := intelem
	F90EXEC := $(FC)
endif

export PROJ_SVN_REV := $(shell svn info | awk '/Rev:/{print $$4}')
export PROJ_DATE := $(shell date '+%Y-%m-%d')

all: gfortran 

clean:
	rm -f *.so *.pyf *.mexa64 pymacos_mod.mod

distclean: clean
	$(MAKE) -w -C $(NPSOL_DIR)/blas -f Makefile_Intel clean
	$(MAKE) -w -C $(NPSOL_DIR)/lapack -f Makefile_Intel clean
	$(MAKE) -w -C $(NPSOL_DIR) -f Makefile_Intel clean NPDIR=$(NPSOL_DIR)
	$(MAKE) -w -C $(MACOS_DIR) -f Makefile_Intel clean PWD=$(MACOS_DIR)
	rm -rf gcc

.PHONY: all clean distclean npsol macos_f90 interface f2py

npsol:
ifneq ($(PYTHON_CMD),nomacos)
	$(MAKE) -w -C $(NPSOL_DIR)/blas -f Makefile_Intel
	$(MAKE) -w -C $(NPSOL_DIR)/lapack -f Makefile_Intel
	$(MAKE) -w -C $(NPSOL_DIR) -f Makefile_Intel NPDIR=$(NPSOL_DIR)
endif

macos_f90: npsol
ifneq ($(PYTHON_CMD),nomacos)
	$(MAKE) -w -C $(MACOS_DIR) -f Makefile_Intel PWD=$(MACOS_DIR)
endif

gfortran: make.py
	$(PYTHON) make.py $(PYTHON_CMD)

ifort: npsol
