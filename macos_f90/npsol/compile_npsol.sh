#!/bin/sh

# Make sure we have the intel compiler.
if [ -z "$intel_version" ]
then
  intel_version=13.1.4
fi

# source intel iccvars file.
source /opt/intel-$intel_version/bin/compilervars.sh intel64

# Compile blas and lapack.
(cd blas; make -f Makefile_Intel clean; make -f Makefile_Intel)
(cd lapack; make -f Makefile_Intel clean; make -f Makefile_Intel)

# Compile npsol
make -f Makefile_Intel clean
make -f Makefile_Intel


