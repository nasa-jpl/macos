OS   = $(shell uname -s)

ifeq ($(OS),Darwin)
    ARCH := $(shell uname -p)
else
    ARCH := $(shell uname -i)
endif

ifeq ($(OS),Linux)
  # Autodetect 32-bit v. 64-bit
  ifeq ($(ARCH),x86_64)
    FC=lf95
    FFLAGS = -O -fs -w -Kfast -Kpic -tpp --verbose --ap --wide --fast --wo 
  else
    FC=lf95
    FFLAGS = -O -fs -w -Kfast -tpp -fPIC
  endif
endif

ifeq ($(OS),Darwin)
  FC     = gfortran
  FFLAGS = -O2
endif

OBJ= \
   ddot.o \
   dswap.o \
   dtrmv.o \
   dtrsv.o \
   dscal.o \
   daxpy.o \
   dgemv.o \
   dnrm2.o \
   dcopy.o \
   idamax.o \
   dger.o \
   dsyr.o \
   xerbla.o \
   dgemm.o \
   lsame.o \
   dtrmm.o 

blas.o: $(OBJ)
	ld -r -o $@ $(OBJ)

clean:
	rm -f *.o

# begin fortran depends
daxpy.o: daxpy.f
dcopy.o: dcopy.f
ddot.o: ddot.f
dgemm.o: dgemm.f
dgemv.o: dgemv.f
dger.o: dger.f
dnrm2.o: dnrm2.f
dscal.o: dscal.f
dswap.o: dswap.f
dsyr.o: dsyr.f
dtrmm.o: dtrmm.f
dtrmv.o: dtrmv.f
dtrsv.o: dtrsv.f
idamax.o: idamax.f
lsame.o: lsame.f
xerbla.o: xerbla.f
# end fortran depends
