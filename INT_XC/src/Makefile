Program_name=test_intxc

F90srcs = program.f90 INTXC_MOD.f90 INT_BASIS.f90\
	  XC_MAT.f90 XC_DFT.f90 XC_DFTFC.f90 XC_GRIDGEN.f90 XC_GTOEVAL.f90 XC_RHOAB.f90\

F77srcs = XC_LEBEDEV.F XC_LIB.F

F90objs = $(F90srcs:.f90=.o)
F77objs = $(F77srcs:.F=.o)
LIB = INT_LIBCINT.a



Complier2=ifort
Complier1=gfortran

Link=-o

#192.9.200.120 with intel mkl
#FLIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm  -lblas

#59.72.122.148 with GNU, gfortran only
FLIB = -L /usr/lib64 -m64 -lblas -llapack

#mac OS with intel mkl
#FLIB = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

$(Program_name):$(F90objs) $(F77objs)
	$(Complier2) $(Link) $(Program_name) $(F90objs) $(F77objs) $(LIB) $(FLIB)

.SUFFIXES : .F .f90 .o

.f90.o:
	$(Complier2) -O2 -g -c $<

.F.o:
	$(Complier2) -O2 -g -c $<

.PHONY:clean
clean : 
	rm -rf *.o *.mod
