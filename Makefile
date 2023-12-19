Program_name=eMC_vMSDFT

F90srcs = MSDFT_MAIN.f90 MSDFT_RDINCAR.f90 MSDFT_COMB.f90 MSDFT_MCCONF.f90 MSDFT_BLWCAS.f90 MSDFT_CHECKFILE.f90 MSDFT_MCFILE.f90\
	  MSDFT_PUTINT.f90 MSDFT_METHOD100.f90 MSDFT_EXCR.f90 MSDFT_INT.f90 MSDFT_MCHAM.f90 MSDFT_UPDATE.f90 MSDFT_1-2RDM.f90\
	  MSDFT_ENG.f90 MSDFT_GENDETMO.f90 MSDFT_TrRAM.f90 MSDFT_TQL.f90 MSDFT_FC.f90 MSDFT_READGEOM.f90\
	  INTXC_MOD.f90 INT_BASIS.f90\
	  XC_MAT.f90 XC_DFT.f90 XC_DFTFC.f90 XC_GRIDGEN.f90 XC_GTOEVAL.f90 XC_RHOAB.f90 XC_KMCF.f90 XC_ZETA.f90\
	  Global_control.f90 Date_time.f90 orbopt_header.f90 Coupling_terms.f90\
	  Matrix.f90 SCFORBOPT.f90 Mathlib.f90 Integrals.f90\
	  Input.f90 Grouptable.f90 HFenergy.f90 PreOneSCF.f90 Orbital_Hessian.f90 Derivaties.f90\
	  NonLinear.f90 Initial_value.f90 Iteration.f90 Iteration_ext.f90 Transform_INT.f90\
	  Run_DMRG.f90 Scf_save.f90 Solver_Hessian.f90 Coupling.f90 Coupled_solver.f90\
	  Derivaties_deltaR.f90 Transform_Nindex.f90 Solver_Augmented_Hessian.f90 Davidson.f90\
	  Operator_update.f90 Solver_WMK.f90 MPSci_update.f90 Print_info.f90 Redundancies.f90\
	  PreDMRGLR.f90 ConDMRGLR.f90 Solver_LR.f90 LagrangeRDMs.f90 Operator_update_CPNR.f90\
	  Coupling_update.f90 Operator_update_LR.f90 Fock_gen.f90 Energy_gen.f90 Gmat_gen.f90\
	  Environment.f90 Y_gen.f90 PmultU.f90 Closed_Shell.f90 Precabb.f90 Precaii.f90 Solver_LR_partition.f90\

F77srcs = XC_LEBEDEV.F XC_LIB.F

F90objs = $(F90srcs:.f90=.o)
F77objs = $(F77srcs:.F=.o)
LIB = libcint.a

Complier1=ifort
Complier2=gfortran

Link=-o

#With intel mkl and GNU compilier
FLIB =  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

#59.72.122.148 with GNU, gfortran only
#FLIB = -L /usr/lib64 -m64 -lblas -llapack

#mac OS with intel mkl
#FLIB = -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

$(Program_name):$(F90objs) $(F77objs)
	$(Complier2) $(Link) $(Program_name) $(F90objs) $(F77objs) $(LIB) $(FLIB)

.SUFFIXES : .F .f90 .o

.f90.o:
	$(Complier2) -O2 -g -c -fallow-argument-mismatch $<

.F.o:
	$(Complier2) -O2 -g -c -fallow-argument-mismatch $<

.PHONY:clean
clean : 
	rm -rf *.o *.mod
