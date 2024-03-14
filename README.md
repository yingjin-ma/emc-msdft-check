qc-modft : quantum computing - molecular orbital density functional theory

0. Pre-install
   
   cd INT_XC folder, build the libcint.so library using CMake
   
2. Install

   a. edit the Makefile;
   
       change the LIB parameter to the file of libcint.so, e.g. LIB = /home/yingjin/Quantum_Soft/eMC_vMSDFT/vMSDFT/INT_XC/LIB/v1/lib/libcint.so

   b. make

3. Run

   a. cd example

   b. LD_LIBRARY_PATH ../qc-modft
   
       e.g. LD_LIBRARY_PATH=../INT_XC/LIB/v1/lib ../qc-modft
   
4. More about the inputs that needed for qc-modft

   a. *.bas file

       Basis sets, which can be downloaded/copied from https://www.basissetexchange.org/
       Notice : Change the "Download basis set" format to Turbomole, select the specific atoms and basis sets, download/copy it
                then, modify the "Head" and add “$ END” (e.g. refer the existing 6-31g.bas)

   b. GEOM.xyz

       geometry file

   c. GAMESSMO.dat

       Initial MO orbitals, extracted using the GAMESS package, it looks like

       ‘’‘
        0.32531328199999998       0.27430139599999998       0.32531328199999998       0.27430139599999998
        0.12460875500000000        1.6899041100000001      -0.12460875500000000       -1.6899041100000001     
        0.77074206700000003      -0.68687553199999996       0.77074206700000003      -0.68687553199999996     
       -1.1145297000000001        1.3418075000000000        1.1145297000000001       -1.3418075000000000
       ’‘’

       can be obtained by implementing GAMESS input (input sample as below) 

       ‘’‘
       $CONTRL SCFTYP=RHF RUNTYP=ENERGY ISPHER=0 NPRINT = 4 $END
       $SYSTEM timlim=1000 $END
       $BASIS  GBASIS=N31 NGAUSS=6 $END
       $DATA
        check the data
        C1
        H   1    0.00000000    0.0    0.375
        H   1    0.00000000    0.0   -0.375
       $END
       ‘’‘
   
   d. INCAR file

       Input parameter for HFSCF/CASSCF/vMSDFT, it looks like

       '''
       MSDFT-OPT.log
       6-31G
       METHOD  0    !INT,(0)CASSCF,(1)MSDFT           $$$$  qc-modft set to 1 $$$$$$
       NATOM   2
       NAO     4   !INT,number of atomic orbitals
       FROZEN  0    !INT,frozen orbitals
       CLOSED  0    !INT,closed orbitals
       ACTIVE  2    !INT,active orbitals
       ELE     2    !INT,total electrons
       SPIN    1    !INT,2S+1 eg,singlet:1
       STATES  1    !INT,number of target states
       IMAX    50   !INT,the max step for orb opt
       RESTART 0    !INT,(1)restart or (0)not
       OPTORB  0    !INT,(0)orbital opt or (1)not
       IWEIGHT 0    !INT,(1)read from WEIGHTS or (0)not
       iMETHOD 100    !INT,(0)HH (1)XHF (2)J+K (3)J
       CONFs   4    !INT,(0)All DETs (>0) Number of DETs
       iROOT   4    !INT,Number of roots for TQL
       KAPPA   0
       ZETA    0
       '''

    * HF   reference for H2/6-31G
      
          [- RHF state energy ] :  -1.1265450511979187
      
    * DFT  reference for H2/6-31G
      
          [- RHF state energy ] :  -1.1823368655417272    
   

 
             
