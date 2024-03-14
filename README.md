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

       Basis sets, which can be copied from 

   b. GEOM.xyz

       geometry file

   c. 
   
   
