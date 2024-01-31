In brief, there are 4 inputs for grad-calculation of BH3(4e3o) examples.
1) Normal CI or MCLR                                                                           (Correct)
2) DMRG-MCLR with workaround solution (ground state Hami for excited state)                    (Correct)
3) Same as 2), while M is limited to 1 (i.e. 2 preserved state)                                (Wrong)
4) same as 3), except the ground state Hamiltonian and excited state Hamiltonian are different (Wrong)

The detailed results and some analysis :

1) The BH3_RASSCF_4e3o_CI.input input is used as the reference calculated by normal CI.
   Usage: molcas -f BH3_RASSCF_4e3o_CI.input
   
      energy=     -22.784225 (act-e : -5.83967)
      conf/sym  111     Coeff      Weight
             1  220  0.99624969 0.99251345
             2  2ud  0.00039859 0.00000016
             3  202  -.06223029 0.00387261
             4  u2d  -.00999069 0.00009981
             5  ud2  -.00001634 0.00000000
             6  022  -.05927871 0.00351397

      energy=     -22.587765 (act-e : -5.64321) 
      conf/sym  111     Coeff      Weight
             1  220  -.00869590 0.00007562
             2  2ud  0.00001467 0.00000000
             3  202  0.02192189 0.00048057
             4  u2d  -.99972151 0.99944309
             5  ud2  0.00052708 0.00000028
             6  022  -.00066770 0.00000045

2) The BH3_RASSCF_4e3o.input input is used as the reference calculated by DMRG.
   Usage: molcas -f BH3_RASSCF_4e3o.input MOLCAS_WORKDIR=./LR_inputs_S1_initial
   It can obtain the same result as the normal CI case. 
   The local Hamiltonian (Hami.x.txt in scratch):
  -5.83511   0.00247   0.00002   0.03627   0.00000   0.00000   0.03796
   0.00247  -5.64306   0.00000   0.00025   0.00004   0.00000   0.00786
   0.00002   0.00000  -5.63861  -0.00004  -0.01599   0.00000   0.00011
   0.03627   0.00025  -0.00004  -5.26910   0.00001   0.00000   0.03710
   0.00000   0.00004  -0.01599   0.00001  -5.30784   0.00000  -0.00004
   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
   0.03796   0.00786   0.00011   0.03710  -0.00004   0.00000  -5.26861
   with the lowest eigenvalues are -5.83967 -5.64321, same as normal CI case.
   (all eigenvalues : -5.83967 -5.64321 -5.63938 -5.30706 -5.30587 -5.22714)

3) The BH3_RASSCF_4e3o_M2.input is the testing example for workaround solution (Same code as 2), but the M value change the states) 
   Usage: molcas -f BH3_RASSCF_4e3o_M2.input MOLCAS_WORKDIR=./LR_inputs_S1_initial_M2
   It give the wrong state due to the truncation / state-specific solver
   The local Hamiltonian (Hami.x.txt in scratch): 
  -5.83511   0.03796  -0.00000   0.00000   0.03627
   0.03796  -5.26861  -0.00011   0.00000   0.03710
  -0.00000  -0.00011  -5.30784   0.00000   0.00002
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.03627   0.03710   0.00002   0.00000  -5.26909
   with the lowest eigenvalues are -5.83964 -5.30784
   (all eigenvalues : -5.83964 -5.30784 -5.30594 -5.22722)

4) The BH3_RASSCF_4e3o_M2_no_irreps.input is the testing example for normal QCmaquis, i.e. state-specific
   Usage: molcas -f BH3_RASSCF_4e3o_M2_no_irreps.input MOLCAS_WORKDIR=./LR_inputs_S1_initial_M2_no_irreps
   It give wrong LR results due to the basis is not match
   It can also reflected by the local Hamiltonian (Hami.0.txt and Hami.1.txt, respectively)
   || Hami.0.txt --- ground state
  -5.83511   0.03796  -0.00000   0.00000   0.03627
   0.03796  -5.26861  -0.00011   0.00000   0.03710
  -0.00000  -0.00011  -5.30784   0.00000   0.00002
   0.00000   0.00000   0.00000   0.00000   0.00000
   0.03627   0.03710   0.00002   0.00000  -5.26909
   with the lowest eigenvalues are -5.83964 
   (all eigenvalues : -5.83964 -5.30784 -5.30594 -5.22722)
   || Hami.1.txt --- excited state
  -5.63861   0.00000  -0.00050   0.01599   0.00000   0.00019
   0.00000  -5.64305  -0.00786  -0.00017   0.00000  -0.00025
  -0.00050  -0.00786  -5.26861  -0.00011   0.00000   0.03710
   0.01599  -0.00017  -0.00011  -5.30784   0.00000   0.00002
   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
   0.00019  -0.00025   0.03710   0.00002   0.00000  -5.26909
    with the lowest eigenvalues are -5.64322 
   (all eigenvalues : -5.643225 -5.63938 -5.30708 -5.30585 -5.23167)

   

