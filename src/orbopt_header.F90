module orbopt_header

use orbopt_version
use hostname

implicit none

public print_orbopt_header

contains

subroutine print_orbopt_header()

!"====================================================================="
!"       Orbital optimization program, version 0.2 (alpha)             "
!  Main references are 
!        H.-J. Werner and W. Meyer,      J. Chem. Phys, 73, 2342 (1980)
!        H.-J. Werner and P. J. Knowles, J. Chem. Phys. 82, 5053 (1985)
!        P.G. Szalay, T. Muller, G. Gidofalvi, H. Lischka, 
!        and R. Shepard,                    Chem. Rev. 112, 108  (2012)
!"---------------------------------------------------------------------"
!" Code written by                                                     "
!                              Yingjin Ma                              "
!"                                                                     "
!"====================================================================="
! Original version        around 2009
!                         updated 2013 for combining Block 
!                         updated 2015 for combining QCmaquis 
!                         updated 2016 for Second-order SCF 
!                         updated 2017 for Linear-response 
! =====================================================================
! a prototype code for Second-order DMRG-SCF 
!                    & DMRG SA-gradients (Hessian matrix part). 

! Steps for this prototype code should be:
! 1) combine Maquis                                          (done) 
! 2) refine the solver, i.e. write a new Hx=g solver         (done)
! 3) write the remaining part of Hessian                     (done)
! 4) import RDM derivatives                                  (done)
! 5) check if it is work                                     (checked) 
! *  improve the integrals part                              (improved)

! to Stefan : Done for orb-opt, but unacceptable slow (40+ iteratons) ??!! (although convergenced results are correct) 
!             I did't check it further due to I turned to SA-gradients    
! 6) closed shell approximation                              (done)

! to Stefan : This part not yet done (code done, but have error), which will deduce the cost of transformation during micro-iteratons
! 7) 2-index trans. split into U * matrix, matrix from 1-index transformation. (to do) 
! 8) Re-code the RDM-deri part (to do)
! 9) AVX improvement (to do)

! Current solvers for SCF
!   1) Werner-Meyer-Knowles with or without coupling
!   2) (Linear) Newton-Raphson
!   3) Augmented Hessian Newton-Raphson
!   4) Step-restricted Augmented Hessian
!   5) Werner-Meyer's finite difference 
!   6) direct inversion
! PS: Only the WMK is well tested and designed
!     Davidson or Full diagionlization is used in 1-4

! Current solvers for LR
!   1) (preconditioner) conjugate gradient  

! If it goes well, shift to MOLCAS

print '(/a)', ' ---------------------------------------------------------------------- '
print '(a/)',' This program is a prototype code for                                    '
print '(a/)',' second-order DMRG-SCF                                                   '
print '(a/)',' and                                                                     '
print '(a//)',' DMRG state-averaged gradients (Lagrange/MCLR part)                     '
print '(a )',' Main author:                                                            '
print '(a )',' ------------                                                            '
print '(a/)','             Y. Ma     (ETH Zurich, U. Nanjing, and CNIC)                '
print '(a )',' Contributing author:                                                    '
print '(a )',' --------------------                                                    '
print '(a/)','             S. Knecht (ETH Zurich)                                      '
print '(a/)','                                                                         '
print '(/a,a )','                ORBOPT version: ',trim(orbopt_v)
print '( a,a/)','                       git SHA: ',trim(orbopt_g(1:12))
print '(/a   )','  Nanjing -> Zurich -> Beijing: 2009 - 2017                           '
print '(a/   )', ' -----------------------------------------------------------------------'
print '(a/)', ' Using this code in publications necessitates to cite:                  '
print '(a/)', ' Y. Ma, S. Knecht, S. Keller and M. Reiher, arXiv:1611.05972 (2016).    '
print '(a/)', ' -----------------------------------------------------------------------'

call print_host(6)

end subroutine print_orbopt_header

end module orbopt_header
