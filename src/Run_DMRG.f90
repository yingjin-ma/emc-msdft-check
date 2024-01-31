      Subroutine run_DMRG_LR()

        use global_control
        use matrix

        character*192 run1,run2,run3
        character    ctmp1
        character*2  ctmp2
        character*2  ctmp3
        character*2  ctmp0
        character*192 string0,string1,string2,string3,string4
        character*192 string5,string6
        logical Lexist

        character(len=255) :: runprefix

        run1='';  run2='';  run3=''; ctmp0=""

        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        do iroot=0,dmrg_nstates-1

!          if(iroot.eq.dmrg_nstates-1)stop

          ! loop all the states
          dmrg_input_name=""
          dmrg_input_name=trim(dmrg_inputs(iroot+1))

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if
          !if(iroot.eq.0)then
            ctmp0=ctmp2
          !end if

          iexe=0
          if(LRirreps)then ! When states in the same irreps, only optimizing the lowest state
            if(iroot.eq.0) iexe=1
          else
            iexe=1
          end if

          iexe=0 
          write(6,*)"Always not re-run the dmrg calculation" 

          if(iexe.eq.1)then
            if(ifmpirun.eqv..true.)then
              run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg"&
                   //" "//trim(dmrg_input_name)//" >> "//&
                   trim(dmrg_output_name)
              run2="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg_meas"&
                   //" "//trim(dmrg_input_name)//" >> "//&
                   trim(dmrg_output_name)
            else
              run1="OMP_NUM_THREADS="//"1"//" "//"dmrg"&
                   //" "//trim(dmrg_input_name)//" >> "//&
                   trim(dmrg_output_name)
              run2="OMP_NUM_THREADS="//"1"//" "//"dmrg_meas"&
                   //" "//trim(dmrg_input_name)//" >> "//&
                   trim(dmrg_output_name)
            end if
            write(1,*)"The MPS was optimized for state",iroot
            call flush(1)
          else

            open(unit=200,file="rdm_ckp_operate0")
              ! checkpoint files
              string0=""
              string0="checkpoint_state."//trim(ctmp0)          ! string0  (ori_)
              string1=""
              string1="checkpoint_state."//trim(ctmp2)          ! string1  (ori_)
              string3=""
              string3="ref_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
              string4=""
              string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_)
              ! checkpoint file (used for calculating energy)
              !string2="rm -rf "//trim(string1)
              !write(200,*)trim(string2)
              !string2="cp -rf "//trim(string0)//" "//trim(string1)
              !write(200,*)trim(string2)
              ! ref-checkpoint file (used for RDMs-deri)
              string2="rm -rf "//trim(string3)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string3)
              write(200,*)trim(string2)
              ! initial_checkpoint file (used for back-up)
              string2="rm -rf "//trim(string4)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string4)
              write(200,*)trim(string2)
              ! Notice: no operatoration for results file
            close(200)
            call system("chmod +x rdm_ckp_operate0")
            call system("./rdm_ckp_operate0")

            ! notice that it is build_v2 for this version
!            run1=trim(runprefix)//"/dmrg"&
!               //" "//trim(dmrg_input_name)//"_Hami_gen > "&
!               //trim(dmrg_output_name)//"_Hami_gen"
!            run2=""

          end if

          string1="results_state."//trim(ctmp2)//".h5"
          write(2,*)"results_file for state ",iroot,string1

          run3="./rdmsave_su2.py "//trim(string1)

          open(unit=200,file="run_dmrg")
            write(200,*)trim(run1)
            write(200,*)trim(run2)
            write(200,*)trim(run3)
          close(200)

          if(ckp_LR)then
          else
            if(iroot.eq.0)then  ! only need to delete for starting
              call system("rm -rf *checkpoint_state.*")
              call system("rm -rf *results_state.*")
            end if
          end if
          inquire(file=trim(string1),exist=Lexist)
          if(Lexist)then
            write(2,*)"restart from previous MPS (",trim(string1),")"
          end if
          if(iroot.eq.0)then
            call system("rm *overlap.txt oneRDM*deri* twoRDM*deri*")
          end if
          call system("chmod +x run_dmrg")
          call system("./run_dmrg")

!!        if(iroot.eq.dmrg_nstates-1) stop

          open(unit=200,file="rdm_ckp_operate1")
            if(iexe.eq.1)then
              ! one-body RDMs
              string1=""
              string1="oneRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv oneparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              ! two-body RDMs
              string1=""
              string1="twoRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv twoparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              ! checkpoint files
              string0=""
              string0="checkpoint_state."//trim(ctmp2)          ! string0  (ori_)
              string3=""
              string3="ref_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
              string4=""
              string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_)
              ! ref-checkpoint file (used for RDMs-deri)
              string2="rm -rf "//trim(string3)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string3)
              write(200,*)trim(string2)
              ! initial_checkpoint file (used for back-up)
              string2="rm -rf "//trim(string4)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string4)
              write(200,*)trim(string2)
              ! Notice: no operatoration for results file
            else
              ! one-body RDMs
              string1=""
              string1="oneRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv oneparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              ! two-body RDMs
              string1=""
              string1="twoRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv twoparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              ! checkpoint files
              string0=""
              string0="checkpoint_state."//trim(ctmp2)          ! string0  (ori_)
              string1=""
              string1="checkpoint_state."//trim(ctmp0)          ! string0  (ori_)
              string3=""
              string3="ref_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
              string4=""
              string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_)
              ! ref-checkpoint file (used for RDMs-deri)
              !string2="rm -rf "//trim(string0)
              !write(200,*)trim(string2)
              !string2="cp -rf "//trim(string1)//" "//trim(string0)   ! 1 --> 0
              !write(200,*)trim(string2)
              ! initial_checkpoint file (used for back-up)
              !string2="rm -rf "//trim(string4)
              !write(200,*)trim(string2)
              !string2="cp -rf "//trim(string1)//" "//trim(string4)   ! 1 --> 4
              !write(200,*)trim(string2)
              ! Notice: no operatoration for results file
            end if
          close(200)
          call system("chmod +x rdm_ckp_operate1")
          call system("./rdm_ckp_operate1")

          write(6,*)"before RDM-deri-R/RDM-deri-L "
          call flush(6)

! ====================================================================
!             Calculating RDM-deri-R (attention: edge MPS is used)
! ====================================================================
          run1=""
          if(iexe.eq.0)then
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "&
                 //trim(runprefix)//"/dmrg "//trim(dmrg_input_name)&
                 //"_deriR > "//trim(dmrg_output_name)//"_deriR"
          else
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "&
                 //trim(runprefix)//"/dmrg "//trim(dmrg_input_name)&
                 //"_SP-deriR > "//trim(dmrg_output_name)//"_deriR"
          end if
          write(6,*)"RDM-deri part need to be re-code"

          open(unit=200,file="run_dmrg_deriR")
            write(200,*)trim(run1)
          close(200)
          !call system("rm overlap.txt")
          call system("chmod +x run_dmrg_deriR")
          call system("./run_dmrg_deriR")

          open(unit=200,file="after_deriR")
           ! rename these RDMs-deri to the right-hand-side
            string1=""
            string1="RDM."//trim(ctmp2)//".deriR_"
            !Debian or Redhat
            if(linux.eq.'ubuntu')then
              string2="rename s/RDMderi_/"//trim(string1)//&
                      "/ *RDMderi_*"
            else
              string2="rename RDMderi_ "//trim(string1)//&
                      " *RDMderi_*"
            end if
            write(200,*)trim(string2)
            ! Save MPS-overlap (later for correctting the sign)
            string2=""
            string2="mv overlap.txt Roverlap."//trim(ctmp2)//".txt"
            write(200,*)trim(string2)
            string2=""
            string2="rm -rf "//trim(string0)
            write(200,*)trim(string2)
            write(200,*) "cp "//trim(dmrg_output_name)//"_deriR "//trim(dmrg_output_name)//"_deriR."//trim(ctmp2)
            ! Recover the checkpoint file in case it is changed
            string2="";
            string2="cp -rf "//trim(string4)//" "//trim(string0)
            write(200,*)trim(string2)
          close(200)
          call system("chmod +x after_deriR")
          call system("./after_deriR")

          write(6,*)"After RDM-deri-R/ before RDM-deri-L "
          call flush(6)
          !if(iroot.eq.0)stop
          !if(iroot.eq.1)stop

! ====================================================================
!             Calculating RDM-deri-L (attention: edge MPS is used)
! ====================================================================
         ! Reset to follow the normal Davidson

          if(iexe.eq.0)then
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "&
                 //trim(runprefix)//"/dmrg "//trim(dmrg_input_name)&
                 //"_deriL > "//trim(dmrg_output_name)//"_deriL"
          else
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "&
                 //trim(runprefix)//"/dmrg "//trim(dmrg_input_name)&
                 //"_SP-deriL > "//trim(dmrg_output_name)//"_deriL"
          end if

          open(unit=200,file="run_dmrg_deriL")
            write(200,*)trim(run1)
          close(200)
!          call system("rm overlap.txt")
          call system("chmod +x run_dmrg_deriL")
          call system("./run_dmrg_deriL")

          open(unit=200,file="after_deriL")
            ! rename these RDMs-deri to the left-hand-side
            string1=""
            string1="RDM."//trim(ctmp2)//".deriL_"
            !Debian or Redhat
            if(linux.eq.'ubuntu')then
              string2="rename s/RDMderi_/"//trim(string1)//&
                      "/ *RDMderi_*"
            else
              string2="rename RDMderi_ "//trim(string1)//&
                      " *RDMderi_*"
            end if
            write(200,*)trim(string2)
            write(200,*) "cp "//trim(dmrg_output_name)//"_deriL "//trim(dmrg_output_name)//"_deriL."//trim(ctmp2)
            ! Save MPS-overlap (later for correctting the sign)
            string2=""
            string2="mv overlap.txt Loverlap."//trim(ctmp2)//".txt"
            write(200,*)trim(string2)
          close(200)
          call system("chmod +x after_deriL")
          call system("./after_deriL")

! ====================================================================
!     Generate the full Hamiltonian (attention: edge MPO is used)
! ====================================================================

!          if(iexe.eq.0)then
            ! no need to be re-run
!          else
            write(6,*)"Before full Hamiltonian "
            call flush(6)

            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "&
                 //trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_input_name)//"_Hami_gen > "& !!! TODO: Sure it's _gs ?
                 //trim(dmrg_output_name)//"_Hami_gen"
           write(*,*) run1
            open(unit=200,file="run_dmrg_gen")
              write(200,*)trim(run1)

              ! checkpoint files
              string0=""
              string0="checkpoint_state."//trim(ctmp2)          ! string0  (ori_)
              string1=""
              string1="checkpoint_state."//trim(ctmp2)          ! string1  (ori_)
              string3=""
              string3="ref_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
              string4=""
              string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_)
              ! checkpoint file (used for calculating energy)
              string2="rm -rf "//trim(string1)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string4)//" "//trim(string1)
              write(200,*)trim(string2)
              ! ref-checkpoint file (used for RDMs-deri)
              string2="rm -rf "//trim(string3)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string4)//" "//trim(string3)
              write(200,*)trim(string2)
              ! initial_checkpoint file (used for back-up)
              ! Notice: no operatoration for results file

            close(200)
            call system("chmod +x run_dmrg_gen")
            call system("./run_dmrg_gen")
!          end if

!!          if(iroot.eq.dmrg_nstates-1) stop

! ====================================================================
!          Rename the coefficients from Davidson vectors
! ====================================================================
          open(unit=200,file="rdm_ckp_operate2")
            string2="rm -rf "//trim(string0)
            write(200,*)trim(string2)
            string2="cp -rf "//trim(string3)//" "//trim(string0)
            write(200,*)trim(string2)

            string1="MPSCi."//trim(ctmp2)//".info"
            string2="mv MPSCi.info "//string1
            write(200,*)trim(string2)

            string1="MPSCi_value."//trim(ctmp2)//".info"
            string2="mv MPSCi_value.info "//string1
            write(200,*)trim(string2)

            string1="MPSCi."//trim(ctmp2)//".davidson"
            string2="mv MPSCi.davidson "//string1
            write(200,*)trim(string2)

            string1="Hami."//trim(ctmp2)//".txt"
            string2="mv Hami.txt "//string1
            write(200,*)trim(string2)
          close(200)
          call system("chmod +x rdm_ckp_operate2")
          call system("./rdm_ckp_operate2")
! ====================================================================

          write(6,*)"After all things for this state"
          call flush(6)

        end do  ! do for every states

        write(6,*)"Before start TDM calculations"
        call flush(6)

! Now, start to calculate
        itdm=0
        do iroot=0,dmrg_nstates-1
          do jroot=iroot+1,dmrg_nstates-1

            ! iroot
            ctmp2=""
            if(iroot.ge.10)then
              write(ctmp2,"(I2)")iroot
            else
              write(ctmp1,"(I1)")iroot
              ctmp2(1:1)=ctmp1
            end if
            ! jroot
            ctmp3=""
            if(jroot.ge.10)then
              write(ctmp3,"(I2)")jroot
            else
              write(ctmp1,"(I1)")jroot
              ctmp3(1:1)=ctmp1
            end if

            string3=""
            string3="checkpoint_state."//trim(ctmp2)          ! string3  (ori_)
            string4=""
            string4="ref_checkpoint_state."//trim(ctmp2)      ! string4  (ref_)
            string5=""
            string5="checkpoint_state."//trim(ctmp3)          ! string5  (ori_)
            string6=""
            string6="ref_checkpoint_state."//trim(ctmp3)      ! string6  (ref_)

! ===================================================================
!           Calculate TDMs (attention: edge MPO is used)
! ===================================================================

! 1)        iroot.jroot
            itdm=itdm+1
            run1=trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_inputs_tdm(itdm))//" >> "& ! dmrg ?
                 //"dmrg_tdms_L.out"

            string0=""
            string0="results_state."&
                    //trim(ctmp2)//"."//trim(ctmp3)//".h5"  ! string0  (ori_)
            run2="./rdmsave_su2.py "//trim(string0)

            open(unit=200,file="run_tdm_L")
              ! ensure ref_checkpoint_state.X is not modified
              string2="rm -rf "//trim(string4)
              write(200,*)trim(string2)
              string2="cp -r "//trim(string3)//" "//trim(string4)
              write(200,*)trim(string2)

              ! ensure ref_checkpoint_state.X is not modified
              string2="rm -rf "//trim(string6)
              write(200,*)trim(string2)
              string2="cp -r "//trim(string5)//" "//trim(string6)
              write(200,*)trim(string2)

              write(200,*)trim(run1)
              write(200,*)trim(run2)
              string1=""
              string1="oneRDM."//trim(ctmp2)//"."//trim(ctmp3)
              string2="mv oneparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              string1=""
              string1="twoRDM."//trim(ctmp2)//"."//trim(ctmp3)
              string2="mv twoparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
            close(200)

            call system("chmod +x run_tdm_L")
            call system("./run_tdm_L")

! 2)        jroot.iroot
            itdm=itdm+1

            run1=trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_inputs_tdm(itdm))//" >> "&
                 //"dmrg_tdms_R.out"

            string0=""
            string0="results_state."&
                    //trim(ctmp3)//"."//trim(ctmp2)//".h5"  ! string0  (ori_)
            run2="./rdmsave_su2.py "//trim(string0)

            open(unit=200,file="run_tdm_R")

              ! ensure ref_checkpoint_state.X is not modified
              string2="rm -rf "//trim(string4)
              write(200,*)trim(string2)
              string2="cp -r "//trim(string3)//" "//trim(string4)
              write(200,*)trim(string2)

              ! ensure ref_checkpoint_state.X is not modified
              string2="rm -rf "//trim(string6)
              write(200,*)trim(string2)
              string2="cp -r "//trim(string5)//" "//trim(string6)
              write(200,*)trim(string2)

              write(200,*)trim(run1)
              write(200,*)trim(run2)
              string1=""
              string1="oneRDM."//trim(ctmp3)//"."//trim(ctmp2)
              string2="mv oneparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
              string1=""
              string1="twoRDM."//trim(ctmp3)//"."//trim(ctmp2)
              string2="mv twoparticle.rdm "//trim(string1)
              write(200,*)trim(string2)
            close(200)

            call system("chmod +x run_tdm_R")
            call system("./run_tdm_R")

          end do
        end do

        write(6,*)"Before finishing run_DMRG_LR"
        call flush(6)

      End Subroutine run_DMRG_LR

      Subroutine run_DMRG(icycle,iloop)

        use global_control
        use matrix
        use date_time

        character*192 run1
        character*192 run2
        character*192 run3
        character*192 rdm1
        character*192 rdm2
        character    ctmp1
        character*2  ctmp2
        character*192 string0,string1,string2,string3,string4
        integer::icycle,iloop
        logical lexist

        character(len=255) :: runprefix

        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        run1=''
        run2=''
        run3=''

        dmrg_rdm1_name=""
        dmrg_rdm2_name=""


        if(dmrg_binary_name.ne."maquis")then

          dmrg_rdm1_name="dmrg.conf_rdm1"
          dmrg_rdm2_name="dmrg.conf_rdm2"

          rdm1="cat "//dmrg_input_name//">> "//dmrg_rdm1_name
          rdm2="cat "//dmrg_input_name//">> "//dmrg_rdm2_name

          open(unit=130,file=dmrg_rdm1_name)
            write(130,1300)"restart"
            write(130,1301)"onerdm"
          close(130)
            open(unit=131,file='temp_cat_rdm1')
              write(131,*)rdm1
          close(131)
          call system("chmod +x temp_cat_rdm1")
          call system("./temp_cat_rdm1")

          open(unit=130,file=dmrg_rdm2_name)
            write(130,1300)"restart"
            write(130,1301)"twordm"
          close(130)
            open(unit=131,file='temp_cat_rdm2')
              write(131,*)rdm2
            close(131)
            call system("chmod +x temp_cat_rdm2")
            call system("./temp_cat_rdm2")
1300      format(A7)
1301      format(A6)

          if(ifmpirun.eqv..true.)then
            run1="mpiexec -np "//nproc_in_mpirun//" "//dmrg_binary_name&
                 //" "//dmrg_input_name//" >> "//dmrg_output_name
            run2="mpiexec -np "//nproc_in_mpirun//" "//dmrg_binary_name&
                 //" "//dmrg_rdm1_name//" >> "//dmrg_output_name
            run3="mpiexec -np "//nproc_in_mpirun//" "//dmrg_binary_name&
               //" "//dmrg_rdm2_name//" >> "//dmrg_output_name
          else
            run1=dmrg_binary_name&
                 //" "//dmrg_input_name//" >> "//dmrg_output_name
            run2=dmrg_binary_name&
                 //" "//dmrg_rdm1_name//" >> "//dmrg_output_name
            run3=dmrg_binary_name&
                 //" "//dmrg_rdm2_name//" >> "//dmrg_output_name
          end if

          open(unit=200,file="run_dmrg")
            write(200,*)run1
            write(200,*)run2
            write(200,*)run3
          close(200)

          !write(*,*)"DMRG calculation started"
!          goto 2222
!          write(*,*)"DMRG calculation ignored"
          call system("chmod +x run_dmrg")
          !write(*,*)"DMRG calculation middle1"
          call system("./run_dmrg")
          !write(*,*)"DMRG calculation middle2"
          call system("mv ./tmp/spatial_onepdm.0.0 .")
          !write(*,*)"DMRG calculation middle3"
          call system("mv ./tmp/spatial_twopdm.0.0 .")

        else
          ! write(6,*)"*********dmrg_nstates*********",dmrg_nstates

          do iroot=0,dmrg_nstates-1
            ! write(6,*)"*********dmrg_nstates_in*********",dmrg_nstates

            ! loop all the states
            dmrg_input_name=""
            dmrg_input_name=trim(dmrg_inputs(iroot+1))

            if(iloop.ne.0)then ! if in the micro-iterations

              write(2,*)"In the micro-iter, only 1 sweep"

              if(ifmpirun.eqv..true.)then
                run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg"&
                    //" "//trim(dmrg_input_name)//"_micro >> "//&
                    trim(dmrg_output_name)//"_micro_out"
                run2="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg_meas"&
                    //" "//trim(dmrg_input_name)//"_micro >> "//&
                    trim(dmrg_output_name)//"_micro_out"
              else
                run1="OMP_NUM_THREADS="//"1"//" "//"dmrg"&
                    //" "//trim(dmrg_input_name)//"_micro >> "//&
                    trim(dmrg_output_name)//"_micro_out"
                run2="OMP_NUM_THREADS="//"1"//" "//"dmrg_meas"&
                    //" "//trim(dmrg_input_name)//"_micro >> "//&
                    trim(dmrg_output_name)//"_micro_out"
              end if

            else  ! If in the macro-iterations

              write(2,*)"In the macro-iter, fully sweeped"

              if(ifmpirun.eqv..true.)then
                
                run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg"&
                    //" "//trim(dmrg_input_name)//" >> "//&
                    trim(dmrg_output_name)
                run2="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg_meas"&
                    //" "//trim(dmrg_input_name)//" >> "//&
                    trim(dmrg_output_name)
                ! write(6,*)"run1",run1
              else
                run1="OMP_NUM_THREADS="//"1"//" "//"dmrg"&
                    //" "//trim(dmrg_input_name)//" >> "//&
                    trim(dmrg_output_name)
                run2="OMP_NUM_THREADS="//"1"//" "//"dmrg_meas"&
                    //" "//trim(dmrg_input_name)//" >> "//&
                    trim(dmrg_output_name)
              end if

            end if

            ctmp2=""
            if(iroot.ge.10)then
              write(ctmp2,"(I2)")iroot
            else
              write(ctmp1,"(I1)")iroot
              ctmp2(1:1)=ctmp1
            end if

            string1="results_state."//trim(ctmp2)//".h5"
            write(2,*)"results_file for state ",iroot,string1

            run3="./rdmsave_su2.py "//trim(string1)

            open(unit=200,file="run_dmrg")
              write(200,*)trim(run1)
              write(200,*)trim(run2)
              write(200,*)trim(run3)
            close(200)

            if(iloop.eq.0)then
              if(restart_maquis)then
              else
                if(icycle.eq.1.or.ith_INTE.eq.1)then
                  if(iroot.eq.0)then  ! only need to delete for starting
                    call system("rm -rf *checkpoint_state.*")
                    call system("rm -rf *results_state.*")
                  end if
                end if
              end if
            else
!             inquire(file="initial_checkpoint_state.0",exist=Lexist)
              inquire(file=trim(string1),exist=Lexist)
              if(Lexist)then
!               call system&
!               ("cp -rf initial_checkpoint_state.0 checkpoint_state.0")
!               write(2,*)"copy initial_ckp as ckp"
                write(2,*)"reference w.f. (",trim(string1),") used"
              end if
            end if

            if(method(icycle).eq."cp-hess".or.ith_INTE.eq.1)then
              if(iroot.eq.0)then
                call system("rm *overlap.txt oneRDM*deri* twoRDM*deri*")
              end if
            end if
            call system("chmod +x run_dmrg")
            walltime(21) = wtime()
            call system("./run_dmrg")

            walltime(22) = wtime()
            write(3,"(I4,A)")&
            iroot,"-th root:"
            call timing_reporter(3,"sweep+RDM",walltime(22)-walltime(21))

!            stop

            open(unit=200,file="rdm_ckp_operate1")

              ! one-body RDMs
              string1=""
              string1="oneRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv oneparticle.rdm "//trim(string1)
              write(200,*)trim(string2)

              ! two-body RDMs
              string1=""
              string1="twoRDM."//trim(ctmp2)//"."//trim(ctmp2)
              string2="mv twoparticle.rdm "//trim(string1)
              write(200,*)trim(string2)

              string0=""
              string0="checkpoint_state."//trim(ctmp2)          ! string0  (ori_)
              string3=""
              string3="ref_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
              string4=""
              string4="initial_checkpoint_state."//trim(ctmp2)  ! string4  (ini_)
              ! ref-checkpoint file (used for RDMs-deri)
              string2="rm -rf "//trim(string3)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string3)
              write(200,*)trim(string2)
              ! initial_checkpoint file (used for back-up)
              string2="rm -rf "//trim(string4)
              write(200,*)trim(string2)
              string2="cp -rf "//trim(string0)//" "//trim(string4)
              write(200,*)trim(string2)

              ! Notice: no operatoration for results file

            close(200)
            call system("chmod +x rdm_ckp_operate1")
            call system("./rdm_ckp_operate1")

          ! For the RDM-deri
            if(method(icycle).eq."cp-hess".or.ith_INTE.eq.1)then

!              if(iroot.eq.0)then
!                call system("rm -rf *RDM*deri*")
!              end if

              goto 1133
! ====================================================================
!          Only used for HF case (M=1, i.e. M=0 in QCMaquis)
! ====================================================================
              call system("mv oneRDM.0.0 oneRDM.normal")
              call system("mv twoRDM.0.0 twoRDM.normal")
              run1=""
            run1=trim(runprefix)//"/dmrg"&
                   //" "//trim(dmrg_input_name)//"_mps_gen > "&
                   //trim(dmrg_output_name)//"_mps_gen"
              open(unit=200,file="run_dmrg_gen")
                write(200,*)trim(run1)
              close(200)
              call system("chmod +x run_dmrg_gen")
              call system("./run_dmrg_gen")
              call system("./rdmsave_su2.py results_state.0.h5")

              call system("mv oneparticle.rdm oneRDM.0.0")
              call system("mv twoparticle.rdm twoRDM.0.0")
              call&
              system("cp -rf checkpoint_state.0 ref_checkpoint_state.0")
!             call system("cp results_state.0 ref_results_state.0")
! ====================================================================
1133          continue

! ====================================================================
!             Calculating RDM-deri-R
! ====================================================================
            run1=trim(runprefix)//"/dmrg"& !! ground or excited state???
                   //" "//trim(dmrg_input_name)//"_deriR > "&
                   //trim(dmrg_output_name)//"_deriR"
              open(unit=200,file="run_dmrg_deriR")
                write(200,*)trim(run1)
              close(200)
!              call system("rm overlap.txt")
              call system("chmod +x run_dmrg_deriR")
              call system("./run_dmrg_deriR")

              open(unit=200,file="after_deriR")
                ! rename these RDMs-deri to the right-hand-side
                string1=""
                string1="RDM."//trim(ctmp2)//".deriR_"
                !Debian or Redhat
                if(linux.eq.'ubuntu')then
                  string2="rename s/RDMderi_/"//trim(string1)//&
                          "/ *RDMderi_*"
                else
                  string2="rename RDMderi_ "//trim(string1)//&
                          " *RDMderi_*"
                end if
                write(200,*)trim(string2)
                ! Save MPS-overlap (later for correctting the sign)
                string2=""
                !string2="cp overlap.txt Roverlap.txt"
                string2="mv overlap.txt Roverlap."//trim(ctmp2)//".txt"
                write(200,*)trim(string2)
                string2=""
                string2="rm -rf "//trim(string0)
                write(200,*)trim(string2)
                ! Recover the checkpoint file in case it is changed
                string2=""
                string2="cp -rf "//trim(string3)//" "//trim(string0)
                write(200,*)trim(string2)
              close(200)
              call system("chmod +x after_deriR")
              walltime(21) = wtime()
              call system("./after_deriR")

              walltime(22) = wtime()
              call timing_reporter(3,"deri-R",walltime(22)-walltime(21))

! ====================================================================
!             Calculating RDM-deri-L
! ====================================================================
             ! Reset to follow the normal Davidson

            run1=trim(runprefix)//"/dmrg"&
                   //" "//trim(dmrg_input_name)//"_deriL > "&
                   //trim(dmrg_output_name)//"_deriL"
              open(unit=200,file="run_dmrg_deriL")
                write(200,*)trim(run1)
              close(200)
              call system("rm overlap.txt")
              call system("chmod +x run_dmrg_deriL")
              call system("./run_dmrg_deriL")

              open(unit=200,file="after_deriL")
                ! rename these RDMs-deri to the left-hand-side
                string1=""
                string1="RDM."//trim(ctmp2)//".deriL_"
                !Debian or Redhat
                if(linux.eq.'ubuntu')then
                  string2="rename s/RDMderi_/"//trim(string1)//&
                          "/ *RDMderi_*"
                else
                  string2="rename RDMderi_ "//trim(string1)//&
                          " *RDMderi_*"
                end if
                write(200,*)trim(string2)
                ! Save MPS-overlap (later for correctting the sign)
                string2=""
                !string2="cp overlap.txt Loverlap.txt"
                string2="mv overlap.txt Loverlap."//trim(ctmp2)//".txt"
                write(200,*)trim(string2)
              close(200)
              call system("chmod +x after_deriL")
              walltime(21) = wtime()
              call system("./after_deriL")

              walltime(22) = wtime()
              call timing_reporter(3,"deri-L",walltime(22)-walltime(21))

! ====================================================================
!             Generate the full Hamiltonian (4 time larger one)
! ====================================================================
              run1=trim(runprefix)//"/dmrg"& ! ground or excited state??
                   //" "//trim(dmrg_input_name)//"_Hami_gen > "&
                   //trim(dmrg_output_name)//"_Hami_gen"
              open(unit=200,file="run_dmrg_gen")
                write(200,*)trim(run1)
              close(200)
              call system("chmod +x run_dmrg_gen")
              walltime(21) = wtime()
              call system("./run_dmrg_gen")

              walltime(22) = wtime()
              call timing_reporter(3,"local Hamiltonian",walltime(22)-walltime(21))
! ====================================================================

! ====================================================================
!             Rename the coefficients from Davidson vectors
! ====================================================================
              open(unit=200,file="rdm_ckp_operate2")
                string2="rm -rf "//trim(string0)
                write(200,*)trim(string2)
                string2="cp -rf "//trim(string4)//" "//trim(string0)
                write(200,*)trim(string2)

                string1="MPSCi."//trim(ctmp2)//".info"
                string2="mv MPSCi.info "//string1
                write(200,*)trim(string2)

                string1="MPSCi_value."//trim(ctmp2)//".info"
                string2="mv MPSCi_value.info "//string1
                write(200,*)trim(string2)

                string1="MPSCi."//trim(ctmp2)//".davidson"
                string2="mv MPSCi.davidson "//string1
                write(200,*)trim(string2)

                string1="Hami."//trim(ctmp2)//".txt"
                string2="mv Hami.txt "//string1
                write(200,*)trim(string2)
              close(200)
              call system("chmod +x rdm_ckp_operate2")
              call system("./rdm_ckp_operate2")

            end if  ! if for RDM-deri calculations
          end do  ! do for every states

        end if

!2222    continue

        !call system("rm run_dmrg")
        !call system("rm temp_cat*")
        !write(*,*)icycle
!        if(icycle.eq.2)stop
        write(2,*)"   DMRG calculation finished   "
        write(2,*)"  Currently no TDM calculation  "

      end Subroutine run_DMRG

! ==============================================
      Subroutine run_DMRG_diag(icycle)
        ! later need to change to the only Hdiag version of DMRG
        !  and only MPO on the edge

        use global_control
        use matrix

        character*192 run1
        character*192 run2
        character*192 run3
        integer::icycle

        character(len=255) :: runprefix
        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        run1=''
        run2=''
        run3=''

        goto 1144

          if(ifmpirun.eqv..true.)then
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//"dmrg"&
                 //" "//dmrg_input_name//" >> "//dmrg_output_name
          else
            run1="OMP_NUM_THREADS="//"1"//" "//"dmrg"&
                 //" "//dmrg_input_name//" >> "//dmrg_output_name
          end if

          open(unit=200,file="run_dmrg")
            write(200,*)run1
          close(200)

          call system("chmod +x run_dmrg")
          call system("./run_dmrg")

          call system("rm -rf ref_checkpoint_state.0")
          call system&
          ("cp -rf checkpoint_state.0 ref_checkpoint_state.0")

1144      continue

! ====================================================================
!            generate the full Hamiltonian (4 time larger one)
! ====================================================================
          run1=trim(runprefix)//"/dmrg"&
               //" "//trim(dmrg_input_name)//"_Hami_gen > "&
               //trim(dmrg_output_name)//"_Hami_gen"
          open(unit=200,file="run_dmrg_gen")
            write(200,*)run1
          close(200)
          call system("chmod +x run_dmrg_gen")
          call system("./run_dmrg_gen")

!          call system("rm -rf checkpoint_state.0")
!          call system("rm -rf results_state.0.h5")

!          call system&
!          ("cp -rf initial_results_state.0.h5 results_state.0.h5")
!          call system&
!          ("cp -rf initial_checkpoint_state.0 checkpoint_state.0")

          write(2,*)"finished the run_DMRG_diag"
          call flush(2)
!          stop

      end Subroutine run_DMRG_diag

! ==============================================
! Update the MPSci part for Lagrange  --  L

! ==============================================
! Update the MPSci part for Lagrange  --  R
      Subroutine MPSci_lagrange_update(iroot,LorR)

        use global_control
        use matrix

          integer :: iroot
        character ::  LorR

        character*192 run1,run2,run3
        character    ctmp1
        character*2  ctmp2
        character*2  ctmp3
        character*192 string0,string1,string2,string3,string4
        character*192 string5,string6
        logical Lexist

        character(len=255) :: runprefix

        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        run1='';  run2='';  run3=''

        ctmp2=""
        if(iroot.ge.10)then
          write(ctmp2,"(I2)")iroot
        else
          write(ctmp1,"(I1)")iroot
          ctmp2(1:1)=ctmp1
        end if

        string1="rm -rf L_checkpoint_state."//trim(ctmp2)
        string2="rm -rf R_checkpoint_state."//trim(ctmp2)
        string3="cp -rf ref_checkpoint_state."//trim(ctmp2)//&
                " L_checkpoint_state."//trim(ctmp2)
        string4="cp -rf ref_checkpoint_state."//trim(ctmp2)//&
                " R_checkpoint_state."//trim(ctmp2)

        open(unit=200,file="ckp_operate_LR")
          write(200,*)trim(string1)
          write(200,*)trim(string2)
          write(200,*)trim(string3)
          write(200,*)trim(string4)
        close(200)

        call system("chmod +x ckp_operate_LR")
        call system("./ckp_operate_LR")

        ! loop all the states
        dmrg_input_name=""
        dmrg_input_name=trim(dmrg_inputs(iroot+1))

        if(ifmpirun.eqv..true.)then
          if(LorR.eq."L")then
            write(1,*)"dmrg_L"
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//         & ! ground or excited state?
            trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_input_name)//"_LR-L >> "//          &
                 trim(dmrg_output_name)//"."//trim(ctmp2)//"_LR-L"
          else
            write(1,*)"dmrg_R"
            run1="OMP_NUM_THREADS="//nproc_in_mpirun//" "//         &
            trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_input_name)//"_LR-R >> "//          &
                 trim(dmrg_output_name)//"."//trim(ctmp2)//"_LR-R"
          end if
        else
          if(LorR.eq."L")then
            write(1,*)"dmrg_L"
            run1="OMP_NUM_THREADS="//"1"//" "//                     &
            trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_input_name)//"_LR-L >> "//          &
                 trim(dmrg_output_name)//"."//trim(ctmp2)//"_LR-L"
          else
            write(1,*)"dmrg_R"
            run1="OMP_NUM_THREADS="//"1"//" "//                     &
            trim(runprefix)//"/dmrg"&
                 //" "//trim(dmrg_input_name)//"_LR-R >> "//          &
                 trim(dmrg_output_name)//"."//trim(ctmp2)//"_LR-R"
          end if
        end if

! If during sweeps
        if(LorR.eq."L")then
          string1="mv oneRDM_1PT oneRDM."//trim(ctmp2)//".LR-L"
          string2="mv twoRDM_1PT twoRDM."//trim(ctmp2)//".LR-L"
        else
          string1="mv oneRDM_1PT oneRDM."//trim(ctmp2)//".LR-R"
          string2="mv twoRDM_1PT twoRDM."//trim(ctmp2)//".LR-R"
        end if

        open(unit=200,file="run_dmrg_LR")
          write(200,*)trim(run1)
          write(200,*)trim(string1)
          write(200,*)trim(string2)
        close(200)

        call system("chmod +x run_dmrg_LR")
        call system("./run_dmrg_LR")

      End Subroutine MPSci_lagrange_update

! ==============================================
      Subroutine run_RDMs_update()

        character*192 run1
        character*192 run2
        character*192 run3
        integer::icycle

        run1=''
        run2=''
        run3=''

      end Subroutine run_RDMs_update

! ==============================================
      Subroutine run_mpsci_update(iroot)
        ! update the MPSci

        use global_control
        use matrix

        character*192 run1
        character*192 run2
        character*192 run3
        integer::iroot

        character    ctmp1
        character*2  ctmp2
        character*192 string0,string1,string2,string3,string4

        character(len=255) :: runprefix
        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

          run1=''
          run2=''
          run3=''

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if

          open(unit=200,file="ckp_operate_up")
            string0="rm -rf UP_checkpoint_state."//trim(ctmp2)
            write(200,*)trim(string0)
            string0="cp -rf initial_checkpoint_state."//trim(ctmp2)&
                    //" UP_checkpoint_state."//trim(ctmp2)
            write(200,*)trim(string0)
          close(200)
          call system("chmod +x ckp_operate_up")
          call system("./ckp_operate_up")

          ! loop all the states
          dmrg_input_name=""
          dmrg_input_name=trim(dmrg_inputs(iroot+1))

          run1=""
          run1=trim(runprefix)//"/dmrg"&
               //" "//trim(dmrg_input_name)//"_MPS_update > "&
               //trim(dmrg_output_name)//"_MPS_update"
          open(unit=200,file="run_mpsci_update")
            write(200,*)run1
          close(200)
          call system("chmod +x run_mpsci_update")
          call system("./run_mpsci_update")

!          call system("rm -rf checkpoint_state.*")
!          call system("rm -rf results_state.0.h5")

!          call system&
!          ("cp -rf initial_results_state.0.h5 results_state.0.h5")
!          call system&
!          ("cp -rf initial_checkpoint_state.0 checkpoint_state.0")

! If using the result-file
!          call system("./rdmsave_su2.py results_state.0.h5")
!          call system("mv oneparticle.rdm oneRDM.0.0")
!          call system("mv twoparticle.rdm twoRDM.0.0")
! If during sweeps

      End Subroutine run_mpsci_update

      Subroutine run_DMRG_direct(icycle)

        use global_control
        use matrix
        use date_time

        character*192 run1
        character*192 run2
        character*192 run3
        integer::icycle

        character(len=255) :: runprefix
        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        dmrg_rdm1_name=""
        dmrg_rdm2_name=""

        run1=''
        run2=''
        run3=''

          run1=trim(runprefix)//"/dmrg"&
               //" "//trim(dmrg_input_name)//"_update > "&
               //trim(dmrg_output_name)//"_update"
          run3="./rdmsave_su2.py edge_results_state.0.h5"
          open(unit=200,file="run_DMRG_update")
            write(200,*)run1
            write(200,*)run3
          close(200)
          call system("chmod +x run_DMRG_update")
          walltime(21) = wtime()
          call system("./run_DMRG_update")


          call system("mv oneparticle.rdm oneRDM.0.0")
          call system("mv twoparticle.rdm twoRDM.0.0")

          call system("rm -rf ref*_checkpoint_state.0")
          call system("rm -rf ref*_results_state.0")

          call&
          system("cp -rf checkpoint_state.0 ref_checkpoint_state.0")
          call&
          system("cp -rf results_state.0.h5 ref_results_state.0.h5")

          walltime(22) = wtime()
          call timing_reporter(3,"partial-DMRG-done(micro-sweep+RDM)",walltime(22)-walltime(21))

!          call&
!          system("cp -rf checkpoint_state.0 initial_checkpoint_state.0")
!          call&
!          system("cp -rf results_state.0.h5 initial_results_state.0.h5")

      end Subroutine run_DMRG_direct

! Re-run to get the updated RDMderi-s
      Subroutine RDMderi_update()

        use global_control
        use matrix

        character*192 run1,run2,run3
        character    ctmp1
        character*2  ctmp2
        character*2  ctmp3
        character*192 string0,string1,string2,string3,string4
        character*192 string5,string6
        logical Lexist

        character(len=255) :: runprefix
        call getenv('QCMAQUIS_LR_PATH',runprefix)
        if (trim(runprefix).eq.'') then
          write (*,*) "WARNING: QCMAQUIS_LR_PATH is not set."
        end if

        run1='';  run2='';  run3=''

        do iroot=0,dmrg_nstates-1

          ctmp2=""
          if(iroot.ge.10)then
            write(ctmp2,"(I2)")iroot
          else
            write(ctmp1,"(I1)")iroot
            ctmp2(1:1)=ctmp1
          end if

          open(unit=200,file="rdm_ckp_operate_mi")
            ! checkpoint files
            string0=""
            string0="checkpoint_state."//trim(ctmp2)         ! string0  (ori_)
            string1=""
            string1="UP_checkpoint_state."//trim(ctmp2)      ! string3  (ref_)
            string3=""
            string3="ref_checkpoint_state."//trim(ctmp2)     ! string3  (ref_)
            string4=""
            string4="initial_checkpoint_state."//trim(ctmp2) ! string4  (ini_)

            string2="rm -rf "//trim(string0)
            write(200,*)trim(string2)
            string2="rm -rf "//trim(string3)
            write(200,*)trim(string2)
            string2="cp -r "//trim(string4)//" "//trim(string0)
            write(200,*)trim(string2)
            string2="cp -r "//trim(string1)//" "//trim(string3) ! ref should be updated
            write(200,*)trim(string2)

            ! Notice: no operatoration for results file
          close(200)
          call system("chmod +x rdm_ckp_operate_mi")
          call system("./rdm_ckp_operate_mi")

! to Stefan : Eq.74 Eq.75 (Left)
! ====================================================================
!             Calculating RDM-deri-R (attention: edge MPS is used)
! ====================================================================

          run1=trim(runprefix)//"/dmrg"&
               //" "//trim(dmrg_input_name)//"_deri_R >> "&
               //trim(dmrg_output_name)//"_deriR_mi"
          open(unit=200,file="run_dmrg_deriR")
            write(200,*)trim(run1)
          close(200)
          !call system("rm overlap.txt")
          call system("chmod +x run_dmrg_deriR")
          call system("./run_dmrg_deriR")

          open(unit=200,file="after_deriR")
           ! rename these RDMs-deri to the right-hand-side
            string1=""
            string1="RDM."//trim(ctmp2)//".deriR_"
            if(linux.eq.'ubuntu')then
              string2="rename s/RDMderi_/"//trim(string1)//&
                      "/ *RDMderi_*"
            else
              string2="rename RDMderi_ "//trim(string1)//&
                      " *RDMderi_*"
            end if
            write(200,*)trim(string2)
            ! Save MPS-overlap (later for correctting the sign)
            string2=""
            string2="mv overlap.txt Roverlap."//trim(ctmp2)//".txt"
            write(200,*)trim(string2)
            string2=""
            string2="rm -rf "//trim(string0)
            write(200,*)trim(string2)
            ! Recover the checkpoint file in case it is changed
            string2="";
            string2="cp -rf "//trim(string3)//" "//trim(string0)
            write(200,*)trim(string2)
          close(200)
          call system("chmod +x after_deriR")
          call system("./after_deriR")

! to Stefan : Eq.74 Eq.75 (Right)
! ====================================================================
!             Calculating RDM-deri-L (attention: edge MPS is used)
! ====================================================================
         ! Reset to follow the normal Davidson
          run1=trim(runprefix)//"/dmrg"&
               //" "//trim(dmrg_input_name)//"_deri_L >> "&
               //trim(dmrg_output_name)//"_deriL_mi"
          open(unit=200,file="run_dmrg_deriL")
            write(200,*)trim(run1)
          close(200)
!          call system("rm overlap.txt")
          call system("chmod +x run_dmrg_deriL")
          call system("./run_dmrg_deriL")

          open(unit=200,file="after_deriL")
            ! rename these RDMs-deri to the left-hand-side
            string1=""
            string1="RDM."//trim(ctmp2)//".deriL_"
            if(linux.eq.'ubuntu')then
              string2="rename s/RDMderi_/"//trim(string1)//&
                      "/ *RDMderi_*"
            else
              string2="rename RDMderi_ "//trim(string1)//&
                      " *RDMderi_*"
            end if
            write(200,*)trim(string2)
            ! Save MPS-overlap (later for correctting the sign)
            string2=""
            string2="mv overlap.txt Loverlap."//trim(ctmp2)//".txt"
            write(200,*)trim(string2)
          close(200)
          call system("chmod +x after_deriL")
          call system("./after_deriL")
        end do

      End Subroutine RDMderi_update
