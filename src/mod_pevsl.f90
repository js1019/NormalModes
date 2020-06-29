module pevsl_mod

  use mpi
  use omp_lib
  use para_mod,                 only: rkind,pin
  use utility_mod
  use geometry_mod,             only: unstrM
  !use cg_create_matrix_mod,     only: CGM
  use cg_matvec_mod


  implicit none 
  
 contains

  subroutine pnm_apply_pevsl()
     integer                                      :: I, K, J, IERR
     integer*8                                    :: POL, CHEBSOL
     integer                                      :: LANSTEP,MDEG,NVEC,&
                                                     NSLICES, NFIRST,&
                                                     EVINT,NEV,MLAN,MAXIT,NEVOUT
     integer                                      :: CHEBDEG, CHEBTYPE
     double precision                             :: LMIN, LMAX,THRE_INT,THRE_EXT
     double precision, dimension(1:4)             :: XINTV
     double precision, dimension(:), allocatable  :: SLI, EIGVAL, EIGVEC
     
     real*8                                       :: PI,tmperr,TOL
     real*8, allocatable                          :: inputv(:),outputav(:),outputbv(:) 
     integer,allocatable                          :: myeigid(:)

     real*8                                       :: npw,errd
     real*8, allocatable                          :: pv(:), pw(:),eigerr(:)

     integer*8                                    :: pevslAB

     ! save files
     integer                                      :: bfsiz
     integer(kind=mpi_offset_kind)                :: offset,leight
     real*8, allocatable                          :: vdat(:),eigtmp(:)
     character(len=1024)                          :: fname,fnb


     PI = 3.14159265359

     MDEG    = 3000
     NVEC    = 60
     NSLICES = 1
     NFIRST  = -1
     MLAN = 3000
     LANSTEP = 5000
    
     TOL = 1.0D-5
      
     call pEVSL_Start_F90(mymatvec%comm,pevslAB)
     
     ! SET PROB SIZE: NFIRST IS UNDEFINED
     call PEVSL_SETPROBSIZES_F90(pevslAB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)
     !call PEVSL_PARCSRCREATE_F90(  mymatvec%pbsiz_Glb,mymatvec%pbsiz_Glb,&
     !     mymatvec%B%idx,mymatvec%B%idx,mymatvec%B%rowloc,mymatvec%B%col,&
     !     mymatvec%B%val,mymatvec%comm,mymatvec%sBV) 


     !call PEVSL_SETAMV_F90(sparseBV, mymatvec)
     !CHEBDEG = 20
     !! NEED TO HAVE EIG BOUNDS OF B IN SETUP [DONE ABOVE BUT WILL DO AGAIN IN
     !! SETUP]

     !call PEVSL_SETBSOL_F90(pevslAB,solveBV, mymatvec)
     call PEVSL_SETBMV_F90(pevslAB, sparseBV, mymatvec)
     
     ! Ruipeng
     CHEBTYPE = 2
     call PEVSL_SETBSOL_CHEBITER_F90(pevslAB, CHEBTYPE, MYMATVEC%CHEBB)
     !call PEVSL_SETLTSOL_LSPOL_F90(pevslAB, MYMATVEC%LSPOLBSQRT)
     
     if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then    
        call PEVSL_SETAMV_F90(pevslAB, sparseAV, mymatvec)
     else
        call PEVSL_SETAMV_F90(pevslAB, sparsefsAV, mymatvec)
     endif 

     call PEVSL_SET_GENEIG_F90(pevslAB)
     
     call pEVSL_LANBOUNDS_F90(pevslAB, MLAN, LANSTEP,TOL, LMIN, LMAX)
     if (mymatvec%rank .eq. 0) then 
        print*, "step 0: eigenvalue bounds for B^{-1}A: lmin", LMIN, "lmax", LMAX
     endif 
    
     call report_time(pin%comm)
     !LMIN = 0.0D0
     !LMAX = 0.0D0

     XINTV(1) = (2.0D0*PI*pin%lowfreq)**2*1.0D-6 !-1.0D-5
     XINTV(2) = (2.0D0*PI*pin%upfreq)**2*1.0D-6 !1.0D-4

     !XINTV(1) = LMIN

     if (XINTV(1).lt.1.0D-10) then
        XINTV(1) = LMIN
     endif

     XINTV(3) = LMIN
     XINTV(4) = LMAX
     NSLICES = 1
     allocate (SLI (NSLICES+1))


     THRE_INT = 0.8D0 ! 0.25D0 0.8
     THRE_EXT = 0.7D0 ! 0.6

     if (.true.) then
     !! TODO
     EVINT = 800
     !
     call pEVSL_FINDPOL_F90(XINTV, THRE_INT, THRE_EXT, POL)
     NEV = EVINT + 2
     MLAN = MAX(4*NEV, 2000)
     !MLAN = MIN(MLAN, mymatvec%Gpbsiz)
     MAXIT = 3*MLAN
     !print*, MAXIT

     call pEVSL_CHEBLANNR_F90(pevslAB, XINTV, MAXIT, TOL, POL)

     call pEVSL_GET_NEV_F90(pevslAB,NEVOUT)
     !print *, "NEV FOUND", NEVOUT
 
     if (NEVOUT > 0) then
        allocate (EIGVAL (NEVOUT))
        allocate (EIGVEC (NEVOUT *mymatvec%pbsiz))
        call pEVSL_COPY_RESULT_F90(pevslAB, EIGVAL, EIGVEC, mymatvec%pbsiz)

        call pEVSL_chebiterstatsprint_f90(mymatvec%chebb)
        if (mymatvec%fsexist.or.mymatvec%purefluid) then 
           call pEVSL_chebiterstatsprint_f90(mymatvec%chebAp)
        endif
        ! --------------------------------------------------------------------------
        allocate(  inputv(mymatvec%pbsiz)) 
        allocate(outputav(mymatvec%pbsiz)) 
        allocate(outputbv(mymatvec%pbsiz)) 
 
        allocate(eigerr(NEVOUT)); eigerr = 0.0D0    
        allocate(eigtmp(NEVOUT)); eigtmp = EIGVAL  

        do J = 1,NEVOUT
            !! added 
            inputv = EIGVEC((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz)
            if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
               call sparseAV(inputv,outputav,mymatvec)
            else
               call sparsefsAV(inputv,outputav,mymatvec)
            !   !call pdfsavpdct(inputv,outputav,mymatvec)
            endif
            call sparseBV(inputv,outputbv,mymatvec)
            tmperr = sum((outputav-EIGVAL(J)*outputbv)**2)
            call mpi_reduce(tmperr,errd,1,mpi_real8,mpi_sum,0,mymatvec%comm,ierr)
            !tmperr = tmperr/dabs(EIGVAL(J))
            if (mymatvec%rank .eq. 0) then
               eigerr(J) = dsqrt(errd/real(mymatvec%Gpbsiz,8))/abs(EIGVAL(J))
               !print*, 'eigenvalues', EIGVAL(J), 'relative error', eigerr(J)
               !print*, tmperr, errd
            endif
        enddo
        

        allocate(myeigid(NEVOUT)); myeigid = (/(i,i=1,NEVOUT)/)  
        call ssort_real(EIGVAL,myeigid)
        !print*,eigtmp(myeigid) - eigval

        if (mymatvec%rank == 0) then
           print*, 'Relative errors: from the smallest to largest eigens:'
           do J = 1,NEVOUT
              i = myeigid(j)
              print*, 'Row', int(j,2), ' relative err. ', eigerr(i)
           enddo
           print*, '==============================================================='
           print*, 'Transform to frequencies (mHz), and periods (s)'
           do J = 1,NEVOUT
              if (EIGVAL(J) > 1.D-15) then 
                 print*, 'Row', int(J,2), dsqrt(EIGVAL(J))/PI/2.0*1.D3,&
                          2.0D0*PI/dsqrt(EIGVAL(J))
              else
                 print*, 'Row', int(J,2), EIGVAL(J), '...too small...'
              endif
           enddo
           print*, '================================================================'
        endif 

        allocate(vdat(mymatvec%Gpbsiz))
        leight = 8
        
        if (mymatvec%rank.eq.0) print*,"save eigenvectors"
        do i = 1,nevout
           j = myeigid(i)
           write(fnb,*) i; fnb = adjustl(fnb) 
           fname  = trim(pin%fvdata)//'_'//trim(fnb)//'.dat'   
           fname  = trim(fname) 
           offset = mymatvec%B%sizdist(mymatvec%rank+1)*leight 
           bfsiz  = mymatvec%pbsiz 
           vdat   = EIGVEC((j-1)*bfsiz+1:j*bfsiz)*mymatvec%B%diag 
           call pnm_save_dreal(bfsiz,vdat,offset,fname)
        enddo 

        call writelog("save eigenvectors",pin%comm)
        call report_time(pin%comm)

        DEALLOCATE (EIGVAL)
        DEALLOCATE (EIGVEC)

     endif

     ! save the solution
     call pnm_save_results()
 
     ! ---------------------------------------------------------------------------
     call pEVSL_FREEPOL_F90(POL)
     endif ! if true 

     DEALLOCATE (SLI)

     call pEVSL_FINISH_F90(pevslAB)

  end subroutine pnm_apply_pevsl


  subroutine pnm_save_results()
     integer                        :: i,j,ierr,bufsize,fid

     integer(kind=mpi_offset_kind)  :: disp,lfour,leight
     
     lfour = 4; leight = 8
     if(mymatvec%rank.eq.0) print*,"save the results" 
     ! save unstrM%new%vlist 
     !call pnm_save_int()
     disp = unstrM%new%vtxdist(unstrM%rank+1)*lfour
     bufsize = unstrM%new%nvtx

     !if (mymatvec%rank.eq.0) print*,unstrM%new%vlist
     !if (mymatvec%rank.eq.0) print*,unstrM%new%vtxdist

     call pnm_save_integer(bufsize,unstrM%new%vlist,disp,pin%fvlist)   
 
     if (mymatvec%fsexist) then 
        call pnm_save_integer(bufsize,unstrM%new%vstat,disp,pin%fvstat)   
        !call pnm_save_integer(bufsize,CGM%B%rstt,disp,pin%fvloct)   
     endif

  end subroutine

  subroutine pnm_save_int()
     integer                        :: i,j,ierr,bufsize,fid

     integer(kind=mpi_offset_kind)  :: disp,lfour
    
     lfour = 4

     call mpi_file_open(mymatvec%comm,pin%fvlist,&
           mpi_mode_wronly+mpi_mode_create,mpi_info_null,fid,ierr)

     disp = unstrM%new%vtxdist(unstrM%rank+1)*lfour

     call mpi_file_set_view(fid,disp,mpi_integer,mpi_integer,&
             'native',mpi_info_null,ierr) 

     bufsize = unstrM%new%nvtx
     call mpi_file_write(fid,unstrM%new%vlist,bufsize,mpi_integer,&
              mpi_status_ignore,ierr) 
     
     call mpi_file_close(fid,ierr) 

  end subroutine pnm_save_int


  subroutine pnm_save_integer(bufsize,buf,disp,fname)
     integer, intent(in)                        :: bufsize
     integer(kind=mpi_offset_kind), intent(in)  :: disp
     integer, intent(in)                        :: buf(bufsize) 
     character(len=1024), intent(in)            :: fname 

     integer                                    :: i,j,ierr,fid

     call mpi_file_open(mymatvec%comm,fname,&
           mpi_mode_wronly+mpi_mode_create,mpi_info_null,fid,ierr)


     call mpi_file_set_view(fid,disp,mpi_integer,mpi_integer,&
             'native',mpi_info_null,ierr) 

     call mpi_file_write(fid,buf,bufsize,mpi_integer,&
              mpi_status_ignore,ierr) 
     
     call mpi_file_close(fid,ierr) 

  end subroutine pnm_save_integer



  subroutine pnm_save_dreal(bufsize,buf,disp,fname)
     integer, intent(in)                        :: bufsize
     integer(kind=mpi_offset_kind), intent(in)  :: disp
     real*8, intent(in)                         :: buf(bufsize) 
     character(len=1024), intent(in)            :: fname 

     integer                                    :: i,j,ierr,fid

     call mpi_file_open(mymatvec%comm,fname,&
           mpi_mode_wronly+mpi_mode_create,mpi_info_null,fid,ierr)

     call mpi_file_set_view(fid,disp,mpi_real8,mpi_real8,&
             'native',mpi_info_null,ierr) 

     call mpi_file_write(fid,buf,bufsize,mpi_real8,&
              mpi_status_ignore,ierr) 
     
     call mpi_file_close(fid,ierr) 

  end subroutine pnm_save_dreal


  subroutine pnm_amv_test()
     integer                        :: i,j
     real(kind=rkind), allocatable  :: v0(:),w0(:)
    
     allocate(v0(mymatvec%A%siz)); v0 = 1.0D0 
     allocate(w0(mymatvec%A%siz)); w0 = 1.0D0 
     do i = 1,100000
        call sparseAV(v0,w0,mymatvec)
     enddo

  end subroutine pnm_amv_test

  subroutine pnm_bmv_test()
     integer                        :: i,j
     real(kind=rkind), allocatable  :: v0(:),w0(:)
    
     allocate(v0(mymatvec%B%siz)); v0 = 1.0D0 
     allocate(w0(mymatvec%B%siz)); w0 = 1.0D0 
     do i = 1,100000
        call sparseBV(v0,w0,mymatvec)
     enddo

  end subroutine pnm_bmv_test



end module pevsl_mod 
