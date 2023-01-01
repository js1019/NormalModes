module pevsl_mod

  use mpi
  use omp_lib
  use para_mod,                 only: rkind,pin
  use utility_mod
  use lapack_interface
  use geometry_mod,             only: unstrM
  use cg_matvec_mod


  implicit none 
  

 contains

  subroutine pnm_pevsl_control() 

     if (pin%JOB.le.3) then 
        call pnm_apply_pevsl()
     elseif (pin%JOB.eq.4) then 
        call pnm_pevsl_JOB4() 
     elseif (pin%JOB.ge.5.and.pin%JOB.le.8) then 
        call pnm_pevsl_rot()
     endif

    
  end subroutine pnm_pevsl_control


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
     real*8, allocatable                          :: vdat(:),eigtmp(:),eigs(:)
     character(len=1024)                          :: lab


     PI = 3.14159265359

     MDEG = 1000; NVEC = 60; NSLICES = 1
     NFIRST = -1; MLAN = 1000; LANSTEP = 3000
      
     call pEVSL_Start_F90(mymatvec%comm,pevslAB)
     
     ! SET PROB SIZE: NFIRST IS UNDEFINED
     call PEVSL_SETPROBSIZES_F90(pevslAB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)

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
     EVINT =pin%evint
     !
     call pEVSL_FINDPOL_F90(XINTV, THRE_INT, THRE_EXT, POL)
     NEV = EVINT + 2
     MLAN = MAX(4*NEV, 2000)
     !MLAN = MIN(MLAN, mymatvec%Gpbsiz)
     MAXIT = 3*MLAN
     !print*, MAXIT
     TOL = 1.0D8

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
            !tmperr = sum((outputav-EIGVAL(J)*outputbv)**2)
            !call mpi_reduce(tmperr,errd,1,mpi_real8,mpi_sum,0,mymatvec%comm,ierr)
            errd = 0.0D0; inputv = outputav-eigval(j)*outputbv
            call pnm_pvecnorm(inputv,mymatvec%pbsiz,errd)
            !tmperr = tmperr/dabs(EIGVAL(J))
            if (mymatvec%rank .eq. 0) then
               eigerr(J) = errd/real(mymatvec%Gpbsiz,8)/abs(EIGVAL(J))
               !print*, 'eigenvalues', EIGVAL(J), 'relative error', eigerr(J)
               !print*, tmperr, errd
            endif
        enddo
        

        allocate(myeigid(NEVOUT)); myeigid = (/(i,i=1,NEVOUT)/)  
        call ssort_real(EIGVAL,myeigid)
        !print*,eigtmp(myeigid) - eigval

        if (mymatvec%rank == 0) then
           print*, 'Relative errors: from the smallest to largest eigens:'
           allocate(eigs(NEVOUT)) ! save eigenfreq
           do J = 1,NEVOUT
              i = myeigid(j)
              print*, 'Row', int(j,2), ' relative err. ', eigerr(i)
           enddo
           print*, '==============================================================='
           print*, 'Transform to frequencies (mHz), and periods (s)'
           do J = 1,NEVOUT
              if (EIGVAL(J) > 1.D-15) then
                 eigs(J) =  dsqrt(EIGVAL(J))/PI/2.0D0*1.D3
                 print*, 'Row', int(J,2), dsqrt(EIGVAL(J))/PI/2.0*1.D3,&
                          2.0D0*PI/dsqrt(EIGVAL(J))
              else
                 eigs(J) = 0.0D0
                 print*, 'Row', int(J,2), EIGVAL(J), '...too small...'
              endif
           enddo
           print*, '================================================================'
        endif 

        ! save eigenfrequencies
        if (mymatvec%rank.eq.0) then
           print*,"save eigenvalues"
           !call pnm_save_freq(EIGVAL,NEVOUT)
           call pnm_save_freq(eigs,NEVOUT)
        endif  


        ! JS 06/29/19 
        if (pin%vecsave) then 
           if (mymatvec%rank.eq.0) print*,"save eigenvectors"
          allocate(vdat(mymatvec%Gpbsiz))
          bfsiz = mymatvec%pbsiz
          do j = 1,nevout  
             vdat   = EIGVEC((j-1)*bfsiz+1:j*bfsiz)*mymatvec%B%diag
             EIGVEC((j-1)*bfsiz+1:j*bfsiz) = vdat 
          enddo

          lab = '_';
          call pnm_save_data(lab,nevout,myeigid,eigvec) 
        endif 

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
  
    
  subroutine pnm_pevsl_JOB4() 
     integer                                      :: I, K, J, IERR
     integer*8                                    :: POL, CHEBSOL
     integer                                      :: LANSTEP,MDEG,NVEC,&
                                                     NSLICES, NFIRST,&
                                                     EVINT,NEV,MLAN,MAXIT,NEVOUT
     integer                                      :: CHEBDEG, CHEBTYPE
     double precision                             :: LMIN, LMAX,THRE_INT,THRE_EXT
     double precision, dimension(1:4)             :: XINTV
     double precision, dimension(:), allocatable  :: SLI, EIGVAL, EIGVEC
     double precision, allocatable                :: lanvecs(:), a20(:,:)
     double precision, allocatable                :: a2(:,:), b2(:,:)
     
     real*8                                       :: PI,tmperr,TOL
     real*8, allocatable                          :: inputv(:),outputav(:),outputbv(:) 
     real*8, allocatable                          :: inputv1(:),inputv2(:) 
     integer,allocatable                          :: myeigid(:)
     real*8                                       :: npw,st0,rv0,errd
     real*8, allocatable                          :: pv(:), pw(:),eigerr(:)
     integer*8                                    :: pevslAB

     ! call dsyev 
     double precision, dimension(:), allocatable  :: wev0, wev1, eigs      

     ! save files
     double precision, allocatable                :: vecs0(:,:), eigvec0(:)
     character(len=1024)                          :: lab
     

     PI = 3.14159265359

     MDEG = 1000; NVEC = 60; NSLICES = 1
     NFIRST = -1; MLAN = 1000; LANSTEP = 3000

     call pEVSL_Start_F90(mymatvec%comm,pevslAB)
     
     ! SET PROB SIZE: NFIRST IS UNDEFINED
     call PEVSL_SETPROBSIZES_F90(pevslAB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)
    
     call PEVSL_SETBMV_F90(pevslAB, sparseBV, mymatvec)
     
     ! Ruipeng
     CHEBTYPE = 2
     call PEVSL_SETBSOL_CHEBITER_F90(pevslAB, CHEBTYPE, MYMATVEC%CHEBB)
     !call PEVSL_SETLTSOL_LSPOL_F90(pevslAB, MYMATVEC%LSPOLBSQRT)
     
     if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then    
        call PEVSL_SETAMV_F90(pevslAB, sparseAVnoP, mymatvec)
     else
        call PEVSL_SETAMV_F90(pevslAB, sparsefsAVnoP, mymatvec)
     endif 

     call PEVSL_SET_GENEIG_F90(pevslAB)
     
     call pEVSL_LANBOUNDS_F90(pevslAB, MLAN, LANSTEP,TOL, LMIN, LMAX)
     if (mymatvec%rank .eq. 0) then 
        print*, "step 0: eigenvalue bounds for B^{-1}A: lmin", LMIN, "lmax", LMAX
     endif 
    
     call report_time(pin%comm)

     XINTV(1) = (2.0D0*PI*pin%lowfreq)**2*1.0D-6 !-1.0D-5
     XINTV(2) = (2.0D0*PI*pin%upfreq)**2*1.0D-6 !1.0D-4

     !XINTV(1) = LMIN

     if (XINTV(1).lt.1.0D-10) then
        XINTV(1) = LMIN
     endif

     XINTV(3) = LMIN; XINTV(4) = LMAX
     NSLICES = 1
     allocate (SLI (NSLICES+1))

     THRE_INT = 0.8D0 ! 0.25D0 0.8
     THRE_EXT = 0.7D0 ! 0.6

     !! TODO
     EVINT = pin%evint
     !
     call pEVSL_FINDPOL_F90(XINTV, THRE_INT, THRE_EXT, POL)
     NEV = EVINT + 2
     MLAN = MAX(4*NEV, 2000)
     !MLAN = MIN(MLAN, mymatvec%Gpbsiz)
     MAXIT = 3*MLAN
     !print*, MAXIT
     TOL = 1.0D-8

     call pEVSL_lanvectors_F90(pevslAB, XINTV, MAXIT, TOL, POL)

     call pEVSL_GET_NEV_F90(pevslAB,NEVOUT)
     call pEVSL_chebiterstatsprint_f90(mymatvec%chebb)
     if (mymatvec%fsexist.or.mymatvec%purefluid) then 
        call pEVSL_chebiterstatsprint_f90(mymatvec%chebAp)
     endif

     if (mymatvec%rank.eq.0) print*, "Lanczos vectors done"
     call report_time(mymatvec%comm)

     ! compute the matrix 
     allocate(lanvecs(nevout*mymatvec%pbsiz))
     call pEVSL_copy_vectors_f90(pevslAB,lanvecs,mymatvec%pbsiz)
     allocate(vecs0(mymatvec%pbsiz,nevout))
     do i = 1,nevout
        vecs0(:,i) = lanvecs((i-1)*mymatvec%pbsiz+1:i*mymatvec%pbsiz)
     enddo

     ! allocate matrix
     allocate(a2(nevout,nevout),b2(nevout,nevout),a20(nevout,nevout))

     allocate( inputv1(mymatvec%pbsiz)) 
     allocate( inputv2(mymatvec%pbsiz)) 
     allocate(outputav(mymatvec%pbsiz)) 
     allocate(outputbv(mymatvec%pbsiz)) 

     do i = 1,nevout
        inputv1 = lanvecs((i-1)*mymatvec%pbsiz+1:i*mymatvec%pbsiz)
        if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
           call sparseAV(inputv1,outputav,mymatvec)
        else
           call sparsefsAV(inputv1,outputav,mymatvec)
        endif
        if (mymatvec%rank.eq.0.and. mod(i,30).eq.0.and.pin%phi1) then
            print*, "matrix generation... i = ", i
            call report_time(mymatvec%comm)
        endif

        do j = 1,i
           ! for matrix A2
           inputv2 = lanvecs((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz)
           st0 = sum(inputv2*outputav); rv0 = 0.0D0;
           call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           a2(i,j) = rv0; a2(j,i) = a2(i,j)
           !! check matrix M
           !call sparseBV(inputv1,outputbv,mymatvec)
           !st0 = sum(inputv2*outputbv); rv0 = 0.0D0; 
           !call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           !b2(i,j) = rv0; b2(j,i) = b2(i,j)
        enddo
     enddo
     if (mymatvec%rank.eq.0) then 
        print*, "matrix (of size ", nevout, ") is done"
     endif
     call report_time(mymatvec%comm)
     !if (mymatvec%rank.eq.0) print*, b2

     allocate(wev0(nevout)); wev0 = 0.0D0; 
     a20 = a2
     call intface_dsyev(nevout,a2,wev0)
     !if (mymatvec%rank.eq.0) print*, matmul(a2,transpose(a2))


     !if (mymatvec%rank.eq.0) print*, wev
     allocate(wev1(nevout)); wev1 = max(wev0,0.0D0); 
     NEV = 0; 
     do i = 1,nevout 
        if (wev1(i).ge.XINTV(1).and.wev1(i).le.XINTV(2)) then 
           NEV = NEV + 1
        endif 
     enddo
     allocate(eigval(nev),eigvec0(nev*mymatvec%pbsiz),eigerr(nev))
     allocate(myeigid(nev)); myeigid = (/(i,i=1,nev)/)  
     j = 0
     do i = 1,nevout 
        if (wev1(i).ge.XINTV(1).and.wev1(i).le.XINTV(2)) then 
           j = j + 1
           eigval(j) = wev0(i)
           ! JS 06/30/19 check errors
           inputv1 = matmul(vecs0,a2(:,i))
           call sparseBV(inputv1,outputbv,mymatvec)
           !npw = 0.0D0
           !call pnm_pvecnorm(outputbv,mymatvec%pbsiz,npw) 
           st0 = sum(inputv1*outputbv); rv0 = 0.0D0;
           call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           npw = rv0
           inputv1  = inputv1/npw
           outputbv = outputbv/npw             
           !if(mymatvec%rank.eq.0) print*, npw, j
           !inputv1 = matmul(vecs0,a2(j,:))
           if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
              call sparseAV(inputv1,outputav,mymatvec)
           else
              call sparsefsAV(inputv1,outputav,mymatvec)
           endif
           !tmperr = sum((outputav-eigval(j)*outputbv)**2)
           !call mpi_reduce(tmperr,errd,1,mpi_real8,mpi_sum,0,mymatvec%comm,ierr)
           errd = 0.0D0; inputv2 = outputav-eigval(j)*outputbv
           call pnm_pvecnorm(inputv2,mymatvec%pbsiz,errd)
           if (mymatvec%rank.eq.0) then
              eigerr(j) = errd/dsqrt(real(mymatvec%Gpbsiz,8))/dsqrt(abs(eigval(j)))
           endif
           eigvec0((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz) =&
           matmul(vecs0,a2(:,i))*mymatvec%B%diag/npw

        endif 
     enddo
     !endif

     if (mymatvec%rank == 0) then
        print*, 'Relative errors:'
        do J = 1,NEV
           print*, 'Row', int(j,2), ' relative err. 1/omega^2',&
                eigerr(j)/dsqrt(abs(eigval(j)))
        enddo
        do J = 1,NEV
           print*, 'Row', int(j,2), ' relative err. 1/omega  ', eigerr(j)
        enddo
        allocate(eigs(NEV)) ! save eigenfreq
        print*, '==============================================================='
        print*, 'Transform to frequencies (mHz), and periods (s)'
        do J = 1,NEV
           if (EIGVAL(J).gt.1.D-15) then 
              eigs(J) =  dsqrt(EIGVAL(J))/PI/2.0D0*1.D3
              print*, 'Row', int(J,2), dsqrt(EIGVAL(J))/PI/2.0*1.D3,&
                       2.0D0*PI/dsqrt(EIGVAL(J))
           else
              eigs(J) = 0.0D0
              print*, 'Row', int(J,2), EIGVAL(J), '...too small...'
           endif
        enddo
        print*, '================================================================'
     endif 
    
     ! save data    
     if (mymatvec%rank.eq.0) then 
        print*,"save eigenvalues"
        !call pnm_save_freq(EIGVAL,NEV)
        call pnm_save_freq(eigs,NEV)
     endif  

     ! JS 06/29/19 
     if (pin%vecsave) then 
        if (mymatvec%rank.eq.0) print*,"save eigenvectors"
        lab = '_' 
        call pnm_save_data(lab,nev,myeigid,eigvec0) 
     endif 
 
     ! save the solution
     call pnm_save_results()
 
     DEALLOCATE(EIGVAL,myeigid)
     DEALLOCATE(lanvecs,eigvec0,vecs0)


     call pEVSL_FREEPOL_F90(POL)
     DEALLOCATE (SLI)
     call pEVSL_FINISH_F90(pevslAB)

  end subroutine pnm_pevsl_JOB4



  subroutine pnm_pevsl_rot() 
     integer                                      :: I, K, J, IERR
     integer*8                                    :: POL, CHEBSOL
     integer                                      :: LANSTEP,MDEG,NVEC,&
                                                     NSLICES, NFIRST,&
                                                     EVINT,NEV,MLAN,MAXIT,NEVOUT
     integer                                      :: CHEBDEG, CHEBTYPE
     double precision                             :: LMIN, LMAX,THRE_INT,THRE_EXT
     double precision, dimension(1:4)             :: XINTV
     double precision, dimension(:), allocatable  :: SLI, EIGVAL, EIGVEC
     double precision, allocatable                :: lanvecs(:), dg(:), ofd(:)
     real*8, dimension(:), allocatable            :: wev, wev1, wev2, wev3

     real*8, dimension(:,:), allocatable          :: CC, Co, Co0, Co1, Co2
     complex*16, allocatable                      :: Zo(:,:)
     
     real*8                                       :: PI,tmperr,TOL
     real*8, allocatable                          :: inputv1(:),inputv2(:) 
     real*8, allocatable                          :: outputav(:),outputbv(:) 
     real*8, allocatable                          :: outputrv(:) 
     integer,allocatable                          :: myeigid(:),sel(:),find(:)
     real*8                                       :: npw,npwr,npwi,st0,rv0,errd
     real*8, allocatable                          :: pv(:),pw(:),pv1(:),pw1(:),eigerr(:)
     integer*8                                    :: pevslAB
     integer                                      :: nev0     

     ! save files
     double precision, allocatable                :: vecs0(:,:),vecs1(:,:),&
                                                     eigr(:),eigi(:),eigs(:)
     character(len=1024)                          :: lab

     PI = 3.14159265359

     MDEG = 1000; NVEC = 60
     NSLICES = 1; NFIRST  = -1
     MLAN = 1000; LANSTEP = 3000

     call pEVSL_Start_F90(mymatvec%comm,pevslAB)
     
     ! SET PROB SIZE: NFIRST IS UNDEFINED
     call PEVSL_SETPROBSIZES_F90(pevslAB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)
    
     call PEVSL_SETBMV_F90(pevslAB, sparseBV, mymatvec)
     
     ! Ruipeng
     CHEBTYPE = 2
     call PEVSL_SETBSOL_CHEBITER_F90(pevslAB, CHEBTYPE, MYMATVEC%CHEBB)
     !call PEVSL_SETLTSOL_LSPOL_F90(pevslAB, MYMATVEC%LSPOLBSQRT)
    
     if (pin%JOB.ne.8) then  
        if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then    
           call PEVSL_SETAMV_F90(pevslAB, sparseAV, mymatvec)
        else
           call PEVSL_SETAMV_F90(pevslAB, sparsefsAV, mymatvec)
        endif 
     else
        if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then    
           call PEVSL_SETAMV_F90(pevslAB, sparseAVnoP, mymatvec)
        else
           call PEVSL_SETAMV_F90(pevslAB, sparsefsAVnoP, mymatvec)
        endif 
     endif

     call PEVSL_SET_GENEIG_F90(pevslAB)
     
     call pEVSL_LANBOUNDS_F90(pevslAB, MLAN, LANSTEP,TOL, LMIN, LMAX)
     if (mymatvec%rank .eq. 0) then 
        print*, "step 0: eigenvalue bounds for B^{-1}A: lmin", LMIN, "lmax", LMAX
     endif 
    
     call report_time(pin%comm)

     XINTV(1) = (2.0D0*PI*pin%lowfreq)**2*1.0D-6 !-1.0D-5
     XINTV(2) = (2.0D0*PI*pin%upfreq)**2*1.0D-6 !1.0D-4
     !XINTV(1) = LMIN

     if (XINTV(1).lt.1.0D-10) then
        XINTV(1) = LMIN
     endif

     XINTV(3) = LMIN; XINTV(4) = LMAX
     NSLICES = 1; allocate (SLI (NSLICES+1))

     THRE_INT = 0.8D0 ! 0.25D0 0.8
     THRE_EXT = 0.7D0 ! 0.6

     !! TODO
     EVINT = pin%evint
     !
     call pEVSL_FINDPOL_F90(XINTV, THRE_INT, THRE_EXT, POL)
     NEV = EVINT + 2; MLAN = MAX(4*NEV, 2000); MAXIT = 3*MLAN
     !print*, MAXIT
     TOL = 1.0D-8

     call pEVSL_lanvectors_F90(pevslAB, XINTV, MAXIT, TOL, POL)
     call pEVSL_GET_NEV_F90(pevslAB,nev0)
     call pEVSL_chebiterstatsprint_f90(mymatvec%chebb)
     if (mymatvec%fsexist.and..not.mymatvec%purefluid) then 
        call pEVSL_chebiterstatsprint_f90(mymatvec%chebAp)
     endif
     if (mymatvec%rank.eq.0) print*, "Lanczos vectors are done"
     call report_time(mymatvec%comm)
    
     ! compute the copy vectors
     allocate(lanvecs(nev0*mymatvec%pbsiz))
     call pEVSL_copy_vectors_f90(pevslAB,lanvecs,mymatvec%pbsiz)


     ! compute the dense matrix
     allocate(CC(nev0,nev0))
     allocate( inputv1(mymatvec%pbsiz)) 
     allocate( inputv2(mymatvec%pbsiz)) 
     allocate(outputav(mymatvec%pbsiz)) 
     allocate(outputbv(mymatvec%pbsiz)) 
     allocate(outputrv(mymatvec%pbsiz)) 
     do i = 1,nev0
        inputv1 = lanvecs((i-1)*mymatvec%pbsiz+1:i*mymatvec%pbsiz)
        if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
           call sparseAV(inputv1,outputav,mymatvec)
        else
           call sparsefsAV(inputv1,outputav,mymatvec)
        endif
        if (mymatvec%rank.eq.0.and.mod(i,30).eq.1.and.pin%phi1) then
            print*, "matrix generation... i = ", i
            call report_time(mymatvec%comm)
        endif
        do j = 1,i 
           inputv2 = lanvecs((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz)
           st0 = sum(inputv2*outputav); rv0 = 0.0D0;
           call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           CC(i,j) = rv0; CC(j,i) = rv0
        enddo
     enddo


     if (mymatvec%rank.eq.0) then 
        print*, "matrix of size, ", nev0, "is done"
     endif
     call report_time(mymatvec%comm)
     
     allocate(wev(nev0),Co(nev0,nev0)) 
     Co = CC
     call intface_dsyev(nev0,Co,wev)

     ! JS 07/24/2019 check eigenfreq
     allocate(sel(nev0)); sel = 0; j = 1 
     allocate(vecs0(mymatvec%pbsiz,nev0))
     do i = 1,nev0
        vecs0(:,i) = lanvecs((i-1)*mymatvec%pbsiz+1:i*mymatvec%pbsiz)
     enddo
     do i = 1,nev0
        inputv1 = matmul(vecs0,Co(:,i))
        if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
           call sparseAV(inputv1,outputav,mymatvec)
        else
           call sparsefsAV(inputv1,outputav,mymatvec)
        endif
        call sparseBV(inputv1,outputbv,mymatvec)
        inputv2 = outputav - wev(i)*outputbv 
        call pnm_pvecnorm(inputv2,mymatvec%pbsiz,errd)
        errd = dsqrt(errd**2/real(mymatvec%Gpbsiz,8))
        !if (mymatvec%rank.eq.0) then
        !   if (wev(i).gt.1.D-10) then 
        !      print*, i, dsqrt(dabs(wev(i)))/2.D0/PI*1.D3,'err', errd/dabs(wev(i))
        !   else
        !      print*, i, wev(i),'err', errd/dabs(wev(i))
        !   endif
        !endif
        if (wev(i).gt.1.D-10.and.dsqrt(dabs(wev(i)))/2.D0/PI*1.D3.gt.0.02&
           .and.errd/dabs(wev(i)).lt.1.D-7.and.pin%JOB.ne.8) then
           sel(i) = j; j = j + 1
        elseif (wev(i).gt.1.D-10.and.dsqrt(dabs(wev(i)))/2.D0/PI*1.D3.gt.0.02&
           .and.errd/dabs(wev(i)).lt.1.D-1.and.pin%JOB.eq.8) then
           sel(i) = j; j = j + 1
        endif 
     enddo 
     if (mymatvec%rank.eq.0) print*,"number of Ritz vectors", j-1

     ! build the complex matrix
     nevout = maxval(sel)
     allocate(find(nevout))
     do i = 1,nev0
        if (sel(i).ne.0) then
           find(sel(i)) = i
        endif 
     enddo 
     allocate(wev1(nevout)); 
     do i = 1,nev0
        if (sel(i).ne.0) then 
           wev1(sel(i)) = dsqrt(wev(i))
        endif
     enddo
     !if (mymatvec%rank.eq.0) then 
     !   print*, "check eigenfrequencies", wev1/PI/2.0D0*1.D3
     !endif
     allocate(vecs1(mymatvec%pbsiz,nevout))
     do i = 1,nevout
        vecs1(:,i) = matmul(vecs0,Co(:,find(i)))
     enddo
     
     allocate(Co0(nevout,nevout), Co1(nevout,nevout)); Co0 = 0.0D0
     allocate(Co2(nevout,nevout)); Co2 = 0.0D0
     do i = 1,nevout
        Co1(i,i) = wev1(i)
     enddo 
     !if(mymatvec%rank.eq.0) print*, nevout,Co1(nevout,:)
     allocate(Zo(2*nevout,2*nevout)); Zo = (0.0D0,0.0D0)
     Zo(1:nevout,nevout+1:2*nevout) = cmplx(Co1,Co0,16) 
     Zo(nevout+1:2*nevout,1:nevout) = cmplx(transpose(Co1),Co0,16) 
         
     do i = 1,nevout
        !inputv1 = lanvecs((find(i)-1)*mymatvec%pbsiz+1:find(i)*mymatvec%pbsiz)
        inputv1 = vecs1(:,i)
        call sparseRTV(inputv1,outputrv,mymatvec)
        do j = 1,i
           !inputv2 = lanvecs((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz)
           inputv2 = vecs1(:,j)
           st0 = sum(inputv2*outputrv); rv0 = 0.0D0;
           call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           Co2(j,i) = rv0
           if (j.ne.i) then
              Co2(i,j) = - Co2(j,i)
           endif
        enddo
     enddo
     Zo(nevout+1:2*nevout,nevout+1:2*nevout) = cmplx(Co0,2.0D0*Co2,16)
     
     if (mymatvec%rank.eq.0) then 
        print*, "complex matrix (of size ", 2*nevout, ") is done"
     endif
     call report_time(mymatvec%comm)
    
     ! solve the reduced system
     allocate(wev2(2*nevout),wev3(nevout)) 
     call intface_zheev(2*nevout,Zo,wev2)
     !print*, wev2(nevout+1:2*nevout)
     !if (mymatvec%rank.eq.0) print*,Zo(nevout+1:2*nevout,1) 
     NEV = 0; 
     do i = 1,nevout
        wev3 = wev2(nevout+1:2*nevout)**2 
        if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
           NEV = NEV + 1
        endif 
     enddo
     allocate(eigval(nev),eigerr(nev))
     allocate(eigr(nev*mymatvec%pbsiz),eigi(nev*mymatvec%pbsiz))
     allocate(myeigid(nev)); myeigid = (/(i,i=1,nev)/)  
     j = 0
     do i = 1,nevout 
        if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
           j = j + 1
           eigval(j) = dsqrt(wev3(i))
        endif
     enddo
    
     !if (mymatvec%rank.eq.0) print*,norm2(wev1-wev2(nevout+1:2*nevout))

     wev2 = max(1.0D-30,wev2)  ! make sure it is positive 
     allocate( pv(nevout), pw(nevout))
     allocate(pv1(nevout),pw1(nevout))
     j = 0
     do i = 1,nevout 
        if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
           j = j + 1
           pv =  real(Zo(nevout+1:2*nevout,i+nevout))/wev2(i+nevout)
           pw = aimag(Zo(nevout+1:2*nevout,i+nevout))/wev2(i+nevout)
           !pv1 =  real(Zo(1:nevout,i+nevout))/wev1(1:nevout)
           !pw1 = aimag(Zo(1:nevout,i+nevout))/wev1(1:nevout)
           !if (mymatvec%rank.eq.0) print*, norm2(pv-pv1),norm2(pw-pw1)
            
           ! JS 07/01/19
           ! check norm 
           npwr = 0.0D0; npwi = 0.0D0; npw = 0.0D0
           inputv1 = matmul(vecs1, pv)
           inputv2 = matmul(vecs1, pw)          
           call pnm_pvecnorm(inputv1,mymatvec%pbsiz,npwr) 
           call pnm_pvecnorm(inputv2,mymatvec%pbsiz,npwi) 
           !if (mymatvec%rank.eq.0) print*,i,npwi,npwr           
           !npwr = 0.0D0; npwi = 0.0D0; npw = 0.0D0

           call sparseBV(inputv1,outputbv,mymatvec)
           npwr = npwr+sum(inputv1*outputbv)
           npwi = npwi+sum(inputv2*outputbv) 
           call sparseBV(inputv2,outputbv,mymatvec)
           npwr = npwr+sum(inputv2*outputbv) 
           npwi = npwi-sum(inputv1*outputbv) 
           call sparseRTV(inputv1,outputrv,mymatvec)
           npwr = npwr-eigval(j)*sum(inputv2*outputrv) 
           npwi = npwi-eigval(j)*sum(inputv1*outputrv) 
           call sparseRTV(inputv2,outputrv,mymatvec)
           npwr = npwr-eigval(j)*sum(inputv1*outputrv) 
           npwi = npwi-eigval(j)*sum(inputv2*outputrv) 
           call mpi_allreduce(npwi,npw,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           npwi = npw; npw = 0.0D0;
           call mpi_allreduce(npwr,npw,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
           npwr = npw
           !if (mymatvec%rank.eq.0) print*,i,npwi/npwr           
  
           ! imagery 
           !pw = aimag(Zo(nevout+1:2*nevout,j))/wev2(i)
           inputv1 = matmul(vecs1,pw)           
           call sparseBV(inputv1,outputbv,mymatvec)
           if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
              call sparseAV(inputv1,outputav,mymatvec)
           else
              call sparsefsAV(inputv1,outputav,mymatvec)
           endif
           !pv =  real(Zo(1:nevout,j))/wev1(1:nevout)
           inputv1 = matmul(vecs1,pv)
           call sparseRTV(inputv1,outputrv,mymatvec)
           inputv2 = outputav - eigval(j)**2*outputbv + 2.0D0*eigval(j)*outputrv
           !inputv2 = - eigval(j)**2*outputbv
           npwi = 0.0D0
           call pnm_pvecnorm(inputv2,mymatvec%pbsiz,npwi)
           !if(mymatvec%rank.eq.0) print*,j,npwi**2/npw,npw

           ! real 
           !pv =  real(Zo(1:nevout,j))/wev1(1:nevout)
           inputv1 = matmul(vecs1,pv)
           call sparseBV(inputv1,outputbv,mymatvec)
           if (.not.mymatvec%fsexist.and..not.mymatvec%purefluid) then 
              call sparseAV(inputv1,outputav,mymatvec)
           else
              call sparsefsAV(inputv1,outputav,mymatvec)
           endif
           !pw = aimag(Zo(nevout+1:2*nevout,j))/wev2(i)
           inputv1 = matmul(vecs1,pw)           
           call sparseRTV(inputv1,outputrv,mymatvec)
           inputv2 = outputav - eigval(j)**2*outputbv - 2.0D0*eigval(j)*outputrv 
           npwr = 0.0D0
           call pnm_pvecnorm(inputv2,mymatvec%pbsiz,npwr)
           !if(mymatvec%rank.eq.0) print*,j,npwr**2/npw,npwi**2/npw,npw

           errd = (npwi**2 + npwr**2)/npw 
           !errd = npwr**2/npw 
           if (mymatvec%rank.eq.0) then
              eigerr(j) = dsqrt(errd/real(mymatvec%Gpbsiz,8))/dabs(eigval(j))
           endif
           eigr((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz) =mymatvec%B%diag*&
                          matmul(vecs1,pv)/dsqrt(npw)
           eigi((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz) =mymatvec%B%diag*&
                          matmul(vecs1,pw)/dsqrt(npw)
        endif
     enddo
     
     if (mymatvec%rank == 0) then
        print*, 'Relative errors:'
        do J = 1,NEV
           print*, 'Row', int(j,2), ' relative err. ', eigerr(j)
        enddo
        allocate(eigs(NEV)) ! save eigenfreq
        print*, '==============================================================='
        print*, 'Transform to frequencies (mHz), and periods (s)'
        do J = 1,NEV
           if (EIGVAL(J).gt.1.D-15) then 
              eigs(J) =  EIGVAL(J)/PI/2.0D0*1.D3
              print*, 'Row', int(J,2), EIGVAL(J)/PI/2.0*1.D3,&
                       2.0D0*PI/EIGVAL(J)
           else
              eigs(J) = 0.0D0
              print*, 'Row', int(J,2), EIGVAL(J), '...too small...'
           endif
        enddo
        print*, '================================================================'
     endif 
     
     ! save data    
     if (mymatvec%rank.eq.0) then 
        print*,"save eigenvalues"
        !call pnm_save_freq(EIGVAL,NEV)
        call pnm_save_freq(eigs,NEV)
     endif  

     ! JS 06/29/19 
     if (pin%vecsave) then 
        if (mymatvec%rank.eq.0) print*,"save eigenvectors"
        lab = '_real_' 
        call pnm_save_data(lab,nev,myeigid,eigr) 
        lab = '_imag_' 
        call pnm_save_data(lab,nev,myeigid,eigi) 
     endif 
 
     ! save the solution
     call pnm_save_results()
 
     DEALLOCATE(EIGVAL,myeigid)
     DEALLOCATE(lanvecs,eigr,eigi,vecs0)

     call pEVSL_FREEPOL_F90(POL)
     DEALLOCATE (SLI)
     call pEVSL_FINISH_F90(pevslAB)

  end subroutine pnm_pevsl_rot


  !subroutine pnm_pevsl_lanrot() 
  !   integer                                      :: I, K, J, IERR
  !   integer*8                                    :: POL, CHEBSOL
  !   integer                                      :: LANSTEP,MDEG,NVEC,&
  !                                                   NSLICES, NFIRST,&
  !                                                   EVINT,NEV,MLAN,&
  !                                                   MAXIT,NEVOUT0,NEVOUT
  !   integer                                      :: CHEBDEG, CHEBTYPE
  !   double precision                             :: LMIN, LMAX,THRE_INT,THRE_EXT
  !   double precision, dimension(1:4)             :: XINTV
  !   double precision, dimension(:), allocatable  :: SLI, EIGVAL, EIGVEC
  !   double precision, allocatable                :: lanvecs(:), dg(:), ofd(:)
  !   real*8, dimension(:), allocatable            :: wev, wev1, wev2, wev3

  !   real*8, dimension(:,:), allocatable          :: CC, Co, Co0, Co1, Co2
  !   complex*16, allocatable                      :: Zo(:,:)
  !   
  !   real*8                                       :: PI,tmperr,TOL,eigcut
  !   real*8, allocatable                          :: inputv1(:),inputv2(:) 
  !   real*8, allocatable                          :: outputav(:),outputbv(:) 
  !   real*8, allocatable                          :: outputrv(:)
  !   integer                                      :: eigs 
  !   integer,allocatable                          :: myeigid(:)
  !   real*8                                       :: npw,npwr,npwi,st0,rv0,errd
  !   real*8, allocatable                          :: pv(:),pw(:),eigerr(:)
  !   integer*8                                    :: pevslAB

  !   ! save files
  !   double precision, allocatable                :: vecs0(:,:),eigr(:),eigi(:)
  !   character(len=1024)                          :: lab

  !   PI = 3.14159265359

  !   MDEG = 1000; NVEC = 60
  !   NSLICES = 1; NFIRST  = -1
  !   MLAN = 1000; LANSTEP = 3000

  !   call pEVSL_Start_F90(mymatvec%comm,pevslAB)
  !   
  !   ! SET PROB SIZE: NFIRST IS UNDEFINED
  !   call PEVSL_SETPROBSIZES_F90(pevslAB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)
  !  
  !   call PEVSL_SETBMV_F90(pevslAB, sparseBV, mymatvec)
  !   
  !   ! Ruipeng
  !   CHEBTYPE = 2
  !   call PEVSL_SETBSOL_CHEBITER_F90(pevslAB, CHEBTYPE, MYMATVEC%CHEBB)
  !   !call PEVSL_SETLTSOL_LSPOL_F90(pevslAB, MYMATVEC%LSPOLBSQRT)
  !  
  !   if (pin%JOB.ne.12) then  
  !      if (.not.mymatvec%fsexist) then    
  !         call PEVSL_SETAMV_F90(pevslAB, sparseAV, mymatvec)
  !      else
  !         call PEVSL_SETAMV_F90(pevslAB, sparsefsAV, mymatvec)
  !      endif 
  !   else
  !      if (.not.mymatvec%fsexist) then    
  !         call PEVSL_SETAMV_F90(pevslAB, sparseAVnoP, mymatvec)
  !      else
  !         call PEVSL_SETAMV_F90(pevslAB, sparsefsAVnoP, mymatvec)
  !      endif 
  !   endif

  !   call PEVSL_SET_GENEIG_F90(pevslAB)
  !   
  !   call pEVSL_LANBOUNDS_F90(pevslAB, MLAN, LANSTEP,TOL, LMIN, LMAX)
  !   if (mymatvec%rank .eq. 0) then 
  !      print*, "step 0: eigenvalue bounds for B^{-1}A: lmin", LMIN, "lmax", LMAX
  !   endif 
  !  
  !   call report_time(pin%comm)

  !   XINTV(1) = (2.0D0*PI*pin%lowfreq)**2*1.0D-6 !-1.0D-5
  !   XINTV(2) = (2.0D0*PI*pin%upfreq)**2*1.0D-6 !1.0D-4
  !   !XINTV(1) = LMIN

  !   if (XINTV(1).lt.1.0D-10) then
  !      XINTV(1) = LMIN
  !   endif

  !   XINTV(3) = LMIN; XINTV(4) = LMAX
  !   NSLICES = 1; allocate (SLI (NSLICES+1))

  !   THRE_INT = 0.8D0 ! 0.25D0 0.8
  !   THRE_EXT = 0.7D0 ! 0.6

  !   !! TODO
  !   EVINT = 500
  !   !
  !   call pEVSL_FINDPOL_F90(XINTV, THRE_INT, THRE_EXT, POL)
  !   NEV = EVINT + 2; MLAN = MAX(4*NEV, 2000); MAXIT = 3*MLAN
  !   !print*, MAXIT
  !   TOL = 1.0D8

  !   call pEVSL_lanvectors_F90(pevslAB, XINTV, MAXIT, TOL, POL)
  !   call pEVSL_GET_NEV_F90(pevslAB,NEVOUT0)
  !   call pEVSL_chebiterstatsprint_f90(mymatvec%chebb)
  !   if (mymatvec%fsexist) then 
  !      call pEVSL_chebiterstatsprint_f90(mymatvec%chebAp)
  !   endif
  !   if (mymatvec%rank.eq.0) print*, "Lanczos vectors are done"
  !   call report_time(mymatvec%comm)
  !  
  !   ! compute the copy vectors
  !   allocate(lanvecs(nevout0*mymatvec%pbsiz))
  !   call pEVSL_copy_vectors_f90(pevslAB,lanvecs,mymatvec%pbsiz)


  !   ! compute the dense matrix
  !   allocate(CC(nevout0,nevout0))
  !   allocate( inputv1(mymatvec%pbsiz)) 
  !   allocate( inputv2(mymatvec%pbsiz)) 
  !   allocate(outputav(mymatvec%pbsiz)) 
  !   allocate(outputbv(mymatvec%pbsiz)) 
  !   allocate(outputrv(mymatvec%pbsiz)) 
  !   do i = 1,nevout0 
  !      inputv1 = lanvecs((i-1)*mymatvec%pbsiz+1:i*mymatvec%pbsiz)
  !      if (.not.mymatvec%fsexist) then 
  !         call sparseAV(inputv1,outputav,mymatvec)
  !      else
  !         call sparsefsAV(inputv1,outputav,mymatvec)
  !      endif
  !      if (mymatvec%rank.eq.0.and.mod(i,30).eq.1.and.pin%phi1) then
  !          print*, "matrix generation... i = ", i
  !          call report_time(mymatvec%comm)
  !      endif
  !      do j = 1,i 
  !         inputv2 = lanvecs((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz)
  !         st0 = sum(inputv2*outputav); rv0 = 0.0D0;
  !         call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
  !         CC(i,j) = rv0; CC(j,i) = rv0
  !      enddo
  !   enddo


  !   if (mymatvec%rank.eq.0) then 
  !      print*, "matrix is done with size", nevout0 
  !   endif
  !   call report_time(mymatvec%comm)

  !   
  !   allocate(wev(nevout0),Co(nevout0,nevout0)) 
  !   Co = CC
  !   call intface_dsyev(nevout0,Co,wev)

  !   ! build the complex matrix 
  !   allocate(wev1(nevout0)); wev1 = dsqrt(max(1.0D-50,wev));
  !   if (mymatvec%rank.eq.0) then 
  !      print*, "check eigenfrequencies", wev1/PI/2.0D0*1.D3
  !   endif
  !  
  !   !! debug JS 07/22/2019
  !   !eigcut = (2.0D0*PI*1.D-2)**2*1.D-6 
  !   !nevout = 0
  !   !do i = 1,nevout0
  !   !   if (wev(i).gt.eigcut) then 
  !   !      nevout = nevout + 1
  !   !      if (nevout.eq.1) eigs = i - 1
  !   !   endif
  !   !enddo
  !   !!if (mymatvec%rank.eq.0) print*,nevout,nevout0
  !   nevout = nevout0; eigs = 0
  !
 
  !   allocate(Co0(nevout,nevout), Co1(nevout,nevout)); Co0 = 0.0D0
  !   allocate(Co2(nevout,nevout)); Co2 = 0.0D0
  !   do i = 1,nevout
  !      !Co1(i,:) = Co(i+eigs,eigs+1:eigs+nevout)*wev1(i+eigs)
  !      Co1(i,:) = Co(eigs+1:eigs+nevout,i+eigs)*wev1(i+eigs)
  !   enddo 
  !   !if(mymatvec%rank.eq.0) print*, nevout,Co1(nevout,:)
  !   allocate(Zo(2*nevout,2*nevout)); Zo = (0.0D0,0.0D0)
  !   Zo(1:nevout,nevout+1:2*nevout) = cmplx(Co1,Co0,16) 
  !   Zo(nevout+1:2*nevout,1:nevout) = cmplx(transpose(Co1),Co0,16) 
  !       
  !   do i = 1,nevout
  !      inputv1 = lanvecs((i+eigs-1)*mymatvec%pbsiz+1:(i+eigs)*mymatvec%pbsiz)
  !      call sparseRTV(inputv1,outputrv,mymatvec)
  !      do j = 1,i
  !         inputv2 = lanvecs((j+eigs-1)*mymatvec%pbsiz+1:(j+eigs)*mymatvec%pbsiz)
  !         st0 = sum(inputv2*outputrv); rv0 = 0.0D0;
  !         call mpi_allreduce(st0,rv0,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
  !         Co2(j,i) = rv0
  !         if (j.ne.i) then
  !            Co2(i,j) = - Co2(j,i)
  !         endif
  !      enddo
  !   enddo
  !   Zo(nevout+1:2*nevout,nevout+1:2*nevout) = cmplx(Co0,2.0D0*Co2,16)
  !   
  !   if (mymatvec%rank.eq.0) then 
  !      print*, "complex matrix (of size ", 2*nevout, ") is done"
  !   endif
  !   call report_time(mymatvec%comm)
  !  
  !   ! solve the reduced system
  !   allocate(wev2(2*nevout),wev3(nevout)) 
  !   call intface_zheev(2*nevout,Zo,wev2)
  !   !print*, wev2(nevout+1:2*nevout)
  !   !if (mymatvec%rank.eq.0) print*,Zo(nevout+1:2*nevout,1) 
  !   NEV = 0; 
  !   do i = 1,nevout
  !      wev3 = wev2(nevout+1:2*nevout)**2 
  !      if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
  !         NEV = NEV + 1
  !      endif 
  !   enddo
  !   allocate(eigval(nev),eigerr(nev))
  !   allocate(eigr(nev*mymatvec%pbsiz),eigi(nev*mymatvec%pbsiz))
  !   allocate(myeigid(nev)); myeigid = (/(i,i=1,nev)/)  
  !   j = 0
  !   do i = 1,nevout 
  !      if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
  !         j = j + 1
  !         eigval(j) = dsqrt(wev3(i))
  !      endif
  !   enddo

  !   allocate(vecs0(mymatvec%pbsiz,nevout))
  !   do i = 1,nevout
  !      vecs0(:,i) = lanvecs((i+eigs-1)*mymatvec%pbsiz+1:(i+eigs)*mymatvec%pbsiz)
  !   enddo
  !   wev2 = max(1.0D-30,wev2)  ! make sure it is positive 
  !   allocate(pv(nevout),pw(nevout))
  !   j = 0
  !   do i = 1,nevout 
  !      if (wev3(i).ge.XINTV(1).and.wev3(i).le.XINTV(2)) then 
  !         j = j + 1
  !         pv =  real(Zo(nevout+1:2*nevout,i+nevout))/wev2(i+nevout)
  !         pw = aimag(Zo(nevout+1:2*nevout,i+nevout))/wev2(i+nevout)
  !         ! JS 07/01/19
  !         ! check norm 
  !         npwr = 0.0D0; npwi = 0.0D0; npw = 0.0D0
  !         inputv1 = matmul(vecs0, pv)
  !         inputv2 = matmul(vecs0, pw)           
  !         call sparseBV(inputv1,outputbv,mymatvec)
  !         npwr = npwr+sum(inputv1*outputbv)
  !         npwi = npwi+sum(inputv2*outputbv) 
  !         call sparseBV(inputv2,outputbv,mymatvec)
  !         npwr = npwr+sum(inputv2*outputbv) 
  !         npwi = npwi-sum(inputv1*outputbv) 
  !         call sparseRTV(inputv1,outputrv,mymatvec)
  !         npwr = npwr-eigval(j)*sum(inputv2*outputrv) 
  !         npwi = npwi-eigval(j)*sum(inputv1*outputrv) 
  !         call sparseRTV(inputv2,outputrv,mymatvec)
  !         npwr = npwr-eigval(j)*sum(inputv1*outputrv) 
  !         npwi = npwi-eigval(j)*sum(inputv2*outputrv) 
  !         call mpi_allreduce(npwi,npw,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
  !         npwi = npw; npw = 0.0D0;
  !         call mpi_allreduce(npwr,npw,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
  !         npwr = npw
  !         !if (mymatvec%rank.eq.0) print*, npwi/npwr           
  !
  !         ! imagery 
  !         inputv1 = matmul(vecs0,pw)           
  !         call sparseBV(inputv1,outputbv,mymatvec)
  !         if (.not.mymatvec%fsexist) then 
  !            call sparseAV(inputv1,outputav,mymatvec)
  !         else
  !            call sparsefsAV(inputv1,outputav,mymatvec)
  !         endif
  !         inputv1 = matmul(vecs0,pv)
  !         call sparseRTV(inputv1,outputrv,mymatvec)
  !         inputv2 = outputav - eigval(j)**2*outputbv + 2.0D0*eigval(j)*outputrv
  !         npwi = 0.0D0
  !         call pnm_pvecnorm(inputv2,mymatvec%pbsiz,npwi)
  !         !if(mymatvec%rank.eq.0) print*,j,npwi**2/npw

  !         ! real 
  !         inputv1 = matmul(vecs0,pv)
  !         call sparseBV(inputv1,outputbv,mymatvec)
  !         if (.not.mymatvec%fsexist) then 
  !            call sparseAV(inputv1,outputav,mymatvec)
  !         else
  !            call sparsefsAV(inputv1,outputav,mymatvec)
  !         endif
  !         inputv1 = matmul(vecs0,pw)           
  !         call sparseRTV(inputv1,outputrv,mymatvec)
  !         inputv2 = outputav - 2.0D0*eigval(j)*outputrv - eigval(j)**2*outputbv
  !         npwr = 0.0D0
  !         call pnm_pvecnorm(inputv2,mymatvec%pbsiz,npwr)
  !         !if(mymatvec%rank.eq.0) print*,j,npwr**2/npw
  !         errd = (npwi**2 + npwr**2)/npw 
  !         !errd = npwr**2/npw 
  !         if (mymatvec%rank.eq.0) then
  !            eigerr(j) = dsqrt(errd/real(mymatvec%Gpbsiz,8))/dabs(eigval(j))
  !         endif
  !         eigr((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz) =mymatvec%B%diag*&
  !                        matmul(vecs0,pv)/dsqrt(npw)
  !         eigi((j-1)*mymatvec%pbsiz+1:j*mymatvec%pbsiz) =mymatvec%B%diag*&
  !                        matmul(vecs0,pw)/dsqrt(npw)
  !      endif
  !   enddo
  !   
  !   if (mymatvec%rank == 0) then
  !      print*, 'Relative errors:'
  !      do J = 1,NEV
  !         print*, 'Row', int(j,2), ' relative err. ', eigerr(j)
  !      enddo
  !      print*, '==============================================================='
  !      print*, 'Transform to frequencies (mHz), and periods (s)'
  !      do J = 1,NEV
  !         if (EIGVAL(J).gt.1.D-15) then 
  !            print*, 'Row', int(J,2), EIGVAL(J)/PI/2.0*1.D3,&
  !                     2.0D0*PI/EIGVAL(J)
  !         else
  !            print*, 'Row', int(J,2), EIGVAL(J), '...too small...'
  !         endif
  !      enddo
  !      print*, '================================================================'
  !   endif 
  !   
  !   ! save data    
  !   if (mymatvec%rank.eq.0) then 
  !      print*,"save eigenvalues"
  !      call pnm_save_freq(EIGVAL,NEV)
  !   endif  

  !   ! JS 06/29/19 
  !   if (pin%vecsave) then 
  !      if (mymatvec%rank.eq.0) print*,"save eigenvectors"
  !      lab = '_real_' 
  !      call pnm_save_data(lab,nev,myeigid,eigr) 
  !      lab = '_imag_' 
  !      call pnm_save_data(lab,nev,myeigid,eigi) 
  !   endif  
 
  !   DEALLOCATE(EIGVAL,myeigid)
  !   DEALLOCATE(lanvecs,eigr,eigi,vecs0)

  !   call pEVSL_FREEPOL_F90(POL)
  !   DEALLOCATE (SLI)
  !   call pEVSL_FINISH_F90(pevslAB)

  !end subroutine pnm_pevsl_lanrot



  subroutine pnm_save_freq(vals,lfqs)
     integer, intent(in)            :: lfqs
     real(kind=rkind), intent(in)   :: vals(lfqs)

     integer                        :: id0,i
     logical                        :: alive   
 
     id0 = 1961
     inquire(file=trim(pin%ffreq),exist=alive)
     if (alive) then
        open(id0, file=trim(pin%ffreq), status="old", action="write")
     else
        open(id0, file=trim(pin%ffreq), status="new", action="write")
     endif
     do i = 1,lfqs
        write(id0,*) vals(i)
     enddo
    
     close(id0)
  end subroutine pnm_save_freq


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
 
     call pnm_save_integer(bufsize,unstrM%new%vstat,disp,pin%fvstat)   

  end subroutine pnm_save_results

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

  subroutine pnm_save_data(lab,nev,ord,vecs)
     integer, intent(in)                          :: nev
     character(len=1024), intent(in)              :: lab
     integer, intent(in)                          :: ord(nev)
     real*8, intent(in)                           :: vecs(nev*mymatvec%pbsiz)
     ! save files
     integer                                      :: i,j,bfsiz
     integer(kind=mpi_offset_kind)                :: offset,leight
     real*8, allocatable                          :: vdat(:)
     character(len=1024)                          :: fname,fnb


     leight = 8
     
     do i = 1,nev
        j = ord(i)
        write(fnb,*) i; fnb = adjustl(fnb) 
        fname  = trim(pin%fvdata)//trim(adjustl(lab))//trim(fnb)//'.dat'   
        fname  = trim(fname) 
        offset = mymatvec%B%sizdist(mymatvec%rank+1)*leight 
        bfsiz  = mymatvec%pbsiz 
        vdat   = vecs((j-1)*bfsiz+1:j*bfsiz) 
        call pnm_save_dreal(bfsiz,vdat,offset,fname)
     enddo 

     call writelog("save eigenvectors",pin%comm)
     call report_time(pin%comm)

 
  end subroutine pnm_save_data


  subroutine pnm_pvecnorm(inputv,lz,nm2)
     integer, intent(in)     :: lz
     real*8, intent(in)      :: inputv(lz)
     real*8, intent(out)     :: nm2

     integer                 :: ierr
     real*8                  :: val0

     val0 = sum(inputv**2)
     call mpi_allreduce(val0,nm2,1,mpi_real8,mpi_sum,mymatvec%comm,ierr)
     nm2 = dsqrt(nm2)      

  end subroutine pnm_pvecnorm

end module pevsl_mod 
