module cg_matvec_mod

  use mpi
  use omp_lib
  use ISO_C_BINDING
  use para_mod,                 only: rkind,pin 
  use cg_datatype_mod,          only: mvparameters
  use geometry_mod,             only: unstrM,refs 
  use cg_create_matrix_mod,     only: CGM
  use cg_models_mod,            only: models
  use utility_mod            
  use exafmm_mod


  implicit none

  public                           :: setupmatvec
  public                           :: mymatvec
  
  type(mvparameters), save         :: mymatvec

contains
   
  subroutine setupmatvec()
    integer                        :: i,j,ierr,nthd,nthd0
    integer                        :: NFIRST,MLAN,LANSTEP,CHEBTYPE,CHEBDEG
    real*8                         :: LMIN,LMAX,TOL
    integer*8                      :: pevslap,pevslB
    integer                        :: LSPOLMAXDEG
    real*8                         :: LSPOLTOL
    logical                        :: ifOMP


!$OMP PARALLEL private(nthd)
  mymatvec%nthd = omp_get_num_threads()
!$OMP END PARALLEL
    if(unstrM%rank.eq.0) print*,'number of threads', mymatvec%nthd
    
    mymatvec%rank              = unstrM%rank
    mymatvec%comm              = unstrM%comm
    mymatvec%nproc             = unstrM%nproc
    mymatvec%fsexist           = unstrM%fsexist
    mymatvec%purefluid         = unstrM%purefluid
    mymatvec%sG                = pin%phi1

    ! B matrix information
    mymatvec%B%siz             = CGM%B%siz  
    mymatvec%B%Gsiz            = CGM%B%Gsiz
    mymatvec%B%NNZ             = CGM%B%NNZ
    mymatvec%B%GNNZ            = CGM%B%GNNZ

    allocate(mymatvec%B%sizdist(mymatvec%nproc+1))
    allocate(mymatvec%B%rowdist(mymatvec%B%siz+1))
    allocate(mymatvec%B%col(mymatvec%B%NNZ))
    allocate(mymatvec%B%val(mymatvec%B%NNZ))

    mymatvec%B%sizdist         = CGM%B%sizdist
    mymatvec%B%rowdist         = CGM%B%rowdist
    mymatvec%B%col             = CGM%B%col-1
    !mymatvec%B%val             = CGM%val

    mymatvec%Gpbsiz            = CGM%B%Gsiz
    mymatvec%pbsiz             = CGM%B%siz

    call Bdiagscaling()
   
    !print*, maxval(mymatvec%B%diag),minval(mymatvec%B%diag),mymatvec%rank
 
    call PEVSL_PARCSRCREATE_F90(  mymatvec%Gpbsiz,mymatvec%Gpbsiz,&
         mymatvec%B%sizdist,mymatvec%B%sizdist,mymatvec%B%rowdist,&
         mymatvec%B%col,mymatvec%B%val,mymatvec%comm,mymatvec%sBV) 

    if (mymatvec%rank.eq.0) print*, 'set up B matvec' 
    ! set up B solve function
    call pEVSL_Start_F90(mymatvec%comm,pevslB)
    ! SET PROB SIZE: NFIRST IS UNDEFINED
    NFIRST = -1
    call PEVSL_SETPROBSIZES_F90(pevslB,mymatvec%Gpbsiz, mymatvec%pbsiz, NFIRST)
 
    call PEVSL_SETAMV_F90(pevslB,sparseBV, mymatvec)
    ! NEED TO HAVE EIG BOUNDS OF B IN SETUP [DONE ABOVE BUT WILL DO AGAIN IN
    ! SETUP]
    MLAN = 2000; LANSTEP = 3000
    TOL = 1.0D-16; CHEBTYPE = 2
    call pEVSL_LANBOUNDS_F90(pevslB, MLAN, LANSTEP,TOL, LMIN, LMAX)
    if (mymatvec%rank.eq.0) print*, 'B cond. number', LMAX/LMIN,LMIN,LMAX

    if (pin%s%porder.eq.2) then 
       CHEBDEG = 60 
    else 
       CHEBDEG = 30 !25 !30
    endif
    call pEVSL_SETUP_CHEBITER_F90(LMIN,LMAX,CHEBDEG,mymatvec%sBV,mymatvec%chebB)
    
    call PEVSL_FINISH_F90(pevslB)

    if (pin%phi1) then
       !call pnm_check_exafmm(mymatvec%nthd,mymatvec%comm)
       call pnm_setup_exafmm() 
    endif

    ! added JS 052319 for RT mat-vec
    if (pin%rotation) then 
       ! RT matrix information
       mymatvec%RT%siz            = CGM%RT%siz  
       mymatvec%RT%Gsiz           = CGM%RT%Gsiz
       mymatvec%RT%NNZ            = CGM%RT%NNZ
       mymatvec%RT%GNNZ           = CGM%RT%GNNZ

       allocate(mymatvec%RT%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%RT%rowdist(mymatvec%RT%siz+1))
       allocate(mymatvec%RT%col(mymatvec%RT%NNZ))
       allocate(mymatvec%RT%val(mymatvec%RT%NNZ))

       mymatvec%RT%sizdist        = CGM%RT%sizdist
       mymatvec%RT%rowdist        = CGM%RT%rowdist
       mymatvec%RT%col            = CGM%RT%col-1
       mymatvec%RT%val            = CGM%RT%val
       ! set up parallel mat-vec
       call PEVSL_PARCSRCREATE_F90(     mymatvec%Gpbsiz,mymatvec%Gpbsiz,&
            mymatvec%RT%sizdist,mymatvec%RT%sizdist,mymatvec%RT%rowdist,&
            mymatvec%RT%col,mymatvec%RT%val,mymatvec%comm,mymatvec%sRTV) 
    endif

    if (unstrM%fsexist.or.unstrM%purefluid) then 
       ! Ad matrix information
       mymatvec%Ad%siz            = CGM%Ad%siz  
       mymatvec%Ad%Gsiz           = CGM%Ad%Gsiz
       mymatvec%Ad%NNZ            = CGM%Ad%NNZ
       mymatvec%Ad%GNNZ           = CGM%Ad%GNNZ

       allocate(mymatvec%Ad%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%Ad%rowdist(mymatvec%Ad%siz+1))
       allocate(mymatvec%Ad%col(mymatvec%Ad%NNZ))
       allocate(mymatvec%Ad%val(mymatvec%Ad%NNZ))

       mymatvec%Ad%sizdist        = CGM%Ad%sizdist
       mymatvec%Ad%rowdist        = CGM%Ad%rowdist
       mymatvec%Ad%col            = CGM%Ad%col-1
       mymatvec%Ad%val            = CGM%Ad%val

       !print*,maxval(CGM%Ad%val),minval(CGM%Ad%val) 
 
       call PEVSL_PARCSRCREATE_F90(     mymatvec%Gpbsiz,mymatvec%Gpbsiz,&
            mymatvec%Ad%sizdist,mymatvec%Ad%sizdist,mymatvec%Ad%rowdist,&
            mymatvec%Ad%col,mymatvec%Ad%val,mymatvec%comm,mymatvec%sAdV) 

       if (mymatvec%rank.eq.0) print*, 'set up Ad matvec'
       
       ! Ap matrix information
       mymatvec%Ap%siz            = CGM%Ap%siz  
       mymatvec%Ap%Gsiz           = CGM%Ap%Gsiz
       mymatvec%Ap%NNZ            = CGM%Ap%NNZ
       mymatvec%Ap%GNNZ           = CGM%Ap%GNNZ

       allocate(mymatvec%Ap%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%Ap%rowdist(mymatvec%Ap%siz+1))
       allocate(mymatvec%Ap%col(mymatvec%Ap%NNZ))
       allocate(mymatvec%Ap%val(mymatvec%Ap%NNZ))

       mymatvec%Ap%sizdist        = CGM%Ap%sizdist
       mymatvec%Ap%rowdist        = CGM%Ap%rowdist
       mymatvec%Ap%col            = CGM%Ap%col-1
       mymatvec%Ap%val            = - CGM%Ap%val 
      
       !if(mymatvec%rank.eq.0) print*,CGM%Ap%sizdist
       !if(mymatvec%rank.eq.0) print*,CGM%Ap%rowdist
 
       !print*, minval(-CGM%Ap%val),unstrM%rank 
       call Apdiagscaling()
       !print*, maxval(mymatvec%Ap%diag),minval(mymatvec%Ap%diag),mymatvec%rank
     
       call PEVSL_PARCSRCREATE_F90(   mymatvec%Ap%Gsiz,mymatvec%Ap%Gsiz,&
            mymatvec%Ap%sizdist,mymatvec%Ap%sizdist,mymatvec%Ap%rowdist,&
            mymatvec%Ap%col,mymatvec%Ap%val,mymatvec%comm,mymatvec%sApV) 

       if (mymatvec%rank.eq.0) print*, 'set up Ap matvec'
       ! set up Ap solve function
       call pEVSL_Start_F90(mymatvec%comm,pevslAp)
       ! SET PROB SIZE: NFIRST IS UNDEFINED
       NFIRST = -1
       call PEVSL_SETPROBSIZES_F90(pevslAp,mymatvec%Ap%Gsiz, mymatvec%Ap%siz, NFIRST)
 
       call PEVSL_SETAMV_F90(pevslAp, sparseApV, mymatvec)
       ! NEED TO HAVE EIG BOUNDS OF B IN SETUP [DONE ABOVE BUT WILL DO AGAIN IN
       ! SETUP]
       MLAN = 2000; LANSTEP = 3000
       TOL = 1.0D-16; CHEBTYPE = 2
       call pEVSL_LANBOUNDS_F90(pevslAp, MLAN, LANSTEP,TOL, LMIN, LMAX)

       call mpi_barrier(mymatvec%comm,ierr)
       if (mymatvec%rank.eq.0) print*, 'Ap cond. number ', LMAX/LMIN, LMIN, LMAX

       if (pin%s%porder.eq.2) then 
          CHEBDEG = 60 
       else 
          CHEBDEG = 30 !25 !30
       endif
       !call pEVSL_SETUP_CHEBITER_F90(LMIN,LMAX,CHEBDEG,&
       !              mymatvec%sApV,MYMATVEC%COMM,mymatvec%chebAp)
       call pEVSL_SETUP_CHEBITER_F90(LMIN,LMAX,CHEBDEG, mymatvec%sApV,mymatvec%chebAp)
       
       ! E matrix information
       mymatvec%E%siz             = CGM%E%siz  
       mymatvec%E%Gsiz            = CGM%E%Gsiz
       mymatvec%E%NNZ             = CGM%E%NNZ
       mymatvec%E%GNNZ            = CGM%E%GNNZ

       allocate(mymatvec%E%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%E%rowdist(mymatvec%E%siz+1))
       allocate(mymatvec%E%col(mymatvec%E%NNZ))
       allocate(mymatvec%E%val(mymatvec%E%NNZ))

       mymatvec%E%sizdist         = CGM%E%sizdist
       mymatvec%E%rowdist         = CGM%E%rowdist
       mymatvec%E%col             = CGM%E%col-1
       mymatvec%E%val             = CGM%E%val
       
       !print*, minval(CGM%ET%val),minval(CGM%E%val),unstrM%rank
       !print*, maxval(CGM%ET%val),maxval(CGM%E%val),unstrM%rank


       call PEVSL_PARCSRCREATE_F90(  mymatvec%Ad%Gsiz,mymatvec%Ap%Gsiz,&
            mymatvec%Ad%sizdist,mymatvec%Ap%sizdist,mymatvec%E%rowdist,&
            mymatvec%E%col,mymatvec%E%val,mymatvec%comm,mymatvec%sEV) 

       if (mymatvec%rank.eq.0) print*, 'set up E matvec'

       ! ET matrix information
       mymatvec%ET%siz            = CGM%ET%siz  
       mymatvec%ET%Gsiz           = CGM%ET%Gsiz
       mymatvec%ET%NNZ            = CGM%ET%NNZ
       mymatvec%ET%GNNZ           = CGM%ET%GNNZ

       allocate(mymatvec%ET%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%ET%rowdist(mymatvec%ET%siz+1))
       allocate(mymatvec%ET%col(mymatvec%ET%NNZ))
       allocate(mymatvec%ET%val(mymatvec%ET%NNZ))

       mymatvec%ET%sizdist        = CGM%ET%sizdist
       mymatvec%ET%rowdist        = CGM%ET%rowdist
       mymatvec%ET%col            = CGM%ET%col-1
       mymatvec%ET%val            = CGM%ET%val
       
       call PEVSL_PARCSRCREATE_F90(   mymatvec%Ap%Gsiz,mymatvec%Ad%Gsiz,&
            mymatvec%Ap%sizdist,mymatvec%Ad%sizdist,mymatvec%ET%rowdist,&
            mymatvec%ET%col,mymatvec%ET%val,mymatvec%comm,mymatvec%sETV) 

       if (mymatvec%rank.eq.0) print*, 'set up ET matvec'

 
    else
       ! A matrix information
       mymatvec%A%siz             = CGM%A%siz  
       mymatvec%A%Gsiz            = CGM%A%Gsiz
       mymatvec%A%NNZ             = CGM%A%NNZ
       mymatvec%A%GNNZ            = CGM%A%GNNZ

       allocate(mymatvec%A%sizdist(mymatvec%nproc+1))
       allocate(mymatvec%A%rowdist(mymatvec%A%siz+1))
       allocate(mymatvec%A%col(mymatvec%A%NNZ))
       allocate(mymatvec%A%val(mymatvec%A%NNZ))

       mymatvec%A%sizdist         = CGM%A%sizdist
       mymatvec%A%rowdist         = CGM%A%rowdist
       mymatvec%A%col             = CGM%A%col-1
       mymatvec%A%val             = CGM%A%val
       
       call PEVSL_PARCSRCREATE_F90(  mymatvec%Gpbsiz,mymatvec%Gpbsiz,&
            mymatvec%A%sizdist,mymatvec%A%sizdist,mymatvec%A%rowdist,&
            mymatvec%A%col,mymatvec%A%val,mymatvec%comm,mymatvec%sAV) 

       if (mymatvec%rank.eq.0) print*, 'set up A matvec' 
    endif 


  end subroutine setupmatvec

  subroutine Bdiagscaling()
    integer                                     :: i,j,k,l,m,ierr
    integer                                     :: ll,lsiz,lsiz0,lsiz1
    integer, allocatable, dimension(:)          :: vtmp0,vtmp1
    integer, allocatable, dimension(:)          :: cnt0,cnt1
    integer, allocatable, dimension(:)          :: itmp,itmp0,itmp1,idtmp
    real(kind=rkind), allocatable, dimension(:) :: dtmp0,dtmp1

    allocate(mymatvec%B%diag(mymatvec%B%siz))
    lsiz = mymatvec%B%sizdist(mymatvec%rank+1) 
    j = lsiz + 1
    !do i = 1,mymatvec%B%NNZ
    !   if (CGM%B%col(i).eq.j) then
    !      mymatvec%B%diag(j-lsiz) = 1.0D0/dsqrt(CGM%B%val(i))
    !      j = j + 1
    !   endif
    !enddo

    do i = 1,mymatvec%B%siz
       ll = mymatvec%B%rowdist(i+1) - mymatvec%B%rowdist(i)
       do l = 1,ll
          k = l + mymatvec%B%rowdist(i)
          if (CGM%B%col(k).eq.j) then
             ! JS add 02132021
             if (CGM%B%val(k).gt.0.0D0) then 
                mymatvec%B%diag(j-lsiz) = 1.0D0/dsqrt(CGM%B%val(k))
             else 
                print*, "Error: B(j,j) with j=",j,"is less than/equal to 0", &
                        CGM%B%val(k), " at rank", mymatvec%rank 
                stop
             endif
          endif
       enddo
       j = j + 1
    enddo

    allocate(itmp(mymatvec%B%NNZ)) 
    itmp = CGM%B%col; ll = int(mymatvec%B%NNZ,4)
    call simplessort(itmp,idtmp)
    deallocate(idtmp)
    call simplepickunique(itmp,ll,itmp0,lsiz0)
    ! itmp1: colid list
    
    allocate(vtmp0(mymatvec%nproc),vtmp1(mymatvec%nproc))
    vtmp0 = 0
    i = 1; j = 1
    do while (i.le.lsiz0)
       if (itmp0(i).gt.mymatvec%B%sizdist(j).and.&
           itmp0(i).le.mymatvec%B%sizdist(j+1)) then
           vtmp0(j) = vtmp0(j) + 1
           i = i + 1
       else
          if (j.lt.mymatvec%nproc) then
             j = j + 1
          endif
       endif
    enddo   
    
    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,mymatvec%comm,ierr)

    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + vtmp0(i)
       cnt1(i) = cnt1(i-1) + vtmp1(i)
    enddo
    lsiz1 = sum(vtmp1)
    allocate(itmp1(lsiz1)) 
     
    call MPI_ALLTOALLV(itmp0,vtmp0,cnt0,mpi_integer,&
                       itmp1,vtmp1,cnt1,mpi_integer,mymatvec%comm,ierr)

    allocate(dtmp1(lsiz1),dtmp0(lsiz0)) 
    do i = 1,lsiz1
       j = itmp1(i) - mymatvec%B%sizdist(mymatvec%rank+1)
       dtmp1(i) = mymatvec%B%diag(j)
    enddo

    call MPI_ALLTOALLV(dtmp1,vtmp1,cnt1,mpi_real8,&
                       dtmp0,vtmp0,cnt0,mpi_real8,mymatvec%comm,ierr)

    !print*,minval(dtmp0),maxval(dtmp0)

    lsiz = mymatvec%B%sizdist(mymatvec%rank+1) 
    do i = 1,mymatvec%B%siz
       ll = mymatvec%B%rowdist(i+1)-mymatvec%B%rowdist(i)
       do j = 1,ll
          l = mymatvec%B%rowdist(i)+j
          m = CGM%B%col(l)
          call findorder(m,itmp0,k)
          mymatvec%B%val(l) = CGM%B%val(l)*dtmp0(k)*mymatvec%B%diag(i)
       enddo
    enddo  
    !print*, minval(abs(mymatvec%B%val)), &
    !        maxval(abs(mymatvec%B%val)), unstrM%rank 
 

  end subroutine Bdiagscaling

  !----------------------------------------------------------------------
  subroutine Apdiagscaling()
    integer                                     :: i,j,k,l,m,ierr
    integer                                     :: ll,lsiz,lsiz0,lsiz1
    integer, allocatable, dimension(:)          :: vtmp0,vtmp1
    integer, allocatable, dimension(:)          :: cnt0,cnt1
    integer, allocatable, dimension(:)          :: itmp,itmp0,itmp1,idtmp
    real(kind=rkind), allocatable, dimension(:) :: dtmp0,dtmp1,dat0

    allocate(mymatvec%Ap%diag(mymatvec%Ap%siz))
    lsiz = mymatvec%Ap%sizdist(mymatvec%rank+1) 
    j = lsiz + 1
    do i = 1,mymatvec%Ap%siz
       ll = mymatvec%Ap%rowdist(i+1) - mymatvec%Ap%rowdist(i)
       do l = 1,ll
          k = l + mymatvec%Ap%rowdist(i)
          if (CGM%Ap%col(k).eq.j) then
             mymatvec%Ap%diag(j-lsiz) = 1.0D0/dsqrt(-CGM%Ap%val(k))
          endif
       enddo
       j = j + 1
    enddo

    allocate(itmp(mymatvec%Ap%NNZ)) 
    itmp = CGM%Ap%col; ll = int(mymatvec%Ap%NNZ,4)
    call simplessort(itmp,idtmp)
    deallocate(idtmp)
    call simplepickunique(itmp,ll,itmp0,lsiz0)
    ! itmp1: colid list
    
    allocate(vtmp0(mymatvec%nproc),vtmp1(mymatvec%nproc))
    vtmp0 = 0
    i = 1; j = 1
    do while (i.le.lsiz0)
       if (itmp0(i).gt.mymatvec%Ap%sizdist(j).and.&
           itmp0(i).le.mymatvec%Ap%sizdist(j+1)) then
           vtmp0(j) = vtmp0(j) + 1
           i = i + 1
       else
          if (j.lt.mymatvec%nproc) then
             j = j + 1
          endif
       endif
    enddo   
    
    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,mymatvec%comm,ierr)

    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + vtmp0(i)
       cnt1(i) = cnt1(i-1) + vtmp1(i)
    enddo
    lsiz1 = sum(vtmp1)
    allocate(itmp1(lsiz1)) 
     
    call MPI_ALLTOALLV(itmp0,vtmp0,cnt0,mpi_integer,&
                       itmp1,vtmp1,cnt1,mpi_integer,mymatvec%comm,ierr)

    allocate(dtmp1(lsiz1),dtmp0(lsiz0)) 
    do i = 1,lsiz1
       j = itmp1(i) - mymatvec%Ap%sizdist(mymatvec%rank+1)
       dtmp1(i) = mymatvec%Ap%diag(j)
    enddo

    call MPI_ALLTOALLV(dtmp1,vtmp1,cnt1,mpi_real8,&
                       dtmp0,vtmp0,cnt0,mpi_real8,mymatvec%comm,ierr)

    !print*,minval(dtmp0),maxval(dtmp0)

    lsiz = mymatvec%Ap%sizdist(mymatvec%rank+1) 
    do i = 1,mymatvec%Ap%siz
       ll = mymatvec%Ap%rowdist(i+1)-mymatvec%Ap%rowdist(i)
       do j = 1,ll
          l = mymatvec%Ap%rowdist(i)+j
          m = CGM%Ap%col(l)
          call findorder(m,itmp0,k)
          mymatvec%Ap%val(l) = - CGM%Ap%val(l)*dtmp0(k)*mymatvec%Ap%diag(i)
          !if (mymatvec%rank.eq.1.and.i.lt.2) then
          !if (mymatvec%rank.eq.0.and.i.lt.12) then
          !   !print*,mymatvec%Ap%val(l),-CGM%Ap%Val(l), m, i+lsiz
          !   print*, i+lsiz,m,CGM%Ap%val(l),mymatvec%Ap%val(l)
          !endif
       enddo
    enddo 
    !print*, minval(CGM%Ap%col),unstrM%rank 
    !print*,mymatvec%Ap%val
    !print*,minval(mymatvec%Ap%val), maxval(mymatvec%Ap%val),unstrM%rank 
    !allocate(dat0(mymatvec%Ap%siz)); dat0 = mymatvec%Ap%diag
    !allocate(idtmp(mymatvec%Ap%siz)); idtmp = (/(i,i=1,mymatvec%Ap%siz)/)
    !call ssort_real(dat0,idtmp); deallocate(idtmp)
    !do i = 1,mymatvec%Ap%siz
    !  print*, i, dat0(i)
    !enddo
    !print*,minval(mymatvec%Ap%diag), maxval(mymatvec%Ap%diag),unstrM%rank 

  end subroutine Apdiagscaling


  !----------------------------------------------------------------------
  subroutine sparseAV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%A%siz)    :: v,v0,w,w0,w2 
  
    w = 0.0D0; w2 = 0.0D0

    v0 = v*mymv%B%diag
    !print*, norm2(v)
    call pevsl_parcsrmatvec_f90(v0,w0,mymv%sAV)
    !print*, norm2(w)
    !todo
    w = w0*mymv%B%diag
    if (mymv%sG) then 
       call pnm_add_su(v0,w2,mymv)  
       !print*,'1',maxval(w2),maxval(w0)      
       w = w - w2*mymv%B%diag 
       !print*,'2',maxval(w)
    endif

  end subroutine sparseAV


  ! added JS 052019
  subroutine sparseAVnoP(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%A%siz)    :: v,v0,w,w0,w2 
  
    w = 0.0D0; w2 = 0.0D0

    v0 = v*mymv%B%diag
    !print*, norm2(v)
    call pevsl_parcsrmatvec_f90(v0,w0,mymv%sAV)
    !print*, norm2(w)
    !todo
    w = w0*mymv%B%diag

  end subroutine sparseAVnoP



  subroutine fmmAV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%dout%nv)  :: v,v0,w,w0,w2 
    
    v0 = v 
    call fmm_v2b(w0,v0,mymv%verbose) 
    w  = w0

  end subroutine fmmAV

  subroutine sparseBV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%B%siz)    :: v,w

    !v0 = v*mymv%Bdiag
    !print*, norm2(v)
    !call pevsl_parcsrmatvec_f90(v0,w0,mymv%sBV)
    !print*, norm2(w)
    !w = w0*mymv%Bdiag

    call pevsl_parcsrmatvec_f90(v,w,mymv%sBV)
  end subroutine sparseBV


  subroutine solveBV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%B%siz)    :: v,v0,w,w0 

    !v0 = v
    call pevsl_chebiter_f90(2,v,w,mymv%chebB)
    !w = w0

  end subroutine solveBV

  subroutine sparseApV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%Ap%siz)   :: v,w

    !v0 = v*mymv%Bdiag
    !print*, norm2(v)
    !call pevsl_parcsrmatvec_f90(v0,w0,mymv%sBV)
    !print*, norm2(w)
    !w = w0*mymv%Bdiag

    call pevsl_parcsrmatvec_f90(v,w,mymv%sApV)
  end subroutine sparseApV

  subroutine sparsefsAV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%Ad%siz)   :: v,v0,w,w0,w1,w2 
    real(kind=rkind), dimension(mymv%Ap%siz)   :: x0,x1,y0,y1

    w = 0.0D0; w0 = 0.0D0;  w1 = 0.0D0;  w2 = 0.0D0

    v0 = v*mymv%B%diag
    !print*, norm2(v)
    call pevsl_parcsrmatvec_f90(v0,w0,mymv%sAdV)

    call pevsl_parcsrmatvec_f90(v0,x0,mymv%sETV)

    x1 = x0*mymv%Ap%diag
    call pevsl_chebiter_f90(2,x1,y0,mymv%chebAp) 
    y1 = y0*mymv%Ap%diag

    call pevsl_parcsrmatvec_f90(y1,w1,mymv%sEV)
 
    if (mymv%sG) then 
       call pnm_add_su(v0,w2,mymv)        
    endif
    !!print*, norm2(w)
    w = w + (w0+w1)*mymv%B%diag
    !w = w0*mymv%Bdiag
    if (mymv%sG) then
       w = w - w2*mymv%B%diag 
    endif

  end subroutine sparsefsAV

  ! added JS 052019 
  subroutine sparsefsAVnoP(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%Ad%siz)   :: v,v0,w,w0,w1,w2 
    real(kind=rkind), dimension(mymv%Ap%siz)   :: x0,x1,y0,y1

    w = 0.0D0; w0 = 0.0D0;  w1 = 0.0D0;  w2 = 0.0D0

    v0 = v*mymv%B%diag
    !print*, norm2(v)
    call pevsl_parcsrmatvec_f90(v0,w0,mymv%sAdV)

    call pevsl_parcsrmatvec_f90(v0,x0,mymv%sETV)

    x1 = x0*mymv%Ap%diag
    call pevsl_chebiter_f90(2,x1,y0,mymv%chebAp) 
    y1 = y0*mymv%Ap%diag

    call pevsl_parcsrmatvec_f90(y1,w1,mymv%sEV)
 
    !!print*, norm2(w)
    w = w + (w0+w1)*mymv%B%diag

  end subroutine sparsefsAVnoP

  ! added JS 052319
  subroutine sparseRTV(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%RT%siz)   :: v,v0,w,w0,w2 
  
    w = 0.0D0; w2 = 0.0D0

    v0 = v*mymv%B%diag
    !print*, norm2(v)
    call pevsl_parcsrmatvec_f90(v0,w0,mymv%sRTV)
    !print*, norm2(w)
    !todo
    w = w0*mymv%B%diag

  end subroutine sparseRTV



  ! construct operation via fmm 
  subroutine pnm_add_su(v,w,mymv)
    type(mvparameters)                         :: mymv
    real(kind=rkind), dimension(mymv%B%siz)    :: v,w
    real(kind=rkind), dimension(mymv%B%siz)    :: v0,v1,w0,w1,w2

    integer(C_INT)                             :: ni,nj 
    integer                                    :: i,j,k,l
    integer                                    :: ii,jj,kk,ll 
    integer, dimension(pin%s%pNp)              :: e0,e1
    integer, dimension(pin%s%Nfp)              :: f0,f1,f2,f3 
    real(kind=rkind), dimension(pin%s%Nfp)     :: uf0,uf1
    real(kind=rkind), dimension(pin%s%pNp)     :: ual0,ual1

    real(kind=rkind)                           :: uc(pin%s%pNp,3),uall(pin%s%pNp)
    real(kind=rkind)                           :: Derv(pin%s%pNp,pin%s%pNp,3)
    
    real(kind=rkind), dimension(pin%s%Nfp,pin%s%Nfp)  :: ufall
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)  :: OP,OPt,OPtmp,rhod,davg
    
   
    real*8, allocatable, dimension(:)          :: vi,vb
    real*8, allocatable, dimension(:)          :: xb,yb,zb
    real*8, allocatable, dimension(:)          :: xv,yv,zv

    real*8                                     :: t0,t1,t2,t3

    davg =  1.D0/real(pin%s%pNp,8)
    do i = 1,pin%s%pNp 
       davg(i,i) = - 1.0D0 + davg(i,i)
    enddo

    call report_loc_time(mymv%comm,t0)

    w = 0.0D0; w0 = 0.0D0; w1 = 0.0D0

    !print*,'v', maxval(v),minval(abs(v)),mymv%rank
    mymv%vin%val = v 
    mymv%vcv%val = 0.0D0
    mymv%din%val = 0.0D0 

    call pnm_find_data(mymv%vin%lth,mymv%vin%ids,mymv%vin%val,&
      mymv%vcv%lth,mymv%vcv%pid,mymv%vcv%ids,mymv%vcv%val,mymv%comm)

    do i = 1,mymv%elm%ClNele
       e0 = mymv%elm%loc_ids(:,i)
       k = unstrM%lelist(i)
       if (models%coeff_loc(k*pin%s%pNp,models%p_vs).lt.pin%TOL) then
          e1 = mymv%vcv%nds(e0) + CGM%vnum(e0) - 3  
       else
          e1 = mymv%vcv%nds(e0) 
       endif

       uc = 0.0D0
       do j = 1,3
          !uc(:,j) = mymv%elm%loc_rho(:,i)*mymv%vcv%val(e1+j) 
          uc(:,j) = mymv%vcv%val(e1+j) 
       enddo
       uall = 0.0D0

       !Derv(:,:,1) = mymv%elm%loc_invJ(1,1,i)*refs%Drst(:,:,1) +&
       !              mymv%elm%loc_invJ(1,2,i)*refs%Drst(:,:,2) +&
       !              mymv%elm%loc_invJ(1,3,i)*refs%Drst(:,:,3)
       !
       !Derv(:,:,2) = mymv%elm%loc_invJ(2,1,i)*refs%Drst(:,:,1) +&
       !              mymv%elm%loc_invJ(2,2,i)*refs%Drst(:,:,2) +&
       !              mymv%elm%loc_invJ(2,3,i)*refs%Drst(:,:,3)
       !
       !Derv(:,:,3) = mymv%elm%loc_invJ(3,1,i)*refs%Drst(:,:,1) +&
       !              mymv%elm%loc_invJ(3,2,i)*refs%Drst(:,:,2) +&
       !              mymv%elm%loc_invJ(3,3,i)*refs%Drst(:,:,3)
     
       !! construct the operator
       !rhod = 0.0D0
       !do j = 1,pin%s%pNp
       !   rhod(j,j) = mymv%elm%loc_rho(j,i)
       !enddo
       !do j = 1,3
       !   OP  = matmul(refs%MassM,Derv(:,:,j))
       !   !OP  = matmul(OP,rhod)
       !   OPt = matmul(transpose(Derv(:,:,j)),refs%MassM) 
       !   !OPt = matmul(rhod,OPt)
       !   OP  = (OP+transpose(OPt))/2.0D0
       !   uall = uall + matmul(OP,mymv%elm%loc_rho(:,i)*uc(:,j))
       !   !uall = uall + matmul(OP,uc(:,j))
       !enddo 
       !!uall = uall*sum(mymv%elm%loc_rho(:,i))/real(pin%s%pNp,8)
       !uall = uall*mymv%elm%loc_dtJ(i)

       ! change JS 03092018
       uall = uall + mymv%elm%Dx(:,i)*mymv%elm%loc_rho(:,i)*uc(:,1)
       uall = uall + mymv%elm%Dy(:,i)*mymv%elm%loc_rho(:,i)*uc(:,2)
       uall = uall + mymv%elm%Dz(:,i)*mymv%elm%loc_rho(:,i)*uc(:,3)
       uall = uall*4.0D0/3.0D0*mymv%elm%loc_dtJ(i)

       !! change JS 02062018
       !ual0 = mymv%elm%loc_rho(:,i)*uc(:,1)
       !ual1 = matmul(davg,ual0)
       !uall = uall + mymv%elm%Dx(:,i)*ual1
       !ual0 = mymv%elm%loc_rho(:,i)*uc(:,2)
       !ual1 = matmul(davg,ual0)
       !uall = uall + mymv%elm%Dy(:,i)*ual1
       !ual0 = mymv%elm%loc_rho(:,i)*uc(:,3)
       !ual1 = matmul(davg,ual0)
       !uall = uall + mymv%elm%Dz(:,i)*ual1
       !
       !uall = uall*mymv%elm%loc_dtJ(i)*4.0D0/3.0D0

       !todo
       mymv%din%val(i) = sum(uall)
       !mymv%din%val(i) = 0.0D0
    enddo

    !print*,'1',maxval(mymv%din%val),minval(mymv%din%val),mymv%rank
    if (mymv%ncf.gt.0) then
       do kk = 1,mymv%ncf
          i = mymv%flist(kk)
       !do i = 1,mymv%srf%lsf
          l  = mymv%srf%inx(i)
          f0 = refs%Fmask(:,l)
          f1 = mymv%celm%loc_ids(f0,mymv%srf%s2e(1,i)) 

          uf0 = 0.0D0; uf1 = 0.0D0
          if (mymv%srf%sta(i).eq.0) then
             f2 = mymv%vcv%nds(f1)
             do j = 1,3
                uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                   mymv%srf%drt(j,i) 
                !uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
                !uf0 = mymv%vcv%val(f2+j)
             enddo
             uf0 = uf0*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))*mymv%srf%rho(i)

             !print*,maxval(uf0)
          elseif (mymv%srf%sta(i).eq.1) then 
             f2 = mymv%vcv%nds(f1) + CGM%vnum(f1) - 3
             do j = 1,3
                uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                   mymv%srf%drt(j,i) 
                !uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
             enddo
             uf0 = uf0*mymv%srf%rho(i)*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))
  
          else
            if (mymv%srf%sta(i).eq.2) then
               !print*,CGM%vnum(f1)
               ! first one is solid
               f2 = mymv%vcv%nds(f1) !+ 3
               do j = 1,3
                  uf0 = uf0 - matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                     mymv%srf%drt(j,i) 
                  !uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
               enddo
               uf0 = uf0*models%srho(1,i)*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))
               !uf0 = uf0

               f2 = mymv%vcv%nds(f1) + 3 
               do j = 1,3
                  uf1 = uf1 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                     mymv%srf%drt(j,i) 
                  !uf1 = uf1 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
               enddo
               uf1 = uf1*models%srho(2,i)*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))
               !uf1 = uf1 
               !print*,models%srho(2,i)-models%srho(1,i)
               uf0 = uf0 + uf1
               !uf0 = 0.0D0 

            elseif (mymv%srf%sta(i).eq.3) then
               ! first one is fluid
               f2 = mymv%vcv%nds(f1) !+ 3
               do j = 1,3
                  uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                     mymv%srf%drt(j,i) 
                  !uf0 = uf0 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
               enddo
               uf0 = uf0*models%srho(2,i)*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))
               !uf0 = uf0

               f2 = mymv%vcv%nds(f1) + 3
               do j = 1,3
                  uf1 = uf1 - matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))*&
                     mymv%srf%drt(j,i) 
                  !uf1 = uf1 + matmul(refs%MassF(:,:,l),mymv%vcv%val(f2+j))
               enddo
               uf1 = uf1*models%srho(1,i)*unstrM%loc_sJac(l,mymv%srf%s2e(1,i))
               !uf1 = uf1     

               uf0 = uf0 + uf1
               !uf0 = 0.0D0 

            endif

          endif 
          ! todo 
          mymv%din%val(kk+mymv%elm%ClNele) = sum(uf0)
          !mymv%din%val(kk+mymv%elm%ClNele) = 0.0D0
       enddo
    endif
    !print*,maxval(mymv%din%val),maxval(v),mymv%rank
    !mymv%din%val = 1.0D0
 
    allocate(vi(mymv%dout%nv),vb(mymv%din%nv))
    ! from din to dout
    vi = 0.0D0; vb = 0.0D0
    mymv%dout%val = 0.0D0
    !call report_loc_time(mymv%comm,t0)
 
    call pnm_pass_data(mymv%din%nv,mymv%din%val,mymv%dout%nv,&
            mymv%dout%pid,mymv%dout%ids,mymv%dout%val,mymv%comm)

    call report_loc_time(mymv%comm,t1)
 
    call fmm_v2b(vi,mymv%dout%val,mymv%verbose)
    !print*, vv   
    vi = vi + mymv%dout%lcR*mymv%dout%val 

    !vi  = mymv%dout%val
    call report_loc_time(mymv%comm,t2)
 
    !ni = mymv%dout%nv; nj = ni
    !allocate(xv(mymv%dout%nv)); xv = mymv%dout%cen(1,:)
    !allocate(yv(mymv%dout%nv)); yv = mymv%dout%cen(2,:)
    !allocate(zv(mymv%dout%nv)); zv = mymv%dout%cen(3,:)
    !allocate(xb(mymv%dout%nv)); xb = mymv%dout%cen(1,:)
    !allocate(yb(mymv%dout%nv)); yb = mymv%dout%cen(2,:)
    !allocate(zb(mymv%dout%nv)); zb = mymv%dout%cen(3,:)
    !call fmm_direct(ni,xv,yv,zv,vi,nj,xv,yv,zv,mymv%dout%val)
    !deallocate(xv,yv,zv,xb,yb,zb)
 
    !print*,"fmm",maxval(vi),maxval(mymv%dout%val),unstrM%rank
    !print*,"fmm_b2b",minval(vi),minval(mymv%dout%val),unstrM%rank
    !vi = mymv%dout%val 
    ! maps it back
    call pnm_pass_data(mymv%dout%nv,vi,mymv%din%nv,&
             mymv%din%pid,mymv%din%ids,vb,mymv%comm)

    !print*,"fmm",maxval(vb),minval(vb),mymv%rank
    !print*,"size",size(vb),mymv%elm%ClNele,mymv%rank
    !vb = vb*mymv%Gcon

    !do i = 1,mymv%elm%ClNele
    !   mymv%elm%val(i) = vb(i)  !*mymv%Gcon
    !enddo

    mymv%elm%val(:) = vb(1:mymv%elm%ClNele)

    !print*,maxval(mymv%elm%val),minval(mymv%elm%val),unstrM%rank

    mymv%celm%val = 0.0D0
    call pnm_find_data(mymv%elm%ClNele,mymv%elm%Clelist,mymv%elm%val,&
      mymv%celm%ClNele,unstrM%ClNids,mymv%celm%Clelist,mymv%celm%val,mymv%comm)

    !mymv%celm%val = mymv%celm%val*mymv%Gcon
    !print*,maxval(mymv%celm%val)

    do i = 1,mymv%celm%ClNele

       uall = 1.0D0*mymv%celm%val(i) 
       uc = 0.0D0
       !print*,maxval(uall)

       !Derv(:,:,1) = mymv%celm%loc_invJ(1,1,i)*refs%Drst(:,:,1) +&
       !              mymv%celm%loc_invJ(1,2,i)*refs%Drst(:,:,2) +&
       !              mymv%celm%loc_invJ(1,3,i)*refs%Drst(:,:,3)
       !
       !Derv(:,:,2) = mymv%celm%loc_invJ(2,1,i)*refs%Drst(:,:,1) +&
       !              mymv%celm%loc_invJ(2,2,i)*refs%Drst(:,:,2) +&
       !              mymv%celm%loc_invJ(2,3,i)*refs%Drst(:,:,3)
       !
       !Derv(:,:,3) = mymv%celm%loc_invJ(3,1,i)*refs%Drst(:,:,1) +&
       !              mymv%celm%loc_invJ(3,2,i)*refs%Drst(:,:,2) +&
       !              mymv%celm%loc_invJ(3,3,i)*refs%Drst(:,:,3)

       !rhod = 0.0D0
       !do j = 1,pin%s%pNp
       !   rhod(j,j) = mymv%celm%loc_rho(j,i)
       !enddo
       !do j = 1,3
       !   !uc(:,j) = matmul(refs%MassM,uall)
       !   !uc(:,j) = matmul(transpose(Derv(:,:,j)),uc(:,j))
       !   OP  = matmul(transpose(Derv(:,:,j)),refs%MassM)
       !   !OP  = matmul(rhod,OP)
       !   OPt = matmul(refs%MassM,Derv(:,:,j))
       !   !OPt = matmul(OPt,rhod) 
       !   OP  = (OP+transpose(OPt))/2.0D0
       !   uc(:,j) = matmul(OP,uall)
       !   uc(:,j) = mymv%celm%loc_rho(:,i)*uc(:,j)
       !   uc(:,j) = uc(:,j)*mymv%celm%loc_dtJ(i)
       !   !uc(:,j) = uc(:,j)*sum(mymv%celm%loc_rho(:,i))/real(pin%s%pNp,8)
       !enddo

       !do j = 1,3
       !   OP  = matmul(refs%MassM,Derv(:,:,j))
       !   OP  = matmul(OP,rhod)
       !   OPt = matmul(Derv(:,:,j),refs%MassM) 
       !   OPt = matmul(rhod,OPt)
       !   OPt  = (OPt+transpose(OP))/2.0D0
       !   uc(:,j) = matmul(OPt,uall)
       !enddo 
       
       ! change JS 03092018
       uc(:,1) = mymv%celm%loc_rho(:,i)*mymv%celm%Dx(:,i)*uall
       uc(:,2) = mymv%celm%loc_rho(:,i)*mymv%celm%Dy(:,i)*uall
       uc(:,3) = mymv%celm%loc_rho(:,i)*mymv%celm%Dz(:,i)*uall
       uc = uc*4.0D0/3.0*mymv%celm%loc_dtJ(i)
       
       !! change JS 02062018
       !ual0 = mymv%celm%Dx(:,i)*uall
       !ual1 = matmul(davg,ual0)
       !uc(:,1) = mymv%celm%loc_rho(:,i)*ual1 
       !ual0 = mymv%celm%Dy(:,i)*uall
       !ual1 = matmul(davg,ual0)
       !uc(:,2) = mymv%celm%loc_rho(:,i)*ual1 
       !ual0 = mymv%celm%Dz(:,i)*uall
       !ual1 = matmul(davg,ual0)
       !uc(:,3) = mymv%celm%loc_rho(:,i)*ual1 
       !uc = uc*mymv%celm%loc_dtJ(i)*4.0D0/3.0D0
       !print*,maxval(uc),unstrM%rank

       do j = 1,pin%s%pNp 
          k = mymv%celm%loc_lab(j,i)
          if (k.gt.0) then
             if (models%coeff_loc(i*pin%s%pNp,models%p_vs).lt.pin%TOL) then
                ll = mymv%vin%nds(k) + CGM%vnum(mymv%vin%lds(k)) - 3
             else
                ll = mymv%vin%nds(k)
             endif
             do l = 1,3
                w0(ll+l) = w0(ll+l) + uc(j,l)
             enddo 
          endif
       enddo
    enddo

    !print*,'1', maxval(w0),maxval(v),mymv%rank
    !todo
    w = w + w0*mymv%Gcon

    mymv%srf%val = 0.0D0
    if (mymv%ncf.gt.0) then
    do i = 1,mymv%ncf
       j = mymv%flist(i)
       mymv%srf%val(j) = vb(i+mymv%elm%ClNele)!*mymv%Gcon
    enddo
    endif
    !print*,vb
    !print*,minval(abs(mymv%srf%val)),mymv%rank

    mymv%esrf%val = 0.0D0

    call pnm_find_data(mymv%srf%lsf,mymv%srf%ids,mymv%srf%val,&
      mymv%esrf%lsf,mymv%esrf%pid,mymv%esrf%ids,mymv%esrf%val,mymv%comm)

    !print*,maxval(mymv%esrf%val),minval(mymv%esrf%val),mymv%rank
    !print*,mymv%esrf%val

    do i = 1,mymv%esrf%lsf
       uf0(:) = mymv%esrf%val(i)
       uf1 = 0.0D0; ufall = 0.0D0
       l  = mymv%esrf%inx(i)
       kk = mymv%esrf%s2e(1,i)

       if (l.lt.0) then
          l = - l
       endif
       f0 = refs%Fmask(:,l)
       do j = 1,3
          uf1 = matmul(refs%MassF(:,:,l),uf0)*mymv%esrf%rho(i) 
          ufall(:,j) = uf1*mymv%esrf%drt(j,i)*&
              unstrM%loc_sJac(l,kk)
          !uf1 = matmul(refs%MassF(:,:,l),uf0)
          !uf1 = uf0 
          !ufall(:,j) = uf1!*unstrM%loc_sJac(l,kk)!*mymv%esrf%rho(i)
       enddo

       ! node list on v
       f1 = mymv%celm%loc_lab(f0,kk)
 
       do j = 1,pin%s%Nfp
          if (f1(j).gt.0) then
             if (mymv%esrf%sta(i).eq.1.or.&
                 mymv%esrf%sta(i).eq.3) then
                ll = mymv%vin%nds(f1(j)) + &
                     CGM%vnum(mymv%vin%lds(f1(j))) - 3   
                !ll = mymv%vin%nds(f1(j))
             else
                ll = mymv%vin%nds(f1(j))
             endif
             !debug
             !if (mymv%esrf%sta(i).lt.2) then
             !print*,CGM%vnum(mymv%vin%lds(f1(j)))
             
             if (mymatvec%esrf%sta(i).ge.2.or.&
            dabs(mymatvec%esrf%rho(i)).gt.mymatvec%rhocut) then 
             do k = 1,3
                w1(ll+k) = w1(ll+k) + ufall(j,k) 
             enddo
             endif
          endif
       enddo
    enddo

    !print*,'2', maxval(w1),minval(abs(w1)),mymv%rank
    w = w + w1*mymv%Gcon

    call report_loc_time(mymv%comm,t3) 
    !print*,'3', maxval(w),minval(abs(w)),mymv%rank

    deallocate(vi,vb)
    
    !if(mymv%rank.eq.0) print*,t1-t0,t2-t1,t3-t2
 
  end subroutine pnm_add_su

  ! prepare for fmm
  !----------------------------------------------------------------------
  subroutine pnm_setup_exafmm()
     integer                             :: i,j,k,l,m,n,ierr
     integer(C_INT)                      :: nb,nv,ncrit,threads
     real(C_DOUBLE)                      :: eps0,kreal,kimag     
     character                           :: path
     integer                             :: ni,nj,ngbl,n0
     integer                             :: sf(3),lct,kk,ll
     real*8, allocatable, dimension(:)   :: vb,vv,vi
     real*8, allocatable, dimension(:)   :: xb,yb,zb
     real*8, allocatable, dimension(:)   :: xv,yv,zv

     ! diagonal correction
     real(kind=rkind)                    :: t0(3,3),t1(3,4),d0,d1
     real(kind=rkind), allocatable       :: davg(:,:),dtp(:),dre(:,:),dism(:,:)
     real(kind=rkind)                    :: ndism(3,3),invd(3,3)

     mymatvec%Gcon = 6.67408D-8!*1.D8
     mymatvec%rhocut = 1.0D-10 !0.1D0 

     allocate(dtp(pin%s%pNp))
     allocate(davg(pin%s%pNp,pin%s%pNp))
     allocate(dre(3,pin%s%pNp)) 
     allocate(dism(pin%s%pNp,3)) 


     davg = 1.D0/real(pin%s%pNp,8)
     do i = 1,pin%s%pNp 
        davg(i,i) = - 1.0D0 + davg(i,i)
     enddo


     ! transfer srf esrf lNele ClNele to mymatvec
     mymatvec%srf%lsf = unstrM%srf%lsf
     allocate(mymatvec%srf%s2e(2,mymatvec%srf%lsf))
     allocate(mymatvec%srf%pid(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%ids(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%inx(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%sta(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%rho(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%cen(3,mymatvec%srf%lsf))
     allocate(mymatvec%srf%drt(3,mymatvec%srf%lsf))
     allocate(mymatvec%srf%vtx(3,mymatvec%srf%lsf))
     allocate(mymatvec%srf%val(  mymatvec%srf%lsf))
     allocate(mymatvec%srf%lcR(  mymatvec%srf%lsf))
 
     mymatvec%srf%s2e = unstrM%srf%s2e
     do i = 1,mymatvec%srf%lsf
        ll = unstrM%srf%s2e(1,i)
        !print*, ll,unstrM%rank
        call findorder(ll,unstrM%Clelist,lct)
        mymatvec%srf%s2e(1,i) = lct
     enddo
     mymatvec%srf%pid = unstrM%srf%pid
     mymatvec%srf%ids = unstrM%srf%ids
     mymatvec%srf%inx = unstrM%srf%inx
     mymatvec%srf%cen = unstrM%srf%cen
     mymatvec%srf%drt = unstrM%srf%drt
     mymatvec%srf%rho = models%srho(2,:)-models%srho(1,:)

     !print*,mymatvec%srf%rho   
 
     do i = 1,mymatvec%srf%lsf
        !call findorder(mymatvec%srf%s2e(1,i),unstrM%Clelist,l)
        l = mymatvec%srf%s2e(1,i)
        !sf = refs%Fmask(:,mymatvec%srf%inx(i))
        sf = refs%FtoV(mymatvec%srf%inx(i),:)
        ll = 0
        do j = 1,3
           !k = unstrM%loc_t2v(sf(j),l)
           !call findorder(k,unstrM%Cvlist,lct)
           lct = unstrM%lt2vid((l-1)*pin%s%pNp+refs%vord(sf(j)))
           if (CGM%vnum(lct).eq.6) then  
              ll = ll + 1
           endif  
        enddo
        if (ll.eq.3) then
           !print*,CGM%vnum(unstrM%lt2vid((l-1)*pin%s%pNp+&
           !        refs%Fmask(:,mymatvec%srf%inx(i)))) 
           if (models%svs(1,i).lt.pin%TOL) then     
              mymatvec%srf%sta(i) = 3
           else
              mymatvec%srf%sta(i) = 2
           endif
           !print*,models%srho(2,i),mymatvec%srf%rho(i),mymatvec%srf%ids(i)
           !print*,mymatvec%srf%rho(i),mymatvec%srf%ids(i)
           !print*,mymatvec%srf%rho(i),mymatvec%srf%ids(i)
           !print*,i 
           !print*,models%svs(1,i)-models%svs(2,i),ll
           !mymatvec%srf%rho(i) = 0.0D0
        elseif (ll.gt.0.and.ll.lt.3) then
           if (models%svs(1,i).gt.pin%TOL) then
              mymatvec%srf%sta(i) = 0
           else
              mymatvec%srf%sta(i) = 1
           endif
           ! debug
           !mymatvec%srf%rho(i) = 0.0D0
        elseif (ll.eq.0) then
            mymatvec%srf%sta(i) = 0
            !print*,mymatvec%srf%rho(i),mymatvec%srf%inx(i),mymatvec%srf%ids(i)
            !mymatvec%srf%rho(i) = 0.0D0
        endif 
        !if (mymatvec%srf%ids(i).eq.53) then
        !   print*,models%svs(1,i),models%svs(2,i),ll
        !endif
        ! todo compute local contribution
        do j = 1,3
           t0(:,j) = unstrM%loc_vtx(:,(l-1)*4+sf(j))
        enddo

        call pnm_srf_lcR(l,mymatvec%srf%inx(i),t0,d0)
        mymatvec%srf%lcR(i) = d0
        !print*, d0,i 
     enddo 
 
     ! add info about node ids
     ! for esrf
     ll = 0
     do i = 1,unstrM%esrf%lsf
        !call findorder(mymatvec%srf%s2e(1,i),unstrM%Clelist,l)
        if (models%evs(1,i).gt.pin%TOL.and.&
            models%evs(2,i).lt.pin%TOL.and.&
           models%erho(2,i).gt.pin%TOL) then
           ll = ll + 1
        elseif (models%evs(1,i).lt.pin%TOL.and.&
               models%erho(1,i).gt.pin%TOL.and.&
                models%evs(2,i).gt.pin%TOL) then
           ll = ll + 1
        endif
     enddo

     ! separate FS
     mymatvec%esrf%lsf = unstrM%esrf%lsf + ll
     !print*,'1',unstrM%esrf%lsf,mymatvec%esrf%lsf,ll,mymatvec%rank

     allocate(mymatvec%esrf%s2e(2,mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%pid(  mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%ids(  mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%inx(  mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%sta(  mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%rho(  mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%drt(3,mymatvec%esrf%lsf))
     allocate(mymatvec%esrf%val(  mymatvec%esrf%lsf))

     !mymatvec%esrf%s2e = unstrM%esrf%s2e
     !mymatvec%esrf%pid = unstrM%esrf%pid
     !mymatvec%esrf%ids = unstrM%esrf%ids
     !mymatvec%esrf%drt = unstrM%esrf%drt
     ll = 1
     do i = 1,unstrM%esrf%lsf
        l = unstrM%esrf%inx(i)
        if (l.gt.0) then
           k = unstrM%esrf%s2e(1,i)
        else
           k = unstrM%esrf%s2e(2,i)
        endif
        call findorder(k,unstrM%Clelist,kk)
 
        if (models%evs(1,i).gt.pin%TOL.and.&
            models%evs(2,i).lt.pin%TOL.and.&
           models%erho(2,i).gt.pin%TOL) then
           mymatvec%esrf%s2e(1,ll:ll+1) = kk ! relative one
           mymatvec%esrf%pid(ll:ll+1)   = unstrM%esrf%pid(i)
           mymatvec%esrf%ids(ll:ll+1)   = unstrM%esrf%ids(i)
           mymatvec%esrf%inx(ll:ll+1)   = unstrM%esrf%inx(i)
           mymatvec%esrf%rho(ll)        = models%erho(1,i) ! solid
           mymatvec%esrf%rho(ll+1)      = models%erho(2,i) ! fluid
           mymatvec%esrf%sta(ll)        = 2 ! solid
           mymatvec%esrf%sta(ll+1)      = 3 ! fluid
           mymatvec%esrf%drt(:,ll)      = - unstrM%esrf%drt(:,i)
           mymatvec%esrf%drt(:,ll+1)    = unstrM%esrf%drt(:,i)
          
           !print*,'esrf',unstrM%esrf%ids(i) 
           !print*, models%evs(1,i)-models%evs(2,i) 
           ll = ll + 2
        elseif (models%evs(1,i).lt.pin%TOL.and.&
               models%erho(1,i).gt.pin%TOL.and.&
                models%evs(2,i).gt.pin%TOL) then
           mymatvec%esrf%s2e(1,ll:ll+1) = kk ! relative one
           mymatvec%esrf%pid(ll:ll+1)   = unstrM%esrf%pid(i)
           mymatvec%esrf%ids(ll:ll+1)   = unstrM%esrf%ids(i)
           mymatvec%esrf%inx(ll:ll+1)   = unstrM%esrf%inx(i)
           mymatvec%esrf%rho(ll)        = models%erho(2,i) ! solid
           mymatvec%esrf%rho(ll+1)      = models%erho(1,i) ! fluid
           mymatvec%esrf%sta(ll)        = 2 ! solid
           mymatvec%esrf%sta(ll+1)      = 3 ! fluid
           mymatvec%esrf%drt(:,ll)      = unstrM%esrf%drt(:,i)
           mymatvec%esrf%drt(:,ll+1)    = - unstrM%esrf%drt(:,i)

           !print*,'esrf',unstrM%esrf%ids(i),models%erho(2,i)-models%erho(1,i)
           !print*,'esrf',unstrM%esrf%ids(i),kk
           !print*, models%evs(1,i)-models%evs(2,i) 
           !print*, models%erho(1,i)-models%erho(2,i) 
           ll = ll + 2
        elseif ((models%evs(1,i).lt.pin%TOL.and.&
                 models%evs(2,i).lt.pin%TOL).or.&
                (models%evs(1,i).gt.pin%TOL.and.&
                 models%evs(2,i).gt.pin%TOL).or.&
                (models%evs(2,i).lt.pin%TOL.and.&
                models%erho(2,i).lt.pin%TOL)) then
           mymatvec%esrf%s2e(1,ll) = kk ! relative one
           mymatvec%esrf%pid(ll)   = unstrM%esrf%pid(i)
           mymatvec%esrf%ids(ll)   = unstrM%esrf%ids(i)
           mymatvec%esrf%inx(ll)   = unstrM%esrf%inx(i)
           mymatvec%esrf%rho(ll)   = models%erho(2,i) - models%erho(1,i)
           !mymatvec%esrf%rho(ll)   = 0.0D0
           mymatvec%esrf%drt(:,ll) = unstrM%esrf%drt(:,i)
           
           if (models%evs(1,i).lt.pin%TOL) then
              mymatvec%esrf%sta(ll) = 1
           else
              mymatvec%esrf%sta(ll) = 0
           endif

           ll = ll + 1
        endif
        !if (i.eq.53) print*,models%evs(:,i),models%erho(:,i),mymatvec%esrf%sta(ll-1)

     enddo
     !print*,'2', ll-1,mymatvec%esrf%lsf,mymatvec%rank
     !print*,minval(models%srho),mymatvec%rank
     !mymatvec%esrf%rho = 1.0D0

 
     ! for lNele
     mymatvec%elm%ClNele = unstrM%lNele
     allocate(mymatvec%elm%Clelist(mymatvec%elm%ClNele))
     allocate(mymatvec%elm%loc_dtJ(mymatvec%elm%ClNele))
     allocate(mymatvec%elm%loc_invJ(3,3,mymatvec%elm%ClNele))
     allocate(mymatvec%elm%loc_ids(pin%s%pNp,mymatvec%elm%ClNele))
     allocate(mymatvec%elm%loc_rho(pin%s%pNp,mymatvec%elm%ClNele))

     allocate(mymatvec%elm%val(mymatvec%elm%ClNele))
     allocate(mymatvec%elm%lcR(mymatvec%elm%ClNele))
    
     allocate(mymatvec%elm%Dx(pin%s%pNp,mymatvec%elm%ClNele))
     allocate(mymatvec%elm%Dy(pin%s%pNp,mymatvec%elm%ClNele))
     allocate(mymatvec%elm%Dz(pin%s%pNp,mymatvec%elm%ClNele))
    
 
     mymatvec%elm%Clelist  = unstrM%Clelist(unstrM%lelist)
     do i = 1,unstrM%lNele
        j = unstrM%lelist(i)
        mymatvec%elm%loc_dtJ(i) = unstrM%loc_detJ(j)
        mymatvec%elm%loc_invJ(:,:,i) = unstrM%loc_invJ(:,:,j)
        do k = 1,pin%s%pNp
           call findorder(unstrM%loc_t2v(k,j),unstrM%Cvlist,l)
           mymatvec%elm%loc_ids(k,i) = l
        enddo
        mymatvec%elm%loc_rho(:,i) = &
         models%coeff_loc((j-1)*pin%s%pNp+1:j*pin%s%pNp,models%p_rho)

        ! todo compute local contributions
        t1 = unstrM%loc_vtx(:,(j-1)*4+1:j*4)
        call pnm_elm_lcR(j,d1,t1)
        mymatvec%elm%lcR(i) = d1
        !print*,d1,i

        ! change JS 03092018
        ndism = 0.0D0; invd = 0.0D0
        dism = 0.0D0; dre = 0.0D0
        do k = 1,3
           !dtp       = unstrM%loc_vtx(k,(j-1)*4+1:j*4)
           dtp       = unstrM%loc_nods(k,(j-1)*pin%s%pNp+1:j*pin%s%pNp)
           dism(:,k) = sum(dtp)/real(pin%s%pNp,8) - dtp
        enddo
        ndism = matmul(transpose(dism),dism)
        call matinv(ndism,invd,3)
        
        dre = matmul(invd,transpose(dism))
        dre = matmul(dre,davg)

        mymatvec%elm%Dx(:,i) = dre(1,:)
        mymatvec%elm%Dy(:,i) = dre(2,:)
        mymatvec%elm%Dz(:,i) = dre(3,:)

        !print*,mymatvec%elm%Dx(:,i)
     enddo
      
     ! for ClNele
     mymatvec%celm%ClNele = unstrM%ClNele  
     allocate(mymatvec%celm%Clelist(mymatvec%celm%ClNele))
     allocate(mymatvec%celm%loc_dtJ(mymatvec%celm%ClNele))
     allocate(mymatvec%celm%loc_invJ(3,3,mymatvec%celm%ClNele))
     allocate(mymatvec%celm%loc_ids(pin%s%pNp,mymatvec%celm%ClNele))
     allocate(mymatvec%celm%loc_lab(pin%s%pNp,mymatvec%celm%ClNele))
     allocate(mymatvec%celm%loc_rho(pin%s%pNp,mymatvec%celm%ClNele))

     allocate(mymatvec%celm%val(mymatvec%celm%ClNele))

     allocate(mymatvec%celm%Dx(pin%s%pNp,mymatvec%celm%ClNele))
     allocate(mymatvec%celm%Dy(pin%s%pNp,mymatvec%celm%ClNele))
     allocate(mymatvec%celm%Dz(pin%s%pNp,mymatvec%celm%ClNele))

     mymatvec%celm%Clelist  = unstrM%Clelist
     mymatvec%celm%loc_dtJ  = unstrM%loc_detJ
     mymatvec%celm%loc_invJ = unstrM%loc_invJ
      
     do i = 1,unstrM%ClNele
        do k = 1,pin%s%pNp
           call findorder(unstrM%loc_t2v(k,i),unstrM%Cvlist,l)
           mymatvec%celm%loc_ids(k,i) = l
        enddo
        mymatvec%celm%loc_rho(:,i) = &
         models%coeff_loc((i-1)*pin%s%pNp+1:i*pin%s%pNp,models%p_rho)

        ! change JS 03092018
        ndism = 0.0D0; invd = 0.0D0
        dism = 0.0D0; dre = 0.0D0
        do k = 1,3
           !dtp       = unstrM%loc_vtx(k,(i-1)*4+1:i*4)
           dtp       = unstrM%loc_nods(k,(i-1)*pin%s%pNp+1:i*pin%s%pNp)
           dism(:,k) = sum(dtp)/real(pin%s%pNp,8) - dtp
        enddo
        !print*,dism(1,1),i
        dre = transpose(dism)
        ndism = matmul(dre,dism)
        call matinv(ndism,invd,3)
        !print*, ndism(1,1),i
        
        dre = matmul(invd,transpose(dism))
        dre = matmul(dre,davg)


        mymatvec%celm%Dx(:,i) = dre(1,:)
        mymatvec%celm%Dy(:,i) = dre(2,:)
        mymatvec%celm%Dz(:,i) = dre(3,:)

     enddo
     
     mymatvec%celm%loc_lab = 0
     do i = 1,unstrM%ClNele
        do j = 1,pin%s%pNp
           k = unstrM%loc_t2v(j,i)
           call findorder(k,unstrM%Cvlist,lct)
           if (unstrM%Cvpid(lct).eq.unstrM%rank) then
              call findorder(k,unstrM%new%vlist,l)
              mymatvec%celm%loc_lab(j,i) = l 
           endif
        enddo
     enddo
 
     ! prepare for exafmm 
     eps0 = 0.0D0; kreal = 0.0; kimag = 0.0 
     path = './'; ncrit = 1000; threads = mymatvec%nthd

     
     k = mod(unstrM%Ntet+unstrM%Nsrf,unstrM%nproc)
     if (k.eq.0) then
        n0 = (unstrM%Ntet+unstrM%Nsrf-k)/unstrM%nproc
     else
        n0 = (unstrM%Ntet+unstrM%Nsrf-k)/unstrM%nproc + 1
     endif

     l = unstrM%lNele + unstrM%srf%lsf
     ngbl = int(max(n0,unstrM%lNele+unstrM%srf%lsf)*1.1)

     allocate(vb(ngbl),vv(ngbl))
     allocate(xb(ngbl),yb(ngbl),zb(ngbl))
     
     allocate(xv(ngbl),yv(ngbl),zv(ngbl))

     xb(1:unstrM%lNele)   = unstrM%sloc(1,:)
     yb(1:unstrM%lNele)   = unstrM%sloc(2,:)
     zb(1:unstrM%lNele)   = unstrM%sloc(3,:)
     !xb(1+unstrM%lNele:l) = mymatvec%srf%cen(1,:)
     !yb(1+unstrM%lNele:l) = mymatvec%srf%cen(2,:)
     !zb(1+unstrM%lNele:l) = mymatvec%srf%cen(3,:)

     kk = 0
     do i = 1,mymatvec%srf%lsf
        if (mymatvec%srf%sta(i).ge.2.or.&
       dabs(mymatvec%srf%rho(i)).gt.mymatvec%rhocut) then 
           kk = kk + 1
           xb(kk+unstrM%lNele) = mymatvec%srf%cen(1,i)
           yb(kk+unstrM%lNele) = mymatvec%srf%cen(2,i)
           zb(kk+unstrM%lNele) = mymatvec%srf%cen(3,i)
        endif
     enddo

     mymatvec%ncf = kk     
     allocate(mymatvec%flist(kk))
     kk = 0
     do i = 1,mymatvec%srf%lsf
        if (mymatvec%srf%sta(i).ge.2.or.&
       dabs(mymatvec%srf%rho(i)).gt.mymatvec%rhocut) then 
           kk = kk + 1
           mymatvec%flist(kk) = i
           !print*,mymatvec%srf%rho(i)
        endif
     enddo
     print*,'reduced # of surfaces',&
      mymatvec%ncf,mymatvec%srf%lsf,mymatvec%rank 

     call mpi_reduce(kk,ll,1,mpi_integer,mpi_sum,0,mymatvec%comm,ierr)
     if (mymatvec%rank.eq.0) then
        print*, "total # of (jumped) surfaces", ll
     endif

 
     xv = xb; yv = yb; zv = zb

     mymatvec%din%nv = mymatvec%ncf + unstrM%lNele
 
     allocate(mymatvec%din%pid(  mymatvec%din%nv))
     allocate(mymatvec%din%ids(  mymatvec%din%nv))
     allocate(mymatvec%din%val(  mymatvec%din%nv))
     allocate(mymatvec%din%cen(3,mymatvec%din%nv))
  
     allocate(mymatvec%din%lcR(  mymatvec%din%nv))

     mymatvec%din%pid = -1
     mymatvec%din%ids = 0

     do i = 1,unstrM%lNele
        mymatvec%din%cen(:,i) = unstrM%sloc(:,i)
        mymatvec%din%lcR(i) = mymatvec%elm%lcR(i)
     enddo
     if (mymatvec%ncf.gt.0) then
        do i = 1,mymatvec%ncf
           j = mymatvec%flist(i)
           mymatvec%din%cen(:,i+unstrM%lNele) = mymatvec%srf%cen(:,j)
           mymatvec%din%lcR(i+unstrM%lNele) = mymatvec%srf%lcR(j)
        enddo
     endif
     !print*,xb(1)
     nv = mymatvec%din%nv; nb = nv
     !print*,'before',nv,nb,ngbl,unstrM%rank
     call fmm_init(eps0,kreal,kimag,ncrit,threads,path,&
                          nb,xb,yb,zb,vb,nv,xv,yv,zv,vv)
     if (unstrM%rank.eq.0) print*, "initialize exafmm"

     ni = mymatvec%din%nv; nj = ni
     !print*,"before", ni,mymatvec%rank
     call fmm_partition(ni,xb,yb,zb,vb,nj,xv,yv,zv,vv) 
     !print*,'after',ni,nj,ngbl,unstrM%rank
     !print*,"after", ni,mymatvec%rank
     if (unstrM%rank.eq.0) print*, "exafmm partition"

     mymatvec%dout%nv = ni 
     allocate(mymatvec%dout%pid(  mymatvec%dout%nv))
     allocate(mymatvec%dout%ids(  mymatvec%dout%nv))
     allocate(mymatvec%dout%val(  mymatvec%dout%nv))
     allocate(mymatvec%dout%cen(3,mymatvec%dout%nv))

     allocate(mymatvec%dout%lcR(  mymatvec%dout%nv))

     mymatvec%dout%cen(1,:) = xb(1:ni) 
     mymatvec%dout%cen(2,:) = yb(1:ni)
     mymatvec%dout%cen(3,:) = zb(1:ni)

     call fmm_buildTree()
     if (unstrM%rank.eq.0) print*, "exafmm build tree"

     ! check the spectra

     ! find out fmm parition
     call pnm_find_fmmpart()
     
     ! find dout%lcR
     call pnm_pass_data(mymatvec%din%nv,mymatvec%din%lcR,mymatvec%dout%nv,&
            mymatvec%dout%pid,mymatvec%dout%ids,mymatvec%dout%lcR,mymatvec%comm)

     ll = 0
     do i = 1,unstrM%cnvtx
        ll = ll + CGM%vnum(i)
     enddo
     mymatvec%vcv%lth = ll !unstrM%cnvtx*3 
     allocate(mymatvec%vcv%ids(mymatvec%vcv%lth))
     allocate(mymatvec%vcv%pid(mymatvec%vcv%lth))
     allocate(mymatvec%vcv%val(mymatvec%vcv%lth))
     allocate(mymatvec%vcv%nds(unstrM%cnvtx))
     ll = 0
     do i = 1,unstrM%cnvtx
        k = CGM%vstt(i) 
        do j = 1,CGM%vnum(i)
           mymatvec%vcv%ids(ll+j) = k + j
           mymatvec%vcv%pid(ll+j) = unstrM%Cvpid(i)
        enddo
        mymatvec%vcv%nds(i) = ll
        ll = ll + CGM%vnum(i)
     enddo

     ll = 0
     do i = 1,unstrM%new%nvtx
        call findorder(unstrM%new%vlist(i),unstrM%Cvlist,l)
        ll = ll + CGM%vnum(l) 
     enddo 

     mymatvec%vin%lth = ll !unstrM%new%nvtx*3  
     mymatvec%vin%stt = CGM%B%sizdist(unstrM%rank+1)
     allocate(mymatvec%vin%ids(mymatvec%vin%lth))
     allocate(mymatvec%vin%val(mymatvec%vin%lth))
     allocate(mymatvec%vin%nds(unstrM%new%nvtx))
     allocate(mymatvec%vin%lds(unstrM%new%nvtx))
     ll = 0
     do i = 1,unstrM%new%nvtx
        call findorder(unstrM%new%vlist(i),unstrM%Cvlist,l)
        mymatvec%vin%lds(i) = l
        k = CGM%vstt(l) 
        do j = 1,CGM%vnum(l)
           mymatvec%vin%ids(ll+j) = k + j
        enddo
        mymatvec%vin%nds(i) = ll
        ll = ll + CGM%vnum(l)
     enddo



  end subroutine pnm_setup_exafmm

  subroutine pnm_find_fmmpart()
     integer                                       :: i,j,k,l,ierr
     integer                                       :: ii,jj,kk,ll
     integer, allocatable, dimension(:)            :: it0,it1
     real(kind=rkind)                              :: tsm
     real(kind=rkind), allocatable, dimension(:)   :: dt0,dt1
     real(kind=rkind), allocatable, dimension(:,:) :: lcen


     do i = 0,mymatvec%nproc-1
        if(mymatvec%rank.eq.i) ll = mymatvec%dout%nv
        call mpi_bcast(ll,1,mpi_integer,i,mymatvec%comm,ierr) 
        allocate(dt0(ll),lcen(3,ll))     
        do j = 1,3
           if(mymatvec%rank.eq.i) dt0 = mymatvec%dout%cen(j,:)
           call mpi_bcast(dt0,ll,mpi_real8,i,mymatvec%comm,ierr)
           lcen(j,:) = dt0
        enddo

        allocate(it0(ll),it1(ll)); it0 = -1; it1 = 0
        do j = 1,ll
           do k = 1,mymatvec%din%nv
              tsm = dsqrt(sum((mymatvec%din%cen(:,k)-lcen(:,j))**2))
              !print*,tsm,k,j
              if (tsm.eq.0.D0) then
                 !print*,tsm,k,j
                 it0(j) = unstrM%rank
                 it1(j) = k
                 mymatvec%din%pid(k) = i 
                 mymatvec%din%ids(k) = j
              endif
           enddo
        enddo
        
        call mpi_reduce(it0,mymatvec%dout%pid,ll,mpi_integer,&
             mpi_sum,i,mymatvec%comm,ierr)
        call mpi_reduce(it1,mymatvec%dout%ids,ll,mpi_integer,&
             mpi_sum,i,mymatvec%comm,ierr)

        if (mymatvec%rank.eq.i) then    
           mymatvec%dout%pid = mymatvec%dout%pid + mymatvec%nproc - 1
           !print*,minval(mymatvec%dout%pid),minval(mymatvec%dout%ids),i
           if (minval(mymatvec%dout%pid).lt.0) then
              print*,"Error, dout pid", i
              stop
           endif
           if (minval(mymatvec%dout%ids).lt.0) then
              print*,"Error, dout ids", i
              stop
           endif
        endif
        deallocate(dt0,lcen,it0,it1)
     enddo

     !print*, minval(mymatvec%din%pid),minval(mymatvec%din%ids),unstrM%rank
     if (minval(mymatvec%din%pid).lt.0) then
        print*,"Error, din pid", i
        stop
     endif
     if (minval(mymatvec%din%ids).lt.0) then
        print*,"Error, din ids", i
        stop
     endif

  end subroutine pnm_find_fmmpart

  subroutine pnm_elm_lcR(inx,dout,t0)
     integer, intent(in)                     :: inx
     real(kind=rkind), intent(in)            :: t0(3,4)
     real(kind=rkind), intent(out)           :: dout

     integer                                 :: i,j,k
     real(kind=rkind)                        :: dtmp,dout0
     real(kind=rkind)                        :: r0,r1,r2,r3,r4
     real(kind=rkind)                        :: c0,c1,c2,c3,c4     
     real(kind=rkind)                        :: tm(3),tc(3)
 
     dtmp = sum(unstrM%loc_sJac(:,inx))*2.D0

     dout0 = dtmp/unstrM%loc_detJ(inx)/4.0D0!*1.5D0 
  
     do i = 1,3
        tc(i) = sum(t0(i,:))/4.0D0 
     enddo
     r0 = 0.D0
     ! vtx 1
     tm = (t0(:,2)+t0(:,3))/2.D0
     call tri_rest(tc,t0(:,1),tm,r1)  
     call theta_sol(tc,t0(:,2),t0(:,3),c1)
     r0 = r0 + dsqrt(sum((tc-tm)**2))*dacos(c1)*r1

     ! vtx 2
     tm = (t0(:,3)+t0(:,4))/2.D0
     call tri_rest(tc,t0(:,2),tm,r1)  
     call theta_sol(tc,t0(:,3),t0(:,4),c1)
     r0 = r0 + dsqrt(sum((tc-tm)**2))*dacos(c1)*r1

     ! vtx 3
     tm = (t0(:,4)+t0(:,1))/2.D0
     call tri_rest(tc,t0(:,3),tm,r1)  
     call theta_sol(tc,t0(:,4),t0(:,1),c1)
     r0 = r0 + dsqrt(sum((tc-tm)**2))*dacos(c1)*r1

     ! vtx 4
     tm = (t0(:,1)+t0(:,2))/2.D0
     call tri_rest(tc,t0(:,4),tm,r1)  
     call theta_sol(tc,t0(:,1),t0(:,2),c1)
     r0 = r0 + dsqrt(sum((tc-tm)**2))*dacos(c1)*r1
    
     r0 = r0/(unstrM%loc_detJ(inx)*8.0D0/3.0D0)     
     !print*, r0/dout0

     dout = r0*2.0D0 
     !dout = dout0
 
  end subroutine pnm_elm_lcR

  subroutine pnm_srf_lcR(inx,inxf,t0,dout)
     integer, intent(in)                     :: inx,inxf
     real(kind=rkind), intent(in)            :: t0(3,3)
     real(kind=rkind), intent(out)           :: dout
     
     real(kind=rkind)                        :: s0,s1(3),sc(3),ds(3)
     real(kind=rkind)                        :: tha0(3),tha1(3),rtmp,rtmp1
     real(kind=rkind)                        :: dout0,dout1
     real(kind=rkind)                        :: cout0,cout1
     integer                                 :: i,j,k,l
 
     s0 = unstrM%loc_sJac(inxf,inx)*2.D0
 
     s1(1) = dsqrt(sum((t0(:,2)-t0(:,3))**2))
     s1(2) = dsqrt(sum((t0(:,3)-t0(:,1))**2))
     s1(3) = dsqrt(sum((t0(:,1)-t0(:,2))**2))

     dout1 = 0.5D0*sum(s1)/s0!*1.54D0 !*2.0D0

     ! center
     do i = 1,3
        sc(i) = sum(t0(i,:))/3.0D0 
     enddo
 
     !do i = 1,3
     !   ds(i) = dsqrt(sum((t0(:,i)-sc)**2))
     !enddo

     !call theta_sol(t0(:,1),sc,t0(:,2),rtmp)
     !tha0(1) = rtmp
     !call theta_sol(t0(:,2),sc,t0(:,3),rtmp)
     !tha0(2) = rtmp
     !call theta_sol(t0(:,3),sc,t0(:,1),rtmp)
     !tha0(3) = rtmp

     !call theta_sol(t0(:,2),sc,t0(:,1),rtmp)
     !tha1(1) = rtmp
     !call theta_sol(t0(:,3),sc,t0(:,2),rtmp)
     !tha1(2) = rtmp
     !call theta_sol(t0(:,1),sc,t0(:,3),rtmp)
     !tha1(3) = rtmp

     !dout0 = 0.0D0
     !do i = 1,3
     !   !dout0 = dout0 + ds(i)*sin(tha0(i))*&
     !   !(dlog(1.D0/dtan(tha1(i)/2.D0))-dlog(dtan(tha0(i)/2.D0)))
     !   ! tan (theta/2) 
     !   rtmp  = dsqrt((1.D0-tha0(i))/(1.D0+tha0(i)))
     !   rtmp1 = dsqrt((1.D0-tha1(i))/(1.D0+tha1(i)))
     !   dout0 = dout0+ds(i)*dsqrt(1.D0-tha0(i)**2)*&
     !         (dlog(1.D0/rtmp1)-dlog(rtmp))
     !enddo 

     cout1 = 0.D0
     call tri_rest(sc,t0(:,1),t0(:,2),cout0)
     cout1 = cout1 + cout0
     call tri_rest(sc,t0(:,2),t0(:,3),cout0)
     cout1 = cout1 + cout0
     call tri_rest(sc,t0(:,3),t0(:,1),cout0)
     cout1 = cout1 + cout0
     !print*,dout0, cout1    
     dout0 = cout1 

     dout = dout0/s0
 
     if (dout/dout1.gt.2.0.or.&
         dout/dout1.le.1.0) then 
        print*,"Warning: surface 1/R",inxf,inx,mymatvec%rank
     endif

     !dout = dout1
     !print*, dout/dout1

  end subroutine pnm_srf_lcR

  subroutine theta_sol(v0,v1,v2,theta)
     real(kind=rkind), intent(in)            :: v0(3),v1(3),v2(3)
     real(kind=rkind), intent(out)           :: theta
     
     real(kind=rkind)                        :: d1(3),d2(3)
     real(kind=rkind)                        :: r1,r2,r0

     d1 = v1 - v0; d2 = v2 - v0 

     r1 = dsqrt(sum(d1**2))
     r2 = dsqrt(sum(d2**2))
 
     r0 = dsqrt(sum((d1*d2)**2))
 
     !theta = dacos(r0/(r1*r2)) 
     ! cos theta
     theta = r0/(r1*r2)

  end subroutine theta_sol


  subroutine tri_rest(vc,v1,v2,rout)
     real(kind=rkind), intent(in)            :: vc(3),v1(3),v2(3)
     real(kind=rkind), intent(out)           :: rout

     real(kind=rkind)                        :: c1,c2
     real(kind=rkind)                        :: r1,r2


     call theta_sol(v1,vc,v2,c1)

     call theta_sol(v2,vc,v1,c2)

     r1  = dsqrt((1.D0-c1)/(1.D0+c1))
     r2  = dsqrt((1.D0-c2)/(1.D0+c2))
      
     rout = dsqrt(sum((v1-vc)**2))*dsqrt(1.D0-c1**2)&
           *(dlog(1.D0/r2)-dlog(r1))

 
  end subroutine tri_rest


end module cg_matvec_mod
