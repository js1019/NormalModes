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
    if(unstrM%rank.eq.0) print*,'check number of threads', mymatvec%nthd
    
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
    MLAN = 1000; LANSTEP = 2000
    TOL = 1.0D-12; CHEBTYPE = 2
    call pEVSL_LANBOUNDS_F90(pevslB, MLAN, LANSTEP,TOL, LMIN, LMAX)
    if (mymatvec%rank.eq.0) print*, 'B cond. number', LMAX/LMIN,LMIN,LMAX

    if (pin%s%porder.eq.2) then 
       CHEBDEG = 45 
    else 
       CHEBDEG = 25 !25 !30
    endif
    call pEVSL_SETUP_CHEBITER_F90(LMIN,LMAX,CHEBDEG,mymatvec%sBV,mymatvec%chebB)
    
    call PEVSL_FINISH_F90(pevslB)


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
       TOL = 1.0D-12; CHEBTYPE = 2
       call pEVSL_LANBOUNDS_F90(pevslAp, MLAN, LANSTEP,TOL, LMIN, LMAX)

       call mpi_barrier(mymatvec%comm,ierr)
       if (mymatvec%rank.eq.0) print*, 'Ap cond. number ', LMAX/LMIN, LMIN, LMAX

       if (pin%s%porder.eq.2) then 
          CHEBDEG = max(min(int(LMAX/LMIN),15)*3,100) 
       else 
          CHEBDEG = 25 !25 !30
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
             mymatvec%B%diag(j-lsiz) = 1.0D0/dsqrt(CGM%B%val(k))
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
    !print*,minval(mymatvec%B%val), maxval(mymatvec%B%val),unstrM%rank 
 

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

  end subroutine sparseAV


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
 
    !!print*, norm2(w)
    w = w + (w0+w1)*mymv%B%diag
    !w = w0*mymv%Bdiag
  end subroutine sparsefsAV





end module cg_matvec_mod
