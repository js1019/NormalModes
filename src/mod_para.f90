!***********************************************************************!
!*  This module declares and reads in parameters variables             *!
!*  Read parameters from file ../bin/global_conf                       *!
!*  Organized by Jia Shi, Nov. 18, 2015                                *!
!***********************************************************************!  
               
!***********************************************************************! 
!*  Since we focus on the eigenvalue problems,                         *!
!*  it is better to keep rkind = 8                                     *!
!***********************************************************************!

!***********************************************************************! 
!*  Change I/O; it is better to use a single process to read           *! 
!*  basic parameters and sends to other processes                      *!
!*  JS 04/25/2018                                                      *!
!***********************************************************************!


!------------------------------------------------------------------------
module para_mod
!------------------------------------------------------------------------
    use mpi
    use string_mod,         only : string_conf

    implicit none
    
    integer, parameter           :: rkind          = 8    ! 8
!-----------------------------------------------------------------------
    type :: singleelm
       ! pOrder: polynomial order
       integer                   :: pOrder         = 1,pNp,Nfp
       ! number of point on one element
       integer                   :: Npfp,Npep,Nip
    end type singleelm

 
    type :: parameterinput
       character(len=1024)       :: outputdir      = "~/"
       character(len=1024)       :: inputdir       = "~/"
       character(len=300)        :: basename
       character(len=1024)       :: fvpt,fvst,frhot,fgrav
       character(len=1024)       :: fhd,fele,fnode,fneigh
       ! add for output files
       character(len=1024)       :: fvlist,fvstat,fvloct,fvdata
       integer                   :: filefid        = 1503
       integer                   :: logfid         = 109
       integer                   :: JOB            = 1
       ! new added  
       integer                   :: comm 
       real(kind=rkind)          :: TOL            = 1.0D-8
       real                      :: lowfreq        = 1.0D-2
       real                      :: upfreq         = 2.0D0
       real                      :: startTime, finalTime, dt

       logical                   :: CGMethod       = .TRUE.      
       logical                   :: buildEdgeOn    = .FALSE.
       logical                   :: buildtetJac    = .FALSE.
       logical                   :: buildtetn      = .FALSE.
       logical                   :: buildBasnodeOn = .FALSE.
       logical                   :: DEBUG
       logical                   :: selfG          = .FALSE.
       logical                   :: phi1           = .FALSE.
       ! separate fluid & solid elements
       type(singleelm)           :: s,f
       logical                   :: mix            = .FALSE. !.TRUE.
    end type parameterinput


    type(parameterinput), save   :: pin

    ! please be aware that these parameters might be used somewhere else
    private
    public                       :: rkind,pin,para_init
    
!------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------

  subroutine para_init(fnm_conf)
    implicit none

    character(len=*)             :: fnm_conf
    !character(len=1024)          :: finput,foutput
    !character(len=300)           :: fbase
    integer                      :: fid,ierr,rank,mpi_size
    logical                      :: alive

    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, mpi_size, ierr)

    pin%comm = mpi_comm_world
    !print*, pin%comm, mpi_comm_world

    ! open file to read input parameter
    fid = 1001
    inquire(file=trim(fnm_conf),exist=alive)
    if(.not. alive)then
       write(*,*) 'Error 1001: File "',trim(fnm_conf),'" does not exist.'
       stop
    endif

    if (rank.eq.0) then  
       open(fid,file=trim(fnm_conf),status="old")
       call string_conf(fid,1,'JOB',      2,      pin%JOB)
       call string_conf(fid,1,'basename', 2, pin%basename)
       call string_conf(fid,1,'inputdir', 2, pin%inputdir)
       call string_conf(fid,1,'outputdir',2,pin%outputdir)
       !call string_conf(fid,1,'CGMethod', 2, pin%CGMethod)
       call string_conf(fid,1,'lowfreq',  2,  pin%lowfreq) 
       call string_conf(fid,1,'upfreq',   2,   pin%upfreq)
       call string_conf(fid,1,'pOrder',   2, pin%s%pOrder)
       close(fid)
    endif 

    ! pass the information
    call mpi_bcast(pin%JOB,1,mpi_integer,0,pin%comm,ierr)
    call mpi_bcast(pin%basename,300,mpi_char,0,pin%comm,ierr)
    !print*,trim(pin%basename),rank
    call mpi_bcast(pin%inputdir, 1024,mpi_char,0,pin%comm,ierr)
    call mpi_bcast(pin%outputdir,1024,mpi_char,0,pin%comm,ierr)
    call mpi_bcast(pin%lowfreq,1,mpi_real,0,pin%comm,ierr)
    call mpi_bcast(pin%upfreq, 1,mpi_real,0,pin%comm,ierr)
    call mpi_bcast(pin%s%pOrder,1,mpi_integer,0,pin%comm,ierr)

    !print*, basename 
    if (pin%JOB.eq.20) then 
       !todo future work mixed with different orders
       pin%mix = .true.
    else
       pin%mix = .false. 
    endif

    if (pin%JOB.eq.2) then 
       pin%selfG = .true.
    elseif (pin%JOB.eq.3) then
       pin%selfG = .true.
       pin%phi1  = .true.
    endif


    if (pin%CGMethod) then
       pin%buildEdgeOn    = .FALSE.
       pin%buildtetJac    = .TRUE.
       pin%buildBasnodeOn = .TRUE.
       pin%buildtetn      = .TRUE.
       call fwd_modelfilesname()
       if (pin%selfG) then 
          call fwd_gravity()
       endif
    endif 

    ! number of nodes on one tet
    pin%s%pNp = (pin%s%pOrder+1)*(pin%s%pOrder+2)*(pin%s%pOrder+3)/6
    ! number of nodes on the face
    pin%s%Nfp = (pin%s%pOrder+1)*(pin%s%pOrder+2)/2
    
    ! number of nodes purely on edges
    pin%s%Npep = 6*(pin%s%pOrder-1)
    ! number of node purely on faces
    pin%s%Npfp = 2*(pin%s%pOrder-1)*(pin%s%pOrder-2)
    ! number of interior nodes
    pin%s%Nip  = (pin%s%pOrder-3)*(pin%s%pOrder-2)*(pin%s%pOrder-1)/6


    ! Todo check fluid part
    if (pin%mix.and.pin%s%pOrder.gt.1) then
       pin%f%pOrder = pin%s%pOrder - 1
    else
       pin%f%pOrder = pin%s%pOrder
    endif
    pin%f%pNp = (pin%f%pOrder+1)*(pin%f%pOrder+2)*(pin%f%pOrder+3)/6
    pin%f%Nfp = (pin%f%pOrder+1)*(pin%f%pOrder+2)/2
    pin%f%Npep = 6*(pin%f%pOrder-1)
    pin%f%Npfp = 2*(pin%f%pOrder-1)*(pin%f%pOrder-2)
    pin%f%Nip  = (pin%f%pOrder-3)*(pin%f%pOrder-2)*(pin%f%pOrder-1)/6

    !print*,pin%s%pOrder,pin%f%pOrder,rank


    ! different TOL with different accuracy 
    if (rkind <= 4) then 
       pin%TOL = 1.0E-6
    else
       pin%TOL = 1.0D-8
    endif

    !if (abs(omega) < 1.0D-3) omega = 1.0D-3 
        
    !call mpi_barrier(pin%comm,ierr)

    if (rank == 0) then 
       print*,'========================================================================='  
       print*, 'Read in parameters from global_conf'
       print*, 'JOB = ', int(pin%JOB,2)
       if (pin%CGmethod) then 
          print*, 'apply the Continuous Galerkin finite element method'
       endif
       print*, 'mesh name:   ', trim(pin%basename) 
       print*, 'input directory: ',  trim(pin%inputdir)
       print*, 'output directory: ', trim(pin%outputdir)
       print*, 'polynomial order', int(pin%s%pOrder,2)
       print*, 'mpi_size', mpi_size
       print*, 'lower frequency in mHz', pin%lowfreq
       print*, 'upper frequency in mHz', pin%upfreq
       print*, 'header file name: ', trim(pin%fhd) 
       print*, 'ele file name: ', trim(pin%fele) 
       print*, 'node file name: ', trim(pin%fnode) 
       print*, 'neigh file name: ', trim(pin%fneigh) 
       print*, 'Vp file name: ', trim(pin%fvpt) 
       print*, 'Vs file name: ', trim(pin%fvst)
       print*, 'rho file name: ', trim(pin%frhot)
       if (pin%selfG) then 
          print*, 'apply gravitational acceleration'
          print*, 'reference gravitational acceleration file name: ',&
                  trim(pin%fgrav)
       endif
       print*,'=========================================================================='  
    endif 

    if (pin%lowfreq .ge. pin%upfreq) then 
       print*, 'error: please check input frequencies' 
       stop
    endif 

  end subroutine para_init

  subroutine fwd_modelfilesname()
     ! construct the model file names 
     character(len=1024)                      :: fvpt,fvst,frhot,str_p
     character(len=1024)                      :: fupf,flowf,np,job

     integer                                  :: rank,mpisize,ierr    

     write(str_p,*) pin%s%pOrder; str_p = adjustl(str_p)
     
     ! mesh file information
     pin%fhd = trim(adjustl(pin%inputdir))//&
               trim(adjustl(pin%basename))//'_mesh.header'
     pin%fhd = trim(pin%fhd)

     pin%fele = trim(adjustl(pin%inputdir))//&
                trim(adjustl(pin%basename))//'_ele.dat'
     pin%fele = trim(pin%fele)

     pin%fnode = trim(adjustl(pin%inputdir))//&
                 trim(adjustl(pin%basename))//'_node.dat'
     pin%fnode = trim(pin%fnode)

     pin%fneigh = trim(adjustl(pin%inputdir))//&
                  trim(adjustl(pin%basename))//'_neigh.dat'
     pin%fneigh = trim(pin%fneigh)

     ! models
     fvpt  = trim(adjustl(pin%inputdir))//trim(adjustl(pin%basename))&
             //'_vp_'//'pod_'//trim(str_p)//'_true.dat'
     pin%fvpt  = trim(fvpt)

     fvst  = trim(adjustl(pin%inputdir))//trim(adjustl(pin%basename))&
             //'_vs_'//'pod_'//trim(str_p)//'_true.dat'
     pin%fvst  = trim(fvst)

     frhot = trim(adjustl(pin%inputdir))//trim(adjustl(pin%basename))&
             //'_rho_'//'pod_'//trim(str_p)//'_true.dat'
     pin%frhot = trim(frhot)
     
     ! output files
     call mpi_comm_rank(pin%comm, rank, ierr)
     call mpi_comm_size(pin%comm, mpisize, ierr)
     write(np,*) mpisize; np = adjustl(np)

     pin%fvlist = trim(pin%outputdir)//trim(pin%basename)&
                  //'_pod'//trim(str_p)//'_np'//trim(np)//'_vlist.dat'
     pin%fvlist = trim(pin%fvlist)
     pin%fvstat = trim(pin%outputdir)//trim(pin%basename)&
                  //'_pod'//trim(str_p)//'_np'//trim(np)//'_vstat.dat'
     pin%fvstat = trim(pin%fvstat)
     pin%fvloct = trim(pin%outputdir)//trim(pin%basename)&
                  //'_pod'//trim(str_p)//'_np'//trim(np)//'_vloct.dat'
     pin%fvloct = trim(pin%fvloct)


     write(fupf,*)  real(pin%upfreq,4);  fupf  = adjustl(fupf) 
     write(flowf,*) real(pin%lowfreq,4); flowf = adjustl(flowf) 
     write(job,*) pin%JOB; job = adjustl(job)

     pin%fvdata = trim(pin%outputdir)//trim(pin%basename)//'_JOB'//trim(job)&
              //'_pod'//trim(str_p)//'_np'//trim(np)//'_'&
              //trim(flowf)//'_'//trim(fupf)
     pin%fvdata = trim(pin%fvdata)

     !if (rank.eq.0.and..false.) then
     if (rank.eq.0.and..false.) then
        print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        print*, 'output files'
        print*, trim(pin%fvlist)
        print*, trim(pin%fvstat)
        print*, trim(pin%fvloct)
        print*, trim(pin%fvdata)
        print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     endif
 
  end subroutine fwd_modelfilesname 

  subroutine fwd_gravity()
    ! construct file name of gravitational acceleration
    character(len=1024)                       :: fgrav,str_p 

    write(str_p,*) pin%s%pOrder; str_p = adjustl(str_p)
    !fgrav = trim(adjustl(inputdir))//trim(adjustl(basename))&
    !       //'_pod_'//trim(str_p)//'_potential_acceleration_true'//'.dat'
    pin%fgrav = trim(adjustl(pin%inputdir))//trim(adjustl(pin%basename))&
           //'_pod_'//trim(str_p)//'_potential_acceleration_true'//'.dat'

    pin%fgrav = trim(pin%fgrav)
  end subroutine fwd_gravity


end module para_mod

