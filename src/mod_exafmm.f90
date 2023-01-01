
module exafmm_mod
  use mpi
  use omp_lib
  use ISO_C_BINDING
  use para_mod,                 only: rkind,pin 
  use geometry_mod,             only: unstrM 
  use cg_models_mod,            only: models
  !use exafmm_interface 

  implicit none

  interface
    subroutine fmm_init(eps,kreal,kimag,ncrit,threads,path,&
              nb,xb,yb,zb,vb,nv,xv,yv,zv,vv) bind(C,name="FMM_Init")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none 
       integer(C_INT), value             :: ncrit,threads,nb,nv 
       character                         :: path     
       real(C_DOUBLE), value             :: eps,kreal,kimag
       real*8, dimension(*)              :: xb,yb,zb,vb,xv,yv,zv,vv
    end subroutine fmm_init
  end interface

  interface 
     subroutine fmm_finalize() bind(C,name="FMM_Finalize") 
       use, intrinsic                        :: ISO_C_BINDING
       implicit none 
     end subroutine fmm_finalize
  end interface

  interface 
     subroutine fmm_partition(nb,xb,yb,zb,vb,&
        nv,xv,yv,zv,vv) bind(C,name="FMM_Partition")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none 
       integer                           :: nb,nv 
       real*8, dimension(*)              :: xb,yb,zb,vb,xv,yv,zv,vv
     end subroutine fmm_partition
  end interface
 
  interface 
     subroutine fmm_buildTree() bind(C,name="FMM_BuildTree")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none 
     end subroutine fmm_buildTree
  end interface

  interface  
     subroutine fmm_b2b(vi,vb,verbose) bind(C,name="FMM_B2B")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none
       !integer(C_INT)                    :: verbose 
       logical*1                         :: verbose
       real*8, dimension(*)              :: vi,vb
     end subroutine fmm_b2b
  end interface

  interface  
     subroutine fmm_v2b(vi,vb,verbose) bind(C,name="FMM_V2B")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none
       !integer(C_INT)                    :: verbose 
       logical*1                         :: verbose
       real*8, dimension(*)              :: vi,vb
     end subroutine fmm_v2b
  end interface

  interface 
     subroutine fmm_direct(ni,xi,yi,zi,vi,nj,xj,yj,zj,vj) bind(C,name="Direct")
       use, intrinsic                    :: ISO_C_BINDING
       implicit none
       integer(C_INT), value             :: ni,nj
       real*8, dimension(*)              :: vi,xi,yi,zi
       real*8, dimension(*)              :: vj,xj,yj,zj
     end subroutine fmm_direct
  end interface

contains

  subroutine pnm_check_exafmm(ntrd,comm)
     integer, intent(in)                 :: ntrd,comm
     logical*1                           :: verbose
     integer(C_INT)                      :: nb,nv,ncrit,threads
     real(C_DOUBLE)                      :: eps0,kreal,kimag     
     character                           :: path
     integer                             :: ni,nj,ngbl,n0,ierr
     real*8, allocatable, dimension(:)   :: vb,vv,vi
     real*8, allocatable, dimension(:)   :: xb,yb,zb
     real*8, allocatable, dimension(:)   :: xv,yv,zv

     integer                             :: i,j,k,l

     eps0 = 0.0D0; kreal = 0.0; kimag = 0.0 
     path = ''; ncrit = 1000; threads = ntrd
     !print*,'ntrd',threads   
     !ngbl = unstrM%Ntet 
     k = mod(unstrM%Ntet,unstrM%nproc)
     if (k.eq.0) then
        n0 = (unstrM%Ntet-k)/unstrM%nproc
     else
        n0 = (unstrM%Ntet-k)/unstrM%nproc + 1
     endif

     ngbl = int(max(n0,unstrM%lNele)*1.05)

     allocate(vb(ngbl),vv(ngbl))
     allocate(xb(ngbl),yb(ngbl),zb(ngbl))
     xb(1:unstrM%lNele) = unstrM%sloc(1,:)
     yb(1:unstrM%lNele) = unstrM%sloc(2,:)
     zb(1:unstrM%lNele) = unstrM%sloc(3,:)

     xv = xb; yv = yb; zv = zb

     !print*,xb(1)
     nv = unstrM%lNele; nb = nv
     call fmm_init(eps0,kreal,kimag,ncrit,threads,path,&
                          nb,xb,yb,zb,vb,nv,xv,yv,zv,vv)
     print*, 'init',nb,nv,ngbl,unstrM%rank
     
     ni = unstrM%lNele; nj = ni
     call fmm_partition(ni,xb,yb,zb,vb,nj,xv,yv,zv,vv) 
     print*, 'partition',ni,nj,ngbl,unstrM%rank


     call fmm_buildTree()

     ! compute sources ! NOT right
     vb = 0.0D0
     do i = 1,unstrM%lNele
        j = unstrM%lelist(i)
        vb(i) = unstrM%loc_detJ(j)*1.0D-6
     enddo
     
     allocate(vi(ngbl)); vi(1:unstrM%lNele) = 0.0D0
     call fmm_b2b(vi,vb,verbose)
     !call fmm_v2b(vi,vb,verbose)

     print*, maxval(vi(1:ni)),minval(vi(1:ni)),unstrM%rank
     call mpi_barrier(comm,ierr)
     !print*, vi(unstrM%lNele),vi(1+unstrM%lNele),unstrM%rank
     !if (unstrM%rank.eq.unstrM%nproc-1) then
     !  do i = 1,unstrM%lNele
     !     print*, i,vi(i) 
     !  enddo  
     !endif
     call fmm_finalize()


  end subroutine pnm_check_exafmm

  


end module exafmm_mod
