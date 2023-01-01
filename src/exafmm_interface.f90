!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!file : exafmm_interface.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module exafmm_interface
  use, intrinsic                   :: ISO_C_BINDING
  use mpi
  use omp_lib
  !use para_mod,                 only: rkind,pin 
  !use geometry_mod,             only: unstrM 

!contains

  !subroutine pnm_setup_exafmm()
  !   integer                                 :: i,j
  !   integer                                 :: ncrit,threads
  !   real(kind=rkind)                        :: eps0,kreal,kimag     
  !   character(len=1024)                     :: path 
  !   real(kind=rkind), allocatable           :: ws(:),wr(:)

  !   eps0 = 0.0D0; kreal = 1.0; kimag = 0.1 
  !   path = './'; ncrit = 10; threads = 1
  !  
  !   allocate(ws(unstrM%lNele),wr(unstrM%lNele))
  !   
  !   call FMM_Init(eps0,kreal,kimag,ncrit,threads,path,&
  !        ws,unstrM%sloc(1,:),unstrM%sloc(2,:),unstrM%sloc(3,:),&
  !           unstrM%sloc(1,:),unstrM%sloc(2,:),unstrM%sloc(3,:),wr)


  !end subroutine pnm_setup_exafmm

  interface
    integer(C_INT) function fmm_init(eps,kreal,kimag,ncrit,threads,path,&
              nb,xb,yb,zb,vb,nv,xv,yv,zv,vv) bind(C,name="fmm_init_")
     use, intrinsic                        :: ISO_C_BINDING
     implicit none 
     integer(C_INT), value                 :: ncrit,threads,nb,nv 
     character                             :: path     
     real(C_double)                        :: eps,kreal,kimag
     real(C_double), dimension(*)          :: xb,yb,zb,vb,xv,yv,zv,vv
    end function fmm_init
  end interface

end module exafmm_interface
