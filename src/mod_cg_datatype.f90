!*****************************************************************!
!* this module contains data type for CG models, matrix and RHS   !
!* By Jia Shi                                                     !  
!*****************************************************************!
module cg_datatype_mod
   
   use ISO_C_BINDING  
   use para_mod,                    only: rkind
   use datatype_mod,                only: surface,elements 
   implicit none

  
   type :: modelcoeff
       ! general information
       integer                               :: p_vp,p_vs,p_rho
       integer                               :: siz,Gsiz,nmps
       ! some other consideration
       logical                               :: FreeBC = .true. ! BC
       ! local information 
       real(kind=rkind), allocatable         :: coeff_loc(:,:)
       !!---------------------- for the fluid-solid --------------------------!
       ! num of f-s interfaces 
       integer                               :: nsf    = 0             
       !--------------------- for self gravitation --------------------------!
       ! local gravitational acceleration
       real(kind=rkind), allocatable         :: g0(:,:,:)       ! (pNp,3,ClNele)
       ! density jumps
       integer                               :: lsf
       real(kind=rkind), allocatable         :: srho(:,:)       ! (2,lsf)
       real(kind=rkind), allocatable         :: svs(:,:)        ! (2,lsf)
       real(kind=rkind), allocatable         :: erho(:,:)       ! (2,lsf)
       real(kind=rkind), allocatable         :: evs(:,:)        ! (2,lsf)
   end type modelcoeff

   type COOmat
       integer                               :: siz,Gsiz 
       integer*8                             :: NNZ,GNNZ
       integer, allocatable                  :: sizdist(:)    ! nproc+1
       ! number of columns
       integer, allocatable                  :: rowdist(:)    ! siz+1
       !integer, allocatable                  :: row(:)        ! NNZ 
       integer, allocatable                  :: col(:)        ! NNZ
       real(kind=rkind), allocatable         :: val(:)        ! NNZ 
       real(kind=rkind), allocatable         :: diag(:)       ! siz
       ! for vertices
       integer, allocatable                  :: rstt(:)       ! nvtx 
       ! number of variables 
       integer, allocatable                  :: rnum(:)       ! nvtx 
   end type COOmat

 
   type :: CGMatrix
       ! stat = .false.  we need to change the structure of the matrix
       logical                               :: stat   = .false. 
       ! set up for communication
       integer, allocatable                  :: vnum(:)       ! cnvtx
       integer, allocatable                  :: vstt(:)       ! cnvtx
       integer, allocatable                  :: pnum(:)       ! cnvtx
       integer, allocatable                  :: pstt(:)       ! cnvtx
       !------------------- separate A matrix 
       type(COOmat)                          :: Ad 
       type(COOmat)                          :: Ap 
       type(COOmat)                          :: E 
       type(COOmat)                          :: ET
       type(COOmat)                          :: A 
       type(COOmat)                          :: B 
   end type CGMatrix

   type :: fmmpart
       integer                               :: nv
       integer, allocatable                  :: pid(:)         ! nv
       integer, allocatable                  :: ids(:)         ! nv
       real(kind=rkind), allocatable         :: cen(:,:)       ! (3,nv)
       real(kind=rkind), allocatable         :: val(:)         ! nv
       real(kind=rkind), allocatable         :: lcR(:)         ! lth
   end type fmmpart

   type :: vector
       integer                               :: lth,stt 
       integer, allocatable                  :: pid(:)         ! lth
       integer, allocatable                  :: ids(:)         ! lth
       integer, allocatable                  :: nds(:)         ! nv/3
       integer, allocatable                  :: lds(:)         ! nv/3
       real(kind=rkind), allocatable         :: val(:)         ! lth
   end type vector  

   ! everything we need for mat-vec
   type :: mvparameters
       integer                               :: rank,comm
       integer                               :: nproc,nthd
       integer                               :: pbsiz,Gpbsiz,nsf
       integer*8                             :: sAV,sBV         ! matvec
       integer*8                             :: sAdV,sApV       ! matvec
       integer*8                             :: sEV,sETV        ! matvec
       integer*8                             :: chebB,chebAp    ! solve
       integer*8                             :: lspolBsqrt
       logical                               :: fsexist,sG,purefluid
       logical*1                             :: verbose
       integer                               :: ncf
       integer, allocatable                  :: flist(:)        ! (ncf)
       real(kind=rkind)                      :: Gcon
       real(kind=rkind)                      :: rhocut
       ! matrix information
       type(COOmat)                          :: A
       type(COOmat)                          :: B
       type(COOmat)                          :: Ad 
       type(COOmat)                          :: Ap 
       type(COOmat)                          :: E 
       type(COOmat)                          :: ET
       ! for fmm
       type(surface)                         :: srf,esrf
       type(elements)                        :: elm,celm
       type(fmmpart)                         :: din,dout
       type(vector)                          :: vin,vcv 
   end type mvparameters
 


end module cg_datatype_mod
