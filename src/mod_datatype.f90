!************************************************************************!
!*  This module declares datatypes and data structures                  *!
!*                                                                      *!
!*  by Jia Shi, shijia1019@gmail.com                                    *!
!*                                                                      *!
!*  This module mainly deals with geometry of the model, and            *!
!*  also local basis, CG to DG, DG to CG maps                           *!
!************************************************************************!


!-------------------------------------------------------------------------
module datatype_mod
    
    use para_mod,                    only: rkind

    implicit none
    
    !----------------------------------------------------------------------
    ! geometry for reference
    type :: tetrahedron
        integer                          :: vert(4)
        integer                          :: edge(6)
        integer                          :: face(4)
        ! defined in the geometry mod
        integer                          :: EtoV(6,2)
        integer                          :: FtoE(4,3)
        integer                          :: FtoV(4,3)
        ! vertice ids
        integer                          :: vord(4)   
        real(kind=rkind)                 :: vert_crs(3,4)
        real(kind=rkind)                 :: sJac(4)
        real(kind=rkind)                 :: n(3,4)            
        ! for nodes
        ! nodes status: 1 verts, 2 edge, 3, face, 4 interior nodes 
        integer, allocatable             :: ns(:)              ! pNp
        real(kind=rkind), allocatable    :: nodes(:,:)         ! (3,pNp)
        real(kind=rkind), allocatable    :: V3D(:,:)           ! (pNp,pNp)
        real(kind=rkind), allocatable    :: invV(:,:)          ! (pNp,pNp)
        real(kind=rkind), allocatable    :: MassM(:,:)         ! (pNp,pNp)
        real(kind=rkind), allocatable    :: Drst(:,:,:)        ! (pNp,pNp,3)
        real(kind=rkind), allocatable    :: D1(:,:)            ! (pNp,pNp)
        real(kind=rkind), allocatable    :: D2(:,:)            ! (pNp,pNp)
        real(kind=rkind), allocatable    :: D3(:,:)            ! (pNp,pNp)
        ! for reference surface
        integer, allocatable             :: Fmask(:,:)         ! (Nfp,4)
        real(kind=rkind), allocatable    :: MassF(:,:,:)       ! (Nfp,Nfp,4)
        real(kind=rkind)                 :: Fscale(4)      
    end type tetrahedron

    ! separate fuild and solid
    type :: elements
        integer                          :: ClNele
        integer, allocatable             :: Clelist(:)        ! (ClNele) 
        integer, allocatable             :: loc_lab(:,:)      ! (4,ClNele)
        integer, allocatable             :: pid(:)            ! (ClNele) 
        integer, allocatable             :: loc_ids(:,:)      ! (pNp,ClNele)
        real(kind=rkind), allocatable    :: loc_nods(:,:)     ! (3,pNp*ClNele)
        real(kind=rkind), allocatable    :: loc_dtJ(:)        ! (ClNele)
        real(kind=rkind), allocatable    :: loc_invJ(:,:,:)   ! (3,3,ClNele)
        real(kind=rkind), allocatable    :: loc_rho(:,:)      ! (4,ClNele)
        real(kind=rkind), allocatable    :: loc_cen(:,:)      ! (3,ClNele)
        real(kind=rkind), allocatable    :: val(:)            ! (ClNele)
        real(kind=rkind), allocatable    :: lcR(:)            ! (ClNele)
        real(kind=rkind), allocatable    :: Dx(:,:)           ! (pNp,ClNele)
        real(kind=rkind), allocatable    :: Dy(:,:)           ! (pNp,ClNele)
        real(kind=rkind), allocatable    :: Dz(:,:)           ! (pNp,ClNele)
    end type elements

    type :: surface
        integer                          :: lsf
        integer, allocatable             :: dist(:)           ! (nproc+1)
        integer, allocatable             :: s2e(:,:)          ! (2,lsf) (-,+)
        integer, allocatable             :: snum(:)           ! (lNele)
        integer, allocatable             :: sdst(:)           ! (lNele+1)
        integer, allocatable             :: pid(:)            ! (lsf) (-) 
        integer, allocatable             :: ids(:)            ! (lsf)  
        integer, allocatable             :: inx(:)            ! (lsf) (-)
        integer, allocatable             :: vtx(:,:)          ! (3,lsf) (-)
        ! stat of the surface F/S: 0 normal; 1 Fluid at FS; 2 FS, S-; 3 FS, F-
        integer, allocatable             :: sta(:)            ! (lsf)
        real(kind=rkind), allocatable    :: cen(:,:)          ! (3,lsf)
        real(kind=rkind), allocatable    :: drt(:,:)          ! (3,lsf) - to +
        real(kind=rkind), allocatable    :: rho(:)            ! (lsf)
        real(kind=rkind), allocatable    :: val(:)            ! (lsf)
        real(kind=rkind), allocatable    :: lcR(:)            ! (lsf)
    end type surface 

    type :: edinfo
        integer                          :: nedg
        integer, allocatable             :: eddist(:)      ! (nporc+1)
        ! edge to its two vertices
        integer, allocatable             :: ed2v(:,:)      ! (nedg,2)
        integer, allocatable             :: ed2vstat(:,:)  ! (nedg,2)
        integer, allocatable             :: ed2vpid(:,:)   ! (nedg,2)
        ! edge list
        integer, allocatable             :: list(:)        ! (nedg) 
        ! edge loc elm
        integer, allocatable             :: locelm(:)      ! (ClNele) 
        ! edge rank
        integer, allocatable             :: pid(:)         ! (nedg)
        ! total size
        integer                          :: nsiz   
        ! edge to all nodes
        integer, allocatable             :: ed2n(:)        ! nsiz      
        integer, allocatable             :: ed2npid(:)     ! nsiz      
        integer, allocatable             :: ed2nstat(:)    ! nsiz      
        ! edge to elements 
        integer                          :: ed2esiz
        integer, allocatable             :: edist(:)       ! (nedg+1)
        integer, allocatable             :: ed2e(:)        ! (ed2esiz) 
    end type edinfo

    type :: v2vmap
        integer                          :: nvtx
        integer, allocatable             :: vtxdist(:)     ! (nproc+1)
        integer, allocatable             :: vlist(:)       ! (nvtx)
        integer, allocatable             :: part(:)        ! (nvtx)
        integer, allocatable             :: vstat(:)       ! (nvtx)
        ! covered vlist for newd & newp
        integer                          :: cnvtx
        integer, allocatable             :: Cvlist(:)      ! (cnvtx)
        integer, allocatable             :: Cvpid(:)       ! (cnvtx)
        ! elem info (add JS 0430/2018)
        integer                          :: nele
        integer, allocatable             :: edist(:)       ! (nproc) 
        integer, allocatable             :: elist(:)       ! (nele) 
        integer, allocatable             :: e2v(:)         ! (nele*4) 
        ! v2e info
        integer                          :: cnt 
        integer, allocatable             :: v2edist(:)     ! (nvtx+1)
        integer, allocatable             :: v2e(:)         ! (cnt)
        ! v2v info
        integer                          :: v2vsiz
        integer, allocatable             :: nv(:)          ! (nvtx)
        integer, allocatable             :: v2vdist(:)     ! (nvtx+1)
        integer, allocatable             :: v2v(:)         ! (v2vsiz)       
        integer, allocatable             :: v2vpid(:)      ! (v2vsiz)
        integer, allocatable             :: v2vstat(:)     ! (v2vsiz)
     end type v2vmap


    ! information for CG 
    type :: mesh
        ! mpi information
        integer                          :: rank,nproc,comm 
        ! global variables 
        ! read from meshfile
        integer                          :: Nvert,Ntet,Nfedge,Nedge,Nsrf
        real(kind=rkind)                 :: xmin,xmax,ymin,ymax,zmin,zmax
        ! analyze the mesh information
        real(kind=rkind)                 :: edmin,edmax,edavg
        ! local input data: build vtx-vtx map
        type(v2vmap)                     :: org,edtmp,orgn,new!,newd,newp
        ! reorder vlist
        integer, allocatable             :: vord(:)           ! (new%nvtx)
        integer, allocatable             :: rord(:)           ! (new%nvtx) 
        ! total number of local vertices
        integer                          :: cnvtx     
        integer, allocatable             :: Cvlist(:)         ! (cnvtx)
        integer, allocatable             :: Cvpid(:)          ! (cnvtx)
        ! locations
        real(kind=rkind), allocatable    :: Cv_crs(:,:)       ! (3,cnvtx)
        !---------------------- for the fluid-solid ---------------------------!
        ! status fluid-solid 
        logical                          :: fsexist
        !--------------------------local variables-----------------------------!
        ! from domain decomposition
        integer                          :: lNele
        integer, allocatable             :: lelist(:)         ! lNele
        real(kind=rkind), allocatable    :: sloc(:,:)         ! (3,lNele)
        type(surface)                    :: srf,csrf,esrf
        ! for parallelism in order to reduce communications
        ! elements contain all nodes in lNele (for CG) 
        integer                          :: ClNele
        integer, allocatable             :: Clelist(:)        ! ClNele
        integer, allocatable             :: ClNeigh(:,:)      ! (4,ClNele)
        integer, allocatable             :: ClNids(:)         ! (ClNele)
        integer, allocatable             :: CNOids(:)         ! (ClNele)
        ! original data
        integer, allocatable             :: loc_t2v(:,:)      ! (pNp,ClNele)
        ! local vtx ids
        integer, allocatable             :: lt2vid(:)         ! (pNp*ClNele)
        real(kind=rkind), allocatable    :: loc_vtx(:,:)      ! (3,4*ClNele)          
        real(kind=rkind), allocatable    :: loc_nods(:,:)     ! (3,pNp*ClNele)
        ! local information
        type(elements)                   :: s,f,efmm
        real(kind=rkind), allocatable    :: loc_Jac(:,:,:)    ! (3,3,ClNele)
        real(kind=rkind), allocatable    :: loc_invJ(:,:,:)   ! (3,3,ClNele)
        real(kind=rkind), allocatable    :: loc_detJ(:)       ! (ClNele)
        ! for internal jumps
        real(kind=rkind), allocatable    :: loc_n(:,:,:)      ! (3,4,ClNele) normal
        real(kind=rkind), allocatable    :: loc_sJac(:,:)     ! (4,ClNele)
        real(kind=rkind), allocatable    :: loc_Fscale(:,:)   ! (4,ClNele)
    end type mesh

  


end module datatype_mod
