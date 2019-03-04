!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!file : parmetis_interface.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parmetis_interface
  use, intrinsic                             :: ISO_C_BINDING

 !interface
 !  integer(C_INT) function METIS_SetDefaultOptions(opts)&
 !             bind(C,name="METIS_SetDefaultOptions")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_SetDefaultOptions
 !end interface

 !interface
 !  integer(C_INT) function METIS_NodeND(nvtxs,xadj,adjncy,vwgt,&
 !                          opts,perm,iperm) bind(C,name="METIS_NodeND")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none
 !  integer(C_INT)                            :: nvtxs
 !  integer(C_INT), dimension(*)              :: xadj,adjncy,perm,iperm
 !  type(C_PTR), value                        :: vwgt
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_NodeND
 !end interface

 !interface
 !  integer(C_INT) function METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsiz,&
 !     ncommon,nparts,tpwgts,opts,objval,epart,npart) bind(C,name="METIS_PartMeshDual")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: ne,nn,ncommon,nparts,objval
 !  integer(C_INT), dimension(*)              :: eptr,eind,epart,npart
 !  type(C_PTR), value                        :: vwgt,vsiz,tpwgts
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_PartMeshDual
 !end interface

 !interface
 !  integer(C_INT) function METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsiz,&
 !            nparts,tpwgts,opts,objval,epart,npart) bind(C,name="METIS_PartMeshNodal")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: ne,nn,nparts,objval
 !  integer(C_INT), dimension(*)              :: eptr,eind,epart,npart
 !  type(C_PTR), value                        :: vwgt,vsiz,tpwgts
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_PartMeshNodal
 !end interface

 !interface
 !  integer(C_INT) function METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt,vsiz,adjwgt,&
 !            nparts,tpwgts,ubvec,opts,objval,part) bind(C,name="METIS_PartGraphKway")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: nvtxs,ncon,nparts,objval
 !  integer(C_INT), dimension(*)              :: xadj,adjncy,part
 !  type(C_PTR), value                        :: vwgt,vsiz,adjwgt,tpwgts,ubvec
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_PartGraphKway
 !end interface


 !interface
 !  integer(C_INT) function METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsiz,adjwgt,&
 !            nparts,tpwgts,ubvec,opts,objval,part) bind(C,name="METIS_PartGraphRecursive")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: nvtxs,ncon,nparts,objval
 !  integer(C_INT), dimension(*)              :: xadj,adjncy,part
 !  type(C_PTR), value                        :: vwgt,vsiz,adjwgt,tpwgts,ubvec
 !  integer(C_INT)                            :: opts(0:40)
 !  end function METIS_PartGraphRecursive
 !end interface

 !interface
 !  integer(C_INT) function METIS_MeshToDual(ne,nn,eptr,eind,&
 !                      ncommon,numflag,xadj,adjncy) bind(C,name="METIS_MeshToDual")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: ne,nn,ncommon,numflag
 !  integer(C_INT), dimension(*)              :: eptr,eind
 !  type(C_PTR), intent(out)                  :: xadj,adjncy
 !  end function METIS_MeshToDual
 !end interface

 !interface
 !  integer(C_INT) function METIS_MeshToNodal(ne,nn,eptr,eind,&
 !                              numflag,xadj,adjncy) bind(C,name="METIS_MeshToNodal")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: ne,nn,numflag
 !  integer(C_INT), dimension(*)              :: eptr,eind
 !  type(C_PTR), intent(out)                  :: xadj,adjncy
 !  end function METIS_MeshToNodal
 !end interface


 interface
   integer(C_INT) function ParMETIS_V3_PartMeshKway(elmdist,eptr,eind,&
         elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,tpwgts,ubvec,&
         opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartMeshKway")
   use, intrinsic                            :: ISO_C_BINDING
   implicit none 
   integer(C_INT)                            :: wgtflag,numflag,ncon,comm,&
                                                nparts,ncommonnodes,edgecut
   integer(C_INT), dimension(*)              :: elmdist,eptr,eind,elmwgt
   !type(C_PTR), value                        :: elmwgt
   real(C_FLOAT), dimension(*)               :: tpwgts,ubvec
   integer(C_INT)                            :: opts(0:2)
   !type(C_PTR), intent(out)                  :: part
   integer(C_INT), dimension(*)              :: part
   end function ParMETIS_V3_PartMeshKway
 end interface

 !interface
 !  integer(C_INT) function ParMETIS_V3now_PartMeshKway(elmdist,eptr,eind,&
 !        elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,tpwgts,ubvec,&
 !        opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartMeshKway")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: wgtflag,numflag,ncon,comm,&
 !                                               nparts,ncommonnodes,edgecut
 !  integer(C_INT), dimension(*)              :: elmdist,eptr,eind
 !  type(C_PTR), value                        :: elmwgt
 !  real(C_FLOAT), dimension(*)               :: tpwgts,ubvec!,elmwgt
 !  integer(C_INT)                            :: opts(0:2)
 !  !type(C_PTR), intent(out)                  :: part
 !  integer(C_INT), dimension(*)              :: part
 !  end function ParMETIS_V3now_PartMeshKway
 !end interface

 !interface
 !  integer(C_INT) function ParMETIS_V3now_PartKway(vtxdist,xadj,adjncy,&
 !                 vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,&
 !            opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartKway")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: wgtflag,numflag,ncon,comm,&
 !                                               nparts,edgecut
 !  integer(C_INT), dimension(*)              :: vtxdist,xadj,adjncy
 !  type(C_PTR), value                        :: vwgt,adjwgt
 !  real(C_FLOAT), dimension(*)               :: tpwgts,ubvec!,elmwgt
 !  integer(C_INT)                            :: opts(0:2)
 !  !type(C_PTR), intent(out)                  :: part
 !  integer(C_INT), dimension(*)              :: part
 !  end function ParMETIS_V3now_PartKway
 !end interface

 !interface
 !  integer(C_INT) function ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,&
 !                 vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,&
 !            opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartKway")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_INT)                            :: wgtflag,numflag,ncon,comm,&
 !                                               nparts,edgecut
 !  integer(C_INT), dimension(*)              :: vtxdist,xadj,adjncy,vwgt
 !  type(C_PTR), value                        :: adjwgt
 !  real(C_FLOAT), dimension(*)               :: tpwgts,ubvec
 !  integer(C_INT)                            :: opts(0:2)
 !  !type(C_PTR), intent(out)                  :: part
 !  integer(C_INT), dimension(*)              :: part
 !  end function ParMETIS_V3_PartKway
 !end interface


 interface
   integer(C_int64_t) function ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,&
                  vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,&
             opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartKway")
   use, intrinsic                             :: ISO_C_BINDING
   implicit none 
   integer(C_int64_t)                         :: wgtflag,numflag,ncon,comm,&
                                                nparts,edgecut
   integer(C_int64_t), dimension(*)           :: vtxdist,xadj,adjncy,vwgt
   type(C_PTR), value                         :: adjwgt
   real(C_float), dimension(*)                :: tpwgts,ubvec
   integer(C_int64_t)                         :: opts(0:2)
   integer(C_int64_t), dimension(*)           :: part
   end function ParMETIS_V3_PartKway
 end interface

 !interface
 !  integer(C_int64_t) function ParMETIS_V3now_PartKway(vtxdist,xadj,adjncy,&
 !                 vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,&
 !            opts,edgecut,part,comm) bind(C,name="ParMETIS_V3_PartKway")
 !  use, intrinsic                            :: ISO_C_BINDING
 !  implicit none 
 !  integer(C_int64_t)                        :: wgtflag,numflag,ncon,comm,&
 !                                               nparts,edgecut
 !  integer(C_int64_t), dimension(*)          :: vtxdist,xadj,adjncy
 !  type(C_PTR), value                        :: vwgt,adjwgt
 !  real(C_FLOAT), dimension(*)               :: tpwgts,ubvec!,elmwgt
 !  integer(C_int64_t)                        :: opts(0:2)
 !  !type(C_PTR), intent(out)                  :: part
 !  integer(C_int64_t), dimension(*)          :: part
 !  end function ParMETIS_V3now_PartKway
 !end interface

end module parmetis_interface
