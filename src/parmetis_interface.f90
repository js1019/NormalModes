!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!file : parmetis_interface.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parmetis_interface
  use, intrinsic                             :: ISO_C_BINDING

  ! This is a bit of a hack, but it is necessary to invoke the Fortran
  ! API entry point for ParMETIS so that MPI handles are converted
  ! properly.  Hence, use the C symbol name parmetis_v3_partmeshkway
  ! (not the C ParMETIS API entry point ParMETIS_V3_PartMeshKway).


 interface
   integer(C_INT) function ParMETIS_V3_PartMeshKway(elmdist,eptr,eind,&
         elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,tpwgts,ubvec,&
         opts,edgecut,part,comm) bind(C,name="parmetis_v3_partmeshkway")
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



 interface
   integer(C_int64_t) function ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,&
                  vwgt,adjwgt,wgtflag,numflag,ncon,nparts,tpwgts,ubvec,&
             opts,edgecut,part,comm) bind(C,name="parmetis_v3_partkway")
   use, intrinsic                             :: ISO_C_BINDING
   implicit none 
   integer(C_int64_t)                         :: wgtflag,numflag,ncon,&
                                                nparts,edgecut
   integer(C_INT)                             :: comm
   integer(C_int64_t), dimension(*)           :: vtxdist,xadj,adjncy,vwgt
   type(C_PTR), value                         :: adjwgt
   real(C_float), dimension(*)                :: tpwgts,ubvec
   integer(C_int64_t)                         :: opts(0:2)
   integer(C_int64_t), dimension(*)           :: part
   end function ParMETIS_V3_PartKway
 end interface


end module parmetis_interface
