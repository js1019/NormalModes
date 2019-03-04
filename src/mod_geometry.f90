!*****************************************************************!
!*  This module converts inputed mesh file to CG mesh structure  *!
!*  Jia Shi, is still working on it                              *!
!*  This module constructs the geometry information needed       *!
!*  It is very important to get everything right                 *!
!*  TODO load balance                                            *!  
!*****************************************************************!
module geometry_mod
    
  use MPI
  use para_mod                   
  use utility_mod 
  use datatype_mod
  use cg_datatype_mod
 
  use ISO_C_BINDING
  use parmetis_interface 

  implicit none
  
  ! variable privacy
  private
  public                           :: unstrM,refs,reff,ednew
  public                           :: build_geometry
        
  type (mesh), save                :: unstrM
  type (tetrahedron), save         :: refs,reff
  type (edinfo), save              :: edge0,edge1,edall,ednew
  
  logical                          :: alive
  

 !------------------------------------------------------------------
 contains

 !--------------------------------------------------------------------
  subroutine Build_geometry()
    ! build up geometry information
    !character (len=1024)                :: filename
    integer                             :: i,j,k,l
    integer                             :: ierr,myrank,mynprocs
    
    ! set up mpi information for mesh
    unstrM%comm = pin%comm !MPI_COMM_WORLD  
    call MPI_COMM_RANK(unstrM%comm,myrank,ierr)
    call MPI_COMM_SIZE(unstrM%comm,mynprocs,ierr)
    unstrM%rank  = myrank 
    unstrM%nproc = mynprocs

    !-------------------------------------
    ! build reference tet
    call writelog("..Building reference tetrahedra", unstrM%comm)
    if (unstrM%rank == 0) print*,"Building reference tetrahedron"
    !call Build_refsolid()
    call Build_reference(pin%s%pOrder,pin%s%Nfp,pin%s%pNp,refs)
    call Build_reference(pin%f%pOrder,pin%f%Nfp,pin%f%pNp,reff)
  
    !do i = 1,4 
    !print*,refs%Fmask(:,i)
    !print*,refs%FtoV(i,:)
    !enddo
    !print*, sum(refs%MassM(:,2))

    !--------------------------------------
    ! new 04302018 JS
    ! read header
    call pnm_read_header()
    ! read element and organize
    call pnm_read_ele() 
    !--------------------------------------

    if (pin%s%pOrder.eq.1) then
       call pnm_apply_parmetis()
       call report_time(pin%comm)
       call pnm_setup_loc_info()
       call pnm_loc_data_build()
    elseif (pin%s%pOrder.eq.2) then
       call pnm_add_edges()
       !if (unstrM%rank.eq.0) print*,"check 0"
       call pnm_p2_parmetis()
       call report_time(pin%comm)
       call pnm_p2_data_build()
    endif 
  

  end subroutine Build_geometry

  !---------------------------------------------------------------------
  subroutine pnm_read_header()
    character (len=1024)                 :: filename
    integer                              :: i,j,k,l
    integer                              :: fid,ierr,error
    integer                              :: tmp1,tmp3 

    if (unstrM%rank == 0) then
       filename = trim(pin%fhd)
       print*,'Opening file ',trim(filename)
       inquire(file=trim(filename),exist=alive)
       if (.not. alive) then
           print*,'Error: File "',trim(filename),'" does not exist.'
           stop
       endif
       fid = 1001
       open(fid,file=trim(filename),status='old',&
                action='read',position='REWIND')
       error = 0
       read(fid,*,iostat=error) unstrM%Ntet,unstrM%Nvert
       close(fid)
       print*,'# of elements:', unstrM%Ntet
       print*,'# of vertices:', unstrM%Nvert
    endif
 
    call mpi_bcast(unstrM%Ntet, 1,mpi_integer,0,unstrM%comm,ierr)
    call mpi_bcast(unstrM%Nvert,1,mpi_integer,0,unstrM%comm,ierr)

    tmp1 = mod(unstrM%Nvert,unstrM%nproc)
    tmp3 = (unstrM%Nvert-tmp1)/unstrM%nproc
    allocate(unstrM%org%vtxdist(unstrM%nproc+1)); unstrM%org%vtxdist = 0
    do i = 1,unstrM%nproc
       if (i.le.tmp1) then
          unstrM%org%vtxdist(i+1) = unstrM%org%vtxdist(i) + tmp3 + 1
       else
          unstrM%org%vtxdist(i+1) = unstrM%org%vtxdist(i) + tmp3 
       endif
       if (i.eq.unstrM%nproc) unstrM%org%vtxdist(i+1) = unstrM%Nvert 
    enddo
    if (unstrM%rank.eq.0) then 
       print*,'initial vertex distribution'
       print*, unstrM%org%vtxdist
       !print*, 'total # of vertices', unstrM%Nvert
    endif

    unstrM%org%nvtx = unstrM%org%vtxdist(unstrM%rank+2)&
                    - unstrM%org%vtxdist(unstrM%rank+1)
     
    allocate(unstrM%org%vlist(unstrM%org%nvtx))
    do i = 1,unstrM%org%nvtx 
       unstrM%org%vlist(i) = unstrM%org%vtxdist(unstrM%rank+1) + i
    enddo

  end subroutine pnm_read_header



  !---------------------------------------------------------------------
  subroutine pnm_read_ele()
    character (len=1024)                 :: filename
    integer                              :: i,j,k,l,nints,l1,lsiz,kk,ll
    integer                              :: fid,ierr,error
    integer                              :: tmp1,tmp3 
    ! mpi I/O 
    integer                              :: stat(mpi_status_size)
    integer(kind=mpi_offset_kind)        :: offset
    integer*8                            :: leight,lfour
 
    integer, allocatable, dimension(:)   :: edist,e2v0,idtmp,vdist 
    integer, allocatable, dimension(:)   :: pid,eid,vt0,vt1,vtp,et,et0,et1
    real(kind=rkind), allocatable        :: vs0(:)
    
    integer, allocatable, dimension(:)   :: idtmp1,idt0,idat0,est0,est1
    integer, allocatable, dimension(:)   :: chk0,chk1
  
    integer                              :: sl(3),sl0(4) 
    integer, allocatable, dimension(:)   :: idv2v0,idv2v1,v2vl,v2vh,v2v0,v2v1
    integer, allocatable, dimension(:)   :: pidv,v2vd,ctv,idat1
    
    leight = 8; lfour = 4
 
    tmp1 = mod(unstrM%Ntet,unstrM%nproc)
    tmp3 = (unstrM%Ntet-tmp1)/unstrM%nproc
    allocate(edist(unstrM%nproc+1)); edist = 0 
    do i = 1,unstrM%nproc
       if (i.le.tmp1) then
          edist(i+1) = edist(i) + tmp3 + 1
       else
          edist(i+1) = edist(i) + tmp3 
       endif
       if (i.eq.unstrM%nproc) edist(i+1) = unstrM%Ntet
    enddo
    if (unstrM%rank.eq.0) then 
       print*,'initial elem. distribution'
       print*, edist
    endif
    
    ! read elements  
    filename = trim(pin%fele)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)

    ! for the length 
    nints = edist(unstrM%rank+2)-edist(unstrM%rank+1)
    nints = nints*4

    ! starting point
    offset = edist(unstrM%rank+1)*4*lfour

    allocate(e2v0(nints)); e2v0 = 0
   
    call mpi_file_read_at(fid,offset,e2v0,nints,mpi_int,stat,ierr)
    !print*,maxval(e2v0),minval(e2v0),unstrM%rank

    call mpi_file_close(fid,ierr)
    ! save it
    unstrM%org%nele = nints/4
    allocate(unstrM%org%edist(unstrM%nproc+1))
    unstrM%org%edist = edist
    allocate(unstrM%org%elist(unstrM%org%nele))
    do i = 1,unstrM%org%nele
       unstrM%org%elist(i) = unstrM%org%edist(unstrM%rank+1) + i
    enddo
    allocate(unstrM%org%e2v(nints)); unstrM%org%e2v = e2v0


 
    ! read fvst
    filename = trim(pin%fvst)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)
 
    offset = edist(unstrM%rank+1)*pin%s%pNp*leight
    allocate(vs0(unstrM%org%nele*pin%s%pNp)); vs0 = 0.0D0
    l = unstrM%org%nele*pin%s%pNp
    call mpi_file_read_at(fid,offset,vs0,l,mpi_real8,stat,ierr)
    !print*,maxval(vs0),minval(vs0),unstrM%rank
    !if (unstrM%rank.eq.0) print*,vs0
    call mpi_file_close(fid,ierr)
   

    ! construct v2e
    allocate(et(nints),et0(nints),vt0(nints))
    do i = 1,nints/4; do j = 1,4
       k = (i-1)*4+j
       et(k) = edist(unstrM%rank+1) + i 
    enddo; enddo     
    vt0 = e2v0
    call simplessort(vt0,idtmp) 
    et0 = et(idtmp)
    
    ! figure out process ids
    allocate(pid(nints),vdist(nints+1)); vdist(1) = 0
    do i = 1,nints
       vdist(i+1) = vdist(i)+1
    enddo
    call pnm_orgpid(nints,vt0,pid)
    ! change data
    call pnm_data_exchange(nints,nints,pid,vdist,vt0,l1,vtp,vt1,unstrM%comm)
    deallocate(vtp)
    !print*,'1',l1,unstrM%rank
    ! change et0
    call pnm_data_exchange(nints,nints,pid,vdist,et0,l1,vtp,et1,unstrM%comm)
    deallocate(vtp)
    !if (unstrM%rank.eq.0) print*, vt1
    !print*,'2',l1,unstrM%rank

    ! about v2e 
    unstrM%org%cnt = l1 
    allocate(unstrM%org%v2edist(unstrM%org%nvtx+1))
    unstrM%org%v2edist = 0
    allocate(unstrM%org%v2e(unstrM%org%cnt))
    call simplessort(vt1,idtmp1) 
    et1 = et1(idtmp1) 
    i = 1; j = 1
    do while (i.lt.unstrM%org%cnt) 
       if (vt1(i).eq.unstrM%org%vlist(j)) then
          unstrM%org%v2edist(j+1) = i
          i = i + 1
       else
          if (j.lt.unstrM%org%nvtx) then
             j = j + 1
          endif
       endif
    enddo
    unstrM%org%v2edist(unstrM%org%nvtx+1)=unstrM%org%cnt   
    !print*,unstrM%org%v2edist(unstrM%org%nvtx+1),unstrM%org%cnt,unstrM%rank
 
    do i = 1,unstrM%org%nvtx
       l = unstrM%org%v2edist(i)
       k = unstrM%org%v2edist(i+1)-unstrM%org%v2edist(i)
       allocate(idat0(k)) 
       idat0 = et1(l+1:l+k)
       call simplessort(idat0,idt0) 
       unstrM%org%v2e(l+1:l+k) = idat0 
       deallocate(idat0,idt0)       
    enddo
    !print*,unstrM%org%cnt,unstrM%rank 
    !if (unstrM%rank.eq.1) print*, unstrM%org%v2e  

    ! pass estat0  
    ! for elem. states
    allocate(est0(nints)); est0 = 0
    do i = 1,nints/4 
       if (maxval(vs0((i-1)*pin%s%pNp+1:i*pin%s%pNp)).lt.1.D-6) then 
          est0((i-1)*4+1:i*4) = 1
       endif
    enddo
    est0 = est0(idtmp) 
    call pnm_data_exchange(nints,nints,pid,vdist,est0,l1,vtp,est1,unstrM%comm)
    deallocate(vtp)
    est1 = est1(idtmp1)
    ! for org%vstat
    allocate(unstrM%org%vstat(unstrM%org%nvtx))
    do i = 1,unstrM%org%nvtx
       l = unstrM%org%v2edist(i)
       k = unstrM%org%v2edist(i+1)-unstrM%org%v2edist(i)
       allocate(idat0(k)) 
       idat0 = est1(l+1:l+k)
       if (maxval(idat0).eq.0) then
          ! pure solid 
          unstrM%org%vstat(i) = 0
       elseif (minval(idat0).eq.1) then
          ! pure fluid
          unstrM%org%vstat(i) = 1
       else
          ! fluid-solid points
          unstrM%org%vstat(i) = 2
       endif
       deallocate(idat0)       
    enddo
    !print*, sum(unstrM%org%vstat),unstrM%rank
 
    allocate(chk0(unstrM%nproc),chk1(unstrM%nproc)); chk0 = 0
    chk0(unstrM%rank+1) = maxval(unstrM%org%vstat(:))
    call mpi_allreduce(chk0,chk1,unstrM%nproc,mpi_integer,&
                 mpi_sum,unstrM%comm,ierr)

    if (maxval(chk1) .gt. 1) then
       unstrM%fsexist = .true.
       if (unstrM%rank.eq.0) print*, 'fluid-solid exists'
    else
       unstrM%fsexist = .false.
       if (unstrM%rank.eq.0) print*, 'NO fluid-solid'
    endif

    ! about v2v
    lsiz = nints*3
    allocate(v2vl(lsiz),v2vh(lsiz))
    do i = 1,nints/4
       sl0 = e2v0((i-1)*4+1:i*4)
       do j = 1,4
          k = (i-1)*12+(j-1)*3
          v2vh(k+1:k+3) = e2v0((i-1)*4+j)
          v2vl(k+1:k+3) = sl0(refs%FtoV(j,:)) 
       enddo
    enddo
    !print*, refs%FtoV(1,:) 
 
    call simplessort(v2vh,idv2v0) 
    v2vl = v2vl(idv2v0)

    allocate(pidv(lsiz),v2vd(lsiz+1)); v2vd = 0
    do i = 1,lsiz
       v2vd(i+1) = v2vd(i) + 1
    enddo
    call pnm_orgpid(lsiz,v2vh,pidv)
    ! change data
    call pnm_data_exchange(lsiz,lsiz,pidv,v2vd,v2vh,l1,vtp,v2v0,unstrM%comm)
    deallocate(vtp) 
    call pnm_data_exchange(lsiz,lsiz,pidv,v2vd,v2vl,l1,vtp,v2v1,unstrM%comm)
    deallocate(vtp)
    
    call simplessort(v2v0,idv2v1)
    v2v1 = v2v1(idv2v1)
    !print*,v2v1
    allocate(ctv(unstrM%org%nvtx+1)); ctv = 0
    i = 1; j = 1
    do while (i.lt.l1) 
       if (v2v0(i).eq.unstrM%org%vlist(j)) then
          ctv(j+1) = i
          i = i + 1
       else
          if (j.lt.unstrM%org%nvtx) then
             j = j + 1
          endif
       endif
    enddo
    ctv(unstrM%org%nvtx+1) = l1
 
    allocate(unstrM%org%v2vdist(unstrM%org%nvtx+1))
    allocate(unstrM%org%nv(unstrM%org%nvtx))
    unstrM%org%v2vdist = 0 
    l = 0   
    do i = 1,unstrM%org%nvtx
       l = ctv(i+1)
       k = ctv(i+1)-ctv(i)
       allocate(idat0(k)) 
       idat0 = v2v1(l-k+1:l)
       call simplessort(idat0,idt0) 
       call simplepickunique(idat0,k,idat1,j)
       !print*, k,j 
       unstrM%org%v2vdist(i+1) = unstrM%org%v2vdist(i) + j
       unstrM%org%nv(i) = j 
       deallocate(idat0,idt0,idat1)       
    enddo
    !print*, unstrM%org%v2vdist(unstrM%org%nvtx+1),unstrM%rank 
    unstrM%org%v2vsiz = unstrM%org%v2vdist(unstrM%org%nvtx+1)
    allocate(unstrM%org%v2v(unstrM%org%v2vsiz))
    do i = 1,unstrM%org%nvtx
       l = ctv(i+1)
       k = ctv(i+1)-ctv(i)
       allocate(idat0(k)) 
       idat0 = v2v1(l-k+1:l)
       call simplessort(idat0,idt0) 
       call simplepickunique(idat0,k,idat1,j)
       ll = unstrM%org%v2vdist(i+1) 
       kk = ll - unstrM%org%v2vdist(i) 
       unstrM%org%v2v(ll-kk+1:ll) = idat1
       deallocate(idat0,idt0,idat1)       
    enddo

    !print*, "check", unstrM%rank

  end subroutine pnm_read_ele 

  !--------------------------------------------------------------------
  subroutine pnm_add_edges()
    integer                                   :: i,j,k,l,ierr,l0
    integer                                   :: lvts,lvts0,lvts1,lvts2,lvts3
   
    integer, allocatable, dimension(:)        :: vts0,v2vpid0,v2vpid1,idtmp
    integer, allocatable, dimension(:)        :: dist0,dist1,ed0,ed1,tmp0,idt
    integer, allocatable, dimension(:)        :: edist,cdist,idt0

    ! establish v2vstat
    allocate(vts0(unstrM%org%v2vsiz))
    allocate(v2vpid0(unstrM%org%v2vsiz),v2vpid1(unstrM%org%v2vsiz))
    vts0 = unstrM%org%v2v
    call simplessort(vts0,idtmp)
    call pnm_orgpid(unstrM%org%v2vsiz,vts0,v2vpid0)
    v2vpid1(idtmp) = v2vpid0

    ! establish v2vstat
    allocate(unstrM%org%v2vstat(unstrM%org%v2vsiz))
    call pnm_find_info(unstrM%org%nvtx,unstrM%org%vlist,unstrM%org%vstat,&
         unstrM%org%v2vsiz,v2vpid1,unstrM%org%v2v,unstrM%org%v2vstat,unstrM%comm) 

    !print*,minval(unstrM%org%v2vstat),maxval(unstrM%org%v2vstat),unstrM%rank

    ! for edges
    allocate(ed0(unstrM%org%nvtx)); ed0 = 0
    lvts0 = 0; lvts1 = 0; lvts2 = 0
    do i = 1,unstrM%org%nvtx
       lvts = unstrM%org%v2vdist(i+1)-unstrM%org%v2vdist(i)
       !if (unstrM%org%vstat(i).gt.0) then
          do j = 1,lvts 
             k = unstrM%org%v2vdist(i) + j 
             !if (unstrM%org%v2vstat(k).gt.0) then
                lvts0 = lvts0 + 1
                ed0(i) = ed0(i) + 1
                ! criterion
                if (v2vpid1(k).lt.unstrM%rank) then 
                   lvts1 = lvts1 + 1
                elseif (v2vpid1(k).eq.unstrM%rank) then
                   lvts2 = lvts2 + 1
                endif
             !endif
          enddo
       !endif
       !unstrM%edtmp%v2vdist(i+1) = unstrM%org%v2vdist(i+1) !+ lvts0
    enddo

    allocate(dist0(unstrM%nproc),dist1(unstrM%nproc))
    dist0 = 0; dist0(unstrM%rank+1) = lvts1 + lvts2/2
 
    call mpi_allreduce(dist0,dist1,unstrM%nproc,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr)

    ! edge0 only for the vertices in this process
    edge0%nedg = lvts1 + lvts2/2 
    call mpi_allreduce(edge0%nedg,unstrM%nedge,1,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr)
    if (unstrM%rank.eq.0) print*,'num of edges', unstrM%nedge

    if (edge0%nedg.le.0) print*,edge0%nedg, unstrM%rank

    ! add together
    unstrM%edtmp%nvtx = unstrM%org%nvtx + edge0%nedg  
    allocate(unstrM%edtmp%v2vdist(unstrM%edtmp%nvtx+1))
    unstrM%edtmp%v2vdist = 0
    do i = 1,unstrM%org%nvtx
       unstrM%edtmp%v2vdist(i+1) = unstrM%org%v2vdist(i+1) !+ lvts0
    enddo

    if (edge0%nedg.gt.0) then
       j = unstrM%org%nvtx
       do i = 1,edge0%nedg
          unstrM%edtmp%v2vdist(j+i+1) = unstrM%edtmp%v2vdist(j+i) + 2
       enddo
    endif 

    allocate(edge0%eddist(unstrM%nproc+1))
    edge0%eddist = 0 
    do i = 1,unstrM%nproc
       edge0%eddist(i+1) = edge0%eddist(i)+dist1(i)
    enddo

    if (edge0%nedg.gt.0) then
       allocate(edge0%list(edge0%nedg)) 
       do i = 1,edge0%nedg
          edge0%list(i) = edge0%eddist(unstrM%rank+1) + i + unstrM%Nvert
       enddo
    endif

    unstrM%edtmp%v2vsiz = unstrM%org%v2vsiz + edge0%nedg*2
    allocate(unstrM%edtmp%v2v(unstrM%edtmp%v2vsiz))
    allocate(unstrM%edtmp%v2vstat(unstrM%edtmp%v2vsiz))
    allocate(unstrM%edtmp%v2vpid(unstrM%edtmp%v2vsiz))

    if (edge0%nedg.gt.0) then
       allocate(edge0%pid(edge0%nedg)); edge0%pid  = unstrM%rank
       allocate(edge0%ed2v(edge0%nedg,2))
       allocate(edge0%ed2vstat(edge0%nedg,2))
       allocate(edge0%ed2vpid(edge0%nedg,2))
    endif 

    ! all the edges
    edall%nedg = lvts0
    allocate(edall%list(edall%nedg),edall%pid(edall%nedg))
    allocate(edall%ed2v(edall%nedg,2))

    lvts0 = 0; lvts3 = 0
    do i = 1,unstrM%org%nvtx
       lvts = unstrM%org%v2vdist(i+1)-unstrM%org%v2vdist(i)
       do j = 1,lvts 
          k = unstrM%org%v2vdist(i) + j 
          lvts0 = lvts0 + 1
          edall%ed2v(lvts0,1) = &
              min(unstrM%org%vlist(i),unstrM%org%v2v(k)) 
          edall%ed2v(lvts0,2) = &
              max(unstrM%org%vlist(i),unstrM%org%v2v(k)) 
          ! criterion
          if (v2vpid1(k).lt.unstrM%rank) then 
             !lvts1 = lvts1 + 1
             lvts3 = lvts3 + 1 
             if (unstrM%org%vlist(i) .lt. unstrM%org%v2v(k)) then
                edge0%ed2v(lvts3,1) = unstrM%org%vlist(i) 
                edge0%ed2v(lvts3,2) = unstrM%org%v2v(k)
                edge0%ed2vstat(lvts3,1) = unstrM%org%vstat(i)
                edge0%ed2vstat(lvts3,2) = unstrM%org%v2vstat(k)
                edge0%ed2vpid(lvts3,1) = unstrM%rank
                edge0%ed2vpid(lvts3,2) = v2vpid1(k)
             else
                edge0%ed2v(lvts3,2) = unstrM%org%vlist(i) 
                edge0%ed2v(lvts3,1) = unstrM%org%v2v(k)
                edge0%ed2vstat(lvts3,2) = unstrM%org%vstat(i)
                edge0%ed2vstat(lvts3,1) = unstrM%org%v2vstat(k)
                edge0%ed2vpid(lvts3,2) = unstrM%rank
                edge0%ed2vpid(lvts3,1) = v2vpid1(k)
             endif

             edall%pid(lvts0) = unstrM%rank 
          elseif (v2vpid1(k).eq.unstrM%rank) then
             !lvts2 = lvts2 + 1
             if (unstrM%org%vlist(i) .lt. unstrM%org%v2v(k)) then
                lvts3 = lvts3 + 1
                edge0%ed2v(lvts3,1) = unstrM%org%vlist(i) 
                edge0%ed2v(lvts3,2) = unstrM%org%v2v(k)
                edge0%ed2vstat(lvts3,1) = unstrM%org%vstat(i)
                edge0%ed2vstat(lvts3,2) = unstrM%org%v2vstat(k)
                edge0%ed2vpid(lvts3,1) = unstrM%rank
                edge0%ed2vpid(lvts3,2) = v2vpid1(k)
             endif
             edall%pid(lvts0) = unstrM%rank
          else
             edall%pid(lvts0) = v2vpid1(k)
          endif
       enddo
    enddo

    ! find local edge ids
    call pnm_find_edges(edge0,edall)

    ! reorganize
    lvts0 = 0
    do i = 1,unstrM%edtmp%nvtx - edge0%nedg 
       lvts = unstrM%org%v2vdist(i+1) - unstrM%org%v2vdist(i) 
       lvts1 = unstrM%org%v2vdist(i); lvts2 = unstrM%edtmp%v2vdist(i)
       unstrM%edtmp%v2v(lvts2+1:lvts2+lvts) = unstrM%org%v2v(lvts1+1:lvts1+lvts)
       unstrM%edtmp%v2vpid(lvts2+1:lvts2+lvts)  =&
                          v2vpid1(lvts1+1:lvts1+lvts)
       unstrM%edtmp%v2vstat(lvts2+1:lvts2+lvts) =&
               unstrM%org%v2vstat(lvts1+1:lvts1+lvts)
       lvts3 = 0
       do j = 1,lvts 
          k = unstrM%org%v2vdist(i) + j 
             lvts0 = lvts0 + 1; lvts3 = lvts3 + 1
             unstrM%edtmp%v2v(lvts2+j) = edall%list(lvts0)
             unstrM%edtmp%v2vstat(lvts2+j) = 3
             if (v2vpid1(k).le.unstrM%rank) then
                unstrM%edtmp%v2vpid(lvts2+j) = unstrM%rank
             else
                unstrM%edtmp%v2vpid(lvts2+j) = v2vpid1(k)
             endif 
       enddo
       if (lvts3.gt.0) then
          allocate(tmp0(lvts)) 
          tmp0 = unstrM%edtmp%v2v(lvts2+1:lvts2+lvts)
          call simplessort(tmp0,idt)
          !print*, '1', edall%list(lvts0-lvts3+1:lvts0)
          unstrM%edtmp%v2v(lvts2+1:lvts2+lvts) = tmp0
          tmp0 = unstrM%edtmp%v2vpid(lvts2+1:lvts2+lvts)
          unstrM%edtmp%v2vpid(lvts2+1:lvts2+lvts) = tmp0(idt)
          tmp0 = unstrM%edtmp%v2vstat(lvts2+1:lvts2+lvts)
          unstrM%edtmp%v2vstat(lvts2+1:lvts2+lvts) = tmp0(idt)
          deallocate(tmp0,idt)
       endif
    enddo

    allocate(unstrM%edtmp%vlist(unstrM%edtmp%nvtx))
    unstrM%edtmp%vlist(1:unstrM%org%nvtx) = unstrM%org%vlist

    if (edge0%nedg.gt.0) then
       do i = 1,edge0%nedg
          j = unstrM%org%v2vsiz 
          unstrM%edtmp%v2v(j+(i-1)*2+1) = edge0%ed2v(i,1)
          unstrM%edtmp%v2v(j+(i-1)*2+2) = edge0%ed2v(i,2)
         
          unstrM%edtmp%v2vstat(j+(i-1)*2+1) = edge0%ed2vstat(i,1)
          unstrM%edtmp%v2vstat(j+(i-1)*2+2) = edge0%ed2vstat(i,2)

          unstrM%edtmp%v2vpid(j+(i-1)*2+1) = edge0%ed2vpid(i,1)
          unstrM%edtmp%v2vpid(j+(i-1)*2+2) = edge0%ed2vpid(i,2)
       enddo
       unstrM%edtmp%vlist(unstrM%org%nvtx+1:unstrM%edtmp%nvtx) = edge0%list
    endif

    !! establish orgn
    unstrM%orgn%nvtx = unstrM%edtmp%nvtx
    allocate(unstrM%orgn%vtxdist(unstrM%nproc+1)) 
    unstrM%orgn%vtxdist = 0 
    allocate(edist(unstrM%nproc),cdist(unstrM%nproc)) 
    edist = 0; cdist = 0
    edist(unstrM%rank+1) = unstrM%orgn%nvtx
    call mpi_allreduce(edist,cdist,unstrM%nproc,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 
    do i = 1,unstrM%nproc
       unstrM%orgn%vtxdist(i+1) = unstrM%orgn%vtxdist(i) + cdist(i) 
    enddo

    !print*, cdist 
    allocate(unstrM%orgn%vlist(unstrM%orgn%nvtx)) 
    do i = 1,unstrM%orgn%nvtx
       unstrM%orgn%vlist(i) = unstrM%orgn%vtxdist(unstrM%rank+1) + i
    enddo
    !!print*,maxval(unstrM%orgn%vlist),int(unstrM%rank,2)

    unstrM%orgn%v2vsiz = unstrM%edtmp%v2vsiz
    allocate(unstrM%orgn%v2v(unstrM%orgn%v2vsiz)) 

    call pnm_find_info(unstrM%edtmp%nvtx,unstrM%edtmp%vlist,unstrM%orgn%vlist,&
                     unstrM%edtmp%v2vsiz,unstrM%edtmp%v2vpid,unstrM%edtmp%v2v,& 
                      unstrM%orgn%v2v,unstrM%comm)  
    !print*, unstrM%orgn%v2v   
 
    allocate(unstrM%orgn%vstat(unstrM%orgn%nvtx)) 
    unstrM%orgn%vstat(1:unstrM%org%nvtx) = unstrM%org%vstat

    !todo
    !unstrM%orgn%vstat(unstrM%org%nvtx+1:unstrM%orgn%nvtx) = 3 
    do i = 1,edge0%nedg  
       if(minval(edge0%ed2vstat(i,:)).eq.0) then
         unstrM%orgn%vstat(unstrM%org%nvtx+i) = 3
       elseif (minval(edge0%ed2vstat(i,:)).eq.1) then
         unstrM%orgn%vstat(unstrM%org%nvtx+i) = 4
       elseif (minval(edge0%ed2vstat(i,:)).eq.2) then
         unstrM%orgn%vstat(unstrM%org%nvtx+i) = 5
       endif      
    enddo
    !print*, unstrM%org%nvtx,unstrM%orgn%nvtx,int(unstrM%rank,2)

    allocate(unstrM%orgn%v2vdist(unstrM%orgn%nvtx+1))
    unstrM%orgn%v2vdist = unstrM%edtmp%v2vdist

    ! find the element list
    call pnm_find_elements()

  end subroutine pnm_add_edges



  !--------------------------------------------------------------------
  subroutine pnm_apply_parmetis()
    integer(C_int64_t)                            :: parmetis_call_status 
    integer(C_int64_t)                            :: options(0:2)
    integer(C_int64_t)                            :: nparts,ncon,edgecut
    integer(C_int64_t)                            :: wgtflag,numflag
    integer(C_int64_t), allocatable, dimension(:) :: xadj,adjncy,opart
    integer(C_int64_t), allocatable, dimension(:) :: vdist,myvwgt
    integer(C_int64_t), allocatable, dimension(:) :: vtxdist,v2vdist,part
    type(C_PTR)                                   :: vwgt,adjwgt
    real(C_float), allocatable, dimension(:)      :: tpwgts,ubvec  

    integer                                       :: i,j
    integer(C_int64_t)                            :: comm

    if (unstrM%nproc.ge.2) then
       if (unstrM%fsexist) then
          allocate(myvwgt(0:3*unstrM%org%nvtx-1))
          myvwgt = 1
          do i = 0,unstrM%org%nvtx-1
             myvwgt(i*3+1) = unstrM%org%vstat(i+1)
             myvwgt(i*3+2) = unstrM%org%nv(i+1)
          enddo
          !print*, "check", unstrM%rank


          allocate(adjncy(0:unstrM%org%v2vsiz-1))
          do i = 0,unstrM%org%v2vsiz-1 
             adjncy(i) = unstrM%org%v2v(i+1) - 1
          enddo

          allocate(opart(0:unstrM%org%nvtx-1)); opart = -1
          allocate(unstrM%org%part(unstrM%org%nvtx))

          adjwgt = c_null_ptr; numflag = 0; nparts = unstrM%nproc
          options = 0; options(0) = 1; options(1) = 1
          
          wgtflag = 2;  ncon = 3
          allocate(tpwgts(ncon*nparts),ubvec(ncon))
          tpwgts = 1.0D0/real(nparts,8); ubvec = 1.0001D0

          ! add JS 10032018
          allocate(vtxdist(unstrM%nproc+1)); vtxdist = unstrM%org%vtxdist
          allocate(v2vdist(unstrM%org%nvtx+1))
          v2vdist = unstrM%org%v2vdist
          comm = unstrM%comm

          parmetis_call_status = ParMETIS_V3_PartKway(vtxdist,&
                v2vdist,adjncy,myvwgt,adjwgt,wgtflag,numflag,&
          ncon,nparts,tpwgts,ubvec,options,edgecut,opart,comm)  
  

          do i = 1,unstrM%org%nvtx
             unstrM%org%part(i) = opart(i-1)
          enddo
          if (unstrM%rank.eq.unstrM%nproc-1) print*, 'partition f-s'
       else 
          allocate(adjncy(unstrM%org%v2vsiz)); adjncy = unstrM%org%v2v - 1
          adjwgt = c_null_ptr
          allocate(myvwgt(0:2*unstrM%org%nvtx-1))
          myvwgt = 1
          do i = 0,unstrM%org%nvtx-1
             myvwgt(i*2+1) = unstrM%org%nv(i+1)
          enddo
       
          wgtflag = 2; numflag = 0; ncon = 2; nparts = unstrM%nproc 
          allocate(tpwgts(ncon*nparts),ubvec(ncon))
          tpwgts = 1.0D0/real(nparts,8); ubvec = 1.00001D0; options = 0
          options(0) = 1; options(1) = 1

          ! add JS 10032018
          allocate(vtxdist(unstrM%nproc+1)); vtxdist = unstrM%org%vtxdist
          allocate(v2vdist(unstrM%org%nvtx+1))
          v2vdist = unstrM%org%v2vdist
          comm = unstrM%comm
          allocate(part(0:unstrM%org%nvtx-1)); part=-1
 
          allocate(unstrM%org%part(unstrM%org%nvtx))
          print*, 'check', unstrM%rank
          parmetis_call_status = ParMETIS_V3_PartKway(vtxdist,&
            v2vdist,adjncy,myvwgt,adjwgt,wgtflag,numflag,ncon,&
          nparts,tpwgts,ubvec,options,edgecut,part,comm)  

          !print*, parmetis_call_status
          !print*, maxval(part),minval(part),unstrM%rank
          !print*, size(part)
          do i = 1,unstrM%org%nvtx
             unstrM%org%part(i) = part(i-1)
          enddo
          !if(unstrM%rank.eq.0) print*, opart
       endif 
       !call pnm_after_parmetis()
       call pnm_parmetis_later()
    else
       allocate(unstrM%org%part(unstrM%org%nvtx))
       unstrM%org%part = 0
       allocate(unstrM%org%v2vpid(unstrM%org%v2vsiz))
       unstrM%org%v2vpid = 0 
    endif

  end subroutine pnm_apply_parmetis


  !--------------------------------------------------------------------
  subroutine pnm_p2_parmetis()
    integer(C_int64_t)                            :: parmetis_call_status 
    integer(C_int64_t)                            :: options(0:2)
    integer(C_int64_t)                            :: nparts,ncon,edgecut
    integer(C_int64_t)                            :: wgtflag,numflag
    integer(C_int64_t), allocatable, dimension(:) :: xadj,adjncy,opart
    integer(C_int64_t), allocatable, dimension(:) :: vdist,myvwgt
    integer(C_int64_t), allocatable, dimension(:) :: vtxdist,v2vdist,part
    type(C_PTR)                                   :: vwgt,adjwgt
    real(C_float), allocatable, dimension(:)      :: tpwgts,ubvec  

    integer                                       :: i,j,nct
    integer(c_int64_t)                            :: comm

   
    if (unstrM%nproc.ge.2) then
       if (unstrM%fsexist) then
          !nct = 3
          nct = 2
          allocate(myvwgt(0:nct*unstrM%orgn%nvtx-1))
          myvwgt = 1
          do i = 0,unstrM%orgn%nvtx-1
             if (unstrM%orgn%vstat(i+1).ne.0.and.&
                 unstrM%orgn%vstat(i+1).ne.3) then
                myvwgt(i*nct) = 1
             else
                myvwgt(i*nct) = 0
             endif

             if (unstrM%orgn%vstat(i+1).ne.1.and.&
                 unstrM%orgn%vstat(i+1).ne.4) then
                myvwgt(i*nct+1) = 1
             else
                myvwgt(i*nct+1) = 0
             endif
          enddo
          allocate(adjncy(0:unstrM%orgn%v2vsiz-1))
          do i = 0,unstrM%orgn%v2vsiz-1 
             adjncy(i) = unstrM%orgn%v2v(i+1) - 1
          enddo

          !print*, unstrM%orgn%v2vsiz, unstrM%rank

          allocate(opart(0:unstrM%orgn%nvtx-1)); opart = -1
          allocate(unstrM%orgn%part(unstrM%orgn%nvtx))

          adjwgt = c_null_ptr; numflag = 0; nparts = unstrM%nproc
          options = 0; options(0) = 1; options(1) = 1
          
          wgtflag = 2;  ncon = nct !2
          allocate(tpwgts(ncon*nparts),ubvec(ncon))
          tpwgts = 1.0D0/real(nparts,8); ubvec = 1.0001D0

          ! add JS 10032018
          allocate(vtxdist(unstrM%nproc+1)); vtxdist = unstrM%orgn%vtxdist
          allocate(v2vdist(unstrM%orgn%nvtx+1))
          v2vdist = unstrM%orgn%v2vdist
          comm = unstrM%comm

          parmetis_call_status = ParMETIS_V3_PartKway(vtxdist,&
                 v2vdist,adjncy,myvwgt,adjwgt,wgtflag,numflag,&
                 ncon,nparts,tpwgts,ubvec,options,edgecut,opart,comm)  
          
          do i = 1,unstrM%orgn%nvtx
             unstrM%orgn%part(i) = opart(i-1)
          enddo
          if (unstrM%rank.eq.unstrM%nproc-1) print*, 'partition f-s'
       else
          allocate(adjncy(unstrM%orgn%v2vsiz)); adjncy = unstrM%orgn%v2v - 1
          !adjwgt = c_null_ptr
          allocate(myvwgt(0:unstrM%orgn%nvtx-1)); myvwgt = 1
          vwgt = c_null_ptr

          wgtflag = 2; numflag = 0; ncon = 1; nparts = unstrM%nproc 
          !wgtflag = 0; numflag = 0; ncon = 1; nparts = unstrM%nproc 
          allocate(tpwgts(ncon*nparts),ubvec(ncon))
          tpwgts = 1.0D0/real(nparts,8); ubvec = 1.05D0; options = 0
          options(0) = 1; options(1) = 1

          ! add JS 10032018
          allocate(vtxdist(unstrM%nproc+1)); vtxdist = unstrM%orgn%vtxdist
          allocate(v2vdist(unstrM%orgn%nvtx+1))
          v2vdist = unstrM%orgn%v2vdist
          comm = unstrM%comm
          

          allocate(unstrM%orgn%part(unstrM%orgn%nvtx))
          allocate(part(0:unstrM%orgn%nvtx-1)); part=-1

          parmetis_call_status = ParMETIS_V3_PartKway(vtxdist,&
            v2vdist,adjncy,myvwgt,adjwgt,wgtflag,numflag,ncon,&
            nparts,tpwgts,ubvec,options,edgecut,part,comm) 
          !parmetis_call_status = ParMETIS_V3now_PartKway(vtxdist,&
          !  v2vdist,adjncy,vwgt,adjwgt,wgtflag,numflag,ncon,&
          !nparts,tpwgts,ubvec,options,edgecut,part,comm) 
 
          do i = 1,unstrM%orgn%nvtx
             unstrM%orgn%part(i) = part(i-1)
          enddo
 
       endif      
 
       !call pnm_mix_parmetis_after()
       call pnm_p2_loc_info()
    else
       allocate(unstrM%orgn%part(unstrM%orgn%nvtx))
       unstrM%orgn%part = 0
       allocate(unstrM%orgn%v2vpid(unstrM%orgn%v2vsiz))
       unstrM%orgn%v2vpid = 0 
       call pnm_p2_loc_info()
    endif

  end subroutine pnm_p2_parmetis


  subroutine pnm_parmetis_later()
    integer                                   :: i,j,k,l
    integer                                   :: lvts,lvts0,lvts1,ierr
   
    integer, allocatable, dimension(:)        :: vlist,vts0,v2vpid0,v2vpid1,idtmp
   
    allocate(vlist(unstrM%org%nvtx))
    do i = 1,unstrM%org%nvtx
       vlist(i) = i + unstrM%org%vtxdist(unstrM%rank+1)
    enddo
    !print*, minval(vlist),maxval(vlist),unstrM%rank
    allocate(vts0(unstrM%org%v2vsiz))
    allocate(v2vpid0(unstrM%org%v2vsiz),v2vpid1(unstrM%org%v2vsiz))
    vts0 = unstrM%org%v2v
    call simplessort(vts0,idtmp)
    call pnm_orgpid(unstrM%org%v2vsiz,vts0,v2vpid0)
    v2vpid1(idtmp) = v2vpid0

    call pnm_orgpid(unstrM%org%v2vsiz,unstrM%org%v2v,v2vpid1)

    allocate(unstrM%org%v2vpid(unstrM%org%v2vsiz)); unstrM%org%v2vpid = 0

    call pnm_find_info(unstrM%org%nvtx,vlist,unstrM%org%part,&
         unstrM%org%v2vsiz,v2vpid1,unstrM%org%v2v,unstrM%org%v2vpid,unstrM%comm) 


  end subroutine pnm_parmetis_later

  !--------------------------------------------------------------------
  subroutine pnm_setup_loc_info()
    integer                                       :: i,j,k,l,ierr
    integer                                       :: lsiz0,lsiz1
    integer                                       :: slgth
    integer, allocatable, dimension(:)            :: vtmp1,dist0,dist1
    integer, allocatable, dimension(:)            :: sdat0,sdat1,idtmp
    integer, allocatable, dimension(:)            :: vlist1,vlist2
    integer, allocatable, dimension(:)            :: vlist3,vlist4

    !**************************************************************
    ! group vtx together
    slgth = unstrM%org%nvtx 
    allocate(sdat0(slgth),dist0(unstrM%org%nvtx+1))
    dist0(1) = 0 
    do i = 1,slgth
       sdat0(i) = i + unstrM%org%vtxdist(unstrM%rank+1)
       dist0(i+1) = dist0(i) + 1
    enddo
    !print*, maxval(unstrM%org%part),minval(unstrM%org%part),unstrM%rank
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
                      dist0,sdat0,lsiz1,vtmp1,vlist1,unstrM%comm)
    
    unstrM%new%nvtx = lsiz1; allocate(unstrM%new%vlist(lsiz1))
    unstrM%new%vlist = vlist1
    call mpi_barrier(unstrM%comm,ierr)
    print*, '# of vertices',lsiz1,' at rank', int(unstrM%rank,4)

    ! add new%vtxdist
    allocate(unstrM%new%vtxdist(unstrM%nproc+1))
    unstrM%new%vtxdist = 0
    allocate(dist1(unstrM%nproc)); dist1 = 0
    allocate(sdat1(unstrM%nproc)); sdat1 = 0
    sdat1(unstrM%rank+1) = unstrM%new%nvtx 
    call mpi_allreduce(sdat1,dist1,unstrM%nproc,mpi_integer,&
           mpi_sum,unstrM%comm,ierr) 
    do i = 1,unstrM%nproc
       unstrM%new%vtxdist(i+1) = unstrM%new%vtxdist(i)+dist1(i)
    enddo
    !print*,unstrM%new%vtxdist

    call mpi_barrier(unstrM%comm,ierr)
    deallocate(vtmp1,vlist1)
    !if (unstrm%rank.eq.0) print*, unstrM%vlist

    ! send out vstat
    !print*, maxval(unstrM%org%vtxwgt(:,2)) 
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
        dist0,unstrM%org%vstat,lsiz1,vtmp1,vlist1,unstrM%comm)
    call mpi_barrier(unstrM%comm,ierr)
    print*,'check weight balance:',sum(vlist1),int(unstrM%rank,4)  
    allocate(unstrM%new%vstat(lsiz1))    
    unstrM%new%vstat = vlist1
    deallocate(vtmp1,vlist1)
    
    ! check weight   
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
        dist0,unstrM%org%nv,lsiz1,vtmp1,vlist1,unstrM%comm)
    call mpi_barrier(unstrM%comm,ierr)
    print*,'check weight (v2v) balance:',sum(vlist1),int(unstrM%rank,4)  


    ! check v2v dist
    !call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
    !            dist0,unstrM%org%nv,lsiz1,vtmp1,vlist1,unstrM%comm)
    allocate(unstrM%new%v2vdist(lsiz1+1)); unstrM%new%v2vdist = 0
    do i = 1,lsiz1
       unstrM%new%v2vdist(i+1) = unstrM%new%v2vdist(i) + vlist1(i) 
    enddo
    deallocate(vtmp1,vlist1)
    call mpi_barrier(unstrM%comm,ierr)
    !if(unstrM%rank.eq.0) print*,unstrM%v2vdist
    !****************************************************************

    ! elm covers vtx
    slgth = unstrM%org%cnt
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
      unstrM%org%v2edist,unstrM%org%v2e,lsiz1,vtmp1,vlist1,unstrM%comm)
    !print*, lsiz1, int(unstrM%rank,2) 
    !if(unstrM%rank.eq.0) print*,vlist1
    call simplessort(vlist1,idtmp)
    deallocate(idtmp)
    call simplepickunique(vlist1,lsiz1,unstrM%Clelist,unstrM%ClNele) 
    call mpi_barrier(unstrM%comm,ierr)
    print*,'# of elements',unstrM%ClNele,' at rank', int(unstrM%rank,4)
    deallocate(vtmp1,vlist1)
    
    !*****************************************************************
    ! v2v data
    slgth = unstrM%org%v2vsiz
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
      unstrM%org%v2vdist,unstrM%org%v2v,lsiz1,vtmp1,vlist1,unstrM%comm)
    deallocate(vtmp1)

    ! v2v pid data
    slgth = unstrM%org%v2vsiz
    call pnm_data_exchange(unstrM%org%nvtx,slgth,unstrM%org%part,&
      unstrM%org%v2vdist,unstrM%org%v2vpid,lsiz1,vtmp1,vlist2,unstrM%comm)
   
    ! add vlist into v2v
    unstrM%new%v2vsiz = lsiz1 + unstrM%new%nvtx
    allocate(unstrM%new%v2vpid(unstrM%new%v2vsiz)) 
    allocate(unstrM%new%v2v(unstrM%new%v2vsiz))
    k = 0
    do i = 1,unstrM%new%nvtx 
       j = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
       allocate(vlist3(j+1),vlist4(j+1))
       vlist3(1:j) = vlist1(unstrM%new%v2vdist(i)+1:unstrM%new%v2vdist(i+1))
       vlist4(1:j) = vlist2(unstrM%new%v2vdist(i)+1:unstrM%new%v2vdist(i+1))
       vlist3(j+1) = unstrM%new%vlist(i) 
       vlist4(j+1) = unstrM%rank
       call simplessort(vlist3,idtmp)
       unstrM%new%v2v(k+1:k+j+1)    = vlist3
       unstrM%new%v2vpid(k+1:k+j+1) = vlist4(idtmp)
       k = k + j + 1
       deallocate(idtmp,vlist3,vlist4) 
    enddo
    do i = 1,unstrM%new%nvtx
       unstrM%new%v2vdist(i+1) = unstrM%new%v2vdist(i+1) + i
    enddo 
    deallocate(vtmp1,vlist1,vlist2)
    !print*, maxval(unstrM%v2vpid),minval(unstrM%v2vpid)

    !***************************************************************** 
    ! for Cvlist
    allocate(vlist1(unstrM%new%v2vsiz),vlist2(unstrM%new%v2vsiz))
    vlist1 = unstrM%new%v2v
    call simplessort(vlist1,idtmp) 
    call simplepickunique(vlist1,unstrM%new%v2vsiz,unstrM%Cvlist,unstrM%cnvtx)
    call mpi_barrier(unstrM%comm,ierr)
    print*,'total # of vertices',unstrM%cnvtx,' at rank',int(unstrM%rank,2)
    !if (unstrM%rank.eq.0) print*,unstrM%Cvlist
    do i = 1,unstrM%new%v2vsiz
       j = idtmp(i)
       vlist2(i) = unstrM%new%v2vpid(j)
    enddo
    allocate(unstrM%Cvpid(unstrM%cnvtx)); unstrM%Cvpid = -1
    i = 1; j = 1; unstrM%Cvpid(1) = vlist2(1)
    do while (i.le.unstrM%new%v2vsiz-1) 
       if (vlist1(i).eq.vlist1(i+1)) then
          i = i + 1
       else
          j = j + 1; i = i + 1
          unstrM%Cvpid(j) = vlist2(i) 
       endif
    enddo
    !if (unstrM%rank.eq.0) print*,unstrM%Cvpid 
    if (j.ne.unstrM%cnvtx) then
       print*, j,unstrM%cnvtx,&
             'Error, local vtx at rank', int(unstrM%rank,4)
       stop 
    endif

    ! check a few things
    j = 0
    do i = 1,unstrM%cnvtx 
       if (unstrM%Cvpid(i).eq.unstrM%rank) then
          j = j + 1
       endif
       !if (unstrM%Cvlist(i).eq.106) print*,unstrM%Cvpid(i),unstrM%rank
    enddo
    !print*, j, unstrM%nvtx,unstrM%cnvtx,'rank',int(unstrM%rank,2)
    if (j.ne.unstrM%new%nvtx) then
       print*, 'Error, Cvpid at rank', int(unstrM%rank,4),j,unstrM%new%nvtx
       stop
    endif

  end subroutine pnm_setup_loc_info

  subroutine pnm_p2_loc_info()
    integer                                       :: i,j,k,l,ll,ierr
    integer                                       :: lsiz0,lsiz1
    integer                                       :: slgth
    integer, allocatable, dimension(:)            :: vtmp1,dist0,dist1
    integer, allocatable, dimension(:)            :: sdat0,sdat1,idtmp,pidt
    integer, allocatable, dimension(:)            :: vlist1,vlist2,vtmp0
    integer, allocatable, dimension(:)            :: vlist3,vlist4,edist

    allocate(dist0(unstrM%orgn%nvtx+1)) 
    dist0(1) = 0
    do i = 1,unstrM%orgn%nvtx
       dist0(i+1) = dist0(i) + 1
    enddo

    call pnm_data_exchange(unstrM%orgn%nvtx,unstrM%orgn%nvtx,&
                   unstrM%orgn%part,dist0,unstrM%edtmp%vlist,&
                               lsiz1,vtmp1,vlist1,unstrM%comm)

    call simplessort(vlist1,unstrM%vord)
    allocate(unstrM%rord(lsiz1))
    do i = 1,lsiz1
       unstrM%rord(unstrM%vord(i)) = i
    enddo

    unstrM%new%nvtx = lsiz1; allocate(unstrM%new%vlist(lsiz1))
    unstrM%new%vlist = vlist1
    !if(unstrM%rank.eq.0) print*,vlist1
    !print*,vlist1
    call mpi_barrier(unstrM%comm,ierr)
    print*, '# of vertices',lsiz1,' at rank', int(unstrM%rank,4)
    call mpi_barrier(unstrM%comm,ierr)
    deallocate(vtmp1,vlist1)
      
 
    ! check weights
    call pnm_data_exchange(unstrM%orgn%nvtx,unstrM%orgn%nvtx,&
      unstrM%orgn%part,dist0,unstrM%orgn%vstat,lsiz1,vtmp1,vlist1,unstrM%comm)
    call mpi_barrier(unstrM%comm,ierr)
    j = 0; k = 0; l = 0; ll = 0
    do i = 1,lsiz1
       if (vlist1(i).eq.0.or.vlist1(i).eq.3) j = j + 1
       if (vlist1(i).eq.1.or.vlist1(i).eq.4) k = k + 1
       if (vlist1(i).eq.2.or.vlist1(i).eq.5) l = l + 1
       if (vlist1(i).ge.3) ll = ll + 1
    enddo
    !print*,'check weight balance:',sum(min(max(vlist1,0),1)),int(unstrM%rank,2)  
    print*,'check balance: fluid',int(k+l,4),'solid',int(j+l,4),int(unstrM%rank,4)  
    allocate(unstrM%new%vstat(lsiz1)); unstrM%new%vstat = vlist1(unstrM%vord)
    call mpi_barrier(unstrM%comm,ierr)
    deallocate(vtmp1,vlist1)

    ! set up edge info
    edge1%nedg = ll
    allocate(edge1%ed2v(edge1%nedg,2),edge1%list(edge1%nedg))
    j = 0 
    do i = 1,lsiz1
       !if (unstrM%new%vlist(i).le.unstrM%Nvert) then
       !   print*, unstrM%new%vlist(i),unstrM%new%vstat(i)
       !endif
       if (unstrM%new%vstat(i).ge.3) then 
          j = j + 1
          edge1%list(j) = unstrM%new%vlist(i)
       endif
    enddo

    allocate(pidt(edge0%nedg),dist1(edge0%nedg+1)); dist1(1) = 0
    allocate(vtmp0(edge0%nedg))
    if (edge0%nedg.gt.0) then
       do i = 1,edge0%nedg
          j = unstrM%org%nvtx + i
          pidt(i) = unstrM%orgn%part(j)
          dist1(i+1) = dist1(i) + 1
       enddo
       vtmp0 = edge0%ed2v(:,1)
    endif

    call pnm_data_exchange(edge0%nedg,edge0%nedg,pidt,dist1,vtmp0,&
                           lsiz1,vtmp1,vlist1,unstrM%comm)
    edge1%ed2v(:,1) = vlist1(:)
    deallocate(vtmp1,vlist1)

    !print*, edge1%nedg,lsiz1
    if (edge0%nedg.gt.0) vtmp0 = edge0%ed2v(:,2)
    call pnm_data_exchange(edge0%nedg,edge0%nedg,pidt,dist1,vtmp0,&
                           lsiz1,vtmp1,vlist1,unstrM%comm)
    edge1%ed2v(:,2) = vlist1(:)
    deallocate(vtmp1,vlist1)

    ! element info
    ! send v2edist
    allocate(sdat0(unstrM%orgn%nvtx))
    do i = 1,unstrM%orgn%nvtx
       sdat0(i) = unstrM%orgn%v2edist(i+1) - unstrM%orgn%v2edist(i) 
    enddo
    call pnm_data_exchange(unstrM%orgn%nvtx,unstrM%orgn%nvtx,unstrM%orgn%part,&
                           dist0,sdat0,lsiz1,vtmp1,vlist1,unstrM%comm)
    allocate(unstrM%new%v2edist(unstrM%new%nvtx+1))
    allocate(edist(unstrM%new%nvtx+1))
    unstrM%new%v2edist = 0; edist = 0
    do i = 1,unstrM%new%nvtx
       edist(i+1) = edist(i) + vlist1(i)
       unstrM%new%v2edist(i+1) = unstrM%new%v2edist(i) + vlist1(unstrM%vord(i))
    enddo
    !print*,i,lsiz1,unstrM%new%nvtx,unstrM%new%v2edist(i),int(unstrM%rank,2) 
    !print*,unstrM%new%v2edist 
    deallocate(vtmp1,vlist1)

    call pnm_data_exchange(unstrM%orgn%nvtx,unstrM%orgn%cnt,unstrM%orgn%part,&
      unstrM%orgn%v2edist,unstrM%orgn%v2e,lsiz1,vtmp1,vlist1,unstrM%comm)
    ! TODO vlist1 important
    unstrM%new%cnt = lsiz1
    allocate(unstrM%new%v2e(unstrM%new%cnt))
    do i = 1,unstrM%new%nvtx
       k = edist(unstrM%vord(i)) 
       j = unstrM%new%v2edist(i+1)-unstrM%new%v2edist(i)
       !print*, j
       l = unstrM%new%v2edist(i)
       unstrM%new%v2e(l+1:l+j) = vlist1(k+1:k+j) 
    enddo 
    !unstrM%new%v2e = vlist1

    call simplessort(vlist1,idtmp); deallocate(idtmp)
    call simplepickunique(vlist1,lsiz1,unstrM%Clelist,unstrM%ClNele) 
    call mpi_barrier(unstrM%comm,ierr)
    print*,'# of elements covered',unstrM%ClNele,' at rank', int(unstrM%rank,4)
    deallocate(vtmp1,vlist1)


  end subroutine pnm_p2_loc_info




  !--------------------------------------------------------------------
  subroutine pnm_loc_data_build()
    integer                              :: i,j,k,l,l0,l1,lc
    integer, allocatable, dimension(:)   :: pid0,idt0,idt1

    integer, allocatable, dimension(:)   :: tmpd0,tmpd1,idtmp

    character (len=1024)                 :: filename
    integer                              :: fid,ierr,error,nints
    ! mpi I/O 
    integer                              :: stat(mpi_status_size)
    integer(kind=mpi_offset_kind)        :: offset
    real(kind=rkind), allocatable        :: buf(:)
    integer*8                            :: leight,lfour

    integer, allocatable, dimension(:)   :: pidv(:),e2e0(:)
    real(kind=rkind), allocatable        :: rd0(:),rd1(:) 

    integer                              :: eleverts(4),location 
    real(kind=rkind)                     :: B(3,3),invB(3,3),detB
    real(kind=rkind)                     :: normal(3,4),pts(3,4),Jacs(4)

    leight = 8; lfour = 4

    ! read elements  
    filename = trim(pin%fneigh)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)
 
    nints = unstrM%org%edist(unstrM%rank+2)-unstrM%org%edist(unstrM%rank+1)
    nints = nints*4

    offset = unstrM%org%edist(unstrM%rank+1)*4*lfour
    allocate(e2e0(nints)); e2e0 = 0
    call mpi_file_read_at(fid,offset,e2e0,nints,mpi_int,stat,ierr)
    !print*,maxval(e2v0),minval(e2v0),unstrM%rank
    call mpi_file_close(fid,ierr)

    if (unstrM%rank.eq.0) print*,'read in neigh info'     

    ! for loc_t2v and ClNeigh
    allocate(pid0(unstrM%ClNele))
    call pnm_org_epid(unstrM%ClNele,unstrM%Clelist,pid0)  
    allocate(unstrM%CNOids(unstrM%ClNele))
    unstrM%CNOids = pid0 

 
    allocate(idt0(unstrM%org%nele),idt1(unstrM%ClNele))
    allocate(unstrM%loc_t2v(4,unstrM%ClNele))
    allocate(unstrM%ClNeigh(4,unstrM%ClNele)) 
    do j = 1,4
       do i = 1,unstrM%org%nele
          idt0(i) = unstrM%org%e2v((i-1)*4+j)
       enddo    

       call pnm_find_info(unstrM%org%nele,unstrM%org%elist,idt0,&
             unstrM%ClNele,pid0,unstrM%Clelist,idt1,unstrM%comm)
        
       unstrM%loc_t2v(j,:) = idt1

       do i = 1,unstrM%org%nele
          idt0(i) = e2e0((i-1)*4+j)
       enddo    

       call pnm_find_info(unstrM%org%nele,unstrM%org%elist,idt0,&
             unstrM%ClNele,pid0,unstrM%Clelist,idt1,unstrM%comm)

       unstrM%ClNeigh(j,:) = idt1
    enddo
    !print*, unstrM%loc_t2v(:,1)

    !if (pin%mix.and.unstrM%fsexist) then 
    !   l0 = unstrM%ClNele*4
    !   allocate(tmpd0(l0)) 
    !   do i = 1,unstrM%ClNele; do j = 1,4
    !      k = (i-1)*4+j
    !      tmpd0(k) = unstrM%loc_t2v(j,i)
    !   enddo; enddo 
    !   call simplessort(tmpd0,idtmp)
    !   call simplepickunique(tmpd0,l0,tmpd1,l1)
    !   unstrM%cnvtx = l1 
    !   allocate(unstrM%Cvlist(unstrM%cnvtx))
    !   unstrM%Cvlist = tmpd1 
    !   !print*, l1, int(unstrM%rank,2)
    !   deallocate(tmpd0,tmpd1,idtmp) 
    !endif 

    ! read vertex info
    filename = trim(pin%fnode)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)

    nints = unstrM%org%nvtx*3 

    offset = unstrM%org%vtxdist(unstrM%rank+1)*3*leight
    !print*, offset,nints,unstrM%rank

    allocate(buf(nints)); buf = 0.0D0
     
    call mpi_file_read_at(fid,offset,buf,nints,mpi_real8,stat,ierr)
    !call mpi_file_seek(fid,offset,mpi_seek_set,ierr)
    !call mpi_file_read(fid,buf,nints,mpi_real8,stat,ierr)
    !print*, maxval(buf),minval(buf),unstrM%rank,ierr
    !print*, stat, unstrM%rank

    call mpi_file_close(fid,ierr)

    ! for Cv_crs(3,cnvtx)
    allocate(unstrM%Cv_crs(3,unstrM%cnvtx))
    allocate(pidv(unstrM%cnvtx))
    call pnm_orgpid(unstrM%cnvtx,unstrM%Cvlist,pidv)  
    allocate(rd0(unstrM%org%nvtx),rd1(unstrM%cnvtx))

    do j = 1,3
       do i = 1,unstrM%org%nvtx
          rd0(i) = buf((i-1)*3+j)
       enddo   

       call pnm_find_data(unstrM%org%nvtx,unstrM%org%vlist,rd0,&
             unstrM%cnvtx,pidv,unstrM%Cvlist,rd1,unstrM%comm)

       unstrM%Cv_crs(j,:) = rd1
    enddo

    allocate(unstrM%lt2vid(4*unstrM%ClNele))
    do i = 1,unstrM%ClNele
       do j = 1,4
          k = (i-1)*4+j
          l = unstrM%loc_t2v(j,i)
          call findorder(l,unstrM%Cvlist,location)
          unstrM%lt2vid(k) = location
       enddo
    enddo

    allocate(unstrM%loc_vtx(3,4*unstrM%ClNele))
   
    allocate(unstrM%loc_nods(3,4*unstrM%ClNele))

    do i = 1,unstrM%ClNele
       do j = 1,4
          k = (i-1)*4+j
          l = unstrM%lt2vid(k) 
          unstrM%loc_vtx(:,k) = unstrM%Cv_crs(:,l)
          !l = unstrM%s%loc_label(k)
          unstrM%loc_nods(:,k) = unstrM%Cv_crs(:,l)
       enddo
    enddo 

    ! for local Jacobian
    if (pin%buildtetJac) then 
       ! -------------------
       ! build tets
       call writelog("..Building local Jacobians",unstrM%comm)
       if (unstrM%rank == 0) print*, "Building local Jacobians" 
       allocate(unstrM%loc_Jac(3,3,unstrM%ClNele))
       allocate(unstrM%loc_invJ(3,3,unstrM%ClNele))
       allocate(unstrM%loc_detJ(unstrM%ClNele))
       do i = 1, unstrM%ClNele
          j = unstrM%Clelist(i)
          call transfMat(i,B,invB,detB)
          unstrM%loc_Jac(:,:,i)  = B/2.0D0
          unstrM%loc_invJ(:,:,i) = invB*2.0D0
          unstrM%loc_detJ(i)     = detB/8.0D0!/1.D9 
       enddo
         
    endif
    
    ! for gravitational acceleration parts
    if (pin%buildtetn) then
       call writelog("..Building local normal vectors",unstrM%comm)
       if (unstrM%rank == 0) print*, "Build local normal vectors"
       allocate(   unstrM%loc_n(3,4,unstrM%ClNele))
       allocate(  unstrM%loc_sJac(4,unstrM%ClNele))
       !allocate(unstrM%loc_Fscale(4,unstrM%ClNele))
       do i = 1, unstrM%ClNele
          eleverts(:) = unstrM%lt2vid((i-1)*4+1:i*4)
          do j = 1,4
             pts(:,j) = unstrM%Cv_crs(:,eleverts(j))
          enddo
          call face_normal(pts,normal,Jacs)
          
          unstrM%loc_n(:,:,i)    = normal
          unstrM%loc_sJac(:,i)   = Jacs/4.0D0
          !unstrM%loc_Fscale(:,i) = Jacs/unstrM%loc_detJ(i)
       enddo
       !print*, refs%Fmask

    endif

    if (pin%phi1) then 
       ! prepare for S(u) 
       call pnm_find_lNele()
       call pnm_find_lsurf()
    endif 
 
  end subroutine pnm_loc_data_build

  !--------------------------------------------------------------------
  subroutine pnm_p2_data_build()
    integer                              :: i,j,k,l,l0,l1
    integer, allocatable, dimension(:)   :: pid0,idt0,idt1,Cvlist

    integer, allocatable, dimension(:)   :: tmpd0,tmpd1,idtmp

    character (len=1024)                 :: filename
    integer                              :: fid,ierr,error,nints
    ! mpi I/O 
    integer                              :: stat(mpi_status_size)
    integer(kind=mpi_offset_kind)        :: offset
    real(kind=rkind), allocatable        :: buf(:)

    integer, allocatable, dimension(:)   :: pidv(:),e2e0(:),e2v0(:,:)
    real(kind=rkind), allocatable        :: rd0(:),rd1(:) 
 
    integer*8                            :: leight,lfour

    integer                              :: eleverts(4),location 
    real(kind=rkind)                     :: B(3,3),invB(3,3),detB
    real(kind=rkind)                     :: normal(3,4),pts(3,4),Jacs(4)
    real(kind=rkind), allocatable        :: tetnodes(:,:)

    ! for edges
    integer                              :: k1,k2,ll,lc,ord(10),org(10)
    integer                              :: lsiz0,lsiz1,vtx4(4),lnew
    integer, allocatable, dimension(:)   :: tmp0,tmp1,v2v,v2vpid
    integer, allocatable, dimension(:)   :: tmpv0,tmpid,oldid,vid1 
    integer, allocatable, dimension(:)   :: cdist,edist 
    integer, allocatable, dimension(:,:) :: ClNids

    leight = 8; lfour = 4

    ! read elements  
    filename = trim(pin%fneigh)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)
 
    nints = unstrM%org%edist(unstrM%rank+2)-unstrM%org%edist(unstrM%rank+1)
    nints = nints*4

    offset = unstrM%org%edist(unstrM%rank+1)*4*lfour
    allocate(e2e0(nints)); e2e0 = 0
    call mpi_file_read_at(fid,offset,e2e0,nints,mpi_int,stat,ierr)
    !print*,maxval(e2v0),minval(e2v0),unstrM%rank
    call mpi_file_close(fid,ierr)

    ! for loc_t2v and ClNeigh
    allocate(pid0(unstrM%ClNele))
    call pnm_org_epid(unstrM%ClNele,unstrM%Clelist,pid0)  
    allocate(unstrM%CNOids(unstrM%ClNele))
    unstrM%CNOids = pid0 

 
    allocate(idt0(unstrM%org%nele),idt1(unstrM%ClNele))
    allocate(e2v0(4,unstrM%ClNele))
    allocate(unstrM%ClNeigh(4,unstrM%ClNele)) 
    do j = 1,4
       do i = 1,unstrM%org%nele
          idt0(i) = unstrM%org%e2v((i-1)*4+j)
       enddo    

       call pnm_find_info(unstrM%org%nele,unstrM%org%elist,idt0,&
             unstrM%ClNele,pid0,unstrM%Clelist,idt1,unstrM%comm)
       
       e2v0(j,:) = idt1
       do i = 1,unstrM%org%nele
          idt0(i) = e2e0((i-1)*4+j)
       enddo    

       call pnm_find_info(unstrM%org%nele,unstrM%org%elist,idt0,&
             unstrM%ClNele,pid0,unstrM%Clelist,idt1,unstrM%comm)

       unstrM%ClNeigh(j,:) = idt1
    enddo
    ! construct Cvlist
    l0 = unstrM%ClNele*4
    allocate(tmpd0(l0))
    do i = 1,unstrM%ClNele; do j = 1,4
       k = (i-1)*4+j
       tmpd0(k) = e2v0(j,i)
    enddo; enddo 
    call simplessort(tmpd0,idtmp)
    call simplepickunique(tmpd0,l0,tmpd1,l1)
    allocate(Cvlist(l1))
    Cvlist = tmpd1 
    deallocate(tmpd0,tmpd1,idtmp) 


    ! read vertex info
    filename = trim(pin%fnode)
    call mpi_file_open(unstrM%comm,filename,mpi_mode_rdonly,&
         mpi_info_null,fid,ierr)
    nints = unstrM%org%nvtx*3 
    offset = unstrM%org%vtxdist(unstrM%rank+1)*3*leight
    !print*, offset,nints,unstrM%rank
    allocate(buf(nints)); buf = 0.0D0
    call mpi_file_read_at(fid,offset,buf,nints,mpi_real8,stat,ierr)
    call mpi_file_close(fid,ierr)

    ! l1 size of all vertices
    allocate(unstrM%Cv_crs(3,l1)); allocate(pidv(l1))
    call pnm_orgpid(l1,Cvlist,pidv)  
    allocate(rd0(unstrM%org%nvtx),rd1(l1))
    do j = 1,3
       do i = 1,unstrM%org%nvtx
          rd0(i) = buf((i-1)*3+j)
       enddo   

       call pnm_find_data(unstrM%org%nvtx,unstrM%org%vlist,rd0,&
                   l1,pidv,Cvlist,rd1,unstrM%comm)

       unstrM%Cv_crs(j,:) = rd1
    enddo
   
    allocate(unstrM%loc_vtx(3,4*unstrM%ClNele))
    do i = 1,unstrM%ClNele
       do j = 1,4
          k = (i-1)*4+j
          call findorder(e2v0(j,i),Cvlist,l) 
          unstrM%loc_vtx(:,k) = unstrM%Cv_crs(:,l)
       enddo
    enddo 


    ! find out Cvlist
    ednew%nedg = unstrM%ClNele*6
    allocate(ednew%ed2v(ednew%nedg,2)) 
    allocate(ednew%list(ednew%nedg),ednew%pid(ednew%nedg))

    do i = 1,unstrM%ClNele
       vtx4 = e2v0(:,i)
       !call simplessort(vtx4,idtmp); deallocate(idtmp)
       ednew%ed2v((i-1)*6+1,1) = min(vtx4(1),vtx4(2)) 
       ednew%ed2v((i-1)*6+1,2) = max(vtx4(1),vtx4(2)) 
       ednew%ed2v((i-1)*6+2,1) = min(vtx4(1),vtx4(3)) 
       ednew%ed2v((i-1)*6+2,2) = max(vtx4(1),vtx4(3)) 
       ednew%ed2v((i-1)*6+3,1) = min(vtx4(2),vtx4(3)) 
       ednew%ed2v((i-1)*6+3,2) = max(vtx4(2),vtx4(3))
       ednew%ed2v((i-1)*6+4,1) = min(vtx4(1),vtx4(4)) 
       ednew%ed2v((i-1)*6+4,2) = max(vtx4(1),vtx4(4))
       ednew%ed2v((i-1)*6+5,1) = min(vtx4(2),vtx4(4)) 
       ednew%ed2v((i-1)*6+5,2) = max(vtx4(2),vtx4(4))
       ednew%ed2v((i-1)*6+6,1) = min(vtx4(3),vtx4(4))
       ednew%ed2v((i-1)*6+6,2) = max(vtx4(3),vtx4(4))
    enddo

    call pnm_orgpid(ednew%nedg,ednew%ed2v(:,2),ednew%pid)
   
    call pnm_find_edges(edge0,ednew)
    allocate(oldid(ednew%nedg)); oldid = ednew%pid
    call pnm_find_info(unstrM%orgn%nvtx,unstrM%edtmp%vlist,unstrM%orgn%part,&
                         ednew%nedg,oldid,ednew%list,ednew%pid,unstrM%comm)
 
    ! check process ids
    allocate(tmpv0(unstrM%ClNele*4),tmpid(unstrM%ClNele*4))
    do i = 1,unstrM%ClNele; do j = 1,4 
       k = (i-1)*4+j
       tmpv0(k) = e2v0(j,i) 
    enddo; enddo
    lsiz0 = unstrM%ClNele*4
    call pnm_orgpid(lsiz0,tmpv0,tmpid)
    allocate(vid1(unstrM%ClNele*4))
    call pnm_find_info(unstrM%orgn%nvtx,unstrM%edtmp%vlist,unstrM%orgn%part,&
                         lsiz0,tmpid,tmpv0,vid1,unstrM%comm) 
    allocate(ClNids(4,unstrM%ClNele)) 
    do i = 1,unstrM%ClNele; do j = 1,4 
       k = (i-1)*4+j
       ClNids(j,i) = vid1(k)
    enddo; enddo

    ! build new%v2v
    allocate(unstrM%new%vtxdist(unstrM%nproc+1))
    unstrM%new%vtxdist(1) = 0
    allocate(edist(unstrM%nproc),cdist(unstrM%nproc)) 
    edist = 0; cdist = 0
    edist(unstrM%rank+1) = unstrM%new%nvtx
    call mpi_allreduce(edist,cdist,unstrM%nproc,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 
    do i = 1,unstrM%nproc
       unstrM%new%vtxdist(i+1) = unstrM%new%vtxdist(i) + cdist(i) 
    enddo
    !if(unstrM%rank.eq.0) print*,unstrM%new%vtxdist  
 
    allocate(unstrM%new%nv(unstrM%new%nvtx))    
    do i = 1,unstrM%new%nvtx
       lc = unstrM%new%vlist(i)
       k = unstrM%new%v2edist(i)
       lsiz0 = unstrM%new%v2edist(i+1) - unstrM%new%v2edist(i)
       k1 = lsiz0*10
       !print*,k1
       allocate(tmp0(k1)); tmp0 = 0
       k1 = 0
       do j = 1,lsiz0
          l = unstrM%new%v2e(k+j)
          call findorder(l,unstrM%Clelist,k2) 
          k1 = k1 + 1
          tmp0((k1-1)*10+1:(k1-1)*10+4)  = e2v0(1:4,k2) 
          tmp0((k1-1)*10+5:k1*10) = ednew%list((k2-1)*6+1:k2*6)
       enddo
       lsiz0 = lsiz0*10
       call pnm_setup_lv2v(lc,lsiz0,10,tmp0,lsiz1)
       unstrM%new%nv(i) = lsiz1
       deallocate(tmp0)  
    enddo

    allocate(unstrM%new%v2vdist(unstrM%new%nvtx+1)) 
    unstrM%new%v2vdist(1) = 0
    do i = 1,unstrM%new%nvtx
       unstrM%new%v2vdist(i+1) = unstrM%new%v2vdist(i) + unstrM%new%nv(i) 
    enddo

    lsiz0 = unstrM%new%v2vdist(unstrM%new%nvtx+1)
    unstrM%new%v2vsiz = lsiz0
    allocate(unstrM%new%v2v(lsiz0))
    allocate(unstrM%new%v2vpid(lsiz0))
    allocate(unstrM%new%v2vstat(lsiz0))

    do i = 1,unstrM%new%nvtx
       lc = unstrM%new%vlist(i)
       k = unstrM%new%v2edist(i)
       lsiz0 = unstrM%new%v2edist(i+1) - unstrM%new%v2edist(i)
       k1 = lsiz0*10
       !print*,k1
       allocate(tmp0(k1)); tmp0 = 0
       allocate(tmp1(k1)); tmp1 = 0
       k1 = 0
       do j = 1,lsiz0
          l = unstrM%new%v2e(k+j)
          call findorder(l,unstrM%Clelist,k2) 
          k1 = k1 + 1
          tmp0((k1-1)*10+1:(k1-1)*10+4)  = e2v0(1:4,k2) 
          tmp0((k1-1)*10+5:k1*10) = ednew%list((k2-1)*6+1:k2*6)
          tmp1((k1-1)*10+1:(k1-1)*10+4)  = ClNids(1:4,k2) 
          tmp1((k1-1)*10+5:k1*10) = ednew%pid((k2-1)*6+1:k2*6)
       enddo
       lsiz0 = lsiz0*10
       call pnm_find_lv2v(lsiz0,tmp0,tmp1,lnew,v2v,v2vpid)

       lc = unstrM%new%v2vdist(i) 
       unstrM%new%v2v(lc+1:lc+lnew)    = v2v
       unstrM%new%v2vpid(lc+1:lc+lnew) = v2vpid
       deallocate(tmp0,tmp1,v2v,v2vpid)  
    enddo

    call pnm_find_info(unstrM%new%nvtx,unstrM%new%vlist,&
           unstrM%new%vstat,unstrM%new%v2vsiz,unstrM%new%v2vpid,&
           unstrM%new%v2v,unstrM%new%v2vstat,unstrM%comm) 

    call mpi_allreduce(unstrM%new%v2vsiz,lsiz0,1,&
                   mpi_integer,mpi_sum,unstrM%comm,ierr)

    ! for Cvlist
    call pnm_constct_cvpid(unstrM%new%v2vsiz,unstrM%new%v2v,unstrM%new%v2vpid,&
                    unstrM%new%cnvtx,unstrM%new%Cvlist,unstrM%new%Cvpid)

    unstrM%cnvtx = unstrM%new%cnvtx
    allocate(unstrM%Cvlist(unstrM%cnvtx))
    unstrM%Cvlist = unstrM%new%Cvlist
    allocate(unstrM%Cvpid(unstrM%cnvtx))
    unstrM%Cvpid = unstrM%new%Cvpid

    ! construct loc_t2v and lt2vid (todo)
    allocate(unstrM%loc_t2v(pin%s%pNp,unstrM%ClNele))
    allocate(unstrM%lt2vid(pin%s%pNp*unstrM%ClNele))
    ord = (/1,3,6,10,2,4,5,7,8,9/)
    do i = 1,unstrM%ClNele
       org(1:4)  = e2v0(:,i)
       org(5:10) = ednew%list((i-1)*6+1:i*6)
       unstrM%loc_t2v(ord,i) = org
    enddo 

    do i = 1,unstrM%ClNele; do j = 1,pin%s%pNp
       call findorder(unstrM%loc_t2v(j,i),unstrM%new%Cvlist,l)
       unstrM%lt2vid((i-1)*pin%s%pNp+j) = l
    enddo; enddo


    ! for local bases nodes
    if (pin%buildBasnodeOn) then 
       call writelog("..Building Base nodes",unstrM%comm)
       if (unstrM%rank == 0) print*, "Building base nodes"
       allocate(unstrM%loc_nods(3,pin%s%pNp*unstrM%ClNele))
       allocate(tetnodes(3,pin%s%pNp))  
       do i = 1,unstrM%ClNele
          !j = unstrM%Clelist(i)
          call node_coo(tetnodes,i,pin%s%pNp) 
          unstrM%loc_nods(:,(i-1)*pin%s%pNp+1:i*pin%s%pNp) = tetnodes
       enddo
       !print*, minval(unstrM%loc_nods)
       !print*,unstrM%loc_nods(1,1:20)
    endif 

    ! for local Jacobian
    if (pin%buildtetJac) then 
       ! -------------------
       ! build tets
       call writelog("..Building local Jacobians",unstrM%comm)
       if (unstrM%rank == 0) print*, "Building local Jacobians" 
       allocate(unstrM%loc_Jac(3,3,unstrM%ClNele))
       allocate(unstrM%loc_invJ(3,3,unstrM%ClNele))
       allocate(unstrM%loc_detJ(unstrM%ClNele))
       do i = 1, unstrM%ClNele
          !j = unstrM%Clelist(i)
          call transfMat(i,B,invB,detB)
          unstrM%loc_Jac(:,:,i)  = B/2.0D0
          unstrM%loc_invJ(:,:,i) = invB*2.0D0
          unstrM%loc_detJ(i)     = detB/8.0D0!/1.D9 
       enddo
         
    endif
 
    ! for gravitational acceleration parts
    if (pin%buildtetn) then
       call writelog("..Building local normal vectors",unstrM%comm)
       if (unstrM%rank == 0) print*, "Build local normal vectors"
       allocate(   unstrM%loc_n(3,4,unstrM%ClNele))
       allocate(  unstrM%loc_sJac(4,unstrM%ClNele))
       do i = 1, unstrM%ClNele
          do j = 1,4
             k = (i-1)*4+j
             call findorder(e2v0(j,i),Cvlist,l) 
             !pts(:,j) = unstrM%loc_vtx(:,(i-1)*4+j)
             pts(:,j) = unstrM%Cv_crs(:,l)
          enddo
          call face_normal(pts,normal,Jacs)
          
          unstrM%loc_n(:,:,i)    = normal
          unstrM%loc_sJac(:,i)   = Jacs/4.0D0
       enddo
    endif

    if (pin%phi1) then 
       ! prepare for S(u) 
       call pnm_find_lNele()
       call pnm_find_lsurf()
    endif 

  end subroutine pnm_p2_data_build


  !--------------------------------------------------------------------
  subroutine pnm_find_elements()
    integer                                   :: i,j,k,l,ierr,l0,l1,lc
    integer                                   :: lvts0,lvts1,lct0,lct1
    integer, allocatable, dimension(:)        :: idtmp,dat0,dat1,npid0,npid1
    integer, allocatable, dimension(:)        :: cnt0,cnt1,tmpid,tmpd0,tmpd1
    integer, allocatable, dimension(:)        :: vts0,vts1,ad0,ad1
    integer, allocatable, dimension(:)        :: np0,np1,ct0,ct1,ciod0
    integer, allocatable, dimension(:)        :: tmp0,tmp1,cnum0,dist,com,cout
 

    allocate(npid0(unstrM%nproc),npid1(unstrM%nproc))
    npid0 = 0; npid1 = 0
 
    if (edge0%nedg.gt.0) then
       allocate(tmpid(edge0%nedg*2),tmpd0(edge0%nedg*2))
       do i = 1,edge0%nedg; do j = 1,2
          tmpid((i-1)*2+j) = edge0%ed2vpid(i,j)
          tmpd0((i-1)*2+j) = edge0%ed2v(i,j)
       enddo; enddo
       call simplessort(tmpid,idtmp)

       allocate(ciod0(2*edge0%nedg))
       do i = 1,edge0%nedg; do j = 1,2
          k = (i-1)*2+j; ciod0(idtmp(k)) = k
       enddo; enddo
       !print*, maxval(tmpd0-edall%pid(idtmp))
       allocate(vts0(2*edge0%nedg))
       vts0 = tmpd0(idtmp) 

       do i = 1,edge0%nedg*2
          j = tmpid(i) + 1; npid0(j) = npid0(j) + 1         
       enddo      
    endif
   
    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,unstrM%comm,ierr) 

    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    lvts0 = sum(npid0); lvts1 = sum(npid1)

    allocate(vts1(lvts1))
    if (.not.allocated(vts0)) allocate(vts0(lvts0)) 

    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,unstrM%comm,ierr)

    allocate(dat0(lvts0),dat1(lvts1))
    !print*, lvts1

    if (lvts1.gt.0) then
       do i = 1,lvts1
          j = vts1(i) - unstrM%org%vtxdist(unstrM%rank+1) 
          dat1(i) = unstrM%org%v2edist(j+1)-unstrM%org%v2edist(j)      
       enddo
    endif

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_integer,&
                       dat0,npid0,cnt0,mpi_integer,unstrM%comm,ierr)

    l0 = sum(dat0); l1 = sum(dat1) 
  
    allocate(ad0(l0),ad1(l1)) 
    allocate(np0(unstrM%nproc),np1(unstrM%nproc))
    allocate(ct0(0:unstrM%nproc-1),ct1(0:unstrM%nproc-1))

    do i = 1,unstrM%nproc-1
       np0(i) = sum(dat0(cnt0(i-1)+1:cnt0(i))) 
    enddo
    if (cnt0(unstrM%nproc-1).lt. lvts0) then
       np0(unstrM%nproc) = sum(dat0(cnt0(unstrM%nproc-1)+1:lvts0))
    else
       np0(unstrM%nproc) = 0
    endif
       
    do i = 1,unstrM%nproc-1
       np1(i) = sum(dat1(cnt1(i-1)+1:cnt1(i))) 
    enddo
    if (cnt1(unstrM%nproc-1).lt. lvts1) then
       np1(unstrM%nproc) = sum(dat1(cnt1(unstrM%nproc-1)+1:lvts1))
    else
       np1(unstrM%nproc) = 0
    endif
    
    ct0 = 0; ct1 = 0
    do i = 1,unstrM%nproc-1
       ct0(i) = ct0(i-1) + np0(i)
       ct1(i) = ct1(i-1) + np1(i)
    enddo

    l = 0
    do i = 1,lvts1
       j = vts1(i) - unstrM%org%vtxdist(unstrM%rank+1)
       k = unstrM%org%v2edist(j)      
       lct1 = unstrM%org%v2edist(j+1)-unstrM%org%v2edist(j)     
       ad1(l+1:l+lct1) = unstrM%org%v2e(k+1:k+lct1)
       l = l + lct1 
    enddo  

    call MPI_ALLTOALLV(ad1,np1,ct1,mpi_integer,&
                       ad0,np0,ct0,mpi_integer,unstrM%comm,ierr)

    !print*, ad0

    if (edge0%nedg.gt.0) then
       allocate(dist(2*edge0%nedg+1)); dist = 0  
       do i = 1,edge0%nedg; do j = 1,2
          k = (i-1)*2+j
          dist(k+1) = dist(k) + dat0(k)
       enddo; enddo

       allocate(cnum0(edge0%nedg))
       do i = 1,edge0%nedg
          lct0 = ciod0((i-1)*2+1); lct1 = ciod0((i-1)*2+2)
          l0 = dat0(lct0); l1  = dat0(lct1)
          !print*, l0,l1
          allocate(tmp0(l0),tmp1(l1))
          tmp0 = ad0(dist(lct0)+1:dist(lct0+1)) 
          tmp1 = ad0(dist(lct1)+1:dist(lct1+1)) 

          call pnm_find_common(l0,tmp0,l1,tmp1,lc,com,cout)
          !if (unstrM%rank.eq.1) print*, 'a',tmp0,'b',tmp1,'c'

          cnum0(i) = lc
          !print*,'a', lc,'b', cout, 'c' 
          !print*, lc
          deallocate(tmp0,tmp1,com,cout) 
       enddo

       allocate(edge0%edist(edge0%nedg+1)); edge0%edist = 0
       do i = 1,edge0%nedg
          edge0%edist(i+1) = edge0%edist(i) + cnum0(i)
       enddo
       edge0%ed2esiz = edge0%edist(edge0%nedg+1)
       allocate(edge0%ed2e(edge0%ed2esiz))  

       do i = 1,edge0%nedg
          lct0 = ciod0((i-1)*2+1); lct1 = ciod0((i-1)*2+2)
          l0 = dat0(lct0); l1  = dat0(lct1)
          !print*, l0,l1
          allocate(tmp0(l0),tmp1(l1))
          tmp0 = ad0(dist(lct0)+1:dist(lct0+1)) 
          tmp1 = ad0(dist(lct1)+1:dist(lct1+1)) 

          call pnm_find_common(l0,tmp0,l1,tmp1,lc,com,cout)
          !if (unstrM%rank.eq.1) print*, 'a',tmp0,'b',tmp1,'c'

          k = edge0%edist(i) 
          edge0%ed2e(k+1:k+lc) = cout

          deallocate(tmp0,tmp1,com,cout) 
       enddo
       !print*, edge0%ed2e
    endif

    ! establish orgn
    if (edge0%nedg.gt.0) then
       unstrM%orgn%cnt = unstrM%org%cnt + edge0%ed2esiz
    else
       unstrM%orgn%cnt = unstrM%org%cnt
    endif

    allocate(unstrM%orgn%v2edist(unstrM%orgn%nvtx+1))
    allocate(unstrM%orgn%v2e(unstrM%orgn%cnt)) 
     
    unstrM%orgn%v2edist(1:unstrM%org%nvtx+1) = unstrM%org%v2edist
    unstrM%orgn%v2e(1:unstrM%org%cnt) = unstrM%org%v2e

    if (edge0%nedg.gt.0) then
       unstrM%orgn%v2edist(unstrM%org%nvtx+1:unstrM%orgn%nvtx+1) =&
        edge0%edist + unstrM%orgn%v2edist(unstrM%org%nvtx+1)
       unstrM%orgn%v2e(unstrM%org%cnt+1:unstrM%orgn%cnt) = edge0%ed2e
    endif


  end subroutine pnm_find_elements

  !--------------------------------------------------------------------
  subroutine pnm_find_edges(ed0,ed1)
    type(edinfo), intent(inout)               :: ed0,ed1

    integer                                   :: i,j,k,l,ierr,check
    integer                                   :: lvts0,lvts1,lct0,lct1
    integer, allocatable, dimension(:)        :: idtmp,dat0,dat1,npid0,npid1
    integer, allocatable, dimension(:)        :: cnt0,cnt1,tmpd0,tmpd1,vts0,vts1
    

    allocate(npid0(unstrM%nproc),npid1(unstrM%nproc))
    npid0 = 0; npid1 = 0
 
    if (ed1%nedg.gt.0) then
       allocate(tmpd0(ed1%nedg)); tmpd0 = ed1%pid
       call simplessort(tmpd0,idtmp)
       !print*, maxval(tmpd0-edall%pid(idtmp))
       allocate(vts0(2*ed1%nedg))
       vts0(1:2*ed1%nedg-1:2) = ed1%ed2v(idtmp,1)
       vts0(2:2*ed1%nedg:2)   = ed1%ed2v(idtmp,2)

       do i = 1,ed1%nedg
          j = tmpd0(i) + 1
          npid0(j) = npid0(j) + 1         
       enddo      
    endif
   
    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,unstrM%comm,ierr) 
    !print*,npid1,unstrM%rank
    npid0 = npid0*2; npid1 = npid1*2 
    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    lvts0 = sum(npid0); lvts1 = sum(npid1)

    allocate(vts1(lvts1))
    if (.not.allocated(vts0)) allocate(vts0(lvts0)) 
 
    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,unstrM%comm,ierr)
  
    !print*, vts0 
    allocate(dat0(lvts0/2),dat1(lvts1/2))
    
    do i = 1,lvts1/2
       lct0 = vts1((i-1)*2+1); lct1 = vts1(i*2) 
       if (ed0%nedg.gt.0) then
          check = 0
          do j = 1,ed0%nedg 
             if (ed0%ed2v(j,1) .eq. lct0 .and. &
                 ed0%ed2v(j,2) .eq. lct1) then
                dat1(i) = ed0%list(j)
                check   = 1
             endif
          enddo
          if (check.ne.1) then
             print*, "Error: find edge",lct0,lct1,int(unstrM%rank,4)
          endif
       endif
    enddo 
    
    cnt0  = cnt0/2;  cnt1  = cnt1/2
    npid0 = npid0/2; npid1 = npid1/2

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_integer,&
                       dat0,npid0,cnt0,mpi_integer,unstrM%comm,ierr)

    if (ed1%nedg.gt.0) then
       ed1%list(idtmp) = dat0
       !print*,edall%list,int(unstrM%rank,2) 
       !print*,minval(edall%list),int(unstrM%rank,2) 
       !print*, size(edall%list),edall%nedg
    endif


  end subroutine pnm_find_edges
  

 

  ! find the original proc. ids
  subroutine pnm_orgpid(lsiz0,dat0,pdat)
    integer, intent(in)                       :: lsiz0
    integer, intent(in)                       :: dat0(lsiz0)
    integer, intent(out)                      :: pdat(lsiz0)

    integer                                   :: i,j,k,l
    integer                                   :: lvts,lvts0,lvts1,ierr
   
    integer, allocatable, dimension(:)        :: vlist,vts0,v2vpid0,v2vpid1,idtmp

    allocate(vts0(lsiz0),v2vpid0(lsiz0),v2vpid1(lsiz0))
    vts0 = dat0
    call simplessort(vts0,idtmp)
    i = 1; j = 1; k = 1
    do while (i.le.lsiz0)    
       if (vts0(i).gt.unstrM%org%vtxdist(j).and.&
           vts0(i).le.unstrM%org%vtxdist(j+1)) then
          v2vpid0(k) = j - 1
          k = k + 1
          i = i + 1
       else
          if (j < unstrM%nproc) then
             j = j + 1
          endif 
       endif
    enddo
    v2vpid1(idtmp) = v2vpid0
    pdat = v2vpid1

    deallocate(vts0,v2vpid0,v2vpid1,idtmp)
 
  end subroutine pnm_orgpid

  ! find the original proc. ids
  subroutine pnm_org_epid(lsiz0,dat0,pdat)
    integer, intent(in)                       :: lsiz0
    integer, intent(in)                       :: dat0(lsiz0)
    integer, intent(out)                      :: pdat(lsiz0)

    integer                                   :: i,j,k,l
    integer                                   :: lvts,lvts0,lvts1,ierr
   
    integer, allocatable, dimension(:)        :: vlist,vts0,v2vpid0,v2vpid1,idtmp

    allocate(vts0(lsiz0),v2vpid0(lsiz0),v2vpid1(lsiz0))
    vts0 = dat0
    call simplessort(vts0,idtmp)
    i = 1; j = 1; k = 1
    do while (i.le.lsiz0)    
       if (vts0(i).gt.unstrM%org%edist(j).and.&
           vts0(i).le.unstrM%org%edist(j+1)) then
          v2vpid0(k) = j - 1
          k = k + 1
          i = i + 1
       else
          if (j < unstrM%nproc) then
             j = j + 1
          endif 
       endif
    enddo
    v2vpid1(idtmp) = v2vpid0
    pdat = v2vpid1

  end subroutine pnm_org_epid



  subroutine pnm_constct_cvpid(v2vsiz,v2v,v2vpid,cnvtx,Cvlist,Cvpid)
    integer, intent(in)                  :: v2vsiz
    integer, intent(in)                  :: v2v(v2vsiz),v2vpid(v2vsiz)

    integer, intent(out)                 :: cnvtx
    integer, allocatable, intent(out)    :: Cvlist(:),Cvpid(:)

    integer                              :: i,j,k,ierr
    integer, allocatable, dimension(:)   :: vlist1,vlist2,idtmp

    allocate(vlist1(v2vsiz),vlist2(v2vsiz))
    vlist1 = v2v
    call simplessort(vlist1,idtmp) 
    call simplepickunique(vlist1,v2vsiz,Cvlist,cnvtx)
    call mpi_barrier(unstrM%comm,ierr)

    do i = 1,v2vsiz
       j = idtmp(i); vlist2(i) = v2vpid(j)
    enddo
    allocate(Cvpid(cnvtx))

    i = 1; j = 1; Cvpid(1) = vlist2(1)
    do while (i.le.v2vsiz-1) 
       if (vlist1(i).eq.vlist1(i+1)) then
          i = i + 1
       else
          j = j + 1; i = i + 1
          Cvpid(j) = vlist2(i) 
       endif
    enddo

    if (j.ne.cnvtx) then
       print*, j,cnvtx,'Error, local vtx at rank', int(unstrM%rank,4)
       stop 
    endif

  end subroutine pnm_constct_cvpid

  subroutine pnm_setup_lv2v(tgt,ll,vnum,sdat,lsiz)
    integer, intent(in)                  :: tgt,ll,vnum
    integer, intent(in)                  :: sdat(ll)
    integer, intent(out)                 :: lsiz
    !integer, allocatable, intent(out)    :: v2v(:)

    integer                              :: i,j,k,l,check
    integer                              :: lsiz0,lsiz1
    integer, allocatable, dimension(:)   :: tmp0,tmp1,idt

    lsiz0 = ll/vnum

    do i = 1,lsiz0
       check = 0 
       do j = 1,vnum
          k = (i-1)*vnum+j 
          if (sdat(k).eq.tgt) then
             check = 1
          endif
       enddo
       if (check.eq.0) then
          print*, 'Error: setup lv2v',tgt,sdat(k-vnum+1:k)
       endif 
    enddo 
 
    allocate(tmp0(ll)); tmp0 = sdat
    call simplessort(tmp0,idt) 

    call simplepickunique(tmp0,ll,tmp1,lsiz)

  end subroutine pnm_setup_lv2v

  subroutine pnm_find_lv2v(ll,sdat,spid,lsiz,v2v,v2vpid)
    integer, intent(in)                  :: ll
    integer, intent(in)                  :: sdat(ll),spid(ll)
    integer, intent(out)                 :: lsiz
    integer, allocatable, intent(out)    :: v2v(:),v2vpid(:)

    integer                              :: i,j,k,l,check
    integer                              :: lsiz0,lsiz1
    integer, allocatable, dimension(:)   :: tmp0,tmp1,tmp2,idt

    allocate(tmp0(ll)); tmp0 = sdat
    call simplessort(tmp0,idt) 
    allocate(tmp2(ll)); tmp2 = spid(idt)

    call simplepickunique(tmp0,ll,v2v,lsiz)

    call simplepickuniord(tmp0,ll,tmp1,lsiz1)
 
    if (lsiz.ne.lsiz1) print*,'error: find lv2v'

    allocate(v2vpid(lsiz))
    v2vpid = tmp2(tmp1)

  end subroutine pnm_find_lv2v


  !---------------------------------------------------------------------
  subroutine transfMat(tetidx,B,Binv,detB)
    ! compute Jacobian, inverse Jacobian and the determinent of the Jacobian
    ! input:  local tetidx 
    ! output: Jacobian information 
    integer,intent(in)              :: tetidx
    real(kind=rkind),intent(out)    :: B(3,3), Binv(3,3),detB
    integer                         :: i,locvert(4)

    !locvert = unstrM%lt2vid((tetidx-1)*4+1:tetidx*4)
    !B(1,:)  = unstrM%Cv_crs(:,locvert(2)) - unstrM%Cv_crs(:,locvert(1))
    !B(2,:)  = unstrM%Cv_crs(:,locvert(3)) - unstrM%Cv_crs(:,locvert(1))
    !B(3,:)  = unstrM%Cv_crs(:,locvert(4)) - unstrM%Cv_crs(:,locvert(1))
    i = tetidx - 1
    B(1,:)  = unstrM%loc_vtx(:,i*4+2) - unstrM%loc_vtx(:,i*4+1)
    B(2,:)  = unstrM%loc_vtx(:,i*4+3) - unstrM%loc_vtx(:,i*4+1)
    B(3,:)  = unstrM%loc_vtx(:,i*4+4) - unstrM%loc_vtx(:,i*4+1)

    call matdet(B,detB)

    if (detB.eq.0) then
      print*, 'Error: the ',tetidx,'-th tetrahedron is flat'
      stop
    end if
    call matinv(B,Binv,3)
  end subroutine transfMat
  
  !--------------------------------------------------------------------
  subroutine node_coo(tetnodes,tetid,pNp)
    ! compute node coordinates for tet id: tetid
    ! input:  tetid
    ! output: tetnodes all nodes coordinates 
    integer, intent(in)              :: tetid,pNp
    real(kind=rkind), intent(out)    :: tetnodes(3,pNp)
     
    integer                          :: k,j
    real(kind=rkind)                 :: a(3),B(3,3),invB(3,3),detB
    !k = unstrM%lt2vid((tetid-1)*4+1)
    !a = unstrM%Cv_crs(:,k)
    a = unstrM%loc_vtx(:,(tetid-1)*4+1)
    call transfMat(tetid,B,invB,detB)
    B = transpose(B/2.0D0)
    do j=1,pNp
       ! X=A+B(X_r+1)/2 : X is a pt in tet; X_r is a pt in ref tet.
       tetnodes(:,j)   = a + matmul(B,refs%nodes(:,j)+1.0D0)
    enddo 

  end subroutine node_coo

  !---------------------------------------------------------------------
  subroutine face_normal(pt,n,sJac)
    ! imput: 4 points of a tet 
    ! output: normal vectors areas of 4 faces
    ! see FtoV their relation
    implicit none
    real(kind=rkind), intent(in)     :: pt(3,4)
    real(kind=rkind), intent(out)    :: n(3,4), sJac(4)

    integer                          :: i,j,k
    real(kind=rkind)                 :: v12(3),v13(3),vref(3),v(3,3)

    do i = 1,4 
       vref(:) = pt(:,i)
       k       = 1
       do j = 1,4
          if (j /= i) then
             v(:,k) = pt(:,j)
             k      = k + 1
          endif
       enddo
       v12(:)  = v(:,2) - v(:,1)
       v13(:)  = v(:,3) - v(:,1)
    
       n(1,i)  = v12(2)*v13(3) - v12(3)*v13(2)
       n(2,i)  = v12(3)*v13(1) - v12(1)*v13(3)
       n(3,i)  = v12(1)*v13(2) - v12(2)*v13(1)
       !sJac(i) = norm2(n(:,i))    
       sJac(i) = sqrt(sum((n(:,i))**2))    
       n(:,i)  = n(:,i)/sJac(i)
       
       if (abs(sJac(i)) < pin%TOL*1.0D-3) then
          print*, "ERROR: Face area close to zero"
          stop
       endif
       vref = vref - v(:,1)
       
       if (sum(vref*n(:,i)) > 0 ) n(:,i) = -n(:,i) 
    enddo

  end subroutine face_normal

  !---------------------------------------------------------------------
  subroutine Build_reference(pOrder,Nfp,pNp,reference)
    ! build reference information for refs
    integer, intent(in)                               :: pOrder,Nfp,pNp
    type(tetrahedron), intent(inout)                  :: reference    
 
    integer                                           :: i,j,k,l,n
    integer, dimension(Nfp)                           :: fm
    real(kind=rkind), dimension(pNp)                  :: r,s,t
    real(kind=rkind), dimension(pNp)                  :: x1,x2,x3
    real(kind=rkind), dimension(Nfp)                  :: face1,face2
    real(kind=rkind)                                  :: Vinv(pNp,pNp)  
    real(kind=rkind), dimension(Nfp,Nfp)              :: V2D,Vtmp
    real(kind=rkind), dimension(pNp,pNp)              :: D1,D2,D3

    ! edges defined by verts
    reference%EtoV(:,1)  =    (/1, 1, 1, 2, 2, 3/)
    reference%EtoV(:,2)  =    (/2, 3, 4, 3, 4, 4/)
    ! Edge 1 -> Nodes 1 2
    ! Edge 2 -> Nodes 1 3
    ! Edge 3 -> Nodes 1 4
    ! Edge 4 -> Nodes 2 3
    ! Edge 5 -> Nodes 2 4
    ! Edge 6 -> Nodes 3 4

    ! faces defined by edges
    reference%FtoE(:,1)  =   (/4, 2, 1, 1/)  ! each col defines a face by edge
    reference%FtoE(:,2)  =   (/5, 3, 3, 2/) 
    reference%FtoE(:,3)  =   (/6, 6, 5, 4/) 
    ! Face 1 -> edges 4, 5, 6 verts: 2, 3, 4
    ! Face 2 -> edges 2, 3, 6 verts: 1, 3, 4
    ! Face 3 -> edges 1, 3, 5 verts: 1, 2, 4
    ! Face 4 -> edges 1. 2, 4 verts: 1, 2, 3

    ! faces defined by verts
    reference%FtoV(:,1)  =   (/2, 1, 1, 1/) ! each col defines a face by vert
    reference%FtoV(:,2)  =   (/3, 3, 2, 2/)
    reference%FtoV(:,3)  =   (/4, 4, 4, 3/)
    ! Face 1 -> verts 2, 3, 4: lack of 1
    ! Face 2 -> verts 1, 3, 4: lack of 2
    ! Face 3 -> verts 1, 2, 4: lack of 3
    ! Face 4 -> verts 2, 3, 4: lack of 4
                    
    !print*, refs%EtoV,refs%FtoE,refs%FtoV

    !-------------------------
    ! build reference vert list
    reference%vert(:)        = (/1,2,3,4/)
    reference%vert_crs(:,1)  = (/-1.0D0,-1.00D0,-1.0D0/) 
    reference%vert_crs(:,2)  = (/+1.0D0,-1.00D0,-1.0D0/) 
    reference%vert_crs(:,3)  = (/-1.0D0,+1.00D0,-1.0D0/) 
    reference%vert_crs(:,4)  = (/-1.0D0,-1.00D0,+1.0D0/) 
    
    !-------------------------
    ! build reference edge list
    reference%edge(:)        = (/1,2,3,4,5,6/)
    
    !--------------e----------
    ! build reference face list
    reference%face(:)        = (/1,2,3,4/)
    ! the relation between face, edge & vert see EtoV, FtoV, FtoE  
    
    ! sJac
    reference%sJac(:)        = (/2.0D0,2.0D0,2.0D0,2.0D0*sqrt(3.0D0)/)    
    
    ! normal vectors outwards
    reference%n(:,1)         = (/ 1.0D0, 1.0D0, 1.0D0/)/sqrt(3.0D0)
    reference%n(:,2)         = (/-1.0D0, 0.0D0, 0.0D0/)
    reference%n(:,3)         = (/ 0.0D0,-1.0D0, 0.0D0/)
    reference%n(:,4)         = (/ 0.0D0, 0.0D0,-1.0D0/)
    
    allocate(reference%nodes(3,pNp))
    !-------------------------
    ! build basis nodal points
    call blend_nodes(r,s,t,pNp,pOrder)
    !print*, "blend_nodes done"
    do i=1,pNp
        reference%nodes(:,i) = (/r(i),s(i),t(i)/)
    end do

    allocate(reference%Fmask(Nfp,4))
    reference%Fmask = -1
    ! Fmask -------------------
    ! build basis nodal points on faces
    ! face_id FtoV
    i=1; j=1; k=1; l=1;
    do n=1,pNp
        if(dabs(1.0D0+r(n)).le.pin%TOL)then
           reference%Fmask(i,2)=n; i=i+1
        endif
        if(dabs(1.0D0+s(n)).le.pin%TOL)then
           reference%Fmask(j,3)=n; j=j+1
        endif
        if(dabs(1.0D0+t(n)).le.pin%TOL)then
           reference%Fmask(l,4)=n; l=l+1
        endif
        if(dabs(1.0D0+r(n)+s(n)+t(n)).le.pin%TOL)then
           reference%Fmask(k,1)=n; k=k+1
        endif
    enddo

    if (i /= Nfp+1 .or. j/= Nfp+1 .or.&
        k /= Nfp+1 .or. l/= Nfp+1) then
       print*, "number of face nodes is wrong"
       stop
    endif

    ! build nodes status
    allocate(reference%ns(pNp)); refs%ns = 4 
    do i = 1,4
       reference%ns(reference%Fmask(:,i)) = 3
    enddo

    do i = 1,pNp
       if (abs(1.0+r(i))<pin%TOL .and. abs(1.0+s(i))<pin%TOL) then
          reference%ns(i) = 2
       endif
       if (abs(1.0+r(i))<pin%TOL .and. abs(1.0+t(i))<pin%TOL) then
          reference%ns(i) = 2
       endif
       if (abs(1.0+t(i))<pin%TOL .and. abs(1.0+s(i))<pin%TOL) then
          reference%ns(i) = 2
       endif 
       if (abs(1.0+r(i))<pin%TOL .and. abs(s(i)+t(i))<pin%TOL) then
          reference%ns(i) = 2
       endif 
       if (abs(1.0+s(i))<pin%TOL .and. abs(r(i)+t(i))<pin%TOL) then
          reference%ns(i) = 2
       endif
       if (abs(1.0+t(i))<pin%TOL .and. abs(s(i)+r(i))<pin%TOL) then
          reference%ns(i) = 2
       endif
       if (i == 1 .or. i == 1 + pOrder .or.&
           i == Nfp .or. i == pNp) then 
          reference%ns(i) = 1
       endif
    enddo
    ! check node status
    j = 0; k = 0; l = 0
    do i = 1,pNp
       if (reference%ns(i) == 2) j = j+1
       if (reference%ns(i) == 3) k = k+1
       if (reference%ns(i) == 4) l = l+1
    enddo
    if (j /= 6*(pOrder-1) .and. pOrder > 1) then 
       print*, "error: nodes status about edges is wrong"
    elseif (pOrder > 1 .and.&
            k /= 2*(pOrder-2)*(pOrder-1)) then 
       print*, "error: nodes status about faces is wrong"
    elseif (pOrder > 3 .and.&
           l /= (pOrder-3)*(pOrder-2)*(pOrder-1)/6) then
       print*, "error: nodes status about interior pts is wrong"
    endif
    !print*, refs%ns
    ! build on Mass matrix on reference
    allocate(reference%V3D(pNp,pNp),reference%MassM(pNp,pNp))
    allocate(reference%invV(pNp,pNp))

    x1 = reference%nodes(1,:) 
    x2 = reference%nodes(2,:) 
    x3 = reference%nodes(3,:) 

    call Vandermonde3D(reference%V3D,pOrder,x1,x2,x3,pNp)       
    call matinv(reference%V3D,Vinv,pNp)


    reference%invV  = Vinv
    reference%MassM = matmul(transpose(Vinv), Vinv)
    !print*, refs%MassM       
    !print*,reference%V3D

    ! build on spatial derivatives
    allocate(reference%Drst(pNp,pNp,3))
    call GradVandermonde3D(D1,D2,D3,pOrder,x1,x2,x3,pNp)

    allocate(reference%D1(pNp,pNp)); reference%D1 = D1
    allocate(reference%D2(pNp,pNp)); reference%D2 = D2
    allocate(reference%D3(pNp,pNp)); reference%D3 = D3

    reference%Drst(:,:,1) = matmul(D1,Vinv)
    reference%Drst(:,:,2) = matmul(D2,Vinv)
    reference%Drst(:,:,3) = matmul(D3,Vinv)
    !print*, refs%Drst(:,:,1)    

    ! build surface mass matrix
    allocate(reference%MassF(Nfp,Nfp,4)) 
    fm = reference%Fmask(:,1); face1 = s(fm); face2 = t(fm)
    call Vandermonde2D(V2D,pOrder,face1,face2,Nfp)
    vtmp  = matmul(V2D,transpose(V2D))
    call matinv(vtmp,reference%MassF(:,:,1),Nfp) 

    fm = reference%Fmask(:,2); face1 = s(fm); face2 = t(fm)
    call Vandermonde2D(V2D,pOrder,face1,face2,Nfp)
    vtmp  = matmul(V2D,transpose(V2D))
    call matinv(vtmp,reference%MassF(:,:,2),Nfp)

    fm = reference%Fmask(:,3); face1 = r(fm); face2 = t(fm)
    call Vandermonde2D(V2D,pOrder,face1,face2,Nfp)
    vtmp  = matmul(V2D,transpose(V2D))
    call matinv(vtmp,reference%MassF(:,:,3),Nfp)

    fm = reference%Fmask(:,4); face1 = r(fm); face2 = s(fm)
    call Vandermonde2D(V2D,pOrder,face1,face2,Nfp)
    vtmp  = matmul(V2D,transpose(V2D))
    call matinv(vtmp,reference%MassF(:,:,4),Nfp)

    reference%Fscale(:)    = reference%sJac(:)/(4.0D0/3.0D0) 

    if (pin%s%pOrder.eq.1) then
       reference%vord = (/1,2,3,4/)
    elseif (pin%s%pOrder.eq.2) then
       reference%vord = (/1,3,6,10/)
    endif

    !print*,reference%MassM    
    !print*, reference%Fmask
    ! Todo: fluid separation
  end subroutine Build_reference

  subroutine pnm_find_lNele()
    !integer, dimension(:,:), intent(in) :: e2v0
    !integer, dimension(:), intent(in)   :: Cvlist,Cvpid

    integer                           :: i,j,k,l,ll,lct,cout
    integer                           :: ierr,lids(4),l0(4),id0(4) 
    integer, allocatable              :: idtmp(:)
    cout = 0 
    do i = 1,unstrM%ClNele
       do j = 1,4 
          k = unstrM%loc_t2v(refs%vord(j),i) 
          call findorder(k,unstrM%Cvlist,lct)
          lids(j) = unstrM%Cvpid(lct)
          !k = e2v0(j,i) 
          !call findorder(k,Cvlist,lct)
          !lids(j) = Cvpid(lct)
       enddo
       call simplessort(lids,idtmp); deallocate(idtmp)
       ll = 0 
       do l = 1,4
          if (lids(l).eq.unstrM%rank) ll = ll + 1
       enddo

       if (ll.ge.3) then 
          cout = cout + 1
       elseif (ll.eq.2) then
          if (unstrM%rank.gt.lids(1)) then
             cout = cout + 1
          elseif (lids(3).ne.lids(4)) then
             cout = cout + 1
          endif
       elseif (ll.eq.1) then
          if (lids(1).lt.lids(2).and.lids(2).lt.lids(3)&
                                .and.lids(3).lt.lids(4)) then
             if (unstrM%rank.eq.lids(1)) then
                cout = cout + 1
             endif
             !print*,unstrM%Clelist(i),k
          endif
       endif
    enddo
 
    unstrM%lNele = cout
    call mpi_allreduce(unstrM%lNele,cout,1,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 

    print*,'# of elems (FMM)',unstrM%lNele,unstrM%ClNele,unstrM%rank

    if (cout.ne.unstrM%Ntet) stop

    allocate(unstrM%lelist(unstrM%lNele))
    allocate(unstrM%ClNids(unstrM%ClNele))   

    cout = 0 
    do i = 1,unstrM%ClNele
       do j = 1,4 
          k = unstrM%loc_t2v(refs%vord(j),i) 
          call findorder(k,unstrM%Cvlist,lct)
          lids(j) = unstrM%Cvpid(lct)
          !k = e2v0(j,i) 
          !call findorder(k,Cvlist,lct)
          !lids(j) = Cvpid(lct)
       enddo
       call simplessort(lids,idtmp); 
       id0=idtmp; deallocate(idtmp)
       !print*,id0 
       ll = 0 
       do l = 1,4
          if (lids(l).eq.unstrM%rank) ll = ll + 1
       enddo
       !if (unstrM%Clelist(i).eq.9) print*,lids,unstrM%rank


       if (ll.ge.3) then 
          cout = cout + 1
          unstrM%lelist(cout) = i
          unstrM%ClNids(i) = unstrM%rank
       elseif (ll.eq.2) then
          if (unstrM%rank.gt.lids(1)) then
             cout = cout + 1
             unstrM%lelist(cout) = i
             unstrM%ClNids(i) = unstrM%rank
          else
             if (lids(3).ne.lids(4)) then
               cout = cout + 1
               unstrM%lelist(cout) = i
               unstrM%ClNids(i) = unstrM%rank
             elseif (lids(3).eq.lids(4)) then
               unstrM%ClNids(i) = lids(4)
             endif
          endif
       elseif (ll.eq.1) then
          if (lids(1).lt.lids(2).and.lids(2).lt.lids(3)&
                                .and.lids(3).lt.lids(4)) then
             if (unstrM%rank.eq.lids(1)) then
                cout = cout + 1
                unstrM%lelist(cout) = i
                !l0 = lids(id0)
                !print*,unstrM%Clelist(i),unstrM%ClNeigh(:,i),l0
             endif
             unstrM%ClNids(i) = lids(1)
          elseif (lids(1).eq.unstrM%rank.and.&
                 (lids(3).eq.lids(2).or.lids(3).eq.lids(4))) then
             unstrM%ClNids(i) = lids(3)
          elseif (lids(2).eq.unstrM%rank.and.lids(3).eq.lids(4)) then
             unstrM%ClNids(i) = lids(3)
          elseif (lids(3).eq.unstrM%rank.and.lids(1).eq.lids(2)) then
             unstrM%ClNids(i) = lids(2)
          elseif (lids(4).eq.unstrM%rank.and.&
                 (lids(2).eq.lids(3).or.lids(2).eq.lids(1))) then
             unstrM%ClNids(i) = lids(2)
          else
             unstrM%ClNids(i) = lids(1) 
          endif
       elseif (ll.eq.0) then
          if (lids(1).lt.lids(2).and.lids(2).lt.lids(3)&
                                .and.lids(3).lt.lids(4)) then
             unstrM%ClNids(i) = lids(1)
          elseif (lids(1).eq.lids(2).and.lids(3).ne.lids(4)) then 
             unstrM%ClNids(i) = lids(1)
          elseif (lids(1).eq.lids(2).and.lids(3).eq.lids(4)) then 
             unstrM%ClNids(i) = lids(3)
          elseif (lids(1).lt.lids(2).and.lids(2).eq.lids(3)) then
             unstrM%ClNids(i) = lids(2)
          elseif (lids(1).lt.lids(2).and.lids(2).lt.lids(3)&
                                    .and.lids(3).eq.lids(4)) then
             unstrM%ClNids(i) = lids(3)
          endif
       endif


       !if (ll.lt.1) then
       !   print*, "Error: ClNids", unstrM%rank
       !   stop
       !endif
    enddo

    ! find order

    !if(unstrM%rank.eq.0) print*,unstrM%lelist

    ! source locations
    allocate(unstrM%sloc(3,unstrM%lNele))
    do i = 1,unstrM%lNele
       j = unstrM%lelist(i) 
       do k = 1,3
          unstrM%sloc(k,i) = sum(unstrM%loc_vtx(k,(j-1)*4+1:j*4))/4.0D0
       enddo
    enddo
    !print*, unstrM%sloc(:,:)
     
  end subroutine pnm_find_lNele 

  ! find surfaces
  subroutine pnm_find_lsurf()
    integer                            :: i,j,k,l,m,n,kk,ll,ierr
    integer                            :: lsiz0,lsiz1,lct,lctp
    integer                            :: c0,c1,cout,sid(3)
    integer                            :: lct0,lct1,lct2 
    integer, allocatable, dimension(:) :: nid,tpid0,tpid1,ck0
 
    ! create a graph via ClNeigh
    ll = 0
    do i = 1,unstrM%lNele; do j = 1,4
       k = unstrM%Clelist(unstrM%lelist(i))
       call findorder(k,unstrM%Clelist,lct)
       if (unstrM%ClNeigh(j,lct).gt.0) ll = ll + 1 
    enddo; enddo

    l = unstrM%lNele*4 - ll

    call mpi_allreduce(ll,lsiz0,1,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 
   
    call mpi_allreduce( l,lsiz1,1,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 
    unstrM%Nsrf = lsiz0/2 + lsiz1

    if(unstrM%rank.eq.0) print*,'# of surfaces', unstrM%Nsrf 

    ! check loc # of surfaces
    kk = 0
    do i = 1,unstrM%ClNele; do j = 1,4
       k = unstrM%Clelist(i)
       if ((unstrM%ClNeigh(j,i).lt.k.and.&
          unstrM%ClNeigh(j,i).gt.0).or.&
          unstrM%ClNeigh(j,i).eq.-1) then
          ll = 0; m = 0 
          do l = 1,4
             if (l.ne.j) then
                ll = ll + 1
                n = unstrM%loc_t2v(refs%vord(l),i) 
                call findorder(n,unstrM%Cvlist,lct)
                sid(ll) = unstrM%Cvpid(lct)
                if (unstrM%Cvpid(lct).eq.unstrM%rank) then
                   m = m + 1
                endif 
             endif
          enddo
          call simplessort(sid,nid); deallocate(nid) 

          if (m.eq.1.and.sid(2).ne.sid(3).and.&
             sid(1).eq.unstrM%rank) then
             kk = kk + 1
          elseif (m.gt.1) then
             kk = kk + 1
          elseif (m.eq.0) then
             !print*,"Error: vertex pid 0",unstrM%rank,sid
          endif
       endif
    enddo; enddo

    unstrM%srf%lsf = kk  
    call mpi_allreduce(unstrM%srf%lsf,lsiz0,1,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr)
 
    !if(unstrM%rank.eq.0) print*,lsiz0,unstrM%Nsrf

    allocate(tpid0(unstrM%nproc)); tpid0 = 0 
    allocate(tpid1(unstrM%nproc)); tpid1 = 0 
    tpid0(unstrM%rank+1) = kk
    call mpi_allreduce(tpid0,tpid1,unstrM%nproc,mpi_integer,&
                       mpi_sum,unstrM%comm,ierr) 
    
    allocate(unstrM%srf%dist(unstrM%nproc+1))
    unstrM%srf%dist(1) = 0
    do i = 1,unstrM%nproc
       unstrM%srf%dist(i+1) = unstrM%srf%dist(i) + tpid1(i)
    enddo 
    

    if (lsiz0.ne.unstrM%Nsrf) then
       if (unstrM%rank.eq.0) then
          print*, 'Error, # of surface',lsiz0,unstrM%Nsrf
          stop
       endif
    endif


    allocate(unstrM%srf%s2e(2,unstrM%srf%lsf))
    allocate(unstrM%srf%inx(unstrM%srf%lsf))
    allocate(unstrM%srf%pid(unstrM%srf%lsf)); unstrM%srf%pid = unstrM%rank
    allocate(unstrM%srf%ids(unstrM%srf%lsf))
    allocate(unstrM%srf%cen(3,unstrM%srf%lsf))
    allocate(unstrM%srf%drt(3,unstrM%srf%lsf))

    do i = 1,unstrM%srf%lsf
       unstrM%srf%ids(i) = unstrM%srf%dist(unstrM%rank+1)+i
    enddo


    allocate(unstrM%srf%snum(unstrM%ClNele)); unstrM%srf%snum = 0
    allocate(unstrM%srf%sdst(unstrM%ClNele+1)); unstrM%srf%sdst(1) = 0

    ! construct surface info
    kk = 0
    do i = 1,unstrM%ClNele
       c0 = 0
       do j = 1,4
          k = unstrM%Clelist(i)
          if ((unstrM%ClNeigh(j,i).lt.k.and.&
             unstrM%ClNeigh(j,i).gt.0).or.&
             unstrM%ClNeigh(j,i).eq.-1) then
             ll = 0; m = 0 
             do l = 1,4
                if (l.ne.j) then
                   ll = ll + 1
                   n = unstrM%loc_t2v(refs%vord(l),i) 
                   call findorder(n,unstrM%Cvlist,lct)
                   sid(ll) = unstrM%Cvpid(lct)
                   if (unstrM%Cvpid(lct).eq.unstrM%rank) then
                      m = m + 1
                   endif 
                endif
             enddo
             call simplessort(sid,nid); deallocate(nid) 
             if ((m.eq.1.and.sid(2).ne.sid(3).and.sid(1).eq.unstrM%rank)&
                .or.m.gt.1) then
                kk = kk + 1; c0 = c0 + 1
                unstrM%srf%s2e(1,kk) = k
                unstrM%srf%s2e(2,kk) = unstrM%ClNeigh(j,i)
                unstrM%srf%inx(kk)   = j
                unstrM%srf%drt(:,kk) = unstrM%loc_n(:,j,i)
                do m = 1,3
                   unstrM%srf%cen(m,kk) = &
                   sum(unstrM%loc_vtx(m,4*(i-1)+refs%FtoV(j,:)))/3.0D0
                   !print*,unstrM%loc_vtx(m,4*(i-1)+refs%FtoV(j,:)) 
                   !print*,unstrM%loc_vtx(m,4*(i-1)+refs%Fmask(:,j)) 
                enddo
             endif
          endif
       enddo
       unstrM%srf%snum(i) = c0
       !print*,k,c0,unstrM%rank
    enddo

    !print*,unstrM%srf%s2e(1,1),unstrM%rank
    !print*, unstrM%srf%inx(:)

    do i = 1,unstrM%ClNele
       unstrM%srf%sdst(i+1) = unstrM%srf%sdst(i) + unstrM%srf%snum(i)
    enddo

    ! deal with csrf
    unstrM%csrf%lsf = unstrM%ClNele*4 

    allocate(unstrM%csrf%s2e(2,unstrM%csrf%lsf))
    allocate(unstrM%csrf%inx(unstrM%csrf%lsf))
    allocate(unstrM%csrf%pid(unstrM%csrf%lsf))
    allocate(unstrM%csrf%ids(unstrM%csrf%lsf))
    allocate(unstrM%csrf%cen(3,unstrM%csrf%lsf))
    allocate(unstrM%csrf%drt(3,unstrM%csrf%lsf))

    do i = 1,unstrM%ClNele
       k = unstrM%Clelist(i)
       do j = 1,4
          kk = (i-1)*4+j
          if (unstrM%ClNeigh(j,i).lt.k) then
             unstrM%csrf%s2e(1,kk) = k
             unstrM%csrf%s2e(2,kk) = unstrM%ClNeigh(j,i)
             unstrM%csrf%inx(kk)   = j
             unstrM%csrf%drt(:,kk) = unstrM%loc_n(:,j,i)
          else         
             unstrM%csrf%s2e(1,kk) = unstrM%ClNeigh(j,i)
             unstrM%csrf%s2e(2,kk) = k 
             unstrM%csrf%inx(kk)   = -j
             unstrM%csrf%drt(:,kk) = -unstrM%loc_n(:,j,i)
          endif
          do m = 1,3
             unstrM%csrf%cen(m,kk) = &
             sum(unstrM%loc_vtx(m,4*(i-1)+refs%FtoV(j,:)))/3.0D0 
          enddo

          ! find out pid
          ll = 0; m = 0 
          do l = 1,4
             if (l.ne.j) then
                ll = ll + 1
                n = unstrM%loc_t2v(refs%vord(l),i) 
                call findorder(n,unstrM%Cvlist,lct)
                sid(ll) = unstrM%Cvpid(lct)
                if (unstrM%Cvpid(lct).eq.unstrM%rank) then
                   m = m + 1
                endif 
             endif
          enddo
          call simplessort(sid,nid); deallocate(nid) 
          if (sid(2).eq.sid(3)) then
             unstrM%csrf%pid(kk) = sid(3)
          else 
             unstrM%csrf%pid(kk) = sid(1)
          endif
       enddo
    enddo

    call pnm_find_srf(unstrM%srf,unstrM%csrf) 

    call pnm_pre_fmm()


  end subroutine pnm_find_lsurf

  subroutine pnm_find_srf(sf0,sf1)
    type(surface), intent(inout)              :: sf0,sf1

    integer                                   :: i,j,k,l,dm,ierr,lct,check
    integer                                   :: lvts0,lvts1,lct0,lct1,lct2
    integer, allocatable, dimension(:)        :: idtmp,dat0,dat1,npid0,npid1
    integer, allocatable, dimension(:)        :: cnt0,cnt1,tmpd0,tmpd1,vts0,vts1

    allocate(npid0(unstrM%nproc),npid1(unstrM%nproc))
    npid0 = 0; npid1 = 0

    allocate(tmpd0(sf1%lsf)); tmpd0 = sf1%pid
    call simplessort(tmpd0,idtmp)
    !print*, maxval(tmpd0-edall%pid(idtmp))
    dm = 3
    allocate(vts0(dm*sf1%lsf))
    vts0(1:dm*sf1%lsf-2:dm) = sf1%s2e(1,idtmp)
    vts0(2:dm*sf1%lsf-1:dm) = sf1%s2e(2,idtmp)
    vts0(3:dm*sf1%lsf-0:dm) = sf1%inx(idtmp)

    do i = 1,sf1%lsf
       j = tmpd0(i) + 1; npid0(j) = npid0(j) + 1         
    enddo      
    !print*,npid0,unstrM%rank

    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,unstrM%comm,ierr) 

    npid0 = npid0*dm; npid1 = npid1*dm 
    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    lvts0 = sum(npid0); lvts1 = sum(npid1)

    allocate(vts1(lvts1))
    if (.not.allocated(vts0)) allocate(vts0(lvts0)) 
 
    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,unstrM%comm,ierr)

    !print*,lvts0/dm,lvts1/dm,unstrM%rank
    allocate(dat0(lvts0/dm),dat1(lvts1/dm))
    
    do i = 1,lvts1/dm
       lct0 = vts1((i-1)*dm+1); lct1 = vts1((i-1)*dm+2) 
       lct2 = vts1((i-1)*dm+3)   
       check = 0
       call findorder(lct0,unstrM%Clelist,lct)
       if (sf0%snum(lct).gt.0) then
          do j = 1,sf0%snum(lct)
             k = sf0%sdst(lct)+j
             if (sf0%s2e(2,k).eq.lct1.and.lct1.eq.-1 &
                 .and.sf0%inx(k).eq.lct2) then
                dat1(i) = sf0%ids(k)
                check = 1
             elseif (sf0%s2e(2,k).eq.lct1.and.lct1.gt.0) then
                dat1(i) = sf0%ids(k)
                check = 1
             endif
          enddo
       else
          print*, "Error: find element",lct0,sf0%snum(lct)
       endif
       !do j = 1,sf0%lsf 
       !   if (sf0%s2e(1,j).eq.lct0 .and. &
       !       sf0%s2e(2,j).eq.lct1 .and. lct1.gt.0) then
       !      dat1(i) = sf0%ids(j)
       !      check   = 1
       !   elseif (sf0%s2e(1,j).eq.lct0 .and. &
       !       sf0%s2e(2,j).eq.lct1 .and. lct1.eq.-1) then
       !      if (lct2.eq.sf0%inx(j)) then
       !         dat1(i) = sf0%ids(j)
       !         check = 1
       !      endif
       !      !check = 1
       !   endif
       !enddo
       if (check.ne.1) then
          print*, "Error: find surf",lct0,lct1,int(unstrM%rank,4)
       endif
    enddo 

    cnt0  = cnt0/dm;  cnt1  = cnt1/dm
    npid0 = npid0/dm; npid1 = npid1/dm

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_integer,&
                       dat0,npid0,cnt0,mpi_integer,unstrM%comm,ierr)
 
    sf1%ids(idtmp) = dat0

 
  end subroutine pnm_find_srf

  ! prepare for fmm 
  subroutine pnm_pre_fmm()
    integer                                    :: i,j,k,l,m,n,ii,jj,kk,ll
    integer                                    :: lct,cck,check
    integer, allocatable                       :: fm(:)
    integer, allocatable, dimension(:)         :: idtmp,dtmp
    
    ! for esrf
    allocate(fm(pin%s%Nfp))
    allocate(dtmp(unstrM%csrf%lsf))
    dtmp = unstrM%csrf%ids
    call simplessort(dtmp,idtmp)
    ! estimate the size of esrf
    ll = 0; cck = -1 
    do i = 1,unstrM%csrf%lsf
       ii = idtmp(i)
       if (unstrM%csrf%ids(ii).ne.cck) then
          cck = unstrM%csrf%ids(ii)
          kk  = unstrM%csrf%inx(ii)
          check = 0
          fm = refs%Fmask(:,int(abs(kk)))
          !fm = refs%FtoV(int(abs(kk)),:)
          if (kk.lt.0) jj = 2
          if (kk.gt.0) jj = 1
          do j = 1,pin%s%Nfp
             l = unstrM%csrf%s2e(jj,ii)
             call findorder(l,unstrM%Clelist,m)
             n = unstrM%loc_t2v(fm(j),m)
             call findorder(n,unstrM%Cvlist,lct)
             if (unstrM%Cvpid(lct).eq.unstrM%rank) then
                check = 1              
             endif
          enddo 
          if (check.eq.1) ll = ll + 1
       else
          ! do nothing 
       endif
    enddo

    unstrM%esrf%lsf = ll
    !print*, ll
    print*,'# of loc surfaces',unstrM%srf%lsf,unstrM%esrf%lsf,int(unstrM%rank,4)

    allocate(unstrM%esrf%s2e(2,unstrM%esrf%lsf))
    allocate(unstrM%esrf%inx(unstrM%esrf%lsf))
    allocate(unstrM%esrf%pid(unstrM%esrf%lsf))
    allocate(unstrM%esrf%ids(unstrM%esrf%lsf))
    allocate(unstrM%esrf%cen(3,unstrM%esrf%lsf))
    allocate(unstrM%esrf%drt(3,unstrM%esrf%lsf))

    ll = 0; cck = -1 
    do i = 1,unstrM%csrf%lsf
       ii = idtmp(i)
       if (unstrM%csrf%ids(ii).ne.cck) then
          cck = unstrM%csrf%ids(ii)
          kk  = unstrM%csrf%inx(ii)
          check = 0
          fm = refs%Fmask(:,int(abs(kk)))
          !fm = refs%FtoV(int(abs(kk)),:)
          if (kk.lt.0) jj = 2
          if (kk.gt.0) jj = 1
          do j = 1,pin%s%Nfp
             l = unstrM%csrf%s2e(jj,ii)
             call findorder(l,unstrM%Clelist,m)
             n = unstrM%loc_t2v(fm(j),m)
             call findorder(n,unstrM%Cvlist,lct)
             if (unstrM%Cvpid(lct).eq.unstrM%rank) then
                check = 1              
             endif
          enddo 
          if (check.eq.1) then
             ll = ll + 1
             unstrM%esrf%s2e(:,ll) = unstrM%csrf%s2e(:,ii)
             unstrM%esrf%inx(ll)   = unstrM%csrf%inx(ii)
             unstrM%esrf%pid(ll)   = unstrM%csrf%pid(ii)
             unstrM%esrf%ids(ll)   = unstrM%csrf%ids(ii)
             unstrM%esrf%cen(:,ll) = unstrM%csrf%cen(:,ii)
             unstrM%esrf%drt(:,ll) = unstrM%csrf%drt(:,ii)
          endif
       else
          ! do nothing 
       endif
    enddo
 
    ! for unstrM%efmm

  end subroutine pnm_pre_fmm



end module geometry_mod
