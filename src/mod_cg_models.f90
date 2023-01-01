!*********************************************************************!
! This module reads models and papramters                            *!
!                                                                    *!
! By Jia Shi                                                         *!
!                                                                    *!
!*********************************************************************!

module cg_models_mod
   use MPI
   use datatype_mod
   use geometry_mod,                only: unstrM,refs
   use para_mod,                    only: rkind,pin
   use utility_mod 
   use cg_datatype_mod
   
   implicit none
 
   private
   public                               :: cg_load_models
   public                               :: written_nodebase
   public                               :: models
   
   type(modelcoeff), save               :: models
   real(kind=rkind)                     :: PI 
   
   !parameter(   PI = 3.14159265359)

 contains

   !-------------------------------------------------------------
   subroutine cg_load_models()
     ! read in the model parameters 
     integer                            :: i,j,k,ierr
     real(kind=rkind)                   :: cpmin,perwave,vstmp
     logical                            :: existence     
 
     !integer                            :: myrank,ierr
      
     ! set up some models parameters
     models%p_vp = 1; models%p_vs = 2; models%p_rho = 3
     models%nmps = 3; models%Gsiz  = pin%s%pNp*unstrM%Ntet
     models%siz = pin%s%pNp*unstrM%ClNele
     ! input frequency     
     !perwave    = 6.0D0

     ! reading the whole model
     allocate(models%coeff_loc(pin%s%pNp*unstrM%ClNele,models%nmps))
        
     models%coeff_loc(:,models%p_vp)   = 10.0D0
     models%coeff_loc(:,models%p_vs)   = 5.7735D0 !10.0D0/sqrt(3.0D0)
     models%coeff_loc(:,models%p_rho)  = 5.51D0
     
     call pnm_read_model(pin%fvpt,models%p_vp,pin%s%pNp)
     if (unstrM%rank == 0) print*, "read Vp file"
     call pnm_read_model(pin%fvst,models%p_vs,pin%s%pNp)
     if (unstrM%rank == 0) print*, "read Vs file"
     call pnm_read_model(pin%frhot,models%p_rho,pin%s%pNp)
     if (unstrM%rank == 0) print*, "read density file"
     !print*, maxval(models%coeff_loc(:,models%p_rho)), &
     !        minval(models%coeff_loc(:,models%p_rho)), unstrM%rank

     call separate_fs()

     !if (pin%selfG) then
     !   inquire(file=trim(pin%fgrav), exist=existence)
     !   if (existence) then 
     !      call fwd_read_gravaccel(pin%fgrav,pin%s%pNp)        
     !      if (unstrM%rank == 0) print*, "read gravitational acceleration g"
     !   else 
     !      if (unstrM%rank == 0) print*, "no gravitational acceleration file"
     !      ! TODO add fmm library into current code
     !      stop 
     !   endif 
     !endif

     if (pin%selfG) then
        call pnm_read_gravaccel(pin%fgrav,pin%s%pNp)
     endif

     ! prepare for density jumps at all surfaces
     if (pin%phi1) then
        call pnm_rho_jumps()
     endif

     call mpi_barrier(unstrM%comm,ierr)
     
   end subroutine cg_load_models

   subroutine pnm_rho_jumps()
      integer                             :: i,j,k,l,kk,ll,ierr
      integer                             :: lct,lctn
      integer, dimension(pin%s%Nfp)       :: f0,f1
      real(kind=rkind)                    :: rho0,rho1,vs0,vs1
      real(kind=rkind), allocatable       :: sd0(:),rd0(:)


      models%lsf = unstrM%srf%lsf
      allocate(models%srho(2,models%lsf)); models%srho = 0.0D0
      allocate(models%svs( 2,models%lsf)); models%svs = 0.0D0

      do i = 1,models%lsf
         k = unstrM%srf%s2e(1,i)
         f0 = refs%Fmask(:,unstrM%srf%inx(i))
         call findorder(k,unstrM%Clelist,lct)
         !rho0 = sum(models%coeff_loc((lct-1)*pin%s%pNp+1:lct*pin%s%pNp,models%p_rho))
         !rho0 = rho0/real(pin%s%pNp,8)
         rho0 = sum(models%coeff_loc((lct-1)*pin%s%pNp+f0,models%p_rho))
         rho0 = rho0/real(pin%s%Nfp,8)
         vs0 = sum(models%coeff_loc((lct-1)*pin%s%pNp+f0,models%p_vs))
         vs0 = vs0/real(pin%s%Nfp,8)
         l = unstrM%srf%s2e(2,i) 
         if (l.gt.0) then
            call findordernostop(l,unstrM%Clelist,lctn)
            if (lctn.gt.0) then
               f1 = 0
               !print*, unstrM%ClNeigh(:,lctn),k
               do j = 1,4
                  if (unstrM%ClNeigh(j,lctn).eq.k) f1 = refs%Fmask(:,j)
               enddo
               !print*,f1
               if (minval(f1).eq.0) then
                  print*,"Error: neigh", unstrM%rank
               endif
               !rho1 = sum(models%coeff_loc((lctn-1)*pin%s%pNp+1:lctn*pin%s%pNp,&
               !           models%p_rho))
               !vs1  = sum(models%coeff_loc((lctn-1)*pin%s%pNp+1:lctn*pin%s%pNp,&
               !           models%p_vs))
               rho1 = sum(models%coeff_loc((lctn-1)*pin%s%pNp+f1,models%p_rho))
               vs1  = sum(models%coeff_loc((lctn-1)*pin%s%pNp+f1,models%p_vs))
        
               rho1 = rho1/real(pin%s%Nfp,8)
               vs1 = vs1/real(pin%s%Nfp,8)
               
               !models%srho(i) = rho0 - rho1 
               models%srho(1,i) = rho0 
               models%srho(2,i) = rho1 
               models%svs( 1,i) =  vs0 
               models%svs( 2,i) =  vs1 
            else
               print*, "Debug: need to find density", unstrM%rank
               stop
               ! todo find density if it is needed
            endif
         else ! neigh -1
            models%srho(1,i) = rho0
            models%svs( 1,i) =  vs0
            !print*,rho0,i 
         endif
      enddo 
      !print*,models%srho(2,:)-models%srho(1,:)

      allocate(sd0(unstrM%srf%lsf),rd0(unstrM%esrf%lsf))

      !print*,maxval(models%srho),minval(models%srho),unstrM%rank
      allocate(models%erho(2,unstrM%esrf%lsf))
      allocate(models%evs( 2,unstrM%esrf%lsf))

      sd0 = models%srho(1,:)
      call pnm_find_sden(unstrM%srf,unstrM%esrf,sd0,rd0)
      !print*,maxval(rd0),unstrM%esrf%lsf,unstrM%rank
      models%erho(1,:) = rd0

      sd0 = models%srho(2,:)
      call pnm_find_sden(unstrM%srf,unstrM%esrf,sd0,rd0)
      models%erho(2,:) = rd0

      sd0 = models%svs(1,:)
      call pnm_find_sden(unstrM%srf,unstrM%esrf,sd0,rd0)
      models%evs(1,:) = rd0

      sd0 = models%svs(2,:)
      call pnm_find_sden(unstrM%srf,unstrM%esrf,sd0,rd0)
      models%evs(2,:) = rd0
   end subroutine pnm_rho_jumps
 

   subroutine pnm_find_sden(sf0,sf1,dt0,dt1)
      type(surface), intent(inout)                :: sf0,sf1
      real(kind=rkind), allocatable               :: dt0(:),dt1(:)

      integer                                     :: i,j,k,l,dm,ierr,lct,check
      integer                                     :: lvts0,lvts1,lct0,lct1,lct2
      integer, allocatable, dimension(:)          :: idtmp,dat0,dat1,npid0,npid1
      integer, allocatable, dimension(:)          :: cnt0,cnt1,tmpd0,tmpd1,vts0,vts1
      real(kind=rkind), allocatable, dimension(:) :: rdat0,rdat1

      allocate(npid0(unstrM%nproc),npid1(unstrM%nproc))
      npid0 = 0; npid1 = 0

      allocate(tmpd0(sf1%lsf)); tmpd0 = sf1%pid
      call simplessort(tmpd0,idtmp)
      dm = 1
      allocate(vts0(dm*sf1%lsf))
      vts0 = sf1%ids(idtmp)
 
      do i = 1,sf1%lsf
         j = tmpd0(i) + 1; npid0(j) = npid0(j) + 1         
      enddo   
   
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
      allocate(rdat0(lvts0/dm),rdat1(lvts1/dm))

      do i = 1,lvts1/dm
         lct0 = vts1(i) - sf0%dist(unstrM%rank+1)
         rdat1(i) = dt0(lct0) !models%srho(lct0) 
      enddo
      cnt0  = cnt0/dm;  cnt1  = cnt1/dm
      npid0 = npid0/dm; npid1 = npid1/dm

      !print*,maxval(rdat1),unstrM%rank

      call MPI_ALLTOALLV(rdat1,npid1,cnt1,mpi_real8,&
                         rdat0,npid0,cnt0,mpi_real8,unstrM%comm,ierr)
     
      !print*,lvts0,sf1%lsf,unstrM%rank
      !allocate(models%erho(lvts0/dm))
      !models%erho(idtmp) = rdat0
      dt1(idtmp) = rdat0

      !print*,maxval(rdat0),minval(rdat0),unstrM%rank

   end subroutine pnm_find_sden
 
   subroutine separate_fs()
     integer                              :: i,j,k,ierr
     real(kind=rkind)                     :: rtmp    
     real(kind=rkind), allocatable        :: tetstat(:)

     unstrM%f%ClNele = 0
     unstrM%s%ClNele = 0

     allocate(tetstat(unstrM%ClNele)); tetstat = 1 
     do i = 1,unstrM%ClNele
        rtmp = sum(models%coeff_loc((i-1)*pin%s%pNp+1:&
                                        i*pin%s%pNp,models%p_vs))
        if (rtmp.le.pin%TOL) then
           ! fluid
           tetstat(i) = 0
           unstrM%f%ClNele = unstrM%f%ClNele + 1
        else
           unstrM%s%ClNele = unstrM%s%ClNele + 1
        endif
     enddo
     call mpi_barrier(unstrM%comm,ierr)

     print*,'# of s',unstrM%s%ClNele,'# of f',unstrM%f%ClNele,&
            'at rank', int(unstrM%rank) 
     !allocate(unstrM%s%loc_nods(3,4*unstrM%ClNele))

     !if (unstrM%s%ClNele.gt.0) then
     !   !allocate(unstrM%s%loc_label(4*unstrM%s%ClNele))
     !   allocate(unstrM%s%loc_nods(3,4*unstrM%s%ClNele))
     !   j = 1
     !   do i = 1,unstrM%ClNele
     !      if (tetstat(i).eq.1) then
     !         !unstrM%s%loc_label((j-1)*4+1:j*4) = &
     !         !unstrM%lt2vid((i-1)*4+1:i*4) 
     !         unstrM%s%loc_nods(:,(j-1)*4+1:j*4) = &
     !         unstrM%Cv_crs(:,unstrM%lt2vid((i-1)*4+1:i*4))
     !         j = j + 1 
     !      endif
     !   enddo
     !endif


   end subroutine separate_fs

   subroutine written_nodebase(filename,buf,siz)
     integer, intent(in)                  :: siz
     character(len=1024), intent(in)      :: filename
     real(kind=rkind), intent(in)         :: buf(siz) 

     integer                              :: rec_len,i,fid
    
     ! to make the code more portable 
     fid = 2015
     inquire(iolength=rec_len) rec_len
     open( unit = fid, file = filename, form = "unformatted", &
           status = "replace", access = "direct", recl = rec_len*int(rkind/2))
   
     do i = 1,siz
        write(fid,rec=i) buf(i)
     enddo
     
     close(fid)
 
   end subroutine written_nodebase    


   subroutine pnm_read_model(fname,ptr,pNp)
      character(len=1024), intent(in)         :: fname
      integer, intent(in)                     :: ptr,pNp

      integer                                 :: i,j,k,l
      integer                                 :: fid,ierr,error,nints
      ! mpi I/O 
      integer                                 :: stat(mpi_status_size)
      integer(kind=mpi_offset_kind)           :: offset
      integer*8                               :: leight,lfour
      real*8, allocatable                     :: buf(:),bufn(:,:)

      real*8, allocatable, dimension(:)       :: rdt0,rdt1
   
      logical                                 :: alive
   
      ! check 
      if (unstrM%rank == 0) then
         inquire(file=trim(fname),exist=alive)
         if (.not.alive) then
            print*,'Error: File "',trim(fname),'" does not exist.'
            stop
         endif
      endif 


      leight = 8; lfour = 4

      call mpi_file_open(unstrM%comm,fname,mpi_mode_rdonly,&
           mpi_info_null,fid,ierr)

      nints = unstrM%org%nele*pNp 

      offset = unstrM%org%edist(unstrM%rank+1)*pNp*leight

      allocate(buf(nints)); buf = 0.0D0
       
      call mpi_file_read_at(fid,offset,buf,nints,mpi_real8,stat,ierr)
      !print*, maxval(buf),minval(buf),unstrM%rank
 
      call mpi_file_close(fid,ierr)
      
      allocate(rdt0(unstrM%org%nele),rdt1(unstrM%ClNele))
      allocate(bufn(pNp,unstrM%ClNele))
      do j = 1,pNp
         do i = 1,unstrM%org%nele
            rdt0(i) = buf((i-1)*pNp+j)
         enddo   

         call pnm_find_data(unstrM%org%nele,unstrM%org%elist,rdt0,&
               unstrM%ClNele,unstrM%CNOids,unstrM%Clelist,rdt1,unstrM%comm)

         bufn(j,:) = rdt1
      enddo
       
      do i = 1,unstrM%ClNele; do j = 1,pNp
         k = (i-1)*pNp+j
         models%coeff_loc(k,ptr) = bufn(j,i)
         !models%coeff_loc(k,ptr) = sum(bufn(:,i))/real(pNp,8)
      enddo; enddo

      deallocate(buf,rdt0,rdt1,bufn)

   end subroutine pnm_read_model

   subroutine pnm_read_gravaccel(fname,pNp)
      character(len=1024), intent(in)         :: fname
      integer, intent(in)                     :: pNp

      integer                                 :: i,j,k,l,nd
      integer                                 :: fid,ierr,error,nints
      ! mpi I/O 
      integer                                 :: stat(mpi_status_size)
      integer(kind=mpi_offset_kind)           :: offset
      integer*8                               :: leight,lfour
      real*8, allocatable                     :: buf(:),bufn(:,:)

      real*8, allocatable, dimension(:)       :: rdt0,rdt1

      leight = 8; lfour = 4

      allocate(models%g0(pNp,3,unstrM%ClNele)); models%g0 = 0.0D0;
      nints = unstrM%org%nele*pNp*3 

      allocate(buf(nints)); buf = 0.0D0
      allocate(rdt0(unstrM%org%nele),rdt1(unstrM%ClNele))
      !allocate(bufn(pNp,unstrM%ClNele))

      call mpi_file_open(unstrM%comm,fname,mpi_mode_rdonly,&
           mpi_info_null,fid,ierr)

      offset = unstrM%org%edist(unstrM%rank+1)*pNp*3*leight
      !print*,offset,nints,unstrM%rank
      call mpi_file_read_at(fid,offset,buf,nints,mpi_real8,stat,ierr)
      !call mpi_file_seek(fid,offset,mpi_seek_set,ierr)
      !call mpi_file_read(fid,buf,nints,mpi_real8,stat,ierr)
      !print*,maxval(buf),minval(buf),unstrM%rank

      do j = 1,pNp; do nd = 1,3
         do i = 1,unstrM%org%nele
            rdt0(i) = buf((i-1)*pNp*3+(j-1)*3+nd)
         enddo   
         call pnm_find_data(unstrM%org%nele,unstrM%org%elist,rdt0,&
               unstrM%ClNele,unstrM%CNOids,unstrM%Clelist,rdt1,unstrM%comm)
         !bufn(j,:) = rdt1
         models%g0(j,nd,:) = rdt1
      enddo; enddo
      
      call mpi_file_close(fid,ierr)
      !print*, maxval(models%g0),minval(models%g0),unstrM%rank

   end subroutine pnm_read_gravaccel


   subroutine fwd_read_model(fname,ptr,pNp)
      character(len=1024), intent(in)         :: fname
      integer, intent(in)                     :: ptr,pNp
      
      integer                                 :: i,j,k
      integer                                 :: fid,rec_len
      real(kind=8), allocatable               :: buf(:)        
      

      allocate(buf(unstrM%ClNele*pNp))
      
      fid = 10

      inquire(iolength = rec_len) rec_len
      open(unit = fid, file = fname, form = "unformatted", status = "old",&
           access = "direct", recl = rec_len*2)
      !do i = 1,models%siz
      !   read(10, rec=i), buf(i)
      !enddo

      do i = 1,unstrM%ClNele 
         k = (unstrM%Clelist(i)-1)*pNp
         do j = 1,pNp
            read(fid, rec=k+j) buf((i-1)*pNp+j)
         enddo
      enddo 

      close(fid)     

      models%coeff_loc(:,ptr) = buf

      !print*, maxval(buf), minval(buf)      

   end subroutine fwd_read_model


   !TODO in the future generate it by itself
   subroutine fwd_read_gravaccel(fname,pNp)
      character(len=1024), intent(in)         :: fname
      integer, intent(in)                     :: pNp

      integer                                 :: i,j,k,l
      integer                                 :: fid,rec_len
      real(kind=8)                            :: buf

      allocate(models%g0(pNp,3,unstrM%ClNele)); models%g0 = 0.0D0;
     
      !allocate(buf(3,unstrM%nnodes))
      fid = 10
      inquire(iolength = rec_len) rec_len
      open(unit = fid, file = fname, form = "unformatted", status = "old",&
           access = "direct", recl = rec_len*2)

      !do i = 1,3  
      !   do j = 1,unstrM%nnodes
      !      read(10, rec=(i-1)*unstrM%nnodes+j) buf(i,j) 
      !   enddo
      !enddo
   
      do j = 1,3
         do k = 1,unstrM%ClNele
            do i = 1,pNp
               l = unstrM%Clelist(k)
               ! compute node id
               read(fid, rec=(j-1)*unstrM%Ntet*pNp+(l-1)*pNp+i) buf 
               models%g0(i,j,k) = buf
            enddo
         enddo
      enddo
    
      close(fid) 
      !l = 50; print*, models%g0(:,3,l), unstrM%Clelist(l)
  
   end subroutine fwd_read_gravaccel   







end module cg_models_mod
