!*********************************************************************!
!* This is the module how to build up a CG matrix in parallel         !
!* By Jia Shi                                                         !
!*********************************************************************!
module cg_create_matrix_mod
  use MPI
  use geometry_mod,                      only: unstrM,ednew,refs,reff
  use para_mod,                          only: rkind,pin
  use cg_models_mod,                     only: models
  use utility_mod                 
  use cg_datatype_mod
 
  implicit none
  ! local CG matrix 
  
  private
  public                                    :: cg_create_matrix
  public                                    :: CGM
  public                                    :: PI,LONG1,RT

  type (CGMatrix), save                     :: CGM

  integer*8                                 :: LONG1 
  real(kind=rkind)                          :: PI,eps0 
  complex(kind=rkind)                       :: RT 

  parameter(   RT = (0.0D0,1.0D0))
  parameter(   PI = 3.14159265359)
  parameter(LONG1 = 1)
  parameter(eps0 = 1.0D-15)
  !purpose of LONG1 is to make it more robust when we compute long integer

 contains
!-------------------------------------------------------------------
  subroutine cg_create_matrix()
    ! construct CG matrix
    ! first we build CG COO matrix structure
    ! form the DG part of matrix to CG structure 
    integer                                 :: ierr

    if (unstrM%fsexist.or.unstrM%purefluid) then
       if (.not. CGM%stat) call matrixstruct_general()
       if (unstrM%rank == 0 .and. .not. CGM%stat) then
          print*, "the matrix structure (with fluid or fluid-solid) is done"
       endif
       call writelog("the matrix (f or f-s) structure is done",pin%comm)
      
       call CGFSE3D_ISO()
    else
       if (.not. CGM%stat) call matrixstruct()
       if (unstrM%rank == 0 .and. .not. CGM%stat) then
          print*, "the matrix structure (no fluid-solid) is done"
       endif
       call writelog("the matrix structure is done",pin%comm)
       call report_time(pin%comm)
 
       call CGE3D_ISO()

    endif

  end subroutine cg_create_matrix


  subroutine cg_mat_write()
    character(len=1024)                                 :: str_p,str_d,fname

    integer                                             :: i,j,k,l,fid
    logical                                             :: alive

    write(str_p,*) pin%s%pOrder; str_p = trim(adjustl(str_p))

    fname = trim(pin%outputdir)//trim(pin%basename)//'_pod_'&
            //trim(str_p)//'_Ap'//'.txt'
    fname = trim(fname)
    print*, trim(fname)
    fid = 2333
    inquire(file=trim(fname),exist=alive)
    if (alive) then
       open(fid,file=trim(fname),status='old',action='write')
    else
       open(fid,file=trim(fname),status='new',action='write')
    endif

    !write(fid,2000) CGM%Ap%Gsiz,CGM%Ap%GNNZ
    do i = 1,CGM%Ap%siz
       l = CGM%Ap%rowdist(i+1)-CGM%Ap%rowdist(i) 
       k = CGM%Ap%rowdist(i)
       do j = 1,l       
          write(fid,2001) i,CGM%Ap%col(k+j),CGM%Ap%val(k+j)
       enddo
    enddo
    close(fid)


2000 format (2I16)
2001 format (2I16,1D20.10)

  end subroutine cg_mat_write



  !--------------------------------------------------------------------------
  subroutine CGFSE3D_ISO()
    integer, allocatable                                :: trow(:),srow(:)
    real(kind=rkind), allocatable                       :: TM(:,:), MM(:,:)
    
    ! add for fluid 
    real(kind=rkind), allocatable                       :: FT(:,:),FM(:,:)
    real(kind=rkind), allocatable                       :: SCM(:,:)
    
    integer                                             :: i,j,k,l,m,n,p,q
    integer                                             :: i1,j1,cout,fcid
    integer                                             :: vid,vid0,vid1,nn,m0
    integer                                             :: dm,les,lef,colid
    integer, dimension(pin%s%pNp)                       :: rows1,cols1
    integer, allocatable                                :: ntmp(:)    
 
    real(kind=rkind), allocatable, dimension(:)         :: C33,C44
    real(kind=rkind), dimension(pin%s%pNp)              :: lam,mu,rhot,&
                                                           rinl,lamiv,rhod
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: OPt,OP1,OP2
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: Derv,TDv,Dv2
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPij,OPs,OPtmp
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPtrs,Mrho
    real(kind=rkind)                                    :: lamavg,muavg,rhoavg

    real(kind=rkind), dimension(pin%f%pNp,pin%f%pNp,3)  :: fOPt,fOP1,fOP2
    real(kind=rkind), dimension(pin%f%pNp,pin%f%pNp,3)  :: fDerv,fTDv,fDv2
    real(kind=rkind), dimension(pin%f%pNp,pin%f%pNp)    :: fOPij,fOPs,fOPtmp
    real(kind=rkind), dimension(pin%f%pNp,pin%f%pNp)    :: fOPtrs,fMrho,flam

    ! for reference gravity
    integer, dimension(pin%s%Nfp)                       :: Fms
    real(kind=rkind)                                    :: N2avg,rhof
    real(kind=rkind), dimension(pin%s%pNp,3)            :: gk1,gk0,gradrho,normalg
    real(kind=rkind), dimension(pin%s%pNp)              :: normg,N2k1,gtmpij
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: rhomtx,rhotmp,&
                                                           dgtmp,dgtmp0,dgtmp1
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: diaggk1,diagnormalg
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPrhodi,OPrhodj,&
                                                           OPrhodij,OPrho,&
                                                           OPrhodit,OPrhodjt
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPrho1,OPrho1t,OPrho2,OPrho2t
    real(kind=rkind), dimension(pin%s%Nfp,pin%s%Nfp)    :: smassi,smassit,smassj,&
                                                           smassjt,smassij,digtmp
    ! for a more accurate derivatives
    real(kind=rkind), dimension(pin%s%pNp,3)            :: dism 
    real(kind=rkind), dimension(3,3)                    :: ndism,invdism,dgk1
    real(kind=rkind), dimension(3)                      :: drho0
   
    ! for fluid surface integral
    real(kind=rkind)                                    :: surfrho,sgn
    real(kind=rkind), dimension(pin%f%Nfp)              :: surfgn
    real(kind=rkind), dimension(pin%f%Nfp,pin%f%Nfp)    :: surfp,difg     


    dm   = 3; les = pin%s%pNp*dm; 
    allocate(trow(pin%s%pNp),TM(les,les),MM(les,les))
    
    ! TODO more careful about mix element method 
    lef = pin%f%pNp+les
    allocate(FT(lef,lef),FM(lef,lef)) 
    l = pin%s%pNp*3
    allocate(SCM(l,l))


    allocate(C44(models%siz),C33(models%siz))
    C44  = models%coeff_loc(:,models%p_rho)*models%coeff_loc(:,models%p_vs)**2
    C33  = models%coeff_loc(:,models%p_rho)*models%coeff_loc(:,models%p_vp)**2 - 2*C44

    do k = 1,unstrM%ClNele
       ! build up global ids 
       do l = 1,pin%s%pNp 
          i = (k-1)*pin%s%pNp + l
          !j = (unstrM%Clelist(k)-1)*pin%s%pNp+l
          rows1(l) = i ! for coeff
          cols1(l) = i ! for global ids
       enddo
      
       ! take the avarage 
       !C33(rows1) = sum(C33(rows1))/real(pin%s%pNp,8)     
       !C44(rows1) = sum(C44(rows1))/real(pin%s%pNp,8)     

       lam   = C33(rows1) 
       mu    = C44(rows1) 
       !rhot  = sum(models%coeff_loc(rows1,models%p_rho))/real(pin%s%pNp,8) 
       rhot  = models%coeff_loc(rows1,models%p_rho) 

       rhoavg = sum(rhot)/real(pin%s%pNp,8)
       lamavg = sum(lam)/real(pin%s%pNp,8)
       ! for reference gravity
       if (pin%selfG) then
          ! get g infomation
          gk1 = models%g0(:,:,k)/1.D3
          gk0 = gk1       
          ! construct for a accurate derivatives
          do i = 1,3
             dism(:,i) = unstrM%loc_nods(i,(k-1)*pin%s%pNp+1:k*pin%s%pNp) -& 
                     sum(unstrM%loc_nods(i,(k-1)*pin%s%pNp+1:k*pin%s%pNp))/real(pin%s%pNp,8)
             gk0(:,i) = gk1(:,i) - sum(gk1(:,i))/real(pin%s%pNp,8)
          enddo
          !print*,dism(1,1),k

          ndism = matmul(transpose(dism),dism) 
          call matinv(ndism,invdism,3)
          !print*,ndism(1,1),k
 
          !dgk1 = matmul(transpose(dism),gk1)      
          !dgk1 = matmul(invdism,dgk1)
          ! change JS 03/092018
          dgk1 = matmul(transpose(dism),gk0)      
          dgk1 = matmul(invdism,dgk1)

          
          ! change JS 03/092018
          rhod = rhot - rhoavg
          drho0 = matmul(transpose(dism),rhod)
          drho0 = matmul(invdism,drho0)

          !print*, '1',rhod,rhoavg
          !print*, '2',drho0
          ! derivatives of density     
          gradrho = 0.0; diaggk1 = 0.0 
          do i = 1,3
             gradrho(:,i) = drho0(i)
             call diagmatrix(diaggk1(:,:,i),gk1(:,i),pin%s%pNp)
          enddo
          !print*,maxval(diaggk1),minval(diaggk1)
          ! normalize g 
          normg = 0.0D0
          do i = 1,3      
             normg = normg + gk1(:,i)**2
          enddo
          normg = dsqrt(normg)
          !print*, normg, k

          ! build N^2 and normal g
          N2k1 = 0.0D0; diagnormalg = 0.0D0 
          do i = 1,3
             N2k1 = N2k1 + gradrho(:,i)*gk1(:,i)
             normalg(:,i) = gk1(:,i)/max(normg,eps0)
             call diagmatrix(diagnormalg(:,:,i),normalg(:,i),pin%s%pNp)
          enddo
          !print*,normalg,k

          N2k1 = N2k1/rhot - normg**2/lam*rhot
          !N2k1 = N2k1/rhot - normg**2/lamavg*rhoavg

          ! JS 05252020 update
          do i = 1,pin%s%pNp 
             if (normg(i).lt.eps0) then 
                N2k1(i) = 0.0D0
             endif 
          enddo
          ! add JS 05512020 
          ! for pure fluid planet, set N2 to zero
          if (unstrM%purefluid) then 
             N2k1 = 0.0D0
          endif  
          N2avg = sum(N2k1)/real(pin%s%pNp,8)    
 
          !N2k1 = N2avg 
          !print*,N2k1 
          !N2avg = 0.0D0 
          !print*, N2avg,k,unstrM%rank 
          
       endif

       ! JS 05252020 an important update
       ! take the avarage 
       !C33(rows1) = sum(C33(rows1))/real(pin%s%pNp,8)     
       !C44(rows1) = sum(C44(rows1))/real(pin%s%pNp,8)     
       !lam   = C33(rows1) 
       !mu    = C44(rows1) 
       !rhot  = sum(models%coeff_loc(rows1,models%p_rho))/real(pin%s%pNp,8) 

       
       ! for the solid region 
       if (maxval(mu) >= pin%TOL) then
          ! compute the spatial derivatives  
          Derv(:,:,1) = unstrM%loc_invJ(1,1,k)*refs%Drst(:,:,1) +&
                        unstrM%loc_invJ(1,2,k)*refs%Drst(:,:,2) +&
                        unstrM%loc_invJ(1,3,k)*refs%Drst(:,:,3)
          
          Derv(:,:,2) = unstrM%loc_invJ(2,1,k)*refs%Drst(:,:,1) +&
                        unstrM%loc_invJ(2,2,k)*refs%Drst(:,:,2) +&
                        unstrM%loc_invJ(2,3,k)*refs%Drst(:,:,3)
          
          Derv(:,:,3) = unstrM%loc_invJ(3,1,k)*refs%Drst(:,:,1) +&
                        unstrM%loc_invJ(3,2,k)*refs%Drst(:,:,2) +&
                        unstrM%loc_invJ(3,3,k)*refs%Drst(:,:,3)
      
          Mrho  = refs%MassM*sum(rhot)/real(pin%s%pNp)
          !call diagmatrix(Mrho,rhot,pin%s%pNp)
          !Mrho = matmul(Mrho,refs%MassM)+matmul(refs%MassM,Mrho)
          !Mrho = Mrho/2.0D0

 
          !call diagmatrix(Mrho,dsqrt(rhot),pin%s%pNp)
          !Mrho = matmul((matmul(Mrho,refs%MassM)),Mrho)      
 
          OPs  = 0.0D0

          ! save derivative matrix 
          Dv2  = Derv
          ! build up derivative matricies
          do i = 1,dm
             ! lambda*D_i
             !call realmcolupdate(Derv(:,:,i),lam,OP1(:,:,i),pin%s%pNp)
             call diagmatrix(dgtmp,lam,pin%s%pNp)
             dgtmp = matmul(refs%MassM,dgtmp)
             dgtmp = (dgtmp+transpose(dgtmp))/2.0D0
             OP1(:,:,i) = matmul(dgtmp,Derv(:,:,i))
             ! mu*D_i
             !call realmcolupdate(Derv(:,:,i),mu,OP2(:,:,i),pin%s%pNp) 
             call diagmatrix(dgtmp,mu,pin%s%pNp)
             dgtmp = matmul(refs%MassM,dgtmp)
             dgtmp = (dgtmp+transpose(dgtmp))/2.0D0
             OP2(:,:,i) = matmul(dgtmp,Derv(:,:,i))
 
             ! D_i'*mu*D_i
             !Derv(:,:,i) = matmul(refs%MassM,Dv2(:,:,i))
             !Dv2(:,:,i)  = Derv(:,:,i)
             TDv(:,:,i)  = transpose(Derv(:,:,i))
             OPt(:,:,i)  = matmul(TDv(:,:,i),OP2(:,:,i))
             OPs         = OPs + OPt(:,:,i)
          enddo


          do i = 1,dm; do j = 1,dm
             if (pin%selfG) then
                ! prepare for g related terms  
                ! for OPrho
                OPrho = 0.0D0
                !call diagmatrix(rhotmp,gk1(:,j),pNp)
                rhomtx   = matmul(refs%MassM,diaggk1(:,:,j))
                OPrho1   = matmul(transpose(Dv2(:,:,i)),rhomtx)
       
                rhomtx   = matmul(refs%MassM,Dv2(:,:,i))
                OPrho1t  = matmul(diaggk1(:,:,j),rhomtx)

                ! symmetrize it 
                !call diagmatrix(rhotmp,rhot*gk1(:,i),pNp)
                rhomtx   = matmul(refs%MassM,Dv2(:,:,j))
                OPrho2   = matmul(diaggk1(:,:,i),rhomtx)
                

                rhomtx   = matmul(refs%MassM,diaggk1(:,:,i))
                OPrho2t  = matmul(transpose(Dv2(:,:,j)),rhomtx) 

                OPrho    = (OPrho1+transpose(OPrho1t))/4.0D0 +&
                           (OPrho2+transpose(OPrho2t))/4.0D0

                rhomtx   = refs%MassM*(dgk1(i,j)+dgk1(j,i))/2.0D0
                OPrho    = OPrho - rhomtx 
                !if (mod(k,1000)==1.and. i==3.and. j==1) print*,'2',maxval(rhomtx)

                rhotmp  = matmul(diaggk1(:,:,i),refs%MassM)
                OPrho1  = matmul(transpose(Dv2(:,:,j)),rhotmp) 
                rhotmp  = matmul(diaggk1(:,:,i),Dv2(:,:,j))
                OPrho1t = matmul(refs%MassM,rhotmp)
                OPrho1  = (OPrho1 + transpose(OPrho1t))/2.D0 

                ! symmetrize it              
                rhotmp  = matmul(diaggk1(:,:,j),Dv2(:,:,i))
                OPrho2  = matmul(refs%MassM,rhotmp)
                rhotmp  = matmul(diaggk1(:,:,j),refs%MassM)
                OPrho2t = matmul(transpose(Dv2(:,:,i)),rhotmp) 
                OPrho2  = (OPrho2 + transpose(OPrho2t))/2.0D0
                
                OPrho   = OPrho - (OPrho1+OPrho2)/2.0D0

                !OPrho   = OPrho*rhoavg
                call diagmatrix(dgtmp,rhot,pin%s%pNp)
                dgtmp   = matmul(OPrho,dgtmp)+matmul(dgtmp,OPrho)
                OPrho   = dgtmp/2.0D0
             endif


             if (i == j) then
                ! diagonal component
                OPtmp         = matmul(TDv(:,:,i),OP1(:,:,i))
                OPij          = OPs + OPt(:,:,i) + OPtmp                             
                
                ! make it exactly symetric
                OPij          = (OPij + transpose(OPij))/2.0D0 

                if(pin%selfG) OPij = OPij + OPrho 

                TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
                MM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = Mrho*unstrM%loc_detJ(k)
             else
                ! off diagonal component
                OPtmp         = matmul(TDv(:,:,i),OP1(:,:,j)) + &
                                matmul(TDv(:,:,j),OP2(:,:,i))
                OPtrs         = matmul(TDv(:,:,j),OP1(:,:,i)) + &
                                matmul(TDv(:,:,i),OP2(:,:,j))
                OPij          = (transpose(OPtrs) + OPtmp)/2.0D0
                
                if(pin%selfG) OPij = OPij + OPrho 
                
                TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
                MM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = 0.0D0
             endif
          enddo; enddo

          ! add the values into the matrices
          do m = 1,pin%s%pNp
             l   = unstrM%lt2vid((k-1)*pin%s%pNp+m)
             if (unstrM%Cvpid(l).eq.unstrM%rank) then
                ! vtx id
                i = unstrM%loc_t2v(m,k)  
                call findorder(i,unstrM%new%vlist,vid)
                !print*,vid
                do p = 1,3
                   ! for Ad 
                   j = sum(CGM%Ad%rnum(1:vid)) - CGM%Ad%rnum(vid) + p
                   les = CGM%Ad%rowdist(j+1)-CGM%Ad%rowdist(j)
                   allocate(ntmp(les))
                   ntmp = CGM%Ad%col(CGM%Ad%rowdist(j)+1:CGM%Ad%rowdist(j+1))
                   do n = 1,pin%s%pNp
                      vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                      do q = 1,3
                         vid1 = CGM%vstt(vid0) + q
                         call findorder(vid1,ntmp,colid)
                         CGM%Ad%val(CGM%Ad%rowdist(j)+colid) = &
                         CGM%Ad%val(CGM%Ad%rowdist(j)+colid) + &
                         TM((p-1)*pin%s%pNp+m,(q-1)*pin%s%pNp+n)  
                      enddo
                   enddo
                   deallocate(ntmp)
                   ! for B matrix
                   j = sum(CGM%B%rnum(1:vid)) - CGM%B%rnum(vid) + p
                   les = CGM%B%rowdist(j+1)-CGM%B%rowdist(j)
                   allocate(ntmp(les))
                   ntmp = CGM%B%col(CGM%B%rowdist(j)+1:CGM%B%rowdist(j+1))
                   do n = 1,pin%s%pNp
                      vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                      vid1 = CGM%vstt(vid0) + p
                      call findorder(vid1,ntmp,colid)
                      CGM%B%val(CGM%B%rowdist(j)+colid) = &
                      CGM%B%val(CGM%B%rowdist(j)+colid) + &
                      MM((p-1)*pin%s%pNp+m,(p-1)*pin%s%pNp+n)  
                   enddo
                   deallocate(ntmp)
                enddo
             endif
          enddo 
       else
          ! for the fluid region
          ! compute the spatial derivatives  
          fDerv(:,:,1) = unstrM%loc_invJ(1,1,k)*reff%Drst(:,:,1) +&
                         unstrM%loc_invJ(1,2,k)*reff%Drst(:,:,2) +&
                         unstrM%loc_invJ(1,3,k)*reff%Drst(:,:,3)
          
          fDerv(:,:,2) = unstrM%loc_invJ(2,1,k)*reff%Drst(:,:,1) +&
                         unstrM%loc_invJ(2,2,k)*reff%Drst(:,:,2) +&
                         unstrM%loc_invJ(2,3,k)*reff%Drst(:,:,3)
          
          fDerv(:,:,3) = unstrM%loc_invJ(3,1,k)*reff%Drst(:,:,1) +&
                         unstrM%loc_invJ(3,2,k)*reff%Drst(:,:,2) +&
                         unstrM%loc_invJ(3,3,k)*reff%Drst(:,:,3)
          
        
          !fMrho  = reff%MassM*sum(rhot)/real(pin%s%pNp,8)
          lamiv = 1.0D0/dsqrt(lam)
          call diagmatrix(flam,lamiv,pin%s%pNp)
          flam = matmul(matmul(flam,refs%MassM),flam)         
 
          !lamiv = 1.0D0/(sum(lam)/real(pin%s%pNp,8))
          !call diagmatrix(flam,lamiv,pin%s%pNp)
          !flam = matmul(flam,refs%MassM)+matmul(refs%MassM,flam)
          !flam = flam/2.0D0

          fOPs  = 0.0D0
          Mrho  = refs%MassM*sum(rhot)/real(pin%s%pNp,8)
          !call diagmatrix(Mrho,rhot,pin%s%pNp)
          !Mrho = matmul(Mrho,refs%MassM)+matmul(refs%MassM,Mrho)
          !Mrho = Mrho/2.0D0
       
          !call diagmatrix(Mrho,dsqrt(rhot),pin%s%pNp)
          !Mrho = matmul(matmul(Mrho,refs%MassM),Mrho)      

          ! save derivative matrix 
          fDv2  = fDerv

          FT = 0.0D0; FM = 0.0D0
          if (pin%selfG) then
             !print*, N2avg,k
             do i = 1,3; do j = 1,3
                if (i .eq. j) then
                   ! N**2 part
                   gtmpij = normalg(:,j)*N2k1*rhot!/N2avg
                   !gtmpij = normalg(:,j)!*rhot!*N2K1!/N2avg
                   call diagmatrix(rhotmp,gtmpij,pin%s%pNp) 
                   rhomtx  = matmul(refs%MassM,rhotmp) 
                   OPij    = matmul(diagnormalg(:,:,i),rhomtx) 

                   !print*,maxval(OPij)

                   rhomtx  = matmul(refs%MassM,diagnormalg(:,:,i))
                   OPtmp   = matmul(rhotmp,rhomtx) 

                   !print*,maxval(OPtmp)

                   !OPij    = (OPij + transpose(OPtmp))/2.0*rhoavg
                   OPij    = (OPij + OPtmp)/2.0!*N2avg!*rhoavg
                   !OPij = 0.0 
                   !print*,maxval(OPij)
                   
                   ! mass part & diagonal part
                   FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                      (i-1)*pin%s%pNp+1:i*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
                   FM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                      (i-1)*pin%s%pNp+1:i*pin%s%pNp) = Mrho*unstrM%loc_detJ(k)

                   fDerv(:,:,i) = matmul(reff%MassM,fDv2(:,:,i))
                   fTDv(:,:,i)  = transpose(fDv2(:,:,i))
                   fOPt(:,:,i)  = matmul(fTDv(:,:,i),reff%MassM) 
                   ! make it symmetric
                   fOPij        = (fDerv(:,:,i) + transpose(fOPt(:,:,i)))/2.0D0
                   ! JS 0413/2018
                   ! additional terms
                   !OPrhodi  = matmul(refs%MassM,diaggk1(:,:,i)) 
                   !OPrhodit = matmul(diaggk1(:,:,i),refs%MassM) 
                   !OPrhodi  = (OPrhodi + transpose(OPrhodit))/2.0D0
                   !OPrhodi = 0.0
                   !print*,'rho', maxval(OPrhodi),minval(OPrhodi)
                   !print*, 'deriv', maxval(fOPij), minval(fOPij)
                   !OPrhodi  = 0.0D0 
                   ! \partial_i p
                   rinl = rhot/lam
                   !print*,maxval(rinl),maxval(rhot),maxval(lam)
                   !print*,minval(rinl),minval(rhot),minval(lam)

                   dgtmp0 = 0.0D0; dgtmp = 0.0D0
                   call diagmatrix(dgtmp0,rinl,pin%s%pNp)
                   ! add
                   dgtmp1 = diaggk1(:,:,i)
                   dgtmp = matmul(dgtmp0,dgtmp1)
                   !print*,'1',maxval(dgtmp),maxval(dgtmp1),maxval(dgtmp0)
                   !print*,'2',minval(dgtmp),minval(dgtmp1),minval(dgtmp0)
                   
                   !OPrhodi = matmul(dgtmp,OPrhodi)+matmul(OPrhodi,dgtmp)
                   !OPrhodi = OPrhodi/2.0D0 
                   OPrhodi  = matmul(refs%MassM,dgtmp) 
                   OPrhodit = matmul(dgtmp,refs%MassM) 
                   OPrhodi  = (OPrhodi + transpose(OPrhodit))/2.0D0
                

                   !FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   !       3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = &
                   !   (fOPij-OPrhodi*rhoavg/lamavg)*unstrM%loc_detJ(k)
                   FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                          3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = &
                      (fOPij-OPrhodi)*unstrM%loc_detJ(k)
                   FM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                          3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = 0.0D0
                   ! \div u
                   !FT(3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp,&
                   !  (i-1)*pin%s%pNp+1:i*pin%s%pNp) = &
                   !  transpose(fOPij-OPrhodi*rhoavg/lamavg)*unstrM%loc_detJ(k)
                   FT(3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp,&
                     (i-1)*pin%s%pNp+1:i*pin%s%pNp) = &
                     transpose(fOPij-OPrhodi)*unstrM%loc_detJ(k)
                   FM(3*pin%f%pNp+1:3*pin%s%pNp+pin%f%pNp,&
                     (i-1)*pin%s%pNp+1:i*pin%s%pNp) = 0.0D0

                   !print*,maxval(OPrhodi)
                else
                   ! N**2 part
                   gtmpij = normalg(:,j)*N2k1*rhot!/N2avg
                   !gtmpij = normalg(:,j)!*rhot*N2K1!/N2avg
                   call diagmatrix(rhotmp,gtmpij,pin%s%pNp) 
                   rhomtx  = matmul(refs%MassM,rhotmp) 
                   OPij    = matmul(diagnormalg(:,:,i),rhomtx) 

                   rhomtx  = matmul(refs%MassM,diagnormalg(:,:,i))
                   OPtmp   = matmul(rhotmp,rhomtx) 

                   OPij    = (OPij + transpose(OPtmp))/2.0D0!*N2avg

                   ! symmetrize it 
                   gtmpij = normalg(:,i)*N2K1*rhot!/N2avg
                   !gtmpij = normalg(:,i)!*rhot*N2K1!/N2avg
                   call diagmatrix(rhotmp,gtmpij,pin%s%pNp) 
                   rhomtx  = matmul(refs%MassM,diagnormalg(:,:,j))
                   OPtrs   = matmul(rhotmp,rhomtx) 
                   
                   rhomtx  = matmul(refs%MassM,rhotmp) 
                   OPtmp   = matmul(diagnormalg(:,:,j),rhomtx) 
       
                   OPtrs   = (OPtrs + transpose(OPtmp))/2.0D0!*N2avg             
                
                   OPij    = (OPij + OPtrs)/2.0!*N2avg!*rhoavg
                   !OPij    = 0.0 
                   ! mass part & off-diagonal part
                   FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                      (j-1)*pin%s%pNp+1:j*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
                   FM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                      (j-1)*pin%s%pNp+1:j*pin%s%pNp) = 0.0D0 

                endif
             enddo; enddo
          else
             do i = 1,3 
                ! mass part
                FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (i-1)*pin%s%pNp+1:i*pin%s%pNp) = 0.0
                FM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                   (i-1)*pin%s%pNp+1:i*pin%s%pNp) = Mrho*unstrM%loc_detJ(k)

                ! \partial_i p
                fDerv(:,:,i) = matmul(reff%MassM,fDv2(:,:,i))
                fTDv(:,:,i)  = transpose(fDv2(:,:,i))
                fOPt(:,:,i)  = matmul(fTDv(:,:,i),reff%MassM) 
                ! make it symmetric
                fOPij        = (fDerv(:,:,i) + transpose(fOPt(:,:,i)))/2.0D0
                
                FT((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                       3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = fOPij*unstrM%loc_detJ(k)
                FM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                       3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = 0.0D0
                ! \div u
                FT(3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp,&
                  (i-1)*pin%s%pNp+1:i*pin%s%pNp) = transpose(fOPij)*unstrM%loc_detJ(k)
                FM(3*pin%f%pNp+1:3*pin%s%pNp+pin%f%pNp,&
                  (i-1)*pin%s%pNp+1:i*pin%s%pNp) = 0.0D0
             enddo
          endif
          ! p/kappa
          !FT(3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp,3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = &
          !           - unstrM%loc_detJ(k)*reff%MassM/(sum(lam)/real(pin%s%pNp,8))
          FT(3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp,3*pin%s%pNp+1:3*pin%s%pNp+pin%f%pNp) = &
                     - unstrM%loc_detJ(k)*flam


          ! JS 05/31/2020 deal with fluid boundary condition
          if (pin%JOB.ne.1.or.pin%JOB.ne.4) then  
             do i = 1,4 
                if (unstrM%ClNeigh(i,k).eq.-1) then
                   ! the surface is at the boundary !!
                   surfrho = sum(rhot(refs%Fmask(:,i)))/real(pin%f%Nfp,8)
                   surfgn  = 0.0D0
                   do j = 1,3
                      !surfgn = surfgn + gk1(refs%Fmask(:,i),j)**2
                      surfgn = surfgn-gk1(refs%Fmask(:,i),j)*unstrM%loc_n(j,i,k)
                   enddo
                   surfgn = surfgn**2 
                   !print*, surfrho,surfgn,i,k
                   !print*, surfrho,sgn*1.0D3,i,k
                   sgn = sum(dsqrt(surfgn))/real(pin%f%Nfp,8) 
                   surfp =  refs%MassF(:,:,i)/sgn/surfrho!*1.D3
                   !do j = 1,3
                   !   surfgn = surfgn-gk1(refs%Fmask(:,i),j)*unstrM%loc_n(j,i,k)
                   !enddo
                   !surfgn = 1.0D0/dsqrt(dabs(surfgn)*rhot(refs%Fmask(:,i)))
                   !call diagmatrix(difg,surfgn,pin%f%Nfp)
                   !surfp =  matmul(matmul(difg,refs%MassF(:,:,i)),difg)  
                   FT(3*pin%s%pNp+refs%Fmask(:,i),3*pin%s%pNp+refs%Fmask(:,i)) =&
                   FT(3*pin%s%pNp+refs%Fmask(:,i),3*pin%s%pNp+refs%Fmask(:,i)) +&
                           - surfp*unstrM%loc_sJac(i,k) 
                else ! todo check fluid-solid boundary
                   ! the surface is at the boundary !!
                   !surfrho = sum(models%coeff_loc((k-1)*pin%s%pNp+refs%Fmask(:,i),&
                   !               models%p_rho))/real(pin%f%Nfp,8)
                   surfrho = sum(rhot(refs%Fmask(:,i)))/real(pin%f%Nfp,8)
                   !surfrho = sum(rhot)/real(pin%f%pNp,8)
                   !print*, sgn-surfrho,sgn,k
                   surfgn  = 0.0D0
                   do j = 1,3
                      surfgn =surfgn + gk1(refs%Fmask(:,i),j)*unstrM%loc_n(j,i,k)
                   enddo 
                   sgn = sum(surfgn)/real(pin%f%Nfp,8) 
                   !print*, surfrho,surfgn,i,k
                   !print*, surfrho,sgn*1.0D3,i,k
                   ! check f-s boundary
                   cout = 0
                   do m = 1,4
                      if (m.ne.i) then
                         l = unstrM%lt2vid((k-1)*pin%s%pNp+refs%vord(m))
                         if (CGM%vnum(l).eq.6) then 
                            cout = cout + 1
                         endif
                      endif
                   enddo 
          
                   if (cout.lt.3) then  
                      do j = 1,3
                         !surfp =  -refs%MassF(:,:,i)*sgn*surfrho                
                         surfp =  refs%MassF(:,:,i)*sgn*surfrho                
                         FT((j-1)*pin%s%pNp+refs%Fmask(:,i),&
                            (j-1)*pin%s%pNp+refs%Fmask(:,i)) =&
                         FT((j-1)*pin%s%pNp+refs%Fmask(:,i),&
                            (j-1)*pin%s%pNp+refs%Fmask(:,i)) +&
                             surfp*unstrM%loc_sJac(i,k)*unstrM%loc_n(j,i,k)**2
                      enddo 
                   endif
                endif
             enddo

          endif ! end pin%JOB
          
          ! add the values into the matrices
          do m = 1,pin%f%pNp
             l = unstrM%lt2vid((k-1)*pin%s%pNp+m)
             if (unstrM%Cvpid(l).eq.unstrM%rank) then
                ! vtx id
                i = unstrM%loc_t2v(m,k)  
                call findorder(i,unstrM%new%vlist,vid)
                ! for Ap 
                j = sum(CGM%Ap%rnum(1:vid)) - CGM%Ap%rnum(vid) + 1
                les = CGM%Ap%rowdist(j+1)-CGM%Ap%rowdist(j)
                allocate(ntmp(les))
                ntmp = CGM%Ap%col(CGM%Ap%rowdist(j)+1:CGM%Ap%rowdist(j+1))
                do n = 1,pin%f%pNp
                   vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                   vid1 = CGM%pstt(vid0) + 1
                   call findorder(vid1,ntmp,colid)
                   CGM%Ap%val(CGM%Ap%rowdist(j)+colid) = &
                   CGM%Ap%val(CGM%Ap%rowdist(j)+colid) + &
                   FT(3*pin%s%pNp+m,3*pin%s%pNp+n)  
                enddo
                deallocate(ntmp)
                
                ! for ET 
                j = sum(CGM%ET%rnum(1:vid)) - CGM%ET%rnum(vid) + 1
                les = CGM%ET%rowdist(j+1)-CGM%ET%rowdist(j)
                allocate(ntmp(les))
                ntmp = CGM%ET%col(CGM%ET%rowdist(j)+1:CGM%ET%rowdist(j+1))
                do n = 1,pin%s%pNp
                   vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n)
                   do q = 1,3 
                      vid1 = CGM%vstt(vid0) + q + CGM%vnum(vid0) - 3
                      call findorder(vid1,ntmp,colid)
                      CGM%ET%val(CGM%ET%rowdist(j)+colid) = &
                      CGM%ET%val(CGM%ET%rowdist(j)+colid) + &
                      FT(3*pin%s%pNp+m,(q-1)*pin%s%pNp+n) 
                   enddo
                enddo
                 
                deallocate(ntmp)

                ! for Ad 
                do p = 1,3
                   ! fluid pts + CGM%Ad%rnum(vid) - 3 
                   j = sum(CGM%Ad%rnum(1:vid)) - 3 + p
                   les = CGM%Ad%rowdist(j+1)-CGM%Ad%rowdist(j)
                   allocate(ntmp(les))
                   ntmp = CGM%Ad%col(CGM%Ad%rowdist(j)+1:CGM%Ad%rowdist(j+1))
                   do n = 1,pin%s%pNp
                      vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                      do q = 1,3
                         vid1 = CGM%vstt(vid0) + q + CGM%vnum(vid0) - 3
                         call findorder(vid1,ntmp,colid)
                         CGM%Ad%val(CGM%Ad%rowdist(j)+colid) = &
                         CGM%Ad%val(CGM%Ad%rowdist(j)+colid) + &
                         FT((p-1)*pin%s%pNp+m,(q-1)*pin%s%pNp+n)  
                      enddo
                   enddo
                   deallocate(ntmp)
                   ! for B matrix
                   j = sum(CGM%B%rnum(1:vid)) - 3 + p
                   les = CGM%B%rowdist(j+1)-CGM%B%rowdist(j)
                   allocate(ntmp(les))
                   ntmp = CGM%B%col(CGM%B%rowdist(j)+1:CGM%B%rowdist(j+1))
                   do n = 1,pin%s%pNp
                      vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                      vid1 = CGM%vstt(vid0) + p + CGM%vnum(vid0) - 3 
                      call findorder(vid1,ntmp,colid)
                      CGM%B%val(CGM%B%rowdist(j)+colid) = &
                      CGM%B%val(CGM%B%rowdist(j)+colid) + &
                      FM((p-1)*pin%s%pNp+m,(p-1)*pin%s%pNp+n)  
                   enddo
                   deallocate(ntmp)
                   
                   ! for E matrix
                   j = sum(CGM%E%rnum(1:vid)) - 3 + p
                   les = CGM%E%rowdist(j+1)-CGM%E%rowdist(j)
                   allocate(ntmp(les))
                   ntmp = CGM%E%col(CGM%E%rowdist(j)+1:CGM%E%rowdist(j+1))
                   do n = 1,pin%f%pNp
                      vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                      vid1 = CGM%pstt(vid0) + 1
                      call findorder(vid1,ntmp,colid)
                      CGM%E%val(CGM%E%rowdist(j)+colid) = &
                      CGM%E%val(CGM%E%rowdist(j)+colid) + &
                      FT((p-1)*pin%s%pNp+m,3*pin%s%pNp+n)  
                   enddo
                   deallocate(ntmp)
 
                enddo
                

 
             endif
          enddo 
          
          ! fluid-soid
          if (.true.) then
             cout = 0
             do m = 1,4
                l = unstrM%lt2vid((k-1)*pin%s%pNp+refs%vord(m))
                if (CGM%vnum(l).eq.6) then 
                   cout = cout + 1
                endif
             enddo
 
             ! todo what if cout = 4 
             if (cout.eq.4) then
                print*,'warning: a tet',k, 'has probably two f-s',unstrM%rank
             endif

             if (cout.eq.3) then
                do m = 1,4
                   l = unstrM%lt2vid((k-1)*pin%s%pNp+refs%vord(m))
                   if (CGM%vnum(l).ne.6) then 
                      ! face id
                      fcid = m 
                   endif

                enddo
   
                ! compute the surface contribution term          
                if (pin%selfG) then
                   Fms = refs%Fmask(:,fcid)
                   rhof = sum(models%coeff_loc(Fms+(k-1)*pin%s%pNp,&
                           models%p_rho))/real(pin%s%Nfp)
                   SCM = 0.0D0
                   do i = 1,3; do j = 1,3
                      smassi  = unstrM%loc_sJac(fcid,k)*unstrM%loc_n(i,fcid,k)*&
                          matmul(refs%MassF(:,:,fcid),diaggk1(Fms,Fms,j))
                      smassit = unstrM%loc_sJac(fcid,k)*unstrM%loc_n(i,fcid,k)*&
                          matmul(diaggk1(Fms,Fms,j),refs%MassF(:,:,fcid))
                   
                      smassi = (smassi + transpose(smassit))/2.0D0
                      ! symmetrize it
                      smassj  = unstrM%loc_sJac(fcid,k)*unstrM%loc_n(j,fcid,k)*&
                          matmul(refs%MassF(:,:,fcid),diaggk1(Fms,Fms,i))
                      smassjt = unstrM%loc_sJac(fcid,k)*unstrM%loc_n(j,fcid,k)*&
                          matmul(diaggk1(Fms,Fms,i),refs%MassF(:,:,fcid))
                   
                      smassj = (smassj + transpose(smassjt))/2.0D0

                      smassij = (smassi + transpose(smassj))/2.0D0 

                      SCM((i-1)*pin%s%Nfp+1:i*pin%s%Nfp,&
                          (j-1)*pin%s%Nfp+1:j*pin%s%Nfp) = smassij*rhof

                   enddo; enddo
                   !print*, maxval(smassij)              
                   !print*, unstrM%loc_sJac(fcid,k),fcid,k

                   ! update Ad 
                   do m = 1,pin%s%Nfp
                      m0 = refs%Fmask(m,fcid) 
                      l   = unstrM%lt2vid((k-1)*pin%s%pNp+m0)
                      if (unstrM%Cvpid(l).eq.unstrM%rank) then
                         i = unstrM%loc_t2v(m0,k)  
                         call findorder(i,unstrM%new%vlist,vid)
                         do p = 1,3
                            ! fluid pts + CGM%Ad%rnum(vid) - 3 
                            j = sum(CGM%Ad%rnum(1:vid)) - CGM%Ad%rnum(vid) + p
                            les = CGM%Ad%rowdist(j+1)-CGM%Ad%rowdist(j)
                            allocate(ntmp(les))
                            ntmp = CGM%Ad%col(CGM%Ad%rowdist(j)+1:CGM%Ad%rowdist(j+1))
                            !nn = 0
                            do n = 1,pin%s%Nfp
                               !if (n.ne.fcid) then 
                                  !nn = nn + 1
                                  nn = refs%Fmask(n,fcid)
                                  vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+nn) 
                                  do q = 1,3
                                     vid1 = CGM%vstt(vid0) + q
                                     !print*,vid1
                                     call findorder(vid1,ntmp,colid)
                                     CGM%Ad%val(CGM%Ad%rowdist(j)+colid) = &
                                     CGM%Ad%val(CGM%Ad%rowdist(j)+colid) - &
                                     SCM((p-1)*pin%s%Nfp+m,(q-1)*pin%s%Nfp+n)  
                                  enddo
                               !endif
                            enddo
                            deallocate(ntmp)
                         enddo
                      endif
                   !endif
                   enddo
                endif
               
                ! for ET and E 
                do m = 1,pin%s%Nfp
                   m0 = refs%Fmask(m,fcid)
                   l   = unstrM%lt2vid((k-1)*pin%s%pNp+m0)
                   if (unstrM%Cvpid(l).eq.unstrM%rank) then
                      i = unstrM%loc_t2v(m0,k)  
                      call findorder(i,unstrM%new%vlist,vid)
                      ! for ET 
                      j = sum(CGM%ET%rnum(1:vid)) - CGM%ET%rnum(vid) + 1
                      les = CGM%ET%rowdist(j+1)-CGM%ET%rowdist(j)
                      allocate(ntmp(les))
                      ntmp = CGM%ET%col(CGM%ET%rowdist(j)+1:CGM%ET%rowdist(j+1))
                      !nn = 0
                      !do n = 1,4
                         !if (n.ne.fcid) then
                            !nn = nn + 1
                         do n = 1,pin%s%Nfp
                            nn = refs%Fmask(n,fcid)  
                            vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+nn)
                            !print*, norm2(unstrM%Cv_crs(:,vid0))
                            do q = 1,3 
                               vid1 = CGM%vstt(vid0) + q
                               call findorder(vid1,ntmp,colid)
                               CGM%ET%val(CGM%ET%rowdist(j)+colid) = &
                               CGM%ET%val(CGM%ET%rowdist(j)+colid) - &
                               unstrM%loc_n(q,fcid,k)*unstrM%loc_sJac(fcid,k)*&
                               refs%MassF(n,m,fcid)
                            enddo 
                         !endif
                      enddo
                      deallocate(ntmp) 
                      ! for E matrix
                      do p = 1,3
                         j = sum(CGM%E%rnum(1:vid)) - 6 + p
                         les = CGM%E%rowdist(j+1)-CGM%E%rowdist(j)
                         allocate(ntmp(les))
                         ntmp = CGM%E%col(CGM%E%rowdist(j)+1:CGM%E%rowdist(j+1))
                         !nn = 0
                         !do n = 1,4
                            !if (n.ne.fcid) then
                               !nn = nn + 1
                            do n = 1,pin%s%Nfp
                               nn = refs%Fmask(n,fcid)  
                            
                               vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+nn)
                               vid1 = CGM%pstt(vid0) + 1
                               call findorder(vid1,ntmp,colid)
                               CGM%E%val(CGM%E%rowdist(j)+colid) = &
                               CGM%E%val(CGM%E%rowdist(j)+colid) - &
                               unstrM%loc_n(p,fcid,k)*unstrM%loc_sJac(fcid,k)*&
                               refs%MassF(m,n,fcid)
                            enddo
                         !enddo
                         deallocate(ntmp) 
                      enddo
                   endif 
                enddo

             endif ! find f-s triangle
          endif ! if false

       endif ! mu > pin%TOL

    enddo 

    !print*, CGM%pstt

    !do i = 1,unstrM%nvtx
    !   if (CGM%Ad%rnum(i).eq.6) then
    !      k  = CGM%Ad%rstt(i) - CGM%Ad%sizdist(unstrM%rank+1) + 1
    !      nn = CGM%E%rowdist(k+1)-CGM%E%rowdist(k)
    !      print*, '#', nn, unstrM%v2vdist(i+1)-unstrM%v2vdist(i)
    !      do l = 1,nn
    !         j = CGM%E%rowdist(k) + l
    !         print*, CGM%E%sizdist(unstrM%rank+1)+i,CGM%E%col(j),CGM%E%val(j)
    !      enddo
    !   endif
    !enddo

    !do i = 1,CGM%E%NNZ
    !   print*,CGM%E%col(i),CGM%E%val(i)
    !enddo

  end subroutine CGFSE3D_ISO

  !--------------------------------------------------------------------------
  subroutine CGE3D_ISO()
    ! construct a local CG matrix for elastic isotropic case
    ! from a local DG matrix 
    integer, allocatable                                :: trow(:)
    real(kind=rkind), allocatable                       :: TM(:,:), MM(:,:)
    
    integer                                             :: i,j,k,l,m,n,p,q
    integer                                             :: dm,les,colid
    integer                                             :: vid,vid0,vid1
    integer, dimension(pin%s%pNp)                       :: rows1,cols1
    integer, allocatable                                :: ntmp(:)    

    real(kind=rkind), allocatable, dimension(:)         :: C33,C44
    
    real(kind=rkind), dimension(pin%s%pNp)              :: lam,mu,rhot
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: OPt,OP1,OP2
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: Derv,TDv,Dv2
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPij,OPs,OPtmp
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPtrs,Mrho,dgtmp

    real(kind=rkind)                                    :: lamavg,muavg,rhoavg

    ! needed for reference gravity
    real(kind=rkind), dimension(pin%s%pNp,3)            :: gk1,gk0,gradrho,normalg
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: rhomtx,rhotmp
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp,3)  :: diaggk1,diagnormalg
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPrhodi,OPrhodj,OPrhodij
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPrho,OPrhodit,OPrhodjt
    real(kind=rkind), dimension(pin%s%pNp,pin%s%pNp)    :: OPrho1,OPrho1t,OPrho2,OPrho2t

    ! for a more accurate derivatives
    real(kind=rkind), dimension(pin%s%pNp,3)            :: dism 
    real(kind=rkind), dimension(3,3)                    :: ndism,invdism,dgk1


    dm = 3; les = pin%s%pNp*dm; 
    allocate(trow(pin%s%pNp),TM(les,les),MM(les,les))

    allocate(C44(models%siz),C33(models%siz))
    C44 = models%coeff_loc(:,models%p_rho)*models%coeff_loc(:,models%p_vs)**2
    C33 = models%coeff_loc(:,models%p_rho)*models%coeff_loc(:,models%p_vp)**2 - 2*C44

    do k = 1,unstrM%ClNele
       ! build up global ids 
       do l = 1,pin%s%pNp 
          i = (k-1)*pin%s%pNp + l
          !j = (unstrM%Clelist(k)-1)*pin%s%pNp+l
          rows1(l) = i ! for coeff
          cols1(l) = i ! for global ids
       enddo
       ! for the coeffients
       lam   = C33(rows1) 
       mu    = C44(rows1) 
       rhot  = models%coeff_loc(rows1,models%p_rho)

       ! take the average 
       lamavg = sum(lam)/real(pin%s%pNp,8)
       muavg  = sum(mu)/real(pin%s%pNp,8)
       rhoavg = sum(rhot)/real(pin%s%pNp,8)

       
 
       ! compute the spatial derivatives  
       Derv(:,:,1) = unstrM%loc_invJ(1,1,k)*refs%Drst(:,:,1) +&
                     unstrM%loc_invJ(1,2,k)*refs%Drst(:,:,2) +&
                     unstrM%loc_invJ(1,3,k)*refs%Drst(:,:,3)
       
       Derv(:,:,2) = unstrM%loc_invJ(2,1,k)*refs%Drst(:,:,1) +&
                     unstrM%loc_invJ(2,2,k)*refs%Drst(:,:,2) +&
                     unstrM%loc_invJ(2,3,k)*refs%Drst(:,:,3)
       
       Derv(:,:,3) = unstrM%loc_invJ(3,1,k)*refs%Drst(:,:,1) +&
                     unstrM%loc_invJ(3,2,k)*refs%Drst(:,:,2) +&
                     unstrM%loc_invJ(3,3,k)*refs%Drst(:,:,3)


       ! prepare for reference gravity
       if (pin%selfG) then
          ! get g infomation
          gk1 = models%g0(:,:,k)/1.0D3
          !gk0 = gk1       
      
          ! construct for a accurate derivatives
          do i = 1,3
             dism(:,i) = unstrM%loc_nods(i,(k-1)*pin%s%pNp+1:k*pin%s%pNp) & 
                    -sum(unstrM%loc_nods(i,(k-1)*pin%s%pNp+1:k*pin%s%pNp))/real(pin%s%pNp,8)
             gk0(:,i) = gk1(:,i) - sum(gk1(:,i))/real(pin%s%pNp,8)
          enddo
          
          ndism = matmul(transpose(dism),dism) 
          call matinv(ndism,invdism,3)
          !print*, matmul(invdism,transpose(dism)),Derv(:,1,1)
          !print*, dism(:,1), Derv(:,1,1)

          !dgk1 = matmul(transpose(dism),gk1)      
          !dgk1 = matmul(invdism,dgk1)
          ! change JS 03/092018
          dgk1 = matmul(transpose(dism),gk0)      
          dgk1 = matmul(invdism,dgk1)
 
          diaggk1 = 0.0D0 
          do i = 1,3
             call diagmatrix(diaggk1(:,:,i),gk1(:,i),pin%s%pNp)
          enddo
       endif


       ! build the mass matrix part
       !Mrho  = (refs%MassM*rhoavg + rhoavg*refs%MassM)/2.0
       call diagmatrix(Mrho,rhot,pin%s%pNp)
       Mrho = matmul(Mrho,refs%MassM)+matmul(refs%MassM,Mrho)
       Mrho = Mrho/2.0D0

       OPs  = 0.0D0
       ! build up derivative matricies
       Dv2  = Derv
       !do i = 1,dm
       !   ! lambda*D_i
       !   OP1(:,:,i) = Derv(:,:,i)*lamavg
       !   ! mu*D_i
       !   OP2(:,:,i) = Derv(:,:,i)*muavg
       !   ! D_i'*mu*D_i
       !   Derv(:,:,i) = matmul(refs%MassM,Dv2(:,:,i))

       !   TDv(:,:,i)  = transpose(Derv(:,:,i))
       !   OPt(:,:,i)  = matmul(TDv(:,:,i),OP2(:,:,i))
       !   OPs         = OPs + OPt(:,:,i)
       !enddo
       ! build up derivative matricies
       do i = 1,dm
          ! lambda*D_i
          !call realmcolupdate(Derv(:,:,i),lam,OP1(:,:,i),pin%s%pNp)
          call diagmatrix(dgtmp,lam,pin%s%pNp)
          dgtmp = matmul(refs%MassM,dgtmp)
          dgtmp = (dgtmp+transpose(dgtmp))/2.0D0
          OP1(:,:,i) = matmul(dgtmp,Derv(:,:,i))
          ! mu*D_i
          !call realmcolupdate(Derv(:,:,i),mu,OP2(:,:,i),pin%s%pNp) 
          call diagmatrix(dgtmp,mu,pin%s%pNp)
          dgtmp = matmul(refs%MassM,dgtmp)
          dgtmp = (dgtmp+transpose(dgtmp))/2.0D0
          OP2(:,:,i) = matmul(dgtmp,Derv(:,:,i))
 
          ! D_i'*mu*D_i
          !Derv(:,:,i) = matmul(refs%MassM,Dv2(:,:,i))
          !Dv2(:,:,i)  = Derv(:,:,i)
          TDv(:,:,i)  = transpose(Derv(:,:,i))
          OPt(:,:,i)  = matmul(TDv(:,:,i),OP2(:,:,i))
          OPs         = OPs + OPt(:,:,i)
       enddo

       do i = 1,dm; do j = 1,dm
          if (pin%selfG) then
             ! prepare for g related terms  
             ! for OPrho
             OPrho = 0.0D0
             !call diagmatrix(rhotmp,gk1(:,j),pNp)
             rhomtx   = matmul(refs%MassM,diaggk1(:,:,j))
             OPrho1   = matmul(transpose(Dv2(:,:,i)),rhomtx)
             rhomtx   = matmul(refs%MassM,Dv2(:,:,i))
             OPrho1t  = matmul(diaggk1(:,:,j),rhomtx)

             OPrho1   = (OPrho1 + transpose(OPrho1t))/2.0D0

             ! symmetrize it 
             !call diagmatrix(rhotmp,rhot*gk1(:,i),pNp)
             rhomtx   = matmul(refs%MassM,Dv2(:,:,j))
             OPrho2   = matmul(diaggk1(:,:,i),rhomtx)
             rhomtx   = matmul(refs%MassM,diaggk1(:,:,i))
             OPrho2t  = matmul(transpose(Dv2(:,:,j)),rhomtx) 
             OPrho2   = (OPrho2 + transpose(OPrho2t))/2.0D0

             OPrho    = (OPrho1 + OPrho2)/2.0D0 
         
             rhomtx   = refs%MassM*(dgk1(i,j)+dgk1(j,i))/2.0D0
             OPrho    = OPrho - rhomtx 

             !if (mod(k,1000)==1.and. i==3.and. j==1) print*,'2',maxval(rhomtx)

             rhotmp  = matmul(diaggk1(:,:,i),refs%MassM)
             OPrho1  = matmul(transpose(Dv2(:,:,j)),rhotmp) 
             rhotmp  = matmul(diaggk1(:,:,i),Dv2(:,:,j))
             OPrho1t = matmul(refs%MassM,rhotmp)
             OPrho1  = (OPrho1 + transpose(OPrho1t))/2.0D0 

             ! symmetrize it              
             rhotmp  = matmul(diaggk1(:,:,j),Dv2(:,:,i))
             OPrho2  = matmul(refs%MassM,rhotmp)
             rhotmp  = matmul(diaggk1(:,:,j),refs%MassM)
             OPrho2t = matmul(transpose(Dv2(:,:,i)),rhotmp) 
             OPrho2  = (OPrho2 + transpose(OPrho2t))/2.0D0
             
             OPrho   = OPrho - (OPrho1+OPrho2)/2.0D0

             !OPrho   = OPrho*rhoavg
             call diagmatrix(dgtmp,rhot,pin%s%pNp)
             dgtmp   = matmul(OPrho,dgtmp)+matmul(dgtmp,OPrho)
             OPrho   = dgtmp/2.0D0
          else 
             OPrho   = 0.0D0
          endif 




          if (i == j) then
             ! diagonal component
             OPtmp = matmul(TDv(:,:,i),OP1(:,:,i))
             OPij  = OPs + OPt(:,:,i) + OPtmp                             
             
             ! make it exactly symetric
             OPij = (OPij + transpose(OPij))/2.0D0 

             !TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
             !   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
             TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                (j-1)*pin%s%pNp+1:j*pin%s%pNp) =  unstrM%loc_detJ(k)*(OPij+OPrho)
             MM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                (j-1)*pin%s%pNp+1:j*pin%s%pNp) = Mrho*unstrM%loc_detJ(k)
          else
             ! off diagonal component
             OPtmp = matmul(TDv(:,:,i),OP1(:,:,j)) + &
                     matmul(TDv(:,:,j),OP2(:,:,i))
             OPtrs = matmul(TDv(:,:,j),OP1(:,:,i)) + &
                     matmul(TDv(:,:,i),OP2(:,:,j))
             OPij  = (transpose(OPtrs) + OPtmp)/2.0D0
            
             !if (k==1) print*, OPij
 
             !TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
             !   (j-1)*pin%s%pNp+1:j*pin%s%pNp) = OPij*unstrM%loc_detJ(k)
             TM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                (j-1)*pin%s%pNp+1:j*pin%s%pNp) = unstrM%loc_detJ(k)*(OPij+OPrho)

             !if (k==3) print*, unstrM%loc_detJ(k)*(OPij+OPrho)

             MM((i-1)*pin%s%pNp+1:i*pin%s%pNp,&
                (j-1)*pin%s%pNp+1:j*pin%s%pNp) = 0.0D0
          endif
       enddo; enddo


       ! finds locations
       do m = 1,pin%s%pNp
          ! get locat information
          l   = unstrM%lt2vid((k-1)*pin%s%pNp+m)
          if (unstrM%Cvpid(l).eq.unstrM%rank) then
             ! vtx id
             i = unstrM%loc_t2v(m,k)  
             call findorder(i,unstrM%new%vlist,vid)
             do p = 1,3
                ! for A
                j = (vid-1)*3+p
                les = CGM%A%rowdist(j+1)-CGM%A%rowdist(j)
                allocate(ntmp(les))
                ntmp = CGM%A%col(CGM%A%rowdist(j)+1:CGM%A%rowdist(j+1))
                do n = 1,pin%s%pNp
                   vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                   do q = 1,3
                      vid1 = CGM%vstt(vid0)+q
                      call findorder(vid1,ntmp,colid)
                      CGM%A%val(CGM%A%rowdist(j)+colid) = &
                      CGM%A%val(CGM%A%rowdist(j)+colid) + &
                      TM((p-1)*pin%s%pNp+m,(q-1)*pin%s%pNp+n)  
                   enddo
                enddo
                deallocate(ntmp)
                ! for B
                j = (vid-1)*3+p
                les = CGM%B%rowdist(j+1)-CGM%B%rowdist(j)
                allocate(ntmp(les))
                ntmp = CGM%B%col(CGM%B%rowdist(j)+1:CGM%B%rowdist(j+1))
                do n = 1,pin%s%pNp
                   vid0 = unstrM%lt2vid((k-1)*pin%s%pNp+n) 
                   vid1 = CGM%vstt(vid0)+p
                   call findorder(vid1,ntmp,colid)  
                   CGM%B%val(CGM%B%rowdist(j)+colid) = &
                   CGM%B%val(CGM%B%rowdist(j)+colid) + &
                   MM((p-1)*pin%s%pNp+m,(p-1)*pin%s%pNp+n)  
                enddo
                deallocate(ntmp)
             enddo
          endif
       enddo

    enddo
  end subroutine CGE3D_ISO

  !--------------------------------------------------------------------
  subroutine matrixstruct()
    integer                                 :: i,j,k,l,m,n,ierr
    integer                                 :: lsiz0,lsiz1,location
    integer, allocatable, dimension(:)      :: vtmp0,vtmp1,cnt0,cnt1
    integer, allocatable, dimension(:)      :: itmp0,itmp1,dtmp0,dtmp1
    integer, allocatable, dimension(:)      :: col0,idtmp
  

    ! A information 
    allocate(CGM%A%rnum(unstrM%new%nvtx)); 
    CGM%A%rnum = 3

    CGM%A%siz = sum(CGM%A%rnum)

    allocate(vtmp0(unstrM%nproc)); vtmp0 = 0
    allocate(vtmp1(unstrM%nproc)); vtmp1 = 0

    vtmp0(unstrM%rank+1) = CGM%A%siz; 

    call mpi_allreduce(vtmp0,vtmp1,unstrM%nproc,&
                     mpi_integer,mpi_sum,unstrM%comm,ierr)

    allocate(CGM%A%sizdist(unstrM%nproc+1)); CGM%A%sizdist = 0
    do i = 1,unstrM%nproc
       CGM%A%sizdist(i+1) = CGM%A%sizdist(i) + vtmp1(i)
    enddo
    CGM%A%Gsiz = CGM%A%sizdist(unstrM%nproc+1)
    !print*, 'local A size',CGM%A%siz,'at rank',int(unstrM%rank) 
    !print*, 'total A size',CGM%A%Gsiz          

    allocate(CGM%A%rstt(unstrM%new%nvtx))
    CGM%A%rstt(1) = CGM%A%sizdist(unstrM%rank+1)
    do i = 2,unstrM%new%nvtx
       CGM%A%rstt(i) = CGM%A%rstt(i-1) + CGM%A%rnum(i-1)
    enddo

    ! set up conmmunication
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
    enddo 
    
    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,unstrM%comm,ierr)
    
    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + vtmp0(i)
       cnt1(i) = cnt1(i-1) + vtmp1(i)
    enddo

    allocate(itmp0(unstrM%cnvtx))
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       itmp0(k) = unstrM%Cvlist(i) 
    enddo
    lsiz1 = sum(vtmp1); allocate(itmp1(lsiz1)) 
 
    call MPI_ALLTOALLV(itmp0,vtmp0,cnt0,mpi_integer,&
                       itmp1,vtmp1,cnt1,mpi_integer,unstrM%comm,ierr)

    allocate(dtmp1(lsiz1))
    do i = 1,lsiz1
       j = itmp1(i)
       call findorder(j,unstrM%new%vlist,location)
       ! offer start locations
       dtmp1(i) = CGM%A%rstt(location) 
    enddo
    allocate(dtmp0(unstrM%cnvtx)) 
    
    call MPI_ALLTOALLV(dtmp1,vtmp1,cnt1,mpi_integer,&
                       dtmp0,vtmp0,cnt0,mpi_integer,unstrM%comm,ierr)
    allocate(CGM%vstt(unstrM%cnvtx)) 
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       CGM%vstt(i) = dtmp0(k)
    enddo
 
    allocate(CGM%vnum(unstrM%cnvtx)); CGM%vnum = 3
  
    ! later
    allocate(CGM%A%rowdist(CGM%A%siz+1))
    CGM%A%rowdist(1) = 0; k = 1
    do i = 1,unstrM%new%nvtx
       do j = 1,CGM%A%rnum(i)
          k = k + 1
          l = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
          CGM%A%rowdist(k) = CGM%A%rowdist(k-1) + l*3
       enddo
    enddo
    CGM%A%NNZ = CGM%A%rowdist(CGM%A%siz+1)
    !print*,CGM%A%NNZ,unstrM%rank

    allocate(CGM%A%col(CGM%A%NNZ),CGM%A%val(CGM%A%NNZ))
    CGM%A%col = -1; CGM%A%val = 0.0D0
    k = 1
    do i = 1,unstrM%new%nvtx
       do j = 1,CGM%A%rnum(i)
          k = k + 1
          lsiz0 = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
          allocate(col0(lsiz0*3))
          do l = 1,lsiz0
             m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
             call findorder(m,unstrM%Cvlist,location)
             do n = 1,CGM%vnum(location)
                col0((l-1)*3+n) = CGM%vstt(location) + n
             enddo
          enddo
          call simplessort(col0,idtmp)
          CGM%A%col(CGM%A%rowdist(k-1)+1:CGM%A%rowdist(k))=col0
          deallocate(col0,idtmp)
       enddo
    enddo
    !print*, minval(CGM%A%col),maxval(CGM%A%col),unstrM%rank

    ! B information
    allocate(CGM%B%rnum(unstrM%new%nvtx))
    CGM%B%rnum = 3
    CGM%B%siz  = CGM%A%siz 
    CGM%B%Gsiz = CGM%A%Gsiz 
    allocate(CGM%B%sizdist(unstrM%nproc+1))
    CGM%B%sizdist = CGM%A%sizdist
    allocate(CGM%B%rstt(unstrM%new%nvtx))
    CGM%B%rstt = CGM%A%rstt

    ! later
    allocate(CGM%B%rowdist(CGM%B%siz+1))
    CGM%B%rowdist(1) = 0; k = 1
    do i = 1,unstrM%new%nvtx
       do j = 1,CGM%B%rnum(i)
          k = k + 1
          l = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
          CGM%B%rowdist(k) = CGM%B%rowdist(k-1) + l
       enddo
    enddo
    CGM%B%NNZ = CGM%B%rowdist(CGM%A%siz+1)

    allocate(CGM%B%col(CGM%B%NNZ),CGM%B%val(CGM%B%NNZ))
    CGM%B%col = -1; CGM%B%val = 0.0D0
    k = 1
    do i = 1,unstrM%new%nvtx
       do j = 1,CGM%B%rnum(i)
          k = k + 1
          lsiz0 = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
          allocate(col0(lsiz0))
          do l = 1,lsiz0
             m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
             call findorder(m,unstrM%Cvlist,location)
             !do n = 1,CGM%vnum(location)
             !CGM%B%col(CGM%B%rowdist(k-1)+l) =&
             !CGM%vstt(location) + j
             !enddo
             col0(l) = CGM%vstt(location) + j
          enddo
          call simplessort(col0,idtmp)
          CGM%B%col(CGM%B%rowdist(k-1)+1:CGM%B%rowdist(k))=col0
          deallocate(col0,idtmp)
       enddo
    enddo
    !print*, minval(CGM%B%col)

    call mpi_allreduce(CGM%A%NNZ,CGM%A%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)

    print*, 'local A NNZ', CGM%A%NNZ, 'at rank', int(unstrM%rank,2)
    call mpi_barrier(unstrM%comm,ierr)

    call mpi_allreduce(CGM%B%NNZ,CGM%B%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
    
    print*, 'local B NNZ', CGM%B%NNZ, 'at rank', int(unstrM%rank,2)
    call mpi_barrier(unstrM%comm,ierr)

    if (unstrM%rank.eq.unstrM%nproc-1) then
       print*, 'total A NNZ', CGM%A%GNNZ
       print*, 'total B NNZ', CGM%B%GNNZ
    endif

  end subroutine matrixstruct

  !-----------------------------------------------------------------
  subroutine matrixstruct_general()
    integer                                 :: i,j,k,l,kk,ll,m,n,mm,nn,ierr
    integer                                 :: lsiz0,lsiz1,location
    integer, allocatable, dimension(:)      :: vtmp0,vtmp1,cnt0,cnt1
    integer, allocatable, dimension(:)      :: itmp0,itmp1,dtmp0,dtmp1,dtmp2
    integer, allocatable, dimension(:)      :: col0,col1,col2,idtmp

    ! matrices information 
    ! find rnum and rstt 
    allocate(CGM%Ad%rnum(unstrM%new%nvtx)); CGM%Ad%rnum = 3
    do i = 1,unstrM%new%nvtx
       if (unstrM%new%vstat(i).eq.2.or.&
           unstrM%new%vstat(i).eq.5) then
          CGM%Ad%rnum(i) = 6
       endif
    enddo

    CGM%Ad%siz = sum(CGM%Ad%rnum)

    allocate(vtmp0(unstrM%nproc)); vtmp0 = 0
    allocate(vtmp1(unstrM%nproc)); vtmp1 = 0

    vtmp0(unstrM%rank+1) = CGM%Ad%siz; 

    call mpi_allreduce(vtmp0,vtmp1,unstrM%nproc,&
                     mpi_integer,mpi_sum,unstrM%comm,ierr)

    allocate(CGM%Ad%sizdist(unstrM%nproc+1)); CGM%Ad%sizdist = 0
    do i = 1,unstrM%nproc
       CGM%Ad%sizdist(i+1) = CGM%Ad%sizdist(i) + vtmp1(i)
    enddo
    CGM%Ad%Gsiz = CGM%Ad%sizdist(unstrM%nproc+1)
    call mpi_barrier(unstrM%comm,ierr)
    print*, 'local Ad size',CGM%Ad%siz,'at rank',int(unstrM%rank,2) 
    if (unstrM%rank.eq.unstrM%nproc-1) print*, 'total Ad size',CGM%Ad%Gsiz         
    !if (unstrM%rank.eq.0) print*, CGM%Ad%sizdist

    allocate(CGM%Ad%rstt(unstrM%new%nvtx))
    CGM%Ad%rstt(1) = CGM%Ad%sizdist(unstrM%rank+1)
    do i = 2,unstrM%new%nvtx
       CGM%Ad%rstt(i) = CGM%Ad%rstt(i-1) + CGM%Ad%rnum(i-1)
    enddo

    ! set up conmmunication
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
    enddo 
    
    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,unstrM%comm,ierr)
    
    allocate(cnt0(0:unstrM%nproc-1),cnt1(0:unstrM%nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,unstrM%nproc-1
       cnt0(i) = cnt0(i-1) + vtmp0(i)
       cnt1(i) = cnt1(i-1) + vtmp1(i)
    enddo

    allocate(itmp0(unstrM%cnvtx))
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       itmp0(k) = unstrM%Cvlist(i) 
    enddo
    lsiz1 = sum(vtmp1); allocate(itmp1(lsiz1)) 

    call MPI_ALLTOALLV(itmp0,vtmp0,cnt0,mpi_integer,&
                       itmp1,vtmp1,cnt1,mpi_integer,unstrM%comm,ierr)

    allocate(dtmp1(lsiz1),dtmp2(lsiz1))
    do i = 1,lsiz1
       j = itmp1(i)
       call findorder(j,unstrM%new%vlist,location)
       ! offer start locations
       dtmp1(i) = CGM%Ad%rstt(location) 
       dtmp2(i) = CGM%Ad%rnum(location)
    enddo
    allocate(dtmp0(unstrM%cnvtx)) 
    
    call MPI_ALLTOALLV(dtmp1,vtmp1,cnt1,mpi_integer,&
                       dtmp0,vtmp0,cnt0,mpi_integer,unstrM%comm,ierr)
    allocate(CGM%vstt(unstrM%cnvtx)) 
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       CGM%vstt(i) = dtmp0(k)
    enddo
    deallocate(dtmp0); allocate(dtmp0(unstrM%cnvtx))
    call MPI_ALLTOALLV(dtmp2,vtmp1,cnt1,mpi_integer,&
                       dtmp0,vtmp0,cnt0,mpi_integer,unstrM%comm,ierr)
   
    allocate(CGM%vnum(unstrM%cnvtx)) 
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       CGM%vnum(i) = dtmp0(k)
    enddo
    deallocate(dtmp0)
    !print*, maxval(CGM%vnum), minval(CGM%vnum)

    ! B information
    allocate(CGM%B%rnum(unstrM%new%nvtx))
    CGM%B%rnum = CGM%Ad%rnum
    allocate(CGM%B%sizdist(unstrM%nproc+1))
    CGM%B%sizdist = CGM%Ad%sizdist
    allocate(CGM%B%rstt(unstrM%new%nvtx))
    CGM%B%rstt = CGM%Ad%rstt

    CGM%B%siz  = CGM%Ad%siz; CGM%B%Gsiz = CGM%Ad%Gsiz
    call mpi_barrier(unstrM%comm,ierr)
    print*, 'local B size',CGM%B%siz,'at rank',int(unstrM%rank,2) 
    if (unstrM%rank.eq.unstrM%nproc-1) print*, 'total B size',CGM%B%Gsiz          

    !**************************************************************
    ! Ap information
    allocate(CGM%Ap%rnum(unstrM%new%nvtx)); CGM%Ap%rnum = 0
    do i = 1,unstrM%new%nvtx
       !if (unstrM%new%vstat(i).ge.1) then
       if (unstrM%new%vstat(i).ne.0.and.&
           unstrM%new%vstat(i).ne.3) then
          CGM%Ap%rnum(i) = 1
       endif
    enddo

    CGM%Ap%siz = sum(CGM%Ap%rnum)
    vtmp0 = 0; vtmp1 = 0
    vtmp0(unstrM%rank+1) = CGM%Ap%siz; 
    call mpi_allreduce(vtmp0,vtmp1,unstrM%nproc,&
                     mpi_integer,mpi_sum,unstrM%comm,ierr)

    allocate(CGM%Ap%sizdist(unstrM%nproc+1)); CGM%Ap%sizdist = 0
    do i = 1,unstrM%nproc
       CGM%Ap%sizdist(i+1) = CGM%Ap%sizdist(i) + vtmp1(i)
    enddo
    CGM%Ap%Gsiz = CGM%Ap%sizdist(unstrM%nproc+1)
    call mpi_barrier(unstrM%comm,ierr)
    print*, 'local Ap size',CGM%Ap%siz,'at rank',int(unstrM%rank,2) 
    if (unstrM%rank.eq.unstrM%nproc-1) print*, 'total Ap size',CGM%Ap%Gsiz         

    allocate(CGM%Ap%rstt(unstrM%new%nvtx)); CGM%Ap%rstt = -1
    !CGM%Ap%rstt(1) = CGM%Ap%sizdist(unstrM%rank+1)
    j = 0
    do i = 1,unstrM%new%nvtx
       if (CGM%Ap%rnum(i) .eq. 1) then
          CGM%Ap%rstt(i) = j + CGM%Ap%sizdist(unstrM%rank+1)
          j = j + 1
       endif
    enddo
    !print*, CGM%Ap%rstt 
    ! set up conmmunication
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
    enddo 
    
    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,unstrM%comm,ierr)

    do i = 1,lsiz1
       j = itmp1(i)
       call findorder(j,unstrM%new%vlist,location)
       ! offer start locations
       dtmp1(i) = CGM%Ap%rstt(location) 
       dtmp2(i) = CGM%Ap%rnum(location)
    enddo
    allocate(dtmp0(unstrM%cnvtx))
    call MPI_ALLTOALLV(dtmp1,vtmp1,cnt1,mpi_integer,&
                       dtmp0,vtmp0,cnt0,mpi_integer,unstrM%comm,ierr)
   
    allocate(CGM%pstt(unstrM%cnvtx)); CGM%pstt = -1 
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       CGM%pstt(i) = dtmp0(k)
    enddo
    deallocate(dtmp0); allocate(dtmp0(unstrM%cnvtx))
    !print*, CGM%pstt

    call MPI_ALLTOALLV(dtmp2,vtmp1,cnt1,mpi_integer,&
                       dtmp0,vtmp0,cnt0,mpi_integer,unstrM%comm,ierr)
   
    allocate(CGM%pnum(unstrM%cnvtx)); CGM%pnum = 0 
    vtmp0 = 0
    do i = 1,unstrM%cnvtx
       j = unstrM%Cvpid(i) + 1
       vtmp0(j) = vtmp0(j) + 1
       k = cnt0(j-1) + vtmp0(j) 
       CGM%pnum(i) = dtmp0(k)
    enddo
    deallocate(dtmp0)
    !print*, CGM%pnum

    ! E and ET
    allocate(CGM%ET%rnum(unstrM%new%nvtx)); CGM%ET%rnum = 0
    CGM%ET%rnum = CGM%Ap%rnum
    CGM%ET%siz  = CGM%Ap%siz; CGM%ET%Gsiz = CGM%Ap%Gsiz
    allocate(CGM%ET%sizdist(unstrM%nproc+1))
    CGM%ET%sizdist = CGM%Ap%sizdist
 
    allocate(CGM%E%rnum(unstrM%new%nvtx)); CGM%E%rnum = 0
    CGM%E%rnum = CGM%Ad%rnum
    CGM%E%siz  = CGM%Ad%siz; CGM%E%Gsiz = CGM%Ad%Gsiz
    allocate(CGM%E%sizdist(unstrM%nproc+1))
    CGM%E%sizdist = CGM%Ad%sizdist


    ! Todo check num of nonzeros
    allocate(CGM%Ad%rowdist(CGM%Ad%siz+1)); CGM%Ad%rowdist(1) = 0
    allocate( CGM%B%rowdist(CGM%Ad%siz+1));  CGM%B%rowdist(1) = 0
    allocate(CGM%Ap%rowdist(CGM%Ap%siz+1)); CGM%Ap%rowdist(1) = 0
    allocate(CGM%ET%rowdist(CGM%ET%siz+1)); CGM%ET%rowdist(1) = 0
    allocate( CGM%E%rowdist( CGM%E%siz+1));  CGM%E%rowdist(1) = 0

    k = 1; kk = 1
    do i = 1,unstrM%new%nvtx
       lsiz0 = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
       if (CGM%Ap%rnum(i).eq.1) then
          if (CGM%Ad%rnum(i).eq.3) then
             ! pure fluid
             ll = 0
             do j = 1,lsiz0 
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.1) ll = ll + 1
             enddo
             if (ll.ne.lsiz0) then
                print*, 'Error: pure fluid', ll, lsiz0
                stop
             endif
             CGM%Ap%rowdist(kk+1) = CGM%Ap%rowdist(kk) + ll
             CGM%ET%rowdist(kk+1) = CGM%ET%rowdist(kk) + ll*3
             kk = kk + 1

             do j = 1,3
                CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + lsiz0*3
                CGM%B%rowdist(k+1)  = CGM%B%rowdist(k)  + lsiz0
                CGM%E%rowdist(k+1)  = CGM%E%rowdist(k)  + ll
                k = k + 1
             enddo 
          else
             ! fluid-solid pts
             ll = 0; nn = 0
             do j = 1,lsiz0 
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.0) then
                   ll = ll + 1
                endif
                if (CGM%vnum(location).eq.6) nn = nn + 1
             enddo

             do j = 1,3 ! solid part
                CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + (nn+ll)*3
                CGM%B%rowdist(k+1)  =  CGM%B%rowdist(k) + (nn+ll)
                CGM%E%rowdist(k+1)  =  CGM%E%rowdist(k) + nn
                k = k + 1
             enddo
             
             CGM%Ap%rowdist(kk+1) = CGM%Ap%rowdist(kk) + (lsiz0-ll)
             CGM%ET%rowdist(kk+1) = CGM%ET%rowdist(kk) + (lsiz0-ll+nn)*3
             kk = kk + 1

             do j = 4,6 ! fluid part
                CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + (lsiz0-ll)*3
                CGM%B%rowdist(k+1)  =  CGM%B%rowdist(k) + (lsiz0-ll)
                CGM%E%rowdist(k+1)  =  CGM%E%rowdist(k) + (lsiz0-ll) 
                k = k + 1
             enddo

          endif
       else
          if (CGM%Ad%rnum(i).eq.3) then
             ! pure solid
             ll = 0
             do j = 1,lsiz0 
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.1) ll = ll + 1
             enddo

             do j = 1,3
                CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + lsiz0*3
                CGM%B%rowdist(k+1)  = CGM%B%rowdist(k)  + lsiz0
                CGM%E%rowdist(k+1)  = CGM%E%rowdist(k)  + 0
                k = k + 1
             enddo 

          else
             ! fluid-solid pts
             print*, 'Error: pure solid pts',i,int(unstrM%rank,2)
          endif
       endif
    enddo 
    if (kk.ne.CGM%Ap%siz+1) then
       print*, 'Error: Ap size', kk, CGM%Ap%siz+1
    endif

    CGM%Ad%NNZ = CGM%Ad%rowdist(CGM%Ad%siz+1)
    CGM%B%NNZ  = CGM%B%rowdist(CGM%B%siz+1)
    CGM%E%NNZ  = CGM%E%rowdist(CGM%E%siz+1)
    CGM%ET%NNZ = CGM%ET%rowdist(CGM%ET%siz+1)
    CGM%Ap%NNZ = CGM%Ap%rowdist(CGM%Ap%siz+1)

    call mpi_barrier(unstrM%comm,ierr)
    print*,'local Ad NNZ',int(CGM%Ad%NNZ,4),'at rank',int(unstrM%rank,2)
    
    call mpi_allreduce(CGM%Ad%NNZ,CGM%Ad%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
   
    call mpi_allreduce(CGM%B%NNZ,CGM%B%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
    
    call mpi_allreduce(CGM%E%NNZ,CGM%E%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
   
    call mpi_allreduce(CGM%ET%NNZ,CGM%ET%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
   
    call mpi_allreduce(CGM%Ap%NNZ,CGM%Ap%GNNZ,1,mpi_integer8,&
           mpi_sum,unstrM%comm,ierr)
   

    if (unstrM%rank.eq.unstrM%nproc-1) then
       print*, 'total NNZ: Ad, B, E, ET, Ap'
       print*, CGM%Ad%GNNZ,CGM%B%GNNZ,CGM%E%GNNZ,CGM%ET%GNNZ,CGM%Ap%GNNZ
    endif
 
    
    ! for the col information
    allocate(CGM%Ad%col(CGM%Ad%NNZ),CGM%Ad%val(CGM%Ad%NNZ))
    CGM%Ad%col = 0; CGM%Ad%val = 0.0D0
    allocate(CGM%Ap%col(CGM%Ap%NNZ),CGM%Ap%val(CGM%Ap%NNZ))
    CGM%Ap%col = 0; CGM%Ap%val = 0.0D0
    allocate(CGM%B%col(CGM%B%NNZ),CGM%B%val(CGM%B%NNZ))
    CGM%B%col = 0; CGM%B%val = 0.0D0
    allocate(CGM%E%col(CGM%E%NNZ),CGM%E%val(CGM%E%NNZ))
    CGM%E%col = 0; CGM%E%val = 0.0D0
    !call mpi_barrier(unstrM%comm,ierr)
    !print*,'local E NNZ ',int(CGM%E%NNZ,4), 'at rank',int(unstrM%rank,2)
    allocate(CGM%ET%col(CGM%ET%NNZ),CGM%ET%val(CGM%ET%NNZ))
    CGM%ET%col = 0; CGM%ET%val = 0.0D0
    !call mpi_barrier(unstrM%comm,ierr)
    !print*,'local ET NNZ',int(CGM%ET%NNZ,4),'at rank',int(unstrM%rank,2)

    call mpi_barrier(unstrM%comm,ierr)
    print*,'local E NNZ ',int(CGM%E%NNZ,4), 'at rank',int(unstrM%rank,2)

    k = 1; kk = 1
    do i = 1,unstrM%new%nvtx
       lsiz0 = unstrM%new%v2vdist(i+1)-unstrM%new%v2vdist(i)
       if (CGM%Ap%rnum(i).eq.1) then
          if (CGM%Ad%rnum(i).eq.3) then
             ! pure fluid
             !CGM%Ap%rowdist(kk+1) = CGM%Ap%rowdist(kk) + lsiz0
             !CGM%ET%rowdist(kk+1) = CGM%ET%rowdist(kk) + lsiz0*3
             allocate(col0(lsiz0*3),col1(lsiz0),col2(lsiz0))
             do l = 1,lsiz0
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
                call findorder(m,unstrM%Cvlist,location)
                col2(l) = CGM%pstt(location) + 1
                if (CGM%pstt(location).lt.0) then
                   print*, 'Error'
                   stop
                endif  
                do n = 1,3
                   col0((l-1)*3+n) = CGM%vstt(location) + n &
                                   + CGM%vnum(location) - 3
                enddo
             enddo
             call simplessort(col2,idtmp); deallocate(idtmp) 
             CGM%Ap%col(CGM%Ap%rowdist(kk)+1:CGM%Ap%rowdist(kk+1)) = col2
             call simplessort(col0,idtmp); deallocate(idtmp)
             CGM%ET%col(CGM%ET%rowdist(kk)+1:CGM%ET%rowdist(kk+1)) = col0
             kk = kk + 1

             do j = 1,3
                !CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + lsiz0*3
                !CGM%B%rowdist(k+1)  = CGM%B%rowdist(k)  + lsiz0
                !CGM%E%rowdist(k+1)  = CGM%E%rowdist(k)  + lsiz0
                do l = 1,lsiz0
                   m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
                   call findorder(m,unstrM%Cvlist,location)
                   col1(l) = CGM%vstt(location) + j + CGM%vnum(location) - 3
                enddo
                call simplessort(col1,idtmp); deallocate(idtmp)
                CGM%Ad%col(CGM%Ad%rowdist(k)+1:CGM%Ad%rowdist(k+1)) = col0
                CGM%B%col(CGM%B%rowdist(k)+1:CGM%B%rowdist(k+1)) = col1 
                CGM%E%col(CGM%E%rowdist(k)+1:CGM%E%rowdist(k+1)) = col2 
                k = k + 1
             enddo 
             deallocate(col0,col1,col2)
          else
             ! fluid-solid pts
             ll = 0; nn = 0
             do j = 1,lsiz0 
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.0) then
                   ll = ll + 1
                endif
                if (CGM%vnum(location).eq.6) nn = nn + 1
             enddo
                
             allocate(col0(nn*3+ll*3),col2(nn))
             l = 0; nn = 0
             do j = 1,lsiz0
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                do n = 1,3
                   if (CGM%pnum(location).eq.0.or.&
                       CGM%vnum(location).eq.6) then
                      l = l + 1
                      col0(l) = CGM%vstt(location) + n
                   endif
                enddo
                if (CGM%vnum(location).eq.6) then
                   nn = nn + 1
                   col2(nn) = CGM%pstt(location) + 1
                   if (CGM%pstt(location).lt.0) then
                      print*, 'Error'
                      stop
                   endif  
                endif
             enddo
             call simplessort(col0,idtmp); deallocate(idtmp)
             call simplessort(col2,idtmp); deallocate(idtmp)
             !print*, col2            
 
             allocate(col1(nn+ll))
             do n = 1,3 ! solid part
                l = 0
                do j = 1,lsiz0
                   m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                   call findorder(m,unstrM%Cvlist,location)
                   if (CGM%pnum(location).eq.0.or.&
                       CGM%vnum(location).eq.6) then
                      l = l + 1
                      col1(l) = CGM%vstt(location) + n
                   endif
                enddo
                call simplessort(col1,idtmp); deallocate(idtmp)
                CGM%Ad%col(CGM%Ad%rowdist(k)+1:CGM%Ad%rowdist(k+1)) = col0 
                CGM%B%col(CGM%B%rowdist(k)+1:CGM%B%rowdist(k+1)) = col1 
                CGM%E%col(CGM%E%rowdist(k)+1:CGM%E%rowdist(k+1)) = col2
                !CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + (nn+ll)*3
                !CGM%B%rowdist(k+1)  =  CGM%B%rowdist(k) + (nn+ll)
                !CGM%E%rowdist(k+1)  =  CGM%E%rowdist(k) + nn
                k = k + 1
             enddo
             deallocate(col0,col1,col2)
            
             allocate(col0(lsiz0-ll),col1((lsiz0-ll+nn)*3))
             l = 0; n = 0
             do j = 1,lsiz0 
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.1) then
                   l = l + 1
                   col0(l) = CGM%pstt(location) + 1
                   do mm = 1,CGM%vnum(location)
                      n = n + 1
                      col1(n) = CGM%vstt(location) + mm
                   enddo 
                endif
             enddo
             call simplessort(col0,idtmp); deallocate(idtmp) 
             CGM%Ap%col(CGM%Ap%rowdist(kk)+1:CGM%Ap%rowdist(kk+1)) = col0
             call simplessort(col1,idtmp); deallocate(idtmp)
             CGM%ET%col(CGM%ET%rowdist(kk)+1:CGM%ET%rowdist(kk+1)) = col1
             !print*, col1
             !CGM%Ap%rowdist(kk+1) = CGM%Ap%rowdist(kk) + (lsiz0-ll)
             !CGM%ET%rowdist(kk+1) = CGM%ET%rowdist(kk) + (lsiz0-ll+nn)*3
             kk = kk + 1
             deallocate(col0,col1)

             allocate(col0((lsiz0-ll)*3),col2(lsiz0-ll))
             l = 0; n = 0
             do j = 1,lsiz0
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                call findorder(m,unstrM%Cvlist,location)
                if (CGM%pnum(location).eq.1) then
                   do mm = 1,3
                      n = n + 1
                      col0(n) = CGM%vstt(location) + mm &
                              + CGM%vnum(location) - 3
                   enddo 
                   l = l + 1
                   col2(l) = CGM%pstt(location) + 1
                endif
             enddo
             call simplessort(col0,idtmp); deallocate(idtmp)
             call simplessort(col2,idtmp); deallocate(idtmp)

             allocate(col1(lsiz0-ll))
             do n = 4,6 ! fluid part
                l = 0
                do j = 1,lsiz0
                   m = unstrM%new%v2v(unstrM%new%v2vdist(i)+j)
                   call findorder(m,unstrM%Cvlist,location)
                   if (CGM%pnum(location).eq.1) then
                      l = l + 1
                      col1(l) = CGM%vstt(location) + n - 3&
                              + CGM%vnum(location) - 3
                   endif
                enddo
                call simplessort(col1,idtmp); deallocate(idtmp)

                CGM%Ad%col(CGM%Ad%rowdist(k)+1:CGM%Ad%rowdist(k+1)) = col0 
                CGM%B%col(CGM%B%rowdist(k)+1:CGM%B%rowdist(k+1)) = col1 
                CGM%E%col(CGM%E%rowdist(k)+1:CGM%E%rowdist(k+1)) = col2
                !CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + (lsiz0-ll)*3
                !CGM%B%rowdist(k+1)  =  CGM%B%rowdist(k) + (lsiz0-ll)
                !CGM%E%rowdist(k+1)  =  CGM%E%rowdist(k) + (lsiz0-ll) 
                k = k + 1
             enddo
             deallocate(col0,col1,col2)
          endif
       else
          if (CGM%Ad%rnum(i).eq.3) then
             ! pure solid
             !ll = 0
             !do j = 1,lsiz0 
             !   m = unstrM%v2v(unstrM%v2vdist(i)+j)
             !   call findorder(m,unstrM%Cvlist,location)
             !   if (CGM%pnum(location).eq.1) ll = ll + 1
             !enddo
             
             allocate(col0(lsiz0*3),col1(lsiz0))
             do l = 1,lsiz0
                m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
                call findorder(m,unstrM%Cvlist,location)
                do n = 1,3
                   col0((l-1)*3+n) = CGM%vstt(location) + n
                enddo
             enddo
             call simplessort(col0,idtmp); deallocate(idtmp) 

             do j = 1,3
                do l = 1,lsiz0
                   m = unstrM%new%v2v(unstrM%new%v2vdist(i)+l)
                   call findorder(m,unstrM%Cvlist,location)
                   col1(l) = CGM%vstt(location) + j
                enddo
                call simplessort(col1,idtmp); deallocate(idtmp) 

                CGM%Ad%col(CGM%Ad%rowdist(k)+1:CGM%Ad%rowdist(k+1)) = col0 
                CGM%B%col(CGM%B%rowdist(k)+1:CGM%B%rowdist(k+1)) = col1 

                !CGM%Ad%rowdist(k+1) = CGM%Ad%rowdist(k) + lsiz0*3
                !CGM%B%rowdist(k+1)  = CGM%B%rowdist(k)  + lsiz0
                !CGM%E%rowdist(k+1)  = CGM%E%rowdist(k)  + 0
                k = k + 1
             enddo 
             deallocate(col0,col1)
          else
             ! fluid-solid pts
             print*, 'Error: pure solid pts',i,int(unstrM%rank,2)
          endif
       endif
    enddo 

    call mpi_barrier(unstrM%comm,ierr)
    print*,'local Ap NNZ',int(CGM%Ap%NNZ,4),'at rank',int(unstrM%rank,2)

  end subroutine matrixstruct_general



end module cg_create_matrix_mod

