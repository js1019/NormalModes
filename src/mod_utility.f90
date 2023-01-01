!************************************************************************!
!* This module contains some utility subroutines                         !
!* It was builded, reorganized and developed by Jia Shi                  ! 
!************************************************************************!

!************************************************************************!
!* most the subroutines here have complete input/output                  !
!* general functions can be added here                                   !
!************************************************************************!
 
module utility_mod
  use mpi
  use omp_lib
  use IFPORT
  use para_mod,     only: rkind,pin
  

  implicit none
  !general tools
  public               :: writelog,openlogfile,closelogfile,report_time
  !matrix tools
  public               :: matinv,matdet, init_random_seed
  !basis functions
  public               :: JacobiGL,combination,JacobiP,Basis3D,&
                          Vandermonde3D,dJacobiP,GradBasis3D,GradVandermonde3D,&
                          Basis2D,Vandermonde2D,blend_nodes

  public               :: pickunique,ssort,findorder,simplepickuniord,&
                          ssort_real,simplepickunique,simplessort,findordernostop
 
  public               :: circumcircle,checkifintriangle,realmcolupdate 
                         
  integer              :: logid
  real*8               :: prog_start_time, current_time


 contains
!-------------------------------------------------------------------
  subroutine report_time(comm)
    implicit none
    character (len=8)   :: msg
    integer             :: time_array(8)
    real*8              :: time_elapsed
    real*4              :: time0,time1,dt(2)
    integer             :: IERR,rank,comm

    call mpi_comm_rank(comm, rank, ierr)
    if (rank.eq.0) then
       if (logid.eq.1.or.logid.eq.0) then
         !call date_and_time(values=time_array)
         !prog_start_time = real(time_array(3),8)*3600.0D0*24.0D0 &
         !                + real(time_array(5),8)*3600.0D0        &
         !                + real(time_array(6),8)*60.0D0          &
         !                + real(time_array(7),8)                 &
         !                + real(time_array(8),8)*0.001D0
         !call cpu_time(prog_start_time)
         !time0 = dtime(dt); 
         prog_start_time = omp_get_wtime() 
         write(pin%logfid,*) " Start program at time elapsed = 0.0 second."
         print*, "Start program at time elapsed = 0.0 second."
       else
         !call date_and_time(values=time_array)
         !current_time    = real(time_array(3),8)*3600.0D0*24.0D0 &
         !                + real(time_array(5),8)*3600.0D0        &
         !                + real(time_array(6),8)*60.0D0          &
         !                + real(time_array(7),8)                 &
         !                + real(time_array(8),8)*0.001D0
         !call cpu_time(current_time)
         !time_elapsed = current_time - prog_start_time
         !write(pin%logfid,*) " Time elapsed = ",time_elapsed,"seconds."
         !print*, "Time elapsed = ",time_elapsed,"seconds."
         !time1 = dtime(dt)
         time_elapsed = omp_get_wtime()
         current_time = time_elapsed - prog_start_time         
         write(pin%logfid,*) " Time elapsed = ", current_time,"seconds."
         print*, "Time elapsed = ", current_time, "seconds."
       endif
    endif
  end subroutine report_time

!-------------------------------------------------------------------
  subroutine report_loc_time(comm,time_elapsed)
    implicit none
    integer             :: time_array(8)
    real*8              :: time_elapsed
    integer             :: IERR,rank,comm

    call mpi_comm_rank(comm, rank, ierr)
    if (rank.eq.0) then
       if (logid.eq.1.or.logid.eq.0) then
         call date_and_time(values=time_array)
         prog_start_time = real(time_array(3),8)*3600.0D0*24.0D0 &
                         + real(time_array(5),8)*3600.0D0        &
                         + real(time_array(6),8)*60.0D0          &
                         + real(time_array(7),8)                 &
                         + real(time_array(8),8)*0.001D0
         !print*, "Start program at time elapsed = 0.0 second."
       else
         call date_and_time(values=time_array)
         current_time    = real(time_array(3),8)*3600.0D0*24.0D0 &
                         + real(time_array(5),8)*3600.0D0        &
                         + real(time_array(6),8)*60.0D0          &
                         + real(time_array(7),8)                 &
                         + real(time_array(8),8)*0.001D0
         time_elapsed = current_time - prog_start_time
         !print*, current_time, time_elapsed
         !print*, "Time elapsed = ",time_elapsed,"seconds."
       endif
    endif
  end subroutine report_loc_time


!-------------------------------------------------------------------
  subroutine openlogfile(comm)
    implicit none
    character(len=1024) :: logfname,str
    logical             :: isopen,alive
    integer             :: ierr,rank,comm,nprc
    
    !print*,comm,mpi_comm_world
    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_rank(comm, nprc, ierr)

    write(str,*) nprc

    if (rank.eq.0) then
       logfname=trim(adjustl(pin%outputdir))//trim(adjustl(pin%basename))&
                //'_'//trim(adjustl(str))//'.log'
       inquire(file=trim(logfname),opened=isopen)
       if (.not. isopen) then
         logid = 0
         inquire(file=trim(logfname),exist=alive)
         if (alive) then
           open(pin%logfid, file=trim(logfname), status="old", action="write")
         else
           open(pin%logfid, file=trim(logfname), status="new", action="write")
         endif
       endif
    endif
  end subroutine openlogfile

!-------------------------------------------------------------------
  subroutine writelog(msg,comm)
    implicit none
    character (len=*)   :: msg
    integer             :: IERR,rank,comm
    call mpi_comm_rank(comm, rank, ierr)
    if (rank.eq.0) then
      logid=logid+1
      write(pin%logfid,'(I4". ", A)') logid, msg
    endif
  end subroutine writelog

!-------------------------------------------------------------------
  subroutine closelogfile(comm)
    integer             :: IERR,rank,comm
    call mpi_comm_rank(comm, rank, ierr)
    if (rank.eq.0) then
      close(pin%logfid)
    endif
  end subroutine closelogfile

!-------------------------------------------------------------------
  subroutine matdet(a,det)
    implicit none
    real(kind=rkind)        :: a(3,3), det
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end subroutine matdet

!-------------------------------------------------------------------
  !============================================================
  ! Downloaded from
  ! http://www.webpages.uidaho.edu/~gabrielp/ME549-CE546/matrix-
  ! inverse.pdf
  !============================================================
  subroutine matinv(a,b,n)
  ! subroutine to calculate the inverse of a matrix using Gauss-Jordan 
  ! elimination the inverse of matrix a(n,n) is calculated and stored 
  ! in the matrix b(n,n)
    implicit none
    integer           :: i,j,k,l,m,n,irow
    real(kind=rkind)  :: big,a(n,n),b(n,n),dum,tmp(n,n)
    tmp=a;
    !build the identity matrix
    do i = 1,n
       do j = 1,n
          b(i,j) = 0.0
       end do
       b(i,i) = 1.0
    end do
    do i = 1,n 
    ! this is the big loop over all the columns of a(n,n)
    ! in case the entry a(i,i) is zero, we need to find a
    ! good pivot; this pivot is chosen as the largest value on 
    ! the column i from a(j,i) with j = 1,n
        big = abs(a(i,i))
        do j = i,n
            if (abs(a(j,i)).gt.big) then
                big = abs(a(j,i))
                irow = j
            end if
        end do
        ! interchange lines i with irow for both a() and 
        ! b() matrices
        if (big.gt.abs(a(i,i))) then
            do k = 1,n
                dum = a(i,k) ! matrix a()
                a(i,k) = a(irow,k)
                a(irow,k) = dum
                dum = b(i,k) ! matrix b()
                b(i,k) = b(irow,k)
                b(irow,k) = dum
            end do
        end if
        ! divide all entries in line i from a(i,j) by the 
        ! value a(i,i); same operation for the identity 
        ! matrix
        dum = a(i,i)
        do j = 1,n
            a(i,j) = a(i,j)/dum
            b(i,j) = b(i,j)/dum
        end do
        ! make zero all entries in the column a(j,i); same 
        ! operation for indent()
        do j = i+1,n
            dum = a(j,i)
            do k = 1,n
                a(j,k) = a(j,k) - dum*a(i,k)
                b(j,k) = b(j,k) - dum*b(i,k)
            end do
        end do
    end do
    ! substract appropiate multiple of row j from row j-1 
    do i = 1,n-1
        do j = i+1,n
            dum = a(i,j)
            do l = 1,n
                a(i,l) = a(i,l)-dum*a(j,l)
                b(i,l) = b(i,l)-dum*b(j,l)
            end do
        end do
    end do
    a=tmp
  end subroutine matinv


  subroutine JacobiGL(gaussX,alpha,beta,N)
    real(kind=rkind)                :: alpha,beta,gaussX(N+1),alphat,betat
    real(kind=rkind),allocatable    :: J(:,:)
    integer                         :: i, i1, N, LDA, INFO, LWORK
    double precision                :: WORK(512)
    real(kind=rkind)                :: swork(512)

    gaussX=0.0D0
    gaussX(1)=-1.0D0; gaussX(N+1)=1.0D0
    alphat=alpha+1.D0; betat=beta+1.D0
    if(N.eq.2) gaussX(2)=-(alphat-betat)/(alphat+betat+2.D0)
    if(N.gt.2) then
       allocate(J(N-1,N-1))
       J = 0.0D0
       J(1,1)=(alphat**2-betat**2)/(2.0D0+alphat+betat)/(alphat+betat)
       if (J(1,1).lt.1D-10) J(1,1)=0.0D0
       do i=1,N-2,1
           J(i+1,i+1)=-0.5D0*(alphat**2-betat**2) &
             /(2.D0*real(i+1,8)+alphat+betat)/(2.D0*real(i,8)+alphat+betat)*2.D0
           J(i,i+1)=2.D0/(2.D0*real(i,8)+alphat+betat)*dsqrt(real(i,8)*&
             (real(i,8)+alphat+betat)*(real(i,8)+alphat)*(real(i,8)+betat) &
             /(2.D0*real(i,8)-1.D0+alphat+betat)/(2.D0*real(i,8)+1.D0+alphat+betat))
       enddo
       LDA   = N-1
       LWORK = -1
       
       ! Computing eigenvalue of matrix J
       if (rkind <= 4) then 
          call SSYEV('N','U',N-1,J,LDA,gaussX(2:N),SWORK,LWORK,INFO)
          LWORK = min(512,int(SWORK(1)))
          call SSYEV('N','U',N-1,J,LDA,gaussX(2:N),SWORK,LWORK,INFO)
       else
          call DSYEV('N','U',N-1,J,LDA,gaussX(2:N),WORK,LWORK,INFO)
          LWORK = min(512,int(WORK(1)))
          call DSYEV('N','U',N-1,J,LDA,gaussX(2:N),WORK,LWORK,INFO)
       endif
    endif
  end subroutine JacobiGL
  
  !--------------------------------------------------------------------------
  function combination(n,alpha)
      integer                           :: n,alpha,beta,i
      real(kind=rkind)                  :: combination
      real(kind=rkind)                  :: tmp
      beta = n-alpha
      tmp  = 1.0D0
      do i=1,beta
         tmp=tmp*real(alpha+i,8)/real(i,8)
      enddo
      combination=tmp
      return
  end function combination
  
  !--------------------------------------------------------------------------
  function JacobiP(x,alpha,beta,N)
      integer             :: alpha,beta,N,i,hh
      real(kind=rkind)    :: JacobiP,x,gamma0,gamma1,a1,a2,a3,p1,p2,ir8,hh8

      gamma0=2.0D0**(alpha+beta+1)/real(alpha+beta+1,8) &
          /combination(alpha+beta,alpha)
      !*gamma(alpha+1.0)*gamma(beta+1.0)/gamma(alpha+beta+1.0)
      !/combination(alpha+beta,alpha)
      gamma1=real(alpha+1,8)*real(beta+1,8)/real(alpha+beta+3,8)*gamma0
      if(N.eq.0)then
          JacobiP=1.0D0/sqrt(gamma0)
      elseif(N.eq.1)then
          JacobiP=(real(alpha+beta+2,8)*x/2.D0+real(alpha-beta,8)/2.D0)&
              /sqrt(gamma1)
      else
          p1=1.0D0/sqrt(gamma0)
          p2=(real(alpha+beta+2,8)*x/2.D0+real(alpha-beta,8)/2.D0)/sqrt(gamma1)
          a1=2.0D0/real(2+alpha+beta,8)*sqrt(real(alpha+1,8)*real(beta+1,8) &
              /real(alpha+beta+3,8))
          do i=1,N-1,1
              hh=2*i+alpha+beta
              ir8 = real(i,8); hh8 = real(hh,8)
              a2=2.0D0/(hh8+2.D0)*sqrt((ir8+1.0D0)*(ir8+1.D0+real(alpha+beta,8)) &
                *real(i+1+alpha,8)*real(i+1+beta,8)/(hh8+1.D0)/(hh8+3.D0))
              a3=-real(alpha**2-beta**2,8)/hh8/(hh8+2.D0)
              JacobiP=1.0D0/a2*(-a1*p1+(x-a3)*p2)
              a1=a2;p1=p2;p2=JacobiP
          enddo
      endif
      return
  end function JacobiP
  
  !--------------------------------------------------------------------------
  subroutine Basis3D(p,r,s,t,xdim,i,j,k)
    ! input: xdim, r, s, t
    integer             :: i,j,k,xdim,n
    real(kind=rkind)    :: p(xdim),r(xdim),s(xdim),t(xdim),a,b,c,h1,h2,h3
    do n=1,xdim
       if(abs(s(n)+t(n)).le.pin%TOL)then
           a=-1.0D0
       else
           a=2D0*(1.0D0+r(n))/(-s(n)-t(n))-1.D0
       endif
       if(abs(t(n)-1).le.pin%TOL)then
           b=-1.D0
       else
           b=2.D0*(1.0D0+s(n))/(1.0D0-t(n))-1.D0
       endif
       c  = t(n)
       h1 = JacobiP(a,0,0,i)
       h2 = JacobiP(b,2*i+1,0,j)
       h3 = JacobiP(c,2*(i+j)+2,0,k)
       p(n)=2.D0*sqrt(2.0D0)*h1*((1.0D0-b)**i)*h2*((1.0D0-c)**(i+j))*h3
  !      print*,'p=',p(n)
    enddo
  end subroutine Basis3D
  
  !--------------------------------------------------------------------------
  subroutine Vandermonde3D(V3D,N,r,s,t,xdim)
    ! input: N: pOrder
    integer          :: N,xdim,i,j,k,sk
    real(kind=rkind) :: V3D(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind) :: r(xdim),s(xdim),t(xdim)
    sk=1;
    do i=0,N,1
        do j=0,N-i,1
            do k=0,N-i-j,1
                call Basis3D(V3D(:,sk),r,s,t,xdim,i,j,k)
  !              print*,V3D(:,sk)
                sk=sk+1
            enddo
        enddo
    enddo
  end subroutine Vandermonde3D
  
  !--------------------------------------------------------------------------
  function dJacobiP(x,alpha,beta,N)
    integer             :: alpha,beta,N
    real(kind=rkind)    :: dJacobiP,x
    if(N.eq.0)then
        dJacobiP=0.0D0
    else
        dJacobiP=sqrt(real(N,8)*(real(N+alpha+beta+1,8)))&
                 *JacobiP(x,alpha+1,beta+1,N-1)
    endif
    return
  end function dJacobiP
  
  !--------------------------------------------------------------------------
  subroutine GradBasis3D(pr,ps,pt,r,s,t,xdim,i,j,k)
    integer            :: i,j,k,xdim,n
    real(kind=rkind)   :: pr(xdim),ps(xdim),pt(xdim),r(xdim),s(xdim),t(xdim)
    real(kind=rkind)   :: a,b,c,h1,h2,h3,dh1,dh2,dh3,tmp

    do n=1,xdim
      if(abs(s(n)+t(n)).le.pin%TOL)then
          a=-1.0D0
      else
          a=2.D0*(1.0D0+r(n))/(-s(n)-t(n))-1.D0
      endif
      if(abs(t(n)-1).le.pin%TOL)then
          b=-1.D0
      else
          b=2.D0*(1.0D0+s(n))/(1.0D0-t(n))-1.D0
      endif
      c=t(n)

      h1=JacobiP(a,0,0,i);dh1=dJacobiP(a,0,0,i)
      h2=JacobiP(b,2*i+1,0,j);dh2=dJacobiP(b,2*i+1,0,j)
      h3=JacobiP(c,2*(i+j)+2,0,k);dh3=dJacobiP(c,2*(i+j)+2,0,k)

      pr(n)=dh1*h2*h3

      if(i.gt.1)pr(n)=pr(n)*((0.5D0*(1.0D0-b))**(i-1))

      if((i+j).gt.1) pr(n)=pr(n)*((0.5d0*(1.0d0-c))**(i+j-1))
      ps(n)=0.5D0*(1.0D0+a)*pr(n)
      tmp=dh2*((0.5D0*(1.0D0-b))**i)

      if(i.gt.0) tmp=tmp+(-0.5D0*real(i,8))*(h2*(0.5D0*(1.0D0-b))**(i-1))

      if((i+j).gt.1) tmp=tmp*((0.5D0*(1.0D0-c))**(i+j-1))
      tmp=tmp*h1*h3
      ps(n)=ps(n)+tmp
      pt(n)=0.5D0*(1.0D0+a)*pr(n)+0.5D0*(1.0D0+b)*tmp
      tmp=dh3*((0.5D0*(1.0D0-c))**(i+j))

      if((i+j).gt.0) then
         tmp=tmp-0.5D0*real(i+j,8)*(h3*((0.5D0*(1.0D0-c))**(i+j-1)))
      endif 

      tmp=h1*h2*tmp*((0.5D0*(1.0D0-b))**i)
      pt(n)=pt(n)+tmp
    enddo
    tmp=2.0D0**(real(2*i+j,8)+1.5D0)
    pr=pr*tmp; ps=ps*tmp; pt=pt*tmp

  end subroutine GradBasis3D
  
  !--------------------------------------------------------------------------
  subroutine GradVandermonde3D(V3Dr,V3Ds,V3Dt,N,r,s,t,xdim)
    integer             :: N,xdim,i,j,k,sk
    real(kind=rkind)    :: V3Dr(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind)    :: V3Ds(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind)    :: V3Dt(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind)    :: r(xdim),s(xdim),t(xdim)
    sk=1;
    do i=0,N,1
        do j=0,N-i,1
            do k=0,N-i-j,1
                call GradBasis3D(V3Dr(:,sk),V3Ds(:,sk),V3Dt(:,sk),&
                                 r,s,t,xdim,i,j,k)
                sk=sk+1
            enddo
        enddo
    enddo
  end subroutine GradVandermonde3D
  
  !--------------------------------------------------------------------------
  subroutine Basis2D(p,r,s,xdim,i,j)
    integer             :: i,j,xdim,n
    real(kind=rkind)    :: p(xdim),r(xdim),s(xdim),a,b,h1,h2
    do n=1,xdim,1
       if(abs(s(n)-1.0).le.pin%TOL)then
           a=-1.0D0
       else
           a=2.0D0*(1.0D0+r(n))/(1.0D0-s(n))-1.0D0
       endif
       b=s(n)
       h1=JacobiP(a,0,0,i)
       h2=JacobiP(b,2*i+1,0,j)
       p(n)=sqrt(2.0D0)*h1*h2*(1.0D0-b)**i
    enddo
  end subroutine Basis2D
  
  !--------------------------------------------------------------------------
  subroutine Vandermonde2D(V2D,N,r,s,xdim)
    integer          :: N,xdim,i,j,sk
    real(kind=rkind) :: V2D(xdim,(N+1)*(N+2)/2)
    real(kind=rkind) :: r(xdim),s(xdim)
    sk=1
    do i=0,N,1
        do j=0,N-i,1
            call Basis2D(V2D(:,sk),r,s,xdim,i,j)
            sk=sk+1
        enddo
    enddo
  end subroutine Vandermonde2D
  

!*********************************************************************************** 
!----------------------------------------------------------------------------------
! This subroutine computes pNp nodal points in reference tetrahedron.
! It is written by Ruichao Ye
!----------------------------------------------------------------------------------
  subroutine blend_nodes(x,y,z,pNp,pOrder)
    integer, intent(in)             :: pNp,pOrder
    real(kind=rkind), intent(out)   :: x(pNp),y(pNp),z(pNp)

    integer                         :: i,j,k,sk
    real(kind=rkind), allocatable   :: r(:),s(:),t(:)
    real(kind=rkind)                :: alphastore(15),alpha
    real(kind=rkind), allocatable   :: L1(:),L2(:),L3(:),L4(:)
    real(kind=rkind), allocatable   :: La(:),Lb(:),Lc(:),Ld(:),XYZ(:,:),shift(:,:)
    real(kind=rkind), allocatable   :: warp1(:),warp2(:),blend(:),denom(:)
    real(kind=rkind)                :: v1(3),v2(3),v3(3),v4(3),t1(3,4),t2(3,4)
    real(kind=rkind)                :: rtmp
    
    data alphastore /0.0D0,   0.0D0,   0.0D0,   0.1002D0,1.1332D0, &
                     1.5608D0,1.3413D0,1.2577D0,1.1603D0,1.10153D0,&
                     0.6080D0,0.4523D0,0.8856D0,0.8717D0,0.9655D0  /

    v1(1)=-1.D0; v1(2)=-1.D0/sqrt(3.D0); v1(3)=-1.D0/sqrt(6.D0)                                                 
    v2(1)= 1.D0; v2(2)=-1.D0/sqrt(3.D0); v2(3)=-1.D0/sqrt(6.D0)
    v3(1)= 0.D0; v3(2)= 2.D0/sqrt(3.D0); v3(3)=-1.D0/sqrt(6.D0)                  
    v4(1)= 0.D0; v4(2)= 0.D0           ; v4(3)= 3.D0/sqrt(6.D0)

    allocate(r(pNp),s(pNp),t(pNp))
    allocate(L1(pNp),L2(pNp),L3(pNp),L4(pNp))
    allocate(La(pNp),Lb(pNp),Lc(pNp),Ld(pNp)) 
    allocate(warp1(pNp),warp2(pNp))
    allocate(blend(pNp),denom(pNp))
    allocate(XYZ(3,pNp),shift(3,pNp))

    if(pOrder .le. 15)then
       alpha=alphastore(pOrder)
    else
       alpha=1.0D0
    endif

    sk=1
    do i=1,pOrder+1,1
       do j=1,pOrder+2-i,1
          do k=1,pOrder+3-i-j,1
             r(sk) = -1.D0 + real(k-1,8)*2.0D0/real(pOrder,8)
             s(sk) = -1.D0 + real(j-1,8)*2.0D0/real(pOrder,8)
             t(sk) = -1.D0 + real(i-1,8)*2.0D0/real(pOrder,8)
             sk    = sk+1
          enddo
       enddo
    enddo
    !print*,r  
 
    L1 =  (1.D0+t)/2.D0;      L2 = (1.D0+s)/2.D0
    L3 = -(1.D0+r+s+t)/2.D0;  L4 = (1.D0+r)/2.D0
    t1(:,1) = v2-v1;            t1(:,2) = v2-v1;
    t1(:,3) = v3-v2;            t1(:,4) = v3-v1;   
    t2(:,1) = v3-0.5D0*(v1+v2); t2(:,2) = v4-0.5D0*(v1+v2);
    t2(:,3) = v4-0.5D0*(v2+v3); t2(:,4) = v4-0.5D0*(v1+v3);  

    ! normalize tangents
    do i=1,4,1
       rtmp=sqrt(sum((t1(:,i)*t1(:,i))))
       if(rtmp.gt.pin%TOL) t1(:,i)=t1(:,i)/rtmp
       rtmp=sqrt(sum((t2(:,i)*t2(:,i))))
       if(rtmp.gt.pin%TOL) t2(:,i)=t2(:,i)/rtmp
    enddo

    do i=1,3,1
       do j=1,pNp,1
          XYZ(i,j)=L3(j)*v1(i)+L4(j)*v2(i) &
                  +L2(j)*v3(i)+L1(j)*v4(i)
       enddo
    enddo
 
    shift=0.0d0

    do i=1,4,1
       if(i .eq. 1)then
         La = L1; Lb = L2; Lc = L3; Ld = L4
       elseif(i .eq. 2)then
         La = L2; Lb = L1; Lc = L3; Ld = L4
       elseif(i .eq. 3)then
         La = L3; Lb = L1; Lc = L4; Ld = L2
       else
         La = L4; Lb = L1; Lc = L3; Ld = L2
       endif

       call evalshift(warp1,warp2,alpha,Lb,Lc,Ld,pNp,pOrder)

       blend=Lb*Lc*Ld
       denom=(Lb+0.5D0*La)*(Lc+0.5D0*La)*(Ld+0.5D0*La)

       do j=1,pNp,1
          if(denom(j)>pin%TOL) then
             blend(j)=(1.D0+(alpha*La(j))**2)*blend(j)/denom(j)
          endif
       enddo

       do j=1,pNp,1
          shift(:,j)=shift(:,j)+blend(j)*warp1(j)*t1(:,i)+&
                     blend(j)*warp2(j)*t2(:,i)
       enddo

       do j=1,pNp,1
          if((La(j).lt.pin%TOL) .and. (.not. ((Lb(j).gt.pin%TOL) .and. &
             (Lc(j).gt.pin%TOL) .and. (Ld(j).gt.pin%TOL))))then
             shift(:,j)=warp1(j)*t1(:,i)+warp2(j)*t2(:,i)
          endif
       enddo
    enddo

    XYZ=XYZ+shift;

    r=XYZ(1,:); s=XYZ(2,:); t=XYZ(3,:)

    x = (v1(1)*(v4(2)*v3(3)-v3(2)*v4(3))+v3(1)*(v1(2)*v4(3)-v4(2)*v1(3)) &
                                        +v4(1)*(v3(2)*v1(3)-v1(2)*v3(3)) &
        +(v1(2)*(v3(3)-v4(3))+v3(2)*(v4(3)-v1(3))+v4(2)*(v1(3)-v3(3)))*r &
        +(v1(1)*(v4(3)-v3(3))+v3(1)*(v1(3)-v4(3))+v4(1)*(v3(3)-v1(3)))*s &
        +(v1(1)*(v3(2)-v4(2))+v3(1)*(v4(2)-v1(2))+v4(1)*(v1(2)-v3(2)))*t)&
        /4.0D0/sqrt(2.0D0);

    y = (v1(2)*(v4(1)*v2(3)-v2(1)*v4(3))+v2(2)*(v1(1)*v4(3)-v4(1)*v1(3)) &
                                        +v4(2)*(v2(1)*v1(3)-v1(1)*v2(3)) &
        +(v1(2)*(v4(3)-v2(3))+v2(2)*(v1(3)-v4(3))+v4(2)*(v2(3)-v1(3)))*r &
        +(v1(1)*(v2(3)-v4(3))+v2(1)*(v4(3)-v1(3))+v4(1)*(v1(3)-v2(3)))*s &
        +(v1(1)*(v4(2)-v2(2))+v2(1)*(v1(2)-v4(2))+v4(1)*(v2(2)-v1(2)))*t)&
         /4.0D0/sqrt(2.0D0);

    z = (v1(3)*(v3(1)*v2(2)-v2(1)*v3(2))+v2(3)*(v1(1)*v3(2)-v3(1)*v1(2)) &
                                        +v3(3)*(v2(1)*v1(2)-v1(1)*v2(2)) &
        +(v1(2)*(v2(3)-v3(3))+v2(2)*(v3(3)-v1(3))+v3(2)*(v1(3)-v2(3)))*r &
        +(v1(1)*(v3(3)-v2(3))+v2(1)*(v1(3)-v3(3))+v3(1)*(v2(3)-v1(3)))*s &
        +(v1(1)*(v2(2)-v3(2))+v2(1)*(v3(2)-v1(2))+v3(1)*(v1(2)-v2(2)))*t)&
        /4.0D0/sqrt(2.0D0);

    x = x*2.D0-1.0D0; y = y*2.D0-1.0D0; z = z*2.D0-1.0D0

    deallocate(XYZ)
    deallocate(La,Lb,Lc,Ld)
    deallocate(L1,L2,L3,L4)
    deallocate(r,s,t)
    deallocate(warp1,warp2)
    deallocate(shift,blend,denom)

  end subroutine blend_nodes

!-----------------------------------------------------------------------------
  subroutine evalshift(w1,w2,alpha,L1,L2,L3,pNp,pOrder)
    integer, intent(in)              :: pNp,pOrder
    real(kind=rkind)                 :: w1(pNp),w2(pNp),alpha
    real(kind=rkind)                 :: L1(pNp),L2(pNp),L3(pNp)
    real(kind=rkind)                 :: T1(pNp),T2(pNp),T3(pNp)
    real(kind=rkind), allocatable    :: gaussX(:),warp1(:),warp2(:),warp3(:)
    real(kind=rkind)                 :: tmp1,tmp2

    allocate(gaussX(pOrder+1))
    allocate(warp1(pNp),warp2(pNp))
    allocate(warp3(pNp))
    tmp1 = 0.0D0; tmp2 = 0.0D0
    call JacobiGL(gaussX,tmp1,tmp2,pOrder)
    gaussX = - gaussX
 
    T1 = L3 - L2
    T2 = L1 - L3
    T3 = L2 - L1

    call evalwarp(warp1,gaussX,T1,pNp,pOrder)
    call evalwarp(warp2,gaussX,T2,pNp,pOrder)
    call evalwarp(warp3,gaussX,T3,pNp,pOrder)

    warp1 = L2*L3*4.D0*warp1*(1.D0+(alpha*L1)**2)
    warp2 = L1*L3*4.D0*warp2*(1.D0+(alpha*L2)**2)
    warp3 = L1*L2*4.D0*warp3*(1.D0+(alpha*L3)**2)
    w1    = warp1-0.5D0*(warp2+warp3)
    w2    = sqrt(3.0D0)/2.D0*(warp2-warp3)
 
    deallocate(gaussX,warp1,warp2,warp3)
  end subroutine evalshift


!-----------------------------------------------------------------------------
  subroutine evalwarp(warp,xnodes,xout,pNp,pOrder)
    integer, intent(in)            :: pNp,pOrder
    real(kind=rkind), intent(in)   :: xnodes(pOrder+1),xout(pNp)
    real(kind=rkind), intent(out)  :: warp(pNp)
    real(kind=rkind), allocatable  :: xeq(:),d(:)
    integer i,j
    
    warp = 0.0D0;
    allocate(xeq(pOrder+1),d(pNp))

    do i = 1,pOrder+1,1
       xeq(i)=-1.D0+2.D0*real(pOrder+1-i,8)/real(pOrder,8)
    enddo

    do i = 1,pOrder+1,1
      d = xnodes(i) - xeq(i)
      do j=2,pOrder,1
        if (i .ne. j) then
          d = d*(xout-xeq(j))/(xeq(i)-xeq(j))
        endif
      enddo

      if (i .ne. 1) then
        d = - d/(xeq(i)-xeq(1))
      endif

      if (i .ne. pOrder+1) then
        d = d/(xeq(i)-xeq(pOrder+1))
      endif
      warp = warp + d
    enddo

    deallocate(xeq,d)
  end subroutine evalwarp

!*******************************************************************************

   
!---------------------------------------------------------------------------
!**********************************************!
! By Jia Shi
! the subroutines below are for CG formulation
!**********************************************!  
!--------------------------------------------------------------  
  subroutine comparetwovects(ids,vec1,vec2,len)
      ! compare two integer vectors for subroutine trans_label
      ! if vec1(i) is less then vec2(i) return 1; else return 0
      integer,intent(in)                  :: len
      integer,intent(in)                  :: vec1(len), vec2(len)
      integer,intent(inout)               :: ids(len) 

      integer                             :: i
      
      ids = 0
      do i = 1, len
         if (vec1(i) < vec2(i)) then
            ids(i) = 1
         endif
      enddo
  
  end subroutine comparetwovects

!--------------------------------------------------------------  
  subroutine pickunique(pos,siz,acol,Nentry)
      integer, intent(in)                  :: siz
      integer*8, intent(in)                :: pos(siz) 
      integer*8, intent(out)               :: Nentry
      integer*8, allocatable, intent(out)  :: acol(:)
      
      integer                              :: i,kk
      integer*8                            :: Ntmp,ntmp2
      integer                              :: nc(siz)


      nc = 0; kk = 1; i = 1; Ntmp = pos(1)
      do while(i <= siz)
         nc(kk) = pos(i)
         Ntmp   = pos(i)
         do while (i < siz) 
            if (pos(i+1) .eq. ntmp) then
               i = i + 1
            else
               exit
            endif
         enddo
         i = i + 1; kk = kk + 1
      enddo
            
      Nentry = kk - 1
      allocate(acol(Nentry))
      acol   = nc(1:Nentry)

  end subroutine pickunique

!--------------------------------------------------------------  
  subroutine simplepickunique(pos,siz,acol,Nentry)
      integer, intent(in)                  :: siz
      integer, intent(in)                  :: pos(siz) 
      integer, intent(out)                 :: Nentry
      integer, allocatable, intent(out)    :: acol(:)
      
      integer                              :: i,kk
      integer                              :: Ntmp
      integer                              :: nc(siz)
     
      nc = 0; kk = 1; i = 1; Ntmp = pos(1)
      do while(i <= siz)
         nc(kk) = pos(i)
         Ntmp   = pos(i)
         do while (i < siz) 
            if (pos(i+1) .eq. ntmp) then
               i = i + 1
            else
               exit
            endif
         enddo
         i = i + 1; kk = kk + 1
      enddo
            
      Nentry = kk - 1
      allocate(acol(Nentry))
      acol   = nc(1:Nentry)

  end subroutine simplepickunique

!--------------------------------------------------------------  
  subroutine simplepickuniord(pos,siz,acol,Nentry)
      integer, intent(in)                  :: siz
      integer, intent(in)                  :: pos(siz) 
      integer, intent(out)                 :: Nentry
      integer, allocatable, intent(out)    :: acol(:)
      
      integer                              :: i,kk
      integer                              :: Ntmp
      integer                              :: nc(siz)
     
      nc = 0; kk = 1; i = 1; Ntmp = pos(1)
      do while(i <= siz)
         nc(kk) = i !pos(i)
         Ntmp   = pos(i)
         do while (i < siz) 
            if (pos(i+1) .eq. ntmp) then
               i = i + 1
            else
               exit
            endif
         enddo
         i = i + 1; kk = kk + 1
      enddo
            
      Nentry = kk - 1
      allocate(acol(Nentry))
      acol   = nc(1:Nentry)

  end subroutine simplepickuniord

!--------------------------------------------------------------  
  subroutine simplessort(a,ID)
      ! reorder array a(:) 
      integer, intent(inout)                 :: a(:)
      integer, allocatable                   :: ID(:)
      integer                                :: N,k,gap,Na,Nb,i,j,itmp
      integer                                :: tem
      
      N  = size(a)
      allocate(ID(N))
      ID = (/(i,i=1,N)/)
      k = 1; gap = 0

      do while(gap .lt. N)
         gap=gap*3+1
      enddo

      do while(gap .gt. 0)
         do i = gap,N
            tem = a(i); itmp = ID(i)
            j   = i - gap
            do while(j .gt. 0)
               if (a(j) .le. tem) exit
               a(j+gap)=a(j); ID(j+gap)=ID(j); j=j-gap
            enddo
            a(j+gap)  = tem
            ID(j+gap) = itmp
         enddo
         gap = (gap-1)/3
      enddo
      !print*, a(:)

  end subroutine simplessort


!--------------------------------------------------------------  
  subroutine ssort(a,ID)
      ! reorder a(:) here a is long integer 
      integer*8, intent(inout)               :: a(:)
      integer, allocatable                   :: ID(:)
      integer                                :: N,k,gap,Na,Nb,i,j,itmp
      integer*8                              :: tem
      
      N  = size(a)
      allocate(ID(N))
      ID = (/(i,i=1,N)/)
      k = 1; gap = 0

      do while(gap .lt. N)
         gap=gap*3+1
      enddo

      do while(gap .gt. 0)
         do i = gap,N
            tem = a(i); itmp = ID(i)
            j   = i - gap
            do while(j .gt. 0)
               if (a(j) .le. tem) exit
               a(j+gap)=a(j); ID(j+gap)=ID(j); j=j-gap
            enddo
            a(j+gap)  = tem
            ID(j+gap) = itmp
         enddo
         gap = (gap-1)/3
      enddo
  end subroutine ssort


!--------------------------------------------------------------  
  subroutine ssort_real(a,id)
      ! reorder real a(:)  
      real(kind=rkind), intent(inout)        :: a(:)
      integer, intent(inout)                 :: ID(:)
      integer                                :: N,k,gap,Na,Nb,i,j,itmp
      real(kind=rkind)                       :: tem
      
      N  = size(a)
      k = 1; gap = 0

      do while(gap .lt. N)
         gap=gap*3+1
      enddo

      do while(gap .gt. 0)
         do i = gap,N
            tem = a(i); itmp = ID(i)
            j   = i - gap
            do while(j .gt. 0)
               if (a(j) .le. tem) exit
               a(j+gap)=a(j); ID(j+gap)=ID(j); j=j-gap
            enddo
            a(j+gap)  = tem
            ID(j+gap) = itmp
         enddo
         gap = (gap-1)/3
      enddo
 
  end subroutine ssort_real

!--------------------------------------------------------------  
  subroutine colupdate(mat1,col,mat2,siz)
    ! compute mat2 = diag(col)*mat1 
    integer, intent(in)                     :: siz
    complex(kind=rkind), intent(in)         :: mat1(siz,siz)
    complex(kind=rkind), intent(in)         :: col(siz) 
    complex(kind=rkind), intent(out)        :: mat2(siz,siz)

    integer                                 :: i

    do i = 1,siz
       mat2(:,i) = mat1(:,i)*col
    enddo

  end subroutine colupdate 
!-------------------------------------------------------------
  subroutine diagmatrix(mat,vct,siz)
    ! make vct into a diagonal matrix 
    integer, intent(in)                    :: siz
    real(kind=rkind), intent(in)           :: vct(siz)
    real(kind=rkind), intent(out)          :: mat(siz,siz)

    integer                                :: i
   
    mat = 0.0D0 
    do i = 1,siz
       mat(i,i) = vct(i)
    enddo
      
  end subroutine diagmatrix
!--------------------------------------------------------------  
  subroutine realmcolupdate(mat1,col,mat2,siz)
    ! compute mat2 = diag(col)*mat1 
    integer, intent(in)                     :: siz
    real(kind=rkind), intent(in)            :: mat1(siz,siz)
    real(kind=rkind), intent(in)            :: col(siz) 
    real(kind=rkind), intent(out)           :: mat2(siz,siz)

    integer                                 :: i

    do i = 1, siz
       mat2(:,i) = mat1(:,i)*col
    enddo

  end subroutine realmcolupdate 
!--------------------------------------------------------------  
  subroutine findorder(aim,ntmp,loc)
     ! find id of aim from array ntmp
     integer, intent(in)                    :: aim,ntmp(:)
     integer, intent(out)                   :: loc
    
     integer                                :: siz,up,down 
   
 
     loc = - 1

     siz = size(ntmp) 
     if (siz.le.0) then 
        print*, 'error siz',siz  
     endif
 
     up = 1; down = siz; siz = siz/2  
     do while (aim /= ntmp(siz) .and. down > up)
        if (aim < ntmp(siz)) then  
           if (siz == down) then
              siz = siz -1; down = down -1
           else
              down = siz; siz = up+(down-up+1)/2
           endif
        elseif (aim > ntmp(siz)) then 
           if (siz == up) then 
              siz = siz + 1; up = up + 1
           else
              up = siz; siz = down-(down-up+1)/2
           endif  
        endif
        !print*, ntmp(up),aim,ntmp(down)
     enddo
   
     if (ntmp(siz) == aim) then
        loc = siz
     elseif (ntmp(up) == aim .and. up == down) then 
        loc = up
     else
        loc = -1 ! aim is not in ntmp
        print*, "error: can not find the id",aim!,ntmp
        stop
        !pause
     endif
     
  end subroutine findorder

  ! JS 01302018
  subroutine findordernostop(aim,ntmp,loc)
     ! find id of aim from array ntmp
     integer, intent(in)                    :: aim,ntmp(:)
     integer, intent(out)                   :: loc
    
     integer                                :: siz,up,down 
    
     loc = - 1

     siz = size(ntmp) 
     if (siz.le.0) then 
        print*, 'error siz',siz  
     endif
 
     up = 1; down = siz; siz = siz/2  
     do while (aim /= ntmp(siz) .and. down > up)
        if (aim < ntmp(siz)) then  
           if (siz == down) then
              siz = siz -1; down = down -1
           else
              down = siz; siz = up+(down-up+1)/2
           endif
        elseif (aim > ntmp(siz)) then 
           if (siz == up) then 
              siz = siz + 1; up = up + 1
           else
              up = siz; siz = down-(down-up+1)/2
           endif  
        endif
        !print*, ntmp(up),aim,ntmp(down)
     enddo
   
     if (ntmp(siz) == aim) then
        loc = siz
     elseif (ntmp(up) == aim .and. up == down) then 
        loc = up
     else
        loc = -1 ! aim is not in ntmp
     endif
     
  end subroutine findordernostop


!--------------------------------------------------------------  
  subroutine circumcircle(abc,r,cen)
     ! compute the radius and center of circumcircle
     ! this is used for box boundary only 
     real(kind=rkind), intent(in)          :: abc(2,3)
     real(kind=rkind), intent(out)         :: r,cen(2)
  
     real(kind=rkind)                      :: r1,r2,r3,vc(2),B(2,2),Binv(2,2)
  
     !vc(1) = norm2(abc(:,1))**2-norm2(abc(:,2))**2
     !vc(2) = norm2(abc(:,2))**2-norm2(abc(:,3))**2   
     vc(1) = sum((abc(:,1))**2)-sum((abc(:,2))**2)
     vc(2) = sum((abc(:,2))**2)-sum((abc(:,3))**2)   

     B(1,1)   = 2.0D0*(abc(1,1)-abc(1,2))
     B(1,2)   = 2.0D0*(abc(2,1)-abc(2,2))
     B(2,1)   = 2.0D0*(abc(1,2)-abc(1,3))
     B(2,2)   = 2.0D0*(abc(2,2)-abc(2,3))
 
     call matinv(B,Binv,2)
     cen      = matmul(Binv,vc)
     
     r1       = sqrt(sum((cen-abc(:,1))**2))
     r2       = sqrt(sum((cen-abc(:,2))**2))
     r3       = sqrt(sum((cen-abc(:,3))**2))

     r        = max(r1,r2,r3)
     ! check 
     if ((r-min(r1,r2,r3))/r > 5.0D-3) then
        print*, "warning: circumcircle is not accurate",(r-min(r1,r2,r3))/r,&
                "suggestion: use double precision"
     endif

  end subroutine circumcircle
  

!--------------------------------------------------------------  
  subroutine checkifintriangle(check,pts,abc)
     ! check if a pt is in a triangle or not
     ! this is used for box boundary only 
     real(kind=rkind), intent(in)          :: pts(2),abc(2,3)

     integer, intent(out)                  :: check
 
     integer                               :: i,j,k
     real(kind=rkind)                      :: PI,vet(2,3),theta(3)

     check = 0; PI = 3.14159265359

     do i = 1,3
        vet(:,i) = abc(:,i) - pts(:) 
     enddo

     do i = 1,3
        !if (norm2(vet(:,i)) <= pin%TOL) then
        if (sqrt(sum((vet(:,i))**2)) <= pin%TOL) then
           check = 1
           return
        else
           !vet(:,i) = vet(:,i)/norm2(vet(:,i))
           vet(:,i) = vet(:,i)/sqrt(sum((vet(:,i))**2))
        endif
     enddo
     
     theta(1) = dot_product(vet(:,1),vet(:,2))
     theta(2) = dot_product(vet(:,2),vet(:,3))
     theta(3) = dot_product(vet(:,3),vet(:,1))
 
     
     do i = 1,3
        if (rkind == 8) then
           theta(i) = dacos(theta(i))
        else
           theta(i) = acos(theta(i))
        endif
     enddo

                 

     if (dabs(sum(theta)-2.0*PI)<1.D-4.and.rkind==8) then
        check = 1
     elseif (abs(sum(theta)/2.0/PI-1.0)<1.0E-4.and.rkind==4) then
        check = 1
     elseif (minval(PI-theta(:)) < 1.0E-4 .and. minval(PI -theta(:)) > 0.0E0) then 
        check = 1 
     else
        check = 0
     endif 
           
     !print*, sum(theta)/2/PI, check

  end subroutine checkifintriangle
 
  !------------------------------------------------
  subroutine pnm_data_exchange(nvtx,slgth,vpids0,dist0,sdat,&
                               lsiz1,vtmp1,vlist1,comm)
    ! apply mpi_alltoall, mpi_alltoallv
    ! input: nvtx,vpids0 -- number of vtx and their proc ids
    !        dist0 -- their distributation
    !        sdat  -- send out information
    integer, intent(in)                           :: nvtx,slgth,comm
    integer, intent(in)                           :: vpids0(nvtx)
    integer, intent(in)                           :: dist0(nvtx+1)
    integer, intent(in)                           :: sdat(slgth)
    integer, intent(out)                          :: lsiz1
    integer, intent(out), allocatable             :: vtmp1(:),vlist1(:)

    integer                                       :: i,j,k,l,ierr
    integer                                       :: rank,nproc,lsiz0
    integer, allocatable, dimension(:)            :: vtmp0
    integer, allocatable, dimension(:)            :: cnt0,cnt1
    integer, allocatable, dimension(:)            :: vlist0

    call mpi_comm_size(comm,nproc,ierr)
    call mpi_comm_rank(comm, rank,ierr)

    allocate(vtmp0(nproc),vtmp1(nproc)); vtmp0 = 0
    do i = 1,nvtx
       j = vpids0(i) + 1
       vtmp0(j) = vtmp0(j) + dist0(i+1)-dist0(i) 
    enddo
    !print*,vtmp0

    if (sum(vtmp0).ne.slgth) then
       print*,'Error, dist0 and slgth do not match',&
             sum(vtmp0),slgth,' at rank',int(rank,4)
    endif

    call MPI_ALLTOALL(vtmp0,1,mpi_integer,&
                      vtmp1,1,mpi_integer,comm,ierr)
    !print*, vtmp1
    allocate(cnt0(0:nproc-1),cnt1(0:nproc-1))

    cnt0(0) = 0; cnt1(0) = 0
    do i = 1,nproc-1
       cnt0(i) = cnt0(i-1) + vtmp0(i)
       cnt1(i) = cnt1(i-1) + vtmp1(i)
    enddo

    lsiz0 = sum(vtmp0); allocate(vlist0(lsiz0))
    vlist0 = 0
    vtmp0 = 0
    do i = 1,nvtx
       j = vpids0(i) + 1
       l = dist0(i+1) - dist0(i)
       vtmp0(j) = vtmp0(j) + l
       k = cnt0(j-1) + vtmp0(j)
       if (l.eq.1) then 
          vlist0(k) = sdat(dist0(i)+1)
       elseif (l.gt.1) then
          vlist0(k-l+1:k) = sdat(dist0(i)+1:dist0(i+1))
       else
          ! do nothing
          !print*,'?',int(rank,2)
       endif
       !vlist0(k-l+1:k) = sdat(dist0(i)+1:dist0(i+1))
    enddo
    !print*,vlist0

    lsiz1 = sum(vtmp1); allocate(vlist1(lsiz1))

    call MPI_ALLTOALLV(vlist0,vtmp0,cnt0,mpi_integer,&
                       vlist1,vtmp1,cnt1,mpi_integer,comm,ierr)

    deallocate(vtmp0,cnt0,cnt1,vlist0)

  end subroutine pnm_data_exchange
  
  subroutine pnm_find_info(nvtx,vlist,vpass,v2vsiz,v2vpid,v2vsnt,v2vrcv,comm)
    ! find information about v2v map
    ! send v2vsct to v2vrec data
    integer, intent(in)                       :: nvtx,v2vsiz,comm
    integer, intent(in)                       :: vlist(nvtx),vpass(nvtx)
    integer, intent(in)                       :: v2vsnt(v2vsiz),v2vpid(v2vsiz)
    integer, intent(out)                      :: v2vrcv(v2vsiz)                 

    integer                                   :: i,j,k,l,lct,nproc,rank
    integer                                   :: lvts1,lvts0,ierr
    integer, allocatable, dimension(:)        :: idtmp,npid0,npid1
    integer, allocatable, dimension(:)        :: vts0,vts1,vts2
    integer, allocatable, dimension(:)        :: cnt0,cnt1
    integer, allocatable, dimension(:)        :: dat0,dat1,pid

    call mpi_comm_size(comm,nproc,ierr)
    call mpi_comm_rank(comm,rank, ierr)

    allocate(vts2(v2vsiz)); vts2 = v2vpid
    call simplessort(vts2,idtmp)

    if (allocated(npid0)) deallocate(npid0)
    allocate(npid0(nproc)); npid0 = 0
    do i = 1,v2vsiz
       j = v2vpid(i)
       npid0(j+1) = npid0(j+1) + 1
    enddo

    if (sum(npid0).ne.v2vsiz) then
       print*,'Error: v2vsiz count at rank',int(rank,4)
       stop
    endif
    !print*, nproc,rank,sum(npid0),v2vsiz 
    !print*, allocated(npid1)
    !if (.not. allocated(npid1)) allocate(npid1(nproc))
    allocate(npid1(nproc)); npid1 = 0
    
    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,comm,ierr)

    lvts0 = sum(npid0); lvts1 = sum(npid1)
    allocate(cnt0(0:nproc-1),cnt1(0:nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    allocate(vts0(lvts0),vts1(lvts1)) 
    vts0 = v2vsnt(idtmp)
    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,comm,ierr)
    
    allocate(dat0(lvts0),dat1(lvts1))

    do i = 1,lvts1
       j = vts1(i)
       call findorder(j,vlist,lct)
       dat1(i) = vpass(lct)
    enddo

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_integer,&
                       dat0,npid0,cnt0,mpi_integer,comm,ierr)

    v2vrcv(idtmp) = dat0(:)
  
    deallocate(vts0,vts1,vts2,dat0,dat1,idtmp)
    deallocate(npid0,npid1,cnt0,cnt1)

  end subroutine pnm_find_info

  subroutine pnm_find_data(nvtx,vlist,vpass,v2vsiz,v2vpid,v2vsnt,v2vrcv,comm)
    ! find information about new data
    ! send v2vsct to v2vrec data
    integer, intent(in)                         :: nvtx,v2vsiz,comm
    integer, intent(in)                         :: vlist(nvtx)
    real(kind=rkind), intent(in)                :: vpass(nvtx)
    integer, intent(in)                         :: v2vsnt(v2vsiz),v2vpid(v2vsiz)
    real(kind=rkind), intent(out)               :: v2vrcv(v2vsiz)                 

    integer                                     :: i,j,k,l,lct,nproc,rank
    integer                                     :: lvts1,lvts0,ierr
    integer, allocatable, dimension(:)          :: idtmp,npid0,npid1
    integer, allocatable, dimension(:)          :: vts0,vts1,vts2
    integer, allocatable, dimension(:)          :: cnt0,cnt1
    real(kind=rkind), allocatable, dimension(:) :: dat0,dat1

    call mpi_comm_size(comm,nproc,ierr)
    call mpi_comm_rank(comm,rank, ierr)

    allocate(vts2(v2vsiz)); vts2 = v2vpid
    call simplessort(vts2,idtmp)

    allocate(npid0(nproc)); npid0 = 0
    do i = 1,v2vsiz
       j = v2vpid(i)
       npid0(j+1) = npid0(j+1) + 1
    enddo

    if (sum(npid0).ne.v2vsiz) then
       print*,'Error: v2vsiz count at rank',int(rank,4)
       stop
    endif

    allocate(npid1(nproc)); npid1 = 0
    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,comm,ierr)

    lvts0 = sum(npid0); lvts1 = sum(npid1)
    allocate(cnt0(0:nproc-1),cnt1(0:nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    allocate(vts0(lvts0),vts1(lvts1)) 
    vts0 = v2vsnt(idtmp)

    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,comm,ierr)
    
    allocate(dat0(lvts0),dat1(lvts1))

    do i = 1,lvts1
       j = vts1(i)
       call findorder(j,vlist,lct)
       dat1(i) = vpass(lct)
    enddo

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_real8,&
                       dat0,npid0,cnt0,mpi_real8,comm,ierr)

    v2vrcv(idtmp) = dat0(:)
  
    deallocate(vts0,vts1,vts2,dat0,dat1,idtmp)
    deallocate(cnt0,cnt1,npid0,npid1)


  end subroutine pnm_find_data

  subroutine pnm_pass_data(nvtx,vpass,v2vsiz,v2vpid,v2vsnt,v2vrcv,comm)
    ! find information about new data
    ! send v2vsct to v2vrec data
    integer, intent(in)                         :: nvtx,v2vsiz,comm
    real(kind=rkind), intent(in)                :: vpass(nvtx)
    integer, intent(in)                         :: v2vsnt(v2vsiz),v2vpid(v2vsiz)
    real(kind=rkind), intent(out)               :: v2vrcv(v2vsiz)                 

    integer                                     :: i,j,k,l,lct,nproc,rank
    integer                                     :: lvts1,lvts0,ierr
    integer, allocatable, dimension(:)          :: idtmp,npid0,npid1
    integer, allocatable, dimension(:)          :: vts0,vts1,vts2
    integer, allocatable, dimension(:)          :: cnt0,cnt1
    real(kind=rkind), allocatable, dimension(:) :: dat0,dat1

    call mpi_comm_size(comm,nproc,ierr)
    call mpi_comm_rank(comm,rank, ierr)

    allocate(vts2(v2vsiz)); vts2 = v2vpid
    call simplessort(vts2,idtmp)

    allocate(npid0(nproc)); npid0 = 0
    do i = 1,v2vsiz
       j = v2vpid(i)
       npid0(j+1) = npid0(j+1) + 1
    enddo

    if (sum(npid0).ne.v2vsiz) then
       print*,'Error: v2vsiz count at rank',int(rank,4)
       stop
    endif

    allocate(npid1(nproc)); npid1 = 0
    call MPI_ALLTOALL(npid0,1,mpi_integer,&
                      npid1,1,mpi_integer,comm,ierr)

    lvts0 = sum(npid0); lvts1 = sum(npid1)
    allocate(cnt0(0:nproc-1),cnt1(0:nproc-1))
    cnt0 = 0; cnt1 = 0
    do i = 1,nproc-1
       cnt0(i) = cnt0(i-1) + npid0(i)
       cnt1(i) = cnt1(i-1) + npid1(i)
    enddo

    allocate(vts0(lvts0),vts1(lvts1)) 
    vts0 = v2vsnt(idtmp)
    call MPI_ALLTOALLV(vts0,npid0,cnt0,mpi_integer,&
                       vts1,npid1,cnt1,mpi_integer,comm,ierr)
    
    allocate(dat0(lvts0),dat1(lvts1))

    do i = 1,lvts1
       j = vts1(i)
       dat1(i) = vpass(j)
    enddo

    call MPI_ALLTOALLV(dat1,npid1,cnt1,mpi_real8,&
                       dat0,npid0,cnt0,mpi_real8,comm,ierr)

    v2vrcv(idtmp) = dat0(:)
  
    deallocate(vts0,vts1,vts2,dat0,dat1,idtmp)
    deallocate(npid0,npid1,cnt0,cnt1)

  end subroutine pnm_pass_data



  subroutine pnm_find_common(l0,dat0,l1,dat1,lc,com0,cout)
     integer, intent(in)                  :: l0,l1
     integer, intent(in)                  :: dat0(l0),dat1(l1)
     integer, intent(out)                 :: lc 
     integer, allocatable, intent(out)    :: com0(:),cout(:)
     
     integer                              :: i,j,k,k0,k1
     integer, allocatable                 :: idtmp(:),dtmp(:) 

     k = l0 + l1
     allocate(dtmp(k)) 
     dtmp(1:l0) = dat0; dtmp(l0+1:l0+l1) = dat1

     call simplessort(dtmp,idtmp)
    
     lc = 0
     do i = 2,l0+l1
        if (dtmp(i-1).eq.dtmp(i)) then 
           lc = lc + 1
        endif
     enddo 
 
     allocate(com0(lc),cout(lc))
     cout = 0; com0 = 0
     k = 0 
     do i = 2,l0+l1
        if (dtmp(i-1).eq.dtmp(i)) then 
           k = k + 1
           cout(k) = dtmp(i)
        endif
     enddo 

     k = 1; i = 1
     do k = 1,lc 
        do while (cout(k) .gt. dat0(i))
           if (i.lt.l0) then 
              i = i + 1
           endif
        enddo
        if (cout(k).eq.dat0(i)) then
           com0(k) = i
        else
           print*,'Error', dat0(i),cout(k) 
        endif    

     enddo

     !print*, maxval(cout - dat0(com0))
     deallocate(idtmp,dtmp)

  end subroutine pnm_find_common

  SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
     
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
     
      CALL SYSTEM_CLOCK(COUNT=clock)
     
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)
     
      DEALLOCATE(seed)
  END SUBROUTINE
     

 
end module utility_mod
