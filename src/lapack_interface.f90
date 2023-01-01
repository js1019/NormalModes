! file: interface for lapack 
module lapack_interface 

  use mpi
  use omp_lib

 implicit none 
 contains
 
   ! solve a dense symmetric eigenvalue problem
   subroutine intface_dsyev(ls,a2,wev0)
      integer, intent(in)                          :: ls
      double precision, intent(out)                :: wev0(ls)
      double precision, intent(inout)              :: a2(ls,ls)
      ! call dsyev 
      character*1                                  :: JOBZ, UPLO
      integer                                      :: LDA, a2siz, lwork, info
      double precision, dimension(:), allocatable  :: work      

      JOBZ = 'V'; UPLO = 'L'; a2siz = ls; LDA = ls;
      lwork = 3*a2siz-1; allocate(work(lwork)) 
      !if (mymatvec%rank.eq.0) then 
      call dsyev(JOBZ, UPLO, a2siz, a2, LDA, wev0, work, lwork, info)


   end subroutine intface_dsyev

   ! solve a symmetric tridiagonal eigenvalue problam
   subroutine intface_dstev(ls,DEV,offd,Zo)
      integer, intent(in)                          :: ls
      double precision, intent(inout)              :: DEV(ls) 
      double precision, intent(in)                 :: offd(ls-1) 
      double precision, intent(out)                :: Zo(ls,ls) 
 
      character*1                                  :: JOBZ
      integer                                      :: LDZ, N, lwork, info
      double precision, dimension(:), allocatable  :: work      

      JOBZ = 'V'; N = ls; LDZ = ls
      lwork = 2*N-2; allocate(work(lwork))

      call dstev(JOBZ, N, DEV, offd, Zo, LDZ, work, info)

   end subroutine intface_dstev

   ! solve a dense hermitian eigenvalue problem
   subroutine intface_zheev(ls,MatD,wev)
      integer, intent(in)                          :: ls
      double precision, intent(out)                :: wev(ls)
      complex*16, intent(inout)                    :: MatD(ls,ls) 
       
      character*1                                  :: JOBZ, UPLO
      integer                                      :: LDA, N, lwork, lrw, info
      double precision, dimension(:), allocatable  :: rwork   
      complex*16, dimension(:), allocatable        :: work   

      JOBZ = 'V'; UPLO = 'L'; N = ls; LDA = ls;
      lwork = 2*N-1; allocate(work(lwork))
      lrw   = 3*N-2; allocate(rwork(lrw)) 

      call zheev(JOBZ, UPLO, N, MatD, LDA, wev, work, lwork, rwork, info)

      !print*, info
   end subroutine intface_zheev





end module lapack_interface


