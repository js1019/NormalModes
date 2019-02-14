!******************************************************************************!
!*  This module deals with text files and string converting.                  *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangwei.zw@gmail.com                        *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $LastChangedDate: 2006-12-10 15:20:12 +0800 (Sun, 10 Dec 2006) $
! $Revision: 62 $
! $LastChangedBy: zhangw $

!-----------------------------------------------------------------------------
module string_mod
!-----------------------------------------------------------------------------

    implicit none
    private
    public             :: string_conf
    
    interface string_conf
       module procedure get_conf_real,    &
                        get_conf_integer, &
                        get_conf_logical, &
                        get_conf_character
    end interface
    
    integer, parameter :: END_OF_LINE=-2
    integer, parameter :: END_OF_FILE=-1
    integer, parameter :: INF_IN_FILE=-3
    integer, parameter :: FD_STR_LEN =132

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! analyse string

  subroutine blank_punct(text)
  !!! remove the tail of text after "#" and replace "= " by " "
  
      character (len=*),intent(inout) :: text
      character (len=*),parameter     ::  &
                comment="#",   &
                punct="=	"
      integer                         :: i
  
      text=adjustl(text)
      i=scan(text,comment)
      if (i>0) then
         text(i:)=" "
      end if
  
      do i=1,len_trim(text)
         if (index(punct,text(i:i))>0) then
            text(i:i)=" "
         end if
      end do
  
  end subroutine blank_punct
  
  subroutine compress_bb(text)
  !!! replace "  " by " "
      character (len=*), intent(inout) :: text
      integer                          :: i

      do
          i=index(trim(text(1:len_trim(text))),"  ")
          !i=index(trim(text),"  ")
          if (i==0) exit
          text(i:)=text(i+1:)
      end do

  end subroutine compress_bb
  
  function count_word(text) result(n)
      character (len=*),intent(inout) :: text
      character (len=FD_STR_LEN)      :: str
      integer                         ::  n,end_of_word
      
      call blank_punct(text)
      call compress_bb(text)
      str=text
      n  = 0
      if (len_trim(str)==len(str)) then
         n=1
      else
         do
            if (len_trim(str)==0) exit
            end_of_word=index(str," ")-1
            n=n+1
            str=str(end_of_word+2:)
         end do
      end if
  end function count_word
  
  function get_part_str(text,col,flag) result(str)
      character (len=*),intent(in) :: text
      integer col
      character (len=FD_STR_LEN)   :: str,tmpstr
      logical,optional             :: flag
      logical stick
      integer i,n,end_of_word
  
      stick=.false.
      if (present(flag)) stick=flag
      tmpstr = text
      n      = count_word(tmpstr)
      if (n<col .and. stick) then
         print *, "there is only",n," column, but you want to read",col
         print *, "the str is ",text
         stop
      elseif (n<col) then
         str=" "
      else
         do i=1,col
            end_of_word=index(tmpstr," ")-1
            str=tmpstr(1:end_of_word)
            tmpstr=tmpstr(end_of_word+2:)
         end do
      end if
  end function get_part_str
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! convert nth part of str to datatype !!!
  subroutine convert_str_integer(text,n,gdata)
      character (len=*),intent(in) :: text
      integer                      :: gdata
      integer                      ::  n
      character (len=FD_STR_LEN)   :: str
      str=get_part_str(text,n,.true.)
      read(str,*) gdata
  end subroutine convert_str_integer

  subroutine convert_str_real(text,n,gdata)
      character (len=*),intent(in) :: text
      real                         ::gdata
      integer                      :: n
      character (len=FD_STR_LEN)   :: str
      str=get_part_str(text,n,.true.)
      read(str,*) gdata
  end subroutine convert_str_real

  subroutine convert_str_logical(text,n,gdata)
      character (len=*),intent(in) :: text
      logical                      :: gdata
      integer                      :: n
      character (len=FD_STR_LEN)   :: str
      str=get_part_str(text,n,.true.)
      read(str,*) gdata
  end subroutine convert_str_logical

  subroutine convert_str_character(text,n,gdata)
      character (len=*),intent(in) :: text
      character (len=*)            :: gdata
      integer                      :: n
      gdata=get_part_str(text,n,.true.)
  end subroutine convert_str_character
  
  !!! read record from conf file !!!
  subroutine get_conf_line(fid,gstr)
      character (len=*) :: gstr
      character (len=8) :: fmt_str
      integer           :: err,fid
      write(fmt_str,"(a2,i5,a1)") "(a",FD_STR_LEN,")"
      do
          read(fid,fmt_str,iostat=err) gstr
          if (err<0) then
             print *, "end of file, can't get_conf_line"
             stop
          end if
          call blank_punct(gstr)
          call compress_bb(gstr)
          if ( len_trim(gstr)/=0 ) exit
      end do
  end subroutine get_conf_line
  
  subroutine strip_conf_file(fid,kcol,kword,gcol,gstr)
      integer ::  &
          fid,    &
          kcol,   &   !key column
          gcol        !get column
      character (len=*) ::  & 
          kword,   &
          gstr
      character (len=Fd_STR_LEN) ::  & 
          text,  &
          prestr
      character (len=8) :: fmt_str
      integer err
      write(fmt_str,"(a2,i5,a1)") "(a",FD_STR_LEN,")"
      rewind(fid)
      do
          read(fid,fmt_str,iostat=err) text
          if (err<0) then
             print *, "there is no ",kword," on column",kcol
             stop
          end if
          call blank_punct(text)
          call compress_bb(text)
          !prestr=get_part_str(text,kcol)
          !if (trim(prestr)==trim(kword)) then
          if (get_part_str(text,kcol)==trim(kword)) then
             gstr=get_part_str(text,gcol,.true.)
             exit
          end if
      end do
  end subroutine strip_conf_file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! shell of strip_conf_file, to convert the gstr to datatype !!!
  subroutine get_conf_real(fid,kcol,kword,gcol,gdata)
      integer ::  &
          fid,    &
          kcol,   &   !key column
          gcol        !get column
      real gdata
      character (len=*) ::  & 
          kword
      character (len=Fd_STR_LEN) ::  & 
          gstr
      call strip_conf_file(fid,kcol,kword,gcol,gstr)
      read(gstr,*) gdata
  end subroutine get_conf_real

  subroutine get_conf_integer(fid,kcol,kword,gcol,gdata)
      integer ::  &
          fid,    &
          kcol,   &   !key column
          gcol        !get column
      integer gdata
      character (len=*) ::  & 
          kword
      character (len=Fd_STR_LEN) ::  & 
          gstr
      call strip_conf_file(fid,kcol,kword,gcol,gstr)
      read(gstr,*) gdata
  end subroutine get_conf_integer

  subroutine get_conf_logical(fid,kcol,kword,gcol,gdata)
      integer ::  &
          fid,    &
          kcol,   &   !key column
          gcol        !get column
      logical gdata
      character (len=*) ::  & 
          kword
      character (len=Fd_STR_LEN) ::  & 
          gstr
      call strip_conf_file(fid,kcol,kword,gcol,gstr)
      read(gstr,*) gdata
  end subroutine get_conf_logical

  subroutine get_conf_character(fid,kcol,kword,gcol,gstr)
      integer ::  &
          fid,    &
          kcol,   &   !key column
          gcol        !get column
      character (len=*) ::  & 
          kword,  &
          gstr
      call strip_conf_file(fid,kcol,kword,gcol,gstr)
  end subroutine get_conf_character
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module string_mod

