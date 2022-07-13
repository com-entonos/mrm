
module io_m
 implicit none
 private
 
 public::output,set_output,flush_output,find_file,add_file,delete_file,read_comments,show_file,fatalio
 interface delete_file
   module procedure delete_filef,delete_filei
 end interface
 interface find_file
   module procedure find_filef,find_filei
 end interface

 type::file_s
   integer::ionum
   character(200)::filename
   type(file_s),pointer::next=>null()
 end type file_s

 integer,save::ionum=-1
 logical,save::file_p=.false.,scrn_p=.false.,fatalio_p=.false.
 type(file_s),pointer,save::file_list => null()
contains

 function fatalio() result (ok_p)
  logical::ok_p
  ok_p=fatalio_p
 end function

 subroutine read_comments(ionum)
 !silly routine to skip comments while reading an ascii file

  integer,intent(in)::ionum
  integer::io,i
  character(1)::a1
  character(3)::com='#!"' 
  
  io=0; i=1
  do while (io.eq.0.and.i.eq.1)
    read(ionum,"(a)",iostat=io) a1    !read first character on the line
    i = scan(adjustl(a1),com)         !returns the index of one of the comment characters if present (i.e. 1), otherwise 0
  enddo
  if (io.eq.0) backspace(ionum)
 end subroutine read_comments


 subroutine output(str,force_p)
  character(*)::str
  logical,intent(in),optional::force_p
  if (scrn_p) then
    print "(a)",trim(str)
  elseif (present(force_p)) then
    if (force_p) print "(a)",trim(str)
  endif
  if (ionum.gt.0) write(ionum,"(a)") trim(str)
 end subroutine

 subroutine flush_output()
   if (ionum.gt.0) then
     if (file_p) then
       close(ionum)
     else
       close(ionum,status='delete')
     endif
     call delete_file(ionum); ionum=-1
   endif
 end subroutine

 subroutine set_output(screen_file,save_screen_p,append_p,quiet_p,ifatalio_p)
  character(*),intent(in)::screen_file
  logical,intent(in)::save_screen_p,append_p,quiet_p,ifatalio_p

  fatalio_p=ifatalio_p
  scrn_p=(.not.quiet_p)
  if (scrn_p.and.len_trim(adjustl(screen_file)).gt.0) then
    ionum=add_file(trim(adjustl(screen_file)))
    if (ionum.gt.0) then
      inquire(file=trim(adjustl(screen_file)),exist=file_p)
      if (append_p.and.file_p) then
        open(ionum,file=trim(adjustl(screen_file)),status='old',position='append')
      else
        open(ionum,file=trim(adjustl(screen_file)),status='unknown')
      endif
      file_p=save_screen_p
    endif
  endif
 end subroutine

 function find_filei(ionum) result(filename)
   integer,intent(in)::ionum
   character(200)::filename
   type(file_s),pointer::f1

   filename=''
   f1=>file_list
   do while (associated(f1))
     if (ionum.eq.f1%ionum) then
       filename=trim(f1%filename)
       exit
     endif
     f1=>f1%next
   enddo
 end function

 function find_filef(filename) result(ionum)
   character(*)::filename
   integer::ionum
   type(file_s),pointer::f1

   ionum=-1
   f1=>file_list
   do while (associated(f1))
     if (trim(filename).eq.trim(f1%filename)) then
       ionum=f1%ionum
       exit
     endif
     f1=>f1%next
   enddo
 end function

 function add_file(filename) result(ionum)
   character(*)::filename
   integer::ionum
   type(file_s),pointer::f1,f2,f3
   logical::ionum_used(10:99)

   ionum=-1; ionum_used=.false.
   f1=>file_list; nullify(f2); nullify(f3)
   do while (associated(f1))
     if (trim(filename).eq.trim(f1%filename)) f2=>f1
     ionum_used(f1%ionum)=.true.
     f3=>f1;f1=>f1%next
   enddo
   if (associated(f2)) return
   do ionum=lbound(ionum_used,DIM=1),ubound(ionum_used,DIM=1)
     if (.not.ionum_used(ionum)) exit
   enddo
   if (ionum.gt.ubound(ionum_used,DIM=1)) then
     ionum=-1
   else
     allocate(f2); f2%filename=trim(filename); f2%ionum=ionum; nullify(f2%next)
     if (associated(file_list)) then
       f3%next=>f2
     else
       file_list=>f2
     endif
   endif
 end function

 subroutine delete_filef(filename)
   character(*),intent(in)::filename
   type(file_s),pointer::f1,f2=>null()

   f1=>file_list
   do while (associated(f1))
     if (trim(filename).eq.trim(f1%filename)) then
       if (associated(f1,file_list)) then
         file_list=>f1%next; deallocate(f1); exit
       else
         f2%next=>f1%next; deallocate(f1); exit
       endif
     endif
     f2=>f1; f1=>f1%next
   enddo
 end subroutine

 subroutine delete_filei(ionum)
   integer,intent(in)::ionum
   type(file_s),pointer::f1,f2

   f1=>file_list
   do while (associated(f1))
     if (f1%ionum.eq.ionum) then
       if (associated(f1,file_list)) then
         file_list=>f1%next; deallocate(f1); exit
       else
         f2%next=>f1%next; deallocate(f1); exit
       endif
     endif
     f2=>f1; f1=>f1%next
   enddo
 end subroutine

 subroutine show_file()
   type(file_s),pointer::f1
   integer::i
   character(300)::os

   call output('list of files:')
   f1=>file_list; i=0
   do while (associated(f1))
     i=i+1
     write(os,"(' ',i2,':',i3,' ',a)") i,f1%ionum,trim(f1%filename); call output(trim(os))
     f1=>f1%next
   enddo
 end subroutine

end module io_m
