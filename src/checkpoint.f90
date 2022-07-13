module checkpoint_m
 implicit none
 private
 public :: checkpoint

 logical,save::save_grain_p=.true.
 real(4),save::checkpoint_time=15.*60.  !time (seconds) between checkpoint saves
 character(200),save::chkfiles(4)=(/'mrm.chk      ', 'mrm_grain.chk', 'mrm_0.chk    ', 'mrm_1.chk    '/)
 real(4),save::start_time=-1.0
 integer,save::idd=1

contains

 function checkpoint(read_p,clean_p,time,iline,isave_grain_p,icheckpoint_time,present_p) result (exit_p)
  use io_m, only:output, add_file, delete_file
  use mag_m, only:file_mag
  use grain_m, only:file_media
  use random_m, only:save_random, restore_random
  logical,intent(in),optional::read_p,isave_grain_p,clean_p,present_p
  real(8),intent(inout),optional::time
  real(4),intent(in),optional::icheckpoint_time
  integer,intent(inout),optional::iline
  logical::exit_p,ok_p
  integer::i,id
  character(200)::os
  real(4)::current_time

  if (present(present_p)) then
     if (present_p) then
        exit_p=.true.
        do id=1,2; os=chkfiles(id); inquire(file=trim(os),exist=ok_p); exit_p=(exit_p.and.ok_p); enddo !is mrm.chk AND mrm_grain.chk present
        if (.not.exit_p) return
        exit_p=.false.
        do id=3,4; os=chkfiles(id); inquire(file=trim(os),exist=ok_p); exit_p=(exit_p.or.ok_p); enddo !is mrm_0.chk OR mrm_1.chk present
        return
     endif
   endif

  if (present(isave_grain_p)) save_grain_p=isave_grain_p
  if (present(icheckpoint_time)) checkpoint_time=icheckpoint_time
  if (start_time.lt.0.0)  call cpu_time(start_time)

  exit_p=.true.
  if (present(isave_grain_p).or.present(icheckpoint_time)) return

  exit_p=.true.; os='mrm.chk'; i=add_file(trim(os)); call delete_file(i)

  if (present(clean_p)) then; if (clean_p) then
    call output(' CHK: attempting to clean up')
    do id = 1, size(chkfiles)
      os=chkfiles(id); inquire(file=trim(os),exist=ok_p)
      if (ok_p) then; open(i,file=trim(os),form='unformatted',status='old'); close(i,status='delete'); endif
    enddo
    call output(' CHK: clean up')
    return
  endif; endif

  if (read_p) then  !read checkpoint
    call cpu_time(start_time)
    do id=1,2; os=chkfiles(id); inquire(file=trim(os),exist=ok_p); if (.not.ok_p) return; enddo !give up if mrm.chk or mrm_grain.chk not present

    call output(' CHK: attempting to read checkpoint...')

!   os=chkfiles(2); open(i,file=trim(os),form='unformatted',status='old')
    os=chkfiles(2); open(i,file=trim(os),form='unformatted',status='old',action='read',position='rewind')
    exit_p= file_media(.true.,i)  !read in media 

    if (exit_p) then
      save_grain_p=.false.

      os=chkfiles(1)
!     open(i,file=trim(os),form='unformatted',status='old')
!     read(i) time,iline,idd; id=1-idd
      open(i,file=trim(os),form='formatted',status='old',action='read',position='rewind')   !save time
!     open(i,file=trim(os),form='formatted',status='old')   !save time
      read(i,*) time,iline,idd; id=1-idd
      close(i)

      ok_p=.true.
      do while (ok_p); ok_p=.false.
        os=chkfiles(3+idd); inquire(file=trim(os),exist=exit_p)
        if (exit_p) then
          open(i,file=trim(os),form='unformatted',status='old',action='read',position='rewind')
!         open(i,file=trim(os),form='unformatted',status='old')
          read(i) time,iline
          exit_p=file_mag(i,.true.)
          if (exit_p) exit_p=restore_random(i)
          close(i)
        endif
        if (.not.exit_p) then; ok_p=(idd.ne.id); idd=id; endif
      enddo
    endif
    if (exit_p) then 
      write(os,"(a,a)") ' CHK: read checkpoint!',' ('//trim(os)//')'; call output(trim(os))
!     call output(' CHK: running through cli...')
    else
      call output(' CHK: failure reading checkpoint!')
      call output(' CHK: aborting...')
    endif
  else
    call cpu_time(current_time)
    ok_p=.false.
    if (abs(current_time-start_time).gt.checkpoint_time.or.ok_p) then
      start_time=current_time

      if (save_grain_p) then
        call output(' CHK: attempting to save microstructure...')
        os=chkfiles(2); open(i,file=trim(os),form='unformatted',status='unknown')
        exit_p= file_media(.false.,i)
        if (.not.exit_p) then
          close(i,status='delete'); return
          call output(' CHK: failure to save microstructure!')
        else
          close(i)
          call output(' CHK: saved microstructure!')
        endif
        save_grain_p=.false.
      endif

      call output(' CHK: attempting to save checkpoint...')
      idd=1-idd; os=chkfiles(3+idd)
      open(i,file=trim(os),form='unformatted',status='unknown')
      write(i) time,iline
      exit_p=file_mag(i,.false.)
      if (exit_p) then
        exit_p=save_random(i)
        if (exit_p) then
          close(i)
          os=chkfiles(1)
!         open(i,file=trim(os),form='unformatted',status='unknown')
!         write(i) time,iline,idd
          open(i,file=trim(os),form='formatted',status='unknown')
          write(i,*) time,iline,idd
          close(i)
          os=chkfiles(3+idd)
          write(os,"(a,a)") ' CHK: saved checkpoint!',' ('//trim(os)//')'; call output(trim(os))
        else
          close(i,status='delete')
          write(os,"(a,a)") ' CHK: failure to save checkpoint (random)!',' ('//trim(os)//')'; call output(trim(os))
      endif
      else
        close(i,status='delete')
        write(os,"(a,a)") ' CHK: failure to save checkpoint!',' ('//trim(os)//')'; call output(trim(os))
      endif
    endif

  endif

 end function
end module checkpoint_m
