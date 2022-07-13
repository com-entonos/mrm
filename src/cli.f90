module cli_m
 implicit none
 private
 public :: cli, checkpoint

 interface checkpoint
   module procedure my_checkpoint
 end interface

 integer,parameter::MAXCHAR=500,MAXTKN=30

 integer,save::line_read=-1,line_target=-2
 real(8),save::time=0.d0,time_target=-1.d300
 integer,save::iol(100)=0,iol_num=0,ioc(100)=0,ioc_num=0

contains

 subroutine close_all_files()
  use io_m, only:output, delete_file, find_file
  integer::i
  character(200)::os

  do i=1,iol_num
    write(os,"(a)") '> CLOSE '//trim(find_file(iol(i)))//' KEEP'; call output(os)
    close(iol(i)); call delete_file(iol(i))
  enddo
  do i=1,ioc_num
    if (ioc(i).gt.0) then
      write(os,"(a)") '> CLOSE '//trim(find_file(ioc(i)))//' KEEP'; call output(os)
      close(ioc(i)); call delete_file(ioc(i))
    endif
  enddo
 end subroutine

 function my_checkpoint(read_p,quit_p,clean_p) result (exit_p)
  use io_m, only:output, add_file, delete_file
  use checkpoint_m, only: checkpoint
  logical,intent(in)::read_p
  logical,intent(in),optional::quit_p,clean_p
  logical::exit_p

  exit_p=.true.

  if (present(quit_p)) then  !caught signal, just close files and get out!
    if (quit_p) then
      call output(' CHK: attempting to quit...'); call close_all_files(); call output(' CHK: quit!')
      exit_p=.false.
      return
    endif
  endif

  if (read_p) then
    exit_p=checkpoint(read_p=.true.,clean_p=clean_p,time=time_target,iline=line_target)
  else
    exit_p=checkpoint(read_p=.false.,clean_p=clean_p,time=time,iline=line_read)
  endif


 end function

 subroutine cli(cli_file)
  use data_m, only:get_size,set_data,update_do_grain,show_window
  use grain_m, only:file_media,generate_hk,generate_area,generate_media,show_media
  use decay_m, only:decay_mag,set_decay_mag
  use mag_m, only:file_mag,generate_mag
  use demag_m, only:init_demag,demag_window
  use applied_m, only:file_temp,file_fld,uniform_fld,suprgaus_temp,set_temp,set_fld,show_fld_set,show_fld,show_temp_set,show_temp,ref_fld, &
         ref_temp,dup_fld,dup_temp
  use temp_scaling_m, only:file_temp_scaling,show_temp_scaling
  use field_m, only: plot_hamr, plot_field,set_source,show_grain,find_grain
  use off_field_m, only: plot_off_field
  use llg_m, only: set_llg,llg,show_llg,plot_llg
  use plot_m, only: set_plot, plot_mh, plot_temp_scaling, plot_wf_scaling, plot_mz_scaling
  use io_m, only:output, find_file, add_file, delete_file, fatalio
  use random_m, only:show_random, save_random, restore_random
  use checkpoint_m, only: checkpoint

  character(*),intent(in)::cli_file
  logical::ok_p,file_p,key_p,x_p,y_p,z_p,v_p,fatalio_p
  integer::i,j,k,n_lay,cliionum,num_token,layer_n
  real(8)::h(MAXTKN)
  character(MAXCHAR)::cl,token(MAXTKN),layer_s,os,s1,s2,s3
  real(4)::start_time,stop_time

  iol_num=0
  call get_size(i,j,n_lay,h(1),h(2))   !get number of layers
  fatalio_p=fatalio()

  inquire(file=trim(adjustl(cli_file)),exist=file_p)
  key_p=(trim(adjustl(cli_file)).eq.'-'); ioc_num=1; ioc(ioc_num)=0
  if (.not.key_p.and..not.file_p) return
  if (file_p.and..not.key_p) then
    cliionum=add_file(trim(adjustl(cli_file))); ioc_num=1; ioc(ioc_num)=cliionum
    open(cliionum,file=trim(adjustl(cli_file)),status='old',action='read')
  endif

  if (key_p) then
    call output('reading cli from keyboard/pipe')
    print "(' command< ',$)"
    read(*,"(a)",iostat=i) cl
  else
    call output('reading cli from file: '//trim(adjustl(cli_file)))
    read(cliionum,"(a)",iostat=i) cl
  endif
  call cpu_time(start_time)

  if (.not.my_checkpoint(.true.)) return
  if (line_read.lt.line_target) call output(' CHK: running through cli...')

  do while (i.eq.0); line_read=line_read+1
    call output('< '//trim(adjustl(cl)))
    call get_token(cl,token,num_token,layer_n,layer_s)

!     print "(i3,x,i3,:,50(x,a))",num_token,layer_n,(trim(token(i)),i=1,num_token),trim(layer_s)
    if (num_token.gt.0) then
      h=0.d0; do i=1,num_token; h(i)=getnum(token(i));enddo
     ! BYE or QUIT or EXIT tells us to quit right now
      if (token(1)(1:3).eq.'BYE'.or.token(1)(1:4).eq.'QUIT'.or.token(1)(1:4).eq.'EXIT') then
        exit
     ! SAVE <filename> !save random number generators
      elseif (token(1)(1:4).eq.'SAVE'.and.num_token.gt.1) then
        write(os,"(a)")'> SAVE '//trim(token(2)); call output(os)
        if (line_read.gt.line_target) then
          i=add_file(trim(token(2)))
          if (i.gt.0) then
            open(i,file=trim(token(2)),form='unformatted',status='unknown')
!           if (.not.save_random(i))  call output('!---- failed to save ---!')
            if (.not.save_random(i,.false.))  call output('!---- failed to save ---!')
            close(i)
            call delete_file(i)
          endif
        endif
     ! RESTORE <filename> !restore random number generators
      elseif (token(1)(1:7).eq.'RESTORE'.and.num_token.gt.1) then
        write(os,"(a)")'> RESTORE '//trim(token(2)); call output(os)
        inquire(file=trim(token(2)),exist=ok_p)
        if (ok_p) then
          i=add_file(trim(token(2)))
          if (i.gt.0) then
            open(i,file=trim(token(2)),form='unformatted',action='read',position='rewind',status='old')
!           if (.not.restore_random(i))  call output('!---- no such file exists ---!')
            if (.not.restore_random(i,.false.))  call output('!---- no such file exists ---!')
            close(i)
            call delete_file(i)
          endif
        endif
     ! OPENR[ead] <filename> [FORMATTED|unformatted]
     ! OPENW[rite] <filename> [FORMATTED|unformatted] [REWIND|append]
      elseif (token(1)(1:5).eq.'OPENR'.or.token(1)(1:5).eq.'OPENW') then
       if (num_token.gt.1) then
        s1='FORMATTED'; s2='REWIND'
        do i=3,num_token
          if (token(i)(1:6).eq.'UNFORM') s1='UNFORMATTED'
          if (token(i)(1:3).eq.'APP') s2='APPEND'
        enddo
        inquire(file=trim(token(2)),exist=ok_p)
        if (token(1)(1:5).eq.'OPENR') then
          if (ok_p) then
            i=add_file(trim(token(2)))
            if (i.gt.0) then
              iol_num=mod(iol_num,size(iol))+1; iol(iol_num)=i
              write(os,"(a)")'> OPENREAD '//trim(token(2))//' '//trim(s1); call output(os)
!             open(i,file=trim(token(2)),form=trim(s1),status='old')
!             open(i,file=trim(token(2)),form=trim(s1),action='read',shared,position='rewind',readonly,status='old')
              open(i,file=trim(token(2)),form=trim(s1),action='read',position='rewind',status='old')
            else
              write(os,"('!--- file ',a,' already open ---')")trim(token(2)); call output(os)
            endif
          else
            write(os,"('!--- no such file ',a,' exists ---')") trim(token(2)); call output(os)
          endif
        else
          i=add_file(trim(token(2)))
          if (i.gt.0) then
            iol_num=mod(iol_num,size(iol))+1; iol(iol_num)=i
            if (ok_p.and.(s2(1:6).eq.'APPEND'.or.line_read.lt.line_target)) then
              write(os,"(a)") '> OPENWRITE '//trim(token(2))//' '//trim(s1)//' '//trim(s2); call output(os)
!             open(i,file=trim(token(2)),position=trim(s2),form=trim(s1),status='old')
              open(i,file=trim(token(2)),position='APPEND',form=trim(s1),status='old')
            else
              write(os,"(a)") '> OPENWRITE '//trim(token(2))//' '//trim(s1); call output(os)
              open(i,file=trim(token(2)),form=trim(s1),status='unknown')
            endif
          else
            write(os,"('!--- file ',a,' already open ---')")trim(token(2)); call output(os)
          endif
        endif
       else
         call output('!--- OPENWRITE/OPENREAD requires at last a file name: filename [APPEND] [UNFORMATTED] ---')
       endif
     ! CLOSE <filename> [KEEP|delete]
      elseif (token(1)(1:5).eq.'CLOSE') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          s1='KEEP'
          do j=3,num_token; if (token(j)(1:3).eq.'DEL') s1='DELETE'; enddo
          write(os,"(a)")'> CLOSE '//trim(token(2))//' '//trim(s1); call output(os)
          close(i,status=trim(s1),iostat=j)
          call delete_file(i)
          do j=1,size(iol); if (iol(j).eq.i) exit; enddo
          do i = j+1,size(iol); iol(i-1)=iol(i); enddo
          if (j+1.le.size(iol)) iol_num=iol_num-1
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
        endif
     ! READ_C[omment] <filename>
      elseif (token(1)(1:6).eq.'READ_C') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          write(os,"(a)") '> READ_COMMENT '//trim(token(2)); call output(os)
          inquire(i,form=s1); s1=adjustl(s1)
          if (s1(1:1).eq.'F') then; read(i,*); else; read(i); endif
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
        endif
     ! WRITE_C[omment] <filename> [token1] [token2] ...
      elseif (token(1)(1:7).eq.'WRITE_C') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
              if (line_read.gt.line_target) then
          write(os,"(a,:,50(x,a))") '> WRITE_COMMENT '//trim(token(2)),(trim(token(j)),j=3,num_token); call output(os)
          inquire(i,form=s1); s1=adjustl(s1)
          if (s1(1:1).eq.'F') then
            write(i,"(:,50(a,x))") (trim(token(j)),j=3,num_token)
          else 
            write(i) (trim(token(j)),j=3,num_token)
          endif
              endif
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
        endif
     ! WRITE_G[rain] <filename>
      elseif (token(1)(1:7).eq.'WRITE_G') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          inquire(i,form=s1); s1=adjustl(s1)
          if (s1(1:1).eq.'F') then
            write(os,"(a)") '!--- file '//trim(token(2))//' is formatted, reopening as unformatted ---'; call output(os)
              inquire(file=trim(token(2)),exist=ok_p)
              if (line_read.gt.line_target.or..not.ok_p) then
            close(i); open(i,file=trim(token(2)),form='unformatted',status='unknown')
              else
            close(i); open(i,file=trim(token(2)),form='unformatted',status='old',position='append')
              endif
          endif
              if (line_read.gt.line_target) then
          write(os,"(a)")'> WRITE_GRAIN '//trim(token(2))//' '//trim(layer_s);call output(os)
          ok_p= file_media(.false.,i,layer=layer_n)
          if (.not.ok_p) call output('!--- WRITE_GRAIN failed ---')
              endif
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
        endif
     ! READ_G[rain] <filename>
      elseif (token(1)(1:6).eq.'READ_G') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          inquire(i,form=s1); s1=adjustl(s1)
          if (s1(1:1).eq.'F') then
            write(os,"(a)") '!--- file '//trim(token(2))//' is formatted, reopening as unformatted ---'; call output(os)
            close(i,status='keep'); open(i,file=trim(token(2)),form='unformatted',status='old',action='read',position='rewind')
          endif
          ok_p=.false.
          do j=num_token,3,-1
            ok_p=(ok_p.or.token(j)(1:3).eq.'OLD')
          enddo
          if (ok_p) then
            write(os,"(a)")'> READ_GRAIN '//trim(token(2))//' OLD '//trim(layer_s);call output(os)
          else
            write(os,"(a)")'> READ_GRAIN '//trim(token(2))//' '//trim(layer_s);call output(os)
          endif
            ok_p= file_media(.true.,i,layer=layer_n,old_p=ok_p)
          if (.not.ok_p) call output('!--- READ_GRAIN failed ---')
          if (.not.ok_p.and.fatalio_p) stop
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2));call output(os)
        endif
     ! WRITE_M[ag] <filename> [transition#]
      elseif (token(1)(1:7).eq.'WRITE_M') then
        i=find_file(trim(token(2)))
              if (line_read.gt.line_target) then
        if (i.gt.0) then
          ok_p=.false.; j=-1
          do k=num_token,3,-1
            if (token(k)(1:3).ne.'OLD'.and.token(k)(1:3).ne.'NEW') j=nint(h(k))
            ok_p=(ok_p.or.token(k)(1:3).eq.'OLD')
          enddo
          if (j.lt.1) then; j=1; os=''; else; write(os,"(x,i4)") j;endif
          if (ok_p) then; os=trim(os)//' OLD'; else; os=trim(os)//' NEW'; endif
          write(os,"(a)")'> WRITE_MAG '//trim(token(2))//trim(os)//' '//trim(layer_s);call output(os)
          ok_p= file_mag(i,.false.,old=ok_p,old_tran=j,layer=layer_n)
          if (.not.ok_p) call output('!--- WRITE_MAG failed ---')
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2));call output(os)
        endif
              endif
     ! READ_M[ag] <filename> 
      elseif (token(1)(1:6).eq.'READ_M') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          ok_p=.false.
          do j=num_token,3,-1
            ok_p=(ok_p.or.token(j)(1:3).eq.'OLD')
          enddo
          if (ok_p) then
            write(os,"(a)")'> READ_MAG '//trim(token(2))//' OLD '//trim(layer_s);call output(os)
          else
            write(os,"(a)")'> READ_MAG '//trim(token(2))//' '//trim(layer_s);call output(os)
          endif
            ok_p= file_mag(i,.true.,layer=layer_n,old=ok_p)
          if (.not.ok_p) call output('!--- READ_MAG failed ---')
          if (.not.ok_p.and.fatalio_p) stop
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2));call output(os)
        endif
     ! READ_S[cale] <filename> [MS|alpha|Hk|Hx]
      elseif (token(1)(1:6).eq.'READ_S') then
        i=find_file(trim(token(2)))
        if (i.gt.0) then
          inquire(i,form=s1); s1=adjustl(s1)
          if (s1(1:1).ne.'F') then
            write(os,"(a)") '!--- file '//trim(token(2))//' is unformatted, reopening as formatted ---';call output(os)
            close(i,status='keep'); open(i,file=trim(token(2)),status='old',action='read',position='rewind')
          endif
          j=0
          do k=num_token,3,-1
            select case (token(k)(1:1))
              case ("M") !Ms
                j=1; s1=token(k)
              case ("A","D") !alpha
                j=2; s1=token(k)
              case ("H","K") !Hk
                j=3; s1=token(k)
              case ("E",'X') !Hx
                j=4; s1=token(k)
            end select
          enddo
          if (j.ne.0) then
            write(os,"(a)")'> READ_SCALE '//trim(token(2))//' '//trim(s1)//' '//trim(layer_s);call output(os)
            ok_p= file_temp_scaling(i,j,ilayer=layer_n)
            if (.not.ok_p) call output('!--- READ_SCALE failed ---')
            if (.not.ok_p.and.fatalio_p) stop
          else
            call output('!--- READ_SCALE requires two tokens: filename and Ms|Alpha|Hk|exchange ---')
          endif
        else
          write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
        endif
     ! READ_T[emp] <filename> <index>
      elseif (token(1)(1:6).eq.'READ_T') then
        if (num_token.gt.2) then
          i=find_file(trim(token(2)))
          if (i.gt.0) then
            inquire(i,form=s1); s1=adjustl(s1)
            if (s1(1:1).ne.'F') then
              write(os,"(a)") '!--- file '//trim(token(2))//' is unformatted, reopening as formatted ---'; call output(os)
              close(i,status='keep'); open(i,file=trim(token(2)),status='old',action='read',position='rewind')
            endif
            do j=num_token,5,-1; if (h(j).ne.0.d0.and.h(4).eq.0.d0) h(4)=h(j); enddo; if (h(4).eq.0.d0.or.num_token.lt.4) h(4)=1.d0; j=nint(h(3))
            write(os,"(a,i4,x,1p,e12.5)")'> READ_TEMP '//trim(token(2)),j,h(4); call output(os)
            ok_p=file_temp(j,i,h(4))
            if (.not.ok_p) call output('!--- READ_TEMP failed ---')
            if (.not.ok_p.and.fatalio_p) stop
          else
            write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
          endif
        else
          call output('!--- READ_TEMP requires two tokens: filename and index [scale] ---')
        endif
     ! READ_F[ield] <filename> <index>
      elseif (token(1)(1:6).eq.'READ_F') then
        if (num_token.gt.2) then
          i=find_file(trim(token(2)))
          if (i.gt.0) then
            inquire(i,form=s1); s1=adjustl(s1)
            if (s1(1:1).ne.'F') then
              write(os,"(a)") '!--- file '//trim(token(2))//' is unformatted, reopening as formatted ---'; call output(os)
              close(i,status='keep'); open(i,file=trim(token(2)),status='old',action='read',position='rewind')
            endif
            do j=num_token,5,-1; if (h(j).ne.0.d0.and.h(4).eq.0.d0) h(4)=h(j); enddo; if (h(4).eq.0.d0.or.num_token.lt.4) h(4)=1.d0; j=nint(h(3))
            write(os,"(a,i4,x,1p,e12.5)")'> READ_FIELD '//trim(token(2)),j,h(4); call output(os)
            ok_p=file_fld(j,i,h(4))
            if (.not.ok_p) call output('!--- READ_FIELD failed ---')
            if (.not.ok_p.and.fatalio_p) stop
          else
            write(os,"('!--- file ',a,' is NOT open ---')")trim(token(2)); call output(os)
          endif
        else
          call output('!--- READ_FIELD requires two tokens: filename and index [scale] ---')
        endif
     ! READ_F[ield] <filename> <index>
      elseif (token(1)(1:10).eq.'SET_LOOP_F') then
        if (num_token.gt.1) then
          write(os,"(a)")'> SET_LOOP_FILE '//trim(token(2))//' '//trim(layer_s); call output(os)
          call set_data(layer_n,trim(token(2)))
        else
          call output('!--- SET_LOOP_FILE requires filename ---')
        endif
     ! CREATE_F[ield] <index> <hx> <hy> <hz>
      elseif (token(1)(1:8).eq.'CREATE_F') then
        if (num_token.gt.4) then
          write(os,"(a,i3,1p,3(x,e19.12))")'> CREATE_FIELD',nint(h(2)),h(3),h(4),h(5); call output(os)
          ok_p=uniform_fld(nint(h(2)),(/h(3),h(4),h(5)/))
        else
          call output('!--- CREATE_FIELD requires four tokens: index hx hy hz ---')
        endif
     ! CREATE_T[emp] <index> <T0>
     ! CREATE_T[emp] <index> <dTmax> <sigma_x> <sigma_y> <exp n>
      elseif (token(1)(1:8).eq.'CREATE_T') then
        if (num_token.gt.5) then
          write(os,"(a,i3,1p,3(x,e19.12),x,0p,i3)")'> CREATE_TEMP',nint(h(2)),h(3),h(4),h(5),nint(h(6)); call output(os)
          ok_p=suprgaus_temp(nint(h(2)),h(3),h(4),h(5),nint(h(6)))
        elseif (num_token.gt.2) then
          write(os,"(a,i3,1p,4(x,e19.12))")'> CREATE_TEMP',nint(h(2)),h(3); call output(os)
          ok_p=suprgaus_temp(nint(h(2)),h(3),0.d0,0.d0,0)
        else
          call output('!--- CREATE_TEMP requires five tokens: index dTmax sigma_x sigma_y exp_n ---')
        endif
     ! LLG_MAXM <component index> <value>
     ! LLG_MAXM <mx> <my> <mz>
      elseif (token(1)(1:8).eq.'LLG_MAXM') then
        if (num_token.eq.3) then
          write(os,"(a,i3,1p,x,e19.12,x,a)") '> LLGMAXM',nint(h(2)),h(3),trim(layer_s); call output(os)
          call set_llg(layer_n,indx=nint(h(2)),max_ave_m=(/h(3),h(3),h(3)/))
        elseif (num_token.eq.4) then
          write(os,"(a,1p,3(x,e19.12),x,a)") '> LLGMAXM',h(2),h(3),h(4),trim(layer_s); call output(os)
          call set_llg(layer_n,max_ave_m=(/h(2),h(3),h(4)/))
        else
          call output('!--- LLG_MAXM requires additional tokens: component_index value OR mx my mz ---')
        endif
     ! LLG_MINM <component index> <value>
     ! LLG_MINM <mx> <my> <mz>
      elseif (token(1)(1:8).eq.'LLG_MINM') then
        if (num_token.eq.3) then
          write(os,"(a,i3,1p,x,e19.12,x,a)") '> LLGMINM',nint(h(2)),h(3),trim(layer_s); call output(os)
          call set_llg(layer_n,indx=nint(h(2)),min_ave_m=(/h(3),h(3),h(3)/))
        elseif (num_token.eq.4) then
          write(os,"(a,1p,3(x,e19.12),x,a)") '> LLGMINM',h(2),h(3),h(4),trim(layer_s); call output(os)
          call set_llg(layer_n,min_ave_m=(/h(2),h(3),h(4)/))
        else
          call output('!--- LLG_MINM requires additional tokens: component_index value OR mx my mz ---')
        endif
     ! SET_TE[mp] <index>
      elseif (token(1)(1:6).eq.'SET_TE') then
        if (num_token.gt.1) then
          write(os,"(a,i3,x,a)") '> SET_TEMP',nint(h(2)),trim(layer_s); call output(os)
          ok_p=set_temp(nint(h(2)),ilayer=layer_n)
          if (.not.ok_p) call output('!--- SET_TEMP failed ---')
        else
          call output('!--- SET_TEMP requires additional token: index ---')
        endif
     ! SET_1T[emp] <index> [scale]
     ! SET_1T[emp] <index> <start_time> <end_time> <start_scale> <end_scale> [scale]
      elseif (token(1)(1:5).eq.'SET_1') then
        v_p=(token(1)(1:6).eq.'SET_1S'); ok_p=(token(1)(1:6).eq.'SET_1T')
        if (ok_p) then; s2='SET_1TEMP'; elseif (v_p) then; s2='SET_1SOURCE';else;s2='SET_1FIELD'; endif
        if (num_token.eq.2.or.num_token.eq.3) then
          if (num_token.eq.2) h(3)=1.d0
          write(os,"(a,i3,'  -1.d22 1.d22 1.d0 1.d0 ',1p,e19.12,x,a)") '> '//trim(s2),nint(h(2)),h(3),trim(layer_s); call output(os)
          if (ok_p.or.v_p) then
            v_p=set_source(v_p,ilayer=layer_n)
            ok_p=set_temp(nint(h(2)),scale=h(3),ilayer=layer_n)
          else
            ok_p=set_fld(nint(h(2)),scale=h(3),ilayer=layer_n)
          endif
          if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
        elseif (num_token.gt.5) then
          if (token(3)(1:3).eq.'STE') then  !step (w/ ramp)
            if (num_token.lt.7) h(7)=0.d0
            if (num_token.lt.8) h(8)=0.d0
            if (num_token.lt.9) h(9)=1.d0
            if (num_token.lt.13) h(13)=h(9)
            if (num_token.lt.12) h(9:12)=1.d0
            if (num_token.lt.12) h(9:10)=(/0.d0,1.d0/)

            write(os,"(a,i3,x,a,1p,6(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),trim(token(3)),(h(i),i=4,9),trim(layer_s); call output(os)
            if (ok_p.or.v_p) then
              v_p=set_source(v_p,ilayer=layer_n)
              ok_p=set_temp(nint(h(2)),func1='STE', gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), param1=h(8), &
                      t_start=h(9),t_end=h(10),s1_start=h(11),s1_end=h(12),scale=h(13),ilayer=layer_n)
            else
              ok_p=set_fld(nint(h(2)),func1='STE', gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), param1=h(8), &
                      t_start=h(9),t_end=h(10),s1_start=h(11),s1_end=h(12),scale=h(13),ilayer=layer_n)
            endif
            if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
          elseif (token(3)(1:3).eq.'NEW') then  !newton cooling
            if (num_token.lt.7) h(7)=0.d0
            if (num_token.lt.8) h(8)=1.d-9
            if (num_token.lt.9) h(9)=0.5d0
            if (num_token.lt.10) h(10)=1.d0
            if (num_token.lt.11) h(11)=1.d0
            if (num_token.lt.15) h(15)=h(11)
            if (num_token.lt.14) h(11:14)=1.d0
            if (num_token.lt.12) h(11:12)=(/0.d0,1.d0/)
            write(os,"(a,i3,x,a,1p,8(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),trim(token(3)),(h(i),i=4,11),trim(layer_s); call output(os)
            h(10)=min(1.d0,max(0.d0,abs(h(10)))); h(9)=min(1.d0,max(0.d0,abs(h(9)))); h(8)=abs(h(8)); if (h(8).le.0.d0) h(8)=1.d-12
!           h(12)=exp(-(1.d0-h(9))*h(4)/h(8)); h(13)=exp(h(9)*h(4)/h(8))
!           h(6)=h(6)*(h(13)-h(12))/((1.d0+h(10))*h(13)-2.d0*h(10)-(1.d0-h(10))*h(12))

            h(6)=h(6)*(exp(h(4)/h(8))-1.d0)/(exp(h(4)/h(8))-1.d0+h(10)*(1.d0-exp(h(4)*(1.d0-h(9))/h(8))))
            if (ok_p.or.v_p) then
              v_p=set_source(v_p,ilayer=layer_n)
              ok_p= set_temp( nint(h(2)), func1='NEW', gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), param1=h(8), duty1=h(9), &
                      moddepth1=h(10),t_start=h(11),t_end=h(12),s1_start=h(13),s1_end=h(14),scale=h(15),ilayer=layer_n)
            else
              ok_p= set_fld( nint(h(2)), func1='NEW', gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), param1=h(8), duty1=h(9), &
                      moddepth1=h(10),t_start=h(11),t_end=h(12),s1_start=h(13),s1_end=h(14),scale=h(15),ilayer=layer_n)
            endif
            if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
          elseif (token(3)(1:3).eq.'SIN'.or.token(3)(1:3).eq.'COS'.or.token(3)(1:3).eq.'EXP'.or.token(3)(1:3).eq.'ROT') then  !newton cooling
            if (num_token.lt.7) h(7)=0.d0
            if (num_token.lt.8) h(8)=1.d0
            if (num_token.lt.12) h(12)=h(8)
            if (num_token.lt.11) h(8:11)=1.d0
            if (num_token.lt.11) h(8:9)=(/0.d0,1.d0/)

            write(os,"(a,i3,x,a,1p,5(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),trim(token(3)),(h(i),i=4,8),trim(layer_s); call output(os)
            if (ok_p.or.v_p) then
              v_p=set_source(v_p,ilayer=layer_n)
              ok_p=set_temp(nint(h(2)),func1=token(3)(1:3), gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), &
                      t_start=h(8),t_end=h(9),s1_start=h(10),s1_end=h(11),scale=h(12),ilayer=layer_n)
            else
              ok_p=set_fld(nint(h(2)),func1=token(3)(1:3), gamma1=h(4), phase1=h(5), amp1=h(6), dc1=h(7), &
                      t_start=h(8),t_end=h(9),s1_start=h(10),s1_end=h(11),scale=h(12),ilayer=layer_n)
            endif
            if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
          else
            if (num_token.lt.7) h(7)=1.d0
            write(os,"(a,i3,1p,5(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),(h(i),i=3,7),trim(layer_s); call output(os)
            if (ok_p.or.v_p) then
              v_p=set_source(v_p,ilayer=layer_n)
              ok_p=set_temp(nint(h(2)),t_start=h(3),t_end=h(4),s1_start=h(5),s1_end=h(6),scale=h(7),ilayer=layer_n)
            else
              ok_p=set_fld(nint(h(2)),t_start=h(3),t_end=h(4),s1_start=h(5),s1_end=h(6),scale=h(7),ilayer=layer_n)
            endif
            if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
          endif
        else
          write(os,"('!--- ',a,' requires at least one additional token: index ---')") trim(s2);call output(os)
        endif
     ! SET_2T[emp] <index> [scale]
     ! SET_2T[emp] <index> <start_time> <end_time> <start_scale> <end_scale> [scale]
      elseif (token(1)(1:5).eq.'SET_2') then
        ok_p=(token(1)(1:6).eq.'SET_2T'); if (ok_p) then; s2='SET_2TEMP'; else; s2='SET_2FIELD'; endif
       ! set_2* index1 index2 [scale] num_token=3 or 4
        if (num_token.eq.3.or.num_token.eq.4) then
          if (num_token.eq.3) h(4)=1.d0
          write(os,"(a,i3,i3,'  -1.d22 1.d22 1.d0 1.d0 ',1p,e19.12,x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),h(4),trim(layer_s); call output(os)
          if (ok_p) then
            ok_p=set_temp(nint(h(2)),ifld2=nint(h(3)),scale=h(4),ilayer=layer_n)
          else
            ok_p=set_fld(nint(h(2)),ifld2=nint(h(3)),scale=h(4),ilayer=layer_n)
          endif
          if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
        else
          i=0
          do j=3,num_token
            s3=token(j)(1:3)
            if (s3(1:3).eq.'NEW'.or.s3(1:3).eq.'STE'.or.s3(1:3).eq.'SIN'.or.s3(1:3).eq.'COS'.or.s3(1:3).eq.'EXP'.or.s3(1:3).eq.'ROT')i=i+1
          enddo
          if (i.eq.0) then
       ! set_2* index1 index2 start_time end_time start_scale1 end_scale1 [scale], num_token=7 or 8
       ! set_2* index1 index2 start_time end_time start_scale1 end_scale1 start_scale2 end_scale2 [scale], num_token=9 or 10
            if (num_token.lt.9) then
              if (num_token.lt.8) h(8)=1.d0
              write(os,"(a,i3,i3,1p,5(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),(h(i),i=4,8),trim(layer_s); call output(os)
              if (ok_p) then
                ok_p=set_temp(nint(h(2)),ifld2=nint(h(3)),t_start=h(4),t_end=h(5),s1_start=h(6),s1_end=h(7),scale=h(8),ilayer=layer_n)
              else
                ok_p=set_fld(nint(h(2)),ifld2=nint(h(3)),t_start=h(4),t_end=h(5),s1_start=h(6),s1_end=h(7),scale=h(8),ilayer=layer_n)
              endif
              if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
            else
              if (num_token.lt.10) h(10)=1.d0
              write(os,"(a,i3,i3,1p,7(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),(h(i),i=4,10),trim(layer_s); call output(os)
              if (ok_p) then
                ok_p=set_temp(nint(h(2)),ifld2=nint(h(3)),t_start=h(4),t_end=h(5),s1_start=h(6),s1_end=h(7), &
                                s2_start=h(8),s2_end=h(9),scale=h(10),ilayer=layer_n)
              else
                ok_p=set_fld(nint(h(2)),ifld2=nint(h(3)),t_start=h(4),t_end=h(5),s1_start=h(6),s1_end=h(7), &
                                s2_start=h(8),s2_end=h(9),scale=h(10),ilayer=layer_n)
              endif
              if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
            endif
          elseif (i.eq.1) then
            do j=3,num_token
              s3=token(j)(1:3)
              if (s3(1:3).eq.'NEW'.or.s3(1:3).eq.'STE'.or.s3(1:3).eq.'SIN'.or.s3(1:3).eq.'COS'.or.s3(1:3).eq.'EXP'.or.s3(1:3).eq.'ROT')i=j
            enddo
           !set_2* index1 index2 STE gamma1 phase1 amp1 [dc1] [param1] [scale], num_token=7-10
           !set_2* index1 index2 STE gamma1 phase1 amp1 [dc1] [param1] start_time end_time start_scale2 end_scale2 [scale], num_token=11-14
           !set_2* index1 index2 SCE gamma1 phase1 amp1 [dc1] [scale], num_token=7-9
           !set_2* index1 index2 SCE gamma1 phase1 amp1 [dc1] start_time end_time start_scale2 end_scale2 [scale], num_token=11-13
!           if (i.ne.4) then
            if (i.eq.4) then
             !set_2* index1 index2 NEW gamma1 phase1 amp1 [dc1] [param1] [duty1] [mod] [scale], num_token=7-12
             !set_2* index1 index2 NEW gamma1 phase1 amp1 dc1 param1 duty1 mod1 start_time end_time start_scale2 end_scale2 scale, num_token=16
             !set_2* index1 index2 NEW gamma1 phase1 amp1 dc1 param1 duty1 mod1 start_time end_time start_scale2 end_scale2      , num_token=15
             !set_2* index1 index2 NEW gamma1 phase1 amp1 dc1 param1 duty1      start_time end_time start_scale2 end_scale2      , num_token=14
             !set_2* index1 index2 NEW gamma1 phase1 amp1 dc1 param1            start_time end_time start_scale2 end_scale2      , num_token=13
              if (token(4)(1:3).eq.'NEW') then
                if (num_token.lt.13) then
                  if (num_token.lt.8) h(8)=0.d0
                  if (num_token.lt.9) h(9)=1.d-9
                  if (num_token.lt.10) h(10)=0.5d0
                  if (num_token.lt.11) h(11)=1.d0
                  if (num_token.lt.12) h(12)=1.d0
                  write(os,"(a,i3,i3,x,a,1p,8(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,12),trim(layer_s); call output(os)
                  h(11)=min(1.d0,max(0.d0,abs(h(11)))); h(10)=min(1.d0,max(0.d0,abs(h(10)))); h(9)=abs(h(9)); if (h(9).le.0.d0) h(9)=1.d-12
!                 h(13)=exp(-(1.d0-h(10))*h(5)/h(9)); h(14)=exp(h(10)*h(5)/h(9))
!                 h(7)=h(7)*(h(14)-h(13))/((1.d0+h(11))*h(14)-2.d0*h(11)-(1.d0-h(11))*h(13))

                  h(7)=h(7)*(exp(h(5)/h(9))-1.d0)/(exp(h(5)/h(9))-1.d0+h(11)*(1.d0-exp(h(5)*(1.d0-h(10))/h(9))))
                  if (ok_p) then
                    ok_p= set_temp( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                         moddepth1=h(11), scale=h(12) ,ilayer=layer_n)
                  else
                    ok_p= set_fld( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                         moddepth1=h(11), scale=h(12) ,ilayer=layer_n)
                  endif
                  if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
                else
                  if (num_token.lt.14) then
                    do j=13,10,-1; h(j+2)=h(j); enddo; h(10)=0.5d0; h(11)=1.d0; h(16)=1.d0
                  elseif (num_token.lt.15) then
                    do j=14,11,-1; h(j+1)=h(j); enddo; h(11)=1.d0; h(16)=1.d0
                  elseif (num_token.lt.16) then
                    h(16)=1.d0
                  endif
                  h(11)=min(1.d0,max(0.d0,abs(h(11)))); h(10)=min(1.d0,max(0.d0,abs(h(10)))); h(9)=abs(h(9)); if (h(9).le.0.d0) h(9)=1.d-12
!                 h(17)=exp(-(1.d0-h(10))*h(5)/h(9)); h(18)=exp(h(10)*h(5)/h(9))
!                 h(7)=h(7)*(h(18)-h(17))/((1.d0+h(11))*h(18)-2.d0*h(11)-(1.d0-h(11))*h(17))
                  h(7)=h(7)*(exp(h(5)/h(9))-1.d0)/(exp(h(5)/h(9))-1.d0+h(11)*(1.d0-exp(h(5)*(1.d0-h(10))/h(9))))
                  write(os,"(a,i3,i3,x,a,1p,12(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,16),trim(layer_s); call output(os)
                  if (ok_p) then
                    ok_p= set_temp( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                         moddepth1=h(11), t_start=h(12), t_end=h(13), s2_start=h(14), s2_end=h(15), scale=h(16) ,ilayer=layer_n)
                  else
                    ok_p= set_fld( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                         moddepth1=h(11), t_start=h(12), t_end=h(13), s2_start=h(14), s2_end=h(15), scale=h(16) ,ilayer=layer_n)
                  endif
                  if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
                endif
             !set_2* index1 index2 STE gamma1 phase1 amp1 [dc1] [param1]                                            [scale], num_token=7-10
             !set_2* index1 index2 STE gamma1 phase1 amp1  dc1   param1 start_time end_time start_scale2 end_scale2  scale , num_token=14
             !set_2* index1 index2 STE gamma1 phase1 amp1  dc1   param1 start_time end_time start_scale2 end_scale2        , num_token=13
             !set_2* index1 index2 STE gamma1 phase1 amp1  dc1          start_time end_time start_scale2 end_scale2        , num_token=12
             !set_2* index1 index2 STE gamma1 phase1 amp1               start_time end_time start_scale2 end_scale2        , num_token=11
              elseif (token(4)(1:3).eq.'STE') then
                if (num_token.lt.11) then
                  if (num_token.lt.8) h(8)=0.d0
                  if (num_token.lt.9) h(9)=0.d0
                  if (num_token.lt.10) h(10)=1.d0
                  write(os,"(a,i3,i3,x,a,1p,6(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,10),trim(layer_s); call output(os)
                  if (ok_p) then
                    ok_p= set_temp( nint(h(2)), nint(h(3)), func1='STE', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), scale=h(10) ,ilayer=layer_n)
                  else
                    ok_p= set_fld( nint(h(2)), nint(h(3)), func1='STE', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), scale=h(10) ,ilayer=layer_n)
                  endif
                  if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
                else
                  if (num_token.lt.12) then
                    do j=11,8,-1; h(j+2)=h(j); enddo; h(8)=0.0d0; h(9)=0.d0; h(14)=1.d0
                  elseif (num_token.lt.13) then
                    do j=12,9,-1; h(j+1)=h(j); enddo; h(9)=0.d0; h(14)=1.d0
                  elseif (num_token.lt.14) then
                    h(14)=1.d0
                  endif
                  write(os,"(a,i3,i3,x,a,1p,10(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,14),trim(layer_s); call output(os)
                  if (ok_p) then
                    ok_p= set_temp( nint(h(2)), nint(h(3)), func1='STE', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), &
                         t_start=h(10), t_end=h(11), s2_start=h(12), s2_end=h(13), scale=h(14) ,ilayer=layer_n)
                  else
                    ok_p= set_fld( nint(h(2)), nint(h(3)), func1='STE', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), &
                         t_start=h(10), t_end=h(11), s2_start=h(12), s2_end=h(13), scale=h(14) ,ilayer=layer_n)
                  endif
                  if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
                endif
              elseif (token(4)(1:3).eq.'COS'.or.token(4)(1:3).eq.'SIN'.or.token(4)(1:3).eq.'EXP'.or.token(4)(1:3).eq.'ROT') then

           !set_2* index1 index2 SCE gamma1 phase1 amp1 [dc1] [scale], num_token=7-9
           !set_2* index1 index2 SCE gamma1 phase1 amp1  dc1  start_time end_time start_scale2 end_scale2 [scale], num_token=12-13
           !set_2* index1 index2 SCE gamma1 phase1 amp1  dc1  start_time end_time start_scale1 end_scale1 start_scale2 end_scale2 [scale], num_token=14-15

                if (num_token.lt.8) h(8)=0.d0
                if (num_token.lt.9) h(9)=1.d0
                if (num_token.lt.13) h(13)=h(9)
                if (num_token.lt.12) h(9:12)=1.d0
                if (num_token.lt.12) h(9:10)=(/0.d0,1.d0/)
                if (num_token.lt.15) h(15)=h(13)
                if (num_token.lt.14) h(13:14)=1.d0
                if (num_token.lt.14.and.num_token.gt.11) then; h(13:14)=h(11:12); h(11:12)=1.d0; endif


                  write(os,"(a,i3,i3,x,a,1p,11(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,15),trim(layer_s); call output(os)
                  if (ok_p) then
                    ok_p= set_temp( nint(h(2)), nint(h(3)), func1=token(4)(1:3), gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), &
                         t_start=h(9),t_end=h(10),s1_start=h(11),s1_end=h(12),s2_start=h(13),s2_end=h(14),scale=h(15),ilayer=layer_n)
                  else
                    ok_p= set_fld( nint(h(2)), nint(h(3)), func1=token(4)(1:3), gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), &
                         t_start=h(9),t_end=h(10),s1_start=h(11),s1_end=h(12),s2_start=h(13),s2_end=h(14),scale=h(15),ilayer=layer_n)
                  endif
                  if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
              else
                write(os,"('!--- this ',a,' requires fourth token to be: NEW|STE|COS|SIN|EXP|ROT ---')") trim(s2);call output(os)
              endif
            else
              write(os,"('!--- this ',a,' requires fourth token to be: NEW|STE|COS|SIN|EXP|ROT ---')") trim(s2);call output(os)
            endif
          elseif (i.eq.2) then
            i=0
            do j=3,num_token
              s3=token(j)
              if (s3(1:3).eq.'NEW'.or.s3(1:3).eq.'STE'.or.s3(1:3).eq.'SIN'.or.s3(1:3).eq.'COS'.or.s3(1:3).eq.'EXP'.or.s3(1:3).eq.'ROT') then
                if (i.ne.0) exit; i=j
              endif
            enddo
            os=token(i)(1:3); s3=token(j)
            if ((os(1:3).eq.'COS'.or.os(1:3).eq.'SIN'.or.os(1:3).eq.'EXP'.or.os(1:3).eq.'ROT').and. &
                (s3(1:3).eq.'COS'.or.s3(1:3).eq.'SIN'.or.s3(1:3).eq.'EXP'.or.s3(1:3).eq.'ROT').and.i.eq.4.and.(j.eq.9.or.j.eq.8)) then
             !set_2* index1 index2 SCE gamma1 phase1 amp1 [dc1] SCE gamma1 phase1 amp1 [dc1] [scale], num_token=11-14
              if (j.eq.8) then
                do k=num_token,8,-1; h(j+1)=h(j); enddo; h(8)=0.d0; num_token=num_token+1
              endif
              if (num_token.lt.13) h(13)=0.d0
              if (num_token.lt.14) h(14)=1.d0
              write(os,"(a,i3,i3,x,a,1p,4(x,e19.12),x,a,5(xe19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,8), &
                    trim(token(9)),(h(i),i=10,14),trim(layer_s); call output(os)
              if (ok_p) then
                ok_p= set_temp( nint(h(2)), nint(h(3)), func1=token(4)(1:3), gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), &
                     func2=token(9)(1:3), gamma2=h(10), phase2=h(11), amp2=h(12), dc2=h(13), scale=h(14) ,ilayer=layer_n)
              else
                ok_p= set_fld( nint(h(2)), nint(h(3)), func1=token(4)(1:3), gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), &
                     func2=token(9)(1:3), gamma2=h(10), phase2=h(11), amp2=h(12), dc2=h(13), scale=h(14) ,ilayer=layer_n)
              endif
              if (.not.ok_p) then; write(os,"('!--- ',a,' failed ---')") trim(s2);call output(os); endif
            elseif (token(i)(1:3).eq.'NEW'.and.token(j)(1:3).eq.'NEW'.and.(num_token.eq.19.or.num_token.eq.20)) then
     !set_2* index1 index2 NEW gamma1 phase1 amp1 dc1 param1 duty1 mod1 NEW gamma2 phase2 amp2 dc2 param2 duty2 mod2 [scale], num_token=7-12
     !  1      2      3     4    5      6     7    8    9     10    11  12   13     14     15   16   17    18    19    20
              if (num_token.lt.20) h(20)=1.d0
              write(os,"(a,i3,i3,x,a,1p,7(x,e19.12),x,a,8(x,e19.12),x,a)") '> '//trim(s2),nint(h(2)),nint(h(3)),trim(token(4)),(h(i),i=5,11),trim(token(12)),&
                       (h(i),i=13,20),trim(layer_s); call output(os)
              h(11)=min(1.d0,max(0.d0,abs(h(11)))); h(10)=min(1.d0,max(0.d0,abs(h(10)))); h(9)=abs(h(9)); if (h(9).le.0.d0) h(9)=1.d-12
              h(19)=min(1.d0,max(0.d0,abs(h(19)))); h(18)=min(1.d0,max(0.d0,abs(h(18)))); h(17)=abs(h(17)); if (h(17).le.0.d0) h(17)=1.d-12
              h(15)=h(15)/h(7)
              h(7)=h(7)*(exp(h(5)/h(9))-1.d0)/(exp(h(5)/h(9))-1.d0+h(11)*(1.d0-exp(h(5)*(1.d0-h(10))/h(9))))
              h(15)=h(15)*h(7)
              if (ok_p) then
                ok_p= set_temp( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                     moddepth1=h(11),                   func2='NEW', gamma2=h(13), phase2=h(14), amp2=h(15), dc2=h(16), param2=h(17), duty2=h(18), &
                     moddepth2=h(19), scale=h(20) ,ilayer=layer_n)
              else
                ok_p= set_fld( nint(h(2)), nint(h(3)), func1='NEW', gamma1=h(5), phase1=h(6), amp1=h(7), dc1=h(8), param1=h(9), duty1=h(10), &
                     moddepth1=h(11),                   func2='NEW', gamma2=h(13), phase2=h(14), amp2=h(15), dc2=h(16), param2=h(17), duty2=h(18), &
                     moddepth2=h(19), scale=h(20) ,ilayer=layer_n)
              endif
            else
              write(os,"('!--- this two func ',a,' form not implemented ---')") trim(s2);call output(os)
            endif
          else
            write(os,"('!--- ',a,' requires at least one additional token: index ---')") trim(s2);call output(os)
          endif
        endif
      elseif (token(1)(1:5).eq.'SET_D') then
        if (num_token.eq.3) then
          ok_p=(token(3)(1:1).eq.'T'.or.token(3)(1:2).eq.'.T')
          if (token(2)(1:1).eq.'N') then; write(os,"('NTRAJ ',i0)") nint(h(3)); ok_p=set_decay_mag(intraj=nint(h(3))); endif
          if (token(2)(1:1).eq.'T') then; write(os,"('TRAJ ',i0)") nint(h(3)); ok_p=set_decay_mag(itraj=nint(h(3))); endif
          if (token(2)(1:1).eq.'O') then; write(os,"('ORBIT ',1p,e12.5)") h(3); ok_p=set_decay_mag(itorbit=h(3)); endif
          if (token(2)(1:1).eq.'M') then; write(os,"('MAX_EB ',1p,e12.5)") h(3); ok_p=set_decay_mag(imax_eb=h(3)); endif
          if (token(2)(1:1).eq.'I') then; write(os,"('ITER ',i0)") nint(h(3)); ok_p=set_decay_mag(iiter=nint(h(3))); endif
          if (token(2)(1:1).eq.'V') then; write(os,"('VOL_FRAD ',1p,e12.5)") h(3); ok_p=set_decay_mag(ivol_frac=h(3)); endif
          if (token(2)(1:1).eq.'D') then; write(os,"('DIFF_M ',l)") ok_p; ok_p=set_decay_mag(idiff_m=ok_p); endif
          if (token(2)(1:1).eq.'R') then; write(os,"('ROT ',l)") ok_p; ok_p=set_decay_mag(itest_rot_p=ok_p); endif
          if (token(2)(1:1).eq.'B') then; write(os,"('BOLTZ ',l)") ok_p; ok_p=set_decay_mag(iboltz=ok_p); endif
          call output('SET_DECAY '//trim(os))
        else
          write(os,"('!--- SET_DECAY requires at least two additional tokens: ntraj|traj|orbit|max_eb|iter|vol_frac|diff_m|rot|boltz and value ---')") ;call output(os)
        endif
      elseif (token(1)(1:5).eq.'DECAY') then
! subroutine decay_mag(t0,dt,temperature,dh_dti*,happliedi*,mkti*,old_gamma_pi*,multispin_pi*,iseedi*), *-optional parameters
        if (num_token.gt.2) then
          if (num_token.gt.8) then
            ok_p=(token(7)(1:1).eq.'T'.or.token(7)(1:2).eq.'.T'); v_p=(token(8)(1:1).eq.'T'.or.token(8)(1:2).eq.'.T')
            write(os,"('> DECAY ',1p,5(e19.12,x),2(1l,x),0p,i0)")h(2),h(3),h(4),h(5),h(6),ok_p,v_p,nint(h(9)); call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4),happliedi=h(5),mkti=h(6),old_gamma_pi=ok_p,multispin_pi=v_p,iseedi=nint(h(9)))
          elseif (num_token.gt.7) then
            ok_p=(token(7)(1:1).eq.'T'.or.token(7)(1:2).eq.'.T'); v_p=(token(8)(1:1).eq.'T'.or.token(8)(1:2).eq.'.T') 
            write(os,"('> DECAY ',1p,5(e19.12,x),2(1l,x))")h(2),h(3),h(4),h(5),h(6),ok_p,v_p; call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4),happliedi=h(5),mkti=h(6),old_gamma_pi=ok_p,multispin_pi=v_p)
          elseif (num_token.gt.6) then 
            ok_p=(token(7)(1:1).eq.'T'.or.token(7)(1:2).eq.'.T')
            write(os,"('> DECAY ',1p,5(e19.12,x),1(1l,x))")h(2),h(3),h(4),h(5),h(6),ok_p; call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4),happliedi=h(5),mkti=h(6),old_gamma_pi=ok_p)
          elseif (num_token.gt.5) then
            write(os,"('> DECAY ',1p,5(e19.12,x))")h(2),h(3),h(4),h(5),h(6); call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4),happliedi=h(5),mkti=h(6))
          elseif (num_token.gt.4) then
            write(os,"('> DECAY ',1p,4(e19.12,x))")h(2),h(3),h(4),h(5); call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4),happliedi=h(5))
          elseif (num_token.gt.3) then
            write(os,"('> DECAY ',1p,3(e19.12,x))")h(2),h(3),h(4); call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3),dh_dti=h(4))
          else
            write(os,"('> DECAY ',1p,2(e19.12,x))")h(2),h(3); call output(os)
            if (line_read.gt.line_target) call decay_mag(time,h(2),h(3))
          endif
          if (line_read.gt.line_target) then; if (.not.my_checkpoint(.false.)) exit; endif
          if (line_read.eq.line_target) then
            if (.not.my_checkpoint(.true.)) return
            time=time_target; call update_do_grain(time)
            call output(' CHK: found last executed cli command, restarted!')
          endif
        else
          write(os,"('!--- DECAY requires at least two additional tokens: dt and temperature ---')") ;call output(os)
        endif

      elseif (token(1)(1:6).eq.'INIT_T') then
        call output('> INIT_TENSOR')
        ok_p=init_demag()
        if (.not.ok_p) call output('!--- INIT_TENSOR failed ---')
      elseif (token(1)(1:7).eq.'START_C') then
        call output('> START_CLOCK'); call cpu_time(start_time)
      elseif (token(1)(1:6).eq.'STOP_C') then
        call cpu_time(stop_time); i=int((stop_time-start_time)/3600); j=int(((stop_time-start_time)-i*3600)/60)
        write(os,"('> STOP_CLOCK: ',i3,' hours, ',i2,' minutes, ',f5.2,' seconds')") i,j,(stop_time-start_time)-60*(j+60*i); call output(os)
      elseif (token(1)(1:9).eq.'SET_CHK_T') then
        if (num_token.gt.1) then
          write(os,"(a,x,1p,e19.12)") '> SET_CHK_TIME',h(2); call output(os); if (h(2).ge.1) ok_p=checkpoint(icheckpoint_time=real(abs(h(2))))
        else
          call output('!--- SET_CHK_TIME requires another token: <seconds> ---')
        endif
      elseif (token(1)(1:8).eq.'SET_TIME') then
        if (num_token.gt.1) then
          write(os,"(a,x,1p,e19.12)") '> SET_TIME',h(2); call output(os); time=h(2)
          if (line_read.gt.line_target) call update_do_grain(time)
        else
          call output('!--- SET_TIME requires another token: time ---')
        endif
      elseif (token(1)(1:7).eq.'INPUT_K') then
        if (key_p) then
          call output('!--- INPUT_KEYBOARD is currently active---')
        else
          call output('> INPUT_KEYBOARD'); key_p=.true.; ioc_num=ioc_num+1; ioc(ioc_num)=0
        endif
      elseif (token(1)(1:7).eq.'INPUT_F') then
        if (num_token.eq.1)then  !backup control file
          if (ioc(ioc_num).gt.0) then; close(ioc(ioc_num)); call delete_file(ioc(ioc_num)); endif
          ioc_num=ioc_num-1; if (ioc_num.lt.1) exit
          cliionum=ioc(ioc_num)
          if (cliionum.lt.1) then
            call output('> INPUT_KEYBOARD'); key_p=.true.
          else
            write(os,"(a)") '> INPUT_FILE '//trim(adjustl(find_file(cliionum)))//' (continue)'; call output(os)
            key_p=.false.;
          endif
        else
          if (trim(adjustl(token(2))).eq.'-') then
            call output('> INPUT_KEYBOARD'); key_p=.true.; ioc_num=ioc_num+1; ioc(ioc_num)=0
          elseif (find_file(trim(adjustl(token(2)))).gt.0) then
            cliionum=find_file(trim(adjustl(token(2))))
            write(os,"(a)") '> INPUT_FILE '//trim(adjustl(token(2)))//' (continue)'; call output(os)
            key_p=.false.
          else
            inquire(file=trim(adjustl(token(2))),exist=ok_p)
            if (ok_p) then
              write(os,"(a)") '> INPUT_FILE '//trim(adjustl(token(2))); call output(os)
              cliionum=add_file(trim(token(2)))
              open(cliionum,file=trim(token(2)),status='old',action='read')
              key_p=.false.; 
              if (ioc_num.eq.size(ioc)) then; do j=2,size(ioc); ioc(j-1)=ioc(j); enddo; endif
              ioc_num=min(ioc_num+1,size(ioc)); ioc(ioc_num)=cliionum
            else
              write(os,"(a)") '!--- no such file exists for INPUT_FILE: '//trim(adjustl(token(2))); call output(os)
            endif
          endif
        endif
      elseif (token(1)(1:5).eq.'GEN_G') then
        call output('> GEN_GRAIN'); ok_p=.true.;if (line_read.gt.line_target) ok_p= generate_media(); if (.not.ok_p) call output('!--- GEN_GRAIN failed ---')
        ok_p=checkpoint(isave_grain_p=.true.)
      elseif (token(1)(1:6).eq.'GEN_HK') then
        call output('> GEN_HK'); ok_p=.true.;if (line_read.gt.line_target) ok_p= generate_hk(); if (.not.ok_p) call output('!--- GEN_HK failed ---')
        ok_p=checkpoint(isave_grain_p=.true.)
      elseif (token(1)(1:6).eq.'GEN_AR') then
        call output('> GEN_AREA'); ok_p=.true.;if (line_read.gt.line_target) ok_p= generate_area(); if (.not.ok_p) call output('!--- GEN_AREA failed ---')
        ok_p=checkpoint(isave_grain_p=.true.)
      elseif (token(1)(1:6).eq.'INIT_M') then
        call output('> INIT_MAG'); ok_p=.true.;if (line_read.gt.line_target) ok_p= generate_mag(); if (.not.ok_p) call output('!--- INIT_MAG failed ---')
      elseif (token(1)(1:4).eq.'LLG ') then
        call output('> LLG'); if (line_read.gt.line_target) call llg(time)
        if (line_read.gt.line_target) then; if (.not.my_checkpoint(.false.)) exit; endif
        if (line_read.eq.line_target) then
          if (.not.my_checkpoint(.true.)) return
          time=time_target; call update_do_grain(time)
          call output(' CHK: found last executed cli command, restarted!')
        endif
      elseif (token(1)(1:6).eq.'FIND_G') then
       if (num_token.gt.2) then
         write(os,"(a,x,i0,x,i0,x,i0)") '> FIND_GRAIN',nint(h(2)),nint(h(3)),nint(h(4)); call output(os)
         call find_grain(nint(h(4)),nint(h(2)),nint(h(3)))
       elseif (num_token.eq.2) then
         write(os,"(a,x,i0,x,i0)") '> FIND_GRAIN',nint(h(2)),nint(h(3)); call output(os)
         call find_grain(0,nint(h(2)),nint(h(3)))
       else
         write(os,"(a)") '!--- FIND_Grain needs grain i j [and layer] ---'; call output(os)
       endif
      elseif (token(1)(1:5).eq.'SHOW_') then
        if (token(1)(1:12).eq.'SHOW_FIELD_S') then
          call output('> SHOW_FIELD_SET'); call show_fld_set()
        elseif (token(1)(1:6).eq.'SHOW_F') then
          call output('> SHOW_FIELD'); call show_fld()
        elseif (token(1)(1:11).eq.'SHOW_TEMP_S') then
          call output('> SHOW_TEMP_SET'); call show_temp_set()
        elseif (token(1)(1:6).eq.'SHOW_T') then
          call output('> SHOW_TEMP'); call show_temp()
        elseif (token(1)(1:6).eq.'SHOW_S') then
          call output('> SHOW_SCALINGS'); call show_temp_scaling()
        elseif (token(1)(1:7).eq.'SHOW_LL') then
          call output('> SHOW_LLG'); call show_llg(time)
        elseif (token(1)(1:6).eq.'SHOW_L') then
          call output('> SHOW_LAYERS'); call show_media()
        elseif (token(1)(1:6).eq.'SHOW_W') then
          call output('> SHOW_WINDOW'); call show_window()
        elseif (token(1)(1:6).eq.'SHOW_R') then
          call output('> SHOW_RANDOM'); call show_random()
        elseif (token(1)(1:6).eq.'SHOW_G') then
          if (num_token.gt.2) then
            write(os,"(a,x,i0,x,i0)") '> SHOW GRAIN',nint(h(2)),nint(h(3)); call output(os)
            call show_grain(nint(h(3)),nint(h(2)))
          elseif (num_token.gt.1) then
            write(os,"(a,x,i0,x,i0)") '> SHOW GRAIN',nint(h(2)); call output(os)
            call show_grain(0,nint(h(2)))
          else
            write(os,"(a)") '!--- SHOW GRAIN needs grain # [and layer] ---'; call output(os)
          endif
        else
          write(os,"(a)") '!--- '//trim(token(1))//' not implemented yet ---'; call output(os)
        endif
      elseif (token(1)(1:4).eq.'FLD_'.or.token(1)(1:5).eq.'TEMP_') then
        k=5; if (token(1)(1:1).eq.'T') k=6
        s2='FLD_'; if (k.eq.6) s2='TEMP_'
        ok_p=.false. !add to value?
        do i=num_token,4,-1
          if (token(i)(1:2).eq.'AD') ok_p=.true.
        enddo; s1=' '; if (ok_p) s1=' ADD'
        if (token(1)(k:k).eq.'V') then
          write(os,"(a,i3,2(x,e19.12),a)") '> '//trim(s2)//'VELOCITY',nint(h(2)),h(3),h(4),trim(s1); call output(os)
          if (k.lt.6) then
            ok_p=ref_fld(ifld=nint(h(2)),vel_x=h(3),vel_y=h(4),add=ok_p)
          else
            ok_p=ref_temp(ifld=nint(h(2)),vel_x=h(3),vel_y=h(4),add=ok_p)
          endif
          if (.not.ok_p) call output('!--- '//trim(s2)//'VELOCITY failed ---')
          call demag_window(v=(/h(3),h(4)/))
        elseif (token(1)(k:k).eq.'O') then
          write(os,"(a,i3,2(x,e19.12),a)") '> '//trim(s2)//'OFFSET',nint(h(2)),h(3),h(4),trim(s1); call output(os)
          if (k.lt.6) then
            ok_p=ref_fld(ifld=nint(h(2)),offset_x=h(3),offset_y=h(4),add=ok_p)
          else
            ok_p=ref_temp(ifld=nint(h(2)),offset_x=h(3),offset_y=h(4),add=ok_p)
          endif
          if (.not.ok_p) call output('!--- '//trim(s2)//'OFFSET failed ---')
        elseif (token(1)(k:k).eq.'R') then
          write(os,"(a,i3,2(x,e19.12),a)") '> '//trim(s2)//'REF',nint(h(2)),h(3),h(4),trim(s1); call output(os)
          if (k.lt.6) then
            ok_p=ref_fld(ifld=nint(h(2)),ref_x=h(3),ref_y=h(4),add=ok_p)
          else
            ok_p=ref_temp(ifld=nint(h(2)),ref_x=h(3),ref_y=h(4),add=ok_p)
          endif
          if (.not.ok_p) call output('!--- '//trim(s2)//'REF failed ---')
          call demag_window(offset=(/-h(3),-h(4)/))
        elseif (token(1)(k:k).eq.'D') then
          write(os,"(a,i3,i3,3(x,e12.5))") '> '//trim(s2)//'DUPLICATE',nint(h(2)),nint(h(3)),h(4),h(5),h(6); call output(os)
          if (k.lt.6) then
            ok_p=dup_fld(nint(h(2)),nint(h(3)),skew=h(4)*acos(-1.d0)/180.d0,rc=(/ h(5),h(6) /))
          else
            ok_p=dup_temp(nint(h(2)),nint(h(3)),skew=h(4)*acos(-1.d0)/180.d0,rc=(/ h(5),h(6) /))
          endif
          if (.not.ok_p) call output('!--- '//trim(s2)//'DUPLICATE failed ---')
        elseif (token(1)(k:k).eq.'T') then
          write(os,"(a,i3,1(x,e19.12),a)") '> '//trim(s2)//'TORIGIN',nint(h(2)),h(3),trim(s1); call output(os)
          if (k.lt.6) then
            ok_p=ref_fld(ifld=nint(h(2)),t_origin=h(3),add=ok_p)
          else
            ok_p=ref_temp(ifld=nint(h(2)),t_origin=h(3),add=ok_p)
          endif
          if (.not.ok_p) call output('!--- '//trim(s2)//'TORIGIN failed ---')
        else
          write(os,"(a)") '!--- '//trim(token(1))//' not implemented yet ---'; call output(os)
        endif
      elseif (token(1)(1:4).eq.'LLG_') then
        if (token(1)(1:7).eq.'LLG_FIX') then
          if (num_token.gt.1) then
            ok_p=(token(2)(1:2).eq.'.T'.or.token(2)(1:1).eq.'T')
            if (ok_p) then
              call output('> LLG_FIXED_DT TRUE')
            else
              call output('> LLG_FIXED_DT FALSE')
            endif
            call set_llg(layer_n,fixed_dt_p=ok_p)
          else
            call output('!--- LLG_FIXED_DT needs additional boolean paramemter: TRUE|FALSE ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_THE') then
          if (num_token.gt.1) then
            ok_p=(token(2)(1:2).eq.'.T'.or.token(2)(1:1).eq.'T')
            write(os,"(a,x,l,x,a)")'> LLG_THERMAL',ok_p,trim(layer_s); call output(os)
            if (layer_n.eq.0) then
              do k=1,n_lay
                call set_data(k,ok_p,thermal_p=.true.)
              enddo
            else
              call set_data(layer_n,ok_p,thermal_p=.true.)
            endif
          else
            call output('!--- LLG_THERMAL needs additional boolean paramemter: TRUE|FALSE ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_TYP') then
          i=nint(h(2))
          if (i.gt.-5.and.i.lt.3) then
            write(os,"(a,i2,x,a)") '> LLG_TYPE',i,trim(layer_s); call output(os)
            call set_llg(layer_n,llg_type=i)
          else
            call output('!--- unknown type for LLG_TYPE ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_FLO') then
          if (num_token.gt.1) then
            h(2)=min(1.d0,max(0.d0,h(2)))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_FLOORMS',h(2),trim(layer_s); call output(os)
            if (layer_n.eq.0) then
              do k=1,n_lay
                call set_data(k,h(2),floor_ms_p=.true.)
              enddo
            else
              call set_data(layer_n,h(2),floor_ms_p=.true.)
            endif
          else
            call output('!--- LLG_FLOORMS needs additional paramemter: 0<=floor_ms<=1 ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_DTF') then
          if (num_token.gt.1) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12)") '> LLG_DTFRAC',h(2); call output(os)
            call set_llg(layer_n,dt_frac=h(2))
          else
            call output('!--- LLG_DTFRAC needs additional paramemter: 0<=dt_frac ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_TNU') then
          if (num_token.gt.1) then
            i=max(1,nint(h(2)))
            write(os,"(a,x,i6)") '> LLG_TNUMAVE',i; call output(os)
            call set_llg(layer_n,thermal_num_ave=i)
          else
            call output('!--- LLG_TNUMAVE needs additional integer paramemter: 0<num ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_GAM') then
          if (num_token.gt.1) then
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_GAMMA',h(2),trim(layer_s); call output(os)
            if (layer_n.eq.0) then
              do k=1,n_lay
                call set_data(k,h(2),gamma_p=.true.)
              enddo
            else
              call set_data(layer_n,h(2),gamma_p=.true.)
            endif
          else
            call output('!--- LLG_GAMMA needs additional paramemter: gamma ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_TEM') then
          if (num_token.gt.1) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_TEMPERATURE',h(2),trim(layer_s); call output(os)
            if (layer_n.eq.0) then
              do k=1,n_lay
                call set_data(k,h(2),background_t_p=.true.)
              enddo
            else
              call set_data(layer_n,h(2),background_t_p=.true.)
            endif
          else
            call output('!--- LLG_TEMPERATURE needs additional paramemter: gamma ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_DTMA') then
          if (num_token.gt.1.and.h(2).gt.0.d0) then
            write(os,"(a,x,1p,e19.12)") '> LLG_DTMAX',h(2); call output(os)
            call set_llg(layer_n,dt_max=h(2))
          else
            call output('!--- LLG_DTMAX needs additional paramemter: 0<dt_max ---')
          endif
        elseif (token(1)(1:6).eq.'LLG_DT') then
          if (num_token.gt.1.and.h(2).gt.0.d0) then
            write(os,"(a,x,1p,e19.12)") '> LLG_DT',h(2); call output(os)
            call set_llg(layer_n,dt=h(2))
          else
            call output('!--- LLG_DT needs additional paramemter: 0<dt_max ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_ERR1') then
          if (num_token.gt.1.and.h(2).ge.0.d0) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_ERR1',h(2),trim(layer_s); call output(os)
            call set_llg(layer_n,conv_angle=h(2))
          else
            call output('!--- LLG_ERR1 needs additional paramemter: 0<=convergence_angle ---')
          endif
        elseif (token(1)(1:7).eq.'LLG_ERR') then
          if (num_token.gt.1.and.h(2).ge.0.d0) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_ERR',h(2),trim(layer_s); call output(os)
            call set_llg(layer_n,conv_torque=h(2))
          else
            call output('!--- LLG_ERR needs additional paramemter: 0<=convergence_torque/Ms ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_MAXD') then
          if (num_token.gt.1.and.h(2).gt.0.d0) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_MAXDM',h(2),trim(layer_s); call output(os)
            call set_llg(layer_n,max_dm=h(2))
          else
            call output('!--- LLG_MAXDM needs additional paramemter: 0<max_rotation ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_MIND') then
          if (num_token.gt.1.and.h(2).ge.0.d0) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_MINDM',h(2),trim(layer_s); call output(os)
            call set_llg(layer_n,min_dm=h(2))
          else
            call output('!--- LLG_MINDM needs additional paramemter: 0<=min_rotation ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_TMIN') then
          if (num_token.gt.1) then
            write(os,"(a,x,1p,e19.12)") '> LLG_TMIN',h(2); call output(os)
            call set_llg(layer_n,t_min=h(2))
          else
            call output('!--- LLG_TMIN needs additional paramemter: min_time ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_TMAX') then
          if (num_token.gt.1) then
            write(os,"(a,x,1p,e19.12)") '> LLG_TMAX',h(2); call output(os)
            call set_llg(layer_n,t_max=h(2))
          else
            call output('!--- LLG_TMAX needs additional paramemter: max_time ---')
          endif
        elseif (token(1)(1:8).eq.'LLG_DELT') then
          if (num_token.gt.1.and.h(2).gt.0.d0) then
            h(2)=max(0.d0,h(2))
            write(os,"(a,x,1p,e19.12)") '> LLG_DELT',h(2); call output(os)
            call set_llg(layer_n,del_t=h(2))
          else
            call output('!--- LLG_DELT needs additional paramemter: 0<elapsed time ---')
          endif
        elseif (token(1)(1:6).eq.'LLG_AL') then
          if (num_token.gt.1.and.h(2).ge.0.d0) then
            write(os,"(a,x,1p,e19.12,x,a)") '> LLG_ALPHA',h(2),trim(layer_s); call output(os)
            if (layer_n.eq.0) then
              do k=1,n_lay
                call set_data(k,h(2),alpha_p=.true.)
              enddo
            else
              call set_data(layer_n,h(2),alpha_p=.true.)
            endif
          else
            call output('!--- LLG_ALPHA needs additional paramemter: 0<=alpha ---')
          endif
!       elseif (token(1)(1:7).eq.'LLG_FLO') then
        else
          write(os,"(a)") '!--- '//trim(token(1))//' not implemented yet ---'; call output(os)
        endif
!     elseif (token(1)().eq.'') then
!     elseif (token(1)().eq.'') then
      elseif (token(1)(1:6).eq.'PRINT ') then
        if (num_token.eq.3.and.(token(3)(1:1).eq.'T'.or.token(3)(1:1).eq.'W'.or.token(3)(1:1).eq.'M')) then
          i=find_file(trim(token(2)))
          if (i.gt.0) then
            if (token(3)(1:1).eq.'T') then
              do j=1,n_lay; call plot_temp_scaling(j,0.d0,(/0.d0,0.d0/),i); enddo
            elseif (token(3)(1:1).eq.'W') then
              do j=1,n_lay; call plot_wf_scaling(j,0.d0,(/0.d0,0.d0/),i); enddo
            elseif (token(3)(1:1).eq.'M') then
              do j=1,n_lay; call plot_mz_scaling(j,0.d0,(/0.d0,0.d0/),i); enddo
            endif
          else
            call output('!--- no such file is open ---')
          endif
        else
          call output('!--- PRINT needs additional paramemters: <filename> T|W|M ---')
        endif
      elseif (token(1)(1:7).eq.'PLOT_MH') then
        ok_p=(token(2)(1:1).eq.'T'.or.token(2)(1:2).eq.'.T')
        if (num_token.eq.2) then
          write(os,"(a,x,l,x,a)") '> PLOT_MH',ok_p,trim(layer_s); call output(os)
!         call set_plot(layer_n,mhloop=ok_p)
        elseif (num_token.eq.5) then
          write(os,"(a,x,l,x,1p,3(e19.12,x),a)") '> PLOT_MH',ok_p,h(3:5),trim(layer_s);call output(os)
!         call set_plot(layer_n,mhloop=ok_p,mh=h(3:5))
          call set_plot(layer_n,mh=h(3:5))
        else
         call output('!--- PLOT_MH requires: TRUE|FALSE [Hx Hy Hz] ---')
        endif
      elseif (token(1)(1:5).eq.'PLOT ') then
        if (token(2)(1:4).eq.'HAMR') then
          call output('> PLOT HAMR'); call plot_hamr()
        elseif (token(2)(1:2).eq.'MH') then
          if (num_token.eq.3) then
            write(os,"(a,x,1p,e19.12)") '> PLOT MH ',h(3); call output(os)
            if (line_read.gt.line_target) call plot_mh(h(3),time)
          elseif (num_token.eq.5) then
            write(os,"(a,1p,3(x,e19.12))") '> PLOT MH ',h(3:5); call output(os)
            call set_plot(0,mh=h(3:5))
            if (line_read.gt.line_target) call plot_mh(sqrt(dot_product(h(3:5),h(3:5))),time)
          else
            call output('!--- PLOT MH requires either: |H| or Hx Hy Hz ---')
          endif
        elseif (token(2)(1:2).eq.'ON') then
          write(os,"(a)")'> PLOT ON '//trim(layer_s);call output(os); call set_plot(plotting=.true.,il=layer_n)
        elseif (token(2)(1:3).eq.'OFF') then
          write(os,"(a)")'> PLOT OFF '//trim(layer_s);call output(os); call set_plot(plotting=.false.,il=layer_n)
        elseif (token(2)(1:3).eq.'MAG') then
          call output('> PLOT MAG'); call plot_llg(time,done_p=.true.)
        elseif (token(2)(1:1).eq.'W') then
          x_p=.false.; y_p=.false.; z_p=.false.; ok_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            ok_p=(ok_p.or.token(i)(1:1).eq.'E');v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT WRITEFIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,wf=(.not.ok_p.and..not.x_p.and..not.y_p.and..not.z_p),wf_eff=ok_p,wf_x=x_p,wf_y=y_p,wf_z=z_p)
        elseif (token(2)(1:1).eq.'A') then
          x_p=.false.; y_p=.false.; z_p=.false.; ok_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            ok_p=(ok_p.or.token(i)(1:1).eq.'E');v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT APPLIED_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,ha=(v_p.or..not.ok_p.and..not.x_p.and..not.y_p.and..not.z_p),ha_eff=ok_p,ha_x=x_p,ha_y=y_p,ha_z=z_p)
        elseif (token(2)(1:1).eq.'B') then
          if (num_token.gt.11) then
            x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
            do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
              v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
            write(os,"(a,:,50(x,a))")'> PLOT BLOCK_FIELD ',(trim(token(i)),i=3,num_token),trim(layer_s);call output(os)
            if (num_token.ge.14) then
              call plot_field(il=layer_n,hb=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hb_x=x_p,hb_y=y_p,hb_z=z_p,&
                  hb_ds=(/h(3),h(4),h(5)/),hb_r=(/h(6),h(7),h(8)/),hb_m=(/h(9),h(10),h(11)/),hb_tilt=(/h(12),h(13),h(14)/))
            else
              call plot_field(il=layer_n,hb=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hb_x=x_p,hb_y=y_p,hb_z=z_p,&
                  hb_ds=(/h(3),h(4),h(5)/),hb_r=(/h(6),h(7),h(8)/),hb_m=(/h(9),h(10),h(11)/))
            endif
          else
            call output("!--- PLOT BLOCK_FIELD dxs dys dzs rx ry gap mx my mz  X|Y|Z|M ---") 
            call output("!--- PLOT BLOCK_FIELD thick width length rx ry gap mx my mz acos(z) X|Y|Z|M ---") 
          endif
        elseif (token(2)(1:1).eq.'D') then
          x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT DEMAG_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,hd=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hd_x=x_p,hd_y=y_p,hd_z=z_p)
        elseif (token(2)(1:1).eq.'E') then
          x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT EXCHANGE_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,hx=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hx_x=x_p,hx_y=y_p,hx_z=z_p)
        elseif (token(2)(1:1).eq.'V') then
          x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT VERTICAL_EXCHANGE_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,hv=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hv_x=x_p,hv_y=y_p,hv_z=z_p)
        elseif (token(2)(1:2).eq.'HK') then
          x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT HK_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,hk=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),hk_x=x_p,hk_y=y_p,hk_z=z_p)
        elseif (token(2)(1:2).eq.'TO') then
          x_p=.false.; y_p=.false.; z_p=.false.; v_p=.false.
          do i=num_token,3,-1; x_p=(x_p.or.token(i)(1:1).eq.'X');y_p=(y_p.or.token(i)(1:1).eq.'Y');z_p=(z_p.or.token(i)(1:1).eq.'Z');
            v_p=(v_p.or.token(i)(1:1).eq.'M');enddo
          write(os,"(a,:,50(x,a))")'> PLOT TOTAL_FIELD ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,htot=(v_p.or..not.x_p.and..not.y_p.and..not.z_p),htot_x=x_p,htot_y=y_p,htot_z=z_p)
        elseif (token(2)(1:1).eq.'T') then
          write(os,"(a,:,50(x,a))")'> PLOT TEMPERATURE ',(token(i)(1:1),i=3,num_token),trim(layer_s);call output(os)
          call plot_field(il=layer_n,temp=.true.,tempf=.true.)
        elseif (token(2)(1:1).eq.'F') then
          if (num_token.eq.4.or.num_token.eq.6) then
            write(os,"(a,:,50(x,a))")'> PLOT FIELD ',(trim(token(i)),i=3,num_token);call output(os)
            if (num_token.eq.4) call plot_off_field(getnum(token(3)),getnum(token(4)))
            if (num_token.eq.6) call plot_off_field(getnum(token(3)),getnum(token(4)),getnum(token(5)),getnum(token(6)))
          endif
        else
          write(os,"(a)") '!--- '//trim(token(1))//' '//trim(token(2))//' not implemented yet ---'; call output(os)
        endif
      else
        write(os,"(a)") '!--- unknown/not implemented command: '//trim(token(1)); call output(os)
      endif
    endif
!           read *

    i=99
    do while (i.ne.0) 
      cliionum = ioc(ioc_num)
      if (key_p) then
        print "(' command< ',$)"
        read(*,"(a)",iostat=i) cl
      else
        read(cliionum,"(a)",iostat=i) cl
      endif
      if (i.ne.0) then
        if (cliionum.gt.0) then; close(cliionum); call delete_file(cliionum); endif
        ioc_num=ioc_num-1
        if (ioc_num.lt.1) exit
        cliionum=ioc(ioc_num); key_p=(cliionum.lt.1)
      endif
    enddo
    if (i.ne.0) exit
  enddo

  call close_all_files()
  ok_p=my_checkpoint(.false.,clean_p=.true.)
  call output('done executing cli')

 contains

  subroutine get_token(clin,token,num_token,layer_n,layer_s)
    character(*)::clin
    character(MAXCHAR)::cl,token(MAXTKN),layer_s
    integer,intent(out)::layer_n,num_token

    integer::i,j,k
    logical::noupper_p

    cl=' '//clin
    j=index(cl,'!'); if (j.gt.0) cl=cl(1:j-1)
    cl=adjustl(cl); num_token=0; token(:)=' '


    noupper_p=.false.
    do while (len_trim(cl).gt.0.and.num_token.lt.MAXTKN)
      if (cl(1:1).eq.'"') then
        j=index(cl(2:),'"'); if (j.eq.0) j=len_trim(cl)
        num_token=num_token+1; token(num_token)=cl(2:j); cl=cl(j+2:)
      elseif (cl(1:1).eq."'") then
        j=index(cl(2:),"'"); if (j.eq.0) j=len_trim(cl)+1
        num_token=num_token+1; token(num_token)=cl(2:j); cl=cl(j+2:)
      else
        j=index(cl(2:),' ')
        k=index(cl(2:),'"'); if (j.eq.0.or.(k.gt.0.and.k.lt.j)) j=k
        k=index(cl(2:),"'"); if (j.eq.0.or.(k.gt.0.and.k.lt.j)) j=k
        if (j.eq.0) j=len_trim(cl)
        num_token=num_token+1; 
        if (noupper_p) then
          token(num_token)=cl(1:j); cl=cl(j+1:)
        else
          token(num_token)=toupper(cl(1:j)); cl=cl(j+1:)
        endif
      endif
      cl=adjustl(cl)
      noupper_p=(num_token.eq.1.and.(token(1)(1:4).eq.'OPEN'.or.token(1)(1:5).eq.'CLOSE'.or.token(1)(1:4).eq.'READ'.or.token(1)(1:4).eq.'SAVE'.or.token(1)(1:4).eq.'REST'.or. &
                                     token(1)(1:5).eq.'WRITE'.or.token(1)(1:7).eq.'INPUT_F'.or.token(1)(1:5).eq.'PRINT'))
    enddo
   !look for token either 'LN' or 'ALL' which gives layer information, if found eliminate from token list
    layer_n=0; layer_s='ALL'  !assume command is for all layers
    do i=num_token,2,-1; j=0
      if (token(i)(1:1).eq.'L') then
        layer_n=nint(getnum(trim(token(i)(2:)))); j=i+1
        layer_s=trim(token(i))
      elseif (token(i)(1:3).eq.'ALL'.and.token(1)(1:5).ne.'PLOT ') then
        layer_n=0; j=i+1
      endif
      if (j.ne.0) then
        do j=i+1,num_token; token(j-1)=token(j); enddo
        token(num_token)=' ';num_token=num_token-1
      endif
    enddo
  end subroutine
  function getnum(str) result (x)
    character(*)::str
    character(len(str))::s1
    real(8)::x
    integer::i
    s1=trim(adjustl(str))
    read(s1,*,iostat=i) x
    if (i.ne.0) x=0.d0
  end function
  function toupper(str) result(s)
    character(*)::str
    character(len(str))::s
    integer::i
    s=''
    do i=1,len_trim(str)
      if (str(i:i).ge.'a'.and.str(i:i).le.'z') then
        s(i:i)=char(ichar(str(i:i))+ichar('A')-ichar('a'))
      else
        s(i:i)=str(i:i)
      endif
    enddo
  end function
 end subroutine cli
end module cli_m
