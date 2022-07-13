
program mrm
 use io_m, only: output, flush_output, add_file, delete_file
 use signal_m, only: set_signal
 use input_m, only: get_input
 use cli_m, only: cli
 implicit none

 integer::i,j
 real(4)::start_time,end_time

! character(200)::in_file='som.inp',cli_file='som.cli'
 character(200)::in_file='',cli_file=''
 character(200)::os

 call cpu_time(start_time)

 call clean_run_files()
 i=set_signal()

 call output(' ')
 call output(' mrm- micromagnetic recording model',.true.)
 call output(' ')

 i=command_argument_count()
 if (i.gt.0) call get_command_argument(1,in_file); in_file=adjustl(in_file)
 if (i.gt.1) call get_command_argument(2,cli_file); cli_file=adjustl(cli_file)

 if (get_input(trim(in_file),cli_file)) call cli(trim(cli_file))

 call cpu_time(end_time)
 i=int((end_time-start_time)/3600); j=int(((end_time-start_time)-i*3600)/60)
 write(os,"('simulation took:',i3,' hours, ',i2,' minutes, ',f5.2,' seconds')") i,j,(end_time-start_time)-60*(j+60*i)
 call output(trim(os),.true.)

 call final_run_files()

 call output('successful completion of mrm')
 call flush_output()

contains

   subroutine clean_run_files()  !delete RUNTIME.LOG and WARNING.LOG; create FAILURE.LOG
     integer::i

     i=add_file('RUNTIME.LOG')
     open(i,file='RUNTIME.LOG',status='unknown'); close(i,status='delete')
     open(i,file='WARNING.LOG',status='unknown'); close(i,status='delete')
     open(i,file='FAILURE.LOG',status='unknown'); write(i,"('mrm failed (blame gjp, perhaps)')"); close(i)
     call delete_file(i)
   end subroutine

   subroutine final_run_files()  !create RUNTIME.LOG w/ elapsed cpu_time; delete FAILURE.LOG
     integer::i

     i=add_file('RUNTIME.LOG')
     open(i,file='RUNTIME.LOG',status='unknown')
     write(i,*) end_time-start_time
     close(i)
     open(i,file='FAILURE.LOG',status='unknown')
     close(i,status='delete')
   end subroutine
end program mrm
