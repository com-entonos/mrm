
module signal_m
 implicit none
 private
 public :: set_signal, catch_sigterm

contains

 function set_signal() result (iret)

#ifdef __INTEL_COMPILER
!for intel fortran
  use ifport

  external signal_m_mp_catch_sigterm
  integer::signal_m_mp_catch_sigterm

  integer::iret
  iret = signal(SIGTERM, signal_m_mp_catch_sigterm, -1)
  iret = signal(SIGINT, signal_m_mp_catch_sigterm, -1)
  iret = signal(SIGABRT, signal_m_mp_catch_sigterm, -1)

#else
!for gfortran
  integer,parameter::SIGTERM=15,SIGINT=2
  integer::iret
  iret = signal(SIGTERM, catch_sigterm)
  iret = signal( SIGINT, catch_sigterm)
#endif

 end function


 function catch_sigterm(sig_num) result (iret)
  use io_m, only: output
  use cli_m, only: checkpoint
  integer::sig_num,iret
  logical::ok_p
  character(200)::os
  write(os,"('we caught a signal, sig_num=',i2)") sig_num; call output(trim(os),.true.)
  ok_p=checkpoint(read_p=.false.,quit_p=.true.)
  iret=1
  stop
 end function

end module signal_m


