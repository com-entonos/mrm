module fftw3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
end module fftw3
! this calls fftw for real fft
!  supports single and double precision
! on windwos, you can compile fftw (good luck!) or use intel's MKL libraries which have it

module fftw_m
  use, intrinsic :: iso_c_binding
  use fftw3
  implicit none
  private
  public :: fft,ffttr,ffttc,fft_quick_plan

! this is how you control the precision of fft. default is single precision (float)
#ifdef FFTD
 integer,parameter::ffttr=C_DOUBLE   !this is double precision fft
 integer,parameter::ffttc=C_DOUBLE_COMPLEX   !this is double precision fft
#else
 integer,parameter::ffttr=C_FLOAT    !this is single precision fft
 integer,parameter::ffttc=C_FLOAT_COMPLEX    !this is single precision fft
#endif

  integer(C_INT),save::planner=FFTW_EXHAUSTIVE !FFTW_ESTIMATE

contains
  subroutine fft_quick_plan()
    use fftw3
    if (planner.eq.FFTW_ESTIMATE) then
      planner=FFTW_EXHAUSTIVE
    else
      planner=FFTW_ESTIMATE
    endif
  end subroutine

  subroutine fft(r,c,nx,ny,inv,p)
   use fftw3
   integer,intent(in)::nx,ny                    !size of transform
   real(ffttr),intent(inout)::r(nx,ny)              !real data
   complex(ffttc),intent(inout)::c(nx/2+1,ny)       !complex data
   logical,intent(in)::inv                      !inverse or forward transformation
   type(C_PTR)::p                               !plans for forward and backward transform

   complex(ffttc),allocatable,dimension(:,:)::d
   real(ffttr),allocatable,dimension(:,:)::s

   if (c_associated(p).and.inv) then
#ifdef FFTD
     call fftw_execute_dft_c2r(p,c,r); return
#else
     call fftwf_execute_dft_c2r(p,c,r); return
#endif
   elseif (c_associated(p).and..not.inv) then
#ifdef FFTD
     call fftw_execute_dft_r2c(p,r,c); return
#else
     call fftwf_execute_dft_r2c(p,r,c); return
#endif
   endif

   if (inv) then
     allocate(d(nx/2+1,ny)); d=c;
     do while (.not.c_associated(p))
#ifdef FFTD
       p = fftw_plan_dft_c2r_2d(ny,nx,c,r,ior(planner,FFTW_DESTROY_INPUT))
#else
       p = fftwf_plan_dft_c2r_2d(ny,nx,c,r,ior(planner,FFTW_DESTROY_INPUT))
#endif
     enddo; c=d
     deallocate(d)
#ifdef FFTD
     call fftw_execute_dft_c2r(p,c,r)
#else
     call fftwf_execute_dft_c2r(p,c,r)
#endif
   else
     allocate(s(nx,ny)); s=r;
     do while (.not.c_associated(p))
#ifdef FFTD
       p = fftw_plan_dft_r2c_2d(ny,nx,r,c,ior(planner,FFTW_DESTROY_INPUT))
#else
       p = fftwf_plan_dft_r2c_2d(ny,nx,r,c,ior(planner,FFTW_DESTROY_INPUT))
#endif
     enddo; r=s
     deallocate(s)
#ifdef FFTD
     call fftw_execute_dft_r2c(p,r,c)
#else
     call fftwf_execute_dft_r2c(p,r,c)
#endif
   endif
  end subroutine

end module fftw_m
