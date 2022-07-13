module demag_m
 use fftw3
 use fftw_m, only: ffttr,ffttc
 implicit none
 private

 public::start_demag,demag_window,get_demag,init_demag,get_demag_window

 integer,save::nx=-1,ny=-1,n_lay=-1
 logical,save::calc_demag_p=.false.,tensor_calculated_p=.false.
 integer,save::num_demags=0

 real(ffttr),pointer,save::rx2d(:,:,:)=>null(),ry2d(:,:,:)=>null(),rz2d(:,:,:)=>null()
 complex(ffttc),pointer,save::x2d(:,:,:),y2d(:,:,:),z2d(:,:,:)
 real(ffttr),pointer,save::xi2d(:,:,:),yi2d(:,:,:),zi2d(:,:,:)
 complex(ffttc),pointer,save::n2d(:,:,:,:,:),xo2d(:,:,:),yo2d(:,:,:),zo2d(:,:,:)
 type(C_PTR),pointer,save::xp(:),yp(:),zp(:),xip(:),yip(:),zip(:)

 type::demagwindow_s
   logical::selfdemag_p
 end type demagwindow_s
 type(demagwindow_s),save::demagwindow
 logical,save::demagwin_init=.false.


contains

  function demag_count() result (i)
   integer::i
   i=num_demags
  end function

  subroutine get_demag(hx,hy,hz)
   use fftw_m, only: fft
   use data_m, only: get_mag
   real(ffttr),pointer::hx(:,:,:),hy(:,:,:),hz(:,:,:)
   integer::k1,k2

   hx=>rx2d; hy=>ry2d; hz=>rz2d
   if (.not.calc_demag_p.or..not.associated(rx2d)) return

   num_demags=num_demags+1

   if (get_mag(xi2d,yi2d,zi2d)) then
     x2d=0.0; y2d=0.0; z2d=0.0
#ifdef OMP
!$omp parallel private(k1,k2)
!$omp sections
!$omp section
     do k1=1,n_lay
       call fft(xi2d(:,:,k1),xo2d(:,:,k1),nx,ny,.false.,xip(k1))
!       do k2=1,n_lay
!         x2d(:,:,k2)=x2d(:,:,k2)+n2d(:,:,1,k2,k1)*xo2d(:,:,k1)
!         y2d(:,:,k2)=y2d(:,:,k2)+n2d(:,:,4,k2,k1)*xo2d(:,:,k1)
!         z2d(:,:,k2)=z2d(:,:,k2)+n2d(:,:,7,k2,k1)*xo2d(:,:,k1)
!       enddo
     enddo
!$omp section
     do k1=1,n_lay
       call fft(yi2d(:,:,k1),yo2d(:,:,k1),nx,ny,.false.,yip(k1))
!       do k2=1,n_lay
!         x2d(:,:,k2)=x2d(:,:,k2)+n2d(:,:,2,k2,k1)*yo2d(:,:,k1)
!         y2d(:,:,k2)=y2d(:,:,k2)+n2d(:,:,5,k2,k1)*yo2d(:,:,k1)
!         z2d(:,:,k2)=z2d(:,:,k2)+n2d(:,:,8,k2,k1)*yo2d(:,:,k1)
!       enddo
     enddo
!$omp section
     do k1=1,n_lay
       call fft(zi2d(:,:,k1),zo2d(:,:,k1),nx,ny,.false.,zip(k1))
!       do k2=1,n_lay
!         x2d(:,:,k2)=x2d(:,:,k2)+n2d(:,:,3,k2,k1)*zo2d(:,:,k1)
!         y2d(:,:,k2)=y2d(:,:,k2)+n2d(:,:,6,k2,k1)*zo2d(:,:,k1)
!         z2d(:,:,k2)=z2d(:,:,k2)+n2d(:,:,9,k2,k1)*zo2d(:,:,k1)
!       enddo
     enddo
!$omp end sections
!$omp sections
!$omp section
     do k1=1,n_lay
       do k2=1,n_lay
         x2d(:,:,k2)=x2d(:,:,k2)+n2d(:,:,1,k2,k1)*xo2d(:,:,k1)+n2d(:,:,2,k2,k1)*yo2d(:,:,k1)+n2d(:,:,3,k2,k1)*zo2d(:,:,k1)
       enddo
     enddo
!$omp section
     do k1=1,n_lay
       do k2=1,n_lay
         y2d(:,:,k2)=y2d(:,:,k2)+n2d(:,:,4,k2,k1)*xo2d(:,:,k1)+n2d(:,:,5,k2,k1)*yo2d(:,:,k1)+n2d(:,:,6,k2,k1)*zo2d(:,:,k1)
       enddo
     enddo
!$omp section
     do k1=1,n_lay
       do k2=1,n_lay
         z2d(:,:,k2)=z2d(:,:,k2)+n2d(:,:,7,k2,k1)*xo2d(:,:,k1)+n2d(:,:,8,k2,k1)*yo2d(:,:,k1)+n2d(:,:,9,k2,k1)*zo2d(:,:,k1)
       enddo
     enddo
!$omp end sections
!$omp sections
!$omp section
     do k1=1,n_lay
       call fft(rx2d(:,:,k1),x2d(:,:,k1),nx,ny,.true.,xp(k1))
     enddo
!$omp section
     do k1=1,n_lay
       call fft(ry2d(:,:,k1),y2d(:,:,k1),nx,ny,.true.,yp(k1))
     enddo
!$omp section
     do k1=1,n_lay
       call fft(rz2d(:,:,k1),z2d(:,:,k1),nx,ny,.true.,zp(k1))
     enddo
!$omp end sections
!$omp end parallel
#else
     do k1=1,n_lay
       call fft(xi2d(:,:,k1),xo2d(:,:,k1),nx,ny,.false.,xip(k1))
       call fft(yi2d(:,:,k1),yo2d(:,:,k1),nx,ny,.false.,yip(k1))
       call fft(zi2d(:,:,k1),zo2d(:,:,k1),nx,ny,.false.,zip(k1))
       do k2=1,n_lay
         x2d(:,:,k2)=x2d(:,:,k2)+n2d(:,:,1,k2,k1)*xo2d(:,:,k1)+n2d(:,:,2,k2,k1)*yo2d(:,:,k1)+n2d(:,:,3,k2,k1)*zo2d(:,:,k1)
         y2d(:,:,k2)=y2d(:,:,k2)+n2d(:,:,4,k2,k1)*xo2d(:,:,k1)+n2d(:,:,5,k2,k1)*yo2d(:,:,k1)+n2d(:,:,6,k2,k1)*zo2d(:,:,k1)
         z2d(:,:,k2)=z2d(:,:,k2)+n2d(:,:,7,k2,k1)*xo2d(:,:,k1)+n2d(:,:,8,k2,k1)*yo2d(:,:,k1)+n2d(:,:,9,k2,k1)*zo2d(:,:,k1)
       enddo
     enddo
     do k1=1,n_lay
       call fft(rx2d(:,:,k1),x2d(:,:,k1),nx,ny,.true.,xp(k1))
       call fft(ry2d(:,:,k1),y2d(:,:,k1),nx,ny,.true.,yp(k1))
       call fft(rz2d(:,:,k1),z2d(:,:,k1),nx,ny,.true.,zp(k1))
     enddo
#endif
!open(66,file='hdemag.dat',status='unknown')
!do k1=1,n_lay;do j=1,ny; do i=1,nx
!  write(66,"(3(i3,x),1p,3(e12.5,x))") i,j,k1,rx2d(i,j,k1),ry2d(i,j,k1),rz2d(i,j,k1)
!enddo; enddo; enddo
!close(66)
   endif
  end subroutine

  subroutine demag_window(demag_p,v,offset,i0,j0,nx,ny)
   use fftw_m
   use window_m, only: set_window,get_window
   logical,intent(in),optional::demag_p
   integer,intent(in),optional::i0,j0,nx,ny
   real(8),intent(in),optional::v(2),offset(2)
   real(8)::vv(2),off(2)
   integer::n2(2)

   if (.not.demagwin_init) then
     demagwin_init=.true.
     demagwindow%selfdemag_p=.false.
     if (present(offset)) then !is this necessary? FIXME
       call set_window(offset=offset)
     else
       call set_window()
     endif
   endif
   if (present(demag_p)) demagwindow%selfdemag_p=demag_p
   if (present(v)) call set_window(v=v)
   if (present(offset)) then
     call get_window(v=vv,offset=off)
     if (vv(1).eq.0.d0) off(1)=offset(1)
     if (vv(2).eq.0.d0) off(2)=offset(2)
     call set_window(offset=off)
   endif
   if (present(i0)) call set_window(i0=i0)
   if (present(j0)) call set_window(j0=j0)
   if (present(nx).or.present(ny)) then
     call get_window(n=n2)
     if (present(nx)) n2(1)=nx
     if (present(ny)) n2(2)=ny
     call set_window(nx=n2(1),ny=n2(2))
   endif
  end subroutine

  function get_demag_window(dx,dy) result (r)
    use window_m, only: get_window
    real(8),intent(in)::dx,dy
    real(8)::r(5,2)
    integer::i0(2)
!   i0=nint(get_demag(i0_p=.true.))
    call get_window(i0=i0)
    r(:,1)=(/ i0(1), i0(1)+nx,i0(1)+nx,i0(1)   ,i0(1) /)*dx
    r(:,2)=(/ i0(2), i0(2),   i0(2)+ny,i0(2)+ny,i0(2) /)*dy
  end function

  function start_demag(nxi,nyi,n_layi,do_p,demag_p) result (ok_p)
   use io_m, only: output
   integer,intent(in)::nxi,nyi,n_layi
   logical,intent(in)::demag_p,do_p
   logical::ok_p
   nx=nxi; ny=nyi; n_lay=n_layi
   ok_p=.true.
   calc_demag_p=demag_p; num_demags=0; call demag_window()
   if (.not.demag_p) then
     call output('  *** magnetostatic field will not be computed! ***')
!    call demag_window(i0=-2*nx,j0=-2*ny)
     return
   endif
   if (do_p) ok_p=init_demag()
  end function
  
  function init_demag(self_only_p) result (ok_p)
   use io_m, only: output
   use tensor_m, only: calculate_tensor
   use fftw_m, only: fft
   logical,intent(in),optional::self_only_p
   logical::ok_p
   character(200)::os

   integer::k1,k2,i
   real(4)::start_time,end_time
   real(8),allocatable::rn2d(:,:,:,:,:)
   type(C_PTR)::p

   ok_p=.true.
   if (tensor_calculated_p.or..not.calc_demag_p) return

   call cpu_time(start_time)

   allocate(rn2d(nx,ny,9,n_lay,n_lay))
   tensor_calculated_p=.true.
   call calculate_tensor(rn2d)
   if (present(self_only_p).and.self_only_p) return

   ok_p=.false.

   allocate(xp(n_lay),yp(n_lay),zp(n_lay),xip(n_lay),yip(n_lay),zip(n_lay))
   xp = C_NULL_PTR; yp = C_NULL_PTR; zp = C_NULL_PTR
   xip= C_NULL_PTR; yip= C_NULL_PTR; zip= C_NULL_PTR
#ifdef FFTD
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,xi2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,xo2d,[nx/2+1,ny,n_lay])
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,yi2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,yo2d,[nx/2+1,ny,n_lay])
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,zi2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,zo2d,[nx/2+1,ny,n_lay])
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,rx2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, x2d,[nx/2+1,ny,n_lay])
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,ry2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, y2d,[nx/2+1,ny,n_lay])
     p=fftw_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,rz2d,[nx    ,ny,n_lay])
     p=fftw_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, z2d,[nx/2+1,ny,n_lay])
#else
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,xi2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,xo2d,[nx/2+1,ny,n_lay])
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,yi2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,yo2d,[nx/2+1,ny,n_lay])
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,zi2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p,zo2d,[nx/2+1,ny,n_lay])
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,rx2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, x2d,[nx/2+1,ny,n_lay])
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,ry2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, y2d,[nx/2+1,ny,n_lay])
     p=fftwf_alloc_real(   int( nx     *ny*n_lay,C_SIZE_T)); call c_f_pointer(p,rz2d,[nx    ,ny,n_lay])
     p=fftwf_alloc_complex(int((nx/2+1)*ny*n_lay,C_SIZE_T)); call c_f_pointer(p, z2d,[nx/2+1,ny,n_lay])
#endif
   allocate(n2d(nx/2+1,ny,9,n_lay,n_lay))
! fftw plans are not thread safe so we can not use openMP here. boo.
   call output('  transforming magnetostatic tensor')
   print "('    ',$)"
   do k1=1,n_lay
     do k2=1,n_lay
       do i=1,9,3
         xi2d(:,:,k1)=real(rn2d(:,:,i  ,k2,k1),ffttr); call fft(xi2d(:,:,k1),xo2d(:,:,k1),nx,ny,.false.,xip(k1)); n2d(:,:,i  ,k2,k1)=xo2d(:,:,k1); print "('.',$)"
         yi2d(:,:,k1)=real(rn2d(:,:,i+1,k2,k1),ffttr); call fft(yi2d(:,:,k1),yo2d(:,:,k1),nx,ny,.false.,yip(k1)); n2d(:,:,i+1,k2,k1)=yo2d(:,:,k1); print "('.',$)"
         zi2d(:,:,k1)=real(rn2d(:,:,i+2,k2,k1),ffttr); call fft(zi2d(:,:,k1),zo2d(:,:,k1),nx,ny,.false.,zip(k1)); n2d(:,:,i+2,k2,k1)=zo2d(:,:,k1); print "('.',$)"
       enddo
     enddo
   enddo
   call output('   -done')
   call output('  setting up fft plans for inverse transform')
   print "('    ',$)"
   do k1=1,n_lay
       x2d=cmplx(0.d0,0.d0,ffttc); y2d=cmplx(0.d0,0.d0,ffttc); z2d=cmplx(0.d0,0.d0,ffttc)
       call fft(rx2d(:,:,k1),x2d(:,:,k1),nx,ny,.true.,xp(k1)); print "('.',$)"
       call fft(ry2d(:,:,k1),y2d(:,:,k1),nx,ny,.true.,yp(k1)); print "('.',$)"
       call fft(rz2d(:,:,k1),z2d(:,:,k1),nx,ny,.true.,zp(k1)); print "('.',$)"
   enddo
   call output('   -done')

   deallocate(rn2d)
   call cpu_time(end_time)
   i=int((end_time-start_time)/3600); k1=int(((end_time-start_time)-i*3600)/60)
   write(os,"('  tensor took:',i3,' hours, ',i2,' minutes, ',f5.2,' seconds')") i,k1,(end_time-start_time)-60*(k1+60*i)
   call output(trim(os),.true.)
   ok_p=.true.; tensor_calculated_p=.true.

  end function


end module demag_m
