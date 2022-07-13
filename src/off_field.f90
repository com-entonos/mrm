
module off_field_m
 use fftw3
 use fftw_m, only: ffttr,ffttc

 implicit none
 private
 public::plot_off_field,set_off_field_tensor

 logical,save::periodic_x=.true.,periodic_y=.true.,sul_image=.false.,head_image=.false.
 character(200),save::filename='mrm_N.dat'
 integer,save::num_images=2
 real(8),save::sul_mu=1.d0,fly_height=8.d-7,space_to_sul=40.d-7,head_mu=1.d0


 integer,save::nx=-1,ny=-1,n_lay=-1

 real(ffttr),pointer,save::rx2d(:,:),ry2d(:,:),rz2d(:,:)
 complex(ffttc),pointer,save::x2d(:,:),y2d(:,:),z2d(:,:)
 real(ffttr),pointer,save::xi2d(:,:),yi2d(:,:),zi2d(:,:)
 complex(ffttc),pointer,save::n2d(:,:,:,:),xo2d(:,:),yo2d(:,:),zo2d(:,:)
 type(C_PTR),save::xp=C_NULL_PTR,yp=C_NULL_PTR,zp=C_NULL_PTR,xip=C_NULL_PTR,yip=C_NULL_PTR,zip=C_NULL_PTR


contains

  subroutine plot_off_field(zoff,dzobs,mu_head,fly_height)
   use io_m, only: output
   use fftw_m, only: fft
   use data_m, only: get_size
   use plot_m, only: get_plot, plot_contour
   real(8),intent(in)::zoff,dzobs
   real(8),intent(in),optional::mu_head,fly_height
   logical::ok_p
   real(ffttr),dimension(:,:),pointer::hx=>null(),hy,hz
   real(8)::z,dz
   character(200)::os
   if (get_plot()) then
      call get_size(nx,ny,n_lay,z,dz)
      z=zoff*1.d-7; dz=dzobs*1.d-7
      if (init_demag(z,dz,mu_head,fly_height)) then
         call get_off_field(hx,hy,hz)
         write(os,"('magnetostatic Hx at ',f8.5,' dz=',f8.5)") zoff,dzobs; call plot_contour(trim(os),dble(hx),.false.,1)
         write(os,"('magnetostatic Hy at ',f8.5,' dz=',f8.5)") zoff,dzobs; call plot_contour(trim(os),dble(hy),.false.,1)
         write(os,"('magnetostatic Hz at ',f8.5,' dz=',f8.5)") zoff,dzobs; call plot_contour(trim(os),dble(hz),.false.,1)
      else
         call output("error generating tensor")
      endif
   else
      call output("plotting needs to be turned on")
   endif
  end subroutine

  subroutine get_off_field(hx,hy,hz)
   use fftw_m, only: fft
   use data_m, only: get_mag
   real(ffttr),pointer::hx(:,:),hy(:,:),hz(:,:)
   integer::k1,k2
   real(ffttr),allocatable::xi(:,:,:),yi(:,:,:),zi(:,:,:)

   hx=>rx2d; hy=>ry2d; hz=>rz2d

   allocate(xi(nx,ny,n_lay),yi(nx,ny,n_lay),zi(nx,ny,n_lay))

   if (get_mag(xi,yi,zi)) then
     x2d=0.0_ffttr; y2d=0.0_ffttr; z2d=0.0_ffttr
!$omp parallel private(k1,k2)
     do k1=1,n_lay
       xi2d=xi(:,:,k1); yi2d=yi(:,:,k1); zi2d=zi(:,:,k1)
!$omp sections
!$omp section
       call fft(xi2d(:,:),xo2d(:,:),nx,ny,.false.,xip)
!$omp section
       call fft(yi2d(:,:),yo2d(:,:),nx,ny,.false.,yip)
!$omp section
       call fft(zi2d(:,:),zo2d(:,:),nx,ny,.false.,zip)
!$omp end sections
#ifdef OMP
!$omp sections
!$omp section
         x2d(:,:)=x2d(:,:)+n2d(:,:,1,k1)*xo2d(:,:)+n2d(:,:,2,k1)*yo2d(:,:)+n2d(:,:,3,k1)*zo2d(:,:)
!$omp section
         y2d(:,:)=y2d(:,:)+n2d(:,:,4,k1)*xo2d(:,:)+n2d(:,:,5,k1)*yo2d(:,:)+n2d(:,:,6,k1)*zo2d(:,:)
!$omp section
         z2d(:,:)=z2d(:,:)+n2d(:,:,7,k1)*xo2d(:,:)+n2d(:,:,8,k1)*yo2d(:,:)+n2d(:,:,9,k1)*zo2d(:,:)
!$omp end sections
#else
         x2d(:,:)=x2d(:,:)+n2d(:,:,1,k1)*xo2d(:,:)+n2d(:,:,2,k1)*yo2d(:,:)+n2d(:,:,3,k1)*zo2d(:,:)
         y2d(:,:)=y2d(:,:)+n2d(:,:,4,k1)*xo2d(:,:)+n2d(:,:,5,k1)*yo2d(:,:)+n2d(:,:,6,k1)*zo2d(:,:)
         z2d(:,:)=z2d(:,:)+n2d(:,:,7,k1)*xo2d(:,:)+n2d(:,:,8,k1)*yo2d(:,:)+n2d(:,:,9,k1)*zo2d(:,:)
#endif
     enddo
     do k1=1,n_lay
!$omp sections
!$omp section
       call fft(rx2d(:,:),x2d(:,:),nx,ny,.true.,xp)
!$omp section
       call fft(ry2d(:,:),y2d(:,:),nx,ny,.true.,yp)
!$omp section
       call fft(rz2d(:,:),z2d(:,:),nx,ny,.true.,zp)
!$omp end sections
     enddo
!$omp end parallel
!open(66,file='hdemag.dat',status='unknown')
!do k1=1,n_lay;do j=1,ny; do i=1,nx
!  write(66,"(3(i3,x),1p,3(e12.5,x))") i,j,k1,rx2d(i,j,k1),ry2d(i,j,k1),rz2d(i,j,k1)
!enddo; enddo; enddo
!close(66)
   endif
   deallocate(xi,yi,zi)
  end subroutine

  function init_demag(zoff,dzobs,mu_head,fly_height) result (ok_p)
   use io_m, only: output
   use fftw_m, only: fft
   real(8),intent(in),optional::mu_head,fly_height
   real(8),intent(in)::zoff,dzobs
   logical::ok_p
   character(200)::os

   integer::k1,k2,i
   real(4)::start_time,end_time
   real(8),allocatable::rn2d(:,:,:,:)
   type(C_PTR)::p

   ok_p=.false.

   call cpu_time(start_time)

   if (.not.c_associated(xp)) then
#ifdef FFTD
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,xi2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,xo2d,[nx/2+1,ny      ])
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,yi2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,yo2d,[nx/2+1,ny      ])
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,zi2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,zo2d,[nx/2+1,ny      ])
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,rx2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, x2d,[nx/2+1,ny      ])
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,ry2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, y2d,[nx/2+1,ny      ])
       p=fftw_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,rz2d,[nx    ,ny      ])
       p=fftw_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, z2d,[nx/2+1,ny      ])
#else
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,xi2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,xo2d,[nx/2+1,ny      ])
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,yi2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,yo2d,[nx/2+1,ny      ])
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,zi2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p,zo2d,[nx/2+1,ny      ])
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,rx2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, x2d,[nx/2+1,ny      ])
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,ry2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, y2d,[nx/2+1,ny      ])
       p=fftwf_alloc_real(   int( nx     *ny      ,C_SIZE_T)); call c_f_pointer(p,rz2d,[nx    ,ny      ])
       p=fftwf_alloc_complex(int((nx/2+1)*ny      ,C_SIZE_T)); call c_f_pointer(p, z2d,[nx/2+1,ny      ])
#endif
     call output('  setting up fft plans for inverse transform')
     print "('    ',$)"
       x2d=cmplx(0.,0.,ffttc); y2d=cmplx(0.,0.,ffttc); z2d=cmplx(0.,0.,ffttc)
       call fft(rx2d(:,:),x2d(:,:),nx,ny,.true.,xp); print "('.',$)"
       call fft(ry2d(:,:),y2d(:,:),nx,ny,.true.,yp); print "('.',$)"
       call fft(rz2d(:,:),z2d(:,:),nx,ny,.true.,zp); print "('.',$)"
     call output('   -done')
   endif
   allocate(n2d(nx/2+1,ny,9,n_lay))
   allocate(rn2d(nx,ny,9,n_lay))

   call calculate_off_field_tensor(rn2d,zoff,dzobs,mu_head,fly_height)
! fftw plans are not thread safe so we can not use openMP here. boo.
   call output('  transforming magnetostatic tensor')
   print "('    ',$)"
   do k1=1,n_lay
       do i=1,9,3
         xi2d(:,:)=real(rn2d(:,:,i  ,k1),ffttr); call fft(xi2d(:,:),xo2d(:,:),nx,ny,.false.,xip); n2d(:,:,i  ,k1)=xo2d(:,:); print "('.',$)"
         yi2d(:,:)=real(rn2d(:,:,i+1,k1),ffttr); call fft(yi2d(:,:),yo2d(:,:),nx,ny,.false.,yip); n2d(:,:,i+1,k1)=yo2d(:,:); print "('.',$)"
         zi2d(:,:)=real(rn2d(:,:,i+2,k1),ffttr); call fft(zi2d(:,:),zo2d(:,:),nx,ny,.false.,zip); n2d(:,:,i+2,k1)=zo2d(:,:); print "('.',$)"
       enddo
   enddo
   call output('   -done')

   deallocate(rn2d)
   call cpu_time(end_time)
   i=int((end_time-start_time)/3600); k1=int(((end_time-start_time)-i*3600)/60)
   write(os,"('  tensor took:',i3,' hours, ',i2,' minutes, ',f5.2,' seconds')") i,k1,(end_time-start_time)-60*(k1+60*i)
   call output(trim(os),.true.)
   ok_p=.true.

  end function


 subroutine set_off_field_tensor(ihead_mu,isul_mu,ifly_height,ispace_to_sul,ifilename,inum_images)
  real(8),intent(in)::ihead_mu,isul_mu,ifly_height,ispace_to_sul
  integer,intent(in)::inum_images
  character(*),intent(in)::ifilename
  sul_mu=isul_mu; fly_height=ifly_height; space_to_sul=ispace_to_sul
  filename=trim(adjustl(ifilename))
  num_images=inum_images; head_mu=ihead_mu
  sul_image=(sul_mu.gt.1.d0); head_image=(head_mu.gt.1.d0)
 end subroutine

 subroutine calculate_off_field_tensor(n2d,zoff,dzobs,mu_head,height_fly)
  use io_m, only: output
  use mselem_m, only: plane_to_plane !,volume_averaged_n
  use data_m, only: get_geometry
  real(8),intent(in)::zoff,dzobs
  real(8),intent(out)::n2d(:,:,:,:)
  real(8),intent(in),optional::mu_head,height_fly

  logical::head_p,sul_p,ok_p
  integer::nx,ny,n_lay,imx,imy,ks,ko
  character(200)::os
  real(8)::dx,dy,zsb,zst,zob,zot,mus,mush
  real(8),allocatable::dz(:),il_space(:),nn(:,:,:,:),hdz(:)
  real(8)::head_mu0,fly_height0

  head_mu0=head_mu; if (present(mu_head)) head_mu0=mu_head
  head_image=(head_mu0.gt.1.d0)
  fly_height0=fly_height; if (present(height_fly)) fly_height0=height_fly*1.d-7
  nx=size(n2d,1); ny=size(n2d,2); n_lay=size(n2d,4)

  imx=0; if (periodic_x) imx=num_images
  imy=0; if (periodic_y) imy=num_images
  write(os,"(' generating magnetostatic tensor: nx=',i4,', ny=',i4,', layers=',i3,' (',i3,',',i3,')')") nx,ny,n_lay,imx,imy

  allocate(dz(n_lay),il_space(n_lay),nn(nx,ny,3,3),hdz(n_lay))

  call get_geometry(dx,dy,dz,il_space)

  mus=0.d0; mush=0.d0
  sul_p=(sul_image.and.sul_mu.ne.1.d0)
  if (sul_p) then
    write(os,"(a,', image in sul')") trim(os)
  else
    write(os,"(a,', no image in sul')") trim(os)
  endif
  head_p=(head_image.and.head_mu0.ne.1.d0)
  if (head_p) then
    write(os,"(a,', image in (infinite) head')") trim(os)
  else
    write(os,"(a,', no image in head')") trim(os)
  endif
  if (sul_p.and.head_p) then
    call output('')
    call output('**** WARNING **** imaging in both head and sul- code does not do imaging of images!')
    call output('')
  endif
  if (sul_p) mus=max(0.d0,(max(1.d0,sul_mu)-1.d0)/(max(1.d0,sul_mu)+1.d0))
  if (head_p) mush=max(0.d0,(max(1.d0,head_mu0)-1.d0)/(max(1.d0,head_mu0)+1.d0))

  if (imx.gt.0.and.imy.gt.0) then
    write(os,"(a,', periodic in plane')") trim(os)
  elseif (imx.gt.0) then
    write(os,"(a,', periodic in x, not y')") trim(os)
  elseif (imy.gt.0) then
    write(os,"(a,', periodic in y, not x')") trim(os)
  else
    write(os,"(a,', not periodic in plane')") trim(os)
  endif
  call output(trim(os))

    n2d=0.d0

    hdz=fly_height0
    do ks=2,n_lay
      hdz(ks)=hdz(ks-1)+dz(ks-1)+il_space(ks)
    enddo
    hdz=hdz+0.5d0*dz  !hdz(k) = distance from head to center of layer k
!print *,hdz

    zob=space_to_sul+zoff
    do ks=1,n_lay
      zob=zob+dz(ks)
      if (ks.ne.n_lay) zob=zob+il_space(ks)
    enddo
    zot=zob+dzobs
    

    zsb=space_to_sul
    do ks=n_lay,1,-1       !source layer
      zst=zsb+dz(ks)
!      zob=space_to_sul
!      do ko=n_lay,1,-1     !observation layer
!        zot=zob+dz(ko)

!print *,zsb,zst,zob,zot,0.5d0*(zst-zot+zsb-zob)

        print "('   computing layer #',i3,' contribution:',$)",ks
 !function plane_to_plane(nx,ny,dx,dy,dz,dzo,rz,n_imgxi,n_imgyi) result(n)
        nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dzobs,0.5d0*(zst-zot+zsb-zob),imx,imy)
        n2d(:,:,1:3,ks)=nn(:,:,:,1)
        n2d(:,:,4:6,ks)=nn(:,:,:,2)
        n2d(:,:,7:9,ks)=nn(:,:,:,3)

        if (sul_p) then
          call output('  -done')
          print "('                   and its SUL image:',$)"
!print *,ks,ko,mus,0.5d0*(-zst-zot-zsb-zob); read *
          nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dzobs,0.5d0*(-zst-zot-zsb-zob),imx,imy)
          nn(:,:,1:2,:)=-mus*nn(:,:,1:2,:); nn(:,:,3,:)=mus*nn(:,:,3,:) !xy component flip sign
          n2d(:,:,1:3,ks)=n2d(:,:,1:3,ks)+nn(:,:,:,1)
          n2d(:,:,4:6,ks)=n2d(:,:,4:6,ks)+nn(:,:,:,2)
          n2d(:,:,7:9,ks)=n2d(:,:,7:9,ks)+nn(:,:,:,3)
        endif
        if (head_p) then
          call output('  -done')
          print "('   computing layer #123 contribution:',$)",ks
          print "('       and its (infinite) head image:',$)"
!print *,hdz(ks),fly_height0-(zoff+0.5d0*dzobs),hdz(ks)+fly_height0-(zoff+0.5d0*dzobs)
          nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dzobs,hdz(ks)+fly_height0-(zoff+0.5d0*dzobs),imx,imy)
          nn(:,:,1:2,:)=-mush*nn(:,:,1:2,:); nn(:,:,3,:)=mush*nn(:,:,3,:) !xy component flip sign
          n2d(:,:,1:3,ks)=n2d(:,:,1:3,ks)+nn(:,:,:,1)
          n2d(:,:,4:6,ks)=n2d(:,:,4:6,ks)+nn(:,:,:,2)
          n2d(:,:,7:9,ks)=n2d(:,:,7:9,ks)+nn(:,:,:,3)
        endif

        n2d(:,:,:,ks)=n2d(:,:,:,ks)/(dx*dy*dzobs*nx*ny)

        call output('  -done')
      zsb=zst+il_space(ks)
    enddo

  deallocate(dz,il_space,nn,hdz)

 end subroutine

end module off_field_m
