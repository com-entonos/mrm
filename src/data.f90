module window_m
 implicit none
 private 
 public::set_window,get_window,show_windw

 type::window_s
   real(8)::v(2)
   real(8)::offset(2),off(2)
   integer::i0,j0
   integer::nx,ny
 end type window_s
 type(window_s),save::window
 logical,save::win_init=.false.

contains

  subroutine show_windw()
   use io_m, only:output
   character(200)::os

   write(os,"(' window: win_init=',l)") win_init; call output(os)
   write(os,"('  nx=',i4,', ny=',i4)") window%nx,window%ny; call output(os)
   write(os,"('  i0=',i4,', j0=',i4)") window%i0,window%j0; call output(os)
   write(os,"(1p,'  vx=',e12.5,', vy=',e12.5)") window%v(1),window%v(2); call output(os)
   write(os,"(1p,'  offsetx=',e12.5,', offsety=',e12.5)") window%offset(1),window%offset(2); call output(os)
   write(os,"(1p,'  offx=',e12.5,', offy=',e12.5)") window%off(1),window%off(2); call output(os)
  end subroutine

  subroutine set_window(v,offset,i0,j0,nx,ny)
   use fftw_m
   integer,intent(in),optional::i0,j0,nx,ny
   real(8),intent(in),optional::v(2),offset(2)

   if (.not.win_init) then
     win_init=.true.
     window%v=0.d0
     window%offset=0.d0
     window%i0=0; window%j0=0
     window%nx=0; window%ny=0
     if (present(offset)) window%offset=offset
   endif
   if (present(v)) window%v=v
   if (present(offset)) then
     if (window%v(1).eq.0.d0) window%offset(1)=offset(1)
     if (window%v(2).eq.0.d0) window%offset(2)=offset(2)
   endif
   if (present(i0)) window%i0=i0
   if (present(j0)) window%j0=j0
   if (present(nx)) window%nx=nx
   if (present(ny)) window%ny=ny
  end subroutine
  subroutine get_window(t,v,off,i0,nxr,nyr,n,offset)
    real(8),intent(in),optional::t
    real(8),intent(out),optional::v(2),off(2),offset(2)
    integer,intent(out),optional::i0(2),nxr(2),nyr(2),n(2)
    if (present(t)) window%off=window%v*t+window%offset
    if (present(off)) off=window%off
    if (present(offset)) offset=window%offset
    if (present(v)) v=window%v
    if (present(i0)) i0=(/window%i0,window%j0/)
    if (present(nxr)) nxr=(/window%i0+1,window%i0+window%nx/)
    if (present(nyr)) nyr=(/window%j0+1,window%j0+window%ny/)
    if (present(n)) n=(/window%nx, window%ny/)
  end subroutine
end module window_m


module data_m
 use FP_m, only: FP_s
 implicit none
 private
 public :: set_data,get_data,set_geometry,get_geometry,get_mag,init_data,get_size,file_grain,show_data,data_s,get_layer,update_do_grain,show_window

   interface get_geometry
     module procedure get_geometrya,get_geometry1
   end interface
   interface get_data
     module procedure get_datai,get_data1i,get_data2i,get_datar,get_data1r,get_data2r,get_datas,get_datal,get_datafp,get_data1l
   end interface
   interface set_data
     module procedure set_datai,set_data1i,set_data2i,set_datar,set_data1r,set_data2r,set_datas,set_datal,set_datafp,set_data1l
   end interface
   interface init_data
     module procedure init_data_one,init_data_all
   end interface

 integer,save::nx,ny,n_lay
 real(8),save::dx,dy
 real(8),allocatable,target,save::dz(:),il_space(:)

 type::data_s
   logical::thermal_p
   integer::in=0,min_n=0,max_n=0
   real(8)::alpha=1.d0,spin=0.5d0,gamma=1.76d7,floor_ms=0.01d0,background_t=0.d0,tau_hs_t=0.d0, thm_alpha=1.d0
   logical,dimension(:),pointer::do_grain=>null()
   real(8),dimension(:),pointer::tc=>null(),ms=>null(),grain_vol=>null(),grain_area=>null(),lambdaMs=>null(),tau=>null(),hthm_m=>null()
   real(8),dimension(:),pointer::ms_scale=>null(),alpha_scale=>null(),hk_scale=>null(),hx_scale=>null(),temp=>null(),energy=>null()
   real(8),dimension(:,:),pointer::mag(:,:)=>null(),h(:,:)=>null(),hthm(:,:)=>null()
   real(8),dimension(:,:),pointer::m=>null(),dm=>null()
   integer,dimension(:),pointer::grain2gridindex=>null()
   integer,dimension(:,:),pointer::grain2grid=>null(),grid2grain=>null()
   type(FP_s),dimension(:),pointer::fp=>null()
   character(80)::loopfile='h.dat'
 end type data_s

 type(data_s),allocatable,target,save::layer(:)
 real(8),dimension(:),pointer::ms_energy=>null()

contains

  subroutine show_window()
   use io_m, only:output
   use window_m, only: get_window, show_windw
   integer::n2(2)
   call get_window(n=n2)
   call show_windw()
   if (n2(1).eq.nx.and.n2(2).eq.ny) call output(' windowing is not used')
  end subroutine

  subroutine update_do_grain(t)
   use window_m, only: get_window, set_window
   real(8),intent(in)::t
   integer::i0(2),n2(2),k,i,j
   real(8)::offset(2)
   type(data_s),pointer::p

   call get_window(n=n2); if (n2(1).eq.nx.and.n2(2).eq.ny) return

   call get_window(t=t,off=offset)
   i=max(1,min(nx-n2(1)+1,nint((nx-n2(1))*0.5d0+offset(1)/dx)))-1
   j=max(1,min(ny-n2(2)+1,nint((ny-n2(2))*0.5d0+offset(2)/dy)))-1
   if (i.eq.i0(1).and.j.eq.i0(2)) return
   i0(1)=i; i0(2)=j
   call set_window(i0=i0(1),j0=i0(2))

   do k=1,n_lay; p=>layer(k)
     p%do_grain=.false.
     p%min_n=size(p%temp)+1; p%max_n=0
     p%m(4-p%in:6-p%in,:)=p%m(1+p%in:3+p%in,:)

 !   do n=1,size(p%temp)
 !     do n1=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); i=p%grain2grid(1,n1)-i0(1); j=p%grain2grid(2,n1)-i0(2)
 !       if (min(i,j).gt.0.and.i.le.n2(1).and.j.le.n2(2)) p%do_grain(n)=.true.
 !     enddo
 !   enddo

     do j=i0(2)+1,i0(2)+n2(2)
       do i=i0(1)+1,i0(1)+n2(1)
         if (p%grid2grain(i,j).gt.0) then
           p%do_grain(p%grid2grain(i,j))=.true.; p%min_n=min(p%min_n,p%grid2grain(i,j)); p%max_n=max(p%max_n,p%grid2grain(i,j))
         endif
       enddo
     enddo
   enddo


  end subroutine

  function get_layer(i) result (s)
   integer::i
   type(data_s),pointer::s
   nullify(s)
   if (i.gt.0.and.i.le.size(layer)) s=>layer(i)
  end function

  subroutine show_data(i)
   use io_m, only:output
   integer::i,n
   type(data_s),pointer::p
   character(200)::os

     p=>layer(i); n=0; if (associated(p%tc)) n=size(p%tc)
     write(os,"(' for layer #',i3)") i; call output(os)
     write(os,"('  nx=',i4,', ny=',i4)") nx,ny; call output(os)
     write(os,"(1p,'  dx=',e12.5,', dy=',e12.5,', dz=',e12.5)") dx,dy,dz(i); call output(os)
     if (i.gt.1) then
       write(os,"(1p,'  spacing to layer above=',e12.5)") il_space(i); call output(os)
     endif
     write(os,"('  m(h) loop file=',a)") trim(p%loopfile); call output(os)
     write(os,"(1p,'  alpha=',e12.5,' (',e12.5,'), gamma=',e12.5,', spin=',e12.5,', HS T=',e12.5)") p%alpha,p%thm_alpha,p%gamma,p%spin,p%tau_hs_t; call output(os)
     write(os,"(1p,'  background T=',e12.5,', floor_ms=',e12.5,', num_grain=',0p,i8)") p%background_t,p%floor_ms,n; call output(os)

  end subroutine

  function init_data_one(il,num_grain) result (ok_p)
   integer,intent(in)::il,num_grain
   type(data_s),pointer::p
   logical::ok_p

   p=>layer(il)
   if (associated(p%tc)) then
     if (size(p%tc).ne.num_grain) then
       deallocate(p%tc,p%ms,p%grain_vol,p%lambdaMs,p%temp,p%mag,p%energy,p%h,p%hthm,p%ms_scale,p%alpha_scale,p%hk_scale,p%hx_scale,p%m,p%fp,p%tau,p%do_grain,p%hthm_m)
       nullify(p%tc)
     endif
   endif
   if (.not.associated(p%tc)) then
     allocate(p%do_grain(num_grain),p%tc(num_grain),p%ms(num_grain),p%grain_vol(num_grain),p%lambdaMs(num_grain),p%temp(num_grain),p%mag(2,num_grain), &
             p%h(3,num_grain),p%hthm(3,num_grain),p%energy(num_grain),p%tau(num_grain), &
             p%ms_scale(num_grain),p%alpha_scale(num_grain),p%hk_scale(num_grain),p%hx_scale(num_grain),p%m(6,num_grain),p%dm(3,num_grain),p%fp(num_grain),p%hthm_m(num_grain))
     p%tc=1.d22; p%ms=0.d0; p%lambdaMs=0.d0; p%temp=0.d0; p%ms_scale=1.d0; p%alpha_scale=1.d0; p%hk_scale=1.d0; p%hx_scale=1.d0; p%m=0.d0
     p%m(3,:)=1.d0; p%m(6,:)=1.d0; p%mag=1.d0; p%h=0.d0; p%hthm=0.d0; p%energy=0.d0; p%dm=0.d0; p%tau=1.d-9; p%hthm_m=0.d0
     p%min_n=1; p%max_n=num_grain
   endif
   p%in=0; p%grain_vol=p%grain_area*dz(il); ok_p=.true.; p%do_grain=.true.
  end function

  function init_data_all(num_grain) result (ok_p)
   integer,intent(in)::num_grain(:)
   integer::i
   logical::ok_p
   
   ok_p=.true.
   do i=1,size(layer)
     ok_p=(ok_p.and.init_data_one(i,num_grain(i)))
   enddo
  end function

  function file_grain(read_p,ionum,il,ng,old_p) result (ok_p)
   logical,intent(in)::read_p,old_p
   integer,intent(in)::ionum,il
   integer,intent(out)::ng

   logical::ok_p
   integer::i1,i2,i3,n,n2
   real(8)::x,y
   type(data_s),pointer::p

   logical::gp
   integer::gi,gdi(2:n_lay)
   integer::gia(nx*ny)
   real(8)::g8
   real(8),allocatable::g8a1(:,:),g8a2(:,:),g8a3(:,:,:),g8a4(:)
   character(31)::gs

   ok_p=.false.
   if (il.lt.1.or.il.gt.size(layer)) return

   p=>layer(il)
   if (read_p) then
     if (old_p) then
!     read(ionum,err=99,end=99) &
!              ms%nx,ms%ny,ms%dx,ms%dy,ms%nrand,ms%nseed,ms%k_ang,ms%kx,ms%ky, &
!              ms%kz,ms%hk_cone,ms%k_cone,ms%hk_3d_random,ms%hk_2d_random,ms%hk_uniform,ms%k_sigma,ms%k1,ms%k2, &
!              ms%intergrain_hx,ms%intergrain_hx_sigma,ms%packing_fraction, ms%grain, &
!              ms%num_grain,ms%vol_sigma,ms%initial_grain_area,ms%dr_grain, &
!              ms%Tc_sigma,ms%Tc_nom,ms%small_packing_fraction,ms%num_small_grain,ms%init_vol_sigma, &
!              ms%initial_ext_grain_area,ms%small_grain_sig,ms%scale_hk,ms%kmax,ms%frac_in_plane,ms%frac_random,&
!              i,dummya
       read(ionum,err=99) i1,i2,x,y,gi,gi,g8,g8,g8, &
                g8,gp,g8,gp,gp,gp,g8,g8,g8, &
                g8,g8,g8,gp, &
                ng,g8,g8,g8, &
                g8,g8,g8,gi,g8, &
                g8,g8,gp,g8,g8,g8,&
                gi,gdi !,dummya
       if (i1.ne.nx.or.i2.ne.ny) then
          ok_p=.false.; go to 99
       endif
       read(ionum,err=99) gs
       if (associated(p%grain_area)) then; if (size(p%grain_area).ne.ng) then; deallocate(p%grain_area); nullify(p%grain_area);endif; endif
       if (associated(p%grain2gridindex)) then; if (size(p%grain2gridindex).ne.ng+1) then; deallocate(p%grain2gridindex); nullify(p%grain2gridindex);endif; endif
       if (.not.associated(p%grain2gridindex)) allocate(p%grain2gridindex(0:ng));
       if (.not.associated(p%grain_area)) allocate(p%grain_area(ng))
       if (.not.associated(p%grid2grain)) allocate(p%grid2grain(nx,ny));

       allocate(g8a1(3,ng),g8a2(2,ng),g8a3(3,3,ng),g8a4(ng))
       read(ionum,err=99) g8a1,g8a2,gia,g8a3,g8a4,g8a1
       backspace(ionum)
       p%grid2grain=0; do i2=1,ny; do i1=1,nx; p%grid2grain(i1,i2)=gia(i1+nx*(i2-1)); enddo; enddo

       n2=count(p%grid2grain.gt.0)
       if (associated(p%grain2grid)) then; if (size(p%grain2grid,2).ne.n2)then; deallocate(p%grain2grid); nullify(p%grain2grid);endif; endif
       if (.not.associated(p%grain2grid)) allocate(p%grain2grid(2,n2));

       p%grain2gridindex=0; p%grain2grid=0; i3=0
       do n=1,ng
         do i2=1,ny; do i1=1,nx
           if (p%grid2grain(i1,i2).eq.n) then
             i3=i3+1; p%grain2grid(1,i3)=i1; p%grain2grid(2,i3)=i2
           endif
         enddo; enddo
         p%grain2gridindex(n)=i3; p%grain_area(n)=dx*dy*(p%grain2gridindex(n)-p%grain2gridindex(n-1))
       enddo
       ok_p=init_data_one(il,ng)
       p%tc=g8a4
       deallocate(g8a1,g8a2,g8a3,g8a4)
     else
       read(ionum,err=99) i1,i2,i3,x,y,n,n2; ng=n
       if (i1.ne.nx.or.i2.ne.ny) then
          ok_p=.false.; go to 99
       endif
       if (associated(p%grain_area)) then; if (size(p%grain_area).ne.n) then; deallocate(p%grain_area); nullify(p%grain_area);endif; endif
       if (associated(p%grain2grid)) then; if (size(p%grain2grid,2).ne.n2)then; deallocate(p%grain2grid); nullify(p%grain2grid);endif; endif
       if (associated(p%grain2gridindex)) then; if (size(p%grain2gridindex).ne.n+1) then; deallocate(p%grain2gridindex); nullify(p%grain2gridindex);endif; endif
       if (.not.associated(p%grain2gridindex)) allocate(p%grain2gridindex(0:n));
       if (.not.associated(p%grain2grid)) allocate(p%grain2grid(2,n2));
       if (.not.associated(p%grain_area)) allocate(p%grain_area(n))
       if (.not.associated(p%grid2grain)) allocate(p%grid2grain(nx,ny));
       read(ionum,err=99) p%grid2grain,p%grain2grid,p%grain2gridindex,p%grain_area
       ok_p=init_data_one(il,n)
       read(ionum,err=99) p%ms,p%tc,p%tau
     endif
     if (associated(p%do_grain)) then; if (size(p%do_grain).ne.ng) then; deallocate(p%do_grain); nullify(p%do_grain); endif; endif
     if (.not.associated(p%do_grain)) allocate(p%do_grain(ng))
     p%do_grain=.true.
     ok_p=.true.
   else
     ng=size(p%tc)
     write(ionum,err=99) nx,ny,il,dx,dy,size(p%tc),size(p%grain2grid,2)
     write(ionum,err=99) p%grid2grain,p%grain2grid,p%grain2gridindex,p%grain_area
     write(ionum,err=99) p%ms,p%tc,p%tau
     ok_p=.true.
   endif
   return

  99 continue
   ok_p=.false.

  end function

  function get_mag(mx,my,mz) result (ok_p)
   use, intrinsic :: iso_c_binding
   use window_m, only: get_window, set_window
   use fftw_m, only: ffttr
   real(ffttr),intent(out)::mx(:,:,:),my(:,:,:),mz(:,:,:)
!  real(8),intent(inout)::offset(:)
   logical::ok_p

   integer::n,n1,i,j,k,id,nxf,nyf,i0(2)
   real(ffttr)::m(3)!,x,offset(2)
   type(data_s),pointer::p

   ok_p=.false.
   nxf=size(mx,1); nyf=size(mx,2)
   if (nxf.eq.nx.and.nyf.eq.ny) then
!    do k=1,n_lay; p=>layer(k); id=p%in
!      do j=1,ny; do i=1,nx
!        n=p%grid2grain(i,j)
!        m=min(1,n)*real(p%ms(max(1,n))*p%ms_scale(max(1,n))*p%m(1+id:3+id,max(1,n)))
!        mx(i,j,k)=m(1); my(i,j,k)=m(2); mz(i,j,k)=m(3)
!      enddo; enddo
!    enddo
     mx=0.0_ffttr; my=0.0_ffttr; mz=0.0_ffttr
     do k=1,n_lay; p=>layer(k); id=p%in
!$omp parallel
!$omp do private(n,m,n1,i,j) firstprivate(p,k,id)
       do n=1,size(p%temp)
         m=real(p%ms(n)*p%ms_scale(n)*p%m(1+id:3+id,n),ffttr)
         do n1=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); i=p%grain2grid(1,n1); j=p%grain2grid(2,n1)
           mx(i,j,k)=m(1); my(i,j,k)=m(2); mz(i,j,k)=m(3)
         enddo
       enddo
!$omp end do nowait
!$omp end parallel
     enddo
   else
     
     call get_window(i0=i0)
!    offset(1)=offset(1)/dx; offset(2)=offset(2)/dy; 
!    i0=max(1,min(nx-nxf+1,nint((nx-nxf)*0.5d0+offset(1))))-1
!    j0=max(1,min(ny-nyf+1,nint((ny-nyf)*0.5d0+offset(2))))-1
!    do k=1,n_lay; p=>layer(k); id=p%in
!      do j=1,nyf; do i=1,nxf
!        n=p%grid2grain(i+i0,j+j0)
!        m=min(1,n)*real(p%ms(max(1,n))*p%ms_scale(max(1,n))*p%m(1+id:3+id,max(1,n)),ffttr)
!        mx(i,j,k)=m(1); my(i,j,k)=m(2); mz(i,j,k)=m(3)
!      enddo; enddo
!    enddo
     mx=0.0_ffttr; my=0.0_ffttr; mz=0.0_ffttr
     do k=1,n_lay; p=>layer(k); id=p%in
!$omp parallel
!$omp do private(n,m,n1,i,j) firstprivate(p,k,id)
!      do n=1,size(p%temp)
       do n=p%min_n,p%max_n
         if (p%do_grain(n)) then
           m=real(p%ms(n)*p%ms_scale(n)*p%m(1+id:3+id,n),ffttr)
           do n1=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); i=p%grain2grid(1,n1)-i0(1); j=p%grain2grid(2,n1)-i0(2)
             if (min(i,j).gt.0.and.i.le.nxf.and.j.le.nyf) then
               mx(i,j,k)=m(1); my(i,j,k)=m(2); mz(i,j,k)=m(3)
             endif
           enddo
         endif
       enddo
!$omp end do nowait
!$omp end parallel
     enddo
   endif
   ok_p=.true.
  end function

  subroutine get_size(nxi,nyi,n_layi,dxi,dyi)
   integer,intent(out)::nxi,nyi,n_layi
   real(8),intent(out)::dxi,dyi
   nxi=nx; nyi=ny; n_layi=n_lay; dxi=dx; dyi=dy
  end subroutine

  subroutine get_geometrya(dxi,dyi,dzi,il_spacei)
   real(8),intent(out)::dxi,dyi,dzi(n_lay),il_spacei(n_lay)
   dxi=dx; dyi=dy; dzi=dz; il_spacei=il_space
  end subroutine
  subroutine get_geometry1(n,dxi,dyi,dzi)
   integer,intent(in)::n
   real(8),intent(out)::dxi,dyi,dzi
   dxi=dx; dyi=dy; dzi=0.d0; if (n.gt.0.and.n.le.n_lay) dzi=dz(n)
  end subroutine
  subroutine set_geometry(nxi,nyi,dxi,dyi,n_layi,dzi,il_spacei)
   integer,intent(in)::nxi,nyi,n_layi
   real(8),intent(in)::dxi,dyi,dzi(n_layi),il_spacei(n_layi)
   nx=nxi; ny=nyi; n_lay=n_layi; dx=dxi; dy=dyi
   allocate(dz(n_lay),il_space(n_lay)); dz=dzi; il_space=il_spacei
   allocate(layer(n_lay))
  end subroutine

  subroutine set_datas(il,str1) !,dz,il_space)
   integer,intent(in)::il
   character(*)::str1
   integer::i
!  logical,intent(in),optional::in_p
   if (il.lt.0.or.il.gt.n_lay) then
     return
   elseif (il.eq.0) then
     do i=1,n_lay
       layer(i)%loopfile=trim(str1)
     enddo
   else
     layer(il)%loopfile=trim(str1)
   endif
  end subroutine
  subroutine get_datas(il,str1) !,dz,il_space)
   integer,intent(in)::il
   character(*)::str1
!  logical,intent(in),optional::in_p
   if (il.lt.1.or.il.gt.n_lay) then
     return
   else
     str1=trim(layer(il)%loopfile)
   endif
  end subroutine
  subroutine get_datai(il,data,in_p,max_n_p,min_n_p) !,dz,il_space)
   integer,intent(in)::il
   integer,pointer::data
   logical,intent(in),optional::in_p,max_n_p,min_n_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(in_p)) then
     data=>layer(il)%in
   elseif (present(min_n_p)) then
     data=>layer(il)%min_n
   elseif (present(max_n_p)) then
     data=>layer(il)%max_n
   endif
  end subroutine
  subroutine set_datai(il,data,in_p,max_n_p,min_n_p) !,dz,il_space)
   integer,intent(in)::il
   integer::data
   logical,intent(in),optional::in_p,max_n_p,min_n_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(in_p)) then
     layer(il)%in=data
   elseif (present(min_n_p)) then
     layer(il)%min_n=data
   elseif (present(max_n_p)) then
     layer(il)%max_n=data
   endif
  end subroutine
  subroutine get_datal(il,data,thermal_p) !,dz,il_space)
   integer,intent(in)::il
   logical,pointer::data
   logical,intent(in),optional::thermal_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(thermal_p)) then
     data=>layer(il)%thermal_p
   endif
  end subroutine
  subroutine get_data1l(il,data,do_grain_p) !,dz,il_space)
   integer,intent(in)::il
   logical,pointer::data(:)
   logical,intent(in),optional::do_grain_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(do_grain_p)) then
      data=>layer(il)%do_grain
   endif
  end subroutine
  subroutine set_datal(il,data,thermal_p) !,dz,il_space)
   integer,intent(in)::il
   logical::data
   logical,intent(in),optional::thermal_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(thermal_p)) then
     layer(il)%thermal_p=data
   endif
  end subroutine
  subroutine set_data1l(il,data,do_grain_p) !,dz,il_space)
   integer,intent(in)::il
   logical,pointer::data(:)
   logical,intent(in),optional::do_grain_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(do_grain_p)) then
      layer(il)%do_grain=>data
   endif
  end subroutine

  subroutine get_data1i(il,data,grain2gridindex_p)
   integer,intent(in)::il
   integer,pointer::data(:)
   logical,intent(in),optional::grain2gridindex_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(grain2gridindex_p)) then
     data=>layer(il)%grain2gridindex
   endif
  end subroutine
  subroutine set_data1i(il,data,grain2gridindex_p)
   integer,intent(in)::il
   integer,pointer::data(:)
   logical,intent(in),optional::grain2gridindex_p

   if (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(grain2gridindex_p)) then
     layer(il)%grain2gridindex=>data
   endif
  end subroutine

  subroutine get_data2i(il,data,grain2grid_p,grid2grain_p)
   integer,intent(in)::il
   integer,pointer::data(:,:)
   logical,intent(in),optional::grain2grid_p,grid2grain_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(grain2grid_p)) then
     data=>layer(il)%grain2grid
   elseif (present(grid2grain_p)) then
     data=>layer(il)%grid2grain
   endif
  end subroutine
  subroutine set_data2i(il,data,grain2grid_p,grid2grain_p)
   integer,intent(in)::il
   integer,pointer::data(:,:)
   logical,intent(in),optional::grain2grid_p,grid2grain_p

   if (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(grain2grid_p)) then
     layer(il)%grain2grid=>data
   elseif (present(grid2grain_p)) then
     layer(il)%grid2grain=>data
   endif
  end subroutine

  subroutine get_datar(il,data,thm_alpha_p,alpha_p,spin_p,gamma_p,floor_ms_p,dz_p,il_space_p,background_t_p,tau_hs_t_p)
   integer,intent(in)::il
   real(8),pointer::data
   logical,intent(in),optional::alpha_p,spin_p,gamma_p,floor_ms_p,dz_p,il_space_p,background_t_p,tau_hs_t_p,thm_alpha_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(thm_alpha_p)) then
     data=>layer(il)%thm_alpha
   elseif (present(alpha_p)) then
     data=>layer(il)%alpha
   elseif (present(spin_p)) then
     data=>layer(il)%spin
   elseif (present(gamma_p)) then
     data=>layer(il)%gamma
   elseif (present(floor_ms_p)) then
     data=>layer(il)%floor_ms
   elseif (present(background_t_p)) then
     data=>layer(il)%background_t
   elseif (present(tau_hs_t_p)) then
     data=>layer(il)%tau_hs_t
   elseif (present(dz_p)) then
     data=>dz(il)
   elseif (present(il_space_p)) then
     data=>il_space(il)
   endif
  end subroutine
  subroutine set_datar(il,data,thm_alpha_p,alpha_p,spin_p,gamma_p,floor_ms_p,dz_p,il_space_p,tau_hs_t_p,background_t_p)
   integer,intent(in)::il
   real(8)::data
   logical,intent(in),optional::alpha_p,spin_p,gamma_p,floor_ms_p,dz_p,il_space_p,tau_hs_t_p,background_t_p,thm_alpha_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(thm_alpha_p)) then
     layer(il)%thm_alpha=data
   elseif (present(alpha_p)) then
     layer(il)%alpha=data
     layer(il)%thm_alpha=data
   elseif (present(spin_p)) then
     layer(il)%spin=data
   elseif (present(gamma_p)) then
     layer(il)%gamma=data
   elseif (present(floor_ms_p)) then
     layer(il)%floor_ms=data
   elseif (present(background_t_p)) then
     layer(il)%background_t=data
   elseif (present(tau_hs_t_p)) then
     layer(il)%tau_hs_t=data
   elseif (present(dz_p)) then
     dz(il)=data
   elseif (present(il_space_p)) then
     il_space(il)=data
   endif
  end subroutine


  subroutine get_data1r(il,data,tau_p,tc_p,ms_p,grain_vol_p,grain_area_p,energy_p,lambdaMs_p,temp_p,ms_scale_p,alpha_scale_p,hk_scale_p,hx_scale_p,ms_energy_p,hthm_p)
   integer,intent(in)::il
   real(8),pointer::data(:)
   logical,intent(in),optional::tau_p,tc_p,ms_p,grain_vol_p,grain_area_p,energy_p,lambdaMs_p,temp_p,ms_scale_p,alpha_scale_p,hk_scale_p,hx_scale_p,ms_energy_p
   logical,intent(in),optional::hthm_p

   nullify(data)
   if (present(ms_energy_p)) then
     data=>ms_energy
   elseif (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(temp_p)) then
     data=>layer(il)%temp
   elseif (present(ms_p)) then
     data=>layer(il)%ms
   elseif (present(tc_p)) then
     data=>layer(il)%tc
   elseif (present(tau_p)) then
     data=>layer(il)%tau
   elseif (present(ms_scale_p)) then
     data=>layer(il)%ms_scale
   elseif (present(alpha_scale_p)) then
     data=>layer(il)%alpha_scale
   elseif (present(hk_scale_p)) then
     data=>layer(il)%hk_scale
   elseif (present(hx_scale_p)) then
     data=>layer(il)%hx_scale
   elseif (present(grain_vol_p)) then
     data=>layer(il)%grain_vol
   elseif (present(grain_area_p)) then
     data=>layer(il)%grain_area
   elseif (present(energy_p)) then
     data=>layer(il)%energy
   elseif (present(lambdaMs_p)) then
     data=>layer(il)%lambdaMs
   elseif (present(hthm_p)) then
     data=>layer(il)%hthm_m
   endif
  end subroutine
  subroutine set_data1r(il,data,tau_p,tc_p,ms_p,grain_vol_p,grain_area_p,energy_p,lambdaMs_p,temp_p,ms_scale_p,alpha_scale_p,hk_scale_p,hx_scale_p,ms_energy_p,hthm_p)
   integer,intent(in)::il
   real(8),pointer::data(:)
   logical,intent(in),optional::tau_p,tc_p,ms_p,grain_vol_p,grain_area_p,energy_p,lambdaMs_p,temp_p,ms_scale_p,alpha_scale_p,hk_scale_p,hx_scale_p,ms_energy_p
   logical,intent(in),optional::hthm_p

   if (present(ms_energy_p)) then
     ms_energy=>data
   elseif (il.lt.1.or.il.gt.n_lay) then
     return
   elseif (present(temp_p)) then
     layer(il)%temp=>data
   elseif (present(ms_p)) then
     layer(il)%ms=>data
   elseif (present(tc_p)) then
     layer(il)%tc=>data
   elseif (present(tau_p)) then
     layer(il)%tau=>data
   elseif (present(ms_scale_p)) then
     layer(il)%ms_scale=>data
   elseif (present(alpha_scale_p)) then
     layer(il)%alpha_scale=>data
   elseif (present(hk_scale_p)) then
     layer(il)%hk_scale=>data
   elseif (present(hx_scale_p)) then
     layer(il)%hx_scale=>data
   elseif (present(grain_vol_p)) then
     layer(il)%grain_vol=>data
   elseif (present(grain_area_p)) then
     layer(il)%grain_area=>data
   elseif (present(energy_p)) then
     layer(il)%energy=>data
   elseif (present(lambdaMs_p)) then
     layer(il)%lambdaMs=>data
   elseif (present(hthm_p)) then
     layer(il)%hthm_m=>data
   endif
  end subroutine

  subroutine get_data2r(il,data,m_p,dm_p,h_p,mag_p,hthm_p) !,dz,il_space)
   integer,intent(in)::il
   real(8),pointer::data(:,:)
   logical,intent(in),optional::m_p,dm_p,h_p,mag_p,hthm_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(m_p)) then
     data=>layer(il)%m
   elseif (present(h_p)) then
     data=>layer(il)%h
   elseif (present(hthm_p)) then
     data=>layer(il)%hthm
   elseif (present(dm_p)) then
     data=>layer(il)%dm
   elseif (present(mag_p)) then
     data=>layer(il)%mag
   endif
  end subroutine
  subroutine set_data2r(il,data,m_p,dm_p,h_p,mag_p,hthm_p) !,dz,il_space)
   integer,intent(in)::il
   real(8),pointer::data(:,:)
   logical,intent(in),optional::m_p,dm_p,h_p,mag_p,hthm_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(m_p)) then
     layer(il)%m=>data
   elseif (present(h_p)) then
     layer(il)%h=>data
   elseif (present(hthm_p)) then
     layer(il)%hthm=>data
   elseif (present(dm_p)) then
     layer(il)%dm=>data
   elseif (present(mag_p)) then
     layer(il)%mag=>data
   endif
  end subroutine
  subroutine get_datafp(il,data,fp_p)
   integer,intent(in)::il
   type(FP_s),pointer::data(:)
   logical,intent(in),optional::fp_p

   nullify(data)
   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(fp_p)) then
      data=>layer(il)%fp
   endif
  end subroutine
  subroutine set_datafp(il,data,fp_p)
   integer,intent(in)::il
   type(FP_s),pointer::data(:)
   logical,intent(in),optional::fp_p

   if (il.lt.1.or.il.gt.n_lay) then
      return
   elseif (present(fp_p)) then
     layer(il)%fp=>data
   endif
  end subroutine

end module data_m
