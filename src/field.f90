!this module finds all of the intrinsic fields acting on a grain
! it then adds on the applied magnetic fields and demag fields
! it also gets the applied temperature field and updates the temperature scaling


module field_m
 use, intrinsic :: iso_c_binding
 use random_m, only: random_s
 use fftw_m, only: ffttr

 implicit none
 private
 public :: construct_field,init_applied,start_field,get_applied,get_field,get_easy_axis,file_field,plot_hamr,set_mode_field,freeze_time,plot_field,set_source,show_grain,find_grain

 type field_ss
   logical::source_p=.false.
   real(8)::last_t=-1.d22,avetau=0.d0
   real(8),pointer::tau_hs_t=>null()
   logical,pointer::thermal_p=>null()
   logical,pointer,dimension(:)::do_grain=>null()
   integer,pointer::in=>null(),max_n=>null(),min_n=>null()
   integer,dimension(:),pointer::grain2gridindex=>null()
   integer,dimension(:,:),pointer::grain2grid=>null(),grid2grain=>null()
   real(8),pointer::dz=>null(),background_t=>null(),floor_ms=>null(),alpha=>null(),gamma=>null()
   real(8),dimension(:),pointer::grain_vol=>null(),grain_area=>null(),tc=>null(),ms=>null(),temp=>null(),energy=>null(),tau=>null()
   real(8),dimension(:,:),pointer::m=>null(),mag(:,:)=>null()
   real(8),dimension(:),pointer::ms_scale=>null(),alpha_scale=>null(),hk_scale=>null(),hx_scale=>null()
   integer,dimension(:),pointer::hxgrainindex=>null(),hxgrain=>null(),lhxgrainindex=>null(),lhxgrain=>null(),uhxgrainindex=>null(),uhxgrain=>null()
   real(8),dimension(:),pointer::hx=>null(),lhx=>null(),uhx=>null(),temp_app=>null(),hthm_m=>null()
   real(8),dimension(:,:),pointer::uni_axis=>null(),hk=>null(),hthm=>null(),h=>null(),h_demag=>null(),h_app=>null(),h_selfdemag=>null(),h_plot=>null()
   real(8),dimension(:,:),pointer::h_ext=>null()
   real(8),dimension(:,:,:),pointer::cubic_axis=>null()
   real(8),dimension(:,:,:,:),pointer::Ndemag=>null()
   type(random_s),pointer::random_thml
 end type field_ss

 type(field_ss),allocatable,target,save::field(:)

 real(8),dimension(:),pointer,save::utemp=>null(),ms_energy=>null()
 real(8),dimension(:,:),pointer,save::ufield=>null()
 real(8),dimension(:,:,:),pointer,save::app_temp=>null(),app_plot=>null()
 real(8),dimension(:,:,:,:),pointer,save::app_field=>null()

 logical,dimension(:,:),allocatable,save::changed
 logical,save::compute_selfdemag_p=.true.,selfdemag_p=.false.,grain_selfdemag_p=.true.

 real(8),save::frozen_time=-1.d300
 integer,save::field_mode=0           ! 0- means compute and collect all fields to grains
                                      ! 1- compute intrinsic fields and add onto whatever is there (intrinsic- self demag + Hk + Hex + ...)
                                      !-1- compute demag and collect it and applied fields to grains subtracting off self-demag

 real(ffttr),dimension(:,:,:),pointer,save::demag_xs=>null(),demag_ys,demag_zs

contains
  subroutine find_grain(i,ig,jg)
   use io_m, only:output,add_file,delete_file
   integer,intent(in)::i,ig,jg
   type(field_ss),pointer::p
   character(200)::os
   integer::i0,i1,il,q,r,ionum

     i0=i; i1=i
     if (i.lt.1.or.i.gt.size(field)) then; i0=1; i1=size(field); endif
     do il=i0,i1
       p=>field(il)
       if (ig.gt.0.and.ig.le.size(p%grid2grain,1).and.jg.gt.0.and.jg.le.size(p%grain2grid,2)) then
         write(os,"(' grain at ',i0, ', ',i0,' in layer # ',i0,' is # ',i0)") ig,jg,il,p%grid2grain(ig,jg); call output(os)
       endif
     enddo

     if (ig.eq.0.and.jg.eq.0) then
       call output ('  dumping out grid with grain index...`grid_grain.dat`')
       ionum=add_file('grid_grain.dat'); call delete_file(ionum)
       open(ionum,file='grid_grain.dat',status='unknown')
       do il=i0,i1
         p=>field(il)
         do r=1,size(p%grid2grain,2)
           do q=1,size(p%grid2grain,1)
             write(ionum,"(4(x,i0))") q,r,il,p%grid2grain(q,r)
           enddo
         enddo
       enddo
       close(ionum)
     endif
  end subroutine

  subroutine show_grain(i,n)
   use io_m, only:output
   integer,intent(in)::i,n
   type(field_ss),pointer::p,q
   character(3000)::os
   integer::i0,i1,il,j

     i0=i; i1=i
     if (i.lt.1.or.i.gt.size(field)) then; i0=1; i1=size(field); endif
     do il=i0,i1
       p=>field(il)
       if (n.gt.0.and.n.le.size(p%tc)) then
       write(os,"(' grain #',i0,' in layer #',i0)") n,il; call output(os)
       write(os,"('   dz=',1p,e12.5,' cm, area=',e12.5,' cm^2, vol=',e12.5,' cm^3')") p%dz,p%grain_area(n),p%grain_vol(n); call output(os)
       write(os,"('   Tc=',1p,e12.5,' K, |Hk|=',e12.5,' (',e12.5,') Oe, Ms=',e12.5,' emu/cc')") p%tc(n),p%hk(:,n),p%ms(n); call output(os)
       write(os,"('   T=',1p,e12.5,' K, |m|=',e12.5,' (',e12.5,'), Hk(T)/Hk(0)=',e12.5)") p%temp(n),p%mag(p%in/3+1,n),p%ms_scale(n),p%hk_scale(n); call output(os)
       write(os,"('   m=',1p,3(x,e12.5))") p%m(1+p%in:3+p%in,n); call output(os)
       write(os,"('   k=',1p,3(x,e12.5))") p%uni_axis(:,n); call output(trim(os))
!      call output('   nodes:'); do j=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); write(os,"(' ',i0,' ',i0)") p%grain2grid(1:2,j); call output(trim(os)); enddo
       write(os,"('   number of nodes: ',i0)") p%grain2gridindex(n)-p%grain2gridindex(n-1); call output(os)
       write(os,"('   nodes:')"); do j=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); write(os,"(a,' (',i0,',',i0,')')") trim(os),p%grain2grid(1:2,j); enddo; call output(os)
!      write(os,"('   nodes:',:,5000(' (',i0,',',i0,')'))") (p%grain2grid(1:2,j),j=p%grain2gridindex(n-1)+1,p%grain2gridindex(n)); call output(os)
       write(os,"('   Happ=',1p,3(x,e12.5))") p%h_app(:,n); call output(trim(os))
       write(os,"('   Hext=',1p,3(x,e12.5))") p%h_ext(:,n); call output(trim(os))
       write(os,"('   Hkcubic=',1p,3(x,e12.5))") get_cubic_hk(n,p,p%in); call output(trim(os))
       write(os,"('   Hkuniax=',1p,3(x,e12.5))") get_uni_hk(n,p,p%in); call output(trim(os))
       if (associated(p%hxgrain)) then; write(os,"('   Hlat=',1p,3(x,e12.5))") get_intralayer_hx(n,p,p%in); call output(trim(os)); endif
       if (associated(p%uhx)) then; q=>field(il-1)
         write(os,"('   Hvabove=',1p,3(x,e12.5))") get_interlayer_hx(n,p,q,p%uhxgrainindex,p%uhxgrain,p%uhx); call output(trim(os))
       endif
       if (associated(p%lhx)) then; q=>field(il+1)
         write(os,"('   Hvbelow=',1p,3(x,e12.5))") get_interlayer_hx(n,p,q,p%lhxgrainindex,p%lhxgrain,p%lhx); call output(trim(os))
       endif
       write(os,"('   Hdemag=',1p,3(x,e12.5))") p%h_demag(:,n); call output(trim(os))
       write(os,"('   Hselfdemag=',1p,3(x,e12.5))") p%h_selfdemag(:,n); call output(trim(os))
       write(os,"('   H=',1p,3(x,e12.5))") p%h(:,n); call output(trim(os))
       else
       write(os,"(' no such grain: #',i0,' in layer #',i0)") n,il; call output(os)
       endif
     enddo
 end subroutine

 function set_source(source_p,ilayer) result (ok_p)
  logical,intent(in)::source_p
  integer,optional,intent(in)::ilayer
  logical::ok_p
  integer::i,k1,k2

  ok_p=.false.
  k1=1; k2=size(field)
  if (present(ilayer)) then
    if (ilayer.gt.0.and.ilayer.le.k2) then; k1=ilayer; k2=ilayer; else; endif
  endif
  do i=k1,k2
    field(i)%source_p=source_p
  enddo
 end function

 function set_mode_field(mode) result (ok_p)
  integer,intent(in)::mode
  logical::ok_p
  ok_p=(mode.lt.-1.or.mode.gt.1)
  if (.not.ok_p) field_mode=mode
 end function

 subroutine plot_hamr()
  use temp_scaling_m, only:plot_temp_scaling
  integer::i
  do i=1,size(field)
    call plot_temp_scaling(i,field(i)%ms,field(i)%alpha,field(i)%hk,field(i)%hx,field(i)%tc,field(i)%grain_area)
  enddo
 end subroutine

 subroutine get_field(t,h_var,thermal_p,varying_p,no_ms_p,plot_p,diff_m,update_do_grain_p,rot_phi)
  use data_m,only: get_data
  real(8),intent(in)::t
  real(8),intent(out),optional::h_var
  logical,intent(in),optional::thermal_p,plot_p,no_ms_p,diff_m,update_do_grain_p
  logical,intent(out),optional::varying_p(:)
  real(8),intent(in),optional::rot_phi(2)

  logical::thermal,v_p(2),p_p,nm_p,diffm
  integer::i,n,ng,id
  real(8)::var
  type(field_ss),pointer::p,q

  if (present(update_do_grain_p)) then; if (update_do_grain_p) then  !support for decay
    do i=1,size(field); call get_data(i,field(i)%do_grain,do_grain_p=.true.); enddo
  endif; endif
  nm_p=.false.; if (present(no_ms_p)) nm_p=no_ms_p
  thermal=.false.; if (present(thermal_p)) thermal=thermal_p
  p_p=.false.; if (present(plot_p)) p_p=plot_p
  diffm=.false.; if (present(diff_m)) diffm=diff_m
! call get_data(1,ms_energy,ms_energy_p=.true.)

  if (field_mode.ne.1) then
    v_p=get_applied(t,plot_p=p_p,rot_phi=rot_phi)
    if (selfdemag_p) then
      call get_selfdemag()
    else
      call get_magnetostatic()
    endif
    if (field_mode.eq.-1) call get_selfdemag()
  else
    v_p=.false.
    call get_selfdemag()
  endif
  var=0.d0; if (field_mode.ne.0) var=0.5d0

  if (associated(ms_energy)) ms_energy=0.d0
  do i=1,size(field); p=>field(i);ng=size(p%temp); id=p%in; p%energy=0.d0
    if (associated(p%hxgrain).and.field_mode.ne.1) then
      do n=p%min_n,p%max_n; if (p%do_grain(n)) then
        p%energy(n) = -dot_product(p%m(1+id:3+id,n),p%h_demag(1:3,n)+p%h_app(1:3,n)-var*p%h_selfdemag(1:3,n))*p%ms(n)*p%ms_scale(n)
        p%h(:,n)=p%h_app(:,n)+p%h_demag(:,n)+get_cubic_hk(n,p,id)+get_uni_hk(n,p,id)+get_intralayer_hx(n,p,id); endif
      enddo
      if (associated(ms_energy)) ms_energy=ms_energy+p%energy*p%grain_vol
    else
      do n=p%min_n,p%max_n; if (p%do_grain(n)) then
        p%energy(n) = -dot_product(p%m(1+id:3+id,n),p%h_demag(1:3,n)+p%h_app(1:3,n)-var*p%h_selfdemag(1:3,n))*p%ms(n)*p%ms_scale(n)
        p%h(:,n)=p%h_app(:,n)+p%h_demag(:,n)+get_cubic_hk(n,p,id)+get_uni_hk(n,p,id); endif
      enddo
      if (associated(ms_energy)) ms_energy=ms_energy+p%energy*p%grain_vol
    endif
    if (associated(p%uhx)) then; q=>field(i-1)
      do n=p%min_n,p%max_n
        if (p%do_grain(n)) p%h(:,n)=p%h(:,n)+get_interlayer_hx(n,p,q,p%uhxgrainindex,p%uhxgrain,p%uhx,associated(ms_energy))
      enddo
    endif
    if (associated(p%lhx)) then; q=>field(i+1)
      do n=p%min_n,p%max_n
        if (p%do_grain(n)) p%h(:,n)=p%h(:,n)+get_interlayer_hx(n,p,q,p%lhxgrainindex,p%lhxgrain,p%lhx)
      enddo
    endif
! n=3;if (field_mode.eq.1.and.present(h_var)) print "(i3,x,i3,1p,:50(x,e12.5))",n,ng,p%energy(n)*p%grain_vol(n)/(1.380662d-16*300.d0),p%h_app(:,n),p%h_demag(:,n),p%h(:,n),p%h_selfdemag(:,n)*var
!call plot_contour('total hx',p%h(1,:),.true.,i,no_write=.true.)
!call plot_contour('total hy',p%h(2,:),.true.,i,no_write=.true.)
!if (p_p) call plot_contour('total hz',p%h(3,:),.true.,i,no_write=.true.)
  enddo
  if (present(h_var)) h_var=0.d0
  if (thermal) then
     var=get_thermal(nm_p,diffm); if (present(h_var)) h_var=var
! else
!    do i=1,size(field); field(i)%hthm=0.d0; enddo
  elseif (diffm) then
    do i=1,size(field); p=>field(i); do n=p%min_n,p%max_n
      if (p%do_grain(n)) p%hthm(:,n)=p%h(:,n)*p%hthm_m(n)/sqrt(1.d-22+dot_product(p%h(:,n),p%h(:,n)))
    enddo; enddo
  endif
  if (present(varying_p)) varying_p=v_p


!if (field_mode.ne.0) then
!n=643
!print "(a,1p,:,50(x,e12.5))",'field:',i,ms_energy(n)
!do i=1,size(field);p=>field(i); id=p%in
!ng=field_mode; field_mode=0
!print "(i,1p,:,50(x,e12.5))",i,p%h_ext(:,n),p%h_app(:,n),p%h_demag(:,n),p%h_selfdemag(:,n),get_uni_hk(n,p,id),get_intralayer_hx(n,p,id),p%energy(n),p%h(:,n)
!if (associated(p%uhx)) then; q=>field(i-1); print "(i,1p,:,50(x,e12.5))",i,get_interlayer_hx(n,p,q,p%uhxgrainindex,p%uhxgrain,p%uhx);endif
!if (associated(p%lhx)) then; q=>field(i+1); print "(i,1p,:,50(x,e12.5))",i,get_interlayer_hx(n,p,q,p%lhxgrainindex,p%lhxgrain,p%lhx);endif
!field_mode=ng
!enddo
!if (field_mode.eq.-1) read *
!endif

!open(68,file='htol.dat',status='unknown')
!open(67,file='grain_happ.dat',status='unknown')
!open(66,file='grain_hdemag1.dat',status='unknown')
!open(65,file='hcubic.dat',status='unknown')
!open(64,file='huni.dat',status='unknown')
!open(63,file='hx.dat',status='unknown')
!open(62,file='hxl.dat',status='unknown')
!open(61,file='hxu.dat',status='unknown')
!do i=1,size(field);p=>field(i);ng=size(p%temp); id=p%in; do n=1,size(p%temp); 
!write(68,"(i3,1p,3(x,e12.5))")n,p%h(:,n);
!write(67,"(i3,1p,3(x,e12.5))")n,p%h_app(:,n);
!write(66,"(i3,1p,3(x,e12.5))")n,p%h_demag(:,n);
!write(65,"(i3,1p,3(x,e12.5))")n,get_cubic_hk(n,p,id);
!write(64,"(i3,1p,3(x,e12.5))")n,get_uni_hk(n,p,id);
!write(63,"(i3,1p,3(x,e12.5))")n,get_intralayer_hx(n,p,id)
!if (associated(p%uhx)) then;q=>field(i-1); write(61,"(i3,x,i3,1p,3(x,e12.5))")i,n,get_interlayer_hx(n,p,q,p%uhxgrainindex,p%uhxgrain,p%uhx);endif
!if (associated(p%lhx)) then;q=>field(i+1); write(62,"(i3,x,i3,1p,3(x,e12.5))")i,n,get_interlayer_hx(n,p,q,p%lhxgrainindex,p%lhxgrain,p%lhx);endif
!enddo;enddo
!close(68); close(67); close(66); close(65); close(64); close(63); close(62); close(61)
!read *

 end subroutine

 subroutine get_selfdemag()
  integer::il,n,is
  real(8)::mm(3)
  type(field_ss),pointer::p

  if (compute_selfdemag_p) call find_grain_N(); compute_selfdemag_p=.false.

! if ((selfdemag_p.or.associated(ms_energy)).and.size(field).gt.1) then
  if ((grain_selfdemag_p.or.associated(ms_energy)).and.size(field).gt.1) then
! if (associated(ms_energy).and.size(field).gt.1) then
!print *,'in get_selfdemag ms_energy',count(field(1)%do_grain),field_mode
    do is=1,size(field); field(is)%h_selfdemag=0.d0; enddo
    do is=1,size(field); p=>field(is)
      do n=p%min_n,p%max_n
        if (p%do_grain(n).or.field_mode.eq.1) then
          mm = p%ms(n)*p%ms_scale(n)*p%m(1+p%in:3+p%in,n)
          do il=1,size(field)
           if (field(il)%do_grain(n)) &
            field(il)%h_selfdemag(:,n) = field(il)%h_selfdemag(:,n)+&
                              (/ p%Ndemag(1,1,il,n)*mm(1) + p%Ndemag(2,1,il,n)*mm(2) + p%Ndemag(3,1,il,n)*mm(3), &
                                 p%Ndemag(1,2,il,n)*mm(1) + p%Ndemag(2,2,il,n)*mm(2) + p%Ndemag(3,2,il,n)*mm(3), &
                                 p%Ndemag(1,3,il,n)*mm(1) + p%Ndemag(2,3,il,n)*mm(2) + p%Ndemag(3,3,il,n)*mm(3) /)
          enddo
        endif
      enddo
    enddo
  else
    do il=1,size(field); p=>field(il)
      do n=p%min_n,p%max_n
        if (p%do_grain(n)) then
          mm = p%ms(n)*p%ms_scale(n)*p%m(1+p%in:3+p%in,n)
          p%h_selfdemag(:,n)= (/ p%Ndemag(1,1,il,n)*mm(1) + p%Ndemag(2,1,il,n)*mm(2) + p%Ndemag(3,1,il,n)*mm(3), &
                                 p%Ndemag(1,2,il,n)*mm(1) + p%Ndemag(2,2,il,n)*mm(2) + p%Ndemag(3,2,il,n)*mm(3), &
                                 p%Ndemag(1,3,il,n)*mm(1) + p%Ndemag(2,3,il,n)*mm(2) + p%Ndemag(3,3,il,n)*mm(3) /)
!if (n.eq.1) print *,mm,hh,p%Ndemag(:,:,:,n),p%h_selfdemag(:,n)
        else
          p%h_selfdemag(:,n)=0.d0
        endif
      enddo
!     if (field_mode.eq.-1) p%h_demag=p%h_demag-p%h_selfdemag
    enddo
  endif
  if (field_mode.eq.-1) then
    do il=1,size(field); field(il)%h_ext=field(il)%h_demag-field(il)%h_selfdemag; enddo
  elseif (field_mode.eq.1) then
    do il=1,size(field); field(il)%h_demag=field(il)%h_ext+field(il)%h_selfdemag; enddo
  else
    do il=1,size(field); field(il)%h_demag=field(il)%h_selfdemag; enddo
  endif
 end subroutine

 function get_interlayer_hx(n,p,q,hxgrainindex,hxgrain,hx,ms_p) result (h)
  integer,intent(in)::n,hxgrainindex(0:),hxgrain(:)
  real(8),intent(in)::hx(:)
  logical,intent(in),optional::ms_p

  integer::i
  real(8)::h(3),yy,xx
  type(field_ss),pointer::p,q
  h=0.d0; yy=0.d0
  do i=hxgrainindex(n-1)+1, hxgrainindex(n)
    h=h+hx(i)*q%ms_scale(hxgrain(i))*q%m(1+q%in:3+q%in,hxgrain(i))
  enddo
  if (field_mode.eq.-1.and..not.associated(ms_energy)) p%h_ext(:,n)=p%h_ext(:,n)+h
  xx=-dot_product(h,p%m(1+p%in:3+p%in,n))*p%ms_scale(n)*p%ms(n)
! if (present(ms_p).and.ms_p) then
  if (present(ms_p)) then; if (ms_p) then
    xx=xx*p%grain_vol(n); ms_energy(n)=ms_energy(n)+xx
    xx=xx/(q%grain_vol(n)+p%grain_vol(n))
    q%energy(n)=q%energy(n)+xx !FIXME? assumes columnar grain
! endif
  endif; endif
  p%energy(n)=p%energy(n)+xx
 end function

 function get_intralayer_hx(n,p,id) result (h)
  integer,intent(in)::id,n
  type(field_ss),pointer::p
  real(8)::h(3),yy
  integer::j
!if (lbound(p%hxgrainindex,dim=1).ne.0) then; print *,'get_intralayer';read *;endif
  h=0.d0
  do j=p%hxgrainindex(n-1)+1, p%hxgrainindex(n)
    h=h+p%hx(j)*p%ms_scale(p%hxgrain(j))*p%m(1+id:3+id,p%hxgrain(j))
  enddo
  yy=-dot_product(h,p%m(1+id:3+id,n))*p%ms(n)*p%ms_scale(n)
  p%energy(n)=p%energy(n)+yy;if (associated(ms_energy)) ms_energy(n)=ms_energy(n)+yy
  if (field_mode.eq.-1) p%h_ext(:,n)=p%h_ext(:,n)+h
 end function


 function get_uni_hk(n,p,id) result (h)
  integer,intent(in)::id,n
  type(field_ss),pointer::p
  real(8)::h(3),hks,yy
  h=0.d0; if (.not.associated(p%uni_axis).or.p%hk_scale(n).le.0.d0) return
! hks = p%m(1+id,n)*p%uni_axis(1,n)+p%m(2+id,n)*p%uni_axis(2,n)+p%m(3+id,n)*p%uni_axis(3,n)
  hks = dot_product(p%m(1+id:3+id,n),p%uni_axis(1:3,n)); yy=1.d0-hks*hks
  h = ( hks * (p%hk(1,n)+p%hk(2,n)*yy) * p%hk_scale(n) ) * p%uni_axis(:,n)
  p%energy(n)=p%energy(n)+0.5d0*yy*(p%hk(1,n)+0.5d0*p%hk(2,n)*yy)*p%hk_scale(n)*p%ms(n)*p%ms_scale(n)
 end function

 function get_cubic_hk(n,p,id) result (h)
  integer,intent(in)::id,n
  type(field_ss),pointer::p
  real(8)::h(3),dcs(3,3)
  h=0.d0; if (.not.associated(p%cubic_axis).or.p%hk_scale(n).le.0.d0) return
  dcs(:,1)=(/ p%cubic_axis(1,1,n)*p%m(1+id,n)+p%cubic_axis(2,1,n)*p%m(2+id,n)+p%cubic_axis(3,1,n)*p%m(3+id,n), &
              p%cubic_axis(1,2,n)*p%m(1+id,n)+p%cubic_axis(2,2,n)*p%m(2+id,n)+p%cubic_axis(3,2,n)*p%m(3+id,n), &
              p%cubic_axis(1,3,n)*p%m(1+id,n)+p%cubic_axis(2,3,n)*p%m(2+id,n)+p%cubic_axis(3,3,n)*p%m(3+id,n) /)
  dcs(:,2) = dcs(:,1)*dcs(:,1)
  dcs(:,3) = (/ dcs(1,1)*(p%hk(1,n)*(dcs(2,2)+dcs(3,2))+p%hk(2,n)*dcs(2,2)*dcs(3,2)), &
                dcs(2,1)*(p%hk(1,n)*(dcs(1,2)+dcs(3,2))+p%hk(2,n)*dcs(1,2)*dcs(3,2)), &
                dcs(3,1)*(p%hk(1,n)*(dcs(1,2)+dcs(2,2))+p%hk(2,n)*dcs(1,2)*dcs(2,2)) /)
  h = - (/ p%cubic_axis(1,1,n)*dcs(1,3)+p%cubic_axis(1,2,n)*dcs(2,3)+p%cubic_axis(1,3,n)*dcs(3,3) , &
           p%cubic_axis(2,1,n)*dcs(1,3)+p%cubic_axis(2,2,n)*dcs(2,3)+p%cubic_axis(2,3,n)*dcs(3,3) , &
           p%cubic_axis(3,1,n)*dcs(1,3)+p%cubic_axis(3,2,n)*dcs(2,3)+p%cubic_axis(3,3,n)*dcs(3,3) /) * p%hk_scale(n)
  p%energy(n)=p%energy(n)+0.5d0*p%hk_scale(n)*p%ms(n)* &
              (p%hk(1,n)*(dcs(1,2)*dcs(2,2)+dcs(1,2)*dcs(3,2)+dcs(2,2)*dcs(3,2)) + &
               p%hk(2,n)* dcs(1,2)*dcs(2,2)*dcs(3,2))*p%hk_scale(n)*p%ms(n)*p%ms_scale(n)
 end function

 function get_thermal(no_ms_p,diff_m) result (var)
  use random_m, only:random_normal3d,random_normal
  logical,intent(in)::no_ms_p,diff_m
  real(8)::var,varl,tsigma,x

  integer::i,n,n0

  var=0.d0
  do i=1,size(field)
    if (field(i)%thermal_p) then
      tsigma=sqrt(2.d0*1.380662d-16*field(i)%alpha/field(i)%gamma)  !2k alpha / gamma
        !var = 2 k T alpha / (gamma Ms V)
      varl=0.d0;n0=0
!     do n=1,size(field(i)%temp)
      do n=field(i)%min_n,field(i)%max_n
       if (field(i)%do_grain(n)) then
        n0=n0+1
        if (no_ms_p) then
          x=field(i)%alpha_scale(n)*field(i)%temp(n)/(field(i)%ms(n)*field(i)%grain_vol(n))
        elseif (field(i)%mag(field(i)%in/3+1,n).ne.1.d0) then
          x=field(i)%alpha_scale(n)*field(i)%temp(n)/(field(i)%ms(n)*field(i)%grain_vol(n)*max(field(i)%floor_ms,abs(field(i)%mag(field(i)%in/3+1,n))))
        elseif (field(i)%ms_scale(n).gt.0.d0) then
          x=field(i)%alpha_scale(n)*field(i)%temp(n)/(field(i)%ms(n)*field(i)%grain_vol(n)*max(field(i)%floor_ms,field(i)%ms_scale(n)))
        else
          x=0.d0
        endif
        x=sqrt(x)*tsigma; varl=varl+x
        if (diff_m) then
          field(i)%hthm_m(n)=x*random_normal(field(i)%random_thml)
          field(i)%hthm(:,n)=field(i)%h(:,n)*field(i)%hthm_m(n)/sqrt(1.d-100+dot_product(field(i)%h(:,n),field(i)%h(:,n)))
        else
          field(i)%hthm(:,n)=x*random_normal3d(field(i)%random_thml)
        endif
       endif
      enddo
      var=var+varl/max(1,n0)
    endif
  enddo
  var=var/max(1,i-1)
 end function

 subroutine get_magnetostatic()
  use demag_m, only: get_demag, demag_window
  use window_m, only: get_window
! real(8)::offset(2)

  type(field_ss),pointer::p
  integer::nx,ny,nxf,nyf,n,n1,i,j,k,i0(2)
  real(8)::y(3)

    call get_demag(demag_xs,demag_ys,demag_zs); if (.not.associated(demag_xs)) return
    nx=size(app_temp,1); ny=size(app_temp,2); 
    nxf=size(demag_xs,1); nyf=size(demag_xs,2);
    if (nx.eq.nxf.and.ny.eq.nyf) then
      do k=1,size(field); p=>field(k)
!if (lbound(p%grain2gridindex,dim=1).ne.0) then; print *,'get_magnetostatic';read *;endif
!       p%h_demag=0.0d0
!       do j=1,ny; do i=1,nx
!         n=p%grid2grain(i,j)
!         if (n.gt.0) p%h_demag(:,n)=p%h_demag(:,n)+(/ dble(demag_xs(i,j,k)), dble(demag_ys(i,j,k)), dble(demag_zs(i,j,k)) /)
!       enddo; enddo
!       do n=1,size(p%temp)
!         p%h_demag(:,n)=p%h_demag(:,n)/(p%grain2gridindex(n)-p%grain2gridindex(n-1))
!       enddo
!$omp parallel
!$omp do private(n,y,n1,i,j) firstprivate(p,k)
        do n=1,size(p%temp); y=0.d0
          do n1=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); i=p%grain2grid(1,n1); j=p%grain2grid(2,n1)
            y=y+(/ dble(demag_xs(i,j,k)), dble(demag_ys(i,j,k)), dble(demag_zs(i,j,k)) /)
          enddo
          p%h_demag(:,n)=y/(p%grain2gridindex(n)-p%grain2gridindex(n-1))
        enddo
!$omp end do nowait
!$omp end parallel
!  call plot_contour('demag_full',dble(demag_zs(:,:,k)),.false.,k,no_write=.true.)
!  call plot_contour('demag_grain',p%h_demag(3,:),.true.,k,no_write=.true.)
!  call plot_contour('mz_grain',p%m(3+p%in,:),.true.,k,no_write=.true.)
      enddo
    else
!     i0=max(1,min(nx-nxf+1,nint((nx-nxf)*0.5d0+offset(1))))-1
!     j0=max(1,min(ny-nyf+1,nint((ny-nyf)*0.5d0+offset(2))))-1
!     call demag_window(i0=i0,j0=j0)
      call get_window(i0=i0)
      do k=1,size(field); p=>field(k)
!if (lbound(p%grain2gridindex,dim=1).ne.0) then; print *,'get_magnetostatic';read *;endif
!       p%h_demag=0.0d0
!       do j=1,nyf; do i=1,nxf
!         n=p%grid2grain(i+i0,j+j0)
!         if (n.gt.0) p%h_demag(:,n)=p%h_demag(:,n)+(/ dble(demag_xs(i,j,k)), dble(demag_ys(i,j,k)), dble(demag_zs(i,j,k)) /)
!       enddo; enddo;
!       do n=1,size(p%temp)
!         p%h_demag(:,n)=p%h_demag(:,n)/(p%grain2gridindex(n)-p%grain2gridindex(n-1))
!       enddo
!$omp parallel
!$omp do private(n,y,n1,i,j) firstprivate(p,k)
!       do n=1,size(p%temp); y=0.d0
        do n=p%min_n,p%max_n; y=0.d0
          if (p%do_grain(n)) then
          do n1=p%grain2gridindex(n-1)+1,p%grain2gridindex(n); i=p%grain2grid(1,n1)-i0(1); j=p%grain2grid(2,n1)-i0(2)
            if (min(i,j).gt.0.and.i.le.nxf.and.j.le.nyf) y=y+(/ dble(demag_xs(i,j,k)), dble(demag_ys(i,j,k)), dble(demag_zs(i,j,k)) /)
          enddo
          endif
          p%h_demag(:,n)=y/(p%grain2gridindex(n)-p%grain2gridindex(n-1))
!write(67,"(i3,1p,3(x,e12.5))") n,p%h_demag(:,n)
        enddo
!$omp end do nowait
!$omp end parallel
      enddo
!open(67,file='grain_hdemag.dat',status='unknown')
!do k=1,size(field); p=>field(k); do n=1,size(p%temp); write(67,"(i4,1p,3(x,e12.5))") n,p%h_demag(:,n); enddo; enddo
!close(67)
    endif
!read *

! k=1
! write(os,"('x Hmag for layer',i3)")k; call plot_contour(trim(os),field(k)%h_demag(1,:),.false.,k,no_write=.true.)
! write(os,"('y Hmag for layer',i3)")k; call plot_contour(trim(os),field(k)%h_demag(2,:),.false.,k,no_write=.true.)
! write(os,"('perp Hmag for layer',i3)")k; call plot_contour(trim(os),field(k)%h_demag(3,:),.false.,k,no_write=.true.)
! read *

 end subroutine

 subroutine freeze_time(t)
  real(8),intent(in)::t
  frozen_time=t
 end subroutine
 subroutine unfreeze_time()
  frozen_time=-1.d300
 end subroutine

 function get_applied(t,silent_p,plot_p,rot_phi) result (ok_p)
  use applied_m, only: get_appl, get_applied_scale
  use temp_scaling_m, only: find_temp_scaling
  use plot_m, only: get_plot, plot_contour, get_plot, plot_wf_scaling, plot_temp_scaling, plot_mk, plot_heff, plot_strip
  real(8),intent(in)::t
  logical,intent(in),optional::silent_p,plot_p
  real(8),intent(in),optional::rot_phi(2)

  logical::ok_p(2),silent,p_p
  integer::n,k,i
  real(8)::x,y(3)
  real(8),dimension(:,:),pointer::tempscale,fldscale
  character(200)::os

  silent=.false.; if (present(silent_p)) silent=silent_p
  p_p=.false.; if (present(plot_p)) p_p=plot_p
  ok_p(1)=get_appl(t,changed,silent,rot_phi)

  do k=1,size(field)
!if (lbound(field(k)%grain2gridindex,dim=1).ne.0) then; print *,'get_applied';read *;endif
    if (changed(k,1).and.changed(k,2)) then
      do n=field(k)%min_n,field(k)%max_n; x=0.d0; y=0.d0
        if (field(k)%do_grain(n)) then
          do i=field(k)%grain2gridindex(n-1)+1,field(k)%grain2gridindex(n)
            x=x+app_temp(field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
            y=y+app_field(:,field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
          enddo
        endif
        field(k)%temp_app(n)=x/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+utemp(k)
        field(k)%h_app(:,n)=y/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+ufield(:,k)
      enddo
      if (field(k)%source_p) then
        call get_grain_temp(t,k)
      else
        field(k)%temp=field(k)%temp_app+field(k)%background_t
      endif
!print *,'get_applied, calling temp_scaling',k,minval(field(k)%ms_scale)
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc)
    elseif (changed(k,1)) then
      do n=field(k)%min_n,field(k)%max_n; x=0.d0;
        if (field(k)%do_grain(n)) then
          do i=field(k)%grain2gridindex(n-1)+1,field(k)%grain2gridindex(n)
            x=x+app_temp(field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
          enddo
        endif
        field(k)%temp_app(n)=x/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+utemp(k)
      enddo
      if (field(k)%source_p) then
        call get_grain_temp(t,k)
      else
        field(k)%temp=field(k)%temp_app+field(k)%background_t
      endif
!print *,'get_applied, calling temp_scaling',k,minval(field(k)%ms_scale)
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc)
    elseif (changed(k,2)) then
      do n=field(k)%min_n,field(k)%max_n; y=0.d0;
        if (field(k)%do_grain(n)) then
          do i=field(k)%grain2gridindex(n-1)+1,field(k)%grain2gridindex(n)
            y=y+app_field(:,field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
          enddo
        endif
        field(k)%h_app(:,n)=y/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+ufield(:,k)
      enddo
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc,.false.)
    else
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc,.false.)
    endif
  enddo
!print *,changed(:,2),ufield(:,:),field(1)%h_app(:,1),app_field(:,1,1,1)
!do k=1,size(field); call plot_strip('Hu',k,t*1.d9,(/sqrt(dot_product(ufield(:,k),ufield(:,k))),ufield(:,k)/)); enddo
  if (p_p.and.get_plot()) then
    call get_applied_scale(tempscale,fldscale)
    do k=1,size(field)
      if (get_plot(k,temp_scale=.true.)) call plot_temp_scaling(k,t*1.d9,tempscale(:,k))
      if (get_plot(k,wf_scale=.true.)) call plot_wf_scaling(k,t*1.d9,fldscale(:,k))
!      if (get_plot(k,wf_scale=.true.)) call plot_strip('Hu',k,t*1.d9,(/sqrt(dot_product(ufield(:,k),ufield(:,k))),ufield(:,k)/))
      if (get_plot(k,wf_eff=.true.)) then; write(os,"('Hsw for layer',i3)")k; call plot_heff(trim(os),field(k)%h_app(:,:),k); endif
      if (get_plot(k,temp=.true.)) then; write(os,"('temperature for layer',i3)")k; call plot_contour(trim(os),field(k)%temp,.true.,k); endif
      if (associated(app_temp).and.get_plot(k,tempf=.true.)) then
              write(os,"('full temperature for layer',i3)")k; call plot_contour(trim(os),app_temp(:,:,k)+field(k)%background_t,.true.,k); endif
      if (get_plot(k,wf_x=.true.)) then; write(os,"('Happ,x for layer',i3)")k; call plot_contour(trim(os),field(k)%h_app(1,:),.true.,k); endif
      if (get_plot(k,wf_y=.true.)) then; write(os,"('Happ,y for layer',i3)")k; call plot_contour(trim(os),field(k)%h_app(2,:),.true.,k); endif
      if (get_plot(k,wf_z=.true.)) then; write(os,"('Happ,z for layer',i3)")k; call plot_contour(trim(os),field(k)%h_app(3,:),.true.,k); endif
      if (get_plot(k,wf_xf=.true.)) then; write(os,"('full Happ,x for layer',i3)")k; call plot_contour(trim(os),app_field(1,:,:,k)+ufield(1,k),.true.,k); endif
      if (get_plot(k,wf_yf=.true.)) then; write(os,"('full Happ,y for layer',i3)")k; call plot_contour(trim(os),app_field(2,:,:,k)+ufield(2,k),.true.,k); endif
      if (get_plot(k,wf_zf=.true.)) then; write(os,"('full Happ,z for layer',i3)")k; call plot_contour(trim(os),app_field(3,:,:,k)+ufield(3,k),.true.,k); endif
! call plot_contour('mz_grain',field(k)%m(3+field(k)%in,:),.true.,k)
    enddo
  endif

  ok_p(1)=any(changed(:,1))
  ok_p(2)=any(changed(:,2))

 end function

 subroutine get_grain_temp(t,k)
  real(8),intent(in)::t
  integer,intent(in)::k

  type(field_ss),pointer::p
  integer::n
  real(8)::dt


  p=>field(k)

  dt=max(0.d0,t-p%last_t); p%last_t=t
!print *,t,k,p%avetau,dt,p%tau(1),p%temp_app(1),p%background_t;read *
! do n=1,size(p%ms)
  do n=p%min_n,p%max_n
    if (p%do_grain(n)) then
      p%temp(n)=p%temp_app(n)*p%tau(n)/p%avetau+p%tau_hs_t+p%background_t+(p%temp(n)-p%temp_app(n)*p%tau(n)/p%avetau-p%tau_hs_t-p%background_t)*exp(-dt/p%tau(n))
    else
      p%temp(n)=p%tau_hs_t+p%background_t
    endif
  enddo
 end subroutine

 function init_applied(t,silent_p,plot_p) result (ok_p)
  use applied_m, only: share_applied, init_appl
  use demag_m, only: demag_window, init_demag
  use temp_scaling_m, only: find_temp_scaling
  real(8),intent(in)::t
  logical,intent(in),optional::silent_p,plot_p

  logical::ok_p,silent
  integer::n,k,i
  real(8)::vel(2),x,y(3)

  call update_data(0)

  ok_p=.false.
  silent=.false.; if (present(silent_p)) silent=silent_p
  if (.not.allocated(changed)) then
    allocate(changed(size(field),2))
    call share_applied(app_temp,app_field,utemp,ufield)
  end if

  if (frozen_time.gt.-1d20) then
    ok_p=init_appl(frozen_time,changed,silent,vel)
  else
    ok_p=init_appl(t,changed,silent,vel)
  endif
  call demag_window(v=vel)

  do k=1,size(field)
!if (lbound(field(k)%grain2gridindex,dim=1).ne.0) then; print *,'init_applied';read *;endif
    if (changed(k,1)) then
      do n=field(k)%min_n,field(k)%max_n; x=0.d0
        if (field(k)%do_grain(n)) then
          do i=field(k)%grain2gridindex(n-1)+1,field(k)%grain2gridindex(n)
            x=x+app_temp(field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
          enddo
        endif
        field(k)%temp_app(n)=x/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+utemp(k)
      enddo
      if (field(k)%source_p) then
        call get_grain_temp(t,k)
      else
        field(k)%temp=field(k)%temp_app+field(k)%background_t
      endif
!print *,'init_applied, calling temp_scaling',k,minval(field(k)%ms_scale)
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc)
    else
      ok_p=find_temp_scaling(k,field(k)%temp,field(k)%tc,.false.)
    endif
    if (changed(k,2)) then
      do n=field(k)%min_n,field(k)%max_n; y=0.d0
        if (field(k)%do_grain(n)) then
          do i=field(k)%grain2gridindex(n-1)+1,field(k)%grain2gridindex(n)
            y=y+app_field(:,field(k)%grain2grid(1,i),field(k)%grain2grid(2,i),k)
          enddo
        endif
        field(k)%h_app(:,n)=y/(field(k)%grain2gridindex(n)-field(k)%grain2gridindex(n-1))+ufield(:,k)
      enddo
    endif
    field(k)%hthm=0.d0; field(k)%hthm_m=0.d0
  enddo
!print *,changed(:,2),ufield(:,:),frozen_time,vel,field(1)%h_app(:,1),app_field(:,1,1,1)

  ok_p=init_demag(selfdemag_p)  !has the demag tensor been calculated? make sure it has
! if (compute_selfdemag_p) call find_grain_N(); compute_selfdemag_p=.false.

 end function

 function start_field(seed,self_fld_only_p,self_fld_grain_only_p) result(ok_p)
   use random_m, only: random_create
   integer,intent(in),optional::seed(:)
   logical,intent(in),optional::self_fld_only_p,self_fld_grain_only_p
   
   integer::i
   logical::ok_p
   character(200)::str

   ok_p=.false.

   if (present(seed)) then
     if (.not.allocated(field)) allocate(field(size(seed)))
     do i=1,size(seed)
       write(str,"(i5)") i
       field(i)%random_thml=>random_create('random_thml'//trim(adjustl(str)),seed(i))
     enddo
   endif
   if (present(self_fld_only_p)) selfdemag_p=self_fld_only_p
   if (present(self_fld_grain_only_p)) grain_selfdemag_p=self_fld_grain_only_p
   ok_p=.true.
 end function

 subroutine update_data(il)
  use data_m, only:get_data
  integer,intent(in),optional::il
  integer::i,k1,k2,j
  type(field_ss),pointer::p
  real(8)::w

  k1=1; k2=size(field)
  if (present(il)) then
    if (il.gt.0.and.il.le.k2) then; k1=il; k2=il; endif
  endif
  do i=k1,k2
    p=>field(i)
    call get_data(i,p%grid2grain,grid2grain_p=.true.)
    call get_data(i,p%grain2grid,grain2grid_p=.true.)
    call get_data(i,p%grain2gridindex,grain2gridindex_p=.true.)
    call get_data(i,p%grain_vol,grain_vol_p=.true.)
    call get_data(i,p%grain_area,grain_area_p=.true.)
    call get_data(i,p%do_grain,do_grain_p=.true.)
    call get_data(i,p%min_n,min_n_p=.true.)
    call get_data(i,p%max_n,max_n_p=.true.)
    call get_data(i,p%tc,tc_p=.true.)
    call get_data(i,p%tau,tau_p=.true.)
    call get_data(i,p%ms,ms_p=.true.)
    call get_data(i,p%mag,mag_p=.true.)
    call get_data(i,p%energy,energy_p=.true.)
    call get_data(i,p%temp,temp_p=.true.)
    call get_data(i,p%ms_scale,ms_scale_p=.true.)
    call get_data(i,p%alpha_scale,alpha_scale_p=.true.)
    call get_data(i,p%hk_scale,hk_scale_p=.true.)
    call get_data(i,p%hx_scale,hx_scale_p=.true.)
    call get_data(i,p%dz,dz_p=.true.)
    call get_data(i,p%in,in_p=.true.)
    call get_data(i,p%floor_ms,floor_ms_p=.true.)
    call get_data(i,p%alpha,thm_alpha_p=.true.)
    call get_data(i,p%gamma,gamma_p=.true.)
    call get_data(i,p%m,m_p=.true.)
    call get_data(i,p%h,h_p=.true.)
    call get_data(i,p%hthm,hthm_p=.true.)
    call get_data(i,p%hthm_m,hthm_p=.true.)
    call get_data(i,p%background_t,background_t_p=.true.)
    call get_data(i,p%thermal_p,thermal_p=.true.)
    call get_data(i,p%tau_hs_t,tau_hs_t_p=.true.)
!   p%avetau=0.d0; do j=1,size(p%grain_vol); p%avetau=p%avetau+p%tau(j)*p%grain_vol(j); enddo; p%avetau=p%avetau/sum(p%grain_vol)
    p%avetau=0.d0; w=0.d0; do j=1,size(p%grain_vol); w=w+p%grain_vol(j); p%avetau=p%avetau+(p%tau(j)-p%avetau)*p%grain_vol(j)/w; enddo

!if (lbound(p%grain2gridindex,dim=1).ne.0) then; print *,'update_data';read *;endif
  enddo
  call get_data(1,ms_energy,ms_energy_p=.true.)

 end subroutine

 function construct_field(nx,ny,dx,dy,i,ihk_p,ihx_p) result (ok_p)
   use construct_grain_m, only: construct_hk_magnitude, construct_hk_direction, construct_hx_transverse, construct_hx_vertical

   integer,intent(in)::i,nx,ny
   real(8),intent(in)::dx,dy
   logical,intent(in),optional::ihk_p,ihx_p

   type(field_ss),pointer::p,q
   logical::ok_p,hk_p,hx_p
   integer::n

   hk_p=.true.; hx_p=.true.
   if (present(ihk_p)) hk_p=ihk_p
   if (present(ihx_p)) hx_p=ihx_p

   ok_p=.true.

   p=>field(i)
     call update_data(i)
     p%last_t=-1.d22


!if (lbound(p%grain2gridindex,dim=1).ne.0) then; print *,'construct_field';read *;endif
     n=size(p%tc)

     if (associated(p%h_app)) then; if (size(p%h_app,2).ne.n) then; 
       deallocate(p%h_demag,p%h_app,p%temp_app,p%h_selfdemag,p%h_plot,p%h_ext); 
       nullify(p%h_demag,p%h_app,p%temp_app,p%h_selfdemag,p%h_plot,p%h_ext); 
     endif; endif
     if (.not.associated(p%h_app)) then
       allocate(p%h_demag(3,n),p%h_app(3,n),p%temp_app(n),p%h_selfdemag(3,n),p%h_plot(3,n),p%h_ext(3,n)); 
       p%h_demag=0.d0; p%h_app=0.d0; p%temp_app=0.d0; p%h_selfdemag=0.d0; p%h_plot=0.d0; p%h_ext=0.d0
     endif

     if (hk_p) then
       ok_p=construct_hk_magnitude(n,p%grain_area,p%ms,p%dz,i,.false.,p%hk)
       if (ok_p) ok_p=construct_hk_direction(n,i,.false.,p%uni_axis,p%cubic_axis)
     endif
     if (hx_p) then
       if (ok_p) ok_p=construct_hx_transverse(nx,ny,dx,dy,n,p%grain_area,p%ms, &
            p%grid2grain,p%grain2gridindex,p%grain2grid,i,.false.,p%hx,p%hxgrainindex,p%hxgrain)
       if (ok_p.and.i.gt.1) then
         q=>field(i-1)
         ok_p=construct_hx_vertical(nx,ny,dx,dy,i,.false., size(q%tc),&
              q%grain_area,q%ms,q%grid2grain,q%grain2gridindex,q%dz,q%grain2grid,q%lhx,q%lhxgrainindex,q%lhxgrain, &
            n,p%grain_area,p%ms,p%grid2grain,p%grain2gridindex,p%dz,p%grain2grid,p%uhx,p%uhxgrainindex,p%uhxgrain)
       endif
     endif
     compute_selfdemag_p=.true.

 end function

 subroutine find_grain_N()
   use io_m, only: output
   use tensor_m, only: get_N_layer

   integer::nx,ny,i,n,nxf,nyf,k1,dj,di,ks,ko,idm,idm0
   type(field_ss),pointer::p
   real(8),pointer::N_layer(:,:,:,:,:)

   call output(' need to calculate grains demag tensor elements!')
   do i=1,size(field); p=>field(i)
     if (associated(p%Ndemag)) then; 
       if (size(p%Ndemag,4).ne.size(p%temp).or.size(p%Ndemag,3).ne.size(field)) then; deallocate(p%Ndemag); nullify(p%Ndemag); endif
     endif
     if (.not.associated(p%Ndemag)) allocate(p%Ndemag(3,3,size(field),size(p%temp))); p%Ndemag=0.d0
   enddo

   call get_N_layer(N_layer)
   if (.not.associated(N_layer)) return

   nx =size(p%grid2grain,1); ny =size(p%grid2grain,2); 
   nxf=size(N_layer,1);    nyf=size(N_layer,2);
   do i=1,size(field); p=>field(i)
!if (lbound(p%grain2gridindex,dim=1).ne.0) then; print *,'find_grain_N';read *;endif
     print "('   computing layer #',i3,' on other layers:',$)",i
     idm0=size(p%temp)/10; idm=0
     do n=1,size(p%temp)
       idm=idm+1; if (mod(idm-1,idm0).gt.mod(idm,idm0)) print "(' ',i2,$)",10-idm/idm0
       do ks=p%grain2gridindex(n-1)+1,p%grain2gridindex(n)
         do ko=p%grain2gridindex(n-1)+1,p%grain2gridindex(n)
           dj=p%grain2grid(2,ks)-p%grain2grid(2,ko); if (abs(dj).gt.ny/2) dj=dj-sign(ny,dj)
           if (dj.gt.0) then; dj=nyf+1-dj; else; dj=1-dj; endif
           di=p%grain2grid(1,ks)-p%grain2grid(1,ko); if (abs(di).gt.nx/2) di=di-sign(nx,di)
           if (di.gt.0) then; di=nxf+1-di; else; di=1-di; endif
           do k1=1,size(field)
             if (di.gt.0.and.di.le.nxf.and.dj.gt.0.and.dj.le.nyf) then
               p%Ndemag(:,1,k1,n) = p%Ndemag(:,1,k1,n) + N_layer(di,dj,1:3,k1,i)
               p%Ndemag(:,2,k1,n) = p%Ndemag(:,2,k1,n) + N_layer(di,dj,4:6,k1,i)
               p%Ndemag(:,3,k1,n) = p%Ndemag(:,3,k1,n) + N_layer(di,dj,7:9,k1,i)
             else
               print *,'oops',di,dj,ks,ko,n,nx,ny,nxf,nyf
               stop
             endif
           enddo
         enddo
       enddo
       p%Ndemag(:,:,:,n)=p%Ndemag(:,:,:,n)*nx*ny/(p%grain2gridindex(n)-p%grain2gridindex(n-1))
!p%Ndemag(1:2,1:2,i,n)=0.d0  !RW
!           field(il)%h_selfdemag(:,n) = field(il)%h_selfdemag(:,n)+&
!                             (/ p%Ndemag(1,1,il,n)*mm(1) + p%Ndemag(2,1,il,n)*mm(2) + p%Ndemag(3,1,il,n)*mm(3), &
!                                p%Ndemag(1,2,il,n)*mm(1) + p%Ndemag(2,2,il,n)*mm(2) + p%Ndemag(3,2,il,n)*mm(3), &
!                                p%Ndemag(1,3,il,n)*mm(1) + p%Ndemag(2,3,il,n)*mm(2) + p%Ndemag(3,3,il,n)*mm(3) /)
     enddo
     call output('  -done')
   enddo
!p%Ndemag(1:2,1:2,:,n)=0.d0  !RW
!print *,field(1)%Ndemag(:,:,1,1)
!print *,field(1)%Ndemag(:,:,1,size(field(1)%temp))
 end subroutine

 function file_field(read_p,ionum,layer,n,old) result (ok_p)
   use data_m, only: get_size
   logical,intent(in)::read_p,old
   integer,intent(in)::ionum,layer,n
   logical::ok_p,u_p,c_p,x_p,ux_p,lx_p
   integer::n2,nx,ny,n_lay
   real(8)::dx,dy
   type(field_ss),pointer::p

   integer::i1,i2
   integer,allocatable::gia(:)
   real(8)::g8,g81
   real(8),allocatable::g8a1(:),g8a2(:,:)
   
   ok_p=.false.; call get_size(nx,ny,n_lay,dx,dy)
   if (layer.lt.1.or.layer.gt.n_lay) return
   p=>field(layer)
   if (read_p) then
     ok_p=construct_field(nx,ny,dx,dy,layer,.false.,.false.)
     if (old) then
       if (associated(p%hk)) then; if (size(p%hk,2).ne.n) then; deallocate(p%hk); nullify(p%hk); endif; endif
       if (.not.associated(p%hk)) allocate(p%hk(2,n))
       if (associated(p%uni_axis)) then; if (size(p%uni_axis,2).ne.n) then; deallocate(p%uni_axis); nullify(p%uni_axis); endif; endif
       if (.not.associated(p%uni_axis)) allocate(p%uni_axis(3,n))
       if (associated(p%cubic_axis)) then; if (size(p%cubic_axis,3).ne.n) then; deallocate(p%cubic_axis); nullify(p%cubic_axis); endif; endif
       if (.not.associated(p%cubic_axis)) allocate(p%cubic_axis(3,3,n))
       allocate(gia(nx*ny),g8a1(n),g8a2(3,n))
       read(ionum,err=99) p%uni_axis,p%hk,gia,p%cubic_axis,g8a1,g8a2
       if (maxval(abs(p%uni_axis(:,:))).eq.0.d0) then; deallocate(p%uni_axis); nullify(p%uni_axis); endif
       if (maxval(abs(p%cubic_axis(:,1,:))).eq.0.d0) then; deallocate(p%cubic_axis); nullify(p%cubic_axis); endif
       read(ionum,err=99) g81,g8,g8,g8,i1,i2
       deallocate(gia,g8a1,g8a2)
       if (g81.ne.0.d0) read(ionum,err=99)
       if (i1.ne.0) read(ionum,err=99)
       if (i2.ne.0) read(ionum,err=99)
       ok_p=construct_field(nx,ny,dx,dy,layer,.false.,.true.)
     else
       if (associated(p%hk)) then; if (size(p%hk,2).ne.n) then; deallocate(p%hk); nullify(p%hk); endif; endif
       if (.not.associated(p%hk)) allocate(p%hk(2,n))
       read(ionum,err=99) u_p,c_p,x_p,lx_p,ux_p,p%hk
       if (u_p) then
         if (associated(p%uni_axis)) then; if (size(p%uni_axis,2).ne.n) then; deallocate(p%uni_axis); nullify(p%uni_axis); endif; endif
         if (.not.associated(p%uni_axis)) allocate(p%uni_axis(3,n))
         read(ionum,err=99) p%uni_axis
       else
         if (associated(p%uni_axis)) then; deallocate(p%uni_axis); nullify(p%uni_axis); endif
       endif
       if (c_p) then
         if (associated(p%cubic_axis)) then; if (size(p%cubic_axis,3).ne.n) then; deallocate(p%cubic_axis); nullify(p%cubic_axis); endif; endif
         if (.not.associated(p%cubic_axis)) allocate(p%cubic_axis(3,3,n))
         read(ionum,err=99) p%cubic_axis
       else
         if (associated(p%cubic_axis)) then; deallocate(p%cubic_axis); nullify(p%cubic_axis); endif
       endif
       if (x_p) then
         if (associated(p%hxgrainindex)) then; if (size(p%hxgrainindex).ne.n+1) then; deallocate(p%hxgrainindex); nullify(p%hxgrainindex); endif; endif
         if (.not.associated(p%hxgrainindex)) allocate(p%hxgrainindex(0:n))
         read(ionum,err=99) p%hxgrainindex,n2
         if (associated(p%hx)) then; if (size(p%hx).ne.n2) then; deallocate(p%hx,p%hxgrain); nullify(p%hx,p%hxgrain); endif; endif
         if (.not.associated(p%hx)) allocate(p%hxgrain(n2),p%hx(n2))
         read(ionum,err=99) p%hxgrain,p%hx
       else
         if (associated(p%hxgrainindex)) then; deallocate(p%hxgrainindex); nullify(p%hxgrainindex); endif
         if (associated(p%hx)) then; deallocate(p%hx,p%hxgrain); nullify(p%hx,p%hxgrain); endif
       endif
       if (lx_p) then
         if (associated(p%lhxgrainindex)) then; if (size(p%lhxgrainindex).ne.n+1) then; deallocate(p%lhxgrainindex); nullify(p%lhxgrainindex); endif; endif
         if (.not.associated(p%lhxgrainindex)) allocate(p%lhxgrainindex(0:n))
         read(ionum,err=99) p%lhxgrainindex,n2
         if (associated(p%lhx)) then; if (size(p%lhx).ne.n2) then; deallocate(p%lhx,p%lhxgrain); nullify(p%lhx,p%lhxgrain); endif; endif
         if (.not.associated(p%lhx)) allocate(p%lhxgrain(n2),p%lhx(n2))
         read(ionum,err=99) p%lhxgrain,p%lhx
       else
         if (associated(p%lhxgrainindex)) then; deallocate(p%lhxgrainindex); nullify(p%lhxgrainindex); endif
         if (associated(p%lhx)) then; deallocate(p%lhx,p%lhxgrain); nullify(p%hx,p%lhxgrain); endif
       endif
       if (ux_p) then
         if (associated(p%uhxgrainindex)) then; if (size(p%uhxgrainindex).ne.n+1) then; deallocate(p%uhxgrainindex); nullify(p%uhxgrainindex); endif; endif
         if (.not.associated(p%uhxgrainindex)) allocate(p%uhxgrainindex(0:n))
         read(ionum,err=99) p%uhxgrainindex,n2
         if (associated(p%uhx)) then; if (size(p%uhx).ne.n2) then; deallocate(p%uhx,p%uhxgrain); nullify(p%uhx,p%uhxgrain); endif; endif
         if (.not.associated(p%uhx)) allocate(p%uhxgrain(n2),p%uhx(n2))
         read(ionum,err=99) p%uhxgrain,p%uhx
       else
         if (associated(p%uhxgrainindex)) then; deallocate(p%uhxgrainindex); nullify(p%uhxgrainindex); endif
         if (associated(p%uhx)) then; deallocate(p%uhx,p%uhxgrain); nullify(p%hx,p%uhxgrain); endif
       endif
     endif
   else
     write(ionum,err=99) associated(p%uni_axis),associated(p%cubic_axis),associated(p%hx),associated(p%lhx),associated(p%uhx),p%hk
     if (associated(p%uni_axis)) write(ionum,err=99) p%uni_axis
     if (associated(p%cubic_axis)) write(ionum,err=99) p%cubic_axis
     if (associated(p%hx)) then
       write(ionum,err=99) p%hxgrainindex,size(p%hxgrain)
       write(ionum,err=99) p%hxgrain,p%hx
     endif
     if (associated(p%lhx)) then
       write(ionum,err=99) p%lhxgrainindex,size(p%lhxgrain)
       write(ionum,err=99) p%lhxgrain,p%lhx
     endif
     if (associated(p%uhx)) then
       write(ionum,err=99) p%uhxgrainindex,size(p%uhxgrain)
       write(ionum,err=99) p%uhxgrain,p%uhx
     endif
     ok_p=.true.
   endif
   return
 99 continue
   ok_p=.false.
   return
 end function

 subroutine get_easy_axis(il,data)
   integer,intent(in)::il
   real(8),pointer::data(:,:)

   if (il.lt.1.or.il.gt.size(field)) then
     return
   elseif (associated(field(il)%uni_axis)) then
     data=>field(il)%uni_axis
   else
     data=>field(il)%cubic_axis(:,1,:)  !return the a-axis direction...
   endif
 end subroutine

 subroutine plot_field(il,hd,hd_x,hd_y,hd_z,hx,hx_x,hx_y,hx_z,hk,hk_x,hk_y,hk_z,hv,hv_x,hv_y,hv_z,htot,htot_x,htot_y,htot_z,wf_eff, &
      wf,wf_x,wf_y,wf_z,ha,ha_x,ha_y,ha_z,ha_eff,temp,tempf,hb,hb_x,hb_y,hb_z,hb_ds,hb_r,hb_m,hb_tilt)
  use plot_m, only: get_plot, plot_contour, plot_heff
  use data_m, only: get_geometry
!  use mselem_m, only: block_to_plane_h
  use sto_m, only: sto
  use applied_m, only: share_applied
  integer,intent(in),optional::il
  logical,intent(in),optional::hd,hd_x,hd_y,hd_z,hx,hx_x,hx_y,hx_z,hk,hk_x,hk_y,hk_z,hv,htot,htot_x,htot_y,htot_z,hv_x,hv_y,hv_z, &
      wf_eff,wf,wf_x,wf_y,wf_z,ha,ha_eff,ha_x,ha_y,ha_z,temp,tempf,hb,hb_x,hb_y,hb_z
  real(8),intent(in),optional::hb_ds(3),hb_r(3),hb_m(3),hb_tilt(3)

  integer::i,n,ng,id,k,j,i1,i2
  real(8)::dx,dy,dz
  type(field_ss),pointer::p,q
  character(200)::os

  if (get_plot()) then
    i1=1; i2=size(field)
    if (present(il)) then
      if (il.lt.0.or.il.gt.i2) return
      if (il.gt.0) then; i1=il; i2=il; endif
    endif
    do i=i1,i2; p=>field(i); ng=size(p%temp); id=p%in
      if (get_plot(i)) then
        if (present(wf_x).and.wf_x) then; write(os,"('Happ,x for layer',i3)")i; call plot_contour(trim(os),p%h_app(1,:),.true.,i,no_write=.true.); endif
        if (present(wf_y).and.wf_y) then; write(os,"('Happ,y for layer',i3)")i; call plot_contour(trim(os),p%h_app(2,:),.true.,i,no_write=.true.); endif
        if (present(wf_z).and.wf_y) then; write(os,"('Happ,z for layer',i3)")i; call plot_contour(trim(os),p%h_app(3,:),.true.,i,no_write=.true.); endif
        if (present(wf).and.wf) then
          do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_app(:,n),p%h_app(:,n))); enddo
          write(os,"('|Happ| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
        endif
        if (present(wf_eff).and.wf_eff) then; write(os,"('Hsw for layer',i3)")k; call plot_heff(trim(os),p%h_app(:,:),i,no_write=.true.); endif

        if (associated(app_field)) then
          if (.not.associated(app_plot)) allocate(app_plot(3,size(app_field,2),size(app_field,3)))
          if (present(ha_x).and.ha_x) then; write(os,"('full Happ,x for layer',i3)")i; call plot_contour(trim(os),app_field(1,:,:,i)+ufield(1,i),.true.,i,no_write=.true.); endif
          if (present(ha_y).and.ha_y) then; write(os,"('full Happ,y for layer',i3)")i; call plot_contour(trim(os),app_field(2,:,:,i)+ufield(2,i),.true.,i,no_write=.true.); endif
          if (present(ha_z).and.ha_z) then; write(os,"('full Happ,z for layer',i3)")i; call plot_contour(trim(os),app_field(3,:,:,i)+ufield(3,i),.true.,i,no_write=.true.); endif
          if (present(ha).and.ha) then
            do j=1,size(app_field,3); do k=1,size(app_field,2); app_plot(1,k,j)= &
                sqrt(dot_product(app_field(:,k,j,i)+ufield(:,i),app_field(:,k,j,i)+ufield(:,i))); enddo; enddo
            write(os,"('full |Happ| for layer',i3)")i; call plot_contour(trim(os),app_plot(1,:,:),.true.,i,no_write=.true.)
          endif
!         if (present(ha_eff).and.ha_eff) then; write(os,"('full Hsw for layer',i3)")k; call plot_heff(trim(os),app_field(:,:,:,i),i,no_write=.true.); endif
        endif

        if (present(temp).and.temp) then; write(os,"('temperature for layer',i3)")i; call plot_contour(trim(os),p%temp,.true.,i,no_write=.true.); endif
        if (associated(app_temp).and.present(tempf).and.tempf) then;
          write(os,"('full temperature for layer',i3)")i; call plot_contour(trim(os),app_temp(:,:,i)+utemp(i)+field(i)%background_t,.true.,i,no_write=.true.)
        endif

        if (selfdemag_p) then
          if (present(hd_x).and.hd_x) then; write(os,"('Hdemag,x for layer',i3)")i; call plot_contour(trim(os),p%h_selfdemag(1,:),.true.,i,no_write=.true.); endif
          if (present(hd_y).and.hd_y) then; write(os,"('Hdemag,y for layer',i3)")i; call plot_contour(trim(os),p%h_selfdemag(2,:),.true.,i,no_write=.true.); endif
          if (present(hd_z).and.hd_z) then; write(os,"('Hdemag,z for layer',i3)")i; call plot_contour(trim(os),p%h_selfdemag(3,:),.true.,i,no_write=.true.); endif
          if (present(hd).and.hd) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_selfdemag(:,n),p%h_selfdemag(:,n))); enddo
            write(os,"('|Hdemag| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        else
          if (present(hd_x).and.hd_x) then; write(os,"('Hdemag,x for layer',i3)")i; call plot_contour(trim(os),p%h_demag(1,:),.true.,i,no_write=.true.); endif
          if (present(hd_y).and.hd_y) then; write(os,"('Hdemag,y for layer',i3)")i; call plot_contour(trim(os),p%h_demag(2,:),.true.,i,no_write=.true.); endif
          if (present(hd_z).and.hd_z) then; write(os,"('Hdemag,z for layer',i3)")i; call plot_contour(trim(os),p%h_demag(3,:),.true.,i,no_write=.true.); endif
          if (present(hd).and.hd) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_demag(:,n),p%h_demag(:,n))); enddo
            write(os,"('|Hdemag| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        endif
        if ((present(hx).or.present(hx_x).or.present(hx_y).or.present(hx_z)).and.associated(p%hxgrain)) then
          do n=1,ng; p%h_plot(:,n)=get_intralayer_hx(n,p,id); enddo
          if (present(hx_x).and.hx_x) then; write(os,"('Hexch,x for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.); endif
          if (present(hx_y).and.hx_y) then; write(os,"('Hexch,y for layer',i3)")i; call plot_contour(trim(os),p%h_plot(2,:),.true.,i,no_write=.true.); endif
          if (present(hx_z).and.hx_z) then; write(os,"('Hexch,z for layer',i3)")i; call plot_contour(trim(os),p%h_plot(3,:),.true.,i,no_write=.true.); endif
          if (present(hx).and.hx) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_plot(:,n),p%h_plot(:,n))); enddo
            write(os,"('|Hexch| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        endif
        if (present(hb).or.present(hb_x).or.present(hb_y).or.present(hb_z)) then
  print *,'hb',present(hb_tilt)
          call get_geometry(i,dx,dy,dz)
          if (.not.associated(app_field)) call share_applied(app_temp,app_field,utemp,ufield)
          if (.not.associated(app_plot)) allocate(app_plot(3,size(app_field,2),size(app_field,3)))
          if (present(hb_tilt)) then
            call sto(hb_ds(1),hb_ds(2),hb_ds(3),hb_tilt(1),size(app_field,2),size(app_field,3),dx,dy,dz,hb_r,hb_m,app_plot)
          else
           call sto(hb_ds(1),hb_ds(2),hb_ds(3),size(app_field,2),size(app_field,3),dx,dy,dz, &
                   (/hb_r(1)+dx*size(app_field,2)*0.5d0,hb_r(2)+dy*size(app_field,3)*0.5d0,hb_r(3)+dz*0.5d0+hb_ds(3)*0.5d0/),hb_m,app_plot)
!            call block_to_plane_h(hb_ds(1),hb_ds(2),hb_ds(3),size(app_field,2),size(app_field,3),dx,dy,dz, &
!                    (/hb_r(1)+dx*size(app_field,2)*0.5d0,hb_r(2)+dy*size(app_field,3)*0.5d0,hb_r(3)+dz*0.5d0+hb_ds(3)*0.5d0/),hb_m,app_plot)
          endif
          if (present(hb_x).and.hb_x) then; write(os,"('Hb,x for layer',i3)")i; call plot_contour(trim(os),app_plot(1,:,:),.true.,i,no_write=.true.); endif
          if (present(hb_y).and.hb_y) then; write(os,"('Hb,y for layer',i3)")i; call plot_contour(trim(os),app_plot(2,:,:),.true.,i,no_write=.true.); endif
          if (present(hb_z).and.hb_z) then; write(os,"('Hb,z for layer',i3)")i; call plot_contour(trim(os),app_plot(3,:,:),.true.,i,no_write=.true.); endif
          if (present(hb).and.hb) then
            do j=1,size(app_field,3); do k=1,size(app_field,2); app_plot(1,k,j)=sqrt(dot_product(app_plot(:,k,j),app_plot(:,k,j))); enddo; enddo
            write(os,"('|Hb| for layer',i3)")i; call plot_contour(trim(os),app_plot(1,:,:),.true.,i,no_write=.true.)
            do j=1,size(app_field,3); do k=1,size(app_field,2); app_plot(1,k,j)=sqrt(dot_product(app_plot(1:2,k,j),app_plot(1:2,k,j))); enddo; enddo
            write(os,"('Hb_perp for layer',i3)")i; call plot_contour(trim(os),app_plot(1,:,:),.true.,i,no_write=.true.)
          endif
        endif
        if ((present(hv).or.present(hv_x).or.present(hv_y).or.present(hv_z)).and.associated(p%uhx)) then
          q=>field(i-1); do n=1,ng; p%h_plot(:,n)=get_interlayer_hx(n,p,q,p%uhxgrainindex,p%uhxgrain,p%uhx); enddo
          if (present(hv_x).and.hv_x) then; write(os,"('Hexch,x from above for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.); endif
          if (present(hv_y).and.hv_y) then; write(os,"('Hexch,y from above for layer',i3)")i; call plot_contour(trim(os),p%h_plot(2,:),.true.,i,no_write=.true.); endif
          if (present(hv_z).and.hv_z) then; write(os,"('Hexch,z from above for layer',i3)")i; call plot_contour(trim(os),p%h_plot(3,:),.true.,i,no_write=.true.); endif
          if (present(hv).and.hv) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_plot(:,n),p%h_plot(:,n))); enddo
            write(os,"('|Hexch| from above for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        endif
        if ((present(hv).or.present(hv_x).or.present(hv_y).or.present(hv_z)).and.associated(p%lhx)) then
          q=>field(i+1); do n=1,ng; p%h_plot(:,n)=get_interlayer_hx(n,p,q,p%lhxgrainindex,p%lhxgrain,p%lhx); enddo
          if (present(hv_x).and.hv_x) then; write(os,"('Hexch,x from below for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.); endif
          if (present(hv_y).and.hv_y) then; write(os,"('Hexch,y from below for layer',i3)")i; call plot_contour(trim(os),p%h_plot(2,:),.true.,i,no_write=.true.); endif
          if (present(hv_z).and.hv_z) then; write(os,"('Hexch,z from below for layer',i3)")i; call plot_contour(trim(os),p%h_plot(3,:),.true.,i,no_write=.true.); endif
          if (present(hv).and.hv) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_plot(:,n),p%h_plot(:,n))); enddo
            write(os,"('|Hexch| from below for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        endif
        if (present(hk).or.present(hk_x).or.present(hk_y).or.present(hk_z)) then
          do n=1,ng; p%h_plot(:,n)=get_cubic_hk(n,p,id)+get_uni_hk(n,p,id); enddo
          if (present(hk_x).and.hk_x) then; write(os,"('Hk,x for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.); endif
          if (present(hk_y).and.hk_y) then; write(os,"('Hk,y for layer',i3)")i; call plot_contour(trim(os),p%h_plot(2,:),.true.,i,no_write=.true.); endif
          if (present(hk_z).and.hk_z) then; write(os,"('Hk,z for layer',i3)")i; call plot_contour(trim(os),p%h_plot(3,:),.true.,i,no_write=.true.); endif
          if (present(hk).and.hk) then
            do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h_plot(:,n),p%h_plot(:,n))); enddo
            write(os,"('|Hk| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
          endif
        endif
        if (present(htot_x).and.htot_x) then; write(os,"('Htotal,x for layer',i3)")i; call plot_contour(trim(os),p%h(1,:),.true.,i,no_write=.true.); endif
        if (present(htot_y).and.htot_y) then; write(os,"('Htotal,y for layer',i3)")i; call plot_contour(trim(os),p%h(2,:),.true.,i,no_write=.true.); endif
        if (present(htot_z).and.htot_z) then; write(os,"('Htotal,z for layer',i3)")i; call plot_contour(trim(os),p%h(3,:),.true.,i,no_write=.true.); endif
        if (present(htot).and.htot) then
          do n=1,ng; p%h_plot(1,n)=sqrt(dot_product(p%h(:,n),p%h(:,n))); enddo
          write(os,"('|Htotal| for layer',i3)")i; call plot_contour(trim(os),p%h_plot(1,:),.true.,i,no_write=.true.)
        endif
      endif
    enddo
  endif
 end subroutine

end module field_m
