
!this module plots and keeps track of what needs to be plotted

module plot_m
  use graphix, only: draw

  implicit none
  private

  public set_plot,plot_histogram,plot_scatter,plot_contour,get_plot,plot_temp_scaling,plot_wf_scaling,plot_mk,plot_mh,plot_heff,plot_m3,plot_mz_scaling,plot_strip

  interface plot_contour
    module procedure plot_contour1,plot_contour2
  end interface
  interface plot_mk
    module procedure plot_mks,plot_mka
  end interface

  type plot_s
    logical::plot=.true.!.false.
    logical::grain=.true.     !collapse contours on to grain structure?
    logical::grain_hist=.true.      !plot grain area histograms
    logical::hk=.true.             !plot hk contour
    logical::hk_ang=.true.          !plot hk angular distribution
    logical::hk_hist=.true.         !plot hk histogram
    logical::hk_vol=.true.         !plot hk histogram
    logical::tc=.true.             !plot tc contour
    logical::tc_hist=.true.        !plot tc histogram
    logical::tc_vol=.true.         !plot tc vs vol histogram
    logical::tau=.true.             !plot tau contour
    logical::tau_hist=.true.        !plot tau histogram
    logical::tau_vol=.true.         !plot tau vs vol histogram
    logical::hx=.true.        !plot hx contour
    logical::hx_hist=.true.        !plot hx histogram
    logical::hx_vol=.true.         !plot hx vs vol histogram
    logical::hv=.true.        !plot interlayer hx contour
    logical::hv_hist=.true.        !plot interlayer hx histogram
    logical::hv_vol=.true.         !plot interlayer hx vs vol histogram
    logical::wf_x=.false.      !plot write field x component volume averaged
    logical::wf_y=.false.      !plot write field y component volume averaged
    logical::wf_z=.true.      !plot write field z component volume averaged
    logical::wf_xf=.false.      !plot write field x component
    logical::wf_yf=.false.      !plot write field y component
    logical::wf_zf=.false.      !plot write field z component
    logical::wf_eff=.true.    !plot write field effective
    logical::temp=.true.      !plot temperature contour
    logical::tempf=.false.      !plot temperature field contour
    logical::mag=.true.       !plot |m| from llb and/or Xu&Zhang
    logical::mx=.false.            !plot mx contour
    logical::my=.false.            !plot my contour
    logical::mz=.false.            !plot mz contour
    logical::mzs=.false.            !plot mz scaling
    logical::ms=.true.             !plot ms contour
    logical::mk=.true.             !plot m.k contour
    logical::scatter=.false.       !plot m scatter plot on sphere
    logical::scatter_len=.false.   !plot m scatter plot using length
    logical::hamr=.true.           !plot parameters vs temperature
    logical::mhloop=.true.         !plot m(h) loop
    logical::wf_scale=.true.       !plot temporal write field scaling
    logical::temp_scale=.true.     !plot temporal temperature scaling
    real(8),pointer::wfscale(:,:)=>null(),tempscale(:,:)=>null(),m_h(:,:)=>null(),mzscale(:,:)=>null(),huscale(:,:)=>null()
    character(2)::strip(20)="  "
  end type plot_s

  type(plot_s),allocatable,target,save::plotl(:)
  real(8),pointer::m_h(:,:)=>null()
  real(8),save::h_d(3)=(/0.d0, 0.d0, 1.d0/) !direction of h for m(h) loop

  logical,save::plot_p=.false.   !global plotting- must be true for any plotting to take effect
  logical,save::plot_to_file=.false.

  integer,save::nx,ny,n_lay=-1
  real(8),save::dx,dy
  real(8),allocatable,save::dz(:)
  real(8),allocatable,save::x(:),y(:)

contains

  subroutine plot_strip(os1,il,t,x,ionum)
   integer,intent(in)::il!,n
!  real(8),intent(in)::t,x(n)
   real(8),intent(in)::t,x(:)
   character(2),intent(in)::os1
   integer,intent(in),optional::ionum
   character(200)::os
   integer::i,j
   real(8)::tt(size(x))

   write(os,"(a,' for layer',i3)") trim(os1),il
   if (os1.eq.'Hu') then
     os='Huniform'
     write(os,"(a,' for layer',i3)") trim(os),il
   endif

   j=0
   do i=1,size(plotl(il)%strip)
     if (plotl(il)%strip(i).eq.os1) exit
     j=j+1
     if (plotl(il)%strip(i).eq.'  ') exit
   enddo
   if (i.le.size(plotl(il)%strip)) then
     if (i.eq.j) then
       call draw(trim(os),strip=.true.,legend=os1,x=(/t/))
     endif
     plotl(il)%strip(i)=os1
     tt=t;
     call draw(trim(os),strip=.true.,add=.true.,x=tt,y=x)
   endif

  end subroutine
  subroutine plot_temp_scaling(il,t,x,ionum)
   integer,intent(in)::il
   real(8),intent(in)::t,x(2)
   integer,intent(in),optional::ionum
   character(200)::os,os2
!  real(8),pointer::tmp(:,:)

   os2='temperature scaling'
   write(os,"(a,' for layer',i3)") trim(os2),il
   call plot_scaling(plotl(il)%tempscale,'temp',plotl(il)%temp_scale); return
  entry plot_wf_scaling(il,t,x,ionum)
   os2='write field scaling'
   write(os,"(a,' for layer',i3)") trim(os2),il
   call plot_scaling(plotl(il)%wfscale,'fld',plotl(il)%wf_scale); return
  entry plot_mz_scaling(il,t,x,ionum)
   os2='mz scaling'
   write(os,"(a,' for layer',i3)") trim(os2),il
   call plot_scaling(plotl(il)%mzscale,'mz',plotl(il)%mzs); return

  contains
    subroutine plot_scaling(scale,os1,plot_p)!,os2)
     real(8),pointer::tmp(:,:),scale(:,:)
     logical,intent(in)::plot_p
     character(*)::os1!,os2
     integer::i

     if (present(ionum)) then
!      if (associated(scale)) then
!        do i=1,size(scale,2); if (scale(1,i).lt.-1.d200) exit; enddo; i=i-1
!        if (scale(3,i).gt.-1.d22) then
!          write(ionum,"(1p,3(e12.5,x))") scale(1:3,1:i)
!        else
!          write(ionum,"(1p,2(e12.5,x))") scale(1:2,1:i)
!        endif
!      endif
       write(ionum,"(1p,:,30(e12.5,x))") t,x
       return
     endif
     if (.not.plot_p) return
     if (.not.associated(scale).or.size(x)+1.gt.size(scale,1)) then
       allocate(tmp(size(x)+1,1)); tmp(:,1)=(/t,x/)
       if (.not.associated(scale)) call draw(title=trim(os),strip=.true.,legend=os(1:2),x=(/t/))
       if (associated(scale)) deallocate(scale)
       scale=>tmp
     endif
     scale(:,1)=(/t,x/)
!    print *,size(scale),scale
     if (size(x).eq.1) then
       call draw(title=trim(os),strip=.true.,add=.true.,x=(/t/),y=x)
     elseif (size(x).eq.2) then
       call draw(title=trim(os),strip=.true.,add=.true.,x=(/t,t,t/),y=(/x(1),x(2),x(1)+x(2)/))
     elseif (size(x).eq.3) then
       call draw(title=trim(os),strip=.true.,add=.true.,x=(/t,t,t/),y=x)
     else
       call draw(title=trim(os),strip=.true.,add=.true.,x=(/t,t,t,t/),y=x)
     endif


!     if (.not.associated(scale)) then; allocate(scale(3,100)); scale=-1.d300; endif
!     do i=1,size(scale,2); if (scale(1,i).lt.-1.d200) exit; enddo; 
!     if (t.lt.scale(1,max(1,i-1))) then; i=1; scale=-1.d300; endif
!     if (i.gt.size(scale,2)) then
!       allocate(tmp(3,i-1+100)); tmp=-1.d300; tmp(1:3,1:i-1)=scale;
!       deallocate(scale); scale=>tmp
!     endif
!   
!     scale(1,i)=t; scale(2:3,i)=x; 
!     if (i.gt.1.and.plot_p) then
!       if (x(2).gt.-1.d22) then
!         call draw(title=trim(os),x=pack(scale(1,1:i),mask=(scale(1,1:i).gt.-1.d200)),&
!                                  y=pack(scale(2,1:i)+scale(3,1:i),mask=(scale(1,1:i).gt.-1.d200)), legend='total '//trim(os1),color=3, &
!            ymin=minval(scale(2:3,1:i)),ymax=maxval(scale(2,1:i)+scale(3,1:i)))
!         call draw(title=trim(os),x=scale(1,1:i),y=scale(2,1:i),legend=trim(os1)//' 1',add=.true.,color=1)
!         call draw(title=trim(os),x=scale(1,1:i),y=scale(3,1:i),legend=trim(os1)//' 2',add=.true.,color=2)
!       else
!         call draw(title=trim(os),x=scale(1,1:i),y=scale(2,1:i),color=1)
!       endif
!     endif
    end subroutine
  end subroutine

  subroutine plot_mh(hmag,time)
   use data_m, only:get_data
   use io_m, only:add_file, delete_file
   use transform_m, only:cross
   real(8),intent(in)::hmag,time
   integer::k,i,n,io1,io2,nup,nup1,ndn,ndn1
   integer,pointer::id
   real(8),pointer::m(:,:),ms(:),ms_scale(:),grain_vol(:),htot(:,:)
   real(8)::x(3),x1(3),norm,norm1,t(3),xs(3),torq2,torq(3),torq2t
   logical::ok_p
   character(80)::loopfile
   character(200)::os
!print *,'in plot_mh'
   call get_data(1,loopfile)
   io1=add_file(trim(loopfile))
   io2=add_file('ave_'//trim(loopfile))

   x=0.d0; norm=0.d0; t=0.d0; nup=0; ndn=0; xs=0.d0; torq2t=0.d0; ok_p=.false.
   do k=1,n_lay
     call get_data(k,id,in_p=.true.)
     call get_data(k,m,m_p=.true.)
     call get_data(k,ms,ms_p=.true.)
     call get_data(k,ms_scale,ms_scale_p=.true.)
     call get_data(k,grain_vol,grain_vol_p=.true.)
     call get_data(k,htot,h_p=.true.)
     x1=0.d0; norm1=0.d0; nup1=0; ndn1=0; torq2=0.d0
     do n=1,size(ms)
       x1=x1+m(1+id:3+id,n)*ms(n)*ms_scale(n)*grain_vol(n)
       norm1=norm1+ms(n)*ms_scale(n)*grain_vol(n)
       torq=cross(m(1+id:3+id,n),htot(:,n)); torq2=torq2+dot_product(torq,torq)*ms(n)*ms_scale(n)*grain_vol(n)
       if (dot_product(m(1+id:3+id,n),h_d).gt.0.d0) then
         nup1=nup1+1
       else
         ndn1=ndn1+1
       endif
     enddo
     write(os,"('m(h) for layer',i2)")k
     ok_p=(ok_p.or.(plotl(k)%plot.and.plotl(k)%mhloop))
     call add_data(plotl(k)%m_h,x1*norm1/(norm1*norm1+1.d-300),x1/sum(grain_vol(1:size(ms))),cross(x1,h_d),torq2,trim(os),io1, &
                trim(loopfile),k,nup1,ndn1,plotl(k)%mhloop.and.(plotl(k)%plot.or.plot_p))
     x=x+x1; norm=norm+norm1; t=t+cross(x1,h_d); nup=nup+nup1; ndn=ndn+ndn1; torq2t=torq2t+torq2
     xs=xs+x1/sum(grain_vol(1:size(ms)))
   enddo

   if (n_lay.gt.1) call add_data(m_h,x*norm/(norm*norm+1.d-300),xs/n_lay,t,torq2t,'m(h)',io2,'ave_'//trim(loopfile),n_lay,nup,ndn,ok_p.or.plot_p) !plot_p)
   call delete_file(io2); call delete_file(io1)
!read *

  contains
    subroutine add_data(m_h,y,ys,t,tq2,os,io2,os2,k,nup,nud,plot_p)
     real(8),intent(in)::y(3),t(3),ys(3),tq2
     real(8),pointer::m_h(:,:),tmp(:,:)
     integer,intent(in)::io2,k,nup,nud
     character(*)::os,os2
     logical::ok_p,plot_p

     if (.not.associated(m_h)) then; allocate(m_h(5,100)); m_h=-1.d300; endif
     do i=1,size(m_h,2); if (m_h(1,i).lt.-1.d200) exit; enddo; 
     if (i.gt.size(m_h,2)) then
       allocate(tmp(5,i-1+100)); tmp=-1.d300; tmp(1:5,1:i-1)=m_h;
       deallocate(m_h); m_h=>tmp
     endif
     m_h(:,i)=(/ hmag, y(1), y(2), y(3), dot_product(y,h_d) /)
     if (plot_p.and.i.gt.1) then
       call draw(title=trim(os),x=m_h(1,1:i),y=m_h(2,1:i),legend='mx',symbol=.true.,nps=4,color=9,ymin=minval(m_h(2:,:i)),ymax=maxval(m_h(2:,:i)))
       call draw(title=trim(os),x=m_h(1,1:i),y=m_h(3,1:i),legend='my',add=.true.,symbol=.true.,nps=4,color=3)
       call draw(title=trim(os),x=m_h(1,1:i),y=m_h(4,1:i),legend='mz',add=.true.,symbol=.true.,nps=4,color=1)
       call draw(title=trim(os),x=m_h(1,1:i),y=m_h(5,1:i),legend='m.h',add=.true.,unsort=.true.,color=15)
     endif

     inquire(file=trim(os2),exist=ok_p)
     if (ok_p) then
       open(io2,file=trim(os2),status='old',position='append')
     else
       open(io2,file=trim(os2),status='new')
       write(io2,"('layer    H (T)      m_x      m_y     m_z     m.h      t      N+      N-  torque (emu T)  Ms|mxH|^2 Hx (T)     Hy(T)    Hz(T)   Mx  My  Mz  M.h')")
     endif
     write(io2,"(:,0p,i3,1p,77(1x,e12.5))") k,hmag,y,dot_product(y,h_d),time,dble(nup),dble(nud),hmag*sqrt(dot_product(t,t)),tq2,hmag*h_d,ys,dot_product(ys,h_d)
     close(io2)
    end subroutine
  end subroutine


  subroutine set_plot(il,plotting,plot,grain,grain_hist,hk,hk_ang,hk_hist,hk_vol,tc,tc_hist,tc_vol,hx,hx_hist,hx_vol,hv,hv_hist,hv_vol, &
           tempf,wf_xf,wf_yf,wf_zf, &
           wf_x,wf_y,wf_z,wf_eff,temp,mag,mx,my,mz,ms,mk,scatter,scatter_len,hamr,mhloop,wf_scale,temp_scale,mh,tau,tau_hist,tau_vol,mzs,write_plot)
   use data_m, only:get_size, get_geometry
   use graphix, only:can_draw
   integer,intent(in),optional::il
   logical,intent(in),optional::plot,grain,grain_hist,hk,hk_ang,hk_hist,hk_vol,tau,tau_hist,tau_vol,tc,tc_hist,tc_vol,hx,hx_hist,hx_vol,hv,hv_hist,hv_vol, &
                      wf_x,wf_y,wf_z,wf_eff,temp,mag,mx,my,mz,ms,mk,scatter,scatter_len,hamr,mhloop,wf_scale,temp_scale,plotting,mzs,write_plot, &
                      tempf,wf_xf,wf_yf,wf_zf
   real(8),intent(in),optional::mh(3)
   real(8),allocatable::xx(:)
   integer::i,i1,i2
   real(8)::z
   type(plot_s),pointer::p
   i=0; if (present(il)) i=il
   if (present(mh)) then; z=dot_product(mh,mh); if (z.gt.0.d0) h_d=mh/sqrt(z); endif
   if (n_lay.lt.0) then
     call get_size(nx,ny,n_lay,dx,dy)
     if (.not.allocated(dz)) allocate(dz(n_lay),xx(n_lay)); call get_geometry(dx,dy,dz,xx); deallocate(xx)
     if (.not.allocated(plotl)) allocate(plotl(n_lay))
     if (.not.allocated(x)) then
       allocate(x(nx),y(ny))
       do i1=1,nx; x(i1)=(i1-0.5d0)*dx; enddo; do i1=1,ny; y(i1)=(i1-0.5d0)*dy; enddo; x=x*1.d7; y=y*1.d7
     endif
   endif
   if (present(plotting)) plot_p=(plotting.and.can_draw())
   if (present(write_plot)) plot_to_file=write_plot
   i1=1; i2=n_lay
   if (i.ne.0) then; i1=i; i2=i; endif
   do i=i1,i2
     p=>plotl(i)
     if (present(plot)) p%plot=plot;                              if (present(grain)) p%grain=grain;
     if (present(grain_hist)) p%grain_hist=grain_hist;            if (present(hk)) p%hk=hk;
     if (present(hk_ang)) p%hk_ang=hk_ang;                        if (present(hk_hist)) p%hk_hist=hk_hist;
     if (present(hk_vol)) p%hk_vol=hk_vol;                        if (present(tc)) p%tc=tc;
     if (present(tc_hist)) p%tc_hist=tc_hist;                     if (present(tc_vol)) p%tc_vol=tc_vol;
     if (present(tau)) p%tau=tau; if (present(tau_hist)) p%tau_hist=tau_hist;                     if (present(tau_vol)) p%tau_vol=tau_vol;
     if (present(hx)) p%hx=hx;                                    if (present(hx_hist)) p%hx_hist=hx_hist;
     if (present(hx_vol)) p%hx_vol=hx_vol;                        if (present(hv)) p%hv=hv;
     if (present(hv_hist)) p%hv_hist=hv_hist;                     if (present(hv_vol)) p%hv_vol=hv_vol;
     if (present(wf_x)) p%wf_x=wf_x;                              if (present(wf_y)) p%wf_y=wf_y;
     if (present(wf_z)) p%wf_z=wf_z;                              if (present(wf_eff)) p%wf_eff=wf_eff;
     if (present(wf_xf)) p%wf_xf=wf_xf;                           if (present(wf_yf)) p%wf_yf=wf_yf;
     if (present(wf_zf)) p%wf_zf=wf_zf;                           if (present(tempf)) p%tempf=tempf;
     if (present(temp)) p%temp=temp;                              if (present(mag)) p%mag=mag;
     if (present(mx)) p%mx=mx;                                    if (present(my)) p%my=my;
     if (present(mz)) p%mz=mz;                                    if (present(ms)) p%ms=ms;
     if (present(mk)) p%mk=mk;                                    if (present(scatter)) p%scatter=scatter;
     if (present(scatter_len)) p%scatter_len=scatter_len;         if (present(hamr)) p%hamr=hamr;
     if (present(mhloop)) p%mhloop=mhloop;                        if (present(wf_scale)) p%wf_scale=wf_scale;
     if (present(temp_scale)) p%temp_scale=temp_scale;            if (present(mzs)) p%mzs=mzs
!    if (present(mhloop)) p%plot=mhloop
!if (p%mhloop) read *
   enddo
  end subroutine


  function get_plot(il,plot,grain,grain_hist,hk,hk_ang,hk_hist,hk_vol,tc,tc_hist,tc_vol,hx,hx_hist,hx_vol,hv,hv_hist,hv_vol, &
           wf_x,wf_y,wf_z,wf_eff,temp,mag,mx,my,mz,ms,mk,scatter,scatter_len,hamr,mhloop,wf_scale,temp_scale, &
           wf_xf,wf_yf,wf_zf,tempf,tau,tau_hist,tau_vol,write_plot) result (ok_p)
   integer,intent(in),optional::il
   logical,intent(in),optional::plot,grain,grain_hist,hk,hk_ang,hk_hist,hk_vol,tc,tc_hist,tc_vol,hx,hx_hist,hx_vol,hv,hv_hist,hv_vol, &
           wf_x,wf_y,wf_z,wf_eff,temp,mag,mx,my,mz,ms,mk,scatter,scatter_len,hamr,mhloop,wf_scale,temp_scale,tau,tau_hist,tau_vol, &
           wf_xf,wf_yf,wf_zf,tempf,write_plot
   logical::ok_p
   integer::i
   type(plot_s),pointer::p

   i=0; if (present(il)) i=il
   ok_p=plot_p; if (present(write_plot)) ok_p=(ok_p.and.plot_to_file)
   if (.not.ok_p.or.i.lt.1.or.i.gt.n_lay) return
   p=>plotl(i); ok_p=p%plot; if (.not.ok_p) return
   if (present(plot)) ok_p=(ok_p.and.p%plot)
   if (present(grain)) ok_p=(ok_p.and.p%grain)
   if (present(grain_hist)) ok_p=(ok_p.and.p%grain_hist)
   if (present(hk)) ok_p=(ok_p.and.p%hk)
   if (present(hk_ang)) ok_p=(ok_p.and.p%hk_ang.and.p%hk)
   if (present(hk_hist)) ok_p=(ok_p.and.p%hk_hist.and.p%hk)
   if (present(hk_vol)) ok_p=(ok_p.and.p%hk_vol.and.p%hk)
   if (present(tc)) ok_p=(ok_p.and.p%tc)
   if (present(tc_hist)) ok_p=(ok_p.and.p%tc_hist.and.p%tc)
   if (present(tc_vol)) ok_p=(ok_p.and.p%tc_vol.and.p%tc)
   if (present(tau)) ok_p=(ok_p.and.p%tau)
   if (present(tau_hist)) ok_p=(ok_p.and.p%tau_hist.and.p%tau)
   if (present(tau_vol)) ok_p=(ok_p.and.p%tau_vol.and.p%tau)
   if (present(hx)) ok_p=(ok_p.and.p%hx)
   if (present(hx_hist)) ok_p=(ok_p.and.p%hx_hist.and.p%hx)
   if (present(hx_vol)) ok_p=(ok_p.and.p%hx_vol.and.p%hx)
   if (present(hv)) ok_p=(ok_p.and.p%hv)
   if (present(hv_hist)) ok_p=(ok_p.and.p%hv_hist.and.p%hv)
   if (present(hv_vol)) ok_p=(ok_p.and.p%hv_vol.and.p%hv)
   if (present(wf_x)) ok_p=(ok_p.and.p%wf_x)
   if (present(wf_y)) ok_p=(ok_p.and.p%wf_y)
   if (present(wf_z)) ok_p=(ok_p.and.p%wf_z)
   if (present(wf_xf)) ok_p=(ok_p.and.p%wf_xf)
   if (present(wf_yf)) ok_p=(ok_p.and.p%wf_yf)
   if (present(wf_zf)) ok_p=(ok_p.and.p%wf_zf)
   if (present(wf_eff)) ok_p=(ok_p.and.p%wf_eff)
   if (present(temp)) ok_p=(ok_p.and.p%temp)
   if (present(tempf)) ok_p=(ok_p.and.p%tempf)
   if (present(mag)) ok_p=(ok_p.and.p%mag)
   if (present(mx)) ok_p=(ok_p.and.p%mx)
   if (present(my)) ok_p=(ok_p.and.p%my)
   if (present(mz)) ok_p=(ok_p.and.p%mz)
   if (present(ms)) ok_p=(ok_p.and.p%ms)
   if (present(mk)) ok_p=(ok_p.and.p%mk)
   if (present(scatter)) ok_p=(ok_p.and.p%scatter)
   if (present(scatter_len)) ok_p=(ok_p.and.p%scatter_len)
   if (present(hamr)) ok_p=(ok_p.and.p%hamr)
   if (present(mhloop)) ok_p=(ok_p.and.p%mhloop)
   if (present(wf_scale)) ok_p=(ok_p.and.p%wf_scale)
   if (present(temp_scale)) ok_p=(ok_p.and.p%temp_scale)
  end function

  subroutine plot_contour1(tit,z,grain_p,i,n0,zmin,zmax,box,no_write)
   use data_m, only:get_data
   use window_m, only:get_window
   character(*),intent(in)::tit
   real(8),intent(in)::z(:)
   logical,intent(in)::grain_p
   integer,intent(in)::i
   integer,intent(in),optional::n0
   real(8),intent(in),optional::zmin,zmax,box(:,:)
   logical,intent(in),optional::no_write
   integer,pointer::g2gindex(:),g2g(:,:),gg(:,:)

   real(8)::z1(nx,ny)
   integer::n,j,nc
   integer::nxr(2),nyr(2)

   nc=16; if (present(n0)) nc=n0

   call get_data(i,g2gindex,grain2gridindex_p=.true.)
   call get_data(i,g2g,grain2grid_p=.true.)
!if (lbound(g2gindex,dim=1).ne.0) then; print *,'plot_contour1';read *;endif

   call get_window(nxr=nxr,nyr=nyr)
   j=nxr(2)-nxr(1)+1; if (3*j.lt.nx) then; nxr(1)=max(1,nxr(1)-j); nxr(2)=min(nx,nxr(1)+3*j-1); nxr(1)=max(1,nxr(2)-3*j+1); endif
   j=nyr(2)-nyr(1)+1; if (3*j.lt.ny) then; nyr(1)=max(1,nyr(1)-j); nyr(2)=min(ny,nyr(1)+3*j-1); nyr(1)=max(1,nyr(2)-3*j+1); endif

   z1=-1.d22
   do n=1,size(z)
     do j=g2gindex(n-1)+1,g2gindex(n)
       z1(g2g(1,j),g2g(2,j))=z(n)
     enddo
   enddo
   if (grain_p) then
     call get_data(i,gg,grid2grain_p=.true.)
     if (present(zmin).and.present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,zmax=zmax,box=box)
     elseif (present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,box=box)
     elseif (present(zmax)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),zmin=minval(z),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmax=zmax,box=box)
     else
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),zmin=minval(z),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,box=box)
     endif
   else
     if (present(zmin).and.present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,zmax=zmax,box=box)
     elseif (present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,box=box)
     elseif (present(zmax)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),zmin=minval(z),num_color=nc,zmax=zmax,box=box)
     else
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z1(nxr(1):nxr(2),nyr(1):nyr(2)),zmin=minval(z),num_color=nc,box=box)
     endif
   endif
   if (plot_to_file) call write_plot(tit,z1,no_write)
  end subroutine

  subroutine plot_contour2(tit,z,grain_p,i,n0,zmin,zmax,box,no_write)
   use data_m, only:get_data 
   use window_m, only: get_window
   character(*),intent(in)::tit
   real(8),intent(in)::z(:,:)
   logical,intent(in)::grain_p
   integer,intent(in)::i
   integer,intent(in),optional::n0
   real(8),intent(in),optional::zmin,zmax,box(:,:)
   logical,intent(in),optional::no_write
   integer,pointer::gg(:,:)
   integer::nc
   integer::nxr(2),nyr(2),j

   nc=16; if (present(n0)) nc=n0

   call get_window(nxr=nxr,nyr=nyr)
   j=nxr(2)-nxr(1)+1; if (3*j.lt.nx) then; nxr(1)=max(1,nxr(1)-j); nxr(2)=min(nx,nxr(1)+3*j-1); nxr(1)=max(1,nxr(2)-3*j+1); endif
   j=nyr(2)-nyr(1)+1; if (3*j.lt.ny) then; nyr(1)=max(1,nyr(1)-j); nyr(2)=min(ny,nyr(1)+3*j-1); nyr(1)=max(1,nyr(2)-3*j+1); endif

   if (grain_p) then
     call get_data(i,gg,grid2grain_p=.true.)
     if (present(zmin).and.present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,zmax=zmax,box=box)
     elseif (present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,box=box)
     elseif (present(zmax)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmax=zmax,box=box)
     else
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),grain=gg(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,box=box)
     endif
   else
     if (present(zmin).and.present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,zmax=zmax,box=box)
     elseif (present(zmin)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmin=zmin,box=box)
     elseif (present(zmax)) then
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,zmax=zmax,box=box)
     else
       call draw(title=trim(tit),x=x(nxr(1):nxr(2)),y=y(nyr(1):nyr(2)),z=z(nxr(1):nxr(2),nyr(1):nyr(2)),num_color=nc,box=box)
     endif
   endif
   if (plot_to_file) call write_plot(tit,z,no_write)
  end subroutine

  subroutine plot_m3(tit,m3,mag,mags,id,il)
   character(*),intent(in)::tit
   integer,intent(in)::il,id
   real(8),intent(in)::m3(:,:),mag(:,:),mags(:)
   if (get_plot(il,scatter_len=.true.)) then
     if (any(mag(1+id/3,:).ne.0.d0)) then
       call draw(title=trim(tit),scatter=m3(1+id:3+id,:),scatter2=m3(4-id:6-id,:),scatter3=mag(1+id/3,:)/maxval(mag(1+id/3,:)))
     else
       call draw(title=trim(tit),scatter=m3(1+id:3+id,:),scatter2=m3(4-id:6-id,:),scatter3=mags/maxval(mags))
     endif
   else
     call draw(title=trim(tit),scatter=m3(1+id:3+id,:),scatter2=m3(4-id:6-id,:))
   endif
  end subroutine

  subroutine plot_heff(tit,heff,i,no_write)
   character(*),intent(in)::tit
   integer,intent(in)::i
   real(8),intent(in)::heff(:,:)
   logical,intent(in),optional::no_write
   real(8)::h(size(heff,2))
   integer::n
   do n=1,size(heff,2)
     h(n)=((heff(3,n)*heff(3,n))**(1.d0/3.d0)+(heff(2,n)*heff(2,n)+heff(1,n)*heff(1,n))**(1.d0/3.d0))**1.5d0
   enddo
   call plot_contour(trim(tit),h,.true.,i,no_write=no_write)
  end subroutine
  subroutine plot_mka(tit,m,k,i,scl,no_write)
   use demag_m, only:get_demag_window
   character(*),intent(in)::tit
   integer,intent(in)::i
   real(8),intent(in)::k(:,:),m(:,:)
   real(8),intent(in),optional::scl(:)
   logical,intent(in),optional::no_write
   real(8)::box(5,2),mk(size(m,2))
   integer::n
   box=get_demag_window(dx,dy)
   do n=1,size(m,2)
     mk(n)=dot_product(k(:,n),m(:,n))
   enddo
   if (present(scl))mk=mk*scl
   call plot_contour(trim(tit),mk,.true.,i,-6,zmin=-1.d0,zmax=1.d0,box=box*1d7,no_write=no_write)
  end subroutine
  subroutine plot_mks(tit,m,k,i,scl,no_write)
   use demag_m, only:get_demag_window
   character(*),intent(in)::tit
   integer,intent(in)::i
   real(8),intent(in)::k(:),m(:,:)
   real(8),intent(in),optional::scl(:)
   logical,intent(in),optional::no_write
   real(8)::box(5,2),mk(size(m,2))
   integer::n
   box=get_demag_window(dx,dy)
   do n=1,size(m,2)
     mk(n)=dot_product(k(:),m(:,n))
   enddo
   if (present(scl))mk=mk*scl
   call plot_contour(trim(tit),mk,.true.,i,-6,zmin=-1.d0,zmax=1.d0,box=box*1d7,no_write=no_write)
  end subroutine

  subroutine plot_scatter(tit,x,y,ave_p)
   character(*),intent(in)::tit
   real(8),intent(in)::x(:),y(:)
   logical,intent(in)::ave_p
   real(8)::xave,xd(150,2),dx,ymin
   integer::i,k
   character(200)::os
   call draw(title=trim(tit),x=x,y=y,symbol=.true.)
   if (ave_p) then
     xave=0.d0;do i=1,size(x); xave=xave+x(i)*y(i); enddo; xave=xave/sum(x)
     write(os,"('average=',f7.2)") xave
     xave=minval(x);dx=(maxval(x)-xave)/size(xd,1); xd=0.d0; ymin=minval(y)
     if (dx.gt.0.d0) then; dx=1.d0/dx
       do i=1,size(x)
         k=max(1,min(size(xd,1),int((x(i)-xave)*dx)+1)); xd(k,1)=xd(k,1)+y(i); xd(k,2)=xd(k,2)+1.d0
       enddo; dx=1.d0/dx
       do i=1,size(xd,1)
         if (xd(i,2).gt.0.d0) then
           xd(i,1)=xd(i,1)/xd(i,2)
         else
           xd(i,1)=ymin-maxval(x)
         endif
         xd(i,2)=xave+(i-0.5d0)*dx
       enddo
       call draw(title=trim(tit),x=xd(:,2),y=xd(:,1),ymin=ymin,symbol=.true.,add=.true.,legend=trim(os),nps=4)
     endif
   endif
  end subroutine

  subroutine plot_histogram(tit,n,x,ionum,xmin,xmax)
   character(*),intent(in)::tit
   integer,intent(in)::n
   real(8),intent(in)::x(:)
   real(8),intent(in),optional::xmin,xmax
   integer,intent(in),optional::ionum
   real(8)::xave,y
   integer::i
   character(200)::os

   if (present(ionum)) then; write(ionum,"(1p,e12.5)") x; return; endif

   call draw(title=trim(tit),histogram=n,x=x,yrange=y,xmin=xmin,xmax=xmax)
   xave=0.d0;do i=1,size(x); xave=xave+x(i); enddo; xave=xave/max(1,size(x))
   write(os,"('average=',f9.4)") xave
   call draw(title=trim(tit),x=(/xave,xave/),y=(/0.d0,1.5d0*y/),add=.true.,legend=trim(os),color=1)
  end subroutine

  subroutine write_plot(tit,z,no_write)
    use io_m, only:add_file,delete_file
    character(*),intent(in)::tit
    real(8),intent(in)::z(:,:)
    logical,intent(in),optional::no_write

    character(400)::fname
    integer::i,j,ionum


    if (present(no_write)) then; if (no_write) return; endif
    fname=find_filename(tit)
    ionum=add_file(trim(fname))
    open(ionum,file=trim(fname),status='unknown')
    do j=1,size(z,2)
      write(ionum,"(1p,:,2048(x,e12.5))") z(:,j)
    enddo
    close(ionum)
    call delete_file(ionum)
    
  end subroutine

  function find_filename(tit) result (fname)
    character(*)::tit
    character(400)::fname
    character(6)::s
    logical::ok_p
    integer::i

    fname=trim(adjustl(tit))//'_'
    do i=1,len_trim(fname)
      if (fname(i:i).eq.'|'.or.fname(i:i).eq.'/'.or.fname(i:i).eq.' '.or.fname(i:i).eq.',') fname(i:i)='_'
    enddo
    i=0; ok_p=.true.
    do while (ok_p)
      write(s,"(i6.6)")i
      inquire(file=trim(fname)//trim(adjustl(s))//'.dat',exist=ok_p)
      if (ok_p) i=i+1
!print *,trim(fname)//trim(adjustl(s))//'.dat?',ok_p
    enddo
    fname=trim(fname)//trim(adjustl(s))//'.dat'
  end function

end module plot_m

