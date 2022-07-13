!this module deals w/ applied fields/temperatures

module applied_m

  implicit none
  private
  public::start_applied,file_fld,file_temp,ref_fld,ref_temp,dup_fld,dup_temp,suprgaus_temp,rot_fld,uniform_fld,show_fld,show_temp, &
      set_fld,set_temp,show_fld_set,show_temp_set, &
      init_appl,get_appl,plot_applied_fld,plot_applied_temp,share_applied,get_applied_scale

   interface share_applied
     module procedure share_applied_all,share_applied_layer
   end interface

  type::applied_s
   integer::index=0       !field index
   integer::nx=0,ny=0       !size of field data
   real(8)::offset_x=-1.d300,offset_y=-1.d300,t_origin=0.d0,ref_x=0.d0,ref_y=0.d0   !offset, reference and time origin
   real(8)::vel_x=0.d0,vel_y=0.d0 !velocity the field is moving (cm/s)
   real(8),dimension(:),pointer::xpos=>null(),ypos=>null()           !coordinate values of field (cm)
   real(8),dimension(:,:),pointer::sfld=>null()                        !scalar field
   real(8),dimension(:,:,:),pointer::vfld=>null()                      !vector field
   real(8)::uniform_h(3)=(/0.d0, 0.d0, 0.d0/)
   real(8)::uniform_t=0.d0
   real(8)::dtemp=0.d0,sigmax=0.d0,sigmay=0.d0                              !supragaussian (only for temperature profiles)
   logical::analytic=.false.
   integer::expn=0
   integer,dimension(:),pointer::ix=>null(),iy=>null()                              !interpolation arrays
   real(8),dimension(:),pointer::fx=>null(),fy=>null()
   type(applied_s),pointer::next=>null()
  end type applied_s

  type(applied_s),pointer,save::all_temp=>null(),all_fld=>null()    !structure that holds all the temperature profils/fields

  type::current_applied_s
   type(applied_s),pointer::fld=>null(),fld2=>null()
   real(8)::t_start=-1.d22,t_end=1.d22,scale1_start=1.d0,scale1_end=1.d0,scale2_start=1.d0,scale2_end=1.d0,scale1=1.d0,scale2=1.d0
   character(3)::func1='   ',func2='   '
   real(8)::gamma1=0.d0,phase1=0.d0,amp1=0.d0,dc1=0.d0,param1=0.d0,duty1=0.d0,moddepth1=0.d0
   real(8)::gamma2=0.d0,phase2=0.d0,amp2=0.d0,dc2=0.d0,param2=0.d0,duty2=0.d0,moddepth2=0.d0
   logical::same_p=.false.,done_p=.false.,const1_p=.false.,const2_p=.false.
   logical::uniform1_p=.false.,uniform2_p=.false.
   real(8)::last_time=-1.d300                                !time we last did interpolation arrays
   real(8),dimension(:,:),pointer::snonmov_fld1=>null(),snonmov_fld2=>null() !optimization for non-moving fields- fields interpolated onto simulation mesh
   real(8),dimension(:,:,:),pointer::vnonmov_fld1=>null(),vnonmov_fld2=>null() !optimization for non-moving fields- fields interpolated onto simulation mesh

  end type current_applied_s

  integer,save::nxs=0,nys=0,nxs0=1,nxs1=0,nys0=1,nys1=0
  real(8),save::dxs=0.d0,dys=0.d0

  type(current_applied_s),allocatable,dimension(:),target,save::current_temp,current_fld

  real(8),allocatable,dimension(:),target,save::last_utemp(:)
  real(8),allocatable,dimension(:,:),target,save::tempscale(:,:),fldscale(:,:),last_ufld(:,:)
  real(8),allocatable,dimension(:,:,:),target,save::last_temp(:,:,:)
  real(8),allocatable,dimension(:,:,:,:),target,save::last_fld(:,:,:,:)

contains

  subroutine get_applied_scale(tscale,hscale)
    real(8),pointer::tscale(:,:),hscale(:,:)
    tscale=>tempscale; hscale=>fldscale
  end subroutine

  subroutine share_applied_all(temperature,field,utemperature,ufield)
    real(8),pointer::temperature(:,:,:),field(:,:,:,:),utemperature(:),ufield(:,:)
    temperature=>last_temp; field=>last_fld; utemperature=>last_utemp; ufield=>last_ufld
  end subroutine

  subroutine share_applied_layer(temperature,field,utemperature,ufield,layer)
    integer,intent(in)::layer
    real(8),pointer::temperature(:,:),field(:,:,:),utemperature,ufield(:)
    temperature=>last_temp(:,:,layer); field=>last_fld(:,:,:,layer); utemperature=>last_utemp(layer); ufield=>last_ufld(:,layer)
  end subroutine

  subroutine plot_applied_temp(layer)
   use plot_m, only: get_plot,plot_contour
   integer,intent(in)::layer
   character(200)::os

   if (get_plot(layer,temp=.true.)) then
     write(os,"('temperature for layer #',i3)") layer
     call plot_contour(trim(os),last_temp(:,:,layer)+last_utemp(layer),get_plot(layer,grain=.true.),layer,no_write=.true.)
   endif
  end subroutine
  subroutine plot_applied_fld(indx,layer)
   use plot_m, only: get_plot, plot_contour
   integer,intent(in)::indx,layer
   character(200)::os
   integer::i,j
   real(8)::gtmp(nxs,nys)

   select case (indx)
    case (1)
     if (get_plot(layer,wf_x=.true.)) then
       write(os,"('Hx for layer #',i3)") layer
       call plot_contour(trim(os),last_fld(1,:,:,layer)+last_ufld(1,layer),get_plot(layer,grain=.true.),layer,no_write=.true.)
     endif
    case (2)
     if (get_plot(layer,wf_y=.true.)) then
       write(os,"('Hy for layer #',i3)") layer
       call plot_contour(trim(os),last_fld(2,:,:,layer)+last_ufld(2,layer),get_plot(layer,grain=.true.),layer,no_write=.true.)
     endif
    case (3)
     if (get_plot(layer,wf_z=.true.)) then
       write(os,"('Hz for layer #',i3)") layer
       call plot_contour(trim(os),last_fld(3,:,:,layer)+last_ufld(3,layer),get_plot(layer,grain=.true.),layer,no_write=.true.)
     endif
    case (4)
     if (get_plot(layer,wf_eff=.true.)) then
       write(os,"('Hsw for layer #',i3)") layer
       do j=1,nys; do i=1,nxs
         gtmp(i,j)=((last_fld(3,i,j,layer)+last_ufld(3,layer))**(2./3.)+((last_fld(1,i,j,layer)+last_ufld(1,layer))**2+ &
                  (last_fld(2,i,j,layer)+last_ufld(2,layer))**2)**(1./3.))**(1.5)
       enddo; enddo
       call plot_contour(trim(os),gtmp,get_plot(layer,grain=.true.),layer,no_write=.true.)
     endif
   end select
  end subroutine

  subroutine show_temp_set()
   use io_m

   type(current_applied_s),dimension(:),pointer::current=>null()
   integer::i
   character(200)::os,tpe
   real(8)::x,y

   current=>current_temp; tpe='temperature'

  entry show_fld_set()


   if (.not.associated(current)) then
     current=>current_fld; tpe='field'
   endif

   do i=1,size(current)
     write(os,"(' ',a,' profile for layer #',i3)") trim(tpe),i
     call output(trim(os))
     if (.not.associated(current(i)%fld)) then
       write(os,"('   warning- first ',a,' profile not set')")trim(tpe)
       call output(trim(os))
     else
       write(os,"(1p,'   t_start=',e12.5,', t_end=',e12.5)") current(i)%t_start,current(i)%t_end
       call output(trim(os))
       write(os,"('   first ',a,' profile index:',i3)") trim(tpe),current(i)%fld%index
       call output(trim(os))
       if (current(i)%func1.eq.'   ') then
         write(os,"(1p,'    scale_start=',e12.5,', scale_end=',e12.5,', constant scale=',e12.5)") &
           current(i)%scale1_start,current(i)%scale1_end,current(i)%scale1
       else
         write(os,"(1p,'    function=',a,', amplitude=',e12.5,', gamma=',e12.5,', phase=',e12.5,', dc=',e12.5,', param=',e12.5,', duty1=',e12.5,', mod depth=',e12.5)") &
           current(i)%func1,current(i)%amp1*current(i)%scale1,current(i)%gamma1,current(i)%phase1, &
           current(i)%dc1*current(i)%scale1,current(i)%param1,current(i)%duty1,current(i)%moddepth1
       endif
       call output(trim(os))
       if (.not.associated(current(i)%fld%sfld).and..not.associated(current(i)%fld%vfld)) then
         if (associated(current,current_fld)) then
           write(os,"(1p,'    uniform field: (',2(e12.5,','),e12.5,')')") current(i)%fld%uniform_h
           call output(trim(os))
         elseif (.not.current(i)%fld%analytic) then
           write(os,"(1p,'    uniform temperature: ',e12.5)") current(i)%fld%uniform_t
           call output(trim(os))
         else
           x=current(i)%fld%sigmax; y=current(i)%fld%sigmay
           write(os,"(1p,'    super-gaussian, sigmax=',e12.5,', sigmay=',e12.5,', dT=',e12.5,', expn=',i3)") x/(x*x+1.d-300), &
                 y/(y*y+1.d-300),current(i)%fld%dtemp,current(i)%fld%expn
           call output(trim(os))
         endif
       endif
     endif

     if (.not.associated(current(i)%fld2)) then
       write(os,"('   no second ',a,' profile')") trim(tpe); call output(trim(os))
     else
       write(os,"('   second ',a,' profile index:',i3)") trim(tpe),current(i)%fld2%index
       call output(trim(os))
       if (current(i)%func2.eq.'   ') then
         write(os,"(1p,'    scale_start=',e12.5,', scale_end=',e12.5,', constant scale=',e12.5)") &
           current(i)%scale2_start,current(i)%scale2_end,current(i)%scale2
       else
         write(os,"(1p,'    function=',a,', amplitude=',e12.5,', gamma=',e12.5,', phase=',e12.5,', dc=',e12.5,', param2=',e12.5,', duty2=',e12.5,', mod depth=',e12.5)") &
           current(i)%func2,current(i)%amp2*current(i)%scale2,current(i)%gamma2,current(i)%phase2, &
           current(i)%dc2*current(i)%scale2,current(i)%param2,current(i)%duty2,current(i)%moddepth2
       endif
       call output(trim(os))
       if (.not.associated(current(i)%fld2%sfld).and..not.associated(current(i)%fld2%vfld)) then
         if (associated(current,current_fld)) then
           write(os,"(1p,'    uniform field: (',2(e12.5,','),e12.5,')')") current(i)%fld2%uniform_h
           call output(trim(os))
         elseif (.not.current(i)%fld2%analytic) then
           write(os,"(1p,'    uniform temperature: ',e12.5)") current(i)%fld2%uniform_t
           call output(trim(os))
         else
           x=current(i)%fld2%sigmax; y=current(i)%fld2%sigmay
           write(os,"(1p,'    super-gaussian, sigmax=',e12.5,', sigmay=',e12.5,', dT=',e12.5,', expn=',i3)") x/(x*x+1.d-300), &
                 y/(y*y+1.d-300),current(i)%fld2%dtemp,current(i)%fld2%expn
           call output(trim(os))
         endif
       endif
     endif
   enddo

   nullify(current)

  end subroutine

  function set_temp(ifld,ifld2,t_start,t_end,s1_start,s1_end,s2_start,s2_end,scale1,scale2,scale,ilayer, &
                    func1,gamma1,phase1,amp1,dc1,param1,duty1,moddepth1,func2,gamma2,phase2,amp2,dc2,param2,duty2,moddepth2) result(ok_p)

   use io_m, only: output
   
   integer,intent(in)::ifld
   integer,intent(in),optional::ifld2,ilayer
   real(8),intent(in),optional::t_start,t_end,s1_start,s1_end,s2_start,s2_end,scale,scale1,scale2
   character(3),intent(in),optional::func1,func2
   real(8),intent(in),optional::gamma1,phase1,amp1,gamma2,phase2,amp2,dc1,dc2,param1,param2,duty1,duty2,moddepth1,moddepth2

   logical::ok_p 
   
   type(applied_s),pointer::fld_list=>null(),fld,fld2,fld3
   type(current_applied_s),dimension(:),pointer::current
   integer::k0,k1,k2


   fld_list=>all_temp; current=>current_temp

  entry  set_fld(ifld,ifld2,t_start,t_end,s1_start,s1_end,s2_start,s2_end,scale1,scale2,scale,ilayer, &
                    func1,gamma1,phase1,amp1,dc1,param1,duty1,moddepth1,func2,gamma2,phase2,amp2,dc2,param2,duty2,moddepth2) result(ok_p)

   if (.not.associated(fld_list)) then
     fld_list=>all_fld; current=>current_fld
   endif

   ok_p=.false.

   if (.not.associated(fld_list)) return

   !find the fields, fld and fld2 will point to them
   nullify(fld,fld2)
   fld3=>fld_list
   do while (associated(fld3))
     if (fld3%index.eq.ifld) fld=>fld3
     if (present(ifld2)) then; if (fld3%index.eq.ifld2) fld2=>fld3; endif
     fld3=>fld3%next
   enddo
   if (.not.associated(fld)) then; nullify(fld_list); return; endif
   if (present(ifld2)) then; if (.not.associated(fld2).or.ifld2.eq.ifld) then; nullify(fld_list); return; endif; endif

   !which layer(s)?
   k1=1; k2=size(current)
   if (present(ilayer)) then
     if (ilayer.gt.0.and.ilayer.le.size(current)) then; k1=ilayer; k2=ilayer; endif
   endif

   do k0=k1,k2   !loop over layers (either just one or all)
     current(k0)%fld=>fld
     current(k0)%fld2=>fld2

  !default values, override if we're passed something else
     current(k0)%t_start=-1.d22;      if (present(t_start )) current(k0)%t_start=t_start
     current(k0)%t_end  = 1.d22;      if (present(t_end   )) current(k0)%t_end  =t_end
     current(k0)%scale1_start= 1.d0;  if (present(s1_start)) current(k0)%scale1_start=s1_start
     current(k0)%scale1_end  = 1.d0;  if (present(s1_end  )) current(k0)%scale1_end  =s1_end
     current(k0)%scale2_start= 1.d0;  if (present(s2_start)) current(k0)%scale2_start=s2_start
     current(k0)%scale2_end  = 1.d0;  if (present(s2_end  )) current(k0)%scale2_end  =s2_end
     current(k0)%scale1 = 1.d0;       if (present(scale1  )) current(k0)%scale1 =scale1
     current(k0)%scale2 = 1.d0;       if (present(scale2  )) current(k0)%scale2 =scale2
     current(k0)%func1 = '   ';       if (present(func1   )) current(k0)%func1=func1
     current(k0)%gamma1 = 0.d0;       if (present(gamma1  )) current(k0)%gamma1=gamma1
     current(k0)%phase1 = 0.d0;       if (present(phase1  )) current(k0)%phase1=phase1
     current(k0)%amp1 = 0.d0;         if (present(amp1    )) current(k0)%amp1=amp1
     current(k0)%dc1 = 0.d0;          if (present(dc1     )) current(k0)%dc1=dc1
     current(k0)%duty1 = 0.d0;        if (present(duty1   )) current(k0)%duty1=max(0.d0,min(1.d0,duty1))
     current(k0)%moddepth1 = 0.d0;    if (present(moddepth1)) current(k0)%moddepth1=max(0.d0,min(1.d0,moddepth1))
     current(k0)%param1 = 0.d0;       if (present(param1  )) current(k0)%param1=param1
     current(k0)%func2 = '   ';       if (present(func2   )) current(k0)%func2=func2
     current(k0)%gamma2 = 0.d0;       if (present(gamma2  )) current(k0)%gamma2=gamma2
     current(k0)%phase2 = 0.d0;       if (present(phase2  )) current(k0)%phase2=phase2
     current(k0)%amp2 = 0.d0;         if (present(amp2    )) current(k0)%amp2=amp2
     current(k0)%dc2 = 0.d0;          if (present(dc2     )) current(k0)%dc2=dc2
     current(k0)%param2 = 0.d0;       if (present(param2  )) current(k0)%param2=param2
     current(k0)%duty2 = 0.d0;        if (present(duty2   )) current(k0)%duty2=max(0.d0,min(1.d0,duty2))
     current(k0)%moddepth2 = 0.d0;    if (present(moddepth2)) current(k0)%moddepth2=max(0.d0,min(1.d0,moddepth2))
     current(k0)%last_time=-1.300

     if ((current(k0)%amp1.eq.0.d0.or.current(k0)%gamma1.eq.0.d0).and.current(k0)%func1.ne.'   ') then
       current(k0)%scale1_start=current(k0)%amp1*funct(current(k0)%func1,0.d0,current(k0)%gamma1,current(k0)%phase1,current(k0)%param1, &
          current(k0)%duty1,current(k0)%moddepth1) +current(k0)%dc1
       current(k0)%func1='   '; current(k0)%scale1_end = current(k0)%scale1_start
     endif
     if ((current(k0)%amp2.eq.0.d0.or.current(k0)%gamma2.eq.0.d0).and.current(k0)%func2.ne.'   ') then
       current(k0)%scale2_start=current(k0)%amp2*funct(current(k0)%func2,0.d0,current(k0)%gamma2,current(k0)%phase2,current(k0)%param2, &
          current(k0)%duty2,current(k0)%moddepth2) +current(k0)%dc2
       current(k0)%func2='   '; current(k0)%scale2_end = current(k0)%scale2_start
     endif

     if (present(scale)) then   !also have scale passed?
       current(k0)%scale1=current(k0)%scale1*scale
       current(k0)%scale2=current(k0)%scale2*scale
!      current(k0)%amp1=current(k0)%amp1*scale
!      current(k0)%amp2=current(k0)%amp2*scale
!      current(k0)%dc1=current(k0)%dc1*scale
!      current(k0)%dc2=current(k0)%dc2*scale
     endif
!    if (current(k0)%scale1_start.eq.current(k0)%scale1_end.and.current(k0)%func1.eq.'   ') then   !make get_appl easier
     if (current(k0)%scale1_start.eq.current(k0)%scale1_end) then   !make get_appl easier
       current(k0)%scale1=current(k0)%scale1*current(k0)%scale1_start
       current(k0)%scale1_end=1.d0; current(k0)%scale1_start=1.d0
     endif
!    if (current(k0)%scale2_start.eq.current(k0)%scale2_end.and.current(k0)%func2.eq.'   ') then   !make get_appl easier
     if (current(k0)%scale2_start.eq.current(k0)%scale2_end) then   !make get_appl easier
       current(k0)%scale2=current(k0)%scale2*current(k0)%scale2_start
       current(k0)%scale2_end=1.d0; current(k0)%scale2_start=1.d0
     endif

   enddo

   ok_p=.true.
   nullify(fld_list)
   return
  end function



  function funct(func,t,gamm,phase,param,duty,moddepth,rot_phi) result (x)
    character(3)::func
    real(8),intent(in)::t,gamm,phase,param,duty,moddepth
    real(8),intent(in),optional::rot_phi(2)
    real(8)::x,y,w
    if (func.eq.'SIN') then
      x=sin(t*gamm+phase)
    elseif (func.eq.'EXP') then
      x=exp(t*gamm+phase)
    elseif (func.eq.'COS') then
      x=cos(t*gamm+phase)
    elseif (func.eq.'NEW') then
!     x=gamm/param
!     y=exp(duty*x); z=exp((duty-1.d0)*x); w=exp(x)
!     x=(t+phase); if (x.lt.0.d0) x = x + gamm * (1+int(abs(x/gamm)))
!     x=mod(x,gamm)
!     if (x.gt.duty*gamm) then  !low level
!       x=exp(-(x-duty*gamm)/param)
!       x=(1.d0-moddepth+w*(moddepth*(2.d0/y-1.d0)-1.d0))*x/(1.d0-w)+(1.d0-moddepth)*(1.d0-x)
!     else  !high level
!       x=exp(-x/param)
!       x=(1.d0+w*(moddepth-1.d0)+moddepth-2.d0*moddepth*y)*x/(1.d0-w)+(1.d0+moddepth)*(1.d0-x)
!     endif

      x=gamm/param
      y=exp(duty*x); w=exp(x)
      x=(t+phase); if (x.lt.0.d0) x = x + gamm * (1+int(abs(x/gamm)))
      x=mod(x,gamm)
      if (x.gt.duty*gamm) then  !low level
        x=exp((gamm-x)/param)
        x=1.d0-moddepth+x*(y-1.d0)*moddepth/(w-1.d0)
      else  !high level
        x=exp(-x/param)
        x=1.d0+x*(y-w)*moddepth/(w-1.d0)
      endif
    elseif (func.eq.'STE') then
      x=t+phase
      if (x.lt.0.d0) x = x + gamm * (1+int(abs(x/gamm)))
      x=mod(x,gamm)
      if (x.lt.0.5d0*param) then
        x=2.d0*x/param
      elseif (x.lt.0.5d0*(gamm-param)) then
        x=1.d0
      elseif (x.lt.0.5d0*(gamm+param)) then
        x=(gamm-2.d0*x)/param
      elseif (x.lt.gamm-0.5d0*param) then
        x=-1.d0
      else
        x=2.d0*(x-gamm)/param
      endif
    elseif (func.eq.'ROT') then
        x=0.d0
        if (present(rot_phi)) then
!print *,t,(/cos(rot_phi+phase*acos(-1.d0)),sin(rot_phi+phase*acos(-1.d0))/)
          if (abs(gamm).ne.1.d0) then
            if (rot_fld((/cos((2.d0*gamm*t+phase)*acos(-1.d0)),sin((2.d0*gamm*t+phase)*acos(-1.d0)),0.d0/))) x=1.d0
          else !if (rot_phi(2).eq.0.d0) then
            y=sqrt(1.d0-rot_phi(2)*rot_phi(2))
            if (rot_fld((/y*cos(rot_phi(1)+phase*acos(-1.d0)),y*sin(rot_phi(1)+phase*acos(-1.d0)),rot_phi(2)/))) x=gamm
!       print *,t,y*cos(rot_phi(1)+phase*acos(-1.d0)),y*sin(rot_phi(1)+phase*acos(-1.d0)),rot_phi(2)
 !          else
 !           y=sqrt(1.d0-rot_phi(2)*rot_phi(2))
 !           if (rot_fld((/y*cos(rot_phi(1)),y*sin(rot_phi(1)),rot_phi(2)/))) x=gamm
          endif
        endif
    else
      x=0.d0
    endif
  end function

  function rot_fld(h) result (ok_p)
   real(8),intent(in)::h(3)
   type(current_applied_s),dimension(:),pointer::current=>null()
   type(applied_s),pointer::f2
   integer::i
   logical::ok_p

   ok_p=.false.

   current=>current_fld
   if (.not.associated(current)) return

   do i=1,size(current)      !loop over layers
     f2=>current(i)%fld2
     if (current(i)%func1.eq.'ROT') current(i)%fld%uniform_h=h 
     if (associated(f2).and.current(i)%func2.eq.'ROT') f2%uniform_h=h
     ok_p=(ok_p.or.(current(i)%func1.eq.'ROT').or.(associated(f2).and.current(i)%func2.eq.'ROT')) !is there a 'ROT' field?
   enddo

  end function
  function uniform_fld(ifld,h) result (ok_p)
   integer,intent(in)::ifld
   real(8),intent(in)::h(3)
   type(applied_s),pointer::fld,fld2,fld3
   logical::ok_p

   ok_p=.false.

   if (.not.associated(all_fld)) then
     allocate(all_fld); fld=>all_fld
   else
     fld2=>all_fld; nullify(fld,fld3)
     do while (associated(fld2))
       if (fld2%index.eq.ifld) fld=>fld2
       fld3=>fld2; fld2=>fld2%next
     enddo
     if (.not.associated(fld)) then
       allocate(fld3%next); fld=>fld3%next
     else
       if (associated(fld%vfld)) deallocate(fld%vfld); nullify(fld%vfld)
     endif
   endif
   fld%index=ifld; fld%nx=1; fld%ny=1 !; fld%last_time=-1.d300
   fld%uniform_h=h
   if (.not.associated(fld%ix)) allocate(fld%xpos(1),fld%ypos(1),fld%ix(nxs),fld%fx(nxs),fld%iy(nys),fld%fy(nys))
   ok_p=.true.
  end function

  function suprgaus_temp(ifld,dtemp,sigmax,sigmay,expn) result (ok_p)
   integer,intent(in)::ifld,expn
   real(8),intent(in)::dtemp,sigmax,sigmay
   type(applied_s),pointer::fld,fld2,fld3
   logical::ok_p

   ok_p=.false.

   if (.not.associated(all_temp)) then
     allocate(all_temp); fld=>all_temp
     fld%offset_x = dxs*nxs*0.5d0
     fld%offset_y = dys*nys*0.5d0
   else
     fld2=>all_temp; nullify(fld,fld3)
     do while (associated(fld2))
       if (fld2%index.eq.ifld) fld=>fld2
       fld3=>fld2; fld2=>fld2%next
     enddo
     if (.not.associated(fld)) then
       allocate(fld3%next); fld=>fld3%next
       fld%offset_x = dxs*nxs*0.5d0
       fld%offset_y = dys*nys*0.5d0
     else
       if (associated(fld%sfld)) deallocate(fld%sfld); nullify(fld%sfld)
     endif
   endif
   fld%index=ifld; fld%analytic=.true.
   fld%nx=1; fld%ny=1; fld%sigmax=0.d0; fld%sigmay=0.d0
   fld%dtemp=dtemp; fld%expn=2*expn !; fld%last_time=-1.d300
   if (sigmax.ne.0.d0) fld%sigmax=1.d0/sigmax
   if (sigmay.ne.0.d0) fld%sigmay=1.d0/sigmay
   if (expn.eq.0) then
     fld%uniform_t=dtemp;fld%nx=0; fld%ny=0; fld%analytic=.false.
   else
     if (.not.associated(fld%ix)) allocate(fld%xpos(1),fld%ypos(1),fld%ix(nxs),fld%fx(nxs),fld%iy(nys),fld%fy(nys))
   endif
   ok_p=.true.
  end function

  function dup_temp(ifld,ifld2,skew,rc) result (ok_p)
   use io_m
   integer,intent(in)::ifld,ifld2
   real(8),intent(in),optional::skew,rc(2)

   real(8)::cskew,sskew,rr(2),dr(2),hw(3),fx,fy
   integer::i0,j0,i,j
   character(200)::os

   type(applied_s),pointer::fld_list=>null(),fld,fld2,fld3,fld4
   logical::ok_p

   ok_p=.false.; if (ifld.eq.ifld2) return
   fld_list=>all_temp

  entry dup_fld(ifld,ifld2,skew,rc) result (ok_p)

   ok_p=.false.; if (ifld.eq.ifld2) return
   if (.not.associated(fld_list)) fld_list=>all_fld

   if (.not.associated(fld_list)) return

   nullify(fld,fld2,fld3); fld4=>fld_list
   do while (associated(fld4))
     if (fld4%index.eq.ifld) fld=>fld4
     if (fld4%index.eq.ifld2) fld2=>fld4
     fld3=>fld4;fld4=>fld4%next
   enddo
   if (.not.associated(fld)) then; nullify(fld_list); return; endif   !source field doesn't exist!

   if (.not.associated(fld2)) then     !target field doesn't exist, create it
     allocate(fld3%next); fld2=>fld3%next
   endif

   fld2%index=ifld2 !; fld2%last_time=fld%last_time
   fld2%nx=fld%nx; fld2%ny=fld%ny; fld2%offset_x=fld%offset_x; fld2%offset_y=fld%offset_y;
   fld2%t_origin=fld%t_origin; fld2%ref_x=fld%ref_x; fld2%ref_y=fld%ref_y; fld2%vel_x=fld%vel_x; fld2%vel_y=fld%vel_y
   fld2%dtemp=fld%dtemp; fld2%sigmax=fld%sigmax; fld2%sigmay=fld%sigmay; fld2%expn=fld%expn; 
   fld2%uniform_h=fld%uniform_h; fld2%uniform_t=fld%uniform_t; fld2%analytic=fld%analytic

   if (associated(fld2%xpos)) then; if (size(fld2%xpos).ne.fld2%nx) then; deallocate(fld2%xpos); nullify(fld2%xpos); endif; endif
   if (.not.associated(fld2%xpos)) allocate(fld2%xpos(fld2%nx)); fld2%xpos=fld%xpos

   if (associated(fld2%ypos)) then; if (size(fld2%ypos).ne.fld2%ny) then; deallocate(fld2%ypos); nullify(fld2%ypos); endif; endif
   if (.not.associated(fld2%ypos)) allocate(fld2%ypos(fld2%ny)); fld2%ypos=fld%ypos

   if (associated(fld%sfld)) then
     if (associated(fld2%sfld)) then
       if (size(fld2%sfld,1).ne.fld2%nx.or.size(fld2%sfld,2).ne.fld2%ny) then; deallocate(fld2%sfld); nullify(fld2%sfld); endif
     endif
     if (.not.associated(fld2%sfld)) allocate(fld2%sfld(fld2%nx,fld2%ny))
     fld2%sfld=fld%sfld
   else; if (associated(fld2%sfld)) then; deallocate(fld2%sfld); nullify(fld2%sfld); endif
   endif

   if (associated(fld%vfld)) then
     if (associated(fld2%vfld)) then
       if (size(fld2%vfld,2).ne.fld2%nx.or.size(fld2%vfld,3).ne.fld2%ny) then; deallocate(fld2%vfld); nullify(fld2%vfld); endif
     endif
     if (.not.associated(fld2%vfld)) allocate(fld2%vfld(3,fld2%nx,fld2%ny))
     fld2%vfld=fld%vfld
   else; if (associated(fld2%vfld)) then; deallocate(fld2%vfld); nullify(fld2%vfld); endif
   endif

   if (associated(fld2%ix)) then; if (size(fld2%ix).ne.nxs) then; deallocate(fld2%ix); nullify(fld2%ix); endif; endif
   if (.not.associated(fld2%ix)) allocate(fld2%ix(nxs)); fld2%ix=fld%ix

   if (associated(fld2%fx)) then; if (size(fld2%fx).ne.nxs) then; deallocate(fld2%fx); nullify(fld2%fx); endif; endif
   if (.not.associated(fld2%fx)) allocate(fld2%fx(nxs)); fld2%fx=fld%fx

   if (associated(fld2%iy)) then; if (size(fld2%iy).ne.nys) then; deallocate(fld2%iy); nullify(fld2%iy); endif; endif
   if (.not.associated(fld2%iy)) allocate(fld2%iy(nys)); fld2%iy=fld%iy

   if (associated(fld2%fy)) then; if (size(fld2%fy).ne.nys) then; deallocate(fld2%fy); nullify(fld2%fy); endif; endif
   if (.not.associated(fld2%fy)) allocate(fld2%fy(nys)); fld2%fy=fld%fy

   cskew=1.d0
   if (present(skew)) cskew=cos(skew)
   if (cskew.ne.1.d0) then
     write(os,"(' duplicating ',i2,' profile to ',i2,' with skew=',1p,e12.5,' and origin=(',e12.5,x,e12.5,')')") ifld,ifld2,skew*180.d0/acos(-1.d0),rc(1),rc(2)
     call output(trim(os))

     sskew=sin(skew)
     rr=(/0.d0, 0.d0/); if (present(rc)) rr=rc
     do j=1,fld%ny
       do i=1,fld%nx  !coordinates in new coordinate: (xpos(i),ypos(j)), in old: (cos(skew)xpos(i)+sin(skew)ypos(J),-sin(skew)xpos(i)+cos(skew)ypos(j))
         dr=(/ (fld%xpos(i)-rr(1))*cskew+(fld%ypos(j)-rr(2))*sskew, -(fld%xpos(i)-rr(1))*sskew+(fld%ypos(j)-rr(2))*cskew /)
         dr=(/ max(fld%xpos(1),min(fld%xpos(fld%nx),dr(1))),max(fld%ypos(1),min(fld%ypos(fld%ny),dr(2))) /)
!        i0=1; do while (fld%xpos(i0+1).lt.dr(1).or.i0.lt.size(fld%ix)-1); i0=i0+1; enddo; fx = (fld%xpos(i0+1)-dr(1))/(fld%xpos(i0+1)-fld%xpos(i0))
!        j0=1; do while (fld%ypos(j0+1).lt.dr(2).or.j0.lt.size(fld%iy)-1); j0=j0+1; enddo; fy = (fld%ypos(j0+1)-dr(2))/(fld%ypos(j0+1)-fld%ypos(j0))
         i0=1; do while (fld%xpos(i0+1).lt.dr(1)); i0=i0+1; enddo; fx = (fld%xpos(i0+1)-dr(1))/(fld%xpos(i0+1)-fld%xpos(i0))
         j0=1; do while (fld%ypos(j0+1).lt.dr(2)); j0=j0+1; enddo; fy = (fld%ypos(j0+1)-dr(2))/(fld%ypos(j0+1)-fld%ypos(j0))
         if (associated(fld2%sfld)) fld2%sfld(i,j)= &
              fy*(fx*fld%sfld(i0,j0)  +(1.d0-fx)*fld%sfld(i0+1,j0))  +(1.d0-fy)*(fx*fld%sfld(i0,j0+1)  +(1.d0-fx)*fld%sfld(i0+1,j0+1))
         if (associated(fld2%vfld)) then
           hw=fy*(fx*fld%vfld(:,i0,j0)+(1.d0-fx)*fld%vfld(:,i0+1,j0))+(1.d0-fy)*(fx*fld%vfld(:,i0,j0+1)+(1.d0-fx)*fld%vfld(:,i0+1,j0+1))
           fld2%vfld(:,i,j)=(/ hw(1)*cskew-hw(2)*sskew, hw(1)*sskew+hw(2)*cskew, hw(3) /)
         endif
       enddo
     enddo

   endif

   ok_p=.true.
   nullify(fld_list)

   return

  end function

  subroutine show_temp()
   use io_m

   type(applied_s),pointer::fld_list=>null(),fld
   character(200)::os,tpe

   fld_list=>all_temp; tpe='temperature'

  entry show_fld()

   if (.not.associated(fld_list)) then
      fld_list=>all_fld; tpe='field'
   endif

   fld=>fld_list
   do while (associated(fld))
     write(os,"(' ',a,' profile index #',i3)") trim(tpe),fld%index
     call output(trim(os))
     write(os,"('   nx=',i5,', ny=',i5,', t_origin=',1p,e12.5,', Vx=',e12.5,', Vy=',e12.5)") fld%nx,fld%ny,fld%t_origin,fld%vel_x,fld%vel_y
     call output(trim(os))
     write(os,"(1p,'   offset_x=',e12.5,', offset_y=',e12.5,', ref_x=',e12.5,', ref_y=',e12.5)") fld%offset_x,fld%offset_y,fld%ref_x,fld%ref_y;
     call output(trim(os))
     if (.not.associated(fld%sfld).and..not.associated(fld%vfld)) then
       if (fld%analytic) then
         write(os,"(1p,'   supra-Gaussian, sigma_x=',e12.5,', sigmay_y=',e12.5,', dT=',e12.5,', exp=',0p,i2)") fld%sigmax,fld%sigmay,fld%dtemp,fld%expn
         call output(trim(os))
       else
         if (associated(fld_list,all_temp)) then
           write(os,"('   spatially uniform ',a,':',1p,e12.5)")trim(tpe),fld%uniform_t
         else
           write(os,"('   spatially uniform ',a,': (',1p,2(e12.5,x),e12.5,')')")trim(tpe),fld%uniform_h
         endif
         call output(trim(os))
       endif
     else
       if (fld%vel_x.eq.0.d0.and.fld%vel_y.eq.0.d0) then
         write(os,"('   non-moving ',a)")trim(tpe)
         call output(trim(os))
       else
         write(os,"('   general (moving and non-uniform) ',a)")trim(tpe)
         call output(trim(os))
       endif
     endif
     fld=>fld%next
   enddo

   nullify(fld_list)

  end subroutine


  function ref_temp(ifld,offset_x,offset_y,t_origin,ref_x,ref_y,vel_x,vel_y,add) result (ok_p)
   integer,intent(in),optional::ifld
   real(8),intent(in),optional::offset_x,offset_y,t_origin,ref_x,ref_y,vel_x,vel_y
   logical,intent(in),optional::add

   type(applied_s),pointer::fld_list=>null(),fld,fld2
   integer::k
   logical::ok_p

   fld_list=>all_temp

  entry ref_fld(ifld,offset_x,offset_y,t_origin,ref_x,ref_y,vel_x,vel_y,add) result (ok_p)

   if (.not.associated(fld_list)) fld_list=>all_fld

   ok_p=.false.

   if (.not.associated(fld_list)) return

   nullify(fld)
   if (present(ifld)) then
     fld2=>fld_list
     do while (associated(fld2))
       if (fld2%index.eq.ifld) fld=>fld2; fld2=>fld2%next
     enddo
     if (.not.associated(fld)) then; nullify(fld_list); return; endif
   endif

   k=0; if (present(add)) then; if (add) k=1; endif

   fld2=>fld_list
   do while (associated(fld2))
     if (.not.associated(fld).or.associated(fld,fld2)) then
       if (present(offset_x)) fld2%offset_x=offset_x+k*fld2%offset_x
       if (present(offset_y)) fld2%offset_y=offset_y+k*fld2%offset_y
       if (present(t_origin)) fld2%t_origin=t_origin+k*fld2%t_origin
       if (present(ref_x)   ) fld2%ref_x   =ref_x   +k*fld2%ref_x
       if (present(ref_y)   ) fld2%ref_y   =ref_y   +k*fld2%ref_y
       if (present(vel_x)   ) fld2%vel_x   =vel_x   +k*fld2%vel_x
       if (present(vel_y)   ) fld2%vel_y   =vel_y   +k*fld2%vel_y
 !       fld2%last_time=-1.d200
     endif
     fld2=>fld2%next
   enddo
   ok_p=.true.
   nullify(fld_list)

   return

  end function


  function file_temp(ifld,ionum,sc) result (ok_p)
   integer,intent(in)::ifld,ionum
   real(8),intent(in)::sc
   logical::ok_p
   ok_p=file_applied(ifld,ionum,all_temp,1,sc)
  end function
  function file_fld(ifld,ionum,sc) result (ok_p)
   integer,intent(in)::ifld,ionum
   real(8),intent(in)::sc
   logical::ok_p
   ok_p=file_applied(ifld,ionum,all_fld,3,sc)
  end function

  function file_applied(ifld,ionum,fld_list,nval,sc) result (ok_p)
!   use graphix, only:draw
   integer,intent(in)::ifld,ionum,nval
   real(8),intent(in)::sc
   type(applied_s),pointer::fld_list,fld,fld2,fld3
   logical::ok_p
   character(300)::os
   integer::i,j
   real(8)::x

   ok_p=.false.

   !let's look to see if this index exists, otherwise create it, fld will point to it
   if (.not.associated(fld_list)) then
     allocate(fld_list); fld=>fld_list
   else
     fld2=>fld_list; nullify(fld,fld3)
     do while (associated(fld2))
       if (fld2%index.eq.ifld) fld=>fld2
       fld3=>fld2; fld2=>fld2%next
     enddo
     if (.not.associated(fld)) then
       allocate(fld3%next); fld=>fld3%next
     else
       fld%analytic=.false.
     endif
   endif
   fld%index=ifld

   !find the size of the data file (nx,ny)
   os='gjp'   !search for the line w/ 'I=' and 'J=', i.e. tecplot file
   do while (index(os,'I=')+index(os,'i=')+index(os,'J=')+index(os,'j=').eq.0)
     read(ionum,'(a)',err=99,end=99) os
   enddo
   os=trim(adjustl(os))
   j=index(os,'I='); if (j.eq.0) j=index(os,'i=')
   do i=j+2,len_trim(os); if (os(i:i).ne.' '.and.(os(i:i).lt.'0'.or.os(i:i).gt.'9')) exit;enddo
   read(os(j+2:min(i-1,len_trim(os))),*) fld%nx
   j=index(os,'J='); if (j.eq.0) j=index(os,'j=')
   do i=j+2,len_trim(os); if (os(i:i).ne.' '.and.(os(i:i).lt.'0'.or.os(i:i).gt.'9')) exit;enddo
   read(os(j+2:min(i-1,len_trim(os))),*) fld%ny

   !now let's find where the data actually starts
   os='gjp';i=1
   do while (i.le.len_trim(os))
     read(ionum,'(a)',err=99,end=99) os
     os=trim(adjustl(os))
     do i=1,len_trim(os); if (os(i:i).ne.' '.and.os(i:i).ne.'e'.and.os(i:i).ne.'E'.and.os(i:).ne.'d'.and.os(i:i).ne.'D' &
            .and.os(i:i).ne.'-'.and.os(i:i).ne.'+'.and.os(i:i).ne.'.'.and.(os(i:i).lt.'0'.or.os(i:i).gt.'9')) exit
     enddo
   enddo
   backspace(ionum)

   if (associated(fld%xpos)) then
     if (size(fld%xpos).ne.fld%nx) then; deallocate(fld%xpos); nullify(fld%xpos); endif
   endif; if (.not.associated(fld%xpos)) allocate(fld%xpos(fld%nx))
   if (associated(fld%ypos)) then
     if (size(fld%ypos).ne.fld%ny) then; deallocate(fld%ypos); nullify(fld%ypos); endif
   endif; if (.not.associated(fld%ypos)) allocate(fld%ypos(fld%ny))
   if (nval.eq.1) then
     if (associated(fld%sfld)) then
       if (size(fld%sfld,1).ne.fld%nx.or.size(fld%sfld,2).ne.fld%ny) then; deallocate(fld%sfld); nullify(fld%sfld); endif
     endif; if (.not.associated(fld%sfld)) allocate(fld%sfld(fld%nx,fld%ny))
   else
     if (associated(fld%vfld)) then
       if (size(fld%vfld,1).ne.nval.or.size(fld%vfld,2).ne.fld%nx.or.size(fld%vfld,3).ne.fld%ny) then; deallocate(fld%vfld); nullify(fld%vfld); endif
     endif; if (.not.associated(fld%vfld)) allocate(fld%vfld(nval,fld%nx,fld%ny))
   endif

   do j=1,fld%ny
     do i=1,fld%nx
       if (nval.eq.1) then
         read(ionum,*,err=99,end=99) fld%xpos(i),fld%ypos(j),fld%sfld(i,j)
       else
!         read(ionum,*,err=99,end=99) fld%xpos(i),x,fld%ypos(j),fld%vfld(:,i,j)
!!         x=fld%vfld(2,i,j); fld%vfld(2,i,j)=fld%vfld(3,i,j); fld%vfld(3,i,j)=-x
!         x=fld%vfld(2,i,j); fld%vfld(2,i,j)=fld%vfld(3,i,j); fld%vfld(3,i,j)=x
         read(ionum,*,err=99,end=99) fld%xpos(fld%nx+1-i),x,fld%ypos(j),fld%vfld(1,fld%nx+1-i,j),fld%vfld(3,fld%nx+1-i,j),fld%vfld(2,fld%nx+1-i,j)
       endif
     enddo
   enddo
   if (nval.ne.1) then
     fld%xpos=-fld%xpos; fld%vfld(1,:,:)=-fld%vfld(1,:,:) !since fem->mrm: x->-x, y->z, z->y
!    do i=1,fld%nx/2
!      x=fld%xpos(i); fld%xpos(i)=-fld%xpos(fld%nx+1-i);fld%xpos(fld%nx+1-i)=-x
!    enddo
!    fld%vfld(1,:,:)=fld%vfld(1,:,:)
     if (sc.ne.1.d0) fld%vfld=fld%vfld*sc
   else
     if (sc.ne.1.d0) fld%sfld=fld%sfld*sc
   endif

 !   fld%last_time=-1.d300
   if (fld%offset_x.lt.-1.d200) fld%offset_x=dxs*nxs*0.5d0
   if (fld%offset_y.lt.-1.d200) fld%offset_y=dys*nys*0.5d0
   if (.not.associated(fld%ix)) allocate(fld%ix(nxs),fld%fx(nxs),fld%iy(nys),fld%fy(nys))

! if (nval.eq.1) call draw(title='temp',x=fld%xpos*1.d7,y=fld%ypos*1.d7,z=fld%sfld)

   ok_p=.true.
  99 continue
   return
  end function

  function start_applied(nx,ny,dx,dy,n_lay) result (ok_p)
   integer,intent(in)::nx,ny,n_lay
   real(8),intent(in)::dx,dy
   logical::ok_p
   nxs=nx; nys=ny; dxs=dx; dys=dy
   nxs1=nxs; nys1=nys
   allocate(current_temp(n_lay),current_fld(n_lay),tempscale(2,n_lay),fldscale(2,n_lay),last_temp(nxs,nys,n_lay),last_fld(3,nxs,nys,n_lay), &
      last_utemp(n_lay),last_ufld(3,n_lay))
   tempscale=0.d0; fldscale=0.d0; last_temp=0.d0; last_fld=0.d0
   ok_p=.true.
  end function

  function get_appl(t,changed,silent_p,rot_phi) result (ok_p)
   use io_m, only: output
   use window_m, only: get_window
   real(8),intent(in)::t
   logical,intent(inout)::changed(:,:)
   logical,intent(in),optional::silent_p
   real(8),intent(in),optional::rot_phi(2)

   logical::s_p,ok_p
   type(current_applied_s),dimension(:),pointer::current=>null()
   integer::k1,k0,i,i00(2)
   real(8)::s1,s2
   real(8),pointer::scalel(:,:)
   real(8),target::gtmp(nxs,nys)
   type(applied_s),pointer::f,f2,f3
   real(8),pointer::f3sfld(:,:)
   real(8),pointer::f3vfld(:,:,:)


   s_p=.false.; if (present(silent_p)) s_p=silent_p
   ok_p=.false.
   call get_window(nxr=i00); nxs0=i00(1); nxs1=i00(2)
   call get_window(nyr=i00); nys0=i00(1); nys1=i00(2)
   do k0=1,2  !temperature and then field 
     if (k0.eq.1) then
       current=>current_temp; scalel=>tempscale
     else
       current=>current_fld; scalel=>fldscale
     endif
     do k1=1,size(current)      !loop over layers
       changed(k1,k0)=(.not.current(k1)%done_p.and.t.ne.current(k1)%last_time)
!print *,k1,k0,t,current(k1)%last_time
       if (changed(k1,k0)) then        !don't do anything if the fields were static in space and time
         f=>current(k1)%fld; f2=>current(k1)%fld2; current(k1)%last_time=t
         last_utemp(k1)=0.d0; last_ufld(:,k1)=0.d0

         if (.not.associated(f2)) then   !one field

          !find temporal scaling
           if (.not.current(k1)%const1_p) then
             if (current(k1)%func1.ne.'   ') then
               scalel(1,k1)=(current(k1)%amp1*funct(current(k1)%func1,t,current(k1)%gamma1,current(k1)%phase1,current(k1)%param1, &
                       current(k1)%duty1,current(k1)%moddepth1,rot_phi)+current(k1)%dc1)*current(k1)%scale1
             else
               scalel(1,k1)=((current(k1)%scale1_end-current(k1)%scale1_start)*max(0.d0,min(1.d0, &
                  (t-current(k1)%t_start)/(current(k1)%t_end-current(k1)%t_start)))+current(k1)%scale1_start)*current(k1)%scale1
             endif
           else
             scalel(1,k1)=current(k1)%scale1
           endif
          !uniform field
           if (current(k1)%uniform1_p) then
             if (k0.eq.1) then
               last_utemp(k1)=scalel(1,k1)*f%uniform_t
             else
               last_ufld(:,k1)=scalel(1,k1)*f%uniform_h
             endif
          !non-moving, time varying
           elseif (associated(current(k1)%snonmov_fld1).or.associated(current(k1)%vnonmov_fld1)) then  !we have a static spatial profile (but time varying)
             if (k0.eq.1) then
               last_temp(:,:,k1)=scalel(1,k1)*current(k1)%snonmov_fld1
             else
               last_fld(:,:,:,k1)=scalel(1,k1)*current(k1)%vnonmov_fld1
             endif
          !moving and (possibly) time varying
           elseif (f%analytic) then !analytic expression for field
!            if (k0.eq.1) then   !we only have analytic expressions for temp
               if (f%vel_x.ne.0.d0) call set_sg(f%expn,f%sigmax,f%ref_x-f%vel_x*(t-f%t_origin)-f%offset_x,nxs,dxs,f%fx)
               if (f%vel_y.ne.0.d0) call set_sg(f%expn,f%sigmay,f%ref_y-f%vel_y*(t-f%t_origin)-f%offset_y,nys,dys,f%fy)
               call get_sg(f%fx,f%fy,last_temp(:,:,k1),f%dtemp*scalel(1,k1))
!            else
!            endif
           else

             if (f%vel_x.ne.0.d0) call set(f%offset_x+f%vel_x*(t-f%t_origin)-f%ref_x,nxs,dxs,f%xpos,f%nx,f%ix,f%fx)
             if (f%vel_y.ne.0.d0) call set(f%offset_y+f%vel_y*(t-f%t_origin)-f%ref_y,nys,dys,f%ypos,f%ny,f%iy,f%fy)
             call interpolate(f%ix,f%fx,f%iy,f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f%sfld,f%vfld,scalel(1,k1))
           endif !type of field

         else  !two fields

          !find temporal scaling
           if (.not.current(k1)%const1_p) then !find temporal scaling for first field
             if (current(k1)%func1.ne.'   ') then   !find temporal scaling
               s1=(current(k1)%amp1*funct(current(k1)%func1,t,current(k1)%gamma1,current(k1)%phase1,current(k1)%param1, &
                       current(k1)%duty1,current(k1)%moddepth1,rot_phi)+current(k1)%dc1)*current(k1)%scale1 &
                       * ((current(k1)%scale1_end-current(k1)%scale1_start)*max(0.d0,min(1.d0, &
                  (t-current(k1)%t_start)/(current(k1)%t_end-current(k1)%t_start)))+current(k1)%scale1_start)
             else
               s1=((current(k1)%scale1_end-current(k1)%scale1_start)*max(0.d0,min(1.d0, &
                  (t-current(k1)%t_start)/(current(k1)%t_end-current(k1)%t_start)))+current(k1)%scale1_start)*current(k1)%scale1
             endif
           else
             s1=current(k1)%scale1
           endif
           if (.not.current(k1)%const2_p) then !find temporal scaling for second field
             if (current(k1)%func2.ne.'   ') then
               s2=(current(k1)%amp2*funct(current(k1)%func2,t,current(k1)%gamma2,current(k1)%phase2,current(k1)%param2, &
                       current(k1)%duty2,current(k1)%moddepth2,rot_phi)+current(k1)%dc2)*current(k1)%scale2 &
                       * ((current(k1)%scale1_end-current(k1)%scale1_start)*max(0.d0,min(1.d0, &
                  (t-current(k1)%t_start)/(current(k1)%t_end-current(k1)%t_start)))+current(k1)%scale1_start)
             else
               s2=((current(k1)%scale2_end-current(k1)%scale2_start)*max(0.d0,min(1.d0, &
                  (t-current(k1)%t_start)/(current(k1)%t_end-current(k1)%t_start)))+current(k1)%scale2_start)*current(k1)%scale2
             endif
           else
             s2=current(k1)%scale2
           endif
           scalel(1,k1)=s1; scalel(2,k1)=s2

           if (current(k1)%uniform1_p.and.current(k1)%uniform2_p) then  ! both uniform?
             if (k0.eq.1) then
               last_utemp(k1)=s1*f%uniform_t+s2*f2%uniform_t
             else
               last_ufld(:,k1)=s1*f%uniform_h+s2*f2%uniform_h
             endif
           elseif (current(k1)%same_p) then !two fields, same mesh or at least one (uniform or analytic) field
             nullify(f3)
             if (current(k1)%uniform1_p) then
               if (k0.eq.1) then; last_utemp(k1)=s1* f%uniform_t; else; last_ufld(:,k1)=s1* f%uniform_h; endif 
               f3=>f2; f3sfld=>current(k1)%snonmov_fld2; f3vfld=>current(k1)%vnonmov_fld2;i=2
             elseif (current(k1)%uniform2_p) then
               if (k0.eq.1) then; last_utemp(k1)=s2*f2%uniform_t; else; last_ufld(:,k1)=s2*f2%uniform_h; endif
               f3=>f;  f3sfld=>current(k1)%snonmov_fld1; f3vfld=>current(k1)%vnonmov_fld1;i=1
             endif
             if (associated(f3)) then ! we had one uniform field and one (non-uniform or analytic) field
               if (associated(f3sfld).or.associated(f3vfld)) then  !we have a static spatial profile (but time varying)
                 if (k0.eq.1) then
                   last_temp(:,:,k1)=scalel(1,k1)*f3sfld
                 else
                   last_fld(:,:,:,k1)=scalel(1,k1)*f3vfld
                 endif
               elseif (f3%analytic) then !analytic expression for field
!                if (k0.eq.1) then   !we only have analytic expressions for temp
                   if (f3%vel_x.ne.0.d0) call set_sg(f3%expn,f3%sigmax,f3%ref_x-f3%vel_x*(t-f3%t_origin)-f3%offset_x,nxs,dxs,f3%fx)
                   if (f3%vel_y.ne.0.d0) call set_sg(f3%expn,f3%sigmay,f3%ref_y-f3%vel_y*(t-f3%t_origin)-f3%offset_y,nys,dys,f3%fy)
                   call get_sg(f3%fx,f3%fy,last_temp(:,:,k1),f3%dtemp*scalel(i,k1))
!                else
!                endif
               else  !moving field
                 if (f3%vel_x.ne.0.d0) call set(f3%offset_x+f3%vel_x*(t-f3%t_origin)-f3%ref_x,nxs,dxs,f3%xpos,f3%nx,f3%ix,f3%fx)
                 if (f3%vel_y.ne.0.d0) call set(f3%offset_y+f3%vel_y*(t-f3%t_origin)-f3%ref_y,nys,dys,f3%ypos,f3%ny,f3%iy,f3%fy)
                 call interpolate(f3%ix,f3%fx,f3%iy,f3%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f3%sfld,f3%vfld,scalel(i,k1))
               endif

            !two fields, neither uniform
             elseif ((associated(current(k1)%snonmov_fld1).and.associated(current(k1)%snonmov_fld2)).or. &
                     (associated(current(k1)%vnonmov_fld1).and.associated(current(k1)%vnonmov_fld2))) then  !we have a static spatial profile (but time varying)
               if (k0.eq.1) then
                 last_temp(:,:,k1)=scalel(1,k1)*current(k1)%snonmov_fld1+scalel(2,k1)*current(k1)%snonmov_fld2
               else
                 last_fld(:,:,:,k1)=scalel(1,k1)*current(k1)%vnonmov_fld1+scalel(2,k1)*current(k1)%vnonmov_fld2
               endif
             else!neither field is uniform, but have same mesh
               if (f%analytic.and.f2%analytic) then
                 if ( f%vel_x.ne.0.d0) call set_sg( f%expn, f%sigmax, f%ref_x- f%vel_x*(t- f%t_origin)- f%offset_x,nxs,dxs, f%fx)
                 if ( f%vel_y.ne.0.d0) call set_sg( f%expn, f%sigmay, f%ref_y- f%vel_y*(t- f%t_origin)- f%offset_y,nxs,dys, f%fy)
                 call get_sg( f%fx, f%fy,gtmp, f%dtemp*scalel(1,k1))
                 if (f2%vel_x.ne.0.d0) call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%vel_x*(t-f2%t_origin)-f2%offset_x,nxs,dxs,f2%fx)
                 if (f2%vel_y.ne.0.d0) call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%vel_y*(t-f2%t_origin)-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,last_temp(:,:,k1),f2%dtemp*scalel(2,k1))
                 last_temp(:,:,k1)=last_temp(:,:,k1)+gtmp
               elseif (f%analytic) then
                 if ( f%vel_x.ne.0.d0) call set_sg( f%expn, f%sigmax, f%ref_x- f%vel_x*(t- f%t_origin)- f%offset_x,nxs,dxs, f%fx)
                 if ( f%vel_y.ne.0.d0) call set_sg( f%expn, f%sigmay, f%ref_y- f%vel_y*(t- f%t_origin)- f%offset_y,nys,dys, f%fy)
                 call get_sg( f%fx, f%fy,gtmp, f%dtemp*scalel(1,k1))
                 if (f2%vel_x.ne.0.d0) call set(f2%offset_x+f2%vel_x*(t-f2%t_origin)-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                 if (f2%vel_y.ne.0.d0) call set(f2%offset_y+f2%vel_y*(t-f2%t_origin)-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f2%sfld,f2%vfld,scalel(2,k1))
                 last_temp(:,:,k1)=last_temp(:,:,k1)+gtmp
               elseif (f2%analytic) then
                 if (f2%vel_x.ne.0.d0) call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%vel_x*(t-f2%t_origin)-f2%offset_x,nxs,dxs,f2%fx)
                 if (f2%vel_y.ne.0.d0) call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%vel_y*(t-f2%t_origin)-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,gtmp,f2%dtemp*scalel(2,k1))
                 if ( f%vel_x.ne.0.d0) call set( f%offset_x+ f%vel_x*(t- f%t_origin)- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 if ( f%vel_y.ne.0.d0) call set( f%offset_y+ f%vel_y*(t- f%t_origin)- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate(f%ix, f%fx, f%iy, f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1), f%sfld, f%vfld,scalel(1,k1))
                 last_temp(:,:,k1)=last_temp(:,:,k1)+gtmp
               else
                 if ( f%vel_x.ne.0.d0) call set( f%offset_x+ f%vel_x*(t- f%t_origin)- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 if ( f%vel_y.ne.0.d0) call set( f%offset_y+ f%vel_y*(t- f%t_origin)- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate2(f%ix,f%fx,f%iy,f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f%sfld,f2%sfld,f%vfld,f2%vfld,s1,s2)
               endif !type of field
             endif !type of field
          !two fields, neither uniform, not same mesh
           else !two fields, neither uniform, not same mesh
             if ((associated(current(k1)%snonmov_fld1).and.associated(current(k1)%snonmov_fld2)).or. &
                     (associated(current(k1)%vnonmov_fld1).and.associated(current(k1)%vnonmov_fld2))) then  !we have a static spatial profile (but time varying)
               if (k0.eq.1) then
                 last_temp(:,:,k1)=scalel(1,k1)*current(k1)%snonmov_fld1+scalel(2,k1)*current(k1)%snonmov_fld2
               else
                 last_fld(:,:,:,k1)=scalel(1,k1)*current(k1)%vnonmov_fld1+scalel(2,k1)*current(k1)%vnonmov_fld2
               endif
             else
               if ( f%vel_x.ne.0.d0) call set( f%offset_x- f%vel_x*(t- f%t_origin)- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
               if ( f%vel_y.ne.0.d0) call set( f%offset_y- f%vel_y*(t- f%t_origin)- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
               if (f2%vel_x.ne.0.d0) call set(f2%offset_x-f2%vel_x*(t-f2%t_origin)-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
               if (f2%vel_y.ne.0.d0) call set(f2%offset_y-f2%vel_y*(t-f2%t_origin)-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
               call interpolate2d(f%ix,f%fx,f%iy,f%fy,f2%ix,f2%fx,f2%iy,f2%fy,last_temp(:,:,k1),last_fld(:,:,:,k1), &
                       f%sfld,f2%sfld,f%vfld,f2%vfld,s1,s2)
             endif
           endif
         endif
       endif  !didn't have to do anything for this layer
     enddo
   enddo
   ok_p=.true.
  end function

  function init_appl(t,changed,silent_p,vel) result (ok_p)
   use io_m, only: output
   use window_m, only: get_window
   real(8),intent(in)::t
   logical,intent(inout)::changed(:,:)
   logical,intent(in),optional::silent_p
   real(8),intent(out),optional::vel(2)

   logical::s_p,ok_p
   type(current_applied_s),dimension(:),pointer::current=>null()
   integer::k1,k0,i00(2)
   real(8)::s1,s2
   character(300)::os,tpe
   real(8),pointer::scalel(:,:)
   real(8),target::gtmp(nxs,nys),g2tmp(3,nxs,nys)
   type(applied_s),pointer::f,f2

   s_p=.false.; if (present(silent_p)) s_p=silent_p
   ok_p=.false.; last_utemp=0.d0; last_ufld=0.d0
   last_temp=0.d0; last_fld=0.d0; if (present(vel)) vel=0.d0
   call get_window(nxr=i00); nxs0=i00(1); nxs1=i00(2)
   call get_window(nyr=i00); nys0=i00(1); nys1=i00(2)
   do k0=1,2
     if (k0.eq.1) then
       current=>current_temp; tpe='temperature'; scalel=>tempscale
     else
       current=>current_fld; tpe='field'; scalel=>fldscale
     endif; scalel=0.d0; !scalel(2,:)=-1.d300
     do k1=1,size(current)
!      if (k0.eq.1.and..not.s_p) then; write(os,"(' for layer #',i2,':')") k1; call output(trim(os)); endif
       write(os,"(' for layer #',i2,':')") k1
       f=>current(k1)%fld; f2=>current(k1)%fld2; current(k1)%last_time=-1.d300
       changed(k1,k0)=(.not.associated(f))
       current(k1)%done_p=changed(k1,k0)
       if (current(k1)%done_p) then
          write(os,"(a,' no ',a,' profile')")trim(os),trim(tpe)
          if (.not.s_p) call output(trim(os))
       elseif (associated(f2)) then
          write(os,"(a,' applied ',a,' profile indices',i3,' and',i3)")trim(os),trim(tpe),f%index,f2%index
!         if (.not.s_p) call output(trim(os))
          if (k0.eq.1) then
            current(k1)%uniform1_p=(.not.associated( f%sfld).and. f%expn.eq.0)
            current(k1)%uniform2_p=(.not.associated(f2%sfld).and.f2%expn.eq.0)
          else
            current(k1)%uniform1_p=(.not.associated( f%vfld))
            current(k1)%uniform2_p=(.not.associated(f2%vfld))
          endif
          if (current(k1)%uniform1_p.or.current(k1)%uniform2_p) then
            current(k1)%same_p=.true.
          elseif (f%analytic.or.f2%analytic) then
            current(k1)%same_p= &
                    (f2%offset_x.eq.f%offset_x.and.f2%offset_y.eq.f%offset_y.and. &
                     f2%t_origin.eq.f%t_origin.and.f2%ref_x.eq.f%ref_x.and.f2%ref_y.eq.f%ref_y.and.f2%vel_x.eq.f%vel_x.and. &
                       f2%vel_y.eq.f%vel_y)
          else
            current(k1)%same_p= & !both uniform fields or same mesh?
                    (f2%nx.eq.f%nx.and.f2%ny.eq.f%ny.and.f2%offset_x.eq.f%offset_x.and.f2%offset_y.eq.f%offset_y.and. &
                     f2%t_origin.eq.f%t_origin.and.f2%ref_x.eq.f%ref_x.and.f2%ref_y.eq.f%ref_y.and.f2%vel_x.eq.f%vel_x.and. &
                       f2%vel_y.eq.f%vel_y.and.all(f2%xpos.eq.f%xpos).and.all(f2%ypos.eq.f%ypos))
!           if (current(k1)%same_p) then
!              write(os,"(' ',a,' indices ',i2,' and ',i2,' use the same mesh (good job!)')") trim(tpe),f%index,f2%index
!              if (.not.s_p) call output(trim(os))
!           endif
          endif
          if (.not.current(k1)%uniform1_p) then; vel=(/f%vel_x,f%vel_y/)
          elseif (.not.current(k1)%uniform2_p) then; vel=(/f2%vel_x,f2%vel_y/)
          endif
       else
          write(os,"(a,' applied ',a,' profile index',i3)")trim(os),trim(tpe),f%index
!         if (.not.s_p) call output(trim(os))
          if (k0.eq.1) then
            current(k1)%uniform1_p =(.not.associated( f%sfld).and. f%expn.eq.0)
          else
            current(k1)%uniform1_p =(.not.associated( f%vfld))
          endif
          if (.not.current(k1)%uniform1_p.and.present(vel)) vel=(/f%vel_x,f%vel_y/)
       endif

       if (.not.current(k1)%done_p) then
         current(k1)%const1_p=(current(k1)%scale1_start.eq.current(k1)%scale1_end.and.current(k1)%func1.eq.'   ')
         current(k1)%const2_p=(current(k1)%scale2_start.eq.current(k1)%scale2_end.and.current(k1)%func2.eq.'   ')

         if (.not.associated(f2)) then   !one field
           s1=current(k1)%scale1*current(k1)%scale1_start
           if (((f%vel_x.eq.0.d0.and.f%vel_y.eq.0.d0).or.current(k1)%uniform1_p.or.s1.eq.0.d0).and.current(k1)%const1_p) then
          !(non-moving or uniform field) and non-time varying
             current(k1)%done_p=.true.; changed(k1,k0)=.true.
             scalel(1,k1)=s1
             if (s1.eq.0.d0) then
               os=trim(os)//', zero scaling (no need to update)'; if (.not.s_p) call output(trim(os))
             elseif (current(k1)%uniform1_p) then  !uniform field and non-time varying
               if (k0.eq.1) then
                 last_utemp(k1)=s1*f%uniform_t
               else
                 last_ufld(:,k1)=s1*f%uniform_h
               endif
               os=trim(os)//', constant and uniform (no need to update)'; if (.not.s_p) call output(trim(os))
             elseif (f%analytic) then !non-moving, non-time varying analytic field
!              if (k0.eq.1) then  !only analytic temperature
                 call set_sg(f%expn,f%sigmax,f%ref_x-f%offset_x,nxs,dxs,f%fx)
                 call set_sg(f%expn,f%sigmay,f%ref_y-f%offset_y,nys,dys,f%fy)
                 call get_sg(f%fx,f%fy,last_temp(:,:,k1),f%dtemp*s1)
!              else
!              endif
               os=trim(os)//', constant and analytic (no need to update)'; if (.not.s_p) call output(trim(os))
             else  !non-moving, non-time varying profile
               call set(f%offset_x-f%ref_x,nxs,dxs,f%xpos,f%nx,f%ix,f%fx)
               call set(f%offset_y-f%ref_y,nys,dys,f%ypos,f%ny,f%iy,f%fy)
               call interpolate(f%ix,f%fx,f%iy,f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f%sfld,f%vfld,s1)
               os=trim(os)//', constant and fixed profile (no need to update)'; if (.not.s_p) call output(trim(os))
             endif
           elseif ((f%vel_x.eq.0.d0.and.f%vel_y.eq.0.d0).or.current(k1)%uniform1_p) then
          !(non-moving or uniform field) and time-varying
             if (.not.current(k1)%uniform1_p) then 
            !non-moving profile that is time varying
               if (f%analytic) then  !then non-moving analytic field
!                if (k0.eq.1) then  !only analytic temperature
                   if (.not.associated(current(k1)%snonmov_fld1)) allocate(current(k1)%snonmov_fld1(nxs,nys))
                   call set_sg(f%expn,f%sigmax,f%ref_x-f%offset_x,nxs,dxs,f%fx)
                   call set_sg(f%expn,f%sigmay,f%ref_y-f%offset_y,nys,dys,f%fy)
                   call get_sg(f%fx,f%fy,current(k1)%snonmov_fld1,f%dtemp)
!                else
!                endif
                 os=trim(os)//', varying and non-moving analytic profile (need to update)'; if (.not.s_p) call output(trim(os))
               else !non-moving and time_varying profile
                 call set(f%offset_x-f%ref_x,nxs,dxs,f%xpos,f%nx,f%ix,f%fx)
                 call set(f%offset_y-f%ref_y,nys,dys,f%ypos,f%ny,f%iy,f%fy)
                 call interpolate(f%ix,f%fx,f%iy,f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f%sfld,f%vfld,1.d0)
                 if (k0.eq.1) then
                   if (.not.associated(current(k1)%snonmov_fld1)) allocate(current(k1)%snonmov_fld1(nxs,nys))
                   current(k1)%snonmov_fld1=last_temp(:,:,k1); last_temp(:,:,k1)=0.d0
                 else
                   if (.not.associated(current(k1)%vnonmov_fld1)) allocate(current(k1)%vnonmov_fld1(3,nxs,nys))
                   current(k1)%vnonmov_fld1=last_fld(:,:,:,k1); last_fld(:,:,:,k1)=0.d0
                 endif
                 os=trim(os)//', varying and non-moving profile (need to update)'; if (.not.s_p) call output(trim(os))
               endif
             else
               os=trim(os)//', varying and uniform (need to update)'; if (.not.s_p) call output(trim(os))
             endif
           else
             if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1)
             if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1)
             nullify(current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1)
          !moving (and possibly time-varying)
             if (f%analytic) then
               call set_sg(f%expn,f%sigmax,f%ref_x-f%offset_x,nxs,dxs,f%fx)
               call set_sg(f%expn,f%sigmay,f%ref_y-f%offset_y,nys,dys,f%fy)
               if (current(k1)%const1_p) then
                 os=trim(os)//', constant but moving analytic profile (need to update)'
               else
                 os=trim(os)//', varying and moving analytic profile (need to update)' 
               endif
               if (.not.s_p) call output(trim(os))
             else
               call set(f%offset_x-f%ref_x,nxs,dxs,f%xpos,f%nx,f%ix,f%fx)
               call set(f%offset_y-f%ref_y,nys,dys,f%ypos,f%ny,f%iy,f%fy)
               if (current(k1)%const1_p) then
                 os=trim(os)//', constant but moving profile (need to update)'
               else
                 os=trim(os)//', varying and moving profile (need to update)' 
               endif
               if (.not.s_p) call output(trim(os))
             endif
           endif
         else  !two fields
           s1=current(k1)%scale1*current(k1)%scale1_start
           s2=current(k1)%scale2*current(k1)%scale2_start
           if (current(k1)%same_p) then
             if (current(k1)%const1_p.and.current(k1)%const2_p) then     !constant in time
               if (current(k1)%uniform1_p.and.current(k1)%uniform2_p) then  !both uniform
                 current(k1)%done_p=.true.; changed(k1,k0)=.true.
                 scalel(1,k1)=s1; scalel(2,k1)=s2
                 if (k0.eq.1) then
                   last_utemp(k1)=s1*f%uniform_t+s2*f2%uniform_t
                 else
                   last_ufld(:,k1)=s1*f%uniform_h+s2*f2%uniform_h
                 endif
                 os=trim(os)//', both constant and uniform profile (no need to update)'; if (.not.s_p) call output(trim(os))
               elseif (current(k1)%uniform1_p.and.f2%vel_x.eq.0.d0.and.f2%vel_y.eq.0.d0) then !constant in time, non-moving
                 current(k1)%done_p=.true.; changed(k1,k0)=.true.
                 scalel(1,k1)=s1; scalel(2,k1)=s2
                 if (k0.eq.1) then
                   last_utemp(k1)=s1*f%uniform_t
                 else
                   last_ufld(:,k1)=s1*f%uniform_h
                 endif
                 if (f2%analytic) then
                   call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                   call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                   call get_sg(f2%fx,f2%fy,last_temp(:,:,k1),f2%dtemp*scalel(2,k1))
                   os=trim(os)//', both constant and one uniform and one fixed analytic profile (no need to update)'; if (.not.s_p) call output(trim(os))
                 else
                   call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                   call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                   call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,last_temp(:,:,k1),last_fld(:,:,:,k1),f2%sfld,f2%vfld,scalel(2,k1))
                   os=trim(os)//', both constant and one uniform and one fixed (no need to update)'; if (.not.s_p) call output(trim(os))
                 endif
               elseif (current(k1)%uniform2_p.and. f%vel_x.eq.0.d0.and. f%vel_y.eq.0.d0) then !constant in time, non-moving
                 current(k1)%done_p=.true.; changed(k1,k0)=.true.
                 scalel(1,k1)=s1; scalel(2,k1)=s2
                 if (k0.eq.1) then
                   last_utemp(k1)=s2*f2%uniform_t
                 else
                   last_ufld(:,k1)=s1*f%uniform_h+s2*f2%uniform_h
                 endif
                 if ( f%analytic) then
                   call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                   call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                   call get_sg( f%fx, f%fy,last_temp(:,:,k1), f%dtemp*scalel(1,k1))
                   os=trim(os)//', both constant and one fixed analytic and one uniform profile (no need to update)'; if (.not.s_p) call output(trim(os))
                 else
                   call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                   call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                   call interpolate(f%ix, f%fx, f%iy, f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1), f%sfld, f%vfld,scalel(1,k1))
                   os=trim(os)//', both constant and one fixed and one uniform (no need to update)'; if (.not.s_p) call output(trim(os))
                 endif
               elseif (f%vel_x.eq.0.d0.and. f%vel_y.eq.0.d0) then !constant in time, non-moving
                 current(k1)%done_p=.true.; changed(k1,k0)=.true.
                 scalel(1,k1)=s1; scalel(2,k1)=s2
                 if ( f%analytic) then
                   call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                   call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                   call get_sg( f%fx, f%fy,last_temp(:,:,k1), f%dtemp*scalel(1,k1))
                 else
                   call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                   call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                   call interpolate( f%ix, f%fx, f%iy, f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1), f%sfld, f%vfld,scalel(1,k1))
                 endif
                 if (f2%analytic) then
                   call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                   call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                   call get_sg(f2%fx,f2%fy,gtmp,f2%dtemp*scalel(2,k1))
                 else
                   call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                   call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                   call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,gtmp,g2tmp,f2%sfld,f2%vfld,scalel(2,k1))
                 endif
                 if (k0.eq.1) then
                    last_temp(:,:,k1)=last_temp(:,:,k1)+gtmp
                 else
                    last_fld(:,:,:,k1)=last_fld(:,:,:,k1)+g2tmp
                 endif
                 os=trim(os)//', both constant and fixed pofiles (no need to update)'; if (.not.s_p) call output(trim(os))
               else
                 if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1)
                 if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1)
                 if (associated(current(k1)%snonmov_fld2)) deallocate(current(k1)%snonmov_fld2)
                 if (associated(current(k1)%vnonmov_fld2)) deallocate(current(k1)%vnonmov_fld2)
                 nullify(current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2)
                 if (f%analytic) then
                   call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx); call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                 else
                   call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx); call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 endif
                 if (f2%analytic) then
                   call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx); call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                 else
                   call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx); call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 endif
                 os=trim(os)//', both constant and at least one moving pofile (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
               endif
             elseif (current(k1)%uniform1_p.and.current(k1)%uniform2_p) then !varying in time, uniform
               if (k0.eq.1) then
                 if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1); nullify(current(k1)%snonmov_fld1)
                 if (associated(current(k1)%snonmov_fld2)) deallocate(current(k1)%snonmov_fld2); nullify(current(k1)%snonmov_fld2)
               else
                 if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1); nullify(current(k1)%vnonmov_fld1)
                 if (associated(current(k1)%vnonmov_fld2)) deallocate(current(k1)%vnonmov_fld2); nullify(current(k1)%vnonmov_fld2)
               endif
               os=trim(os)//', both uniform and at least one varying (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
             elseif (current(k1)%uniform1_p.and.f2%vel_x.eq.0.d0.and.f2%vel_y.eq.0.d0) then !varying in time, non-moving
               if (k0.eq.1) then
                 if (.not.associated(current(k1)%snonmov_fld2)) allocate(current(k1)%snonmov_fld2(nxs,nys))
                 if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1); nullify(current(k1)%snonmov_fld1)
               else
                 if (.not.associated(current(k1)%vnonmov_fld2)) allocate(current(k1)%vnonmov_fld2(3,nxs,nys))
                 if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1); nullify(current(k1)%vnonmov_fld1)
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                 call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,current(k1)%snonmov_fld2,f2%dtemp)
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                 call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2,f2%sfld,f2%vfld,1.d0)
               endif
               os=trim(os)//', at least one varying (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
             elseif (current(k1)%uniform2_p.and. f%vel_x.eq.0.d0.and. f%vel_y.eq.0.d0) then !varying in time, non-moving
               if (k0.eq.1) then
                 if (.not.associated(current(k1)%snonmov_fld1)) allocate(current(k1)%snonmov_fld1(nxs,nys))
                 if (associated(current(k1)%snonmov_fld2)) deallocate(current(k1)%snonmov_fld2); nullify(current(k1)%snonmov_fld2)
               else
                 if (.not.associated(current(k1)%vnonmov_fld1)) allocate(current(k1)%vnonmov_fld1(3,nxs,nys))
                 if (associated(current(k1)%vnonmov_fld2)) deallocate(current(k1)%vnonmov_fld2); nullify(current(k1)%vnonmov_fld2)
               endif
               if ( f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                 call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                 call get_sg( f%fx, f%fy,current(k1)%snonmov_fld1, f%dtemp)
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate( f%ix, f%fx, f%iy, f%fy,current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1, f%sfld, f%vfld,1.d0)
               endif
               os=trim(os)//', at least one varying (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
             elseif (f%vel_x.eq.0.d0.and. f%vel_y.eq.0.d0) then !varying in time, non-moving
               if (k0.eq.1) then
                 if (.not.associated(current(k1)%snonmov_fld1)) allocate(current(k1)%snonmov_fld1(nxs,nys))
                 if (.not.associated(current(k1)%snonmov_fld2)) allocate(current(k1)%snonmov_fld2(nxs,nys))
               else
                 if (.not.associated(current(k1)%vnonmov_fld1)) allocate(current(k1)%vnonmov_fld1(3,nxs,nys))
                 if (.not.associated(current(k1)%vnonmov_fld2)) allocate(current(k1)%vnonmov_fld2(3,nxs,nys))
               endif
               if ( f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                 call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                 call get_sg( f%fx, f%fy,current(k1)%snonmov_fld1, f%dtemp)
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate( f%ix, f%fx, f%iy, f%fy,current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1, f%sfld, f%vfld,1.d0)
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                 call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,current(k1)%snonmov_fld2,f2%dtemp)
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                 call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2,f2%sfld,f2%vfld,1.d0)
               endif
               os=trim(os)//', at least one varying and fixed (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
             else
               if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1)
               if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1)
               if (associated(current(k1)%snonmov_fld2)) deallocate(current(k1)%snonmov_fld2)
               if (associated(current(k1)%vnonmov_fld2)) deallocate(current(k1)%vnonmov_fld2)
               nullify(current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2)
               if (f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx); call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx); call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx); call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx); call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
               endif
               os=trim(os)//', both varying and moving (need to update, same mesh)'; if (.not.s_p) call output(trim(os))
             endif
           else !not same mesh
             if (f%vel_x.eq.0.d0.and.f%vel_y.eq.0.d0.and.f2%vel_x.eq.0.d0.and.f2%vel_y.eq.0.d0.and.current(k1)%const1_p.and.current(k1)%const2_p) then
               current(k1)%done_p=.true.; changed(k1,k0)=.true.
               scalel(1,k1)=s1; scalel(2,k1)=s2
               if ( f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                 call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                 call get_sg( f%fx, f%fy,last_temp(:,:,k1), f%dtemp*scalel(1,k1))
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate( f%ix, f%fx, f%iy, f%fy,last_temp(:,:,k1),last_fld(:,:,:,k1), f%sfld, f%vfld,scalel(1,k1))
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                 call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,gtmp,f2%dtemp*scalel(2,k1))
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                 call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,gtmp,g2tmp,f2%sfld,f2%vfld,scalel(2,k1))
               endif
               if (k0.eq.1) then
                  last_temp(:,:,k1)=last_temp(:,:,k1)+gtmp
               else
                  last_fld(:,:,:,k1)=last_fld(:,:,:,k1)+g2tmp
               endif
               os=trim(os)//', both constant and fixed (no need to update, not same mesh)'; if (.not.s_p) call output(trim(os))
             elseif (f%vel_x.eq.0.d0.and.f%vel_y.eq.0.d0.and.f2%vel_x.eq.0.d0.and.f2%vel_y.eq.0.d0) then
               if (k0.eq.1) then
                 if (.not.associated(current(k1)%snonmov_fld1)) allocate(current(k1)%snonmov_fld1(nxs,nys))
                 if (.not.associated(current(k1)%snonmov_fld2)) allocate(current(k1)%snonmov_fld2(nxs,nys))
               else
                 if (.not.associated(current(k1)%vnonmov_fld1)) allocate(current(k1)%vnonmov_fld1(3,nxs,nys))
                 if (.not.associated(current(k1)%vnonmov_fld2)) allocate(current(k1)%vnonmov_fld2(3,nxs,nys))
               endif
               if ( f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx)
                 call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
                 call get_sg( f%fx, f%fy,current(k1)%snonmov_fld1, f%dtemp)
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx)
                 call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
                 call interpolate( f%ix, f%fx, f%iy, f%fy,current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1, f%sfld, f%vfld,1.d0)
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx)
                 call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
                 call get_sg(f2%fx,f2%fy,current(k1)%snonmov_fld2,f2%dtemp)
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx)
                 call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
                 call interpolate(f2%ix,f2%fx,f2%iy,f2%fy,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2,f2%sfld,f2%vfld,1.d0)
               endif
               os=trim(os)//', at least one varying and both fixed (need to update, not same mesh)'; if (.not.s_p) call output(trim(os))
             else
               if (associated(current(k1)%snonmov_fld1)) deallocate(current(k1)%snonmov_fld1)
               if (associated(current(k1)%vnonmov_fld1)) deallocate(current(k1)%vnonmov_fld1)
               if (associated(current(k1)%snonmov_fld2)) deallocate(current(k1)%snonmov_fld2)
               if (associated(current(k1)%vnonmov_fld2)) deallocate(current(k1)%vnonmov_fld2)
               nullify(current(k1)%snonmov_fld1,current(k1)%vnonmov_fld1,current(k1)%snonmov_fld2,current(k1)%vnonmov_fld2)
               if (f%analytic) then
                 call set_sg( f%expn, f%sigmax, f%ref_x- f%offset_x,nxs,dxs, f%fx); call set_sg( f%expn, f%sigmay, f%ref_y- f%offset_y,nys,dys, f%fy)
               else
                 call set( f%offset_x- f%ref_x,nxs,dxs, f%xpos, f%nx, f%ix, f%fx); call set( f%offset_y- f%ref_y,nys,dys, f%ypos, f%ny, f%iy, f%fy)
               endif
               if (f2%analytic) then
                 call set_sg(f2%expn,f2%sigmax,f2%ref_x-f2%offset_x,nxs,dxs,f2%fx); call set_sg(f2%expn,f2%sigmay,f2%ref_y-f2%offset_y,nys,dys,f2%fy)
               else
                 call set(f2%offset_x-f2%ref_x,nxs,dxs,f2%xpos,f2%nx,f2%ix,f2%fx); call set(f2%offset_y-f2%ref_y,nys,dys,f2%ypos,f2%ny,f2%iy,f2%fy)
               endif
               os=trim(os)//', at least one varying and moving (need to update, not same mesh)'; if (.not.s_p) call output(trim(os))
             endif
           endif
         endif !one or two fields
       endif !not done_p
     enddo ! layer
   enddo !temp and field

   ok_p=.true.
   
  end function

  subroutine interpolate(ix,fx,iy,fy,sf,vf,sfi,vfi,s)
   integer,intent(in)::ix(nxs),iy(nys)
   real(8),intent(in)::fx(nxs),fy(nys),s
   real(8),intent(out)::sf(nxs,nys),vf(3,nxs,nys)
   real(8),dimension(:,:),pointer::sfi
   real(8),dimension(:,:,:),pointer::vfi

   integer::i,j
   real(8)::yf,yfc

!print *,'inter:',associated(sfi),associated(vfi)
   if (associated(sfi)) then
!    do j=1,nys      !loop over simulation mesh
!$omp parallel
!$omp do private(j,i,yf,yfc)
     do j=nys0,nys1
       yf = s*fy(j); yfc = s*(1.d0-fy(j))   !fraction from nearest nodes of write field
!      do i=1,nxs
       do i=nxs0,nxs1
         sf(i,j) = fx(i)*yf *sfi(ix(i),iy(j)  ) + (1.d0-fx(i))*yf *sfi(ix(i)+1,iy(j)  ) + &
                   fx(i)*yfc*sfi(ix(i),iy(j)+1) + (1.d0-fx(i))*yfc*sfi(ix(i)+1,iy(j)+1) !+ t0
       enddo
     enddo
!$omp end do nowait
!$omp end parallel
     if (nys0.gt.1) sf(:,nys0-1)=0.d0; if (nys1.lt.nys) sf(:,nys1+1)=0.d0
     if (nxs0.gt.1) sf(nxs0-1,:)=0.d0; if (nxs1.lt.nxs) sf(nxs1+1,:)=0.d0
   else
!    do j=1,nys      !loop over simulation mesh
!$omp parallel
!$omp do private(j,i,yf,yfc)
     do j=nys0,nys1      !loop over simulation mesh
       yf = s*fy(j); yfc = s*(1.d0-fy(j))   !fraction from nearest nodes of write field
!      do i=1,nxs
       do i=nxs0,nxs1
!        f1 = fx(i)*yf; f2=(1.d0-fx(i))*yf; f3=fx(i)*yfc; f4=(1.d0-fx(i))*yfc
!        vf(1,i,j) = f1*vfi(1,ix(i),iy(j)  ) + f2*vfi(1,ix(i)+1,iy(j)  ) + &
!                    f3*vfi(1,ix(i),iy(j)+1) + f4*vfi(1,ix(i)+1,iy(j)+1)
!        vf(2,i,j) = f1*vfi(2,ix(i),iy(j)  ) + f2*vfi(2,ix(i)+1,iy(j)) + &
!                    f3*vfi(2,ix(i),iy(j)+1) + f4*vfi(2,ix(i)+1,iy(j)+1)
!        vf(3,i,j) = f1*vfi(3,ix(i),iy(j)  ) + f2*vfi(3,ix(i)+1,iy(j)) + &
!                    f3*vfi(3,ix(i),iy(j)+1) + f4*vfi(3,ix(i)+1,iy(j)+1)
         vf(:,i,j) = fx(i)*yf *vfi(:,ix(i),iy(j)  ) + (1.d0-fx(i))*yf *vfi(:,ix(i)+1,iy(j)  ) + &
                     fx(i)*yfc*vfi(:,ix(i),iy(j)+1) + (1.d0-fx(i))*yfc*vfi(:,ix(i)+1,iy(j)+1)
       enddo
     enddo
!$omp end do nowait
!$omp end parallel
     if (nys0.gt.1) vf(:,:,nys0-1)=0.d0; if (nys1.lt.nys) vf(:,:,nys1+1)=0.d0
     if (nxs0.gt.1) vf(:,nxs0-1,:)=0.d0; if (nxs1.lt.nxs) vf(:,nxs1+1,:)=0.d0
   endif
  end subroutine
!  subroutine interpolate_g(k,npx,npy,ix,fx,iy,fy,sf,vf,sfi,vfi,s)
!   use data_m, only: get_data
!   integer,intent(in)::npx,npy,ix(nxs),iy(nys),k
!   real(8),intent(in)::fx(nxs),fy(nys),s
!   real(8),intent(out)::sf(:),vf(:,:)
!   real(8),dimension(:,:),pointer::sfi
!   real(8),dimension(:,:,:),pointer::vfi
!
!   integer::n,nn,i,j
!   integer,pointer::min_n,max_n,grain2gridindex(:),grain2grid(:,:)
!   logical,pointer::do_grain(:)
!
!   call get_data(k,min_n,min_n_p=.true.)
!   call get_data(k,max_n,max_n_p=.true.)
!   call get_data(k,grain2gridindex,grain2gridindex_p=.true.)
!   call get_data(k,grain2grid,grain2grid_p=.true.)
!
!   if (associated(sfi)) then
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         sf(n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              sf(n)=s*(      fy(j) *(fx(i)*sfi(ix(i),iy(j)  )+(1.d0-fx(i))*sfi(ix(i)+1,iy(j)  ))+ &
!                       (1.d0-fy(j))*(fx(i)*sfi(ix(i),iy(j)+1)+(1.d0-fx(i))*sfi(ix(i)+1,iy(j)+1)))
!         enddo
!         sf(n)=sf(n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   else
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         vf(:,n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              vf(:,n)=s*(      fy(j) *(fx(i)*vfi(:,ix(i),iy(j)  )+(1.d0-fx(i))*vfi(:,ix(i)+1,iy(j)  ))+ &
!                         (1.d0-fy(j))*(fx(i)*vfi(:,ix(i),iy(j)+1)+(1.d0-fx(i))*vfi(:,ix(i)+1,iy(j)+1)))
!         enddo
!         vf(:,n)=vf(:,n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   endif
!  end subroutine interpolate_g


  subroutine interpolate2(ix,fx,iy,fy,sf,vf,sfi1,sfi2,vfi1,vfi2,s1,s2)
   integer,intent(in)::ix(nxs),iy(nys)
   real(8),intent(in)::fx(nxs),fy(nys),s1,s2
   real(8),intent(out)::sf(nxs,nys),vf(3,nxs,nys)
   real(8),dimension(:,:),pointer::sfi1(:,:),sfi2(:,:)
   real(8),dimension(:,:,:),pointer::vfi1(:,:,:),vfi2(:,:,:)

   integer::i,j
   real(8)::yf,yfc,xf,xfc

   if (associated(sfi1)) then
!$omp parallel
!$omp do private(j,i,yf,yfc,xf,xfc)
!    do j=1,nys
     do j=nys0,nys1
       yf = fy(j); yfc=1.d0-yf
!      do i=1,nxs
       do i=nxs0,nxs1
         xf = fx(i); xfc=1.d0-xf
         sf(i,j) = xf *yf *(s1*sfi1(ix(i)  ,iy(j)  )+s2*sfi2(ix(i)  ,iy(j)  )) + &
                   xfc*yf *(s1*sfi1(ix(i)+1,iy(j)  )+s2*sfi2(ix(i)+1,iy(j)  )) + &
                   xf *yfc*(s1*sfi1(ix(i)  ,iy(j)+1)+s2*sfi2(ix(i)  ,iy(j)+1)) + &
                   xfc*yfc*(s1*sfi1(ix(i)+1,iy(j)+1)+s2*sfi2(ix(i)+1,iy(j)+1)) !+ t0
       enddo
     enddo
!$omp end do nowait
!$omp end parallel
     if (nys0.gt.1) sf(:,nys0-1)=0.d0; if (nys1.lt.nys) sf(:,nys1+1)=0.d0
     if (nxs0.gt.1) sf(nxs0-1,:)=0.d0; if (nxs1.lt.nxs) sf(nxs1+1,:)=0.d0
   else
!$omp parallel
!$omp do private(j,i,yf,yfc,xf,xfc)
!    do j=1,nys
     do j=nys0,nys1
       yf = fy(j); yfc=1.d0-yf
!      do i=1,nxs
       do i=nxs0,nxs1
         xf = fx(i); xfc=1.d0-xf
         vf(:,i,j) = xf *yf *(s1*vfi1(:,ix(i)  ,iy(j)  )+s2*vfi2(:,ix(i)  ,iy(j)  )) + &
                     xfc*yf *(s1*vfi1(:,ix(i)+1,iy(j)  )+s2*vfi2(:,ix(i)+1,iy(j)  )) + &
                     xf *yfc*(s1*vfi1(:,ix(i)  ,iy(j)+1)+s2*vfi2(:,ix(i)  ,iy(j)+1)) + &
                     xfc*yfc*(s1*vfi1(:,ix(i)+1,iy(j)+1)+s2*vfi2(:,ix(i)+1,iy(j)+1))
       enddo
     enddo
!$omp end do nowait
!$omp end parallel
     if (nys0.gt.1) vf(:,:,nys0-1)=0.d0; if (nys1.lt.nys) vf(:,:,nys1+1)=0.d0
     if (nxs0.gt.1) vf(:,nxs0-1,:)=0.d0; if (nxs1.lt.nxs) vf(:,nxs1+1,:)=0.d0
   endif
   return
  end subroutine interpolate2
!  subroutine interpolate2_g(k,npx,npy,ix,fx,iy,fy,sf,vf,sfi1,sfi2,vfi1,vfi2,s1,s2)
!   use data_m, only: get_data
!   integer,intent(in)::npx,npy,ix(nxs),iy(nys),k
!   real(8),intent(in)::fx(nxs),fy(nys),s1,s2
!   real(8),intent(out)::sf(:),vf(:,:)
!   real(8),dimension(:,:),pointer::sfi1(:,:),sfi2(:,:)
!   real(8),dimension(:,:,:),pointer::vfi1(:,:,:),vfi2(:,:,:)
!
!   integer::n,nn,i,j
!   integer,pointer::min_n,max_n,grain2gridindex(:),grain2grid(:,:)
!   logical,pointer::do_grain(:)
!
!   call get_data(k,min_n,min_n_p=.true.)
!   call get_data(k,max_n,max_n_p=.true.)
!   call get_data(k,grain2gridindex,grain2gridindex_p=.true.)
!   call get_data(k,grain2grid,grain2grid_p=.true.)
!
!   if (associated(sfi1)) then
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         sf(n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              sf(n)=         fy(j) *(      fx(i) *(s1*sfi1(ix(i)  ,iy(j)  )+s2*sfi2(ix(i)  ,iy(j)  ))+  &
!                                     (1.d0-fx(i))*(s1*sfi1(ix(i)+1,iy(j)  )+s2*sfi2(ix(i)+1,iy(j)  )))+ &
!                       (1.d0-fy(j))*(      fx(i) *(s1*sfi1(ix(i),  iy(j)+1)+s2*sfi2(ix(i)  ,iy(j)+1))+  &
!                                     (1.d0-fx(i))*(s1*sfi1(ix(i)+1,iy(j)+1)+s2*sfi2(ix(i)+1,iy(j)+1)))  
!         enddo
!         sf(n)=sf(n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   else
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         vf(:,n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              vf(:,n)=         fy(j) *(      fx(i) *(s1*vfi1(:,ix(i)  ,iy(j)  )+s2*vfi2(:,ix(i)  ,iy(j)  ))+  &
!                                       (1.d0-fx(i))*(s1*vfi1(:,ix(i)+1,iy(j)  )+s2*vfi2(:,ix(i)+1,iy(j)  )))+ &
!                         (1.d0-fy(j))*(      fx(i) *(s1*vfi1(:,ix(i),  iy(j)+1)+s2*vfi2(:,ix(i)  ,iy(j)+1))+  &
!                                       (1.d0-fx(i))*(s1*vfi1(:,ix(i)+1,iy(j)+1)+s2*vfi2(:,ix(i)+1,iy(j)+1)))
!         enddo
!         vf(:,n)=vf(:,n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   endif
!  end subroutine interpolate2_g

  subroutine interpolate2d(ix,fx,iy,fy,ix2,fx2,iy2,fy2,sf,vf,sfi1,sfi2,vfi1,vfi2,s1,s2)
   integer,intent(in)::ix(nxs),iy(nys),ix2(nxs),iy2(nys)
   real(8),intent(in)::fx(nxs),fy(nys),s1,s2,fx2(nxs),fy2(nys)
   real(8),intent(out)::sf(nxs,nys),vf(3,nxs,nys)
   real(8),pointer::sfi1(:,:),sfi2(:,:),vfi1(:,:,:),vfi2(:,:,:)

   integer::i,j
   real(8)::yf,yfc,xf,xfc,yf2,yfc2,xf2,xfc2

   if (associated(sfi1)) then
!    do j=1,nys
     do j=nys0,nys1
       yf = s1*fy(j); yfc=s1*(1.d0-fy(j)); yf2 = s2*fy2(j); yfc2=s2*(1.d0-fy2(j))
!      do i=1,nxs
       do i=nxs0,nxs1
         xf = fx(i); xfc=1.d0-xf; xf2 = fx2(i); xfc2=1.d0-xf
         sf(i,j) = xf *yf *sfi1(ix(i)  ,iy(j)  )+xf2 *yf2 *sfi2(ix2(i)  ,iy2(j)  ) + &
                   xfc*yf *sfi1(ix(i)+1,iy(j)  )+xfc2*yf2 *sfi2(ix2(i)+1,iy2(j)  ) + &
                   xf *yfc*sfi1(ix(i)  ,iy(j)+1)+xf2 *yfc2*sfi2(ix2(i)  ,iy2(j)+1) + &
                   xfc*yfc*sfi1(ix(i)+1,iy(j)+1)+xfc2*yfc2*sfi2(ix2(i)+1,iy2(j)+1) !+ t0
       enddo
     enddo
     if (nys0.gt.1) sf(:,nys0-1)=0.d0; if (nys1.lt.nys) sf(:,nys1+1)=0.d0
     if (nxs0.gt.1) sf(nxs0-1,:)=0.d0; if (nxs1.lt.nxs) sf(nxs1+1,:)=0.d0
   else
!    do j=1,nys
     do j=nys0,nys1
       yf = s1*fy(j); yfc=s1*(1.d0-fy(j)); yf2 = s2*fy2(j); yfc2=s2*(1.d0-fy2(j))
!      do i=1,nxs
       do i=nxs0,nxs1
         xf = fx(i); xfc=1.d0-xf; xf2 = fx2(i); xfc2=1.d0-xf
         vf(:,i,j) = xf *yf *vfi1(:,ix(i)  ,iy(j)  )+xf2 *yf2 *vfi2(:,ix2(i)  ,iy2(j)  ) + &
                     xfc*yf *vfi1(:,ix(i)+1,iy(j)  )+xfc2*yf2 *vfi2(:,ix2(i)+1,iy2(j)  ) + &
                     xf *yfc*vfi1(:,ix(i)  ,iy(j)+1)+xf2 *yfc2*vfi2(:,ix2(i)  ,iy2(j)+1) + &
                     xfc*yfc*vfi1(:,ix(i)+1,iy(j)+1)+xfc2*yfc2*vfi2(:,ix2(i)+1,iy2(j)+1)
       enddo
     enddo
     if (nys0.gt.1) vf(:,:,nys0-1)=0.d0; if (nys1.lt.nys) vf(:,:,nys1+1)=0.d0
     if (nxs0.gt.1) vf(:,nxs0-1,:)=0.d0; if (nxs1.lt.nxs) vf(:,nxs1+1,:)=0.d0
   endif
   return
  end subroutine interpolate2d

!  subroutine interpolate2d_g(k,npx,npy,npx2,npy2,ix,fx,iy,fy,ix2,fx2,iy2,fy2,sf,vf,sfi1,sfi2,vfi1,vfi2,s1,s2)
!   use data_m, only: get_data
!   integer,intent(in)::npx,npy,npx2,npy2,ix(nxs),iy(nys),ix2(nxs),iy2(nys),k
!   real(8),intent(in)::fx(nxs),fy(nys),s1,s2,fx2(nxs),fy2(nys)
!   real(8),intent(out)::sf(:),vf(:,:)
!   real(8),pointer::sfi1(:,:),sfi2(:,:),vfi1(:,:,:),vfi2(:,:,:)
!
!   integer::n,nn,i,j
!   integer,pointer::min_n,max_n,grain2gridindex(:),grain2grid(:,:)
!   logical,pointer::do_grain(:)
!
!   call get_data(k,min_n,min_n_p=.true.)
!   call get_data(k,max_n,max_n_p=.true.)
!   call get_data(k,grain2gridindex,grain2gridindex_p=.true.)
!   call get_data(k,grain2grid,grain2grid_p=.true.)
!
!   if (associated(sfi1)) then
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         sf(n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              sf(n)=         fy(j) *(      fx(i) *(s1*sfi1(ix(i)  ,iy(j)  )+s2*sfi2(ix2(i)  ,iy2(j)  ))+  &
!                                     (1.d0-fx(i))*(s1*sfi1(ix(i)+1,iy(j)  )+s2*sfi2(ix2(i)+1,iy2(j)  )))+ &
!                       (1.d0-fy(j))*(      fx(i) *(s1*sfi1(ix(i),  iy(j)+1)+s2*sfi2(ix2(i)  ,iy2(j)+1))+  &
!                                     (1.d0-fx(i))*(s1*sfi1(ix(i)+1,iy(j)+1)+s2*sfi2(ix2(i)+1,iy2(j)+1)))
!         enddo
!         sf(n)=sf(n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   else
!     do n=min_n,max_n
!       if (do_grain(n)) then
!         vf(:,n)=0.d0
!         do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!           if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) &
!              vf(:,n)=         fy(j) *(      fx(i) *(s1*vfi1(:,ix(i)  ,iy(j)  )+s2*vfi2(:,ix2(i)  ,iy2(j)  ))+  &
!                                       (1.d0-fx(i))*(s1*vfi1(:,ix(i)+1,iy(j)  )+s2*vfi2(:,ix2(i)+1,iy2(j)  )))+ &
!                         (1.d0-fy(j))*(      fx(i) *(s1*vfi1(:,ix(i),  iy(j)+1)+s2*vfi2(:,ix2(i)  ,iy2(j)+1))+  &
!                                       (1.d0-fx(i))*(s1*vfi1(:,ix(i)+1,iy(j)+1)+s2*vfi2(:,ix2(i)+1,iy2(j)+1)))
!         enddo
!         vf(:,n)=vf(:,n)/(grain2gridindex(n)-grain2gridindex(n-1))
!       endif
!     enddo
!   endif
!  end subroutine interpolate2d_g
!
  subroutine set(offset,n,dd,pos,npos,ind,f)
   integer,intent(in)::n,npos
   real(8),intent(in)::offset,dd,pos(npos)
   integer,intent(inout)::ind(n)
   real(8),intent(inout)::f(n)

   integer::i,i0
   real(8)::x

   i0 = 2
   do i=1,n
     x = (i-.5d0)*dd - offset
     do while (pos(i0).le.x.and.i0.lt.npos)
       i0=i0+1
     enddo
    !interpolation: A(i) = f(i)*Ai(ind(i)) + (1.d0-f(i))*Ai(ind(i)+1)
     ind(i) = min(npos-1,i0-1)
     f(i) = max(0.d0,min(1.d0,(pos(i0)-x)/(pos(i0)-pos(i0-1)) ))
   enddo
  end subroutine

  subroutine set_sg(expn,sig,del,n,dd,fx)
 !subroutine to set super-gaussian temperature profile onto simulation sub-grain
 !mesh. only do nx+ny exponetials...
   integer,intent(in)::n,expn
   real(8),intent(in)::del,dd,sig
   real(8),intent(out)::fx(n)

   integer::i

   do i=1,n
    fx(i)=exp( - (((dd*(i-0.5d0)+del)*sig)**expn) )
   enddo
  end subroutine set_sg

  subroutine get_sg(fx,fy,happ,s)
 !subroutine to set super-gaussian temperature profile onto simulation sub-grain
 !mesh. only do nx+ny exponetials...
   real(8),intent(in)::fx(nxs),fy(nys),s
   real(8),intent(out)::happ(nxs,nys)

   integer::i,j
   real(8)::x

!  do j=1,nys
   do j=nys0,nys1
     x=s*fy(j)
!    do i=1,nxs
     do i=nxs0,nxs1
       happ(i,j)=x*fx(i) !+t0
     enddo
   enddo
   if (nys0.gt.1) happ(:,nys0-1)=0.d0; if (nys1.lt.nys) happ(:,nys1+1)=0.d0
   if (nxs0.gt.1) happ(nxs0-1,:)=0.d0; if (nxs1.lt.nxs) happ(nxs1+1,:)=0.d0
  end subroutine get_sg

!  subroutine get_sg_g(k,fx,fy,happ,s)
!   use data_m, only: get_data
!   integer,intent(in)::k
!   real(8),intent(in)::fx(nxs),fy(nys),s
!   real(8),intent(out)::happ(:)
!
!   integer::n,nn,i,j
!   integer,pointer::min_n,max_n,grain2gridindex(:),grain2grid(:,:)
!   logical,pointer::do_grain(:)
!
!   call get_data(k,min_n,min_n_p=.true.)
!   call get_data(k,max_n,max_n_p=.true.)
!   call get_data(k,grain2gridindex,grain2gridindex_p=.true.)
!   call get_data(k,grain2grid,grain2grid_p=.true.)
!   do n=min_n,max_n
!     if (do_grain(n)) then
!       happ(n)=0.d0
!       do nn=grain2gridindex(n-1)+1,grain2gridindex(n); i=grain2grid(1,nn); j=grain2grid(2,nn)
!         if (i.ge.nxs0.and.i.le.nxs1.and.j.ge.nys0.and.j.le.nys1) happ(n)=happ(n)+s*fy(j)*fx(i)
!       enddo
!       happ(n)=happ(n)/(grain2gridindex(n)-grain2gridindex(n-1))
!     endif
!   enddo
!  end subroutine get_sg_g


end module applied_m
