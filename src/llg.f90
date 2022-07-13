
module llg_m
 use FP_m, only: FP_s
 implicit none
 private
 public :: start_llg,set_llg,init_llg,llg,show_llg,backup_llg,restore_llg,need_to_initialize_llg,plot_llg

 type::llgp_s
   logical:: mt=.false.,meq_bs=.false.,zhang=.false.,llb=.false.,fp=.false.,fixed_dt_p=.false.,garanin_tau=.true.
   integer::iter=100, plot_iter=100, thermal_num_ave=5001, type=2, min_it=10
   real(8)::dt=1.d-14, dt_max=1.d-8, t_max=1.d22, t_min=-1.d22, del_t=1.d22, max_dmdt=1.d0, dt_frac=0.2d0, max_dmdtnext=1.d0, errmx=0.d0, errmx1=0.d0
   real(8)::mave(3)=0.d0,mvar(3)=0.d0,h_var=0.d0,mmag=0.d0,have(3)=0.d0
 end type llgp_s
 type(llgp_s),pointer,save::llgp=>null(),llgpb=>null()

 type::llg_s
   integer,pointer::in,min_n=>null(),max_n=>null()
   logical,pointer::do_grain(:)=>null()
   real(8),pointer::alpha=>null(),gamma=>null(),m(:,:)=>null(),h(:,:)=>null(),hthm(:,:)=>null(),mag(:,:)=>null(),dm(:,:)=>null(),grain_vol(:)=>null(),msc(:,:)=>null()
   real(8),pointer::alpha_scale(:)=>null(),ms_scale(:)=>null(),ms(:)=>null(),lambdaMs(:)=>null(),spin=>null(),temp(:)=>null(),tc(:)=>null(),hthm_m(:)=>null()
   type(FP_s),pointer::fp(:)=>null()
   real(8)::max_dm=3.1415925d0,min_dm=0.d0,err=0.005d0,err1=0.d0,minm(3)=-1.d22*(/1.d0,1.d0,1.d0/),maxm(3)=1.d22*(/1.d0,1.d0,1.d0/),scale=0.d0,floor_ms=0.01d0
   real(8)::mave(3)=0.d0,mvar(3)=0.d0,mmag=0.d0,have(3)=0.d0
   logical::errl=.true.,err1l=.true.,mave_min,mave_max
 end type llg_s

 type(llg_s),pointer,save::llgl(:)=>null(),llglb(:)=>null()
 logical,save::need_to_initialize=.true.

 type::llg_sc
   type(llg_s),pointer::llgl(:)=>null()
   type(llgp_s),pointer::llgp=>null()
   type(llg_sc),pointer::prev=>null()
 end type llg_sc

 type(llg_sc),pointer,save::llgb=>null()


 
 integer,save::NTRAJ=0
 integer,save::iter_plot=0   !hack to keep track of time steps for correct timing of plotting/dumping 
!integer,save::NTRAJ=400 !313 !400 !329 !592
!integer,save::NTRAJ=26 !124 !0 !63 !8 !568 !63 !362 !326 !169 !307 !28 !52 !224 !643

contains

  subroutine need_to_initialize_llg()
   need_to_initialize=.true.
  end subroutine
  subroutine backup_llg()
    type(llg_sc),pointer::llgbb
    if (.not.associated(llgb)) then   !first time? create linked list of parameters
            allocate(llgb); llgb%llgl=>llgl; llgb%llgp=>llgp; llgb%prev=>null()
    endif
    allocate(llgbb); llgbb%prev=>llgb; llgb=>llgbb    !allocate new parameter block, link this block to first
    allocate(llgb%llgp); llgb%llgp=llgp; llgp=>llgb%llgp  !copy over llgp and link llgp to new copy
    if (associated(llgl)) then
      allocate(llgb%llgl(size(llgl))); llgb%llgl=llgl; llgl=>llgb%llgl  !copy over llgl if it exists, and link llgl to new copy
    endif
!   llgpb=llgp
!   if (allocated(llgl)) then
!     if (.not.allocated(llglb)) allocate(llglb(size(llgl)))
!     llglb=llgl
!   endif
  end subroutine
  subroutine restore_llg()
    type(llg_sc),pointer::llgbb
    if (associated(llgb)) then    !we can only restore if we have something to restore from!
      llgbb=>llgb; llgb=>llgb%prev   !get link to previous parameter block we're going to restore to
      if (associated(llgb)) then  !we can only restore if previous block exists
        llgp=>llgb%llgp; deallocate(llgbb%llgp); llgl=>null()  !link llgp to old, can free up current llgp
        if (associated(llgb%llgl)) llgl=>llgb%llgl             !link llgl to old if present
        if (associated(llgbb%llgl)) deallocate(llgbb%llgl)     !free up llgl if exists
        deallocate(llgbb)                                      !free up this parameter block
      endif
    endif
!   llgp=llgpb
!   if (allocated(llglb)) then
!     llgl=llglb
!   endif
  end subroutine

  subroutine show_llg(time)
   use io_m, only: output
   real(8),intent(in)::time
   logical::ok_p
   integer::k
   character(200)::os

   write(os,"(' mt=',l1,'zhang=',l1,', llb=',l1,', fp=',l1)") llgp%mt,llgp%zhang,llgp%llb,llgp%fp; call output(os)
   write(os,"(' garanin_tau=',l1,', fixed_dt=',l1,', use Bs for meq=',l1)") llgp%garanin_tau,llgp%fixed_dt_p,llgp%meq_bs; call output(os)
   write(os,"(' type=',i2,', iter=',i6,', plot_iter=',i6,', thermal_num_ave=',i6)") llgp%type,llgp%iter,llgp%plot_iter,llgp%thermal_num_ave; call output(os)
   write(os,"(1p,' dt=',e12.5,', dt_max=',e12.5,', dt_frac=',e12.5,', max_dmdt=',e12.5)") llgp%dt,llgp%dt_max,llgp%dt_frac,llgp%max_dmdt; call output(os)
   write(os,"(1p,' t_max=',e12.5,', t_min=',e12.5,', del_t=',e12.5,', time=',e12.5)") llgp%t_max,llgp%t_min,llgp%del_t,time; call output(os)
!  if (allocated(llgl)) then
   if (associated(llgl)) then
     ok_p=init_llg()
     do k=1,size(llgl)
       write(os,"('  for layer #',i2)")k; call output(os)
       write(os,"(1p,'   min_dm=',e12.5,', max_dm=',e12.5,', floor_ms=',e12.5)")llgl(k)%min_dm*180.d0/acos(-1.d0),llgl(k)%max_dm*180.d0/acos(-1.d0), &
                  llgl(k)%floor_ms; call output(os)
       write(os,"(1p,'   max|mxH/Ms|=',e12.5,', max dm=',e12.5,', alpha=',e12.5)")llgl(k)%err,llgl(k)%err1,llgl(k)%alpha; call output(os)
       write(os,"(1p,'   minm=(',e12.5,', ',e12.5,', ',e12.5,')')")llgl(k)%minm; call output(os)
       write(os,"(1p,'   maxm=(',e12.5,', ',e12.5,', ',e12.5,')')")llgl(k)%maxm; call output(os)
     enddo
   endif

  end subroutine


  subroutine llg(t,silent_p,track,trace_p,fe_p,diff_m)
   use io_m, only: output
   use data_m, only: update_do_grain
   use plot_m, only: get_plot
   use applied_m, only: get_applied_scale
   use field_m, only: init_applied
   use track_m, only: track_s, do_tracking

   use graphix, only: draw

   real(8),intent(inout)::t
   logical,intent(in),optional::silent_p,trace_p,fe_p,diff_m
   type(track_s),pointer,optional::track
   
   logical::silent,ok_p,cont_p,vary_p(2)
   integer::n_lay,num_ave,num_ave_i,i,it,n,k
   real(8)::pi,del_t,x,x1,x3
   real(8),allocatable::ave_tq(:,:,:)
   character(500)::os,os1
   real(8),pointer::hscale(:,:),tscale(:,:)
!  integer::ng
   real(8)::msdata(3,size(llgl)+1,2)

   if (.not.get_plot(write_plot=.true.)) iter_plot=0
   silent=init_llg()
   silent=.false.; if (present(silent_p)) silent=silent_p
   pi=acos(-1.d0); n_lay=size(llgl)
   call update_do_grain(t)
   ok_p=init_applied(t,silent,.true.)
   allocate(ave_tq(3,0:llgp%thermal_num_ave,n_lay)); ave_tq=0.d0; num_ave=0; num_ave_i=0
   call get_applied_scale(tscale,hscale)

!  do i=1,size(llgl)
!    call plot_contour('mz_grain',llgl(i)%m(3+llgl(i)%in,:),.true.,i,no_write=.true.)
!    ng=minloc(llgl(i)%m(3+llgl(i)%in,:))
!    print *,ng,llgl(i)%m(3+llgl(i)%in,ng)
!  enddo; read *

   if (.not.llgp%zhang.and..not.llgp%llb.and..not.llgp%fp.and..not.llgp%mt) then
     write(os,"(a)") 'using llg'
   elseif (llgp%llb) then
     write(os,"(a)") 'using llb'
     if (llgp%type.lt.-3) llgp%type=-3
   elseif (llgp%fp) then
     write(os,"(a)") 'using fp'
   else
     if (llgp%zhang) then
       write(os,"(a)") 'using zhang'
       if (llgp%meq_bs) then
         write(os,"(a)") trim(os)//' with Bs(x,s)'
       else
         write(os,"(a)") trim(os)//' with L(x)'
       endif
     elseif (llgp%mt) then
       write(os,"(a)") 'using mt'
       if (llgp%type.lt.-3) llgp%type=-3
       if (llgp%meq_bs) then
         write(os,"(a)") trim(os)//' with Bs(x,s)'
       else
         write(os,"(a)") trim(os)//' with L(x)'
       endif
     endif
     if (llgp%garanin_tau) then
       write(os,"(a)") trim(os)//' and garanin''s tau'
     else
       write(os,"(a)") trim(os)//' and zhang''s tau'
     endif
   endif
   if (abs(llgp%type).eq.1) then
     write(os,"(a)") trim(os)//', rotations with quaternions'
     if (llgp%type.lt.0) write(os,"(a)") trim(os)//' mid-point'
   elseif (abs(llgp%type).eq.2) then
     if (llgp%zhang.or.llgp%mt) then
       write(os,"(a)") trim(os)//', analytic decay and explicit rotation'
     else
       write(os,"(a)") trim(os)//', rotations with boris'
     endif
     if (llgp%type.lt.0) write(os,"(a)") trim(os)//' mid-point'
   elseif (llgp%type.eq.-3) then
     write(os,"(a)") trim(os)//', Heun'
   elseif (llgp%type.eq.-4) then
     if (llgp%zhang) then
       write(os,"(a)") trim(os)//', old LLG rotation- don''t use!'
     else
       write(os,"(a)") trim(os)//', Heun in LLG'
     endif
   else
           write(os,"(a)") ' *** abort '//trim(os)//': unknown type ***'; call output(trim(os)); return
   endif
   if (.not.silent) call output(trim(os))

   write(os,"('   iter     mxH/Ms      theta          dt         delta t      t            dt_m         H_thml       <|m|>        scale(L1)')")
   do i=2,n_lay; write(os1,"(i3)")i; os=trim(os)//'    scale(L'//trim(adjustl(trim(os1)))//')'; enddo
   if (.not.silent) call output(trim(os))

   llgp%max_dmdt=1.d0 !which one? FIXME
   if (llgp%dt_frac.ne.0.d0) llgp%max_dmdt=1.d0
   it=0; del_t=0.d0

   cont_p=.false.; do i=1,n_lay; cont_p=(cont_p.or.any(llgl(i)%do_grain)); enddo

   do while (cont_p)
     !vary_p=get_fields(t)

     if (need_to_initialize.and.llgp%fp) call initialize_fp()

     llgp%mvar=0.d0; llgp%mave=0.d0; x=0.d0; llgp%have=0.d0
     do i=1,n_lay
!       llgl(i)%mave=0.d0; x=0.d0
        llgl(i)%mave=0.d0; x1=0.d0
        llgl(i)%mvar=0.d0; llgl(i)%have=0.d0
        do n=1,size(llgl(i)%grain_vol)
!          x=x+llgl(i)%grain_vol(n); llgl(i)%mave=llgl(i)%mave+llgl(i)%grain_vol(n)*(llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,n)-llgl(i)%mave)/x
           x3=llgl(i)%grain_vol(n)/(x1+llgl(i)%grain_vol(n))  ! w_k/W_k
           llgl(i)%mvar=llgl(i)%mvar+x3*x1*(llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,n)-llgl(i)%mave)*(llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,n)-llgl(i)%mave)
           x1=x1+llgl(i)%grain_vol(n)
           llgl(i)%mave=llgl(i)%mave+x3*(llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,n)-llgl(i)%mave)
           llgl(i)%have=llgl(i)%have+x3*(llgl(i)%h(:,n)-llgl(i)%have)
        enddo
        llgl(i)%mvar=llgl(i)%mvar/x1
        x=x+x1
        llgp%mave=llgp%mave+x1/x*(llgl(i)%mave-llgp%mave)
        llgp%have=llgp%have+x1/x*(llgl(i)%have-llgp%have)
        llgp%mvar=llgp%mvar+x1/x*(llgl(i)%mvar-llgp%mvar)
!       x1=x1+x; llgp%mave=llgp%mave+x*llgl(i)%mave/x1
        llgl(i)%mave_min=( any(llgl(i)%mave.le.llgl(i)%minm) ); llgl(i)%mave_max=( any(llgl(i)%mave.ge.llgl(i)%maxm) )
     enddo
     if (all(llgl(:)%mave_min.or.llgl(:)%mave_max).and..not.vary_p(1)) then; cont_p=.false.; exit; endif

     vary_p=get_fields(t)
     call push_m(t)
     
!    if (NTRAJ.gt.0.and.llgl(1)%do_grain(max(1,NTRAJ)).and.(present(trace_p).or.present(track)).and.mod(it,10).eq.0) then
     if (NTRAJ.gt.0.and.(present(trace_p).or.present(track)).and.mod(it,10).eq.0) then
       ok_p=.false.; do i=1,n_lay; ok_p=ok_p.or. llgl(i)%do_grain(NTRAJ); msdata(:,i,2)=msdata(:,i,1); msdata(:,i,1)=llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,NTRAJ);  enddo
       msdata(:,n_lay+1,2)=msdata(:,n_lay+1,1); msdata(:,n_lay+1,1)=(/ sum(msdata(1,1:n_lay,1)), sum(msdata(2,1:n_lay,1)), sum(msdata(3,1:n_lay,1)) /)/n_lay

!      msdata(:,n_lay+1,1)=0.d0
!      do i=1,n_lay 
!        msdata(:,n_lay+1,1)=msdata(:,n_lay+1,1)+msdata(:,i,1)*llgl(i)%ms(NTRAJ)*llgl(i)%ms_scale(NTRAJ)*llgl(i)%grain_vol(NTRAJ)
!      enddo
!      msdata(:,n_lay+1,1)=msdata(:,n_lay+1,1)/sqrt(dot_product(msdata(:,n_lay+1,1),msdata(:,n_lay+1,1)))
!      call draw('M trace',scatter=msdata(:,n_lay+1:n_lay+1,1),scatter2=msdata(:,n_lay+1:n_lay+1,2),add=(it.gt.0))

!      i=size(msdata,2); if (n_lay.eq.1) i=i-1
       if (ok_p) then
       i=size(llgl)+min(1,n_lay-1)
       if (it.ne.0) then
         call draw('M trace',scatter=msdata(:,1:i,1),scatter2=msdata(:,1:i,2),add=.true.)
         call draw('M',      scatter=msdata(:,1:i,1),scatter2=msdata(:,1:i,2),add=.false.)
       else
         call draw('M trace',scatter=msdata(:,1:i,1),scatter2=msdata(:,1:i,2))
         call draw('M',      scatter=msdata(:,1:i,1),scatter2=msdata(:,1:i,2))
       endif
       endif
     endif

     llgp%dt = llgp%dt * llgp%max_dmdt; llgp%max_dmdt =llgp%max_dmdtnext
     t = t + llgp%dt; del_t = del_t + llgp%dt; it=it+1; 

     if (.not.silent.and.mod(it,llgp%iter).eq.0) then
       write(os,"(i8,1p,:,1310(' ',e12.5))") it,llgp%errmx,acos(1.d0-.5d0*llgp%errmx1)*180.d0/pi, &
                   llgp%dt,del_t,t,llgp%max_dmdt,llgp%h_var,llgp%mmag/n_lay,hscale(1,1:n_lay)
       call output(trim(os))
     endif

     cont_p= ((((any(llgl(:)%errl).and.any(llgl(:)%err1l)).or.vary_p(1)).and.del_t.lt.llgp%del_t.and.&
          t.lt.max(llgp%t_max, llgp%t_min)) .or. it.lt.llgp%min_it)
     if (present(track)) cont_p=(cont_p.and.do_tracking(track,t,llgp%dt,fe_p=fe_p))
     n=0; do i=1,n_lay; n=n+count(llgl(i)%do_grain); enddo
     cont_p=(cont_p.and.n.gt.0)

     if (mod(iter_plot,llgp%plot_iter).eq.0) call plot_llg(t,it=it-1); iter_plot=iter_plot+1
     call update_do_grain(t)
   enddo

   call plot_llg(t,it=it,done_p=.true.)
   i=0; do k=1,size(llgl); i=i+sum(max(0,int(sign(1.d0,llgl(k)%m(llgl(k)%in+3,:))))); enddo
   write(os,"(2x,'iter:',i7,' dt=',1p,e13.6,' t=',e13.6,' delta t=',e13.6)") it,llgp%dt,t,del_t
   if (.not.present(track))call output(trim(os))
   if (all(llgl(:)%errl.eqv..false.)) then                          !converged!
      write(os,"('CONVERGED |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   elseif (all(llgl(:)%err1l.eqv..false.)) then
      write(os,"('CONVERGED theta:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         acos(1.d0-.5d0*llgp%errmx1)*180.d0/pi,hscale(1,1),llgp%mave(1:3),i
   elseif (del_t.ge.llgp%del_t) then              !elapsed time experied
      write(os,"('TIME ELAPSED, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   elseif (all(llgl(:)%mave_min)) then
      write(os,"('MIN M, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   elseif (all(llgl(:)%mave_max)) then
      write(os,"('MAX M, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   elseif (all(llgl(:)%mave_max.or.llgl(:)%mave_min)) then
      write(os,"('MAX OR MIN M, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   elseif (n.eq.0) then                                             !t>specified time
      write(os,"('NO MORE SPINS, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   else                                             !t>specified time
      write(os,"('INTEGRATED TIME, |mxH|/Ms:',1p,e13.6,', scale(L1):',e12.5,', <m>=(',2(e12.5,','),e12.5,'), N+=',0p,i7)") &
         llgp%errmx,hscale(1,1),llgp%mave(1:3),i
   endif
   call output(trim(os))

   deallocate(ave_tq)
  contains
   function get_fields(t) result (ok_p)
    use field_m, only: get_field
    use transform_m, only: cross
!   use applied_m, only: rot_fld
    real(8),intent(in)::t
    logical::ok_p(2)
    integer::i,n
    real(8)::hmag,dti,fact1,fact0
    logical,save::do_rot_p=.false. !.true.
    real(8),save::r0(2)
    real(8)::ha(2),hav(3)

    if (it.eq.0.or.do_rot_p) then
!      dti=atan2(llgp%mave(2),llgp%mave(1))+pi/2.d0
       dti=atan2(llgp%mave(2),llgp%mave(1))
       fact0=sqrt(dot_product(llgp%mave,llgp%mave)+1.d-300)
       fact1=acos(llgp%mave(3)/fact0)
!      do_rot_p=rot_fld((/cos(dti),sin(dti),0.d0/))
!       if (mod(it+1,llgp%iter).eq.0) print "(1p,6(e12.5,x))",llgp%mave,dti
       if (do_rot_p) then
         if (it.eq.0) then
               call draw('std(M)',strip=.true.,legend='std(M)',x=(/t/)) !,add=.false.)
!              call draw('<M>',strip=.true.,legend='<M>') !,add=.false.)
               call draw('spherical angles',strip=.true.,legend='<M>1',x=(/t/)) !,add=.false.)
!              call draw('<M>2',strip=.true.,legend='<M2>') !,add=.false.)
               call draw('<M>',strip=.true.,legend='<M>3',x=(/t/)) !,add=.false.)
         endif
         if (mod(iter_plot,llgp%plot_iter).eq.0) call draw('std(M)',strip=.true.,add=.true.,x=(/t,t,t/)*1e9,y=sqrt(llgp%mvar))
!        if (mod(iter_plot,llgp%plot_iter).eq.0) call draw('<M>',strip=.true.,add=.true.,x=(/t,t,t,t/)*1e9,y=(/fact0,cos(dti)*sin(fact1),sin(dti)*sin(fact1),cos(fact1)/))
         if (mod(iter_plot,llgp%plot_iter).eq.0) call draw('spherical angles',strip=.true.,add=.true.,x=(/t,t/)*1e9,y=(/dti,fact1/))
!        if (mod(iter_plot,llgp%plot_iter).eq.0) call draw('<M>2',strip=.true.,add=.true.,x=(/t,t,t,t/)*1e9,y=(/fact0,llgp%mave(1),llgp%mave(2),llgp%mave(3)/))
         if (mod(iter_plot,llgp%plot_iter).eq.0) call draw('<M>',strip=.true.,add=.true.,x=(/t,t,t,t/)*1e9,y=(/fact0,llgp%mave(1:3)/))

         if (it.eq.0) r0=llgp%mave(1:2)
!if (it.eq.0) print *,r0

         ha=(/ atan2(llgp%mave(2)-r0(2),llgp%mave(1)-r0(1)), 0.d0 /)
         hav=cross(llgp%have,llgp%mave); dti=dot_product(hav,hav) *0.d0
         if (dti.ne.0.d0) then
            hav=hav/sqrt(dti)
            ha=(/ atan2(hav(2),hav(1)), hav(3) /)
         endif
       endif
       call get_field(t,llgp%h_var,.true.,ok_p,no_ms_p=llgp%llb.or.llgp%mt,plot_p=(it.eq.0.or.mod(iter_plot,llgp%plot_iter).eq.0),diff_m=diff_m,rot_phi=ha)
    else
       call get_field(t,llgp%h_var,.true.,ok_p,no_ms_p=llgp%llb.or.llgp%mt,plot_p=(it.eq.0.or.mod(iter_plot,llgp%plot_iter).eq.0),diff_m=diff_m)
    endif
    

    dti=1.d0/sqrt(llgp%dt*llgp%max_dmdt); llgp%h_var=llgp%h_var*dti
    num_ave_i=mod(num_ave,llgp%thermal_num_ave)+1; num_ave=num_ave+1
    if (llgp%dt_frac.ne.0.d0) llgp%dt=1.d22

    do i=1,n_lay
      hmag=0.d0
      do n=llgl(i)%min_n,llgl(i)%max_n
       if (llgl(i)%do_grain(n)) then
!     do n=1,size(llgl(i)%grain_vol)
        hmag=max(hmag, dot_product(llgl(i)%h(:,n)+llgl(i)%hthm(:,n)*dti,llgl(i)%h(:,n)+llgl(i)%hthm(:,n)*dti))
       endif
      enddo
      hmag=sqrt(hmag)
      ave_tq(1:3,0,i)=ave_tq(1:3,0,i)-ave_tq(1:3,num_ave_i,i)

      fact1=(1.d0+llgl(i)%alpha*llgl(i)%alpha)/((1.d0+llgl(i)%alpha)*llgl(i)%gamma)*llgp%dt_frac
      fact0 = 1.d0; if (fact1.ne.0.d0) fact0=0.d0
      llgp%dt=min( llgp%dt, llgp%dt_max, fact1/max(1.d0,hmag) + fact0*llgp%dt )
    enddo
   end function

   function bs(x,s) result (b)
   !kittel, ch 11, eq 22
   ! Tc = (p mub)^2 lambda N /3kV, where lambda M is the mean field, p = g sqrt(S(s+1))
   ! lambda = 3 kTc V / N (g mub)^2 S (S+1) = 3 kTc /Ms g (S+1) mub 
    real(8),intent(in)::x,s
    real(8)::b 
    b=x/(2.d0*s); if (b.ne.0.d0) b=( (2.d0*s+1.d0)/tanh((2.d0*s+1.d0)*b) - 1.d0/tanh(b) ) / (2.d0*s)
   end function
   function dbs(x,s) result (b)
   !  d(Bs(x,s))/dx
    real(8),intent(in)::x,s
    real(8)::b,b1
    b=(1.d0+s)/(3.d0*s)
    if (x.ne.0.d0) then
      b=x/(2.d0*s)
      b1=sinh(b); b=sinh(x+b)
      b=(1.d0/(b1*b1)-(1.d0+2.d0*s)*(1.d0+2.d0*s)/(b*b))/(4.d0*s*s)
    endif
   end function

   subroutine push_zhang(k1,iter,errmx_l,errmx1_l,mmag1,dt,dt0)
    use transform_m, only: cross, rotate
    integer,intent(in)::k1,iter
    real(8),intent(in)::dt,dt0
    real(8),intent(inout)::errmx_l,errmx1_l,mmag1
    type(llg_s),pointer::p
    integer::n,id
    real(8)::scale1,h(3),m(3),t(3),y,scale11,scale2,o(3),d(3),mag,htot(3),mhtot,meq,lambdam,etar,sig_p,th0,dmag,hthm(3)
    real(8),parameter::gmub_kb = 1.343427768d-4 !=g mu_b / k_B = 2.d0*9.27400968d-21/1.3806488d-16

          p=>llgl(k1); id=p%in
          p%scale=-0.5d0*p%gamma
          scale1=p%scale * dt
!         do n=1,size(p%grain_vol)
          do n=p%min_n,p%max_n
           m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n)
           if (p%do_grain(n)) then
            lambdam=p%lambdaMs(n)
            if (.not.llgp%meq_bs) lambdam=lambdam*(p%spin+1.d0)/p%spin
            y=gmub_kb * p%spin / p%temp(n)

            if (abs(llgp%type).lt.3) then  !analytic decay of (m-meq)/tau, explicit rotation
              if (iter.eq.1.and.llgp%type.lt.0) then   !doing the mid-point method, so load in m and mag at dt/2
                m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
              endif
              h=p%h(:,n)
              htot=h+lambdam*mag*m; mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              if (llgp%meq_bs) then
                meq=bs(y*mhtot, p%spin)
              else
                meq=lv(y*mhtot)
              endif

              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                if (llgp%meq_bs) then
                  sig_p=meq/(y*mhtot*dbs(y*mhtot,p%spin)+1.d-80)
                else
                  sig_p=meq/(y*mhtot*dlv(y*mhtot)+1.d-80)  !sig_perp^2
                endif
                scale2 = 2.d0*sig_p/y !Garinin's 1/(tau_s*gamma*alpha)
              else
                scale2 = lambdam !Xu&Zhang's 1/(tau_s*gamma*alpha)
              endif

              hthm=p%hthm(:,n)*sqrt(abs(scale2*meq/(mhtot*mag*dt0+1.d-80)))
              scale2=p%alpha*p%alpha_scale(n)*scale2   !1/(gamma*tau_s)

              if (iter.eq.1.and.llgp%type.lt.0) then; m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n); endif

              t=cross(m,h+hthm); etar=dot_product(t,t)
              if (etar.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=etar/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif

              if (abs(llgp%type).lt.2) then  !analytic rotation and decay of |m|
               !analytic |m|
                dmag = (meq*dot_product(m,htot)-mag)*(1.d0-exp(scale1*scale2)) !half decay
                etar=sqrt(dot_product(h,h)); scale11=2.d0*scale1*etar   !scale11 = - gammaa dt |Heff|
                o=h/(etar+1.d-80); th0=acos(min(1.d0,max(-1.d0,dot_product(o,m*sign(1.d0,mag+dmag)))))
               !deterministic rotation for damping, deterministic+stochastic precession rotation
                d=rotate(th0-2.d0*atan(exp(scale11*meq*scale2/(mag*mhtot+1.d-66))*tan(0.5d0*th0)),cross(m*sign(1.d0,mag+dmag),o), &
                         abs(2.d0*scale1*sqrt(dot_product(h+hthm,h+hthm))),h+hthm,(mag+dmag)*m,.true.,.true.)-mag*m
              else   !explicit gyro-rotation and analytic decay
                d = (meq*htot - mag*m)*(1.d0-exp(scale1*scale2)) !half decay
                d = d+2.d0*scale1*cross(mag*m+d,h+hthm) !full precession rotation (deterministic and stochastic)
!               d = rotate(abs(2.d0*scale1*sqrt(dot_product(h+hthm,h+hthm))),h+hthm,mag*m+d,.true.)-mag*m
              endif

            !now get new meq after 1/2 decay and full rotation
              o=mag*m+d
              htot=p%h(:,n)+lambdam*o; mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              if (llgp%meq_bs) then
                meq=bs(y*mhtot, p%spin)
              else
                meq=lv(y*mhtot)
              endif
              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                if (llgp%meq_bs) then
                  sig_p=meq/(y*mhtot*dbs(y*mhtot,p%spin)+1.d-80)
                else
                  sig_p=meq/(y*mhtot*dlv(y*mhtot)+1.d-80)  !sig_perp^2
                endif
                scale2 = 2.d0*sig_p/y !Garinin's 1/(tau_s*gamma)
              else
                scale2 = lambdam
              endif
              scale2=p%alpha*p%alpha_scale(n)*scale2   !1/(gamma*tau_s)
              if (abs(llgp%type).lt.2) then  !analytic rotation and decay of |m|
                y=sqrt(dot_product(o,o)); if (y.gt.0.d0) o=o/y
                d = d +(meq*dot_product(o,htot)-y)*(1.d0-exp(scale1*scale2))*o !+ eta !half decay + stochastic
              else
                d = d+ (meq*htot - o)*(1.d0-exp(scale1*scale2)) !+eta  !half decay
              endif

              errmx1_l = max(errmx1_l, dot_product(d,d))

              o=m*mag+d; etar=sqrt(dot_product(o,o))
              if (etar.gt.0.d0) then
                p%m(4-id:6-id,n)=o/etar
              else
                etar=0.d0; p%m(4-id:6-id,n)=-m
              endif
              p%mag(2-id/3,n)=etar
!           if (abs(llgp%type).lt.3) then  !semi-implicit push for dm/dt = - g mx(H+Hth) - (m-meq)/tau
!
!               if (iter.eq.1.and.llgp%type.lt.0) then   !doing the mid-point method, so load in m and mag at dt/2
!                 m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
!               endif
!               h=p%h(:,n)
!               htot=h+lambdam*mag*m
!               mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
!               if (llgp%meq_bs) then
!                 meq=bs( gmub_kb * p%spin *mhtot / p%temp(n), p%spin)
!               else
!                 meq = lv(gmub_kb * p%spin *mhtot / p%temp(n))
!               endif
!               scale11 = 2.d0 * scale1
!               scale2 = p%alpha*lambdam
!
!               y=abs(scale2*meq/(mhtot*mag+1.d-66))
!               h=h+p%hthm(:,n)*sqrt(abs(y/(p%alpha*p%alpha_scale(n)*dt0))) !*0.d0
!
!               m=p%m(1+id:3+id,n); mag=p%mag(1+id/3,n)
!
!               t=cross(m,h); y=dot_product(t,t)
!               if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
!
!     ! trying to solve dm/dt = - g m x H - (m-meq)/tau
!     ! discretize as dm = - g dt (m + dm/2) x H - (m + dm/2 - meq)/tau, where dm = m(new)-m(old) = m(new) - m and H is Heff + Hthermal
!     ! so (1+dt/2tau) dm = B - (g dt/2) dm x H   (i.e. B = constant = - g dt m x H - (m - meq)dt/tau), 1/tau = alpha gamma lambdam 
!     !                   = B - (g dt/2) / (1+dt/2tau) [ B     - (g dt/2) dm x H] x H
!     !                   = B - (g dt/2) / (1+dt/2tau) [ B x H - (g dt/2) ( - |H|^2 dm + dm.H H )]
!     ! (1+dt/2tau + (g dt |H|/2)^2/(1+dt/2tau)) dm = B - (g dt/2) / (1+dt/2tau) B x H + [(g dt/2) / (1+dt/2tau)]^2 B.H
!
!               lambdam = 1.d0 - scale2 * scale1  ! = 1+ alpha lambdaMs gamma dt / 2 = 1 + dt / 2 tau
!               y = scale1 / lambdam  ! = - gamma dt / (2(1 + dt/2tau))
!
!               t = scale11 * ( mag*t + scale2 * (mag*m - meq*htot) )  ! - gamma dt ( m x H + a lambdaMs (m-meq))
!               t = t + y * cross(t,h) + y * y * dot_product(t,h) * h
!               d = t / (lambdam * (1.d0 + y * y * dot_product(h,h)))
!
!               errmx1_l = max(errmx1_l, dot_product(d,d))
!
!               o=m*mag+d; y=sqrt(dot_product(o,o))
!               if (y.gt.0.d0) then
!                 p%m(4-id:6-id,n)=o/y
!               else
!                 y=0.d0; p%m(4-id:6-id,n)=-m
!               endif
!               p%mag(2-id/3,n)=y
            elseif (abs(llgp%type).eq.4) then  !rotational, i.e. boris push- for historical reasons... shouldn't be used anymore!
              if (mag.gt.0.d0) then
                scale11 = scale1
                if (iter.eq.2.or.llgp%type.gt.0) then !first or only step
                  htot=p%h(:,n)+lambdam*mag*m
                  mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                  if (llgp%meq_bs) then
                    mag=bs( gmub_kb * p%spin *mhtot / p%temp(n), p%spin)*dot_product(m,htot)
                  else
                    mag = lv(gmub_kb * p%spin *mhtot / p%temp(n))*dot_product(m,htot)
                  endif
                  mag = mag + (p%mag(id/3+1,n)-mag)*exp(scale11*p%alpha*lambdam)
                else
                  m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
                endif
! if (n.eq.1) print *,'<',p%lambdaMs(n),p%mag(id/3+1,n),htot,mhtot,mag,gmub_kb * p%spin,dot_product(m,htot),bs(gmub_kb * p%spin *mhtot / p%temp(n), p%spin)

                htot=p%h(:,n)+lambdam*mag*m
                mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                if (llgp%meq_bs) then
                  meq=bs( gmub_kb * p%spin *mhtot / p%temp(n), p%spin)
                else
                  meq= lv(gmub_kb * p%spin *mhtot / p%temp(n))
                endif
                scale2 = -p%alpha*lambdam*meq/(abs(mhtot*mag)+1.d-66)

                m=p%m(1+id:3+id,n)
                h=p%h(:,n)+p%hthm(:,n)*sqrt(abs(scale2/(p%alpha*p%alpha_scale(n)*dt0)))

                t=cross(m,h)
                y=dot_product(t,t);
                if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif

!                if (abs(llgp%type).ne.2) then
!!                 o=(1.d0+scale2*scale2)*scale11*h+scale2*m !damping goes as -gamma alpha m x m x h
!                  o=scale11*(1.d0+scale2*scale2)*h+scale2*m !damping goes as alpha m x dm/dt, but still LL form, dm/dt=-g(1-a^2) mxH + a mxdm/dt
!                else
                  o=scale11*h+scale2*m                      !damping goes as alpha m x dm/dt
!                endif

                t=m+cross(m,o); y=2.d0/(1.d0+dot_product(o,o))
                d=y*cross(t,o)
!if (n.eq.1) print *,'>',scale2,scale11,m,h,cross(m,h),y,o,t,d
                p%m(4-id:6-id,n)=m+d

                htot=p%h(:,n)+lambdam*p%mag(id/3+1,n)*p%m(1+id:3+id,n)
                mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                if (llgp%meq_bs) then
                  mag=bs( gmub_kb * p%spin *mhtot / p%temp(n), p%spin)*dot_product(m,htot)
                else
                  mag = lv(gmub_kb * p%spin *mhtot / p%temp(n))*dot_product(m,htot)
                endif
                mag = mag + (p%mag(id/3+1,n)-mag)*exp(2.d0*scale11*p%alpha*lambdam)
                p%mag(2-id/3,n)=abs(mag)
                if (mag.lt.0.d0) then
                  o=-(m+0.5d0*d); p%m(4-id:6-id,n)=o/sqrt(dot_product(o,o))
                  d=p%m(4-id:6-id,n)-m
                endif

                y=0.5d0*(p%mag(2-id/3,n)+p%mag(id/3+1,n))
                errmx1_l=max(errmx1_l,dot_product(d*y+m*(p%mag(2-id/3,n)-p%mag(id/3+1,n)),d*y+m*(p%mag(2-id/3,n)-p%mag(id/3+1,n))))
              else   !ms_scale = 0
                h=p%h(:,n); y = sqrt(dot_product(h,h))
                if (y.gt.0.d0) then
                  p%m(4-id:6-id,n)=h/y
                  if (llgp%meq_bs) then
                    mag=bs(gmub_kb*p%spin*y/p%temp(n),p%spin)
                  else
                    mag=lv(gmub_kb*p%spin*y/p%temp(n))
                  endif
                  p%mag(2-id/3,n)=mag*(1.d0-exp(2.d0*scale1*p%alpha*lambdam))
                else
                  p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
                  p%mag(2-id/3,n)=p%mag(1+id/3,n)
                endif
                errmx1_l=max(errmx1_l,(p%mag(2-id/3,n)-p%mag(1+id/3,n))**2)
              endif
            else   !heun
              scale11 = 2.d0 * scale1
              if (iter.eq.1) then
                m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
                if (mag.le.0.d0) then; mag=0.d0; y=sqrt(dot_product(p%h(:,n),p%h(:,n))); if (y.ne.0.d0) m=p%h(:,n)/y; endif
              endif

              htot=p%h(:,n)+lambdam*mag*m
              mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              y=y*mhtot
              if (llgp%meq_bs) then
                meq = bs(y, p%spin)
              else
                meq = lv(y)
              endif
              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                if (llgp%meq_bs) then
                  y=meq/(y*dbs(y,p%spin)+1.d-66)   !(sig_perp/sig_parallel)^2
                else
                  y=meq/(y*dlv(y)+1.d-66)   !(sig_perp/sig_parallel)^2
                endif
                scale2 = 2.d0*p%temp(n)*y/(gmub_kb*p%spin) !Garinin's 1/(tau_s*gamma*alpha)
              else
                scale2 = lambdam !Xu&Zhang's 1/(tau_s*gamma*alpha)
              endif

   ! p%hthm(1:3,n) = sqrt( 2 k T alpha / gamma Ms(0) V |m| )
              y=abs(scale2*meq/(mhtot*mag+1.d-66))
              h=p%h(:,n)+p%hthm(:,n)*sqrt(abs(y/dt0))
              scale2=scale2*p%alpha*p%alpha_scale(n)  !1/(gamma*tau_s)

              if (llgp%type.eq.-3) then
                d=scale11*(cross(m*mag,h)+(m*mag-meq*htot)*scale2)
              else
                d=cross(m,h)
                d=scale11*( (mag*d+scale2*meq/(mhtot+1.d-22)*cross(m,d))/(1.d0+y*y)+scale2*(mag-meq*dot_product(m,htot))*m)
              endif

              if (iter.eq.1) then
                d=0.5d0*(d+p%dm(:,n))
              else
                p%dm(:,n)=d
              endif

              t=cross(m,h); y=dot_product(t,t)
              if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
              errmx1_l = max(errmx1_l, dot_product(d,d))

              m=p%m(1+id:3+id,n); mag=p%mag(1+id/3,n)
              o=m*mag+d; y=sqrt(dot_product(o,o))
              if (y.gt.0.d0) then
                p%m(4-id:6-id,n)=o/y
              else
                y=0.d0; p%m(4-id:6-id,n)=-m
              endif
              p%mag(2-id/3,n)=y

            endif
           else
            p%m(4-id:6-id,n)=m
            p%mag(2-id/3,n)=mag
            p%dm(:,n)=0.d0
           endif
           mmag1 = mmag1 +p%mag(2-id/3,n)*p%grain_vol(n)
          enddo
!    if (k1.eq.1) print *,'!',id,p%m(:,1),p%h(:,1),p%hthm(:,1),p%m(4-id:6-id,1)-p%m(1+id:3+id,1),p%ms_scale(1),p%temp(1)
!    read *
   end subroutine

   function lv(x) result (b)
  !lv(x) = coth(x)-1/x
    real(8),intent(in)::x
    real(8)::b 
    
    b=x/3.d0; if (abs(x).gt.1.d-2) b=1.d0/tanh(x)-1.d0/x
   end function 

   function find_meq(heff,hx) result (b)
    real(8),intent(in)::heff,hx
    real(8)::b0,b1,b

    b0=0.d0; b1=1.d0
    do while (b1-b0.gt.1.d-9)
      b=(b1+b0)*0.5d0
      if (b.gt.lv(heff+hx*b)) then
        b1=b
      else
        b0=b
      endif
    enddo
   end function
   
   function dlv(x) result (b)
  !ld(lv(x))/dx = 1/x^2 - Csch(x)^2
    real(8),intent(in)::x
    real(8)::b 
    
    b=1.d0/3.d0; if (abs(x).gt.4.d-3) then; b=sinh(x); b=1.d0/(x*x)-1.d0/(b*b); endif
   end function

   subroutine push_mt(k1,iter,errmx_l,errmx1_l,mmag1,dt,dt0)
   ! do dm/dt = - gamma m x Heff - (m-meq)/tau + eta
    use transform_m, only: cross, rotate
    integer,intent(in)::k1,iter
    real(8),intent(in)::dt,dt0
    real(8),intent(inout)::errmx_l,errmx1_l,mmag1
    type(llg_s),pointer::p
    integer::n,id
    real(8)::scale1,h(3),m(3),t(3),y,scale11,scale2,o(3),d(3),mag,htot(3),mhtot,meq,lambdam,eta(3),etar,dmag
    real(8)::sig_l,sig_p,th0
    real(8),parameter::gmub_kb = 1.343427768d-4 !=g mu_b / k_B = 2.d0*9.27400968d-21/1.3806488d-16

          p=>llgl(k1); id=p%in
          p%scale=-0.5d0*p%gamma
          scale1=p%scale * dt
          scale11 = 2.d0 * scale1
!         do n=1,size(p%grain_vol)
          do n=p%min_n,p%max_n
           m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n)
           if (p%do_grain(n)) then

            lambdam=p%lambdaMs(n)
            if (.not.llgp%meq_bs) lambdam=lambdam*(p%spin+1.d0)/p%spin
!   print *,1.d0/(p%lambdaMs(n)*p%gamma*p%alpha),p%spin/((p%spin+1.d0)*p%lambdaMs(n)*p%gamma*p%alpha); read *

            y=gmub_kb * p%spin / p%temp(n)

            if (abs(llgp%type).lt.3) then  !analytic decay of (m-meq)/tau, explicit rotation
              if (iter.eq.1.and.llgp%type.lt.0) then   !doing the mid-point method, so load in m and mag at dt/2
                m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
              endif
              h=p%h(:,n)
              htot=h+lambdam*mag*m; mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              if (llgp%meq_bs) then
                meq=bs(y*mhtot, p%spin)
              else
                meq=lv(y*mhtot)
              endif

              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                if (llgp%meq_bs) then
                  sig_p=meq/(y*mhtot*dbs(y*mhtot,p%spin)+1.d-80)  !sig_perp^2/sig_parallel^2
                else
                  sig_p=meq/(y*mhtot*dlv(y*mhtot)+1.d-80)  !sig_perp^2/sig_parallel^2
                endif
                scale2 = 2.d0*p%alpha_scale(n)*p%alpha*sig_p/y !Garinin's 1/(tau_s*gamma)
              else
                scale2 = p%alpha_scale(n)*p%alpha*lambdam !Xu&Zhang's 1/(tau_s*gamma)
              endif

              if (iter.eq.1.and.llgp%type.lt.0) then; m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n); endif

              t=cross(m,h); etar=dot_product(t,t)
              if (etar.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=etar/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif

              if (abs(llgp%type).lt.2) then  !analytic rotation and decay of |m|
                dmag = (meq*dot_product(m,htot)-mag)*(1.d0-exp(scale1*scale2)) !half decay
                etar=dot_product(h,h)
                if (etar.gt.0.d0) then
                  etar=sqrt(etar); scale11=2.d0*scale1*etar   !scale11 = - gammaa dt |Heff|
                  o=h/etar; th0=acos(min(1.d0,max(-1.d0,dot_product(o,m*sign(1.d0,mag+dmag)))))
                  d=rotate(th0-2.d0*atan(exp(scale11*meq*scale2/(mag*mhtot+1.d-66))*tan(0.5d0*th0)), &
                         cross(m*sign(1.d0,mag+dmag),o),abs(scale11),o,(mag+dmag)*m,.true.)-mag*m
                else
                  d=dmag*m
                endif
              else   !explicit gyro-rotation and analytic decay
                d = (meq*htot - mag*m)*(1.d0-exp(scale1*scale2)) !half decay
                d = d+scale11*cross(mag*m+d,h) !full rotation
              endif

            !now get new meq after 1/2 decay and full rotation
              o=mag*m+d
              htot=h+lambdam*o; mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              if (llgp%meq_bs) then
                meq=bs(y*mhtot, p%spin)
                sig_p=meq/(y*mhtot+1.d-80)  !sig_perp^2
                sig_l=dbs(y*mhtot,p%spin)   !sig_parallel^2
              else
                meq=lv(y*mhtot)
                sig_p=meq/(y*mhtot+1.d-80)  !sig_perp^2
                sig_l=dlv(y*mhtot)   !sig_parallel^2
              endif
              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                eta=dt*p%hthm(:,n)*p%gamma*sig_p*sqrt(2.d0/((sig_l+1.d-80)*dt0))    !eta, using Garanin's tau_s
!               scale2 = 2.d0*p%alpha_scale(n)*p%alpha*p%temp(n)*sig_p/(gmub_kb*p%spin*(sig_l+1.e-22)) !Garinin's 1/(tau_s*gamma)
                scale2 = 2.d0*p%alpha_scale(n)*p%alpha*sig_p/(y*(sig_l+1.d-80)) !Garinin's 1/(tau_s*gamma)
              else
                eta=dt*p%hthm(:,n)*p%gamma*sqrt(lambdam*meq/((mhtot+1.d-80)*dt0))   !eta
              endif
              eta=sqrt(sig_l/sig_p)*dot_product(htot,eta)*htot-cross(htot,cross(htot,eta))
              if (abs(llgp%type).lt.2) then  !analytic rotation and decay of |m|
                y=sqrt(dot_product(o,o)); if (y.gt.0.d0) o=o/y
!               d = d +(meq*dot_product(o,htot)-y)*(1.d0-exp(scale1*scale2))*y*o + eta !half decay + stochastic
                d = d +(meq*dot_product(o,htot)-y)*(1.d0-exp(scale1*scale2))*o + eta !half decay + stochastic
              else
                d = d+ (meq*htot - o)*(1.d0-exp(scale1*scale2))+eta  !half decay
              endif

              errmx1_l = max(errmx1_l, dot_product(d,d))

              o=m*mag+d; etar=sqrt(dot_product(o,o))
              if (etar.gt.0.d0) then
                p%m(4-id:6-id,n)=o/etar
              else
                etar=0.d0; p%m(4-id:6-id,n)=-m
              endif
              p%mag(2-id/3,n)=etar
            else   !heun
              if (iter.eq.1) then
                m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
                if (mag.le.0.d0) then; mag=0.d0; etar=sqrt(dot_product(p%h(:,n),p%h(:,n))); if (etar.ne.0.d0) m=p%h(:,n)/etar; endif
              endif

              htot=p%h(:,n)+lambdam*mag*m
              mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
              y=y*mhtot
              if (llgp%meq_bs) then
                meq=bs( y, p%spin)
                sig_l=dbs(y,p%spin)   !sig_parallel^2
                sig_p=meq/(y+1.d-22)  !sig_perp^2
              else
                meq = lv(y)
                sig_l=dlv(y)   !sig_parallel^2
                sig_p=meq/(y+1.d-22)  !sig_perp^2
              endif

   !p%scale = -0.5 gamma
   ! scale1 = -0.5 gamma dt
   ! scale2 = alpha * lambdaMs
   ! scale11= - gamma dt
   ! htot = \hat{Htot}
   ! mhtot = Htot
   ! p%hthm(1:3,n) = sqrt( 2 k T alpha / gamma Ms(0) V )
   ! p%hthm(1:3,n) = sqrt( 2 k T alpha / gamma Ms(0) |m| V )
   ! gmub_kb = 1.343427768d-4 !=g mu_b / k_B = 2.d0*9.27400968d-21/1.3806488d-16

              if (llgp%garanin_tau) then !use Garinin's tau_s (based on Gamma_1)
                eta=dt*p%hthm(:,n)*p%gamma*sig_p*sqrt(2.d0/((sig_l+1.d-22)*dt0))    !eta, using Garanin's tau_s
                eta=sqrt(sig_l/sig_p)*dot_product(htot,eta)*htot-cross(htot,cross(htot,eta))
                scale2 = 2.d0*p%temp(n)*sig_p/(gmub_kb*p%spin*(sig_l+1.e-22)) !Garinin's 1/(tau_s*gamma*alpha)
              else
                eta=dt*p%hthm(:,n)*p%gamma*sqrt(lambdam*meq/((mhtot+1.d-22)*dt0))   !eta
!               eta=etar*dot_product(htot,eta)*htot-cross(htot,cross(htot,eta))
                eta=sqrt(sig_l/sig_p)*dot_product(htot,eta)*htot-cross(htot,cross(htot,eta))
                scale2 = lambdam !Xu&Zhang's 1/(tau_s*gamma*alpha)
              endif

              h=p%h(:,n)
              d=scale11*(cross(m*mag,h)+(m*mag-meq*htot)*scale2*p%alpha*p%alpha_scale(n)) + eta

              if (iter.eq.1) then
                d=0.5d0*(d+p%dm(:,n))
              else
                p%dm(:,n)=d
              endif

              t=cross(m,h); etar=dot_product(t,t)
              if (etar.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=etar/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
              errmx1_l = max(errmx1_l, dot_product(d,d))

              m=p%m(1+id:3+id,n); mag=p%mag(1+id/3,n)
              o=m*mag+d; etar=sqrt(dot_product(o,o))
              if (etar.gt.0.d0) then
                p%m(4-id:6-id,n)=o/etar
              else
                etar=0.d0; p%m(4-id:6-id,n)=-m
              endif
              p%mag(2-id/3,n)=etar
            endif
           else
            p%m(4-id:6-id,n)=m
            p%mag(2-id/3,n)=mag
            p%dm(:,n)=0.d0
           endif
           mmag1 = mmag1 +p%mag(2-id/3,n)*p%grain_vol(n)
          enddo
!    if (k1.eq.1) print *,'!',id,p%m(:,1),p%h(:,1),p%hthm(:,1),p%m(4-id:6-id,1)-p%m(1+id:3+id,1),p%ms_scale(1),p%temp(1)
!    read *
   end subroutine

   subroutine push_llb(k1,iter,errmx_l,errmx1_l,mmag1,dt,dt0)
    use transform_m, only: cross
    integer,intent(in)::k1,iter
    real(8),intent(in)::dt,dt0
    real(8),intent(inout)::errmx_l,errmx1_l,mmag1
    type(llg_s),pointer::p
    integer::n,id
    real(8)::scale1,h(3),m(3),t(3),y,scale11,scale2,o(3),d(3),mag,htot(3),mhtot,meq,aperp,apar,hthmpar(3),hthmperp(3),lambdam
    real(8),parameter::gmub_kb = 1.343427768d-4 !=g mu_b / k_B = 2.d0*9.27400968d-21/1.3806488d-16

          p=>llgl(k1); id=p%in
          p%scale=-0.5d0*p%gamma
          scale1=p%scale * dt
!         do n=1,size(p%grain_vol)
          do n=p%min_n,p%max_n
           if (p%do_grain(n)) then
            lambdam=p%lambdaMs(n)*(p%spin+1)/p%spin
            if (p%mag(id/3+1,n).gt.0.d0) then
              if (abs(llgp%type).lt.3) then  !rotational, i.e. boris push

                scale11=scale1 

                if (iter.eq.2.or.llgp%type.gt.0) then !first or only step
                  m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n)
                  htot=p%h(:,n)+lambdam*mag*m
                  mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                  meq = gmub_kb * p%spin * mhtot / p%temp(n)
                  mag = lv(meq)*dot_product(m,htot)
                  mag=mag+(p%mag(id/3+1,n)-mag) &
                     * exp(scale1*2.d0*p%alpha*lambdam*lambdam*p%temp(n)*p%temp(n)*lv(meq) &
                     / (3.d0*p%tc(n)*p%tc(n)*dlv(meq)*(mhtot+1.d-20)))
                else
                  m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
                endif

                htot=p%h(:,n)+lambdam*mag*m
                mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                meq = gmub_kb * p%spin * mhtot / p%temp(n)

                aperp=p%alpha*(1.d0-lambdam*lv(meq)*p%temp(n)/(3.d0*mhtot*p%tc(n)))
                apar=p%alpha*2.d0*p%temp(n)/(3.d0*p%tc(n))

                scale2 = - aperp / (abs(mag)+1.d-20)

                h=p%h(:,n)+p%hthm(:,n)*sqrt(abs(scale2/(p%alpha*p%alpha_scale(n)*dt0*(abs(mag)+1.d-20))))
                m=p%m(1+id:3+id,n)

                t=cross(m,h); y=dot_product(t,t)
                if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif

                if (abs(llgp%type).ne.2) then
                  o=(1.d0+scale2*scale2)*scale11*h+scale2*m !damping goes as -gamma alpha m x m x h
                else
                  o=scale11*h+scale2*m                      !damping goes as alpha m x dm/dt
                endif

                t=m+cross(m,o); y=2.d0/(1.d0+dot_product(o,o))
                d=y*cross(t,o)
                p%m(4-id:6-id,n)=m+d

                htot=p%h(:,n)+lambdam*mag*(m+0.5d0*d)
                mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                meq = gmub_kb * p%spin * mhtot /p%temp(n)
                mag = lv(meq)*dot_product(sign(1.d0,mag)*(m+0.5d0*d),htot)
                meq=mag+(p%mag(id/3+1,n)-mag) &
                     * exp(scale1*4.d0*p%alpha*lambdam*lambdam*p%temp(n)*p%temp(n)*lv(meq) &
                     / (3.d0*p%tc(n)*p%tc(n)*dlv(meq)*(mhtot+1.d-20)))
                meq=meq-2.d0*scale11*dot_product(p%hthm(:,n)*sqrt(abs(apar/(p%alpha*p%alpha_scale(n)*dt0))),m+0.5d0*d)
                p%mag(2-id/3,n)=abs(meq)
                if (meq.lt.0.d0) then
                  o=-(m+0.5d0*d); p%m(4-id:6-id,n)=o/sqrt(dot_product(o,o))
                  d=p%m(4-id:6-id,n)-m
                endif
                y=0.5d0*(p%mag(2-id/3,n)+p%mag(id/3+1,n))
                errmx1_l=max(errmx1_l,dot_product(d*y+m*(p%mag(2-id/3,n)-p%mag(id/3+1,n)),d*y+m*(p%mag(2-id/3,n)-p%mag(id/3+1,n))))

              else   !heun
                if (iter.eq.1) then
                  m=p%m(4-id:6-id,n); mag=p%mag(2-id/3,n)
                  if (mag.le.0.d0) then; mag=0.d0; y=sqrt(dot_product(p%h(:,n),p%h(:,n))); if (y.ne.0.d0) m=p%h(:,n)/y; endif
                else
                  m=p%m(1+id:3+id,n); mag=p%mag(id/3+1,n)
                endif

                h=p%h(:,n)
                htot=h+lambdam*mag*m
                mhtot=sqrt(dot_product(htot,htot)); if (mhtot.ne.0.d0) htot=htot/mhtot
                meq = gmub_kb * p%spin * mhtot / p%temp(n)
                aperp = p%alpha*(1.d0-lambdam*lv(meq)*p%temp(n)/(3.d0*mhtot*p%tc(n)))
                apar  = p%alpha* 2.d0 * p%temp(n)/(3.d0*p%tc(n))

                h=p%hthm(:,n)/sqrt(p%alpha*p%alpha_scale(n)*dt0)
                hthmpar  = dot_product(h,m)*m; hthmperp = h-hthmpar

!      if (p%temp(n).ge.p%tc(n)) then
!        apar  = p%alpha* 2.d0/3.d0
!        aperp=apar
!      endif
                hthmpar  = hthmpar *sqrt(apar) 
                hthmperp = hthmperp*sqrt(max(0.d0,aperp-apar)/(aperp*aperp))

                scale11 = 2.d0*scale1
                if (llgp%type.ne.-3) scale11=scale11/(1.d0+aperp*aperp/(mag*mag+1.d-66))
                meq = apar*(lambdam/mhtot)*lambdam*2.d0*scale1*p%temp(n)*lv(meq)/(3.d0*dlv(meq)*p%tc(n))*(mag-lv(meq)*dot_product(m,htot))
                t = cross(m,p%h(:,n)+hthmperp)
                d = scale11*(mag*cross(m,p%h(:,n))+aperp*cross(m,t)+hthmpar)+sign(min(0.5d0,abs(meq)),meq)*m

                if (iter.eq.1) then
                  d=0.5d0*(d+p%dm(:,n))
                else
                  p%dm(:,n)=d
                endif

                y=dot_product(t,t); if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
                errmx1_l = max(errmx1_l, dot_product(d,d))

                o=p%m(1+id:3+id,n)*p%mag(id/3+1,n)+d !m(t+dt)=m(t)+dm
                y=sqrt(dot_product(o,o))
                if (y.gt.0.d0) then
                  p%m(4-id:6-id,n)=o/y
                else
                  y=0.d0
                  p%m(4-id:6-id,n)=-p%m(1+id:3+id,n)
                endif
                p%mag(2-id/3,n)=y
              endif
            else   !ms_scale = 0
              h=p%h(:,n)
              y = sqrt(dot_product(h,h))
              if (y.gt.0.d0) then
                p%m(4-id:6-id,n)=h/y
                scale11=2.d0*scale1
                meq=gmub_kb*p%spin*y/p%temp(n)
                scale2=lambdam*p%temp(n)/(3.d0*dlv(meq)*p%tc(n))
                p%mag(2-id/3,n)=(y/scale2+lv(meq))*(1.d0-exp(scale11*scale2*2.d0*p%temp(n)*p%alpha/(3.d0*p%tc(n))))
              else
                p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
                p%mag(2-id/3,n)=p%mag(id/3+1,n)
              endif
              errmx1_l=max(errmx1_l,(p%mag(2-id/3,n)-p%mag(id/3+1,n))**2)
            endif
           else
            p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
            p%mag(2-id/3,n)=p%mag(id/3+1,n)
            p%dm(:,n)=0.d0
           endif
           mmag1 = mmag1 +p%mag(2-id/3,n)*p%grain_vol(n)
          enddo

   end subroutine
   subroutine push_llg(k1,iter,errmx_l,errmx1_l,mmag1,dt,dt0)
    use transform_m, only: cross,rotate
    integer,intent(in)::k1,iter
    real(8),intent(in)::dt,dt0
    real(8),intent(inout)::errmx_l,errmx1_l,mmag1
    type(llg_s),pointer::p
    integer::n,id
    real(8)::fact1,fact0,scale1,h(3),m(3),t(3),y,scale11,scale2,o(3),d(3),th0

          p=>llgl(k1); id=p%in
          fact1=(1.d0+p%alpha*p%alpha)/((1.d0+p%alpha)*p%gamma)*llgp%dt_frac
          fact0 = 1.d0; if (fact1.ne.0.d0) fact0=0.d0
          p%scale=-0.5d0*p%gamma
          scale1=p%scale * dt
!         do n=1,size(p%grain_vol)
          do n=p%min_n,p%max_n
           if (p%do_grain(n)) then
            h=p%hthm(:,n)/sqrt(dt0) + p%h(:,n)
            if (p%ms_scale(n).gt.0.d0) then
              m=p%m(1+id:3+id,n)
              if (abs(llgp%type).lt.3) then  !rotational, i.e. boris push

                t=cross(m,h)
                y=dot_product(t,t)
                if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif

                scale2=p%alpha_scale(n)*p%alpha
                if (abs(llgp%type).ne.2) then
         !do analytic rotations
         ! precession term is dphi = g dt H
         ! damping term is dthet = 2 ArcTan[ exp(-g alpha dt H) tan(theta0/2) ] - theta0
         ! where theta0 = ArcCos(\hat{m}.\hat{h})
!     if (mod(it,50).eq.0.and.n.eq.1) print "(1p,:,10(e12.5,x))",d,sqrt(dot_product(d,d))
                  y=dot_product(h,h)
                  if (y.gt.0.d0) then
                    y=sqrt(y); scale11=2.d0*scale1*y/(1.d0+scale2*scale2)
                    o=h/y; th0=acos(min(1.d0,max(-1.d0,dot_product(o,m))))
                    d=rotate(th0-2.d0*atan(exp(scale11*scale2)*tan(0.5d0*th0)),cross(m,o),abs(scale11),o,m,.true.)-m
!if (mod(it,50).eq.0.and.n.eq.1) print "(1p,:,50(e12.5,x))",d,sqrt(dot_product(d,d)),exp(scale11*scale2),y,scale2,th0,2.d0*atan(exp(scale11*scale2)*tan(0.5d0*th0)),acos(dot_product(o,m+d)),dot_product(m+d,m+d),dot_product(t,t)
                  else
                    d=0.d0
                  endif
                else
         ! for |type| = 2
         !dm/dt = - g m x H + a m x dm/dt, do semi-implicit differencing:
         !-> dm = - g dt (m + dm/2) x H + a (m + dm/2) x dm, where m is the initial m and dm=m(final)-m
         !->    = (2m + dm) x (-g H dt/2 - a dm) = (2m + dm) x O    notice a (2m + dm) x dm / 2 = (2m+dm) x (-a m)
         !->    = 2m x O + dm x O = 2 m x O + ((2m + dm) x O) x O = 2m x O + 2(m x O) x O + O.dm O - |O|^2 dm, but O.dm=0 (i.e. dm = (2m+dm) x O)
         !=> dm (1 + O^2) =2(m + m x O) x O
                  o=scale1*h-scale2*m                      !damping goes as alpha m x dm/dt
                  t=m+cross(m,o)
                  y=2.d0/(1.d0+dot_product(o,o))
                  d=y*cross(t,o)
                endif
                errmx1_l=max(errmx1_l,dot_product(d,d))
                p%m(4-id:6-id,n)=m+d
              else   !heun
                scale2 = p%alpha_scale(n)*p%alpha
                scale11 = 2.d0 * scale1/(1.d0+scale2*scale2)
                if (iter.eq.1) m=p%m(4-id:6-id,n)  !second step
                t = cross(m,h); d=scale11*(t+scale2*cross(m,t))
                if (iter.eq.1) then
                  d=0.5d0*(d+p%dm(:,n))
                else
                  p%dm(:,n)=d
                endif
                y=dot_product(t,t)
                if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
                errmx1_l = max(errmx1_l, dot_product(d,d))
                p%m(4-id:6-id,n)=p%m(1+id:3+id,n)+d

                y=dot_product(p%m(4-id:6-id,n),p%m(4-id:6-id,n))
                if (y.gt.0.d0) p%m(4-id:6-id,n)=p%m(4-id:6-id,n)/sqrt(y)
              endif
            else   !ms_scale = 0
              h=p%h(:,n)
              y = sqrt(dot_product(h,h))
              if (y.gt.0.d0) then
                p%m(4-id:6-id,n)=h/y
              else
                p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
              endif
            endif
           else
            p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
            p%dm(:,n)=0.d0
           endif
           mmag1 = mmag1 + sqrt(dot_product(p%m(4-id:6-id,n),p%m(4-id:6-id,n)))*p%ms_scale(n)*p%grain_vol(n)
          enddo
   end subroutine

   subroutine push_fp(k1,errmx_l,errmx1_l,mmag1,dt,dt0)
    use FP_m, only: advance_FP !,FP_s,new_FP
    use transform_m, only: cross
    integer,intent(in)::k1
    real(8),intent(in)::dt,dt0
    real(8),intent(inout)::errmx_l,errmx1_l,mmag1
    type(llg_s),pointer::p
    integer::n,id
    real(8)::h(3),m(3),t(3),y,d(3)
    real(8),parameter::gmub_kb = 1.343427768d-4 !=g mu_b / k_B = 2.d0*9.27400968d-21/1.3806488d-16

          p=>llgl(k1); id=p%in
!         do n=1,size(p%grain_vol)
          do n=p%min_n,p%max_n
           if (p%do_grain(n)) then
            if (abs(llgp%type).lt.3) then  !rotational, i.e. boris push
              m=p%m(1+id:3+id,n)
              h=p%hthm(:,n)/sqrt(dt0) + p%h(:,n)
              t=cross(m,h); y=dot_product(t,t)
              if (y.gt.errmx_l*p%ms(n)*p%ms(n)) then; errmx_l=y/(p%ms(n)*p%ms(n)); ave_tq(1:3,num_ave_i,k1)=t(1:3)/p%ms(n); endif
              h=h/1.d4  !need it Telsa
              call advance_FP(p%fp(n),h(1),h(2),h(3),p%temp(n) / (1.d4*gmub_kb * p%spin) ,p%tc(n)/(1.d4*gmub_kb*p%spin),dt)
              p%mag(2-id/3,n)=sqrt(h(1)*h(1)+h(2)*h(2)+h(3)*h(3))
              p%m(4-id:6-id,n) = h/(p%mag(2-id/3,n)+1.d-66)
              d=p%mag(2-id/3,n)*p%m(4-id:6-id,n)-p%mag(1+id/3,n)*p%m(1+id:3+id,n)
              errmx1_l=max(errmx1_l,dot_product(d,d))
            else   !heun
              print *,'fp Heun not implemented!'; stop
            endif
           else
            p%m(4-id:6-id,n)=p%m(1+id:3+id,n)
            p%mag(2-id/3,n)=p%mag(id/3+1,n)
            p%dm(:,n)=0.d0
           endif
           mmag1 = mmag1 +p%mag(2-id/3,n)*p%grain_vol(n)
          enddo
   end subroutine

   subroutine push_m(t0)
    use field_m, only: get_field
    real(8),intent(in)::t0

    logical::too_far_p
    real(8)::dt,dt_s,mmag1,rot_ratio,errmx_l,errmx1_l
    integer::k1,iter,n
    type(llg_s),pointer::p

    too_far_p=.true.  !did any spin rotate too much?
    do while (too_far_p)
      dt=llgp%dt * llgp%max_dmdt   !guess for time step
      dt_s=1.d0; iter=1; if (llgp%type.lt.0) then; iter=2; if (llgp%type.gt.-3) dt_s=0.5d0; endif
      do while (iter.gt.0)  !iter starts at 1 or 2 and decreases
        llgp%errmx=0.d0; llgp%errmx1=0.d0; llgp%mmag=0.d0

        llgl(:)%errl=.false.; llgl(:)%err1l=.false.; too_far_p=.false.; llgp%max_dmdtnext=1.d22
        do k1=1,size(llgl); p=>llgl(k1)

          errmx_l=0.d0; errmx1_l=0.d0; mmag1=0.d0
          if (llgp%zhang) then
            call push_zhang(k1,iter,errmx_l,errmx1_l,mmag1,dt*dt_s,dt)
          elseif (llgp%llb) then
            call push_llb(k1,iter,errmx_l,errmx1_l,mmag1,dt*dt_s,dt)
          elseif (llgp%fp) then
            call push_fp(k1,errmx_l,errmx1_l,mmag1,dt*dt_s,dt)
          elseif (llgp%mt) then
            call push_mt(k1,iter,errmx_l,errmx1_l,mmag1,dt*dt_s,dt)
          else
            call push_llg(k1,iter,errmx_l,errmx1_l,mmag1,dt*dt_s,dt)
          endif
          if (iter.eq.1) then  !last (only) step
            if (llgp%fixed_dt_p) then
              rot_ratio=1.d0; too_far_p=.false.
            else
              rot_ratio = acos(1.d0-0.5d0*errmx1_l)/p%max_dm;
!             rot_ratio = max(0.1d0*p%max_dm,acos(1.d0-0.5d0*errmx1_l))/p%max_dm;
              too_far_p = (rot_ratio.gt.1.d0)
            endif
            if (too_far_p) then
              llgp%max_dmdt = 0.9d0 * llgp%max_dmdt / rot_ratio
              exit
            else
!             llgp%mmag=llgp%mmag+mmag1/sum(p%grain_vol)
              llgp%mmag=llgp%mmag+mmag1/sum(p%grain_vol,mask=p%do_grain(:))
              if (any(p%hthm.ne.0.d0)) then !thermal fields present
                ave_tq(1:3,0,k1)=ave_tq(1:3,0,k1)+ave_tq(1:3,num_ave_i,k1)
                errmx_l=dot_product(ave_tq(:,0,k1),ave_tq(:,0,k1))/min(num_ave,llgp%thermal_num_ave)/min(num_ave,llgp%thermal_num_ave)
              endif
              if (rot_ratio.lt.2.d0/3.d0) then
!               llgp%max_dmdtnext=min(llgp%max_dmdtnext,llgp%max_dmdt*min(1.5d0,1.d0/(rot_ratio+1.d-22)))
!               llgp%max_dmdtnext=min(llgp%max_dmdtnext,min(1.d0,llgp%max_dmdt/(rot_ratio+1.d-200)))
!               llgp%max_dmdtnext=min(llgp%max_dmdtnext,min(1.5d0,llgp%max_dmdt/(rot_ratio+1.d-200)))
                llgp%max_dmdtnext=min(llgp%max_dmdtnext,min(2.d0,llgp%max_dmdt/(rot_ratio+1.d-200)))
!             if (rot_ratio.lt.0.5d0) then
!               llgp%max_dmdtnext=min(llgp%max_dmdtnext,min(1.d0,2.d0*llgp%max_dmdt))
              else
                llgp%max_dmdtnext=min(llgp%max_dmdtnext,llgp%max_dmdt)
              endif
!             p%errl=(errmx_l.gt.p%err*p%err); p%err1l=(errmx1_l.gt.p%err1)
              p%errl=(errmx_l.ge.p%err*p%err); p%err1l=(errmx1_l.ge.p%err1)
              llgp%errmx=max(llgp%errmx,sqrt(errmx_l)); llgp%errmx1=max(llgp%errmx1,errmx1_l)
            endif
          endif
        enddo  !loop over layers
        if (too_far_p) exit
        iter=iter-1
        if (iter.gt.0) then !going to do another step
          do k1=1,size(llgl); llgl(k1)%in=3-llgl(k1)%in; enddo
          call get_field(t0+dt*dt_s,diff_m=diff_m); dt_s=1.d0
          do k1=1,size(llgl); llgl(k1)%in=3-llgl(k1)%in; enddo
        endif
      enddo !loop over iterations (1 or 2)
    enddo  !loop until we have a valid timestep
    do k1=1,size(llgl); llgl(k1)%in=3-llgl(k1)%in; enddo
    if(llgp%max_dmdtnext.gt.1.d20) llgp%max_dmdtnext=llgp%max_dmdt
   end subroutine

  end subroutine

  subroutine plot_llg(t,it,done_p)
   use graphix, only: draw
   use plot_m, only: get_plot,plot_mk,plot_contour,plot_m3,plot_mz_scaling,plot_histogram
   use construct_grain_m, only: get_k_nom
   real(8),intent(in)::t
   integer,optional::it
   logical,optional::done_p

   integer::k,j,i
   character(200)::os
   real(8)::x(20),x1,x0,dx,y(20),x2,s,sx,sy,sxx,sxy,d,stt,a,b
     if (get_plot()) then
       do k=1,size(llgl)
         if (get_plot(k,mk=.true.)) then
           write(os,"('m.k for layer',i3)")k
           if (present(done_p)) then
             call plot_mk(trim(os),llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),get_k_nom(k),k)
           else
             call plot_mk(trim(os),llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),get_k_nom(k),k,llgl(k)%ms_scale)
           endif
         endif
         if (get_plot(k,ms=.true.)) then; write(os,"('Ms for layer',i3)")k; call plot_contour(trim(os),llgl(k)%ms*llgl(k)%ms_scale,.true.,k,no_write=done_p); endif
         if (get_plot(k,mag=.true.)) then; write(os,"('|m| for layer',i3)")k; call plot_contour(trim(os),llgl(k)%mag(llgl(k)%in/3+1,:),.true.,k,no_write=done_p); endif
         if (get_plot(k,mx=.true.)) then; write(os,"('mx for layer',i3)")k; call plot_contour(trim(os),llgl(k)%m(1+llgl(k)%in,:),.true.,k,no_write=done_p); endif
         if (get_plot(k,my=.true.)) then; write(os,"('my for layer',i3)")k; call plot_contour(trim(os),llgl(k)%m(2+llgl(k)%in,:),.true.,k,no_write=done_p); endif
         if (get_plot(k,mz=.true.)) then; write(os,"('mz for layer',i3)")k; call plot_contour(trim(os),llgl(k)%m(3+llgl(k)%in,:),.true.,k,no_write=done_p); endif
         if (get_plot(k,scatter=.true.).or.get_plot(k,scatter_len=.true.)) then
           write(os,"('scatter for layer',i3)")k
           if (.not.associated(llgl(k)%msc)) then; allocate(llgl(k)%msc(3,size(llgl(k)%m,2))); llgl(k)%msc=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:); endif
!          if (it.eq.0)  &
!                  call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:))
           if (present(it).and.it.eq.0) then
         call draw(title=trim(os)//' trace',scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%msc)
         call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%msc)
!        call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:))
           else
         call draw(title=trim(os)//' trace',scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%msc,add=.true.)
         call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%msc,add=.false.)
           endif
         llgl(k)%msc=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:)
!        call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:),add=.false.)
!          call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:),add=(it.gt.0))
!           if (it.eq.0) then
!           call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:))
!                   print *,it
!           else
!           call draw(title=trim(os),scatter=llgl(k)%m(1+llgl(k)%in:3+llgl(k)%in,:),scatter2=llgl(k)%m(4-llgl(k)%in:6-llgl(k)%in,:),add=.true.)
!           endif
!          write(os,"('scatter for layer',i3)")k; call plot_m3(trim(os),llgl(k)%m,llgl(k)%mag,llgl(k)%ms_scale,llgl(k)%in,k)
!          write(os,"('scatter for layer',i3)")k; call plot_m3(trim(os),llgl(k)%m,llgl(k)%mag,llgl(k)%ms_scale,llgl(k)%in,k,add=.true.)
!          write(os,"('trajectory for layer',i3)")k; call plot_m3(trim(os),llgl(k)%m,llgl(k)%mag,llgl(k)%ms_scale,llgl(k)%in,k,add=(it.ne.1))
         endif

         if (.false.) then
           write(os,"('cos(theta) histogram for layer',i3)")k;
           call plot_histogram(trim(os),20,llgl(k)%m(3+llgl(k)%in,:))
           write(os,"('theta histogram for layer',i3)")k;
           call plot_histogram(trim(os),20,acos(llgl(k)%m(3+llgl(k)%in,:))*180/acos(-1.d0))
           x=0.d0; x0 = minval(llgl(k)%m(3+llgl(k)%in,:)); x1 = maxval(llgl(k)%m(3+llgl(k)%in,:)); dx = (x1-x0)/size(x)
           do i=1,size(llgl(k)%grain_vol)
             j=min(dble(size(x)),max(1.d0,1.d0+(llgl(k)%m(3+llgl(k)%in,i)-x0)/dx))
             x(j)=x(j)+1.d0!/size(llgl(k)%grain_vol)
           enddo
           do i =1,size(x)
             if (x(i).ne.0.d0) x(i)=log(x(i))
             y(i)=x0+(i-0.5d0)*dx
           enddo
           call draw('ln(N) vs cos(theta)',x=y,y=x,ymin=minval(x),ymax=maxval(x),symbol=.true.,nps=10)
           x=0.d0; y=0.d0
           s=112.1979487572d0
           x1=s*(1-x0*x0); x0=0.d0; dx = (x1-x0)/size(x)
           do i=1,size(llgl(k)%grain_vol)
             x2 = s*(1.d0-llgl(k)%m(3+llgl(k)%in,i)*llgl(k)%m(3+llgl(k)%in,i))
             j=min(dble(size(x)),max(1.d0,1.d0+(x2-x0)/dx))
             x(j)=x(j)+1.d0!/size(llgl(k)%grain_vol)
           enddo
           s=0.d0; sx=0.d0; sy=0.d0; sxx=0.d0; sxy=0.d0; stt=0.d0
           do i =1,size(x)
             x2=x(i)
             if (x(i).ne.0.d0) x(i)=log(x(i))
             y(i)=x0+(i-0.5d0)*dx
             s=s+x2; sx=sx+y(i)*x2; sy=sy+x(i)*x2; sxx=sxx+y(i)*y(i)*x2; sxy=sxy+x(i)*y(i)*x2
           enddo
           a=0.d0; b=0.d0
           do i =1,size(x)
             x2=exp(x(i)/2)*(y(i)-sx/s)
             stt=stt+x2*x2
             b=b+x2*x(i)*exp(x(i)/2)
           enddo
           call draw('ln(N) vs energy/kT',x=y,y=x,xmin=0.d0,xmax=x1,ymin=minval(x),ymax=maxval(x)*1.1d0,symbol=.true.,nps=10)
           b=b/stt
           a=(sy-sx*b)/s
           write(os,"('m =',f10.6,' +- ',f8.6)")b,sqrt(1.d0/stt)
           call draw('ln(N) vs energy/kT',add=.true.,x=(/x0,x1/),y=(/a+b*x0,a+b*x1/),legend=trim(os))
         endif

       !    d=s*sxx-sx*sx
       !    b=(s*sxy-sx*sy)/d !b
       !    a=(sxx*sy-sx*sxy)/d ! a; y=a-bx
       !    write(os,"('m=',1p,e12.5,' +- ',e12.5)")b,sqrt(s/d)
!          print *,x

!      call plot_mz_scaling(k,t,(/sum(llgl(k)%mag(llgl(k)%in/3+1,:)*llgl(k)%m(3+llgl(k)%in,:))/size(llgl(k)%m(3,:)),-1.d3002/))
!      call plot_mz_scaling(k,t,(/sum(llgl(k)%mag(llgl(k)%in/3+1,:)*llgl(k)%m(3+llgl(k)%in,:))/size(llgl(k)%m(3,:)),0.d0/))
       enddo
     endif
   end subroutine


  function init_llg(il) result (ok_p)
   use data_m, only: get_data
   integer,intent(in),optional::il
   logical::ok_p
   integer::i1,i2,i

   i1=1; i2=size(llgl)
   if (present(il)) then
     if (il.gt.0.and.il.le.i2) then; i1=il; i2=il; endif
   endif
   do i=i1,i2
     call get_data(i,llgl(i)%do_grain,do_grain_p=.true.)
     call get_data(i,llgl(i)%min_n,min_n_p=.true.)
     call get_data(i,llgl(i)%max_n,max_n_p=.true.)
     call get_data(i,llgl(i)%alpha,alpha_p=.true.)
     call get_data(i,llgl(i)%gamma,gamma_p=.true.)
     call get_data(i,llgl(i)%in,in_p=.true.)
     call get_data(i,llgl(i)%m,m_p=.true.)
     call get_data(i,llgl(i)%ms,ms_p=.true.)
     call get_data(i,llgl(i)%lambdaMs,lambdaMs_p=.true.)
     call get_data(i,llgl(i)%spin,spin_p=.true.)
     call get_data(i,llgl(i)%h,h_p=.true.)
     call get_data(i,llgl(i)%hthm,hthm_p=.true.)
     call get_data(i,llgl(i)%hthm_m,hthm_p=.true.)
     call get_data(i,llgl(i)%mag,mag_p=.true.)
     call get_data(i,llgl(i)%temp,temp_p=.true.)
     call get_data(i,llgl(i)%tc,tc_p=.true.)
     call get_data(i,llgl(i)%dm,dm_p=.true.)
     call get_data(i,llgl(i)%alpha_scale,alpha_scale_p=.true.)
     call get_data(i,llgl(i)%ms_scale,ms_scale_p=.true.)
     call get_data(i,llgl(i)%grain_vol,grain_vol_p=.true.)
     call get_data(i,llgl(i)%fp,fp_p=.true.)
     llgl(i)%scale=-llgl(i)%gamma/(2.d0*(1.d0+llgl(i)%alpha*llgl(i)%alpha))
   enddo
   ok_p=.true.
  end function

  subroutine initialize_fp()
    use FP_m, only: FP_s, new_FP, advance_FP, null_fp
    use data_m, only:set_data
    integer::n,i
    real(8)::dt,mfp(3),lambda

    dt=1.d-11
    do i=1,size(llgl)
      do n=1,size(llgl(i)%grain_vol)
        if(null_FP(llgl(i)%fp(n))) call new_fp(llgl(i)%fp(n))
        mfp=llgl(i)%h(:,n)/1.d4; lambda=llgl(i)%lambdaMs(n)*(llgl(i)%spin+1.d0)/llgl(i)%spin/1.d4;
        call advance_FP(llgl(i)%fp(n),mfp(1),mfp(2),mfp(3),llgl(i)%temp(n),lambda,dt)
        llgl(i)%mag(1+llgl(i)%in/3,n)=sqrt(dot_product(mfp,mfp))
        llgl(i)%m(1+llgl(i)%in:3+llgl(i)%in,n)=mfp/llgl(i)%mag(1+llgl(i)%in/3,n)
      enddo
     call set_data(i,llgl(i)%fp,fp_p=.true.)
    enddo
    need_to_initialize=.false.
  end subroutine

  function start_llg(num_lay) result (ok_p)
   integer,intent(in)::num_lay
   logical::ok_p
   if (.not.associated(llgp)) allocate(llgp)
   if (.not.associated(llgl)) allocate(llgl(num_lay))
!  if (.not.allocated(llgp)) allocate(llgp)
!  if (.not.allocated(llgl)) allocate(llgl(num_lay))
   ok_p=.true.
  end function


  subroutine set_llg (il,dt,max_dm,min_dm,dt_frac,dt_max, t_max,t_min,del_t,conv_torque,conv_angle, &
           fixed_dt_p,iter,min_ave_m,max_ave_m,plot_iter,thermal_num_ave,llg_type,mt,meq_bs,zhang,llb,fp,indx,garanin_tau,min_it,traj)
   integer,intent(in)::il
   logical,intent(in),optional::fixed_dt_p,mt,meq_bs,zhang,llb,fp,garanin_tau
   integer,intent(in),optional::iter,plot_iter,thermal_num_ave,llg_type,indx,min_it,traj
   real(8),intent(in),optional::dt,max_dm,min_dm,dt_frac,dt_max,t_max,t_min,del_t,conv_torque,conv_angle,min_ave_m(3),max_ave_m(3)
   integer::i1,i2,i
   real(8),save::dt_frac_old=0.2d0
   if (present(traj)) NTRAJ=traj
   if (present(llb)) llgp%llb=llb;                      if (present(zhang)) llgp%zhang=zhang
   if (present(fp)) llgp%fp=fp;                         if (present(garanin_tau)) llgp%garanin_tau=garanin_tau
   if (present(iter)) llgp%iter=iter;                   if (present(mt)) llgp%mt=mt; if (present(meq_bs)) llgp%meq_bs=meq_bs
   if (present(plot_iter)) llgp%plot_iter=plot_iter;    if (present(thermal_num_ave)) llgp%thermal_num_ave=max(1,thermal_num_ave)
   if (present(dt)) llgp%dt=dt;                         if (present(dt_max)) llgp%dt_max=dt_max
   if (present(t_max)) llgp%t_max=t_max;                if (present(t_min)) llgp%t_min=t_min
   if (present(del_t)) llgp%del_t=del_t;                if (present(llg_type)) llgp%type=llg_type       
   if (present(dt_frac)) then; if (llgp%dt_frac.ne.0.d0) dt_frac_old=llgp%dt_frac; llgp%dt_frac=dt_frac; endif
   if (present(fixed_dt_p)) then; llgp%fixed_dt_p=fixed_dt_p; if (.not.fixed_dt_p.and.llgp%dt_frac.eq.0.d0) llgp%dt_frac=dt_frac_old; endif
   if (present(min_it)) llgp%min_it=min_it
   i1=1; i2=size(llgl)
   if (il.gt.0.and.il.le.i2) then; i1=il; i2=il; endif
   do i=i1,i2
     if (present(max_dm)) llgl(i)%max_dm=max_dm*acos(-1.d0)/180.d0
     if (present(min_dm)) llgl(i)%min_dm=min_dm*acos(-1.d0)/180.d0
     if (present(conv_torque)) llgl(i)%err=abs(conv_torque); if (present(conv_angle)) llgl(i)%err1=abs(conv_angle)
     if (present(indx)) then
       if (present(min_ave_m)) llgl(i)%minm(indx)=min_ave_m(indx)
       if (present(max_ave_m)) llgl(i)%maxm(indx)=max_ave_m(indx)
     else
       if (present(min_ave_m)) llgl(i)%minm=min_ave_m;       if (present(max_ave_m)) llgl(i)%maxm=max_ave_m
     endif
   enddo
   if (llgp%fixed_dt_p) llgp%dt_frac=0.d0
   if (llgp%dt_frac.eq.0.d0) llgp%fixed_dt_p=.true.
   if (llgp%fixed_dt_p) llgp%max_dmdt=1.d0
  end subroutine
end module llg_m
