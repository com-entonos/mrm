module decay_m
 use data_m, only: data_s
 use random_m, only: random_s
 use graphix, only:draw
 implicit none
 private

 public::decay_mag,set_decay_mag

! real(8)::energy0

 integer,parameter::nitersaddle=2
 integer,save::NTRAJ=20,TRAJ=0
 real(8),save::torbit=1.d-9
 logical,save::boltz_p=.true.,diff_m=.true.,test_rot_p=.true.
!integer,parameter::NTRAJ=20,TRAJ=0
!integer,parameter::NTRAJ=20,TRAJ=313 !400 !313 !400 !592
!integer,parameter::NTRAJ=20,TRAJ=256 !340 !290 !63 !124 !0 !63 !8 !568 !63 !362 !326 !169 !307 !28 !52 !224 !643

 type::decay_s
   logical,dimension(:),pointer::traj_dump=>null(),do_grain=>null()
   integer,dimension(:),pointer::traj_lch=>null(),num_flip=>null()
   real(8),dimension(:),pointer::gamma=>null(),gamma1=>null()
   real(8),dimension(:,:),pointer::m=>null(),energy=>null(),trj_energy=>null(),trj_mxh2=>null(),trj_lm=>null()
   real(8),dimension(:,:,:),pointer::m_min=>null(),bltz=>null(),trj_mstart=>null(),trj_m0=>null(),bltz_coef=>null()
   real(8),dimension(:,:,:,:),pointer::bltz_m0=>null()
   type(data_s),pointer::layer=>null()
 end type decay_s

 type(decay_s),allocatable,dimension(:),target::decay
 type(random_s),pointer,save::random_decay=>null(),random_orbit=>null()

 real(8)::last_flip

 real(8),parameter::kb=1.380662d-16

 real(8),save::VOL_FRAC=0.01d0     !during one MC, fraction of grain volume to flip before recalculating transition rates
 real(8),save::max_eb=50.d0  !max energy barrier in units of kT
 real(8),save::min_eb=1.d0 !0.1d0 !1.d0 !min energy barrier to worry about   1/gmmaa ~ (exp(e/kT) - 1 - e/kT) / (a/kT)
 real(8),save::last_cum_prob=0.d0
 integer,save::iter=0

contains
  function set_decay_mag(itorbit,intraj,itraj,iiter,imax_eb,iboltz,ivol_frac,ilast_cum_prob,idiff_m,itest_rot_p) result (ok_p)
   real(8),intent(in),optional::itorbit,imax_eb,ivol_frac,ilast_cum_prob
   integer,intent(in),optional::intraj,itraj,iiter
   logical,intent(in),optional::iboltz,idiff_m,itest_rot_p
   logical::ok_p
   if (present(itorbit)) torbit=itorbit
   if (present(iboltz)) boltz_p=iboltz
   if (present(itest_rot_p)) test_rot_p=itest_rot_p
   if (present(idiff_m)) diff_m=idiff_m
   if (present(intraj)) NTRAJ=intraj
   if (present(itraj)) TRAJ=itraj
   if (present(ilast_cum_prob)) last_cum_prob=ilast_cum_prob
   if (present(imax_eb)) max_eb=imax_eb
   if (present(ivol_frac)) VOL_FRAC=ivol_frac
   if (present(iiter)) iter=iiter
   ok_p=.true.
  end function
!  function catch_falls(t0,ntrys,alpha,grain_do,n_lay,ng) result (ok_p)
!   use llg_m, only: llg, backup_llg, restore_llg, set_llg, show_llg
!
!   real(8),intent(in)::t0
!   integer,intent(in)::ntrys,n_lay,ng
!   logical,intent(in)::grain_do(ng)
!   real(8),intent(in)::alpha(n_lay)
!   logical::ok_p,torig(n_lay)
!   integer::i,n,it,mflip(ng)
!   real(8)::morig(3,ng,n_lay),m(3,ng,n_lay),t,aorig(n_lay)
!
!   do i=1,n_lay
!     morig(:,:,n_lay)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)
!     torig(n_lay)=decay(i)%layer%thermal_p
!     aorig(n_lay)=decay(i)%layer%alpha
!   enddo
!call show_llg(t0); read *
!   call backup_llg(); t=t0
!
!   do it=1,ntrys
!
!   !load up m
!     do i=1,n_lay
!       do n=1,ng
!         if (grain_do(n)) then
!           m(:,n,n_lay)=morig(:,n,i)
!         else
!           m(:,n,n_lay)=decay(i)%m_min(:,n,3)
!         endif
!       enddo
!       decay(i)%layer%thermal_p=.true.
!     enddo
!
!   !set up llg w/ thermals and do it
!     call set_llg(0,conv_torque=0.d0,t_max=t0+2.d-9,mt=.false.,zhang=.false.,llb=.false.,llg_type=1,dt=20.d-14,del_t=1.d22,fixed_dt_p=.true.,min_dm=0.d0)   !only use LLG/LL
!call show_llg(t0); read *
!     t=t0; call llg(t)
!
!   !turn off thermal and just let the grains settle
!     do i=1,n_lay
!       decay(i)%layer%thermal_p=.false.
!     enddo
!     call set_llg(0,conv_torque=0.0003d0,t_max=t0+1d22,llg_type=2,dt=20.d-14,del_t=1.d22,fixed_dt_p=.true.,min_dm=0.d0)   !only use LLG/LL
!call show_llg(t0); read *
!     t=t0;call llg(t)
!
!    !collect stats on who flipped
!     do i=1,n_lay
!       do n=1,ng
!         if (grain_do(n).and.dot_product(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n),morig(:,n,i)).lt.0.d0) mflip(n)=mflip(n)+1
!       enddo
!     enddo
!   enddo
!
!   do i=1,n_lay
!     decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=morig(:,:,n_lay)
!     decay(i)%layer%thermal_p=torig(n_lay)
!     decay(i)%layer%alpha=aorig(n_lay)
!   enddo
!   call restore_llg()
!call show_llg(t0); read *
!
!   ok_p=.true.
!  end function

  subroutine decay_mag(t0,dt,temperature,dh_dti,happliedi,mkti,old_gamma_pi,multispin_pi,iseedi)
   use io_m, only: output,delete_file,add_file,show_file
   use random_m, only: random_create, random_destroy
   use data_m, only: get_layer, get_data, set_data
   use llg_m, only: llg, backup_llg, restore_llg, set_llg
   use field_m, only: freeze_time, set_mode_field, get_field
   use track_m, only: track_s, init_track, destroy_track, set_track
   real(8),intent(in)::t0,dt,temperature
   real(8),intent(in),optional::dh_dti,happliedi,mkti
   logical,intent(in),optional::old_gamma_pi,multispin_pi
   integer,intent(in),optional::iseedi
   real(8),dimension(:),pointer::ms_energy=>null()
   real(8),allocatable,dimension(:,:),target::msdata

   type(data_s),pointer::layer
   type(track_s),pointer::track

   logical::err_p, old_gamma_p, multispin_p, flip_ms
   integer::i,n,id,n_lay,flipped

   real(8)::total_time,x,total_rate,mkt,dh_dt,happlied,t
   real(8),allocatable,dimension(:)::old_alpha,old_temperature
   logical,allocatable,dimension(:)::old_thermal
   character(400)::os,dirs

!  print "('t0,dt,temperature: ',1p,:,50(e12.5,x))",t0,dt,temperature
!  if (present(dh_dti))print "('dh_dti=',1p,e12.5)",dh_dti
!  if (present(happliedi))print "('happliedi=',1p,e12.5)",happliedi
!  if (present(mkti))print "('mkti=',1p,e12.5)",mkti
!  if (present(old_gamma_pi))print "('old_gamma_pi=',l)",old_gamma_pi
!  if (present(multispin_pi))print "('mputispin_pi=',l)",multispin_pi
!  if (present(iseedi))print "('iseedi=',i0)",iseedi
!  print *,associated(random_orbit),associated(random_decay)

!   do i=1,10
!     write(os,"(i)") i; os='file'//trim(adjustl(os))//'.dat'
!     n=add_file(trim(os)); !call delete_file(n)
!     print *,n,' ',trim(os)
!   enddo
!   call show_file()
!   read *

   if (.not.associated(random_decay)) then
     i=4357; if (present(iseedi)) i=i+iseedi; random_decay=>random_create('random_decay',i)
   endif

   old_gamma_p=.false.; if (present(old_gamma_pi)) old_gamma_p=old_gamma_pi
   multispin_p=.false.; if (present(multispin_pi)) multispin_p=multispin_pi
   mkt=0.d0; if (present(mkti)) mkt=mkti
   dh_dt=0.d0; if (present(dh_dti)) dh_dt=dh_dti
   happlied=0.d0; if (present(happliedi)) happlied=happliedi

   i=1; layer=>get_layer(i)
   do while (associated(layer))
     i=i+1; layer=>get_layer(i)
   enddo; n_lay=i-1; if (n_lay.lt.1) return


   allocate(decay(n_lay))
   do i=1,n_lay
     layer=>get_layer(i); if (i.eq.1) n=size(layer%grain_vol)
     multispin_p=(multispin_p.and.n.eq.size(layer%grain_vol))
     decay(i)%layer=>layer; n=size(layer%grain_vol); id=layer%in
     allocate(decay(i)%m(3,n),decay(i)%num_flip(n), &
               decay(i)%m_min(3,n,6), decay(i)%energy(4,n),decay(i)%traj_dump(n), &
               decay(i)%trj_energy(n,NTRAJ),decay(i)%trj_m0(3,n,NTRAJ),decay(i)%trj_mxh2(n,NTRAJ), &
               decay(i)%trj_lm(n,NTRAJ),decay(i)%traj_lch(n), decay(i)%gamma1(n),decay(i)%gamma(n),decay(i)%bltz_coef(nitersaddle*NTRAJ,2,n), &
               decay(i)%trj_mstart(6,n,NTRAJ),decay(i)%bltz(0:nitersaddle*NTRAJ,4,n),decay(i)%bltz_m0(3,2,0:nitersaddle*NTRAJ+1,n))
     decay(i)%m=layer%m(1+id:3+id,:); decay(i)%num_flip=0
   enddo

   call get_data(1,ms_energy,ms_energy_p=.true.); layer=>get_layer(1)
!  if (n_lay.lt.2.or..not.multispin_p)then
   allocate(msdata(size(layer%grain_vol),4)); msdata=0.d0
   if (.not.multispin_p)then
     multispin_p=.false.; if (associated(ms_energy)) deallocate(ms_energy); nullify(ms_energy)
   else
     if (n_lay.lt.2) call output('  *** warning ***  doing multispin with less than 2 layers')
     if (.not.associated(ms_energy)) allocate(ms_energy(size(layer%grain_vol)))
!    allocate(msdata(size(layer%grain_vol),4)); msdata=0.d0
   endif
   call set_data(1,ms_energy,ms_energy_p=.true.)
   if (multispin_p .and..not.associated(random_orbit)) then
     i=4357; if (present(iseedi)) i=2*(i+iseedi)+1; random_orbit=>random_create('random_orbit',i)
   endif
   allocate(old_alpha(n_lay),old_temperature(n_lay),old_thermal(n_lay))
   do i=1,n_lay; old_alpha(i)=decay(i)%layer%alpha; old_temperature(i)=decay(i)%layer%background_t; old_thermal(i)=decay(i)%layer%thermal_p;enddo
   total_time=0.d0; last_flip=0.d0

   call backup_llg()
   call freeze_time(t0)
   nullify(track)

   if (multispin_p) then
     do i=1,n_lay-1; id=0; do n=1,size(decay(i)%m,2); if (dot_product(decay(i)%m(:,n),decay(i+1)%m(:,n)).lt.0.d0) id=id+1; enddo
       write(os,"('!<broken spin between ',i0,' -',2(x,i0),' out of ',i0)") i,i+1,id,size(decay(1)%m,2); call output(trim(os))
     enddo
     id=0; do n=1,size(decay(i)%m,2); do i=1,n_lay-1; if (dot_product(decay(i)%m(:,n),decay(i+1)%m(:,n)).lt.0.d0) exit; enddo; id=id+1-i/n_lay; enddo
     write(os,"('!<broken grains ',i0,' out of ',i0)") id,size(decay(1)%m,2); call output(trim(os))
!n=TRAJ;if (TRAJ.ne.0) then; do i=1,n_lay; print "(1p,3(e12.5,x))",decay(i)%m(:,n);enddo; endif
   endif

   do while (.true.)  !loop forever
     call backup_llg()
!    call set_llg(0,llg_type=2,zhang=.false.,llb=.false.,plot_iter=10000000)
     call set_llg(0,plot_iter=10000000,t_max=t0+1.d-9)
     err_p=set_mode_field(0); dirs=''
     do i=1,n_lay
       decay(i)%layer%alpha=old_alpha(i); decay(i)%layer%background_t=old_temperature(i); decay(i)%layer%thermal_p=old_thermal(i)
       write(dirs,"(a,l)") trim(dirs),decay(i)%layer%m(3+decay(i)%layer%in,max(1,TRAJ)).gt.0.d0
     enddo
n=TRAJ;if (n.ne.0) then; write(os,"(' ###  initial state for grain ',i0,x,a,x,l)")n,trim(dirs),decay(1)%layer%do_grain(n); call output(trim(os)); endif

     call output(' '); call output('*** relaxing magnetization collective...'); call output(' ')

    !run as set up
     t=t0; call llg(t) !,trace_p=.true.) 
     if (iter.lt.0) then; do i=1,n_lay; decay(i)%m=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:); enddo; endif
!    total_time=total_time+t-t0; if (total_time.ge.dt) then; call restore_llg(); exit; endif
     if (total_time.ge.dt) then; call restore_llg(); exit; endif

n=TRAJ;if (n.ne.0) then
     dirs=''; do i=1,n_lay; write(dirs,"(a,l)") trim(dirs),decay(i)%layer%m(3+decay(i)%layer%in,n).gt.0.d0; enddo
     write(os,"(' ###      new state for grain ',i0,x,a,x,l)")n,trim(dirs),decay(1)%layer%do_grain(n); call output(trim(os))
endif
     

     do i=1,n_lay; decay(i)%m(:,:)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:); enddo

!    call set_llg(0,t_max=t0+20.d-14,mt=.false.,zhang=.false.,llb=.false.,llg_type=2,dt=20.d-14,del_t=1.d22,fixed_dt_p=.false.,max_dm=2.5d0,min_dm=0.d0)   !only use LLG/LL
     call set_llg(0,t_max=t0+20.d-14,mt=.false.,zhang=.false.,llb=.false.,llg_type=2,dt=20.d-14,del_t=1.d22,fixed_dt_p=.false.,             min_dm=0.d0)   !only use LLG/LL
     do i=1,n_lay; decay(i)%layer%background_t=temperature; decay(i)%layer%thermal_p=.false.; enddo
     t=t0; call llg(t,.true.,trace_p=.true.);

     t=t0
     call output(' '); call output('*** running magnetization to equilibrium position...')
     err_p=set_mode_field(-1); call get_field(t); err_p=set_mode_field(1)
     call output('**** w/ frozen spins...'); call output(' ')
     call set_llg(0,conv_torque=0.00003d0,t_max=1.d22)
     do i=1,n_lay; decay(i)%layer%alpha=100.d0; enddo
     call llg(t,.true.) !,trace_p=.true.)

     dirs=''
     do i=1,n_lay
       decay(i)%energy=1.d220; decay(i)%bltz=max_eb; decay(i)%bltz_m0=0.d0
       do n=1,size(decay(i)%layer%grain_vol)
         decay(i)%m_min(:,n,1)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n)
if (n.eq.TRAJ) write(dirs,"(a,l)") trim(dirs),decay(i)%m_min(3,n,1).gt.0.d0
       enddo
       decay(i)%m_min(:,:,3)=decay(i)%m_min(:,:,1)
       decay(i)%bltz(0,1,:)=decay(i)%layer%energy
       decay(i)%bltz(0,2,:)=0.d0
       decay(i)%energy(1,:)=decay(i)%layer%energy
       decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=-decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)
     enddo
     if (multispin_p) msdata(:,1)=ms_energy

n=TRAJ;if (n.ne.0) then; write(os,"(' ###    newer state for grain ',i0,x,a,x,l)")n,trim(dirs),decay(1)%layer%do_grain(n); call output(trim(os)); endif

n=TRAJ;if (multispin_p.and.TRAJ.ne.0) then; do i=1,n_lay; print "(1p,3(e12.5,x))",decay(i)%m_min(:,n,1);enddo; endif

     call output(' '); call output('*** running inverted magnetization to equilibrium position...'); call output(' ')
     t=t0
!    call set_llg(0,del_t=1.d22,conv_torque=0.00003d0); do i=1,n_lay; decay(i)%layer%alpha=100.d0; enddo
     call llg(t,.true.) !,trace_p=.true.)


     flipped=0; x=1.d22; dirs=''
     do i=1,n_lay
       do n=1,size(decay(i)%layer%grain_vol)
         decay(i)%m_min(:,n,2)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n)
         decay(i)%layer%m(4-decay(i)%layer%in:6-decay(i)%layer%in,n)=decay(i)%m_min(:,n,1)
         if (acos(dot_product(decay(i)%m_min(:,n,1),decay(i)%m_min(:,n,2))).gt.5.d0*acos(-1.d0)/180.d0) &
             x=min(x,(decay(i)%layer%energy(n)-decay(i)%energy(1,n))*decay(i)%layer%grain_vol(n)/(kb*temperature))
if (n.eq.TRAJ) write(dirs,"(a,l)") trim(dirs),decay(i)%m_min(3,n,2).gt.0.d0
       enddo
       decay(i)%m_min(:,:,4)=decay(i)%m_min(:,:,2)
       decay(i)%energy(2,:)=decay(i)%layer%energy
       flipped = flipped + count(acos(decay(i)%m_min(1,:,1)*decay(i)%m_min(1,:,2)+decay(i)%m_min(2,:,1)*decay(i)%m_min(2,:,2)+ &
            decay(i)%m_min(3,:,1)*decay(i)%m_min(3,:,2)).gt.5.d0*acos(-1.d0)/180.d0)
     enddo
n=TRAJ;if (n.ne.0) then; write(os,"(' ### inverted state for grain ',i0,x,a,x,l)")n,trim(dirs),decay(1)%layer%do_grain(n); call output(trim(os));endif

n=TRAJ;if (multispin_p.and.TRAJ.ne.0) then; do i=1,n_lay; print "(1p,3(e12.5,x))",decay(i)%m_min(:,n,2);enddo; endif

     if (multispin_p) then
       msdata(:,2)=ms_energy
       flipped=0; x=1.d300
       do n=1,size(decay(1)%layer%grain_vol)
         flip_ms=.false.
         do i=1,n_lay
           flip_ms=(flip_ms.or.dot_product(decay(i)%m_min(:,n,1),decay(i)%m_min(:,n,2)).lt.0.d0)
         enddo
         if (flip_ms) then
           flipped=flipped+1
           x=min(x,(msdata(n,2)-msdata(n,1))/(kb*temperature))
         endif
         decay(1)%traj_dump(n)=(flip_ms)
       enddo
     endif

     write(os,"(i5,' flipped, min((e0p-e0)/kT)=',1p,e12.5)") flipped,x; call output(os); call output(' ')
     if (flipped.le.0.or.x.gt.max_eb) then
       total_time=dt
       do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=decay(i)%m; enddo
       if (flipped.le.0) call output('  no grain has more than one minimum')
       if (x.gt.300.d0) call output('  minimum energy barrier > 300 kT')
       call restore_llg()
       exit
     endif
     
     call set_llg(0,traj=TRAJ)
     if (.not.associated(track)) track=>init_track(TRAJ)
     call set_track(track)  !requires local min to be in m(4-id:6-id,:) and other minimum in m(1+id:3+id,:)
     
     call set_llg(0,llg_type=-2); do i=1,n_lay; decay(i)%layer%alpha=0.d0; enddo

     if (.not.multispin_p) then; do i=1,n_lay
       do n=1,size(decay(i)%layer%grain_vol)
         decay(i)%traj_dump(n)=(acos(dot_product(decay(i)%m_min(:,n,1),decay(i)%m_min(:,n,2))).gt.5.d0*acos(-1.d0)/180.d0)
       enddo; enddo
     endif

     flipped=0; total_rate=0.d0
     if (find_energy_barrier(iter)) then
       do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=decay(i)%m; enddo
       call restore_llg()
       exit
     endif

     total_rate=find_rates()
     if (iter.lt.0) then
       total_time=0.d0
       if (multispin_p) then
         inquire(file='eb.dat',exist=err_p)
         i=add_file(trim(os)); call delete_file(i)
         if (err_p) then
           open(i,file='eb.dat',status='old',position='append')
         else
           open(i,file='eb.dat',status='new')
         endif
         do n=1,size(msdata,1)
           write(i,"(i4,1p,:,40(x,e12.5))")n,decay(1)%bltz(NTRAJ,3,n),decay(1)%gamma(n),msdata(n,1)/(kb*temperature),msdata(n,2)/(kb*temperature)
         enddo
         close(i)
         iter=iter+1; if (iter.ge.0) stop
!read *
       else
               print *,'energy barrier distribution not implemented for .not. multispin'; stop
       endif
     else
print *,'total_rate: ',total_rate
       if (total_rate.le.0.d0) then
         total_time=dt
       else
         call kinetic_monte_carlo(total_rate)
       endif
     endif

     do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=decay(i)%m; enddo
     write(os,"('total elasped time:',1p,e12.5,'(s), ',e12.5)")total_time,total_time/dt; call output(os)
     call restore_llg()
 !read *
   enddo
!write(os,"('cum_prob ',4(e12.5,x))")last_cum_prob,last_cum_prob+1.d0-exp(-total_rate*(dt-last_flip)); call output(os)
!   last_cum_prob=last_cum_prob+1.d0-exp(-total_rate*(dt-last_flip))
if (total_rate.ne.0.d0) then
write(os,"('cum_prob ',:,4(e12.5,x))") last_cum_prob,last_cum_prob*exp(total_rate*(dt-last_flip)); call output(os)
   last_cum_prob=last_cum_prob*exp(total_rate*(dt-last_flip))
else
write(os,"('cum_prob ',:,4(e12.5,x))") last_cum_prob; call output(os)
endif

   if (multispin_p) then
     do i=1,n_lay-1; id=0; do n=1,size(decay(i)%m,2); if (dot_product(decay(i)%m(:,n),decay(i+1)%m(:,n)).lt.0.d0) id=id+1; enddo
       write(os,"('!>broken spin between ',i0,' -',2(x,i0),' out of ',i0)") i,i+1,id,size(decay(1)%m,2); call output(trim(os))
     enddo
     id=0; do n=1,size(decay(i)%m,2); do i=1,n_lay-1; if (dot_product(decay(i)%m(:,n),decay(i+1)%m(:,n)).lt.0.d0) exit; enddo; id=id+1-i/n_lay; enddo
     write(os,"('!>broken grains ',i0,' out of ',i0)") id,size(decay(1)%m,2); call output(trim(os))
   endif

   if (associated(track)) track=>destroy_track(track)
   call restore_llg()
   call freeze_time(-1.d300); err_p=set_mode_field(0)
   do i=1,n_lay
     decay(i)%layer%alpha=old_alpha(i); decay(i)%layer%background_t=old_temperature(i); decay(i)%layer%thermal_p=old_thermal(i)
     deallocate(decay(i)%m,decay(i)%num_flip,decay(i)%m_min,decay(i)%energy, &
             decay(i)%traj_dump,decay(i)%trj_energy,decay(i)%trj_m0,decay(i)%trj_mxh2,decay(i)%trj_lm,decay(i)%traj_lch, &
             decay(i)%gamma1,decay(i)%gamma,decay(i)%bltz_coef,decay(i)%trj_mstart,decay(i)%bltz,decay(i)%bltz_m0)
   enddo
   deallocate(decay,old_alpha, old_thermal, old_temperature); if (allocated(msdata)) deallocate(msdata)
   if (associated(ms_energy)) then; deallocate(ms_energy); nullify(ms_energy); call set_data(1,ms_energy,ms_energy_p=.true.); endif
   if (associated(random_decay)) call random_destroy(random_decay)
   if (associated(random_orbit))  call random_destroy(random_orbit)

  contains
    function damp(ng,dotraj,alpha,move,max_dth,conv_dth) result (it)
     use graphix, only:draw
     use transform_m, only: cross, rotate
     real(8),intent(in)::max_dth,alpha(n_lay),conv_dth
     integer,intent(in)::ng
     logical,intent(in)::dotraj(ng),move(n_lay,ng)
     integer::k,n,it
     logical::dog(ng)
     real(8)::hmag(n_lay),th0(n_lay),dth,min_dt,x
     real(8)::m(3,n_lay),h(3,n_lay) !,mm(3)
     real(8)::msdata(3,n_lay+1,2)

!if (TRAJ.ne.0) print "(a,1p,4(e12.5,x))",'in damp() ',ms_energy(TRAJ),decay(4)%layer%m(1+decay(4)%layer%in:3+decay(4)%layer%in,TRAJ)
     it=-1
     dog=dotraj
!    do while (count(dog).gt.0.or.it.lt.20) ; it=it+1
do while ((count(dog).gt.0.or.it.lt.20).and.it.lt.1000) ; it=it+1
       do n=1,ng
         if (dog(n)) then
           min_dt=1.d300; th0=0.d0
           do k=1,n_lay
             m(:,k)=decay(k)%layer%m(1+decay(k)%layer%in:3+decay(k)%layer%in,n)
if (n.eq.TRAJ) msdata(:,k,2)=m(:,k)
if (n.eq.TRAJ) msdata(:,k,1)=m(:,k)
             if (move(k,n)) then
               h(:,k)=decay(k)%layer%h(:,n)
               hmag(k)=dot_product(h(:,k),h(:,k))
               if (hmag(k).gt.0.d0) then 
                 hmag(k)=sqrt(hmag(k)); h(:,k)=h(:,k)/hmag(k)
                 th0(k)=acos(max(-1.d0,min(1.d0,dot_product(m(:,k),h(:,k)))))
                 if (th0(k).gt.max_dth) then
                   min_dt = min(min_dt,log(max(1.d0,tan(0.5d0*th0(k))/tan(0.5d0*(th0(k)-max_dth))))*(1.d0+alpha(k)*alpha(k))/(hmag(k)*alpha(k)+1.d-300))
                 elseif (th0(k).gt.0.d0) then
                   min_dt = min(min_dt,log(max(1.d0,tan(0.5d0*th0(k))/tan(0.25d0*th0(k))))*(1.d0+alpha(k)*alpha(k))/(hmag(k)*alpha(k)+1.d-300))
                 endif
               endif
!print *,n,k,move(k),h(:,k),hmag(k),th0(k),min_dt
             endif
           enddo
!if (n.eq.TRAJ.and.all(move)) print "(1p,:,40(e12.5,x))",min_dt,th0*180.d0/acos(-1.d0),hmag
           if (min_dt.lt.1.d295) then
             x=0.d0
             do k=1,n_lay
               if (move(k,n).and.hmag(k).ne.0.d0) then
                 dth=th0(k)-2.d0*atan(exp(-min_dt*alpha(k)*hmag(k)/(1.d0+alpha(k)*alpha(k)))*tan(0.5d0*th0(k)))
                 x=max(x,dth)
                 decay(k)%layer%m(1+decay(k)%layer%in:3+decay(k)%layer%in,n) = rotate(dth,cross(m(:,k),h(:,k)),m(:,k),.true.)
if (n.eq.TRAJ) msdata(:,k,1)= decay(k)%layer%m(1+decay(k)%layer%in:3+decay(k)%layer%in,n)
               endif
             enddo
             if (it.gt.19) dog(n)=(maxval(th0).gt.conv_dth*0.99d0)
if (n.eq.TRAJ) then
 msdata(:,n_lay+1,2)=msdata(:,n_lay+1,1); msdata(:,n_lay+1,1)=(/ sum(msdata(1,1:n_lay,1)), sum(msdata(2,1:n_lay,1)), sum(msdata(3,1:n_lay,1)) /)/n_lay
 k=n_lay+min(1,n_lay-1)
 if (it.ne.0) then
 call draw(title='M trace',scatter=msdata(:,1:k,1),scatter2=msdata(:,1:k,2),add=.true.)
 call draw(title='M', scatter=msdata(:,1:k,1),scatter2=msdata(:,1:k,2),add=.false.)
 else
 call draw(title='M trace',scatter=msdata(:,1:k,1),scatter2=msdata(:,1:k,2))
 call draw(title='M', scatter=msdata(:,1:k,1),scatter2=msdata(:,1:k,2))
 endif
endif
           endif
         endif
       enddo
       call get_field(t0) !get energy & field
     enddo
    end function

    function find_energy_barrier_ms(ng,iter) result (err_p)
     use io_m, only: output
     use random_m, only: random_uniform
     use field_m, only:get_field
     use transform_m, only:cross,rotate
     use track_m, only: reset_track, dump_track 
     use data_m, only: get_data
     use graphix, only:draw
     use plot_m, only: plot_contour,plot_histogram,get_plot,set_plot
     integer,intent(in)::ng,iter

     logical::err_p,move(n_lay,ng),broke_grain(n_lay),flipd(ng),plotting
     integer::im(ng)
     logical,pointer::grain_do(:),grain_do_false(:),grain_do_copy(:)
     real(8)::mo(3,ng,n_lay,0:5,0:1),me(ng,0:5,0:1),dangle(ng),u(3),m0(3),x
     integer::i,j,n,n1,n2,n3,olay,olay1,it
     type(track_s),pointer::tracki
     real(4)::time1,time2,time3

     integer,pointer,save::layer_order(:)=>null()

     type::mo_s
!      real(8)::me,mo(3,15)
       real(8)::me
       real(8),pointer::mo(:,:)=>null()
       type(mo_s),pointer::next
     end type mo_s
     type::mo_ss
       type(mo_s),pointer::p=>null()
     end type mo_ss

     type(mo_ss)::mogood(ng,0:1)  !effectively an array of pointers to mo_s. e.g. mogood(n,0)%p%me
     type(mo_s),pointer::moi
!    type(mo_s),pointer::mogood(:,:)=>null(),moi


!    if (n_lay.gt.15) then
!      call output(' **** fatal error in find_energy_barrier_ms ****')
!      stop
!    endif

if (TRAJ.ne.0.and..not.associated(layer_order)) then; write(os,"('e vs angle for grain',i5)") TRAJ; call draw(title=trim(os),strip=.true.); endif

     if (.not.associated(layer_order)) then; allocate(layer_order(n_lay)); do i=1,n_lay;layer_order(i)=n_lay-i+1; enddo; endif

     plotting=get_plot(); call set_plot(plotting=.false.)

     nullify(tracki)
     allocate(grain_do(ng),grain_do_copy(ng),grain_do_false(ng)); grain_do=decay(1)%traj_dump; grain_do_false=.false.
     do i=1,n_lay; decay(i)%do_grain=>decay(i)%layer%do_grain;call set_data(i,grain_do_copy,do_grain_p=.true.); enddo
     err_p=.false.

     do i=1,n_lay; decay(i)%layer%alpha=0.d0; enddo

 call cpu_time(time1)
     me(:,1,1)=1.d300
     do olay1=1,n_lay  !for each spin
       olay=layer_order(olay1)


      !initialize
       do i=1,n_lay; mo(:,:,i,0,0)=decay(i)%m_min(:,:,1); mo(:,:,i,5,0)=decay(i)%m_min(:,:,2); enddo; me(:,0,0)=msdata(:,1); me(:,5,0)=msdata(:,2) 
       mo(:,:,:,1,0)=mo(:,:,:,0,0); mo(:,:,:,4,0)=mo(:,:,:,5,0); me(:,1,0)=me(:,0,0); me(:,4,0)=me(:,5,0);
!print "('olay: ',:,50(i0,x))",olay,count(grain_do),ng;read *

       it=0; flipd=.false.
!      dangle=acos(-1.d0)/NTRAJ/2.d0
       dangle=4.5d0*acos(-1.d0)/180.d0
!      dangle=acos(-1.d0)/4.d0
!      dangle=acos(-1.d0)/NTRAJ!/2.d0
       do while (any(grain_do)) !something?
         it=it+1
print "('iter: ',:,50(i0,x))",it,olay,count(grain_do),ng
!do n=1,ng;if (grain_do(n)) exit;enddo; print *,'still spinning- ',n

       !load largest energy(?) initial orientation that goes to local minimum
         do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:) = mo(:,:,i,1,0); enddo
!print *,'init'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

       !rotate spin by dangle
         do n=1,ng
           if (grain_do(n)) then
             m0=mo(:,n,olay,1,0) !decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n) !initial m0 orientation
             if (abs(dot_product(m0,mo(:,n,olay,0,0))).lt.0.999998d0) then  !current m and min energy orientation are not the same
               u=cross(mo(:,n,olay,0,0),m0)            ! vector perpendicular to m and min energy orientation
             else
               i=maxloc(abs(m0),DIM=1); j=max(1,mod(i+1,4))  ! some vector perpendicular to m
               u=0.d0; u(j)=-m0(i); u(i)=m0(j)
             endif
             mo(:,n,olay,2,0)=rotate(dangle(n),u,m0,.true.)
             decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n)=rotate(dangle(n),u,m0,.true.) !rotate!
             grain_do(n)=dot_product(mo(:,n,olay,0,0),m0).gt.dot_product(mo(:,n,olay,0,0),decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n))
!if (n.eq.TRAJ) print "(' rot: ',l,x,:,50(f12.5,x))",grain_do(n),(me(n,1,0)-me(n,0,0))/(kb*temperature),(/dangle(n),acos(min(1.d0,max(-1.d0,dot_product(m0,rotate(dangle(n),u,m0,.true.))))),acos(min(1.d0,max(-1.d0,dot_product(m0,mo(:,n,olay,0,0))))),acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,0,0),rotate(dangle(n),u,m0,.true.)))))/)*180.d0/acos(-1.d0),rotate(dangle(n),u,m0,.true.)

!if (abs(1.d0-dot_product(rotate(dangle(n),u,m0,.true.),rotate(dangle(n),u,m0,.true.))).gt.1.d-8) print *,n,rotate(dangle(n),u,m0,.true.),dot_product(rotate(dangle(n),u,m0,.true.),rotate(dangle(n),u,m0,.true.))

           endif
         enddo

         call get_field(t0); me(:,2,0)=ms_energy !get energy & field

       !run single spin for one orbit
         do i=1,n_lay; if (i.ne.olay) call set_data(i,grain_do_false,do_grain_p=.true.); enddo; grain_do_copy=grain_do !only this spin spins
!print *,'after rot'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)
!if (olay.eq.2) read *
!        t=t0; call reset_track(track,max_dti=1.d0,one_layer=olay); call set_llg(0,llg_type=-2); call llg(t,.true.,track,trace_p=.true.)
         t=t0; call reset_track(track,max_dti=1.d0,one_layer=olay); call set_llg(0,llg_type=-1); call llg(t,.true.,track,trace_p=.true.)
!         t=t0; call reset_track(track,max_dti=1.d22,min_dti=torbit,test_rot_pi=test_rot_p)
!         call set_llg(0,llg_type=-1,t_max=1.d22,del_t=torbit*1000.d0,conv_torque=0.d0,conv_angle=0.d0)
!         call llg(t,.true.,track,trace_p=.true.)

!print *,'after single spin'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!call get_field(t0,update_do_grain_p=.true.)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

       !get m such min(|mxH|) and -mxdm/dt.m0 <= 0 or NOT
         tracki=>track; do i=1,olay-1; tracki=>tracki%next; enddo
         call dump_track(tracki,decay(olay)%trj_energy(:,1),decay(olay)%trj_mxh2(:,1),decay(olay)%trj_lm(:,1),decay(olay)%trj_mstart(:,:,1),grain_do)
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
         do i=1,n_lay; call set_data(i,grain_do_copy,do_grain_p=.true.); enddo; grain_do_copy=grain_do !all multispins
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)

!print *,grain_do(TRAJ),(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print *,associated(grain_do_copy),(associated(grain_do_copy,decay(i)%layer%do_grain),i=1,n_lay)
         call get_field(t0,update_do_grain_p=.true.)
!print *,'after single spin'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)
!if (olay.eq.2) read *

!n=TRAJ;if (n.gt.0.and.grain_do(n)) print "('   ec? ',:,55(f18.13))",ms_energy(n)/me(n,2,0),(/(me(n,2,0)-ms_energy(n)),ms_energy(n),me(n,2,0),ms_energy(n)-me(n,0,0),me(n,2,0)-me(n,0,0)/)/(kb*temperature)
!print *,ms_energy(TRAJ)/(kb*temperature),' orbit'

      !now relax other spins, fixing this spin. just do damping (small steps?)
         do n=1,ng    ! try to avoid a trap where we don't increase the angle from minimum energy orientation
           if (grain_do(n)) then
             m0=decay(olay)%trj_mstart(4:6,n,1)
!if (n.eq.NTRAJ.and.dot_product(m0,mo(:,n,olay,0,0)).gt.dot_product(m0,mo(:,n,olay,1,0))) print *,'force increase',dot_product(m0,mo(:,n,olay,0,0)),dot_product(m0,mo(:,n,olay,1,0)),dot_product(mo(:,n,olay,2,0),mo(:,n,olay,1,0))
!if (n.eq.NTRAJ.and.dot_product(m0,mo(:,n,olay,0,0)).gt.dot_product(m0,mo(:,n,olay,1,0))) read *
             if (dot_product(m0,mo(:,n,olay,0,0)).gt.dot_product(mo(:,n,olay,1,0),mo(:,n,olay,0,0))) m0=mo(:,n,olay,2,0)
             decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n)=m0
             mo(:,n,olay,2,0)=m0
!if (n.eq.TRAJ) print *,(/acos(min(1.d0,max(-1.d0,dot_product(m0,mo(:,n,olay,0,0))))),acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,1,0),mo(:,n,olay,0,0)))))/)*180.d0/acos(-1.d0)
           endif
         enddo
         call get_field(t0); me(:,2,0)=ms_energy !get energy & field
!print *,'min mxH'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

!print *,ms_energy(TRAJ)/(kb*temperature),' min mxH'
!print *,'min mxH',(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
         move=.true.; move(olay,:)=.false.
         i=damp(ng,grain_do_copy,old_alpha,move,5.d0*acos(-1.d0)/180.d0,0.1d0*acos(-1.d0)/180.d0)

!print *,'relaxed other m'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

!grain_do_copy=grain_do;print *,'relaxed',(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print *,2,(decay(i)%layer%do_grain(TRAJ),i=1,n_lay),(ms_energy(n)-me(n,0,0))/(kb*temperature);read *
      !store positions
         do i=1,n_lay; mo(:,:,i,3,0)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:); enddo  !equilibrium w/ olay frazen
         me(:,3,0)=ms_energy

!print *,(/ms_energy(TRAJ),me(TRAJ,0,0)/)/(kb*temperature),' relaxed'

!n=TRAJ;if (n.gt.0.and.grain_do(n))print "(' new? ',:,50(f12.5,x))",(me(n,3,0)-me(n,0,0))/(kb*temperature),(/acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,0,0),mo(:,n,olay,3,0))))),dangle(n),acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,1,0),mo(:,n,olay,3,0)))))/)*180.d0/acos(-1.d0)

!print *,3,(ms_energy(n)-me(n,0,0))/(kb*temperature),mo(:,TRAJ,olay,3,0);read *

      !and now let them all go!
         grain_do_copy=grain_do; move=.true.
         i=damp(ng,grain_do_copy,old_alpha,move,5.d0*acos(-1.d0)/180.d0,0.1d0*acos(-1.d0)/180.d0)
!print *,'freed'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

!print *,'freed'
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!n=TRAJ;if (n.gt.0.and.grain_do(n)) then
!write(os,"(:,500l)") mo(3,n,:,0,0).gt.0.d0; os=trim(os)//'->';do i=1,n_lay; write(os,"(a,l)") trim(os),decay(i)%layer%m(3+decay(i)%layer%in,n).gt.0.d0; enddo
!print "(' fnl: ',a,:,x,50(f12.5,x))",trim(os),(ms_energy(n)-me(n,0,0))/(kb*temperature),(/acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,0,0),decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n))))),acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,1,0),decay(olay)%layer%m(1+decay(olay)%layer%in:3+decay(olay)%layer%in,n)))))/)*180.d0/acos(-1.d0)
!endif

!do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,:)=mo(:,:,i,1,0); enddo
!call get_field(t0)
!print *,'reset!'
!print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay); grain_do_copy=grain_do; print "(:,50l)",(decay(i)%layer%do_grain(TRAJ),i=1,n_lay)
!print "(:,3(x,f12.5))",(/ms_energy(TRAJ),me(TRAJ,0,0),me(TRAJ,1,0)/)/(kb*temperature)
!print "(3(x,f12.5))",(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ),i=1,n_lay)
!print "(3(x,f12.5))",(cross(decay(i)%layer%h(:,TRAJ),decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,TRAJ)),i=1,n_lay)

!read *

      !where did they go?
         do n=1,ng
           if (grain_do(n)) then
             do i=1,n_lay  !is any spin more than 60 degrees from local energy minimum orientation?
               broke_grain(i)=(dot_product(decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n),mo(:,n,i,0,0)).lt..5d0)
             enddo
             if (any(broke_grain)) then !at least one spin is 60 degrees or more away from the local minimum orientation
         ! we didn't return to local energy orientation
!              grain_do(n)=(dangle(n).gt.0.1d0*acos(-1.d0)/180.d0)
!              grain_do(n)=(dangle(n).gt.0.3d0*acos(-1.d0)/180.d0)
               grain_do(n)=(dangle(n).gt.0.29d0*acos(-1.d0)/180.d0)
               dangle(n)=dangle(n)*3.d0/10.d0 !/3.d0 !2.5d0 !/2.d0
               do i=1,n_lay; mo(:,n,i,4,0)=decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n); enddo; me(n,4,0)=ms_energy(n)
if (TRAJ.eq.n) then
write(os,"('broke/flip grain: ',i0,x,l,x,:,500l)")count(broke_grain),grain_do(n),broke_grain
write(os,"(a,' ',:,500l)") trim(os),mo(3,n,:,0,0).gt.0.d0; write(os,"(a,' -> ',:,500l)") trim(os),mo(3,n,:,4,0).gt.0.d0;
print *,trim(os)
print "(l,:,33(x,f12.5))",grain_do(n),(me(n,:,0)-me(n,0,0))/(kb*temperature)
!read *
endif

            !lets' back up two if we can...
               moi=>mogood(n,0)%p; if (associated(moi)) moi=>moi%next; if (.not.flipd(n).and.associated(moi)) moi=>moi%next
               if (associated(moi)) then
                 dangle(n)=dangle(n)*10.d0/3.d0/4.5d0 !/1.6d0 !/2.d0
                 deallocate(mogood(n,0)%p%mo); deallocate(mogood(n,0)%p)
                 mogood(n,0)%p=>moi
                 mo(:,n,:,1,0)=moi%mo  !equilibrium w/ olay frazen
                 me(n,1,0)=moi%me
               endif

               flipd(n)=.true.
             else
         ! we fell back to local energy orientation
               grain_do(n)=((me(n,3,0)-me(n,0,0).lt.max_eb*kb*temperature).and.(me(n,3,0).lt.me(n,1,1)))  !Eb must be less than max_eb and any previous Eb
if (n.eq.TRAJ) print *,'did not break',grain_do(n),acos(min(1.d0,max(-1.d0,(/dot_product(mo(:,n,olay,3,0),mo(:,n,olay,0,0)),dot_product(mo(:,n,olay,1,0),mo(:,n,olay,0,0))/))))*180.d0/acos(-1.d0)
if (n.eq.TRAJ.and.dot_product(mo(:,n,olay,3,0),mo(:,n,olay,0,0)).gt.dot_product(mo(:,n,olay,1,0),mo(:,n,olay,0,0))) read *
               if (me(n,3,0).lt.me(n,1,0)) then
!print "('damn energy:',f12.5,' ->',f12.5,x,f12.5)",(me(n,1,0)-me(n,0,0))/(kb*temperature),(me(n,3,0)-me(n,0,0))/(kb*temperature),(me(n,1,0)-me(n,3,0))/(kb*temperature)
!write(os,"('damn energy: ',i0,x,f12.5,' ->',f12.5,x,f12.5)")n,(me(n,1,0)-me(n,0,0))/(kb*temperature),(me(n,3,0)-me(n,0,0))/(kb*temperature),(me(n,1,0)-me(n,3,0))/(kb*temperature);call output(os)
if (n.eq.TRAJ)then;write(os,"('damn energy: ',i0,x,f12.5,' ->',f12.5,x,f12.5)")n,(me(n,1,0)-me(n,0,0))/(kb*temperature),(me(n,3,0)-me(n,0,0))/(kb*temperature),(me(n,1,0)-me(n,3,0))/(kb*temperature);call output(os);endif
if (n.eq.TRAJ) print "(l,:,33(x,f12.5))",grain_do(n),(me(n,:,0)-me(n,0,0))/(kb*temperature)
if (n.eq.TRAJ) read *
               endif
          !new good starting position
               mo(:,n,:,1,0)=mo(:,n,:,3,0)  !equilibrium w/ olay frazen
               me(n,1,0)=me(n,3,0)

               if (.not.flipd(n).and.dot_product(mo(:,n,olay,3,0),mo(:,n,olay,0,0)).lt.-0.866d0) grain_do(n)=.false.
               grain_do(n)=grain_do(n).and.(flipd(n).or.dot_product(mo(:,n,olay,3,0),mo(:,n,olay,0,0)).gt.-0.866d0)

               allocate(moi); allocate(moi%mo(3,n_lay)); moi%me=me(n,3,0); moi%mo=mo(:,n,:,3,0); moi%next=>mogood(n,0)%p; mogood(n,0)%p=>moi

if (TRAJ.eq.n) then; write(os,"('e vs angle for grain',i5)") n
if (it.eq.1) call draw(title=trim(os),strip=.true.,add=.false.)
call draw(title=trim(os),strip=.true.,add=.true.,x=(/ acos(min(1.d0,max(-1.d0,dot_product(mo(:,n,olay,1,0),mo(:,n,olay,0,0))))) /)*180.0/acos(-1.d0),y=(/me(n,3,0)-me(n,0,0)/)/(kb*temperature))
endif

             endif
           endif
         enddo

       enddo  !iteration
!print *,4;read *

       grain_do=decay(1)%traj_dump
       do n=1,ng
         !did we find an alternative minimum and the barrier is less than the current minimum value?
         if (grain_do(n).and.flipd(n).and.me(n,1,0)-me(n,0,0).lt.min(max_eb*kb*temperature,me(n,1,1)-me(n,0,0))) then 
!print *,'oh my!',n,olay,(me(n,1,0)-me(n,0,0))/(kb*temperature)
           im(n)=olay; mo(:,n,:,:,1)=mo(:,n,:,:,0); me(n,:,1)=me(n,:,0);
           do while (associated(mogood(n,1)%p)); moi=>mogood(n,1)%p; mogood(n,1)%p=>moi%next; deallocate(moi%mo); deallocate(moi); enddo; mogood(n,1)%p=>mogood(n,0)%p
         endif
         nullify(mogood(n,0)%p)
       enddo
     enddo  !layers
!print *,4;read *

  !in range of interest?
     decay(1)%traj_dump=grain_do.and.me(:,1,1)-me(:,0,1).lt.max_eb*kb*temperature
     grain_do=decay(1)%traj_dump

if (TRAJ.ne.0) then
write(os,"('msenergy barrier, e, final')")
call plot_histogram(trim(os),20,pack(me(:,1,1)-me(:,0,1),mask=grain_do)/(kb*temperature))
endif

     err_p=(count(grain_do).lt.1)


call cpu_time(time2);time1=time2-time1
     if (.not.err_p) then
       write(os,"('finding transport coefficients for ',i4,' out of ',i4,' multispins')")count(grain_do),ng; call output(os)
n=TRAJ;if (TRAJ.ne.0) print "(' fixed layer ',i0,' for spin ',i0,x,l)",im(n),n,grain_do(n)
     do i=1,n_lay;layer_order(i)=count(grain_do.and.im.eq.i);write(os,"(i0,' fixed in layer',:,50(x,i0))")layer_order(i),i; call output(trim(os)); enddo
     n=0; do i=1,n_lay;n=n*10+maxloc(layer_order,dim=1); layer_order(maxloc(layer_order,dim=1))=-1; enddo
     do i=n_lay,1,-1; layer_order(i)=n-(n/10)*10;n=n/10; enddo
print *,layer_order

       if (diff_m) then
         do i=1,n_lay; decay(i)%layer%background_t=temperature; decay(i)%layer%thermal_p=.true.; enddo  !turn on thermal kisk along H
       endif

  ! we know which spin to advance (im(n) stored in mo(:,n,:,1,1)) and approximately the energy barrier (me(n,1,1))
  !  let's just run on the orbits we found


  !initialize
       it=0; decay(1)%traj_lch=0  !iteration and running number of orbits for each multispin
       do i=1,n_lay; call set_data(i,grain_do_copy,do_grain_p=.true.); enddo !all multispins
       do n=1,ng; grain_do(n)=grain_do(n).and.associated(mogood(n,1)%p); enddo !only do spins that have orientations
       decay(1)%traj_dump=grain_do  !final list of multispins we are following 

  !do all or NTRAJ trajectories for each valid spin
       do while (any(grain_do))
         it=it+1


    ! load up a good m orientation
         do n=1,ng
           if (grain_do(n)) then
             moi=>mogood(n,1)%p%next; n1=decay(1)%traj_lch(n); x=decay(1)%bltz(1,1,n)+log(1.d0-n1*(1.d0-exp(-5.d0))/(NTRAJ-1))
             do while (n1.gt.0.and.associated(moi))  !try to sample in energy via exponential spacing
               if (moi%me-msdata(n,1).lt.x*(kb*temperature).and.mogood(n,1)%p%me-msdata(n,1).lt.(kb*temperature)*decay(1)%bltz(1,1,n)) exit
!print "('skipping ',i0,x,i0,:,50(x,f12.5))",n,n1,x,(/mogood(n,1)%p%me-msdata(n,1),moi%me-msdata(n,1)/)/(kb*temperature)
               moi=>mogood(n,1)%p; mogood(n,1)%p=>moi%next; deallocate(moi%mo); deallocate(moi); moi=>mogood(n,1)%p%next
             enddo
             do i=1,n_lay; decay(i)%layer%m(1+decay(i)%layer%in:3+decay(i)%layer%in,n)=mogood(n,1)%p%mo(:,i); enddo; ms_energy(n)=mogood(n,1)%p%me
           endif
         enddo

print "('iter: ',:,3(i0,x),2l)",it,count(grain_do),ng,diff_m,test_rot_p

    ! track |mxH|^2 in dampless llg orbit
         t=t0; grain_do_copy=grain_do; call reset_track(track,max_dti=1.d22,min_dti=torbit,test_rot_pi=test_rot_p)
!if (TRAJ.ne.0) then; n=TRAJ
!if (associated(mogood(n,1)%p)) then; call get_field(t0); print *,(/mogood(n,1)%p%me,ms_energy(n),mogood(n,1)%p%me-ms_energy(n)/)/(kb*temperature); endif
!endif
!        call set_llg(0,llg_type=-2,t_max=1.d22,del_t=torbit*1000.d0,conv_torque=0.d0,conv_angle=0.d0)
         call set_llg(0,llg_type=-1,t_max=1.d22,del_t=torbit*1000.d0,conv_torque=0.d0,conv_angle=0.d0)
         call llg(t,.true.,track,trace_p=.true.,diff_m=diff_m);   !diffuse m along mxH

    !grab |mxH|^2 and 
         tracki=>track
         do i=1,n_lay
           call dump_track(tracki,decay(i)%trj_energy(:,1),decay(i)%trj_mxh2(:,1),decay(i)%trj_lm(:,1),decay(i)%trj_mstart(:,:,1),grain_do)
         enddo

    !grab |mxH|^2 and 
         do n=1,ng
           if (grain_do(n)) then
             moi=>mogood(n,1)%p; mogood(n,1)%p=>moi%next; grain_do(n)=associated(mogood(n,1)%p); n1=decay(1)%traj_lch(n)+1
!if (n.eq.TRAJ) print *,(/mogood(n,1)%p%me,moi%me,moi%me-msdata(n,1),msdata(n,1)/)/(kb*temperature),decay(1)%bltz(1,1,n)
             if (.not.test_rot_p.or.decay(1)%trj_lm(n,1).gt.0.d0) then
!            if ((n1.lt.2.or.moi%me-msdata(n,1).le.decay(1)%bltz(1,1,n)*(kb*temperature)).and.(.not.test_rot_p.or.decay(1)%trj_lm(n,1).gt.0.d0)) then
               decay(1)%traj_lch(n)=n1
               grain_do(n)=(grain_do(n).and.n1.lt.NTRAJ) !don't need more than NTRAJ orbits...
!if (n.eq.TRAJ) read *
!              decay(1)%bltz(n1:NTRAJ,1,n)=(ms_energy(n)-msdata(n,1))/(kb*temperature) !decay(1)%bltz(1:nint(),1,1:n) = de/kT
               decay(1)%bltz(n1:NTRAJ,1,n)=(moi%me-msdata(n,1))/(kb*temperature) !decay(1)%bltz(1:nint(),1,1:n) = de/kT
               decay(1)%bltz(n1,2,n)=0.d0
               do i=1,n_lay
                 decay(1)%bltz(n1,2,n)=decay(1)%bltz(n1,2,n)+decay(i)%layer%gamma*old_alpha(i)* &
                                decay(i)%layer%ms(n)*decay(i)%layer%ms_scale(n)* &
                                decay(i)%layer%grain_vol(n)*decay(i)%trj_mxh2(n,1)/(1.d0+old_alpha(i)**2)
               enddo
           grain_do(n)=(grain_do(n).and.maxval(decay(1)%bltz(1:n1,1,n))-decay(1)%bltz(n1,1,n).lt.5.d0)
!if (n.eq.TRAJ) print "(i0,x,l,:,55(x,f18.13))",n1,grain_do(n),(/moi%me,moi%me-msdata(n,1)/)/(kb*temperature),decay(1)%bltz(1:n1,1,n),decay(1)%bltz(n1,2,n)
if (n.eq.TRAJ) then; write(os,"('a, ',i4)") n; call draw(trim(os),x=decay(1)%bltz(1:n1,1,n),y=decay(1)%bltz(1:n1,2,n)); endif
             endif
             deallocate(moi%mo); deallocate(moi)
           endif
         enddo

       enddo

       grain_do=decay(1)%traj_dump
       decay(1)%traj_dump=(grain_do.and.decay(1)%traj_lch.gt.1)
       write(os,"('found transport coefficients for ',i0,' out of ',i0,' with ',i0,' multispins')")count(decay(1)%traj_dump),count(grain_do),ng
       call output(os)
       grain_do=decay(1)%traj_dump

        !sort via energy
       do n=1,ng
         if (decay(1)%traj_dump(n)) then; x=decay(1)%bltz(1,1,n)
           decay(1)%bltz(1,3,n)=decay(1)%bltz(1,1,n)
           decay(1)%bltz(1,4,n)=decay(1)%bltz(1,2,n)
           do n2=2,decay(1)%traj_lch(n)
             do n3=1,n2-1; if (decay(1)%bltz(n2,1,n).lt.decay(1)%bltz(n3,3,n)) exit; enddo
             do n1=n2-1,n3,-1
               decay(1)%bltz(n1+1,3,n)=decay(1)%bltz(n1,3,n)
               decay(1)%bltz(n1+1,4,n)=decay(1)%bltz(n1,4,n)
             enddo
             decay(1)%bltz(n3,3,n)=decay(1)%bltz(n2,1,n)
             decay(1)%bltz(n3,4,n)=decay(1)%bltz(n2,2,n)
           enddo
!          do n2=1,decay(1)%traj_lch(n); if (decay(1)%bltz(n3,3,n).gt.x) exit;enddo; decay(1)%traj_lch(n)=min(decay(1)%traj_lch(n),n2)
!      write(os,"('found transport coefficients for ',i0,' out of ',i0,' with ',i0,' multispins')")count(decay(1)%traj_dump),count(grain_do),ng
!      call output(os)

!        elseif (iter.ge.0.and.decay(1)%bltz(decay(1)%traj_lch(n),1,n).lt.max_eb) then  !grain w/ barrier less than kT! don't pollute the gammas
!          decay(1)%bltz(:,3,n)=decay(1)%bltz(:,1,n); dirs=''
!          do n3=1,n_lay;write(dirs,"(a,l)") trim(dirs),decay(n3)%m(3,n).gt.0.d0; enddo
!          if (random_uniform(random_decay).gt.1.d0/(1.d0+exp((msdata(n,1)-msdata(n,2))/(kb*temperature)))) then
!            decay(1)%num_flip(n)=decay(1)%num_flip(n)+1
!            dirs=''; do n3=1,n_lay;decay(n3)%m(:,n)=decay(n3)%m_min(:,n,4); write(dirs,"(a,l)") trim(dirs),decay(n3)%m(3,n).gt.0.d0; enddo
!            write(os,"(i0,a,1p,5(x,e12.5),x,a)")n,' flopped',decay(1)%bltz(decay(1)%traj_lch(n),1,n),msdata(n,1)/(kb*temperature),msdata(n,2)/(kb*temperature), &
!                    (msdata(n,1)-msdata(n,2))/(kb*temperature),1.d0/(1.d0+exp((msdata(n,1)-msdata(n,2))/(kb*temperature))),trim(dirs)
!          else
!            write(os,"(i0,a,1p,5(x,e12.5),x,a)")n,' did not flopped',decay(1)%bltz(decay(1)%traj_lch(n),1,n),msdata(n,1)/(kb*temperature),msdata(n,2)/(kb*temperature), &
!                    (msdata(n,1)-msdata(n,2))/(kb*temperature),1.d0/(1.d0+exp((msdata(n,1)-msdata(n,2))/(kb*temperature))),trim(dirs)
!          endif
!          call output(trim(os))
         endif
       enddo
       do n=1,ng; decay(1)%bltz(NTRAJ,3,n)=decay(1)%bltz(decay(1)%traj_lch(n),3,n); enddo
       decay(1)%bltz(:,1,:)=decay(1)%bltz(:,3,:)
       decay(1)%bltz(:,2,:)=decay(1)%bltz(:,4,:)

       do n=1,ng
         if (grain_do(n)) then
           msdata(n,2)=me(n,4,1)
           do i=1,n_lay; decay(i)%m_min(:,n,4)=mo(:,n,i,4,1);enddo
!print "(i4,1p,:,40(x,e12.5))",n,decay(1)%bltz(NTRAJ,1,n),msdata(n,1)/(kb*temperature),msdata(n,2)/(kb*temperature)
         endif
       enddo

       write(os,"(' final number of gammas to be computed: ',i0)") count(decay(1)%traj_dump); call output(trim(os))

!if (get_plot()) call plot_contour('grain energy barrier (kT)b',decay(1)%bltz(NTRAJ,1,:),.true.,1,n0=5)
!if (get_plot()) call plot_contour('grain energy barrier (kT)',decay(1)%bltz(NTRAJ,1,:),.true.,1)
!if (get_plot()) call plot_histogram('energy barrier (kT)',20,pack(decay(1)%bltz(NTRAJ,1,:),mask=decay(1)%bltz(NTRAJ,1,:).lt.max_eb))
!if (get_plot()) call plot_histogram('small energy barrier (kT)',20,pack(decay(1)%bltz(NTRAJ,1,:),mask=decay(1)%bltz(NTRAJ,1,:).lt.10.d0))
!if (get_plot()) call plot_histogram('energy barrier (kT)',20,pack(decay(1)%bltz(NTRAJ,1,:),mask=decay(1)%bltz(NTRAJ,1,:).lt.max_eb),xmin=20.d0,xmax=40.d0)

       err_p=.false.
       call cpu_time(time3);time2=time3-time2
!      print "('timing of search: ',f8.3,', timing of coefficients: ',3(f8.3,x))",time1,time2,time1*1.d2/(time2+time1),time2*1.d2/(time2+time1)
       write(os,"('timing of search: ',f8.3,', timing of coefficients: ',3(f8.3,x))")time1,time2,time1*1.d2/(time2+time1),time2*1.d2/(time2+time1); call output(os)
     else
       call output('no transport coefficients to find!')
!      print "('timing of search: ',f8.3)",time1
       write(os,"('timing of search: ',f8.3)")time1; call output(os)
     endif

     do i=1,n_lay; call set_data(i,decay(i)%do_grain,do_grain_p=.true.); enddo
     deallocate(grain_do,grain_do_copy,grain_do_false); nullify(grain_do,grain_do_copy,grain_do_false)
!read *

     call set_plot(plotting=plotting)
    end function


    function find_energy_barrier(iter) result (err_p)
     use track_m, only: reset_track, dump_track
     integer,intent(in)::iter
     logical::err_p
     type(track_s),pointer::tracki
     integer::itersaddle,n,il,i,n1,n3,n4,n2
     real(8)::x


     flipped=0
     if (multispin_p) then; err_p=find_energy_barrier_ms(size(decay(1)%layer%grain_vol),iter); return; endif

     do i=1,n_lay; decay(i)%do_grain=>decay(i)%layer%do_grain;call set_data(i,decay(i)%traj_dump,do_grain_p=.true.); enddo
     do itersaddle=nitersaddle,0,-1
       do n=1,NTRAJ
         do il=1,n_lay
           decay(il)%traj_lch=n
           call fan_trajectory(decay(il)%m_min(:,:,1),decay(il)%m_min(:,:,2),decay(il)%layer%m,decay(il)%layer%in,decay(il)%traj_lch, &
                   NTRAJ+2,il,decay(il)%trj_m0(:,:,n))
         enddo
         t=t0; call reset_track(track); 
!  read *
         call llg(t,.true.,track)
!        call llg(t,.false.,track)
         tracki=>track
         do il=1,n_lay
           call dump_track(tracki,decay(il)%trj_energy(:,n),decay(il)%trj_mxh2(:,n),decay(il)%trj_lm(:,n),decay(il)%trj_mstart(:,:,n),decay(il)%traj_dump)
!          tracki=>tracki%next
         enddo
       enddo

       i=0
       do il=1,n_lay
         do n=1,size(decay(il)%layer%grain_vol)
           decay(il)%bltz_m0(:,1,0,n)=decay(il)%m_min(:,n,1)
           decay(il)%bltz_m0(:,2,0,n)=decay(il)%m_min(:,n,1)
           decay(il)%bltz_m0(:,1,NTRAJ+1,n)=decay(il)%m_min(:,n,2)
           decay(il)%bltz_m0(:,2,NTRAJ+1,n)=decay(il)%m_min(:,n,2)

           if (decay(il)%traj_dump(n)) then
             do n1=1,NTRAJ
               if (decay(il)%trj_lm(n,n1).le.0.d0) exit !orbit escaped!
             enddo
             n3=min(n1,NTRAJ)
             x=decay(il)%layer%grain_vol(n)/(kb*temperature)
             if (itersaddle.eq.nitersaddle.and.n1.le.NTRAJ)then
               if ((decay(il)%trj_energy(n,max(1,n1-1))-decay(il)%energy(1,n))*x.lt.max_eb)flipped=flipped+1 !FIXME eb/kT < 50
             endif
             if (n.eq.TRAJ) print *,itersaddle,n1,TRAJ,decay(il)%energy(1,n)*x,decay(il)%trj_energy(n,1:n3)*x,decay(il)%trj_lm(n,1:n3)
             if (n.eq.TRAJ) print *,'n1.gt.NTRAJ:',(n1.gt.NTRAJ)
             if (n.eq.TRAJ) print *,'old m_min:',decay(il)%m_min(:,n,1),decay(il)%m_min(:,n,2)
             if (n1.ne.n3) then
               do n3=size(decay(il)%bltz,1)-1,1,-1; if (decay(il)%bltz(n3,1,n).lt.1.d200) exit;enddo
               if (n3.gt.0) then
                 decay(il)%bltz(n1,1:2,n)=decay(il)%bltz(n3,1:2,n)
                 decay(il)%bltz_m0(:,1:2,n1,n)=decay(il)%bltz_m0(:,1:2,n3,n)
                 decay(il)%energy(3,n)=decay(il)%trj_energy(n,n1-1) !- best guess for saddle energy
                 decay(il)%energy(4,n)=decay(il)%trj_mxh2(n,n1-1)
               else
                 decay(il)%bltz(n1,1,n)=decay(il)%energy(2,n)
                 decay(il)%bltz(n1,2,n)=0.d0
                 decay(il)%bltz_m0(:,2,n1,n)=decay(il)%m_min(:,n,2)
                 decay(il)%bltz_m0(:,1,n1,n)=decay(il)%m_min(:,n,2)
               endif
               n3=min(n1,NTRAJ)
             endif
             decay(il)%bltz(1:n3,1,n)=decay(il)%trj_energy(n,1:n3)      !energy of orbit
             decay(il)%bltz(1:n3,2,n)=decay(il)%trj_mxh2(n,1:n3)        !time averaged mxh2
             decay(il)%bltz_m0(:,2,1:n3,n)=decay(il)%trj_mstart(1:3,n,1:n3) !m for max mxH on orbit (?)
             decay(il)%bltz_m0(:,1,1:n3,n)=decay(il)%trj_mstart(4:6,n,1:n3) !m for min mxH on orbit (?)
             if (n1.eq.n3.and.n1.gt.1) then !n1 holds the first escaped trajectory
               decay(il)%energy(3,n)=decay(il)%bltz(n1-1,1,n) !decay%trj_energy(n,n1-1) !- best !guess for !saddle !energy
               decay(il)%energy(4,n)=decay(il)%bltz(n1-1,2,n) !decay%trj_mxh2(n,n1-1)
             endif
             decay(il)%bltz(n1+1:,:,n)=1.d220

           !FIXME- do these 2 lines do anything
             decay(il)%m_min(:,n,5)=decay(il)%trj_mstart(1:3,n,max(1,n1-1))       !point of maximum |mxH| on orbit 
             decay(il)%m_min(:,n,6)=decay(il)%trj_mstart(4:6,n,max(1,n1-1))       !point of minimum |mxH| on orbit 

             decay(il)%m_min(:,n,1)=decay(il)%bltz_m0(:,1,max(0,n1-2),n) !something close to good trajectory
             decay(il)%m_min(:,n,2)=decay(il)%bltz_m0(:,1,n1,n) !first bad trajectory
             if (n.eq.TRAJ) print *,'new m_min:',decay(il)%m_min(:,n,1),decay(il)%m_min(:,n,2)
           endif
         enddo
         i=i+count(decay(il)%traj_dump)
       enddo
       write(os,"('(',i2,'/',i2,') number of saddle points:',i4,', number flipped:',i4)") nitersaddle-itersaddle+1,nitersaddle+1,i,flipped; call output(os)
!  read *
       if (flipped.eq.0) exit
     enddo
     err_p=(flipped.eq.0)!; if (err_p) return
     if (.not.err_p) then

     n4=0
     do il=1,n_lay
       decay(il)%traj_dump=.false.
       do n=1,size(decay(il)%layer%grain_vol)
         if (acos(max(-1.d0,min(1.d0,dot_product(decay(il)%m_min(:,n,3),decay(il)%m_min(:,n,4))))).gt.0.01d0*acos(-1.d0)/180.d0) then
           if (decay(il)%energy(3,n).gt.1.d200) print "('>>> ',i4,' two minimum didn''t find saddle!')",n
           if (decay(il)%energy(3,n).lt.1.d200) then; i=i+1; decay(il)%traj_dump(n)=.true.; endif
         else
           if (decay(il)%energy(3,n).lt.1.d200) print "('>>> ',i4,' found saddle w/ only one min!')",n
         endif
         x=decay(il)%layer%grain_vol(n)/(kb*temperature)
         if (TRAJ.gt.0) write(89,"(i4,1p,22(x,e17.10))")n,decay(il)%energy(:,n),(decay(il)%energy(3,n)-decay(il)%energy(1,n))*x, &
                   (decay(il)%energy(3,n)-decay(il)%energy(2,n))*x
       enddo
       do n=1,size(decay(il)%layer%grain_vol)
         if (decay(il)%traj_dump(n)) then
           do n3=size(decay(il)%bltz,1)-1,1,-1; if (decay(il)%bltz(n3,1,n).lt.1.d200) exit;enddo
           if (n3.gt.0) then
             decay(il)%m_min(:,n,1)=decay(il)%m_min(:,n,3)               !local minimum
             decay(il)%m_min(:,n,2)=decay(il)%bltz_m0(:,2,n3-1,n) !last good orbit
             decay(il)%traj_dump(n)=(decay(il)%energy(3,n).lt.1.d200)
           else
             decay(il)%traj_dump(n)=.false.
           endif
         endif
       enddo
       n4=n4+count(decay(il)%traj_dump)
     enddo

     write(os,"('number of saddle points:',i4,', number flipped:',i4)") i,flipped; call output(os)
     write(os,"('now find transport coefficients near ',i4,' critical orbits...')")n4; call output(os)

     do n2=1,nitersaddle*NTRAJ
       do il=1,n_lay
         decay(il)%traj_lch=n2
         call fan_trajectory(decay(il)%m_min(:,:,1),decay(il)%m_min(:,:,2), &
             decay(il)%layer%m,decay(il)%layer%in,decay(il)%traj_lch,nitersaddle*NTRAJ+1,il,onesided_p=.true.)
       enddo
       t=t0; call reset_track(track,.true.); call llg(t,.true.,track)
       tracki=>track
       do il=1,n_lay
         call dump_track(tracki,decay(il)%trj_energy(:,1),decay(il)%trj_mxh2(:,1),decay(il)%trj_lm(:,1),decay(il)%trj_mstart(:,:,1),decay(il)%traj_dump, &
            1.d0*old_alpha(il)*kb*temperature*decay(il)%layer%gamma/(1.d0+old_alpha(il)*old_alpha(il)))
!           1.d0*old_alpha(il)*kb*temperature*               1.76d7/(1.d0+old_alpha(il)*old_alpha(il)))
         do n=1,size(decay(il)%layer%grain_vol)
           decay(il)%gamma(n)=0.d0; decay(il)%gamma1(n)=0.d0
           if (decay(il)%traj_dump(n)) then
             decay(il)%bltz(n2,3,n)=decay(il)%trj_energy(n,1)
             decay(il)%bltz(n2,4,n)=decay(il)%trj_mxh2(n,1)
             x=decay(il)%layer%grain_vol(n)/(kb*temperature)
!            decay(il)%gamma(n)=1.76d7*old_alpha(il)*decay(il)%layer%ms(n)*decay(il)%layer%ms_scale(n)*decay(il)%bltz(n2,4,n)/(1.d0+old_alpha(il)**2)* &
             decay(il)%gamma(n)=decay(il)%layer%gamma*old_alpha(il)*decay(il)%layer%ms(n)*decay(il)%layer%ms_scale(n)*decay(il)%bltz(n2,4,n)/(1.d0+old_alpha(il)**2)* &
                   x*exp((decay(il)%energy(1,n)-decay(il)%energy(3,n))*x)
           endif
         enddo
       enddo
     enddo
     endif

     do i=1,n_lay; call set_data(i,decay(i)%do_grain,do_grain_p=.true.); enddo
    end function

    function find_rates() result(total_rate)
     use graphix, only:draw
     use plot_m, only: plot_contour, get_plot, plot_histogram, plot_scatter
      real(8)::total_rate,es,x,y,old_gamma,inty(NTRAJ,2)
      integer::il,n,n3,n2,n1!,n0
      logical::err_p

      total_rate=0.d0
      if (multispin_p) then
    print *,'in find_rates',count(decay(1)%traj_dump),NTRAJ
        x=kb*temperature
        decay(1)%gamma=0.d0; decay(1)%gamma1=0.d0
        do n=1,size(decay(1)%layer%grain_vol)
          if (decay(1)%traj_dump(n)) then
            decay(1)%bltz(1:,1,n)=max(0.d0,decay(1)%bltz(1:,1,n)-mkt)*x
            do n3=1,decay(1)%traj_lch(n) !NTRAJ !nint(decay(1)%bltz(0,1,1)) !nitersaddle*NTRAJ
              if (decay(1)%bltz(n3,1,n).ge.0.d0) exit
            enddo
!FIXME hard limit on when to reject computed a and replace by average with previous a
!y=maxval(decay(1)%bltz(n3:NTRAJ,2,n)); decay(1)%bltz(n3:NTRAJ,2,n)=max(0.01d0*y,decay(1)%bltz(n3:NTRAJ,2,n))
!            decay(1)%bltz(NTRAJ,2,n)=decay(1)%bltz(NTRAJ-1,2,n)
! decay(1)%bltz(n3:NTRAJ,2,n)=1.d-3
!           do n2=n3+1,NTRAJ
!             if (decay(1)%bltz(n2,2,n)/decay(1)%bltz(n2-1,2,n).lt.10.d0) decay(1)%bltz(n2,2,n)=(decay(1)%bltz(n2,2,n)+decay(1)%bltz(n2-1,2,n))/2.d0
!             if (decay(1)%bltz(n2,2,n)/decay(1)%bltz(n2-1,2,n).lt.20.d0) decay(1)%bltz(n2,2,n)=decay(1)%bltz(n2-1,2,n)
!           enddo
  

            decay(1)%traj_dump(n)=(n3.lt.decay(1)%traj_lch(n))
            if (decay(1)%traj_dump(n)) then !NTRAJ) then
              inty=0.d0

if (n.eq.TRAJ) then
 do n2=n3,decay(1)%traj_lch(n) !NTRAJ !nint(decay(1)%bltz(0,1,1)) !nitersaddle*NTRAJ
   print "(i5,i3,1p,2(x,e12.5))",n,n2,decay(1)%bltz(n2,1,n)/x,decay(1)%bltz(n2,2,n)
 enddo
 print *,mkt,NTRAJ,n3 !nint(decay(1)%bltz(0,1,1)),n3
! read *
endif
              y=(exp(decay(1)%bltz(n3,1,n)/x)-1.d0)*0.5d0*(decay(1)%bltz(n3+1,1,n)-decay(1)%bltz(n3,1,n))/decay(1)%bltz(n3,2,n)
  if (y.lt.0.d0) print *,n,n3,y
              inty(n3,2)=y; inty(n3,1)=decay(1)%bltz(n3,3,n)
if (n.eq.TRAJ) print "(i5,i3,1p,:,4(x,e12.5))",n,n3,decay(1)%bltz(n3,1:2,n)/x,y,1.d0/y
              do n2=n3+1,decay(1)%traj_lch(n)-1 !NTRAJ-1 !nint(decay(1)%bltz(0,1,1))-1
                y=y+(exp(decay(1)%bltz(n2,1,n)/x)-1.d0)*0.5d0*(decay(1)%bltz(n2+1,1,n)-decay(1)%bltz(n2-1,1,n))/decay(1)%bltz(n2,2,n)
                inty(n2,2)=y; inty(n2,1)=decay(1)%bltz(n2,3,n)
  if (y.lt.0.d0) print *,n,n2,y
if (n.eq.TRAJ) print "(i5,i3,1p,:,4(x,e12.5))",n,n2,decay(1)%bltz(n2,1:2,n)/x,y,1.d0/y
              enddo
              n2=decay(1)%traj_lch(n) !NTRAJ !nint(decay(1)%bltz(0,1,1)) !nitersaddle*NTRAJ
              y=y+(exp(decay(1)%bltz(n2,1,n)/x)-1.d0)*0.5d0*(decay(1)%bltz(n2,1,n)-decay(1)%bltz(n2-1,1,n))/decay(1)%bltz(n2,2,n)
              inty(n2,2)=y; inty(n2,1)=decay(1)%bltz(n2,3,n)
  if (y.lt.0.d0) print *,n,n2,y
if (n.eq.TRAJ) print "(i5,i3,1p,:,4(x,e12.5))",n,n2,decay(1)%bltz(n2,1:2,n)/x,y,1.d0/y
              if (y.gt.0.d0) decay(1)%gamma1(n)=1.d0/y
              decay(1)%gamma(n)=decay(1)%gamma1(n)
if (n.eq.TRAJ) then
do n1=n2-1,n3+1,-1; if (inty(n2,1)-inty(n1,1).gt.3.d0) exit;enddo
write(os,"('1/gamma(e) vs e, for grain ',i0)") n !; call draw(title=trim(os),strip=.true.); call draw(title=trim(os),strip=.true.,add=.false.)
call draw(title=trim(os),x=inty(n1:n2,1),y=inty(n1:n2,2))
call draw(title=trim(os),x=inty(n1:n2,1),y=inty(n1:n2,2),add=.true.,symbol=.true.,nps=5,color=1)
call draw(title=trim(os),x=inty(n1:n2,1),y=inty(n2,2)*exp(-inty(n2,1)+inty(n1:n2,1)),add=.true.,color=3,unsort=.true.)
write(os,"('log10(1/gamma(e)) vs e, for grain ',i0)") n !; call draw(title=trim(os),strip=.true.); call draw(title=trim(os),strip=.true.,add=.false.)
call draw(title=trim(os),x=inty(n1:n2,1),y=log10(inty(n1:n2,2)))
call draw(title=trim(os),x=inty(n1:n2,1),y=log10(inty(n1:n2,2)),add=.true.,symbol=.true.,nps=5,color=1)
call draw(title=trim(os),x=inty(n1:n2,1),y=log10(inty(n2,2)*exp(-inty(n2,1)+inty(n1:n2,1))),add=.true.,color=3,unsort=.true.)
endif
            endif
          endif
        enddo

!BS
!       do n=1,size(decay(1)%layer%grain_vol)
!         if (decay(1)%traj_dump(n)) then
!           decay(1)%gamma(n)=1.d11/(exp(decay(1)%bltz(NTRAJ,3,n))-1.d0-decay(1)%bltz(NTRAJ,3,n))
!         endif
!       enddo
!end of BS

        do il=2,n_lay
          decay(il)%gamma=decay(1)%gamma
          decay(il)%gamma1=decay(1)%gamma1
        enddo
        total_rate=sum(decay(1)%gamma)
        write(os,"('total rate:',1p,3(e17.10,x),i4)")total_rate,sum(decay(1)%gamma),maxval(decay(1)%gamma),maxloc(decay(1)%gamma); call output(os)
!       os='energy_rates.dat'; 
!       os='er.dat'; il=add_file(trim(os)); call delete_file(il)
!       open(il,file=trim(os),form='formatted',status='unknown')
!       do n=1,size(decay(1)%layer%grain_vol)
!         if (decay(1)%traj_dump(n)) then
!           write(il,"(i5,:,1p,50(x,e12.5))") n,happlied,decay(1)%bltz(nint(decay(1)%bltz(0,1,1)),1,n),decay(1)%gamma(n)
!         endif
!       enddo
!       close(il)
      else
      do il=1,n_lay
        do n=1,size(decay(il)%layer%grain_vol)
          es=0.d0; x=1.d0/(kb*temperature); n3=1
          if (decay(il)%traj_dump(n)) then
            x=decay(il)%energy(1,n)
            decay(il)%bltz(1:,3,n)=(decay(il)%bltz(1:,3,n)-x)*decay(il)%layer%grain_vol(n)
            es=(decay(il)%energy(3,n)-decay(il)%energy(1,n))*decay(il)%layer%grain_vol(n)
!           x=decay(il)%layer%grain_vol(n)*1.76d7*old_alpha(il)*decay(il)%layer%ms(n)*decay(il)%layer%ms_scale(n)/(1.d0+old_alpha(il)*old_alpha(il))
            x=decay(il)%layer%grain_vol(n)*decay(il)%layer%gamma*old_alpha(il)*decay(il)%layer%ms(n)*decay(il)%layer%ms_scale(n)/(1.d0+old_alpha(il)*old_alpha(il))
            decay(il)%bltz(1:,4,n)=decay(il)%bltz(1:,4,n)*x

            x=1.d0/(kb*temperature); n3=1
            err_p=(es*x.lt.300.d0)
            if (.not.err_p) decay(il)%gamma(n)=0.d0
            do while (err_p); err_p=.false.
              if (n.eq.TRAJ) print *,'a:',decay(il)%bltz(:,4,n),(es.lt.300.d0),es
              do n2=1,20
        !B's
                decay(il)%bltz_coef(n3,2,n)=-1.d0
                do n1=n3+1,nitersaddle*NTRAJ
                  y=0.5d0*(decay(il)%bltz(n1-1,3,n)+decay(il)%bltz(n1,3,n))
                  decay(il)%bltz_coef(n1,2,n)=decay(il)%bltz_coef(n1-1,2,n)*(2.d0*decay(il)%gamma(n)-decay(il)%bltz(n1-1,4,n)*x)/ &
                         (2.d0*decay(il)%gamma(n)-decay(il)%bltz(n1,4,n)*x)* &
                         exp(y*decay(il)%gamma(n)*(1.d0/decay(il)%bltz(n1,4,n)-1.d0/decay(il)%bltz(n1-1,4,n)))
                enddo
        !A(N)
                n1=nitersaddle*NTRAJ
                decay(il)%bltz_coef(n1,1,n)=-decay(il)%bltz_coef(n1,2,n)*exp(es*(x-2.d0*decay(il)%gamma(n)/decay(il)%bltz(n1,4,n)))

        !A's
                do n1=nitersaddle*NTRAJ-1,n3,-1
                  y=0.5d0*(decay(il)%bltz(n1+1,3,n)+decay(il)%bltz(n1,3,n)) !energy at boundary
                  decay(il)%bltz_coef(n1,1,n)=decay(il)%bltz_coef(n1+1,1,n)*exp(y*decay(il)%gamma(n)*(1.d0/decay(il)%bltz(n1+1,4,n)-1.d0/decay(il)%bltz(n1,4,n)))+&
                     exp(y*(x-2.d0*decay(il)%gamma(n)/decay(il)%bltz(n1,4,n)))*(-decay(il)%bltz_coef(n1,2,n)+ &
                       decay(il)%bltz_coef(n1+1,2,n)*exp(y*decay(il)%gamma(n)*(1.d0/decay(il)%bltz(n1,4,n)-1.d0/decay(il)%bltz(n1+1,4,n))))
                enddo
                if (n.eq.TRAJ) print *,decay(il)%gamma(n),decay(il)%bltz(n3,4,n)*decay(il)%bltz_coef(n3,2,n)*x/ &
                      (decay(il)%bltz_coef(n3,2,n)-decay(il)%bltz_coef(n3,1,n))
                old_gamma=decay(il)%gamma(n)
                decay(il)%gamma(n)=0.5d0*decay(il)%gamma(n)+0.5d0*decay(il)%bltz(n3,4,n)*decay(il)%bltz_coef(n3,2,n)*x/ &
                        (decay(il)%bltz_coef(n3,2,n)-decay(il)%bltz_coef(n3,1,n))
              enddo
              do n2=n3,nitersaddle*NTRAJ-1
                if (decay(il)%bltz_coef(n2,2,n)*decay(il)%bltz_coef(n3+1,2,n).lt.0.d0) then
                  err_p=.true.; exit
                endif
              enddo
              if (err_p) n3=n3+1
              if (n3.ge.nitersaddle*NTRAJ-1.or.(.not.err_p.and.abs(decay(il)%gamma(n)).lt.1.d-15)) then
                decay(il)%gamma(n)=0.d0; err_p=.false.
              endif
            enddo
            
            decay(il)%bltz(1:,3,n)=max(0.d0,decay(il)%bltz(1:,3,n)-mkt*kb*temperature)
            y=(exp(decay(il)%bltz(n3,3,n)*x)-1.d0)*0.5d0*(decay(il)%bltz(n3+1,3,n)-decay(il)%bltz(n3,3,n))/decay(il)%bltz(n3,4,n)

            if (n.eq.TRAJ) print *,'gamma(e):',n3,decay(il)%bltz(n3,3,n)*x,y

            do n2=n3+1,nitersaddle*NTRAJ-2
              y=y+(exp(decay(il)%bltz(n2,3,n)*x)-1.d0)*0.5d0*(decay(il)%bltz(n2+1,3,n)-decay(il)%bltz(n2-1,3,n))/decay(il)%bltz(n2,4,n)
              if (n.eq.TRAJ) print *,n2,decay(il)%bltz(n2,3,n)*x,y
            enddo
            n2=nitersaddle*NTRAJ-1
            y=y+(exp(decay(il)%bltz(n2,3,n)*x)-1.d0)*0.5d0*(decay(il)%bltz(n2,3,n)-decay(il)%bltz(n2-1,3,n))/decay(il)%bltz(n2,4,n)
decay(il)%bltz(NTRAJ,3,n)=decay(il)%bltz(n2,3,n)*x
            if (n.eq.TRAJ) print *,n2,decay(il)%bltz(n2,3,n)*x,y
            decay(il)%gamma1(n)=0.d0; if (y.gt.0.d0) decay(il)%gamma1(n)=1.d0/y
          else
            decay(il)%gamma(n)=0.d0; decay(il)%gamma1(n)=0.d0
          endif
          n1=nitersaddle*NTRAJ
          if (TRAJ.ne.0) then
!         if (.true.) then
            print "(i4,1p,:,22(x,e17.10))",n,es*x,decay(il)%gamma(n),decay(il)%gamma1(n),decay(il)%energy(1,n)*decay(il)%layer%grain_vol(n) &
                  /(kb*temperature),decay(il)%bltz(n1,4,n)*exp(-es*x)*x,decay(il)%bltz(1,4,n), &
                  decay(il)%bltz(nitersaddle*NTRAJ,4,n),decay(il)%gamma(n)/(x*decay(il)%bltz(nitersaddle*NTRAJ,4,n)*exp(-es*x))
!           write(85, "(i4,1p,:,22(x,e17.10))")n,es*x,decay(il)%gamma(n),decay(il)%gamma1(n),decay(il)%energy(1,n)*decay(il)%layer%grain_vol(n) &
!                 /(kb*temperature),decay(il)%bltz(n1,4,n)*exp(-es*x)*x,decay(il)%bltz(1,4,n), &
!                 decay(il)%bltz(nitersaddle*NTRAJ,4,n),decay(il)%gamma(n)/(x*decay(il)%bltz(nitersaddle*NTRAJ,4,n)*exp(-es*x))
          endif
          if (.not.old_gamma_p) then; y=decay(il)%gamma(n); decay(il)%gamma(n)=decay(il)%gamma1(n); decay(il)%gamma1(n)=y; endif
        enddo
        total_rate=total_rate+sum(decay(il)%gamma)
      enddo
      write(os,"('total rate:',1p,3(e17.10,x),i4)")total_rate,sum(decay(1)%gamma),maxval(decay(1)%gamma),maxloc(decay(1)%gamma); call output(os)
      endif


!if (get_plot()) then
!       x=count(decay(1)%bltz(NTRAJ,3,:).lt.max_eb.and.decay(1)%gamma.gt.0.d0)  !number of grains >0 and < max_eb
!       y=maxval(decay(1)%gamma,mask=decay(1)%bltz(NTRAJ,3,:).lt.max_eb)
!       do while (x.gt.0.d0.and.count(decay(1)%gamma.ge.y.and.decay(1)%bltz(NTRAJ,3,:).lt.max_eb)/x.lt.0.9d0)
!         y=y-(maxval(decay(1)%gamma,mask=decay(1)%bltz(NTRAJ,3,:).lt.max_eb)-minval(decay(1)%gamma,mask=decay(1)%bltz(NTRAJ,3,:).lt.max_eb))/50.d0
!        read *
!       enddo
!       x=maxval(decay(1)%gamma,mask=decay(1)%bltz(NTRAJ,3,:).lt.max_eb.and.decay(1)%gamma.gt.0.d0)
!print *,x,y
!if (x.gt.0.d0) call plot_contour('grain log10(life time)',log10(1.d0/max(y/10.d0,decay(1)%gamma)),.true.,1,n0=5,zmax=log10(1.d0/y),zmin=log10(1.d0/x))
!if (x.gt.0.d0) call plot_histogram('log10(life time)', 20,pack(log10(1.d0/decay(1)%gamma),mask=decay(1)%gamma.ge.y.and.decay(1)%bltz(NTRAJ,3,:).lt.max_eb))
!endif
!      x=x/minval(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0)
!       x=maxval(decay(1)%gamma,mask=decay(1)%bltz(NTRAJ,3,:).lt.max_eb)
!print *,x
!if (x.gt.0.d0.and.get_plot()) call plot_contour('grain log10(gamma)',log10(max(x/1.d5,decay(1)%gamma)),.true.,1,n0=5,zmin=log10(x/1.d5))
!if (x.gt.0.d0.and.get_plot()) call plot_contour('grain log10(gamma) a',log10(max(x/1.d5,decay(1)%gamma)),.true.,1,zmin=log10(x/1.d5))
!if (x.gt.0.d0.and.get_plot()) call plot_contour('grain log10(life time).',log10(1.d0/max(x/1.d5,decay(1)%gamma)),.true.,1,n0=5,zmin=log10(1.d0/x))
!if (x.gt.0.d0.and.get_plot()) call plot_histogram('log10(life time).', 20,pack(log10(1.d0/max(x/1.d5,decay(1)%gamma)),mask=decay(1)%gamma.gt.x/1.d4))
!if (x.gt.0.d0.and.get_plot()) call plot_histogram('log10(life time)', 20,pack(log10(1.d0/max(x/1.d5,decay(1)%gamma)),mask=decay(1)%gamma.gt.x/1.d4),xmin=-2.d0,xmax=2.d0)
!if (x.gt.0.d0.and.get_plot()) call plot_contour('grain log10(life time) a',log10(1.d0/max(x/1.d5,decay(1)%gamma)),.true.,1,zmin=log10(1.d0/x))

    end function

    subroutine kinetic_monte_carlo(total_rate)
      use random_m, only: random_uniform
      use plot_m, only: plot_scatter
      use graphix, only:draw
      real(8),intent(inout)::total_rate
      real(8)::vol_flip,max_vol_flip,es,x,tau !,es0
      integer::i,il,n,n_lay1,j !,ion
      character(100)::dirs,dirs0
!     logical::ok_p

      decay(1)%bltz(1,3,:)=0.d0
      max_vol_flip=0.d0
      if (multispin_p) then 
        max_vol_flip=sum(decay(1)%layer%grain_vol)
        decay(1)%energy(1,:)=msdata(:,1)
        decay(1)%energy(2,:)=msdata(:,2)
        n_lay1=1
        total_rate=sum(decay(1)%gamma)
!n=TRAJ;if (n.ne.0.and.associated(ms_energy)) then
!call draw(title=trim(os),20,pack(me(:,1,1)-me(:,0,1),mask=grain_do)/(kb*temperature))
!endif
        do n=1,size(msdata,1)
!          if (decay(1)%traj_lch(n).gt.1) then
           if (decay(1)%gamma(n).gt.0.d0) then
             write(os,"(i4,1p,:,40(x,e12.5))")n,decay(1)%bltz(NTRAJ,3,n),msdata(n,1)/(kb*temperature),msdata(n,2)/(kb*temperature),decay(1)%gamma(n)
!if(decay(1)%bltz(NTRAJ,3,n).ne.max_eb) &
! if(decay(1)%bltz(decay(1)%traj_lch(n),3,n).ne.max_eb) &
!print *,trim(os)
call output(trim(os))
           endif
        enddo
!       write(os,"('total rate:',1p,1(e17.10,x),i0)")maxval(decay(1)%gamma),maxloc(decay(1)%gamma); call output(os)
        write(os,"('total rate:',1p,3(e17.10,x),i4)")total_rate,sum(decay(1)%gamma),maxval(decay(1)%gamma),maxloc(decay(1)%gamma); call output(os)
!if (TRAJ.ne.0.and.count(decay(1)%gamma.gt.0.d0).gt.1) call draw(title='gamma vs Eb',x=pack(decay(1)%bltz(NTRAJ,3,:),mask=decay(1)%gamma.gt.0.d0),y=log10(pack(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0)))
print *,TRAJ,count(decay(1)%gamma.gt.0.d0),count(decay(1)%bltz(NTRAJ,3,:).ne.-1),size(msdata,1),minval(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0),maxval(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0)
print *,size(pack(decay(1)%bltz(NTRAJ,3,:),mask=decay(1)%gamma.gt.0.d0)),size(pack(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0))
!read *
if (TRAJ.ne.0.and.count(decay(1)%gamma.gt.0.d0).gt.1) call draw(title='gamma vs Eb',x=pack(decay(1)%bltz(NTRAJ,3,:),mask=decay(1)%gamma.gt.0.d0),y=log10(pack(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0)),symbol=.true.,nps=5,color=1)
if (TRAJ.ne.0) then;if (decay(1)%gamma(TRAJ).gt.0.d0) call draw(title='gamma vs Eb',x=(/decay(1)%bltz(NTRAJ,3,TRAJ)/),y=(/log10(decay(1)%gamma(TRAJ))/),symbol=.true.,nps=7,color=3,add=.true.);endif
!call plot_scatter('gamma vs Eb',pack(decay(1)%bltz(NTRAJ,3,:),mask=decay(1)%gamma.gt.0.d0),log10(pack(decay(1)%gamma,mask=decay(1)%gamma.gt.0.d0)),.true.)
      else
        n_lay1=n_lay; total_rate=0.d0
        do i=1,n_lay
          max_vol_flip=max_vol_flip+sum(decay(i)%layer%grain_vol)
          decay(i)%energy(1,:)=decay(i)%energy(1,:)*decay(i)%layer%grain_vol
          decay(i)%energy(2,:)=decay(i)%energy(2,:)*decay(i)%layer%grain_vol
          total_rate=total_rate+sum(decay(i)%gamma)
        enddo
      endif
      max_vol_flip=max_vol_flip*VOL_FRAC; vol_flip=0.d0
      do while (vol_flip.lt.max_vol_flip.and.total_time.lt.dt.and.total_rate.gt.0.d0)
!write(os,*)last_flip,last_cum_prob,(last_flip.gt.0.d0.or.last_cum_prob.eq.0.d0); call output(os)
        if (last_flip.gt.0.d0.or.last_cum_prob.eq.0.d0) last_cum_prob=1.d0-random_uniform(random_decay)
!write(os,*)last_flip,last_cum_prob,(last_flip.gt.0.d0.or.last_cum_prob.eq.0.d0),total_rate; call output(os)
        tau=-log(last_cum_prob)/total_rate !time to flip
!es0=random_uniform(random_decay)*total_rate
!last_cum_prob=1.d0-random_uniform(random_decay)
!x=0.d0; do n=1,size(decay(1)%layer%grain_vol); x=x+decay(1)%gamma(n); if (es0-x.le.0.d0) exit; enddo; n=min(n,size(decay(1)%layer%grain_vol))
!tau=-log(last_cum_prob)/decay(1)%gamma(n)
        if (total_time+tau.gt.dt) then
          total_time=dt
        else
          total_time=total_time+tau; last_flip=total_time
          es=random_uniform(random_decay)*total_rate; i=0
!es=es0; i=0
          do il=1,n_lay1; x=sum(decay(il)%gamma)   !find the layer the grain sits in...
            if (es-x.gt.0.d0) then; es=es-x  !grain is not in this layer, move to next (not valid for multi-spin)
            else; exit                       !this is the layer the grain is in
            endif
          enddo; x=0.d0
          do n=1,size(decay(il)%layer%grain_vol)   !now find the actual grain
            x=x+decay(il)%gamma(n); if (es-x.le.0.d0) exit
          enddo; n=min(n,size(decay(il)%layer%grain_vol)); dirs0=''; dirs=''; 
          if (multispin_p) then
            do j=1,n_lay; write(dirs,"(a,l)") trim(dirs),decay(j)%m(3,n).gt.0.d0; write(dirs0,"(a,l)") trim(dirs0),decay(j)%m_min(3,n,4).gt.0.d0;enddo
          else; write(dirs,"(l)") decay(il)%m(3,n).gt.0.d0; write(dirs0,"(l)") decay(il)%m_min(3,n,4).gt.0.d0
          endif

!print *,n,il,tau,NTRAJ
!print *,decay(il)%bltz(NTRAJ,3,n)
!print *,temperature
!print *,kb
!       if (decay(il)%gamma(n)*dt.gt.1.d0.and. &   !if a quick flipper, then allow it to flip back...
          if (boltz_p.and. &                                     !allow the grain to flip back...
              random_uniform(random_decay).le.1.d0/(1.d0+exp((decay(il)%energy(1,n)-decay(il)%energy(2,n))/(kb*temperature)))) then; j=NTRAJ; if (multispin_p) j=decay(1)%traj_lch(n)
                  write(os,"('grain/layer, time to flip, gamma, 1/gamma:',i4,'/',i1,1p,:,7(e12.5,x),0p,i5,' but flipped back! ',a)") &
                      n,il,tau,decay(il)%gamma(n),1.d0/decay(il)%gamma(n),decay(il)%bltz(j,3,n),decay(il)%energy(1,n)/(kb*temperature),decay(il)%energy(2,n)/(kb*temperature), &
                      1.d0/(1.d0+exp((decay(il)%energy(1,n)-decay(il)%energy(2,n))/(kb*temperature))),decay(il)%num_flip(n),trim(dirs)//'<-'//trim(dirs0)
!           if (multispin_p) decay(1)%bltz(1,3,n)=-1.d0
!     write(82)
          else
            if (multispin_p) then
              do j=1,n_lay; decay(j)%m(:,n)=decay(j)%m_min(:,n,4); enddo
!             decay(1)%bltz(1,3,n)=1.d0
            else
              decay(il)%m(:,n)=decay(il)%m_min(:,n,4)
            endif
            decay(il)%num_flip(n)=decay(il)%num_flip(n)+1
            vol_flip=vol_flip+decay(il)%layer%grain_vol(n)*(2*mod(decay(il)%num_flip(n),2)-1)
            j=NTRAJ; if (multispin_p) j=decay(1)%traj_lch(n)
            write(os,"('grain/layer, time to flip, gamma, 1/gamma:',i4,'/',i1,1p,:,7(e12.5,x),0p,i5,x,a)") &
                    n,il,tau,decay(il)%gamma(n),1.d0/decay(il)%gamma(n),decay(il)%bltz(j,3,n),decay(il)%energy(1,n)/(kb*temperature),decay(il)%energy(2,n)/(kb*temperature), &
                    1.d0/(1.d0+exp((decay(il)%energy(1,n)-decay(il)%energy(2,n))/(kb*temperature))),decay(il)%num_flip(n),trim(dirs)//'->'//trim(dirs0)
!     write(83)
          endif
          call output(os)
          decay(il)%gamma(n)=0.d0; total_rate=0.d0; do il=1,n_lay1; total_rate=total_rate+sum(decay(il)%gamma); enddo
        endif
      enddo
      if (vol_flip.ge.max_vol_flip) call output(' maximum volume flipped')
      if (total_time.ge.dt) call output(' total time elapsed')
      if (total_rate.le.0.d0) total_time=dt
      if (total_rate.le.0.d0) call output(' no more spins to flip')
    end subroutine

     subroutine fan_trajectory_ms(m_min0,m_min1,grain_do,mlayer,ith,num,il,m_initial,onesided_p)  !start, end, layer for this m, %in, traj, ntraj, layer
      use transform_m, only:cross
      real(8),intent(in)::m_min0(:,:,:),m_min1(:,:,:)
      real(8),intent(out),optional::m_initial(:,:)
      logical,intent(in),optional::onesided_p
      logical,intent(in)::grain_do(:)
      integer,intent(in)::ith(:),num,il,mlayer(:)
      integer::i
      real(8)::m(3)
      real(8)::m0(3),m1(3),zp(3),x
      logical::oneside_p
      
      oneside_p=.false.; if (present(onesided_p)) oneside_p=onesided_p
      do i=1,size(mlayer)
       if (grain_do(i)) then
        m0=m_min0(:,i,mlayer(i)); m1=m_min1(:,i,mlayer(i))
        zp = cross(cross(m0,m1),m0)
        if (maxval(abs(zp)).eq.0.d0) then
          zp(minloc(abs(m0),DIM=1))=1.d0        
          zp=cross(cross(m0,zp),m0)
        endif
        x=sqrt(dot_product(zp,zp)); zp=zp/x
        x=acos(min(1.d0,max(-1.d0,dot_product(m0,m1))))
        if (x.eq.0.d0) x=acos(-1.d0)
        if (oneside_p) then
          x=x*cos((num-1-ith(i))*acos(-1.d0)/2.d0/(num-1))
        else
          if (ith(i).le.(n-1)/2) then
            x=x*0.5d0*(1.d0-cos(ith(i)*acos(-1.d0)/(num-1)))
          else 
            x=x*(1.d0-0.5d0*(1.d0-cos(((num-1)-ith(i))*acos(-1.d0)/(num-1))))
          endif
        endif
        m=m0*cos(x)+zp*sin(x)
        x=dot_product(m,m)
        if (abs(1.d0-x).gt.1.d-12) print "('fan: |m|=',4(x,i0),:,1p,52(x,e12.5))",i,mlayer(i),size(mlayer),size(m_min0),abs(1.d0-x),abs(1.d0-sqrt(x)),dot_product(m0,m0),dot_product(m1,m1),dot_product(cross(m0,m1),cross(m0,m1)),m0,m1,m
!if (TRAJ.eq.i) print *,i,ith(i),num,m/sqrt(x),dot_product(m/sqrt(x),m0)
        decay(mlayer(i))%layer%m(1+decay(mlayer(i))%layer%in:3+decay(mlayer(i))%layer%in,i)=m/sqrt(x)
!if (TRAJ.eq.i) print "('f> ',5(i0,x),:,1p,66(e12.5,x))",i,mlayer(i),decay(mlayer(i))%layer%in,ith(i),num,m0,m1,m/sqrt(x),decay(mlayer(i))%layer%m(1+decay(mlayer(i))%layer%in:3+decay(mlayer(i))%layer%in,i)
        if (present(m_initial)) m_initial(:,i)=m/sqrt(x)
       endif
      enddo
      i=TRAJ; if (i.le.0) i=1
!     if (il.ne.0.and.TRAJ.ne.0) print "('fan #',i4,' out of',i4,1p,2(x,e12.5))",ith(i),num,acos(min(1.d0,max(-1.d0, & 
!           dot_product(m_min0(:,i),m(1+id:3+id,i)))))*180.d0/acos(-1.d0), & 
!           acos(min(1.d0,max(-1.d0,dot_product(m_min1(:,i),m(1+id:3+id,i)))))*180.d0/acos(-1.d0)
!    read *
     end subroutine

    subroutine fan_trajectory(m_min0,m_min1,m,id,ith,num,il,m_initial,onesided_p)
     use transform_m, only:cross
     real(8),intent(in)::m_min0(:,:),m_min1(:,:)
     real(8),intent(inout)::m(:,:)
     real(8),intent(out),optional::m_initial(:,:)
     logical,intent(in),optional::onesided_p
     integer,intent(in)::ith(:),num,id,il

     integer::i,j
     real(8)::m0(3),m1(3),zp(3),yp(3),m2(3),x,xp(3)
     logical::colnr_p,oneside_p

     oneside_p=.false.; if (present(onesided_p)) oneside_p=onesided_p
     do i=1,size(m,2)
       m0=m_min0(:,i); m1=m_min1(:,i); zp = cross(m0,m1)
       colnr_p=(abs(1.d0-dot_product(m0,m1)).lt.1.d-7.or.maxval(abs(zp)).eq.0.d0)
       if (colnr_p) then
         j=minloc(abs(m0),DIM=1); zp=0.d0; zp(j)=1.d0; zp=cross(m0,zp)
       endif
       x=sqrt(dot_product(zp,zp))
       zp=zp/x; xp=m0; yp=cross(zp,xp)
       x=acos(min(1.d0,max(-1.d0,dot_product(m0,m1))))
       if (colnr_p) x= (1.d0*acos(-1.d0)-acos(min(1.d0,max(-1.d0,dot_product(m0,m1)))))
       if (oneside_p) then
         x=x*cos((num-1-ith(i))*acos(-1.d0)/2.d0/(num-1))
       else
         if (ith(i).le.(n-1)/2) then
           x=x*0.5d0*(1.d0-cos(ith(i)*acos(-1.d0)/(num-1)))
         else
           x=x*(1.d0-0.5d0*(1.d0-cos(((num-1)-ith(i))*acos(-1.d0)/(num-1))))
         endif
       endif
       m2=(/ cos(x), sin(x), 0.d0 /)
  x=dot_product(m(1+id:3+id,i),m(1+id:3+id,i))
  if (abs(1.d0-x).gt.1.d-12) print "('before fan: |m|=',i3,1p,2(x,e12.5))",i,abs(1.d0-x),abs(1.d0-sqrt(x))
       m(1+id:3+id,i)= &
          (/ xp(1)*m2(1)+yp(1)*m2(2)+zp(1)*m2(3), &
             xp(2)*m2(1)+yp(2)*m2(2)+zp(2)*m2(3), &
             xp(3)*m2(1)+yp(3)*m2(2)+zp(3)*m2(3) /)
       x=dot_product(m(1+id:3+id,i),m(1+id:3+id,i))
       if (abs(1.d0-x).gt.1.d-12) print "('fan: |m|=',i3,1p,2(x,e12.5))",i,abs(1.d0-x),abs(1.d0-sqrt(x))
       m(1+id:3+id,i)=m(1+id:3+id,i)/sqrt(x)
       if (present(m_initial)) m_initial(:,i)=m(1+id:3+id,i)
     enddo
     i=TRAJ; if (i.le.0) i=1
     if (il.eq.1) print "('fan #',i4,' out of',i4,1p,2(x,e12.5))",ith(max(1,i)),num,acos(min(1.d0,max(-1.d0, &
         dot_product(m_min0(:,max(1,i)),m(1+id:3+id,max(1,i))))))*180.d0/acos(-1.d0), &
         acos(min(1.d0,max(-1.d0,dot_product(m_min1(:,max(1,i)),m(1+id:3+id,max(1,i))))))*180.d0/acos(-1.d0)
!    read *
    end subroutine
!    subroutine e_surface(t,n0,x,emin)
!     use field_m, only:get_field
!     real(8)::t,emin
!     integer,intent(in),optional::n0
!
!     real(8)::x,phi,theta,m0(3)
!     integer::i,j,n
!
!     n=1; if (present(n0)) n=n0
!     x=decay(1)%layer%grain_vol(n)/(kb*temperature)
!     m0=decay(1)%layer%m(1+decay(1)%layer%in:3+decay(1)%layer%in,n)
!     write(os,"('e_',i4.4,'.dat')") n
!     open(23,file=trim(adjustl(os)),status='unknown')
!     write(23,"('VARIABLES = ""phi"" ""theta"" ""m_x"" ""m_y"" ""m_z"" ""h_x"" ""h_y"" ""h_z"" ""energy/kT"" ""|mxH|""')")
!     write(23,"('ZONE T=""cubic""')")
!     write(23,"('I=',i4,', J=',i4)") 401,201
!     do j=0,200
!       theta = j*acos(-1.d0)/200
!       do i=0,400
!         phi = i*acos(-1.d0)/200
!         decay(1)%layer%m(1+decay(1)%layer%in:3+decay(1)%layer%in,n) = (/ cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta) /)
!         call get_field(t)
!         write(23,"(1p,10(e12.5,x))") phi,theta,decay(1)%layer%m(1+decay(1)%layer%in:3+decay(1)%layer%in,n), &
!             decay(1)%layer%h(:,n),(decay(1)%layer%energy(n)-emin)*x,&
!             sqrt( (decay(1)%layer%m(2+decay(1)%layer%in,n)*decay(1)%layer%h(3,n)-decay(1)%layer%m(3+decay(1)%layer%in,n)*decay(1)%layer%h(2,n))**2+ &
!             (decay(1)%layer%m(3+decay(1)%layer%in,n)*decay(1)%layer%h(1,n)-decay(1)%layer%m(1+decay(1)%layer%in,n)*decay(1)%layer%h(3,n))**2+ &
!             (decay(1)%layer%m(1+decay(1)%layer%in,n)*decay(1)%layer%h(2,n)-decay(1)%layer%m(2+decay(1)%layer%in,n)*decay(1)%layer%h(1,n))**2)
!       enddo
!     enddo
!     decay(1)%layer%m(1+decay(1)%layer%in:3+decay(1)%layer%in,n)=m0
!     close(23)
!    end subroutine
  end subroutine
end module decay_m
