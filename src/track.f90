module track_m
 use data_m, only: data_s

  use graphix, only: draw


 implicit none
 private

 public::init_track,set_track,reset_track,destroy_track,do_tracking,track_s,dump_track

 integer,save::TRAJ,one_lay
 real(8),save::min_dt,max_dt,dt0
 logical,save::calc_density_p,rot_p

 type::orbit_s
   integer::n_orbit
   real(8)::min_dm
   real(8),dimension(:,:),pointer::orbit_param=>null(),orbit_m=>null()
 end type orbit_s

 type::track_s
   real(8),dimension(:),pointer::period=>null(),angle=>null(),dm_dt2=>null(),ddm=>null(),heff2=>null(),shxm2=>null(),hxm2=>null(),energy=>null()
   real(8),dimension(:,:),pointer::dm_dt=>null(),heff=>null(),hxm=>null(),start=>null(),hxm0=>null(),hxm1=>null(),m_min0=>null(),m_min1=>null(),fe=>null()
   type(orbit_s),dimension(:),pointer::orbit=>null()
   logical,dimension(:),pointer::orbit_done=>null()
   type(data_s),pointer::layer=>null()
   type(track_s),pointer::next=>null()
 end type track_s
 real(8),dimension(:),pointer::ms_energy=>null(),ms_lm=>null()
 real(8),dimension(:,:),pointer::ms_fe=>null()
 real(8),parameter::dfe=25.d0/50.d0*1.380662d-16*300.d0


 real(8),dimension(:),pointer::ms_energy0=>null()
  logical,save::xgidp=.false.
  character(50),dimension(0+2),save::xgid !#layer + 2

contains

 subroutine set_track(trackk)
  type(data_s),pointer::layer
  type(track_s),pointer::track,trackk

  track=>trackk
  do while (associated(track))
    layer=>track%layer
    track%m_min0(:,:)=layer%m(4-layer%in:6-layer%in,:)
    track%m_min1(:,:)=layer%m(1+layer%in:3+layer%in,:)
    track=>track%next
  enddo
  call reset_track(trackk)
 end subroutine

 function init_track(traji) result (trackk)
  use data_m, only: get_layer, get_data
  integer,intent(in)::traji
  type(track_s),pointer::track,trackk
  type(data_s),pointer::layer
  integer::i,n

  nullify(trackk,track)
  TRAJ=traji;
  i=1; layer=>get_layer(i); call get_data(i,ms_energy,ms_energy_p=.true.)
  if (associated(ms_energy)) allocate(ms_fe(0:100,size(layer%grain_vol)),ms_lm(size(layer%grain_vol)),ms_energy0(size(layer%grain_vol)))
  do while (associated(layer))
    if (associated(trackk)) then
      allocate(track%next); track=>track%next
    else
      allocate(track); trackk=>track
    endif
    track%layer=>layer; n=size(layer%grain_vol)
    allocate(track%period(n),track%dm_dt(3,n),track%heff(3,n),track%hxm(3,n), &
             track%dm_dt2(n),track%ddm(n),track%heff2(n),track%hxm2(n),track%shxm2(n),track%angle(n), &
             track%start(3,n),track%orbit_done(n),track%energy(n), &
             track%hxm0(3,n),track%hxm1(6,n),track%m_min0(3,n),track%m_min1(3,n),track%fe(0:100,n))
    track%m_min0(:,:)=layer%m(4-layer%in:6-layer%in,:); track%fe=0.d0
    track%m_min1(:,:)=layer%m(1+layer%in:3+layer%in,:)
    i=i+1; layer=>get_layer(i)
  enddo
  call reset_track(trackk)
 end function
 function destroy_track(trackk) result (track)
  type(track_s),pointer::trackk,track
  integer::i
  track=>trackk
  do while (associated(track))
    if (associated(track%orbit)) then
      do i=1,size(track%orbit)
        if (associated(track%orbit(i)%orbit_param)) deallocate(track%orbit(i)%orbit_param,track%orbit(i)%orbit_m)
      enddo
      deallocate(track%orbit)
    endif
    deallocate(track%period,track%dm_dt,track%heff,track%hxm, &
               track%dm_dt2,track%ddm,track%heff2,track%hxm2,track%shxm2,track%angle, &
               track%start,track%orbit_done,track%energy, &
               track%hxm0,track%hxm1,track%m_min0,track%m_min1)
    track=>track%next; deallocate(trackk); trackk=>track
  enddo
  nullify(trackk,track)
 end function
 subroutine reset_track(trackk,calc_density_pi,test_rot_pi,max_dti,min_dti,one_layer)
  use data_m, only: get_data
  type(track_s),pointer::trackk,track
  real(8),intent(in),optional::min_dti,max_dti
  logical,intent(in),optional::calc_density_pi,test_rot_pi
  integer,intent(in),optional::one_layer

  type(data_s),pointer::layer
  integer::n

 if (TRAJ.ne.0) then
 if (.not.xgidp) then
   do n=1,size(xgid)-2; write(xgid(n),"('|mxH|^2 layer ',i0)") n; call draw(trim(xgid(n)),strip=.true.); enddo; xgidp=.true.
   n=size(xgid)-1; xgid(n)='weighted m0.(mxdm/dt)'; call draw(trim(xgid(n)),strip=.true.)
   n=size(xgid);   xgid(n)='weighted |mxH|^2';      call draw(trim(xgid(n)),strip=.true.)
   call draw('ms_energy',strip=.true.)
!else
!  do n=1,size(xgid); call draw(trim(xgid(n)),strip=.true.,add=.false.); enddo
!  call draw('ms_energy',strip=.true.,add=.false.)
 endif
!read *
endif

  call get_data(1,ms_energy,ms_energy_p=.true.)

!print "('associated(ms_energy)=',L)",associated(ms_energy)
  one_lay=0; if (present(one_layer)) one_lay=one_layer
  rot_p=.false.; if (present(test_rot_pi)) rot_p=test_rot_pi
  calc_density_p=.false.; if (present(calc_density_pi)) calc_density_p=calc_density_pi
  min_dt=0.d0; if (present(min_dti)) min_dt=min_dti
  max_dt=1.d-9+min_dt; if (present(max_dti)) max_dt=max_dti
  if (associated(ms_lm)) ms_lm=1.d22
  if (associated(ms_energy0)) ms_energy0=1.d22
  track=>trackk
  do while (associated(track))
    layer=>track%layer
    track%period=-1.d0; track%dm_dt=0.d0; track%dm_dt2=0.d0; track%heff=0.d0; track%ddm=0.d0
    track%heff2=0.d0; track%hxm=0.d0; track%hxm2=0.d0; track%shxm2=0.d0; track%angle=0.d0; track%start=0.d0; track%orbit_done=.false.;
    track%hxm0=0.d0; track%hxm1=0.d0; track%hxm0(1,:)=1.d22; track%hxm0(2,:)=-1.d22; track%hxm0(3,:)=1.d22
    if (calc_density_p .and. .not. associated(track%orbit)) then
      allocate(track%orbit(size(layer%grain_vol)))
      do n=1,size(track%orbit)
        allocate(track%orbit(n)%orbit_param(3,300),track%orbit(n)%orbit_m(3,300))
      enddo
    endif
    track=>track%next
  enddo
 end subroutine
 subroutine dump_track(trackk,energy,mxh2,lm,mstart,dmp,diff_m,msfe)
  use solve1d_m, only: solve1d
  type(track_s),pointer::trackk,track
  real(8),intent(out)::energy(:),mxh2(:),lm(:),mstart(:,:)
  real(8),intent(in),optional::diff_m
  real(8),intent(out),optional::msfe(:,:)
  logical,intent(in)::dmp(:)

  integer::n,i,j

  real(8),allocatable,dimension(:)::a,b,c,r,x,dxi
  real(8)::zz,yy,xx,ww,p
  real(8)::diffm,diffm0

!print "('associated(ms_energy)=',L,L)",associated(ms_energy),associated(ms_lm)

  diffm0=1.d8 !FIXME
  if (present(diff_m)) diffm0=diff_m

  if (associated(ms_energy).and.present(msfe)) msfe=ms_fe

! if (TRAJ.ne.0) then
!   track=>trackk
!   if (associated(track)) then
!     n=TRAJ
!  print *,n,track%layer%do_grain(n)
!     if (track%layer%do_grain(n)) then
!       p=1.d0/track%period(n)
!       write(92,"(i4,1p,20(x,e12.5))")n,track%period(n),track%energy(n),track%dm_dt2(n)*p,track%heff2(n),track%hxm2(n), &!track%hxm2(n)*p, &
!               track%hxm0(:,n),track%hxm1(:,n),track%ddm(n)
!       print "(i4,1p,20(x,e12.5))",n,track%period(n),track%energy(n),track%dm_dt2(n)*p,track%heff2(n),track%hxm2(n), &!track%hxm2(n)*p, &
!               track%hxm0(:,n),track%hxm1(:,n),track%ddm(n)
!     endif
!   endif
! endif

  !FIXME- this assumes just one layer, otherwise we get only the last layer of spins
  track=>trackk
! do while (associated(track))
  if (associated(track)) then
    do n=1,size(track%layer%grain_vol)
      if (dmp(n)) then
        energy(n)=track%energy(n)
!       mxh2(n)=track%hxm2(n)/track%period(n)
        mxh2(n)=track%hxm2(n)!/track%period(n)
        mstart(:,n)=track%hxm1(:,n)
        if (associated(ms_lm)) then
          lm(n)=ms_lm(n)!/track%period(n)
        else
          lm(n)=track%hxm0(3,n)/track%period(n)
!if (ms_lm(n).le.0.d0) mstart(1:3,n)=track%hxm(:,n)
        endif
        if (calc_density_p) then
          i=track%orbit(n)%n_orbit
          if (.not.allocated(a)) then
            allocate(a(i*3/2),b(i*3/2),c(i*3/2),r(i*3/2),x(i*3/2),dxi(i*3/2))
          else if (size(a).lt.i) then
            deallocate(a,b,c,r,x,dxi)
            allocate(a(i*3/2),b(i*3/2),c(i*3/2),r(i*3/2),x(i*3/2),dxi(i*3/2))
          endif

          diffm=diffm0/(track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)) !(track%layer%grain2gridindex(n)-track%layer%grain2gridindex(n-1))

          dxi(1)=(track%orbit(n)%orbit_param(1,1)-track%orbit(n)%orbit_param(1,2))
          c(1)=diffm/dxi(1)
          do j=2,i
            a(j)=c(j-1); dxi(j)=(track%orbit(n)%orbit_param(1,j+1)-track%orbit(n)%orbit_param(1,j))
            c(j)=diffm/dxi(j); b(j)=-(c(j)+a(j))
          enddo
          do j=2,i
            a(j)=a(j)+track%orbit(n)%orbit_param(2,j-1)
            c(j)=c(j)-track%orbit(n)%orbit_param(2,mod(j,i)+1)
          enddo
          r=0.d0; x=1.d0; a(1)=0.d0; c(1)=0.d0; b(1)=1.d0
          r(1)=diffm/(track%orbit(n)%orbit_param(2,1)*track%orbit(n)%orbit_param(1,i+1))
          call solve1d(a,b,c,r,x,i,n.eq.TRAJ.and..false.)
          zz=0.d0; yy=0.d0; ww=0.d0
          xx=0.5d0*(dxi(1)+dxi(i))
          do j=1,i
            ww=ww+x(j)*xx; zz=zz+x(j)*xx; yy=yy+x(j)*xx*track%orbit(n)%orbit_param(3,j)
            xx=0.5d0*(dxi(mod(j,i)+1)+dxi(j))
          enddo
          mxh2(n)=yy/zz

!         if (n.eq.TRAJ) then
!           write(71,"('ZONE T = ""trajectory density"", I=',i4)")i
!           do j=1,i
!             write(71,"(i4,1p,:,40(x,e17.10))")j,track%orbit(n)%orbit_m(:,j),sum(dxi(1:j)),x(j),x(j)/ww,sum(dxi(1:i))*x(j)/ww,sum(dxi(1:i))
!           enddo
!         endif
        endif
      endif
    enddo
!   track=>track%next
! enddo
    trackk=>track%next
  endif
  if (allocated(a)) deallocate(a,b,c,r,x,dxi)
 end subroutine
 function do_tracking(trackk,t,dt,fe_p) result(go_p)
  use transform_m, only: cross
  logical::go_p
  type(track_s),pointer::track,trackk
  real(8),intent(in)::t,dt
  logical,intent(in),optional::fe_p

  real(8)::m(3),m0(3),dm(3),mh(3),msv,dm0(3)
  real(8),dimension(:,:),pointer::new_orbit_param
  integer::n,id,i,nla
  real(8)::lmx,ms_lmx(size(trackk%layer%grain_vol),5)

!real(8)::ms_lmi(size(trackk%layer%grain_vol))
!if (associated(ms_lm)) ms_lmi=ms_lm

!print "('associated(ms_energy)=',L)",associated(ms_energy)
  dt0=dt; go_p=.false.; track=>trackk; ms_lmx=0.d0; nla=0
  do while (associated(track))
    id=track%layer%in; nla=nla+1
    do n=1,size(track%layer%grain_vol)
     if (track%layer%do_grain(n)) then
      m=track%layer%m(1+id:3+id,n); m0=track%layer%m(4-id:6-id,n)
      dm=m-m0
      msv=sqrt(dot_product(track%layer%h(:,n),track%layer%h(:,n)));dm0=dm*msv/(msv+track%layer%hthm_m(n))
      if (track%period(n).lt.0.d0) then !start tracking
        track%start(:,n)=m0; track%angle(n)=dot_product(m,m0); track%period(n)=t-dt
        track%orbit_done(n)=.false.; track%energy(n)=track%layer%energy(n)
        if (calc_density_p) then
          track%orbit(n)%min_dm=sqrt(dot_product(dm,dm))*0.95d0/2.d0
          track%orbit(n)%n_orbit=0; track%orbit(n)%orbit_param(1,1)=-1.d22
        endif
        if (associated(ms_energy)) then
          ms_fe=0.d0
          ms_energy0(n)=ms_energy(n)
!          print *,n
!          print *,TRAJ
!          print *,associated(track%next)
if (TRAJ.eq.n.and.(.not.associated(track%next).or.nla.eq.one_lay)) then
  do i=1,size(xgid); call draw(trim(xgid(i)),strip=.true.,add=.false.); enddo
  call draw('ms_energy',strip=.true.,add=.false.)
endif
        endif

      endif
      i=min(size(track%fe,1)-1,max(0,int(track%layer%energy(n)/dfe))); track%fe(i,n)=track%fe(i,n)+dt
      if (associated(ms_energy).and.associated(track,trackk)) then; i=min(size(ms_fe,1)-1,max(0,int(ms_energy(n)/dfe))); ms_fe(i,n)=ms_fe(i,n)+dt; endif
      if (present(fe_p)) then; if (fe_p) cycle; endif

      if (.not.track%orbit_done(n)) then
        msv=track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
        ms_lmx(n,5)=ms_lmx(n,5)+msv
        track%dm_dt(:,n)=track%dm_dt(:,n)+dm  !sum of dm
!       track%dm_dt2(n)=track%dm_dt2(n)+dot_product(dm,dm) ! sum of |dm|MsV ! sum of |dm|^2
!       track%dm_dt2(n)=track%dm_dt2(n)+sqrt(dot_product(dm,dm))*msv ! sum of |dm|MsV ! sum of |dm|^2
        track%dm_dt2(n)=track%dm_dt2(n)+sqrt(dot_product(dm0,dm0))*msv ! sum of |dm|MsV ! sum of |dm|^2
        track%heff(:,n)=track%heff(:,n)+(track%layer%h(:,n)-track%heff(:,n))*dt/(t-track%period(n))  !time average of Heff
        track%heff2(n)=track%heff2(n)+(dot_product(track%layer%h(:,n),track%layer%h(:,n))-track%heff2(n))*dt/(t-track%period(n)) !time average of |Heff|^2*dt

        mh = cross(m0,track%layer%h(:,n))

!if (n.eq.TRAJ) track%shxm2(n)=track%shxm2(n)+(dot_product(mh,mh)-track%hxm2(n))**2*dt*max(0.d0,t-track%period(n)-dt)/(t-track%period(n))
                track%hxm2( n)=track%hxm2( n)+(dot_product(mh,mh)-track%hxm2(n))   *dt                               /(t-track%period(n))  !time average of |mxH|^2

!if (n.eq.TRAJ) call draw(trim(xgid(nla)),strip=.true.,add=.true.,x=(/ 1.d0, 1.d0, 1.d0, 1.d0 /)*(t-track%period(n))*1.d9,y=track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)*1.d21*(/ track%hxm2(n), track%hxm2(n)+sqrt(track%shxm2(n)/(t-track%period(n))),track%hxm2(n)-sqrt(track%shxm2(n)/(t-track%period(n))),dot_product(mh,mh)/))

!       lmx=dot_product( cross(m0,dm),track%m_min0(:,n))*track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
        lmx=dot_product( cross(m0,dm0),track%m_min0(:,n))*track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
        ms_lmx(n,1)=ms_lmx(n,1)+lmx
if (n.eq.TRAJ) ms_lmx(n,2)=ms_lmx(n,2)+     track%hxm2( n)*track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
if (n.eq.TRAJ) ms_lmx(n,3)=ms_lmx(n,3)+sqrt(track%shxm2(n)/(t-track%period(n)))*track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
if (n.eq.TRAJ) ms_lmx(n,4)=ms_lmx(n,4)+dot_product(mh,mh)*track%layer%ms(n)*track%layer%ms_scale(n)*track%layer%grain_vol(n)
        if (.not.associated(ms_energy)) then
!if (n.eq.TRAJ) print *,track%layer%energy(n),track%layer%alpha,m,dot_product(m,m)
          track%hxm0(3,n)=min(track%hxm0(3,n),lmx)
        elseif (.not.associated(track%next).or.nla.eq.one_lay) then   !last layer
          ms_lm(n)=min(ms_lm(n),ms_lmx(n,1))
 !        track%layer%do_grain(n)=(.not.(ms_lm(n).le.0.d0))  !strop tracking if we think we can flip...
!         track%layer%do_grain(n)=(ms_lm(n).gt.0.d0.and.(t-track%period(n).lt.min_dt.or.track%dm_dt2(n)/(ms_lmx(n,5)+1.d-300).lt.6.28d0))  !strop tracking if we think we can flip...
          track%layer%do_grain(n)=(t-track%period(n).lt.min_dt.or.track%dm_dt2(n)/(ms_lmx(n,5)+1.d-300).lt.6.28d0)  !strop tracking once we go 2 pi 
          if (track%layer%do_grain(n).and.rot_p) track%layer%do_grain(n)=(ms_lm(n).gt.0.d0)  !strop tracking if we think we can flip...

if (n.eq.TRAJ) call draw(trim(xgid(size(xgid)-1)),strip=.true.,add=.true.,x=(/ 1.d0, 1.d0 /)*(t-track%period(n))*1.d9,y=1.d21*(/ ms_lmx(n,1), ms_lm(n) /))
if (n.eq.TRAJ) track%shxm2(n)=track%shxm2(n)+(ms_lmx(n,4)-ms_lmx(n,2))**2*dt*max(0.d0,t-track%period(n)-dt)/(t-track%period(n))
if (n.eq.TRAJ) call draw(trim(xgid(size(xgid))),strip=.true.,add=.true.,x=(/ 1.d0, 1.d0, 1.d0, 1.d0 /)*(t-track%period(n))*1.d9,y=1.d21*(/ ms_lmx(n,2), ms_lmx(n,2)+sqrt(track%shxm2(n)/(t-track%period(n))),ms_lmx(n,2)-sqrt(track%shxm2(n)/(t-track%period(n))),ms_lmx(n,4) /))
if (n.eq.TRAJ) call draw('ms_energy',strip=.true.,add=.true.,x=(/ 1.d0 /)*(t-track%period(n))*1.d9,y=(/ (ms_energy(n)-ms_energy0(n))/4.14d-14 /))
        endif


!       track%hxm(:,n)=track%hxm(:,n)+mh*dt
!gp!       track%hxm(:,n)=track%hxm(:,n)+(mh-track%hxm(:,n))*dt/(t-track%period(n)) ! time average of mxH
!       track%hxm2(n)=track%hxm2(n)+dot_product(mh,mh)*dt


        if (track%hxm0(1,n).gt.dot_product( mh,mh )) track%hxm1(4:6,n)=m0
        track%hxm0(1,n)=min(track%hxm0(1,n),dot_product( mh,mh ))
        if (track%hxm0(2,n).lt.dot_product( mh,mh )) track%hxm1(1:3,n)=m0
        track%hxm0(2,n)=max(track%hxm0(2,n),dot_product( mh,mh ))

!       if (n.eq.TRAJ) write(70,"(i4,1p,20(x,e12.5))")n,m,dot_product(mh,mh),sqrt(dot_product(mh,mh)), &
!              dot_product( (/ m0(2)*dm(3)-m0(3)*dm(2), m0(3)*dm(1)-m0(1)*dm(3), m0(1)*dm(2)-m0(2)*dm(1) /),track%m_min0(:,n)), &
!              dot_product(dm,dm),track%ddm(n)+0.5d0*sqrt(dot_product(dm,dm))

        track%ddm(n)=track%ddm(n)+sqrt(dot_product(dm,dm)) !sum of |dm|

        if (calc_density_p) then
          i=track%orbit(n)%n_orbit
          if (i+1.gt.size(track%orbit(n)%orbit_param,2)) then
            allocate(new_orbit_param(3,2*(i+1))); new_orbit_param(:,1:i)=track%orbit(n)%orbit_param
            deallocate(track%orbit(n)%orbit_param); track%orbit(n)%orbit_param=>new_orbit_param;nullify(new_orbit_param)
            allocate(new_orbit_param(3,2*(i+1))); new_orbit_param(:,1:i)=track%orbit(n)%orbit_m
            deallocate(track%orbit(n)%orbit_m); track%orbit(n)%orbit_m=>new_orbit_param;nullify(new_orbit_param)
          endif
        endif

        if ((dot_product(track%start(:,n),track%layer%m(1+id:3+id,n)).ge.track%angle(n).and.track%period(n).ne.t-dt.and.  &
            t-track%period(n).ge.min_dt.and..not.associated(ms_energy)).or.t-track%period(n).gt.max_dt.or..not.track%layer%do_grain(n)) then
          track%orbit_done(n)=.true.
          if (t-track%period(n).gt.max_dt.and..not.associated(ms_energy)) then
            track%period(n)=1.d22
            track%hxm0(3,n)=0.d0
          else
            track%period(n)=t-track%period(n)
          endif
          if (calc_density_p) then
            i=track%orbit(n)%n_orbit
            i=i+1; track%orbit(n)%orbit_param(1,i)=track%ddm(n)
          endif
          if (n.eq.TRAJ) print "('orbit#',i4,' done',1p,5(x,e12.5),0p,x,i4,x,l,x,l,x,l)",n,track%layer%energy(n), &
              1.d0-track%energy(n)/track%layer%energy(n),track%period(n)/dt, &
              sqrt(dot_product(track%dm_dt(:,n),track%dm_dt(:,n))),t,count(track%orbit_done(:)),track%layer%do_grain(n), &
              associated(ms_energy),(track%period(n).gt.max_dt)
        else
          go_p=.true.
          if (calc_density_p) then
            i=track%orbit(n)%n_orbit
            if (track%ddm(n)-track%orbit(n)%orbit_param(1,max(1,i)).gt.track%orbit(n)%min_dm) then
              i=i+1; track%orbit(n)%orbit_param(1,i)=track%ddm(n)
              track%orbit(n)%orbit_param(2,i)=sqrt(dot_product(dm,dm))/dt;
              track%orbit(n)%orbit_param(3,i)=dot_product( mh,mh)
              track%orbit(n)%orbit_m(:,i)=m0; track%orbit(n)%n_orbit=i
            endif
          endif
        endif
      endif
     endif
    enddo
    track=>track%next
  enddo

!if (associated(ms_lm)) then
!do n=1,size(trackk%layer%grain_vol)
!  if (trackk%layer%do_grain(n).and.ms_lm(n).le.0.d0.and.ms_lmi(n).ne.ms_lm(n)) then
!    track=>trackk
!    do while (associated(track))
!      track%hxm(:,n)=track%layer%m(1+track%layer%in:3+track%layer%in,n)
!      track=>track%next
!    enddo
!  endif
!enddo
!endif

  !multispin
! go_p=.false.; gg_p=.true.
! track=>trackk
! do while (associated(track))
!   do n=1,size(track%layer%grain_vol)
!     if (gg_p(n)) then
!       if (track%angle(n).gt.0.d0.and.track%period(n).lt.1.d20) then
!         print *,n,acos(dot_product(track%layer%m(1+track%layer%in:3+track%layer%in,n),track%start(:,n)))*180.d0/acos(-1.d0)
!         if (track%period(n).lt.1.d20) track%angle(n)=-track%angle(n)
!       endif
!     else
!       go_p=.true.; track%orbit_done(n)=.false.; track%period(n)=0.d0
!     endif
!   enddo
!   track=>track%next
! enddo
! if (go_p) go_p=(t.lt.max_dt)
 end function
end module track_m
