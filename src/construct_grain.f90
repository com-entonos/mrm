
!this module creates layers of grains' microstructure, Hk, Hx, Tc, etc.

module construct_grain_m
  use random_m, only: random_s
  use m_indmed, only: indmed
  implicit none
  private
  public set_construct_grain,construct_microstructure,construct_ms,construct_tc,construct_hk_magnitude, &
         construct_hk_direction,construct_hx_transverse,get_ms_nom,construct_hx_vertical,show_grain,get_k_nom,construct_tau

  type::grain_s
    !grain micostructure
    logical::grain=.false.,old_small_grain=.true.,grow_small_grain=.false.      !are grains collections of nodes?
    real(8)::final_packing_fraction=1.d0, packing_fraction=1.d0,small_packing_fraction=0.d0,initial_grain_area=0.d0,initial_grain_ext_area=0.d0, &
             small_grain_area=0.d0,grain_sigma=0.d0,small_grain_sigma=0.d0,average_grain_area=0.d0, dr_grain=0.2d0, initial_small_grain_area=0.d0
    !Hk
    real(8)::kx=0.d0,ky=0.d0,kz=0.d0,k1=0.d0,k2=0.d0   !nominal direction of Hk and magnitudes
    real(8)::k_sigma_small=0.d0,k1_small=0.d0,k2_small=0.d0,frac_hk_small=0.d0
    real(8)::k_ang=0.d0,k_cone_sigma=0.d0,kmax=1.d300,k_frac_inplane=0.d0,k_frac_random=0.d0,k_sigma=0.d0,k_ac=0.37d-7, hkmin=-1.d22
    logical::k_3d=.false.,k_2d=.false.,k_uniform=.false.,k_cone=.false.,k_cubic=.false.,k_scale=.false.,k_rescale=.false.,k_scale_al=.false.
    integer::num_end=2
    !Tc
    real(8)::tc_nom=1.d22,tc_sigma=0.d0,tc_max=1.d22,tc_d0=0.9d0,tc_eta=0.d0,tau_nom=1.d-9,tau_sigma=0.d0,vol_scale=1.d0
    logical::tc_rescale=.false.
    !Ms
    real(8)::ms_nom=0.d0,ms_sigma=0.d0,ms_max=1.d300
    !hx_lateral
    real(8)::ex=0.d0,ex_sigma=0.d0,ex_upper=0.d0,ex_upper_sigma=0.d0
    !hx_vertical
    !internal
    type(random_s),pointer::random_grain=>null(),random_grains=>null(),random_hk_mag=>null(),random_hk_uni_angles=>null(), &
           random_hk_frac_random=>null(), random_hk_frac_inplane=>null(),random_hk_cubic_angles=>null(),random_tc=>null(), &
           random_ms=>null(),random_ex=>null(),random_exv=>null(),random_tau=>null(),random_hk_mag_small=>null()
  end type grain_s

  type::grain_ss
    type(grain_s),pointer::p
  end type grain_ss

  type(grain_ss),pointer,save::grain(:)=>null()
  integer,save::dist_io=0
contains

  subroutine show_grain(i)
   use io_m, only:output
   integer::i
   type(grain_s),pointer::p
   character(200)::os

     p=>grain(i)%p
     write(os,"('  grains are nodes=',l1)") (.not.p%grain);call output(os)
     write(os,"(1p,'  Ms=',e12.5,', sigma Ms=',e12.5,', max Ms=',e12.5)") p%ms_nom,p%ms_sigma,p%ms_max; call output(os)
     write(os,"(1p,'  Tc=',e12.5,', sigma Tc=',e12.5,', max Tc=',e12.5)") p%tc_nom,p%tc_sigma,p%tc_max; call output(os)
     write(os,"(1p,'  Tc d0=',e12.5,', Tc eta=',e12.5,', rescale Tc=',0p,l1)") p%tc_d0,p%tc_eta,p%tc_rescale; call output(os)
     write(os,"(1p,'  <tau>=',e12.5,', sigma tau=',e12.5)") p%tau_nom,p%tau_sigma; call output(os)
     write(os,"(1p,'  Hk angle=',e12.5,', Hk cone=',e12.5)") p%k_ang,p%k_cone_sigma; call output(os)
     write(os,"(1p,'  Hk sigma=',e12.5,', Hk_1=',e12.5,', Hk_2=',e12.5)") p%k_sigma,p%k1,p%k2; call output(os)
     write(os,"(1p,'  vol_scale=',e12.5,0p,' num_end=',i0)") p%vol_scale,p%num_end; call output(os)
     write(os,"(1p,'  small Hk frac=',e12.5,', Hk sigma=',e12.5,', Hk_1=',e12.5,', Hk_2=',e12.5)") p%frac_hk_small,p%k_sigma_small, &
       p%k1_small,p%k2_small; call output(os)
     write(os,"(1p,'  Hkx=',e12.5,', Hky=',e12.5,', Hkz=',e12.5)") p%kx,p%ky,p%kz; call output(os)
     write(os,"('  Hk 3d random=',l1,', Hk 2d random=',l1,', Hk uniform=',l1,', Hk cone=',l1)") &
       p%k_3d,p%k_2d,p%k_uniform,p%k_cone; call output(os)
     write(os,"(1p,'  frac_inplane=',e12.5,', frac_random=',e12.5,', kmax=',e12.5,', hk_min=',e12.5,', k_ac=',e12.5)") p%k_frac_inplane, &
       p%k_frac_random,p%kmax,p%hkmin,p%k_ac; call output(os)
     write(os,"('  uniaxial=',l1,', cubic=',l1,', scale Hk=',l1,', A.L.''s scale Hk=',l1,', rescale Hk=',l11)") (.not.p%k_cubic),p%k_cubic, &
       p%k_scale,p%k_scale_al,p%k_rescale; call output(os)
     write(os,"(1p,'  intralayer Ex=',e12.5,', intralayer Ex sigma=',e12.5)") p%ex,p%ex_sigma; call output(os)
     if (i.gt.1) then
       write(os,"(1p,'  to layer above, interlayer Ex=',e12.5,', interlayer Ex sigma=',e12.5)") p%ex_upper,p%ex_upper_sigma; call output(os)
     endif
     write(os,"(1p,'  packing fraction=',e12.5,', small grain packing fraction=',e12.5)") p%packing_fraction,p%small_packing_fraction; call output(os)
     write(os,"(1p,'  final packing fraction=',e12.5)") p%final_packing_fraction; call output(os)
     write(os,"(1p,'  average grain area=',e12.5,', small grain area=',e12.5)") p%average_grain_area,p%small_grain_area; call output(os)
     write(os,"(1p,'  init grain area sigma=',e12.5,', small grain area sigma=',e12.5)") p%grain_sigma,p%small_grain_sigma; call output(os)
     write(os,"(1p,'  grow speed=',e12.5,', initial grain area=',e12.5)") p%dr_grain,p%initial_grain_area; call output(os)

  end subroutine

  function get_k_nom(il) result(k)
   integer,intent(in)::il
   real(8)::k(3)
   k=(/0.d0, 0.d0, 1.d0/)
   if (il.gt.0.and.il.le.size(grain)) k=(/grain(il)%p%kx,grain(il)%p%ky,grain(il)%p%kz/)
  end function
  function get_ms_nom(il) result(ms)
   integer,intent(in)::il
   real(8)::ms
   ms=0.d0
   if (il.gt.0.and.il.le.size(grain)) ms=grain(il)%p%ms_nom
  end function

  function construct_microstructure(nx,ny,dx,dy,ilayer,isilent_p,num_grain,grid2grain,grain2grid,grain2gridindex,grain_area) result (ok_p)

   use io_m, only: output
   use random_m, only: random_uniform, random_normal
   use transform_m, only: rotate
   use plot_m, only: get_plot, plot_histogram

   integer,intent(in)::nx,ny
   real(8),intent(in)::dx,dy
   integer,intent(in)::ilayer
   logical,intent(in)::isilent_p

   integer,intent(out)::num_grain
   integer,pointer::grid2grain(:,:),grain2grid(:,:),grain2gridindex(:)
   real(8),pointer::grain_area(:)

   logical::silent_p,ok_p,look_p
   integer::layer
   type(grain_s),pointer::p
   integer::i,j,i1,j1,num_small_grain,ng
   real(8)::pi,pf,pff,pfs,r0,x,y,small_grain_list(3,300),cdfm,r0max
   integer::grid1(nx*ny),grid(nx,ny),gridfree(nx,ny),num_free,n
   character(200)::os
   real(8),pointer::pos_seed(:,:),rad_seed(:)
   integer,allocatable::sort_grain(:)
   integer::nn(50),i0,j0
   real(8)::pp(50)
   real(8),pointer::pitch(:)=>null(),pitchb(:)=>null(),tp(:)=>null()

   integer::current_ns,current_n,n1,n2
   real(8)::rs,small_grain_sigma,cdfms,dn,fs

   ok_p=.false.; num_grain=-1
!  silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
!  layer=1;if (present(ilayer)) layer=max(1,ilayer)
   silent_p=isilent_p
   layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   if (associated(p)) then

     pi=acos(-1.d0)
     pf=p%packing_fraction; if (pf.gt.1.d0) pf=pf/1.d2
     pff=p%final_packing_fraction; if (pff.gt.1.d0) pff=pff/1.d2
     pfs=p%small_packing_fraction; if (pfs.gt.1.d0) pfs=pfs/1.d2

     if (.not.p%grain) then  !grains are nodes
       pf=min(1.d0,max(0.d0,pf+pfs))
       num_small_grain=0
       num_grain=nint(nx*ny*pf)
       num_grain=int(nx*ny*pf)
     else                    !grains are collections of nodes
       pf=max(0.d0, min(1.d0, pf-pfs))
       num_small_grain=0; if (p%small_grain_area.ne.0.d0) num_small_grain=nint(pfs*nx*ny*dx*dy/p%small_grain_area)
       num_grain=nint(pf*nx*ny*dx*dy/p%average_grain_area)+num_small_grain
       if (p%initial_small_grain_area.le.0.d0) p%initial_small_grain_area=p%small_grain_area
     endif

     if (.not.silent_p) then
       write(os,"(' constructing media microstructure for layer #',0p,i3)")layer
       call output(trim(os))
       write(os,"('  constructing ',i7' grains (',i5,' small) for layer #',i3)")num_grain,num_small_grain,layer
       call output(trim(os))
     endif


     if (.not.p%grain) then  !grains are nodes
       write(os,"('   grains are equal to elements for layer #',0p,i3)")layer; if (.not.silent_p) call output(trim(os))
       i=0; grid1=1
       do while (i.lt.nx*ny-num_grain)
         j=min(nx*ny,max(1,1+int(nx*ny*random_uniform(p%random_grain))))
         if (grid1(j).gt.0) then; i=i+1; grid1(j)=0;endif
       enddo
       if (num_grain.ne.count(grid1.gt.0)) then
               print *,'error- num_grain != occupied grid elements!'
               stop
       endif
       if (associated(grid2grain)) then; if (size(grid2grain,1).ne.nx.or.size(grid2grain,2).ne.ny)then; deallocate(grid2grain); nullify(grid2grain);endif; endif
       if (associated(grain2grid)) then; if (size(grain2grid,2).ne.num_grain)then; deallocate(grain2grid); nullify(grain2grid);endif; endif
       if (associated(grain2gridindex)) then; if (size(grain2gridindex).ne.num_grain+1) then; deallocate(grain2grid); nullify(grain2grid);endif; endif
       if (associated(grain_area)) then; if (size(grain_area).ne.num_grain) then; deallocate(grain_area); nullify(grain_area);endif; endif
       if (.not.associated(grain2gridindex)) allocate(grain2gridindex(0:num_grain));
       if (.not.associated(grid2grain)) allocate(grid2grain(nx,ny))
       if (.not.associated(grain2grid)) allocate(grain2grid(2,num_grain));
       if (.not.associated(grain_area)) allocate(grain_area(num_grain)); grain_area=dx*dy
!if (lbound(grain2gridindex,dim=1).ne.0) then; print *,'construct_microstructure';read *;endif


       i=0; grid2grain=0; grain2gridindex=0
       do j=1,nx*ny
         if (grid1(j).gt.0) then !grain!
            i=i+1; j1=(j-1)/nx+1; i1=j-nx*(j1-1)
            grid2grain(i1,j1)=i
            grain2grid(:,i)=(/i1,j1/); grain2gridindex(i)=i
         endif
       enddo
     else
       allocate(pos_seed(2,num_grain),rad_seed(num_grain)); r0=cdf_erf(0.d0,0.d0,0.d0,.true.)
       r0=max(1.d0,sqrt(p%initial_grain_area/(pi*dx*dy))); r0max=-1.d300
       cdfm=0.5d0*(p%grain_sigma*p%grain_sigma+2.d0*log(r0))
       grid=0; gridfree=0; num_free=nx*ny
       if (p%old_small_grain) then
       if (.not.silent_p) then
         write(os,"('   trying to place ',i5,' grains with initial area',1p,e12.5,' nodes (',e12.5,' cm^2) and sigma ',e12.5,', for layer #',0p,i3)") &
           num_grain-num_small_grain,p%initial_grain_area/dx/dy,p%initial_grain_area,p%grain_sigma,layer
         call output(trim(os))
         write(os,"('    average grain will have',1p,e12.5,' node (',1p,e12.5,' cm^2) for layer #',0p,i3)")p%average_grain_area/dx/dy,p%average_grain_area,layer
         call output(trim(os))
       endif
       do n=num_small_grain+1,num_grain

         if (p%grain_sigma.ne.0.d0) then
           rad_seed(n)=max(0.75d0,cdf_erf((num_grain-n+0.5d0)/(num_grain-num_small_grain),cdfm,p%grain_sigma))
         else
           rad_seed(n)=max(0.75d0,r0)
         endif
         gridfree=0
         do i=num_small_grain+1,n-1
           x=pos_seed(1,i); y=pos_seed(2,i)
           do j1=int(y-rad_seed(i)-rad_seed(n)-1),int(y+rad_seed(i)+rad_seed(n)+1)
             do i1=int(x-rad_seed(i)-rad_seed(n)-1),int(x+rad_seed(i)+rad_seed(n)+1)
               if (radius(i1-x-0.5d0,j1-y-0.5d0).lt.(rad_seed(i)+rad_seed(n))*(rad_seed(i)+rad_seed(n))) gridfree(iwrpx(i1),iwrpy(j1))=i
             enddo
           enddo
         enddo
         num_free=count(gridfree.eq.0)
         do while (num_free.eq.0)
           r0=r0-0.05d0; gridfree=0
           cdfm=0.5d0*(p%grain_sigma*p%grain_sigma+2.d0*log(r0))
!   print *,n,r0,rad_seed(n),(num_grain-n+0.5d0)/(num_grain-num_small_grain),exp(cdfm),p%grain_sigma
           if (p%grain_sigma.ne.0.d0) then
             rad_seed(n)=max(0.75d0,cdf_erf((num_grain-n+0.5d0)/(num_grain-num_small_grain),cdfm,p%grain_sigma))
           else
             rad_seed(n)=max(0.75d0,r0)
           endif
           do i=num_small_grain+1,n-1
             x=pos_seed(1,i); y=pos_seed(2,i)
             do j1=int(y-rad_seed(i)-rad_seed(n)-1),int(y+rad_seed(i)+rad_seed(n)+1)
               do i1=int(x-rad_seed(i)-rad_seed(n)-1),int(x+rad_seed(i)+rad_seed(n)+1)
                 if (radius(i1-x-0.5d0,j1-y-0.5d0).lt.(rad_seed(i)+rad_seed(n))*(rad_seed(i)+rad_seed(n))) gridfree(iwrpx(i1),iwrpy(j1))=i
               enddo
             enddo
           enddo
           num_free=count(gridfree.eq.0)
!   print *,n,r0,num_free
         enddo
         r0max=max(r0max,rad_seed(n))
         if (.not.silent_p.and.mod(num_grain-n,(num_grain-num_small_grain)/10).eq.0) then
           write(os,"('     trying to place grain #',i5,' with radius',1p,e12.5,' nodes')")n-num_small_grain+1,rad_seed(n)
           call output(trim(os))
         endif

         i1=min(num_free,max(1,1+int(num_free*random_uniform(p%random_grain))))
         do j=1,ny
           do i=1,nx
             if (gridfree(i,j).eq.0) i1=i1-1
             if (i1.eq.0) exit
           enddo
           if (i1.eq.0) exit
         enddo
         x=min(1.d0,max(0.d0,random_uniform(p%random_grain)))+i-1
         y=min(1.d0,max(0.d0,random_uniform(p%random_grain)))+j-1
         pos_seed(:,n)=(/x,y/)
         do j1=int(y-rad_seed(n)-1),int(y+rad_seed(n)+1)
           do i1=int(x-rad_seed(n)-1),int(x+rad_seed(n)+1)
             if (radius(i1-x-0.5d0,j1-y-0.5d0).le.rad_seed(n)*rad_seed(n)) grid(iwrpx(i1),iwrpy(j1))=n
           enddo
         enddo
!   print *,n,count(grid.eq.n)
       enddo
       write(os,"('   trying to place ',i5,' small grains with area ',1p,e12.5,' nodes (',e12.5,' cm^2) and sigma ',e12.5,',for layer #',0p,i3)") &
             num_small_grain,p%small_grain_area/dx/dy,p%small_grain_area,p%small_grain_sigma,layer
       if (.not.silent_p.and.num_small_grain.ne.0) call output(trim(os))
       num_free=count(grid.eq.0); gridfree=grid
       r0=sqrt(p%initial_small_grain_area/(pi*dx*dy))
       do n=1,num_small_grain
         look_p=.true.
         rad_seed(n)=max(0.75d0,random_normal(p%random_grains)*max(0.d0,p%small_grain_sigma/2.d0/dx)+r0)
         do while (look_p)
           i1=min(num_free,max(1,1+int(num_free*random_uniform(p%random_grains))))
!print *,i1,num_free
           do j=1,ny
             do i=1,nx
               if (gridfree(i,j).eq.0) i1=i1-1
               if (i1.eq.0) exit
             enddo
             if (i1.eq.0) exit
           enddo
           x=min(1.d0,max(0.d0,random_uniform(p%random_grains)))+i-1
           y=min(1.d0,max(0.d0,random_uniform(p%random_grains)))+j-1
           ng=0
           do j1=int(y-rad_seed(n)-1),int(y+rad_seed(n)+1)
             do i1=int(x-rad_seed(n)-1),int(x+rad_seed(n)+1)
               do i=1,ng
                 if (radius(i1-x-0.5d0,j1-y-0.5d0).lt.small_grain_list(3,i)) exit
               enddo
               do j=ng,i,-1
                 small_grain_list(:,j+1)=small_grain_list(:,j)
               enddo
               small_grain_list(:,i)=(/ dble(i1), dble(j1), radius(i1-x-0.5d0,j1-y-0.5d0) /)
               ng=ng+1
             enddo
           enddo
           do i=1,nint(pi*rad_seed(n)*rad_seed(n))
             look_p=(look_p.and.grid(iwrpx(nint(small_grain_list(1,i))),iwrpy(nint(small_grain_list(2,i)))).eq.0)
           enddo
           if (look_p) then
!    print *,n,nint(pi*rad_seed(n)*rad_seed(n)),r0,rad_seed(n)
             do i=1,nint(pi*rad_seed(n)*rad_seed(n))
               grid(iwrpx(nint(small_grain_list(1,i))),iwrpy(nint(small_grain_list(2,i))))=n
               num_free=num_free-1
             enddo
             pos_seed(:,n)=(/x,y/)
           endif
           look_p=(.not.look_p)
         enddo
       enddo
       n2=-1
       else ! new way of doing small grains
         if (.not.silent_p) then
           write(os,"('   trying to place ',i5,' ln-normal grains with initial area',1p,e12.5,' nodes (',e12.5,' cm^2) and sigma ',e12.5,', for layer #',0p,i3)") &
             num_grain-num_small_grain,p%initial_grain_area/dx/dy,p%initial_grain_area,p%grain_sigma,layer
           call output(trim(os))
           write(os,"('    average ln-normal grain will have',1p,e12.5,' node (',1p,e12.5,' cm^2) for layer #',0p,i3)")p%average_grain_area/dx/dy,p%average_grain_area,layer
           call output(trim(os))
           write(os,"('   trying to place ',i5,' normal grains with initial area',1p,e12.5,' nodes (',e12.5,' cm^2) and sigma ',e12.5,', for layer #',0p,i3)") &
             num_small_grain,p%initial_small_grain_area/dx/dy,p%initial_small_grain_area,p%small_grain_sigma/sqrt(dx*dy),layer
           call output(trim(os))
           write(os,"('    average normal grain will have',1p,e12.5,' node (',1p,e12.5,' cm^2) for layer #',0p,i3)")p%small_grain_area/dx/dy,p%small_grain_area,layer
           call output(trim(os))
         endif
         r0=cdf_serf(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,.true.)
!        r0=max(1.d0,sqrt(p%initial_grain_area/(pi*dx*dy))); r0max=-1.d300
!        cdfm=0.5d0*(p%grain_sigma*p%grain_sigma+2.d0*log(r0))
         r0max=-1.d300
         r0=sqrt(p%initial_grain_area/(pi*dy*dx))
         cdfm=log(r0)-0.5d0*p%grain_sigma*p%grain_sigma
         cdfms=sqrt(p%initial_small_grain_area/(dx*dy*pi)); small_grain_sigma=p%small_grain_sigma*0.5d0/sqrt(dx*dy)
         grid=0; gridfree=0; num_free=nx*ny
         current_ns=0;current_n=num_small_grain
         fs=dble(num_small_grain)/num_grain
!print *,fs,pfs/(pf+pfs)
!print *,pfs/(pf+pfs),cdfm,p%grain_sigma
!print *,cdfms,small_grain_sigma
         n2=num_grain; dn=0.d0
         do n1=1,num_grain
           n2=n2-1; dn=max(0.d0,dn-1.d0)

           if (p%grain_sigma.ne.0.d0.or.small_grain_sigma.ne.0.d0) then
             rs=cdf_serf(max(0.1d0,n2-dn+0.5d0)/num_grain,fs,cdfm,p%grain_sigma,cdfms,small_grain_sigma,random_uniform(p%random_grains))
!            rs=sign(min(10.d0,max(0.75d0,abs(rs))),rs)
             rs=sign(max(0.75d0,abs(rs)),rs)
           else
             rs=max(0.75d0,r0)
           endif
           gridfree=0
           do i=1,num_grain
             if (i.le.current_ns.or.(i.gt.num_small_grain.and.i.le.current_n)) then
               x=pos_seed(1,i); y=pos_seed(2,i)
               do j1=int(y-rad_seed(i)-abs(rs)-1),int(y+rad_seed(i)+abs(rs)+1)
                 do i1=int(x-rad_seed(i)-abs(rs)-1),int(x+rad_seed(i)+abs(rs)+1)
                   if (radius(i1-x-0.5d0,j1-y-0.5d0).lt.(rad_seed(i)+abs(rs))*(rad_seed(i)+abs(rs))) gridfree(iwrpx(i1),iwrpy(j1))=i
                 enddo
               enddo
             endif
           enddo
           num_free=count(gridfree.eq.0)
!if (num_free.eq.0) then
           do while (num_free.eq.0)
!print *,'shit'
!            r0=r0-0.05d0; gridfree=0
!            n2=n2-1
             gridfree=0
             dn=dn+0.2d0
             if (n2.lt.0) r0=r0-0.05d0
             cdfm=0.5d0*(p%grain_sigma*p%grain_sigma+2.d0*log(r0))
!   print *,n,r0,rad_seed(n),(num_grain-n+0.5d0)/(num_grain-num_small_grain),exp(cdfm),p%grain_sigma
             if (p%grain_sigma.ne.0.d0.or.small_grain_sigma.ne.0.d0) then
               rs=cdf_serf(max(0.1d0,n2-dn+0.5d0)/num_grain,fs,cdfm,p%grain_sigma,cdfms,small_grain_sigma,random_uniform(p%random_grains))
!              rs=sign(min(10.d0,max(0.75d0,abs(rs))),rs)
               rs=sign(max(0.75d0,abs(rs)),rs)
             else
               rs=max(0.75d0,r0)
             endif
             do i=1,num_grain
               if (i.le.current_ns.or.(i.gt.num_small_grain.and.i.le.current_n)) then
                 x=pos_seed(1,i); y=pos_seed(2,i)
                 do j1=int(y-rad_seed(i)-abs(rs)-1),int(y+rad_seed(i)+abs(rs)+1)
                   do i1=int(x-rad_seed(i)-abs(rs)-1),int(x+rad_seed(i)+abs(rs)+1)
                     if (radius(i1-x-0.5d0,j1-y-0.5d0).lt.(rad_seed(i)+abs(rs))*(rad_seed(i)+abs(rs))) gridfree(iwrpx(i1),iwrpy(j1))=i
                   enddo
                 enddo
               endif
             enddo
             num_free=count(gridfree.eq.0)
!     print *,n,r0,num_free
           enddo
!print *,n2,dn
!print *,cdfm,cdfms
!print *,max(0.1d0,n2-dn+0.5d0)/num_grain,fs
!endif
!print *,rs
           if (current_n.ge.num_grain.or.(rs.lt.0.d0.and.current_ns.lt.num_small_grain)) then
             current_ns=current_ns+1; n=current_ns
           else
             current_n=current_n+1; n=current_n
           endif
           rad_seed(n)=abs(rs)
!print *,current_ns,current_n,rs
           
           r0max=max(r0max,rad_seed(n))
           if (.not.silent_p.and.mod(n1,num_grain/15).eq.0) then
             write(os,"('     trying to place grain #',i5,' with radius',1p,e12.5,' nodes')")n1,rs
             call output(trim(os))
           endif

           i1=min(num_free,max(1,1+int(num_free*random_uniform(p%random_grain))))
           do j=1,ny
             do i=1,nx
               if (gridfree(i,j).eq.0) i1=i1-1
               if (i1.eq.0) exit
             enddo
             if (i1.eq.0) exit
           enddo
           x=min(1.d0,max(0.d0,random_uniform(p%random_grain)))+i-1
           y=min(1.d0,max(0.d0,random_uniform(p%random_grain)))+j-1
           pos_seed(:,n)=(/x,y/)
           do j1=int(y-rad_seed(n)-1),int(y+rad_seed(n)+1)
             do i1=int(x-rad_seed(n)-1),int(x+rad_seed(n)+1)
               if (radius(i1-x-0.5d0,j1-y-0.5d0).le.rad_seed(n)*rad_seed(n)) grid(iwrpx(i1),iwrpy(j1))=n
             enddo
           enddo
!   print *,n,count(grid.eq.n)
         enddo
         n=num_small_grain+1; if (p%grow_small_grain) n=1
         n2=1
       endif ! new way of doing small grains

       if (.not.silent_p) call output('    growing grains')
       num_free=(count(grid.eq.0))
       r0=0.d0; n=num_small_grain+1
!      do while (num_free.gt.(1.d0-pfs-pf)*nx*ny)
       do while (num_free.gt.(1.d0-pff)*nx*ny)
         do j1=int(pos_seed(2,n)-r0-rad_seed(n)-1),int(pos_seed(2,n)+r0+rad_seed(n)+1)
           do i1=int(pos_seed(1,n)-r0-rad_seed(n)-1),int(pos_seed(1,n)+r0+rad_seed(n)+1)
             if (radius(i1-pos_seed(1,n)-0.5d0,j1-pos_seed(2,n)-0.5d0).le.(r0+rad_seed(n))*(r0+rad_seed(n)).and.grid(iwrpx(i1),iwrpy(j1)).eq.0.and. &
                 (grid(iwrpx(i1),iwrpy(j1-1)).eq.n.or.grid(iwrpx(i1-1),iwrpy(j1)).eq.n.or. &
                  grid(iwrpx(i1+1),iwrpy(j1)).eq.n.or.grid(iwrpx(i1),iwrpy(j1+1)).eq.n)) then
               grid(iwrpx(i1),iwrpy(j1))=n; num_free=num_free-1
!              if (num_free.lt.(1.d0-pfs-pf)*nx*ny) exit
               if (num_free.lt.(1.d0-pff)*nx*ny) exit
             endif
           enddo
!          if (num_free.lt.(1.d0-pfs-pf)*nx*ny) exit
           if (num_free.lt.(1.d0-pff)*nx*ny) exit
         enddo
         n=n+n2 !n=n-1
         if ((n.le.num_small_grain.and..not.p%grow_small_grain).or.n.lt.1.or.n.gt.num_grain) then
           n=num_grain
           if (n2.gt.0) then; n=num_small_grain+1; if (p%grow_small_grain) n=1; endif
           r0=r0+p%dr_grain
           write(os,"('     additional radius is',1p,e12.5,' nodes, with',0p,i6,' nodes open')")r0,num_free
           if (.not.silent_p) call output(trim(os))
         endif
       enddo
!      r0max=2.d0*(r0max+r0)

!if (associated(grain_area)) then; if (size(grain_area).ne.num_grain) then; deallocate(grain_area); nullify(grain_area);endif; endif
!allocate(grain_area(num_grain))
!grain_area=0.d0
!do j1=1,ny; do i1=1,nx
!if (grid(i1,j1).ne.0) grain_area(grid(i1,j1))=grain_area(grid(i1,j1))+1
!enddo;enddo
!write(os,"('small grain diameter for layer #',0p,i3)")layer
!if (num_small_grain.gt.0) call plot_histogram(trim(os),50,sqrt(grain_area(1:num_small_grain)*dx*dy/pi)*2.d7)
!write(os,"('big grain diameter for layer #',0p,i3)")layer
!call plot_histogram(trim(os),50,sqrt(grain_area(num_small_grain+1:num_grain)*dx*dy/pi)*2.d7)
!write(os,"('all grain diameter for layer #',0p,i3)")layer
!call plot_histogram(trim(os),50,sqrt(grain_area*dx*dy/pi)*2.d7)
!deallocate(grain_area);nullify(grain_area)
!read *
!      do n=1,num_grain
!        print *,n,count(grid.eq.n)
!        if (count(grid.eq.n).eq.0) pause
!      enddo
       write(os,"('   computing pitch for layer #',i3)") layer
       if (.not.silent_p) call output(trim(os))
       do n=1,num_grain   !compute centroid of each grain
         x=0.d0; y=0.d0; j=0
         do j1=int(pos_seed(2,n)-r0max-1),int(pos_seed(2,n)+r0max+1)
           do i1=int(pos_seed(1,n)-r0max-1),int(pos_seed(1,n)+r0max+1)
             if (grid(iwrpx(i1),iwrpy(j1)).eq.n) then
               x=x+i1-0.5d0; y=y+j1-0.5d0; j=j+1
             endif
           enddo
         enddo
         pos_seed(:,n)=(/ x/j, y/j /)
       enddo; gridfree=grid


      !find pitch w/ all grains
       num_free=(count(grid.eq.0)); r0=0.d0; n=num_grain
       do while (num_free.gt.0)    !first grow grains to cover the plane
         do j1=int(pos_seed(2,n)-r0-rad_seed(n)-1),int(pos_seed(2,n)+r0+rad_seed(n)+1)
           do i1=int(pos_seed(1,n)-r0-rad_seed(n)-1),int(pos_seed(1,n)+r0+rad_seed(n)+1)
             if (radius(i1-pos_seed(1,n)-0.5d0,j1-pos_seed(2,n)-0.5d0).le.(r0+rad_seed(n))*(r0+rad_seed(n)).and.grid(iwrpx(i1),iwrpy(j1)).eq.0.and. &
                 (grid(iwrpx(i1),iwrpy(j1-1)).eq.n.or.grid(iwrpx(i1-1),iwrpy(j1)).eq.n.or. &
                  grid(iwrpx(i1+1),iwrpy(j1)).eq.n.or.grid(iwrpx(i1),iwrpy(j1+1)).eq.n)) then
               grid(iwrpx(i1),iwrpy(j1))=n; num_free=num_free-1
               if (num_free.lt.1) exit
             endif
           enddo
           if (num_free.lt.1) exit
         enddo
         n=n-1
         if (n.lt.1) then
           n=num_grain; r0=r0+p%dr_grain
         endif
       enddo
       r0=2.d0*(r0max+r0)
!      allocate(pp(num_grain,num_grain)); pp=-1.d0
       allocate(pitch(num_grain*6)); pitch=0.d0; j0=0
       do n=1,num_grain          !now find the grains which touch (i.e. nearest neighbors) and store the distance^2 between centroids
         nn=0; pp=-1.d0; j=0
         do j1=int(pos_seed(2,n)-r0-1),int(pos_seed(2,n)+r0+1)
           do i1=int(pos_seed(1,n)-r0-1),int(pos_seed(1,n)+r0+1)
             i=grid(iwrpx(i1),iwrpy(j1))
!            if (i.ne.0.and.i.ne.n.and. &
             if (i.gt.n.and. &
                 (grid(iwrpx(i1),iwrpy(j1-1)).eq.n.or.grid(iwrpx(i1-1),iwrpy(j1)).eq.n.or. &
                  grid(iwrpx(i1+1),iwrpy(j1)).eq.n.or.grid(iwrpx(i1),iwrpy(j1+1)).eq.n)) then
!                   pp(min(i,n),max(i,n))=radius(dx*mod(abs(pos_seed(1,n)-pos_seed(1,i)),nx*0.5d0),dy*mod(abs(pos_seed(2,n)-pos_seed(2,i)),ny*0.5d0))
!                   pp(min(i,n),max(i,n))=radius(dx*min(abs(pos_seed(1,n)-pos_seed(1,i)+nx),abs(pos_seed(1,n)-pos_seed(1,i)-nx),abs(pos_seed(1,n)-pos_seed(1,i))),dy*min(abs(pos_seed(2,n)-pos_seed(2,i)+ny),abs(pos_seed(2,n)-pos_seed(2,i)-ny),abs(pos_seed(2,n)-pos_seed(2,i))))
                    do i0=1,j; if (nn(i0).eq.i) exit; enddo; j=max(j,i0); nn(i0)=i
                    pp(i0)=radius(dx*min(abs(pos_seed(1,n)-pos_seed(1,i)+nx),abs(pos_seed(1,n)-pos_seed(1,i)-nx),abs(pos_seed(1,n)-pos_seed(1,i))),dy*min(abs(pos_seed(2,n)-pos_seed(2,i)+ny),abs(pos_seed(2,n)-pos_seed(2,i)-ny),abs(pos_seed(2,n)-pos_seed(2,i))))
             endif
           enddo
         enddo
         if (j0+j.gt.size(pitch)) then
           allocate(tp(2*size(pitch))); tp(1:size(pitch))=pitch
           deallocate(pitch); pitch=>tp; nullify(tp)
         endif
         pitch(j0+1:j0+j)=pp(1:j); j0=j0+j
       enddo;
       allocate(tp(j0)); tp=pitch(1:j0); deallocate(pitch); pitch=>tp; nullify(tp)
!      allocate(pitch(count(pp.gt.0.d0)));j=0
!      do n=2,num_grain
!        do i=1,n-1
!          if (pp(i,n).gt.0.d0) then
!            j=j+1; pitch(j)=pp(i,n)
! if (1.d7*sqrt(pp(i,n)).gt.60.d0) print *,i,n,pos_seed(:,i),pos_seed(:,n),sqrt(pp(i,n)/dx/dy)
!          endif
!        enddo
!      enddo
       pitch=sqrt(pitch)
      !find pitch w/ only big grains
       if (num_small_grain.gt.0) then
         grid=gridfree; where (grid.le.num_small_grain) grid=0
         num_free=(count(grid.eq.0)); r0=0.d0; n=num_grain
         do while (num_free.gt.0)    !first grow grains to cover the plane
           do j1=int(pos_seed(2,n)-r0-rad_seed(n)-1),int(pos_seed(2,n)+r0+rad_seed(n)+1)
             do i1=int(pos_seed(1,n)-r0-rad_seed(n)-1),int(pos_seed(1,n)+r0+rad_seed(n)+1)
               if (radius(i1-pos_seed(1,n)-0.5d0,j1-pos_seed(2,n)-0.5d0).le.(r0+rad_seed(n))*(r0+rad_seed(n)).and.grid(iwrpx(i1),iwrpy(j1)).eq.0.and. &
                   (grid(iwrpx(i1),iwrpy(j1-1)).eq.n.or.grid(iwrpx(i1-1),iwrpy(j1)).eq.n.or. &
                    grid(iwrpx(i1+1),iwrpy(j1)).eq.n.or.grid(iwrpx(i1),iwrpy(j1+1)).eq.n)) then
                 grid(iwrpx(i1),iwrpy(j1))=n; num_free=num_free-1
                 if (num_free.lt.1) exit
               endif
             enddo
             if (num_free.lt.1) exit
           enddo
           n=n-1
           if (n.le.num_small_grain) then
             n=num_grain; r0=r0+p%dr_grain
           endif
         enddo
         r0=2.d0*(r0max+r0)
!        pp=-1.d0
         allocate(pitchb(num_grain*6)); pitchb=0.d0; j0=0
         do n=num_small_grain+1,num_grain          !now find the grains which touch (i.e. nearest neighbors) and store the distance^2 between centroids
           nn=0; pp=-1.d0; j=0
           do j1=int(pos_seed(2,n)-r0-1),int(pos_seed(2,n)+r0+1)
             do i1=int(pos_seed(1,n)-r0-1),int(pos_seed(1,n)+r0+1)
               i=grid(iwrpx(i1),iwrpy(j1))
!              if (i.ne.0.and.i.ne.n.and. &
               if (i.gt.n.and. &
                   (grid(iwrpx(i1),iwrpy(j1-1)).eq.n.or.grid(iwrpx(i1-1),iwrpy(j1)).eq.n.or. &
                    grid(iwrpx(i1+1),iwrpy(j1)).eq.n.or.grid(iwrpx(i1),iwrpy(j1+1)).eq.n)) then
!                     pp(min(i,n),max(i,n))=radius(dx*(pos_seed(1,n)-pos_seed(1,i)),dy*(pos_seed(2,n)-pos_seed(2,i)))
!                     pp(min(i,n),max(i,n))=radius(dx*min(abs(pos_seed(1,n)-pos_seed(1,i)+nx),abs(pos_seed(1,n)-pos_seed(1,i)-nx),abs(pos_seed(1,n)-pos_seed(1,i))),dy*min(abs(pos_seed(2,n)-pos_seed(2,i)+ny),abs(pos_seed(2,n)-pos_seed(2,i)-ny),abs(pos_seed(2,n)-pos_seed(2,i))))
                      do i0=1,j; if (nn(i0).eq.i) exit; enddo; j=max(j,i0); nn(i0)=i
                      pp(i0)=radius(dx*min(abs(pos_seed(1,n)-pos_seed(1,i)+nx),abs(pos_seed(1,n)-pos_seed(1,i)-nx),abs(pos_seed(1,n)-pos_seed(1,i))),dy*min(abs(pos_seed(2,n)-pos_seed(2,i)+ny),abs(pos_seed(2,n)-pos_seed(2,i)-ny),abs(pos_seed(2,n)-pos_seed(2,i))))
               endif
             enddo
           enddo
           if (j0+j.gt.size(pitchb)) then
             allocate(tp(2*size(pitchb))); tp(1:size(pitchb))=pitchb
             deallocate(pitchb); pitchb=>tp; nullify(tp)
           endif
           pitchb(j0+1:j0+j)=pp(1:j); j0=j0+j
         enddo;
         allocate(tp(j0)); tp=pitchb(1:j0); deallocate(pitchb); pitchb=>tp; nullify(tp)
!        allocate(pitchb(count(pp.gt.0.d0)));j=0
!        do n=num_small_grain+2,num_grain
!          do i=num_small_grain+1,n-1
!            if (pp(i,n).gt.0.d0) then
!              j=j+1; pitchb(j)=pp(i,n)
!            endif
!          enddo
!        enddo
         pitchb=sqrt(pitchb)
       endif

       write(os,"('   constructing internal grain structures for layer #',i3)") layer
       if (.not.silent_p) call output(trim(os))
       grid=gridfree
       allocate(sort_grain(num_grain)); sort_grain=0
       n=0
       if (associated(grid2grain)) then; if (size(grid2grain,1).ne.nx.or.size(grid2grain,2).ne.ny)then; deallocate(grid2grain); nullify(grid2grain);endif; endif
       if (.not.associated(grid2grain)) allocate(grid2grain(nx,ny)); grid2grain=0
       do j=1,ny
         do i=1,nx
           if (grid(i,j).gt.0) then
             if (.not.any(sort_grain.eq.grid(i,j))) then
               n=n+1; sort_grain(n)=grid(i,j)
               where (grid.eq.grid(i,j)) grid2grain=n
!       print *,grid(i,j),count(grid.eq.grid(i,j)),n
             endif
           endif
         enddo
       enddo
       num_free=nx*ny-count(grid.eq.0)
!    print *,num_free,count(grid.ne.0),count(grid2grain.ne.0),num_grain

       if (associated(grain_area)) then; if (size(grain_area).ne.num_grain) then; deallocate(grain_area); nullify(grain_area);endif; endif
       if (associated(grain2grid)) then; if (size(grain2grid,2).ne.num_free)then; deallocate(grain2grid); nullify(grain2grid);endif; endif
       if (associated(grain2gridindex)) then; if (size(grain2gridindex).ne.num_grain+1) then; deallocate(grain2gridindex); nullify(grain2gridindex);endif; endif
       if (.not.associated(grain2gridindex)) allocate(grain2gridindex(0:num_grain));
       if (.not.associated(grain2grid)) allocate(grain2grid(2,num_free));
       if (.not.associated(grain_area)) allocate(grain_area(num_grain))
!if (lbound(grain2gridindex,dim=1).ne.0) then; print *,'construct_microstructure';read *;endif
       grain2grid=0; grain2gridindex=0; i1=0
       do n=1,num_grain
         do j=1,ny
           do i=1,nx
             if (grid2grain(i,j).eq.n) then
                i1=i1+1; grain2grid(:,i1)=(/i,j/)
             endif
           enddo
         enddo
         grain2gridindex(n)=i1
         if (grain2gridindex(n).eq.grain2gridindex(n-1)) then
            print *,'grain size of zero! ',n
            stop
         endif
       enddo
       if (i1.ne.num_free) then
           print *,'number elements not equal to num_free!'
           stop
       endif
!      do n=1,num_grain
!        print *,n,grain2gridindex(n)-grain2gridindex(n-1),count(grid2grain.eq.n)
!        if (grain2gridindex(n)-grain2gridindex(n-1).ne.count(grid2grain.eq.n)) pause
!        if (grain2gridindex(n)-grain2gridindex(n-1).eq.0) pause
!      enddo

       write(os,"('    packing fraction:',f6.2,'% for layer #',0p,i3)")1.d2*grain2gridindex(num_grain)/(nx*ny),layer
       if (.not.silent_p) call output(trim(os))

       if (.not.silent_p) then
         call find_mean_sig(pitch,x,y)
         write(os,"('    average pitch',1p,e12.5,' and sigma',e12.5,' cm for layer #',0p,i3)")x,y,layer
         call output(trim(os))
         call find_mean_sig(log(pitch),x,y)
         write(os,"('    ln-normal pitch exp(mu)',1p,e12.5,' (cm), sigma',e12.5,' for layer #',0p,i3)")exp(x),y,layer
         call output(trim(os))
         call indmed(pitch,i)
         write(os,"('    median pitch',1p,e12.5,' cm for layer #',0p,i3)")pitch(i),layer
         call output(trim(os))
       endif


       do n=1,num_grain
         grain_area(n)=grain2gridindex(n)-grain2gridindex(n-1)
       enddo; grain_area=grain_area*dx*dy

!do n=1,size(pitch);write(81,"(i6,x,e12.5)")n,pitch(n)*1.d7;enddo
!do n=1,size(grain_area);write(82,"(i6,2(x,e12.5))")n,sqrt(grain_area(n)*1e14/pi)*2.d0,grain_area(n)*1e14;enddo
!do n=1,num_grain; write(88,*) n,2.d0*sqrt(grain_area(n)/pi);enddo
       if (dist_io.ne.0) then
         open(dist_io,file='d.dat',form='formatted',status='unknown')
         call plot_histogram('d.dat',size(grain_area),2.d0*sqrt(grain_area/pi),ionum=dist_io)
         close(dist_io)
         open(dist_io,file='pitch.dat',form='formatted',status='unknown')
         call plot_histogram('pitch.dat',size(pitch),pitch,ionum=dist_io)
         close(dist_io)
       endif

       if (.not.silent_p) then
         i1=-nx*ny-1; j1=-i1
         do n=1,num_grain
           i1=max(i1,grain2gridindex(n)-grain2gridindex(n-1)); j1=min(j1,grain2gridindex(n)-grain2gridindex(n-1))
         enddo
         call indmed(grain_area,j)
         write(os,"('    minimum, maximum and median grain area:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")dble(j1),dble(i1),grain_area(j)/(dx*dy),layer
         call output(trim(os))
         write(os,"('    minimum, maximum and median grain area:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")j1*dx*dy,i1*dx*dy,grain_area(j),layer
         call output(trim(os))
         call find_mean_sig(grain_area,x,y)
         write(os,"('    average grain area:',1p,e12.5,' and sigma:',e12.5,' cm^2 for layer #',0p,i3)")x,y,layer
         call output(trim(os))
         call find_mean_sig(log(grain_area),x,y)
         write(os,"('    ln-normal grain area exp(mu): '1p,e12.5,' (cm^2), sigma of grain area:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
         call output(trim(os))

         x=-1.d30; y=1.d30
         do n=1,num_grain
           x=max(x,2.d0*sqrt(grain_area(n)/pi))
           y=min(y,2.d0*sqrt(grain_area(n)/pi))
         enddo
         call indmed(2.d0*sqrt(grain_area/pi),i)
         write(os,"('    minimum, maximum and median grain diameter:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")y/sqrt(dx*dy),x/sqrt(dx*dy),2.d0*sqrt(grain_area(i)/(pi*dx*dy)),layer
         call output(trim(os))
         write(os,"('    minimum, maximum and median grain diameter:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")y,x,2.d0*sqrt(grain_area(i)/pi),layer
         call output(trim(os))
         call find_mean_sig(2.d0*sqrt(grain_area/pi),x,y)
         write(os,"('    average grain diameter:',1p,e12.5,' and sigma:',e12.5,' cm for layer #',0p,i3)")x,y,layer
         call output(trim(os))
         call find_mean_sig(log(2.d0*sqrt(grain_area/pi)),x,y)
         write(os,"('    ln-normal grain diameter exp(mu): '1p,e12.5,' (cm), sigma of grain diameter:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
         call output(trim(os))
       endif

       if (num_small_grain.ne.0) then

         if (dist_io.ne.0) then
           open(dist_io,file='pitch_big.dat',form='formatted',status='unknown')
           call plot_histogram('pitch_big.dat',size(pitchb),pitchb,ionum=dist_io)
           close(dist_io)
         endif

         if (.not.silent_p) then
           call find_mean_sig(pitchb,x,y)
           write(os,"('     average big grain pitch',1p,e12.5,' and sigma',e12.5,' cm for layer #',0p,i3)")x,y,layer; call output(trim(os))
           call find_mean_sig(log(pitchb),x,y);
           write(os,"('     ln-normal big grain pitch exp(mu)',1p,e12.5,' (cm), sigma',e12.5,' for layer #',0p,i3)")exp(x),y,layer; call output(trim(os))
           call indmed(pitchb,i)
           write(os,"('     median big grain pitch',1p,e12.5,' cm for layer #',0p,i3)")pitch(i),layer; call output(trim(os))

  !big grains
           i1=-nx*ny-1; j1=-i1
           do n=1,num_grain
             if (sort_grain(n).gt.num_small_grain) then; i1=max(i1,grain2gridindex(n)-grain2gridindex(n-1)); j1=min(j1,grain2gridindex(n)-grain2gridindex(n-1)); endif
           enddo
           call indmed(pack(grain_area,mask=sort_grain.gt.num_small_grain),i); do n=1,num_grain; if (sort_grain(n).gt.num_small_grain) i=i-1; if (i.le.0) exit;enddo
           write(os,"('     minimum, maximum and median grain big area:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")dble(j1),dble(i1),grain_area(n)/(dx*dy),layer; call output(trim(os))
           write(os,"('     minimum, maximum and median grain big area:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")j1*dx*dy,i1*dx*dy,grain_area(n),layer; call output(trim(os))
           call find_mean_sig(pack(grain_area,mask=sort_grain.gt.num_small_grain),x,y)
           write(os,"('     average big grain area:',1p,e12.5,' and sigma:',e12.5,' cm^2 for layer #',0p,i3)")x,y,layer; call output(trim(os))
           call find_mean_sig(log(pack(grain_area,mask=sort_grain.gt.num_small_grain)),x,y)
           write(os,"('     ln-normal big grain area exp(mu): '1p,e12.5,' (cm^2), sigma of big grain area:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer; call output(trim(os))

           x=-1.d30; y=1.d30
           do n=1,num_grain
             if (sort_grain(n).gt.num_small_grain) then; x=max(x,2.d0*sqrt(grain_area(n)/pi)); y=min(y,2.d0*sqrt(grain_area(n)/pi)); endif
           enddo
           call indmed(2.d0*sqrt(pack(grain_area,mask=sort_grain.gt.num_small_grain)/pi),i); do n=1,num_grain; if (sort_grain(n).gt.num_small_grain) i=i-1; if (i.le.0) exit;enddo
           write(os,"('     minimum, maximum and median big grain diameter:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")y/sqrt(dx*dy),x/sqrt(dx*dy),2.d0*sqrt(grain_area(n)/(pi*dx*dy)),layer
           call output(trim(os))
           write(os,"('     minimum, maximum and median big grain diameter:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")y,x,2.d0*sqrt(grain_area(n)/pi),layer; call output(trim(os))
           call find_mean_sig(2.d0*sqrt(pack(grain_area,mask=sort_grain.gt.num_small_grain)/pi),x,y)
           write(os,"('     average big grain diameter:',1p,e12.5,' and sigma:',e12.5,' cm for layer #',0p,i3)")x,y,layer; call output(trim(os))
           call find_mean_sig(log(2.d0*sqrt(pack(grain_area,mask=sort_grain.gt.num_small_grain)/pi)),x,y)
           write(os,"('     ln-normal big grain diameter exp(mu): '1p,e12.5,' cm, (sigma) of big grain diameter:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer; call output(trim(os))

  !small grains
           i1=-nx*ny-1; j1=-i1
           do n=1,num_grain
             if (sort_grain(n).le.num_small_grain) then; i1=max(i1,grain2gridindex(n)-grain2gridindex(n-1)); j1=min(j1,grain2gridindex(n)-grain2gridindex(n-1)); endif
           enddo
           call indmed(pack(grain_area,mask=sort_grain.le.num_small_grain),i); do n=1,num_grain; if (sort_grain(n).le.num_small_grain) i=i-1; if (i.le.0) exit;enddo
           write(os,"('      minimum, maximum and median grain small area:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")dble(j1),dble(i1),grain_area(n)/(dx*dy),layer; call output(trim(os))
           write(os,"('      minimum, maximum and median grain small area:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")j1*dx*dy,i1*dx*dy,grain_area(n),layer; call output(trim(os))
           call find_mean_sig(pack(grain_area,mask=sort_grain.le.num_small_grain),x,y)
           write(os,"('      average small grain area:',1p,e12.5,' and sigma:',e12.5,' cm^2 for layer #',0p,i3)")x,y,layer; call output(trim(os))
           call find_mean_sig(log(pack(grain_area,mask=sort_grain.le.num_small_grain)),x,y)
           write(os,"('      ln-normal small grain area exp(mu): '1p,e12.5,' (cm^2), sigma of small grain area:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer; call output(trim(os))

           x=-1.d30; y=1.d30
           do n=1,num_grain
             if (sort_grain(n).le.num_small_grain) then; x=max(x,2.d0*sqrt(grain_area(n)/pi)); y=min(y,2.d0*sqrt(grain_area(n)/pi)); endif
           enddo
           call indmed(2.d0*sqrt(pack(grain_area,mask=sort_grain.le.num_small_grain)/pi),i); do n=1,num_grain; if (sort_grain(n).le.num_small_grain) i=i-1; if (i.le.0) exit;enddo
           write(os,"('      minimum, maximum and median small grain diameter:',1p,3(e12.5,1x),'nodes for layer #',0p,i3)")y/sqrt(dx*dy),x/sqrt(dx*dy),2.d0*sqrt(grain_area(n)/(pi*dx*dy)),layer
           call output(trim(os))
           write(os,"('      minimum, maximum and median small grain diameter:',1p,3(e12.5,1x),'cm^2 for layer #',0p,i3)")y,x,2.d0*sqrt(grain_area(n)/pi),layer; call output(trim(os))
           call find_mean_sig(2.d0*sqrt(pack(grain_area,mask=sort_grain.le.num_small_grain)/pi),x,y)
           write(os,"('      average small grain diameter:',1p,e12.5,' and sigma:',e12.5,' cm for layer #',0p,i3)")x,y,layer; call output(trim(os))
           call find_mean_sig(log(2.d0*sqrt(pack(grain_area,mask=sort_grain.le.num_small_grain)/pi)),x,y)
           write(os,"('      ln-normal small grain diameter exp(mu): '1p,e12.5,' (cm), sigma of small grain diameter:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer; call output(trim(os))

         endif

       endif


       if (get_plot(layer)) then
         write(os,"('pitch for layer #',0p,i3)")layer
         call plot_histogram(trim(os),20,pitch*1.d7)
         write(os,"('big grain pitch for layer #',0p,i3)")layer
         if (num_small_grain.gt.0) call plot_histogram(trim(os),20,pitchb*1.d7)
         if (get_plot(layer,grain_hist=.true.)) then
           write(os,"('grain area for layer #',0p,i3)")layer
           call plot_histogram(trim(os),20,grain_area*1.d14)
           write(os,"('grain diameter for layer #',0p,i3)")layer
           call plot_histogram(trim(os),20,sqrt(grain_area/pi)*2.d7)
         endif
       endif

       deallocate(pos_seed,rad_seed,sort_grain,pitch) !,pp)
       if (associated(pitchb)) deallocate(pitchb)


     endif   !grain or not
     ok_p=.true.
   endif  !this layer is defined?
  contains

   function radius(x,y) result(r)
    real(8),intent(in)::x,y
    real(8)::r
    r=x*x+y*y
   end function
   function iwrp(i,n) result(r)
    integer,intent(in)::i,n
    integer::r
    r=mod(i-1,n)+1-n*(-1+sign(1,mod(i-1,n)))/2
   end function
   function iwrpx(i) result(r)
    integer,intent(in)::i
    integer::r
    r=iwrp(i,nx)
   end function
   function iwrpy(i) result(r)
    integer,intent(in)::i
    integer::r
    r=iwrp(i,ny)
   end function


   function cdf_serf(f,fs,mean,sig,means,sigs,random,reset) result (s)
   !appears to use normal cdf, returns s corresponding to cdf = f
    real(8),intent(in)::f,fs,mean,sig,means,sigs,random
    logical,intent(in),optional::reset
    real(8),save::x1=10.d0,x0=0.d0,dx=-1.d0
    real(8)::s,sig0,sig0s,x,d,v

    s=0.d0
    if (present(reset)) then
      if (reset) then
        x1=10.d0; x0=0.d0; dx=-1.d0; s=0.d0
      endif
      return
    endif

    sig0=1.d0/(sqrt(2.d0)*sig)
    sig0s=1.d0/(sqrt(2.d0)*sigs)
    d=0.5d0*erfc(means*sig0s)  !missing probability

!print *,f,d,fs
    if (dx.le.0.d0) dx=(x1-x0)/100.d0

    s=(x1+x0)/2.d0
 
    x=fs*(0.5d0*(1.d0+erf((s-means)*sig0s))-d)/(1.d0-d)+(1.d0-fs)*0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
    if (x.eq.0.d0) then
      return
    elseif (x.lt.0.d0) then
      do while (fs*(0.5d0*(1.d0+erf((x1-means)*sig0s))-d)/(1.d0-d)+(1.d0-fs)*0.5d0*(1.d0+erf((log(x1)-mean)*sig0))-f.lt.0.d0)
        x0=x1
        x1=x1+dx
      enddo
    else
      do while (fs*(0.5d0*(1.d0+erf((x0-means)*sig0s))-d)/(1.d0-d)+(1.d0-fs)*0.5d0*(1.d0+erf((log(x0)-mean)*sig0))-f.ge.0.d0)
        x1=x0
        x0=max(0.d0,x0-dx)
      enddo
    endif
    dx=x1-x0;
    s=(x1+x0)/2.d0; x=fs*(0.5d0*(1.d0+erf((s-means)*sig0s))-d)/(1.d0-d)+(1.d0-fs)*0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
!print *,s,x
!print *,fs*(0.5d0*(1.d0+erf((s-means)*sig0s))-d)/(1.d0-d),(1.d0-fs)*0.5d0*(1.d0+erf((log(s)-mean)*sig0))
    do while (abs(x).gt.1.d-6)
!print *,s,x
      if (x.gt.0.d0) then
        x1=s
      else
        x0=s
      endif
      s=(x1+x0)/2.d0; x=fs*(0.5d0*(1.d0+erf((s-means)*sig0s))-d)/(1.d0-d)+(1.d0-fs)*0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
    enddo
    x1=s
    x0=max(0.d0,x-dx)

    v=exp(-0.5d0*((s-means)/sigs)**2)
    if (random.le.fs*v/(fs*v+(1.d0-fs)*exp(-0.5d0*((log(s)-mean)/sig)**2)*sigs/(sig*s))) s=-s

   end function

   function cdf_nerf(f0,mean,sig,reset) result (s)
   !appears to use normal cdf, returns s corresponding to cdf = f
    real(8),intent(in)::f0,mean,sig
    logical,intent(in),optional::reset

    real(8),save::x1=10.d0,x0=0.d0,dx=-1.d0
    real(8)::s,sig0,x,f,d

    s=0.d0
    if (present(reset)) then
      if (reset) then
        x1=10.d0; x0=0.d0; dx=-1.d0; s=0.d0
      endif
      return
    endif

    sig0=1.d0/(sqrt(2.d0)*sig)
    d=0.5d0*erfc(mean*sig0)  !missing probability
    f=min(1.d0,max(d,f0*(1.d0-d)+d))

    if (dx.le.0.d0) dx=(x1-x0)/100.d0

    s=(x1+x0)/2.d0
 
    x=0.5d0*(1.d0+erf((s-mean)*sig0))-f
    if (x.eq.0.d0) then
      return
    elseif (x.lt.0.d0) then
      do while (0.5d0*(1.d0+erf((x1-mean)*sig0))-f.lt.0.d0)
        x0=x1
        x1=x1+dx
      enddo
    else
      do while (0.5d0*(1.d0+erf((x0-mean)*sig0))-f.ge.0.d0)
        x1=x0
        x0=max(0.d0,x0-dx)
      enddo
    endif
    dx=x1-x0;
    s=(x1+x0)/2.d0; x=0.5d0*(1.d0+erf((s-mean)*sig0))-f
    do while (abs(x).gt.1.d-6)
      if (x.gt.0.d0) then
        x1=s
      else
        x0=s
      endif
      s=(x1+x0)/2.d0; x=0.5d0*(1.d0+erf((s-mean)*sig0))-f
    enddo
    x1=s
    x0=max(0.d0,x-dx)
 
   end function
   function cdf_erf(f,mean,sig,reset) result (s)
   !appears to use ln-normal cdf, returns s corresponding to cdf = f
    real(8),intent(in)::f,mean,sig
    logical,intent(in),optional::reset

    real(8),save::x1=10.d0,x0=0.d0,dx=-1.d0
    real(8)::s,sig0,x

    s=0.d0
    if (present(reset)) then
      if (reset) then
        x1=10.d0; x0=0.d0; dx=-1.d0; s=0.d0
      endif
      return
    endif

    if (dx.le.0.d0) dx=(x1-x0)/100.d0
    sig0=1.d0/(sqrt(2.d0)*sig)

    s=(x1+x0)/2.d0
 
    x=0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
    if (x.eq.0.d0) then
      return
    elseif (x.lt.0.d0) then
      do while (0.5d0*(1.d0+erf((log(x1)-mean)*sig0))-f.lt.0.d0)
        x0=x1
        x1=x1+dx
      enddo
    else
      do while (0.5d0*(1.d0+erf((log(x0)-mean)*sig0))-f.ge.0.d0)
        x1=x0
        x0=max(0.d0,x0-dx)
      enddo
    endif
    dx=x1-x0;
    s=(x1+x0)/2.d0; x=0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
    do while (abs(x).gt.1.d-6)
      if (x.gt.0.d0) then
        x1=s
      else
        x0=s
      endif
      s=(x1+x0)/2.d0; x=0.5d0*(1.d0+erf((log(s)-mean)*sig0))-f
    enddo
    x1=s
    x0=max(0.d0,x-dx)
 
   end function
  end function

  function construct_hx_vertical(nx,ny,dx,dy,layer,silent_p, &
                   unum_grain,ugrain_area,ums,ugrid2grain,ugrain2gridindex,udz,ugrain2grid,uhx,uhxgrainindex,uhxgrain, &
                   lnum_grain,lgrain_area,lms,lgrid2grain,lgrain2gridindex,ldz,lgrain2grid,lhx,lhxgrainindex,lhxgrain) result (ok_p)
   use io_m, only: output
   use random_m, only: random_normal, random_uniform

   integer,intent(in)::unum_grain,nx,ny,lnum_grain
   real(8),intent(in)::ugrain_area(unum_grain),ums(unum_grain),dx,dy,lgrain_area(lnum_grain),lms(lnum_grain),udz,ldz
   integer,intent(in)::ugrain2gridindex(0:unum_grain),ugrain2grid(:,:),ugrid2grain(nx,ny),layer
   integer,intent(in)::lgrain2gridindex(0:lnum_grain),lgrain2grid(:,:),lgrid2grain(nx,ny)
   logical,intent(in)::silent_p

   real(8),pointer::uhx(:),lhx(:)
   integer,pointer::uhxgrainindex(:),lhxgrainindex(:),uhxgrain(:),lhxgrain(:)

   logical::ok_p
   integer::n,n1,n2,num_nn,ng,k,i
   real(8)::ungrain(lnum_grain),lngrain(unum_grain)
   real(8)::y
   type(grain_s),pointer::p
   character(200)::os

   ok_p=.false.
   nullify(p); if (layer.gt.1.or.layer.le.size(grain)) p=>grain(layer)%p
   if (associated(p)) then
     if (p%ex_upper.eq.0.d0) then
       write(os,"('  there is no exchange between layer #',0p,i3,' and #',i3)")layer-1,layer
       if (.not.silent_p) call output(trim(os))
       if (associated(uhx)) deallocate(uhx); nullify(uhx)
       if (associated(uhxgrainindex)) deallocate(uhxgrainindex); nullify(uhxgrainindex)
       if (associated(uhxgrain)) deallocate(uhxgrain); nullify(uhxgrain)
       if (associated(lhx)) deallocate(lhx); nullify(lhx)
       if (associated(lhxgrainindex)) deallocate(lhxgrainindex); nullify(lhxgrainindex)
       if (associated(lhxgrain)) deallocate(lhxgrain); nullify(lhxgrain)
     else
       write(os,"('  constructing exchange between layer #',0p,i3,' and #',i3)")layer-1,layer
       if (.not.silent_p) call output(trim(os))
       n=min(unum_grain,lnum_grain)
       if (unum_grain.eq.lnum_grain.and.all(ugrid2grain.eq.lgrid2grain)) then !two layers have identical grain structure
         write(os,"('   grain structure is identical between layers')"); if (.not.silent_p) call output(trim(os))
         if (associated(uhx)) then; if (size(uhx).ne.n) then; deallocate(uhx,uhxgrain); nullify(uhx,uhxgrain); endif; endif
         if (associated(lhx)) then; if (size(lhx).ne.n) then; deallocate(lhx,lhxgrain); nullify(lhx,lhxgrain); endif; endif
         if (associated(uhxgrainindex)) then; if (size(uhxgrainindex).ne.n+1) then; deallocate(uhxgrainindex); nullify(uhxgrainindex); endif; endif
         if (associated(lhxgrainindex)) then; if (size(lhxgrainindex).ne.n+1) then; deallocate(lhxgrainindex); nullify(lhxgrainindex); endif; endif

         if (.not.associated(uhxgrainindex)) allocate(uhxgrainindex(0:n)); uhxgrainindex=0
         if (.not.associated(lhxgrainindex)) allocate(lhxgrainindex(0:n)); lhxgrainindex=0
         if (.not.associated(uhx)) allocate(uhx(n),uhxgrain(n))
         if (.not.associated(lhx)) allocate(lhx(n),lhxgrain(n))

!if (lbound(uhxgrainindex,dim=1).ne.0) then; print *,'construct_hx_vertical u';read *;endif
!if (lbound(lhxgrainindex,dim=1).ne.0) then; print *,'construct_hx_vertical l';read *;endif

         do n=1,n
           uhxgrainindex(n)=uhxgrainindex(n-1)+1
           uhxgrain(n)=n; uhx(n)=ugrain_area(n)
           lhxgrainindex(n)=lhxgrainindex(n-1)+1
           lhxgrain(n)=n; lhx(n)=lgrain_area(n)
         enddo

       else !two layers do not have identical grain structures
         if (associated(uhxgrainindex)) then; if (size(uhxgrainindex).ne.unum_grain+1) then; deallocate(uhxgrainindex); nullify(uhxgrainindex); endif; endif
         if (.not.associated(uhxgrainindex)) allocate(uhxgrainindex(0:unum_grain))
         if (associated(lhxgrainindex)) then; if (size(lhxgrainindex).ne.lnum_grain+1) then; deallocate(lhxgrainindex); nullify(lhxgrainindex); endif; endif
         if (.not.associated(lhxgrainindex)) allocate(lhxgrainindex(0:lnum_grain))
!if (lbound(uhxgrainindex,dim=1).ne.0) then; print *,'construct_hx_vertical u';read *;endif
!if (lbound(lhxgrainindex,dim=1).ne.0) then; print *,'construct_hx_vertical l';read *;endif
         num_nn=0
         do n=1,unum_grain; ungrain=0.d0
           do k=ugrain2gridindex(n-1)+1,ugrain2gridindex(n)
             ng=lgrid2grain(ugrain2grid(1,k),ugrain2grid(2,k)); if (ng.ne.0) ungrain(ng)=1.d0
           enddo
           num_nn=num_nn+nint(sum(ungrain))
         enddo
         if (associated(uhxgrain)) then; if (size(uhxgrain).ne.num_nn) then; deallocate(uhx,uhxgrain); nullify(uhx,uhxgrain); endif; endif
         if (.not.associated(uhx)) allocate(uhx(num_nn),uhxgrain(num_nn))

         num_nn=0; lhxgrainindex=0
         do n=1,lnum_grain; lngrain=0.d0
           do k=lgrain2gridindex(n-1)+1,lgrain2gridindex(n)
             ng=ugrid2grain(lgrain2grid(1,k),lgrain2grid(2,k)); if (ng.ne.0) lngrain(ng)=1.d0
           enddo
           num_nn=num_nn+nint(sum(lngrain))
         enddo
         if (associated(lhxgrain)) then; if (size(lhxgrain).ne.num_nn) then; deallocate(lhx,lhxgrain); nullify(lhx,lhxgrain); endif; endif
         if (.not.associated(lhx)) allocate(lhx(num_nn),lhxgrain(num_nn))

         uhxgrainindex=0; n2=0
         do n=1,unum_grain; ungrain=0.d0
           do k=ugrain2gridindex(n-1)+1,ugrain2gridindex(n)
             ng=lgrid2grain(ugrain2grid(1,k),ugrain2grid(2,k)); if (ng.ne.0) ungrain(ng)=ungrain(ng)+1.d0
           enddo
           do n1=1,lnum_grain; uhxgrainindex(n)=n2
             if (ungrain(n1).ne.0) then
               n2=n2+1; uhxgrainindex(n)=n2; uhxgrain(n2)=n1; uhx(n2)=ungrain(n1)*dx*dy
             endif
           enddo
         enddo

         lhxgrainindex=0; n2=0
         do n=1,lnum_grain; lngrain=0.d0
           do k=lgrain2gridindex(n-1)+1,lgrain2gridindex(n)
             ng=ugrid2grain(lgrain2grid(1,k),lgrain2grid(2,k)); if (ng.ne.0) lngrain(ng)=lngrain(ng)+1.d0
           enddo
           do n1=1,unum_grain; lhxgrainindex(n)=n2
             if (lngrain(n1).ne.0) then
               n2=n2+1; lhxgrainindex(n)=n2; lhxgrain(n2)=n1; lhx(n2)=lngrain(n1)*dx*dy
             endif
           enddo
         enddo
       endif

          call diag1(lnum_grain,lhxgrainindex,lhxgrain,lhx,lgrain_area,layer,layer-1)
          call diag1(unum_grain,uhxgrainindex,uhxgrain,uhx,ugrain_area,layer-1,layer)

      !ok, u(l)hxgrainindex(0:u(l)num_grain) gives the index to u(l)hx and u(l)hxgrain
      !  u(l)hxgrain is the grain index on the lower(upper) layer
      !  u(l)hx is currently the total area shared between grains
       do n=1,lnum_grain     !loop over lower layer
         do i=lhxgrainindex(n-1)+1, lhxgrainindex(n)
           n1=lhxgrain(i)    !grain n1 on upper layer coupling to grain n on lower layer

           do n2=uhxgrainindex(n1-1)+1, uhxgrainindex(n1)  !find grain n on lower layer that couples with upper grain n1
             if (uhxgrain(n2).eq.n) exit
           enddo
           if (n2.gt.uhxgrainindex(n1)) then; call output(' **** internal error in construct_hx_vertical, aborting simulation! ****',.true.); stop; endif
           if (lhx(i).ne.uhx(n2)) then
!                  print *,lhx(i),uhx(n2),1.d0-min(lhx(i),uhx(n2))/max(lhx(i),uhx(n2))
              call output(' **** internal error in construct_hx_vertical (face coupling not consistent), aborting simulation! ****',.true.); stop
           endif

           y=1.d0 
           if (p%ex_upper_sigma.ne.0.d0) then
             y = 1.d0 + p%ex_upper_sigma*random_normal(p%random_exv)   !gaussian distribution
             if (p%ex_upper_sigma.gt.0.d0) y=exp(y-1.d0)               !ln-normal distribution
           endif

         ! l(u)hx has the area in commone between grains
           lhx(i )=p%ex_upper*y*lhx(i )/(lms(n )*ldz*lgrain_area(n ))
           uhx(n2)=p%ex_upper*y*uhx(n2)/(ums(n1)*udz*ugrain_area(n1))
         enddo
       enddo

          call diag(lnum_grain,lhxgrainindex,lhxgrain,lhx,layer,layer-1,lgrain2gridindex)
          call diag(unum_grain,uhxgrainindex,uhxgrain,uhx,layer-1,layer,ugrain2gridindex)
         
     endif
     ok_p=.true.
   endif
  contains 
    subroutine diag1(num_grain,hxgrainindex,hxgrain,hx,grain_area,i1,i2)
     integer,intent(in)::num_grain,hxgrainindex(0:),hxgrain(:),i1,i2
     real(8),intent(in)::hx(:),grain_area(:)
     real(8)::x,y,z,sig

       if (.not.silent_p) then
     x=0.d0; y=1.d22; z=-1.d22
     do n=1,num_grain; sig=0.d0
       sig=hxgrainindex(n)-hxgrainindex(n-1)
       x=x+sig; y=min(y,sig); z=max(z,sig)
     enddo; x=x/num_grain
     write(os, "('   minimum, maximum and average number of interlayer neighboring grains:',1p,3(e12.5,1x),'for layer #',0p,i3,' due to #',i3)")y,z,x,i1,i2
     call output(trim(os))

     x=0.d0; y=1.d22; z=-1.d22
     do n=1,num_grain; sig=0.d0
       do i=hxgrainindex(n-1)+1,hxgrainindex(n)
         sig=sig+hx(i)
       enddo; sig=sig/(dx*dy)
       x=x+sig; y=min(y,sig); z=max(z,sig)
     enddo; x=x/num_grain
     write(os, "('   minimum, maximum and average number of interlayer faces:',1p,3(e12.5,1x),'for layer #',0p,i3,' due to #',i3)")y,z,x,i1,i2
     call output(trim(os))

     x=0.d0; y=1.d22; z=-1.d22
     do n=1,num_grain; sig=0.d0
       do i=hxgrainindex(n-1)+1,hxgrainindex(n)
         sig=sig+hx(i)
       enddo; sig=sig/grain_area(n)
       x=x+sig; y=min(y,sig); z=max(z,sig)
     enddo; x=x/num_grain
     write(os, "('   minimum, maximum and average interlayer ratio of faces to nodes:',1p,3(e12.5,1x),'for layer #',0p,i3,' due to #',i3)")y,z,x,i1,i2
     call output(trim(os))
       endif

    end subroutine
    subroutine diag(num_grain,hxgrainindex,hxgrain,hx,i1,i2,grain2gridindex)
     use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour
     integer,intent(in)::num_grain,hxgrainindex(0:),hxgrain(:),i1,i2,grain2gridindex(0:)
     real(8),intent(in)::hx(:)
     real(8)::x,y,z,sig
     real(8)::dt(num_grain,2)

       if (.not.silent_p) then
     x=0.d0; y=1.d22; z=-1.d22
     do n=1,num_grain; sig=0.d0
       do i=hxgrainindex(n-1)+1,hxgrainindex(n)
         sig=sig+hx(i)
       enddo
       x=x+sig; y=min(y,sig); z=max(z,sig)
     enddo; x=x/num_grain
     write(os, "('   minimum, maximum and average interlayer exchange for a grain:',1p,3(e12.5,1x),'(Oe) for layer #',0p,i3,' from #',i3)")y,z,x,i1,i2
     call output(trim(os))

     z=0.d0
     do n=1,num_grain; sig=0.d0
       do i=hxgrainindex(n-1)+1,hxgrainindex(n)
         sig=sig+hx(i)
       enddo
       z=z+(sig-x)*(sig-x)
     enddo; z=sqrt(z/max(1,num_grain-1))
     write(os, "('   standard deviation for interlayer exchange for a grain:',1p,1(e12.5,1x),'(Oe) for layer #',0p,i3,' from #',i3)")z,i1,i2
     call output(trim(os))
       endif

       if (get_plot(i1,hv=.true.)) then
         dt=0.d0
         do n=1,num_grain
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             dt(n,1)=dt(n,1)+hx(i)
           enddo
           dt(n,2)=grain2gridindex(n)-grain2gridindex(n-1)
         enddo
         write(os,"('Hx histogram for layer #',0p,i3,' from #',i3)")i1,i2
         if (get_plot(i1,hv_hist=.true.)) call plot_histogram(trim(os),20,dt(:,1)/1.d4)
         write(os,"('Hx vs grain area for layer #',0p,i3,' from #',i3)")i1,i2
         if (get_plot(i1,hv_vol=.true.)) call plot_scatter(trim(os),dt(:,2)*dx*dy*1.d14,dt(:,1)/1.d4,.true.)
         write(os,"('max Hx for layer #',0p,i3,' from #',i3)")i1,i2
         call plot_contour(trim(os),dt(:,1)/1.d4,.true.,i1)
       endif
    end subroutine
  end function

  function construct_hx_transverse(nx,ny,dx,dy,num_grain,grain_area,ms,grid2grain,grain2gridindex,grain2grid,ilayer,isilent_p, &
                                  hx,hxgrainindex,hxgrain) result (ok_p)
   use io_m, only: output
   use random_m, only: random_normal, random_uniform
   use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour

   integer,intent(in)::num_grain,nx,ny
   real(8),intent(in)::grain_area(num_grain),ms(num_grain),dx,dy
   integer,intent(in),optional::ilayer,grain2gridindex(0:num_grain),grain2grid(:,:),grid2grain(nx,ny)
   logical,intent(in),optional::isilent_p

   real(8),pointer::hx(:)
   integer,pointer::hxgrainindex(:),hxgrain(:)
   logical::silent_p,ok_p
   integer::layer,ngrain(num_grain)
   type(grain_s),pointer::p
   integer::i,j,n,num_nn,ng(4),k
   real(8)::mean,sigma,x,y,dt(num_grain)
   character(200)::os

   ok_p=.false.
   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   if (associated(p)) then
     if (associated(hx)) then; if (size(hx).ne.num_grain)then; deallocate(hx); nullify(hx);endif; endif
     if (p%ex.eq.0.d0) then
       write(os,"('  there is no intralayer exchange for layer #',0p,i3)")layer
       if (.not.silent_p) call output(trim(os))
       if (associated(hx)) deallocate(hx); nullify(hx)
       if (associated(hxgrainindex)) deallocate(hxgrainindex); nullify(hxgrainindex)
       if (associated(hxgrain)) deallocate(hxgrain); nullify(hxgrain)
     else
       write(os,"('  constructing intralayer exchange for layer #',0p,i3)")layer
       if (.not.silent_p) call output(trim(os))

       num_nn=0
       do n=1,num_grain
         ngrain=0
         do k=grain2gridindex(n-1)+1,grain2gridindex(n)
           i=grain2grid(1,k); j= grain2grid(2,k)
           ng= (/grid2grain(i,iwrpy(j-1)), grid2grain(iwrpx(i-1),j), grid2grain(iwrpx(i+1),j),grid2grain(i,iwrpy(j+1)) /)
           do i=1,4
             if (ng(i).ne.0.and.ng(i).ne.n) ngrain(ng(i))=1
           enddo
         enddo
         num_nn=num_nn+sum(ngrain)
       enddo
       if (associated(hx)) then; if (size(hx).ne.num_nn) then; deallocate(hx); nullify(hx); endif; endif
       if (associated(hxgrainindex)) then; if (size(hxgrainindex).ne.num_grain+1) then; deallocate(hxgrainindex); nullify(hxgrainindex); endif; endif
       if (associated(hxgrain)) then; if (size(hxgrain).ne.num_nn) then; deallocate(hxgrain); nullify(hxgrain); endif; endif
       if (.not.associated(hx)) allocate(hx(num_nn))
       if (.not.associated(hxgrainindex)) allocate(hxgrainindex(0:num_grain))
       if (.not.associated(hxgrain)) allocate(hxgrain(num_nn))
!if (lbound(hxgrainindex,dim=1).ne.0) then; print *,'construct_hx_transverse';read *;endif
       hxgrainindex=0; ngrain=0
       do n=1,num_grain
         dt=0.d0
         do k=grain2gridindex(n-1)+1,grain2gridindex(n)
           i=grain2grid(1,k); j= grain2grid(2,k)
           ng= (/grid2grain(i,iwrpy(j-1)), grid2grain(iwrpx(i-1),j), grid2grain(iwrpx(i+1),j),grid2grain(i,iwrpy(j+1)) /)
           if (ng(1).ne.0.and.ng(1).ne.n) dt(ng(1))=dt(ng(1))+dx
           if (ng(2).ne.0.and.ng(2).ne.n) dt(ng(2))=dt(ng(2))+dy
           if (ng(3).ne.0.and.ng(3).ne.n) dt(ng(3))=dt(ng(3))+dy
           if (ng(4).ne.0.and.ng(4).ne.n) dt(ng(4))=dt(ng(4))+dx
           do i=1,4
             if (ng(i).ne.0.and.ng(i).ne.n) ngrain(n)=ngrain(n)+1
           enddo
         enddo
         hxgrainindex(n)=hxgrainindex(n-1)
         do k=1,num_grain
           if (dt(k).gt.0) then
             hxgrainindex(n)=hxgrainindex(n)+1
             hxgrain(hxgrainindex(n))=k
             hx(hxgrainindex(n))=dt(k)
           endif
         enddo
       enddo

       !ok, hxgrainindex(n-1)+1:hxgrainindex(n) is the range of grains that are exchanged coupled to grain n
       !   hxgrain(hxgrainindex(n-1)+1:hxgrainindex(n)) = are grain numbers that are exchanged coupled to grain n
       !   hx(hxgrainindex(n-1)+1:hxgrainindex(n)) = is currently the total length shared between two grains, soon to be Hx
       if (.not.silent_p) then
         mean=0.d0; x=1.d300; y=-1.d300
         do n=1,num_grain
           mean=mean+hxgrainindex(n)-hxgrainindex(n-1)
           x=min(x,dble(hxgrainindex(n)-hxgrainindex(n-1)))
           y=max(y,dble(hxgrainindex(n)-hxgrainindex(n-1)))
         enddo;mean=mean/num_grain
         write(os,"('   minimum, maximum and average number of neighboring grains:',1p,3(e12.5,1x),' for layer #',0p,i3)")x,y,mean,layer
         call output(trim(os))

         mean=0.d0; x=1.d300; y=-1.d300
         do n=1,num_grain
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             mean=mean+hx(i); x=min(x,hx(i)); y=max(y,hx(i))
           enddo
         enddo;mean=mean/num_grain
         write(os,"('   minimum, maximum and average length of shared boundary:',1p,3(e12.5,1x),' for layer #',0p,i3)")x,y,mean,layer
         call output(trim(os))

         mean=0.d0; x=1.d300; y=-1.d300
         do n=1,num_grain
           sigma=0.d0
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             sigma=sigma+hx(i)
           enddo
           mean=mean+sigma/grain_area(n); x=min(x,sigma/grain_area(n)); y=max(y,sigma/grain_area(n))
         enddo;mean=mean/num_grain
         write(os,"('   minimum, maximum and average length of shared boundary per grain area:',1p,3(e12.5,1x),' for layer #',0p,i3)")x,y,mean,layer
         call output(trim(os))

       endif

       !ok, let's set hx
       do n=1,num_grain
         do i=hxgrainindex(n-1)+1,hxgrainindex(n)
           if (n.lt.hxgrain(i)) then   !have we done this pair yet?
             y=1.d0;x=hx(i)      !total length shared between grain n and grain hxgrain(i)
             if (p%ex_sigma.ne.0.d0) then
               y = 1.d0 + p%ex_sigma*random_normal(p%random_ex)   !gaussian distribution
               if (p%ex_sigma.gt.0.d0) y=exp(y-1.d0)        !ln-normal distribution
             endif
             y=y*p%ex
             hx(i)=x*y/(ms(n)*grain_area(n))
             do j=hxgrainindex(hxgrain(i)-1)+1,hxgrainindex(hxgrain(i))  !must be symmetric between grains
               if (hxgrain(j).eq.n) hx(j)=x*y/(ms(hxgrain(i))*grain_area(hxgrain(i)))
             enddo
          endif
         enddo
       enddo

       if (.not.silent_p) then

         mean=0.d0; x=1.d300; y=-1.d300
         do n=1,num_grain; sigma=0.d0
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             sigma=sigma+hx(i)
           enddo
!write(83,*)n,sigma
           mean=mean+sigma; x=min(x,sigma); y=max(y,sigma)
         enddo;mean=mean/num_grain
         write(os,"('   minimum, maximum and average exchange per grain:',1p,3(e12.5,1x),'Oe for layer #',0p,i3)")x,y,mean,layer
         call output(trim(os))

         sigma=0.d0
         do n=1,num_grain
           x=0.d0
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             x=x+hx(i)
           enddo
           sigma=sigma+(x-mean)**2
         enddo;sigma=sqrt(sigma/(num_grain-1))
         write(os,"('   standard deviation of exchange per grain:',1p,1(e12.5,1x),'Oe for layer #',0p,i3)")sigma,layer
         call output(trim(os))

         mean=0.d0; x=1.d300; y=-1.d300; k=0
         do n=1,num_grain
           sigma=0.d0
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             sigma=sigma+hx(i); k=k+1
             x=min(x,sigma); y=max(y,sigma)
           enddo
           mean=mean+sigma
         enddo;mean=mean/max(1,k)
         write(os,"('   minimum, maximum and average exchange for pair of grains:',1p,3(e12.5,1x),'Oe for layer #',0p,i3)")x,y,mean,layer
         call output(trim(os))

         sigma=0.d0;k=0
         do n=1,num_grain
           x=0.d0
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             x=x+hx(i)
           enddo
           sigma=sigma+(x-mean)**2;k=k+1
         enddo;sigma=sqrt(sigma/max(1,k-1))
         write(os,"('   standard deviation of exchange for pair of grains:',1p,1(e12.5,1x),'Oe for layer #',0p,i3)")sigma,layer
         call output(trim(os))

       endif

       if (get_plot(layer,hx=.true.)) then
         dt=0.d0;
         do n=1,num_grain
           do i=hxgrainindex(n-1)+1,hxgrainindex(n)
             dt(n)=dt(n)+hx(i)
           enddo
         enddo
         write(os,"('Hx histogram for layer #',0p,i3)")layer
         if (get_plot(layer,hx_hist=.true.)) call plot_histogram(trim(os),20,dt/1.d4)
         write(os,"('Hx vs grain area for layer #',0p,i3)")layer
         if (get_plot(layer,hx_vol=.true.)) call plot_scatter(trim(os),grain_area*1.d14,dt/1.d4,.true.)
         write(os,"('max Hx for layer #',i3)")layer
         call plot_contour(trim(os),dt/1.d4,.true.,layer)
       endif


     endif
     ok_p=.true.
   endif
  contains
   function iwrp(i,n) result(r)
    integer,intent(in)::i,n
    integer::r
    r=mod(i-1,n)+1-n*(-1+sign(1,mod(i-1,n)))/2
   end function
   function iwrpx(i) result(r)
    integer,intent(in)::i
    integer::r
    r=iwrp(i,nx)
   end function
   function iwrpy(i) result(r)
    integer,intent(in)::i
    integer::r
    r=iwrp(i,ny)
   end function
  end function

  function construct_hk_magnitude(num_grain,grain_area,grain_ms,grain_height,ilayer,isilent_p,hk) result (ok_p)

   use io_m, only: output
   use random_m, only: random_normal, random_uniform
   use data_m, only: get_data
   use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour

   integer,intent(in)::num_grain
   real(8),intent(in)::grain_area(num_grain),grain_ms(num_grain),grain_height
   integer,intent(in),optional::ilayer
   logical,intent(in),optional::isilent_p

   real(8),pointer::hk(:,:)
   logical::silent_p,ok_p
   integer::layer
   type(grain_s),pointer::p
   integer::i,iter
   real(8)::mean,sigma,t,xx,x,pi,dt(num_grain),k_sigma,k1,k2
   character(200)::os
   real(8),pointer::gv(:),ms(:)

   ok_p=.false.
   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   if (associated(p)) then
     write(os,"('  constructing Hk magnitude for layer #',0p,i3)")layer
     if (.not.silent_p) call output(trim(os))
     if (associated(hk)) then; if (size(hk,2).ne.num_grain)then; deallocate(hk); nullify(hk);endif; endif
     if (.not.associated(hk)) allocate(hk(2,num_grain))
     t=-1.d0; iter=1; hk(1,:)=4.d0*p%kmax/grain_ms; pi=acos(-1.d0)
     k_sigma=p%k_sigma;       k1=p%k1;       k2=p%k2; xx=1.d0
     do while (t.lt.0.d0)
       do i=1,num_grain
         do while (hk(1,i).gt.2.d0*p%kmax/grain_ms(i).or.hk(1,i).lt.p%hkmin)
           x = 1.d0 + k_sigma*random_normal(p%random_hk_mag)   !gaussian distribution
           if (k_sigma.gt.0.d0) x=exp(x-1.d0)        !ln-normal distribution
!          x=max(0.d0,x)
           x=abs(x)
           if (p%k_scale_al) then
             xx=max(0.d0, min(1.d0, 1.d0-(0.24864d-7/grain_height/abs(p%vol_scale)+0.163736d-7/sqrt(grain_area(i)/pi))))
           elseif (p%k_scale) then
!            x=x*max(0.d0,1.d0-2.d0*p%k_ac/grain_height)* &
!            xx=max(0.d0,1.d0-p%num_end*p%k_ac/grain_height)* &
             xx=max(0.d0,1.d0-2.d0*p%k_ac/grain_height/abs(p%vol_scale))* &
                max(0.d0,1.d0-p%k_ac/sqrt(grain_area(i)/pi))**2
           endif
           hk(1,i)=xx*x*k1; hk(2,i)=xx*x*k2
         enddo
         if (p%vol_scale.gt.1.d0.and.p%k_scale) then
           xx=max(0.d0,1.d0-p%num_end*p%k_ac/grain_height)* &
              max(0.d0,1.d0-p%k_ac/sqrt(grain_area(i)/pi))**2
           hk(1,i)=xx*x*k1; hk(2,i)=xx*x*k2
         endif
       enddo
       t=0.d0
       do i=1,num_grain
!        t=t+(hk(1,i)+hk(2,i))*grain_area(i)
         t=t+(hk(1,i)+hk(2,i))*grain_area(i)
       enddo; t=t/sum(grain_area)
       if (p%k_rescale.and.abs(1.d0-t/p%k1).gt.3.d-2.and.iter.lt.4) then
         iter=iter+1; x=p%k1/t; hk(1,:)=hk(1,:)*x; hk(2,:)=hk(2,:)*x
         if ((count(hk(1,:).lt.p%hkmin)).gt.0.or.count(hk(1,:).gt.2.d0*p%kmax/grain_ms(:)).gt.0.and.iter.lt.4)t=-1.d0
       endif
     enddo
     if (p%frac_hk_small.gt.0.d0) then
       k_sigma=p%k_sigma_small; k1=p%k1_small; k2=p%k2_small; xx=1.d0
       do i=1,num_grain; x=4.d0*p%kmax/grain_ms(i)/k1
         do while (xx*x.gt.2.d0*p%kmax/grain_ms(i)/k1.or.xx*x.lt.p%hkmin/k1)
         x = 1.d0 + k_sigma*random_normal(p%random_hk_mag_small)   !gaussian distribution
         if (k_sigma.gt.0.d0) x=exp(x-1.d0)        !ln-normal distribution
!        x=abs(x)
         if (p%k_scale_al) then
           xx=max(0.d0, min(1.d0, 1.d0-(0.24864d-7/grain_height/abs(p%vol_scale)+0.163736d-7/sqrt(grain_area(i)/pi))))
         elseif (p%k_scale) then
!          x=x*max(0.d0,1.d0-2.d0*p%k_ac/grain_height)* &
!          x=x*max(0.d0,1.d0-2.d0*p%k_ac/grain_height/p%vol_scale_tc)* &
!          xx=max(0.d0,1.d0-p%num_end*p%k_ac/grain_height)* &
!             max(0.d0,1.d0-p%k_ac/sqrt(grain_area(i)/pi))**2
           xx=max(0.d0,1.d0-2.d0*p%k_ac/grain_height/abs(p%vol_scale))* &
              max(0.d0,1.d0-p%k_ac/sqrt(grain_area(i)/pi))**2
         endif
         enddo
         if (p%vol_scale.gt.1.d0.and.p%k_scale) &
           xx=max(0.d0,1.d0-p%num_end*p%k_ac/grain_height)* &
              max(0.d0,1.d0-p%k_ac/sqrt(grain_area(i)/pi))**2
         if (random_uniform(p%random_hk_mag_small).lt.p%frac_hk_small) then; hk(1,i)=xx*x*k1; hk(2,i)=xx*x*k2; endif
         xx=1.d0
       enddo
     endif
!do i=1,num_grain;write(84,"(i6,3(x,e12.5))")i,hk(:,i),grain_area(i)*1e14;enddo
     if (.not.silent_p) then
       mean=-1.d30; sigma=1.d30
       do i=1,num_grain; mean=max(mean,hk(1,i)+hk(2,i)); sigma=min(sigma,hk(1,i)+hk(2,i)); enddo
       call indmed(hk(1,:)+hk(2,:),i)
       write(os,"('    minimum, maximum and median Hk:',1p,3(e12.5,1x),'(Oe) for layer #',0p,i3)")sigma,mean,hk(1,i)+hk(2,i),layer; call output(trim(os))
       call find_mean_sig(hk(1,:)+hk(2,:),mean,sigma)
       write(os,"('    average Hk:',1p,e12.5,' and sigma:',e12.5,' (Oe) for layer #',0p,i3)")mean,sigma,layer; call output(trim(os))
       call find_mean_sig(log(hk(1,:)+hk(2,:)),mean,sigma)
       write(os,"('    ln-normal Hk exp(mu): '1p,e12.5,' (Oe), sigma of Hk:',1p,e12.5,' for layer #',0p,i3)")exp(mean),sigma,layer
       call find_mean_sig(    hk(1,:)+hk(2,:) ,mean,sigma,grain_area)
       write(os,"('    volume average Hk:',1p,e12.5,' and sigma:',e12.5,' (Oe) for layer #',0p,i3)")mean,sigma,layer; call output(trim(os))
       call find_mean_sig(log(hk(1,:)+hk(2,:)),mean,sigma,grain_area)
       write(os,"('    ln-normal volume weighted Hk exp(mu): '1p,e12.5,' (Oe), sigma of Hk:',1p,e12.5,' for layer #',0p,i3)")exp(mean),sigma,layer
       call output(trim(os))
     endif

     if (get_plot(layer,hk=.true.)) then
       write(os,"('|Hk| histogram for layer #',0p,i3)")layer
       do i=1,num_grain; dt(i)=hk(1,i)+hk(2,i); enddo
       if (get_plot(layer,hk_hist=.true.)) call plot_histogram(trim(os),20,dt/1.d4)
       write(os,"('|Hk| vs grain area for layer #',0p,i3)")layer
       if (get_plot(layer,hk_vol=.true.)) call plot_scatter(trim(os),grain_area*1.d14,dt/1.d4,.true.)
       write(os,"('|Hk| for layer #',i3)")layer
       call plot_contour(trim(os),dt/1.d4,.true.,layer)
       write(os,"('KuV/kT for layer #',0p,i3)")layer
       call get_data(layer,gv,grain_vol_p=.true.); call get_data(layer,ms,ms_p=.true.)
       call plot_contour(trim(os),0.5d0*dt*ms*gv/(1.380662d-16*300.d0),.true.,layer)
       write(os,"('KuV/kT histogram for layer #',0p,i3)")layer
       if (get_plot(layer,hk_hist=.true.)) call plot_histogram(trim(os),20,0.5d0*dt*ms*gv/(1.380662d-16*300.d0))
     endif

     ok_p=.true.
   endif
  end function

  function construct_hk_direction(num_grain,ilayer,isilent_p,uni_axis,cubic_axis) result (ok_p)

   use io_m, only: output
   use random_m, only: random_uniform, random_normal, random_gaussian_sphere
   use transform_m, only: rotate
   use plot_m, only: get_plot, plot_histogram

   integer,intent(in)::num_grain
   integer,intent(in),optional::ilayer
   logical,intent(in),optional::isilent_p

   real(8),pointer::cubic_axis(:,:,:),uni_axis(:,:)
   logical::silent_p,ok_p
   integer::layer
   type(grain_s),pointer::p
   integer::i
   real(8)::mean,sigma,phi,theta,z,pi,tsig,tcsig,theta_n,phi_n,psi,euler(3,3,2),dt(num_grain)
   character(200)::os

   ok_p=.false.
   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   if (associated(p)) then
     pi=acos(-1.d0)
     z=p%kx*p%kx+p%ky*p%ky+p%kz*p%kz
     if (z.ne.0d0) then
        z=1.d0/z; p%kx=p%kx*z; p%ky=p%ky*z ; p%kz=p%kz*z
     else
        p%kz=1.d0
     endif
     theta_n=acos(p%kz); phi_n=0.d0; if (abs(p%ky)+abs(p%kx).ne.0.d0) phi_n=atan2(p%ky,p%kx)
     tsig=p%k_ang*pi/180.d0; tcsig=p%k_cone_sigma*pi/180.d0

     if (p%k_cubic) then
       write(os,"('  constructing cubic Hk orientation for layer #',0p,i3)")layer
       if (.not.silent_p) call output(trim(os))
       if (associated(cubic_axis)) then; if (size(cubic_axis,3).ne.num_grain)then; deallocate(cubic_axis); nullify(cubic_axis);endif; endif
       if (.not.associated(cubic_axis)) allocate(cubic_axis(3,3,num_grain))
       if (associated(uni_axis)) then; deallocate(uni_axis); endif; nullify(uni_axis)
       psi=0.d0; if (abs(p%ky)+abs(p%kx).ne.0.d0) psi=atan2(p%kx,-p%ky)
       euler(1,:,1) = (/ cos(psi), sin(psi), 0.d0 /)
       euler(2,:,1) = (/ -cos(theta_n)*sin(psi), cos(theta_n)*cos(psi), sin(theta_n) /)
       euler(3,:,1) = (/ sin(psi)*sin(theta_n), -cos(psi)*sin(theta_n), cos(theta_n) /)
       do i=1,num_grain
         phi=2.d0*pi*random_uniform(p%random_hk_cubic_angles)
         psi=2.d0*pi*random_uniform(p%random_hk_cubic_angles)
         if (p%k_3d) then
           theta = acos( max(-1.d0, min(1.d0, 1.d0-2.d0*random_uniform(p%random_hk_cubic_angles))))
         elseif (tsig.eq.0.d0) then !point along nominal
           theta=0.d0
           if (p%k_2d) theta=pi/2.d0
         elseif (p%k_uniform) then  !uniform in cos(theta), not theta
           if (p%k_2d) then
             theta = acos( cos(0.5d0*pi-abs(tsig))*(1.d0-2.d0*random_uniform(p%random_hk_cubic_angles)))
           else
             theta = acos( min(1.d0, max( cos(abs(tsig)), 1.d0-random_uniform(p%random_hk_cubic_angles)*(1.d0-cos(abs(tsig))))))
           endif
         elseif (p%k_cone) then
           theta = abs(tsig)+tcsig*random_normal(p%random_hk_cubic_angles)
         elseif (tsig.gt.0.d0) then
           theta = random_gaussian_sphere(p%random_hk_cubic_angles,tsig)
         else
           theta = tsig*random_normal(p%random_hk_cubic_angles)
         endif
         theta = abs(theta)
         theta = abs(theta-2.d0*pi*int((theta+mod(int(theta/pi),2)*pi)/(2.d0*pi)))
         euler(:,1,2) = (/  cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi), &
                            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi), &
                            sin(psi)*sin(theta) /)
         euler(:,2,2) = (/ -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi), &
                           -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi), &
                            cos(psi)*sin(theta) /)
         euler(:,3,2) = (/  sin(phi)*sin(theta), -cos(phi)*sin(theta), cos(theta) /)
         cubic_axis(:,1,i) = (/ sum(euler(:,1,1)*euler(:,1,2)), sum(euler(:,2,1)*euler(:,1,2)), sum(euler(:,3,1)*euler(:,1,2)) /)  !a-axis
         cubic_axis(:,2,i) = (/ sum(euler(:,1,1)*euler(:,2,2)), sum(euler(:,2,1)*euler(:,2,2)), sum(euler(:,3,1)*euler(:,2,2)) /)  !b-axis
         cubic_axis(:,3,i) = (/ sum(euler(:,1,1)*euler(:,3,2)), sum(euler(:,2,1)*euler(:,3,2)), sum(euler(:,3,1)*euler(:,3,2)) /)  !c-axis
       enddo
     else
       write(os,"('  constructing uniaxial Hk orientation for layer #',0p,i3)")layer
       if (.not.silent_p) call output(trim(os))
       if (associated(uni_axis)) then; if (size(uni_axis,2).ne.num_grain)then; deallocate(uni_axis); nullify(uni_axis);endif; endif
       if (.not.associated(uni_axis)) allocate(uni_axis(3,num_grain))
       if (associated(cubic_axis)) then; deallocate(cubic_axis); endif; nullify(cubic_axis)

       do i=1,num_grain
         phi=2.d0*pi*random_uniform(p%random_hk_uni_angles)
         if (p%k_3d) then
           theta = acos( max(-1.d0, min(1.d0, 1.d0-2.d0*random_uniform(p%random_hk_uni_angles))))
         elseif (tsig.eq.0.d0) then !point along nominal
           theta=0.d0
           if (p%k_2d) theta=pi/2.d0
         elseif (p%k_uniform) then  !uniform in cos(theta), not theta
           if (p%k_2d) then
             theta = acos( cos(0.5d0*pi-abs(tsig))*(1.d0-2.d0*random_uniform(p%random_hk_uni_angles)))
           else
             theta = acos( min(1.d0, max( cos(abs(tsig)), 1.d0-random_uniform(p%random_hk_uni_angles)*(1.d0-cos(abs(tsig))))))
           endif
         elseif (p%k_cone) then
           theta = abs(tsig)+tcsig*random_normal(p%random_hk_uni_angles)
         elseif (tsig.gt.0.d0) then
           theta = random_gaussian_sphere(p%random_hk_uni_angles,tsig)
         else
           theta = tsig*random_normal(p%random_hk_cubic_angles)
         endif
         if (random_uniform(p%random_hk_frac_random).le.p%k_frac_random) &
           theta = acos( max(-1.d0, min(1.d0, 1.d0-2.d0*random_uniform(p%random_hk_uni_angles))))
         if (random_uniform(p%random_hk_frac_inplane).le.p%k_frac_inplane) theta=0.5d0*pi
         theta = abs(theta)
         theta = abs(theta-2.d0*pi*int((theta+mod(int(theta/pi),2)*pi)/(2.d0*pi)))
!        call rotate(phi_n,theta_n,phi,theta,uni_axis(1,i),uni_axis(2,i),uni_axis(3,i))
         uni_axis(1:3,i)=rotate(phi_n,theta_n,phi,theta)
       enddo

!do i=1,num_grain;write(85,"(i6,3(x,e12.5))")i,uni_axis(:,i);enddo
     !compute min, max,  average and std of k.knominal

       if (.not.silent_p) then
         do i=1,num_grain; dt(i)=acos(dot_product((/p%kx,p%ky,p%kz/),uni_axis(:,i)))*180.d0/pi; enddo
         phi=-1.d30; psi=1.d30
         do i=1,num_grain; phi=max(dt(i),phi); psi=min(dt(i),psi); enddo
         call indmed(dt,i)
         write(os,"('    minimum, maximum and median k.knom:',1p,3(e12.5,1x),'(deg) for layer #',0p,i3)")psi,phi,dt(i),layer; call output(trim(os))
         call find_mean_sig(dt,mean,sigma)
         write(os,"('    average k.knom:',1p,e12.5,' and sigma:',e12.5,' (s) for layer #',0p,i3)")mean,sigma,layer; call output(trim(os))
!        call find_mean_sig(dt,mean,sigma,grain_area)
!        write(os,"('    volume average k.knom:',1p,e12.5,' and sigma:',e12.5,' (K) for layer #',0p,i3)")mean,sigma,layer; call output(trim(os))
!        call find_mean_sig(log(dt),mean,sigma,grain_area)
!        write(os,"('    ln-normal k.knom exp(mu): '1p,e12.5,', sigma of k.knom:',1p,e12.5,' for layer #',0p,i3)")exp(mean),sigma,layer; call output(trim(os))
       endif

       if (get_plot(layer,hk_ang=.true.)) then
         write(os,"('Hk angle histogram from nominal for layer #',0p,i3)")layer
         do i=1,num_grain; dt(i)=acos(dot_product((/p%kx,p%ky,p%kz/),uni_axis(:,i)))*180.d0/pi; enddo
         call plot_histogram(trim(os),20,dt)
         write(os,"('Hk polar angle histogram for layer #',0p,i3)")layer
         do i=1,num_grain; dt(i)=acos(uni_axis(3,i))*180.d0/pi; enddo
         call plot_histogram(trim(os),20,dt)
       endif

     endif
     ok_p=.true.
   endif
  end function

  function construct_tau(num_grain,grain_area,ilayer,isilent_p,tau) result (ok_p)

   use io_m, only: output
   use random_m, only: random_normal
   use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour

   integer,intent(in)::num_grain
   real(8),intent(in)::grain_area(num_grain)
   integer,intent(in),optional::ilayer
   logical,intent(in),optional::isilent_p

   real(8),pointer::tau(:)
   logical::silent_p,ok_p
   integer::layer
   type(grain_s),pointer::p
   integer::i
   real(8)::mean,x,y,pi
   character(200)::os

   ok_p=.false.

   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   nullify(tau)
   if (associated(p)) then
     pi=acos(-1.d0)

     write(os,"('  constructing tau distribution for layer #',0p,i3)")layer
     if (.not.silent_p) call output(trim(os))

     if (associated(tau)) then; if (size(tau).ne.num_grain) then; deallocate(tau); nullify(tau); endif; endif
     if (.not.associated(tau)) allocate(tau(num_grain))
     mean=p%tau_nom; if (mean.eq.0.d0) mean=1.d-12; tau=mean

     !add in gaussian distribution
     if (p%tau_sigma.ne.0.d0) then
       do i=1,num_grain
         x=1.d0+p%tau_sigma*random_normal(p%random_tau)
         do while (x.le.0.d0)
           x=1.d0+p%tau_sigma*random_normal(p%random_tau)
         enddo
         tau(i)=tau(i)*x
       enddo
     endif

     !compute min, max and volume average
     if (.not.silent_p) then
       x=-1.d30; y=1.d30
       do i=1,num_grain; x=max(x,tau(i)); y=min(y,tau(i)); enddo
       call indmed(tau,i)
       write(os,"('    minimum, maximum and median tau:',1p,3(e12.5,1x),'(s) for layer #',0p,i3)")y,x,tau(i),layer; call output(trim(os))
       call find_mean_sig(tau,x,y)
       write(os,"('    average tau:',1p,e12.5,' and sigma:',e12.5,' (s) for layer #',0p,i3)")x,y,layer; call output(trim(os))
       call find_mean_sig(tau,x,y,grain_area)
       write(os,"('    volume average tau:',1p,e12.5,' and sigma:',e12.5,' (s) for layer #',0p,i3)")x,y,layer; call output(trim(os))
       call find_mean_sig(log(tau),x,y,grain_area)
       write(os,"('    ln-normal tau exp(mu): '1p,e12.5,' (s), sigma of tau:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
       call output(trim(os))
     endif

     if (get_plot(layer,tau=.true.)) then
       write(os,"('tau histogram layer #',0p,i3)")layer
       if (get_plot(layer,tau_hist=.true.)) call plot_histogram(trim(os),20,tau*1.d9)
       write(os,"('tau vs. area for layer #',0p,i3)")layer
       if (get_plot(layer,tau_vol=.true.)) call plot_scatter(trim(os),grain_area*1.d14,tau*1.d9,.true.)
       write(os,"('tau for layer #',0p,i3)")layer
       call plot_contour(trim(os),tau*1.e9,.true.,layer)
     endif

     ok_p=.true.
   endif

  end function

  function construct_tc(num_grain,grain_vol,ilayer,isilent_p,tc) result (ok_p)
! function construct_tc(num_grain,grain_area,ilayer,isilent_p,tc) result (ok_p)

   use io_m, only: output
   use random_m, only: random_normal
   use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour

   integer,intent(in)::num_grain
   real(8),intent(in)::grain_vol(num_grain)
!  real(8),intent(in)::grain_area(num_grain)
   integer,intent(in),optional::ilayer
   logical,intent(in),optional::isilent_p

   real(8),pointer::tc(:)
   logical::silent_p,ok_p
   integer::layer
   type(grain_s),pointer::p
   integer::i
   real(8)::mean,x,y,z,pi
   character(200)::os
   real(8),parameter::tc_min=200.d0  !TODO fix me

   ok_p=.false.

   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   nullify(tc)
   if (associated(p)) then
     pi=acos(-1.d0)

     write(os,"('  constructing Tc distribution for layer #',0p,i3)")layer
     if (.not.silent_p) call output(trim(os))

     if (associated(tc)) then; if (size(tc).ne.num_grain) then; deallocate(tc); nullify(tc); endif; endif
     if (.not.associated(tc)) allocate(tc(num_grain))
     mean=p%tc_nom; if (mean.eq.0.d0) mean=1.d22; tc=mean

     !do scaling due to size
     if (p%tc_eta.ne.0.d0) then
       y=-1.d0/p%tc_eta
       do i=1,num_grain
!        x=2.d0*sqrt(grain_area(i)/pi)*1.d7
         x=(6.d0*grain_vol(i)*abs(p%vol_scale)/pi)**(1.d0/3.d0)*1.d7
         tc(i)=max( tc_min, min( p%tc_max, mean * max(0.d0, 1.d0 - (x/p%tc_d0)**(y))))
       enddo
     endif

     !add in gaussian distribution
     if (p%tc_sigma.ne.0.d0) then
       do i=1,num_grain
         x=tc(i)*(1.d0+p%tc_sigma*random_normal(p%random_tc))
         do while (x.gt.p%tc_max)
           x=tc(i)*(1.d0+p%tc_sigma*random_normal(p%random_tc))
         enddo
         tc(i)=max(tc_min, min(p%tc_max, x))
       enddo
     endif

     if (mean.lt.1.d20) then
     !compute min, max and volume average
     if (.not.silent_p) then
       x=-1.d30; y=1.d30
       do i=1,num_grain; x=max(x,tc(i)); y=min(y,tc(i)); enddo
       call indmed(tc,i)
       write(os,"('    (before rescaling) minimum, maximum and median Tc:',1p,3(e12.5,1x),'(K) for layer #',0p,i3)")y,x,tc(i),layer; call output(trim(os))
       call find_mean_sig(tc,x,y)
       write(os,"('    (before rescaling) average Tc:',1p,e12.5,' and sigma:',e12.5,' (K) for layer #',0p,i3)")x,y,layer; call output(trim(os))
!      call find_mean_sig(tc,mean,y,grain_area); z=mean*sum(grain_area)
       call find_mean_sig(tc,mean,y,grain_vol); z=mean*sum(grain_vol)
       write(os,"('    (before rescaling) volume average Tc:',1p,e12.5,' and sigma:',e12.5,' (K) (',e12.4,' %) for layer #',0p,i3)")mean,y,y/mean*100.d0,layer; call output(trim(os))
!      call find_mean_sig(log(tc),x,y,grain_area)
       call find_mean_sig(log(tc),x,y,grain_vol)
       write(os,"('    (before rescaling) ln-normal Tc exp(mu): '1p,e12.5,', sigma of Tc:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
       call output(trim(os))
     else
!      call find_mean_sig(tc,mean,y,grain_area); z=mean*sum(grain_area)
       call find_mean_sig(tc,mean,y,grain_vol); z=mean*sum(grain_vol)
     endif

     !do we need to rescale?
     if (z.ne.0.d0) then
!      z = p%tc_nom * sum(grain_area) / z
       z = p%tc_nom * sum(grain_vol) / z
       if (.not. p%tc_rescale) z=1.d0
       do i=1,num_grain
         tc(i)=max( tc_min, min(p%tc_max,tc(i)*z))
       enddo
     endif

!do i=1,num_grain;write(86,"(i6,2(x,e12.5))")i,tc(i),grain_vol(i)*1e21;enddo
       if (dist_io.ne.0) then
         open(dist_io,file='tc.dat',form='formatted',status='unknown')
         call plot_histogram('tc.dat',size(tc),tc,ionum=dist_io)
         close(dist_io)
       endif
     if (.not.silent_p) then
       x=-1.d30; y=1.d30
       do i=1,num_grain; x=max(x,tc(i)); y=min(y,tc(i)); enddo
       call indmed(tc,i)
       write(os,"('    minimum, maximum and median Tc:',1p,3(e12.5,1x),'(K) for layer #',0p,i3)")y,x,tc(i),layer; call output(trim(os))
       call find_mean_sig(tc,x,y)
       write(os,"('    average Tc:',1p,e12.5,' and sigma:',e12.5,' (K) for layer #',0p,i3)")x,y,layer; call output(trim(os))
!      call find_mean_sig(tc,x,y,grain_area); mean = x; z=mean*sum(grain_area)
       call find_mean_sig(tc,x,y,grain_vol); mean = x; z=mean*sum(grain_vol)
       write(os,"('    volume average Tc:',1p,e12.5,' and sigma:',e12.5,' (K) (',e12.4,' %) for layer #',0p,i3)")mean,y,y/mean*100.d0,layer; call output(trim(os))
!      call find_mean_sig(log(tc),x,y,grain_area)
       call find_mean_sig(log(tc),x,y,grain_vol)
       write(os,"('    ln-normal Tc exp(mu): '1p,e12.5,' (K), sigma of Tc:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
       call output(trim(os))
     endif

     if (get_plot(layer,tc=.true.).and.minval(tc).lt.1.d4) then
       write(os,"('Tc histogram layer #',0p,i3)")layer
       if (get_plot(layer,tc_hist=.true.)) call plot_histogram(trim(os),20,tc)
!      write(os,"('Tc vs. area for layer #',0p,i3)")layer
!      if (get_plot(layer,tc_vol=.true.)) call plot_scatter(trim(os),grain_area*1.d14,tc,.true.)
       write(os,"('Tc vs. grain volume for layer #',0p,i3)")layer
       if (get_plot(layer,tc_vol=.true.)) call plot_scatter(trim(os),grain_vol*1.d21,tc,.true.)
       write(os,"('Tc for layer #',0p,i3)")layer
       call plot_contour(trim(os),tc,.true.,layer)
     endif
     endif

     ok_p=.true.
   endif

  end function

  function construct_ms(num_grain,grain_area,ilayer,isilent_p,ms) result (ok_p)

   use io_m, only: output
   use random_m, only: random_normal
   use plot_m, only: get_plot, plot_histogram, plot_scatter, plot_contour

   integer,intent(in)::num_grain
   real(8),intent(in)::grain_area(num_grain)
   integer,intent(in),optional::ilayer
   logical,intent(in),optional::isilent_p

   real(8),pointer::ms(:)
   logical::silent_p,ok_p
   integer::layer
   type(grain_s),pointer::p
   integer::i
   real(8)::mean,x,y,pi
   character(200)::os

   ok_p=.false.

   silent_p=.false.; if (present(isilent_p)) silent_p=isilent_p
   layer=1;if (present(ilayer)) layer=max(1,ilayer)
   nullify(p); if (layer.le.size(grain)) p=>grain(layer)%p

   nullify(ms)
   if (associated(p)) then
     pi=acos(-1.d0)

     write(os,"('  constructing Ms distribution for layer #',0p,i3)")layer
     if (.not.silent_p) call output(trim(os))

     if (associated(ms)) then; if (size(ms).ne.num_grain) then; deallocate(ms); nullify(ms); endif; endif
     if (.not.associated(ms)) allocate(ms(num_grain))
     mean=p%ms_nom; if (mean.eq.0.d0) mean=800.d0; ms=mean

     !add in gaussian distribution
     if (p%ms_sigma.ne.0.d0) then
       do i=1,num_grain
         x=ms(i)*(1.d0+p%ms_sigma*random_normal(p%random_ms))
         do while (x.gt.p%ms_max)
           x=ms(i)*(1.d0+p%ms_sigma*random_normal(p%random_ms))
         enddo
         ms(i)=x
       enddo
     endif

     !compute min, max and volume average
     if (.not.silent_p) then
       x=-1.d30; y=1.d30
       do i=1,num_grain
         x=max(x,ms(i))
         y=min(y,ms(i))
       enddo
       call indmed(ms,i)
       write(os,"('    minimum, maximum and median Ms:',1p,3(e12.5,1x),'(emu/cc) for layer #',0p,i3)")y,x,ms(i),layer; call output(trim(os))
       call find_mean_sig(ms,x,y)
       write(os,"('    average Ms:',1p,e12.5,' and sigma:',e12.5,' (emu/cc) for layer #',0p,i3)")x,y,layer; call output(trim(os))
       call find_mean_sig(ms,x,y,grain_area)
       write(os,"('    volume average Ms:',1p,e12.5,' and sigma:',e12.5,' (emu/cc) for layer #',0p,i3)")x,y,layer; call output(trim(os))
       call find_mean_sig(log(ms),x,y,grain_area)
       write(os,"('    ln-normal exp(mu): '1p,e12.5,' (emu/cc), sigma of Ms:',1p,e12.5,' for layer #',0p,i3)")exp(x),y,layer
       call output(trim(os))
     endif

     if (get_plot(layer,ms=.true.).and.p%ms_sigma.ne.0.d0) then
       write(os,"('Ms histogram layer #',0p,i3)")layer
       call plot_histogram(trim(os),20,ms)
       write(os,"('Ms vs. area for layer #',0p,i3)")layer
       call plot_scatter(trim(os),grain_area*1.d14,ms,.true.)
       write(os,"('Ms(T=0K) for layer #',0p,i3)")layer
       call plot_contour(trim(os),ms,.true.,layer)
     endif

     ok_p=.true.
   endif

  end function

  subroutine find_mean_sig(x,m,s,w)
   real(8),dimension(:),intent(in)::x
   real(8),dimension(:),intent(in),optional::w
   real(8),intent(out)::m,s
   real(8)::ww
   integer n
   m=0.d0; s=0.d0; ww=0.d0
   if (present(w)) then
     do n=1,size(x)
       s=s+ww*w(n)*(x(n)-m)**2/(w(n)+ww)
       ww=ww+w(n)
       m=m+(x(n)-m)*w(n)/ww
     enddo
     s=sqrt(s*size(x)/(ww*max(1,size(x)-1)))
   else
     do n=1,size(x)
       s=s+(n-1)*(x(n)-m)**2/n
       m=m+(x(n)-m)/n
     enddo
     s=sqrt(s/max(1,size(x)-1))
   endif
  end subroutine


  function set_construct_grain(layer,grain_p,final_packing_fraction,packing_fraction,small_packing_fraction,initial_grain_area,initial_grain_ext_area, &
                small_grain_area,grain_sigma,small_grain_sigma,average_grain_area,dr_grain, & 
                kx,ky,kz,k1,k2,k_ang,k_cone_sigma,kmax,hkmin,k_frac_inplane,k_frac_random,k_sigma,k_3d,k_2d,k_uniform,k_cone,k_cubic,k_ac,k_scale,k_rescale, &
                tc_nom,tc_sigma,tc_max,tc_d0,tc_eta,tc_rescale,random_seed,ms_nom,ms_sigma,ms_max,ex,ex_sigma,ex_upper,ex_upper_sigma, &
                same_grain_p,same_hk_p,same_exchange_p,same_tc_p,same_ms_p,k_scale_al,tau_nom,tau_sigma,num_end,vol_scale, &
                k_sigma_small,k1_small,k2_small,frac_hk_small,old_small_grain,grow_small_grain,initial_small_grain_area, &
                rand_mat,print_dist_p) result (ok_p)

   use random_m, only: random_copy, random_destroy
   use io_m
   integer,intent(in)::layer
   logical,intent(in),optional::grain_p,same_grain_p,same_hk_p,same_exchange_p,same_tc_p,same_ms_p,print_dist_p
   real(8),intent(in),optional::packing_fraction,small_packing_fraction,initial_grain_area,initial_grain_ext_area,small_grain_area,grain_sigma,small_grain_sigma,initial_small_grain_area,final_packing_fraction
   real(8),intent(in),optional::kx,ky,kz,k1,k2,average_grain_area,dr_grain,vol_scale
   real(8),intent(in),optional::k_ang,k_cone_sigma,kmax,k_frac_inplane,k_frac_random,k_sigma,k_ac,hkmin
   logical,intent(in),optional::k_3d,k_2d,k_uniform,k_cone,k_cubic,k_scale,k_rescale,k_scale_al
   real(8),intent(in),optional::tc_nom,tc_sigma,tc_max,tc_d0,tc_eta,ms_nom,ms_sigma,ms_max,ex,ex_sigma,ex_upper,ex_upper_sigma,tau_nom,tau_sigma
   logical,intent(in),optional::tc_rescale,old_small_grain,grow_small_grain
   integer,intent(in),optional::random_seed,num_end,rand_mat
   real(8),intent(in),optional::k_sigma_small,k1_small,k2_small,frac_hk_small

   integer::rand_seed,rand_mult

   logical::ok_p
    character(200)::str


   rand_seed=34; if (present(random_seed)) rand_seed=random_seed;
   rand_seed=rand_seed*(layer+1)
   rand_mult=89032211; if (present(rand_mat)) rand_mult=rand_mat;
   ok_p=.true.; call create_grain_ss(layer)
   if (present(grain_p)) grain(layer)%p%grain=grain_p
   if (present(num_end)) grain(layer)%p%num_end=num_end
   if (present(vol_scale)) grain(layer)%p%vol_scale=vol_scale
   if (present(packing_fraction)) grain(layer)%p%packing_fraction=packing_fraction
   if (present(final_packing_fraction)) grain(layer)%p%final_packing_fraction=final_packing_fraction
   if (present(small_packing_fraction)) grain(layer)%p%small_packing_fraction=small_packing_fraction
   if (present(initial_grain_area)) grain(layer)%p%initial_grain_area=initial_grain_area
   if (present(average_grain_area)) grain(layer)%p%average_grain_area=average_grain_area
   if (present(dr_grain)) grain(layer)%p%dr_grain=dr_grain
   if (present(initial_grain_ext_area)) grain(layer)%p%initial_grain_ext_area=initial_grain_ext_area
   if (present(small_grain_area)) grain(layer)%p%small_grain_area=small_grain_area
   if (present(initial_small_grain_area)) grain(layer)%p%initial_small_grain_area=initial_small_grain_area
   if (present(old_small_grain)) grain(layer)%p%old_small_grain=old_small_grain
   if (present(grow_small_grain)) grain(layer)%p%grow_small_grain=grow_small_grain
   if (present(grain_sigma)) grain(layer)%p%grain_sigma=grain_sigma
   if (present(small_grain_sigma)) grain(layer)%p%small_grain_sigma=small_grain_sigma
   if (present(kx)) grain(layer)%p%kx=kx
   if (present(ky)) grain(layer)%p%ky=ky
   if (present(kz)) grain(layer)%p%kz=kz
   if (present(k1)) grain(layer)%p%k1=k1
   if (present(k2)) grain(layer)%p%k2=k2
   if (present(k1_small)) grain(layer)%p%k1_small=k1_small
   if (present(k2_small)) grain(layer)%p%k2_small=k2_small
   if (present(k_sigma_small)) grain(layer)%p%k_sigma_small=k_sigma_small
   if (present(frac_hk_small)) grain(layer)%p%frac_hk_small=frac_hk_small
   if (present(k_ang)) grain(layer)%p%k_ang=k_ang
   if (present(k_cone_sigma)) grain(layer)%p%k_cone_sigma=k_cone_sigma
   if (present(kmax)) grain(layer)%p%kmax=kmax
   if (present(hkmin)) grain(layer)%p%hkmin=hkmin
   if (present(k_frac_inplane)) grain(layer)%p%k_frac_inplane=k_frac_inplane
   if (present(k_frac_random)) grain(layer)%p%k_frac_random=k_frac_random
   if (present(k_sigma)) grain(layer)%p%k_sigma=k_sigma
   if (present(k_3d)) grain(layer)%p%k_3d=k_3d
   if (present(k_2d)) grain(layer)%p%k_2d=k_2d
   if (present(k_uniform)) grain(layer)%p%k_uniform=k_uniform
   if (present(k_cone)) grain(layer)%p%k_cone=k_cone
   if (present(k_cubic)) grain(layer)%p%k_cubic=k_cubic
   if (present(k_scale)) grain(layer)%p%k_scale=k_scale
   if (present(k_scale_al)) grain(layer)%p%k_scale_al=k_scale_al
   if (present(k_ac)) grain(layer)%p%k_ac=k_ac
   if (present(k_rescale)) grain(layer)%p%k_rescale=k_rescale
   if (present(tc_nom)) grain(layer)%p%tc_nom=tc_nom
   if (present(tc_sigma)) grain(layer)%p%tc_sigma=tc_sigma
   if (present(tc_max)) grain(layer)%p%tc_max=tc_max
   if (present(tc_d0)) grain(layer)%p%tc_d0=tc_d0
   if (present(tc_eta)) grain(layer)%p%tc_eta=tc_eta
   if (present(tc_rescale)) grain(layer)%p%tc_rescale=tc_rescale
   if (present(tau_nom)) grain(layer)%p%tau_nom=tau_nom
   if (present(tau_sigma)) grain(layer)%p%tau_sigma=tau_sigma
   if (present(ms_nom)) grain(layer)%p%ms_nom=ms_nom
   if (present(ms_sigma)) grain(layer)%p%ms_sigma=ms_sigma
   if (present(ms_max)) grain(layer)%p%ms_max=ms_max
   if (present(ex)) grain(layer)%p%ex=ex
   if (present(ex_sigma)) grain(layer)%p%ex_sigma=ex_sigma
   if (present(ex_upper)) grain(layer)%p%ex_upper=ex_upper
   if (present(ex_upper_sigma)) grain(layer)%p%ex_upper_sigma=ex_upper_sigma

   write(str,"(i5)") layer;str=adjustl(str)
   if (present(same_grain_p)) then
     if (same_grain_p.and.layer.ne.1) then
        call random_destroy(grain(layer)%p%random_grain);  grain(layer)%p%random_grain=>random_copy('random_grain'//trim(str),grain(layer-1)%p%random_grain)
        call random_destroy(grain(layer)%p%random_grains); grain(layer)%p%random_grains=>random_copy('random_grains'//trim(str),grain(layer-1)%p%random_grains)
     endif
   endif
   if (present(same_hk_p)) then
     if (same_hk_p.and.layer.ne.1) then
        call random_destroy(grain(layer)%p%random_hk_mag); grain(layer)%p%random_hk_mag=>random_copy('random_hk_mag'//trim(str),grain(layer-1)%p%random_hk_mag)
        call random_destroy(grain(layer)%p%random_hk_mag_small); 
                            grain(layer)%p%random_hk_mag_small=>random_copy('random_hk_mag_small'//trim(str),grain(layer-1)%p%random_hk_mag_small)
        call random_destroy(grain(layer)%p%random_hk_uni_angles); 
                            grain(layer)%p%random_hk_uni_angles=>random_copy('random_hk_uni_angles'//trim(str),grain(layer-1)%p%random_hk_uni_angles)
        call random_destroy(grain(layer)%p%random_hk_frac_random); 
                            grain(layer)%p%random_hk_frac_random=>random_copy('random_hk_frac_random'//trim(str),grain(layer-1)%p%random_hk_frac_random)
        call random_destroy(grain(layer)%p%random_hk_frac_inplane);
                            grain(layer)%p%random_hk_frac_inplane=>random_copy('random_hk_frac_inplane'//trim(str),grain(layer-1)%p%random_hk_frac_inplane)
        call random_destroy(grain(layer)%p%random_hk_cubic_angles);
                            grain(layer)%p%random_hk_cubic_angles=>random_copy('random_hk_cubic_angles'//trim(str),grain(layer-1)%p%random_hk_cubic_angles)
     endif
   endif
   if (present(same_exchange_p)) then
     if (same_exchange_p.and.layer.ne.1) then
        call random_destroy(grain(layer)%p%random_ex); grain(layer)%p%random_ex=>random_copy('random_ex'//trim(str),grain(layer-1)%p%random_ex)
        call random_destroy(grain(layer)%p%random_exv); grain(layer)%p%random_exv=>random_copy('random_exv'//trim(str),grain(layer-1)%p%random_exv)
     endif
   endif
   if (present(same_tc_p)) then
     if (same_tc_p.and.layer.ne.1) then
        call random_destroy(grain(layer)%p%random_tc); grain(layer)%p%random_tc=>random_copy('random_tc'//trim(str),grain(layer-1)%p%random_tc)
     endif
   endif
   if (present(same_ms_p)) then
     if (same_ms_p.and.layer.ne.1) then
        call random_destroy(grain(layer)%p%random_ms); grain(layer)%p%random_ms=>random_copy('random_ms'//trim(str),grain(layer-1)%p%random_ms)
     endif
   endif

   if (present(print_dist_p)) then
     if (print_dist_p.and.dist_io.eq.0) then
       dist_io=add_file('temp')
     endif
   endif

  contains

   subroutine create_grain_ss(layer)
    use random_m, only: random_s, random_create, random_copy, random_destroy
    integer,intent(in)::layer
    type(grain_ss),pointer::new_grain(:)
    type(grain_s),pointer::p
    integer::n,n_l

    n_l=0; if (associated(grain)) n_l=size(grain)

    if (layer.gt.n_l) then
      allocate(new_grain(layer))
      do n=1,n_l
        new_grain(n)%p=>grain(n)%p
      enddo
      do n=n_l+1,layer
        allocate(p)
        new_grain(n)%p=>p
        call create_random_generators(new_grain(n)%p,rand_seed,rand_mult,n)
      enddo
      if (associated(grain)) deallocate(grain)
      grain=>new_grain
    endif
   end subroutine

   subroutine destroy_random_generators(p)
    use random_m, only: random_destroy
    type(grain_s),pointer::p

    if (associated(p%random_grain)) call random_destroy(p%random_grain)
    if (associated(p%random_grains)) call random_destroy(p%random_grains)
    if (associated(p%random_hk_mag)) call random_destroy(p%random_hk_mag)
    if (associated(p%random_hk_mag_small)) call random_destroy(p%random_hk_mag_small)
    if (associated(p%random_hk_uni_angles)) call random_destroy(p%random_hk_uni_angles)
    if (associated(p%random_hk_frac_random)) call random_destroy(p%random_hk_frac_random)
    if (associated(p%random_hk_frac_inplane)) call random_destroy(p%random_hk_frac_inplane)
    if (associated(p%random_hk_cubic_angles)) call random_destroy(p%random_hk_cubic_angles)
    if (associated(p%random_tc)) call random_destroy(p%random_tc)
    if (associated(p%random_tau)) call random_destroy(p%random_tau)
    if (associated(p%random_ms)) call random_destroy(p%random_ms)
    if (associated(p%random_ex)) call random_destroy(p%random_ex)
    if (associated(p%random_exv)) call random_destroy(p%random_exv)
   end subroutine

   subroutine create_random_generators(p,r,m,n)
    use random_m, only: random_s, random_create, random_destroy, random_uniform
    integer,intent(in)::r,n,m
    type(grain_s),pointer::p
    type(random_s),pointer::rand
    character(200)::str

    write(str,"(i5)") n;str=adjustl(str)

    call destroy_random_generators(p)
    rand=>random_create('temp',r)
    p%random_grain=>random_create('random_grain'//trim(str),nint(89032211*random_uniform(rand)))
    p%random_grains=>random_create('random_grains'//trim(str),nint(89032211*random_uniform(rand)))
    p%random_hk_mag=>random_create('random_hk_mag'//trim(str),nint(m*random_uniform(rand)))
    p%random_hk_uni_angles=>random_create('random_hk_uni_angles'//trim(str),nint(m*random_uniform(rand)))
    p%random_hk_frac_random=>random_create('random_hk_frac_random'//trim(str),nint(m*random_uniform(rand)))
    p%random_hk_frac_inplane=>random_create('random_hk_frac_inplane'//trim(str),nint(m*random_uniform(rand)))
    p%random_hk_cubic_angles=>random_create('random_hk_cubic_angles'//trim(str),nint(m*random_uniform(rand)))
    p%random_tc=>random_create('random_tc'//trim(str),nint(m*random_uniform(rand)))
    p%random_ms=>random_create('random_ms'//trim(str),nint(m*random_uniform(rand)))
    p%random_ex=>random_create('random_ex'//trim(str),nint(m*random_uniform(rand)))
    p%random_exv=>random_create('random_exv'//trim(str),nint(m*random_uniform(rand)))
    p%random_tau=>random_create('random_tau'//trim(str),nint(m*random_uniform(rand)))
    p%random_hk_mag_small=>random_create('random_hk_mag_small'//trim(str),nint(m*random_uniform(rand)))
    call random_destroy(rand)
   end subroutine
  end function

end module construct_grain_m
