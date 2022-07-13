! this file contains module to read in main input file, do initialization and
! return cli file name to the calling routine

module input_m
  implicit none

  private
  public::get_input

contains


  function get_input(in_file0,clifile) result(ok_p)
 !this function reads in the main input file
  use io_m
  use plot_m, only: set_plot
  use data_m, only: set_geometry, set_data
  use tensor_m, only: set_tensor
  use off_field_m, only: set_off_field_tensor
  use demag_m, only: demag_window, start_demag
  use construct_grain_m, only: set_construct_grain
  use temp_scaling_m, only:set_temp_scaling
  use applied_m, only: start_applied
  use mag_m, only: start_mag, set_mag
  use field_m, only: start_field
  use llg_m, only: start_llg, set_llg
  use decay_m, only: set_decay_mag
  use checkpoint_m, only: checkpoint
  use fftw_m, only: fft_quick_plan

   character(*),intent(in)::in_file0
   character(200)::in_file
   character(200),intent(inout)::clifile
   logical::ok_p

   logical::self_fld_only_p,self_fld_grain_only_p,hk_3d_random_p,hk_2d_random_p,hk_uniform_p,thermal_p,boltz_p, &
            fixed_dt_p,plot_p,write_plot,plot_hk_p,plot_grain_p,plot_fld_p,plot_hx_p,plot_ms_p,plot_mag_p,plot_mz_p,diff_m,test_rot_p, &
            plot_mh_p,plot_heff_p,plot_signed_heff_p, plot_tc_p, plot_tau_p,plot_wfx_p,plot_wfy_p,plot_m_p, &
            plot_wfxf_p,plot_wfyf_p,plot_wfzf_p,plot_tempf_p,plot_mx_p,plot_my_p, &
            save_screen_p,append_screen_p,oldmag_write_p,vertical_same_hk, &
            oldmag_read_p,ac_p,saturate_p,oldgrain_read_p,oldgrain_write_p,old_exch_p, print_dist_p, &
            mt,mt_bs,zhang,llb,fp,plot_hscale_p,new_integrator,win_close_p,plot_temp_p,plot_hamr_scales_p, &
            hk_cone_p, compute_tensor_only_p, demag_p, plot_hx_il_p, hk_cubic_p, scale_hk, rescale_hk, tc_rescale_p, plot_scatter_p, plot_scatter_len_p, &
            fft_selfdemag_p,same_grain_p,same_hk_p,same_exchange_p,same_init_mag_p,same_tc_p, same_ms_p, scale_hk_al, use_mag_p, &
            zhang_l, garanin_tau, old_small_grain,grow_small_grain, quiet_p, fatalio_p, quick_fft_plan
   integer::nx,ny,num_images,rand_seed,thermal_seed,iter,min_iter,plot_iter,thermal_num_ave, &
            num_level_mag,compute_tensor_layer,compute_tensor_layer_obs,llg_type,nxfft,nyfft,media_layer,num_end,rand_mat,vertical_spins
   real(8)::dx,dy,dz,interlayer_space,sul_mu,fly_height,head_mu,ms,ms_sigma,ms_max,k_ang,kx,ky,kz, &
            k_sigma,hk,hk2,hx,hx_sigma,final_packing_fraction,packing_fraction,vol_sigma,dr_grain, &
            average_grain_area,alpha,gamma,dt,max_dm,min_dm,dt_frac,dt_max, &
            t_max,t_min,del_t,conv_torque,conv_angle,temperature,initial_grain_area, &
            mx,my,mz,xmin,xmax,ymin,ymax,min_ave_m(3),max_ave_m(3),hk_cone_sigma,tc,tc_max,tau,tau_sigma,heatsink_dt, &
            max_k, frac_in_plane, frac_random, tc_d0, tc_eta, hk_ac, min_hk, vol_scale, &
            tc_sigma,init_vol_sigma,initial_ext_grain_area,average_small_grain_area,small_packing_fraction,small_grain_sigma,initial_small_grain_area, &
            k_sigma_small,hk_small,hk2_small,frac_hk_small, &
            offset_xfft,offset_yfft,spin,ex,ex_sigma,ex_upper,ex_upper_sigma,space_to_upper,hkexp,msexp,vertical_a
   integer::ntraj,traj
   real(8)::torbit,vol_frac,max_eb

   character(200)::cli_file,screen_file,demag_file,loop_file

   integer::i,num_layers,ionum,prev_rand_seed,j
   real(8)::x,y
   character(200)::os
   real(8),allocatable::dzl(:),il_space(:)
   integer,allocatable::thm_seed(:)

  !new variables needed for multiple layers

  !see below for documentation
   namelist /decay/ ntraj,traj,torbit,vol_frac,max_eb,iter,boltz_p,test_rot_p,diff_m
   namelist /grid/ nx,ny,dx,dy,dz,nxfft,nyfft,offset_xfft,offset_yfft,fft_selfdemag_p
   namelist /demag/ interlayer_space,sul_mu,fly_height,head_mu,num_images,self_fld_only_p,self_fld_grain_only_p, &
      demag_file,compute_tensor_only_p,compute_tensor_layer,demag_p,compute_tensor_layer_obs,quick_fft_plan
   namelist /media/ ms,ms_sigma,ms_max,k_ang,kx,ky,kz,hk_cubic_p, hk_3d_random_p,hk_2d_random_p,hk_uniform_p, &
      k_sigma,hk,hk2,hx,hx_sigma,final_packing_fraction,packing_fraction,vol_sigma,dr_grain, &
      k_sigma_small,hk_small,hk2_small,frac_hk_small,vertical_a,vertical_spins,vertical_same_hk, &
      average_grain_area,rand_seed,initial_grain_area,ac_p,saturate_p,mx,my,mz,old_exch_p, &
      hk_cone_p,hk_cone_sigma,tc,tc_sigma,tc_d0,tc_eta,tc_max,tc_rescale_p,tau,tau_sigma,heatsink_dt,  &
      dz,ex,ex_sigma,ex_upper,ex_upper_sigma,same_grain_p,same_hk_p,same_exchange_p, &
      init_vol_sigma,initial_ext_grain_area,max_k,min_hk,scale_hk,rescale_hk,scale_hk_al,hk_ac,frac_in_plane,frac_random, &
      spin, space_to_upper,same_tc_p,same_ms_p,same_init_mag_p,loop_file, average_small_grain_area,small_packing_fraction,small_grain_sigma, &
      old_small_grain,grow_small_grain,initial_small_grain_area,vol_scale,num_end,rand_mat,print_dist_p
   namelist /llg/ alpha,gamma,dt,thermal_p,thermal_seed,max_dm,min_dm,dt_frac,dt_max, &
      t_max,t_min,del_t,conv_torque,conv_angle,temperature,fixed_dt_p,iter,min_iter,min_ave_m, &
      max_ave_m,plot_iter,thermal_num_ave,new_integrator,llg_type, &
      media_layer,mt,mt_bs,zhang,llb,fp,use_mag_p,hkexp,msexp,zhang_l,garanin_tau
   namelist /plot/plot_p,write_plot,plot_hk_p,plot_grain_p,plot_fld_p,plot_hx_p,plot_ms_p,plot_mag_p,plot_mz_p,plot_mh_p, &
      plot_heff_p,plot_signed_heff_p, plot_hscale_p,num_level_mag,plot_temp_p,plot_hamr_scales_p, &
      plot_wfxf_p,plot_wfyf_p,plot_wfzf_p,plot_tempf_p,plot_mx_p,plot_my_p, &
      plot_tc_p,plot_tau_p,media_layer,plot_hx_il_p,plot_scatter_p ,plot_scatter_len_p,plot_wfx_p,plot_wfy_p,plot_m_p
   namelist /control/cli_file,save_screen_p,append_screen_p,screen_file,oldmag_write_p, &
      oldmag_read_p,oldgrain_read_p,oldgrain_write_p,win_close_p,loop_file,quiet_p,fatalio_p
   namelist /pattern/mx,my,mz,ac_p,saturate_p,xmin,ymin,xmax,ymax,media_layer


   in_file=adjustl(in_file0)
   if (in_file(1:1).eq.' ') in_file='som.inp'
   
   inquire(file=trim(in_file),exist=ok_p)
   if (.not.ok_p) then     !no new input file. does old one exist?
     write(os,"('ERROR- no input files: ',a)") trim(in_file); call output(trim(os)); return
   endif

   write(os,"('reading namelists from ',a)")trim(in_file); call output(trim(os))
   ok_p = .false.

   ionum=add_file(trim(in_file))
   open(ionum,file=trim(in_file),status='old',action='read',position='rewind')

  !control parameters
   call output(' reading `control` namelist')

   clifile=adjustl(clifile); cli_file=trim(clifile) !cli file name {'som.cli'}
   if (cli_file(1:1).eq.' ') cli_file='som.cli'
   screen_file='scrn.out'          !screen file name 
   quiet_p=.true.
   save_screen_p=.true.           !save screen prints to screen_file? 
   append_screen_p=.false.         !append to screen_file?
   oldmag_read_p=.false.           !read media.f magnetization files?
   oldmag_write_p=.false.          !write media.f magnetization files?
   oldgrain_read_p=.false.         !read media.f grain files?
   oldgrain_write_p=.false.        !write media.f grain files? (not recommended)
   win_close_p=.false.             !should we exit directly or pause at the end
   loop_file='h.dat'; loop_file=adjustl(loop_file)
   fatalio_p=.false.
    rewind(ionum)
    read(ionum,control,iostat=i)
    if (i.gt.0) then
      write(os,"(' ERROR reading namelist `control` in ',a)")trim(in_file); call output(trim(os))
      close(ionum); ok_p=.false.
      return
    endif

   if (clifile(1:1).eq.' ') clifile=trim(adjustl(cli_file));    !only update clifile if user did not pass in on command line
   append_screen_p=checkpoint(present_p=.true.)    !append if checkpoint files exist, regardless of setting
   call set_output(screen_file,save_screen_p,append_screen_p,quiet_p,fatalio_p)
   call set_data(0,trim(loop_file))

  !grid parameters
   rewind(ionum)
   call output(' reading `grid` namelist')
   nx=1   !number of sub-grains along x (GM's z) direction
   ny=1   !number of sub-grains along y (GM's x) direction
   dx=0.d0 !size of sub-grains along x (cm)
   dy=0.d0 !size of sub-grains along y (cm)
   dz=0.d0 !media thickness (cm)
   nxfft=0; nyfft=0; fft_selfdemag_p=.false.
   offset_xfft=0.d0; offset_yfft=0.d0; 
    read(ionum,grid,iostat=i)
    if (i.gt.0) then
      write(os,"(' ERROR reading namelist `grid` in ',a)")trim(in_file); call output(trim(os));
      close(ionum); ok_p=.false.
      return
    endif
   nx=abs(nx); ny=abs(ny); dx=abs(dx); dy=abs(dy); dz=abs(dz)

   num_layers=0; rewind(ionum); vertical_spins=1;
   read(ionum,media,iostat=i); do while (i.eq.0); num_layers=num_layers+vertical_spins; vertical_spins=1; read(ionum,media,iostat=i); enddo;num_layers=max(1,num_layers)
   allocate(dzl(num_layers), il_space(num_layers), thm_seed(num_layers)); space_to_upper=0.d0
   
   rewind(ionum); j=0
   do i=1,num_layers
     if (j.lt.1) then
       vertical_spins=1; read(ionum,media); j=vertical_spins
       if (j.ne.1) then
         dzl(i)=dz/j; il_space(i)=0.d0
       else
         dzl(i)=dz;il_space(i)=max(0.d0,space_to_upper)
       endif
     else
        dzl(i)=dzl(i-1); il_space(i)=0.d0
     endif
     j=j-1

   enddo

!  rewind(ionum)
!  thermal_seed=34; thm_seed=thermal_seed
!  read(ionum,llg,iostat=i)
!  do while (i.eq.0)
!    if (media_layer.gt.0.and.media_layer.le.num_layers) then
!      thm_seed(media_layer)=thermal_seed
!    elseif (media_layer.eq.0) then
!      thm_seed=thermal_seed
!    endif
!    read(ionum,llg,iostat=i)
!  enddo

   j=0; rewind(ionum); thermal_seed=-34; thm_seed=0; media_layer=0
   read(ionum,llg,iostat=i); j=j+1
   do while (i.eq.0)
     if (media_layer.gt.0.and.media_layer.le.num_layers) then
       thm_seed(media_layer)=abs(thermal_seed)
     elseif (media_layer.eq.0) then
       thm_seed(min(num_layers,j))=abs(thermal_seed)
     endif
     thermal_seed=-34; media_layer=0; read(ionum,llg,iostat=i); j=j+1
   enddo
   i=maxval(thm_seed); if (i.le.0) i=34
   do j=1,num_layers
     if (thm_seed(j).le.0) thm_seed(j)=i+j*10000
   enddo
!print *,thm_seed

   ok_p=start_mag(num_layers); ok_p=start_field(thm_seed)
   ok_p=start_applied(nx,ny,dx,dy,num_layers);

   call set_geometry(nx,ny,dx,dy,num_layers,dzl,il_space); dz=dzl(1)
   if (nxfft.lt.1.or.nxfft.gt.nx) nxfft=nx
   if (nyfft.lt.1.or.nyfft.gt.ny) nyfft=ny

  !magnetostatic parameters
   call output(' reading `demag` namelist')
   interlayer_space=0.d0   !space between top of SUL and bottom of hard layer (cm)
   sul_mu=1.d0             !relative dc permeability of SUL
   fly_height=0.d0         !space between top of hard layer and bottom of head (cm) - NOT USED
   head_mu=1.d0            !relative dc permeability of head - NOT USED
   num_images=2           !number of image pairs if both SUL and head of permeability > 1 - NOT USED
   self_fld_only_p=.false. !only have self-magnetostatic fields (approximate for grains) - NOT USED
   self_fld_grain_only_p=.true.   !if self demag only, then do full multi-spin grain
   demag_file='mrm_N.dat'  !file name to store magnetostatic tensor - NOT USED
   compute_tensor_only_p=.false. !only compute magentostatic tensor and quit
   compute_tensor_layer=0  !just compute this layer only (if 0 or > number of layers, then all layers) - NOT USED
   compute_tensor_layer_obs=0  !just compute this observation layer only (if 0 or > number of layers, then all layers) - NOT USED
   quick_fft_plan=.false.

   demag_p=.true.          !do we calculate demag field?
   rewind(ionum)
!  open(ionum,file=trim(in_file),status='old')
    read(ionum,demag,iostat=i)
    if (i.gt.0) then
      write(os,"(' ERROR reading namelist `demag` in ',a)")trim(in_file)
      call output(trim(os))
      close(ionum); ok_p=.false.
      return
    endif

   if (quick_fft_plan) call fft_quick_plan()
   call set_tensor(max(1.d0,head_mu),max(1.d0,sul_mu),abs(fly_height),abs(interlayer_space),trim(demag_file),num_images)
   call set_off_field_tensor(max(1.d0,head_mu),max(1.d0,sul_mu),abs(fly_height),abs(interlayer_space),trim(demag_file),num_images)
   call demag_window(offset=(/ offset_xfft, offset_yfft /),nx=nxfft,ny=nyfft)
   if (compute_tensor_only_p) then
     ok_p=start_demag(nxfft,nyfft,num_layers,.true.,demag_p.or.self_fld_only_p); ok_p=.false.; return
   else !if (.not.demag_p) then
     ok_p=start_demag(nxfft,nyfft,num_layers,.false.,demag_p.or.self_fld_only_p)
   endif
   ok_p=start_field(self_fld_only_p=self_fld_only_p,self_fld_grain_only_p=self_fld_grain_only_p)

  !media parameters
       scale_hk=.false.   !scale Hk according to volume
       scale_hk_al=.false. !use Hk scaling according to A. Lyberatos
       vol_scale=1.d0
       num_end=2
       hk_ac = 0.d0
       rescale_hk=.false.
       max_k=6.d10        !maximum K1 for uniaxial
       min_hk=-1.d22
       frac_in_plane=0.d0 !fraction of grains w/ uniaxial Hk in-plane
       frac_random=0.d0   !fraction of grains w/ uniaxial Hk randomly orientated
       ms=800.d0  !nominal Ms of hard layer
       ms_sigma=0.d0 !ln-normal/gaussian ms sigma
       ms_max=1.d22
       spin=0.5d0 !S (or J) of material (i.e. total angular momentum, Ms(T=0) = N g S mu_b
       k_ang=0.d0 !parameter for angular dispersion for Hk (degrees)
       kx =0.d0   !x direction of nominal Hk direction
       ky =0.d0   !y direction of nominal Hk direction
       kz =0.d0   !z direction of nominal Hk direction
       hk_3d_random_p=.false.  !HK should be random on a sphere?
       hk_2d_random_p=.false.  !Hk should be random in 2d perpendicular to vector (kx, ky, kz)?
       hk_uniform_p=.false.    !Hk should be uniformly distributed between 0 and k_ang?
       hk_cone_p=.false.   !Hk is a cone around (kx,ky,kz) at nominal angle k_ang w/ gaussian spread k_cone_sigma
       hk_cone_sigma=0.d0  !gaussian spread for Hk cone distribution (degrees)
       k_sigma=0.d0  !ln-normal sigma for Hk magnitude
       hk_cubic_p=.false. !cubic Hk? otherwise uniaxial
       hk=0.d0       !Hk magnitude (first term) (Oe)
       hk2=0.d0      !Hk magnitude (second term) (Oe) 
       hx=0.d0       !nominal intergrainular exchange an average grain will feel in dc neighborhood (Oe)
       hx_sigma=0.d0 !explicity ln-normal dispersion of intergrain exchange coupling
       k_sigma_small=0.d0
       hk_small=0.d0
       hk2_small=0.d0
       frac_hk_small=0.d0
       packing_fraction=1.d0 !packing fraction (decimal, otherwise a percent if > 1)
       small_packing_fraction=0.d0
       vol_sigma=0.d0 !target ln-normal sigma for grain volume
       old_small_grain=.true.
       grow_small_grain=.false.
       print_dist_p=.false.
       small_grain_sigma=0.d0     !this is a gaussian sigma
       dr_grain=0.2d0 !rate of grain growth (sub-grain/iteration)
       average_grain_area=0.d0 !desired grain area (<0 or 1.d0 or dx*dy, then grains are subgrains
                             !                   (<1.d0, then (cm^2)
                             !                   (otherwise in terms of sub-grains^2)
       average_small_grain_area=0.d0
       initial_small_grain_area=0.d0
       initial_grain_area=0.d0 !used instead of volume_sigma (if > 0 and < 1, in (cm^2) otherwise in sub-grain^2)
       init_vol_sigma=0.d0        !initial ln-normal sigma of seeds
       initial_ext_grain_area=0.d0        !initial grain+boundary area
       ac_p=.false.    !should sample be AC erased?
       mx=0.d0        !nominal initial direction of x component of magnetization
       my=0.d0        !nominal initial direction of y component of magnetization
       mz=1.d0        !nominal initial direction of z component of magnetization
       saturate_p= .false. ! should we initialize along (/ mx,my,mz /) instead of easy axes?
       old_exch_p=.false.  !exchange between grain varied by one random number?  otherwise each face coupling is varied
       tau=1.d-9; tau_sigma=0.d0
       heatsink_dt=0.d0
       tc=-1.d0
       tc_max=1.d40; tc_rescale_p=.true. !max Tc and should we rescale Tc distribution so that volume averarged Tc = tc
       tc_sigma= 0.d0      ! sigma(Tc)/mean(Tc) in HAMR simulations
       tc_d0=0.90d0
       tc_eta=0.d0
       ex=0.d0        !in-plane intergrain exchange as energy/area
       ex_sigma=0.d0 !in-plane intergrain exchange sigma
       rand_mat=89032211   !some large integer to scale random seeds for material properities of grains
       vertical_a=1.13d-6  !typical intragranular exchange stiffness (erg/cm)?
       vertical_same_hk=.true.
   rewind(ionum)
   call output(' reading `media` namelist')
!  do i=1,num_layers
   i=0;do while (i.lt.num_layers); i=i+1
     rand_seed=57; final_packing_fraction=-1.d0
     if (i.eq.1.or.vertical_spins.lt.1) then; vertical_spins=1
       ex_upper=0.d0; ex_upper_sigma=0.d0
       same_grain_p=.false.;same_hk_p=.false.;same_exchange_p=.false.
       same_tc_p=.false.;same_ms_p=.false.;same_init_mag_p=.false.
       read(ionum,media)
       if (vertical_spins.gt.1) then; num_end=1; vol_scale=dble(vertical_spins); if (vertical_same_hk) vol_scale=-vol_scale; ex_upper = 2.d0*vertical_a/dzl(i); endif
     else
       num_end=0; if (vertical_spins.eq.1) num_end=1; 
       same_grain_p =.true.; same_hk_p = .true.; same_exchange_p=.true.; same_init_mag_p =.true.; 
       same_tc_p = .true.; same_ms_p = .true.;
     endif
     write(os,("(a,x,i3)")) '  read media layer #',i; call output(os)
     packing_fraction=abs(packing_fraction);if (packing_fraction.gt.1.d0) packing_fraction=packing_fraction/1.d2
     if (final_packing_fraction.le.0.d0) then
       final_packing_fraction=packing_fraction
     else
       final_packing_fraction=abs(final_packing_fraction);if (final_packing_fraction.gt.1.d0) final_packing_fraction=final_packing_fraction/1.d2
     endif
     small_packing_fraction=abs(small_packing_fraction);if (small_packing_fraction.gt.1.d0) small_packing_fraction=small_packing_fraction/1.d2
     if (i.eq.1) then
       ok_p=set_construct_grain(i,random_seed=rand_seed,rand_mat=rand_mat); ok_p=start_mag(i,rand_seed*2)
       ok_p=set_mag(layer=i,ac_p=ac_p,saturate_p=saturate_p,m=(/mx,my,mz/)); prev_rand_seed=rand_seed*2
     else
       ok_p=set_construct_grain(i,random_seed=rand_seed,rand_mat=rand_mat,same_grain_p=same_grain_p,same_hk_p=same_hk_p,same_exchange_p=same_exchange_p, & 
                same_tc_p=same_tc_p,same_ms_p=same_ms_p)
       if (same_init_mag_p) then; j=prev_rand_seed; else; j=rand_seed*(i+1); endif; ok_p=start_mag(i,j); prev_rand_seed=j
       ok_p=set_mag(layer=i,ac_p=ac_p,saturate_p=saturate_p,m=(/mx,my,mz/))
!      ok_p=set_init_mag(i,random_seed=rand_seed,ac_p=ac_p,same_init_mag_p=same_init_mag_p,saturate_p=saturate_p,m=(/mx,my,mz/))
     endif

     if (initial_grain_area.ge.1.d0) initial_grain_area=initial_grain_area*dx*dy
     if (average_small_grain_area.ge.1.d0) average_small_grain_area=average_small_grain_area*dx*dy
     if (initial_small_grain_area.ge.1.d0) initial_small_grain_area=initial_small_grain_area*dx*dy
     if (average_grain_area.ge.1.d0) average_grain_area=average_grain_area*dx*dy
     if (small_grain_sigma.gt.20.d0*min(dx,dy)) small_grain_sigma=small_grain_sigma*min(dx,dy)
     ok_p=set_construct_grain(i,grain_p=(average_grain_area.ne.1.d0.and.average_grain_area.ne.dx*dy),packing_fraction=packing_fraction, &
       small_packing_fraction=small_packing_fraction,initial_grain_area=initial_grain_area,small_grain_area=average_small_grain_area, &
       initial_small_grain_area=initial_small_grain_area,grain_sigma=init_vol_sigma,small_grain_sigma=small_grain_sigma,average_grain_area=average_grain_area, &
       dr_grain=abs(dr_grain),old_small_grain=old_small_grain,grow_small_grain=grow_small_grain,final_packing_fraction=final_packing_fraction)
     x=kx*kx+ky*ky+kz*kz; if (x.ne.0.d0) x=1.d0/sqrt(x)
     kx=x*kx; ky=x*ky; kz=x*kz
     x=hk_ac; if (scale_hk.and.x.eq.0.d0)x=0.37d-7; if (scale_hk.and.hk_ac.lt.0.d0) x=0.37d-7*abs(x)
     ok_p=set_construct_grain(i,kx=kx,ky=ky,kz=kz,k1=hk,k2=hk2,k_ang=k_ang,k_cone_sigma=hk_cone_sigma,kmax=max_k,k_frac_inplane=frac_in_plane, &
           k_frac_random=frac_random,k_sigma=k_sigma,k_3d=hk_3d_random_p,k_2d=hk_2d_random_p,k_uniform=hk_uniform_p,k_cone=hk_cone_p, &
           k_cubic=hk_cubic_p,k_ac=x,k_scale=scale_hk,k_rescale=rescale_hk,k_scale_al=scale_hk_al,hkmin=min_hk,vol_scale=vol_scale, &
           k_sigma_small=k_sigma_small,k1_small=hk_small,k2_small=hk2_small,frac_hk_small=frac_hk_small,num_end=num_end,print_dist_p=print_dist_p)

     x=tc; if (x.lt.1.d0) x=1.d22
     y=tc_max; if (y.le.1.d0) y=1.d22
!    if (tc_eta.ne.0.d0) tc_eta=-1.d0/tc_eta
     ok_p=set_construct_grain(i,tc_nom=x,tc_sigma=tc_sigma,tc_max=y,tc_d0=tc_d0,tc_eta=tc_eta,tc_rescale=tc_rescale_p,tau_nom=tau,tau_sigma=tau_sigma)
     ok_p=set_temp_scaling(i,tc_nom=x)

     ok_p=set_construct_grain(i,ms_nom=abs(ms),ms_sigma=ms_sigma,ms_max=abs(ms_max),ex=ex,ex_sigma=ex_sigma,ex_upper=ex_upper,ex_upper_sigma=ex_upper_sigma)

     call set_data(i,spin,spin_p=.true.)
     call set_data(i,heatsink_dt,tau_hs_t_p=.true.)
     call set_data(i,trim(loop_file))
  vertical_spins=vertical_spins-1
   enddo

  !LLG parameters
   call output(' reading `llg` namelist')
   alpha=1.d0  !Gilbert damping parameter
   gamma=1.76d7 !gyromagnetic ration (rad/Oe/s)
   dt=.1d-13    !(initial) LLG time step (s) 
   thermal_p=.false. !stochastic thermal fields?
   thermal_seed=34;  !random seed for stochastic fields?
   max_dm=180.d0  !maximum allowed magnetization rotation in a time step (degree)
   min_dm=0.d0    !minimum allowed magnetization rotation in a time step (degree)
   zhang=.false.
   mt=.false.
   mt_bs=.false.
   zhang_l=.false.
   garanin_tau=.true.
   llb=.false.
   fp=.false.
   dt_frac=0.2d0  !fraction of gyroperiod to estimate LLG timestep
   dt_max=1.d22   !maximum allowed LLG time step (s)
   t_max=1.d22    !simulation time to integrate LLG to (s)
   t_min=-1.d22   !must integrate LLG to at least this simulation time (s)
   del_t=1.d22    !maximum time to integrate LLG (s)
   conv_torque=5.d-3 !convergence criteria for |mxH|/Ms
   conv_angle=0.d0 !convergence criteria for magnetization rotation (degree)
   temperature=293.15d0 !temperature (K) 
   fixed_dt_p=.false. !LLG time step is fixed?
   iter=100        !print to screen between this many iterations in LLG
   min_iter=10
   min_ave_m=-1.d22*(/ 1.d0, 1.d0, 1.d0 /)  !convergence criteria on minimum average normalized magnetization
   max_ave_m= 1.d22*(/ 1.d0, 1.d0, 1.d0 /)   !convergence criteria on minimum average normalized magnetization
   plot_iter=-1    !plot magnetization/write field between this many iterations in LLG (<1 means plot_iter=iter)
   thermal_num_ave = 5001  !maximum number of time steps to average max torque for convergence w/ stochastic fields
   new_integrator = .true. !use LLG integrator, or just LL?
   llg_type = -9
   msexp=-1.d0; hkexp=-1.d0; use_mag_p=.false.
   media_layer=0
   do i=1,num_layers
     call set_data(i,alpha,alpha_p=.true.)
     call set_data(i,gamma,gamma_p=.true.)
     call set_data(i,temperature,background_t_p=.true.)
     call set_data(i,thermal_p,thermal_p=.true.)
   enddo

   ok_p=start_llg(num_layers)
   rewind(ionum); i=0
   do while (i.eq.0)    !read multiple llg?
     read(ionum,llg,iostat=i)
     if (i.gt.0) then
        write(os,"(' ERROR reading namelist `llg` in ',a)")trim(in_file)
        call output(trim(os))
        close(ionum); ok_p=.false.
        return
     elseif (i.eq.0) then
       media_layer=max(0,min(num_layers,media_layer))
       write(os,"('  setting default llg values')")
       if (media_layer.ne.0) write(os,"('  setting llg values for layer #',i3)")media_layer
       call output(trim(os))
       if (media_layer.eq.0) then
         do i=1,num_layers
           call set_data(i,alpha,alpha_p=.true.)
           call set_data(i,gamma,gamma_p=.true.)
           call set_data(i,temperature,background_t_p=.true.)
           call set_data(i,thermal_p,thermal_p=.true.)
!          ok_p= set_temp_scaling(i,hkexp=hkexp,use_mag_p=use_mag_p)
           ok_p= set_temp_scaling(i,msexp=msexp,hkexp=hkexp,use_mag_p=use_mag_p)
         enddo
       else
         call set_data(media_layer,alpha,alpha_p=.true.)
         call set_data(media_layer,gamma,gamma_p=.true.)
         call set_data(media_layer,temperature,background_t_p=.true.)
         call set_data(media_layer,thermal_p,thermal_p=.true.)
!        ok_p= set_temp_scaling(media_layer,hkexp=hkexp,use_mag_p=use_mag_p)
         ok_p= set_temp_scaling(media_layer,msexp=msexp,hkexp=hkexp,use_mag_p=use_mag_p)
       endif
       i=iter; if (plot_iter.gt.0) i=plot_iter
       if (llg_type.eq.-9) then; llg_type=1; if (new_integrator) llg_type=2; endif
       if (zhang) then; mt_bs=.not.zhang_l; elseif (llb) then; mt_bs=.false.; endif
       call set_llg(media_layer,dt=dt,max_dm=max_dm,min_dm=min_dm,dt_frac=dt_frac,dt_max=dt_max, &
                    t_max=t_max,t_min=t_min,del_t=del_t,conv_torque=conv_torque,conv_angle=conv_angle, &
                    fixed_dt_p=fixed_dt_p,iter=iter,min_ave_m=min_ave_m,max_ave_m=max_ave_m,plot_iter=i,&
!                   thermal_num_ave=thermal_num_ave,llg_type=llg_type,zhang=zhang,llb=llb,fp=fp)
                    garanin_tau=garanin_tau,min_it=max(1,min_iter), &
                    thermal_num_ave=thermal_num_ave,llg_type=llg_type,meq_bs=mt_bs,mt=mt,zhang=zhang,llb=llb,fp=fp)
       i=0
     endif
   enddo


  !plot parameters
   rewind(ionum)
   call output(' reading `plot` namelist')
!  call init_plot(m1,num_layers)
   i=0; num_layers=-num_layers
   do while (i.eq.0)    !read multiple llg?
    if (num_layers.lt.0) then
     num_layers=abs(num_layers)
     plot_p=.false.  !turn on graphics?
     write_plot=.false. !dump ascii files of contour plots
     plot_hk_p=.true. !plot Hk histograms?
     plot_temp_p=.false.  !plot temperature contour averaged on grain?
     plot_hamr_scales_p=.false. !plot Ms,Hk,gamma,alpha vs T?
     plot_grain_p=.true. !plot grain area histogram?
     plot_fld_p=.true.   !plot (perpendicular) write field contour averaged on grain?
     plot_tc_p=.true.
     plot_tau_p=.true.
     plot_hx_p=.true.  !plot intralayer exchange histogram, scatter and contour?
     plot_hx_il_p=.true.  !plot interlayer exchange histogram, scatter and contour?
     plot_mag_p=.true. !plot magnetization projected along nominal Hk direction?
     plot_ms_p=.true.
     plot_mz_p=.true. !plot mz component
     plot_my_p=.false. !plot my component
     plot_mx_p=.false. !plot mx component
     plot_wfx_p=.false. !plot write Hx averaged on grain
     plot_wfy_p=.false. !plot write Hy averaged on grain
     plot_wfxf_p=.false. !plot write field Hx as field
     plot_wfyf_p=.false. !plot write field Hy as field
     plot_wfzf_p=.false. !plot write field Hz as field
     plot_tempf_p=.false. !plot temperature as field
     plot_m_p=.true.    !plot |m|
     plot_mh_p=.false. !.true.  !plot/print m(h)?
     plot_heff_p=.false.  !plot SW effective field contour?
     plot_signed_heff_p= .false. !plot contours of SW field with sign of Hz ?
     plot_hscale_p=.true.   !plot write field scaling vs. time
     num_level_mag = -6    !number of contour levels for m.k plot
     plot_scatter_p=.false.
     plot_scatter_len_p=.false.
    endif
     media_layer=0
      read(ionum,plot,iostat=i)
      if (i.gt.0) then
        write(os,"(' ERROR reading namelist `plot` in ',a)")trim(in_file)
        call output(trim(os))
        close(10); ok_p=.false.
        return
      elseif (i.eq.0) then
       media_layer=max(0,min(num_layers,media_layer))
       write(os,"('  setting default plotting options')")
       if (media_layer.ne.0) write(os,"('  setting plotting options for layer #',i3)")media_layer
       call output(trim(os))
       call set_plot(il=media_layer,hk=plot_hk_p,hk_hist=plot_hk_p,hk_vol=plot_hk_p,temp=plot_temp_p,hamr=plot_hamr_scales_p, &
          grain_hist=plot_grain_p,wf_z=plot_fld_p,tc=plot_tc_p,tc_hist=plot_tc_p,tc_vol=plot_tc_p,hx=plot_hx_p,hx_hist=plot_hx_p,hx_vol=plot_hx_p, &
          hv=plot_hx_il_p,hv_hist=plot_hx_il_p,hv_vol=plot_hx_il_p,mz=plot_mz_p,wf_x=plot_wfx_p,wf_y=plot_wfy_p,mag=plot_m_p,mk=plot_mag_p, &
          wf_xf=plot_wfxf_p,wf_yf=plot_wfyf_p,wf_zf=plot_wfzf_p,tempf=plot_tempf_p,mx=plot_mx_p,my=plot_my_p, &
          mhloop=plot_mh_p,wf_eff=(plot_heff_p.or.plot_signed_heff_p),ms=plot_ms_p, &
          scatter=plot_scatter_p,scatter_len=plot_scatter_len_p,tau=plot_tau_p,tau_hist=plot_tau_p,tau_vol=plot_tau_p)
       if (media_layer.eq.0) call set_plot(plotting=plot_p,write_plot=write_plot)
       if (media_layer.ne.0) call set_plot(il=media_layer,plot=plot_p)
      endif
   enddo

   call create_checker(dx,dy)

  !pre-written bit patterns
   rewind(ionum)
   call output(' reading `pattern` namelist(s)')
   i=0
   do while (i.eq.0)
     media_layer=0
     mx=0.d0   !nominal x component of magnetization
     my=0.d0   !nominal y component of magnetization
     mz=1.d0   !nominal y component of magnetization
     ac_p=.false. !AC erased? 
     saturate_p= .false. ! init along m instead of Hk?
     xmin=-1.d22  !minimum x position of this region (cm)
     xmax=1.d22   !maximum x position of this region (cm)
     ymin=-1.d22  !minimum y position of this region (cm)
     ymax=1.d22   !maximum y position of this region (cm)
     read(10,pattern,iostat=i)
     if (i.gt.0) then
       write(os,"(' ERROR reading namelist `pattern` in ',a)")trim(in_file)
       call output(trim(os))
       close(ionum);ok_p=.false.
       return
     elseif (i.eq.0) then
       ok_p=set_mag( m=(/ mx, my, mz /), ac_p=ac_p, saturate_p= saturate_p, &
          xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax,layer=media_layer)
     endif
   enddo

  !decay parameters
   call output(' reading `decay` namelist')
   ntraj=-1; traj=-1; torbit=-1.d0; max_eb=-1.d0; iter=0; boltz_p=.true.; vol_frac=-1.d0
   test_rot_p=.true.; diff_m=.true.
   rewind(ionum)
!  open(ionum,file=trim(in_file),status='old')
    read(ionum,decay,iostat=i)
    if (i.gt.0) then
      write(os,"(' ERROR reading namelist `decay` in ',a)")trim(in_file)
      call output(trim(os))
      close(ionum); ok_p=.false.
      return
    endif
   if (ntraj.ne.-1) ok_p=set_decay_mag(intraj=ntraj)
   if (traj.ne.-1) ok_p=set_decay_mag(itraj=traj)
   if (torbit.gt.0.d0) ok_p=set_decay_mag(itorbit=torbit)
   if (max_eb.gt.0.d0) ok_p=set_decay_mag(imax_eb=max_eb)
   if (iter.ne.0) ok_p=set_decay_mag(iiter=iter)
   if (vol_frac.ge.0) ok_p=set_decay_mag(ivol_frac=vol_frac)
   ok_p=set_decay_mag(idiff_m=diff_m)
   ok_p=set_decay_mag(itest_rot_p=test_rot_p)
   ok_p=set_decay_mag(iboltz=boltz_p)

   if (.not.checkpoint(present_p=.true.)) ok_p=create_loop() !create loop cli only if not restarting from checkpoint
   close(ionum)

   inquire(file=trim(adjustl(clifile)),exist=ok_p)
   ok_p = (ok_p .or.trim(adjustl(clifile)).eq.'-')

   deallocate(dzl,il_space, thm_seed)

   return
  contains
  
   subroutine create_checker(dxs,dys)
  !this subroutine creates a checker-board pattern for the initial magnetization
  ! if the namelist 'checker' is in the input file

    real(8),intent(in)::dxs,dys
    integer::j
    real(8)::x0,y0,xcenter,ycenter
    namelist /checker/mx,my,mz,dx,dy,xcenter,ycenter,media_layer

   !checker pattern
    call output(' reading `checker` namelist')
    media_layer=0
    kx=mx; ky=my; kz=mz
    mx=0.d0  !nominal x component of magnetization
    my=0.d0  !nominal y component of magnetization
    mz=0.d0  !nominal z component of magnetization
    dx=.5d0*dx*nx  !size along x of one checker (1/2 the simulation area) (cm)
    dy=.5d0*dy*ny  !size along y of one checker (1/2 the simulation area) (cm)
    xcenter=1.d22  !center position of one checker along x (cm)
    ycenter=1.d22  !center position of one checker along y (cm)
    rewind(ionum)
     read(ionum,checker,iostat=i)
     if (i.gt.0) then
       write(os,"(' ERROR reading namelist `checker` in ',a)")trim(in_file)
       call output(trim(os))
       close(10)
       return
     elseif (i.eq.0) then
       media_layer=max(0,media_layer)
       write(os,"(' creating checker pattern on all layers')")
       if (media_layer.ne.0) write(os,"('  creating checker pattern for layer #',i3)")media_layer
       call output(trim(os))
       dx=abs(dx); if (dx.eq.0.d0) dx=.5d0*dxs*nx
       dy=abs(dy); if (dy.eq.0.d0) dy=.5d0*dys*ny
       if (xcenter.gt.1.d20) xcenter=dx*.5d0
       if (ycenter.gt.1.d20) ycenter=dy*.5d0
       do while (xcenter-.5d0*dx.gt.0.d0);xcenter=xcenter-dx;enddo
       do while (xcenter+.5d0*dx.le.0.d0);xcenter=xcenter+dx;enddo;xcenter=xcenter-.5d0*dx
       do while (ycenter-.5d0*dy.gt.0.d0);ycenter=ycenter-dy;enddo
       do while (ycenter+.5d0*dy.le.0.d0);ycenter=ycenter+dy;enddo;ycenter=ycenter-.5d0*dy
       x=mx*mx+my*my+mz*mz; if (x.eq.0.d0) then; mx=kx; my=ky; mz=kz; endif
       x=mx*mx+my*my+mz*mz; if (x.eq.0.d0) mz=1.d0
       y0=ycenter; i=-1
       do while (y0.lt.dys*ny)
        x0=xcenter; i=-i; j=-i
        do while (x0.lt.dxs*nx)
         ok_p=set_mag( m=(/ mx, my, mz /)*sign(1,j), xmin=x0, xmax=x0+dx, ymin=y0, ymax=y0+dy,layer=media_layer)
         x0=x0+dx;j=-j
        enddo
        y0=y0+dy
       enddo
     endif
   end subroutine create_checker

   function create_loop() result (ok_p)
  !this routine creates a cli file for loops
    use io_m, only: add_file, delete_file
    logical::loop_p,ok_p,continous_p
    real(8):: hx,hy,hz,h_start,dh,hmax,hmin,hmin_fine,hmax_fine,dh_fine,dt_max,dt_max_fine

    logical::end_p,write_grain_p,write_m_p,rotate_p,decay_p,decay_old_gamma_p,decay_multispin_p,replace_m_p
    integer::i,io,decay_seed
    real(8)::x,dhc,time,h_next,last_dt,h_rot,rotate_rate,rotate_axis,rotate_ang1,rotate_ang2,decay_temp,decay_sweep_rate,decay_wait,decay_mkt
    character(80)::grain_file,m_file


    namelist /loop/ hx,hy,hz,h_start,dh,hmax,hmin,hmin_fine,hmax_fine,dh_fine,loop_p,continous_p, &
        dt_max, dt_max_fine,write_grain_p,write_m_p,grain_file,m_file,rotate_p,rotate_rate,rotate_axis,rotate_ang1,rotate_ang2, &
        decay_p,decay_temp,decay_sweep_rate,decay_wait,decay_mkt,decay_old_gamma_p,decay_multispin_p,decay_seed,replace_m_p

    decay_p=.false.; decay_temp=300.d0; decay_sweep_rate=190.d0; decay_wait=0.d0; decay_mkt=0.d0
    ok_p=.false.; decay_old_gamma_p=.false.;decay_multispin_p=.false.;decay_seed=4357

    write_grain_p=.true.  !write grain file?
    write_m_p=.false.     !write magnetization file?
    replace_m_p=.true.   !overwrite magnetization file?
    grain_file='grainloop.dat'  !name of grain file
    m_file='magloop.dat'   !name of magnetization file
    hx=0.d0; hy=0.d0; hz=0.d0   !uniform applied field direction (normalized field)
    hmin=0.d0    !minimum applied field magnitude (Oe)
    hmax=0.d0    !maximum applied field magnitude (Oe)
    dh=0.d0      !field step (Oe)- sign determins the initial field sweep direction
    hmin_fine=1.d22 !minimum applied field amplitude using dh_file field steps (Oe)
    hmax_fine=1.d22 !minimum applied field amplitude using dh_file field steps (Oe)
    dh_fine=0.d0  !field step (Oe)
    h_start=1.d22  !starting field for 'loop'
    loop_p=.false. !want a full loop?
    continous_p=.false. !change the field continously during LLG integration?
    dt_max = del_t  !maximum time in LLG between field steps when steppin w/ dh (s) {default via som.inp's `llg` namelist}
    dt_max_fine = del_t  !maximum time in LLG between field steps when stepping w/ dh_fine (s) {default via som.inp's `llg` namelist}
!   dt_max=get_llg_delt(); dt_max_fine=dt_max
    rotate_p = .false.    !constant magnitude field, just rotate it
    rotate_rate = 0.d0    !degrees per second or degrees per iteration
    rotate_axis = 2       !rotate field about which axis (1-x, 2-y and 3-z)
    rotate_ang1 = 0.d0    !starting angle
    rotate_ang2 = 360.d0  !ending angle

    rewind(ionum)
!   open(10,file=trim(in_file),status='old')
     read(ionum,loop,iostat=i)
     if (i.gt.0) then
       write(os,"(' ERROR reading namelist `loop` in ',a)")trim(in_file)
       call output(trim(os))
!      close(10)
       return
     elseif (i.ne.0) then
!      close(10)
       return
     endif
!   close(10)

    if (decay_p) continous_p=.false.

    rotate_rate=rotate_rate*acos(-1.d0)/180.d0  !rad/s or rad per iteration
    rotate_ang1=rotate_ang1*acos(-1.d0)/180.d0  !starting angle (rad)
    rotate_ang2=rotate_ang2*acos(-1.d0)/180.d0  !ending angle (rad)
    if (abs(rotate_axis).lt.1.or.abs(rotate_axis).gt.3) rotate_axis=2

    grain_file=adjustl(grain_file); m_file=adjustl(m_file)
    if (grain_file(1:1).eq.' ') grain_file='grainloop.dat'
    if (m_file(1:1).eq.' ') m_file='mloop.dat'
    if (hmin.gt.hmax) then
      x=hmin; hmin=hmax; hmax=x
    endif
    if (hmin_fine.gt.hmax_fine) then
      x=hmin_fine; hmin_fine=hmax_fine; hmax_fine=x
    endif
    if (min(hmin_fine,hmax_fine).lt.1.d20) then   !fine stepping
      if (hmin_fine.gt.1d20) hmin_fine=hmax_fine
      if (hmax_fine.gt.1d20) hmax_fine=hmin_fine
      hmax = max(hmax,hmax_fine); hmin = min(hmin,hmin_fine)
    endif
    if (h_start.gt.1.d20) then                    !get starting field
      h_start=hmin; if (dh.lt.0.d0) h_start=hmax
      end_p=.false.
    else
      h_start=min(hmax,max(h_start,hmin))
      end_p=.true.
    endif
    if (dh.eq.0.d0) dh=max(1.d0,hmax-hmin)*1.1d0
    if (dh_fine.eq.0.d0) dh_fine=dh
    dh_fine=sign(dh_fine,dh); dhc=dh
    if (hmin_fine.ne.hmax_fine.and.(hmin_fine.lt.hmax.or.hmax_fine.gt.hmin)) then    !we are doing fine seps
      if ((dhc.gt.0.d0.and.h_start+dhc.gt.hmin_fine.and.h_start.lt.hmax_fine) .or. &
          (dhc.lt.0.d0.and.h_start+dhc.lt.hmax_fine.and.h_start.gt.hmin_fine)) dhc=dh_fine
    else
      hmin_fine=1.d22; hmax_fine=1.d22; dh_fine=1.d22
    endif
    h_next=h_start+dhc; last_dt=del_t

    cli_file=adjustl(cli_file)
    if (cli_file(1:1).eq.' ') cli_file='som.cli'

    write(os,"(' read namelist `loop` in ',a,', constructing cli file ',a)")trim(in_file),trim(cli_file)
    call output(trim(os))

    io=add_file(trim(cli_file))
    open(io,file=trim(cli_file),status='unknown')
    write(io,"('!')")
    write(io,"('!automatically generated file for loop simulation')")
    write(io,"('!')")
    write(io,"('!read in temperature scalings')")
    write(io,"('openr      MsHktemp.dat')")
    write(io,"('read_scale MsHktemp.dat Ms')")
    write(io,"('read_scale MsHktemp.dat Hk')")
    write(io,"('close      MsHktemp.dat')")

    write(io,"('!')")
    write(io,"('start_clock')")
    write(io,"('!')")
    write(io,"('plot_fld .false.')")
    write(io,"('plot_heff .false.')")
    write(io,"('plot_signed_heff .false.')")
    if (rotate_p) then
      h_rot=sqrt(hx*hx+hy*hy+hz*hz)
      if (h_rot.eq.0.d0) return
      hz=h_rot*cos(acos(-1.d0)*rotate_ang1/180.d0)
      hx=h_rot*sin(acos(-1.d0)*rotate_ang1/180.d0)
      if (abs(rotate_axis).eq.1) then     !rotate about x axis, zero is along z (normal to plane) and increasing to y
        write(io,"('plot_mh .true.',1p,3(1x,e19.12))") 0.d0,hx,hz
      elseif (abs(rotate_axis).eq.2) then !rotate about y axis, zero is along z (normal to plane) and increasing to x
        write(io,"('plot_mh .true.',1p,3(1x,e19.12))") hx,0.d0,hz
      else
        write(io,"('plot_mh .true.',1p,3(1x,e19.12))") hz,hx,0.d0
      endif
    else
      x=hx*hx+hy*hy+hz*hz; if (x.eq.0.d0) hz=1.d0
      x=1.d0/sqrt(hx*hx+hy*hy+hz*hz); hx=hx*x; hy=hy*x; hz=hz*x
      write(io,"('plot_mh .true.',1p,3(1x,e19.12))") hx,hy,hz
    endif
   !set applied field
    write(io,"('!')")
    if (rotate_p) then
     if (continous_p) then
      write(io,"('!create two uniform fields, label them `1` and `2`')")
      i=1
      if (abs(rotate_axis).lt.3) then; write(io,"('create_field ',i1,1p,3(1x,e19.12))") i,0.d0,0.d0,1.d0; i=i+1; endif
      if (abs(rotate_axis).eq.2) then; write(io,"('create_field ',i1,1p,3(1x,e19.12))") i,1.d0,0.d0,0.d0; i=i+1; endif
      if (abs(rotate_axis).eq.1) write(io,"('create_field ',i1,1p,3(1x,e19.12))") i,0.d0,1.d0,0.d0
     endif
    else
      write(io,"('!create uniform field, label it `1`')")
      write(io,"('create_field 1',1p,3(1x,e19.12))") hx,hy,hz
    endif
   !generate grains
    write(io,"('!')")
    write(io,"('!create grain microstructure')")
    write(io,"('gen_grain')")
   !save grain structure
    if (write_grain_p) then
      write(io,"('!')")
      write(io,"('!save grain microstructure')")
!     write(io,"('openw ',a,' !append !unformatted')") trim(grain_file)
      write(io,"('openw ',a,' unformatted')") trim(grain_file)
      write(io,"('write_grain ',a)") trim(grain_file)
      write(io,"('close ',a,' !delete')") trim(grain_file)
    endif
    if (write_m_p.and..not.replace_m_p) then
      write(io,"('!')")
!     write(io,"('openw ',a,' !append !unformatted')") trim(m_file)
      write(io,"('openw ',a,' unformatted')") trim(m_file)
    endif
   !initialize magnetization
    write(io,"('!')")
    write(io,"('!initialize magnetization')")
    write(io,"('init_mag')")
    write(io,"('!')")
    write(io,"('plot all')")
    write(io,"('!')")
    write(io,"('!set simulation time')")
    write(io,"('set_time 0.d0')")
    write(io,"('!')")
    if (decay_p) then
       write(io,"('set_1field 1',1x,1p,e19.12)") h_start
       write(io,"('llg')")
    endif
    write(io,"('!now do loop over field')")
   !do field loop
    ok_p=.true.; time=0.d0
    do while (ok_p)
     write(io,"('!')")
     x=dt_max; if (abs(dhc).eq.abs(dh_fine)) x=dt_max_fine
     if (continous_p) then
       if (rotate_p) then
         write(io,"('set_1field 1 ',a,1p,5(1x,e19.12))") 'COS',rotate_rate*sign(1.d0,rotate_axis),rotate_ang1,h_rot,0.d0,1.d0
         write(io,"('set_1field 2 ',a,1p,5(1x,e19.12))") 'SIN',rotate_rate*sign(1.d0,rotate_axis),rotate_ang1,h_rot,0.d0,1.d0
       else
         write(io,"('set_1field 1',1p,4(1x,e19.12))") time,time+x,h_start,h_next
       endif
     else
       if (rotate_p) then
         write(io,"('!create uniform field, label it `1`')")
         hx=sin(rotate_ang1); hz=cos(rotate_ang1); hy=0.d0
         if (abs(rotate_axis).eq.1) then
           write(io,"('create_field 1',1p,3(1x,e19.12))") 0.d0,hx,hz
           hy=hx;hx=0.d0
         elseif (abs(rotate_axis).eq.2) then
           write(io,"('create_field 1',1p,3(1x,e19.12))") hx,0.d0,hz
         else
           write(io,"('create_field 1',1p,3(1x,e19.12))") hz,hx,0.d0
           hy=hx;hx=hz;hz=0.d0
         endif
         write(io,"('set_1field 1',1x,1p,e19.12)") h_rot
       else
         write(io,"('set_1field 1',1x,1p,e19.12)") h_start
       endif
     endif
     if (x.ne.last_dt) write(io,"('LLG_DELT',1x,1p,e19.12)") x; last_dt=x
     time=time+x
     write(io,"('start_clock')")
     if (decay_p) then
       write(io,"('decay ',1p,5(e19.12,x),0p,2(1l,x),i8)") abs(h_start-h_next)/decay_sweep_rate,decay_temp,decay_sweep_rate, &
                          h_start,decay_mkt,decay_old_gamma_p,decay_multispin_p,decay_seed
       decay_seed=decay_seed+10000
!      if (decay_old_gamma_p) then
!        write(io,"('decay ',1p,5(e19.12,x),0p,a,1l)") abs(h_start-h_next)/decay_sweep_rate,decay_temp,decay_sweep_rate,h_start,decay_mkt,'T',decay_multispin_p
!      else
!        write(io,"('decay ',1p,5(e19.12,x),0p,a,1l)") abs(h_start-h_next)/decay_sweep_rate,decay_temp,decay_sweep_rate,h_start,decay_mkt,'F',decay_multispin_p
!      endif
     else
       write(io,"('llg')")
     endif
     write(io,"('stop_clock')")
     if (rotate_p) then
       if (continous_p) then
         hx=h_rot*sin(rotate_rate*sign(1.d0,rotate_axis)*time); hz=h_rot*cos(rotate_rate*sign(1.d0,rotate_axis)*time)
         if (abs(rotate_axis).eq.1) then
           write(io,"('plot mh',4(x,e19.12))") 0.d0,hx,hz
         elseif (abs(rotate_axis).eq.2) then
           write(io,"('plot mh',4(x,e19.12))") hx,0.d0,hz
         else
           write(io,"('plot mh',4(x,e19.12))") hz,hx,0.d0
         endif
       else
         write(io,"('plot mh',4(x,e19.12))") h_rot*hx,h_rot*hy,h_rot*hz
       endif
     else
       if (continous_p) then
         write(io,"('plot mh ',e19.12)") h_next/1.d4    !print field in Telsa
       else
         write(io,"('plot mh ',e19.12)") h_start/1.d4   !print field in Telsa
       endif
     endif
     write(io,"('plot mag')")
     if (write_m_p) then
       if (replace_m_p) then
             write(io,"('openw ',a,' unformatted')") trim(m_file)
             write(io,"('write_mag ',a)") trim(m_file)
             write(io,"('close ',a)") trim(m_file)
       else
             write(io,"('write_mag ',a)") trim(m_file)
       endif
     endif
     if (rotate_p) then
       if (continous_p) then
         hx=rotate_rate*time
         ok_p = ( abs(hx) .lt. abs(rotate_ang1-rotate_ang2) )
       else
         hx=rotate_ang1+rotate_rate*sign(1.d0,rotate_axis)
         ok_p = ( ((rotate_ang1.gt.rotate_ang2.or.hx.le.rotate_ang2).and.rotate_axis.gt.0) .or. &
                  ((rotate_ang1.lt.rotate_ang2.or.hx.ge.rotate_ang2).and.rotate_axis.lt.0) )
         rotate_ang1=hx
       endif
     else
       h_start=h_start+dhc
       if ((h_start.lt.hmin-.1d0*abs(dhc).or.h_start.gt.hmax+.1d0*abs(dhc)).and.(loop_p.or.end_p)) then
         if (end_p) then; 
            end_p=.false.; 
         elseif (loop_p) then 
            loop_p=.false.; 
            if (decay_p.and.decay_wait.gt.0.d0) then
              write(io,"('start_clock')")
!             write(10,"('decay ',1p,5(e19.12,x))") decay_wait,decay_temp,decay_sweep_rate,h_start,decay_mkt
              write(io,"('decay ',1p,5(e19.12,x),0p,2(1l,x),i8)") decay_wait,decay_temp,decay_sweep_rate, &
                          h_start,decay_mkt,decay_old_gamma_p,decay_multispin_p,decay_seed
              decay_seed=decay_seed+10000
              write(io,"('stop_clock')")
            endif
         endif
         dhc=-dhc; dh=-dh; dh_fine=-dh_fine; h_start=h_start+2.d0*dhc
         hmin_fine=-hmin_fine; hmax_fine=-hmax_fine; x=min(hmin_fine, hmax_fine)
         hmax_fine=max(hmin_fine,hmax_fine); hmin_fine=x
       endif
       if ((dhc.gt.0.d0.and.h_start+dhc.gt.hmin_fine+.1d0*abs(dh_fine).and. &
             h_start.lt.hmax_fine-.1d0*abs(dh_fine)).or. &
           (dhc.lt.0.d0.and.h_start+dhc.lt.hmax_fine-.1d0*abs(dh_fine).and. &
             h_start.gt.hmin_fine+1.d0*abs(dh_fine))) then
         dhc=dh_fine
       else
         dhc=dh
       endif
       h_next=h_start+dhc
       if ((h_next.lt.hmin-.1d0*abs(dhc).or.h_start.gt.hmax+.1d0*abs(dhc)).and.(loop_p.or.end_p)) &
         h_next=h_next-2.d0*dhc
!        h_next=-2.d0*dhc
       ok_p=(loop_p.or.end_p.or.(h_start.gt.hmin-0.1d0*abs(dhc).and. &
                                 h_start.lt.hmax+0.1d0*abs(dhc)))
     endif
    enddo
    write(io,"('!')")
    write(io,"('!time for simulation:')")
    write(io,"('stop_clock')")
    close(io)
    call delete_file(io)

    ok_p=.true.
    return
   end function create_loop

 end function get_input

end module input_m

