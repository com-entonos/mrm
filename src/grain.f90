

module grain_m
 implicit none
 private
 public :: generate_area, generate_hk, generate_media, file_media, show_media
contains

  function generate_hk() result (ok_p)
   logical::ok_p
   ok_p=generate_media(iarea_p=.false.)
  end function
  function generate_area() result (ok_p)
   logical::ok_p
   ok_p=generate_media(ihk_p=.false.)
  end function

  function file_media(read_p,ionum,layer,old_p) result (ok_p)
   use data_m, only: file_grain, get_size, get_data, set_data
   use field_m, only: file_field
   use temp_scaling_m, only: init_temp_scaling
   use construct_grain_m, only: construct_ms
   use llg_m, only: init_llg
   use io_m, only: output
   logical,intent(in)::read_p
   integer,intent(in)::ionum
   integer,intent(in),optional::layer
   logical,intent(in),optional::old_p

   logical::ok_p,old
   integer::k,k1,k2,nx,ny,n_lay,i1,n
   real(8)::dx,dy
   real(8),pointer::grain_area(:),ms(:),tc(:),mag(:,:),lambdaMs(:),spin

   old=.false.; if (present(old_p)) old=old_p
   ok_p=.false.
   call get_size(nx,ny,n_lay,dx,dy)
   k1=1; k2=n_lay; if (present(layer)) then; if (layer.gt.0.and.layer.le.n_lay) then; k1=layer; k2=layer; endif; endif

   do k=k1,k2
     if (read_p) then
       if (old) then
         call output(' trying to read old grain microstructure from binary file...')
       else
         call output(' trying to read grain microstructure from binary file...')
       endif
       if (.not.file_grain(.true.,ionum,k,n,old)) go to 99
       if (old) then
         call get_data(k,grain_area,grain_area_p=.true.)
         call get_data(k,ms,ms_p=.true.)
         if (.not.construct_ms(size(grain_area),grain_area,k,.true.,ms)) go to 99
         call set_data(k,ms,ms_p=.true.)
        endif
       if (.not.init_temp_scaling(il=k)) go to 99
       if (.not.init_llg(il=k)) go to 99
       call get_data(k,spin,spin_p=.true.)
       call get_data(k,tc,tc_p=.true.)
       call get_data(k,mag,mag_p=.true.)
       call get_data(k,lambdaMs,lambdaMs_p=.true.)
       do i1=1,n
         mag(:,i1)=1.d0; lambdaMs(i1)=3.d0*tc(i1)/(1.343427768d-4*(spin+1.d0))
       enddo
       if (.not.file_field(.true.,ionum,k,n,old)) go to 99
       call output(' successful read of grain microstructure!')
     else
       call output(' trying to write grain microstructure to binary file...')
       if (.not.file_grain(.false.,ionum,k,n,.false.)) go to 99
       if (.not.file_field(.false.,ionum,k,n,.false.)) go to 99
       call output(' successful write of grain microstructure!')
     endif
   enddo
   ok_p=.true.
   return

  99 continue
   ok_p=.false.
   if (read_p) then; call output(' failure reading grain microstructure file!')
   else; call output(' failure writing grain microstructure file!'); endif
   return
   
  end function

  subroutine show_media()
   use data_m, only:show_data,get_size
   use construct_grain_m, only: show_grain
   integer::nx,ny,n_lay
   real(8)::dx,dy

   call get_size(nx,ny,n_lay,dx,dy)
   do nx=1,n_lay
     call show_data(nx)
     call show_grain(nx)
   enddo

  end subroutine

  function generate_media(iarea_p,ihk_p) result (ok_p)
   use data_m, only: set_data, get_data, get_geometry, init_data, get_size
   use construct_grain_m, only: construct_microstructure,construct_ms,construct_tc,construct_tau
   use temp_scaling_m, only: init_temp_scaling
   use field_m, only: construct_field
   use llg_m, only: init_llg

   logical,intent(in),optional::iarea_p,ihk_p

   integer::nx,ny,n_lay
   real(8)::dx,dy
   integer::i,num_grain,n
   logical::ok_p,area_p,hk_p
   integer,pointer::grid2grain(:,:),grain2grid(:,:),grain2gridindex(:)
   real(8),pointer::grain_area(:),ms(:),tc(:),mag(:,:),lambdaMs(:),spin,tau(:),grain_vol(:)

   area_p=.true.; hk_p=.true.
   if (present(iarea_p)) area_p=iarea_p
   if (present(ihk_p)) hk_p=ihk_p

   ok_p=.true.
   call get_size(nx,ny,n_lay,dx,dy)
   do i=1,n_lay
     if (area_p) then
       call get_data(i,grid2grain,grid2grain_p=.true.)
       call get_data(i,grain2grid,grain2grid_p=.true.)
       call get_data(i,grain2gridindex,grain2gridindex_p=.true.)
       call get_data(i,grain_area,grain_area_p=.true.)
     
       ok_p=construct_microstructure(nx,ny,dx,dy,i,.false.,num_grain,grid2grain,grain2grid,grain2gridindex,grain_area)
       if (ok_p) then
         call set_data(i,grid2grain,grid2grain_p=.true.)
         call set_data(i,grain2grid,grain2grid_p=.true.)
         call set_data(i,grain2gridindex,grain2gridindex_p=.true.)
         call set_data(i,grain_area,grain_area_p=.true.)
! if (lbound(grain2gridindex,dim=1).ne.0) then; print *,'generate_media';read *;endif

         ok_p=init_data(i,num_grain)
         call get_data(i,ms,ms_p=.true.)
         ok_p=construct_ms(num_grain,grain_area,i,.false.,ms)
         if (ok_p) then
           call set_data(i,ms,ms_p=.true.)
           call get_data(i,tc,tc_p=.true.)
           call get_data(i,tau,tau_p=.true.)
           call get_data(i,grain_vol,grain_vol_p=.true.)
           ok_p=(construct_tc(num_grain,grain_vol,i,.false.,tc).and.construct_tau(num_grain,grain_area,i,.false.,tau))
!          ok_p=(construct_tc(num_grain,grain_area,i,.false.,tc).and.construct_tau(num_grain,grain_area,i,.false.,tau))
           if (ok_p) then
             call set_data(i,tc,tc_p=.true.)
             call set_data(i,tau,tau_p=.true.)
             ok_p=init_temp_scaling(il=i)
             ok_p=init_llg(il=i)
             call get_data(i,spin,spin_p=.true.)
             call get_data(i,mag,mag_p=.true.)
             call get_data(i,lambdaMs,lambdaMs_p=.true.)
             do n=1,num_grain
               mag(:,n)=1.d0; lambdaMs(n)=3.d0*tc(n)/(1.343427768d-4*(spin+1.d0))
             enddo
             ok_p=construct_field(nx,ny,dx,dy,i,hk_p,.true.)
           endif
         endif
       endif
     else
       ok_p=construct_field(nx,ny,dx,dy,i,hk_p,.true.)
     endif
     
   enddo
  end function

end module grain_m

