module tensor_m
 implicit none
 private
 public::calculate_tensor,set_tensor,get_N_layer

 logical,save::periodic_x=.true.,periodic_y=.true.,sul_image=.false.,head_image=.false.
 character(200),save::filename='mrm_N.dat'
 integer,save::num_images=2
 real(8),save::sul_mu=1.d0,fly_height=8.d-7,space_to_sul=40.d-7,head_mu=1.d0
 real(8),allocatable,target,save::sn2d(:,:,:,:,:)
contains

 subroutine set_tensor(ihead_mu,isul_mu,ifly_height,ispace_to_sul,ifilename,inum_images)
  real(8),intent(in)::ihead_mu,isul_mu,ifly_height,ispace_to_sul
  integer,intent(in)::inum_images
  character(*),intent(in)::ifilename
  sul_mu=isul_mu; fly_height=ifly_height; space_to_sul=ispace_to_sul
  filename=trim(adjustl(ifilename))
  num_images=inum_images; head_mu=ihead_mu
  sul_image=(sul_mu.gt.1.d0); head_image=(head_mu.gt.1.d0)
 end subroutine

 subroutine get_N_layer(x)
  real(8),pointer::x(:,:,:,:,:)
  nullify(x); if (allocated(sn2d)) x=>sn2d
 end subroutine

 function file_tensor(read_p,n2d,nx,ny,n_lay,head_p,sul_p,imx,imy,dx,dy,dz,il_space) result (ok_p)
  use io_m, only: output,add_file,delete_file

  real(8),intent(inout)::n2d(:,:,:,:,:)
  logical,intent(in)::sul_p,head_p
  integer,intent(in)::nx,ny,n_lay,imx,imy
  real(8),intent(in)::dx,dy,dz(n_lay),il_space(n_lay)
  
  logical::isul_p,ihead_p
  integer::inx,iny,in_lay,iimx,iimy
  real(8)::idx,idy,idz(n_lay),iil_space(n_lay),isul_mu,ifly_height,ispace_to_sul,ihead_mu

  logical::ok_p,read_p
  integer::ionum,i
  character(200)::os

  ok_p=.false.
  if (read_p) then
    inquire(file=trim(filename),exist=ok_p)
    if (.not.ok_p) return
    write(os,"('  trying to read magnetostatic tensor from file: ',a)") trim(filename); call output(trim(os))
    ionum=add_file(trim(filename)); if (ionum.lt.1) return
    call delete_file(ionum)
    open(ionum,file=trim(filename),status='old',form='unformatted',action='read',position='rewind')
    read(ionum,iostat=i) ihead_p,isul_p,inx,iny,in_lay,iimx,iimy,idx,idy,ihead_mu,isul_mu,ifly_height,ispace_to_sul,idz,iil_space
    ok_p=(i.eq.0.and.(ihead_p.eqv.head_p).and.(isul_p.eqv.sul_p).and.inx.eq.nx.and.iny.eq.ny.and.in_lay.eq.n_lay.and.iimx.eq.imx.and.iimy.eq.imy.and.  &
          idx.eq.dx.and.idy.eq.dy.and.isul_mu.eq.sul_mu.and.ihead_mu.eq.head_mu.and.ifly_height.eq.fly_height.and.ispace_to_sul.eq.space_to_sul.and.  &
          all(idz.eq.dz).and.all(iil_space.eq.il_space))
    if (ok_p) read(ionum,iostat=i) n2d
    ok_p=(i.eq.0.and.ok_p)
    close(ionum)
    if (ok_p) then
      call output('   -successful read!')
    else
      call output('   -unsuccessful read, will compute...')
    endif
  else
    write(os,"('  trying to write magnetostatic tensor to file: ',a)") trim(filename); call output(trim(os))
    ionum=add_file(trim(filename)); if (ionum.lt.1) return
    call delete_file(ionum)
    open(ionum,file=trim(filename),status='unknown',form='unformatted')
    write(ionum,iostat=i) head_p,sul_p,nx,ny,n_lay,imx,imy,dx,dy,head_mu,sul_mu,fly_height,space_to_sul,dz,il_space; ok_p=(i.eq.0)
    if (ok_p) write(ionum,iostat=i) n2d
    ok_p=(i.eq.0.and.ok_p)
    close(ionum)
    if (ok_p) then
      call output('   -successful save!')
    else
      call output('   -unsuccessful save, oh no!')
    endif
  endif

 end function

 subroutine calculate_tensor(n2d)
  use io_m, only: output
  use mselem_m, only: plane_to_plane !,volume_averaged_n
  use data_m, only: get_geometry
  real(8),intent(out)::n2d(:,:,:,:,:)

  logical::head_p,sul_p,ok_p
  integer::nx,ny,n_lay,imx,imy,ks,ko
  character(200)::os
  real(8)::dx,dy,zsb,zst,zob,zot,mus,mush
  real(8),allocatable::dz(:),il_space(:),nn(:,:,:,:),hdz(:)

  nx=size(n2d,1); ny=size(n2d,2); n_lay=size(n2d,4)

  imx=0; if (periodic_x) imx=num_images
  imy=0; if (periodic_y) imy=num_images
  write(os,"(' generating magnetostatic tensor: nx=',i4,', ny=',i4,', layers=',i3,' (',i3,',',i3,')')") nx,ny,n_lay,imx,imy

  allocate(dz(n_lay),il_space(n_lay),nn(nx,ny,3,3),hdz(n_lay))

  call get_geometry(dx,dy,dz,il_space)

  mus=0.d0; mush=0.d0
  sul_p=(sul_image.and.sul_mu.ne.1.d0)
  if (sul_p) then
    write(os,"(a,', image in sul')") trim(os)
  else
    write(os,"(a,', no image in sul')") trim(os)
  endif
  head_p=(head_image.and.head_mu.ne.1.d0)
  if (head_p) then
    write(os,"(a,', image in (infinite) head')") trim(os)
  else
    write(os,"(a,', no image in head')") trim(os)
  endif
  if (sul_p.and.head_p) then
    call output('')
    call output('**** WARNING **** imaging in both head and sul- code does not do imaging of images!')
    call output('')
  endif
  if (sul_p) mus=max(0.d0,(max(1.d0,sul_mu)-1.d0)/(max(1.d0,sul_mu)+1.d0))
  if (head_p) mush=max(0.d0,(max(1.d0,head_mu)-1.d0)/(max(1.d0,head_mu)+1.d0))

  if (imx.gt.0.and.imy.gt.0) then
    write(os,"(a,', periodic in plane')") trim(os)
  elseif (imx.gt.0) then
    write(os,"(a,', periodic in x, not y')") trim(os)
  elseif (imy.gt.0) then
    write(os,"(a,', periodic in y, not x')") trim(os)
  else
    write(os,"(a,', not periodic in plane')") trim(os)
  endif
  call output(trim(os))

  if (.not.file_tensor(.true.,n2d,nx,ny,n_lay,head_p,sul_p,imx,imy,dx,dy,dz,il_space)) then !try to read in tensor
    n2d=0.d0

    hdz=fly_height
    do ks=2,n_lay
!     hdz(ks)=hdz(ks-1)+dz(ks)+il_space(ks)
      hdz(ks)=hdz(ks-1)+dz(ks-1)+il_space(ks)
    enddo
!   hdz=0.5d0*dz  !hdz(k) = distance from head to center of layer k
    hdz=hdz+0.5d0*dz  !hdz(k) = distance from head to center of layer k
    

    zsb=space_to_sul
    do ks=n_lay,1,-1       !source layer
      zst=zsb+dz(ks)
      zob=space_to_sul
      do ko=n_lay,1,-1     !observation layer
        zot=zob+dz(ko)

        print "('   computing layer #',i3,' on layer #',i3,':',$)",ko,ks
 !function plane_to_plane(nx,ny,dx,dy,dz,dzo,rz,n_imgxi,n_imgyi) result(n)
        nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dz(ko),0.5d0*(zst-zot+zsb-zob),imx,imy)
        n2d(:,:,1:3,ko,ks)=nn(:,:,:,1)
        n2d(:,:,4:6,ko,ks)=nn(:,:,:,2)
        n2d(:,:,7:9,ko,ks)=nn(:,:,:,3)

        if (sul_p) then
          call output('  -done')
          print "('                    and its SUL image:',$)"
!print *,ks,ko,mus,0.5d0*(-zst-zot-zsb-zob); read *
          nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dz(ko),0.5d0*(-zst-zot-zsb-zob),imx,imy)
          nn(:,:,1:2,:)=-mus*nn(:,:,1:2,:); nn(:,:,3,:)=mus*nn(:,:,3,:) !xy component flip sign
          n2d(:,:,1:3,ko,ks)=n2d(:,:,1:3,ko,ks)+nn(:,:,:,1)
          n2d(:,:,4:6,ko,ks)=n2d(:,:,4:6,ko,ks)+nn(:,:,:,2)
          n2d(:,:,7:9,ko,ks)=n2d(:,:,7:9,ko,ks)+nn(:,:,:,3)
        endif
        if (head_p) then
          call output('  -done')
          print "('        and its (infinite) head image:',$)"
!         nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dz(ko),hdz(ks)-hdz(ko),imx,imy)
          nn=plane_to_plane(nx,ny,dx,dy,dz(ks),dz(ko),hdz(ks)+hdz(ko),imx,imy)
          nn(:,:,1:2,:)=-mush*nn(:,:,1:2,:); nn(:,:,3,:)=mush*nn(:,:,3,:) !xy component flip sign
          n2d(:,:,1:3,ko,ks)=n2d(:,:,1:3,ko,ks)+nn(:,:,:,1)
          n2d(:,:,4:6,ko,ks)=n2d(:,:,4:6,ko,ks)+nn(:,:,:,2)
          n2d(:,:,7:9,ko,ks)=n2d(:,:,7:9,ko,ks)+nn(:,:,:,3)
        endif

        n2d(:,:,:,ko,ks)=n2d(:,:,:,ko,ks)/(dx*dy*dz(ko)*nx*ny)

!do j=1,ny; do i=1,nx
!print "(6(i4),1p,2(x,e19.12))",ks,ko,j,i,imx,imy,dx*dy*dz(ko)*nx*ny,0.5d0*(zst-zot+zsb-zob)
!print "(1p,3(e12.5,x))",n2d(i,j,1:3,ko,ks)
!print "(1p,3(e12.5,x))",n2d(i,j,4:6,ko,ks)
!print "(1p,3(e12.5,x))",n2d(i,j,7:9,ko,ks)
!write(55,"(6(i4),1p,2(x,e19.12))")ks,ko,j,i,imx,imy,dx*dy*dz(ko)*nx*ny,0.5d0*(zst-zot+zsb-zob)
!write(55,"(1p,3(e12.5,x))")n2d(i,j,1:3,ko,ks)
!write(55,"(1p,3(e12.5,x))")n2d(i,j,4:6,ko,ks)
!write(55,"(1p,3(e12.5,x))")n2d(i,j,7:9,ko,ks)
!enddo; enddo;

        zob=zot+il_space(ko)
        call output('  -done')
      enddo
      zsb=zst+il_space(ks)
    enddo
! do ks=1,n_lay; do ko=1,n_lay
! do j=1,ny; do i=1,nx
!!write(55,"(6(i4),1p,2(x,e19.12))")ks,ko,i,j !,imx,imy,dx*dy*dz(ko)*nx*ny,0.5d0*(zst-zot+zsb-zob)
! write(55,"(1p,3(e12.5,x),0p,4(i4))")n2d(i,j,1:3,ko,ks),i,j,ko,ks
! write(55,"(1p,3(e12.5,x),0p,4(i4))")n2d(i,j,4:6,ko,ks),i,j,ko,ks
! write(55,"(1p,3(e12.5,x),0p,4(i4))")n2d(i,j,7:9,ko,ks),i,j,ko,ks
! enddo; enddo;
! enddo; enddo
    ok_p=file_tensor(.false.,n2d,nx,ny,n_lay,head_p,sul_p,imx,imy,dx,dy,dz,il_space)
  endif

  deallocate(dz,il_space,nn,hdz)
  if (allocated(sn2d)) deallocate(sn2d)
  allocate(sn2d(nx,ny,9,n_lay,n_lay)); sn2d=n2d

 end subroutine

end module tensor_m
