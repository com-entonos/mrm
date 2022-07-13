module mselem_m
 !this module calculates magnetostatic tensor elements and/or fields for
 !parallel pipeds or spheres

 implicit none

 private
 public::volume_averaged_n, volume_averaged_h, volume_averaged_sphere_n, volume_averaged_sphere_h, plane_to_plane, block_to_plane, &
         block_to_plane_h

contains

 function volume_averaged_n(rx,ry,rz,xs,ys,zs,xo,yo,zo,n) result (vol)
!calculate the volume average demag tensor element.
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
!  n = tensor element, i.e. hx = [nx(1,1)*Mx + nx(2,1)*My + nx(3,1)*Mz] (Oe cc/emu)
!                           hy = [ny(1,2)*Mx + ny(2,2)*My + ny(3,2)*Mz]
!                           hz = [nz(1,3)*Mx + nz(2,3)*My + nz(3,3)*Mz]
!  vol = 1 / volume of observation cell (cm^-3)

  real(8),intent(in)::rx,ry,rz,xs,ys,zs,xo,yo,zo
  real(8),intent(inout)::n(3,3)
  real(8)::vol
  real(8)::dxo(2),dyo(2),dzo(2)
  real(8)::dxs(2),dys(2),dzs(2),x,d,dx,dy,dz,t(3,3)

  integer::i1,i2,i3,i4,i5,i6,i

  dxo=xo*(/ -0.5d0, 0.5d0 /); dyo=xo*(/ -0.5d0, 0.5d0 /); dzo=zo*(/ -0.5d0, 0.5d0 /)
  dxs=xs*(/ -0.5d0, 0.5d0 /); dys=ys*(/ -0.5d0, 0.5d0 /); dzs=zs*(/ -0.5d0, 0.5d0 /)
  vol=1.d0/(xo*yo*zo); t=0.d0

  do i6=1,2  !source cell, x
   do i5=1,2  !source cell, y
    do i4=1,2  !source cell, z
     do i3=1,2  !observation cell, x
      dx=-rx+dxo(i3)-dxs(i6)
      do i2=1,2  !observation cell, y
       dy=-ry+dyo(i2)-dys(i5)
       do i1=1,2  !observation cell, z
        dz=-rz+dzo(i1)-dzs(i4)
        d=sqrt(dx*dx+dy*dy+dz*dz)
        if (d.ne.0.d0) then
          i=(-1)**(i1+i2+i3+i4+i5+i6-1)

          t(1,1)=t(1,1)+i*g1(dx,dy,dz,d)
          t(2,2)=t(2,2)+i*g1(dy,dz,dx,d)
          t(3,3)=t(3,3)+i*g1(dz,dx,dy,d)
          
          x=i*g2(dx,dy,dz,d)
          t(3,2)=t(3,2)+x
          t(2,3)=t(2,3)+x
          
          x=i*g2(dy,dz,dx,d)
          t(3,1)=t(3,1)+x
          t(1,3)=t(1,3)+x
          
          x=i*g2(dz,dx,dy,d)
          t(2,1)=t(2,1)+x
          t(1,2)=t(1,2)+x
        endif
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
! print *,t(:,1)
! print *,t(:,2)
! print *,t(:,3)
  n=n-t
  return

 contains

   function g1(x,y,z,d) result (r)
     real(8),intent(in)::x,y,z,d
     real(8)::r

     r =  x*y*z*atan(y*z/sign(max(1.d-300,abs(x)*d),x)) &
         +0.5d0*y*(z*z-x*x)*log(max(1.d-300,d-y)) &
         +0.5d0*z*(y*y-x*x)*log(max(1.d-300,d-z)) &
         +(y*y+z*z-2.d0*x*x)*d/6.d0

   end function g1

   function g2(z,x,y,d) result (r)
     real(8),intent(in)::x,y,z,d
     real(8)::r

     r = -x*y*z*log(max(1.d-300,d+z)) &
         +y*(y*y-3.d0*z*z)*log(max(1.d-300,d+x))/6.d0 &
         +x*(x*x-3.d0*z*z)*log(max(1.d-300,d+y))/6.d0 &
         +0.5d0*x*x*z*atan(y*z/sign(max(1.d-300,abs(x)*d),x)) &
         +0.5d0*y*y*z*atan(x*z/sign(max(1.d-300,abs(y)*d),y)) &
         +z*z*z*atan(x*y/sign(max(1.d-300,abs(z)*d),z))/6.d0+x*y*d/3.d0

   end function g2

 end function volume_averaged_n

 function volume_averaged_h(rx,ry,rz,xs,ys,zs,xo,yo,zo,m) result (h)
!calculate the volume average demag tensor element.
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
!  m = (mx,my,mz) = source block magnetization (emu/cc_
!  h = volume average field in observation block (Oe)
  real(8),intent(in)::rx,ry,rz,xs,ys,zs,xo,yo,zo,m(3)
  real(8)::h(3)

  real(8)::n(3,3),vol

  n=0.d0
  vol=volume_averaged_n(rx,ry,rz,xs,ys,zs,xo,yo,zo,n)
  n=n*vol
  h = (/ m(1)*n(1,1) + m(2)*n(2,1) + m(3)*n(3,1), &
         m(1)*n(1,2) + m(2)*n(2,2) + m(3)*n(3,2), &
         m(1)*n(1,3) + m(2)*n(2,3) + m(3)*n(3,3) /)
  return
 end function volume_averaged_h

 function volume_averaged_sphere_n(rx,ry,rz,radius,n) result (vol)
  real(8),intent(in)::rx,ry,rz,radius
  real(8),intent(inout)::n(3,3)
  real(8)::vol

  real(8)::e,nrm,x,t(3,3)

  nrm=4.d0*acos(-1.d0)/3.d0
  vol = nrm*radius*radius*radius
  t=0.d0

  e=sqrt(rx*rx+ry*ry+rz*rz)
  if (e.gt.1.d-20) then !not same sphere
    e=1.d0/e; x=abs(radius*e*radius*e*radius*e)
    t(1,1) = (3.d0*rx*e*rx*e-1.d0)*x
    t(2,1) =  3.d0*rx*e*ry*e      *x
    t(3,1) =  3.d0*rx*e*rz*e      *x

    t(1,2) =  3.d0*rx*e*ry*e      *x
    t(2,2) = (3.d0*ry*e*ry*e-1.d0)*x
    t(3,2) =  3.d0*ry*e*rz*e      *x

    t(1,3) =  3.d0*rx*e*rz*e      *x
    t(2,3) =  3.d0*ry*e*rz*e      *x
    t(3,3) = (3.d0*rz*e*rz*e-1.d0)*x
  else
    t(1,1) = -1.d0
    t(2,2) = -1.d0
    t(3,3) = -1.d0
  endif
  n=n+t*nrm*vol; vol=1.d0/vol
 end function volume_averaged_sphere_n

 function volume_averaged_sphere_h(rx,ry,rz,radius,m) result (h)
  real(8),intent(in)::rx,ry,rz,radius,m(3)
  real(8)::h(3)
  real(8)::vol,n(3,3)

  n=0.d0; vol=volume_averaged_sphere_n(rx,ry,rz,radius,n)
  n=n*vol
  h = (/ m(1)*n(1,1) + m(2)*n(2,1) + m(3)*n(3,1), &
         m(1)*n(1,2) + m(2)*n(2,2) + m(3)*n(3,2), &
         m(1)*n(1,3) + m(2)*n(2,3) + m(3)*n(3,3) /)
  return
  end function volume_averaged_sphere_h

 function block_to_plane(dxs,dys,dzs,nx,ny,dx,dy,dz,r) result(n)
  integer,intent(in)::nx,ny
  real(8),intent(in)::dx,dy,dz,r(3),dxs,dys,dzs
  real(8)::n(3,3,nx,ny)

  integer::i,j
  real(8)::vol

  n=0.d0
  do j=1,ny  !di =  observation
    do i=1,nx  !dj = observation
      vol=volume_averaged_n(r(1)-(i-0.5d0)*dx,r(2)-(j-0.5d0)*dy,r(3),dxs,dys,dzs,dx,dy,dz,n(:,:,i,j))
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
    enddo
  enddo
  n=n*vol
!  print *,nx,ny,vol,n(:,:,:,ny/2)
 end function block_to_plane

 subroutine block_to_plane_h(dxs,dys,dzs,nx,ny,dx,dy,dz,r,m,h)
  integer,intent(in)::nx,ny
  real(8),intent(in)::dx,dy,dz,r(3),dxs,dys,dzs,m(3)
  real(8),intent(out)::h(3,nx,ny)
  real(8)::n(3,3)

  integer::i,j
  real(8)::vol

  n=0.d0
  h=0.d0
  do j=1,ny  !di =  observation
    do i=1,nx  !dj = observation
      n=0.d0
      vol=volume_averaged_n(r(1)-(i-0.5d0)*dx,r(2)-(j-0.5d0)*dy,r(3),dxs,dys,dzs,dx,dy,dz,n)
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
      h(:,i,j) = (/ m(1)*n(1,1) + m(2)*n(2,1) + m(3)*n(3,1), &
                    m(1)*n(1,2) + m(2)*n(2,2) + m(3)*n(3,2), &
                    m(1)*n(1,3) + m(2)*n(2,3) + m(3)*n(3,3) /)
    enddo
  enddo
!  h=h/(vol*nx*ny)
  h=h*vol
 end subroutine block_to_plane_h

 function plane_to_plane(nx,ny,dx,dy,dz,dzo,rz,n_imgxi,n_imgyi) result(n)
  integer,intent(in)::nx,ny
  real(8),intent(in)::dx,dy,dz,rz,dzo
  integer,intent(in),optional::n_imgxi,n_imgyi
  real(8)::n(nx,ny,3,3)

  integer::i,j,nx_img,ny_img,n_imgx,n_imgy,idm,idm0,xp0,xp1,yp0,yp1
  real(8)::vol,t(3,3)

  n_imgx=0; n_imgy=0; idm=0; idm0=max(1,nint(nx*ny/10.))
  if (present(n_imgxi))n_imgx=max(0,n_imgxi)  !if > 0, then this is periodic
  if (present(n_imgyi))n_imgy=max(0,n_imgxi)  !if > 0, then this is periodic
 
  n=0.d0
  do j=-ny/2,ny/2  !di = source - observation
    yp0=-n_imgy*(1-max(0,j/((ny+1)/2))); yp1=n_imgy*(1+min(0,j/((ny+1)/2))) !fix for ny even, otherwise overcounts
    do i=-nx/2,nx/2  !dj = source - observation
      xp0=-n_imgx*(1-max(0,i/((nx+1)/2))); xp1=n_imgx*(1+min(0,i/((nx+1)/2))) !fix for nx even, otherwise overcounts
      idm=idm+1; if (mod(idm-1,idm0).gt.mod(idm,idm0)) print "(' ',i2,$)",10-idm/idm0
      t=0.d0
!write(44,*)i,j,xp0,xp1,yp0,yp1
!write(44,*)i,j,mod(nx-i,nx)+1,mod(ny-j,ny)+1
      do ny_img=yp0,yp1
        do nx_img=xp0,xp1
! print *,i+nx_img*nx,j+ny_img*ny
!if (i.eq.56) write(44,"(4(i3,x),1p,9(e12.5,x),:,0p,10(i3,x))") i,j,mod(nx-i,nx)+1,mod(ny-j,ny)+1,i*dx+nx_img*dx*nx,j*dy+ny_img*dy*ny,rz,dx,dy,dz,dx,dy,dzo,nx_img,ny_img,xp0,xp1,yp0,yp1
          vol=volume_averaged_n(i*dx+nx_img*dx*nx,j*dy+ny_img*dy*ny,rz,dx,dy,dz,dx,dy,dzo,t)
!function volume_averaged_n(rx,ry,rz,xs,ys,zs,xo,yo,zo,n) result (vol)
!calculate the volume average demag tensor element.
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
!  n = tensor element, i.e. hx = [nx(1,1)*Mx + nx(2,1)*My + nx(3,1)*Mz] (Oe cc/emu)
!                           hy = [ny(1,2)*Mx + ny(2,2)*My + ny(3,2)*Mz]
!                           hz = [nz(1,3)*Mx + nz(2,3)*My + nz(3,3)*Mz]
!  vol = 1 / volume of observation cell (cm^-3)
        enddo
      enddo
      n(mod(nx-i,nx)+1,mod(ny-j,ny)+1,:,:)=n(mod(nx-i,nx)+1,mod(ny-j,ny)+1,:,:)+t
! print *,mod(nx-i,nx)+1,mod(ny-j,ny)+1,i,j
    enddo
  enddo
 end function plane_to_plane
end module mselem_m
