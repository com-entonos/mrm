module sto_m
 !this module calculates magnetostatic tensor elements and/or fields for
 !parallel pipeds or spheres

 implicit none

 private
 public:: sto

 interface sto
  module procedure sto_rect,sto_tilt
 end interface

 type sto_ss
   integer::nx,ny
   real(8)::dx,dy,dz,r(3),thickness,width,sheight,angle,t
   real(8)::dxs,dys,dzs
   real(8),dimension(:,:,:,:),allocatable::n!(3,3)
   type(sto_ss),pointer::next=>null()
 end type sto_ss

 type(sto_ss),pointer,save::stod !=>null()
 real(8),save::t=-1.d0
contains

  function calc_h(n,m) result (h)
    use plot_m, only: plot_strip
    real(8),intent(in)::n(:,:,:,:),m(3)
    real(8)::h(3,size(n,3),size(n,4))

    integer::i,j
    integer::nx,ny

! n=0.d0
  nx=size(n,3); ny=size(n,4)
  h=0.d0
  do j=1,ny
    do i=1,nx
!  (rx,ry,rz) = source position - observation position vector, center to center  (cm)
!  (xs,ys,zs) = size of source block  (cm)
!  (xo,yo,zo) = size of observation block (cm)
      h(:,i,j) = (/ m(1)*n(1,1,i,j) + m(2)*n(2,1,i,j) + m(3)*n(3,1,i,j), &
                    m(1)*n(1,2,i,j) + m(2)*n(2,2,i,j) + m(3)*n(3,2,i,j), &
                    m(1)*n(1,3,i,j) + m(2)*n(2,3,i,j) + m(3)*n(3,3,i,j) /)
    enddo
  enddo

  t=t+1.d0
  call plot_strip('1xy1',1,t,(/sqrt(dot_product(h(:,56,48),h(:,56,48))),h(:,56,48)/))
  call plot_strip('2xy2',1,t,(/sqrt(dot_product(h(:,76,48),h(:,76,48))),h(:,76,48)/))
  call plot_strip('3xy3',1,t,(/sqrt(dot_product(h(:,56,68),h(:,56,68))),h(:,56,68)/))
  call plot_strip('4xy4',1,t,(/sqrt(dot_product(h(:,76,68),h(:,76,68))),h(:,76,68)/))
  call plot_strip('cxyc',1,t,(/sqrt(dot_product(h(:,nx/2,ny/2),h(:,nx/2,ny/2))),h(:,nx/2,ny/2)/))
  call plot_strip('clyc',1,t,(/sqrt(dot_product(h(:,nx/2,48),h(:,nx/2,48))),h(:,nx/2,48)/))
 end function

 subroutine sto_rect(dxs,dys,dzs,nx,ny,dx,dy,dz,r,m,h)
  use mselem_m,only: block_to_plane
  integer,intent(in)::nx,ny
  real(8),intent(in)::dx,dy,dz,r(3),dxs,dys,dzs,m(3)
  real(8),intent(out)::h(3,nx,ny)
  type(sto_ss),pointer::isto,nsto=>null(),lsto=>null()

  lsto=>stod
  do while (associated(lsto))
    nsto=>lsto
    if (lsto%dxs.eq.dxs.and.lsto%dys.eq.dys.and.lsto%dzs.eq.dzs.and. &
            lsto%r(1).eq.r(1).and.lsto%r(2).eq.r(2).and.lsto%r(3).eq.r(3)) exit
    lsto=>lsto%next
  enddo
  if (.not.associated(lsto)) then 
    allocate(isto); 
    if (.not.associated(stod)) then
       stod=>isto
    else
       nsto%next=>isto
    endif
    allocate(isto%n(3,3,nx,ny)); isto%n=0.d0
    isto%dxs=dxs; isto%dys=dys; isto%dzs=dzs; isto%r=r
    isto%n = block_to_plane(dxs,dys,dzs,nx,ny,dx,dy,dz,r)
!   print *,isto%n(:,:,:,ny/2)
  else
    isto=>lsto
  endif
  h = calc_h(isto%n,m)
 end subroutine sto_rect
!subroutine sto_tilt(dxs,dys,dzs,nx,ny,dx,dy,dz,r,m,t,h)
 subroutine sto_tilt(thickness, width, sheight, angle,nx,ny,dx,dy,dz,r,m,h)
  integer,intent(in)::nx,ny
  real(8),intent(in)::dx,dy,dz,r(3),thickness,width,sheight,angle,m(3)
  real(8),intent(out)::h(3,nx,ny)

  integer::i,j
  real(8)::vol
! real(8)::n(3,3)

  !angle in radian
  !r(1),r(2)= center of tilted slab in x,y; r(3)=gap between bottom of tilted slab and sampling xy plane

  real(8)::yc,zc,mm,a1,a0,x0,x1,z0,z1,y00,y01,y10,y11,xp,x,y

  type(sto_ss),pointer::isto,nsto=>null(),lsto=>null()

print *,thickness, width, sheight, angle,nx,ny,dx,dy,dz,r,m

  lsto=>stod
  do while (associated(lsto))
    nsto=>lsto
    if (lsto%thickness.eq.thickness.and.lsto%width.eq.width.and.lsto%sheight.eq.sheight.and.lsto%angle.eq.angle.and. &
            lsto%r(1).eq.r(1).and.lsto%r(2).eq.r(2).and.lsto%r(3).eq.r(3)) exit
    lsto=>lsto%next
  enddo
  if (.not.associated(lsto)) then 
    allocate(isto); 
    if (.not.associated(stod)) then
       stod=>isto
    else
       nsto%next=>isto
    endif
    allocate(isto%n(3,3,nx,ny)); isto%n=0.d0
    isto%thickness=thickness; isto%width=width; isto%sheight=sheight; isto%r=r; isto%angle=angle

    ! dxs = width
    ! dys = thickness
    ! dzs = slanted length (i.e. dz = dzs*cos(atan(kkkkkk
    ! the z-planes
    mm = tan(angle)
    z0 = sheight*cos(angle); zc = r(3) + z0*0.5d0; z1=r(3)+z0; z0 = r(3)
    x0 = -width/2.d0 + r(1) + nx * dx * 0.5d0; x1 = width/2.d0 + r(1) + nx * dx * 0.5d0
    yc = ny * dy * 0.5d0 + r(2)
!   a1 = yc - zc*mm; a0 = a1 - thickness/2.d0; a1 = a1 + thickness/2.d0
!   y00 = z0 * mm + a0
!   y01 = z0 * mm + a1
!   y10 = z1 * mm + a0
!   y11 = z1 * mm + a1
    a0 = yc - thickness/2.d0; a1 = yc + thickness/2.d0
    y00 = (z0-zc) * mm + a0
    y01 = (z0-zc) * mm + a1
    y10 = (z1-zc) * mm + a0
    y11 = (z1-zc) * mm + a1
    ! the x-planes
print *,yc,zc,mm
print *,z0,z1
print *,x0,x1
print *,a0,a1
print *,y00,y01
print *,y10,y11
    !need to compute n!
    print *,'calculating tilted slab'
    do j=1,ny  !di =  observation
 print *,ny-j+1
      y=(j-0.5d0)*dy
      do i=1,nx  !dj = observation
        x=(i-0.5d0)*dx
        ! xy sheets, i.e. Mz
        isto%n(3,1,i,j)=  zs(y-y11,x-x1,-z1) -  zs(y-y10,x-x1,-z1) -  zs(y-y11,x-x0,-z1) +  zs(y-y10,x-x0,-z1) &
                        - zs(y-y01,x-x1,-z0) +  zs(y-y00,x-x1,-z0) +  zs(y-y01,x-x0,-z0) -  zs(y-y00,x-x0,-z0)
        isto%n(3,2,i,j)=  zs(x-x1,y-y11,-z1) -  zs(x-x0,y-y11,-z1) -  zs(x-x1,y-y10,-z1) +  zs(x-x0,y-y10,-z1) &
                        - zs(x-x1,y-y01,-z0) +  zs(x-x0,y-y01,-z0) +  zs(x-x1,y-y00,-z0) -  zs(x-x0,y-y00,-z0)
        isto%n(3,3,i,j)= zsz(x-x1,y-y11,-z1) - zsz(x-x0,y-y11,-z1) - zsz(x-x1,y-y10,-z1) + zsz(x-x0,y-y10,-z1) &
                        -zsz(x-x1,y-y01,-z0) + zsz(x-x0,y-y01,-z0) + zsz(x-x1,y-y00,-z0) - zsz(x-x0,y-y00,-z0)
        !parallelogram in yz, i.e. Mx
        isto%n(1,1,i,j)= integrate(xsxi,z0,z1,(x-x1),y,a1)-integrate(xsxi,z0,z1,(x-x1),y,a0) &
                        -integrate(xsxi,z0,z1,(x-x0),y,a1)+integrate(xsxi,z0,z1,(x-x0),y,a0)
!        isto%n(1,2,i,j)= xsy(x-x1,y,a1,z1) - xsy(x-x1,y,a1,z0) - xsy(x-x1,y,a0,z1) + xsy(x-x1,y,a0,z0) &
!                        -xsy(x-x0,y,a1,z1) + xsy(x-x0,y,a1,z0) + xsy(x-x0,y,a0,z1) - xsy(x-x0,y,a0,z0)
!  if (abs(isto%n(1,2,i,j)).gt.1.d0) print *,i,j, &
!                         xsy(x-x1,y,a1,z1),- xsy(x-x1,y,a1,z0),- xsy(x-x1,y,a0,z1),+ xsy(x-x1,y,a0,z0),&
!                        -xsy(x-x0,y,a1,z1),+ xsy(x-x0,y,a1,z0),+ xsy(x-x0,y,a0,z1),- xsy(x-x0,y,a0,z0), &
!                        (a1-y+mm*(z0-zc)+z0/mm)/sqrt(1.d0+1.d0/(mm*mm)),sqrt((x-x0)**2+(a1-y+mm*(z0-zc))**2+z0*z0), &
!                        (a0-y+mm*(z0-zc)+z0/mm)/sqrt(1.d0+1.d0/(mm*mm)),sqrt((x-x0)**2+(a0-y+mm*(z0-zc))**2+z0*z0)
        isto%n(1,2,i,j)= integrate(xsyi,z0,z1,(x-x1),y,a1)-integrate(xsyi,z0,z1,(x-x1),y,a0) &
                        -integrate(xsyi,z0,z1,(x-x0),y,a1)+integrate(xsyi,z0,z1,(x-x0),y,a0)
        isto%n(1,3,i,j)= integrate(xszi,z0,z1,(x-x1),y,a1)-integrate(xszi,z0,z1,(x-x1),y,a0) &
                        -integrate(xszi,z0,z1,(x-x0),y,a1)+integrate(xszi,z0,z1,(x-x0),y,a0)
!        x= integrate(xsyi,z0,z1,(x-x1),y,a1)-integrate(xsyi,z0,z1,(x-x1),y,a0) &
!                        -integrate(xsyi,z0,z1,(x-x0),y,a1)+integrate(xsyi,z0,z1,(x-x0),y,a0)
!        if (min(abs(x),abs(isto%n(1,2,i,j)))/max(abs(isto%n(1,2,i,j)),abs(x)).gt.1.d-3) print *,i,j,x,min(abs(x),abs(isto%n(1,2,i,j)))/max(abs(isto%n(1,2,i,j)),abs(x))
      enddo
    enddo
  else
    isto=>lsto
  endif
  h = calc_h(isto%n,m)

!  h=h/(vol*nx*ny)
! h=h*vol
 contains
   function zs(p,q,r) result (s)  ! Hx: p=(y-yp), q=(x-xp); Hy: p=(x-xp), q=(y-yp) from a rectangle charge sheet peprnedicular to z axis
     real(8),intent(in)::p,q,r
     real(8)::s
     s = -log(p+sqrt(p*p+q*q+r*r))
   end function
   function zsz(p,q,r) result (s)
     real(8),intent(in)::p,q,r
     real(8)::s
     s = atan(p*q/(r*sqrt(p*p+q*q+r*r)))
   end function
   function xsy(r,y,a,zp) result (s)  ! Hy: q=(y-yp), p=(z-zp), r=(x-xp)
     real(8),intent(in)::r,a,y,zp  !amy=a-y
     real(8)::s,x
!    s = log(zp + mm * (amy-mm*zp) / sqrt(1.d0+mm*mm) + &
!        sqrt( r*r + amy*amy + 2.d0*zp*mm*amy + (1.d0+mm*mm)*zp*zp ))
!    s = log(zp + mm * (-amy+mm*zp) / sqrt(1.d0+mm*mm) + &
!        sqrt( r*r + amy*amy - 2.d0*zp*mm*amy + (1.d0+mm*mm)*zp*zp ))
     x = a - y + mm*(zp-zc)
!    print *,(x + zp/mm) / sqrt(1.d0+1.d0/(mm*mm)), sqrt( r*r + x*x + zp*zp )
     s = log( (x + zp/mm) / sqrt(1.d0+1.d0/(mm*mm)) + sqrt( r*r + x*x + zp*zp ) )
!    s = log( (mm*x + zp) / sqrt(1.d0+mm*mm) + sqrt( r*r + x*x + zp*zp ) )
!    s = log( x/sqrt(1.d0+1.d0/(mm*mm)) + zp / sqrt(1.d0+mm*mm) + sqrt( r*r + x*x + zp*zp ) )
!    x = y - a + mm*(zp-zc)
!    s = -log( (mm*x - zp) / sqrt(1.d0+mm*mm) + sqrt( r*r + x*x + zp*zp ) )
!    s = -log( x/sqrt(1.d0+1.d0/(mm*mm)) - zp / sqrt(1.d0+mm*mm) + sqrt( r*r + x*x + zp*zp ) )
   end function

   function xsyi(zp,r,y,a) result (s)
     real(8),intent(in)::zp,r,y,a
     real(8)::s,x
!    s = -r*(yma-mm*zp)/(r*r+zp*zp)/sqrt(r*r+zp*zp+(yma-mm*zp)*(yma-mm*zp))
     x = y - a - mm*(zp-zc)
     s = 1.d0/sqrt(r*r+zp*zp+x*x)
   end function
   function xszi(zp,r,y,a) result (s)
     real(8),intent(in)::zp,r,y,a
     real(8)::s,x
!    s = zp*(yma-mm*zp)/(r*r+zp*zp)/sqrt(r*r+zp*zp+(yma-mm*zp)*(yma-mm*zp))
     x = y - a - mm*(zp-zc)
     s = zp*x/(r*r+zp*zp)/sqrt(r*r+zp*zp+x*x)
   end function
   function xsxi(zp,r,y,a) result (s)
     real(8),intent(in)::zp,r,y,a
     real(8)::s,x
!    s = -r*(yma-mm*zp)/(r*r+zp*zp)/sqrt(r*r+zp*zp+(yma-mm*zp)*(yma-mm*zp))
     x = y - a - mm*(zp-zc)
     s = -r*x/(r*r+zp*zp)/sqrt(r*r+zp*zp+x*x)
   end function

   function integrate(func,a,b,p1,p2,p3) result (s)
    real(8),intent(in)::a,b,p1,p2,p3
    real(8)::func
    external func

    real(8),parameter::eps=1.d-8
    integer,parameter::jmax=30 !20
    integer::j
    real(8)::olds,s
    real(8)::ss

    olds=-1.d30
    do j=1,jmax
      ss=olds
      call trapzd(func,a,b,p1,p2,p3,j,s)
!     if (j.gt.5.and.(abs(s-olds).lt.eps*max(1.d-5,abs(olds)).or.(s.eq.0.d0.and.olds.eq.0.d0))) print *,j,jmax
!     if (j.gt.5.and.(abs(s-olds).lt.eps*max(1.d-5,abs(olds)).or.(s.eq.0.d0.and.olds.eq.0.d0))) return
      if (j.gt.5.and.(abs(s-olds).lt.eps*max(1.d-5,abs(olds)).or.(abs(s).lt.1.d-5.and.abs(olds).lt.1.d-5))) return
      olds=s
    enddo
print *,'many steps in integrate',abs(s-ss),eps*(ss),s,ss

   end function integrate
   subroutine trapzd(func,a,b,p1,p2,p3,n,s)
    real(8),intent(in)::a,b,p1,p2,p3
    integer,intent(in)::n
    real(8),intent(inout)::s
    real(8)::func
    external func

    integer::it,j
    real(8)::del,sum,tnm,x

    if (n.eq.1) then
      s=0.5d0*(b-a)*(func(a,p1,p2,p3)+func(b,p1,p2,p3))
    else
      it = 2**(n-2)
      tnm = it
      del=(b-a)/tnm
      x=a+0.5d0*del
      sum=0.d0
      do j=1,it
        sum=sum+func(x,p1,p2,p3)
        x=x+del
      enddo
      s=0.5d0*(s+(b-a)*sum/tnm)
    endif
   end subroutine trapzd
 end subroutine sto_tilt

end module sto_m
