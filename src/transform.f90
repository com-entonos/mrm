module transform_m
!module transform_m: this module does some simple rotational transformations of
!coordinate systems. it also declares a pure function for a vector cross product

  implicit none

  private
  public::rotate,rotated,cross

  interface rotate
    module procedure rotater,rotate_1,rotate_2
  end interface

contains

  pure function cross(m,h) result(t)
    real(8),intent(in)::m(3),h(3)
    real(8)::t(3)
    t = (/ m(2)*h(3)-m(3)*h(2), m(3)*h(1)-m(1)*h(3), m(1)*h(2)-m(2)*h(1) /)
  end function cross

  function rotater(phi,theta, phi_s, theta_s) result (x)
 ! given a initial vector orientated by (phi,theta). we now scatter that vector
 ! by phi_s, theta_s (measured w.r.t the initial orientation), we end up w/ a
 ! unit vector (x,y,z). phi_s is typically random (uniformly distributed), so it doesn't matter were we
 ! measure from. we assume it's from a vector drawn from the z axis to the
 ! initial orientation.
 ! this routine assumes phi, theta, phi_s, theta_s are in radians!

   implicit none

   real(8),intent(in)::phi, theta, phi_s, theta_s
   real(8)::x(3)

  !local variables
   real(8)::cp,sp,ct,st,xs,ys,zs

   cp = cos(phi);     sp = sin(phi);      ct = cos(theta);    st = sin(theta)
   zs = sin(theta_s); xs = cos(phi_s)*zs; ys = sin(phi_s)*zs; zs = cos(theta_s)

!  z = zs*ct - xs*st
!  x = (xs*ct + zs*st)*cp - ys*sp
!  y = (xs*ct + zs*st)*sp + ys*cp
   x(3) = zs*ct - xs*st
   x(1) = (xs*ct + zs*st)*cp - ys*sp
   x(2) = (xs*ct + zs*st)*sp + ys*cp
  end function rotater

  function rotated(phi,theta, phi_s, theta_s) result (x)
 ! same as routine rotate above, but
 ! this routine assumes phi, theta, phi_s, theta_s are in degrees

   implicit none

   real(8),intent(in)::phi, theta, phi_s, theta_s
   real(8)::x(3)

  !local variables
   real(8)::f

   f=acos(-1.d0)/180.d0
   x=rotate( phi*f, theta*f, phi_s*f, theta_s*f)
  end function rotated

  function rotate_q(w,m) result (t)
 ! rotate vector m using quaterion w, i.e. m-> w ** (0,m(1),m(2),m(3)) ** w*
 ! this is used if more than a single rotation is needed
   real(8),intent(in)::w(4),m(3)
   real(8)::t(3),w1,w2,w3,w4
   w1=w(1)*w(1); w2=w(2)*w(2); w3=w(3)*w(3); w4=w(4)*w(4)
   t = (/ 2.d0*(m(3)*(w(1)*w(3) + w(2)*w(4)) + m(2)*(w(2)*w(3) - w(1)*w(4))) + m(1)*(w1 + w2 - w3 - w4), &
          2.d0*(m(1)*(w(1)*w(4) + w(2)*w(3)) + m(3)*(w(3)*w(4) - w(1)*w(2))) + m(2)*(w1 + w3 - w2 - w4), &
          2.d0*(m(2)*(w(1)*w(2) + w(3)*w(4)) + m(1)*(w(2)*w(4) - w(1)*w(3))) + m(3)*(w1 + w4 - w3 - w2) /)
  end function

  function rotate_1(theta,u,m,norm) result (t)
 ! rotate vector m about vector u by angle of theta; theta>0 means rotation obeys right hand rule (thumb along u)
 ! use quaternions to do the rotation: p = (cos(theta/2),u(1) sin(theta/2), u(2) sin(theta/2), u(3) sin(theta/2))
 !  then m -> p ** (0, m(1), m(2), m(3)) ** p*
   real(8),intent(in)::theta,u(3),m(3)
   logical,intent(in),optional::norm
   real(8)::t(3),ct,st,dp
   t=u; 
   if (present(norm)) then; ct=sqrt(dot_product(t,t)); if (ct.ne.0.d0) then; t=t/ct; else; t=m; return; endif; endif
   ct=cos(theta); st=sin(theta); dp=dot_product(t,m)
   t = (/ t(1) * dp - ((dp - m(1) * t(1)) * t(1) + m(1) * (t(1) * t(1) - 1.d0)) * ct + (m(3) * t(2) - m(2) * t(3)) * st, &
          t(2) * dp - ((dp - m(2) * t(2)) * t(2) + m(2) * (t(2) * t(2) - 1.d0)) * ct + (m(1) * t(3) - m(3) * t(1)) * st, &
          t(3) * dp - ((dp - m(3) * t(3)) * t(3) + m(3) * (t(3) * t(3) - 1.d0)) * ct + (m(2) * t(1) - m(1) * t(2)) * st /)
!  t=t/(sqrt(dot_product(t,t))+1.d-300)
!  t = (/ t(1) * dp - ((m(2) * t(2) + m(3) * t(3)) * t(1) + m(1) * (t(1) * t(1) - 1.d0)) * ct + (m(3) * t(2) - m(2) * t(3)) * st, &
!         t(2) * dp - ((m(1) * t(1) + m(3) * t(3)) * t(2) + m(2) * (t(2) * t(2) - 1.d0)) * ct + (m(1) * t(3) - m(3) * t(1)) * st, &
!         t(3) * dp - ((m(1) * t(1) + m(2) * t(2)) * t(3) + m(3) * (t(3) * t(3) - 1.d0)) * ct + (m(2) * t(1) - m(1) * t(2)) * st /)
! or we could use rotate_q, but above appears to have less round-off and should be faster
!  t=rotate_q( (/ cos(theta/2.d0), t(1)*sin(theta/2.d0), t(2)*sin(theta/2.d0), t(3)*sin(theta/2.d0) /),m)
  end function

  function rotate_2(theta0,u0,theta1,u1,m,norm0,norm1) result (t)
 ! rotate vector m first about vector u0 by theta0, then around u1 by theta1
   real(8),intent(in)::theta0,u0(3),theta1,u1(3),m(3)
   logical,intent(in),optional::norm0,norm1
   real(8)::t(3),t0(3),t1(3),ct0,st0,ct1,st1
   st1=0.5d0*theta1; st0=0.5d0*theta0
   t0=u0; if (present(norm0)) then; ct0=sqrt(dot_product(t0,t0)); if (ct0.ne.0.d0) then; t0=t0/ct0; else; t0=(/0.d0,0.d0,1.d0/); st0=0.d0; endif; endif
   t1=u1; if (present(norm1)) then; ct0=sqrt(dot_product(t1,t1)); if (ct0.ne.0.d0) then; t1=t1/ct0; else; t1=(/0.d0,0.d0,1.d0/); st1=0.d0; endif; endif
   if (st0.eq.0.d0.and.st0.eq.st1) then; t=m; return; endif
   ct1=cos(st1); st1=sin(st1)
   ct0=cos(st0); st0=sin(st0)
   t = rotate_q( (/ ct0*ct1 - dot_product(t0,t1) * st0 * st1, &
                    t0(1)*ct1*st0+(t1(1)*ct0+(t1(2)*t0(3)-t0(2)*t1(3))*st0)*st1, &
                    t0(2)*ct1*st0+(t1(2)*ct0+(t1(3)*t0(1)-t0(3)*t1(1))*st0)*st1, &
                    t0(3)*ct1*st0+(t1(3)*ct0+(t1(1)*t0(2)-t0(1)*t1(2))*st0)*st1 /), m )
!  t=t/(sqrt(dot_product(t,t))+1.d-300)
  ! or we could just call rotate_1 twice
  !t = rotate_1(theta1,t1,rotate_1(theta0,t0,m))
  end function


end module transform_m

