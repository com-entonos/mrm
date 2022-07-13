module solve1d_m
 !this module solves tridiagonal or cyclic tridiagonal systems
 !parallel pipeds or spheres

 implicit none

 private
 public::solve1d
  interface solve1d
     module procedure solve1d_tri_cyc
  end interface

contains
  subroutine solve1d_tri_cyc(a,b,c,r,x,n,print_p)
    integer,intent(in)::n
    real(8),intent(in)::a(n),b(n),c(n),r(n)
    real(8),intent(out)::x(n)
    logical,intent(in)::print_p
    real(8)::alpha,beta,fact,gamma,bb(n),u(n),z(n)
    integer::i

    if (a(1).eq.0.d0.and.c(n).eq.0.d0) then
      call solve1d_tri(a,b,c,r,x,n,print_p)                 !straight forward tridiagonal, not periodic
    else                                                !cyclic tridiagonal
      if (print_p) then
           do i=1,n
             print "('>>',i4,' a,b,c:',3(e17.10,x))",i,a(i),b(i),c(i)
           enddo
      endif
      beta=a(1); alpha=c(n)
      gamma=-b(1)
      bb=b; bb(1)=b(1)-gamma; bb(n)=b(n)-alpha*beta/gamma
      call solve1d_tri(a,bb,c,r,x,n,print_p)
      if (print_p) print *,'# of zero x:',count(x.eq.0.d0),n
      u=0.d0; u(1)=gamma; u(n)=alpha
      call solve1d_tri(a,bb,c,u,z,n,print_p)
      fact=(x(1)+beta*x(n)/gamma)/(1.d0+z(1)+beta*z(n)/gamma)
      if (print_p) print *,'# of zero z:',count(z.eq.0.d0),n,fact
      x=x-fact*z
      if (print_p) then
           do i=1,n
             print "('>>',i4,' a,b,c,x,a+b+c:',10(e17.10,x))",i,a(i),b(i),c(i),x(i),a(i)*x(max(1,i-1))+b(i)*x(i)+c(i)*x(min(n,i+1))
           enddo
      endif
    endif
  end subroutine solve1d_tri_cyc

  subroutine solve1d_tri(a,b,c,r,u,n,print_p)       !straight tridiagonal solve
    integer,intent(in)::n
    real(8),intent(in)::a(n),b(n),c(n),r(n)
    real(8),intent(out)::u(n)
    logical,intent(in)::print_p
    real(8)::gamm(n),bet
    integer::j


            if (b(1).eq.0.d0) stop
    bet=b(1)
    u(1)=r(1)/bet
    do j=2,n
      gamm(j)=c(j-1)/bet
      bet=b(j)-a(j)*gamm(j)
            if (bet.eq.0.d0) stop
      u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
    do j=n-1,1,-1
      u(j)=u(j)-gamm(j+1)*u(j+1)
    enddo
    if (print_p) print *,'# of zero u:',count(u.eq.0.d0),n

  end subroutine solve1d_tri

end module solve1d_m
                       
