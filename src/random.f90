
module random_m
  implicit none
  
  private
  public::random_s,random_create,random_copy,random_uniform,random_normal,random_normal3d,random_gaussian,random_gaussian_sphere,random_destroy, &
      restore_random, save_random, show_random

  !Mersenne twister algorithm
  integer,parameter::N=624,N1=N+1,M=397,MATA=-1727483681,UMASK=-2147483648,LMASK=2147483647,TMASKB=-1658038656,TMASKC=-272236544,mag01(0:1)=(/0,MATA/)
  type::random_s
    integer::mti,mt(0:N-1),iset
    real(8)::gset
    character(40)::name   !unique ID for
  end type random_s

  type::random_list_s
    type(random_s),pointer::generator
    type(random_list_s),pointer::next=>null()
  end type

  type(random_list_s),pointer,save::random_list=>null()

contains

  subroutine show_random()
    use io_m, only: output
    type(random_list_s),pointer::list
    character(200)::os
    list=>random_list
    call output('random number generator list:');
    do while (associated(list))
      write(os,"(' ',a,' (',i0,')')") trim(adjustl(list%generator%name)),list%generator%mti; call output(trim(os))
      list=>list%next
    enddo
  end subroutine

  function save_random(ionum,v_p) result (ok_p)
    use io_m, only: output
    integer,intent(in)::ionum
    logical,intent(in),optional::v_p
    logical::ok_p
    integer::n
    type(random_list_s),pointer::list
    character(200)::os

    ok_p=.true.
    n=0; list=>random_list;do while (associated(list)); n=n+1; list=>list%next; enddo
    write(os,"(' trying to save ',i0,' random generators...')") n; call output(trim(os))
    write(ionum) n
    list=>random_list
    do while (associated(list))
      write(ionum,iostat=n) list%generator%mti,list%generator%mt,list%generator%iset,list%generator%gset,list%generator%name
      ok_p=(ok_p.and.n.eq.0)
      if (n.eq.0) then
        write(os,"('  saved random generator ',a,' (',i0,')')") trim(adjustl(list%generator%name)),list%generator%mti
      else
        write(os,"('  >>> failed to save random generator ',a,' (',i0,')')") trim(adjustl(list%generator%name)),list%generator%mti
      endif
      if (present(v_p)) call output(trim(os))
      list=>list%next
    enddo
    call output(' saved all random generators!')
  end function

  function restore_random(ionum,v_p) result (ok_p)
    use io_m, only: output
    integer,intent(in)::ionum
    logical,intent(in),optional::v_p
    logical::ok_p
    integer::n,i,nr,nn
    type(random_s)::generator
    type(random_list_s),pointer::list
    character(200)::os

    nr=0; nn=0
    read(ionum) n
    write(os,"(' trying to restore ',i0,' random generators...')") n; call output(trim(os))
    do while (n.gt.0)
      read(ionum,iostat=i) generator%mti,generator%mt,generator%iset,generator%gset,generator%name
      nn=nn+1
      ok_p=(i.eq.0)
      if (.not.ok_p) then
        call output('  >>> error reading a random generators! <<<')
        return
      endif
      list=>random_list
      do while (associated(list))
        if (trim(adjustl(list%generator%name)).eq.trim(adjustl(generator%name))) then
          nr=nr+1
          write(os,"('  restored random generator ',a,' (',i0,')')") trim(adjustl(generator%name)),generator%mti; if (present(v_p)) call output(trim(os))
          list%generator%mti  = generator%mti
          list%generator%mt   = generator%mt
          list%generator%iset = generator%iset
          list%generator%gset = generator%gset
          exit
        endif
        list=>list%next
      enddo
      if (.not.associated(list)) then
        write(os,"('  WARNING- random generator (',a,') not in list! ignoring...')") trim(adjustl(generator%name))
        call output(trim(os))
      endif
      n=n-1
    enddo
    write(os,"(' restored ',i0,' out of ',i0,' random generators...')") nr,nn; call output(trim(os))
    ok_p=.true.
  end function

  subroutine add_random(r)
   type(random_s),pointer::r
   type(random_list_s),pointer::list
   if (associated(r)) then
     if (associated(random_list)) then
       list=>random_list
       do while (associated(list%next)); list=>list%next; enddo
       allocate(list%next); list%next%generator=>r; list%next%next=>null()
     else
       allocate(random_list); random_list%generator=>r; random_list%next=>null()
     endif
   endif
  end subroutine
  subroutine remove_random(r)
   type(random_s),pointer::r
   type(random_list_s),pointer::list,last

   if (associated(r).and.associated(random_list)) then
     list=>random_list
     do while (associated(list))
       if (associated(list%generator,r)) then  !this is the generator
         if (associated(list,random_list)) then
           random_list=>list%next
         else
           last%next=>list%next
         endif
         deallocate(list);exit
       endif
       last=>list; list=>list%next
     enddo
   endif
  end subroutine

  subroutine random_destroy(s)
   type(random_s),pointer::s
   call remove_random(s)
   if (associated(s)) deallocate(s)
   nullify(s)
   return
  end subroutine random_destroy

  function random_copy(name,s) result (r)
 !return a random number between [0, 1), based on Mersenne twister alg
   character(*),intent(in)::name
   type(random_s),pointer::s,r

   if (associated(s)) then
     allocate(r)
     r=s
     r%name=trim(adjustl(name))
     call add_random(r)
   else
     nullify(r)
   endif
  end function random_copy

  function random_create(name,iseedi) result (s)
 !create a random number generator and seeds it, based on Mersenne twister alg
    character(*),intent(in)::name
    integer,intent(in),optional::iseedi
    type(random_s),pointer::s
    integer::iseed,i

    iseed=4357; if (present(iseedi)) iseed=iseed+abs(iseedi)

    allocate(s)
    s%mt(0)=iand(iseed,-1)
    do i=1,N-1
      s%mt(i)=iand(69069*s%mt(i-1),-1)
    enddo
    s%iset=0; s%mti=i
    s%name=trim(adjustl(name))
    call add_random(s)
  end function random_create

  function random_uniform(s) result (r)
 !return a random number between [0, 1), based on Mersenne twister alg
   type(random_s),pointer::s
   real(8)::r
   integer::j,i

   if (s%mti.ge.N) then
     do i=0,N-M-1
       j=ior(iand(s%mt(i),UMASK),iand(s%mt(i+1),LMASK))
       s%mt(i)=ieor(ieor(s%mt(i+M),ishft(j,-1)),mag01(iand(j,1)))
     enddo
     do i=N-M,N-2
       j=ior(iand(s%mt(i),UMASK),iand(s%mt(i+1),LMASK))
       s%mt(i)=ieor(ieor(s%mt(i+M-N),ishft(j,-1)),mag01(iand(j,1)))
     enddo
     j=ior(iand(s%mt(N-1),UMASK),iand(s%mt(0),LMASK))
     s%mt(N-1)=ieor(ieor(s%mt(M-1),ishft(j,-1)),mag01(iand(j,1)))
     s%mti=0
   endif

   j=s%mt(s%mti); s%mti=s%mti+1
   j=ieor(j,ishft(j,-11)); j=ieor(j,iand(ishft(j,7),TMASKB))
   j=ieor(j,iand(ishft(j,15),TMASKC)); j=ieor(j,ishft(j,-18))
   if (j.lt.0) then
     r=(dble(j)+2.d0**32)/(2.d0**32)
   else
     r=dble(j)/(2.d0**32)
   endif
  end function random_uniform

  function random_normal(s) result(r)
 !return a normal (gaussian w/ average=0, sigma=1) random number
   type(random_s),pointer::s
   real(8)::r,rsq,v1,v2,f

   if (s%iset.ne.0) then
     r=s%gset; s%iset=0
   else
     v1=2.d0*random_uniform(s)-1.d0; v2=2.d0*random_uniform(s)-1.d0; rsq=v1*v1+v2*v2
     do while (rsq.ge.1.d0 .or. rsq .eq. 0.d0)
       v1=2.d0*random_uniform(s)-1.d0; v2=2.d0*random_uniform(s)-1.d0; rsq=v1*v1+v2*v2
     enddo
     f=sqrt(-2.d0*log(rsq)/rsq); s%gset=v1*f; r=v2*f; s%iset=1
   endif
  end function random_normal

  function random_normal3d(s) result(r3)
 !return a vector of normal (gaussian w/ average=0, sigma=1) random number
   type(random_s),pointer::s
   real(8)::r3(3),r,rsq,v1,v2,f
   integer::i

   do i=1,3
     if (s%iset.ne.0) then
       r=s%gset; s%iset=0
     else
       v1=2.d0*random_uniform(s)-1.d0; v2=2.d0*random_uniform(s)-1.d0; rsq=v1*v1+v2*v2
       do while (rsq.ge.1.d0 .or. rsq .eq. 0.d0)
         v1=2.d0*random_uniform(s)-1.d0; v2=2.d0*random_uniform(s)-1.d0; rsq=v1*v1+v2*v2
       enddo
       f=sqrt(-2.d0*log(rsq)/rsq); s%gset=v1*f; r=v2*f; s%iset=1
     endif
     r3(i)=r
   enddo
  end function random_normal3d

  function random_gaussian(s,sigmai,avei) result(r)
 !return a guassian random number
   type(random_s),pointer::s
   real(8)::r
   real(8),intent(in),optional::sigmai,avei

   real(8)::sigma,ave
   sigma=1.d0; ave=0.d0; if (present(sigmai)) sigma=sigmai; if (present(avei)) ave=avei
   r=random_normal(s)*sigma+ave
  end function random_gaussian

  function random_gaussian_sphere(s,sigma) result(r)
 !return an gaussian random number w/ average=0 and sigma=sigma on a sphere
 !that is, probablility that theta is between theta and theta+dtheat =
 !  A exp(-.5 (theta/sigma)^2) sin(theta) dtheta, where A is a
 !  normalization constant
 ! notice that sin(theta) is from the metric (i.e. we're on a sphere)
   type(random_s),pointer::s
   real(8),intent(in)::sigma
   real(8)::r,ss

   ss=sigma/(sigma*sigma+1.d-300); ss=ss*ss
   r=acos(max(-1.d0,min(1.d0, 2.d0*random_uniform(s)-1.d0)))
   do while (random_uniform(s).gt.exp(-0.5d0*(r*r*ss)))
     r=acos(max(-1.d0,min(1.d0, 2.d0*random_uniform(s)-1.d0)))
   enddo
  end function random_gaussian_sphere

end module random_m
