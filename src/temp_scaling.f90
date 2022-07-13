
!this module creates layers of grains' microstructure, Hk, Hx, Tc, etc.
!this module scales physical parameters according to temperature

module temp_scaling_m
  use io_m
  implicit none
  private
  public find_temp_scaling,init_temp_scaling,file_temp_scaling,set_temp_scaling,show_temp_scaling,plot_temp_scaling

  type::temp_table_s
    real(8),pointer::temp(:)=>null()
    real(8),pointer::scaling(:)=>null()
    real(8),pointer::val(:)=>null()
  end type temp_table_s

  type::temp_table_ss
    type(temp_table_s),pointer::p
  end type temp_table_ss

  type::temp_scaling_s
    real(8)::tc_nom=1.d22,hkexp=-1.d0,msexp=-1.d0
    logical::use_mag_p=.false.
    real(8),pointer::mag(:,:)=>null()
    type(temp_table_ss),pointer::table(:)=>null()
  end type temp_scaling_s

  type::temp_scaling_ss
    type(temp_scaling_s),pointer::p
  end type temp_scaling_ss

  type(temp_scaling_ss),pointer,save::temp_scaling(:)=>null()

  character(5),dimension(4),save::pstr=(/'Ms   ','alpha','Hk   ','Hx   '/)
contains

  function set_temp_scaling(layer,tc_nom,hkexp,msexp,use_mag_p) result (ok_p)
  
   integer,intent(in)::layer
   real(8),optional,intent(in)::tc_nom,hkexp,msexp
   logical,optional,intent(in)::use_mag_p
   logical::ok_p
!  integer::i

   ok_p=.true.
   call create_temp_scaling_ss(layer)
   if (present(tc_nom))    temp_scaling(layer)%p%tc_nom=tc_nom
   if (present(use_mag_p)) temp_scaling(layer)%p%use_mag_p=use_mag_p
   if (present(hkexp))     temp_scaling(layer)%p%hkexp=hkexp
   if (present(msexp))     temp_scaling(layer)%p%msexp=msexp

!    print *,layer,associated(temp_scaling(layer)%p),associated(temp_scaling(layer)%p%table(1)%p),associated(temp_scaling(layer)%p%table(2)%p),associated(temp_scaling(layer)%p%table(3)%p),associated(temp_scaling(layer)%p%table(4)%p),temp_scaling(layer)%p%msexp,temp_scaling(layer)%p%hkexp,temp_scaling(layer)%p%use_mag_p,temp_scaling(layer)%p%tc_nom
!  do i=1,layer
!    print *,i,associated(temp_scaling(i)%p),associated(temp_scaling(i)%p%table(1)%p),associated(temp_scaling(i)%p%table(2)%p),associated(temp_scaling(i)%p%table(3)%p),associated(temp_scaling(i)%p%table(4)%p)
!  enddo

  contains

   subroutine create_temp_scaling_ss(layer)
    integer,intent(in)::layer
    type(temp_scaling_ss),pointer::new_temp_scaling(:)
    integer::n,n_l,i

    n_l=0; if (associated(temp_scaling)) n_l=size(temp_scaling)

    if (layer.gt.n_l) then
      allocate(new_temp_scaling(layer))
      do n=1,n_l
        new_temp_scaling(n)%p=>temp_scaling(n)%p
      enddo
      do n=n_l+1,layer
          allocate(new_temp_scaling(n)%p)
          allocate(new_temp_scaling(n)%p%table(4))
          do i=1,4; allocate(new_temp_scaling(n)%p%table(i)%p); enddo
      enddo
      if (associated(temp_scaling)) deallocate(temp_scaling)
      temp_scaling=>new_temp_scaling
    endif
   end subroutine

  end function

  subroutine show_temp_scaling()
    use io_m, only: output

    integer::k,indx,n
    character(200)::os
    type(temp_table_s),pointer::t
    type(temp_scaling_s),pointer::p

    do k=1,size(temp_scaling)
      p=>temp_scaling(k)%p
      write(os,"(' temperature scaling for layer #',i3,':')") k
      call output(trim(os))
      do indx=1,4
        t=>p%table(indx)%p
        if (pstr(indx)(1:2).eq.'Ms'.and.p%use_mag_p) then
          write(os,"('  using dynamic scaling')")
        elseif (pstr(indx)(1:2).eq.'Ms'.and.p%msexp.ge.0.d0) then
          write(os,"('  Ms ~ (1-T/Tc)^',1p,e12.5)") p%msexp
        elseif (pstr(indx)(1:2).eq.'Hk'.and.p%hkexp.ge.0.d0) then
          write(os,"('  Hk ~ (1-T/Tc)^',1p,e12.5)") p%hkexp
        elseif (.not.associated(t%temp)) then
          write(os,"('  no ',a,' temperature scaling')") trim(pstr(indx))
        else
          write(os,"('  ',a,' temperature scaling:')") trim(pstr(indx))
          call output(trim(os))
          write(os,"(a)") '   temp (K)             '//trim(pstr(indx))//'(T)/'//trim(pstr(indx))//'(T=0K)'
          do n=1,size(t%temp)
            call output(trim(os))
            write(os,"(1p,2(2x,e19.12))") t%temp(n)*p%tc_nom,t%scaling(n)
          enddo
        endif
        call output(trim(os))
      enddo
    enddo
  end subroutine

  function file_temp_scaling(ionum,indx,ilayer) result (ok_p)
    use io_m, only: output, read_comments
    integer,intent(in)::ionum,indx
    integer,intent(in),optional::ilayer

    logical::ok_p
    integer::k1,k2,k,n,i,j
    character(200)::os
    real(8),allocatable::dat(:,:,:)
    type(temp_scaling_s),pointer::p
    type(temp_table_s),pointer::t

    do k=1,size(temp_scaling)
      p=>temp_scaling(k)%p
      do i=2,4
        do j=1,i-1
          if (associated(p%table(i)%p%temp,p%table(j)%p%temp)) then
            allocate(p%table(i)%p%temp(size(p%table(j)%p%temp)))
!           p%table(j)%p%temp=p%table(j)%p%temp
            p%table(i)%p%temp=p%table(j)%p%temp
          endif
        enddo
      enddo
    enddo

    ok_p=.false.
    if (indx.lt.1.or.indx.gt.size(pstr)) return

    call output(' trying to read '//trim(pstr(indx))//'(T) scaling table')

    call read_comments(ionum)
    read(ionum,*,err=99,end=99) n
    if (n.lt.1) then; call output('  no entries in table!'); return; endif
    call read_comments(ionum)
    allocate(dat(n,2,2));dat=0.d0
    do k=1,n
      read(ionum,*,err=99,end=99) dat(k,1,2),dat(k,2,2)
    enddo

   !let's sort the table
    dat(1,1,1)=dat(1,1,2); dat(1,2,1)=dat(1,2,2)
    k=1
    do i=2,n
      j=1
      do while (j.le.k)
        if (dat(i,1,2).lt.dat(j,1,1)) exit
        j=j+1
      enddo
      do k1=k,j,-1
        dat(k1+1,1,1)=dat(k1,1,1); dat(k1+1,2,1)=dat(k1,2,1)
      enddo
      dat(j,1,1)=dat(i,1,2); dat(j,2,1)=dat(i,2,2); 
      k=k+1
    enddo

   !eliminate temperature duplicates
    k=0
    do i=n,2,-1
      if (dat(i-1,1,1).eq.dat(i,1,1)) then
        k=k+1
        dat(i-1,2,1)=0.5d0*(dat(i-1,2,1)+dat(i,2,1))
        do j=i+1,n
          dat(j-1,1,1)=dat(j,1,1); dat(j-1,2,1)=dat(j,2,1)
        enddo
      endif
    enddo
    n=max(1,n-k)
   !get rid of scaling duplicates at beginning
    if (dat(1,2,1).eq.dat(2,2,1)) then
      do i=2,n; dat(i-1,1,1)=dat(i,1,1); dat(i-1,2,1)=dat(i,2,1); enddo
      n=n-1
    endif
   !get rid of scaling duplicates at end
    if (n.gt.1) then
      if (dat(n-1,2,1).eq.dat(n,2,1)) n=n-1
    endif

    if (dat(n,1,1).gt.5.d0) then
       write(os,"(' WARNING: T for ',a,'(T) scaling is normalized by ',1p,e12.5)") trim(pstr(indx)),dat(n,1,1)
       call output(trim(os))
       dat(1:n,1,1)=dat(1:n,1,1)/dat(n,1,1)
    endif

    k1=1; k2=size(temp_scaling);
    if (present(ilayer)) then
      if (ilayer.gt.0.and.ilayer.le.k2) then; k1=ilayer;k2=ilayer; endif
    endif


    p=>temp_scaling(k1)%p
    t=>p%table(indx)%p
    allocate(t%temp(n),t%scaling(n))
    t%temp(1:n)=dat(1:n,1,1); t%scaling(1:n)=dat(1:n,2,1)
    do k=k1+1,k2
      temp_scaling(k)%p%table(indx)%p%temp=>t%temp; temp_scaling(k)%p%table(indx)%p%scaling=>t%scaling
    enddo

    deallocate(dat)

    if (k1.eq.k2.and.size(temp_scaling).ne.1) then
      write(os,"(' success reading ',a,'(T) scaling table for layer #',i3)") trim(pstr(indx)),k1
    else
      write(os,"(' success reading ',a,'(T) scaling table for all',i3,' layer(s)')") trim(pstr(indx)),size(temp_scaling)
    endif
    call output(trim(os))
    call output('  reduced table:')
    call output('    K                    '//trim(pstr(indx))//'(T)/'//trim(pstr(indx))//'(T=0K)')
    do i=1,n
      write(os,"(1p,2(2x,e19.12))") t%temp(i)*p%tc_nom,t%scaling(i)
      call output(trim(os))
    enddo
    ok_p=.true.
    return

 99 continue
    call output(' failure reading '//trim(pstr(indx))//'(T) scaling table!')
    if (allocated(dat)) deallocate(dat)
    return
  end function

  subroutine plot_temp_scaling(i,ms,alpha,hk,hx,tc0,grain_area) 
   use data_m, only:get_data
   use graphix, only:draw
   use plot_m, only: get_plot
   use io_m, only: output
   integer,intent(in)::i
   integer::j,n
   logical::ok_p
   real(8),pointer::tc0(:),hx(:),ms(:),hk(:,:),alpha,grain_area(:)
   real(8),allocatable::tgrain(:),tc(:),yd(:,:)
   real(8)::dx,xmax,x
   character(200)::os
   logical::plot_p(size(pstr))
   type(temp_table_s),pointer::t

   if (.not.associated(tc0).or..not.get_plot(i,hamr=.true.)) return
   ok_p= init_temp_scaling(i)
   n=size(tc0); allocate(tgrain(n),tc(n),yd(n,size(pstr))); yd=0.d0; plot_p=.false.
   xmax=1.05d0*maxval(tc0); dx=xmax/max(1,n-1)
   do j=1,n; tgrain(j)=(j-1)*dx; enddo

   do j=1,n
     tc=tc0(j)
     ok_p=find_temp_scaling(i,tgrain,tc)
     t=>temp_scaling(i)%p%table(1)%p
     if (associated(ms).and.associated(t%temp)) then; x=ms(j)
       write(os,"(a,i3)") trim(pstr(1))//'(T) for layer',i
       if (j.eq.1) then; call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15,xmin=0.d0,xmax=xmax,ymin=0d0,ymax=1.05d0*maxval(ms))
       else;             call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15)
       endif
       yd(:,1)=yd(:,1)+t%val*x*grain_area(j); plot_p(1)=.true.
     endif
     t=>temp_scaling(i)%p%table(2)%p
     if (associated(alpha).and.associated(t%temp)) then; x=alpha
       write(os,"(a,i3)") trim(pstr(2))//'(T) for layer',i
       if (j.eq.1) then; call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15,xmin=0.d0,xmax=xmax,ymin=0d0,ymax=1.05d0*alpha)
       else;             call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15)
       endif
       yd(:,2)=yd(:,2)+t%val*x*grain_area(j); plot_p(2)=.true.
     endif
     t=>temp_scaling(i)%p%table(3)%p
     if (associated(hk).and.associated(t%temp)) then; x=(hk(1,j)+hk(2,j))*1.d-4
       write(os,"(a,i3)") trim(pstr(3))//'(T) for layer',i
       if (j.eq.1) then; call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15,xmin=0.d0,xmax=xmax,ymin=0d0,ymax=1.05d-4*maxval(hk(1,:)+hk(2,:)))
       else;             call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15)
       endif
       yd(:,3)=yd(:,3)+t%val*x*grain_area(j); plot_p(3)=.true.
     endif
     t=>temp_scaling(i)%p%table(4)%p
     if (associated(hx).and.associated(t%temp)) then; x=hx(j)*1.d-4
       write(os,"(a,i3)") trim(pstr(4))//'(T) for layer',i
       if (j.eq.1) then; call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15,xmin=0.d0,xmax=xmax,ymin=0d0,ymax=1.05d-4*maxval(hx))
       else;             call draw(trim(os),x=tgrain,y=t%val*x,add=.true.,color=15)
       endif
       yd(:,4)=yd(:,4)+t%val*x*grain_area(j); plot_p(4)=.true.
     endif
   enddo
   do j=1,size(pstr)
     if (plot_p(j)) then
       write(os,"(a,i3)") trim(pstr(j))//'(T) for layer',i
       call draw(trim(os),x=tgrain,y=yd(:,j)/sum(grain_area),add=.true.,color=1,legend='average')
     else
       write(os,"(a,i3)") ' no '//trim(pstr(j))//'(T) for layer',i; call output(os)
     endif
   enddo
   deallocate(tc,tgrain,yd)
   ok_p= init_temp_scaling(i)

  end subroutine

  function init_temp_scaling(il) result (ok_p)
   use data_m, only:get_data
   use io_m, only:output
   integer,intent(in),optional::il

   type(temp_scaling_s),pointer::p
   real(8),pointer::data(:),data2(:,:)
   logical::ok_p
   integer::i,j,i1,i2,j1
   character(200)::os

   i1=1; i2=size(temp_scaling)
   if (present(il)) then
     if (il.ge.il.and.il.le.i2) then; i1=il; i2=il; endif
   endif
   ok_p=.false.
   do i=i1,i2
     p=>temp_scaling(i)%p
     call get_data(i,data,ms_scale_p=.true.);    p%table(1)%p%val=>data
     call get_data(i,data,alpha_scale_p=.true.); p%table(2)%p%val=>data
     call get_data(i,data,hk_scale_p=.true.);    p%table(3)%p%val=>data
     call get_data(i,data,hx_scale_p=.true.);    p%table(4)%p%val=>data
     call get_data(i,data2,mag_p=.true.);    p%mag=>data2
     do j=1,4
       if (.not.associated(p%table(j)%p%temp).and.associated(p%table(j)%p%val)) p%table(j)%p%val=1.d0
       if (pstr(j)(1:2).eq.'Ms'.and.p%use_mag_p) then
         write(os,"(' Ms using dynamic scaling for layer',i3)") i
         call output(trim(os))
       elseif (pstr(j)(1:2).eq.'Ms'.and.p%msexp.ge.0.d0) then
         write(os,"(' Ms ~ (1-T/Tc)^',1p,e12.5,0p,' for layer',i3)") p%msexp,i
         call output(trim(os))
       elseif (pstr(j)(1:2).eq.'Hk'.and.p%hkexp.ge.0.d0) then
         write(os,"(' Hk ~ (1-T/Tc)^',1p,e12.5,0p,' for layer',i3)") p%hkexp,i
         call output(trim(os))
       endif
     enddo
     do j=2,4
       if (associated(p%table(j)%p%temp)) then
         do j1=1,j-1
           if (associated(p%table(j1)%p%temp)) then
             if (size(p%table(j)%p%temp).eq.size(p%table(j1)%p%temp).and..not.associated(p%table(j)%p%temp,p%table(j1)%p%temp)) then
               if (all(p%table(j)%p%temp.eq.p%table(j1)%p%temp)) then
  !              deallocate(p%table(j)%p%temp)
                 p%table(j)%p%temp=>p%table(j1)%p%temp !possible memory leak if nothing else is pointing to table(j)%p%temp FIXME
                 write(os,"(' temperature scale for layer',i3,' is the same for ',a,' => ',a,' (good job!)')")i,trim(pstr(j)),trim(pstr(j1)); call output(os)
               endif
             endif
           endif
         enddo
       endif
     enddo
   enddo
   ok_p=.true.
!  read *
  end function

  function find_temp_scaling(layer,tgrain,tc,t_changed_p) result (ok_p)
   use data_m, only:get_data
   integer,intent(in)::layer
   real(8),intent(in)::tgrain(:),tc(:)
   logical,intent(in),optional::t_changed_p
   type(temp_scaling_s),pointer::p
   logical::ok_p,do_p
   integer::i,n,k1,k2,k,num_grain
   real(8)::x,y
   logical::found(4)
   integer,pointer::in
   logical,pointer::do_grain(:)=>null()

   ok_p=.false.; do_p=.true.
   if (present(t_changed_p)) do_p=t_changed_p
   if (layer.lt.1.or.layer.gt.size(temp_scaling)) return
   ok_p=.true.; p=>temp_scaling(layer)%p; if (p%tc_nom.gt.1.d10) return
   found=.false.; num_grain=size(tgrain)

   do i=1,4
     if (p%use_mag_p.and.pstr(i)(1:2).eq.'Ms') then
       call get_data(layer,in,in_p=.true.)
       p%table(i)%p%val(:)=min(1.d0,max(0.d0,abs(p%mag(in/3+1,:))))
     elseif (do_p.and.p%msexp.ge.0.d0.and.pstr(i)(1:2).eq.'Ms') then
       p%table(i)%p%val(:)=min(1.d0,max(0.d0,(1.d0-tgrain(:)/tc(:))))**p%msexp
     elseif (do_p.and.p%hkexp.ge.0.d0.and.pstr(i)(1:2).eq.'Hk') then
!      p%table(i)%p%val(:)=min(1.d0,max(0.d0,p%table(1)%p%val(:)**p%hkexp)) !assumes Ms is first table entry FIXME
       p%table(i)%p%val(:)=min(1.d0,max(0.d0,(1.d0-tgrain(:)/tc(:))))**p%hkexp
     elseif (do_p.and..not.found(i).and.associated(p%table(i)%p%temp)) then
       if (size(p%table(i)%p%temp).gt.1) then
         call get_data(layer,do_grain,do_grain_p=.true.)
         do n=1,num_grain
          if (do_grain(n)) then
           x=tgrain(n)/tc(n); k1=1; k2=size(p%table(i)%p%temp)
           if (x.le.p%table(i)%p%temp(1)) then
             y=0.d0
           elseif (x.ge.p%table(i)%p%temp(k2)) then
             y=0.d0; k1=k2
           else
             do while (k2-k1.gt.1)
               k=(k2+k1)/2
               if (x.lt.p%table(i)%p%temp(k)) then
                 k2=k
               else
                 k1=k
               endif
             enddo
             y=max(0.d0,min(1.d0,(x-p%table(i)%p%temp(k1))/(p%table(i)%p%temp(k2)-p%table(i)%p%temp(k1))))
           endif
           p%table(i)%p%val(n)=p%table(i)%p%scaling(k1)+y*(p%table(i)%p%scaling(k2)-p%table(i)%p%scaling(k1))
!  if (n.eq.1) print "(i5,1p,:50(x,e12.5))",i,tgrain(n),tc(n),x,p%table(i)%p%temp(k1),p%table(i)%p%temp(k2),p%table(i)%p%val(n)
!  print "(i1,i5,1p,:50(x,e12.5))",i,n,tgrain(n),tc(n),x,p%table(i)%p%temp(k1),p%table(i)%p%temp(k2),y,p%table(i)%p%val(n)
           do k=i+1,4
             if (associated(p%table(i)%p%temp,p%table(k)%p%temp)) then
                found(k)=.true. 
                p%table(k)%p%val(n)=p%table(k)%p%scaling(k1)+y*(p%table(k)%p%scaling(k2)-p%table(k)%p%scaling(k1))
!  if (n.eq.1) print "(i5,1p,:50(x,e12.5))",k,tgrain(n),tc(n),x,p%table(k)%p%temp(k1),p%table(k)%p%temp(k2),p%table(k)%p%val(n)
!  print "(i1,i5,1p,:50(x,e12.5))",k,n,tgrain(n),tc(n),x,p%table(k)%p%temp(k1),p%table(k)%p%temp(k2),y,p%table(k)%p%val(n)
             endif
           enddo
          endif
         enddo
       endif
     endif
   enddo

  end function

end module temp_scaling_m
