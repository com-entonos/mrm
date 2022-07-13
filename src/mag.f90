module mag_m
  use random_m, only: random_s
  implicit none
  private
  public::set_mag,generate_mag,start_mag,file_mag

  type::mag_s
    integer::layer
    logical::ac_p, saturate_p
    real(8)::m(3),xmin,xmax,ymin,ymax
    type(mag_s),pointer::next=>null()
  end type mag_s

  type(mag_s),pointer,save::mag=>null()

  type::magr_s
    type(random_s),pointer::random_mag=>null()
  end type magr_s
  type(magr_s),allocatable,save::mr(:)
contains

  function start_mag(il,seed) result (ok_p)
   use random_m, only: random_create
   integer,intent(in)::il
   integer,intent(in),optional:: seed
   logical::ok_p
   character(200)::str

   ok_p=.true.
   if (.not.allocated(mr)) then; allocate(mr(il)); return; endif
   if (present(seed)) then
     write(str,"(i5)") il;str=adjustl(str)
     mr(il)%random_mag=>random_create('random_mag'//trim(str),seed)
   endif

  end function

  function generate_mag(il) result (ok_p)
   use random_m, only: random_uniform
   use data_m, only:get_data,get_size
   use field_m, only:get_easy_axis
   use llg_m, only:need_to_initialize_llg
   integer,intent(in),optional::il
   logical::ok_p

   integer::i1,i2,i,nx,ny,n_lay,n,n1,j
   real(8)::dx,dy,x
   real(8),pointer::m(:,:),k(:,:)
   integer,pointer::id,grain2gridindex(:),grain2grid(:,:)
   type(mag_s),pointer::m1

   ok_p=.true.
   call get_size(nx,ny,n_lay,dx,dy)

   i1=1; i2=n_lay
   if (present(il)) then
     if (il.gt.0.and.il.le.n_lay) then; i1=il; i2=il; endif
   endif

   do i=i1,i2
     call get_data(i,m,m_p=.true.)
     call get_data(i,id,in_p=.true.)
     call get_data(i,grain2gridindex,grain2gridindex_p=.true.)
     call get_data(i,grain2grid,grain2grid_p=.true.)
     call get_easy_axis(i,k)
     call need_to_initialize_llg()
!if (lbound(grain2gridindex,dim=1).ne.0) then; print *,'generate_mag';read *;endif
     do n=1,size(k,2)
       m(1+id:3+id,n)=k(:,n)
     enddo

     m1=>mag
     do while (associated(m1))
       if (m1%layer.eq.0.or.m1%layer.eq.i) then
         do n=1,size(k,2); n1=0
           do j=grain2gridindex(n-1)+1, grain2gridindex(n)
             if ( (grain2grid(1,j)-0.5d0)*dx.ge.m1%xmin .and. (grain2grid(1,j)-0.5d0)*dx.le.m1%xmax .and. &
                  (grain2grid(2,j)-0.5d0)*dy.ge.m1%ymin .and. (grain2grid(2,j)-0.5d0)*dy.le.m1%ymax) n1=n1+1
           enddo
           if (random_uniform(mr(i)%random_mag).le.dble(n1)/(grain2gridindex(n)-grain2gridindex(n-1)) ) then
             if (m1%ac_p) then
               x=sign(1.d0,random_uniform(mr(i)%random_mag)-0.5d0)
               if (m1%saturate_p) then
                 m(1+id:3+id,n)=x*m1%m
               else
                 m(1+id:3+id,n)=x*k(:,n)
               endif
             else
               m(1+id:3+id,n)=m1%m
               if (.not.m1%saturate_p) then
                 x=sign(1.d0,dot_product(m(1+id:3+id,n),k(:,n)))
                 m(1+id:3+id,n)=x*k(:,n)
               endif
             endif
           endif
         enddo
       endif
       m1=>m1%next
     enddo
   enddo

  end function

  function set_mag(m,ac_p,saturate_p,xmin,xmax,ymin,ymax,reset_p,layer) result (ok_p)
   real(8),intent(in),optional::m(3),xmin,xmax,ymin,ymax
   integer,intent(in),optional::layer
   logical,intent(in),optional::ac_p,saturate_p,reset_p
   logical::ok_p

   type(mag_s),pointer::m1,m2
   real(8)::x

   ok_p=.true.
   if (present(reset_p)) then; if (reset_p) then
     m1=>mag
     do while (associated(m1)); m2=>m1%next; deallocate(m1); m1=>m2; enddo;
     return
   endif; endif

   m1=>mag; nullify(m2)
   do while (associated(m1)); m2=>m1; m1=>m1%next; enddo
   allocate(m1)
   if (associated(m2)) then
     m2%next=>m1
   else
     mag=>m1
   endif

   m1%layer=0; if (present(layer)) m1%layer=max(0,layer)

   m1%m=0.d0
   if (present(m)) then
     x=dot_product(m,m);if (x.ne.0.d0) x=1.d0/sqrt(x); m1%m=m*x
   endif

   m1%ac_p=.false.; if (present(ac_p)) m1%ac_p=ac_p
   m1%saturate_p=.false.; if (present(saturate_p)) m1%saturate_p=saturate_p

   if (maxval(abs(m1%m)).eq.0.d0) then
     m1%ac_p=.true.; m1%saturate_p=.false.
   endif

   m1%xmin=-1.d22; if (present(xmin)) m1%xmin=xmin
   m1%ymin=-1.d22; if (present(ymin)) m1%ymin=ymin
   m1%xmax=1.d22; if (present(xmax)) m1%xmax=xmax
   m1%ymax=1.d22; if (present(ymax)) m1%ymax=ymax

   if (m1%xmin.gt.m1%xmax) then; x=m1%xmin; m1%xmin=m1%xmax; m1%xmax=x; endif
   if (m1%ymin.gt.m1%ymax) then; x=m1%ymin; m1%ymin=m1%ymax; m1%ymax=x; endif
  end function

  function file_mag(ionum,read_p,old_tran,layer,old) result (ok_p)
   use data_m, only:get_data,get_size
   use construct_grain_m, only: get_ms_nom
   use io_m, only:output
   integer,intent(in)::ionum
   logical,intent(in)::read_p
   integer,intent(in),optional::old_tran,layer
   logical,intent(in),optional::old

   real(8),pointer::m(:,:),ms(:),ms_scale(:)
   integer,pointer::id,grain2gridindex(:),grain2grid(:,:),grid2grain(:,:)
   character(200)::os
   logical::ok_p,ascii_p,old_p
   integer::i,i1,i2,nx,ny,n_lay,k,j,n,ll
   real(8)::dx,dy,mg,mm(3),ms_nom
   real(8),allocatable::mms(:,:)
   real(4),allocatable::mms4(:,:,:)

   ok_p=.false.
   old_p=.false.; if (present(old)) old_p=old
   call get_size(nx,ny,n_lay,dx,dy)
   i1=1; i2=n_lay
   if (present(layer)) then
     if (layer.gt.0.and.layer.le.n_lay) then; i1=layer; i2=layer; endif
   endif

   inquire(ionum,form=os); os=adjustl(os); ascii_p=(os(1:2).eq.'FO')

   do k=i1,i2
     call get_data(k,m,m_p=.true.)
     call get_data(k,ms,ms_p=.true.)
     call get_data(k,ms_scale,ms_scale_p=.true.)
     call get_data(k,id,in_p=.true.)
     call get_data(k,grain2gridindex,grain2gridindex_p=.true.)
     call get_data(k,grain2grid,grain2grid_p=.true.)
     call get_data(k,grid2grain,grid2grain_p=.true.)
!if (lbound(grain2gridindex,dim=1).ne.0) then; print *,'file_mag';read *;endif
     ms_nom=get_ms_nom(k)
     if (read_p) then
       if (ascii_p) then
         call output(' trying to read magnetization from ascii file...')
         read(ionum,"(22x,i5,4x,i5)",err=99) i,j 
         if (i.ne.nx.or.j.ne.ny) then
           if (i.ne.nx) call output('  ERROR- nx is incorrect')
           if (j.ne.ny) call output('  ERROR- ny is incorrect')
           go to 99
         endif
         ok_p=.true.
         do j=1,ny
           do i=1,nx
!            read(ionum,*,err=99,end=99)n,ll,ll,mm
             read(ionum,"((3(i5),3(x,e12.5)))",err=99,end=99)n,ll,ll,mm
             if (grid2grain(i,j).gt.0) then
               mg=dot_product(mm,mm);if (mg.ne.0.d0) mm=mm/sqrt(mg)
               m(1+id:3+id,grid2grain(i,j))=mm
               ok_p=(ok_p.and.mg.ne.0.d0)
             endif
           enddo
         enddo
       else
         if (old_p) then
           call output(' trying to read old magnetization from binary file...')
           allocate(mms4(3,nx,ny))
           read(ionum,err=99) mms4
           ok_p=.true.
           do j=1,ny
             do i=1,nx
               n=grid2grain(i,j)
               if (n.gt.size(m,2)) then; ok_p=.false.
               elseif (n.gt.0) then
                 mg=dot_product(dble(mms4(:,i,j)),dble(mms4(:,i,j)))
                 if (mg.gt.0.d0) then; m(1+id:3+id,n)=dble(mms4(:,i,j))/mg
                 else; ok_p=.false.
                 endif
               endif
             enddo
           enddo
           deallocate(mms4)
           if (.not.ok_p) go to 99
         else
           call output(' trying to read magnetization from binary file...')
           read(ionum,err=99) i,j,n,dx
           if (i.ne.nx.or.j.ne.ny.or.n.ne.size(m,2)) then
             if (i.ne.nx) call output('  ERROR- nx is incorrect')
             if (j.ne.ny) call output('  ERROR- ny is incorrect')
             if (n.ne.size(m,2)) call output('  ERROR- number of grains is incorrect')
             go to 99
           endif
           read(ionum,err=99) m(1+id:3+id,1:n)
           ok_p=.true.; if (ms_nom.eq.0.d0) ms_nom=dx
           do i=1,n
             mg=dot_product(m(1+id:3+id,i),m(1+id:3+id,i))
             if (mg.ne.0.d0) then
               m(1+id:3+id,i)=m(1+id:3+id,i)/sqrt(mg)
               ms(n)=sqrt(mg)*ms_nom/dx
             else
               ok_p=.false.
             endif
           enddo
           if (dx.ne.ms_nom) then
             write(os,"('  Ms changed from',1p,e19.12,' to',e19.12)") dx,ms_nom; call output(trim(os))
           endif
         endif
       endif
     else
       if (ascii_p) then
         call output(' trying to write magnetization to ascii file...')
         j=0; if (present(old_tran)) j=old_tran
         write(ionum,"('ZONE T=""',f9.2,'"", I=',i5,', J=',i5,', K=',i3,' F=POINT')",err=99) real(j),nx,ny,1
         do j=1,ny
           do i=1,nx
             n=grid2grain(i,j)
             if (n.gt.0) then
               write(ionum,"(3(i5),1p,7(1x,e12.5))",err=99)i,j,k,m(1+id:3+id,n)*ms(n)*ms_scale(n)
             else
               write(ionum,"(3(i5),1p,7(1x,e12.5))",err=99)i,j,k,0.d0,0.d0,0.d0
             endif
           enddo
         enddo
         ok_p=.true.
       else
         if (old_p) then
           call output(' trying to write old magnetization to binary file...')
!          os=find_file(ionum) 
           allocate(mms4(3,nx,ny));mms4=0.
           do j=1,ny
             do i=1,nx
               n=grid2grain(i,j); if (n.gt.0) mms4(:,i,j)=real(m(1+id:3+id,n))
             enddo
           enddo
           write(ionum,err=99) mms4
           deallocate(mms4)
         else
           call output(' trying to write magnetization to binary file...')
           n=size(m,2); allocate(mms(3,n))
           do i=1,n; mms(:,i)=m(1+id:3+id,i)*ms(i); enddo
           write(ionum,err=99) nx,ny,n,ms_nom  !'1.d0' is for backward compatability
           write(ionum,err=99) mms(1:3,1:n)
           deallocate(mms)
         endif
         ok_p=.true.
       endif
     endif
   enddo
   if (read_p) then; call output(' successful read of magnetization file!')
   else; call output(' successful write of magnetization file!')
   endif
   return
 99 continue
   if (read_p) then; call output(' failure reading magnetization file!')
   else; call output(' failure writing magnetization file!')
   endif
  end function

end module mag_m

