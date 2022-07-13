!  this simple module (graphix) uses PLPlot
!   to do xy plots and contour plots on Linux. 
!
!   this routine was written so that the API is identical to the windows version
!   which uses Digital/Compaq fortran95
!    g.j.parker (c) 2013
!    ver 0.1.0b1

! to use this, the fist line in your program/routine must be 'USE graphix'

module graphix

#ifdef XWIN
! only use plplot if we were explicitly told w/ a -DXWIN in compile line
 use plplot
#endif

 implicit none       !force explicit typing
 private             !hide everything to the outside except for the main subroutine call
 public :: draw,can_draw
  
contains
 

function can_draw() result (ok_p)
  logical:: ok_p
#ifdef XWIN
  ok_p = .true.
#else
  ok_p = .false.
#endif
end function


subroutine draw(title,x,y,z,nx,ny,symbol,num_color,xmin,xmax,ymin,ymax, &
        zmin,zmax,npx,npy,nps,add,legend,color,grain,unsort,scatter,scatter2,scatter3,alpha,beta,box,histogram,clegend,yrange,strip)

!routine to draw either line plots (xy) or contour plots

!  title: unique character string which is used to identify the window
!   x: one dimensional double precision array holding the x-coordinate
!      for xy plots, does not have to be in ascending order, will be sorted (internally) if symbol=.false.
!      for contour plots, must be in ascending order
!   y: one dimensional double precision array holding the corresponding y-coordinate
!      for contour plots, must be in ascending order
!   z: two dimensional double precision array holding contour variable values
!  nx: logical size of array x; default = size of x if xmin > xmax
!  ny: logical size of array y; default = size of y if ymin > ymax
!      for xy plots the number of points = min( nx, ny, n
! symbol:
!      logical to plot symbols (vs. connecting points w/ line segments) for xy plots.
!      ignored for contours. default is .false. (draw straight lines between points)
! num_color:
!      number of contour levels. for < 18, a walk around the 'unit cube' (i.e. rainbow)
! xmin/xmax/ymin/ymax:
!      specify the range of x/y-axis. default is data dependent (auto-scaling) if x(y)min > x(y)max
!      default is xmin=ymin=0.d0, xmax=ymax=-1.d0 (i.e. auto-scaling)
!      for xy plots, only used if new window or redraw (add=.false.)
! zmin/zmax:
!      specify the range of contour levels
! npx: number of pixels along x direction for entire window
!      for xy plots default is 400; for contour plots default at least 400, more to resolve data
! npy: number of pixels along y direction for entire window
!      for xy plots default is 300; for contour plots default at least 400, more to resolve data
! nps: size, in pixels, of symbols on xy plots (default is to scale size to fit graph)
! add: logical. for xy plots, add a curve to the currently displayed graph if .true.
!      otherwise erase window and draw just this curve. default is .false. does not effect 
!      contour plots
! legend:
!      character string holding legend entry for xy plot
! color:
!      color index (0-15: 0-black; 1-dark blue; 2-dark green; 3-dark cyan; 4-dark red; 5-dark violet; 6-dark yellow; 
!        7-dark white; 8-gray; 9-blue; 10-green; 11-cyan; 12-red; 13-violet; 14-yellow; 15-white) for
!      xy plots' line/symbol and legend. default is 15 or the next color if add=.true.
! grain:
!      special case for drawing borders around 'grains'
! histogram:
!      histogram of unbinned data? this is the number of bins


! for xy plots typical call is
!
!  call draw('title string',x=xdata,y=ydata,nx=n,symbol=.true.)
!
! where 'title string' is the name of the xy plot. it must be unique. 
!  xdata is one dimensional double precision array storing the x-coordinate
!    xdata need not be in ascending order
!  ydata is one dimensional double precision array storing the y-coordinate
!  n is an integer which means points 1-n will be plotted. 
!   if n is not specified, it will be equal to the minimum of length of xdata or ydata array
!  symbol=.true. if plot should be scatter plot 
!   otherwise symbol=.false or symbol not specified will connect the points w/ line segments
!    the xdata will be sorted if the latter case
!  to update the plot, use the same 'title string'

! for contour plots typical call is
!
!  call draw('title string',x=xdata,y=ydata,z=zdata,nx=imax,ny=jmax)
!
! where 'title string' is the name of the contour plot. it must be unique. 
!  xdata is one dimensional double precision array storing the x-coordinate- must be in ascending order
!  ydata is one dimensional double precision array storing the y-coordinate- must be in ascending order
!  zdata is two dimensional double precision array storing the contour variable- 
!      i.e. zdata(i,j) is the value of the contour value at coordinates x(i), y(j)
!  imax is an integer which gives the number of x coordinates to use. 
!   if not present then equal to the the minimum size of the first dimension of zdata or xdata
!  jmax is an integer which gives the number of y coordinates to use. 
!   if not present then equal to the the minimum size of the second dimension of zdata or ydata

! if one wants to 'update' the window, draw() must be called w/ the same 'title string'. all the other 
!   data may change (either values and/or dimensions)

! potential problems
!  on contour plots, routine tries to map 1 pixel to the smallest mesh size.
!  if the length divided by the smallest size is large, windows may not let us create a large
!  enough window. if it doesn't, the contour won't plot all our data.

implicit none

!all the input variables- we don't allow them to be changed!
character(*),intent(in)::title
character(*),intent(in),optional::legend
integer,intent(in),optional::nx,ny,num_color,npx,npy,color,nps
double precision,intent(in),optional::x(:),y(:)
double precision,intent(in),optional::z(:,:),xmin,xmax,ymin,ymax,zmin,zmax
logical,intent(in),optional::symbol,add,unsort,clegend,strip
integer,intent(in),optional::grain(:,:),histogram
double precision,optional::scatter2(:,:),scatter(:,:),alpha,beta,box(:,:),scatter3(:)
double precision,optional,intent(out)::yrange

!local variables
logical::first_p,clear_p,unsort_p
integer::nxs,nys,i1,np,nq,i,j,j2,nc,idum
real(8),save,pointer::colorrgb(:,:) !,color(:)
double precision,allocatable::xs(:),ys(:)

double precision::xx,x1,y1,x2,y2,z1,z2,clevels(50)
!character(9) str
character(70) str

!list of windows we have- we'll keep a linked list so we know what is there and what isn't
type::windo
 character(200)::name
 integer::id,nleg,ncur
 type(windo),pointer::nxt
end type
type(windo),save,pointer::window_list
type(windo),pointer::w1_ptr,w2_ptr,w3_ptr

!PGPlot functions and variables
!integer pgopen
real(8):: pgx1,pgx2,pgy2,pgy1
!integer,save::color_dec(15)=(/ 9,10,11,13,12,8,15,14,4,3,5,2,6,7,1 /)
!integer,save::pgc1=0,pgc2=0
integer,save::cscat(10)=(/1,3,9,10,11,2,8,13,5,6/)

!ok, let's do this.
!  first some sanity/bookkeeping
if ((.not.present(x).or..not.present(y)).and..not.present(scatter).and..not.present(strip).and.(.not.present(x).or..not.present(histogram))) return   !we need at least these for anything!

#ifndef XWIN
if (present(yrange)) yrange=0.d0

return  !outta here!

#else
!926 format(F<max(6,i)>.<max(6,i)-i-1>)  !tricky (and non-standard fortran!) way to move decimal point


;
!(present(x).and.present(y)).or.present(scatter).or.(present(x).and.present(histogram))

!ok, first we find if it's a new window being requested or an old one.
!  for new windows, create it in our linked list and guess at the size we should
!  use
!  for old windows, select it and change the size if the user told us to

! we then either plot a contour (we always erase for a contour) or an xy plot
!  for xyplots we are either adding to the plot or erasing and starting over.
!  xy plots are either symbols or straight lines connecting the points. for the
!  latter, we have to sort the input data.


!first see if the window has been created
w3_ptr=>window_list; nullify(w1_ptr); str=adjustl(trim(title))
do while (associated(w3_ptr))
  if (trim(w3_ptr%name) .eq. trim(str)) w1_ptr=>w3_ptr
  w2_ptr=>w3_ptr; w3_ptr=>w3_ptr%nxt
  first_p=associated(w3_ptr)
enddo 

!print *,trim(str),.not.associated(w1_ptr),present(z),present(strip)

if (.not.associated(w1_ptr)) then  !we have to create the window!
 clear_p=.true.                     !clear a new window, always
 allocate(w1_ptr); 
 if (associated(window_list)) then  !store the window on the window list
   w2_ptr%nxt=>w1_ptr
   w1_ptr%id=w2_ptr%id+1
 else
   window_list => w1_ptr
   w1_ptr%id=0
 endif
 w1_ptr%name=adjustl(trim(title)); w1_ptr%nleg=0; w1_ptr%ncur=0
 nullify(w1_ptr%nxt)
!print *,'new window',w1_ptr%id

!need to compute size of window. for contours, this can be involved. 
!for xy plots, it's fixed at 400x300 (though one can change that right here!)
 nxs=600; nys=420; if (present(npx)) nxs=abs(npx); if (present(npy)) nys=abs(npy)
 if (present(scatter)) then
    nxs=800; nys=800
 endif
!nxs=600; nys=450; if (present(npx)) nxs=abs(npx); if (present(npy)) nys=abs(npx)

 if (present(z)) then  !then creating a new contour window
  np=min(size(x),size(z,1));if (present(nx)) np=min(nx,np) !size along x/i
  nq=min(size(y),size(z,2));if (present(ny)) nq=min(ny,nq) !size along y/j
  x1=x(np)-x(1)+.5d0*(x(2)-x(1)+x(np)-x(np-1)); ! x-range
  y1=y(nq)-y(1)+.5d0*(y(2)-y(1)+y(nq)-y(nq-1)); ! y-range
  x2=1.d22; y2=1.d22
  do i=1,np-1; x2=min(x2,abs(x(i+1)-x(i))); enddo  !find finest mesh
  do i=1,nq-1; y2=min(y2,abs(y(i+1)-y(i))); enddo
  nxs=nint((x1)/min(x2,y2)); nys=nint((y1)/min(x2,y2))          !finest mesh gets at least one pixel
! if (max(nxs,nys).lt.600) then         !we can expand the window to 400 pixels...
    if (nxs.ge.nys) then
      nys=int(1000.*y1/x1)
      nxs=1000
    else
      nxs=int(1.1*1000.*x1/y1)
      nys=1000
    endif
! endif
  if (nys.gt.1500) then
    nxs=int(1.1*nxs*1500/nys)
    nys=1500
! elseif (nxs.lt.nys) then
!   nxs=int(1.1*nxs)
  endif
 
  ! allow caller override
  if (present(npx)) then                !preserve the aspect ratio of the viewport
      nxs=abs(npx);
      if(.not.present(npy)) nys=int(nxs*y1/x1)
  endif
  if (present(npy)) then
      nys=abs(npy);
      if(.not.present(npx)) nxs=int(nys*x1/y1)
  endif
! if (nxs.le.650) nxs=nxs/0.9
! if (nys.le.650) nys=nys/0.9

 endif  !end of contour window size
!
!need to set the size of the view surface. assume we're on 72pixel/inch device (Mac is!)
!if (.not.present(strip)) then
call plsstrm(w1_ptr%id)
call plspage(72.d0,72.d0,nxs,nys,0,0)
call plsdev("xwin") !call plsdev("aqt")
call plinit()
idum=plsetopt('plwindow',trim(title))
!call plinit()
if (present(strip)) then
  call pladv(0) 
  call plvsta()
! call plstripc(w1_ptr%id,'bcnst','bcnst',0.d0,1.d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 1, 2 /),(/ 1, 1, 1, 1 /),(/'average','ave-sig','ave+sig','inst   '/),'t','',trim(title))
  xx=0.d0; if (present(x)) xx=x(1)
  if (.not.present(legend)) then
          print *,'not present'
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'0','1','2','3'/),'t','',trim(title))
  elseif (trim(legend).eq.'Hu') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'m','x','y','z'/),'t','',trim(title))
  elseif (trim(legend).eq.'std(M)') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'sig_x','sig_y','sig_z','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'  acos(<mz>)','cos(p)sin(t)','sin(p)sin(t)','      cos(t)'/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>1') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'phi  ','theta','     ','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>3') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'<mz>.<mz>','<mx>     ','<my>     ','<mz>     '/),'t','',trim(title))
  elseif (trim(legend).eq.'te') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'   T1','   T2','T1+T2','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'wr') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'   H1','   H2','H1+H2','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'mz') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'    mz1','    mz2','mz1+mz2','       '/),'t','',trim(title))
  else
          print *,'else ',trim(legend)
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'0','1','2','3'/),'t','',trim(title))
  endif
endif
!else
!call plsdev("xwin") !call plsdev("aqt")
!call plsetopt('plwindow',trim(title))
!call plinit()
!call pladv(0) 
!call plvsta()
!endif

else  !window exists. don't create it again!

!       print *,'found window:',w1_ptr%id

!if (.not.present(strip)) then
 call plsstrm(w1_ptr%id)  !select window
 clear_p=.true.; if (present(add)) clear_p=(.not.add)  !; if (present(z)) clear_p=.true.
 call plgpage(x1,x2,nxs,nys,i,j)
 if (present(npx)) nxs=-abs(npx); if (present(npy)) nys=-abs(npy)
 if (min(nxs,nys).lt.1) then
   nxs=abs(nxs); nys=abs(nys)
   call plspage(72.d0,72.d0,nxs,nys,0,0)
 endif
!endif

endif


 call plspause(.false.)
!ok, window is created/selected. contour or xy plot?

!do a contour plot?
if (present(z)) then  !we're doing a contour


                
!  do we have to create a new color table?
! nc=17; if (present(num_color)) nc=max(2,abs(num_color))
! nc=min(nc,abs(pgc2-pgc1)+1)
  nc=16; if (present(num_color)) nc=max(2,abs(num_color))
  if (associated(colorrgb)) then    !if there's a color table and it doesn't
    if (size(colorrgb,1).ne.nc) then  ! have the correct number of colors, create a new one
      deallocate(colorrgb); nullify(colorrgb)
    endif
  endif

  if (.not.associated(colorrgb)) then   !yes, do something stupid for the colors and get moving
    allocate(colorrgb(nc,4)); colorrgb=0.
    select case (nc)
     case (2:3)
      colorrgb(nc,3)=1.
      pgx1=2.
     case (4)
      colorrgb(nc,1)=1.
      colorrgb(nc,2)=1.
      colorrgb(nc,3)=1.
      pgx1=2.
     case default
      pgx1=5./(nc-2)
      colorrgb(nc,1)=1.
      colorrgb(nc,2)=1.
      colorrgb(nc,3)=1.
    end select; pgx2=0.; 
    do i=1,nc-1
      if (pgx2.lt.2.) colorrgb(i,1)=min(1.,max(0.,2.-pgx2))  !red
      if (pgx2.gt.4.) colorrgb(i,1)=min(1.,max(0.,pgx2-4.))
      if (pgx2.le.3.) colorrgb(i,2)=min(1.,max(0.,pgx2))       !green
      if (pgx2.gt.3.) colorrgb(i,2)=min(1.,max(0.,4.-pgx2))
      if (pgx2.gt.2.) colorrgb(i,3)=min(1.,max(0.,pgx2-2.))  !blue
      pgx2=pgx2+pgx1
      colorrgb(i,4)=(i-1.d0)/(nc-1.d0)
    enddo
    colorrgb(nc,4)=1.d0

    if (present(num_color)) then
      if (num_color.eq.-6) then
        colorrgb=0.d0
        colorrgb(1,1)=1.d0
        colorrgb(2,1)=1.d0
        colorrgb(2,2)=.5d0
        colorrgb(3,:)=1.d0
        colorrgb(4,1)=1.d0
        colorrgb(4,3)=1.d0
        colorrgb(5,2)=1.d0
        colorrgb(6,3)=1.d0
      endif
    endif

    if (nc.gt.17) then   !too many colors, do a sort of banding thing...
     do i=1,nc/2,2;pgy1=colorrgb(i,1); colorrgb(i,1)=colorrgb(nc/2+i,1); colorrgb(nc/2+i,1)=pgy1
                   pgy1=colorrgb(i,2); colorrgb(i,2)=colorrgb(nc/2+i,2); colorrgb(nc/2+i,2)=pgy1
                   pgy1=colorrgb(i,3); colorrgb(i,3)=colorrgb(nc/2+i,3); colorrgb(nc/2+i,3)=pgy1
     enddo
    endif


  endif   !done creating color table
    call plscmap1n(nc)
    if (present(num_color).and.num_color.lt.0) then
      call plscmap1(nint(255*colorrgb(:,1)),nint(255*colorrgb(:,2)),nint(255*colorrgb(:,3)))!,abs(num_color))
    else
      call plscmap1l(.true.,colorrgb(:,4),colorrgb(:,1),colorrgb(:,2),colorrgb(:,3))
    endif

         if (present(zmin)) then; z1=zmin; else; z1=minval(z); endif
         if (present(zmax)) then; z2=zmax; else; z2=maxval(z); endif
         if (z1.gt.z2) then; x1=z1; z1=z2; z2=x1; endif
         if (z1.eq.z2) z2=max(0.1d0,abs(z1)*1.01d0)*sign(1.d0,z1)
         if (z1.gt.z2) then; x1=z1; z1=z2; z2=x1; endif
         if (z2-z1.lt.z2*1.d-4) then; z2=z1*1.05d0; z1=z1*.95d0; endif
         do i=1,nc+1
           clevels(i)=z1+(z2-z1)*(i-1)/nc
         enddo
         x1=0.5d0*(x(2)-x(1)); y1=0.5d0*(y(2)-y(1))
         x2=maxval(x)+x1;x1=minval(x)-x1
         y2=maxval(y)+y1;y1=minval(y)-y1
!        x1=minval(x);x2=maxval(x)
!        y1=minval(y);y2=maxval(y)

!        call plspal0('cmap0_white_on_black.pal')
!        call plspal1('cmap1_blue_red.pal',1)
!        call plspal1('cmap1_highfreq.pal',1)
!        call plscmap0n(3)

         call plcol0(0)
         call pladv(0)
!        call plvpor(0.15d0,0.95d0,0.2d0,0.9d0)
         call plvpor(0.2d0,0.95d0,0.2d0,0.9d0)
         call plwind(x1,x2,y1,y2)
         call plpsty(0)
!        call plshades(z,' ',x1,x2,y1,y2,clevels(1:nc+1),1.d0,0,1.d0,x,y)
         pgy1=y1;
         np=min(size(x),size(z,1));if (present(nx)) np=min(nx,np) !size along x/i
         nq=min(size(y),size(z,2));if (present(ny)) nq=min(ny,nq) !size along y/j
         do j=1,size(z,2); pgy2=0.5d0*(y(j)+y(min(nq,j+1))); if (j.eq.nq) pgy2=y2; pgx1=x1
           do i=1,size(z,1); pgx2=0.5d0*(x(i)+x(min(np,i+1))); if (i.eq.np) pgx2=x2
             if (z(i,j).ge.z1.and.z(i,j).le.z2) then
               call plcol1(min(1.d0,max(0.d0,(z(i,j)-z1)/(z2-z1))))
               call plfill((/pgx1,pgx2,pgx2,pgx1/),(/pgy1,pgy1,pgy2,pgy2/))
             endif
             pgx1=pgx2
           enddo
           pgy1=pgy2
         enddo
         call plcol0(0)
         call pl_setcontlabelformat(4,3)
         if (present(clegend)) then 
            call pl_setcontlabelparam(0.006, 0.3, 0.1, 1)
            call plcont(z, 1, size(z,1), 1, size(z,2), clevels(1:nc+1), x, y)
         endif
         call plcol0(15)
!        call plbox('bcnst', 0.0_plflt, 0, 'bcnstv', 0.0_plflt, 0)
         call plbox('BCINTS', 0.0_plflt, 0, 'BCINTSV', 0.0_plflt, 0)
         call pllab('x (nm)', 'y (nm)',  trim(title))

       !do legend (newer plplot can do this, but screw it...)
         call plvpor(0.1d0,0.95d0,0.1d0,0.15d0)
         call plwind(z1,z2,0.d0,0.05d0)
!        call plshades(reshape((/z1,z2,z1,z2/),(/2,2/)),' ',z1,z2,0.d0,0.01d0,clevels(1:nc+1),1,0,1,(/z1,z2/),(/0.d0,0.01d0/))
!        call plshades(reshape((/z1,z2,z1,z2/),(/2,2/)),' ',z1,z2,0.d0,0.01d0,clevels(1:nc+1),1.d0,0,1.d0,(/z1,z2/),(/0.d0,0.01d0/))
         call plshades(reshape((/z1,z2,z1,z2/),(/2,2/)),z1,z2,0.0_plflt,0.01_plflt,clevels(1:nc+1),1.d0,0,1.d0,.true.) !,(/z1,z2/),(/0.0,0.01/))
         call plcol0(15)
         call plbox('BINTS', 0.0_plflt, 0, ' ', 0.0_plflt, 0)
!        print *,x1,x2

 !do we want to draw a border around the grains?
!        if (present(grain).and.max(nxs,nys).lt.900) then
         if (present(grain)) then
!          call plvpor(0.15d0,0.95d0,0.2d0,0.9d0)
           call plvpor(0.2d0,0.95d0,0.2d0,0.9d0)
           call plwind(x1,x2,y1,y2)
           call plcol0(0)      !black borders
           y1=(y(2)-y(1))*.5d0; x1=(x(2)-x(1))*.5d0
           do j=1,size(z,2)
            do i=1,size(z,1)-1
              if (grain(i,j).ne.grain(i+1,j)) call pljoin(x(i)+x1,y(j)-y1,x(i)+x1,y(j)+y1)
!            if (grain(i,j).ne.grain(i+1,j)) print *,i,j,'x'
            enddo
           enddo
           do j=1,size(z,2)-1
            do i=1,size(z,1)
              if (grain(i,j).ne.grain(i,j+1)) call pljoin(x(i)-x1,y(j)+y1,x(i)+x1,y(j)+y1)
!             if (grain(i,j).ne.grain(i,j+1)) print *,i,j,'y'
            enddo
           enddo
           call plcol0(15)      !black borders
         endif
         if (present(box)) then
           call plcol0(1)
           call plline(box(:,1),box(:,2))
           call plcol0(15)      !black borders
         endif

! call pgebuf()   !buffer PGPlot commands- faster drawing
! call pgupdt()     !force an update on the display

elseif (present(scatter)) then !doing a 3d scatter plot

  !draw sphere
! call pgsci(1)   !color of sphere

  xx = acos(-1.d0)/180.d0
  z1=180.d0; z2=45.d0
  z1=35.264d0; z2=45.d0
  z1=35.264d0; z2=-135.d0
  if (present(alpha)) z1=alpha
  if (present(beta)) z2=beta
  z1 = z1 * xx; z2 = z2 * xx

  x1=cos(z1); x2=cos(z2); y1=sin(z1); y2=sin(z2)

  !plot.f90: never calls w/ add
  ! if add present and false; xor dots
  ! if add present and true; add dots
  ! if add not present; redraw window and dots (assume they are xor)

  clear_p=.false.; if (present(add)) clear_p=.not.add  !xor dots
  if (.not.present(add)) then !draw grid
!     clear_p=.true.; if (present(add)) clear_p=(.not.add)  !; if (present(z)) clear_p=.true.
!
! if (.not.present(add).or..not.add.or.clear_p) then
!print *,'drawing grid',clear_p,.not.present(add) !,present(add).and..not.add
         call pladv(0)
         call plvpor(0.025d0,0.975d0,0.0d0,0.95d0)
         call plwind(-1.05d0,1.05d0,-1.05d0,1.05d0)
         call plcol0(15)
         call pllab(' ', ' ',  trim(title))
         call plcol0(7)
         do i=0,359,30/2
           z2=1.d0; if (mod(i,90).eq.0) z2=0.25d0
           z1=0.d0; do while (z1.le.180.d0)
!          do z1=0.d0,180.d0,z2
             call plpoin((/x2*cos(i*xx)*sin(z1*xx)-y2*sin(i*xx)*sin(z1*xx)/),(/x1*cos(z1*xx)+y1*sin(z1*xx)*(sin(i*xx)*x2+cos(i*xx)*y2)/),-1)  !plot the symbols
             z1=z1+z2
           enddo
         enddo
         do j=15,180-15,15
           z2=max(1,nint(1/abs(sin(j*xx)))); if (mod(j,90).eq.0) z2=z2*0.25d0
           z1=0.d0; do while (z1.lt.360.d0)
!          do z1=0.d0, 359.99999d0, z2
             call plpoin((/x2*cos(z1*xx)*sin(j*xx)-y2*sin(z1*xx)*sin(j*xx)/),(/x1*cos(j*xx)+y1*sin(j*xx)*(sin(z1*xx)*x2+cos(z1*xx)*y2)/),-1)  !plot the symbols
             z1=z1+z2
           enddo
         enddo
         call plcol0(1)   !color of plot x-axis
         call pljoin(0.d0,0.d0,x2,y1*y2)
         call pllsty(2)
         call pljoin(0.d0,0.d0,-x2,-y1*y2)
         call plcol0(3)   !color of plot y-axis
         call pllsty(1)
         call pljoin(0.d0,0.d0,-y2,y1*x2)
         call pllsty(2)
         call pljoin(0.d0,0.d0,y2,-y1*x2)
         call plcol0(9)   !color of plot  z-axis
         call pllsty(1)
         call pljoin(0.d0,0.d0,0.d0,x1)
         call pllsty(2)
         call pljoin(0.d0,0.d0,0.d0,-x1)
         call pllsty(1)
       call plxormod(.true.,clear_p);!print *,'xormod?',clear_p
       clear_p=.false.
  endif
! call plgchr(z1,z2)  !size in mm, scale
! print *,z1,z2
! print *,z1,z2
  call plcol0(15)  !color of spin
  if (present(scatter3)) then   !magnitude of spins can be other than 1
    xx=0.01d0
    do i=1,size(scatter,2)
!     call plpoin((/scatter3(i)*(x2*scatter(1,i)-y2*scatter(2,i))/),(/scatter3(i)*(x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2))/),-1)  !plot the symbols
      call plarc(scatter3(i)*(x2*scatter(1,i)-y2*scatter(2,i)),scatter3(i)*(x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2)),xx,xx,0.d0,360.d0,0.d0,.true.)
    enddo
    if (present(scatter2)) then
      if (maxval(abs(scatter2)).gt.0.d0) then
        call plcol0(5)  !color of tail
        do i=1,size(scatter,2)
          call pljoin(scatter3(i)*(x2*scatter(1,i)-y2*scatter(2,i)),scatter3(i)*(x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2)), &
                scatter3(i)*(x2*scatter2(1,i)-y2*scatter2(2,i)),scatter3(i)*(x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2)))
        enddo
      endif
    endif
  else
    if (size(scatter,2).gt.10) then
     xx=0.01d0

     clear_p=(clear_p.and.present(scatter2))
     if (clear_p) call plxormod(.true.,first_p)
     do i=1,size(scatter,2)
       if (clear_p) call plarc(x2*scatter2(1,i)-y2*scatter2(2,i),x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
       call plarc(x2*scatter(1,i)-y2*scatter(2,i),x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
       if (present(scatter2)) then
         call plcol0(14)  !color of tail
         call pljoin(x2*scatter(1,i)-y2*scatter(2,i),x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2), &
                x2*scatter2(1,i)-y2*scatter2(2,i),x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2))
         call plcol0(15)  !color of spin
       endif
     enddo

!     if (clear_p.and.present(scatter2)) then
!       call plxormod(.true.,clear_p)
!       do i=1,size(scatter2,2)
!         call plarc(x2*scatter2(1,i)-y2*scatter2(2,i),x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
!       enddo
!     endif
!     do i=1,size(scatter,2)
!!      call plpoin((/x2*scatter(1,i)-y2*scatter(2,i)/),(/x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2)/),-1)  !plot the symbols
!       call plarc(x2*scatter(1,i)-y2*scatter(2,i),x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
!     enddo
!    if (present(scatter2).and..not.clear_p) then
!     if (present(scatter2)) then
!!    if (.false.) then
!       if (maxval(abs(scatter2)).gt.0.d0) then
!         call plcol0(14)  !color of tail
!         do i=1,size(scatter,2)
!          call pljoin(x2*scatter(1,i)-y2*scatter(2,i),x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2), &
!                x2*scatter2(1,i)-y2*scatter2(2,i),x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2))
!         enddo
!       endif
!     endif
     call plxormod(.false.,clear_p)
    else
      if (clear_p) then
        call plxormod(.true.,clear_p)
        xx=0.035d0
        do i=size(scatter,2),1,-1
          call plcol0(cscat(i))
          call plarc(x2*scatter2(1,i)-y2*scatter2(2,i),x1*scatter2(3,i)+y1*(scatter2(2,i)*x2+scatter2(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
          xx=0.05d0
        enddo
!       call plxormod(.false.,clear_p)
      endif
      xx=0.05d0
      do i=1,size(scatter,2)
        call plcol0(cscat(i))
        if (i.eq.size(scatter,2)) xx=0.035d0
        call plarc(x2*scatter(1,i)-y2*scatter(2,i),x1*scatter(3,i)+y1*(scatter(2,i)*x2+scatter(1,i)*y2),xx,xx,0.d0,360.d0,0.d0,.true.)
      enddo
      call plxormod(.false.,clear_p)
    endif
  endif

elseif (present(histogram)) then
    np=size(x);if (present(nx)) np=min(nx,np)
    i=10;if (histogram.ne.0) i=histogram
    x1=0.d0; x2=-1.d0;
    if (present(xmin)) x1=xmin; if (present(xmax)) x2=xmax
    if (x1.gt.x2) then; x1=minval(x(1:np)); x2=maxval(x(1:np)); endif  !get x and y range of values. 
    if (x2-x1.lt.x2*1.d-4) then; x2=x1*1.05d0; x1=x1*.95d0; endif
    if (x1.eq.x2) then; x2=x1*1.d0; x1=x1*.9d0; else; z1=(x2-x1)/i;x2=x2+z1;x1=x1-z1;endif
    call plcol0(15)
    j=pl_hist_ignore_outliers+pl_hist_noexpand
    if (clear_p) then; w1_ptr%nleg=0; w1_ptr%ncur=14; endif
    if (present(add)) then; if (add) j=j+PL_HIST_NOSCALING; endif
    i1=mod(w1_ptr%ncur,15)+1
    if (present(color)) i1=min(15,max(1,abs(color))); 
    call plcol0(i1)
    call plhist(x(1:np),x1,x2,i+2,j)
    call pllab(' ',' ',trim(title))
    w1_ptr%ncur=mod(w1_ptr%ncur+1,15)  !number of curves- for colors if user adds plots

 !done plotting, do we need a legend?
    if (present(legend)) then;    !is there a legend?
     call plmtex('lv',-1.d0,.95d0-0.05d0*w1_ptr%nleg,0.d0,trim(legend))
     w1_ptr%nleg = w1_ptr%nleg+1
    endif
    if (present(yrange)) call plgvpw(x1,x2,y1,yrange)
elseif (present(strip)) then
    if (present(add)) then
      if (add) then
!print*,min(size(y),size(x))
        do i=1,min(size(y),size(x))
          call plstripa(w1_ptr%id, i-1, x(i), y(i))
        enddo
      else
        call plstripd(w1_ptr%id)
        call pladv(0) 
        call plvsta()
!print "(i0,' R**',a,'**R')",w1_ptr%id,trim(title)
!       call plstripc(w1_ptr%id,'bcnst','bcnst',0.d0,1.d-14,0.3d0,-1d-22,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 1, 2 /),(/ 1, 1, 1, 1 /),(/'average','ave-sig','ave+sig','inst   '/),'t','',trim(title))
!        call plstripc(w1_ptr%id,'bcnst','bcnst',0.d0,1.d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 1, 2 /),(/ 1, 1, 1, 1 /),(/'average','ave-sig','ave+sig','inst   '/),'t','',trim(title))
!       call plstripc(w1_ptr%id,'bcnst','bcnst',0.d0,1.d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'0','1','2','3'/),'t','',trim(title))
  xx=0.d0; if (present(x)) xx=x(1)
  if (.not.present(legend)) then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'0','1','2','3'/),'t','',trim(title))
  elseif (trim(legend).eq.'Hu') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'m','x','y','z'/),'t','',trim(title))
  elseif (trim(legend).eq.'std(M)') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'sig_x','sig_y','sig_z','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'  acos(<mz>)','cos(p)sin(t)','sin(p)sin(t)','      cos(t)'/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>1') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'phi  ','theta','     ','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'<M>3') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'<mz>.<mz>','<mx>     ','<my>     ','<mz>     '/),'t','',trim(title))
  elseif (trim(legend).eq.'te') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'   T1','   T2','T1+T2','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'wr') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'   H1','   H2','H1+H2','     '/),'t','',trim(title))
  elseif (trim(legend).eq.'mz') then
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'    mz1','    mz2','mz1+mz2','       '/),'t','',trim(title))
  else
    call plstripc(w1_ptr%id,'bcnst','bcnst',xx,xx+1d-3,0.3d0,0.d0,1.d-22,0.d0,0.25d0,.true.,.true.,15,1,(/ 15, 1, 2, 3 /),(/ 1, 1, 1, 1 /),(/'0','1','2','3'/),'t','',trim(title))
  endif
!print "(i0,' R**',a,'**R')",w1_ptr%id,trim(title)
      endif
    endif
else !we're doing xy plot

 !get number of points to plot = min (size of x, size of y, user specified)
  np=min(size(y),size(x));if (present(nx)) np=min(nx,np); if (present(ny)) np=min(np,ny);  

  if (clear_p) then  !clear device
!         print *,'clear!'
    if (clear_p) call plclear()
    x1=0.d0; x2=-1.d0; y1=0.d0; y2=-1.d0
    if (present(xmin)) x1=xmin; if (present(xmax)) x2=xmax
    if (present(ymin)) y1=ymin; if (present(ymax)) y2=ymax
    if (x1.gt.x2) then; x1=minval(x(1:np)); x2=maxval(x(1:np)); endif  !get x and y range of values. 
    if (y1.gt.y2) then; y1=minval(y(1:np)); y2=maxval(y(1:np)); endif  !get x and y range of values. 
    if (x1.eq.x2) then; x2=x1*1.05d0; x1=x1*.95d0; endif
    if (y1.eq.y2) then; y2=y1*1.05d0; y1=y1*.95d0; endif
    if (x2-x1.lt.x2*1.d-4) then; x2=x1*1.05d0; x1=x1*.95d0; endif
    if (y2-y1.lt.y2*1.d-4) then; y2=y1*1.05d0; y1=y1*.95d0; endif
    if (y2.eq.y1) then; y2=0.05d0; y1=-0.05d0; endif
!         print *,'clear!'
!   call pladv(0)
    call plcol0(15)
!         print *,'plenv0',x1,x2,y1,y2,(x1.eq.x2),(y1.eq.y2)
!   call plinit()
    call plenv0(x1,x2,y1,y2,0,2)
!         print *,'pllab'
!   call plenv(0.d0,1.d0,0.d0,1.d0,0,2)
    call pllab(' ',' ',trim(title))
!         print *,'done'
    w1_ptr%nleg=0; w1_ptr%ncur=14
  endif

 !select the color (first white, then we cycle through them, unless user wants something else)
  i=mod(w1_ptr%ncur,15)+1
  if (present(color)) i=min(15,max(1,abs(color))); 
  call plcol0(i)


 !finally, let's do the xy plot!
  first_p=.false.
  if (present(symbol)) first_p=symbol

 ! do line segments
  unsort_p=.false.; if (present(unsort)) unsort_p=unsort
  if (.not.first_p) then                   !do w/ line segments, we have to sort!
   if (unsort_p) then
     call plline(x(1:np),y(1:np))
   else
     allocate(xs(np),ys(np))                  !temp storage for sorting
     xs(1)=x(1); ys(1)=y(1)                  !simple insertion sort
     do i=2,np; j=1; do j=1,i-1; if (x(i).lt.xs(j)) exit; enddo
       do j2=i-1,j,-1;  xs(j2+1)=xs(j2); ys(j2+1)=ys(j2); enddo; 
       xs(j)=x(i); ys(j)=y(i); enddo
     call plline(xs(1:np),ys(1:np))
     deallocate(xs,ys)                       !free sort array memory
   endif

  else                     !just points

 ! do symbols
   call plgvpd(pgx1,pgx2,pgy1,pgy2)
   xx=min(abs(pgx1-pgx2),abs(pgy1-pgy2))/np   !size of symbols in inches
   xx=1.d0
!  if (present(nps)) xx=abs(nps/72.)
   if (present(nps)) xx=abs(nps)
   call plssym(0.d0,xx)
!  call plstring(x(1:np),y(1:np),'*')
   call plpoin(x(1:np),y(1:np),1)
   call plssym(0.d0,1.d0)

  endif
  w1_ptr%ncur=mod(w1_ptr%ncur+1,15)  !number of curves- for colors if user adds plots

 !done plotting, do we need a legend?
  if (present(legend)) then;    !is there a legend?
   call plmtex('lv',-1.d0,.95d0-0.05d0*w1_ptr%nleg,0.d0,trim(legend))
   w1_ptr%nleg = w1_ptr%nleg+1
  endif

endif  !done w/ contour or xy plot
call plflush()

return  !outta here!

!926 format(F<max(6,i)>.<max(6,i)-i-1>)  !tricky (and non-standard fortran!) way to move decimal point

#endif

end subroutine
end module


