!****************************************************************
! MODULE mod_contour
!****************************************************************


MODULE mod_contour
use f90_kind
implicit none

type contour_level_curve_pos_type
  real(rkind)      :: x_indx                        ! real valued index of x-position
  real(rkind)      :: y_indx                        ! real valued index of y-position
  real(rkind)      :: r                             ! x-dimension (major radius)
  real(rkind)      :: z                             ! y-dimension (z)
end type contour_level_curve_pos_type

type contour_level_curve_type
  integer(ikind)   :: N_pos
  type(contour_level_curve_pos_type), dimension(:), allocatable :: pos
end type contour_level_curve_type

type contour_level_type
  integer(ikind)   :: N_curves                      ! number of curves found
  real(rkind)      :: amp                           ! level height
  type(contour_level_curve_type), dimension(:), allocatable :: curve
end type contour_level_type

type contour_type
  integer(ikind)   :: N_levels
  type(contour_level_type), dimension(:), allocatable :: level
end type contour_type

!private :: 
public  :: contourgq,                   &
           contourgq_non_eqd,           &
           contouring,                  &
           contour_indx2rz,             &
           map_contour_indx_to_rz

CONTAINS


subroutine contourgq(ier, a, mx, nx, nz, lx, idim, &
                        ncsmax, h, xpl, zpl, ncs, npt)
! file=/afs/ipp/u/pem/libpem/Plotsource/contourgq.f

!                                  C O N T O U R G Q
! input:
!   a     - real matrix
!   mx    - leading dimension of matrix a
!   nx    - integer(ikind) ::, number of grid points in r-direction, nx > 3
!   nz    - integer(ikind) ::, number of grid points in z-direction, nz > 3
!   lx    - aux. integer(ikind) :: vector of length  idim
!   idim  - integer(ikind) ::, dimension of vectors  xpl, zpl, lx
!           (suggested: idim = 4*(nx+nz)
!   ncsmax- integer(ikind) ::, dimension of output vector npt
!   h     - real value for which contours in a are searched
! output:
!   ier   - integer(ikind) ::, error code. ier = 0  okay
!           ier = 1 : nx < 4  oder nz < 4
!           ier = 2 : idim insufficient
!           ier = 3 : ncsmax insufficient
!   xpl, zpl
!         - r-values of plot points, dimension idim
!           z-values of plot points, dimension idim
!   ncs   - integer(ikind) :: number of curves found
!   npt   - integer(ikind) :: vector, number of points for curves found
!
!                                                    revised aug.1995
      real(rkind),   intent(in) :: a(*), h
      integer(ikind),intent(in) :: mx, nx, nz, idim, ncsmax
      real(rkind),   intent(out):: xpl(*), zpl(*)
      integer(ikind),intent(out):: lx(*), npt(*)
      integer(ikind),intent(out):: ncs, ier

      integer(ikind) :: ja(6), ix, ixa, mxr, mxrt, mpl, mpllast, jx, jm, in, &
              jabs, kx, ikx, jfor, is, ie, k, i, j, instart, jnb, npl
      real(rkind) :: adn, a1, b1, dp, p1, p2

      if(nx < 4  .or. nz < 4) then
         ier=1
         return
      end if

      ier=0
      ix=0
      ixa=1
      mxr=mx*(nz-1)
      mxrt=mxr+mx
      ncs=0
      mpl=0
      mpllast=0

! find intersects contour line - grid lines: are stored in table lx
! 1st: look to the grid boundaries

      do jx=2,nx-1
         if((ix+(mxr-jx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxr,mx
            j=jm+mx
            if( ((a(j)-h)*(a(jm)-h))  <= 0.0d0) then
               ix=ix+1
!              if(ix > idim) goto 6000
               lx(ix)=-j
            end if
         end do
      end do

      do jx=2,nx
         if((ix+(mxr-jx-mx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx+mx,mxr,mx
            if( ((a(jm)-h)*(a(jm-1)-h))  <= 0.0d0) then
               ix=ix+1
!              if(ix > idim) goto 6000
               lx(ix)=jm
            end if
         end do
      end do

      do jx=2,nx
         if((ix+(mxrt-jx)/mxr+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxrt,mxr
            if( ((a(jm)-h)*(a(jm-1)-h))  <= 0.0d0) then
               ix=ix+1
!              if(ix > idim) goto 6000
               lx(ix)=jm
            end if
         end do
      end do

      do jx=1,nx,nx-1
         if((ix+(mxr-jx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxr,mx
            j=jm+mx
            if( ((a(j)-h)*(a(jm)-h))  <= 0.0d0) then
               ix=ix+1
!              if(ix > idim) goto 6000
               lx(ix)=-j
            end if
         end do
      end do

      if(ix < 2) return

! start a new curve

  do ! 108 continue
      if(ncs == ncsmax) then
         ier=3
         return
      end if
      ncs=ncs+1
      in=ix
      jx=lx(in)
      jfor=0
      instart=in
      npl=0

! continue curve

   30 continue
      jabs=iabs(jx)
      kx=(jabs-1)/mx+1
      ikx= jabs-mx*(kx-1)
      if(jx < 0) jnb=jabs-mx
      if(jx > 0) jnb=jabs-1

      if(mpl == idim)then
        ier=2
        return
      endif

      mpl=mpl+1
      npl=npl+1

! interpolation for point found
! linear interpolation at boundaries, else quadr. interpolation.
! note: for method used, check only for upper boundaries required.
      if(jx < 0) then
         if((jabs+mx) <= mxrt) then
            a1=(a(jabs+mx)-2.0d0*a(jabs)+a(jabs-mx))*2.0d0
            if(a1 /= 0.0d0) then
               b1=a(jabs+mx)-a(jabs-mx)
               dp=sqrt(b1**2 - 4.0d0*a1*(a(jabs)-h))
               p1= (dp-b1)/a1
               p2=-(dp+b1)/a1
               if(abs(p1+0.5d0) > abs(p2+0.5d0)) p1=p2
               zpl(mpl)=dble(kx-1)+p1
            else
               adn=a(jabs)-a(jabs-mx)
               if(adn /= 0.0d0) p1=-(a(jabs)-h)/adn
               if(adn == 0.0d0) p1=-0.5d0
               zpl(mpl)=dble(kx-1)+p1
            end if
         else
            adn=a(jabs)-a(jabs-mx)
            if(adn /= 0.0d0) p1=-(a(jabs)-h)/adn
            if(adn == 0.0d0) p1=-0.5d0
            zpl(mpl)=dble(kx-1)+p1
         end if
         xpl(mpl)=dble(ikx-1)
      else
         if(ikx < nx) then
            a1=(a(jabs+1)-2.0d0*a(jabs)+a(jabs-1))*2.0d0
            if(a1 /= 0.0d0) then
               b1=a(jabs+1)-a(jabs-1)
               dp=sqrt(b1**2 - 4.0d0*a1*(a(jabs)-h))
               p1= (dp-b1)/a1
               p2=-(dp+b1)/a1
               if(abs(p1+0.5d0) > abs(p2+0.5d0)) p1=p2
               xpl(mpl)=dble(ikx-1)+p1
            else
               adn=a(jabs)-a(jabs-1)
               if(adn /= 0.0d0) p1=-(a(jabs)-h)/adn
               if(adn == 0.0d0) p1=-0.5d0
               xpl(mpl)=dble(ikx-1)+p1
            end if
         else
            adn=a(jabs)-a(jabs-1)
            if(adn /= 0.0d0) p1=-(a(jabs)-h)/adn
            if(adn == 0.0d0) p1=-0.5d0
            xpl(mpl)=dble(ikx-1)+p1
         end if
         zpl(mpl)=dble(kx-1)
      end if

! look for next point

      is=1
      ie=6
      if(jx <= 0) then
         ja(1)=iabs(jx)
         ja(2)=-(iabs(jx)-1)
         ja(3)=ja(1)-mx

         ja(4)=-(iabs(jx)+1)
         ja(5)=iabs(jx)+1
         ja(6)=ja(5)-mx
         if(ikx == 1) is=4
         if(ikx == nx) ie=3
      else
         ja(1)=jx+mx
         ja(2)=-(jx+mx)
         ja(3)=-(iabs(ja(2))-1)
         ja(4)=-jx
         ja(5)=jx-mx
         ja(6)=-(jx-1)
         if(kx == 1) ie=3
         if(kx == nz) is=4
      end if
! make this to an internel function
      do  k=is,ie
         if(ja(k) /= jfor) then
            do i=ixa,ix
               if(lx(i) == ja(k)) then
                  if(jfor /= 0) lx(in)=0
                  jfor=jx
                  jx=lx(i)
                  in = i
                  go to 30
               end if
            end do
         end if
      end do
! hier exit
! hier 30 enddo
! curve found

      if(npl > 1) then
         npt(ncs)=npl
         mpllast=mpl
      else
         ncs=ncs-1
         mpl=mpllast
      end if

! remove used intersects from table lx and try to find another
! curve line for contour value

      if(in == instart) lx(instart)=0
      if(in > 0) lx(in)=0

!    6 continue
!      if(lx(ix) == 0) then
!         ix=ix-1
!         if(ix > ixa) goto 6
!      end if
     do
       if(lx(ix) /= 0) exit
       ix=ix-1
       if(ix <= ixa) exit
     enddo

!    7 continue
!      if(lx(ixa) == 0) then
!         ixa=ixa+1
!         if(ixa < ix) goto 7
!      end if
     do
       if(lx(ixa) /= 0) exit  
       ixa=ixa+1
       if(ixa >= ix) exit  
     enddo

!     if((ix-ixa+1) >= 2) goto 108
     if((ix-ixa+1) < 2) exit
   enddo ! 108


end subroutine contourgq

!*****************************************************************

subroutine contourgq_non_eqd(ier, a, xcoord, zcoord, mx, nx, nz, lx, idim, &
                        ncsmax, h, xpl, zpl, ncs, npt)
! file=/afs/ipp/u/pem/libpem/Plotsource/contourgq.f

!                                  C O N T O U R G Q
! input:
!   a     - real matrix -  read in as a vector
!   mx    - leading dimension of matrix a
!   nx    - integer(KIND=4) ::, number of grid points in r-direction, nx > 3
!   nz    - integer(KIND=4) ::, number of grid points in z-direction, nz > 3
!   lx    - aux. integer(KIND=4) :: vector of length  idim - stores the coordinates of the penetration points to be interpolated
!   idim  - integer(KIND=4) ::, dimension of vectors  xpl, zpl, lx
!           (suggested: idim = 4*(nx+nz)
!   ncsmax- integer(KIND=4) ::, dimension of output vector npt
!   h     - real value for which contours in a are searched
! output:
!   ier   - integer(KIND=4) ::, error code. ier = 0  okay
!           ier = 1 : nx < 4  oder nz < 4
!           ier = 2 : idim insufficient
!           ier = 3 : ncsmax insufficient
!   xpl, zpl
!         - r-values of plot points, dimension idim
!           z-values of plot points, dimension idim
!   ncs   - integer(KIND=4) ::, number of curves found
!   npt   - integer(KIND=4) :: vector, number of points for curves found
!
!                                                    revised aug.1995
! Minor revisions 20.12.12: S. Denk
      real(rkind), dimension(*),  intent(in) :: a
      real(rkind), dimension(:), intent(in) :: xcoord, zcoord
      real (rkind), intent(in)              :: h
      integer(ikind),intent(in) :: mx, nx, nz, idim, ncsmax
      real(rkind), dimension(:),  intent(out):: xpl, zpl
      integer(ikind), dimension(:), intent(out):: lx, npt
      integer(ikind),intent(out):: ncs, ier

      integer(ikind) :: pos_neighbour(6), & ! look up table to find the next point on the curve
                         pnt_cnt, & ! penetration points found so far, that are not yet interpolated -
                         !this is counted down when we start interpolating
                         lower_boundary, & ! index for lx, all points below this index are already interpolated
                         mxr, mxrt, pnt_int, & !counts the points interpolated so fa
                         pnt_int_last, & ! saves the pnt_int from the last curve - used to eliminate curves with just 1 point
                         jx, & ! used as coordinates in a - negative coordinates singal horizontal sign changes in a
                         jm, in, &
                         jabs, kx, ikx, jfor, is, ie, k, i, j, instart, jnb, npl, startcoord
      real(rkind) :: adn, &
                      a1, a2, a3, &! Newton interpolation coefficients
                      dp, p1, p2, x1,x2,x3
      logical :: found_neighbour

      if(nx < 4  .or. nz < 4) then
         ier=1
         return
      end if
      ier=0
      pnt_cnt=0
      lower_boundary=1
      mxr=mx*(nz-1)
      mxrt=mxr+mx
      ncs=0
      pnt_int=0
      pnt_int_last=0

! find intersects contour line - grid lines: are stored in table lx
! The order in which we store the points is important, since it is important
! to start all curves at the edges. Since we step backwards through the points found,
! the edges of the matrix are added last.
! First
! Saving points that are right of a horizontal sign change of a - h
! Disregarding the left and the right column
      do jx=2,nx-1 ! here jx is the x coordinate
         if((pnt_cnt+(mxr-jx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxr,mx ! jm is the x + z coordinate
            j=jm+mx
            if( ((a(j)-h)*(a(jm)-h))  <= 0.0d0) then
              !  if( (a(j)-h)*(a(jm)-h) == 0) print*, pnt_cnt + 1
               pnt_cnt=pnt_cnt+1
!              if(pnt_cnt > idim) goto 6000
               lx(pnt_cnt)=-j
            end if
         end do
      end do
! Second
! Saving points that are below a vertical sign change of a - h
! Disregarding the bottom and the top row
      do jx=2,nx
         if((pnt_cnt+(mxr-jx-mx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx+mx,mxr,mx
            if( ((a(jm)-h)*(a(jm-1)-h))  <= 0.0d0) then
              !  if(  ((a(jm)-h)*(a(jm-1)-h)) == 0) print*, pnt_cnt + 1
               pnt_cnt=pnt_cnt+1
!              if(pnt_cnt > idim) goto 6000
               lx(pnt_cnt)=jm
            end if
         end do
      end do
! Third
! Now only tread the right and the left column (again horizontal)
      do jx=2,nx
         if((pnt_cnt+(mxrt-jx)/mxr+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxrt,mxr
            if( ((a(jm)-h)*(a(jm-1)-h))  <= 0.0d0) then
               pnt_cnt=pnt_cnt+1
!              if(pnt_cnt > idim) goto 6000
               lx(pnt_cnt)=jm
            end if
         end do
      end do
! Four
! Now only tread the top and the bottom column (finally vertical)
      do jx=1,nx,nx-1
         if((pnt_cnt+(mxr-jx)/mx+1) > idim)then
           ier=2
           return
         endif
         do jm=jx,mxr,mx
            j=jm+mx
            if( ((a(j)-h)*(a(jm)-h))  <= 0.0d0) then
               pnt_cnt=pnt_cnt+1
!              if(pnt_cnt > idim) goto 6000
               lx(pnt_cnt)=-j
            end if
         end do
      end do

      if(pnt_cnt < 2) return
! start a new curve
    !print*, "found ", pnt_cnt, "points"
! Here we start interpolating
! We start with the last point we found. After that we follow the path of the curve.
! When the curve is finished we add the starting point once more, if the curve is a closed
! loop. Then we start a new curve. The first point of the new curve is the one that has the
! highest index and has not been interpolated yet.
  
   in=pnt_cnt
   jx=lx(in) ! here jx is the complete coordinate in the matrix
   startcoord = jx
   jfor=0
   instart=in
   npl=0
   
   do ! 108 continue
      if(ncs == ncsmax) then
         ier=3
         return
      end if
      ncs=ncs+1
      in=pnt_cnt
      jx=lx(in)
      startcoord = jx
      jfor=0
      instart=in
      npl=0
! continue curve

   30 continue
      jabs=iabs(jx)
      kx=(jabs-1)/mx+1
      ikx= jabs-mx*(kx-1)
      !if(jx < 0) jnb=jabs-mx
      !if(jx > 0) jnb=jabs-1

      if(pnt_int == idim)then
        ier=2
        return
      endif

      pnt_int=pnt_int+1
      npl=npl+1
! interpolation for point found
! linear interpolation at boundaries, else quadr. interpolation.
! note: for method used, check only for upper boundaries required.
      if(jx < 0) then
         x1 = zcoord(kx - 1) ! kx is the z + 1 coordinate
         x2 = zcoord(kx)
         a1 = a(jabs - mx) ! first newton interpolation coeff.
         a2 = (a(jabs) - a1) / (x2 - x1) ! second coeff.
         if((jabs+mx) <= mxrt) then !(only if not at the boarder)
         ! Note that we can never have a penetration point at the left boarder,
         ! since penetration points are always left of the saved point
            x3 = zcoord(kx + 1) 
            a3 = ((a(jabs + mx) - a(jabs)) / (x3 - x2) - a2) &
                / ( x3 - x1)! third coeff.
            if(a3 /= 0.0d0) then
               dp =sqrt(a2 ** 2 + a3 *(-4 * a1 + 4 * h + a3 * (x1 - x2)**2) + &
                2 * a2 * a3 * (x1 - x2)) ! root in the
                  ! quadratic equation
               p1 = -(dp + a2 - a3 * x1 - a3 * x2)/(2 * a3) ! - solution
               p2 = (dp - a2 + a3 * x1 + a3 * x2)/(2 * a3) ! + solution
               if(p1 >= x1 .and. p1 <= x2) then ! in most cases only one of them
               ! is within the interval
               ! If both are, we arbitrarily take the first one
                  zpl(pnt_int) = p1
               else
                  zpl(pnt_int) = p2
               end if
            else
               ! linear interpolation if the curvature is zero
               if(a2 /= 0.0d0) p1 = (h - a1 + a2 * x1) / a2
               if(a2 == 0.0d0) p1 = x1 + (x2 - x1) * 0.5
               ! just the midpoint if the slope is also zero
               if(p1 /= 0.0d0) zpl(pnt_int) = p1
               if(p1 == 0.0d0) zpl(pnt_int) = x1 + (x2 - x1) * 0.5
            end if
         else
            if(a2 /= 0.0d0) p1 = (h - a1 + a2 * x1) / a2
            if(a2 == 0.0d0) p1 = x1 + (x2 - x1) * 0.5
            if(p1 /= 0.0d0) zpl(pnt_int) = p1
            if(p1 == 0.0d0) zpl(pnt_int) = x1 + (x2 - x1) * 0.5
         end if
         xpl(pnt_int) = xcoord(ikx)
      else
         x1 = xcoord(ikx - 1)
         x2 = xcoord(ikx)
         a1 = a(jabs - 1)
         a2 = (a(jabs) - a1) / (x2 - x1)
         if(ikx < nx) then
            x3 = xcoord(ikx + 1)
            a3 = ((a(jabs + 1) - a(jabs)) / (x3 - x2) - a2) &
                / ( x3 - x1)
            if(a3 /= 0.0d0) then
               dp =sqrt(a2 ** 2 + a3 *(-4 * a1 + 4 * h + a3 * (x1 - x2)**2) + &
                2 * a2 * a3 * (x1 - x2)) ! root in the
                  ! quadratic equation
               p1 = -(dp + a2 - a3 * x1 - a3 * x2)/(2 * a3)
               p2 = (dp - a2 + a3 * x1 + a3 * x2)/(2 * a3)
               if(p1 >= x1 .and. p1 <= x2) then
                  xpl(pnt_int) = p1
               else
                  xpl(pnt_int) = p2
               end if
            else
               if(a2 /= 0.0d0) p1 = (h - a1 + a2 * x1) / a2
               if(a2 == 0.0d0) p1 = x1 + (x2 - x1) * 0.5
               if(p1 /= 0.0d0) xpl(pnt_int) = p1
               if(p1 == 0.0d0) xpl(pnt_int) = x1 + (x2 - x1) * 0.5
            end if
         else
            if(a2 /= 0.0d0) p1 = (h - a1 + a2 * x1) / a2
            if(a2 == 0.0d0) p1 = x1 + (x2 - x1) * 0.5
            if(p1 /= 0.0d0) xpl(pnt_int) = p1
            if(p1 == 0.0d0) xpl(pnt_int) = x1 + (x2 - x1) * 0.5
         end if
         zpl(pnt_int) = zcoord(kx)
      end if
    !if( xpl(pnt_int) == 0.0d0 .or. zpl(pnt_int)== 0.0d0) print*, xcoord(ikx - 1), zcoord(kx - 1), &
    !ikx,kx, jx
    ! look for next point
! Here we create a look-up-table and then see if any of the points we found
! qualify as the next point within our curve.
! Here the order is important again.
! Since we have to test every point for the first direction first,
! before we can try another this step is the most expensive in this code
      found_neighbour = .false.
      is=1
      ie=6
      if(jx <= 0) then ! horizontal sign change
         pos_neighbour(1)=iabs(jx) ! base - vertical
         pos_neighbour(2)=-(iabs(jx)-1) ! top- horizontal
         pos_neighbour(3)=pos_neighbour(1)-mx ! left- vertical

         pos_neighbour(4)=-(iabs(jx)+1) ! bottom - horizontal
         pos_neighbour(5)=iabs(jx)+1 ! bottom - vertical
         pos_neighbour(6)=pos_neighbour(5)-mx ! bottom left - vertical
         if(ikx == 1) is=4
         if(ikx == nx) ie=3
      else ! vertical sign change
         pos_neighbour(1)=jx+mx ! right - vertical
         pos_neighbour(2)=-(jx+mx) ! right - horizontal
         pos_neighbour(3)=-(iabs(pos_neighbour(2))-1) ! right - top - horizontal
         pos_neighbour(4)=-jx ! base - horizontal
         pos_neighbour(5)=jx-mx ! left - vertical
         pos_neighbour(6)=-(jx-1) ! top - horizontal
         if(kx == 1) ie=3
         if(kx == nz) is=4
      end if
! make this to an internel function
      do  k=is,ie
        !if(pos_neighbour(k) /= jfor) then
            do i=lower_boundary,pnt_cnt
               if(lx(i) == pos_neighbour(k)) then
                  lx(in)=0
                  if(npl > 2) lx(instart) = startcoord
                  jfor=jx
                  jx=lx(i)
                  in = i
                  go to 30
               end if
            end do
         !end if
      end do
! hier exit
! hier 30 enddo
! curve found
      if(npl > 1) then
         npt(ncs)=npl
         pnt_int_last=pnt_int
      else
         !print*, "erased single point"
         ncs=ncs-1
         pnt_int=pnt_int_last
      end if

! remove used intersects from table lx and try to find another
! curve line for contour value

      if(in == instart) lx(instart)=0
      if(in > 0) lx(in)=0

!    6 continue
!      if(lx(pnt_cnt) == 0) then
!         pnt_cnt=pnt_cnt-1
!         if(pnt_cnt > lower_boundary) goto 6
!      end if
     do
       if(lx(pnt_cnt) /= 0) exit
       pnt_cnt=pnt_cnt-1
       if(pnt_cnt <= lower_boundary) exit
     enddo

!    7 continue
!      if(lx(lower_boundary) == 0) then
!         lower_boundary=lower_boundary+1
!         if(lower_boundary < pnt_cnt) goto 7
!      end if
     do
       if(lx(lower_boundary) /= 0) exit  
       lower_boundary=lower_boundary+1
       if(lower_boundary >= pnt_cnt) exit  
     enddo

!     if((pnt_cnt-lower_boundary+1) >= 2) goto 108
     if((pnt_cnt-lower_boundary+1) < 2) exit
   enddo ! 108


end subroutine contourgq_non_eqd

!*****************************************************************

subroutine contouring(A, h, contour, xcoord, zcoord)
! evaluates contour lines of A at amplitudes h(:)
! if present(xcoord, zcoord) contour lines in x/z coordinates
!   otherwise contour lines in real valued index coordinates (starts at 1 and ends at N)
use f90_kind
implicit none
real(rkind),    dimension(:,:), intent(in)    :: A          ! matrix
real(rkind),    dimension(:),   intent(in)    :: h          ! values for which contours in A are searched
type(contour_type),             intent(inout) :: contour
real(rkind),    dimension(:),   intent(in), optional :: xcoord, zcoord  ! coordinates of 1st(x) and 2nd(y) index of A
integer(ikind), parameter                  :: ncsmax = 20 ! dimension of output vector npt
integer(ikind), dimension(ncsmax)          :: npt         ! number of points for curves found
real(rkind),    dimension(:), allocatable  :: xpl, zpl    ! r- and z-values of plot points
integer(ikind), dimension(:), allocatable  :: lx          ! aux. vector of length idim
integer(ikind) :: mx   ! leading dimension of matrix A
integer(ikind) :: nx   ! number of grid points in r-direction, nx > 3
integer(ikind) :: nz   ! number of grid points in z-direction, nz > 3
integer(ikind) :: ier  ! error code. ier = 0  okay
                       !             ier = 1 : nx < 4  oder nz < 4
                       !             ier = 2 : idim insufficient
                       !             ier = 3 : ncsmax insufficient
integer(ikind) :: idim ! dimension of vectors  xpl, zpl, lx
integer(ikind) :: i, i1, i2, j, k

mx = size(A, dim=1)   ! first dimension of A
nx = mx               ! number of grid points in r-direction
nz = size(A, dim=2)   ! second dimension of A
!idim = 4*(nx+nz)      ! 2*(nx+nz) is not sufficient close to the separatrix
idim = 10*(nx+nz)      ! 2*(nx+nz) is not sufficient close to the separatrix
allocate(lx(idim), xpl(idim), zpl(idim))

! allocate structure contour
!---------------------------
contour%N_levels = size(h)
allocate(contour%level(contour%N_levels))

do k = 1, contour%N_levels
  contour%level(k)%amp = h(k)

  if (.not.present(xcoord)) then     ! contour lines as a function of real valued indices
    call contourgq(ier, A, mx, nx, nz, lx, idim, &
                   ncsmax, contour%level(k)%amp, xpl, zpl, contour%level(k)%N_curves, npt)
    xpl = xpl + 1   ! first indx is xpl=0
    zpl = zpl + 1

  else                               ! contour lines as a function of xcoord and zcoord
    if (.not.present(zcoord)) stop "sub contouring: either both or none of xcoord and zcoord should be present!"
    call contourgq_non_eqd(ier, A, xcoord, zcoord, mx, nx, nz, lx, idim, &
                           ncsmax, contour%level(k)%amp, xpl, zpl, contour%level(k)%N_curves, npt)
  endif

  allocate(contour%level(k)%curve( contour%level(k)%N_curves ))

  ! map: npt(:), xpl(:), zpl(:) to structure contour
  !---------------------------
  i1 = 1
  !print*,"jkasghdfkgahkjdhl",k,contour%level(k)%N_curves,npt(1:contour%level(k)%N_curves)
  do j = 1, contour%level(k)%N_curves
    contour%level(k)%curve(j)%N_pos         = npt(j)
    allocate(contour%level(k)%curve(j)%pos( contour%level(k)%curve(j)%N_pos ))
    i2                                      = i1 + contour%level(k)%curve(j)%N_pos - 1
    contour%level(k)%curve(j)%pos(:)%x_indx = xpl(i1:i2)
    contour%level(k)%curve(j)%pos(:)%y_indx = zpl(i1:i2)
    i1 = i2 + 1               ! next first index
  enddo

  !do i = 1, sum(npt(1:contour%level(k)%N_curves))
  !  write(6,'(i3,2f9.3)')i, xpl(i), zpl(i)
  !enddo
  !print*,h(k), maxval(A), minval(A)

  !write(6,'(a,i3)')"nuber of curves found =", contour%level(k)%N_curves
  !write(6,'(a,<contour%level(k)%N_curves>i5)')"nuber of points found =", npt(1:contour%level(k)%N_curves)
  !write(6,'(a,i3)')"contouring error =", ier

enddo ! k

deallocate(lx, xpl, zpl)

!stop "sub contouring"
end subroutine contouring

!*****************************************************************

subroutine contour_indx2rz(r_grid, z_grid, contour)
! evaluate for each contour point (r,z) from (x,y)_indx
use f90_kind
implicit none
real(rkind), dimension(:), intent(in)    :: r_grid, z_grid
type(contour_type),        intent(inout) :: contour
real(rkind)    :: dr, dz
integer(ikind) :: j, jr, jz, k, l

dr = r_grid(2) - r_grid(1)
dz = z_grid(2) - z_grid(1)

do l = 1, size(contour%level)
  do k = 1, contour%level(l)%N_curves
    do j = 1, size(contour%level(l)%curve(k)%pos)
      jr = INT(contour%level(l)%curve(k)%pos(j)%x_indx)
      jz = INT(contour%level(l)%curve(k)%pos(j)%y_indx)
      if (jr < 1 .or. jz < 1) then
        write(6,'(a,5i4,2f16.8)')"sub contour_indx2rz:",l, k, j, jr, jz, contour%level(l)%curve(k)%pos(j)%x_indx, contour%level(l)%curve(k)%pos(j)%y_indx
      endif

      contour%level(l)%curve(k)%pos(j)%r = r_grid(jr) + dr * (contour%level(l)%curve(k)%pos(j)%x_indx - jr)
      contour%level(l)%curve(k)%pos(j)%z = z_grid(jz) + dz * (contour%level(l)%curve(k)%pos(j)%y_indx - jz)
      !write(91,'(5i4,6f14.6)')l, k, j, jr, jz, contour%level(l)%curve(k)%pos(j)%r, contour%level(l)%curve(k)%pos(j)%z, &
      !    contour%level(l)%curve(k)%pos(j)%x_indx,  contour%level(l)%curve(k)%pos(j)%y_indx, &
      !    contour%level(l)%curve(k)%pos(j)%x_indx - jr,  contour%level(l)%curve(k)%pos(j)%y_indx - jz
    enddo
  enddo !k = 1, contour%level(l)%N_curves
enddo !l = 1, size(contour%level)

!stop "sub contour_indx2rz"
end subroutine contour_indx2rz

!*****************************************************************

subroutine map_contour_indx_to_rz(r_grid, z_grid, r_contour_point, z_contour_point, r, z)
use interpolation_routines,      only: linear_interpolation
implicit none
real(rkind), dimension(:),   intent(in)  :: r_grid, z_grid
real(rkind), dimension(:),   intent(in)  :: r_contour_point  ! r-values of points on contour line [real valued indx]
real(rkind), dimension(:),   intent(in)  :: z_contour_point  ! z-values of points on contour line [real valued indx]
real(rkind), dimension(:),   intent(out) :: r, z             ! [m] (r,z) coordinates of contour points
real(rkind)    :: krr, kzr
integer(ikind) :: i, kr, kz, Nr, Nz
integer(ikind) :: Nc             ! number of points on contour line

Nr = size(r_grid)                ! number of R points on grid
Nz = size(z_grid)                ! number of z points on grid
Nc = size(r_contour_point)       ! number of points on contour line

do i = 1, Nc
  kr  = INT(r_contour_point(i))
  kz  = INT(z_contour_point(i))
  krr = real(kr,kind(1.d0))       ! real version for interpolation
  kzr = real(kz,kind(1.d0))       ! real version for interpolation

  if (kr < 1 .or. kr > Nr-1) then
    write(6,'(a,i4,a,i4)')"kr =", kr, "   Nr-1 =", Nr-1
    stop "sub integration_on_closed_flux_surface: index kr not within [1, Nr-1]!"
  endif
  if (kz < 1 .or. kz > Nz-1) then
    write(6,'(a,i4,a,i4)')"kz =", kz, "   Nz-1 =", Nz-1
    stop "sub integration_on_closed_flux_surface: index kz not within [1, Nz-1]!"
  endif

  call linear_interpolation(krr, krr+1, r_grid(kr), r_grid(kr+1), r_contour_point(i), r(i))
  call linear_interpolation(kzr, kzr+1, z_grid(kz), z_grid(kz+1), z_contour_point(i), z(i))
  !write(6,'(i4,e12.4,i4,4e12.4)')i,r_contour_point(i),INT(r_contour_point(i)),r_grid(INT(r_contour_point(i))),r(i), r_grid(kr), r_grid(kr+1)
  !write(6,'(i4,e12.4,i4,4e12.4)')i,z_contour_point(i),INT(z_contour_point(i)),z_grid(INT(z_contour_point(i))),z(i), z_grid(kz), z_grid(kz+1)
enddo

!stop "sub map_contour_indx_to_rz"
end subroutine map_contour_indx_to_rz

!*****************************************************************




END MODULE mod_contour

!********************************************************************
!********************************************************************
!********************************************************************

