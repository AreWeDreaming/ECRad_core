module mod_ecfm_refr_interpol
    implicit none
#ifdef IDA
  public ::   spline_1d, & ! overloaded function
              make_rect_spline, &
              make_1d_spline, &
              deallocate_rect_spline, &
              deallocate_1d_spline, &
              rect_spline, &
              rect_spline_vec, &
              bispline_1d, &
              bispline_1d_vec, &
              spline_1d_get_roots, &
              spline_1d_integrate, &
              splint_1d, &
              splint_1d_vec
#else
  public ::   spline_1d, &
              make_rect_spline, &
              make_1d_spline, &
              deallocate_rect_spline, &
              deallocate_1d_spline, &
              rect_spline, &
              rect_spline_vec, &
              bispline_1d, &
              bispline_1d_vec, &
              spline_1d_get_roots, &
              spline_1d_integrate
#endif
  private :: print_2d_spline_params

    interface spline_1d
#ifdef IDA
      module procedure bispline_1d, bispline_1d_vec, splint_1d, splint_1d_vec
#else
      module procedure bispline_1d, bispline_1d_vec
#endif
    end interface spline_1d

    contains

  subroutine make_rect_spline(spl, m, n, x, y, mat, iopt, m_max)
    use f90_kind
    use mod_ecfm_refr_types,        only: spl_type_2d
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_2d), intent(inout)   :: spl
    integer*4, intent(in)            :: m, n
    real(rkind), dimension(:), intent(in) :: x,y
    real(rkind), dimension(:,:), intent(in) :: mat
    integer*4, intent(in), optional :: iopt
    integer*4, intent(in), optional :: m_max
    integer*4                                :: q, ier, kx, ky
    real*8                                   :: fp
    real*8, dimension(m*n)                   :: temp_mat
    integer(ikind)                           :: i, j
    kx = 3 ! this is hardcoded in the interpolation -> cubic splines
    ky = 3 ! this is hardcoded in the interpolation -> cubic splines
    if(.not. all(x(1:m - 1) < x(2:m))) then
      print*, "y has to be monotonically increasing for spline interpolation (2D)"
      print*, "Check input!"
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
    if(.not. all(y(1:n - 1) < y(2:n))) then
      print*, "y has to be monotonically increasing for spline interpolation (2D)"
      print*, "Check input!"
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
      stop "Failure in make_rect_spline"
#endif
    end if
    if(present(iopt)) then
      spl%iopt_int = int(iopt,4)
    end if
    if(spl%iopt_int == 0) then !
      spl%nuest= m + 4
      spl%nvest = n + 4
      spl%nu = spl%nuest
      spl%nv = spl%nvest
      q = n
      if(spl%nuest > q) q = spl%nuest
      spl%lwrk = 4+spl%nuest*(n+2*kx+5)+spl%nvest*(2*ky+5)+m*(kx+1)+ n*(ky+1) + q
      spl%kwrk = 3+m+n+spl%nuest+spl%nvest
      ! For derivative calculation
      ! Upper limit
      if(present(m_max)) then
        ! If we interpolate arrays that are larger than the matrix we are interpolating we need to replace m and n by the
        ! maximum length of the array we will be interpolating
        spl%lwrk_in_spl = max(m,m_max) * (kx + 1) +  max(n,m_max) * (ky + 1) + (max(spl%nuest,m_max) - kx - 1) * (max(spl%nvest,m_max) - ky - 1)
        spl%kwrk_in_spl = max(m,m_max) + max(n,m_max)
      else
        spl%lwrk_in_spl = m * (kx + 1) + n * (ky + 1) + (spl%nuest - kx - 1) * (spl%nvest - ky - 1)
        spl%kwrk_in_spl = m + n
      end if
      allocate(spl%tu(spl%nuest), spl%tv(spl%nvest), &
               spl%c(spl%nuest - kx - 1, spl%nvest - ky - 1), &
               spl%wrk(spl%lwrk), spl%iwrk(spl%kwrk))
    end if
    do i = 1, m
      do j = 1, n
        temp_mat(n*(i-1)+j) = mat(i,j)
      end do
    end do
    spl%x_start = x(1)
    spl%x_end = x(m)
    spl%y_start = y(1)
    spl%y_end = y(n)
    call regrid(spl%iopt_int, m, x, n, y, temp_mat, &
                x(1), x(m), y(1), y(n), kx, ky, 0.d0, &
                spl%nuest,  spl%nvest,  spl%nu,  spl%tu, &
                spl%nv, spl%tv,  spl%c, fp,  spl%wrk, &
                spl%lwrk,  spl%iwrk,  spl%kwrk, ier)
    if(spl%iopt_int == 0) spl%iopt_int = 1
    if(ier /= -1) then
      print*, "ier", ier
      print*, "Spline interpolation in make_rect_spline failed"
#ifdef INTEL
      call tracebackqq()
      call abort
#else
      call backtrace()
      stop "Failure in mod_interpolation"
#endif
    end if
  end subroutine make_rect_spline

  subroutine make_1d_spline(spl, m, x, y, iopt, k)
    use f90_kind
    use mod_ecfm_refr_types,        only: spl_type_1d
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_1d), intent(inout)   :: spl
    integer*4, intent(in)            :: m
    real(rkind), dimension(:), intent(in) :: x, y
    integer*4,   intent(in), optional     :: iopt
    integer*4,   intent(in), optional     :: k
    integer*4                                :: ier
    real*8                                   :: fp
    real*8, dimension(m)                     :: w
    if(present(k)) then
      if(k == 1 .or. k == 3) then
        spl%k = k
      else
        print*, "The order of the splines must be either linear or cubic"
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    else
      spl%k = 3
    end if
    if(m <= 1) then
      print*, "Cannot interpolate a single point"
      print*, "1D Interpolation called with singular point"
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
    if(.not. all(x(1:m - 1) < x(2:m))) then
      print*, "x has to be monotonically increasing for spline interpolation 1D"
      print*, "Check input!"
      print*, "x", x
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
    if(m <= spl%k) spl%k = 1
    if(present(iopt)) then
      spl%iopt_int = iopt
    end if
    w(:) = 1.d0
    if(spl%iopt_int == 0) then !
      if(allocated(spl%t)) then
        deallocate(spl%t, spl%c, spl%wrk, spl%iwrk)
      end if
      spl%nest= m + spl%k + 2*spl%k+2
      spl%n = spl%nest
      spl%lwrk = (spl%k + 1) *  m + spl%nest * (7 + 3 * spl%k)
      allocate(spl%t(spl%nest), &
               spl%c(spl%nest), spl%wrk(spl%lwrk), spl%iwrk(spl%nest))
    end if
    spl%x_start = x(1)
    spl%x_end = x(m)
    call curfit(spl%iopt_int, m, x, y, w, x(1), x(m), spl%k, 0.d0, spl%nest, spl%n,&
                spl%t, spl%c, fp, spl%wrk, spl%lwrk, spl%iwrk, ier)
    if(spl%iopt_int == 0) spl%iopt_int = 1
    if(ier /= -1) then
      print*, "x", x
      print*, "y", y
      print*, "ier", ier
      print*, "Spline interpolation in make 1d spline failed"
      print*, "m", m, "k", spl%k
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
  end subroutine make_1d_spline

  subroutine deallocate_rect_spline(spl)
    use f90_kind
    use mod_ecfm_refr_types,        only: spl_type_2d
    type(spl_type_2d), intent(inout)   :: spl
    if(allocated(spl%tu)) deallocate(spl%tu, spl%tv, spl%c, &
                                      spl%wrk, spl%iwrk)
    spl%iopt_int = 0
  end subroutine

  subroutine deallocate_1d_spline(spl)
    use f90_kind
    use mod_ecfm_refr_types,        only: spl_type_1d
    type(spl_type_1d), intent(inout)   :: spl
    if(allocated(spl%t)) deallocate(spl%t, spl%c, spl%wrk, spl%iwrk)
    spl%iopt_int = 0
  end subroutine

#ifdef NAG
  subroutine rect_spline(spl,x,y,f,dfdx,dfdy, nag_spline)
#else
  subroutine rect_spline(spl,x,y,f,dfdx,dfdy)
#endif
  ! Spline evaluation routine for B and Te/ne in 2D mode
  ! WARNING: This routine does not check bounds - out of bounds interpolations are prone to very large errors
    use f90_kind
#ifdef NAG
    USE mod_ecfm_refr_types , only : plasma_params_type, spl_type_2d, double_check_splines
#else
    USE mod_ecfm_refr_types , only : plasma_params_type, spl_type_2d
#endif
#ifdef NAG
    USE nag_spline_2d             , only: nag_spline_2d_eval, &
                                          nag_spline_2d_comm_wp => nag_spline_2d_comm_dp
    USE nag_error_handling
#endif
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_2d), intent(in)                :: spl
    real(rkind),                      intent(in)  :: x, y
    real(rkind),                      intent(out) :: f
    real(rkind),            intent(out), optional :: dfdx,dfdy
#ifdef NAG
    type(nag_spline_2d_comm_wp),    intent(in), optional :: nag_spline
    type(nag_error)                        :: error
    real(rkind)                             :: h_x, nag_val
#endif
    integer*4                               :: nux, nuy, kx, ky, m, ier
    real*8, dimension(spl%lwrk_in_spl)      :: real_work
    integer*4, dimension(spl%kwrk_in_spl)   :: integer_work
    real*8, dimension(1)                    :: x_ar, y_ar, val_ar
    x_ar(1) = x
    y_ar(1) = y
    m = 1
    kx = 3
    ky = 3
    if(present(dfdx)) then
      nux = 1
      nuy = 0
      call pardeu(spl%tu, spl%nu, spl%tv, spl%nv, spl%c, kx, ky, nux, nuy, x_ar, y_ar, val_ar,m,&
                real_work,spl%lwrk_in_spl,integer_work,spl%kwrk_in_spl,ier)
      dfdx = val_ar(1)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdx"
        print*, ier
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_ar
          print*, y_ar
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    if(present(dfdy)) then
      nux = 0
      nuy = 1
      call pardeu(spl%tu,spl%nu, spl%tv, spl%nv, spl%c, kx, ky, nux,nuy,x_ar, y_ar, val_ar,m,&
                real_work,spl%lwrk_in_spl,integer_work,spl%kwrk_in_spl,ier)
      dfdy = val_ar(1)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdy"
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_ar
          print*, y_ar
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    call bispeu(spl%tu,spl%nu,spl%tv,spl%nv,spl%c,kx,ky,x_ar, y_ar, val_ar,m,&
                real_work, spl%lwrk_in_spl, ier)
    f = val_ar(1)
    if(ier /= 0) then
      print*, "Critical error in spline evaluation: f"
      print*, "spline in mod_ecfm_refr_utils failed"
      if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_ar
          print*, y_ar
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
#ifdef NAG
    if(present(nag_spline) .and. double_check_splines .and. output_level) then
      CALL nag_set_error(error, halt_level=4)
      call nag_spline_2d_eval(nag_spline, x, y, nag_val, error=error)
      if(error%level > 1) then
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
      if(abs(nag_val - f) > 1.e-4 .and. abs(nag_val - f)/ abs(nag_val + f) > 1.e-4) then
        print*, "Large deviation between the two splines"
        print*, "nag", nag_val
        print*, "spline", f
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
#endif
  end subroutine rect_spline

  subroutine print_2d_spline_params(spl,m)
    use f90_kind
    USE mod_ecfm_refr_types , only : spl_type_2d
    implicit None
    type(spl_type_2d), intent(in) :: spl
    integer*4, intent(in)         :: m
    print*, "nu, nv, nuest, nvest, lwrk, kwrk, lwrk_in_spl, kwrk_in_spl"
    print*, spl%nu, spl%nv, spl%nuest, spl%nvest, spl%lwrk, spl%kwrk
    print*, "size tu, tv, c, wrk, iwrk"
    print*, size(spl%tu), size(spl%tv),  size(spl%c), size(spl%wrk), size(spl%iwrk)
    print*, "Required work array size (for dx =1 and dy =1):", (spl%nu - 3)  *(spl%nv - 3)+(3-1) * m + (3 - 1) * m
  end subroutine print_2d_spline_params

#ifdef NAG
  subroutine rect_spline_vec(spl,x_vec,y_vec,f,dfdx,dfdy, nag_spline)
#else
  subroutine rect_spline_vec(spl,x_vec,y_vec,f,dfdx,dfdy)
#endif
  ! Spline evaluation routine
    use f90_kind
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
#ifdef NAG
    USE nag_spline_2d             , only: nag_spline_2d_eval, &
                                          nag_spline_2d_comm_wp => nag_spline_2d_comm_dp
    USE nag_error_handling
    USE mod_ecfm_refr_types , only : plasma_params_type, spl_type_2d, double_check_splines, output_level
#else
    USE mod_ecfm_refr_types , only : plasma_params_type, spl_type_2d
#endif
    implicit none
    type(spl_type_2d), intent(in)                 :: spl
    real(rkind), dimension(:),        intent(in)  :: x_vec, y_vec
    real(rkind), dimension(:),        intent(out) :: f
    real(rkind), dimension(:), intent(out), optional :: dfdx,dfdy
#ifdef NAG
    type(nag_spline_2d_comm_wp),    intent(in), optional :: nag_spline
    type(nag_error)                        :: error
#endif
    integer*4                               :: nux, nuy, kx, ky, m, ier
    real*8, dimension(spl%lwrk_in_spl)      :: real_work
    integer*4, dimension(spl%kwrk_in_spl)   :: integer_work
#ifdef NAG
    real(rkind), dimension(size(f))          :: nag_vals
#endif
    m = size(x_vec)
    kx = 3
    ky = 3
    if(present(dfdx)) then
      nux = 1
      nuy = 0
      call pardeu(spl%tu,spl%nu,spl%tv,spl%nv,spl%c,kx,ky,nux,nuy,x_vec, y_vec, dfdx,m,&
                real_work,spl%lwrk_in_spl,integer_work,spl%kwrk_in_spl,ier)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdx"
        print*, ier
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_vec
          print*, y_vec
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
        end if
        call print_2d_spline_params(spl, m)
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    if(present(dfdy)) then
      nux = 0
      nuy = 1
      call pardeu(spl%tu,spl%nu,spl%tv,spl%nv,spl%c,kx,ky,nux,nuy,x_vec, y_vec, dfdy,m,&
                real_work,spl%lwrk_in_spl,integer_work,spl%kwrk_in_spl,ier)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdy"
        print*, ier
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_vec
          print*, y_vec
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
          if(any(x_vec < spl%x_start)) print*, "some x values smaller than the minimum"
          if(any(x_vec > spl%x_end)) print*, "some x values larger than the maximum"
          if(any(y_vec < spl%y_start)) print*, "some y values smaller than the minimum"
          if(any(y_vec > spl%y_end)) print*, "some y values smaller than the maximum"
        end if
        call print_2d_spline_params(spl,m)
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    call bispeu(spl%tu,spl%nu,spl%tv,spl%nv,spl%c,kx,ky,x_vec, y_vec, f,m,&
                real_work,spl%lwrk_in_spl,ier)
    if(ier /= 0) then
      print*, "Critical error in spline evaluation: f"
      print*, ier
      print*, "spline in mod_ecfm_refr_utils failed"
      if(ier == 10) then
          print*, "Array value out of bounds?"
          print*, x_vec
          print*, y_vec
          print*, "x boundary", spl%x_start, spl%x_end
          print*, "y boundary", spl%y_start, spl%y_end
        end if
      call print_2d_spline_params(spl,m)
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
#ifdef NAG
    if(present(nag_spline) .and. double_check_splines .and. output_level) then
      CALL nag_set_error(error, halt_level=4)
      call nag_spline_2d_eval(nag_spline, x_vec, y_vec, nag_vals, error=error)
      if(error%level > 1) then
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
      if(any(abs(nag_vals - f) > 1.e-3) .and. any(abs(nag_vals - f)/ abs(nag_vals + f) > 1.e-3)) then
        print*, "Large deviation between the two splines"
        print*, "nag", nag_vals
        print*, "spline", f
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
#endif
  end subroutine rect_spline_vec

#ifdef IDA
  subroutine splint_1d(knot_pos, val, deriv2, x, y, dydx)
  ! Spline evaluation routine for B and Te/ne in 2D mode
  ! WARNING: This routine does not check bounds - out of bounds interpolations are prone to very large errors
    use f90_kind
    USE nr_spline, only  : splint, splintg
    implicit none
    real(rkind), dimension(:),        intent(in)  :: knot_pos, val, deriv2
    real(rkind),                      intent(in)  :: x
    real(rkind),                      intent(out) :: y
    real(rkind),            intent(out), optional :: dydx
    if(present(dydx)) then
      call splintg(knot_pos, val, deriv2, x, y, dydx)
    else
      y = splint(knot_pos, val, deriv2, x)
    end if
  end subroutine splint_1d

  subroutine splint_1d_vec(knot_pos, val, deriv2, x, y, dydx)
  ! Spline evaluation routine for B and Te/ne in 2D mode
  ! WARNING: This routine does not check bounds - out of bounds interpolations are prone to very large errors
    use f90_kind
    implicit none
    real(rkind), dimension(:),        intent(in)  :: knot_pos, val, deriv2
    real(rkind), dimension(:),        intent(in)  :: x
    real(rkind), dimension(:),        intent(out) :: y
    real(rkind), dimension(:), intent(out), optional :: dydx
    integer(ikind)                                   :: i
    do i= 1, size(x)
      if(present(dydx)) then
        call splint_1d(knot_pos, val, deriv2, x(i), y(i), dydx(i))
      else
        call splint_1d(knot_pos, val, deriv2, x(i), y(i))
      end if
    end do
  end subroutine splint_1d_vec
#endif

#ifdef NAG
  subroutine bispline_1d(spl,x,f,dfdx, nag_spline)
#else
  subroutine bispline_1d(spl,x,f,dfdx)
#endif
  ! Spline evaluation routine for B and Te/ne in 2D mode
  ! WARNING: This routine does not check bounds - out of bounds interpolations are prone to very large errors
    use f90_kind
#ifdef NAG
    USE mod_ecfm_refr_types , only : spl_type_1d, double_check_splines
#else
    USE mod_ecfm_refr_types , only : spl_type_1d
#endif
#ifdef NAG
    USE nag_spline_1d        , only : nag_spline_1d_eval, &
                                      nag_spline_1d_comm_wp => nag_spline_1d_comm_dp
    USE nag_error_handling
#endif
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_1d), intent(in)         :: spl
    real(rkind),                      intent(in)  :: x
    real(rkind),                      intent(out) :: f
    real(rkind),            intent(out), optional :: dfdx
#ifdef NAG
    type(nag_spline_1d_comm_wp),    intent(in), optional :: nag_spline
    type(nag_error)                        :: error
    real(rkind)                             :: nag_val
#endif
    real*8, dimension(spl%n)                :: real_work
    integer*4                               :: nx, m, ier
    real*8, dimension(1)                    :: x_ar, val_ar
    x_ar(1) = x
    m = 1
    if(present(dfdx)) then
      nx = 1
      call splder(spl%t, spl%n, spl%c, spl%k, nx, x_ar, val_ar, m, 2, real_work, ier)
      dfdx = val_ar(1)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdx"
        print*, ier
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds"
          print*, "requested value", x_ar
          print*, "x boundary", spl%x_start, spl%x_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    call splev(spl%t, spl%n, spl%c, spl%k, x_ar, val_ar, m, 2, ier)
    f = val_ar(1)
    if(ier /= 0) then
      print*, "Critical error in spline evaluation: f"
      print*, "spline in mod_ecfm_refr_utils failed"
      if(ier == 10) then
        print*, "Array value out of bounds"
        print*, "requested value", x_ar
        print*, "x boundary", spl%x_start, spl%x_end
      end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
#ifdef NAG
    if(present(nag_spline) .and. double_check_splines .and. output_level) then
      CALL nag_set_error(error, halt_level=4)
      call nag_spline_1d_eval(nag_spline, x, nag_val, error=error)
      if(error%level > 1) then
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
      if(abs(nag_val - f) > 1.e-3 .and. abs(nag_val - f)/ abs(nag_val + f) > 1.e-3) then
        print*, "Large deviation between the two splines"
        print*, "nag", nag_val
        print*, "spline", f
        if(ier == 1) then
          print*, "Array value out of bounds"
          print*, "requested value",  x_ar
          print*, "x boundary", spl%x_start, spl%x_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
#endif
  end subroutine bispline_1d

#ifdef NAG
  subroutine bispline_1d_vec(spl, x, f, dfdx, nag_spline)
#else
  subroutine bispline_1d_vec(spl, x, f, dfdx)
#endif
  ! Spline evaluation routine for B and Te/ne in 2D mode
  ! WARNING: This routine does not check bounds - out of bounds interpolations are prone to very large errors
    use f90_kind
#ifdef NAG
    USE nag_spline_1d        , only : nag_spline_1d_eval, &
                                      nag_spline_1d_comm_wp => nag_spline_1d_comm_dp
    USE nag_error_handling
    USE mod_ecfm_refr_types , only  : plasma_params_type, spl_type_1d, double_check_splines
#else
    USE mod_ecfm_refr_types , only  : plasma_params_type, spl_type_1d
#endif
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_1d), intent(in)         :: spl
    real(rkind), dimension(:),        intent(in)     :: x
    real(rkind), dimension(:),        intent(out)    :: f
    real(rkind), dimension(:), intent(out), optional :: dfdx
#ifdef NAG
    type(nag_spline_1d_comm_wp),    intent(in), optional :: nag_spline
    type(nag_error)                        :: error
    real(rkind), dimension(size(x))         :: nag_val
#endif
    real*8, dimension(spl%n)                :: real_work

    integer*4                               :: nx, m, ier
    m = int(size(x), 4)
    if(present(dfdx)) then
      nx = 1
      call splder(spl%t, spl%n, spl%c, spl%k, nx, x, dfdx, m, 2, real_work, ier)
      if(ier /= 0) then
        print*, "Critical error in spline evaluation: dfdx"
        print*, ier
        print*, "spline in mod_ecfm_refr_utils failed"
        if(ier == 10) then
          print*, "Array value out of bounds"
          print*, "requested value", x
          print*, "x boundary", spl%x_start, spl%x_end
        end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
    call splev(spl%t, spl%n, spl%c, spl%k, x, f, m, 2, ier)
    if(ier /= 0) then
      print*, "Critical error in spline evaluation: f"
      print*, "spline in mod_ecfm_refr_utils failed"
      print*, "error message", ier
      if(ier == 10) then
        print*, "Array value out of bounds"
        print*, "requested value", x
        print*, "x boundary", spl%x_start, spl%x_end
      end if
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
#ifdef NAG
    if(present(nag_spline) .and. double_check_splines .and. output_level) then
      CALL nag_set_error(error, halt_level=4)
      call nag_spline_1d_eval(nag_spline, x, nag_val, error=error)
      if(error%level > 1) call tracebackqq()
      if(any(abs(nag_val - f) > 1.e-3 .and. abs(nag_val - f)/ abs(nag_val + f) > 1.e-3)) then
        print*, "Large deviation between the two splines"
        print*, "nag", nag_val
        print*, "spline", f
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
      end if
    end if
#endif
  end subroutine bispline_1d_vec

  subroutine spline_1d_get_roots(spl, roots, root_cnt)
  ! Finds the root of a 1D spline
    use f90_kind
    USE mod_ecfm_refr_types , only  : spl_type_1d
#ifdef INTEL
    use ifcore,                     only: tracebackqq
#endif
    implicit none
    type(spl_type_1d)                     :: spl
    real(rkind), dimension(:),        intent(out)    :: roots
    integer(ikind),                   intent(out)    :: root_cnt
    integer*4                                        :: m_root, ier
    if(spl%k /= 3) then
      print*, "The root search only works for cubic splines"
      print*, "Order of the given spline", spl%k
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
    call sproot(spl%t,spl%n,spl%c, roots, int(size(roots), 4), m_root, ier)
    !print*, "roots", roots(1:3)
    root_cnt = int(m_root,8)
    if(ier /= 0) then
      print*, "Critical error in evaluation of spline roots"
      print*, ier
      print*, "spline in mod_ecfm_refr_utils failed"
      print*, "t", spl%t
#ifdef INTEL
      call tracebackqq()
#else
      call backtrace()
#endif
      call abort()
    end if
  end subroutine spline_1d_get_roots

  subroutine spline_1d_integrate(spl, a, b, int_val)
  ! Definite integral a, b of spline
    use f90_kind
    USE mod_ecfm_refr_types , only  : spl_type_1d
    implicit none
    type(spl_type_1d)                     :: spl
    real(rkind),           intent(in)     :: a, b
    real(rkind),           intent(out)    :: int_val
    real(rkind), dimension(spl%n)         :: wrk
    real(rkind), external :: splint_fitpack
    int_val = splint_fitpack(spl%t,spl%n,spl%c, spl%k, a, b, wrk)
  end subroutine spline_1d_integrate
end module mod_ecfm_refr_interpol
