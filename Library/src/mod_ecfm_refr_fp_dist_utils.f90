module mod_ecfm_refr_fp_dist_utils
  use f90_kind
  implicit none

contains


subroutine make_rhop_Bmin()
  USE mod_ecfm_refr_types, only: ant, rad, plasma_params, data_folder, OERT, N_freq, N_ray, ffp, output_level
  use constants,           only: mass_e, e0, pi
  use mod_ecfm_refr_interpol,    only: make_1d_spline, rect_spline
  use mod_contour, only: contour_type, contouring, contour_indx2rz
#ifdef NAG
  USE nag_spline_1d,                only: nag_spline_1d_interp
#endif
  implicit none
  type(contour_type)          :: rhop_contour
  integer(ikind)              :: iint, icurve, ipnt
  real(rkind)                 :: cur_B_min
  real(rkind), dimension(3)   :: R_vec, B_vec
  allocate(ffp%rhop_B_min(ffp%N_B_min), ffp%B_min(ffp%N_B_min))
  do iint = 1, ffp%N_B_min
    ffp%rhop_B_min(iint) = (real(iint) - 1) /  (ffp%N_B_min - 1)
  end do
  call contouring(plasma_params%rhop, ffp%rhop_B_min, rhop_contour)
  call contour_indx2rz(plasma_params%R, plasma_params%z, rhop_contour)
  ffp%B_min(:) = 200.d0 ! Large value for Bt for which we'll certainly find smaller ones
  do iint = 1, ffp%N_B_min
    do icurve = 1, rhop_contour%level(iint)%N_curves
      do ipnt = 1, rhop_contour%level(iint)%curve(icurve)%N_pos
        R_vec(1) = rhop_contour%level(iint)%curve(icurve)%pos(ipnt)%R
        R_vec(3) = rhop_contour%level(iint)%curve(icurve)%pos(ipnt)%z
        if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
           R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
              cur_B_min = 200.d0 ! Give positions outside the plasma a large value so they are not considered
        else
#ifdef NAG
          if(plasma_params%debug_level == 0 .and. output_level) then
            call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                        B_vec(1))
            call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                        B_vec(2))
            call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                        B_vec(3))
          else
            call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                        B_vec(1), nag_spline=plasma_params%B_r_spline_nag)
            call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                        B_vec(2), nag_spline=plasma_params%B_t_spline_nag)
            call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                        B_vec(3), nag_spline=plasma_params%B_z_spline_nag)
          end if
#else
          call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                      B_vec(1))
          call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                      B_vec(2))
          call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                      B_vec(3))
#endif
          cur_B_min = sqrt(sum(B_vec**2))
        end if
        if(cur_B_min < ffp%B_min(iint)) ffp%B_min(iint) = cur_B_min
      end do
    end do
  end do
  do iint = ffp%N_B_min, 1, -1
    ! Contouring fails near magnetic axis
    ! To fix this we go from the outer point from B_min and extrapolate inwards for all points were contouring failed
    if(ffp%B_min(iint) > 10.d0) then
      if( iint + 1 > ffp%N_B_min) then
        print*, "Failed to find B_min for rhop = 1.0 - cannot set up B_min"
        print*, ffp%rhop_B_min(iint)
        print*, ffp%B_min(iint)
        call abort()
      else if(ffp%B_min(iint + 1) > 10.d0) then
        print*, "Failed to fix B_min using the next outer neighbor"
        print*, "point", ffp%rhop_B_min(iint), ffp%B_min(iint)
        print*, "neighbor", ffp%rhop_B_min(iint + 1), ffp%B_min(iint + 1)
        call abort()
      else if(iint + 2 > ffp%N_B_min) then
        ffp%B_min(iint) = ffp%B_min(iint + 1) ! directly replace with neighbor if only one good point next
      else
        ! Extrapolate from last two points
        ffp%B_min(iint) = ffp%B_min(iint + 1) - (ffp%B_min(iint + 2) - ffp%B_min(iint + 1)) / (ffp%rhop_B_min(iint + 2) - ffp%rhop_B_min(iint + 1)) * (ffp%rhop_B_min(iint + 1) - ffp%rhop_B_min(iint))
      end if
    end if
  end do
!  do iint = 1, ffp%N_B_min
!    print*, "rhop, B_min", ffp%rhop_B_min(iint), ffp%B_min(iint)
!  end do
  call make_1d_spline(ffp%B_min_spl, int(ffp%N_B_min, kind=4), ffp%rhop_B_min, ffp%B_min)
#ifdef NAG
  call nag_spline_1d_interp(ffp%rhop_B_min, ffp%B_min, ffp%B_min_nag_spl)
#endif
end subroutine make_rhop_Bmin

subroutine setup_f_rhop_splines(ffp)
  use mod_ecfm_refr_types,          only: ffp_type, N_absz_large, double_check_splines
  use mod_ecfm_refr_interpol,       only: make_1d_spline
#ifdef NAG
  USE nag_spline_1d,                only: nag_spline_1d_interp
#endif
  implicit none
  type(ffp_type), intent(inout)   :: ffp
  integer(ikind)                  :: iu, ipitch
  ! N_u x N_pitch splines are set up along rhop
  allocate(ffp%f_rhop_spl(ffp%N_u, ffp%N_pitch))
!#ifdef NAG
!  if(double_check_splines) then
!    allocate(ffp%f_rhop_nag_spl(ffp%N_u, ffp%N_pitch))
!  end if
!#endif
  do iu = 1, ffp%N_u
    do ipitch = 1, ffp%N_pitch
      call make_1d_spline(ffp%f_rhop_spl(iu, ipitch), ffp%N_rhop, ffp%rhop, ffp%f(:, iu, ipitch), k = 1)
!#ifdef NAG
!      if(double_check_splines) then
!        call nag_spline_1d_interp(ffp%rhop, ffp%f(:, iu, ipitch), ffp%f_rhop_nag_spl(iu, ipitch))
!      end if
!#endif
    end do
  end do
end subroutine setup_f_rhop_splines

subroutine make_B_min_and_f_inter(svec, f_spl, B_min)
  use mod_ecfm_refr_types,          only: ffp, rad_diag_ch_mode_ray_freq_svec_type, N_absz_large, double_check_splines, spl_type_2d
  use mod_ecfm_refr_interpol,       only: make_rect_spline, deallocate_rect_spline, rect_spline, rect_spline_vec, spline_1d
#ifdef NAG
  USE nag_spline_2d,                only: nag_spline_2d_interp
#endif
  use constants,                    only: pi, mass_e, e0, c0
  implicit none
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  type(spl_type_2d), intent(inout) :: f_spl
  real(rkind)                      :: B_min
  integer(ikind)                   :: irhop, iu, ipitch, i
  real(rkind)                      :: rhop_best, u, f, dfdu, mu, a, c
  real(rkind), dimension(30)       :: test_u, test_pitch, test_f, test_df_du
  integer*4                        :: ier, kx, ky
  real*8                           :: fp
  real(rkind), dimension(ffp%N_u, ffp%N_pitch) :: f_inter
  rhop_best = svec%rhop
  if(rhop_best > ffp%rhop_max) return ! automatic switch to thermal distributions when the distribution function is evaluated
  if(rhop_best < ffp%rhop_min) rhop_best = ffp%rhop_min ! if rhop grid sensible this should be fine
  do iu = 1, ffp%N_u
    do ipitch = 1, ffp%N_pitch
!#ifdef NAG
!      if(double_check_splines) then
!        call spline_1d(ffp%f_rhop_spl(iu, ipitch), rhop_best, f_inter(iu, ipitch), nag_spline = ffp%f_rhop_nag_spl(iu, ipitch))
!      else
!        call spline_1d(ffp%f_rhop_spl(iu, ipitch), rhop_best, f_inter(iu, ipitch))
!      end if
!#else
      call spline_1d(ffp%f_rhop_spl(iu, ipitch), rhop_best, f_inter(iu, ipitch))
!#endif
    end do
  end do
!  if(ffp%restart_spline) then
  f_spl%iopt_int = 0 ! restart
  call deallocate_rect_spline(f_spl)
  !end if
  call make_rect_spline(f_spl, int(ffp%N_u,4), int(ffp%N_pitch,4), &
       ffp%u, ffp%pitch, f_inter, iopt = 0, m_max=N_absz_large)
!#ifdef NAG
!  if(double_check_splines) then
!    call nag_spline_2d_interp(ffp%u, ffp%pitch, f_inter, ffp%f_nag_spl)
!  end if
!#endif
!  if(ffp%restrict_u_max) then
!    !ffp%restart_spline = .false.
!    u  = sqrt((1.d0 + svec%Te * e0 / (mass_e* c0**2))**2 - 1.d0) * 3.0/ 2.0
!    test_u(:) = u
!    do i = 1, size(test_u)
!      test_pitch(i) = ffp%pitch_min + real(i - 1,8) / real(size(test_u) - 1,8) * ffp%pitch_min
!    end do
!#ifdef NAG
!    if(double_check_splines) then
!      call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f, test_df_du, nag_spline = ffp%f_nag_spl)
!    else
!      call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f, test_df_du)
!    end if
!#else
!    call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f, test_df_du)
!#endif
    ! Determine gradient near origin - a proxy for the temperature
    ! since we use the logarithm of f the f cancels out
!    mu = -sqrt(1.d0 + u**2) / u * minval(test_df_du)
!    !print*, "non-thermal mu", mu
!    a = 1.0 / (1 + 105.0 / (128.0 * mu ** 2) + 15.0 / (8.0 * mu))
!    ! Smallest value of a distribution that can still be interpolated correctly
!    c = ffp%f_min / (a * sqrt(mu / (2.0 * pi)) ** 3)
!    ! Set corresponding maximum u
!    ffp%u_max = min(sqrt((1.d0 - log(c) / mu)**2 - 1.d0), ffp%u_max_data)
!    test_u(:) = ffp%u_max
!    if(ffp%u_max < ffp%u_min) then
!      print*, "u_min > u_max???"
!      print*, "Te", svec%Te
!      print*, "Te from distribution", mass_e* c0**2 / e0 / mu
!      print*, "df/du", dfdu * f
!      stop "Something wrong with the distribution?"
!    end if
!    do i = 1, size(test_u)
!      test_pitch(i) = ffp%pitch_min + real(i - 1,8) / real(size(test_u) - 1,8) * ffp%pitch_min
!    end do
!#ifdef NAG
!    if(double_check_splines) then
!      call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f, nag_spline = ffp%f_nag_spl)
!    else
!      call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f)
!    end if
!#else
!    call rect_spline_vec(ffp%f_spl, test_u, test_pitch, test_f)
!#endif
!    if(any(exp(test_f) > 1.d-8)) then
!      print*, "The estimation criterion for u_max cut of a significant part of the distribution"
!      print*, "test_u", ffp%u_max
!      print*, "test_f", maxval(exp(test_f))
!      print*, "limit", 1.d-8
!      print*, "Te", svec%Te
!      print*, "rhop", svec%rhop
!      print*, "Te from distribution", mass_e* c0**2 / e0 / mu
!      stop "Rework estimation criterion!"
!  !  else
!  !    print*, "Current Te and u_max", svec%Te, sqrt((1.d0 - log(c) / mu)**2 - 1.d0)
!    end if
!    ffp%f_spl%iopt_int = 0 ! restart
!    call deallocate_rect_spline(ffp%f_spl)
!    i = 1
!    do while(ffp%u(i) < ffp%u_max)
!      i = i + 1
!    end do
!    call make_rect_spline(ffp%f_spl, int(size(ffp%u(1:i)),4), int(size(ffp%pitch),4), &
!         ffp%u(1:i), ffp%pitch, f_inter(1:i, :), iopt = ffp%f_spl%iopt_int)
!#ifdef NAG
!    if(double_check_splines) then
!      call nag_spline_2d_interp(ffp%u(1:i), ffp%pitch, ffp%f_inter(1:i, :), ffp%f_nag_spl)
!    end if
!#endif
!  else
!    ffp%u_max = ffp%u_max_data
!  end if

#ifdef NAG
    call spline_1d(ffp%B_min_spl, svec%rhop, B_min, nag_spline = ffp%B_min_nag_spl)
#else
    call spline_1d(ffp%B_min_spl, svec%rhop, B_min)
#endif
end subroutine make_B_min_and_f_inter

subroutine cyl_to_pol(u_par, u_perp, svec, B_min, u, pitch)
  use mod_ecfm_refr_types,          only: ffp,rad_diag_ch_mode_ray_freq_svec_type
  use constants,                    only: mass_e, pi, e0, c0
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  real(rkind), intent(in)                       :: B_min
  real(rkind),  dimension(size(u_par)), intent(out)         :: u, pitch
  real(rkind), dimension(size(u_par))                       :: mu
  real(rkind)                     :: zeta
  zeta = pi*svec%freq_2X*mass_e/ (e0 * B_min)
  if(zeta < 1.d0 .or. B_min == 0.d0 .or. ffp%LUKE) zeta = 1.d0
  u(:) = sqrt(u_par(:)**2 + u_perp(:)**2)
  if(zeta == 1.d0) then
    mu(:) = u_par(:)/u(:)
  else
    mu(:) = sign(sqrt((u_par(:)**2 + u_perp(:)**2*(zeta - 1.0)/zeta))/u(:) , u_par(:))
  end if
  where(u(:) == 0.0) mu(:) = 0.0
  where(mu(:) < -1.0) mu(:) = mu(:) + 2
  where(mu(:) > 1.0) mu(:) = mu(:) - 2
  pitch = acos(mu)
  where(pitch(:)  < ffp%pitch(1)) pitch(:) = pitch(:)  + pi ! should only cause a small error
  where(pitch(:)  > ffp%pitch(ffp%N_pitch)) pitch(:) = pitch(:) - pi ! should only cause a small error
  where(u > ffp%u_max) u(:) = ffp%u_max!
  where(u > ffp%u_max) pitch(:) = pi/2.d0
  if(any(pitch /= pitch)) then
    print*, "Nan in coordinate transform!"
    print*, "u_par, u_perp, zeta,", u_par, u_perp, zeta
    stop "Nan in coordinate transform in mod_ecfm_refr_fp_dist_utils.f90"
  end if
end subroutine cyl_to_pol

subroutine make_f_and_f_grad_along_line(u_par, u_perp, svec, f_spl, B_min, f, df_du_par, df_du_perp, debug)
  use mod_ecfm_refr_types,          only: ffp,rad_diag_ch_mode_ray_freq_svec_type, double_check_splines, spl_type_2d
  use constants,                    only: mass_e, pi, e0, c0
  use mod_ecfm_refr_interpol,       only: rect_spline_vec
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  type(spl_type_2d), intent(in) :: f_spl
  real(rkind), intent(in)       :: B_min
  real(rkind), dimension(:), intent(out) :: f, df_du_par, df_du_perp
  logical, intent(in), optional          :: debug
  logical                                :: dbg
  real(rkind), dimension(size(u_par))    :: sign_array
  real(rkind)                     :: interpolate_f_u
  integer(ikind)                  :: i
  real(rkind), dimension(size(u_par)) :: temp_u, temp_pitch, u, pitch, df_du, df_dpitch,&
                                         du_du_par, du_du_perp, &
                                         dpitch_du_par, dpitch_du_perp, u_step, pitch_step
  real(rkind)                         :: zeta, h, mu, a, norm
  integer*4                           :: ier, kx, ky, nux, nuy, m
  sign_array(:) = 1.d0
  sign_array = sign(sign_array, u_par)
  dbg = .false.
  if(present(debug)) dbg = debug
  m = int(size(u_par), 4)
  call cyl_to_pol(u_par, u_perp, svec, B_min, u, pitch)
  mu = mass_e * c0**2 / (e0 * svec%Te)
  a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
  norm = a * (sqrt(mu / (2 * pi))**3)
  temp_u = u
  temp_pitch = pitch
  ! Using temp_u and temp_pitch to avoid interpolation errors
  ! Values outside the range where the non-thermal distribution is not defined
  ! will be correctly replaced by values from a thermal distribution later
  where(temp_u <= ffp%u_min) temp_u = ffp%u_min + 1.d-11
  where(temp_u >= ffp%u_max) temp_u = ffp%u_max - 1.d-11
  where(temp_pitch <= ffp%pitch_min .or. temp_pitch >= ffp%pitch_max) temp_pitch = 0.d0
!#ifdef NAG
!  if(double_check_splines) then
!    call rect_spline_vec(f_spl, temp_u, temp_pitch, f, df_du, df_dpitch, ffp%f_nag_spl)
!  else
!    call rect_spline_vec(f_spl, temp_u, temp_pitch, f, df_du, df_dpitch)
!  end if
!#else
  call rect_spline_vec(f_spl, temp_u, temp_pitch, f, df_du, df_dpitch)
!#endif
  f = exp(f)
  zeta = pi*svec%freq_2X*mass_e/ (e0 * B_min)
  if(zeta < 1.d0  .or. ffp%LUKE) then
      zeta = 1.d0
  end if
  h = (zeta - 1.d0)/zeta
  if(B_min /= 0.0) then
     !dpitch_du_par = ((u_par*Sqrt((u**2 - u_par**2 - h*u_perp**2)*(u_par**2 + h*u_perp**2)))/ &
     !                 (u**2*(u_par**2 + h*u_perp**2)))
     !dpitch_du_perp =  -((u_perp*(u_par**2 + h*(-u**2 + u_perp**2)))/ &
     !                  (u**2*Sqrt((u**2 - u_par**2 - h*u_perp**2)*(u_par**2 + h*u_perp**2))))
     dpitch_du_par = -((u_par*Sqrt((1.d0 - h)*u_perp**2))/(u**2*Sqrt(u_par**2 + h*u_perp**2)))
     dpitch_du_perp = (u_par**2*Sqrt((1.d0 - h)*u_perp**2))/(u**2*u_perp*Sqrt(u_par**2 + h*u_perp**2))
     dpitch_du_par = dpitch_du_par*sign_array
     dpitch_du_perp = dpitch_du_perp*sign_array
  else
    dpitch_du_par = -u_perp/(u_par**2 + u_perp**2)
    dpitch_du_perp = u_par/(u_par**2 + u_perp**2)
  end if
  du_du_par = u_par/ u
  du_du_perp = u_perp / u
  if(dbg) then
    if(B_min /= 0.0) then
      print*,"B_min not zero"
    else
      print*,"B_min zero"
    end if
    open(93, file="finite_diff_par.txt")
    open(94, file="finite_diff_perp.txt")
    call cyl_to_pol(u_par + 1.d-4, u_perp, svec, B_min, u_step, pitch_step)
    do i = 1, size(u_par)
      !print*,"dpar", u_par(i), dpitch_du_par(i), (pitch_step - pitch(i)) / 1.d-4, du_du_par(i), (u_step - u(i)) / 1.d-4
      write(93,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)"), &
              u_par(i), " ",&
              dpitch_du_par(i),  " ",&
              (pitch_step(i) - pitch(i)) / 1.d-4, " ", &
              du_du_par(i), " ", &
              (u_step(i) - u(i)) / 1.d-4

    end do
    call cyl_to_pol(u_par, u_perp + 1.d-4, svec, B_min, u_step, pitch_step)
    do i = 1, size(u_par)
      !print*, "dperp",u_perp(i), dpitch_du_perp(i), (pitch_step - pitch(i)) / 1.d-4, du_du_perp(i), (u_step - u(i)) / 1.d-4
      write(94,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)"), &
              u_perp(i), " ",&
              dpitch_du_perp(i), " ", &
              (pitch_step(i) - pitch(i)) / 1.d-4, " ", &
              du_du_perp(i), " ", &
              (u_step(i) - u(i)) / 1.d-4
   end do
  close(93)
  close(94)
  end if
  df_du_par = (df_du * du_du_par + dpitch_du_par * df_dpitch) * f
  df_du_perp = (df_du * du_du_perp + dpitch_du_perp * df_dpitch) * f
  do i = 1, int(m,8)
    if(u(i) >= ffp%u_max .or. u(i) <= ffp%u_min .or. pitch(i) <= ffp%pitch_min .or. pitch(i) >= ffp%pitch_max) then
      f(i) = exp( mu* (1.d0 - sqrt(1.d0 + u(i)**2))) * norm
      df_du_par(i) = -((mu * u_par(i))/sqrt(1.d0 + u(i)**2)) * f(i)
      df_du_perp(i) = -((mu * u_perp(i))/sqrt(1.d0 + u(i)**2)) * f(i)
    end if
  end do
end subroutine make_f_and_f_grad_along_line

end module mod_ecfm_refr_fp_dist_utils
