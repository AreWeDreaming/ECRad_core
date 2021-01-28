! module mod_ECRad_rad_transp
!        subroutine benchmark_abs_and_N, calculate_Irad


!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ECRad_rad_transp

  use f90_kind
  implicit none

  public :: calculate_Trad

  private :: evaluate_em_ab_single


contains

subroutine evaluate_em_ab_single(rad_freq, j, omega, mode, ds2, eval_pol_coeff, x_launch, &
                                 em, ab, em_secondary, ab_secondary, pol_coeff, pol_coeff_secondary)
! Wrapper function for em, ab, .... since a lot of if clauses regarding dstf and dstf_comp are required this subroutine cleans up the code below.
use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_type, &
                                      output_level, dstf, ignore_Te, ignore_ne, ne_max
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ECRad_abs_Al,           only: abs_Albajar, abs_Albajar_fast, abs_Al_Fa_abs, func_N_cold, func_rel_N
implicit none
type(rad_diag_ch_mode_ray_freq_type), intent(inout) :: rad_freq
integer(ikind),             intent(in)    :: j
real(rkind),                intent(in)    :: omega, ds2
integer(ikind),             intent(in)    :: mode
logical,                    intent(in)    :: eval_pol_coeff
real(rkind), dimension(:), intent(in)     :: x_launch
real(rkind),                intent(out)   :: em, ab
real(rkind),                intent(out)   :: em_secondary, ab_secondary
real(rkind),                intent(out)   :: pol_coeff,  pol_coeff_secondary
!call benchmark_abs_and_N()
!stop "benchmarking"
em = 0.d0
ab = 0.d0
em_secondary = 0.d0
ab_secondary = 0.d0
! Negative values are possible in the chi^2 calcuation of IDA since we use the nag spline in that instance, which may lead to ringing
! Hence, we only evaluate the absorption coefficient for decent densities, Te > Te_sol or if the polarization coefficient is needed
if((rad_freq%svec(j)%Te <= ignore_Te .or. rad_freq%svec(j)%ne <= ignore_ne .or. rad_freq%svec(j)%ne > ne_max) .and. .not. eval_pol_coeff) return
rad_freq%svec(j)%N_cold = func_N_cold(omega, rad_freq%svec(j), mode)! note the convention for mode
if(output_level .or. eval_pol_coeff .or. trim(dstf) /= "Th") then
  if(output_level) rad_freq%svec_extra_output(j)%N_cor = func_rel_N(omega, rad_freq%svec(j), mode)
  if(eval_pol_coeff) then
    if(trim(dstf) == "Th") then
      call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                       pol_coeff = pol_coeff, x_launch = x_launch)!
    else  !non-thermal
      call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                       pol_coeff = pol_coeff, c_abs_secondary = ab_secondary, &
                       j_secondary = em_secondary, x_launch = x_launch)!
      pol_coeff_secondary = pol_coeff
    end if
  else
    if(trim(dstf) == "Th") then
      call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em)!
    else !non-thermal
      call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                       c_abs_secondary = ab_secondary, &
                       j_secondary = em_secondary)!
    end if
  end if
else
  call abs_Albajar_fast(rad_freq%svec(j), omega, mode, ds2, ab)
  em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
        rad_freq%svec(j)%Te / c0**2
end if
if(output_level .and. dstf == "Th") then
  if(rad_freq%svec_extra_output(j)%N_warm <= 0 .or. &
    rad_freq%svec_extra_output(j)%N_warm /= rad_freq%svec_extra_output(j)%N_warm) then
    rad_freq%svec_extra_output(j)%N_warm = rad_freq%svec(j)%N_cold ! do not reuse last, but start with cold
  end if
  if(eval_pol_coeff) then
    ! last point within the separatrix
    ab_secondary = abs_Al_Fa_abs(rad_freq%svec(j), omega, mode, rad_freq%svec_extra_output(j)%N_warm, &
                pol_coeff_secondary = pol_coeff_secondary, x_launch = x_launch)
  else
    ab_secondary = abs_Al_Fa_abs(rad_freq%svec(j), omega, mode, rad_freq%svec_extra_output(j)%N_warm)
  end if
  em_secondary = ab_secondary * (omega / ( 2.d0 * pi))**2 * e0 * &
        rad_freq%svec(j)%Te / c0**2
end if
if(output_level) then
  rad_freq%svec_extra_output(j)%em = em
  rad_freq%svec_extra_output(j)%em_secondary = em_secondary
  rad_freq%svec_extra_output(j)%ab = ab
  rad_freq%svec_extra_output(j)%ab_secondary = ab_secondary
end if
end subroutine evaluate_em_ab_single

subroutine calculate_Trad(rad_ray_freq, freq, x_vec_launch, mode, Trad, Trad_secondary, tau, tau_secondary, tau_array, tau_secondary_array, error, debug)
! Solves the radiation transport equation for Trad and tau on the svec grid with a RK4 solver.
! The routine checks itself for numerical instability by comparing Trad
! attained with Rk4 to Trad calculated with explicit Euler. This method provides some estimate or the error but the method is entirely empirical.
! If the empircically estimated error is beyond a threshold this routine will abort and the step size in svec will be halfed.
! Due to the purely empirical nature of the current method a replacement method is desireable.
! A good, theorertically motivated replacement that also works with the fixed grid is the doubling method described in
! "Comparing Error Estimators for Runge-Kutta Methods"
! By L. F. Shampine and H. A. Watts (1971 in MATHEMATICS OF COMPUTATION, VOLUME 25, NUMBER 115)
! Tau is calculated with the rapezoid method. It is desireable to replace this method by a more sophisticated technique.


use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_type, output_level, &
                                      dstf, ffp, mode_cnt, tau_thick, ignore_Te, &
                                      static_grid, plasma_vac_boundary, spl_type_1d
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ECRad_abs_Al,       only: func_N_cold, func_rel_N
use mod_ECRad_interpol,     only: make_1d_spline, spline_1d_integrate, deallocate_1d_spline
implicit none
type(rad_diag_ch_mode_ray_freq_type), intent(inout) :: rad_ray_freq
real(rkind),                intent(in)    :: freq
real(rkind), dimension(3),  intent(in)    :: x_vec_launch
integer(ikind),             intent(in)    :: mode
real(rkind),                intent(out)   :: Trad, Trad_secondary
real(rkind),                intent(out)   :: tau, tau_secondary
integer(ikind),             intent(out)   :: error
real(rkind), dimension(:),  intent(inout) :: tau_array, tau_secondary_array ! Used twice - first to store ab then to store tau
real(rkind), dimension(rad_ray_freq%total_LOS_points) :: s_cons ! s stored consequtively in the memory
type(spl_type_1d)                         :: tau_spline
logical,    intent(in), optional          :: debug
real(rkind)    :: ds, ds2
real(rkind), dimension(3) :: em, ab, em_secondary, ab_secondary
real(rkind)    :: k1, k2, k3, k4, smallest_rhop_on_LOS, omega, delta, Trad_0
integer(ikind) :: i, j, k, last_index_inside, i_smallest_rhop_on_LOS
logical        :: first, debug_internal
!call benchmark_abs_and_N()
!stop "benchmark"
error = 0
debug_internal = .false.
if(present(debug)) debug_internal = debug
if(.not. rad_ray_freq%use_external_pol_coeff .and. mode_cnt > 1) then
  rad_ray_freq%pol_coeff = -1.0
  if(output_level) rad_ray_freq%pol_coeff_secondary = -1.0
end if
first = .false.
!-----------------------------------------------------------------
omega = freq * 2 * Pi
j = 1; k = 1
ds2  = abs(rad_ray_freq%svec(j + 1)%s - rad_ray_freq%svec(j)%s)   ! ds for all but last point
if(dstf == "Re") ffp%restart_spline = .True.
if(output_level) rad_ray_freq%svec_extra_output(j)%N_warm = 0.d0 ! force cold plasma eval
call evaluate_em_ab_single(rad_ray_freq, j, omega, mode, ds2, .False., x_vec_launch, &
                           em(k), ab(k), em_secondary(k), ab_secondary(k), rad_ray_freq%pol_coeff, &
                           rad_ray_freq%pol_coeff_secondary)
em(1) = em(1) * c0**2 / (freq**2 * e0)
tau_array(1) = ab(k)
smallest_rhop_on_LOS = 2.d0
last_index_inside = -1
i_smallest_rhop_on_LOS = -1
Trad_0 = Trad
!if(mode > 0) print*, "Solving for X-mode"
!if(mode < 0) print*, "Solving for O-mode"
if(output_level) tau_secondary_array(1) = 0.5d0 * ab_secondary(k) * ds2
! loop over all indices
!----------------------
do i = 1, rad_ray_freq%total_LOS_points-2, 2           ! integration over every second point on LOS
  ! step size
  !----------
  if(.not. all(rad_ray_freq%svec(i + 1:i + 2)%plasma) .or. all(rad_ray_freq%svec(i + 1:i + 2)%Te <= ignore_Te)) then
    em(1) = 0.d0
    ab(1) = 0.d0
    tau_array(i + 1) = 0.d0
    tau_array(i + 2) = 0.d0
    if(output_level) then
      j = i + 1
      rad_ray_freq%svec_extra_output(j)%em = 0.d0
      rad_ray_freq%svec_extra_output(j)%ab = 0.d0
      rad_ray_freq%svec_extra_output(j)%em_secondary = 0.d0
      rad_ray_freq%svec_extra_output(j)%ab_secondary = 0.d0
      rad_ray_freq%svec_extra_output(j)%N_warm = 1.d0
      j = i + 2
      rad_ray_freq%svec_extra_output(j)%em = 0.d0
      rad_ray_freq%svec_extra_output(j)%ab = 0.d0
      rad_ray_freq%svec_extra_output(j)%em_secondary = 0.d0
      rad_ray_freq%svec_extra_output(j)%ab_secondary = 0.d0
      rad_ray_freq%svec_extra_output(j)%N_warm = 1.d0
      tau_secondary_array(i + 1) = 0.d0
      tau_secondary_array(i + 2) = 0.d0
    end if
    cycle
  end if
  !if(output_level .and. mod(i,50) ==0 ) print*, i, rad_ray_freq%svec(i)%rhop, rad_ray_freq%svec(i)%Te, Trad
  ds2  = abs(rad_ray_freq%svec(i+1)%s - rad_ray_freq%svec(i)%s)   ! ds for all but last point
  ds = ds2 * 2.d0
  !------------------------------------------------
  do k = 2, 3             ! 3 Runge-Kutta evaluations; the first comes from previous iteration
    j = i + k - 1         ! corresponding indices on svec
    if(rad_ray_freq%svec(j - 1)%rhop < plasma_vac_boundary) last_index_inside = j
    if(rad_ray_freq%svec(j - 1)%rhop < smallest_rhop_on_LOS .and. rad_ray_freq%svec(j - 1)%rhop > 0.d0) i_smallest_rhop_on_LOS = i
    if(output_level) rad_ray_freq%svec_extra_output(j)%N_warm = rad_ray_freq%svec_extra_output(j - 1)%N_warm
    !print*, j, rad_ray_freq%svec(j)%plasma, rad_ray_freq%svec(j)%rhop, rad_ray_freq%svec(j)%Te, rad_ray_freq%svec(j)%ne
    call evaluate_em_ab_single(rad_ray_freq, j, omega, mode, ds2, .False., x_vec_launch, &
                         em(k), ab(k), em_secondary(k), ab_secondary(k), rad_ray_freq%pol_coeff, &
                         rad_ray_freq%pol_coeff_secondary)
  end do
  em(2:3) = em(2:3) * c0**2 / (freq**2 * e0)
  k1   = em(1) - ab(1) * Trad !dN_abs_2_ds(2) / N_abs_2(2)
  k2   = em(2) - ab(2) * (Trad + k1*ds2)
  k3   = em(2) - ab(2) * (Trad + k2*ds2)
  k4   = em(3) - ab(3) * (Trad + k3*ds)
  if(rad_ray_freq%svec(j)%N_cold <= 0) then  ! Cut_pff
    Trad = Trad_0
    tau_array(1:i + 1) = 0.d0
    if(output_level) then
      rad_ray_freq%svec_extra_output(1:i + 1)%ab = 0.d0
      rad_ray_freq%svec_extra_output(1:i + 1)%em = 0.d0
    end if
!    print*, "found negative refractive index", rad_ray_freq%svec(j)%N_cold
!    print*, "i", i, rad_ray_freq%svec(i)%ne, rad_ray_freq%svec(i)%Te
!    print*, "i", i, rad_ray_freq%svec(i)%freq_2X, rad_ray_freq%svec(i)%theta / pi * 180.d0
  else if((ab(2) + ab(3)) * ds2 > tau_thick) then
      ! For strongly radiating plasmas we get in trouble with numerical stability.
      ! If, however, the optical depth of the last plasma slab is large we can fall back on just using
      ! the black body intensity - for tau = 3 the error relative error is about 5 %
      if(dstf /= "relamax" .and. dstf /= "Hu") then
        if(Trad > 1.d7 .or. Trad < 0) then
          error = 1
          return
        end if
      else
        Trad = rad_ray_freq%svec(j)%Te
      end if
  else
    if(Trad > ignore_Te .and. .not. (rad_ray_freq%max_points_svec_reached .or. static_grid) ) then
      ! If last step not optically think double check if current step size is adequate
      delta = Trad + ds2 * (em(2) - ab(2) * Trad)
      delta = delta + ds2 * (em(3) - ab(3) * delta)
      delta = abs((Trad + ds/6.d0 * (k1 + 2.d0*(k2 + k3) + k4)) / Trad  - delta / Trad)
      if(delta > 2.0) then
        if(output_level) then
          print*, "Optical depth of slab from last step", (ab(2) + ab(3)) * ds2
          print*, "Error between Rk4 and 2 times fwd. Euler = ", delta, " > 200% of Trad - switching to smaller step size"
          print*, "Abs coeff of last step", ab(1), ab(2), ab(3)
          print*, "Trad/s of last step", em(1), em(2), em(3)
          print*, "Te of last steps", rad_ray_freq%svec(i)%Te, rad_ray_freq%svec(i + 1)%Te, &
                                      rad_ray_freq%svec(i + 2)%Te
          print*, "ne of last steps", rad_ray_freq%svec(i)%ne, rad_ray_freq%svec(i + 1)%ne, &
                                      rad_ray_freq%svec(i + 2)%ne
          print*, "ds of last steps", rad_ray_freq%svec(i + 1)%s - rad_ray_freq%svec(i)%s, &
                                      rad_ray_freq%svec(i + 2)%s - rad_ray_freq%svec(i + 1)%s
        end if
        error = 1
        return
      end if
    end if
    if((Trad > max(rad_ray_freq%svec(j)%Te * 3, 1.d5) .or. Trad < 0) .and. &
       (ab(1) + ab(2)) / 2.0 * ds2 > tau_thick * 0.5d0 .and. &
       .not. rad_ray_freq%max_points_svec_reached .and. .not. static_grid) then
      if(dstf /= "relamax" .and. dstf /= "Hu") then
        error = 1
        return
      end if
    end if
    Trad = Trad + ds/6.d0 * (k1 + 2.d0*(k2 + k3) + k4)
  end if
  tau_array(i + 1) = ab(2)
  tau_array(i + 2) = ab(3)
  if(output_level) then
    em_secondary(2:3) = em_secondary(2:3) * c0**2 / (freq**2 * e0)
    rad_ray_freq%svec_extra_output(i: i + 2)%Trad = Trad
    k1   = em_secondary(1) - (ab_secondary(1)) *  Trad_secondary
    k2   = em_secondary(2) - (ab_secondary(2)) * (Trad_secondary + k1*ds2)
    k3   = em_secondary(2) - (ab_secondary(2)) * (Trad_secondary + k2*ds2)
    k4   = em_secondary(3) - (ab_secondary(3)) * (Trad_secondary + k3*ds )
    if(rad_ray_freq%svec(j)%N_cold <= 0) then  ! Cut_off
      Trad_secondary = Trad_0
      tau_secondary_array(1:i + 1) = 0.d0
      rad_ray_freq%svec_extra_output(1:i + 1)%ab_secondary = 0.d0
      rad_ray_freq%svec_extra_output(1:i + 1)%em_secondary = 0.d0
    else if((ab_secondary(2) + ab_secondary(3)) * ds2 > tau_thick) then
        ! For strongly radiating plasmas we get in trouble with numerical stability.
        ! If, however, the optical depth of the last plasma slab is large we can fall back on just using
        ! the black body intensity - for tau = 3 the error relative error is about 5 %
        Trad_secondary = rad_ray_freq%svec(j)%Te
    ! No reinterpolations for secondary model
    else if(Trad_secondary > max(rad_ray_freq%svec(j)%Te * 3, 1.d5) .or. Trad_secondary < 0.d0) then
      print*, "Warning Trad secondary was reset due to divergence or negative values"
      print*, "Optical depth of last slab was below tau_thick", (ab_secondary(1) + ab_secondary(2)) * ds* 0.5d0
      Trad_secondary = rad_ray_freq%svec(j)%Te
    else
      Trad_secondary = Trad_secondary + ds/6.d0 * (k1 + 2.d0*(k2 + k3) + k4)
    end if
    rad_ray_freq%svec_extra_output(i: i + 2)%Trad_secondary = Trad_secondary
  end if
  if(output_level) then
    tau_secondary_array(i + 1) = ab_secondary(2)
    tau_secondary_array(i + 2) = ab_secondary(3)
  end if
  if(debug_internal) then
    j = i
    if(j + 2 > rad_ray_freq%total_LOS_points) j = rad_ray_freq%total_LOS_points - 2
    print*, "-----------------", j, "-th step ---------------"
    print*, "Cold refr. index", rad_ray_freq%svec(j)%N_cold
    print*, "Trad", Trad
    print*, "tau", tau_array(i + 1), tau_array(i + 2)
    print*, "em",em
    print*, "ab",ab
    if(output_level) then
      print*, "N warm",  rad_ray_freq%svec_extra_output(j)%N_warm
      print*, "Trad secondary", Trad_secondary
      print*, "tau secondary", tau_secondary_array(i + 1), tau_secondary_array(i + 2)
      print*, "em_secondary",em_secondary
      print*, "ab_secondary",ab_secondary
    end if
    print*, "rhop", rad_ray_freq%svec(j)%rhop, rad_ray_freq%svec(j+1)%rhop, rad_ray_freq%svec(j+2)%rhop
    print*, "Te", rad_ray_freq%svec(j)%Te, rad_ray_freq%svec(j+1)%Te, rad_ray_freq%svec(j+2)%Te
    print*, "ne", rad_ray_freq%svec(j)%ne, rad_ray_freq%svec(j+1)%ne, rad_ray_freq%svec(j+2)%ne
    print*, "theta", rad_ray_freq%svec(j)%theta / pi * 180, rad_ray_freq%svec(j+1)%theta / pi * 180, rad_ray_freq%svec(j+2)%theta / pi * 180
    print*, "freq_2X", rad_ray_freq%svec(j)%freq_2X, rad_ray_freq%svec(j+1)%freq_2X, rad_ray_freq%svec(j+2)%freq_2X
  end if
  if(Trad /= Trad .or. Trad > max(rad_ray_freq%svec(j)%Te * 10.0, 1.d9)) then !/= Trad ! .or. tau_array(i + 2) < 0.d0
    ! Trad >> 1 MeV => much larer than saturation value of most ECEs
    j = i
    if(j + 2 > rad_ray_freq%total_LOS_points) j = rad_ray_freq%total_LOS_points - 2
    print*, "Trad", Trad
    print*, "I_bb", rad_ray_freq%svec(j)%Te / c0**2 * (freq**2 * e0)
    print*, "em",em
    print*, "ab",ab
    print*, "optical depth of ds", ab * ds2
    print*, "Te", rad_ray_freq%svec(j)%Te, rad_ray_freq%svec(j+1)%Te, rad_ray_freq%svec(j+2)%Te
    print*, "ne", rad_ray_freq%svec(j)%ne, rad_ray_freq%svec(j+1)%ne, rad_ray_freq%svec(j+2)%ne
    print*, "theta", rad_ray_freq%svec(j)%theta / pi * 180, rad_ray_freq%svec(j+1)%theta / pi * 180, rad_ray_freq%svec(j+2)%theta / pi * 180
    print*, "freq_2X", rad_ray_freq%svec(j)%freq_2X, rad_ray_freq%svec(j+1)%freq_2X, rad_ray_freq%svec(j+2)%freq_2X
    print*, "N_abs",  rad_ray_freq%svec(j)%N_cold, rad_ray_freq%svec(j+1)%N_cold, rad_ray_freq%svec(j+2)%N_cold
    print*, "j", j, j + 2
    error = -1
    return
  end if
  if(output_level) then
    if(Trad_secondary /= Trad_secondary .or. Trad_secondary > max(rad_ray_freq%svec(j)%Te * 10.0, 1.d9) .or. Trad_secondary < 0.d0) then
      print*, "Trad_secondary", Trad_secondary
      print*, "em_secondary",em_secondary
      print*, "ab_secondary",ab_secondary
      print*, "Te", rad_ray_freq%svec(j - 1)%Te
      print*, "ne", rad_ray_freq%svec(j - 1)%ne
      print*, "optical depth of ds", ab_secondary * ds2
      print*, "theta", rad_ray_freq%svec(j - 1)%theta / pi * 180
      print*, "freq_2X", rad_ray_freq%svec(j - 1)%freq_2X
      print*, "N_abs",  rad_ray_freq%svec(j - 1)%N_cold
      print*, "j", j
!      error = -1
!      return
    end if
  end if
  ! emissivity and absorption to be used in next iteration
  !-------------------------------------------------------
  em(1) = em(3)
  ab(1) = ab(3)
  if(output_level) then
    em_secondary(1) = em_secondary(3)
    ab_secondary(1) = ab_secondary(3)
  end if
!  print*, "third", i, rad_ray_freq%svec_extra_output(i:i + 1)%ab
enddo !i = 1, rad_ray_freq%total_LOS_points           ! integration over all points on LOS
!print*, "test", i, rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%ab
tau_array(rad_ray_freq%total_LOS_points) = ab(3)
!print*, "pol after radtransp", rad_ray_freq%pol_coeff
if(output_level) then
  tau_secondary_array(rad_ray_freq%total_LOS_points) = ab_secondary(3)
  rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%em = em(3)
  rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%em_secondary = em_secondary(3)
  rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%ab = ab(3)
  rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%ab_secondary = ab_secondary(3)
end if
s_cons(:) = rad_ray_freq%svec(1:rad_ray_freq%total_LOS_points)%s
call make_1d_spline(tau_spline, rad_ray_freq%total_LOS_points, s_cons, \
                    tau_array(1:rad_ray_freq%total_LOS_points))
call spline_1d_integrate(tau_spline, rad_ray_freq%svec(1)%s, rad_ray_freq%svec(rad_ray_freq%total_LOS_points)%s, tau)
if(output_level) then
  tau_array(:) = 0.d0 ! Not necessary but just to emphasize that from now on this will store tau and not ab
  do i = 1, rad_ray_freq%total_LOS_points          ! integration over every second point on LOS
    call spline_1d_integrate(tau_spline, rad_ray_freq%svec(1)%s, rad_ray_freq%svec(i)%s, tau_array(i))
  enddo !i = 1, rad_ray_freq%total_LOS_points
  call deallocate_1d_spline(tau_spline)
  call make_1d_spline(tau_spline, rad_ray_freq%total_LOS_points, s_cons, \
                      tau_secondary_array(1:rad_ray_freq%total_LOS_points))
  tau_secondary_array(:) = 0.d0  ! Not necessary but just to emphasize that from now on this will store tau_secondary and not ab_secondary
  do i = 1, rad_ray_freq%total_LOS_points          ! integration over every second point on LOS
    call spline_1d_integrate(tau_spline, rad_ray_freq%svec(1)%s, rad_ray_freq%svec(i)%s, tau_secondary_array(i))
  enddo !i = 1, rad_ray_freq%total_LOS_points
  tau_secondary = tau_secondary_array(rad_ray_freq%total_LOS_points)
  rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T = &
    exp(-(tau - tau_array(1:rad_ray_freq%total_LOS_points)))
  rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T_secondary = &
    exp(-(tau_secondary - tau_secondary_array(1:rad_ray_freq%total_LOS_points)))
end if
call deallocate_1d_spline(tau_spline)
if(last_index_inside == -1) then
! Did not make into the confined region of the plasma
  last_index_inside = i_smallest_rhop_on_LOS
end if
if( last_index_inside == -1) then
  if(output_level) print*, "WARNING! Not a single point of the LOS lies inside the given Flux matrix!!!!"
  if(output_level) print*, "Most likely the point of launch lies already inside an evanescent reason -> Check input"
  error = 0
  return
end if
if(.not. rad_ray_freq%use_external_pol_coeff .and. mode_cnt > 1) then
  call evaluate_em_ab_single(rad_ray_freq, last_index_inside, omega, mode, ds2, .True.,  x_vec_launch, &
                         em(1), ab(1), em_secondary(1), ab_secondary(1), rad_ray_freq%pol_coeff, &
                         rad_ray_freq%pol_coeff_secondary)
end if
if(debug_internal) then
  print*, "Static settings of radiation transport are:"
  print*, "measured frequency", freq
  print*, "Maximum Te along LOS", maxval(rad_ray_freq%svec(1:rad_ray_freq%total_LOS_points)%Te)
  print*, "Maximum ne along LOS", maxval(rad_ray_freq%svec(1:rad_ray_freq%total_LOS_points)%ne)
  print*, "Rho_pol value of last position in confined region", rad_ray_freq%svec(last_index_inside)%rhop
  print*, "Primary polarization coefficient", rad_ray_freq%pol_coeff
  if(output_level) print*, "Secondary polarization coefficient", rad_ray_freq%pol_coeff_secondary
end if
!print*, "tau, tau max", tau, maxval(tau_array(1:rad_ray_freq%total_LOS_points))
   !print*, "tau_scnd, tau_scnd max", tau_secondary, maxval(tau_secondary_array(1:rad_ray_freq%total_LOS_points))
!print*, "Calculated polarization filter", rad_ray_freq%pol_coeff, rad_ray_freq%use_external_pol_coeff
end subroutine calculate_Trad


subroutine get_em_T_fast(rad_ray_freq, freq, x_vec_launch, mode)
use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_type, spl_type_1d
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ECRad_abs_Al,       only: func_N_cold, func_rel_N
use mod_ECRad_interpol,     only: make_1d_spline, spline_1d_integrate, deallocate_1d_spline
implicit none
type(rad_diag_ch_mode_ray_freq_type), intent(inout) :: rad_ray_freq
real(rkind),                intent(in)    :: freq
real(rkind), dimension(3),  intent(in)    :: x_vec_launch
integer(ikind),             intent(in)    :: mode
type(spl_type_1d)                         :: ab_spline
real(rkind), dimension(rad_ray_freq%total_LOS_points) :: s, ab
real(rkind)    :: ab_dummy, em_dummy, omega, ds2
integer(ikind) :: i

!-----------------------------------------------------------------
omega = freq * 2 * Pi
!----------------------
do i = 1, rad_ray_freq%total_LOS_points          ! integration over every second point on LOS
  if(i > 1) then
    ds2 = rad_ray_freq%svec(i)%s - rad_ray_freq%svec(i - 1)%s
  else
    ds2 = rad_ray_freq%svec(i + 1)%s - rad_ray_freq%svec(i)%s
  end if
  call evaluate_em_ab_single(rad_ray_freq, i, omega, mode, ds2, .False., x_vec_launch, &
                             rad_ray_freq%svec_extra_output(i)%em, rad_ray_freq%svec_extra_output(i)%ab, em_dummy, ab_dummy, \
                             rad_ray_freq%pol_coeff, &
                             rad_ray_freq%pol_coeff_secondary)
enddo !i = 1, rad_ray_freq%total_LOS_points           ! integration over all points on LOS
s = rad_ray_freq%svec(1:rad_ray_freq%total_LOS_points)%s
ab = rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%ab
call make_1d_spline(ab_spline, rad_ray_freq%total_LOS_points, s, ab)
rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%T = 0.d0
do i = 1, rad_ray_freq%total_LOS_points - 1          ! integration over every second point on LOS
  call spline_1d_integrate(ab_spline, rad_ray_freq%svec(i)%s, rad_ray_freq%svec(rad_ray_freq%total_LOS_points)%s, rad_ray_freq%svec_extra_output(i)%T)
enddo !i = 1, rad_ray_freq%total_LOS_points
rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T = exp(-rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T)
call deallocate_1d_spline(ab_spline)
end subroutine get_em_T_fast


end module mod_ECRad_rad_transp
