! module mod_ecfm_refr_rad_transp
!        subroutine benchmark_abs_and_N, calculate_Irad


!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ecfm_refr_rad_transp

  use f90_kind
  implicit none

  public :: calculate_Trad

  private :: benchmark_abs_and_N, evaluate_em_ab_single


contains


subroutine benchmark_abs_and_N()
! produces profiles of both the absorption coefficient and the refractive index
use mod_ecfm_refr_types,         only: rad_diag_ch_mode_ray_freq_svec_type, straight
use f90_kind
use mod_ecfm_refr_abs_Al,             only: abs_Albajar, abs_Al_tor_abs, func_N_cold, func_rel_N, get_upper_limit_tau
use mod_ecfm_refr_em_Hu,                    only: calculate_em
use constants,                    only: pi, e0, mass_e, eps0, c0
#ifdef OMP
use omp_lib
#endif
implicit none
type(rad_diag_ch_mode_ray_freq_svec_type) :: svec
real(rkind)                               :: omega
integer(ikind)                            :: mode
real(rkind)                               :: c_abs_Alb, c_abs_warm_disp, c_abs_Hutch
real(rkind)                               :: N_cold, N_cor, N_gray
real(rkind), dimension(600)               :: Te_prof, ne_prof, freq_2X_prof
integer(ikind)                            :: i, m
real(rkind)                               :: ds, dummy_1, dummy_2, abs_crude_approx!, X, Y,
logical                                   :: warm_plasma
!Te_prof(1) = 200.d0
!Te_prof(200) = 25.d3
!do i  = 2, 199
!  Te_prof(i) = Te_prof(1) + (Te_prof(200) - Te_prof(1)) / 200.d0 * i
!end do
!ne_prof(1) = 0.5d19
!ne_prof(200) = 2.d20
!do i  = 2, 199
!  ne_prof(i) = ne_prof(1) + (ne_prof(200) - ne_prof(1)) / 200.d0 * i
!end do
straight = .true.
m = size(freq_2X_prof)
freq_2X_prof(1) = 65.d9
freq_2X_prof(m) = 220.d9
do i  = 2, m -2
  freq_2X_prof(i) = freq_2X_prof(1) + (freq_2X_prof(m) - freq_2X_prof(1)) / real(m,8) * i
end do
mode = 1
svec%theta = 85.d0/180.d0 * Pi
svec%cos_theta = cos(svec%theta)
svec%sin_theta = sin(svec%theta)
!svec%ne = 1.d19
svec%Te = 8.e3
omega = 140.d9 * 2.d0 * Pi
!X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
!Y = freq_2X_prof(1) * Pi / omega
warm_plasma = .true.
ds = 1.d0
svec%v_g_perp = svec%sin_theta
!print*, "X, Y", X, Y
svec%ne = 10.d19
open(81, file = "v_high_ne_abs_prof.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
svec%ne = 8.d19
open(81, file = "high_ne_abs_prof.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
open(81,  file = "low_ne_abs_prof.dat")
svec%ne = 6.d19
do i = 1, m
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  svec%freq_2X = freq_2X_prof(i)
  !X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  !print*, "X, Y", X, Y
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
omega = 110.d9 * 2.d0 * Pi
svec%ne = 6.d19
open(81, file = "v_high_ne_abs_prof_low_f.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
svec%ne = 4.d19
open(81, file = "high_ne_abs_prof_low_f.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
open(81,  file = "low_ne_abs_prof_low_f.dat")
svec%ne = 2.d19
do i = 1, m
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  svec%freq_2X = freq_2X_prof(i)
  !X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  !print*, "X, Y", X, Y
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp," ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
end subroutine benchmark_abs_and_N



subroutine evaluate_em_ab_single(rad_freq, j, omega, mode, ds2, eval_pol_coeff, x_launch, em, ab, em_secondary, ab_secondary, pol_coeff, pol_coeff_secondary)
! Wrapper function for em, ab, .... since a lot of if clauses regarding dstf and dstf_comp are required this subroutine cleans up the code below.
use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_type, &
                                      output_level, data_folder, Ich_name, &
                                      dstf, dstf_comp, ffp, SOL_Te, SOL_ne, ne_max
use mod_ecfm_refr_em_Hu,        only: calculate_em, simple_in_cutoff
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,           only: abs_Albajar, abs_Albajar_fast, abs_Al_tor_abs, func_N_cold, func_rel_N
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
real(rkind)                               :: dummy
!call benchmark_abs_and_N()
!stop "benchmarking"
em = 0.d0
ab = 0.d0
em_secondary = 0.d0
ab_secondary = 0.d0
if((rad_freq%svec(j)%Te > SOL_Te .and. rad_freq%svec(j)%ne > SOL_ne .and. rad_freq%svec(j)%ne < ne_max) .or. eval_pol_coeff) then
   ! Negative values are possible in the chi^2 calcuation since we use the nag spline in that instance, which may lead to ringing
  rad_freq%svec(j)%N_cold = func_N_cold(omega, rad_freq%svec(j), mode)! note the convention for mode
  if(output_level .or. (dstf /= "relamax" .and. dstf /= "Hu") .or. eval_pol_coeff) then
    if(output_level) rad_freq%svec_extra_output(j)%N_cor = func_rel_N(omega, rad_freq%svec(j), mode)
    if(dstf_comp /= "TO" .and. trim(dstf) /= "Hu_nbes" .and. trim(dstf) /= "Hu_bess" .and. trim(dstf) /= "Hu") then
      if(eval_pol_coeff) then
        ! Take last point within the separatrix
        if(trim(dstf_comp) == "Mx" .or. trim(dstf_comp) == "Al") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           pol_coeff = pol_coeff, c_abs_secondary = ab_secondary, x_launch = x_launch)!
          if(trim(dstf_comp) == "Al") then
            ab = ab_secondary
            em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
          end if
          pol_coeff_secondary = pol_coeff
        else if(trim(dstf) == "gene" .or. trim(dstf) == "gcomp") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           pol_coeff = pol_coeff, c_abs_secondary = ab_secondary, &
                           j_secondary = em_secondary, x_launch = x_launch)!
          pol_coeff_secondary = pol_coeff
        else if(trim(dstf_comp) == "Th") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           pol_coeff = pol_coeff, c_abs_secondary = ab_secondary, x_launch = x_launch)!
          pol_coeff_secondary = pol_coeff
        else
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           pol_coeff = pol_coeff, x_launch = x_launch)!
        end if
      else
        if(trim(dstf_comp) == "Mx" .or. trim(dstf_comp) == "Al") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           c_abs_secondary = ab_secondary)!
          if(trim(dstf_comp) == "Al") then
            ab = ab_secondary
            em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
          end if
        else if(trim(dstf) == "gene" .or. trim(dstf) == "gcomp") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           c_abs_secondary = ab_secondary, &
                           j_secondary = em_secondary)!
        else if(trim(dstf_comp) == "Th") then
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em, &
                           c_abs_secondary = ab_secondary)
        else
          call abs_Albajar(rad_freq%svec(j), omega, mode, ds2, ab, em)!
        end if
        !call abs_Albajar_fast(rad_freq%svec(j), omega, mode, ds2, ab)
        !em = ab * rad_freq%svec(j)%Ibb
      end if
      if(dstf_comp == "maxwell") em_secondary = ab_secondary * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
    else
      if(mode /= -1) then
        call calculate_em(rad_freq%svec(j), omega, em, em_secondary, ab)                        ! out    [W m^-3 sr^-1 Hz^-1]
        if(eval_pol_coeff) pol_coeff = 1.d0
      else
        if(eval_pol_coeff) pol_coeff = 0.d0
      end if
    end if
  else
    if(dstf == "relamax") then
      call abs_Albajar_fast(rad_freq%svec(j), omega, mode, ds2, ab)
      em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
    else if(dstf == "Hu") then
      call calculate_em(rad_freq%svec(j), omega, em, em_secondary, ab)                          ! out    [W m^-3 sr^-1 Hz^-1]
      ab =  em / ((omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2)
    else
      print*, "dstf flag not set correctly - allowed is relamax or Hu"
      print*, "dstf flag is", dstf
      stop "mod_ecfm_rad_int calculate_Trad"
    end if
  end if
  if(output_level) then
    if (dstf_comp == "Al") then
      if(mode /= -1) then
          call calculate_em(rad_freq%svec(j), omega, dummy, em_secondary, ab_secondary)                      ! out    [W m^-3 sr^-1 Hz^-1]
          ab_secondary = em_secondary / ((omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2)
        if(eval_pol_coeff) pol_coeff_secondary = 1.d0
      else
        if(eval_pol_coeff) pol_coeff_secondary = 0.d0
      end if
    else if((dstf_comp == "TB" .or. dstf_comp == "TO")) then
      if(rad_freq%svec_extra_output(j)%N_warm <= 0 .or. &
        rad_freq%svec_extra_output(j)%N_warm /= rad_freq%svec_extra_output(j)%N_warm) then
        rad_freq%svec_extra_output(j)%N_warm = rad_freq%svec(j)%N_cold ! do not reuse last, but start with cold
      end if
      if(eval_pol_coeff) then
        ! last point within the separatrix
        ab_secondary = abs_Al_tor_abs(rad_freq%svec(j), omega, mode, rad_freq%svec_extra_output(j)%N_warm, &
                    pol_coeff_secondary = pol_coeff_secondary, x_launch = x_launch)
      else
        ab_secondary = abs_Al_tor_abs(rad_freq%svec(j), omega, mode, rad_freq%svec_extra_output(j)%N_warm)
      end if
      em_secondary = ab_secondary * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
    else if( dstf_comp == "O1") then
      rad_freq%svec_extra_output(j)%N_warm = rad_freq%svec(j)%N_cold
      ab_secondary = abs_Al_tor_abs(rad_freq%svec(j), omega, mode, rad_freq%svec_extra_output(j)%N_warm)
      em_secondary = ab_secondary * (omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2
    else if(trim(dstf) == "Hu_nbes" .or. trim(dstf) == "Hu_bess") then
      ab_secondary = em_secondary / ((omega / ( 2.d0 * pi))**2 * e0 * &
            rad_freq%svec(j)%Te / c0**2)
      if(mode /= -1) then
        if(eval_pol_coeff) pol_coeff_secondary = 1.d0
      else
        if(eval_pol_coeff) pol_coeff_secondary = 0.d0
      end if
    else if(trim(dstf) /= "gene" .and. trim(dstf) /="gcomp" .and. (dstf_comp == "Th" .or. dstf_comp == "Mx" )) then
      em_secondary = ab_secondary * (omega / ( 2.d0 * pi))**2 * e0 * &
                    rad_freq%svec(j)%Te / c0**2
    end if
    rad_freq%svec_extra_output(j)%em = em
    rad_freq%svec_extra_output(j)%em_secondary = em_secondary
    rad_freq%svec_extra_output(j)%ab = ab
    rad_freq%svec_extra_output(j)%ab_secondary = ab_secondary
  end if
else
  return
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


use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_type, output_level, data_folder, Ich_name, straight, &
                                      dstf, dstf_comp, ffp, mode_cnt, tau_thick, SOL_Te, &
                                      static_grid, plasma_vac_boundary, spl_type_1d
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,       only: func_N_cold, func_rel_N
use mod_ecfm_refr_interpol,     only: make_1d_spline, spline_1d_integrate, deallocate_1d_spline
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
character(120)               :: cur_filename
character(20)                :: ich_str
real(rkind), dimension(3) :: em, ab, em_secondary, ab_secondary
real(rkind)    :: k1, k2, k3, k4, em_diff, dummy, N_cor, smallest_rhop_on_LOS, &
                  ab_diff,trans_coeff,trans_coeff_secondary,  omega, N_est, delta, Trad_0
integer(ikind) :: i, j, k, ich_tot, last_index_inside, i_smallest_rhop_on_LOS
logical        :: first, eval_pol_coeff, debug_internal
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
  if(.not. all(rad_ray_freq%svec(i + 1:i + 2)%plasma) .or. all(rad_ray_freq%svec(i + 1:i + 2)%Te < SOL_Te)) then
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
  if(output_level .and. mod(i,50) ==0 ) print*, i, rad_ray_freq%svec(i)%rhop, rad_ray_freq%svec(i)%Te, Trad
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
    if(Trad > SOL_Te .and. .not. (rad_ray_freq%max_points_svec_reached .or. static_grid) ) then
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


subroutine get_em_T_fast(rad_ray_freq, freq, x_vec_launch, mode, Trad)
! Routine still work in progress - DO NOT USE!
use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_type, output_level, data_folder, Ich_name, straight, &
                                      dstf, dstf_comp, ffp, mode_cnt, tau_thick, SOL_Te, &
                                      static_grid, plasma_vac_boundary, spl_type_1d
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,       only: func_N_cold, func_rel_N
use mod_ecfm_refr_interpol,     only: make_1d_spline, spline_1d_integrate, deallocate_1d_spline
implicit none
type(rad_diag_ch_mode_ray_freq_type), intent(inout) :: rad_ray_freq
real(rkind),                intent(in)    :: freq
real(rkind), dimension(3),  intent(in)    :: x_vec_launch
integer(ikind),             intent(in)    :: mode
real(rkind),                intent(in)    :: Trad
type(spl_type_1d)                         :: ab_spline
real(rkind), dimension(rad_ray_freq%total_LOS_points) :: s, ab
character(120)               :: cur_filename
character(20)                :: ich_str
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


end module mod_ecfm_refr_rad_transp
