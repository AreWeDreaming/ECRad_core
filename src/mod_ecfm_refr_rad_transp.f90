! module mod_ecfm_refr_rad_transp
!        subroutine benchmark_abs_and_N, calculate_Irad


!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ecfm_refr_rad_transp

  use f90_kind
  use mod_ecfm_refr_types,         only: rad_diag_ch_mode_ray_freq_svec_type, &
                                         rad_diag_ch_mode_ray_freq_svec_extra_output_type
  implicit none

  public :: calculate_Trad, calculate_Trad_LSODE

  private :: benchmark_abs_and_N, evaluate_em_ab_single
   real(rkind), dimension(116) :: work_lsode
!    real(rkind), dimension(112) :: work_lsode
  integer(ikind), dimension(20) :: iwork_lsode
  integer(ikind)                :: j_LSODE, idiag_LSODE, ich_LSODE, imode_LSODE, iray_LSODE, ifreq_LSODE
  real(rkind)                   :: glob_ab, glob_ab_secondary
  type(rad_diag_ch_mode_ray_freq_svec_type) :: glob_svec
  type(rad_diag_ch_mode_ray_freq_svec_extra_output_type) :: glob_svec_extra_output

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

subroutine evaluate_em_ab_LSODE(svec, svec_extra_output, s, ds, omega, mode, model, eval_pol_coeff, x_launch, em, ab, pol_coeff)
! Wrapper function for em, ab for the integration of the radiation tranposrt equation with DLSODE
! Since a lot of if clasus are required here this routines cleans up the code below
use mod_ecfm_refr_types,          only: rad_diag_ch_mode_ray_freq_svec_type, rad_diag_ch_mode_ray_freq_svec_extra_output_type, &
                                      output_level, data_folder, Ich_name, &
                                      dstf, dstf_comp, rad, ffp, mode_cnt
use mod_ecfm_refr_em_Hu,                  only: calculate_em, simple_in_cutoff
use constants,                    only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,         only: abs_Albajar, abs_Albajar_fast, abs_Al_tor_abs, func_N_cold, func_rel_N
implicit none
type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout) :: svec
type(rad_diag_ch_mode_ray_freq_svec_extra_output_type), intent(out) :: svec_extra_output
real(rkind),                intent(in)    :: s, omega, ds
integer(ikind),             intent(in)    :: mode
character(*),               intent(in)    :: model
logical,                    intent(in)    :: eval_pol_coeff
real(rkind), dimension(:), intent(in)     :: x_launch
real(rkind),                intent(out)   :: em, ab
real(rkind),                intent(out)   :: pol_coeff
real(rkind)                               :: ab_secondary, dummy
!call benchmark_abs_and_N()
!stop "benchmarking"
em = 0.d0
ab = 0.d0
if(eval_pol_coeff) then
  pol_coeff = 0.d0
end if
! Interpolated svec
if(svec%Te > 1.d0 .and. svec%ne > 1.d15 .and. svec%ne < 1.d21) then
   ! Negative values are possible in the chi^2 calcuation since we use the nag spline in that instance, which may lead to ringing
  svec%N_cold = func_N_cold(omega, svec, mode)! note the convention for mode
  if(trim(model) == "primary") then
    if(output_level .or. (dstf /= "relamax" .and. dstf /= "Hu")) then
      if(output_level) svec_extra_output%N_cor = func_rel_N(omega, svec, mode)
      if(dstf_comp /= "TO") then
       if(eval_pol_coeff .and. mode_cnt == 2) then
          ! Take last point within the separatrix
          call abs_Albajar(svec, omega, mode, ds, ab, em, &
                           pol_coeff = pol_coeff, x_launch = x_launch)
       else
          call abs_Albajar(svec, omega, mode, ds, ab, em)
          !call abs_Albajar_fast(svec, omega, mode, ds2, ab)
          !em = ab * svec%Ibb
        end if
      else
        call calculate_em(svec, omega, em, dummy, ab)                          ! out    [W m^-3 sr^-1 Hz^-1]
      end if
    else
      if(dstf == "relamax") then
        call abs_Albajar_fast(svec, omega, mode, ds, ab)
        em = ab * svec%Ibb
      else if(dstf == "Hu") then
        call calculate_em(svec, omega, em, dummy, ab)                          ! out    [W m^-3 sr^-1 Hz^-1]
        ab =  em / svec%Ibb
      else
        print*, "dstf flag not set correctly - allowed is relamax or Hu"
        print*, "dstf flag is", dstf
        stop "mod_ecfm_rad_int calculate_Trad"
      end if
    end if
  else if(trim(model) == "secondary") then
    if ((dstf_comp == "Al" .or. dstf_comp == "Th") .and. dstf == "relamax") then
      if(mode /= -1) then
        call calculate_em(svec, omega, em, dummy, ab)
      end if
    else if( dstf_comp == "TB" .or. dstf_comp == "TO") then
      if(svec_extra_output%N_warm <= 0 .or. &
        svec_extra_output%N_warm /= svec_extra_output%N_warm) then
        svec_extra_output%N_warm = svec%N_cold ! do not reuse last, but start with cold
      end if
      if(eval_pol_coeff .and. mode_cnt == 2) then
        ! last point within the separatrix
        ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm, &
                    pol_coeff_secondary = pol_coeff, x_launch = x_launch)
      else
        ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm)
      end if
      call abs_Albajar(svec, omega, mode, ds, ab, em, &
                           pol_coeff = pol_coeff, x_launch = x_launch)
      ab = ab_secondary
      em = ab * svec%Ibb
    else if( dstf_comp == "O1") then
      svec_extra_output%N_warm = svec%N_cold
      ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm)
      em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            svec%Te / c0**2
    else
      em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            svec%Te / c0**2
    end if
!    if(ab_secondary > 1.d0) then
!      print*, em_secondary , ab_secondary
!    end if
    svec_extra_output%em_secondary = em
    svec_extra_output%ab_secondary = ab
  else
    print*,"Got model", model
    print*,"Expected either 'primary' or 'secondary'"
    stop "Bad arguments in evaluate_em_ab_single_LSODE in mod_ecfm_refr_rad_stranp.f90"
  end if
else
  return
end if
end subroutine evaluate_em_ab_LSODE

subroutine rad_transp(m, s, Trad, dIds)
! RHS of radiation transport equation to be integrated by DLSODE
  use f90_kind
  USE mod_ecfm_refr_types, only : ant, rad, output_level
  use mod_ecfm_refr_interpol,       only: spline_1d
  use constants,                    only: pi, e0, mass_e, eps0, c0
  implicit none
  integer          , intent(in)           :: m
  real(rkind),  intent(in)                :: s
  real(rkind),  dimension(m), intent(in)  :: Trad
  real(rkind),  dimension(m), intent(inout) :: dIds
  real(rkind)                             :: ds, omega, ab, em
  logical                                 :: eval_pol_coeff

  if(s > rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_max) then
    print*, "s in rad_transp >  s_max",s, rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_max
    call abort()
  end if
  glob_svec%s = s
  omega = ant%diag(idiag_LSODE)%ch(ich_LSODE)%freq(ifreq_LSODE) * 2.0 * pi
  if(j_LSODE > 1) then
    ds = s - rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec(j_LSODE - 1)%s
  else
    ds = s - rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_min
  end if
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%x, s, &
                  glob_svec%x_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%y, s, &
                 glob_svec%x_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%z, s, &
                 glob_svec%x_vec(3))
  glob_svec%R = sqrt(glob_svec%x_vec(1)**2 + glob_svec%x_vec(2)**2)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%rhop, s, &
                 glob_svec%rhop)
  if(glob_svec%rhop > 0.98 .and. glob_svec%rhop < 1.d0 .and. &
        .not. rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%use_external_pol_coeff) then
     eval_pol_coeff = .true.
  end if
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%theta, s, &
                 glob_svec%theta)
   glob_svec%cos_theta = cos(glob_svec%theta)
   glob_svec%sin_theta = sin(glob_svec%theta)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%freq_2X, s, &
                 glob_svec%freq_2X)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Te, s, &
                 glob_svec%Te)
  glob_svec%Ibb = glob_svec%Te * (omega / (2.d0 * pi))**2 * e0 / c0**2
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%ne, s, &
                 glob_svec%ne)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%v_g_perp, s, &
                 glob_svec%v_g_perp)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Nx, s, &
                 glob_svec%N_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Ny, s,  &
                 glob_svec%N_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Nz, s,  &
                 glob_svec%N_vec(3))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Bx, s, &
                 glob_svec%B_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%By, s,  &
                 glob_svec%B_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Bz, s,  &
                 glob_svec%B_vec(3))
  call evaluate_em_ab_LSODE( glob_svec, glob_svec_extra_output, &
                             s, ds, omega, rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%mode, "primary", eval_pol_coeff, &
                             ant%diag(idiag_LSODE)%ch(ich_LSODE)%ray_launch(iray_LSODE)%x_vec, em, ab, &
                             rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%pol_coeff)
  glob_ab = ab
  dIds(1) = em  * c0**2 / (omega**2 * e0) * 4.d0 * pi**2 - ab * Trad(1)
end subroutine rad_transp

subroutine Jac(nq,path, y_vec,ml,mu,dfdy,npowpd)
! Empty routine - gradients evaluated numerically by DLSODE
  use f90_kind
  implicit none
  integer,   intent(in)   :: nq,ml,mu,npowpd
  real(rkind), intent(in) :: path
  real(rkind), intent(inout) :: y_vec(:), dfdy(:,:)
   ! empty subroutine for the ODE-solver - if implicit scheme is used use finite differences
  return
end subroutine Jac

subroutine calculate_Trad_LSODE(rad_ray, idiag, ich, imode, ir, ifreq, mode, Trad, Trad_secondary, tau, tau_secondary, tau_array, tau_secondary_array, error, debug)
! In its current state this method does not work. Neither DLSODE mode 4 nor mode 5 produce reasonable results. Even if the resonance is marked as a
! difficult point, the predictor, corrector method of DLSODE is not capable of solving the radiation transport equation more efficiently than the RK4 method below.
! DO NOT USE THIS ROUTINE!!!!
use mod_ecfm_refr_types,        only: ant, rad_diag_ch_mode_ray_type, output_level, data_folder, Ich_name, straight, max_points_svec, &
                                      tau_array, tau_secondary_array, dstf, dstf_comp, rad, ffp, mode_cnt, tau_thick
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,           only: func_N_cold, func_rel_N
use mod_ecfm_refr_interpol,       only: spline_1d
implicit none
type(rad_diag_ch_mode_ray_type), intent(inout) :: rad_ray
integer(ikind),             intent(in)    :: idiag, ich, imode, ir, ifreq
integer(ikind),             intent(in)    :: mode
real(rkind),                intent(out)   :: Trad, Trad_secondary
real(rkind),                intent(out)   :: tau, tau_secondary
real(rkind), dimension(:), intent(inout)  :: tau_array, tau_secondary_array
integer(ikind),             intent(out)   :: error
logical,    intent(in), optional          :: debug
real(rkind)    :: s, ds
real(rkind), dimension(:), allocatable    :: Trad_vec, dIds
integer(ikind)                            :: nTrad, istate, n_s_important, i_important
character(120)               :: cur_filename
character(20)                :: ich_str
real(rkind)    :: h, h_large, h_small, freq_2X_step, freq, s_in, s_switch
real(rkind), dimension(20) :: s_important
logical         :: first, eval_pol_coeff, debug_internal
error = 0
idiag_LSODE = idiag; ich_LSODE = ich; imode_LSODE = imode; iray_LSODE = ir; ifreq_LSODE = ifreq
! initialization
Trad     =  0.d0
Trad_secondary =  0.d0
tau      =  0.d0
tau_secondary  =  0.d0
debug_internal = .false.
rad_ray%freq(ifreq)%svec(:)%rhop = -1.d0 ! to mark all entries in svec that are not used for the radiation transport
if(present(debug)) debug_internal = debug
if(.not. rad_ray%freq(ifreq)%use_external_pol_coeff) then
  rad_ray%freq(ifreq)%pol_coeff = 1.0
  rad_ray%freq(ifreq)%pol_coeff_secondary = 1.0
end if
first = .false.
tau_array(:) = 0.d0
if(output_level) then
  nTrad = 2
  tau_secondary_array(:) = 0.d0
else
  nTrad = 1
end if
allocate(Trad_vec(nTrad), dIds(nTrad))
Trad_vec(:) = 0.d0
freq =  ant%diag(idiag)%ch(ich)%freq(ifreq)
!-----------------------------------------------------------------
j_LSODE = 1
s = rad_ray%freq(ifreq)%svec_spline%s_min
if(dstf == "Re") ffp%restart_spline = .True.
if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE)%N_warm = 0.d0 ! force cold plasma eval
call rad_transp(nTrad, s, Trad_vec, dIds)
rad_ray%freq(ifreq)%svec(j_LSODE) = glob_svec
if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE) = glob_svec_extra_output
h_large = 5.d-2
h_small = 1.d-2
istate = 1
n_s_important = 3
s_important(1) = rad_ray%freq(ifreq)%svec_spline%s_center
s_important(2) = rad_ray%freq(ifreq_LSODE)%s_res
s_important(3) = rad_ray%freq(ifreq)%svec_spline%s_max
!print*, "s start", s_important(1:n_s_important -1)
!print*, "s end", s_important(2:n_s_important)
do while(any(s_important(1:n_s_important -1) > s_important(2:n_s_important)))
  do i_important = 1, n_s_important - 1
    if(s_important(i_important) > s_important(i_important + 1)) then
      s_switch = s_important(i_important)
      s_important(i_important) = s_important(i_important + 1)
      s_important(i_important + 1) = s_switch
    end if
  end do
end do
work_lsode(1) = s_important(1)
i_important = 1
do while(s < rad_ray%freq(ifreq_LSODE)%svec_spline%s_max)
  if(j_LSODE ==  max_points_svec) then
    print*, "Need more maximum points in svec"
    print*, "Please increase max_points_svec and recompile"
    call abort()
  end if
  if(abs(rad_ray%freq(ifreq)%svec(j_LSODE)%freq_2X / freq - 1.d0) <  0.02) then
    h = h_small
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
  else
    h = h_large
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
    call spline_1d(rad_ray%freq(ifreq)%svec_spline%freq_2X, s, freq_2X_step)
    if(abs(freq_2X_step / freq - 1.d0) <  0.02) h = h_small
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
  end if
  if(h == 0.d0) then
    print*, "h in rad transp last step was zero", s
    call abort
  end if
  s_in = s
  !print*, "s + h, s_important(i_important), i_important", s + h, s_important(i_important), i_important
  do while(s + h >= s_important(i_important) .and. i_important < n_s_important)
  ! Need do loop here in case s < s(center) < s(res) < s + h
    i_important = i_important + 1
    work_lsode(1) = s_important(i_important)
    !print*, "New s_crit, s_max", s_important(i_important), rad_ray%freq(ifreq)%svec_spline%s_max
  end do
  call dlsode(rad_transp, (/nTrad/) , Trad_vec, s, s + h,       &
                             1, (/1d-7/), (/1.d-7/), 5, istate, 1, &
                             work_lsode, size(work_lsode), iwork_lsode, size(iwork_lsode), Jac,10)
  if (istate /= 2) then
    print*, "h" ,h, "istate", istate, "Trad", Trad_vec
    print*, "LSODE failed to solve the radiation transport equation"
    call abort()
  end if
  ds = s - s_in
  !print*, "s, s_in, ds", s, s_in, ds
  if(glob_svec%s - rad_ray%freq(ifreq)%svec(j_LSODE)%s <= 0) cycle
  j_LSODE = j_LSODE + 1
  if(j_LSODE > 1) then
    tau_array(j_LSODE) = tau_array(j_LSODE) + glob_ab * ds
  else
    tau_array(j_LSODE) = glob_ab  * ds
  end if
  if(output_level) then
    if(j_LSODE > 1) then
      tau_secondary_array(j_LSODE) = tau_secondary_array(j_LSODE - 1)  + glob_ab_secondary  * ds
    else
      tau_secondary_array(j_LSODE) = glob_ab_secondary  * ds
    end if
  end if
  rad_ray%freq(ifreq)%svec(j_LSODE) = glob_svec
  if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE) = glob_svec_extra_output
  !print*, "j_LSODE, s", j_LSODE, rad_ray%freq(ifreq)%svec(j_LSODE)%s
  if(Trad_vec(1) /= Trad_vec(1) .or. Trad_vec(1) > max(glob_svec%Te * 10.0, 1.d6) &
                                            .or. tau_array(j_LSODE) < 0.d0) then
    ! Trad >> 1 MeV => much larer than saturation value of most ECEs
    print*, "Trad", Trad_vec(1)
    print*, "I_bb", glob_svec%Te / c0**2 * &
                    (freq**2 * e0)
    print*, "ab", glob_ab
    print*, "tau", tau_array(j_LSODE)
    print*, "optical depth of ds", glob_ab * ds
    print*, "length of ds", ds
    print*, "Te", glob_svec%Te
    print*, "ne", glob_svec%ne
    print*, "theta", glob_svec%theta / pi * 180
    print*, "freq_2X", glob_svec%freq_2X
    print*, "N_abs",  glob_svec%N_cold
    call abort
  end if
  if(output_level) then
    if(Trad_vec(2) /= Trad_vec(2) .or. Trad_vec(2) > max(glob_svec%Te * 10.0, 1.d6) &
                                            .or. tau_secondary_array(j_LSODE) < 0.d0) then
        ! Trad >> 1 MeV => much larer than saturation value of most ECEs
        print*, "Trad secondary", Trad_vec(2)
        print*, "I_bb", glob_svec%Te / c0**2 * &
                        (freq**2 * e0)
        print*, "ab secondary", glob_ab_secondary
        print*, "tau", tau_secondary_array(j_LSODE)
        print*, "optical depth of ds", glob_ab_secondary * ds
        print*, "length of ds", ds
        print*, "Te", glob_svec%Te
        print*, "ne", glob_svec%ne
        print*, "theta", glob_svec%theta / pi * 180
        print*, "freq_2X", glob_svec%freq_2X
        print*, "N_abs",  glob_svec%N_cold
        call abort
      end if
  end if
end do
if(j_LSODE == 1) then
  print*, "LSODE finished after just one step"
  print*, "s", rad_ray%freq(ifreq)%svec_spline%s_min, s, rad_ray%freq(ifreq)%svec_spline%s_max
  call abort()
end if
j_LSODE = j_LSODE - 1
rad_ray%freq(ifreq)%total_LOS_points = j_LSODE
tau = tau_array(j_LSODE)
Trad = Trad_vec(1)
if((ant%diag(idiag)%diag_name == "CTC" .or. ant%diag(idiag)%diag_name == "CTA" .or. &
   ant%diag(idiag)%diag_name == "IEC") .and. .not.  &
        rad_ray%freq(ifreq)%use_external_pol_coeff) then
  if(mode == - 1) then
    rad_ray%freq(ifreq)%pol_coeff = 0.d0
    rad_ray%freq(ifreq)%pol_coeff_secondary = 0.d0
  else
    rad_ray%freq(ifreq)%pol_coeff = 1.d0
    rad_ray%freq(ifreq)%pol_coeff_secondary = 1.d0
  end if
end if
if(output_level) then
  tau_secondary = tau_secondary_array(j_LSODE)
  Trad_secondary = Trad_vec(2)
  tau_secondary = tau_secondary_array(rad_ray%freq(ifreq)%total_LOS_points)
  rad_ray%freq(ifreq)%svec_extra_output(1:rad_ray%freq(ifreq)%total_LOS_points)%T = &
    exp(-(tau - tau_array(1:rad_ray%freq(ifreq)%total_LOS_points)))
  rad_ray%freq(ifreq)%svec_extra_output(1:rad_ray%freq(ifreq)%total_LOS_points)%T_secondary = &
    exp(-(tau_secondary - tau_secondary_array(1:rad_ray%freq(ifreq)%total_LOS_points)))
end if
end subroutine calculate_Trad_LSODE

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
call make_1d_spline(ab_spline, rad_ray_freq%total_LOS_points, rad_ray_freq%svec(1:rad_ray_freq%total_LOS_points)%s, \
                    rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%ab)
rad_ray_freq%svec_extra_output(rad_ray_freq%total_LOS_points)%T = 0.d0
do i = 1, rad_ray_freq%total_LOS_points - 1          ! integration over every second point on LOS
  call spline_1d_integrate(ab_spline, rad_ray_freq%svec(i)%s, rad_ray_freq%svec(rad_ray_freq%total_LOS_points)%s, rad_ray_freq%svec_extra_output(i)%T)
enddo !i = 1, rad_ray_freq%total_LOS_points
rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T = exp(-rad_ray_freq%svec_extra_output(1:rad_ray_freq%total_LOS_points)%T)
call deallocate_1d_spline(ab_spline)
end subroutine get_em_T_fast


end module mod_ecfm_refr_rad_transp
