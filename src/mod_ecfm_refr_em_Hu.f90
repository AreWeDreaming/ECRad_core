! module radiation
!        subroutine calculate_abs
!        subroutine calculate_em_tot
!        subroutine calculate_mode_frac
!        subroutine calculate_sf
!        subroutine calculate_sf_rrf
!        subroutine calculate_int_beta_nume
!        subroutine radiation_df_over_f


!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ecfm_refr_em_Hu

  use f90_kind

  implicit none
  public :: calculate_em, &
            radiation_gauss_init, &
            radiation_gauss_clean_up!, &
            !calculate_N
    real(rkind), dimension(:), allocatable :: Int_weights
    real(rkind), dimension(:), allocatable :: Int_absz
    real(rkind), dimension(:), allocatable :: Int_weights_many
    real(rkind), dimension(:), allocatable :: Int_absz_many
    integer(ikind)                         :: total_cnt, shortcut_cnt
  private :: BesselJ, &
             BesselJ_fast_1, &
             BesselJ_fast_2, &
             BesselJ_fast_3, & ! different syntax due to use of spline data
             calculate_em_tot,        & !
             calculate_mode_frac,     & !
             calculate_sf_rrf,        & ! analytical calculation from rrf
             calculate_int_u_gauss, &
             calculate_int_u_gauss_many, &
             calculate_int_u_gauss_arbitrary, &
             bessel_term, &
             Int_weights, &
             Int_absz, &
             Int_weights_many, &
             Int_absz_many, &
             total_cnt, &
             shortcut_cnt

contains

!  subroutine calculate_N(svec, omega, m, N, mode, old_model)
!    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
!    use constants,                  only: pi, e0 , eps0, mass_e
!    implicit None
!    type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout)  :: svec
!    real(rkind), intent(in)                           :: omega
!    integer(ikind), intent(in)                        :: m, mode
!    real(rkind), intent(out)                          :: N
!    logical, intent(in), optional                     :: old_model
!    real(rkind)                   :: omega_bar_2, alpha, beta,omega_p_c_2, omega_c, omega_p
!    real(rkind)                   :: f,rho
!    real(rkind)                   :: omega_LC, omega_RC
!    logical                       :: old_model_flag
!    old_model_flag = .False.
!    if(present(old_model)) old_model_flag = old_model
!    omega_c =  svec%freq_2X * PI
!    omega_p     = sqrt((svec%ne * e0**2.d0)/(eps0 * mass_e))
!    omega_bar_2 = (omega / omega_c)**2
!    alpha = omega_p**2 / omega**2
!    beta = omega_c**2 / omega**2
!    if(int(abs(mode)) /= 1) then
!      print*, "Only -1 (X-mode) and 1 (O-Mode) allowed for variable mode in"
!      print*, "the routine abs_non_therm_N"
!      stop "abs_non_therm_N in mod_ecfm_abs_Albajar.f90"
!    end if
!    omega_RC = sqrt(0.25 * omega_c**2 + omega_p**2)
!    omega_LC = omega_RC - 0.5 * omega_c
!    omega_RC = omega_RC + 0.5 * omega_c
!    omega_p_c_2 = omega_p / omega_c
!    if(alpha > 1.0) then
!      N = 0.0
!    else if(mode == 1) then
!      if((omega**2 - omega_LC**2)*(omega**2 - omega_RC**2) < 0.0) then
!        N = 0.0
!        return
!      end if
!    end if
!    rho =  svec%sin_theta**4 + 4.d0/real(m,8)**2 *(real(m,8)**2 - omega_p_c_2)**2 * svec%cos_theta**2
!    if(rho < 0.d0) then
!      print*, "root of negative number avoided: variable rho"
!      stop "calculate_N in mod_ecfm_radiation.f90"
!    end if
!    rho = sqrt(rho)
!    f = (2.d0 * (reaL(m,8)**2 - omega_p_c_2))/ (2.d0 * (real(m,8)**2 - omega_p_c_2) - (svec%sin_theta**2 + real(mode,8) * rho))
!    N = 1.d0 - omega_p_c_2/real(m,8)**2 * f
!    if(N < 0.d0) then
!      N = 0.0
!      return
!    end if
!    N = sqrt(N)
!    if(old_model_flag) N =  simple_in_cutoff(N,omega / (2 * Pi), svec%ne)
!  end subroutine calculate_N

  function simple_in_cutoff(N_abs,freq, ne)
  use mod_ecfm_refr_types,        only: old_cutoff
  implicit none
  real(rkind), intent(in)       :: N_abs, freq, ne
  logical                       :: simple_in_cutoff
  if(old_cutoff) then
    simple_in_cutoff = (0.620221303d-2 * freq**2 < ne)
    !print*, "cutoff", .620221303d-2 * freq**2, simple_in_cutoff
  else
    simple_in_cutoff = (N_abs <= 0.d0)
  end if
  end function simple_in_cutoff

  subroutine calculate_em(svec,omega, em, em_secondary, abs)
    ! S. Denk 4. 2013
    use mod_ecfm_refr_types,         only: ant, rad, non_maxwellian, &
                                           rad_diag_ch_mode_ray_freq_svec_type, &
                                           bi_max,drift_m, Spitzer, multi_slope, output_level, old_cutoff, &
                                           spl_type_2d, non_therm_params_type, dstf, dstf_comp
    use constants,                   only: pi, e0, mass_e, eps0, c0
    !use mod_ecfm_refr_fp_dist_utils, only: interpolate_f_rhop, export_Te_rhop
    !use ecfm_non_therm_abs,         only: abs_Albajar, abs_non_therm_tor_abs
    use mod_ecfm_radiation_dist,             only: radiation_dist_f_norm, prepare_dist
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout)  :: svec
    real(rkind), intent(in)                           :: omega
    real(rkind),    intent(out)                       :: em, em_secondary, abs                  ! [W m^-3 sr^-1 Hz^-1]
    character(30)                                     :: secondary_flag
    real(rkind)                                       :: freq, em_tot, em_tot_secondary, eta_2X, eta_2O, sf, sf_abs , &
                                                         sf_secondary,const, ne_frac, f_rel,f_rel_sq, mu, &
                                                         int_u, int_u_abs, f_norm,a, N, temp_norm , &
                                                         em_omega, abs_omega, em_secondary_omega, int_u_second
    integer(ikind)                                    ::  irhop_TDiff, irhop, m, l_omega
    type(spl_type_2d) :: f_spl
    type(non_therm_params_type)        :: dist_params
    em = 0.d0
    abs = 0.d0
    em_secondary = 0.d0
    if(svec%Te < 1.d0) return ! very low Te => absorption can be ignored
    if(svec%ne < 1.e17) return
    if(svec%N_cold <= 0.d0 .or. svec%N_cold > 1.d0) return
    freq = omega / (2.d0 * Pi)
    f_rel = freq /  svec%freq_2X
    if(dstf /= "relamax" .and. dstf /=  "Hu_nbes" .and.  dstf /= "Hu_bess") call prepare_dist(svec,Int_weights_many, Int_absz_many, f_spl, dist_params)
    do m = 2, 2
      sf = 0.d0
      sf_secondary = 0.d0
      sf_abs = 0.d0
      eta_2X = 0.d0
      eta_2O = 0.d0
      em_tot = 0.d0
      em_tot_secondary = 0.d0
      f_rel_sq = f_rel**2.d0
      const = (mass_e * c0**2.d0) / (2.d0 * e0 * svec%Te)
      mu = 2 * const
! Compare em/abs
      if(svec%cos_theta == 0.0 .or. svec%sin_theta == 0.0) then
          print*, "svec",svec
          print*, f_rel, svec%theta, svec%sin_theta, svec%cos_theta
          print*, "Sub_calculate_em"
          stop "Calculation not possible for theta = 90 deg. or theta == 0 deg."
      end if
      sf_secondary = 0.d0
      if(output_level) then
        if(trim(dstf_comp) == "Mx") then !"relamax"!! "maxwell"
          if(m == 2) then
            call calculate_sf_rrf(svec%Te, mu / 2.d0, svec%theta,       &
                  svec%cos_theta, svec%sin_theta, svec%freq_2X,f_rel, f_rel_sq, int_u_second)
!            print*, "doing the fischer", int_u_second
            call calculate_em_tot_secondary(svec%ne, &
              svec%cos_theta, &
              svec%sin_theta, &
              svec%freq_2X, const,m , em_tot_secondary)
            sf_secondary = int_u_second / (const**2 )
           else
            sf_secondary = 0.d0
          end if
        else if( trim(dstf_comp) == "NB") then
          call calculate_int_u_gauss_no_bessel(svec,mu,f_rel,f_rel_sq, int_u_second, int_u_abs,m,svec%N_cold)!
          call calculate_em_tot_secondary(svec%ne, &!
                    svec%cos_theta, &
                    svec%sin_theta, &
                    svec%freq_2X, const,m , em_tot_secondary)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf_secondary = int_u_second * pi/ freq * f_norm
        else if( trim(dstf_comp) == "Al") then
          call calculate_em_tot(svec%ne, &
            svec%cos_theta, &
            svec%sin_theta, &
            svec%freq_2X, const, m, em_tot_secondary)
          call calculate_int_u_gauss_refr(svec, mu,f_rel,f_rel_sq, int_u_second, int_u_abs,m,svec%N_cold)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf_secondary = int_u_second * pi/ freq * f_norm
        else
          print*, "Unknown dstf_comp:", dstf_comp
          stop "Critical ERROR in calcualte_em"
        end if
        sf = 0.0
        sf_abs = 0.0
        if( dstf == "Hu_nbes") then
          call calculate_int_u_gauss_no_bessel(svec,mu,f_rel,f_rel_sq, int_u, int_u_abs,m,svec%N_cold)!
!          print*, "doing no bessels whatsoever", int_u
          call calculate_em_tot_secondary(svec%ne, &!
            svec%cos_theta, &
            svec%sin_theta, &
            svec%freq_2X, const,m , em_tot)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf = int_u * pi/ freq * f_norm
          sf_abs =  - int_u_abs * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
        else if(dstf == "Hu_bess") then
          call calculate_int_u_gauss(svec, mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf = int_u * pi/ freq * f_norm
          sf_abs =  - int_u_abs * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
          call calculate_em_tot(svec%ne, & !
                                svec%cos_theta, &
                                svec%sin_theta, &
                                svec%freq_2X, const,m , em_tot)
        else
          call calculate_int_u_gauss_arbitrary(svec, dist_params,mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)
          call calculate_em_tot(svec%ne, &
                                svec%cos_theta, &
                                svec%sin_theta, &
                                svec%freq_2X, const,m , em_tot)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = radiation_dist_f_norm(svec%Te, dstf)
          sf = int_u * pi/ freq * f_norm
          sf_abs =  int_u * mu * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
        end if
      else
        call calculate_int_u_gauss(svec,mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)!_arbitrary
        a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
        f_norm = a * (sqrt(mu / (2 * pi))**3)
        sf = int_u * pi/ freq * f_norm
        call calculate_em_tot(svec%ne, &
                              svec%cos_theta, &
                              svec%sin_theta, &
                              svec%freq_2X, const,m , &
                              em_tot_secondary)
      end if
      call calculate_mode_frac(svec%cos_theta, svec%sin_theta, 0.5d0 / f_rel , eta_2X, eta_2O)
      em = em + eta_2X * sf * em_tot
      abs = abs + eta_2X * sf_abs * em_tot
      em_secondary = em_secondary + em_tot_secondary * eta_2X * sf_secondary
    end do
    if( em /= em  .or. em < 0.d0 .or. abs /= abs .or. abs < 0.d0) then !c_abs < 0.d0 .or.
        print*, "rhop", svec%rhop
        print*, "Te",svec%Te
        print*, "ne",svec%ne
        print*, "freq", freq
        print*, "freq 2",  svec%freq_2X
        print*, "eta_2X",eta_2X
        print*, "em_tot", em_tot
        print*, "int_u_second", int_u_second
        print*, "int_u", int_u
        print*, "sf", sf
        print*, "sf_abs", sf_abs
        print*, "em_tot_secondary",  em_tot_secondary
        print*, "sf_secondary", sf_secondary
        print*, "em", em
        print*, "abs", abs
        print*, "em_secondary", em_secondary
        stop "Nan in calculate_em"
      end if
  end subroutine calculate_em

  !*******************************************************************************

  subroutine calculate_em_tot(ne, cos_theta, sin_theta, freq_2X, const, m, em_tot)

    use constants,                  only: pi, e0, mass_e, eps0, c0

    implicit none

    real(rkind), intent(in)  :: ne, cos_theta, sin_theta, freq_2X, const
    integer(ikind), intent(in) :: m
    real(rkind), intent(out) :: em_tot
    real(rkind) ::  freq_m
    freq_m = freq_2X / 2 * m
    em_tot = (e0 * freq_m)**2.d0 * &
      ne / (eps0 * c0 )
  end subroutine calculate_em_tot

  subroutine calculate_em_tot_secondary(ne, cos_theta, sin_theta, freq_2X, const, m, em_tot)

    use constants,                  only: pi, e0, mass_e, eps0, c0

    implicit none

    real(rkind), intent(in)  :: ne, cos_theta, sin_theta, freq_2X, const
    integer(ikind), intent(in) :: m ! order of the harmonic
    real(rkind), intent(out) :: em_tot
    real(rkind) ::  freq_m
    integer(ikind) :: i, fac
    fac = 1
    do i = 1,m - 1
      fac = fac * i
    end do
    freq_m = freq_2X / 2 * m
    em_tot = (m**(2*(m-1))/fac**2) / (2.d0**(2 * m))* sin_theta**(2*(m-1)) * (cos_theta**2.d0  + 1.d0)
    em_tot = em_tot * (e0 * freq_m)**2.d0 * &
      ne / (eps0 * c0 )
  end subroutine calculate_em_tot_secondary

  !*******************************************************************************

  subroutine calculate_mode_frac(cos_theta, sin_theta, Y, eta_2X, eta_2O)

    implicit none

    real(rkind), intent(in)  :: cos_theta, sin_theta, Y
    real(rkind), intent(out) :: eta_2X, eta_2O
    real(rkind)  :: a1, a2, a3

    a1 = sin_theta**4.d0 / (4.d0 / Y)
    a2 = cos_theta**2.d0
    a3 = (a1 + a2) / ( (a2+1.d0) * sqrt(a2+a1 * Y) )

    eta_2X = 0.5d0 + a3
    eta_2O = 0.5d0 - a3

  end subroutine calculate_mode_frac

  !*******************************************************************************


  !*******************************************************************************

  subroutine calculate_sf_rrf(Te, const, theta, cos_theta, &
    sin_theta, freq_2X, f_rel, f_rel_sq, sf)

    use mod_ecfm_refr_types,        only: ant
    use constants,                  only: pi, e0, mass_e, c0, sPi

    implicit none

    real(rkind),    intent(in)  :: Te, const
    real(rkind),    intent(in)  :: theta, cos_theta, sin_theta
    real(rkind),    intent(in)  :: freq_2X, f_rel, f_rel_sq
    real(rkind),    intent(out) :: sf

    real(rkind)     :: eps, eta, xh1, xh2, prefact, t1, t2, sin_theta_2, num_eps, freq_new
    real(rkind)     :: beta_par_1, beta_par_2, alpha1, alpha2
    real(rkind)     :: e1, e2, sa1, sa2, ca1, ca2, F1, F2

    real(rkind)     :: sf2
    sin_theta_2 = sin_theta**2
    eps = 1.d0 / (const * f_rel_sq)                    ! a small quantity
    eta = 1.d0 + f_rel_sq * cos_theta**2               ! 1 + (freq/freq_2X * cos(theta))**2
    if(1.d0 - f_rel_sq * sin_theta_2 < 0.d0) then
      sf = 0.d0
      return
    end if
    t1 = f_rel_sq * cos_theta                          ! (freq/freq_2X)**2 * cos(theta)
    t2 = sqrt(1.d0 - f_rel_sq * sin_theta_2)           ! sqrt(1-(freq/freq_2X * sin(theta))**2)
    beta_par_1 = (t1 - t2) / eta                       ! lower integration limit
    beta_par_2 = (t1 + t2) / eta
    num_eps = abs(beta_par_2 - beta_par_1) * 1.d-9
    beta_par_1 = beta_par_1 + num_eps
    beta_par_2 = beta_par_2 - num_eps
    sa1 = 1.d0 - beta_par_1 * cos_theta                ! uses beta_par_1
    sa2 = 1.d0 - beta_par_2 * cos_theta                ! uses beta_par_2
    alpha1 = sa1**2                                    ! uses beta_par_1
    alpha2 = sa2**2                                    ! uses beta_par_2
    e1  = exp(-const*(1.d0 - f_rel_sq * alpha1))
    e2  = exp(-const*(1.d0 - f_rel_sq * alpha2))
    ca1 = sqrt(const*f_rel_sq*alpha1)
    ca2 = sqrt(const*f_rel_sq*alpha2)
    xh1 = eta**2 * eps
    xh2 = (3.d0*eta - 2.d0*sin_theta_2/eps) * sqrt(eps)

    call dawson_integral(ca1, F1)
    call dawson_integral(ca2, F2)
    freq_new = f_rel * freq_2X
    prefact = const**1.5d0 / (sPi * freq_new * f_rel_sq * cos_theta**5)

    sf = prefact * ( e1 * (xh1 - sin_theta_2 / sa1 - xh2 * F1) &
      -e2 * (xh1 - sin_theta_2 / sa2 - xh2 * F2) )
    if(sf < 0.d0) sf = 0.d0
  end subroutine calculate_sf_rrf

  !*******************************************************************************

  subroutine dawson_integral(x, daw)
    ! Dawson integral: daw = exp(-x**2) int_0^x exp(u**2) du

    use nr_special_funcs,           only: dawson_s

    implicit none

    real(rkind), intent(in)  :: x
    real(rkind), intent(out) :: daw   ! value of Dawson integral
    real(rkind)  :: y

    if (x >= 4.8d0) then
      y   = 1.d0/(2.d0*x**2)
      daw = (1.d0 + y*(1.d0 + 3.d0*y*(1.d0 + 5.d0*y*(1.d0 + 7.d0*y*(1.d0 + 9.d0*y*(1.d0 + 11.d0*y)))))) / (2.d0*x)
    else
      daw = dawson_s(x)
    endif

  !stop 'subroutine dawson_integral'
  end subroutine dawson_integral

  !*******************************************************************************

subroutine radiation_gauss_init(N_step_sf_int)
#ifdef NAG
    use nag_quad_util,              only: nag_quad_gs_wt_absc
    USE nag_error_handling
#endif
    use quadrature,                 only: cdgqf
    implicit none
    integer(ikind), intent(in)     :: N_step_sf_int
#ifdef NAG
    type(nag_error)                :: error
#endif
    integer(ikind)                 :: i
    real(rkind)                    :: h
    real(rkind), dimension(:), allocatable :: Int_weights_check, Int_absz_check
    allocate( Int_weights(N_step_sf_int),Int_absz(N_step_sf_int))
    allocate( Int_weights_many(64),Int_absz_many(64))
    call cdgqf( int(N_step_sf_int,kind=4), int(1,kind=4), 0.d0, 0.d0, Int_absz, Int_weights)
    call cdgqf( int(64,kind=4), int(1,kind=4), 0.d0, 0.d0, Int_absz_many, Int_weights_many)
#ifdef NAG
    allocate( Int_weights_check(N_step_sf_int), Int_absz_check(N_step_sf_int))
    call nag_quad_gs_wt_absc( 0, -1.d0, 1.d0, Int_weights_check, Int_absz_check)
    if(sum((Int_weights - Int_weights_check)**2) > 1.d-5) then
      print*, "Weights deviate by more than 1.d-10"
      do i = 1, N_step_sf_int
        print*, Int_weights(i), Int_weights_check(i)
      end do
      call abort()
    end if
    if(sum((Int_absz - Int_absz_check)**2) > 1.d-5) then
      print*, "Abszissae deviate by more than 1.d-10"
      do i = 1, N_step_sf_int
        print*, Int_absz(i), Int_absz_check(i)
      end do
      call abort()
    end if
    deallocate(Int_weights_check, Int_absz_check)
    !call nag_quad_gs_wt_absc( 0, -1.d0, 1.d0, Int_weights, Int_abszz)
    if (error%level >= 1) print*, error%msg
    if (error%level >= 1) stop "Failed to setup "
#endif
    total_cnt = 0
    shortcut_cnt = 0
  end subroutine radiation_gauss_init

  subroutine radiation_gauss_clean_up()
    use mod_ecfm_refr_types,        only: output_level
  implicit none

    if(allocated(Int_weights)) deallocate(Int_weights)
    if(allocated(Int_absz)) deallocate(Int_absz)
    if(allocated(Int_weights_many)) deallocate(Int_weights_many)
    if(allocated(Int_absz_many)) deallocate(Int_absz_many)
    if(output_level) print*, shortcut_cnt, " of ", total_cnt, " were calculated quickly"
  end subroutine radiation_gauss_clean_up




   subroutine calculate_int_u_gauss(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((real(m,8)/ (f_rel *2.0 ))**2  < 1.e0 - svec%cos_theta**2) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2
    cos_theta_2 = (cos_theta)**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta  + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 *  f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_JN(m, zeta(k)/ sin_theta))**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
      (Bessel_JN(m - 1, zeta(k)) - Bessel_JN(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * exp( -mu / 2.d0 * ((u_par(k) / gamma(k))**2 + u_perp(k)**2 / gamma(k)**2))* Int_weights(k)
      !* exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights(k)
      !if(stop) print*, u_par(k), int_u_add(k)
      !if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then

        !print*,"trub",
        !k_int = k
       ! stop = .true.
        !goto 1
        !return

      !end if
    enddo
    !if(stop) stop "done"
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
    end do
    int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss

  subroutine calculate_int_u_gauss_no_bessel(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2! * N**2
    cos_theta_2 = (cos_theta)**2! * N**2
    if(m_omega_bar_2 - sin_theta_2 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta  + m_omega_bar !* N
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      int_u_add(k) = u_perp(k)**(2*m) / gamma(k)**(2*m - 2)
      int_u_add(k) = int_u_add(k) * du_dalpha * exp( -mu / 2.d0 * ((u_par(k) / gamma(k))**2 + u_perp(k)**2 / gamma(k)**2)) * Int_weights(k)
      !exp( mu* (1.d0 - gamma(k))) * Int_weights(k)
      !if(stop) print*, u_par(k), int_u_add(k)
!      if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then
!
!        !print*,"trub",
!        !k_int = k
!        stop = .true.
!        !goto 1
!        !return
!
!      end if
    enddo
    !if(stop) stop "done"
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
    end do
    int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_no_bessel


  subroutine calculate_int_u_gauss_arbitrary(svec, dist_params, mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, k_int, non_therm_params_type
    use mod_ecfm_radiation_dist,    only: radiation_dist_f_u, radiation_dist_Rf,radiation_dist_f_norm
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    type(non_therm_params_type), intent(in) :: dist_params
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2, det!,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights_many))     :: u_par,u_perp, gamma, zeta, int_u_add, int_u_abs_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, N_par, N_perp
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    N_par = svec%N_cold * cos_theta !
    N_perp = svec%N_cold * sin_theta !
    m_omega_bar =  real(m,8) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    if(m_omega_bar_2  + N_par**2  - 1.e0 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    !sin_theta_2 = (sin_theta)**2!  * N**2
    !cos_theta_2 = (cos_theta)**2! * N**2
    !det = sqrt(m_omega_bar_2 - sin_theta_2)
    !a = (m_omega_bar*cos_theta - det) / sin_theta_2
    !b = (m_omega_bar*cos_theta + det) / sin_theta_2
!    if(m_omega_bar**2  + N_par**2  - 1.e0 < 0.d0) then
!      Int_u = 0.d0
!      int_u_abs = 0.d0
!      return
!    end if
    a =  (m_omega_bar * N_par - sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    b = (m_omega_bar * N_par + sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    !1 continue
    do k = 1, size(Int_weights_many)
      u_par(k) = Int_absz_many(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * N_par + m_omega_bar
      !gamma(k) = u_par(k) * cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)

      ! Hutchinson
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_Jn(m, zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
      (Bessel_JN(m - 1, zeta(k)) - Bessel_Jn(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * gamma(k)**2 * Int_weights_many(k)
      !if(stop)  print*,  u_par(k),  radiation_dist_Rf(u_par(k), u_perp(k), gamma(k), m_omega_bar, cos_theta, mu, svec)
      !if(u_perp(k) /= u_perp(k) .or. int_u_add(k) /= int_u_add(k)) then
      !  print*, m_omega_bar**2  + N_par**2  - 1.e0
      !stop "nan"
      !end if
      !if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then

        !print*,"trub",
        !k_int = k
        !stop = .true.
        !goto 1
        !return
      !end if
      !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
        !print*, "u_par",u_par(k)
        !print*, "u_perp", u_perp(k)
        !print*, "f", radiation_dist_f_u(u_par(k), u_perp(k), gamma(k), mu, svec)
        !print*, "int_u",int_u_add(k)
      !  print*,u_par(k),radiation_dist_f_u(u_par(k), u_perp(k), gamma(k), mu, svec)
      !end if
    enddo
    int_u_abs_add = int_u_add * radiation_dist_Rf(u_par, u_perp, gamma, m_omega_bar, cos_theta, mu, svec, dist_params)
    int_u_add = int_u_add * radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params)
    !if(stop) stop "done"
    !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
    !  stop "Sense?"
    !end if
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights_many)
      Int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs  +  int_u_abs_add(k)
    end do
    !int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_arbitrary

  subroutine calculate_int_u_gauss_refr(svec, mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, k_int, non_therm_params_type
    use mod_ecfm_radiation_dist,    only: radiation_dist_f_u, radiation_dist_Rf,radiation_dist_f_norm
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2, det!,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add, int_u_abs_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, N_par, N_perp
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    N_par = svec%N_cold * cos_theta !
    N_perp = svec%N_cold * sin_theta !
    m_omega_bar =  real(m,8) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    if(m_omega_bar_2  + N_par**2  - 1.e0 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    !sin_theta_2 = (sin_theta)**2!  * N**2
    !cos_theta_2 = (cos_theta)**2! * N**2
    !det = sqrt(m_omega_bar_2 - sin_theta_2)
    !a = (m_omega_bar*cos_theta - det) / sin_theta_2
    !b = (m_omega_bar*cos_theta + det) / sin_theta_2
!    if(m_omega_bar**2  + N_par**2  - 1.e0 < 0.d0) then
!      Int_u = 0.d0
!      int_u_abs = 0.d0
!      return
!    end if
    a =  (m_omega_bar * N_par - sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    b = (m_omega_bar * N_par + sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    !1 continue
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * N_par + m_omega_bar
      !gamma(k) = u_par(k) * cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)

      ! Hutchinson
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_Jn(m, zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
                     (Bessel_JN(m - 1, zeta(k)) - Bessel_Jn(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * gamma(k)**2 * Int_weights(k)
    enddo
    int_u_abs_add = -int_u_add * mu * exp( -mu / 2.d0 * ((u_par(:) / gamma(:))**2 + u_perp(:)**2 / gamma(:)**2))
    int_u_add = int_u_add  * exp( -mu / 2.d0 * ((u_par(:) / gamma(:))**2 + u_perp(:)**2 / gamma(:)**2))
    !if(stop) stop "done"
    !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
    !  stop "Sense?"
    !end if
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs  +  int_u_abs_add(k)
    end do
    !int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_refr


  function full_relativistics(cos_theta, sin_theta, freq_2X, freq)
   ! Returns true if the resonance curve is in the relativistic regime of velocity space
   ! the relativistic regime is defined here as beta >= 1 %
   ! Therefore gamma may not be larger than approximately 1.005
   use mod_ecfm_refr_types,             only: rad_diag_ch_mode_ray_freq_svec_type
   implicit none
   real(rkind),                 intent(in) :: cos_theta, sin_theta, freq_2X, freq
   logical                                 :: full_relativistics
   real(rkind)                             :: u_par_max
   real(rkind)                             :: gamma_max
   full_relativistics = .true.
   u_par_max = ((freq_2X/freq)*cos_theta + sqrt((freq_2X/freq)**2 - sin_theta**2)) / sin_theta**2
   gamma_max = (freq_2X/freq) + u_par_max * cos_theta
   total_cnt = total_cnt + 1
   if(u_par_max / gamma_max < 1.d-1) full_relativistics = .false.
   !full_relativistics = .true.
  end function full_relativistics

  subroutine calculate_int_u_gauss_opt(cos_theta, sin_theta,mu, f_rel, f_rel_sq, int_u) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none
    real(rkind),    intent(in)  :: cos_theta,sin_theta,mu,f_rel, f_rel_sq
    real(rkind),    intent(out) :: int_u
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2, ab_mean
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_omega_bar =  1.0 / f_rel
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2
    cos_theta_2 = (cos_theta)**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    ab_mean = (b + a) / 2.d0
    int_u       = 0.d0
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta + m_omega_bar
     if(gamma(k)**2 - u_par(k)**2 - 1.d0 < 0) then
      print*, k
      print*, u_par(k)
      print*, gamma(k)
      stop "imaginary value"
     end if
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* BesselJ_fast_2(zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + &
        ( u_perp(k)/gamma(k))**2  * (BesselJ_fast_1(zeta(k)) - BesselJ_fast_3(zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha *  exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights(k)
    enddo
    do k = 1, size(Int_weights)
      int_u = int_u + int_u_add(k)
    end do
  end subroutine calculate_int_u_gauss_opt

  function BesselJ_fast_1(x)
  ! Evaluates the Besselfunction of the first kind for n = 1
  ! Fast implementation of BesselJ_fast_1
  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_1
    ! BesselJ = 1/ (n !)
    BesselJ_fast_1 = 1.0 - x**2/8.d0 ! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_1 = BesselJ_fast_1 + x**4/192.d0 ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_1 = BesselJ_fast_1 - x**6/(9216.d0) ! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_1 = BesselJ_fast_1 * x / 2.0
  end function BesselJ_fast_1

  function BesselJ_fast_2(x)
  ! Evaluates the Besselfunction of the first kind for n = 1

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_2
    ! BesselJ = 1/ (n !)                :: n_fac
    BesselJ_fast_2 = 0.5 - x**2/24.d0 ! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_2 = BesselJ_fast_2 + x**4/768.d0 ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_2 = BesselJ_fast_2 - x**6/(46080.d0) ! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_2 = BesselJ_fast_2 * x**2 / 4.0
  end function BesselJ_fast_2

  function BesselJ_fast_3(x)
  ! Evaluates the Besselfunction of the first kind for n = 1

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_3
    ! BesselJ = 1/ (n !)
    BesselJ_fast_3 = 1/6.d0 - x**2/(96.d0)! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_3 = BesselJ_fast_3 + x**4/(3840.d0) ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_3 = BesselJ_fast_3 - x**6/(276480.d0)! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_3 = BesselJ_fast_3 * x**3 / (8.0)
  end function BesselJ_fast_3

  function BesselJ(n,x)
  ! Evaluates the Besselfunction of the first kind
  ! For n = 1, the relativ numerical error for x < 1.5 is below 10^-4
  ! For values in between 1.5 < x < 2.5 the numerical error is below 1 %
  ! For values 2.5 > x > 3.0 the numerical error is below 10 %
  ! For values x > 3 this routine should not be used since the truncation error becomes very large
  ! The larger n, the greater the range where this function gives an almost exact result
  ! e.g. for  n = 2 the numerical for x = 3.0 is still below 1 %

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    integer(ikind),    intent(in)  :: n
    real(rkind)                    :: besselj
    integer(ikind)                 :: i
    real(rkind)                    :: n_fac
    n_fac = 1.0
    do i= 1,n
      n_fac = real(i) * n_fac
    end do
    BesselJ = 1.0/n_fac
    n_fac = n_fac * real(n + 1)
    BesselJ = BesselJ - x**2/(n_fac * 2.d0**(2))
    n_fac = n_fac * real(n + 2)
    BesselJ = BesselJ + x**4/(n_fac * 2.d0**(5))
    n_fac = n_fac * real(n + 3)
    BesselJ = BesselJ - x**6/(n_fac * 2.d0**(7)  * 3.d0)
    BesselJ = BesselJ * x**n / (2.0)**n
  end function BesselJ



  subroutine calculate_int_u_gauss_many(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m, N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind), dimension(size(Int_weights_many))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind)     :: eps
    real(rkind)     :: a, b, du_dalpha
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta =svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2! * N**2
    cos_theta_2 = (cos_theta)**2! * N**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    do k = 1, size(Int_weights_many)
      u_par(k) = Int_absz_many(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k)  *  cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta  - u_par(k)/gamma(k))* BesselJ_fast_2(zeta(k))/ sin_theta )**2
      int_u_add(k) = int_u_add(k)  + &
        ( u_perp(k)/gamma(k))**2  * (BesselJ_fast_1(zeta(k)) - BesselJ_fast_3(zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha *  exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights_many(k)
    end do
    do k = 1, size(Int_weights_many)
      int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs - int_u_add(k) * mu
    enddo
    if(int_u < 0 .or. int_u /= int_u) then
      print*, "Negative Shape function"
      print*, Int_absz_many(1) * du_dalpha + (b + a) / 2.d0
      print*,(Int_absz_many(1) * du_dalpha + (b + a) / 2.d0)*cos_theta + m_omega_bar
      stop "Bug"
    end if
  end subroutine calculate_int_u_gauss_many

function bessel_term(u_par, u_perp, gamma,zeta, m, cos_theta, sin_theta)
implicit none

real(rkind),    intent(in)  :: u_par, u_perp, gamma,zeta, cos_theta, sin_theta
integer(ikind), intent(in)  :: m
real(rkind)                 :: bessel_term
 bessel_term = ((cos_theta - u_par/gamma)/ (sin_theta))**2 * Bessel_JN(m, zeta)**2
 bessel_term = bessel_term +  ( u_perp/gamma)**2  * &
   (Bessel_JN(m - 1, zeta) - Bessel_JN(m + 1,zeta))**2 / 4.d0

end function bessel_term
end module mod_ecfm_refr_em_Hu
