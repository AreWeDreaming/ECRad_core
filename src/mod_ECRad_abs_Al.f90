! Possible Upgrades:
! - Add support for first harmonic emission: -> See subroutine abs_Albajar

module mod_ECRad_abs_Al
! This routine calculates the absorption using the formalism from the article [1]:
! Electron-cyclotron absorption in high-temperature plasmas:
! quasi-exact analytical evaluation and comparative numerical analysis"
! by
! F. Albajar, N. Bertelli, M. Bornatici and F. Engelman
! published 2007 in Plasma Phys. Control. Fusion 49 15-29
  Use f90_kind
  Implicit none
  real(rkind), dimension(:), allocatable :: Int_weights
  real(rkind), dimension(:), allocatable :: Int_absz
  real(rkind)                            :: m_omega_glob, cos_theta_glob, mu_glob
  integer, parameter :: r8=selected_real_kind(15,300)
  public :: abs_Al_init, &
            abs_Al_clean_up, &
            func_N_cold, &
            func_rel_N, &
            abs_Al_Fa_abs, &
            abs_Albajar, &
            abs_Albajar_fast
  private :: Int_weights, &
             Int_absz, &
             rotate_vec_around_axis, abs_Al_N_with_pol_vec, &
             abs_Al_integral_nume, &
             abs_Al_integral_nume_fast, &
             abs_Al_pol_fact, &
             get_E_factors, &
             BesselJ, &
             BesselJ_custom, &
             get_filter_transmittance




contains

  subroutine rotate_vec_around_axis(vec, axis, theta, vec_rot)
    use f90_kind
    Implicit none
    complex(r8), dimension(:), intent(in) :: vec
    real(rkind), dimension(:), intent(in) :: axis
    real(rkind)                           :: theta
    complex(r8), dimension(:), intent(out):: vec_rot
    real(rkind)                           :: a, b, c, d, aa, bb, cc, dd, &
                                             ab, ac, ad, bc, bd, cd
    real(rkind), dimension(3, 3)          :: rotation_matrix
    real(rkind), dimension(3)             :: axis_norm
   ! Return the rotation matrix associated with counterclockwise rotation about
   ! the given axis by theta radians. Source: https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
   ! This routine has been double checked and returns identical results to initial calculation
    axis_norm = axis / sqrt(dot_product(axis, axis))
    if(abs(theta) < 1.d-5) then
        vec_rot = vec
    else
        a = cos(theta / 2.0)
        b = -axis_norm(1) * sin(theta / 2.0)
        c = -axis_norm(2) * sin(theta / 2.0)
        d = -axis_norm(3) * sin(theta / 2.0)
        aa = a * a
        bb = b * b
        cc = c * c
        dd = d * d
        bc = b * c
        ad = a * d
        ac = a * c
        ab = a * b
        bd = b * d
        cd = c * d
        rotation_matrix(1,1) = aa + bb - cc - dd
        rotation_matrix(1,2) = 2 * (bc + ad)
        rotation_matrix(1,3) = 2 * (bd - ac)
        rotation_matrix(2,1) = 2 * (bc - ad)
        rotation_matrix(2,2) =  aa + cc - bb - dd
        rotation_matrix(2,3) = 2 * (cd + ab)
        rotation_matrix(3,1) = 2 * (bd + ac)
        rotation_matrix(3,2) = 2 * (cd - ab)
        rotation_matrix(3,3) = aa + dd - bb - cc
        vec_rot(1) = dot_product(rotation_matrix(1,:), vec)
        vec_rot(2) = dot_product(rotation_matrix(2,:), vec)
        vec_rot(3) = dot_product(rotation_matrix(3,:), vec)
    end if
  end subroutine rotate_vec_around_axis

  subroutine abs_Al_init(N_absz)
#ifdef NAG
    use nag_quad_util,              only: nag_quad_gs_wt_absc
    USE nag_error_handling
#endif
    use quadrature,                 only: cdgqf
    use mod_ECRad_abs_Fa,       only: set_extv
    use mod_ECRad_types,        only: dstf_comp
    implicit none
    integer(ikind)                         :: N_absz
#ifdef NAG
    type(nag_error)                        :: error
    real(rkind), dimension(:), allocatable :: Int_weights_check, Int_absz_check
#endif
    allocate( Int_weights(N_absz),Int_absz(N_absz))
    call cdgqf( int(N_absz,kind=4), int(1,kind=4), 0.d0, 0.d0, Int_absz, Int_weights)
#ifdef NAG
    allocate( Int_weights_check(N_absz), Int_absz_check(N_absz))
    call nag_quad_gs_wt_absc( 0, -1.d0, 1.d0, Int_weights_check, Int_absz_check)
    if(sum((Int_weights - Int_weights_check)**2) > 1.d-5) then
      print*, "Weights deviate by more than 1.d-10"
      do i = 1, N_absz
        print*, Int_weights(i), Int_weights_check(i)
      end do
      call abort()
    end if
    if(sum((Int_absz - Int_absz_check)**2) > 1.d-5) then
      print*, "Abszissae deviate by more than 1.d-10"
      do i = 1, N_absz
        print*, Int_absz(i), Int_absz_check(i)
      end do
      call abort()
    end if
    !print*, Int_absz(1)
    if (error%level >= 1) print*, error%msg
    if (error%level >= 1) stop "Nag Screwed up when setting weights for Albajar abs. coeff."
    deallocate(Int_weights_check, Int_absz_check)
#endif
    if(dstf_comp == "DF") call set_extv
  end subroutine abs_Al_init

  subroutine abs_Al_clean_up()
  implicit none
    if(allocated(Int_weights)) deallocate(Int_weights)
    if(allocated(Int_absz)) deallocate(Int_absz)

  end subroutine abs_Al_clean_up

  function func_N_cold( omega, svec,mode) ! Following [2]
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none
    real(rkind), intent(in)       :: omega
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
    integer(ikind), intent(in)    :: mode
    !mode =  +1 -> X | mode = -1 -> O
    real(rkind)                   :: func_N_cold
    real(rkind)                   :: rho, f, omega_c, omega_p_2, X, Y, sin_theta, cos_theta
    omega_c = svec%freq_2X * Pi
    omega_p_2 = svec%ne * e0**2 / (eps0 * mass_e)
    X = omega_p_2 / omega**2
    if(X >= 1.d0) then
      func_N_cold = 0.d0
      return
    end if
    Y = omega_c / omega
    sin_theta = svec%sin_theta
    cos_theta = svec%cos_theta
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    if(rho < 0.d0) then
      func_N_cold = 0.d0
      return
    end if
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 - real(mode,8) *  Y* rho)
    func_N_cold = 1.d0 - X * f
    if(func_N_cold < 0.d0) then
      func_N_cold = 0.d0
      return
    end if
    func_N_cold = sqrt(func_N_cold)
  end function func_N_cold

  function func_rel_N( omega, svec,mode) ! Following [2]
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none
    real(rkind), intent(in)       :: omega
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
    integer(ikind), intent(in)    :: mode
    !mode =  +1 -> X | mode = -1 -> O
    real(rkind)                   :: func_rel_N
    real(rkind)                   :: rho, f, w_mass_e, omega_c, omega_p_2, X, Y, sin_theta, cos_theta, mu
    mu = mass_e * c0**2 / (e0 * svec%Te)
    w_mass_e = mass_e * sqrt(1.d0 + 5.d0 / mu)
    omega_c = svec%freq_2X * Pi * mass_e / w_mass_e
    omega_p_2 = svec%ne * e0**2 / (eps0 * w_mass_e)
    X = omega_p_2 / omega**2
    if(X >= 1.d0) then
      func_rel_N = 0.d0
      return
    end if
    Y = omega_c / omega
    sin_theta = svec%sin_theta
    cos_theta = svec%cos_theta
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    if(rho < 0.d0) then
      func_rel_N = 0.d0
      return
    end if
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 - real(mode,8) *  Y* rho)
    func_rel_N = 1.d0 - X * f
    if(func_rel_N < 0.d0) then
      func_rel_N = 0.d0
      return
    end if
    func_rel_N = sqrt(func_rel_N)
  end function func_rel_N

  subroutine abs_Al_N_with_pol_vec( X, Y, cos_theta, sin_theta, mode, N, e) ! Following [2]
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    real(rkind), intent(in)                :: X, Y, cos_theta, sin_theta
    integer(ikind), intent(in)             :: mode
    real(rkind), intent(out)               :: N
    complex(r8), dimension(3), intent(out) :: e
    real(rkind)                            :: rho, f, a_sq, b_sq
    e(:) = 0.d0
    if(X >= 1.d0) then
      N = 0.d0
      return
    end if
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    if(rho < 0.d0) then
      N = 0.d0
      return
    end if
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 - real(mode,8) *  Y* rho)
    N = 1.d0 - X * f
    if(N < 0.d0) then
      N = 0.d0
      f = 0.d0
      return
    end if
    N = sqrt(N)
    if(cos_theta ** 2 < 1.d-5 .or. 1.0 - sin_theta ** 2 < 1.d-5) then
    ! Quasi-perpendicular
      if(mode > 0) then  ! X-mode
                e(2) = cmplx(0.0, sqrt(1.0 / N))
                e(1) = cmplx(0.0, 1.0 / Y * (1.0 - (1.0 - Y ** 2) * f)) * e(2)
                ! e_z zero for quasi-perpendicular X-mode
      else
          e(3) = sqrt(1.d0 / N)  !
      end if
    else
      a_sq = sin_theta**2 * (1.d0 + (((1.d0 - X) * N**2 * cos_theta**2) / &
        (1.d0 - X - N**2 * sin_theta**2)**2) * &
        1.d0 / Y**2 * (1.d0 - (1.d0 - Y**2)* f)**2)**2
      b_sq = cos_theta**2 * (1.d0 + ((1.d0 - X) / &
        (1.d0 - X - N**2 * sin_theta**2)) * &
        1.d0 / Y**2 * (1.d0 - (1.d0 - Y**2)* f)**2)**2
      if(mode > 0) then
      ! X-mode - e_y positive, purely imaginary
        e(2) = cmplx(0.d0, sqrt(1.d0 / (N * sqrt(a_sq + b_sq))))
      else
      ! O-mode -e_y negative, purely imaginary
        e(2) = cmplx(0.d0, -sqrt(1.d0 / (N * sqrt(a_sq + b_sq))))
      end if
      e(1) = cmplx(0.d0, 1.d0 / Y * (1.d0 - (1.d0 - Y ** 2) * f)) * e(2)
      e(3) = -cmplx((N ** 2 * sin_theta * cos_theta) / (1.d0 - X - N ** 2 * sin_theta ** 2), 0.d0) * e(1)
    end if
  end subroutine abs_Al_N_with_pol_vec


  function abs_Al_Fa_abs(svec, omega, mode, Nr, pol_coeff_secondary, pol_vec, x_launch) ! note that this routine uses angular frequency
  ! Calculates the absorption coefficient using Grays warm_disp routine. Note that this always includes both 2nd and 3rd harmonic.
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type, &
                                          ratio_for_third_harmonic, dstf, straight, Hamil, &
                                          ignore_Te, ignore_ne
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ECRad_abs_Fa,       only: warmdamp
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
    real(rkind), intent(in)       :: omega
    integer(ikind), intent(in)       :: mode
    real(rkind), intent(inout)       :: Nr
    real(rkind), intent(out), optional :: pol_coeff_secondary
    complex(r8), dimension(3), intent(out), optional :: pol_vec
    real(rkind), dimension(:), intent(in), optional :: x_launch
    real(rkind)                   :: abs_Al_Fa_abs
    real(rkind)                   :: omega_p,omega_c, alpha, beta, omega_p_sq
    real(rkind)                   :: ni_perp_sq
    complex(r8)                   :: N_perp_cmplx, scal_prod
    complex(r8), dimension(3)     :: pol_vec_dummy
    real(rkind), dimension(3)     :: abs_pol_vec
    integer(ikind)                :: max_harmonic
      abs_Al_Fa_abs = 0.d0
      if(svec%Te < ignore_Te .and. .not. present(pol_coeff_secondary)) return ! very low Te => absorption can be ignored
      if(svec%ne < ignore_ne .and. .not. present(pol_coeff_secondary)) return
      ni_perp_sq = 0.d0
      omega_c        = svec%freq_2X * Pi
      omega_p_sq     = (svec%ne * e0**2.d0)/(eps0 * mass_e)
      omega_p = sqrt(omega_p_sq)
      alpha = omega_p_sq / omega**2
      beta = omega_c**2 / omega**2
      max_harmonic = 2 ! consider only second harmonic
      if(dstf /= "Th") then
        max_harmonic = 5 ! always also third harmonic for non-thermal discharges
      else
        if(sqrt(beta) < ratio_for_third_harmonic) max_harmonic = 3
      end if
      if(present(pol_coeff_secondary) .or. present(pol_vec)) then
        if(.not. present(pol_vec)) then
          call warmdamp(alpha, beta, Nr, svec%theta, svec%Te, max_harmonic, mode, N_perp_cmplx, pol_vec = pol_vec_dummy)!-1
          abs_pol_vec(1) = abs(pol_vec_dummy(1))
          abs_pol_vec(2) = abs(pol_vec_dummy(2))
          abs_pol_vec(3) = abs(pol_vec_dummy(3))
        else
          call warmdamp(alpha, beta, Nr, svec%theta, svec%Te, max_harmonic, mode, N_perp_cmplx, pol_vec = pol_vec)!-1
          pol_vec_dummy = pol_vec
          abs_pol_vec(1) = abs(pol_vec(1)) ! Absolute value - the sign just gives us the phase shift which is averaged out
          abs_pol_vec(2) = abs(pol_vec(2))
          abs_pol_vec(3) = abs(pol_vec(3))
        end if
        scal_prod = 0.d0
        if(present(pol_coeff_secondary)) pol_coeff_secondary = get_filter_transmittance(omega, &
                                         alpha, sqrt(beta), svec%cos_theta, svec%sin_theta, &
                                         mode, svec%x_vec, svec%N_vec, svec%B_vec, x_launch, &
                                         pol_vec_ext = pol_vec_dummy)
      else
       call warmdamp(alpha, beta, Nr, svec%theta, svec%Te, max_harmonic, mode, N_perp_cmplx)!-1
      end if
      if(straight .or. Hamil /= "Dani") then
        abs_Al_Fa_abs = 2.d0 * aimag(N_perp_cmplx) * omega / c0 * svec%sin_theta
      else
        abs_Al_Fa_abs = 2.d0 * aimag(N_perp_cmplx**2) * omega / c0 * svec%v_g_perp
      end if
      Nr = real(N_perp_cmplx,8) / svec%sin_theta
      if(abs_Al_Fa_abs .ne. abs_Al_Fa_abs .or. Nr <= 0.d0 .or. abs_Al_Fa_abs < 0.d0 .or. Nr .ne. Nr .or. &
        abs_Al_Fa_abs > 1.e6) then  ! .or. Nr > 1.d0
        Nr = 0.d0
        abs_Al_Fa_abs = 0.d0
      end if
  end function abs_Al_Fa_abs

  function get_upper_limit_tau(svec, omega, ds2)
  ! Gets upper limit of absorption coefficient by evaluating the integral
  ! Int u_perp ^(2 * m) Exp[mu(1-gamma)]
  ! This uses the Expansion of the Besselfunctions (see Hutchinson/ S. Denk Msc. Thesis) and neglects the term gamma(-2n + 2) making it an
  ! estimation for the upper limit
  ! The integration was carried out with mathematica
  ! Considers n=2 and n=3
  use mod_ECRad_types,        only: ratio_for_third_harmonic, rad_diag_ch_mode_ray_freq_svec_type
  use constants,                  only: pi, e0, mass_e, eps0, c0
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
  real(rkind), intent(in)       :: omega, ds2
  real(rkind)                   :: get_upper_limit_tau
  real(rkind)                   ::  omega_bar, m_omega_bar,  omega_c, mu, &
                                   norm, f, f_sum, N_par, t1, t2
  integer(ikind)                :: m, m_max
  get_upper_limit_tau = 0.d0
  omega_c = (svec%freq_2X * pi)
  omega_bar =  omega / omega_c
  m_max = 2
  if(1.d0 / omega_bar < ratio_for_third_harmonic) m_max = 3
  mu = mass_e * c0**2 / (e0 * svec%Te)
  f_sum = 0.d0
  N_par = svec%cos_theta ! We do not consider N_abs here -> tenous plasma limit
  ! This causes for much larger values of alpha in the oblique regime
  ! If we use N_par = N_abs * cos_theta this formula underestimates alpha for |cos_theta - 90| > 5 deg which we do not want
  ! Rather evaluate more than necessary instead of omitting important points
  do m =2, m_max
    m_omega_bar = real(m, 8) / omega_bar
    if(m_omega_bar**2 + N_par**2 < 1.d0) cycle
    if(m == 2) then
      t1 =   (-8*Exp((mu*(-1 + N_par**2 + m_omega_bar - N_par*Sqrt(-1 + N_par**2 + m_omega_bar**2)))/ &
            (-1 + N_par**2))*(-3 - N_par**4*(3 + mu**2) +  &
            3*N_par*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2) -  &
            3*svec%cos_theta**3*mu*Sqrt(-1 + svec%cos_theta**2 + m_omega_bar**2) + &
            svec%cos_theta**2*(6 - mu**2*(-1 + m_omega_bar**2))))/(svec%cos_theta**5*mu**5)
      t2 =  (-8*Exp((mu*(-1 + N_par**2 + m_omega_bar + N_par*Sqrt(-1 + N_par**2 + m_omega_bar**2)))/ &
             (-1 + N_par**2))*(3 + N_par**4*(3 + mu**2) + &
            3*N_par*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2) - &
            3*N_par**3*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2) + &
            N_par**2*(-6 + mu**2*(-1 + m_omega_bar**2))))/(N_par**5*mu**5)
    else
      t1 =   (48*Exp((mu*(-1 + N_par**2 + m_omega_bar - N_par*Sqrt(-1 + N_par**2 + m_omega_bar**2)))/ &
            (-1 + N_par**2))*(-15 + 3*N_par**6*(5 + 2*mu**2) + &
            15*N_par*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2) + &
            N_par**5*mu*(15 + mu**2)*Sqrt(-1 + N_par**2 + m_omega_bar**2) + &
            3*N_par**4*(-15 + 2*mu**2*(-2 + m_omega_bar**2)) + &
            N_par**2*(45 - 6*mu**2*(-1 + m_omega_bar**2)) + &
            N_par**3*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2)*(-30 + mu**2*(-1 + m_omega_bar**2))))/ &
            (N_par**7*mu**7)
      t2 =  (48*Exp((mu*(-1 + N_par**2 + m_omega_bar + N_par*Sqrt(-1 + N_par**2 + m_omega_bar**2)))/ &
             (-1 + N_par**2))*(15 - 3*N_par**6*(5 + 2*mu**2) + &
            15*N_par*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2) + &
            N_par**5*mu*(15 + mu**2)*Sqrt(-1 + N_par**2 + m_omega_bar**2) + &
            N_par**4*(45 - 6*mu**2*(-2 + m_omega_bar**2)) + &
            N_par**3*mu*Sqrt(-1 + N_par**2 + m_omega_bar**2)*(-30 + mu**2*(-1 + m_omega_bar**2)) + &
            N_par**2*(-45 + 6*mu**2*(-1 + m_omega_bar**2))))/(N_par**7*mu**7)
    end if
    f = (t1 + t2) * mu * (N_par**2 + 1.d0) * m**2 ! from omega_m**2
    if( m == 3) then
      f = f * 20.25d0 * (1.d0 - N_par**2)**2 * 0.5**6
    else
      f = f * 4.0d0 * (1.d0 - N_par**2) * 0.5**4
    end if
    f_sum = f_sum + f ! * arc_length
  end do
  if(f_sum == 0.d0) return ! no resonance -> no contribution
  norm = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
  get_upper_limit_tau = f_sum * norm * (sqrt(mu / (2 * pi))**3) * ds2
  get_upper_limit_tau = get_upper_limit_tau * (e0 * omega_c)**2.d0 * &
                        svec%ne / (eps0 * c0 ) * pi * 8.d0 / (omega**3 * mass_e)


  end function get_upper_limit_tau

  subroutine abs_Albajar(svec, omega, mode, ds2,  c_abs, j, pol_coeff, c_abs_secondary, j_secondary, x_launch)
    ! Calculates the absorption coefficient and emissivity
    ! This algorithm is developed according to F. Albajar et al (2006) [1]
    ! Note that an arbitrary propagation direction and a polarization filter in phi (tokamak) direction is considered.
    ! Does not support the fundamental
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type, dstf, output_level,&
                                          ratio_for_third_harmonic, not_eval, eval, warm_plasma, &
                                          tau_ignore, spl_type_2d, non_therm_params_type, &
                                          ignore_Te, ignore_ne
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ECRad_radiation_dist,    only: prepare_dist
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
    real(rkind), intent(in)       :: omega, ds2
    integer(ikind), intent(in)    :: mode
    real(rkind), intent(out)      :: c_abs, j
    real(rkind), intent(out), optional      :: pol_coeff
    real(rkind), intent(out), optional      :: c_abs_secondary, j_secondary
    real(rkind), dimension(:), intent(in), optional :: x_launch
    integer(ikind)                :: m_sum, max_harmonic
    real(rkind)                   :: c_abs_m, c_abs_m_secondary, j_m, j_m_secondary
    real(rkind)                   :: w_mass_e, omega_p,omega_p_cold, mu, omega_c, &
                                     X, Y, omega_bar, N_abs, N_par, N_perp, m_0
    complex(r8), dimension(3)     :: e
    type(spl_type_2d) :: f_spl
    type(non_therm_params_type)        :: dist_params
    c_abs = 0.d0
    if(present(c_abs_secondary)) c_abs_secondary = 0.d0
    if(present(j_secondary)) j_secondary = 0.d0
    j = 0.d0
    if(svec%Te < ignore_Te .and. .not. present(pol_coeff)) return ! very low Te => absorption can be ignored
    if(svec%ne < ignore_ne .and. .not. present(pol_coeff)) return
    if(svec%ne < 0.d0) then
      print*, "Negative density!", svec%ne
      call abort()
    end if
    mu = mass_e * c0**2 / (e0 * svec%Te)
    if(warm_plasma) then
      w_mass_e = mass_e * sqrt(1.d0 + 5.d0 / mu)
    else
      w_mass_e = mass_e
    end if
    omega_c = svec%freq_2X * Pi
    Y = svec%freq_2X * Pi * mass_e / w_mass_e / omega
    max_harmonic = 2 ! consider only second harmonic
    omega_p = e0 * sqrt( svec%ne / (eps0 * w_mass_e))
    omega_p_cold = e0 * sqrt( svec%ne / (eps0 * mass_e))
    X = omega_p**2 / omega**2
    omega_bar = omega / omega_c
    call abs_Al_N_with_pol_vec(X, Y, svec%cos_theta, svec%sin_theta, mode, N_abs, e) ! mode = X -> + 1
    if(N_abs /= N_abs .or. N_abs <= 0.0 .or. N_abs > 1.0) return
    if(dstf /= "Th") then
      call prepare_dist(svec, Int_absz, Int_weights, f_spl, dist_params)
      max_harmonic = 5
    else
      if((Y * w_mass_e / mass_e) < ratio_for_third_harmonic) max_harmonic = 3
      if(get_upper_limit_tau(svec, omega, ds2) < tau_ignore .and. .not. present(pol_coeff)) then
        not_eval = not_eval + 1
        return
      else
        eval = eval + 1
      end if
    end if
    if(present(pol_coeff)) then
      pol_coeff = get_filter_transmittance(omega, X, Y, svec%cos_theta, svec%sin_theta, mode, &
                                           svec%x_vec, svec%N_vec, svec%B_vec, x_launch)
    end if
    if(present(pol_coeff) .and. svec%Te < ignore_Te) then
      c_abs = 0.d0
      if(present(c_abs_secondary)) c_abs_secondary = 0.d0
      if(present(j_secondary)) j_secondary = 0.d0
      return
    end if
    N_par = svec%cos_theta * N_abs
    N_perp = abs(svec%sin_theta * N_abs)
    m_0 = sqrt(1.d0 - N_par**2) * omega_bar
    ! Here would be a possible place to add the first harmonic
    ! A few expressions can be found in Bornatici 1983
    ! They could be easily benchmarked against the absoprtion routine by D. Farina (see above)
    do m_sum = 2, max_harmonic ! First harmonic needs to be treated seperately (for now ignored)
      if(real(m_sum,8) < m_0 ) cycle
      if(present(c_abs_secondary) .and. present(j_secondary) ) then
        call abs_Al_integral_nume(svec, f_spl, dist_params, X, Y, omega_bar, m_0, N_abs, svec%cos_theta, svec%sin_theta, e, mode, m_sum, c_abs_m, j_m, &
                                  c_abs_secondary=c_abs_m_secondary, j_secondary=j_m_secondary)
      else if(present(c_abs_secondary)) then
        call abs_Al_integral_nume(svec, f_spl, dist_params, X, Y, omega_bar, m_0, N_abs, svec%cos_theta, svec%sin_theta, e, mode, m_sum, c_abs_m, j_m, &
                                  c_abs_secondary= c_abs_m_secondary)
      else if(present(j_secondary)) then
        call abs_Al_integral_nume(svec, f_spl, dist_params, X, Y, omega_bar, m_0, N_abs, svec%cos_theta, svec%sin_theta, e, mode, m_sum, c_abs_m, j_m, &
                                  j_secondary=j_m_secondary)
      else
        call abs_Al_integral_nume(svec, f_spl, dist_params, X, Y, omega_bar, m_0, N_abs, svec%cos_theta, svec%sin_theta, e, mode, m_sum, c_abs_m, j_m)
      end if
      c_abs_m = -(c_abs_m * 2.d0 * pi**2 / m_0) ! Splitting this is just for overview
      c_abs_m = c_abs_m * omega_p_cold**2 / (omega_c  * c0 ) ! revert the norminalization (w /wp^2)
      if(present(c_abs_secondary)) then
        c_abs_m_secondary = -(c_abs_m_secondary * 2.d0 * pi**2 / m_0) ! Splitting this is just for overview
        c_abs_m_secondary =  c_abs_m_secondary * omega_p_cold**2 / (omega_c  * c0 )
      end if
      j_m = j_m * 2.d0 * pi**2 / m_0
      j_m = j_m *  omega_p_cold**2 / (omega_c  * c0 ) * omega**2 * mass_e / (4.d0 * pi**2)
      if(present(j_secondary)) then
        j_m_secondary = j_m_secondary * 2.d0 * pi**2 / m_0
        j_m_secondary = j_m_secondary *  omega_p_cold**2 / (omega_c  * c0 ) * omega**2 * mass_e / (4.d0 * pi**2)
      end if
      c_abs = c_abs + sqrt((real(m_sum,8) / m_0)**2 - 1.d0) * c_abs_m
      if(present(c_abs_secondary)) c_abs_secondary = c_abs_secondary + sqrt((real(m_sum,8) / m_0)**2 - 1.d0) *  c_abs_m_secondary
      j = j + sqrt((real(m_sum,8) / m_0)**2 - 1.d0) * j_m
      if(present(j_secondary)) j_secondary = j_secondary + sqrt((real(m_sum,8) / m_0)**2 - 1.d0) * j_m_secondary
      if( c_abs_m /= c_abs_m  .or. (c_abs_m < 0.0 .and. .not. (trim(dstf) == 'gene' .or. dstf == "numeric")) &
          .or. j_m < 0.d0) then !c_abs < 0.d0 .or.
        print*, "rhop", svec%rhop
        print*, "Te",svec%Te
        print*, "ne",svec%ne
        print*, "freq", omega / (2.d0 * pi)
        print*, "wp/wc",omega_p / (svec%freq_2X * Pi )
        print*, "int", c_abs
        print*, "N", N_abs
        print*, "pref",  sqrt((real(m_sum,8) / m_0)**2 - 1.d0)
        print*, "m_0",m_0
        print*, "Y", 1.d0 / omega_bar
        print*, "N_par", N_par
        print*, "cos theta/ theta", svec%cos_theta, svec%theta * 180.d0 / pi
        print*, "c_abs_m", c_abs_m
        print*, "harmonic", m_sum
        print*, "j_m*ds/Ibb", j_m * ds2 / svec%Ibb
        print*, "j_m", j_m
        print*, "fully rel. dispersion c_abs", abs_Al_Fa_abs(svec, omega, mode, N_abs)
        stop "Nan in c_abs Albajar"
      end if
    end do
    if(output_level) then
      m_0 = get_upper_limit_tau(svec,  omega, ds2)
      if(m_0 > 1.d-30 .and. c_abs > 1.d-30 .and. dstf == "relamax") then
        if(abs(m_0 / ds2 / c_abs) < 1.d-1) then
          print*, "Estimate", m_0 / ds2
          print*, "Abs", c_abs
          print*, abs(m_0 / ds2 / c_abs)
          print*, "Large underestimation of c_abs"
          c_abs = - 1.d0 ! for debug output
        end if
      end if
    end if
  end subroutine abs_Albajar

  subroutine abs_Albajar_fast(svec, omega, mode, ds2, c_abs)
    ! Calculates the absorption coefficient -> reduced version of abs_Albajar for faster evaluation
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type, ratio_for_third_harmonic, eval, not_eval, tau_ignore
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)    :: svec
    real(rkind), intent(in)       :: omega, ds2
    integer(ikind), intent(in)    :: mode
    real(rkind), intent(out)      :: c_abs
    integer(ikind)                :: m_sum, max_harmonic
    real(rkind)                   :: c_abs_m
    real(rkind)                   :: omega_p, mu, omega_c, &
                                     X, Y, omega_bar, N_abs, N_par, N_perp, m_0, &
                                     c_abs_int, tau_upper_limit
    complex(r8), dimension(3)     :: e
    c_abs = 0.d0
    if(svec%Te < 20.d0) return ! very low Te => absorption can be ignored
    if(svec%ne < 1.e17) return
    mu = mass_e * c0**2 / (e0 * svec%Te)
    omega_c = svec%freq_2X * Pi
    Y = svec%freq_2X * Pi / omega
    max_harmonic = 2 ! consider only second harmonic
    if(Y < ratio_for_third_harmonic) max_harmonic = 3
    omega_p = e0 * sqrt( svec%ne / (eps0 * mass_e))
    X = omega_p**2 / omega**2
    omega_bar = omega / omega_c
    c_abs = 0.d0
    call abs_Al_N_with_pol_vec(X, Y, svec%cos_theta, svec%sin_theta, mode, N_abs, e) ! mode = X -> + 1
    if(N_abs /= N_abs .or. N_abs <= 0.0 .or. N_abs > 1.0) return
    tau_upper_limit = get_upper_limit_tau(svec, omega, ds2)
    if(tau_upper_limit < tau_ignore) then
      not_eval = not_eval + 1
      return
    end if
    N_par = svec%cos_theta * N_abs
    N_perp = abs(svec%sin_theta * N_abs)
    eval = eval + 1
    m_0 = sqrt(1.d0 - N_par**2) * omega_bar
    do m_sum = 2, max_harmonic ! First harmonic needs to be treated seperately (for now ignored)
      if(real(m_sum,8) < m_0 ) cycle
      call abs_Al_integral_nume_fast(svec, X, Y, omega_bar, m_0, N_abs, svec%cos_theta, svec%sin_theta, e, m_sum, c_abs_m)
      c_abs = c_abs + sqrt((real(m_sum,8) / m_0)**2 - 1.d0) * c_abs_m
    end do
    c_abs_int = c_abs
    c_abs = -(c_abs * 2.d0 * pi**2 / m_0) ! Splitting this is just for overview
    c_abs = c_abs * omega_p**2 / (omega_c  * c0 ) ! revert the norminalization (w /wp^2)
    if(c_abs < 0.0 .and. c_abs > -1.d-1) c_abs = 0.d0
    if( c_abs /= c_abs .or. c_abs < 0.d0 ) then !
      print*, "rhop", svec%rhop
      print*, "Te",svec%Te
      print*, "ne",svec%ne
      print*, "freq", omega / (2.d0 * pi)
      print*, "wp/wc",omega_p / (svec%freq_2X * Pi )
      print*, "int", c_abs_int
      print*, "N", N_abs
      print*, "pref",sqrt((real(2,8) / m_0)**2 - 1.d0)
      print*, "m_0",m_0
      print*, "Y", 1.d0 / omega_bar
      print*, "N_par", N_par
      print*, "cos theta/ theta", svec%cos_theta, svec%theta * 180.d0 / pi
      print*,"c_abs", c_abs
      stop "Nan or negative in c_abs Albajar"
    end if
  end subroutine abs_Albajar_fast


  subroutine abs_Al_integral_nume(svec, f_spl, dist_params, X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta, e, mode, m,c_abs, j, c_abs_secondary, j_secondary, debug)
    Use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type, dstf, ffp, spl_type_2d, non_therm_params_type
    use constants,                  only: pi,e0, mass_e ,c0
    use mod_ECRad_radiation_dist,    only: radiation_dist_f_norm, make_f_and_Rf_along_line
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)   :: svec
    type(spl_type_2d), intent(in)                           :: f_spl
    type(non_therm_params_type), intent(in)                   :: dist_params
    real(rkind), intent(in)                            :: X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta
    complex(r8), dimension(:), intent(in)              :: e
    integer(ikind), intent(in)                         :: mode, m
    real(rkind), intent(out)                           :: c_abs, j
    real(rkind), intent(out), optional                 :: c_abs_secondary, j_secondary
    logical, intent(in), optional                      :: debug
    real(rkind), dimension(size(Int_weights)) :: u_par, u_perp_sq, gamma, pol_fact, c_abs_int, j_int, &
                                                 c_abs_secondary_int,  j_secondary_int, f_dist, Rf_dist, &
                                                 f_dist_comp, Rf_dist_comp ! these two belong to the distribution
    real(rkind)     :: m_omega_bar, mu, norm, N_par, N_perp
    integer(ikind)                :: k
    c_abs= 0.d0
    N_par = N_abs * cos_theta
    N_perp = abs(N_abs * sin_theta)
    if(present(c_abs_secondary)) c_abs_secondary= 0.d0
    j = 0.d0
    if(present(j_secondary)) j_secondary= 0.d0
    m_omega_bar = real(m) / omega_bar
    mu = c0**2 * mass_e/(svec%Te * e0)
    u_par(:) = 1.d0 / sqrt(1.d0 - N_par**2) * ( real(m,8)/m_0 * N_par + &
                       sqrt((real(m) /  m_0)**2 - 1.d0) * Int_absz(:))
    if(any(u_par /= u_par)) then
      print*, "Nan in u_par encountered"
      print*, "N_par", N_par
      print*, "m_0", m_0
      stop "Nan in resonance (u_par) in mod_ECRad_abs_Alb.f90"
    end if
    u_perp_sq(:) =    ((real(m,8)/m_0)**2 - 1.d0) *(1.d0 - Int_absz(:)**2)
    if(any(u_perp_sq < 0)) then
      print*, "u_perp_sq smaller zero", u_perp_sq
      stop "Nan in resonance (u_perp)  in mod_ECRad_abs_Alb.f90"
    end if
    gamma(:) = sqrt(1.d0 + u_par(:)**2 + u_perp_sq(:))
    call abs_Al_pol_fact(svec, Int_absz, X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta, e, mode,  m, pol_fact)
    call make_f_and_Rf_along_line(u_par, sqrt(u_perp_sq), gamma, m_omega_bar, N_par, mu, svec, f_spl, dist_params, dstf, f_dist, Rf_dist)
    ! primary distribution
    c_abs_int =  Int_weights * pol_fact * Rf_dist
    j_int =  Int_weights * pol_fact * f_dist
    ! Secondary distribution
    if(present(c_abs_secondary) .and. present(j_secondary)) then
      if(dstf == "gene") then
        call make_f_and_Rf_along_line(u_par, sqrt(u_perp_sq), gamma, m_omega_bar, N_par, mu, svec, f_spl, dist_params, "genef0", f_dist_comp, Rf_dist_comp)
      else if(dstf == "gcomp") then
        call make_f_and_Rf_along_line(u_par, sqrt(u_perp_sq), gamma, m_omega_bar, N_par, mu, svec, f_spl, dist_params, "gcomp0", f_dist_comp, Rf_dist_comp)
      else
      ! Thermal distribution for comparison with non-thermal distributions -> Could also use Kirchhoff's öaw for either secondary j or abs
        call make_f_and_Rf_along_line(u_par, sqrt(u_perp_sq), gamma, m_omega_bar, N_par, mu, svec, f_spl, dist_params, "Th", f_dist_comp, Rf_dist_comp)
      end if
    else if((present(c_abs_secondary) .and. .not. present(j_secondary)) .or. &
            (.not. present(c_abs_secondary) .and. present(j_secondary))) then
      print*, "Currently it is only sensible to call abs_Albajar with both, secondary c_abs and secondary j or neither of the two"
      print*, dstf
      call abort()
    end if
    c_abs_secondary_int =  Int_weights * pol_fact * Rf_dist_comp
    j_secondary_int =  Int_weights * pol_fact * f_dist_comp
    do k = 1, size(Int_weights)
      c_abs = c_abs + c_abs_int(k)
      if(present(c_abs_secondary)) c_abs_secondary = c_abs_secondary + c_abs_secondary_int(k)
      j = j + j_int(k)
      if(present(j_secondary)) j_secondary = j_secondary + j_secondary_int(k)
    end do
    norm = radiation_dist_f_norm(svec%Te, dstf)
    c_abs = c_abs * norm
    j = j * norm
    if(present(c_abs_secondary) .and. (dstf /= "gene" .and. dstf /= "gcomp")) c_abs_secondary = c_abs_secondary * radiation_dist_f_norm(svec%Te, "Th")
    if(present(j_secondary) .and. (dstf /= "gene" .and. dstf /= "gcomp")) j_secondary = j_secondary * radiation_dist_f_norm(svec%Te, "Th")
    if(sum(c_abs_int) > 1.d-8 .and. .not. trim(dstf) == 'gene' .and. .not. all(c_abs_int < 4.d-9)) then
      print*, "Error region of significant negative absorption found"
      print*, "u_par u_perp, c_abs"
      do k = 1, size(Int_weights)
        print*, u_par(k), sqrt(u_perp_sq(k)), c_abs_int(k)
      end do
      print*, "u_max", ffp%u_max
    end if
    if(present(debug) .or. j < 0.d0) then
      if(debug) then
        print*, "u_par u_perp, c_abs"
        do k = 1, size(Int_weights)
          print*, u_par(k), sqrt(u_perp_sq(k)), c_abs_int(k)
        end do
        print*, "u_par u_perp, f"
        do k = 1, size(Int_weights)
          print*, u_par(k), sqrt(u_perp_sq(k)), f_dist(k)
        end do
        print*, "u_par u_perp, Rf"
        do k = 1, size(Int_weights)
          print*, u_par(k), sqrt(u_perp_sq(k)), Rf_dist(k)
        end do
      end if
    end if
  end subroutine abs_Al_integral_nume

  subroutine abs_Al_integral_nume_fast(svec, X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta, e,  m,c_abs)
    Use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    use constants,                  only: pi,e0, mass_e ,c0
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)   :: svec
    real(rkind), intent(in)                            :: X, Y, omega_bar ,m_0, N_abs, cos_theta, sin_theta
    complex(r8), dimension(:), intent(in)              :: e
    integer(ikind), intent(in)                         :: m
    real(rkind), intent(out)                           :: c_abs
    real(rkind), dimension(size(Int_weights)) :: u_par, u_perp_sq, gamma, pol_fact, c_abs_int
    real(rkind)     :: m_omega_bar, mu, a, N_par, N_perp
    integer(ikind)                :: k
    c_abs= 0.d0
    N_par = N_abs * cos_theta
    N_perp = abs(N_abs * sin_theta)
    m_omega_bar = real(m) / omega_bar
    mu = c0**2 * mass_e/(svec%Te * e0)
    u_par(:) = 1.d0 / sqrt(1.d0 - N_par**2) * ( real(m,8)/m_0 * N_par + &
                       sqrt((real(m) /  m_0)**2 - 1.d0) * Int_absz(:))
    u_perp_sq(:) =    ((real(m,8)/m_0)**2 - 1.d0) *(1.d0 - Int_absz(:)**2)
    gamma(:) = sqrt(1.d0 + u_par(:)**2 + u_perp_sq(:))
    call abs_Al_pol_fact(svec, Int_absz, X, Y, omega_bar, m_0, N_abs, N_par, N_perp, e, 1, m, pol_fact)
    do k = 1, size(Int_weights)
      c_abs_int(k) =  Int_weights(k) * pol_fact(k) * (-mu) * exp(mu * (1.0 - gamma(k)))
      c_abs = c_abs + c_abs_int(k)
    end do
    a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
    c_abs = c_abs * a * (sqrt(mu / (2 * pi))**3)
  end subroutine abs_Al_integral_nume_fast

  subroutine abs_Al_pol_fact(svec, t, X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta, e, mode, m, pol_fact)
  ! According to formula (2a and 2c) in [1]
  ! The 1D ECE is slightly oblique and the Inline and imagining systems are very oblique,
  ! hence the approximation of N_par = 0 is not appropriate
    use mod_ECRad_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)   :: svec
    real(rkind), dimension(:), intent(in)              :: t
    real(rkind), intent(in)                            ::  X, Y, omega_bar, m_0, N_abs, cos_theta, sin_theta
    complex(r8), dimension(:), intent(in)              :: e
    integer(ikind), intent(in)                         :: mode, m
    real(rkind), dimension(:), intent(out)             :: pol_fact
    real(rkind)                   :: x_m
    real(rkind)                   :: N_eff, Axz_sq, Re_Axz_ey, Re_Axz_ez, Re_ey_ez, N_gray, ey_sq, ez_sq, abs_c, N_par, N_perp
    real(rkind), dimension(3,3)  :: E_mat
    complex(r8), dimension(3)  :: pol_vect
    complex(r8)                :: Axz
    real(rkind), dimension(size(t))  :: bessel_arg, bessel_n_l, bessel_n_2, bessel_n, bessel_n_u, abs_Al_bessel_sqr_deriv
    logical                      :: cold, debug
    cold = .true. ! to see if Gray the polarzation vector increases accuracy -> cold = .false.
    debug = .false.
    N_par = N_abs * cos_theta
    N_perp = abs(N_abs * sin_theta)
    x_m =  N_perp * omega_bar * sqrt((real(m,8)/ m_0)**2 - 1.d0)
    N_eff = (N_perp * N_par)/ (1.d0 - N_par**2)
    if(cold) then
      pol_vect = e
	    call get_E_factors(X, Y, N_abs, N_perp, N_par, e,  E_mat )
      ! E_mat (1,2), (2,1) , (3,2) and (2,3) have additional factor of i
      Axz_sq = E_mat(1,1)  + N_eff**2 * E_mat(3,3) + 2.d0 * N_eff * E_mat(1,3) ! mixed terms do not vanish!
      Re_Axz_ey = -E_mat(1,2) - N_eff * E_mat(3,2) ! this includes the additional i
      if(debug) print*, "E_mat(1,2), N_eff * E_mat(3,2)", E_mat(1,2), N_eff * E_mat(3,2)
      Re_Axz_ez = E_mat(1,3) + N_eff * E_mat(3,3)
      Re_ey_ez = -E_mat(3,2)
      ey_sq = E_mat(2,2)
      ez_sq = E_mat(3,3)
    else
      N_gray = N_abs
      abs_c =  abs_Al_Fa_abs(svec, svec%freq_2X * Pi / Y, mode, N_gray, pol_vec =  pol_vect)
    end if
    Axz = e(1) + N_eff * e(3)
    if(debug) print*, "Re_Axz^2 cold", Axz_sq
    Axz_sq = abs(Axz**2 )! mixed terms do not vanish!
    if(debug) print*, "Re_Axz^2 warm", Axz_sq
    if(debug) print*, "Re_Axz_ey cold", Re_Axz_ey
    Re_Axz_ey = real(cmplx(0.d0, 1.d0) * Axz * conjg(pol_vect(2)))
    if(debug) print*, "Re_Axz_ey warm", Re_Axz_ey
    if(debug) print*, "Re_Axz_ez cold", Re_Axz_ez
    Re_Axz_ez = real(Axz * conjg(pol_vect(3)))
    if(debug) print*, "Re_Axz_ez warm", Re_Axz_ez
    if(debug) print*, "Re_ey_ez cold", Re_ey_ez
    Re_ey_ez = real(cmplx(0.0, 1.0) * conjg(pol_vect(2)) * pol_vect(3))
    if(debug) print*, "Re_ey_ez warm", Re_ey_ez
    if(debug) print*, "ey_sq cold", ey_sq
    ey_sq = abs(pol_vect(2))**2
    if(debug) print*, "ey_sq warm", ey_sq
    if(debug) print*, "ez_sq cold", ez_sq
    ez_sq = abs(pol_vect(3))**2
    if(debug) print*, "ez_sq warm", ez_sq
    bessel_arg = x_m * sqrt(1.d0 - t(:)**2)
    bessel_n_l = BesselJ(m - 1 , bessel_arg)
    bessel_n = BesselJ(m , bessel_arg)
    bessel_n_2 = BesselJ(m , bessel_arg)**2
    bessel_n_u = BesselJ(m + 1 , bessel_arg)
    abs_Al_bessel_sqr_deriv = bessel_arg / x_m * bessel_n * ( bessel_n_l - bessel_n_u)
    pol_fact(:) = ( Axz_sq  + ey_sq) * bessel_n_2
    pol_fact(:) = pol_fact(:) + Re_Axz_ey * x_m / real(m,8) * abs_Al_bessel_sqr_deriv
    pol_fact(:) = pol_fact(:) - (bessel_arg / real(m,8))**2 * &
      ey_sq * bessel_n_l * bessel_n_u
    pol_fact(:) = pol_fact(:) + (x_m / (real(m,8) * sqrt( 1.0 - N_par**2)))**2 * &
      ez_sq * t(:)**2 * bessel_n_2
    pol_fact(:) = pol_fact(:) + x_m / (real(m,8) * sqrt( 1.0 - N_par**2)) * &
      2.d0 * Re_Axz_ez * t(:) * bessel_n_2
    pol_fact(:) = pol_fact(:) + x_m / (real(m,8) * sqrt( 1.0 - N_par**2)) * &
      Re_ey_ez * t(:) * x_m / real(m,8) * abs_Al_bessel_sqr_deriv !
    pol_fact(:) = pol_fact(:)  * (real(m,8)  / (N_perp * omega_bar))**2
  end subroutine abs_Al_pol_fact

subroutine get_E_factors(X, Y, N_abs, cos_theta, sin_theta, e, E_mat )
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
!    ! Computeses the matrix:
!    ! (ex ex* ex ey* ex ez*)
!    ! (ey ex* ey ey* ey ez*)
!    ! (ez ex* ez ey* ez ez*)
!    ! From Bornatici Review 1983 p. 1198 eq. 3.1.61 and 3.1.62
!    ! No support for first harmonic nor omega_p > omega_c
!    ! ex := A c ey ! c contains all complex factors
!    ! ez := B c' ey ! c' contains all complex factors
!    ! ey ey^* = 1 / (N sqrt(a_sq + b_sq))
    real(rkind), intent(in) :: X, Y, N_abs,  sin_theta, cos_theta
    complex(r8), dimension(:),  intent(in) :: e
    real(rkind), dimension(:,:), intent(out) :: E_mat
    real(rkind)         :: A, B, ey_sq,rho, f
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 - real(1,8) *  Y* rho)
    A = 1.d0 / Y * (1.d0 - (1 - Y**2) * f)
    B = (N_abs**2 * sin_theta * cos_theta) / ( 1.d0 - X -  N_abs**2 * sin_theta**2) * A
    ey_sq = abs(e(2))**2
     E_mat(1,1) = A**2 ! real
     E_mat(1,2) = A ! imaginary (* i)
     E_mat(1,3) = -A * B ! real
     E_mat(2,1) = -A !imaginary (* i)
     E_mat(2,2) = 1.d0 ! real
     E_mat(2,3) = B !imaginary (* i)
     E_mat(3,1) = -A * B ! real
     E_mat(3,2) = -B !imaginary (* i)
     E_mat(3,3) = B**2 ! real
     E_mat(:,:) = E_mat(:,:) * ey_sq
    if(Any(abs(E_mat(:,:)) > 1.d99) .or. Any(E_mat(:,:) /= E_mat(:,:))) then
      print*, A, Y
      stop "E_mat"
    end if
  end subroutine get_E_factors

  function BesselJ(n,x)
  ! Evaluates the Besselfunction of the first kind
  ! For n = 1, the relativ numerical error for x < 1.5 is below 10^-4
  ! For values in between 1.5 < x < 2.5 the numerical error is below 1 %
  ! For values 2.5 > x > 3.0 the numerical error is below 10 %
  ! For values x > 3 this routine should not be used since the truncation error becomes very large
  ! The larger n, the greater the range where this function gives an almost exact result
  ! e.g. for  n = 2 the numerical error for x = 3.0 is still below 1 %

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind), dimension(:), intent(in)  :: x
    integer(ikind),    intent(in)          :: n
    real(rkind), dimension(size(x))        :: besselj
    BesselJ = Bessel_JN(n,x)
  end function BesselJ

  function BesselJ_custom(n,x)
  ! Evaluates the Besselfunction of the first kind
  ! For n = 1, the relativ numerical error for x < 1.5 is below 10^-4
  ! For values in between 1.5 < x < 2.5 the numerical error is below 1 %
  ! For values 2.5 > x > 3.0 the numerical error is below 10 %
  ! For values x > 3 this routine should not be used since the truncation error becomes very large
  ! The larger n, the greater the range where this function gives an almost exact result
  ! e.g. for  n = 2 the numerical error for x = 3.0 is still below 1 %

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind), dimension(:), intent(in)  :: x
    integer(ikind),    intent(in)          :: n
    real(rkind), dimension(size(x))        :: Besselj_custom
    integer(ikind)                         :: i
    real(rkind)                            :: n_fac
    n_fac = 1.0
    do i= 1,n
      n_fac = real(i) * n_fac
    end do
    BesselJ_custom = 1.0/n_fac
    n_fac = n_fac * real(n + 1)
    BesselJ_custom = BesselJ_custom - x**2/(n_fac * 4.d0)
    n_fac = n_fac * real(n + 2)
    BesselJ_custom = BesselJ_custom + x**4/(n_fac * 32.d0)
    n_fac = n_fac * real(n + 3)
    BesselJ_custom = BesselJ_custom - x**6/(n_fac * 384.0)
    BesselJ_custom = BesselJ_custom * x**n / (2.0)**n
  end function BesselJ_custom

  function get_filter_transmittance(omega, X, Y, cos_theta, sin_theta, mode, x_vec, N_vec, B_vec, x_launch, pol_vec_ext, filter_rot)
  ! Calculates the filtered intensity for a polarization filter aligned with the toroidal direction (e_phi).
  ! Seven steps are required:
  !    1. Calculate polarization vector e in Stix coordinate system
  !    2. Express the polarization vector in the carthesian coordinate system of the ray.
  !    3. The opitcal system of the diagnostic rotates the normalized wave vector N and the polarization vector
  !       such that N -> N' with N' perpendicular to the filter. Hence, step 3 is to determine N' and e'.
  !    4. Express e' in the coordinate system spanned by the passing and the blocking direction of the filter, to determine the Jones Vector.
  !    5. Apply the linear polarizer on the Jones vector.
  !    6. Calculate the Jones vector behing the filter.
  !    7. Calculate the passing intensity.
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ECRad_utils,        only: sub_remap_coords
    use mod_ECRad_types,        only: output_level
    implicit none
    real(rkind), intent(in)               :: omega, X, Y, cos_theta, sin_theta
    integer(ikind), intent(in)            :: mode
    real(rkind), dimension(:), intent(in) :: x_vec, N_vec, B_vec, x_launch
    complex(r8), dimension(:), intent(in), optional :: pol_vec_ext
    real(rkind), intent(in), optional :: filter_rot
    real(rkind)                 :: get_filter_transmittance, N_abs_ray, B_abs, cos_phi, &
                                   sin_phi, N_abs, N_par, N_perp, sigma
    real(rkind), dimension(3)   :: R_vec, f_pass, f_block, f_perp, e_x, e_y, e_z, &
                                   N_vec_norm, norm_vec_N_rot_plane, N_vec_perp_test ! unit vectors of the coordinate system where the polarization vector is defined
    complex(r8), dimension(3)   :: pol_vec, pol_vector_lab, pol_vec_perp
    complex(r8), dimension(2)   :: Jones_vector, filtered_Jones_vector
    complex(r8), dimension(2,2) :: filter_mat
    integer(ikind)              :: i
    logical                     :: debug
    debug = .false.
    call sub_remap_coords(x_launch, R_vec)
    cos_phi = cos(R_vec(2))
    sin_phi = sin(R_vec(2))
    f_pass(:) = 0.d0
    f_block(:) = 0.d0
    f_pass(3) = 1.d0
    f_block(1) = sin_phi
    f_block(2) = cos_phi
    f_perp(1) = f_pass(2) * f_block(3) - f_pass(3) * f_block(2)
    f_perp(2) = f_pass(3) * f_block(1) - f_pass(1) * f_block(3)
    f_perp(3) = f_pass(1) * f_block(2) - f_pass(2) * f_block(1)
    !debug = .false.
    if(present(pol_vec_ext)) then
      pol_vec = pol_vec_ext
      if(debug)  print*, "External polarization vector", pol_vec
    else
      call abs_Al_N_with_pol_vec(X, Y, cos_theta, sin_theta, mode, N_abs, pol_vec)
      if(N_abs == 0.d0) then
        get_filter_transmittance = 0.d0 ! cut off - polarization vector undefined
        return
      end if
      N_perp = abs(N_abs * sin_theta)
      N_par = N_abs * cos_theta
      if(debug) print*, "Internal polarization vector", pol_vec
    end if
    ! To calculate the X/O-mode fraction we need a polarization vector normalized to unity
    pol_vec(:) = pol_vec(:) / sqrt(sum(abs(pol_vec) ** 2))
    if(any(pol_vec /= pol_vec) .and. output_level) then
      print*, "Something wrong with the pol_vec"
      print*, pol_vec
    end if
    if(debug) print*, "pol_vec", pol_vec(:)
    if(debug) print*, "norm pol_vec", sum(abs(pol_vec(:))**2)
    ! First we rotate the polarization vector into the carthesian coordinate system of the machine
    ! Now calculate rotation from Cartesian to reference frame of the wave and the polarization vector.
    ! The polarization coefficients are given in a way so that
    ! vec(N) = (vec(e_x) * sin(theta) + vec(e_z) * cos(theta)) * N_abs
    ! vec(B) = B e_z
    ! Hence, the unit vectors are given by:
    ! vec(e_x) = vec(N) / N_abs - cos(theta) * vec(B) / B_abs
    ! vec(e_y) = vec(e_x) x vec(e_z)
    ! vec(e_z) = vec(B) / B_abs
    ! The rotation is then given according to https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotations_and_motions
    N_abs_ray = sqrt(sum(N_vec**2))
    N_vec_norm = N_vec / N_abs_ray
    B_abs = sqrt(sum(B_vec**2))
    if(debug) then
      print*, "X", X
      print*, "Y", Y
      print*, "omega", omega
      print*, "x_vec", x_vec
      print*, "N_vec", N_vec
      print*, "B_vec", B_vec
      print*, "theta", acos(cos_theta) / pi * 180.d0
    end if
    if(mode > 0 .and. debug) print*, "X-mode"
    if(mode < 0 .and. debug) print*, "O-mode"
    ! N_vec already points towards the antenna when its is copied into svec
    e_x = N_vec_norm - cos_theta  * B_vec / B_abs
    e_x = e_x / sqrt(sum(e_x)**2)
    e_z = B_vec/B_abs
    ! e_y = e_z x e_x
    e_y(1) = e_z(2) * e_x(3) - e_z(3) * e_x(2)
    e_y(2) = e_z(3) * e_x(1) - e_z(1) * e_x(3)
    e_y(3) = e_z(1) * e_x(2) - e_z(2) * e_x(1)
    e_y(:) = e_y(:) / sqrt(sum(e_y**2)) ! not necessary because e_x and e_z perpendicular
    if(debug) then
      print*, "e_x", e_x
      print*, "e_y", e_y
      print*, "e_z", e_z
      print*, "e_x . e_y", sum(e_x * e_y)
      print*, "e_x . e_z", sum(e_x * e_y)
      print*, "e_y . e_z", sum(e_y * e_z)
      print*, "e_x.N_vec", sum(e_x * N_vec_norm)
      print*, "e_y.N_vec", sum(e_y * N_vec_norm)
      print*, "e_z.N_vec", sum(e_z * N_vec_norm)
    end if
    pol_vector_lab(:) = pol_vec(1) * e_x + &
                        pol_vec(2) * e_y + &
                        pol_vec(3) * e_z
    ! Now remove portion that points along N_vec
    pol_vector_lab(:) = pol_vector_lab(:) - N_vec_norm(:) * sum(N_vec_norm(:) * pol_vector_lab(:))
    pol_vector_lab(:) = pol_vector_lab(:) / sqrt(sum(abs(pol_vector_lab(:))**2))
    if(debug) print*, "Polvec in laboratory frame", pol_vector_lab(:)
    if(debug) print*, "norm pol_vec in lab frame", sqrt(sum(abs(pol_vec(:))**2))
    if(debug) print*, "Dot product N_vec pol vec in lab frame", sum(N_vec * pol_vector_lab)
    ! Next the rotation of k by the quasi-optical system:
    ! Normalized vector perpendicular to the filter
    if(debug) print*, "N_vec norm", N_vec_norm
    ! Rotation matrix around angle sigma with axis N_filter x N_vec
    sigma = -acos(sum(N_vec_norm * (-f_perp)))
    if(debug) print*, "Sigma [deg.]", sigma / pi * 180.d0
    norm_vec_N_rot_plane(:) = 0.e0
    ! Rotate the polarization vector in the plane spanned by N_filter and N_vec_norm by the angle sigma
    norm_vec_N_rot_plane(1) = (-f_perp(2)) * N_vec_norm(3) - (-f_perp(3)) * N_vec_norm(2)
    norm_vec_N_rot_plane(2) =  (-f_perp(3)) * N_vec_norm(1) - (-f_perp(1)) * N_vec_norm(3)
    norm_vec_N_rot_plane(3) =  (-f_perp(1)) * N_vec_norm(2) - (-f_perp(2)) * N_vec_norm(1)
    norm_vec_N_rot_plane(:) = norm_vec_N_rot_plane / sqrt(sum(norm_vec_N_rot_plane**2))
    if(debug) print*, "Axis of rotation", norm_vec_N_rot_plane
    call rotate_vec_around_axis(pol_vector_lab, norm_vec_N_rot_plane, sigma, pol_vec_perp)
    if(debug) print*, "pol vec perp filter routine", pol_vec_perp
    if(debug) print*, "dot product rotated N_vec and N_filter vec - should be one", sum(N_vec_perp_test * (-f_perp))
    if(debug) print*, "dot product rotated polarization vector and N_filter vec - should be 0", sum(pol_vec_perp * (-f_perp))
    Jones_vector(1) = sum(f_pass * pol_vec_perp)
    Jones_vector(2) = sum(f_block *  pol_vec_perp)
    if(debug) print*, "Jones_vector:" , Jones_vector
    if(present(filter_rot)) then
      filter_mat(1,1) = cos(filter_rot)**2
      filter_mat(1,2) = cos(filter_rot) * sin(filter_rot)
      filter_mat(2,1) = filter_mat(1,2)
      filter_mat(2,2) = sin(filter_rot)**2
    else
      filter_mat(:,:) = 0.d0
      filter_mat(1,1) = 1.d0
    end if
    do i = 1,2
      filtered_Jones_vector(i) = sum(filter_mat(i,:) * Jones_vector(:))
    end do
    if(debug) print*, "filtered Jones_vector:" , filtered_Jones_vector
    get_filter_transmittance = sum(abs(filtered_Jones_vector)**2)
    if(debug) print*, "Transmittance", get_filter_transmittance
    if(get_filter_transmittance /= get_filter_transmittance .or. get_filter_transmittance < 0.0 .or. &
       get_filter_transmittance > 1.d0) then
      print*, "pol_vec in original coordinates", pol_vec
      print*, "pol_vec in carth. coordinates", pol_vec
      print*, "phi", R_vec(2)
      print*, "cos_theta", cos_theta
      print*, "X", X
      print*, "Y", Y
      print*, "N_abs", N_abs
      print*, "Transmittance", get_filter_transmittance
      stop "Bad value in calculation of polarization filter rejection"
    end if
  end function get_filter_transmittance

end module mod_ECRad_abs_Al

