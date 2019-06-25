! Possible Upgrades:
! -> Introduce beam tracing, see function f_Lambda
! -> Introduce 3D equilibria, see sub_spatial_grad_N_par, sub_grad_N_par


module mod_ecfm_refr_raytrace
    use f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type
#ifdef OMP
    use omp_lib
#endif
    implicit none
    type(plasma_params_type)   :: glob_plasma_params
    real(rkind)   :: glob_omega
    real(rkind), dimension(:), allocatable :: x_loc_vec, N_loc_vec !  dimension(3)
    integer(ikind)          :: glob_mode, debug_level
    !$OMP THREADPRIVATE(glob_omega, &
    !$OMP               x_loc_vec, N_loc_vec, &
    !$OMP               glob_mode, debug_level)
    public ::  span_svecs, dealloc_rad, reinterpolate_svec
    private :: glob_plasma_params, glob_omega, x_loc_vec, N_loc_vec, &
               func_Delta, func_N_s_2, func_N_s_star_2, &
               func_Lambda, func_Lambda_star, &
               func_N_perp_pt, func_H, func_N_s_2_stix, func_N, &
               func_dH_dA, func_dH_dB, func_dH_dC, func_A, func_B, &
               func_C, func_dA_dX, func_dA_dY, func_dB_dX, func_dB_dY, &
               func_dB_dN_par, func_dC_dX, func_dC_dY, func_dC_dN_par, &
               func_dNs_sq_dX, func_dNs_sq_dY2, func_dNs_sq_dN_par, &
               func_dLambda_star_dX, func_dLambda_star_dY2, &
               func_dN_s_star_2_dN_par, func_dR_dx, func_dR_dy, &
               func_dphi_dx, func_dphi_dy, sub_spatial_grad_rhop, &
               func_rhop, func_B_abs, make_B_vec, func_flux_norm_vec, &
               sub_N_par, sub_theta, func_B_r, &
               func_B_x, func_B_y, func_B_z, func_B_x_R, func_B_y_R, &
               func_B_z_R, func_B_r_R, func_B_phi_R, sub_grad_N_par, delta, sub_spatial_grad_X, &
               func_X, sub_spatial_grad_N_par, &
               func_Y, sub_spatial_grad_Y2, sub_spatial_grad_Y, &
               sub_grad_H, sub_grad_Lambda, sub_grad_Lambda_star, f_H, f_Lambda, jac, &
               sub_single_step_LSODE, &
               sub_local_params, reflect_off_surface, make_Snells_refraction, make_H, &
               make_H_wrapper, sub_calculate_initial_N, &
               func_within_plasma, find_first_point_in_plasma, make_ray_segment, &
               make_s_grid, interpolate_svec

    contains
! This module is based on the same routines as D. Farina's code GRAY
! Reference: IFP-CNR Internal Report FP 05/01 (2005)
! To omit diffractive effects S has been set to 1 and therefore k' = grad S = 0


  function func_Delta(X, Y, N_par)
  ! Delta in formula for refractive index
    use f90_kind
    implicit none
    real(rkind), intent(in)       :: X, Y, N_par
    real(rkind)                   :: func_Delta
    func_Delta = (1.d0 - N_par**2)**2 + 4.d0 * N_par**2 * (1.d0 - X) / Y**2
    if(func_Delta < 0.d0) then
      print*, "root of negative number avoided: variable delta"
      stop "Delta in mod_raytrace.f90"
    end if
    func_Delta = sqrt(func_Delta)
  end function func_Delta

  function func_N_s_2(X, Y, N_par, mode)
  ! Calculates the refractive index N_s**2 for either X or O-mode
    use f90_kind
    implicit none
    real(rkind), intent(in)       :: X, Y, N_par
    integer(ikind), intent(in)    :: mode
    real(rkind)                   :: func_N_s_2
    if(mode /= 1 .and. mode /= -1) then
      print*, "Only -1 (X-mode) and 1 (O-Mode) allowed for variable mode in"
      print*, "Function: func_N_s_star_2 in mod_raytrace.f90"
      print*, "mode is", mode
      stop "calculate_N in mod_raytrace.f90"
    end if
    if(X > 1.0) then
      print*, "Ray traced into evanescent region"
      print*, "Subroutine: calculate_N in mod_raytrace.f90"
      call abort()
    end if
    func_N_s_2 = 1.d0 - X + (1.d0 +  real(mode) * func_Delta(X, Y, N_par) + N_par**2)/(2.d0 * (-1.d0 + X + Y**2)) * X * Y**2
  end function func_N_s_2

  function func_N_s_star_2(X, Y, N_par, mode)
  ! Calculates the refractive index N_s**2 for either X or O-mode
    use f90_kind
    implicit none
    real(rkind), intent(in)       :: X, Y, N_par
    integer(ikind), intent(in)    :: mode
    real(rkind)                   :: func_N_s_star_2
    if(mode /= 1 .and. mode /= -1) then
      print*, "Only -1 (X-mode) and 1 (O-Mode) allowed for variable mode in"
      print*, "Function: func_N_s_star_2 in mod_raytrace.f90"
      print*, "mode is", mode
      stop "calculate_N in mod_raytrace.f90"
    end if
    if(X > 1.0) then
      print*, "Ray traced into evanescent region"
      print*, "Subroutine: calculate_N in mod_raytrace.f90"
      stop "calculate_N in mod_raytrace.f90"
    end if
    func_N_s_star_2 =  (1.d0 - X) * (2.d0 * (- 1.d0 + X + Y**2))**2 + &
                       (1.d0 + real(mode) * func_Delta(X, Y, N_par) + N_par**2) * &
                       X * Y**2 * (2.d0 * (- 1.d0 + X + Y**2))
  end function func_N_s_star_2

  function func_Lambda(N_abs, X, Y, N_par, mode)
  ! Caclulates Lambda which is also the Hamiltonian
  ! If even somebody wanted to implemented beam tracing this would be the place
  ! By extending the spline routines and using Mathematica to compute the second derivative of N_s_2 it should be relatively little work
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: N_abs, X, Y ,N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_Lambda
    func_Lambda = N_abs**2  - func_N_s_2(X, Y, N_par, mode)
    !func_Lambda = N_abs**2 * (1.d0 - Y**2) - func_N_s_star_2(X, Y, N_par, mode)**2
  end function func_Lambda


  function func_Lambda_star(N_abs, X, Y, N_par, mode)
  ! Caclulates Lambda which is also the Hamiltonian
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: N_abs, X, Y ,N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_Lambda_star
    func_Lambda_star= N_abs**2 * (2.d0 * (-1.d0 + X + Y**2))**2 - func_N_s_star_2(X, Y, N_par, mode)
    !func_Lambda = N_abs**2 * (1.d0 - Y**2) - func_N_s_star_2(X, Y, N_par, mode)**2
  end function func_Lambda_star


  function func_N_perp_pt(N_abs, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: N_abs, N_par
    real(rkind)                           :: func_N_perp_pt
    func_N_perp_pt = 2.d0 * (N_abs**2 - N_par**2)
  end function func_N_perp_pt

  function func_H(N_perp, A, B, C, mode)
    USE f90_kind
    USE mod_ecfm_refr_types , only : eps
    implicit none
    real(rkind),              intent(in)  :: N_perp, A, B, C
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_H
    func_H = 2.d0 * N_perp**2  + (B - real(mode,8) * sqrt(B**2 - 4.d0 * A * C)) * (A / (A**2 + eps**2))
  end function func_H


  function func_N_s_2_stix(A,  B, C, mode)
    USE f90_kind
    USE mod_ecfm_refr_types , only : eps
    implicit none
    real(rkind),              intent(in)  :: A, B, C
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_N_s_2_stix
    func_N_s_2_stix = (B - real(mode,8) * sqrt(B**2 - 4.d0 * A * C)) * (A / (A**2 + eps**2))
  end function func_N_s_2_stix

  function func_N( X, Y, theta, mode) ! Following eq 13(a -e) of ref [1]
    implicit none
    real(rkind), intent(in)       :: X, Y, theta
    integer(ikind), intent(in)    :: mode
    real(rkind)                   :: func_N
    real(rkind)                   :: sin_theta, cos_theta, rho,f
    sin_theta = sin(theta)
    cos_theta = cos(theta)
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    if(rho < 0.d0) then
      func_N = 0.d0
      return
    end if
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 + real(mode,8) * Y* rho)
    func_N = 1.d0 - X * f
    if(func_N < 0.d0) then
      func_N = 0.d0
      return
    end if
    func_N = sqrt(func_N)
  end function func_N

  function func_dH_dA(A, B, C, mode)
    USE f90_kind
    USE mod_ecfm_refr_types , only : eps
    implicit none
    real(rkind),              intent(in)  :: A, B, C
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dH_dA
    func_dH_dA =(-2.d0*A**3*C*real(mode,8) + 6.d0*A*C*eps**2*real(mode,8) + &
                B*eps**2*(Sqrt(B**2 - 4.d0*A*C) - B*real(mode,8)) + &
                A**2*B*(-Sqrt(B**2 - 4.d0*A*C) + B*real(mode,8)))/ &
                (Sqrt(B**2 - 4.d0*A*C)*(A**2 + eps**2)**2)
  end function func_dH_dA

  function func_dH_dB(A, B, C, mode)
    USE f90_kind
    USE mod_ecfm_refr_types , only : eps
    implicit none
    real(rkind),              intent(in)  :: A, B, C
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dH_dB
    func_dH_dB = (A - (A*B*real(mode,8))/Sqrt(B**2 - 4.d0*A*C))/(A**2 + eps**2)
  end function func_dH_dB
!
  function func_dH_dC(A, B, C, mode)
    USE f90_kind
    USE mod_ecfm_refr_types , only : eps
    implicit none
    real(rkind),              intent(in)  :: A, B, C
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dH_dC
    func_dH_dC = (2.d0*A**2*real(mode,8))/(Sqrt(B**2 - 4.d0*A*C)*(A**2 + eps**2))
  end function func_dH_dC

  function func_A(X, Y)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y
    real(rkind)                           :: func_A
    func_A = (-1.d0 + X + Y**2)/(-1.d0 + Y**2)
  end function func_A

  function func_B(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y, N_par
    real(rkind)                           :: func_B
    func_B =  (X**2*Y**2)/(1 - Y**2)**2 - (1 - X)*(1 - X/(1 - Y**2)) - &
              (1 - X/(1 - Y**2))**2 + N_par**2*(2 - X - X/(1 - Y**2))
  end function func_B

  function func_C(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y, N_par
    real(rkind)                           :: func_C
    func_C =  ((-1.d0 + X)*(1.d0 - 2*X + X**2 - Y**2 - N_par**4*(-1.d0 + Y**2) + &
              2.d0*N_par**2*(-1.d0 + X + Y**2)))/(-1.d0 + Y**2)
  end function func_C

  function func_dA_dX(X, Y)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y
    real(rkind)                           :: func_dA_dX
    func_dA_dX = 1.d0/(-1.d0 + Y**2)
  end function func_dA_dX

function func_dA_dY(X, Y)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y
    real(rkind)                           :: func_dA_dY
    func_dA_dY =  (-2.d0*X*Y)/(-1.d0 + Y**2)**2
  end function func_dA_dY

  function func_dB_dX(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  ::  X, Y, N_par
    real(rkind)                           :: func_dB_dX
    func_dB_dX = (-4.d0 + 4.d0*X + Y**2 - N_par**2*(-2.d0 + Y**2))/(-1.d0 + Y**2)
  end function func_dB_dX

  function func_dB_dY(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  ::  X, Y, N_par
    real(rkind)                           :: func_dB_dY
    func_dB_dY = (-2.d0*X*(-3.d0 + N_par**2 + 2.d0*X)*Y)/(-1.d0 + Y**2)**2
  end function func_dB_dY

  function func_dB_dN_par(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y, N_par
    real(rkind)                           :: func_dB_dN_par
    func_dB_dN_par = 2.d0*N_par*(2d0 + X*(-1d0 + 1d0/(-1d0 + Y**2)))

  end function func_dB_dN_par

    function func_dC_dX(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  ::  X, Y, N_par
    real(rkind)                           :: func_dC_dX
    func_dC_dX = -((-3.d0 + 6.d0*X - 3.d0*X**2 + Y**2 + N_par**4*(-1.d0 + Y**2) - &
                  2.d0*N_par**2*(-2.d0 + 2.d0*X + Y**2))/(-1.d0 + Y**2))
  end function func_dC_dX

  function func_dC_dY(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  ::  X, Y, N_par
    real(rkind)                           :: func_dC_dY
    func_dC_dY = (-2.d0*(-1.d0 + X)*X*(-2.d0 + 2.d0*N_par**2 + X)*Y)/(-1.d0 + Y**2)**2
  end function func_dC_dY

  function func_dC_dN_par(X, Y, N_par)
    USE f90_kind
    implicit none
    real(rkind),              intent(in)  :: X, Y, N_par
    real(rkind)                           :: func_dC_dN_par
    func_dC_dN_par = (4.d0*N_par*(-1.d0 + X)*(-1.d0 + X + Y**2 - N_par**2*(-1.d0 + Y**2)))/(-1.d0 + Y**2)
  end function func_dC_dN_par

  function func_dNs_sq_dX(N_abs, X, Y, N_par, mode)
  ! Calculates dN_s/ dX
    !USE mod_ecfm_refr_types , only : h_x_glob
    USE f90_kind
    implicit None
    real(rkind),              intent(in)  :: N_abs, X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dNs_sq_dX
    func_dNs_sq_dX = (-2.d0 - (2.d0*real(mode,8)*N_par**2*X)/((-1.d0 + X + Y**2)*func_delta(X,Y,N_par)) - &
                           (X*Y**2*(1.d0  + N_par**2 + real(mode,8)*func_delta(X,Y,N_par)))/ &
                            (-1.d0  + X + Y**2)**2 + &
                              (Y**2*(1.d0  + N_par**2 + real(mode,8)*func_delta(X,Y,N_par)))/(-1.d0 + X + Y**2.d0))/2.d0
  end function func_dNs_sq_dX

 function func_dNs_sq_dY2(N_abs, X, Y, N_par, mode)
  ! Calculates dNs_sq/dY^2
    USE f90_kind
    !USE mod_ecfm_refr_types , only : h_x_glob
    implicit None
    real(rkind),              intent(in)  :: N_abs, X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dNs_sq_dY2
    func_dNs_sq_dY2 =((-1.d0 + X)*X*(real(mode,8)*(-2.d0*N_par**2*(-1.d0 + X) + Y**2 + N_par**4*Y**2) + &
                            (1.d0 + N_par**2)*Y**2*Sqrt(1.d0 + N_par**4 - (2.d0*N_par**2*(-2.d0 + 2.d0*X + Y**2))/Y**2)))/ &
                            (Y*(-1.d0 + X + Y**2)**2*func_delta(X,Y,N_par))
    func_dNs_sq_dY2 = func_dNs_sq_dY2 / (2.d0 * Y)
  end function func_dNs_sq_dY2

  function func_dNs_sq_dN_par(N_abs, X, Y, N_par, mode)
  ! Calculates dN_s/ d N_par
    USE f90_kind
    implicit None
    real(rkind),              intent(in)  :: N_abs, X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dNs_sq_dN_par
    func_dNs_sq_dN_par =  (N_par*X*Y**2*(1 + (real(mode,8)*(2 - 2*X + (-1 + N_par**2)*Y**2))/ &
                               (Y**2*func_delta(X,Y,N_par))))/(-1.d0 + X + Y**2)
    !func_dNs_sq_dN_par =  -(real(mode,8)*X*Y**2*(-4.d0*N_par*(1.d0 - N_par**2) + (8*N_par*(1.d0 - X))/Y**2))/ &
    !  (4.d0*Sqrt(1.d0 - N_par**2)**2 + (4.d0*N_par**2*(1.d0 - X))/Y**2)* &
    !  (1.d0 - X - Y**2)
  end function func_dNs_sq_dN_par

  function func_dLambda_star_dX(N_abs, X, Y, N_par, mode)
  ! Calculates dN_s/ dX
    !USE mod_ecfm_refr_types , only : h_x_glob
    USE f90_kind
    implicit None
    real(rkind),              intent(in)  :: N_abs,  X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dLambda_star_dX
    !func_dLambda_star_dX = - 2.d0 * N_abs**2  - (2.d0 * (-1.d0 + X)+(2.d0 * real(mode,8) * N_par**2 * X)/func_Delta(X,Y,N_par) - &
    !               (1.d0 + N_par**2 + real(mode,8) * func_Delta(X,Y,N_par)) * Y**2 + 2.d0 * (-1.d0 + X + Y**2))
    func_dLambda_star_dX = 8.d0*N_abs**2*(-1.d0 + X + Y**2) + 8.d0*(-1.d0 + X)*(-1.d0 + X + Y**2) + &
                           4.d0*(-1.d0 + X + Y**2)**2 + (4.d0*real(mode,8)*N_par**2*X*(-1.d0 + X + Y**2))/ &
                           func_delta(X,Y,N_par) - &
                           2.d0*X*Y**2*(1.d0 + N_par**2 + real(mode,8)* &
                           func_delta(X,Y,N_par)) - &
                           2.d0*Y**2*(-1.d0 + X + Y**2)*(1.d0 + N_par**2 + &
                           real(mode,8)*func_delta(X,Y,N_par))
  end function func_dLambda_star_dX


  function func_dLambda_star_dY2(N_abs, X, Y, N_par, mode)
    USE f90_kind
    !USE mod_ecfm_refr_types , only : h_x_glob
    implicit None
    real(rkind),              intent(in)  :: N_abs, X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dLambda_star_dY2
    !func_dLambda_star_dY2 =-((4.d0 * real(mode, 8) * N_par**2 * (-1.d0 + X) * X)/(func_Delta(X,Y,N_par)  * Y)) + &
    !                 4.d0 * (-1.d0 + X) * Y - 2.d0 * X * (1.d0 + N_par**2 + real(mode, 8) * func_Delta(X,Y,N_par) ) * Y
    !func_dLambda_star_dY2 =  - 2.d0 * N_abs**2 -func_dLambda_star_dY2 / (2.d0 * Y)
    func_dLambda_star_dY2 = 8.d0*N_abs**2*(-1.d0 + X + Y**2) + 8.d0*(-1.d0 + X)*(-1.d0 + X + Y**2) - &
                            (4*real(mode,8)*N_par**2*(-1.d0 + X)*X*(-1.d0 + X + Y**2))/ &
                            (Y**2*func_delta(X,Y,N_par)) - &
                            2.d0*X*Y**2*(1.d0 + N_par**2 + real(mode,8)* &
                            func_delta(X,Y,N_par)) - &
                            2.d0*X*(-1.d0 + X + Y**2)*(1.d0 + N_par**2 + &
                           real(mode,8)*func_delta(X,Y,N_par))
  end function func_dLambda_star_dY2

  function func_dN_s_star_2_dN_par(N_abs,X, Y, N_par, mode)
  ! Calculates dN_s/ d N_par
    USE f90_kind
    implicit None
    real(rkind),              intent(in)  :: N_abs, X, Y, N_par
    integer(ikind),           intent(in)  :: mode
    real(rkind)                           :: func_dN_s_star_2_dN_par
!    func_dN_s_star_2_dN_par=  - (-2.d0 * N_par * X  * Y**2 * (1.d0 + real(mode,8) * &
!                       (2.d0 - 2.d0 * X +(N_par**2 -1.d0) * Y**2)/(func_Delta(X,Y,N_par) * Y**2)))
    func_dN_s_star_2_dN_par = 2.d0*X*Y**2*(-1.d0 + X + Y**2)*(2*N_par + &
                              (2.d0*real(mode,8)*N_par*(2.d0 - 2.d0*X + (-1.d0 + N_par**2)*Y**2))/ &
                              (Y**2*func_delta(X,Y,N_par)))
  end function func_dN_s_star_2_dN_par

  function func_dR_dx(x,y)
  use f90_kind
  implicit none
    real(rkind), intent(in)    :: x, y
    real(rkind)                :: func_dR_dx
    func_dR_dx = x / sqrt(x**2 + y**2)
  end function func_dR_dx

  function func_dR_dy(x,y)
    use f90_kind
    implicit none
    real(rkind), intent(in)    :: x, y
    real(rkind)                :: func_dR_dy
    func_dR_dy = y / sqrt(x**2 + y**2)
  end function func_dR_dy

  function func_dphi_dx(x,y)
    use f90_kind
    implicit none
    real(rkind), intent(in)    :: x, y
    real(rkind)                :: func_dphi_dx
    real(rkind)                :: h_x
    h_x = 1.d-4
    if(abs(y) == 1.d-17) then
      func_dphi_dx = 0.d0
    else
      func_dphi_dx = (y*(-1.d0 + x/Sqrt(x**2 + y**2)))/(x**2 + y**2 - x*Sqrt(x**2 + y**2))
    end if
    !print*,"ana dphi_dx",func_dphi_dx
    !print*, "num dphi_dx",(- func_calc_phi(x + 2.d0 * h_x, y) + &
    !                              8.d0 *  func_calc_phi(x + 1.d0 * h_x, y) - &
    !                          8.d0 *  func_calc_phi(x - 1.d0 * h_x, y)  &
    !                          +   func_calc_phi(x - 2.d0 * h_x, y )) / (12.d0 *h_x)
  end function func_dphi_dx

  function func_dphi_dy(x,y)
    use f90_kind
    implicit none
    real(rkind), intent(in)    :: x, y
    real(rkind)                :: func_dphi_dy
    real(rkind)                :: h_x
    h_x = 1.d-4
    if(abs(y) == 1.d-17) then
      func_dphi_dy = 0.d0
    else
      func_dphi_dy =  (x - x**2/Sqrt(x**2 + y**2))/(x**2 + y**2 - x*Sqrt(x**2 + y**2))
    end if
    !print*,"ana dphi_dy",func_dphi_dy
    !print*, "num dphi_dy",(- func_calc_phi(x, y + 2.d0 * h_x) + &
    !                              8.d0 *  func_calc_phi(x, y + 1.d0 * h_x) - &
    !                          8.d0 *  func_calc_phi(x, y - 1.d0 * h_x)  &
    !                          +   func_calc_phi(x, y - 2.d0 * h_x)) / (12.d0 *h_x)
    !stop "derivatives"
  end function func_dphi_dy


  subroutine sub_spatial_grad_rhop(plasma_params, x_vec, rhop, spatial_grad_rhop)
    USE f90_kind
    USE mod_ecfm_refr_types , only: plasma_params_type
    USE ripple3d,                 only: grad_type
    use mod_ecfm_refr_utils,      only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    !USE nag_spline_2d             , only: nag_spline_2d_eval, &
    !                                      nag_error, nag_set_error
    implicit none
    type(plasma_params_type)                   :: plasma_params
    real(rkind), dimension(:)  , intent(in)    :: x_vec
    real(rkind)                , intent(out)   :: rhop
    real(rkind), dimension(3)  , intent(out)   :: spatial_grad_rhop
    type(grad_type)                            :: R_grad_rhop
    real(rkind), dimension(3)                  :: R_vec
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) - plasma_params%R_shift > plasma_params%R_max .or. R_vec(1) - plasma_params%R_shift < plasma_params%R_min .or. &
       R_vec(3) - plasma_params%z_shift > plasma_params%z_max .or. R_vec(3) - plasma_params%z_shift < plasma_params%z_min) then
       rhop = -1.d0
       spatial_grad_rhop(:) = 0.d0
       return
    end if
    call rect_spline(plasma_params%rhop_spline, R_vec(1) - plasma_params%R_shift, &
                     R_vec(3) - plasma_params%z_shift, rhop,R_grad_rhop%dR,R_grad_rhop%dz)
    if(rhop /= rhop) then
      rhop = -1.d0 !plasma_params%rhop_max + h_x_glob
      spatial_grad_rhop(:) = 0.d0
      return
    end if
    spatial_grad_rhop(1) = func_dR_dx(x_vec(1),x_vec(2)) * R_grad_rhop%dR
    spatial_grad_rhop(2) = func_dR_dy(x_vec(1),x_vec(2)) * R_grad_rhop%dR
    spatial_grad_rhop(3) = R_grad_rhop%dz
    !call bilin_inter_regular(plasma_params, R_vec, plasma_params%rhop, rhop)
    !print*,"rhop_bilin", rhop
    !rhop = func_n_e(plasma_params, x_vec)
  end subroutine sub_spatial_grad_rhop


  function func_rhop(plasma_params, x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    USE ripple3d,                 only: grad_type
    USE mod_ecfm_refr_utils,        only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    type(plasma_params_type)                  :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3)                 :: R_vec
    real(rkind)                               :: func_rhop
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) - plasma_params%R_shift > plasma_params%R_max .or. R_vec(1) - plasma_params%R_shift < plasma_params%R_min .or. &
       R_vec(3) - plasma_params%z_shift > plasma_params%z_max .or. R_vec(3) - plasma_params%z_shift < plasma_params%z_min) then
       func_rhop =  -1.d0
       return
    end if
#ifdef NAG
    if(debug_level == 0 .or. .not. output_level) then
      call rect_spline(plasma_params%rhop_spline, R_vec(1) - plasma_params%R_shift, &
                R_vec(3) - plasma_params%z_shift, &
                func_rhop)
    else
      call rect_spline(plasma_params%rhop_spline, R_vec(1) - plasma_params%R_shift, &
              R_vec(3)  - plasma_params%z_shift, &
              func_rhop, nag_spline=plasma_params%rhop_spline_nag)
    end if
#else
    call rect_spline(plasma_params%rhop_spline, R_vec(1) - plasma_params%R_shift, &
                R_vec(3) - plasma_params%z_shift, &
                func_rhop)
#endif
    if(func_rhop /= func_rhop) then
      print*, "warning nan in rhop encountered"
      func_rhop = -1.d0!plasma_params%rhop_max + h_x_glob
    end if
  end function func_rhop

  function func_B_abs(plasma_params, x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    USE ripple3d,                 only: get_ripple, grad_type
    USE mod_ecfm_refr_utils,      only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type)                  :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind)                               :: func_B_abs
    real(rkind), dimension(3)                 :: R_vec, B_r_vec,B_ripple
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       func_B_abs = 0.d0
       return
    else
#ifdef NAG
      if(debug_level == 0 .or. .not. output_level) then
        call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1))
        call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                    B_R_vec(2))
        call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                    B_R_vec(3))
      else
        call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1), nag_spline=plasma_params%B_r_spline_nag)
        call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                    B_R_vec(2), nag_spline=plasma_params%B_t_spline_nag)
        call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                    B_R_vec(3), nag_spline=plasma_params%B_z_spline_nag)
      end if
#else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1))
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3))
#endif
    end if
    if(plasma_params%w_ripple) then
      call get_ripple(R_vec, B_ripple)
      B_R_vec = B_R_vec + B_ripple
    end if
    func_B_abs = sqrt(B_R_vec(1)**2 + B_R_vec(2)**2 + B_R_vec(3)**2)
  end function func_B_abs

  subroutine make_B_vec(plasma_params, x_vec, B_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    USE ripple3d,                 only: get_ripple, grad_type
    USE mod_ecfm_refr_utils,      only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type)                  :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(:)  , intent(out)  :: B_vec
    real(rkind), dimension(3)                 :: R_vec, B_r_vec,B_ripple
    real(rkind)                               :: cos_phi, sin_phi
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       B_R_vec(1) = 0.d0
       B_R_vec(2) = 1.d0 ! purely toroidal field
       B_R_vec(3) = 0.d0
       return
    else
#ifdef NAG
      if(debug_level == 0 .or. .not. output_level)  then
        call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1))
        call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                    B_R_vec(2))
        call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                    B_R_vec(3))
      else
        call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1), nag_spline=plasma_params%B_r_spline_nag)
        call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                    B_R_vec(2), nag_spline=plasma_params%B_t_spline_nag)
        call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                    B_R_vec(3), nag_spline=plasma_params%B_z_spline_nag)
      end if
#else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                    B_R_vec(1))
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3))
#endif
      if(plasma_params%w_ripple) then
        call get_ripple(R_vec, B_ripple)
        B_R_vec = B_R_vec + B_ripple
      end if
    end if
    cos_phi = cos(R_vec(2))
    sin_phi = sin(R_vec(2))
    B_vec(1) = cos_phi * B_R_vec(1) - sin_phi * B_R_vec(2)
    B_vec(2) = sin_phi * B_R_vec(1) + cos_phi * B_R_vec(2)
    B_vec(3) = B_R_vec(3)
  end subroutine make_B_vec

  function func_flux_norm_vec(plasma_params, x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    USE ripple3d,                 only: get_ripple, grad_type
    USE mod_ecfm_refr_utils,      only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:), intent(in)     :: x_vec
    real(rkind), dimension(3)                 :: func_flux_norm_vec
    real(rkind), dimension(3)                 :: R_vec, B_ripple
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       func_flux_norm_vec(1) = 0.d0
       func_flux_norm_vec(2) = 1.d0 ! purely toroidal field
       func_flux_norm_vec(3) = 0.d0
       return
    else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  func_flux_norm_vec(3))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  func_flux_norm_vec(1))
    end if
    func_flux_norm_vec(2) = 0.d0 ! Assume no gradient in phi direction except ripple
    ! This is only correct for tokamaks !
    if(plasma_params%w_ripple) then
      call get_ripple(R_vec, B_ripple)
      func_flux_norm_vec(3) = func_flux_norm_vec(3) + B_ripple(1)
      func_flux_norm_vec(1) = func_flux_norm_vec(1) + B_ripple(3)
      func_flux_norm_vec(2) = B_ripple(2)
    end if
    func_flux_norm_vec(3) = - func_flux_norm_vec(3)
  end function func_flux_norm_vec

  subroutine sub_N_par(plasma_params, x_vec, N_vec, N_par, N_abs, B_abs)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type,get_ripple
    Use mod_ecfm_refr_utils,  only : sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec, N_vec
    real(rkind), intent(out)                  :: N_par, N_abs, B_abs
    real(rkind), dimension(3)                 :: R_vec, B_R_vec, B_x_vec, B_ripple
    real(rkind)                               :: cos_phi, sin_phi, scal_prod
    integer(ikind)                            :: i
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       N_par = 0.01
       return
    end if
#ifdef NAG
    if(debug_level == 0 .or. .not. output_level) then
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3))
    else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1), nag_spline=plasma_params%B_r_spline_nag)
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2), nag_spline=plasma_params%B_t_spline_nag)
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3), nag_spline=plasma_params%B_z_spline_nag)
    end if
#else
    call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
    call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                B_R_vec(2))
    call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                B_R_vec(3))
#endif
    if(plasma_params%w_ripple) then
      call get_ripple(R_vec, B_ripple)
      B_R_vec = B_R_vec + B_ripple
    end if
    cos_phi = cos(R_vec(2))
    sin_phi = sin(R_vec(2))
    B_x_vec(1) = cos_phi * B_R_vec(1) - sin_phi * B_R_vec(2)
    B_x_vec(2) = sin_phi * B_R_vec(1) + cos_phi * B_R_vec(2)
    B_x_vec(3) = B_R_vec(3)
!    B_x_vec = B_r_vec
    N_abs = sqrt(N_vec(1)**2 + N_vec(2)**2 + N_vec(3)**2)
    B_abs = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
    scal_prod = 0.d0
    do i = 1, 3
      scal_prod = scal_prod + N_vec(i) * B_x_vec(i)
    end do
    N_par = scal_prod / B_abs
  end subroutine sub_N_par

  subroutine sub_theta(plasma_params, x_vec, N_vec, theta, N_abs, B_abs, B_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type, get_ripple
    use mod_ecfm_refr_utils, only  : sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    use constants                     , only: pi
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec, N_vec
    real(rkind), intent(out)                  :: theta, N_abs, B_abs
    real(rkind), dimension(3), intent(out)    :: B_vec
    real(rkind), dimension(3)                 :: R_vec, B_R_vec, B_x_vec, B_ripple
    real(rkind)                               :: cos_phi, sin_phi, scal_prod
    integer(ikind)                            :: i
    call sub_remap_coords(x_vec, R_vec)
    cos_phi = cos(R_vec(2))
    sin_phi = sin(R_vec(2))
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       theta = 0.d0
       N_abs = 1.0
       B_abs = 0.d0 !plasma_params%B_ax * plasma_params%R_ax / R_vec(1) ! Not entirely correct since diamagnetic field included
                                                                  ! However, this is just an output and does not affect the computation
!       B_vec(1) = - sin_phi * B_abs
!       B_vec(2) =   cos_phi * B_abs
       B_vec(:) = 0.d0
       return
    end if
#ifdef NAG
    if(debug_level == 0 .or. .not. output_level)  then
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3))
    else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1), nag_spline=plasma_params%B_r_spline_nag)
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2), nag_spline=plasma_params%B_t_spline_nag)
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3), nag_spline=plasma_params%B_z_spline_nag)
    end if
#else
    call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
    call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                B_R_vec(2))
    call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                B_R_vec(3))
#endif
    if(plasma_params%w_ripple) then
      call get_ripple(R_vec, B_ripple)
      B_R_vec = B_R_vec + B_ripple
    end if
    B_x_vec(1) = cos_phi * B_R_vec(1) - sin_phi * B_R_vec(2)
    B_x_vec(2) = sin_phi * B_R_vec(1) + cos_phi * B_R_vec(2)
    B_x_vec(3) = B_R_vec(3)
    !print*, "R_vec(2), B_vec",R_vec(2) / pi * 180.d0, B_x_vec, B_R_vec
    B_vec = B_x_vec
    N_abs = sqrt(N_vec(1)**2 + N_vec(2)**2 + N_vec(3)**2)
    B_abs = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
!    if(abs(sqrt(B_R_vec(1)**2 + B_R_vec(2)**2 + B_R_vec(3)**2) - B_abs) > 1.d-6) then
!      print*, "Large deviation between the Bs", sqrt(B_R_vec(1)**2 + B_R_vec(2)**2 + B_R_vec(3)**2) - B_abs
!      stop "Calculation of B faulty"
!    end if
    scal_prod = 0.d0
    do i = 1, 3
      scal_prod = scal_prod + N_vec(i) * B_x_vec(i)
    end do
    theta = acos(scal_prod / (B_abs * N_abs))
  end subroutine sub_theta

  subroutine func_B_r(plasma_params, R_vec, B_R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type, get_ripple
#ifdef NAG
    USE nag_spline_2d             , only: nag_spline_2d_eval
#endif
   use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(:)  , intent(out)  :: B_R_vec
    real(rkind), dimension(3)  :: B_ripple
#ifdef NAG
    if(debug_level == 0 .or. .not. output_level)  then
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2))
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3))
    else
      call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1), nag_spline=plasma_params%B_r_spline_nag)
      call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                  B_R_vec(2), nag_spline=plasma_params%B_t_spline_nag)
      call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                  B_R_vec(3), nag_spline=plasma_params%B_z_spline_nag)
    end if
#else
    call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                  B_R_vec(1))
    call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                B_R_vec(2))
    call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                B_R_vec(3))
#endif
    if(plasma_params%w_ripple) then
      call get_ripple(R_vec, B_ripple)
      B_R_vec = B_R_vec + B_ripple
    end if
    if(debug_level == 0) return
  end subroutine func_B_r

  function  func_B_x(plasma_params, x_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3) :: R_vec, B_R_vec
    real(rkind)                :: func_B_x
    call sub_remap_coords(x_vec, R_vec)
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_x =  B_R_vec(1) * cos(R_vec(2)) - B_R_vec(2) * sin(R_vec(2))
  end function func_B_x

  function  func_B_y(plasma_params, x_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3) :: R_vec, B_R_vec
    real(rkind)                :: func_B_y
    call sub_remap_coords(x_vec, R_vec)
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_y =  B_R_vec(1) * sin(R_vec(2)) + B_R_vec(2) * cos(R_vec(2))
  end function func_B_y

    function  func_B_z(plasma_params, x_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3) :: R_vec, B_R_vec
    real(rkind)                :: func_B_z
    call sub_remap_coords(x_vec, R_vec)
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_z =  B_R_vec(3)
  end function func_B_z

 function  func_B_x_R(plasma_params, R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3) :: B_R_vec
    real(rkind)               :: func_B_x_R
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_x_R =  B_R_vec(1) * cos(R_vec(2)) - B_R_vec(2)* sin(R_vec(2))
  end function func_B_x_R

  function  func_B_y_R(plasma_params, R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3)  :: B_R_vec
    real(rkind)                :: func_B_y_R
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_y_R =  B_R_vec(1) * sin(R_vec(2)) + B_R_vec(2)* cos(R_vec(2))
  end function func_B_y_R

  function  func_B_z_R(plasma_params, R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3)  :: B_R_vec
    real(rkind)                :: func_B_z_R
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_z_R =  B_R_vec(3)
  end function func_B_z_R

  function  func_B_r_R(plasma_params, R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,         only     : grad_type
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3)  :: B_R_vec
    real(rkind)                :: func_B_r_R
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_r_R =  B_R_vec(1)
  end function func_B_r_R

  function  func_B_phi_R(plasma_params, R_vec)
    use f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3)  :: B_R_vec
    real(rkind)                :: func_B_phi_R
    call func_B_r(plasma_params, R_vec, B_R_vec)
    func_B_phi_R =  B_R_vec(2)
  end function func_B_phi_R


  subroutine sub_grad_N_par(plasma_params, x_vec, N_vec, N_abs, B_abs, spatial_grad_B_abs, &
                           N_par, N_grad_N_par, spatial_grad_N_par) !spatial_grad_B_x, spatial_grad_B_y, spatial_grad_B_z,
  ! Gradient of vec(B) along LOS coordinates x ripple not (yet) included => dB_vec/dphi = 0
  ! Also calculates N_par, d N_par(theta)/dN_i, d N_par/dx_i, d|B|dx_i since all information is readily available
  ! (Calculating these quantities here safes interpolations)
  ! This is the routine that would need to be extended to allow for 3D equilibria
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,         only : grad_type, get_ripple_w_grad
    USE mod_ecfm_refr_utils, only: sub_remap_coords
    use mod_ecfm_refr_interpol,      only: rect_spline
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec, N_vec
    real(rkind)                , intent(in)   :: N_abs
    real(rkind), dimension(:),  intent(out)   :: spatial_grad_B_abs
    real(rkind)                , intent(out)  :: B_abs, N_par
    real(rkind), dimension(:)  , intent(out)  :: N_grad_N_par, spatial_grad_N_par
    real(rkind), dimension(3, 3)              :: spatial_grad_B
    real(rkind)                               :: cos_phi_tok, sin_phi_tok, scal_prod, h_x
    type(grad_type)                           :: dB_r_inter, dB_t_inter, dB_z_inter, dB_ripple_r, &
                                                 dB_ripple_t, dB_ripple_z
    real(rkind), dimension(3)                 :: R_vec, B_R_vec, B_x_vec, B_ripple!, N_vec_norm, B_vec_norm
    real(rkind), dimension(4, 3)              :: aux_x, aux_R
    real(rkind), dimension(3)                 :: dB_x_dR, dB_y_dR, dB_z_dR
    integer(ikind)                            :: i, j,k
    h_x = 1.d-6!h_x_glob
    call sub_remap_coords(x_vec, R_vec)
    if(R_vec(1) > plasma_params%R_max .or. R_vec(1) < plasma_params%R_min .or. &
       R_vec(3) > plasma_params%z_max .or. R_vec(3) < plasma_params%z_min) then
       B_abs = 0.d0
       spatial_grad_B_abs = 0.d0
       N_par = 0.01
       N_grad_N_par = 0.d0
       spatial_grad_N_par = 0.d0
       return
    end if
    call rect_spline(plasma_params%B_r_spline, R_vec(1), R_vec(3), &
                B_R_vec(1), dB_r_inter%dR, dB_r_inter%dz)
    call rect_spline(plasma_params%B_t_spline, R_vec(1), R_vec(3), &
                B_R_vec(2), dB_t_inter%dR, dB_t_inter%dz)
    call rect_spline(plasma_params%B_z_spline, R_vec(1), R_vec(3), &
                B_R_vec(3), dB_z_inter%dR, dB_z_inter%dz)
    if(plasma_params%w_ripple) then
      call get_ripple_w_grad(R_vec, B_ripple, dB_ripple_r, dB_ripple_t, dB_ripple_z)
      B_R_vec = B_R_vec + B_ripple
      dB_r_inter%dR = dB_r_inter%dR  + dB_ripple_r%dR
      dB_r_inter%dphi = 0.d0  + dB_ripple_r%dphi ! dB_r_inter%dphi, replace 0.0 with the derivatives of B in the toroidal direction1
      dB_r_inter%dz = dB_r_inter%dz  + dB_ripple_r%dz
      dB_t_inter%dR = dB_t_inter%dR  + dB_ripple_t%dR
      dB_t_inter%dphi = 0.d0  + dB_ripple_t%dphi ! dB_t_inter%dphi, replace 0.0 with the derivatives of B in the toroidal direction1
      dB_t_inter%dz = dB_t_inter%dz  + dB_ripple_t%dz
      dB_z_inter%dR = dB_z_inter%dR  + dB_ripple_z%dR
      dB_z_inter%dphi = 0.d0 + dB_ripple_z%dphi ! dB_z_inter%dphi, replace 0.0 with the derivatives of B in the toroidal direction1
      dB_z_inter%dz = dB_z_inter%dz  + dB_ripple_z%dz
    else
      dB_r_inter%dphi = 0.d0 ! replace 0.0 with the derivatives of B in the toroidal direction1
      dB_t_inter%dphi = 0.d0 ! replace 0.0 with the derivatives of B in the toroidal direction1
      dB_z_inter%dphi = 0.d0 ! replace 0.0 with the derivatives of B in the toroidal direction1
    end if
    cos_phi_tok = cos(R_vec(2))
    sin_phi_tok = sin(R_vec(2))
    B_x_vec(1) = B_R_vec(1) * cos_phi_tok - B_R_vec(2) * sin_phi_tok
    B_x_vec(2) = B_R_vec(1) * sin_phi_tok + B_R_vec(2) * cos_phi_tok
    B_x_vec(3) = B_R_vec(3)
    ! Next dB_x_j/d_x_i --> j (first index) direction of magnetic field
    !                   --> i (second index) direction of the derivative
    ! To do this analyitaclly the multidimensional chain rule is used
    ! First dvec(B_x)dR and dvec(B_x)dz are computed
    dB_x_dR(1) = dB_r_inter%dR * cos_phi_tok - dB_t_inter%dR * sin_phi_tok
    ! dB_x/dphi
    dB_x_dR(2) = -B_R_vec(1) * sin_phi_tok - B_R_vec(2) * cos_phi_tok + &
                  dB_r_inter%dphi * cos_phi_tok - dB_t_inter%dphi * sin_phi_tok
    ! dB_x/dz
    dB_x_dR(3) = dB_r_inter%dz * cos_phi_tok - dB_t_inter%dz * sin_phi_tok
    ! dB_y/dR
    dB_y_dR(1) = dB_r_inter%dR * sin_phi_tok + dB_t_inter%dR * cos_phi_tok
    ! dB_y/dphi
    dB_y_dR(2) = B_R_vec(1) * cos_phi_tok - B_R_vec(2) * sin_phi_tok + &
                 dB_r_inter%dphi * sin_phi_tok + dB_t_inter%dphi * cos_phi_tok
    ! dB_y/dz
    dB_y_dR(3) = dB_r_inter%dz * sin_phi_tok + dB_t_inter%dz * cos_phi_tok
    ! dB_z/dR
    dB_z_dR(1) = dB_z_inter%dR
    ! dB_z/dphi
    dB_z_dR(2) = dB_z_inter%dphi
    ! dB_z/dz
    dB_z_dR(3) = dB_z_inter%dz
    ! Now with the chain rule compute the entries of dB_i/dx_j
    ! dBx/dx third term is zero since dz/dx = 0
    spatial_grad_B(1,1) = dB_x_dR(1) * func_dR_dx(x_vec(1), x_vec(2)) + dB_x_dR(2) * func_dphi_dx(x_vec(1), x_vec(2))
    ! dBx/dy third term is zero since dz/dx = 0
    spatial_grad_B(1,2) = dB_x_dR(1) * func_dR_dy(x_vec(1), x_vec(2)) + dB_x_dR(2) * func_dphi_dy(x_vec(1), x_vec(2))
    ! dBx/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(1,3) = dB_x_dR(3)
    ! dBy/dx third term is zero since dz/dy = 0
    spatial_grad_B(2,1) = dB_y_dR(1) * func_dR_dx(x_vec(1), x_vec(2)) + dB_y_dR(2) * func_dphi_dx(x_vec(1), x_vec(2))
    ! dBy/dy third term is zero since dz/dy = 0
    spatial_grad_B(2,2) = dB_y_dR(1) * func_dR_dy(x_vec(1), x_vec(2)) + dB_y_dR(2) * func_dphi_dy(x_vec(1), x_vec(2))
    ! dBy/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(2,3) = dB_y_dR(3)
    ! dBz/dx third term is zero since dz/dx = 0
    spatial_grad_B(3,1) = dB_z_dR(1) * func_dR_dx(x_vec(1), x_vec(2))  + dB_z_dR(2) * func_dphi_dx(x_vec(1), x_vec(2))
    ! dBz/dy third term is zero since dz/dy = 0
    spatial_grad_B(3,2) = dB_z_dR(1) * func_dR_dy(x_vec(1), x_vec(2))  + dB_z_dR(2) * func_dphi_dy(x_vec(1), x_vec(2))
    ! dBz/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(3,3) = dB_z_dR(3)
    !stop "correct?"
    !---gradient of total magnetic field
    B_abs = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
    ! to obtain d|B|/dx it is split into d|B|/dB_x_j and dB_x_i/dx_j
    ! -> d|B|/dx_i = B_x_j / |B| * dB_j/dx_i
    spatial_grad_B_abs(:) = 0.d0
    do i = 1, 3
        do j = 1,3
          spatial_grad_B_abs(i) = spatial_grad_B_abs(i) + B_x_vec(j) * spatial_grad_B(j,i)
        end do
    end do
    spatial_grad_B_abs = spatial_grad_B_abs / B_abs
    !B_vec_norm = B_x_vec / B_abs
    !N_vec_norm = N_vec / N_abs
    scal_prod = 0.d0
    do i =1, 3
      scal_prod = scal_prod + B_x_vec(i) * N_vec(i)
    end do
    N_par = scal_prod / B_abs
    !------------------------------ Now that all mangetic field gradients are known compute dNpar/dN_i------!
    !
    N_grad_N_par =  B_x_vec(:) / B_abs
    !----------------------------- Finally the spatial gardiend of N_par -------------------------!
    !----------------------------- Here we resolve dN_par/dx into dN_par/dB_j * dB_j/dx_i  (multidimensional chain rule)!
    spatial_grad_N_par(:) = 0.d0
    do i = 1,3
      do j =1,3
        do k = 1,3
        spatial_grad_N_par(i)  = spatial_grad_N_par(i) + (delta(j,k) * (N_vec(j) * B_abs**2 - B_x_vec(j)**2 * N_vec(j)) - &
                                 (1.d0 - delta(j,k)) * (N_vec(k) * B_x_vec(j) * B_x_vec(k))) * spatial_grad_B(j,i)
        end do
      end do
    end do
    spatial_grad_N_par(:) = spatial_grad_N_par / (B_abs)**3
    if(debug_level >= 3 .or. (debug_level ==  2 .and.  sqrt(x_vec(1)**2 + x_vec(2)**2) <= 2.13)) then
      print*, "/---------------DEBUG OUTPUT------------------\"
      print*,"--------Begin of magnetic field debug output-------"
      print*,"x = ", x_vec(1), "m "
      print*,"y = ", x_vec(2), "m "
      print*,"z = ", x_vec(3), "m "
      print*,"R = ", R_vec(1), " m"
      print*,"z = ", R_vec(3), " m"
      print*, "B_R", B_R_vec
      print*, "N_vec", N_vec
      print*, "B_abs", B_abs
      print*, "N_par", N_par
      print*, "dR/dx", func_dR_dx(x_vec(1), x_vec(2))
      print*, "dphi/dx", func_dphi_dx(x_vec(1), x_vec(2))
      print*, "dR/dy", func_dR_dy(x_vec(1), x_vec(2))
      print*, "dphi/dy", func_dphi_dy(x_vec(1), x_vec(2))
!      print*, "Closest spline knot"
!      print*,  plasma_params%R(minloc(abs(plasma_params%R - R_vec(1))))
!      print*,  plasma_params%z(minloc(abs(plasma_params%z - R_vec(3))))
      do i = 1, 3
         if(i ==1) print*,"--------dR-------"
         if(i ==2) print*,"--------dphi-------"
         if(i ==3) print*,"--------dz-------"
         do j = 1, 4
            aux_x(j,:) = x_vec
            if(j < 3) aux_x(j,i) = x_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
            if(j >= 3) aux_x(j,i) = x_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
            aux_R(j,:) = R_vec
            if(j < 3) aux_R(j,i) = R_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
            if(j >= 3) aux_R(j,i) = R_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
         end do
         print*,"x ana",spatial_grad_B(1,i)
         print*, "x num",(- func_B_x(plasma_params, aux_x(1,:)) + &
                                    8.d0 *  func_B_x(plasma_params, aux_x(2,:)) - &
                                8.d0 *  func_B_x(plasma_params, aux_x(3,:))  &
                                +   func_B_x(plasma_params, aux_x(4,:))) / (12.d0 *h_x)
         print*,"y ana",spatial_grad_B(2,i)
         print*, "y num",(- func_B_y(plasma_params, aux_x(1,:)) + &
                                    8.d0 *  func_B_y(plasma_params, aux_x(2,:)) - &
                                8.d0 *  func_B_y(plasma_params, aux_x(3,:))  &
                                +   func_B_y(plasma_params, aux_x(4,:))) / (12.d0 *h_x)
        print*,"z ana",spatial_grad_B(3,i)
        print*, "z num",(- func_B_z(plasma_params, aux_x(1,:)) + &
                                    8.d0 *  func_B_z(plasma_params, aux_x(2,:)) - &
                                8.d0 *  func_B_z(plasma_params, aux_x(3,:))  &
                                +   func_B_z(plasma_params, aux_x(4,:))) / (12.d0 *h_x)
        print*,"dBx/dR ana",dB_x_dR(i)
        print*,"dBx/dR num",(- func_B_x_R(plasma_params, aux_R(1,:)) + &
                                    8.d0 *  func_B_x_R(plasma_params, aux_R(2,:)) - &
                                8.d0 *  func_B_x_R(plasma_params, aux_R(3,:))  &
                                +   func_B_x_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
        print*,"dBy/dR ana",dB_y_dR(i)
        print*,"dBy/dR num",(- func_B_y_R(plasma_params, aux_R(1,:)) + &
                                    8.d0 *  func_B_y_R(plasma_params, aux_R(2,:)) - &
                                8.d0 *  func_B_y_R(plasma_params, aux_R(3,:))  &
                                +   func_B_y_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
        print*,"dBz/dR ana ", dB_z_dR(i)
        print*,"dBz/dR num",(- func_B_z_R(plasma_params, aux_R(1,:)) + &
                                    8.d0 *  func_B_z_R(plasma_params, aux_R(2,:)) - &
                                8.d0 *  func_B_z_R(plasma_params, aux_R(3,:))  &
                                +   func_B_z_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
       if(i == 1) then
          print*,"dBR/dR ana",dB_r_inter%dR
          print*,"dBR/dR num",(- func_B_r_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_r_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_r_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_r_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBt/dR ana",dB_t_inter%dR
          print*,"dBt/dR num",(- func_B_phi_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_phi_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_phi_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_phi_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBz/dR ana", dB_z_inter%dz
          print*,"dBz/dR num", (- func_B_z_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_z_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_z_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_z_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
      else if(i == 2) then
          print*,"dBR/dphi ana",dB_r_inter%dphi
          print*,"dBR/dphi num",(- func_B_r_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_r_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_r_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_r_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBt/dphi ana",dB_t_inter%dphi
          print*,"dBt/dphi num",(- func_B_phi_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_phi_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_phi_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_phi_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBz/dz ana", dB_z_inter%dz
          print*,"dBz/dz num",(- func_B_z_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_z_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_z_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_z_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
      else
          print*,"dBR/dz ana",dB_r_inter%dz
          print*,"dBR/dz num",(- func_B_r_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_r_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_r_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_r_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBt/dz ana",dB_t_inter%dz
          print*,"dBt/dz num",(- func_B_phi_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_phi_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_phi_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_phi_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
          print*,"dBz/dz ana",dB_z_inter%dz
          print*,"dBz/dz num",(- func_B_z_R(plasma_params, aux_R(1,:)) + &
                                      8.d0 *  func_B_z_R(plasma_params, aux_R(2,:)) - &
                                  8.d0 *  func_B_z_R(plasma_params, aux_R(3,:))  &
                                  +   func_B_z_R(plasma_params, aux_R(4,:))) / (12.d0 *h_x)
      end if
      end do
    end if
  end subroutine sub_grad_N_par



  function delta(i,j) !Kronecker-Delta
    USE f90_kind
    implicit none
    integer(ikind), intent(in) :: i,j
    real(rkind)                :: delta
    delta = 0.5d0 * real(sign(1,i-j) + sign(1,j-i),8)
  end function delta




  subroutine sub_spatial_grad_X(plasma_params, omega, x_vec, X,  spatial_grad_X, rhop_out)
  ! dX/ dx
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma
    USE ripple3d,                 only: grad_type
    use constants,                 only : pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_utils, only: retrieve_n_e, retrieve_T_e, retrieve_n_e_mat_single, retrieve_T_e_mat_single
    ! corresponds to flux coordinates
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind)              ,   intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind)              ,   intent(out)  :: X
    real(rkind), dimension(:),   intent(out)  :: spatial_grad_X
    real(rkind)              ,   intent(out)  :: rhop_out
    real(rkind)                               :: rhop, n_e, grad_n_e, T_e, grad_T_e
    real(rkind), dimension(3)                 :: spatial_grad_rhop_ne, spatial_grad_rhop_Te
    if(.not. plasma_params%Te_ne_mat) then
      call sub_spatial_grad_rhop(plasma_params, x_vec, rhop, spatial_grad_rhop_ne)
      if(rhop == -1.d0) then
        spatial_grad_X = 0.d0
        X= 0.d0
        rhop_out = -1.d0
        return
      end if
      call retrieve_n_e(plasma_params, rhop, n_e, grad_n_e)
      if(warm_plasma) then
        call retrieve_T_e(plasma_params, rhop, T_e, grad_T_e)
        spatial_grad_rhop_Te = spatial_grad_rhop_ne * plasma_params%rhop_scale_Te * grad_T_e
      end if
      spatial_grad_rhop_ne =  spatial_grad_rhop_ne * plasma_params%rhop_scale_ne * grad_n_e
    else
      call retrieve_n_e_mat_single(plasma_params, x_vec, n_e, spatial_grad_rhop_ne) ! grad ne
      if(warm_plasma) then
        call retrieve_T_e_mat_single(plasma_params, x_vec, T_e, spatial_grad_rhop_Te) ! grad Te
      end if
      rhop = func_rhop(plasma_params, x_vec)
    end if
    if(warm_plasma) then
      spatial_grad_X = (e0**2*(2.d0*(c0**2 * mass_e + 5.d0 * e0 * T_e) * &
                       spatial_grad_rhop_ne - &
                       5.d0 * e0 * n_e * spatial_grad_rhop_Te)) / &
                       (2.d0*eps0*mass_e * omega**2 * (c0**2 * mass_e + 5.d0 * e0 * T_e) * &
                       Sqrt(1.d0 + (5.d0 * e0 * T_e)/(c0**2 * mass_e)))
    else
      spatial_grad_X = spatial_grad_rhop_ne(:) * e0**2.d0/(eps0 * mass_e * omega**2)
    end if
    X  =  func_X(plasma_params, omega, n_e, T_e)
    rhop_out = rhop
  end subroutine sub_spatial_grad_X



  function func_X(plasma_params, omega, n_e, T_e)
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma
    use constants,            only: pi, e0, mass_e, eps0, c0
    implicit None
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind),                 intent(in)   :: omega, n_e, T_e
    real(rkind)                               :: func_X
    if(warm_plasma) then
      func_X = n_e * e0**2.d0/(eps0 * mass_e * Sqrt(1.d0 + (5.0 * e0 * T_e)/(c0**2 * mass_e)) * omega**2)
    else
      func_X = n_e * e0**2.d0/(eps0 * mass_e * omega**2)
    end if
  end function func_X

 subroutine sub_spatial_grad_N_par(plasma_params, x_vec, N_vec, N_abs, N_par, spatial_grad_N_par, spatial_grad_B_abs, B_abs, N_grad_N_par)
  ! dN_par / dx
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type
    implicit none
    type(plasma_params_type), intent(in)        :: plasma_params
    real(rkind), dimension(:),      intent(in)  :: x_vec, N_vec
    real(rkind),                    intent(in)  :: N_abs
    real(rkind),                    intent(out) :: N_par, B_abs
    real(rkind), dimension(:),      intent(out) :: spatial_grad_N_par
    real(rkind), dimension(:),      intent(out) :: spatial_grad_B_abs
    real(rkind), dimension(:),      intent(out) :: N_grad_N_par
    call sub_grad_N_par(plasma_params, x_vec, N_vec, N_abs, B_abs, spatial_grad_B_abs, &
                           N_par, N_grad_N_par, spatial_grad_N_par)

  end subroutine sub_spatial_grad_N_par

  function func_Y(plasma_params, omega, B_abs, T_e)
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma
    use constants,            only: pi, e0, mass_e, eps0, c0
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind),                 intent(in)   :: omega, B_abs, T_e
    real(rkind)                               :: func_Y
    if(warm_plasma) then
      func_Y = B_abs * e0 / (mass_e * Sqrt(1.d0 + (5.0 * e0 * T_e)/(c0**2 * mass_e)) * omega)
    else
      func_Y = B_abs * e0 / (mass_e * omega)
    end if
  end function func_Y

  subroutine sub_spatial_grad_Y2(plasma_params, omega, x_vec, B_abs, spatial_grad_B_abs, Y, spatial_grad_Y2)
  ! dY^2 / dx
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type
    ! corresponds to flux coordinates
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind),                 intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind),                 intent(in)   :: B_abs
    real(rkind), dimension(:),   intent(in)   :: spatial_grad_B_abs
    real(rkind),                 intent(out)  :: Y
    real(rkind), dimension(:),   intent(out)  :: spatial_grad_Y2
      call sub_spatial_grad_Y(plasma_params, omega, x_vec, B_abs, spatial_grad_B_abs, Y, spatial_grad_Y2)
      spatial_grad_Y2 = 2.d0 * Y * spatial_grad_Y2
  end subroutine sub_spatial_grad_Y2

  subroutine sub_spatial_grad_Y(plasma_params, omega, x_vec, B_abs, spatial_grad_B_abs, Y, spatial_grad_Y)
  ! dY^2 / dx
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma
    use constants,            only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_utils, only: retrieve_T_e, retrieve_T_e_mat_single
    ! corresponds to flux coordinates
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind),                 intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind),                 intent(in)   :: B_abs
    real(rkind), dimension(:),   intent(in)   :: spatial_grad_B_abs
    real(rkind),                 intent(out)  :: Y
    real(rkind), dimension(:),   intent(out)  :: spatial_grad_Y
    real(rkind), dimension(3)                 :: spatial_grad_rhop ! grad_x(B^2)
    real(rkind)                               :: rhop, T_e, grad_T_e
    if(warm_plasma) then
      if(.not. plasma_params%Te_ne_mat) then
        call sub_spatial_grad_rhop(plasma_params, x_vec, rhop, spatial_grad_rhop)
        if(rhop == -1.d0) then
          spatial_grad_Y = 0.d0
          Y = 0.d0
          return
        end if
        spatial_grad_rhop = spatial_grad_rhop * plasma_params%rhop_scale_Te
        call retrieve_T_e(plasma_params, rhop, T_e, grad_T_e)
        spatial_grad_rhop = grad_T_e * spatial_grad_rhop
      else
        call retrieve_T_e_mat_single(plasma_params, x_vec, T_e, spatial_grad_rhop)
      end if
      spatial_grad_Y(:) = (e0*((0.2*c0**2*mass_e + e0*T_e) * spatial_grad_B_abs(:) - &
                          0.5*e0*B_abs*spatial_grad_rhop)) / &
                          (mass_e*omega*(0.2*c0**2*mass_e + e0*T_e) * &
                          Sqrt(1. + (5.*e0*T_e)/(c0**2*mass_e)))
    else
      spatial_grad_Y(:) = spatial_grad_B_abs(:) * e0 / (mass_e * omega)
      T_e = 0.d0
    end if
    Y = func_Y(plasma_params, omega, B_abs, T_e)
  end subroutine sub_spatial_grad_Y



 subroutine sub_grad_H(plasma_params, omega, mode, x_vec, dx_dsigma)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, h_x_glob
    USE ripple3d,                 only: grad_type
    USE constants,                 only : eps0, mass_e, e0, c0
    use mod_ecfm_refr_utils, only: retrieve_n_e, retrieve_T_e, sub_remap_coords, &
                                   retrieve_n_e_mat_single, &
                                   retrieve_T_e_mat_single
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params ! only slice of the entire ray
    real(rkind),                 intent(in)   :: omega
    integer(ikind)                            :: mode
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind), dimension(:),   intent(out)  :: dx_dsigma
    real(rkind)                               :: N_abs, N_par, N_perp, X, Y, rhop_out, B_abs
    real(rkind), dimension(3)                 :: R_vec, spatial_grad_X, spatial_grad_Y, &
                                                 spatial_grad_N_par, spatial_grad_B_abs, N_grad_N_par
    real(rkind)                               :: dH_dA, dH_dB, dH_dC, H, A, B ,C, dB_dN_par, &
                                                 dC_dN_par, h_x, B_abs_aux, rhop, T_e, n_e
    real(rkind), dimension(4)                 :: X_aux, Y_aux, N_perp_aux, N_perp_aux_2, N_par_aux,N_abs_aux, &
                                                 N_par_aux_2, A_aux,B_aux, C_aux, B_aux_2, C_aux_2
    real(rkind), dimension(4,3)               :: aux_x, aux_N
    integer(ikind)                            :: i, j
    N_abs = sqrt(x_vec(4)**2 + x_vec(5)**2 + x_vec(6)**2)
    call sub_spatial_grad_X(plasma_params, omega, x_vec(1:3), X, spatial_grad_X, rhop_out)
    call sub_spatial_grad_N_par(plasma_params, x_vec(1:3), x_vec(4:6), N_abs, N_par, spatial_grad_N_par, spatial_grad_B_abs, B_abs, N_grad_N_par)
    call sub_spatial_grad_Y(plasma_params, omega, x_vec(1:3), B_abs, spatial_grad_B_abs, Y, spatial_grad_Y)
    A = func_A(X, Y)
    B = func_B(X, Y, N_par)
    C = func_C(X, Y, N_par)
    dH_dA = func_dH_dA(A, B, C, mode)
    dH_dB = func_dH_dB(A, B, C, mode)
    dH_dC = func_dH_dC(A, B, C, mode)
    dB_dN_par = func_dB_dN_par(X, Y, N_par)
    dC_dN_par = func_dC_dN_par(X, Y, N_par)
    ! correction factor of (1 - Y^2) applied to hamiltonian
    ! first term is from d/dx (N^2 ( 1- Y^2)) = -N^2 dY^2/dx
    dx_dsigma(4:6) =  spatial_grad_N_par(:) * (-4.d0 * N_par +  dB_dN_par * dH_dB + &
                          dC_dN_par * dH_dC)
    dx_dsigma(4:6) = dx_dsigma(4:6) + spatial_grad_X(:) * &
                       (func_dA_dX(X, Y) * dH_dA + &
                        func_dB_dX(X, Y, N_par) * dH_dB + &
                       func_dC_dX(X, Y, N_par) * dH_dC)
    dx_dsigma(4:6) = dx_dsigma(4:6) + spatial_grad_Y(:) * &
                       (func_dA_dY(X, Y) * dH_dA + &
                        func_dB_dY(X, Y, N_par) * dH_dB + &
                        func_dC_dY(X, Y, N_par) * dH_dC)
    dx_dsigma(4:6) = - dx_dsigma(4:6)
    dx_dsigma(1:3) = 4.d0 * x_vec(4:6) + N_grad_N_par(:) * (-4.d0 * N_par +  dB_dN_par * dH_dB + &
                     dC_dN_par * dH_dC)
    !! Debug !!

    N_perp = sqrt(N_abs**2 - N_par**2)
    H = func_H(N_perp, A, B, C, mode)
!    if(abs(H / plasma_params%H_last) > plasma_params%trigger_level .and. H > 1.e-6) then
!      debug_level = 3
!      PRINT*, plasma_params%H_last
!    else
!     plasma_params%H_last =  H
!    end if
    if((debug_level >= 3 .or. (debug_level ==  2 .and.  sqrt(x_vec(1)**2 + x_vec(2)**2) <= 2.13)))then! .or. abs(H) > 1.d-2) then
      print*, "/---------------DEBUG OUTPUT------------------\"
    else
      return
    end if
    rhop = func_rhop(plasma_params, x_vec(1:3))
    print*, "H", H
    print*, "omega", omega
    print*, "rhop", rhop
    print*, "x", x_vec(1:3)
    print*, "N", x_vec(4:6)
    print*, "ABC", A, B, C
    print*, "X,Y, N_par", X, Y, N_par
    print*, "N_abs", N_abs
    call sub_N_par(plasma_params, x_vec(1:3), x_vec(4:6), N_par_aux(1), N_abs_aux(1), B_abs_aux)
    if(.not. plasma_params%Te_ne_mat) then
      rhop = func_rhop(plasma_params, x_vec(1:3))
      call retrieve_n_e(plasma_params, rhop, n_e)
      call retrieve_T_e(plasma_params, rhop, T_e)
    else
      call retrieve_n_e_mat_single(plasma_params, x_vec, n_e)
      call retrieve_T_e_mat_single(plasma_params, x_vec, T_e)
    end if
    X_aux(1) = func_X(plasma_params, omega, n_e, T_e)
    Y_aux(1) =  func_Y(plasma_params,omega, B_abs_aux,T_e)
    print*, "N_par defaul", N_par, "N_par debug", N_par_aux(1)
    print*, "X defaul", X, "X debug", X_aux(1)
    print*, "Y defaul", Y, "Y debug", Y_aux(1)
    do i =1,3
      h_x = h_x_glob! * abs(dx_dsigma(i + 3))
      do j = 1 , 4
        aux_x(j,:) = x_vec(1:3)
        aux_N(j,:) = x_vec(4:6)
        if(j < 3) aux_x(j,i) = x_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
        if(j >= 3) aux_x(j,i) = x_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
        if(j < 3) aux_N(j,i) = x_vec(i + 3) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
        if(j >= 3) aux_N(j,i) = x_vec(i + 3) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
        call sub_N_par(plasma_params, aux_x(j,:), x_vec(4:6), N_par_aux(j), N_abs_aux(j),  B_abs_aux)
        if(.not. plasma_params%Te_ne_mat) then
          rhop = func_rhop(plasma_params, x_vec(1:3))
          call retrieve_n_e(plasma_params, rhop, n_e)
          call retrieve_T_e(plasma_params, rhop, T_e)
        else
          call retrieve_n_e_mat_single(plasma_params, x_vec(1:3), n_e)
          call retrieve_T_e_mat_single(plasma_params, x_vec(1:3), T_e)
        end if
        X_aux(j) = func_X(plasma_params, omega, n_e, T_e)
        Y_aux(j) = func_Y(plasma_params, omega, B_abs_aux, T_e)
        N_perp_aux(j) = sqrt(N_abs_aux(j)**2 - N_par_aux(j)**2)
        call sub_N_par(plasma_params, x_vec(1:3), aux_N(j,:), N_par_aux_2(j), N_abs_aux(j),  B_abs_aux)
        N_perp_aux_2(j) = sqrt(N_abs_aux(j)**2 - N_par_aux_2(j)**2)
        A_aux(j) = func_A(X_aux(j), Y_aux(j))
        B_aux(j) = func_B(X_aux(j), Y_aux(j), N_par_aux(j))
        C_aux(j) = func_C(X_aux(j), Y_aux(j), N_par_aux(j))
        B_aux_2(j) = func_B(X, Y, N_par_aux_2(j))
        C_aux_2(j) = func_C(X, Y, N_par_aux_2(j))
      end do
      call sub_remap_coords(x_vec(1:3), R_vec)
      print*, "Current position", R_vec
      print*, "Current rhop", rhop_out
      print*, "X_aux", X_aux
      print*, "Y_aux", Y_aux
      print*, "N_par_aux", N_par_aux
      print*, "dX/dx ana", spatial_grad_X(i)
      print*, "dX/dx num fwd. :", (X_aux(2) - X)/h_x
      print*, "dX/dx num ctr. :" , (- X_aux(1) + 8.d0 *  X_aux(2) - &
                          8.d0 *  X_aux(3) +   X_aux(4)) / (12.d0 *h_x)
      print*, "dY/dx ana",spatial_grad_Y(i)
      print*, "dY/dx num fwd.", (Y_aux(2) - Y)/h_x
      print*, "dY/dx num ctr. :" ,(- Y_aux(1) + 8.d0 *  Y_aux(2) - &
                          8.d0 *  Y_aux(3) +   Y_aux(4)) / (12.d0 *h_x)
      print*, "dN_par/dx ana",spatial_grad_N_par(i)
      print*, "dN_par/dx num fwd.", (N_par_aux(2) - N_par)/h_x
      print*, "dN_par/dx num ctr. :" , (- N_par_aux(1) + 8.d0 *  N_par_aux(2) - &
                          8.d0 *  N_par_aux(3) +   N_par_aux(4)) / (12.d0 *h_x)
      print*, "dN_par/dN ana",N_grad_N_par(i)
      print*, "dN_par/dN num fwd.", (N_par_aux_2(2) - N_par)/h_x
      print*, "dN_par/dN num ctr. :" , (- N_par_aux_2(1) + 8.d0 *  N_par_aux_2(2) - &
                          8.d0 *  N_par_aux_2(3) +   N_par_aux_2(4)) / (12.d0 *h_x)

      print*, "dA/dX ana", func_dA_dX(X, Y)
      print*, "dA/dX num fwd", ( func_A( X + h_x, Y) - A) / h_x
      print*, "dA/dX num ctr", (- func_A( X +  2.d0 * h_x, Y) + &
                              8.d0 *  func_A( X +  h_x, Y) - &
                          8.d0 *  func_A( X - h_x, Y)  &
                          +   func_A(  X - 2.d0 * h_x, Y)) / (12.d0 *h_x)
      print*, "dA/dY ana", func_dA_dY(X, Y)
      print*, "dA/dY num fwd", ( func_A( X , Y + h_x) - A) / h_x
      print*, "dA/dY num ctr", (- func_A( X , Y +  2.d0 * h_x) + &
                              8.d0 *  func_A( X, Y +  h_x) - &
                          8.d0 *  func_A( X, Y -  h_x)  &
                          +   func_A(  X , Y- 2.d0 * h_x)) / (12.d0 *h_x)
      print*, "dB/dX ana", func_dB_dX(X, Y, N_par)
      print*, "dB/dX num fwd", ( func_B( X + h_x, Y, N_par) - B) / h_x
      print*, "dB/dX num ctr", (- func_B( X +  2.d0 * h_x, Y, N_par) + &
                              8.d0 *  func_B( X +  h_x, Y, N_par) - &
                          8.d0 *  func_B( X - h_x, Y, N_par)  &
                          +   func_B(  X - 2.d0 * h_x, Y, N_par)) / (12.d0 *h_x)
      print*, "dB/dY ana", func_dB_dY(X, Y, N_par)
      print*, "dB/dY num fwd", ( func_B( X , Y + h_x, N_par) - B) / h_x
      print*, "dB/dY num ctr", (- func_B( X , Y +  2.d0 * h_x, N_par) + &
                              8.d0 *  func_B( X, Y +  h_x, N_par) - &
                          8.d0 *  func_B( X, Y -  h_x, N_par)  &
                          +   func_B(  X , Y- 2.d0 * h_x, N_par)) / (12.d0 *h_x)
      print*, "dB/dN_par ana", dB_dN_par
      print*, "dB/dN_par num fwd", ( func_B( X , Y, N_par + h_x) - B) / h_x
      print*, "dB/dN_par num ctr", (- func_B( X , Y, N_par +  2.d0 * h_x) + &
                              8.d0 *  func_B( X, Y, N_par +  h_x) - &
                          8.d0 *  func_B( X, Y, N_par -  h_x)  &
                          +   func_B(  X , Y , N_par - 2.d0 * h_x)) / (12.d0 *h_x)
      print*, "dC/dX ana", func_dC_dX(X, Y, N_par)
      print*, "dC/dX num fwd", ( func_C( X + h_x, Y, N_par) - C) / h_x
      print*, "dC/dX num ctr", (- func_C( X +  2.d0 * h_x, Y, N_par) + &
                              8.d0 *  func_C( X +  h_x, Y, N_par) - &
                          8.d0 *  func_C( X - h_x, Y, N_par)  &
                          +   func_C(  X - 2.d0 * h_x, Y, N_par)) / (12.d0 *h_x)
      print*, "dC/dY ana", func_dC_dY(X, Y, N_par)
      print*, "dC/dY num fwd", ( func_C( X , Y + h_x, N_par) - C) / h_x
      print*, "dC/dY num ctr", (- func_C( X , Y +  2.d0 * h_x, N_par) + &
                              8.d0 *  func_C( X, Y +  h_x, N_par) - &
                          8.d0 *  func_C( X, Y -  h_x, N_par)  &
                          +   func_C(  X , Y- 2.d0 * h_x, N_par)) / (12.d0 *h_x)
      print*, "dC/dN_par ana", dC_dN_par
      print*, "dC/dN_par num fwd", ( func_C( X , Y, N_par + h_x) - C) / h_x
      print*, "dC/dN_par num ctr", (- func_C( X , Y, N_par +  2.d0 * h_x) + &
                              8.d0 *  func_C( X, Y, N_par +  h_x) - &
                          8.d0 *  func_C( X, Y, N_par -  h_x)  &
                          +   func_C(  X , Y , N_par - 2.d0 * h_x)) / (12.d0 *h_x)
      print*, "dH/dA ana", dH_dA
      print*, "dH/dA num fwd", ( func_H(N_perp, A + h_x, B, C, mode) - H) / h_x
      print*, "dH/dA num ctr", (- func_H(N_perp,  A + 2.d0 * h_x, B, C, mode) + &
                              8.d0 *  func_H(N_perp,  A + h_x, B, C, mode) - &
                          8.d0 *  func_H(N_perp,  A - h_x, B, C, mode)  &
                          +   func_H(N_perp,  A - 2.d0 * h_x, B, C, mode)) / (12.d0 *h_x)
      print*, "dH/dB ana", dH_dB
      print*, "dH/dB num fwd", ( func_H(N_perp,  A, B + h_x, C, mode) - H) / h_x
      print*, "dH/dB num ctr", (- func_H(N_perp,  A, B + 2.d0 * h_x, C, mode) + &
                              8.d0 *  func_H(N_perp,  A, B + h_x, C, mode) - &
                          8.d0 *  func_H(N_perp,  A, B - h_x, C, mode)  &
                          +   func_H(N_perp,  A,  B - 2.d0 * h_x, C, mode)) / (12.d0 *h_x)
      print*, "dH/dC ana", dH_dC
      print*, "dH/dC num fwd", ( func_H(N_perp,  A, B, C + h_x, mode) - H) / h_x
      print*, "dH/dC num ctr", (- func_H(N_perp,  A, B, C + 2.d0 * h_x, mode) + &
                              8.d0 *  func_H(N_perp,  A, B, C + h_x, mode) - &
                          8.d0 *  func_H(N_perp,  A, B, C - h_x, mode)  &
                          +   func_H(N_perp, A, B, C - 2.d0 * h_x, mode)) / (12.d0 *h_x)
      print*, "d(2* N_perp^2)/dN ana", 4.d0 * x_vec(i + 3) - N_grad_N_par(i) * 4.d0 * N_par
      print*, "d(2* N_perp^2)/dN num fwd", ( func_N_perp_pt(N_abs_aux(2), N_par_aux_2(2)) - N_perp**2 * 2.d0) / h_x
      print*, "d(2* N_perp^2)/dN num ctr", (- func_N_perp_pt(N_abs_aux(1), N_par_aux_2(1)) + &
                              8.d0 *  func_N_perp_pt(N_abs_aux(2), N_par_aux_2(2)) - &
                          8.d0 *  func_N_perp_pt(N_abs_aux(3), N_par_aux_2(3)) &
                          +   func_N_perp_pt(N_abs_aux(4), N_par_aux_2(4))) / (12.d0 *h_x)
      print*,"ana", dx_dsigma(i + 3)
      print*,"num fwd",(func_H(N_perp_aux(2), A_aux(2),  B_aux(2), C_aux(2), mode) - H) / h_x
      print*, "num backward",  -(func_H( N_perp_aux(3), A_aux(3), B_aux(3), C_aux(3), mode) - H) / h_x
      print*, "num ctr:", (- func_H(N_perp_aux(1), A_aux(1),B_aux(1), C_aux(1), mode) + &
                              8.d0 *  func_H(N_perp_aux(2), A_aux(2), B_aux(2), C_aux(2), mode) - &
                          8.d0 *  func_H(N_perp_aux(3), A_aux(3), B_aux(3), C_aux(3), mode)  &
                          +   func_H(N_perp_aux(4), A_aux(4), B_aux(4), C_aux(4), mode)) / (-12.d0 *h_x)
      print*,"ana", dx_dsigma(i)
      print*,"num fwd",(func_H( N_perp_aux_2(2), A, B_aux_2(2), C_aux_2(2), mode) - H) / h_x
      print*, "num backward",  -(func_H(N_perp_aux_2(3), A, B_aux_2(3), C_aux_2(3), mode) - H) / h_x
      print*, "num ctr:", (- func_H(N_perp_aux_2(1), A,  B_aux_2(1), C_aux_2(1), mode) + &
                              8.d0 *  func_H(N_perp_aux_2(2), A,  B_aux_2(2), C_aux_2(2), mode) - &
                          8.d0 *  func_H(N_perp_aux_2(3), A,  B_aux_2(3), C_aux_2(3), mode)  &
                          +   func_H(N_perp_aux_2(4),  A, B_aux_2(4), C_aux_2(4), mode)) / (12.d0 *h_x)
!        dx_dsigma(i) = (- func_H(N_perp_aux_2(1), A,  B_aux_2(1), C_aux_2(1), mode) + &
!                                8.d0 *  func_H(N_perp_aux_2(2), A,  B_aux_2(2), C_aux_2(2), mode) - &
!                            8.d0 *  func_H(N_perp_aux_2(3), A,  B_aux_2(3), C_aux_2(3), mode)  &
!                            +   func_H(N_perp_aux_2(4),  A, B_aux_2(4), C_aux_2(4), mode)) / (12.d0 *h_x)
!        dx_dsigma(i + 3) = (- func_H(N_perp_aux(1), A_aux(1),B_aux(1), C_aux(1), mode) + &
!                                8.d0 *  func_H(N_perp_aux(2), A_aux(2), B_aux(2), C_aux(2), mode) - &
!                           8.d0 *  func_H(N_perp_aux(3), A_aux(3), B_aux(3), C_aux(3), mode)  &
!                            +   func_H(N_perp_aux(4), A_aux(4), B_aux(4), C_aux(4), mode)) / (-12.d0 *h_x)
    end do
    print*, " Critical error in raytracing"
    call abort
  end subroutine sub_grad_H

  subroutine sub_grad_Lambda(plasma_params, omega, mode, x_vec, dxds)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, h_x_glob
    USE ripple3d,                 only: grad_type
    USE constants,                 only : c0, eps0, mass_e, e0
    use mod_ecfm_refr_utils, only: retrieve_n_e, retrieve_T_e, sub_remap_coords, &
                                   retrieve_n_e_mat_single, &
                                   retrieve_T_e_mat_single
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params ! only slice of the entire ray
    real(rkind),                 intent(in)   :: omega
    integer(ikind),              intent(in)   :: mode
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind), dimension(:),   intent(out)  :: dxds
    real(rkind)                               :: N_abs, N_par, X, Y, rhop_out
    real(rkind), dimension(3)                 :: R_vec, N_vec, spatial_grad_X, spatial_grad_Y, spatial_grad_Y2, &
                                                 spatial_grad_N_par, spatial_grad_B_abs, N_grad_N_par
    real(rkind)                               :: Hamil, N_s, B_abs, n_e, T_e, rhop, N_abs_aux, h_x, dNs_sq_dN_par, B_abs_aux
    real(rkind), dimension(4)                 :: X_aux, Y_aux, N_par_aux, N_abs_aux_2
    real(rkind), dimension(4,3)               :: aux_x, aux_N
    integer(ikind)                            :: i, j
    N_abs = sqrt(x_vec(4)**2 + x_vec(5)**2 + x_vec(6)**2)
    N_vec = x_vec(4:6)
    call sub_spatial_grad_X(plasma_params, omega, x_vec(1:3), X, spatial_grad_X, rhop_out)
    call sub_spatial_grad_N_par(plasma_params, x_vec(1:3), x_vec(4:6), N_abs, N_par, spatial_grad_N_par, spatial_grad_B_abs, B_abs, N_grad_N_par)
    call sub_spatial_grad_Y(plasma_params, omega, x_vec(1:3), B_abs, spatial_grad_B_abs, Y, spatial_grad_Y)
    spatial_grad_Y2 = 2.d0 * Y * spatial_grad_Y
    dNs_sq_dN_par = func_dNs_sq_dN_par(N_abs,  X, Y, N_par, mode)
    dxds(4:6) = -spatial_grad_N_par(:) * dNs_sq_dN_par
    dxds(4:6) = dxds(4:6) - spatial_grad_X(:) * func_dNs_sq_dX(N_abs,  X, Y, N_par, mode)
    dxds(4:6) = dxds(4:6) - spatial_grad_Y2(:) * func_dNs_sq_dY2(N_abs,  X, Y, N_par, mode)
    dxds(1:3) = 2.d0  * N_vec(:) - N_grad_N_par(:) * dNs_sq_dN_par
    if((debug_level >= 3 .or. (debug_level ==  2 .and.  sqrt(x_vec(1)**2 + x_vec(2)**2) <= 2.13)))then! .or. abs(H) > 1.d-2) then
      print*, "/---------------DEBUG OUTPUT------------------\"
    else
      return
    end if
    Hamil = func_Lambda(N_abs, X, Y,  N_par, mode)
    print*, "Hamil", Hamil
    N_s = sqrt(func_N_s_2(X, Y, N_par, mode))
    do i =1,3
      h_x = h_x_glob
      do j = 1 , 4
        aux_x(j,:) = x_vec(:)
        if(j < 3) aux_x(j,i) = x_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
        if(j >= 3) aux_x(j,i) = x_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
        call sub_N_par(plasma_params, aux_x(j,:), N_vec, N_par_aux(j), N_abs_aux,  B_abs_aux)
        if(.not. plasma_params%Te_ne_mat) then
          rhop = func_rhop(plasma_params, aux_x(j,:))
          call retrieve_n_e(plasma_params, rhop, n_e)
          call retrieve_T_e(plasma_params, rhop, T_e)
        else
          call retrieve_n_e_mat_single(plasma_params, aux_x(j,:), n_e)
          call retrieve_T_e_mat_single(plasma_params, aux_x(j,:), T_e)
        end if
        X_aux(j) = func_X(plasma_params, omega, n_e, T_e)
        Y_aux(j) = func_Y(plasma_params, omega, B_abs_aux, T_e)
      end do
      call sub_remap_coords(x_vec(1:3), R_vec)
      print*, "Current position", R_vec
      print*, "Current rhop", rhop_out
      print*, "X_aux", X_aux
      print*, "Y_aux", Y_aux
      print*, "N_par_aux", N_par_aux
      print*, "dX/dx ana", spatial_grad_X(i)
      print*, "dX/dx num fwd. :", (X_aux(2) - X)/h_x
      print*, "dX/dx num ctr. :" , (- X_aux(1) + 8.d0 *  X_aux(2) - &
                          8.d0 *  X_aux(3) +   X_aux(4)) / (12.d0 *h_x)
      print*, "dY^2/dx ana",spatial_grad_Y2(i)
      print*, "dY^2/dx num fwd.", (Y_aux(2)**2 - Y**2)/h_x
      print*, "dY^2/dx num ctr. :" ,(- Y_aux(1)**2 + 8.d0 *  Y_aux(2)**2 - &
                          8.d0 *  Y_aux(3)**2 +   Y_aux(4)**2) / (12.d0 *h_x)
      print*, "dN_par/dx ana",spatial_grad_N_par(i)
      print*, "dN_par/dx num fwd.", (N_par_aux(2) - N_par)/h_x
      print*, "dN_par/dx num ctr. :" , (- N_par_aux(1) + 8.d0 *  N_par_aux(2) - &
                          8.d0 *  N_par_aux(3) +   N_par_aux(4)) / (12.d0 *h_x)
      print*, "dLambda/dX ana", func_dNs_sq_dX(N_abs,  X, Y, N_par, mode) ! derivative dN/dX = 0
      print*, "dLambda/dX num fwd", ( func_Lambda(N_abs, X + h_x, Y, N_par, mode) - N_s) / h_x
      print*, "dLambda/dX num ctr", (- func_Lambda(N_abs, X +  2.d0 * h_x, Y, N_par, mode) + &
                              8.d0 *  func_Lambda(N_abs, X +  h_x, Y, N_par, mode) - &
                          8.d0 *  func_Lambda(N_abs, X - h_x, Y, N_par, mode)  &
                          +   func_Lambda(N_abs, X - 2.d0 * h_x, Y, N_par, mode)) / (12.d0 *h_x)
      print*, "dLambda/dY2 ana", func_dNs_sq_dY2(N_abs,  X, Y, N_par, mode)! derivative dN/dY^2 = 0
      print*, "dLambda/dY2 num fwd", ( func_N_s_2( X, Y + h_x, N_par, mode) - N_s) / (h_x * 2.d0 * Y)
      print*, "dLambda/dY2 num ctr",(- func_Lambda(N_abs, X, Y + 2.d0 * h_x, N_par, mode) + &
                              8.d0 *  func_Lambda(N_abs, X, Y + h_x, N_par, mode) - &
                          8.d0 *  func_Lambda(N_abs, X, Y - 1.d0 * h_x, N_par, mode)  &
                          +   func_Lambda(N_abs, X, Y - 2.d0 * h_x, N_par, mode)) / (12.d0 *h_x * 2.d0 * Y)
      print*, "dNs_sq/dN_par ana", func_dNs_sq_dN_par(N_abs,  X, Y, N_par, mode)
      print*, "dNs_sq/dN_par num fwd", ( func_N_s_2(X, Y, N_par + h_x, mode) - N_s) / h_x
      print*, "dNs_sq/dN_par num ctr", (- func_N_s_2(X, Y, N_par + 2.d0 * h_x, mode) + &
                              8.d0 *  func_N_s_2(X, Y, N_par + h_x, mode) - &
                          8.d0 *  func_N_s_2(X, Y,  N_par - h_x, mode)  &
                          +   func_N_s_2(X, Y,  N_par - 2.d0 * h_x, mode)) / (12.d0 *h_x)
      print*,"ana", dxds(i + 3)
      print*,"num fwd",(func_Lambda(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - Hamil) / h_x
      print*, "num backward",  -(func_Lambda(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode) - Hamil) / h_x
      print*, "num ctr:", (- func_Lambda(N_abs, X_aux(1), Y_aux(1), N_par_aux(1), mode) + &
                              8.d0 *  func_Lambda(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - &
                          8.d0 *  func_Lambda(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode)  &
                          +   func_Lambda(N_abs, X_aux(4), Y_aux(4), N_par_aux(4), mode)) / (12.d0 *h_x)
!        spatial_grad_Lambda(i) = (- func_Lambda(N_abs, X_aux(1), Y_aux(1), N_par_aux(1), mode) + &
!                               8.d0 *  func_Lambda(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - &
!                               8.d0 *  func_Lambda(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode) + &
!                               func_Lambda(N_abs, X_aux(4), Y_aux(4), N_par_aux(4), mode)) / (12.d0 *h_x)
    end do
    do i = 1,3
      do j = 1 , 4
        h_x = h_x_glob
        aux_N(j,:) = N_vec(:)
        if(j < 3) aux_N(j,i) = N_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
        if(j >= 3) aux_N(j,i) = N_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
        call sub_N_par(plasma_params,x_vec, aux_N(j,:), N_par_aux(j), N_abs_aux_2(j),  B_abs_aux)
      end do
      print*, "dN_par/dN ana", N_grad_N_par(i)
      print*, "dN_par/dN num forward", (N_par_aux(2) - N_par) / h_x
      print*, "dN_par/dN Numerical central:" , (- N_par_aux(1) + 8.d0 * N_par_aux(2) - &
                          8.d0 * N_par_aux(3) +  N_par_aux(4)) / (12.d0 *h_x)
      print*, "dNs_sq/dN_par ana", func_dNs_sq_dN_par(N_abs,  X, Y, N_par, mode)
      print*,"dNs_sq/dN_par num forward",(func_N_s_2(X, Y,  N_par + h_x, mode) - N_s) / h_x
      print*, "dNs_sq/dN_par num ctr", (- func_N_s_2(X, Y,  N_par + 2.d0 * h_x, mode) + &
                            8.d0 *  func_N_s_2(X, Y,  N_par + h_x, mode) - &
                        8.d0 *  func_N_s_2(X, Y,  N_par - h_x, mode)  &
                        +   func_N_s_2(X, Y,  N_par - 2.d0 * h_x, mode)) / (12.d0 *h_x)
      print*, "ana", dxds(i)
      print*, "num forward",  (func_Lambda(N_abs_aux_2(2), X, Y, N_par_aux(2), mode) - Hamil) / h_x
      print*, "num backward",  -(func_Lambda(N_abs_aux_2(3), X, Y, N_par_aux(3), mode) - Hamil) / h_x
      print*,"num central", (- func_Lambda(N_abs_aux_2(1), X, Y, N_par_aux(1), mode)+ 8.d0 * func_Lambda(N_abs_aux_2(2), X, Y, N_par_aux(2), mode)- &
                          8.d0 * func_Lambda(N_abs_aux_2(3), X, Y, N_par_aux(3),  mode) + func_Lambda(N_abs_aux_2(4), X, Y, N_par_aux(4), mode)) / (12.d0 *h_x)
!      N_grad_Lambda(i) = (- func_Lambda(N_abs_aux(1), X, Y, N_par_aux(1), mode)+ 8.d0 * func_Lambda(N_abs_aux(2), X, Y, N_par_aux(2), mode)- &
!                            8.d0 * func_Lambda(N_abs_aux(3), X, Y, N_par_aux(3),  mode) + func_Lambda(N_abs_aux(4), X, Y, N_par_aux(4), mode)) / (12.d0 *h_x)
    end do
    print*, " Critical error in raytracing"
    call abort
    !end if
  end subroutine sub_grad_Lambda


  subroutine sub_grad_Lambda_star(plasma_params, omega, mode, x_vec, dxds)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, h_x_glob
    USE ripple3d,                 only: grad_type
    USE constants,                 only : c0, eps0, mass_e, e0
    use mod_ecfm_refr_utils, only: retrieve_n_e, retrieve_T_e, &
                                   retrieve_n_e_mat_single, &
                                   retrieve_T_e_mat_single
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params ! only slice of the entire ray
    real(rkind),                 intent(in)   :: omega
    integer(ikind),              intent(in)   :: mode
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind), dimension(:),   intent(out)  :: dxds
    real(rkind)                               :: N_abs, N_par, X, Y, rhop_out
    real(rkind), dimension(3)                 :: N_vec, spatial_grad_X, spatial_grad_Y, spatial_grad_Y2, &
                                                 spatial_grad_N_par, spatial_grad_B_abs, N_grad_N_par
    real(rkind)                               :: Hamil, N_s, B_abs, n_e, T_e, rhop, N_abs_aux, &
                                                 h_x, dN_s_star_2_dN_par, B_abs_aux
    real(rkind), dimension(4)                 :: X_aux, Y_aux, N_par_aux, N_abs_aux_2
    real(rkind), dimension(4,3)               :: aux_x, aux_N
    integer(ikind)                            :: i, j
    N_abs = sqrt(x_vec(4)**2 + x_vec(5)**2 + x_vec(6)**2)
    N_vec = x_vec(4:6)
    call sub_spatial_grad_X(plasma_params, omega, x_vec(1:3), X, spatial_grad_X, rhop_out)
    call sub_spatial_grad_N_par(plasma_params, x_vec(1:3), x_vec(4:6), N_abs, N_par, spatial_grad_N_par, spatial_grad_B_abs, B_abs, N_grad_N_par)
    call sub_spatial_grad_Y(plasma_params, omega, x_vec(1:3), B_abs, spatial_grad_B_abs, Y, spatial_grad_Y)
    Hamil = func_Lambda_star(N_abs, X, Y,  N_par, mode)
    spatial_grad_Y2 = 2.d0 * Y * spatial_grad_Y
    dN_s_star_2_dN_par = func_dN_s_star_2_dN_par(N_abs,  X, Y, N_par, mode)
    dxds(4:6) = -spatial_grad_N_par(:) * dN_s_star_2_dN_par
    dxds(4:6) = dxds(4:6) + spatial_grad_X(:) * func_dLambda_star_dX(N_abs,  X, Y, N_par, mode)
    dxds(4:6) = dxds(4:6) + spatial_grad_Y2(:) * func_dLambda_star_dY2(N_abs,  X, Y, N_par, mode)
    dxds(1:3) = 2.d0  * N_vec(:) * (2.d0 * (-1.d0 + X + Y**2))**2 - N_grad_N_par(:) * dN_s_star_2_dN_par
    if((debug_level >= 3 .or. (debug_level ==  2 .and.  sqrt(x_vec(1)**2 + x_vec(2)**2) <= 2.13)))then! .or. abs(H) > 1.d-2) then
      print*, "/---------------DEBUG OUTPUT------------------\"
    else
      return
    end if
    Hamil = func_Lambda_star(N_abs, X, Y,  N_par, mode)
    print*, "Hamil", Hamil
    print*,"X, Y, N_par, Delta, N_abs,", X, Y, N_par, func_delta(X, Y, N_par), N_abs
    N_s = sqrt(func_N_s_2(X, Y, N_par, mode))
   ! if(sqrt(x_vec(1)**2 + x_vec(2)**2) < 2.3) then
      do i =1,3
        h_x = h_x_glob
        do j = 1 , 4
          aux_x(j,:) = x_vec(:)
          if(j < 3) aux_x(j,i) = x_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
          if(j >= 3) aux_x(j,i) = x_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
          call sub_N_par(plasma_params, aux_x(j,:), N_vec, N_par_aux(j), N_abs_aux, B_abs_aux)
          if(.not. plasma_params%Te_ne_mat) then
            rhop = func_rhop(plasma_params, aux_x(j,:))
            call retrieve_n_e(plasma_params, rhop, n_e)
            call retrieve_T_e(plasma_params, rhop, T_e)
          else
            call retrieve_n_e_mat_single(plasma_params,  aux_x(j,:), n_e)
            call retrieve_T_e_mat_single(plasma_params, aux_x(j,:), T_e)
          end if
          X_aux(j) = func_X(plasma_params,omega, n_e, T_e)
          Y_aux(j) = func_Y(plasma_params,omega, B_abs_aux, T_e)
        end do
        print*, "dX/dx ana", spatial_grad_X(i)
        print*, "dX/dx num fwd. :", (X_aux(2) - X)/h_x
        print*, "dX/dx num ctr. :" , (- X_aux(1) + 8.d0 *  X_aux(2) - &
                            8.d0 *  X_aux(3) +   X_aux(4)) / (12.d0 *h_x)
        print*, "dY^2/dx ana",spatial_grad_Y2(i)
        print*, "dY^2/dx num fwd.", (Y_aux(2)**2 - Y**2)/h_x
        print*, "dY^2/dx num ctr. :" ,(- Y_aux(1)**2 + 8.d0 *  Y_aux(2)**2 - &
                            8.d0 *  Y_aux(3)**2 +   Y_aux(4)**2) / (12.d0 *h_x)
        print*, "dN_par/dx ana",spatial_grad_N_par(i)
        print*, "dN_par/dx num fwd.", (N_par_aux(2) - N_par)/h_x
        print*, "dN_par/dx num ctr. :" , (- N_par_aux(1) + 8.d0 *  N_par_aux(2) - &
                            8.d0 *  N_par_aux(3) +   N_par_aux(4)) / (12.d0 *h_x)
        print*, "dLambda*/dX ana", func_dLambda_star_dX(N_abs,  X, Y, N_par, mode)
        print*, "dLambda*/dX num fwd", ( func_Lambda_star(N_abs, X + h_x, Y, N_par, mode) - N_s) / h_x
        print*, "dLambda*/dX num ctr", (- func_Lambda_star(N_abs, X +  2.d0 * h_x, Y, N_par, mode) + &
                                8.d0 *  func_Lambda_star(N_abs, X +  h_x, Y, N_par, mode) - &
                            8.d0 *  func_Lambda_star(N_abs, X - h_x, Y, N_par, mode)  &
                            +   func_Lambda_star(N_abs, X - 2.d0 * h_x, Y, N_par, mode)) / (12.d0 *h_x)
        print*, "dLambda*/dY2 ana", func_dLambda_star_dY2(N_abs,  X, Y, N_par, mode)
        print*, "dLambda*/dY2 num fwd", ( func_N_s_star_2( X, Y + h_x, N_par, mode) - N_s) / (h_x * 2.d0 * Y)
        print*, "dLambda*/dY2 num ctr",(- func_Lambda_star(N_abs, X, Y + 2.d0 * h_x, N_par, mode) + &
                                8.d0 *  func_Lambda_star(N_abs, X, Y + h_x, N_par, mode) - &
                            8.d0 *  func_Lambda_star(N_abs, X, Y - 1.d0 * h_x, N_par, mode)  &
                            +   func_Lambda_star(N_abs, X, Y - 2.d0 * h_x, N_par, mode)) / (12.d0 *h_x * 2.d0 * Y)
        print*, "dN_s*^2/dN_par ana", func_dN_s_star_2_dN_par(N_abs,  X, Y, N_par, mode)
        print*, "dN_s*^2/dN_par num fwd", ( func_N_s_star_2(X, Y, N_par + h_x, mode) - N_s) / h_x
        print*, "dN_s*^2/dN_par num ctr", (- func_N_s_star_2(X, Y, N_par + 2.d0 * h_x, mode) + &
                                8.d0 *  func_N_s_star_2(X, Y, N_par + h_x, mode) - &
                            8.d0 *  func_N_s_star_2(X, Y, N_par - h_x, mode)  &
                            +   func_N_s_star_2(X, Y, N_par - 2.d0 * h_x, mode)) / (12.d0 *h_x)
        print*,"ana", dxds(i + 3)
        print*,"num fwd",(func_Lambda_star(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - Hamil) / h_x
        print*, "num backward",  -(func_Lambda_star(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode) - Hamil) / h_x
        print*, "num ctr:", (- func_Lambda_star(N_abs, X_aux(1), Y_aux(1), N_par_aux(1), mode) + &
                                8.d0 *  func_Lambda_star(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - &
                            8.d0 *  func_Lambda_star(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode)  &
                            +   func_Lambda_star(N_abs, X_aux(4), Y_aux(4), N_par_aux(4), mode)) / (12.d0 *h_x)
!        spatial_grad_Lambda(i) = (- func_Lambda(N_abs, X_aux(1), Y_aux(1), N_par_aux(1), mode) + &
!                               8.d0 *  func_Lambda(N_abs, X_aux(2), Y_aux(2), N_par_aux(2), mode) - &
!                               8.d0 *  func_Lambda(N_abs, X_aux(3), Y_aux(3), N_par_aux(3), mode) + &
!                               func_Lambda(N_abs, X_aux(4), Y_aux(4), N_par_aux(4), mode)) / (12.d0 *h_x)
      end do
      do i = 1,3
        do j = 1 , 4
          h_x = h_x_glob
          aux_N(j,:) = N_vec(:)
          if(j < 3) aux_N(j,i) = N_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
          if(j >= 3) aux_N(j,i) = N_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
          call sub_N_par(plasma_params,x_vec, aux_N(j,:), N_par_aux(j), N_abs_aux_2(j),  B_abs_aux)
        end do
        print*, "dN_par/dN ana",N_grad_N_par(i)
        print*, "dN_par/dN num forward", (N_par_aux(2) - N_par) / h_x
        print*, "dN_par/dN Numerical central:" , (- N_par_aux(1) + 8.d0 * N_par_aux(2) - &
                            8.d0 * N_par_aux(3) +  N_par_aux(4)) / (12.d0 *h_x)
        print*, "dLambda*/dN_par ana", func_dN_s_star_2_dN_par(N_abs,  X, Y, N_par, mode)
        print*,"dLambda*/dN_par num forward",(func_N_s_star_2(X, Y, N_par + h_x, mode) - N_s) / h_x
        print*, "dN_s*/dN_par num ctr", (- func_N_s_star_2(X, Y, N_par + 2.d0 * h_x, mode) + &
                              8.d0 *  func_N_s_star_2(X, Y, N_par + h_x, mode) - &
                          8.d0 *  func_N_s_star_2( X, Y, N_par - h_x, mode)  &
                          +   func_N_s_star_2(X, Y, N_par - 2.d0 * h_x, mode)) / (12.d0 *h_x)
        print*, "ana", dxds(i)
        print*, "num forward",  (func_Lambda_star(N_abs_aux_2(2), X, Y, N_par_aux(2), mode) - Hamil) / h_x
        print*, "num backward",  -(func_Lambda_star(N_abs_aux_2(3), X, Y, N_par_aux(3), mode) - Hamil) / h_x
        print*,"num central",(- func_Lambda_star(N_abs_aux_2(1), X, Y, N_par_aux(1), mode)+ 8.d0 * func_Lambda_star(N_abs_aux_2(2), X, Y, N_par_aux(2), mode)- &
                            8.d0 * func_Lambda_star(N_abs_aux_2(3), X, Y, N_par_aux(3),  mode) + func_Lambda_star(N_abs_aux_2(4), X, Y, N_par_aux(4), mode)) / (12.d0 *h_x)
!      N_grad_Lambda(i) = (- func_Lambda(N_abs_aux(1), X, Y, N_par_aux(1), mode)+ 8.d0 * func_Lambda(N_abs_aux(2), X, Y, N_par_aux(2), mode)- &
!                            8.d0 * func_Lambda(N_abs_aux(3), X, Y, N_par_aux(3),  mode) + func_Lambda(N_abs_aux(4), X, Y, N_par_aux(4), mode)) / (12.d0 *h_x)
      end do
      stop "ok?"
    !end if
  end subroutine sub_grad_Lambda_star

  subroutine f_H(m, sigma, x_vec, dx_dsigma)
    use f90_kind
    implicit none
    integer          , intent(in)           :: m
    real(rkind),  intent(in)                :: sigma
    real(rkind),  dimension(m), intent(in)  :: x_vec
    real(rkind),  dimension(m), intent(inout) :: dx_dsigma
    call sub_grad_H(glob_plasma_params, glob_omega, glob_mode, x_vec, dx_dsigma)
  end subroutine f_H

  subroutine f_Lambda(m, s, x_vec, dxds)
    use f90_kind
    USE mod_ecfm_refr_types, only : Lambda_star
    implicit none
    integer          , intent(in)           :: m
    real(rkind)                             :: s
    real(rkind),  dimension(m), intent(in)  :: x_vec
    real(rkind),  dimension(m), intent(out) :: dxds
    real(rkind)                             :: dxds_abs
    !print*, "current m", m
    if(Lambda_star) then
      call sub_grad_Lambda_star(glob_plasma_params, glob_omega, glob_mode, x_vec, dxds)
    else
      call sub_grad_Lambda(glob_plasma_params, glob_omega, glob_mode, x_vec, dxds)
    end if
    dxds_abs = sqrt(dxds(1)**2 + dxds(2)**2 + dxds(3)**2)
    dxds(1:3) = dxds(1:3) / dxds_abs
    dxds(4:6) = -dxds(4:6) / dxds_abs
!    if (any(x_vec /= x_vec) .or. any(dxds /= dxds)) then
!        glob_debug_level = 3
!        if(lambda_star) then
!          call sub_grad_Lambda_star(glob_plasma_params, glob_omega, x_vec, dxds)
!        else
!          call sub_grad_Lambda(glob_plasma_params, glob_omega, x_vec, dxds)
!        end if
!     end if
  end subroutine f_Lambda

  subroutine Jac(nq,path, y_vec,ml,mu,dfdy,npowpd)
    use f90_kind
    implicit none
    integer,   intent(in)   :: nq,ml,mu,npowpd
    real(rkind), intent(in) :: path
    real(rkind), intent(inout) :: y_vec(:), dfdy(:,:)
     ! empty subroutine for the ODE-solver - if implicit scheme is used use finite differences
    return
 end subroutine Jac

  subroutine sub_single_step_LSODE(sigma, x_vec, N_vec, h, x_vec_out, N_vec_out, B_vec_out, &
              theta_out, Hamil_out, N_s_out, n_e_out, omega_c_out, T_e_out, rhop_out, v_g_perp, &
              work_lsode, iwork_lsode, istate)
    USE f90_kind
    USE mod_ecfm_refr_types, only : Hamil, Lambda_star
    USE constants,                 only : c0, eps0, mass_e, e0
    implicit None
    real(rkind)                , intent(inout):: sigma
    real(rkind), dimension(:),   intent(in)   :: x_vec, N_vec
    real(rkind),                 intent(inout):: h
    real(rkind), dimension(:),   intent(out)  :: x_vec_out, N_vec_out, B_vec_out
    real(rkind),                 intent(out)  :: theta_out, Hamil_out, N_s_out, n_e_out, omega_c_out, T_e_out, rhop_out, v_g_perp
    real(rkind), dimension(:),   intent(in)   :: work_lsode
    integer(ikind), dimension(:), intent(in)  :: iwork_lsode
    integer(ikind),              intent(inout):: istate
    real(rkind)                               :: X, Y, A, B, C, N_par, dX1, dX2, sigma0
    integer, dimension(1)                     :: neq
    double precision, dimension(1)            :: rtol, atol
    real(rkind), dimension(6)                 :: y_vec, dy_vec_dummy, y_vec_init
    neq = 6
    rtol = 1d-6
    atol = 1d-6
    y_vec(1:3) = x_vec
    y_vec(4:6) = N_vec
    y_vec_init = y_vec
    sigma0 = sigma
    if(Hamil == "Dani") then
      call dlsode(f_Lambda, neq,y_vec,sigma,sigma + h,       &
                  1, rtol, atol,4,istate,1, &
                  work_lsode,size(work_lsode),iwork_lsode,size(iwork_lsode),Jac,10) !22
      if(Lambda_star) then
        call sub_grad_Lambda_star(glob_plasma_params, glob_omega, glob_mode, y_vec_init, dy_vec_dummy)
        dX1 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
        call sub_grad_Lambda_star(glob_plasma_params, glob_omega, glob_mode, y_vec, dy_vec_dummy)
        dX2 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
      else
        call sub_grad_Lambda(glob_plasma_params, glob_omega, glob_mode, y_vec_init, dy_vec_dummy)
        dX1 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
        call sub_grad_Lambda(glob_plasma_params, glob_omega, glob_mode, y_vec, dy_vec_dummy)
        dX2 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
      end if
    else
      call dlsode(f_H,neq,y_vec,sigma,sigma + h,       &
                  1, rtol, atol,4,istate,1, &
                  work_lsode,size(work_lsode),iwork_lsode,size(iwork_lsode),Jac,10) !22
      call sub_grad_H(glob_plasma_params, glob_omega, glob_mode, y_vec_init, dy_vec_dummy)
      dX1 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
      call sub_grad_H(glob_plasma_params, glob_omega, glob_mode, y_vec, dy_vec_dummy)
      dX2 = sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2 + &
                 dy_vec_dummy(4)**2 + dy_vec_dummy(5)**2 + dy_vec_dummy(6)**2)
    end if

    !call sub_spatial_grad_X(glob_plasma_params, glob_omega, x_vec, X,  spatial_grad_X, rhop_out)
    !spatial_grad_X = spatial_grad_X * glob_omega / (2.d0 * sqrt(X)) ! We want just dw_p/dx here
    !dX1 = sqrt(spatial_grad_X(1)**2 + spatial_grad_X(2)**2 + spatial_grad_X(3)**2)
    x_vec_out = y_vec(1:3)
    N_vec_out = y_vec(4:6)
    !call sub_spatial_grad_X(glob_plasma_params, glob_omega, x_vec_out, X,  spatial_grad_X, rhop_out)
    !spatial_grad_X = spatial_grad_X * glob_omega / (2.d0 * sqrt(X)) ! We want just dw_p/dx here
    !dX2 = sqrt(spatial_grad_X(1)**2 + spatial_grad_X(2)**2 + spatial_grad_X(3)**2)
    h = h * abs(dX1 / dX2)
    if(h < glob_plasma_params%h_min) h = glob_plasma_params%h_min
    if(h > glob_plasma_params%h_max) h = glob_plasma_params%h_max
    call sub_local_params(glob_plasma_params, glob_omega, x_vec_out, N_vec_out, B_vec_out, N_s_out, n_e_out, omega_c_out, T_e_out, theta_out, rhop_out)
    X = func_X(glob_plasma_params, glob_omega, n_e_out,T_e_out)
    Y = func_Y(glob_plasma_params, glob_omega, omega_c_out * mass_e / e0, T_e_out)
    N_par = cos(theta_out) * N_s_out
    v_g_perp = 1.d0/sqrt(dy_vec_dummy(1)**2 + dy_vec_dummy(2)**2 + dy_vec_dummy(3)**2)
    if(Hamil == "Dani") then
      if(lambda_star) then
        Hamil_out = func_Lambda_star(N_s_out, X, Y, N_par, glob_mode)
      else
        Hamil_out = func_Lambda(N_s_out, X, Y, N_par, glob_mode)
      end if
    else
      A = func_A(X, Y)
      B = func_B(X, Y, N_par)
      C = func_C(X, Y, N_par)
      Hamil_out = func_H(sqrt(N_s_out**2 - N_par**2), A, B, C, glob_mode)
    end if
    !if(abs(dX1 / dX2  - 1) > 1.d-2) then
    !  print*, "Changed h by",  dX1 / dX2, "h", h, "H", Hamil_out, "rhop", rhop_out
    !end if
    if (istate /= 2 .or. abs(Hamil_out) > 1.0e+1 .or. any(y_vec /= y_vec) .or. Hamil_out /= Hamil_out) then
      if(istate /= 2 ) then
        print'(a,i4,a,2e16.8)','WARNING from DLSODE: istate =',istate, &
                              '   p_got =',sigma0, sigma0 + h
      else
        print*, "Very large Hamiltonian encountered, rays most likely highly inaccurate", Hamil_out
      end if
      print*, "Entering Debug mode"
      print*, "Current step size", h
      print*, "Final position", x_vec_out
      debug_level = 3
      istate= 0
      y_vec(1:3) = x_vec
      y_vec(4:6) = N_vec
      if(Hamil == "Dani") then
        if(lambda_star) then
          call sub_grad_lambda_star(glob_plasma_params,glob_omega, glob_mode, y_vec, dy_vec_dummy)
        else
          call sub_grad_lambda(glob_plasma_params,glob_omega, glob_mode, y_vec, dy_vec_dummy)
        end if
      else
          call sub_grad_H(glob_plasma_params,glob_omega, glob_mode, y_vec, dy_vec_dummy)
      end if
    end if
  end subroutine sub_single_step_LSODE

  subroutine sub_local_params(plasma_params, omega, x_vec, N_vec, B_vec, N_abs, n_e, omega_c, T_e, theta, rhop_out)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, &
                                    straight, SOL_ne, SOL_Te !, h_x_glob
    USE ripple3d,                 only: grad_type
    USE constants,                 only : e0, mass_e
    use mod_ecfm_refr_utils, only: retrieve_n_e, retrieve_T_e, &
                                   retrieve_T_e_mat_single, retrieve_n_e_mat_single
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params ! only slice of the entire ray
    real(rkind),                 intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec, N_vec
    real(rkind), dimension(:),   intent(out)  :: B_vec
    real(rkind),                 intent(out)  :: N_abs, theta, n_e, omega_c, T_e, rhop_out
    real(rkind)                               :: B_abs
    if(any(x_vec /= x_vec)) then
      print*, "nan in x_vec"
      stop "Nan in x_vec in sub_local_params in mod_raytrace.f90"
    end if
    call sub_theta(plasma_params, x_vec, N_vec, theta, N_abs, B_abs, B_vec)
    omega_c = e0 * B_abs /mass_e
    rhop_out = func_rhop(plasma_params, x_vec)
    if(rhop_out > plasma_params%rhop_max .or. &
       rhop_out > plasma_params%rhop_entry .or. &
       rhop_out == -1.d0) then
      n_e = SOL_ne
      T_e = SOL_Te
    else
      if(plasma_params%No_ne_te) then
        if(.not. straight) then
          print*, "Ne and Te must be present if raytracing is considered"
          call abort
        end if
        n_e = 0.d0
        T_e = 0.d0
      else if(.not. plasma_params%Te_ne_mat) then
        call retrieve_n_e(plasma_params, rhop_out, n_e)
        call retrieve_T_e(plasma_params, rhop_out, T_e)
      else
        call retrieve_n_e_mat_single(plasma_params, x_vec, n_e)
        call retrieve_T_e_mat_single(plasma_params, x_vec, T_e)
      end if
    end if
  end subroutine sub_local_params

  function reflect_off_surface(N_vec, n_surf, theta)
  USE f90_kind
  real(rkind), dimension(:),   intent(in)  :: N_vec, n_surf
  real(rkind), intent(in), optional        :: theta
  real(rkind), dimension(3)                :: reflect_off_surface
  real(rkind), dimension(3,3)              :: R, uxu, ux
  integer(ikind)                           :: i,j
  real(rkind)                              :: sin_theta, cos_theta
  reflect_off_surface(:) = 0.d0
  if(present(theta)) then
    sin_theta = sin(theta)
    cos_theta = cos(theta)
  else
    sin_theta = 1.d0
    cos_theta = 0.d0
  end if
  do i =1,3
    ux(:,i) = 0.d0
    R(:,i) = 0.d0
    do j= 1,3
      if(i ==j) R(i,j) = cos_theta
      uxu(j,i) = n_surf(i) * n_surf(j)
    end do
  end do
  ux(1,2) = -n_surf(3)
  ux(2,1) = n_surf(3)
  ux(3,1) = -n_surf(2)
  ux(1,3) = n_surf(2)
  ux(2,3) = -n_surf(1)
  ux(3,2) = n_surf(1)
  do j =1,3
    R(:,j) = R(:,j) + sin_theta * ux(:,j) + (1.d0 - cos_theta) * uxu(:,j)
  end do
  do j =1,3
    reflect_off_surface(j) = sum(N_vec(:) * R(j,:))
  end do
  end function

  function make_Snells_refraction(N_abs, total_reflection)
  USE f90_kind
  use mod_ecfm_refr_utils,  only: sub_remap_coords
  use constants,            only: pi
  implicit none
  real(rkind),                 intent(in)  :: N_abs
  logical, intent(out), optional           :: total_reflection
  real(rkind), dimension(3)                :: make_Snells_refraction
  real(rkind)                              :: ratio, scalar_k_surf, cos_phi_tok, sin_phi_tok, abs_norm_vec
  real(rkind), dimension(3)                :: R_vec, flux_norm_vec, flux_vec_R, N_vec_norm
  !print*, "Before Snell's law", N_loc_vec
  call sub_remap_coords(x_loc_vec, R_vec)
  flux_vec_R = func_flux_norm_vec(glob_plasma_params, x_loc_vec)
  if(flux_vec_R(1) < 0.d0) flux_vec_R = -flux_vec_R
  ! From cylindrical to carthesian coordinates
  cos_phi_tok = cos(R_vec(2))
  sin_phi_tok = sin(R_vec(2))
  flux_norm_vec(1) = flux_vec_R(1) * cos_phi_tok - flux_vec_R(2) * sin_phi_tok
  flux_norm_vec(2) = flux_vec_R(1) * sin_phi_tok + flux_vec_R(2) * cos_phi_tok
  flux_norm_vec(3) = flux_vec_R(3)
  abs_norm_vec = sqrt(flux_norm_vec(1)**2 + flux_norm_vec(2)**2 + flux_norm_vec(3)**2)
  flux_norm_vec = flux_norm_vec / abs_norm_vec
  if(N_abs <= 0.d0) then
  ! If the first point in the plasma is evanescent all we can do is reflect off it
    if(present(total_reflection)) total_reflection = .True.
    make_Snells_refraction = reflect_off_surface(N_loc_vec, flux_norm_vec)
    return
  end if
  !print*, "Flux surface normal", x_loc_vec / sqrt(x_loc_vec(1)**2 + x_loc_vec(2)**2 + x_loc_vec(3)**2), flux_norm_vec

  !flux_norm_vec = -flux_norm_vec
  ratio = 1.d0 / N_abs
  !print*, "Ratio", ratio
  N_vec_norm = N_loc_vec
  !print*, "|N_vec_norm|", sqrt(N_loc_vec(1)**2 + N_loc_vec(2)**2 + N_loc_vec(3)**2)
  scalar_k_surf =  - dot_product(flux_norm_vec, N_vec_norm)
  !print*, "k . n", scalar_k_surf
  !print*, "correction weight", ( ratio * scalar_k_surf - &
  !                        sqrt(1.d0 - ratio**2 * ( 1.d0 - scalar_k_surf**2)))
!  if(scalar_k_surf < 0.d0) then
!    flux_norm_vec = -flux_norm_vec
!    scalar_k_surf =  - dot_product(flux_norm_vec, N_vec_norm)
!  end if
  if(1.d0 - ratio**2 * ( 1.d0 - scalar_k_surf**2) < 0.d0 ) then
    make_Snells_refraction = reflect_off_surface(N_loc_vec, flux_norm_vec)
    if(present(total_reflection)) total_reflection = .True.
    return
  end if
  if(present(total_reflection)) total_reflection = .False.
  make_Snells_refraction = ratio * N_vec_norm + ( ratio * scalar_k_surf - &
                          sqrt(1.d0 - ratio**2 * ( 1.d0 - scalar_k_surf**2))) * flux_norm_vec
  !print*, "|N_norm| after Snell", &
  !  sqrt(make_Snells_refraction(1)**2 + make_Snells_refraction(2)**2 + &
  !       make_Snells_refraction(3)**2)
  !print*, "After Snell's law", make_Snells_refraction
  end function make_Snells_refraction

  function make_H(N_abs)
    USE f90_kind
    USE mod_ecfm_refr_types, only : Hamil, lambda_star
    USE constants,                 only : mass_e, e0
    implicit none
    real(rkind),                 intent(in)  :: N_abs
    real(rkind)                              :: make_H
    real(rkind)                              :: X, Y, theta, A, B, C, rhop, omega_c, N_abs_aux, N_par, n_e, T_e
    real(rkind), dimension(3)                :: N_vec, B_vec
    N_vec = make_Snells_refraction(N_abs) * N_abs
    call sub_local_params(glob_plasma_params, glob_omega, x_loc_vec, N_vec, B_vec, N_abs_aux, n_e, omega_c, T_e, theta, rhop)
    X = func_X(glob_plasma_params,glob_omega, n_e, T_e)
    Y = func_Y(glob_plasma_params,glob_omega, omega_c * mass_e / e0, T_e)
    N_par = cos(theta) * N_abs
    if(Hamil == "Dani") then
      if(lambda_star) then
        make_H = func_Lambda_star(N_abs, X, Y, N_par, glob_mode)
        !print*, N_abs, make_H
      else
        make_H = func_Lambda(N_abs, X, Y, N_par, glob_mode)
      end if
    else
        A = func_A(X, Y)
        B = func_B(X, Y, N_par)
        C = func_C(X, Y, N_par)
        make_H = func_H(sqrt(N_abs**2 - N_par**2), A, B, C, glob_mode)
        if(make_H /= make_H) stop "nan with initial conditions"
        !print*, N_abs, make_H
    end if
    !print*, "H", make_H
  end function make_H

  function make_Hamil_Ns(X, Y, N_abs, theta, mode)
    USE f90_kind
    USE mod_ecfm_refr_types, only : Hamil, lambda_star, eps
    USE constants,                 only : mass_e, e0
    implicit none
    real(rkind),                 intent(in)  :: X, Y, N_abs, theta
    integer(ikind),              intent(in)  :: mode
    real(rkind)                              :: make_Hamil_Ns
    real(rkind)                              :: A, B, C, N_par
    N_par = cos(theta) * N_abs
    if(Hamil == "Dani") then
      if(lambda_star) then
        make_Hamil_Ns = func_N_s_star_2(X, Y, N_par, mode)
      else
        make_Hamil_Ns = func_N_s_2(X, Y, N_par, mode)
      end if
    else
        A = func_A(X, Y)
        B = func_B(X, Y, N_par)
        C = func_C(X, Y, N_par)
        make_Hamil_Ns =((B - real(mode,8) * sqrt(B**2 - 4.d0 * A * C)) * (A / (A**2 + eps**2))) / ( 2.d0 * sin(theta))
        !print*, N_abs, make_H
    end if
    !print*, "H", make_H
  end function make_Hamil_Ns

  function make_H_wrapper(N_abs)
    USE f90_kind
    implicit none
    real(rkind), dimension(:),  intent(in)  :: N_abs
    real(rkind)                :: make_H_wrapper
    make_H_wrapper = make_H(N_abs(1))**2 ! need positive value for minimization algorithm
  end function make_H_wrapper

  subroutine sub_calculate_initial_N(plasma_params, omega, mode, x_vec, N_vec, H_init, total_reflection)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, output_level
    USE mod_ecfm_refr_utils, only : BrentRoots
#ifdef NAG
    USE nag_nlin_eqn,             only : nag_nlin_eqn_sol
#endif
    use nr_mod, only                        : powell
    USE constants,                 only : mass_e, e0, pi
    implicit none
    type(plasma_params_type)          , intent(in)                   :: plasma_params
    real(rkind),                 intent(in)                          :: omega
    integer(ikind),              intent(in)                          :: mode
    real(rkind), dimension(:),   intent(in)                          :: x_vec
    real(rkind), dimension(:),   intent(inout)                       :: N_vec
    real(rkind),                 intent(out)                         :: H_init
    logical, intent(out)                                             :: total_reflection
    real(rkind)                                                      :: N_abs, N_abs_lower, N_abs_upper, N_abs_0, X, Y, n_e, T_e, omega_c, theta,N_par, rhop_out
    real(rkind)                                                      :: H_lower, H_upper, H_lower_last, H_upper_last
    real(rkind)                                                      :: h
    real(rkind), dimension(3)                                        :: B_vec
    integer(ikind)                                                   :: n, iter, error
    call sub_local_params(plasma_params, omega, x_vec, N_vec, B_vec, N_abs_0, n_e, omega_c, T_e, theta, rhop_out)
    if(output_level .and. debug_level > 0) print*, "Rho pol and ne at first point in plasma", rhop_out, n_e
    N_par = cos(theta) * N_abs_0
    N_vec  = N_vec / N_abs_0
    h = 1.d-5
    x_loc_vec= x_vec
    N_loc_vec = N_vec
    X = func_X(plasma_params,omega, n_e, T_e)
    Y = func_Y(plasma_params,omega, omega_c * mass_e / e0, T_e)
    N_abs = func_N( X, Y, theta, mode)
    if(N_abs /= N_abs .or. X > 1.d0) then
      N_vec = make_Snells_refraction(0.d0, total_reflection)
      return
    end if
    N_abs = sqrt(make_Hamil_Ns(X, Y, N_abs, theta, mode))
    N_abs_lower = 0.999d0 * N_abs
    N_abs_upper = 1.001d0 * N_abs
!    print*, "N_abs range", N_abs_lower, N_abs,  N_abs_upper
!    print*, "Hamil range", make_H(N_abs_lower), make_H(N_abs), make_H(N_abs_upper)
    H_lower = make_H(N_abs_lower)
    H_upper = make_H(N_abs_upper)
    H_lower_last = 0.d0
    H_upper_last = 0.d0
    n = 0
    do while(H_lower * H_upper > 0.d0)
      if(mod(n,2) == 0) then
        N_abs_lower = 0.999d0 * N_abs_lower
        H_lower_last = H_lower
        H_lower = make_H(N_abs_lower)
      else
        N_abs_upper = 1.001d0 * N_abs_upper
        H_upper_last = H_upper
        H_upper = make_H(N_abs_upper)
      end if
      if(H_lower /= H_lower .or. H_upper /= H_upper .or. &
         H_lower == H_lower_last .or. H_upper == H_upper_last .or. &
         n > 10000) then
        print*, "X, Y", X, Y
        print*, "N_abs lower, N_abs_ray, N_abs, N_abs upper", N_abs_lower, N_abs_0, N_abs, N_abs_upper
        print*, "Hamiltonian lower, upper", H_lower, H_upper
        print*, "Last Hamiltonian lower, upper", H_lower_last, H_upper_last
        print*, "ne, Te, omega, omega_c, theta", n_e, T_e, omega, omega_c, theta / pi * 180
        print*, "x_vec", x_vec
        print*, "n", n
        stop "error in starting conditions for plasma coupling"
      end if
      n = n + 1
    end do
    N_abs = func_N( X, Y, theta, mode)
    N_abs =  BrentRoots( N_abs_lower, N_abs_upper, 1.d-12,  &
                         150, make_H,  &
                         H_init, iter, error )
    N_vec = make_Snells_refraction(N_abs, total_reflection) * N_abs
    if(total_reflection) return
    if(error == 1) then
      print*, "Could not solve vacuum plasma transition"
      print*, "No N_abs does not seem to be in the range between" , N_abs * 0.7d0 , 1.d0
      call abort()
    end if
#ifdef NAG
    if(debug_level > 0) then
      N_abs_opt(1) = N_abs
      call nag_nlin_eqn_sol(make_H, N_abs_opt(1), N_abs * 0.7d0, 1.d0, x_tol = 1.d-10, f_tol = 0.d0)
      if(output_level) print*, "N_abs Brent vs N_abs Nag", N_abs , N_abs_opt(1)
    end if
#endif
    call sub_local_params(plasma_params, omega, x_vec, N_vec, B_vec, N_abs, n_e, omega_c, T_e, theta, rhop_out)
    H_init = make_H(N_abs)
    if(error == 2 .and. abs(H_init) > 1.e-6) then
      print*, "Brent method did not reach target accuracy due to lack of iterations"
      print*, "H_init", H_init
      call abort()
    end if
    if(output_level) print*, "B_vec", B_vec
    if(output_level) print*,"Initial N_abs, H and cold N_abs", N_abs, H_init, func_N( X, Y, theta, mode)
    !N_vec = make_Snells_refraction(N_abs)
!    print*, "N_root stix", N_abs, make_H(N_abs)
!    Hamil = "Dani"
!    call nag_nlin_eqn_sol(make_H, N_abs, 0.9d0, 1.d0, x_tol = 1.d-13, f_tol = 0.d0)
!    N_vec = N_vec /N_abs_0 * N_abs
!    print*, "N_root Dani", N_abs, make_H(N_abs)
!    stop "compare"
!    call sub_local_params(plasma_params, omega, x_vec, N_vec, B_vec, N_abs, X, Y, N_par, rhop_out)
!    if(Hamil == "Dani") then
!        print*, "N_abs", N_abs, "Ns_sq", func_Lambda(N_abs, X, Y, N_par, mode)
!        stop "ok"
!    else
!        A = func_A(X, Y)
!        B = func_B(X, Y, N_par)
!        C = func_C(X, Y, N_par)
!        print*, "N_abs num,", N_abs, "H num", func_H(sqrt(N_abs**2 - N_par**2), A, B,C, mode)
!        stop "ok"
!    end if
    !stop "sensible launch ?"
  end subroutine sub_calculate_initial_N

  function func_within_plasma(plasma_params,x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, h_check
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type)      , intent(in)      :: plasma_params
    real(rkind), dimension(:)     , intent(in)      :: x_vec
    logical                                         :: func_within_plasma
    real(rkind), dimension(3)                       :: R_vec
    call sub_remap_coords(x_vec, R_vec)
    !print*,plasma_params%R_min, "<", R_tok, ">",plasma_params%R_max
    !print*,plasma_params%z_min, "<", z_tok, ">",plasma_params%z_max
    if(R_vec(1) - h_check > plasma_params%R_min  .and. R_vec(1) + 2.d0 * h_check < plasma_params%R_max) then
      if(R_vec(3) - h_check > plasma_params%z_min .and. R_vec(3) + 2.d0 * h_check < plasma_params%z_max) then
        if(func_rhop(plasma_params, x_vec) < plasma_params%rhop_max) then
          func_within_plasma = .true.
          return
        end if
      end if
    end if
    func_within_plasma = .false.
  end function func_within_plasma

  subroutine find_first_point_in_plasma(plasma_params, omega, mode, ray_segment, last_N, wall_hits, LOS_end)
  ! Straight line until we hit the wall
  ! this approach is brute force and therefore very slow
  ! FIXME :  Use geometry to find the intersection between LOS and first wall
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, ray_element_full_type, &
                                          straight, max_points_svec, output_level, SOL_ne, SOL_Te
    USE mod_ecfm_refr_utils, only : sub_remap_coords, func_in_poly, distance_to_poly
    USE constants,                 only : mass_e, e0, eps0
    implicit none
    type(plasma_params_type)          , intent(in)                                  :: plasma_params
    real(rkind), intent(in)                         :: omega
    integer(ikind), intent(in)                      :: mode
    integer(ikind), intent(out)                     :: last_N, wall_hits
    type(ray_element_full_type), dimension(:), intent(inout)      :: ray_segment !temporary ray
    logical, intent(out)                                             :: LOS_end
    logical                                                          :: plasma_prop
    integer(ikind)                                                   :: N, first_N
    real(rkind), dimension(3)                                        :: R_vec
    real(rkind)                                                      :: h, R_cur, R_last
    wall_hits = 0
    glob_omega = omega
    glob_mode = mode
    plasma_prop = .false.
    LOS_end = .false. ! If ray does not enter vessel or is immediately reflected upon plasma entry
    ray_segment(1)%sigma = 0.d0
    N = 1
    ray_segment(:)%theta = 0.d0
    ray_segment(1)%s = 0.d0
    first_N = 1 ! reset the first point of interpolation, because we are starting a new ray
    h = 1.d-2
    R_last = distance_to_poly(plasma_params%vessel_poly, Sqrt(ray_segment(N)%x_vec(1)**2 + ray_segment(N)%x_vec(1)**2), ray_segment(N)%x_vec(3))
    !print*, "-----------------Ray init-------------"
    do while(.not. plasma_prop)
      !print*, ray_segment(N)%s
      if(any(ray_segment(N)%x_vec /= ray_segment(N)%x_vec)) then
        print*, "nan in x_vec"
        if(N > 1 ) print*,ray_segment(N - 1)%x_vec
        print*, "step", N
        stop "Nan in x_vec in mod_raytrace.f90"
      end if
      if(N + 1 >= max_points_svec) then
        print*,"WARNING ! Maximum step count exceeded!"
        call sub_remap_coords(ray_segment(1)%x_vec, R_vec)
        print*, "R init", R_vec
        call sub_remap_coords(ray_segment(N)%x_vec, R_vec)
        print*, "R final", R_vec
        print*, "travelled distance", ray_segment(N)%s
        if(all(ray_segment(:)%rhop == -1)) print*, "No plasma along the ray - check viewing geometry!"
        last_N = N + 1
        !stop "Increase MAXIT in sub_make_ray in mod_rayrace.f90"
        return
      end if
      call sub_remap_coords(ray_segment(N)%x_vec, R_vec)
      if(wall_hits == 1 .and. .not. func_in_poly(plasma_params%vessel_poly, R_vec(1), R_vec(3))) then
         if(debug_level > 0  .and. output_level) print*, "Passed through port out of the plasma while straight line"
        print*, "Something wrong with rhop spline in mod_ecfm_refr_raytrace.f90"
        print*, ray_segment(1:N)%rhop
        print*, "Critical error during ray tracing - no plasma on LOS"
        call abort
      end if
      if (func_in_poly(plasma_params%vessel_poly, R_vec(1), R_vec(3))) then
        if(debug_level > 0 .and. output_level) print*, "Passed through port into the plasma while straight line"
        !print*, "first point in vessel", R_vec
        !plasma_prop = .true.
        wall_hits = 1
      else if(mod(N,500) == 0) then
        R_cur = distance_to_poly(plasma_params%vessel_poly, Sqrt(ray_segment(N)%x_vec(1)**2 + ray_segment(N)%x_vec(1)**2), ray_segment(N)%x_vec(3))
        if(R_cur > R_last) then
          LOS_end = .true.
          exit
        end if
        R_last = R_cur
      end if
      !print*, "plasma_prop", plasma_prop
      call sub_remap_coords(ray_segment(N)%x_vec, ray_segment(N)%R_vec)
      if((ray_segment(N)%R_vec(1)  > plasma_params%R_min .and.  ray_segment(N)%R_vec(1)  < plasma_params%R_max) .and.  &
         (ray_segment(N)%R_vec(3)  > plasma_params%z_min .and.  ray_segment(N)%R_vec(3)  < plasma_params%z_max)) then
         call sub_local_params(plasma_params, omega, ray_segment(N)%x_vec, ray_segment(N)%N_vec, ray_segment(N)%B_vec, &
          ray_segment(N)%N_s, ray_segment(N)%n_e, ray_segment(N)%omega_c,  ray_segment(N)%T_e, ray_segment(N)%theta, ray_segment(N)%rhop)
      else
          ray_segment(N)%omega_c = 0.d0
          ray_segment(N)%rhop = -1.d0
          ray_segment(N)%n_e = SOL_ne
          ray_segment(N)%T_e = SOL_Te
          ray_segment(N)%N_s = 1.d0
      end if
      ray_segment(N)%Hamil = 0.d0
      ray_segment(N + 1)%x_vec = ray_segment(N)%x_vec + h * ray_segment(N)%N_vec ! Vacuum (straight) propagation when we have no information on the plasma
      ray_segment(N + 1)%s = ray_segment(N)%s + sqrt((ray_segment(N + 1)%x_vec(1) - ray_segment(N)%x_vec(1))**2 + &
                             (ray_segment(N + 1)%x_vec(2) - ray_segment(N)%x_vec(2))**2 + &
                             (ray_segment(N + 1)%x_vec(3) - ray_segment(N)%x_vec(3))**2)
      ray_segment(N + 1)%N_vec = ray_segment(N)%N_vec
      !if(mod(N - 1, 50) == 0) print*, ray_segment(N)%N_s
      if (ray_segment(N)%rhop /= -1) then
        if((ray_segment(N)%rhop < plasma_params%rhop_entry .or. &
          plasma_params%X_entry < func_X(plasma_params, omega, ray_segment(N)%n_e, ray_segment(N)%T_e)) .and. &
          ray_segment(N)%rhop < plasma_params%rhop_max) then
          plasma_prop = .true.
        end if
      end if
      N = N + 1
    end do
    call sub_remap_coords(ray_segment(N)%x_vec, ray_segment(N)%R_vec)
    if(ray_segment(N)%R_vec(1)  > plasma_params%R_min .and.  ray_segment(N)%R_vec(1)  < plasma_params%R_max .and.  &
       ray_segment(N)%R_vec(3)  > plasma_params%z_min .and.  ray_segment(N)%R_vec(3)  < plasma_params%z_max) then
        ray_segment(N)%omega_c = func_B_abs(plasma_params, ray_segment(N)%x_vec) * e0 / mass_e
    else
        ray_segment(N)%omega_c = 0.d0
    end if
    if(.not. ((ray_segment(N)%R_vec(1)  > plasma_params%R_min .and.  ray_segment(N)%R_vec(1)  < plasma_params%R_max) .and.  &
         (ray_segment(N)%R_vec(3)  > plasma_params%z_min .and.  ray_segment(N)%R_vec(3)  < plasma_params%z_max))) then
          print*, "first point in plasma has no equilibrium"
          print*, "R, z",ray_segment(N)%R_vec(1), ray_segment(N)%R_vec(3)
          print*, "R_min, R_max, z_min, z_max given by eq:", plasma_params%R_min, plasma_params%R_max,plasma_params%z_min,plasma_params%z_max
          stop "no equilibrium for within vessel in mod_raytrace.f90"
    end if
    call sub_local_params(plasma_params, omega, ray_segment(N)%x_vec, ray_segment(N)%N_vec, ray_segment(N)%B_vec, &
          ray_segment(N)%N_s, ray_segment(N)%n_e, ray_segment(N)%omega_c,  ray_segment(N)%T_e, ray_segment(N)%theta, ray_segment(N)%rhop)
    last_N = N
    !print*,"Vacuum propagation from", ray_segment(1)%x_vec, "to", ray_segment(N)%x_vec
    if(straight .or. LOS_end) then
      ray_segment(N)%Hamil = 0.d0
      return
    end if
    call sub_calculate_initial_N(plasma_params, omega, mode, ray_segment(N)%x_vec, ray_segment(N)%N_vec, ray_segment(N)%Hamil, LOS_end)
  end subroutine find_first_point_in_plasma

  subroutine make_ray_segment(distance, plasma_params, omega, mode, ray_segment, last_N, wall_hits, N_start)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, ray_element_full_type, Hamil, &
                                    straight, max_points_svec, output_level, UH_stop
    USE mod_ecfm_refr_utils, only : sub_remap_coords, func_in_poly
    USE constants,                 only : pi, mass_e, e0, eps0
    implicit none
    type(plasma_params_type)          , intent(in)                                  :: plasma_params
    real(rkind), intent(in)                         :: omega, distance
    integer(ikind), intent(in)                      :: mode
    integer(ikind), intent(out)                     :: last_N
    type(ray_element_full_type), dimension(:), intent(inout)      :: ray_segment !temporary ray
    integer(ikind), intent(inout)                   :: wall_hits
    integer(ikind), intent(in), optional            :: N_start
    logical                                                          :: propagating
    integer(ikind)                                                   :: N, i, istate, first_N
    real(rkind), dimension(3)                                        :: R_vec
    real(rkind)                                                      :: first_s, X, Y, angle_change
    real(rkind), dimension(116)               :: work_lsode
    integer(ikind), dimension(20)             :: iwork_lsode
    iwork_lsode(:) = 0.d0
    iwork_lsode(6) = 8000
    work_lsode(:) = 0.d0
    propagating = .true.
    ray_segment(1)%sigma = 0.d0
    if(present(N_start)) then
      N = N_start
      first_N = N_start
    else
      N = 1
      first_N = 1
    end if
    first_s = ray_segment(first_N)%s
    if(N == 1) then
      ray_segment(N)%h = plasma_params%h
    else
      ray_segment(N - 1)%h = plasma_params%h ! for plasma propagation we want to start with a small h
    end if
    !call sub_remap_coords(ray_segment(N)%x_vec, ray_segment(N)%R_vec)
    last_N = N
    ray_segment(:)%theta = 0.d0
    istate = 1
    !print*, "Initial wall hits", wall_hits
!    open(96, file= "k_out_ecfm")
    !print*, "-----------------Ray trace init-------------"
    do while(propagating)
      !print*, ray_segment(N)%s
      if(N + 1 > max_points_svec) then
        print*,"WARNING ! Maximum step count exceeded!"
        call sub_remap_coords(ray_segment(first_N)%x_vec, R_vec)
        print*, "R init", R_vec
        call sub_remap_coords(ray_segment(N - 1)%x_vec, R_vec)
        print*, "R final", R_vec
        print*, "travelled distance", ray_segment(N -1)%s
        print*, "last step size and step count", ray_segment(N - 1)%h , N
        print*, "arc length of the last 2000 steps", ray_segment(N - 1)%s - ray_segment(N - 1000)%s , N
        !stop "Increase MAXIT in sub_make_ray in mod_rayrace.f90"
        return
      end if
      if(any(ray_segment(N)%x_vec /= ray_segment(N)%x_vec)) then
        print*, "nan in x_vec"
        if(N > 1 ) print*,ray_segment(N - 1)%x_vec
        print*, "step", N
        stop "Nan in x_vec in mod_raytrace.f90"
      end if
      !print*, "plasma_prop", plasma_prop
      if(.not. straight) then
         ! also retrieve N_par, N_s, X, Y for current step
        if(Hamil == "Dani") then
          if(N > 1) ray_segment(N)%h = ray_segment(N - 1)%h ! retrieve updated h from last step
          ray_segment(N + 1)%s = ray_segment(N)%s
          work_lsode(1) = max(ray_segment(N + 1)%s + 4.d0 * ray_segment(N)%h, work_lsode(13) + 4.d0 * ray_segment(N)%h)! s can be smaller than s from the solver
          call sub_single_step_LSODE(ray_segment(N + 1)%s, ray_segment(N)%x_vec, ray_segment(N)%N_vec, &
                               ray_segment(N)%h, ray_segment(N + 1)%x_vec, &
                               ray_segment(N + 1)%N_vec, ray_segment(N + 1)%B_vec, ray_segment(N + 1)%theta, &
                               ray_segment(N + 1)%Hamil, ray_segment(N + 1)%N_s, ray_segment(N + 1)%n_e, &
                               ray_segment(N + 1)%omega_c,  ray_segment(N + 1)%T_e, ray_segment(N + 1)%rhop, &
                               ray_segment(N + 1)%v_g_perp, work_lsode, iwork_lsode, istate)
          call sub_remap_coords(ray_segment(N + 1)%x_vec, ray_segment(N + 1)%R_vec)
            !print*, "x, R", ray_segment(N + 1)%x_vec, ray_segment(N + 1)%R_vec
        else
          ray_segment(N + 1)%sigma = ray_segment(N)%sigma
          work_lsode(1) = ray_segment(N + 1)%sigma  + 4.d0 * ray_segment(N)%h
          call sub_single_step_LSODE(ray_segment(N + 1)%sigma, ray_segment(N)%x_vec, ray_segment(N)%N_vec, &
                               ray_segment(N)%h, ray_segment(N + 1)%x_vec, ray_segment(N + 1)%N_vec, &
                               ray_segment(N + 1)%B_vec, ray_segment(N + 1)%theta, ray_segment(N + 1)%Hamil, ray_segment(N + 1)%N_s, &
                               ray_segment(N + 1)%n_e, ray_segment(N + 1)%omega_c, ray_segment(N + 1)%T_e, &
                               ray_segment(N + 1)%rhop, ray_segment(N + 1)%v_g_perp, work_lsode, iwork_lsode, istate)
          !ray_segment(N)%N_s = sqrt(ray_segment(N)%N_vec(1)**2 + ray_segment(N)%N_vec(2)**2 + ray_segment(N)%N_vec(3)**2)
          call sub_remap_coords(ray_segment(N + 1)%x_vec, ray_segment(N + 1)%R_vec)
          !print*, ray_segment(N + 1)%R_vec, ray_segment(N)%rhop
          ray_segment(N + 1)%s = ray_segment(N)%s + sqrt((ray_segment(N + 1)%x_vec(1) - ray_segment(N)%x_vec(1))**2 + (ray_segment(N + 1)%x_vec(2) - ray_segment(N)%x_vec(2))**2 + &
                (ray_segment(N + 1)%x_vec(3) - ray_segment(N)%x_vec(3))**2)
        end if
        !if(ray_segment(N + 1)%rhop < 1.2) print*, "rhop, B, n_e", ray_segment(N + 1)%rhop, ray_segment(N + 1)%omega_c / e0 * mass_e, ray_segment(N + 1)%n_e
        X = func_X(plasma_params,omega, ray_segment(N + 1)%n_e, ray_segment(N + 1)%T_e)
        Y = func_Y(plasma_params,omega, ray_segment(N + 1)%omega_c / e0 * mass_e, ray_segment(N + 1)%T_e)
        !N_cold = func_N(X, Y, &
        !  ray_segment(N + 1)%theta,- mode)
!        write(96, fmt="(10E15.7)") ray_segment(N + 1)%x_vec, ray_segment(N + 1)%N_vec,ray_segment(N + 1)%rhop,  X, Y, ray_segment(N + 1)%theta
        if(istate /= 2) then
          propagating = .false.
          wall_hits = 2
          print*, "stopped propagating because of LSODE error"
        end if
      else
        ray_segment(N + 1)%x_vec = ray_segment(N)%x_vec + plasma_params%h * ray_segment(N)%N_vec ! Vacuum (straight) propagation
        ray_segment(N + 1)%s = ray_segment(N)%s + sqrt((ray_segment(N + 1)%x_vec(1) - ray_segment(N)%x_vec(1))**2 + (ray_segment(N + 1)%x_vec(2) - ray_segment(N)%x_vec(2))**2 + &
                (ray_segment(N + 1)%x_vec(3) - ray_segment(N)%x_vec(3))**2)
        ray_segment(N + 1)%N_vec = ray_segment(N)%N_vec
        call sub_local_params(plasma_params, omega, ray_segment(N + 1)%x_vec, ray_segment(N + 1)%N_vec, ray_segment(N + 1)%B_vec, &
                              ray_segment(N + 1)%N_s, ray_segment(N + 1)%n_e, ray_segment(N + 1)%omega_c,  &
                              ray_segment(N + 1)%T_e, ray_segment(N + 1)%theta, ray_segment(N + 1)%rhop)
        call sub_remap_coords(ray_segment(N + 1)%x_vec, ray_segment(N + 1)%R_vec)
        ray_segment(N + 1)%Hamil = 0.d0
      end if
      if(ray_segment(N + 1)%s - first_s > distance ) then
        ! .or. (ray_segment(N + 1)%Y < 0.51 .and. ray_segment(N + 1)%Y > 0.49 .and. .not. fine)
        if(N <= 1) print*, "stopped propagating because of target distance reached: ", ray_segment(N)%s - first_s, distance
        propagating = .false.
      end if
      if(UH_stop  .and. .not. straight .and. ((mode > 0 .and. Hamil == "Dani") .or. (mode < 0 .and. Hamil == "Stix"))) then
        if(abs(1.d0 - sqrt(X + Y**2)) < 0.01) then
          if(output_level) print*, "Close to UH-resonance - finishing ray"
          wall_hits =  2
          propagating = .false.
        end if
      end if
      if(N > 3 .and. .not. straight) then
          angle_change = 0.d0
          do i = 1, 3
            angle_change = angle_change + ray_segment(N + 1)%N_vec(i) * ray_segment(first_N)%N_vec(i)
          end do
          angle_change = angle_change / sqrt(ray_segment(N + 1)%N_vec(1)**2 + ray_segment(N + 1)%N_vec(2)**2 + ray_segment(N + 1)%N_vec(3)**2)
          angle_change = angle_change / sqrt(ray_segment(1)%N_vec(1)**2 + ray_segment(1)%N_vec(2)**2 + ray_segment(1)%N_vec(3)**2)
          if(abs(acos(angle_change)) > plasma_params%angle_threshold) then
              wall_hits =  wall_hits + 1
              if(debug_level > 0 .and. output_level) print*, "Propagation was rotated by", acos(angle_change) / pi * 180.d0, " deg - finishing ray"
              !print*, ray_segment(N - 2)%omega_c * mass_e / e0,ray_segment(N - 1)%omega_c * mass_e / e0,ray_segment(N)%omega_c * mass_e / e0
          end if
          if(wall_hits > 1) propagating = .false.
          ! TODO: For HFS ECE diagnostics the statement above needs to be inverted
          !if(debug_level > 0) print*, acos(angle_change) / pi * 180.d0
      end if
      N = N + 1
      call sub_remap_coords(ray_segment(N)%x_vec, R_vec)
      if(R_vec(1) - plasma_params%h < plasma_Params%R_min .or. &
         R_vec(1) + plasma_params%h > plasma_Params%R_max .or. &
         R_vec(3) - plasma_params%h < plasma_Params%z_min .or. &
         R_vec(3) + plasma_params%h > plasma_Params%z_max) then
         wall_hits =  2
         if(debug_level > 0 .and. output_level) print*, "Left the domain on which the flux matrix is given"
         if(debug_level > 0 .and. output_level) print*, "Position",  R_vec(1), R_vec(3)
         propagating  = .false.
      else if(ray_segment(N)%rhop == -1.d0) then
        if(debug_level > 0 .and. output_level) print*, "Rhop not useful anymore stopping propagation"
        if(debug_level > 0 .and. output_level) print*, "Position",  R_vec(1), R_vec(3)
        wall_hits =  2
        propagating  = .false.
      else if (any(abs(ray_segment(1:N)%rhop) < plasma_params%rhop_inside) .and. ray_segment(N)%rhop > plasma_params%rhop_exit) then
        if(debug_level > 0 .and. output_level) print*, "Rhop now larger than rhop_exit after pass through plasma"
        if(debug_level > 0 .and. output_level) print*, "Position",  R_vec(1), R_vec(3)
        wall_hits =  2
        propagating  = .false.
      else if(.not. func_in_poly(plasma_params%vessel_poly, R_vec(1), R_vec(3))) then
        if(wall_hits > 0) then !left machine
          wall_hits =  wall_hits + 1
          if(debug_level > 0 .and. output_level) print*, "Passed through port out of the plasma"
          if(debug_level > 0 .and. output_level) print*, "Position",  R_vec(1), R_vec(3)
        end if
      else if(wall_hits == 0) then
        wall_hits =  wall_hits + 1
        if(debug_level > 0 .and. output_level) print*, "Passed through port into plasma"
        if(debug_level > 0 .and. output_level) print*, "Position",  R_vec(1), R_vec(3)
      end if
      if(wall_hits > 1) propagating = .false.
    end do
    last_N = N - 1
    !print*,"Plasma propagation from", ray_segment(first_N)%x_vec, "to", ray_segment(last_N)%x_vec
    if(last_N <= N_start + 2) then
      print*, "Raytracing stopped after just 2 iterations"
      print*, "place",ray_segment(N - 1)%x_vec
      print*, "trajectory",ray_segment(N - 1)%N_vec
      stop "Error in raytracing in mod_refr_raytracing.f90"
    end if
!    close(96)
  end subroutine make_ray_segment

  subroutine make_s_grid(plasma_params, omega, Y_res, svec, ray, total_LOS_points, &
                         last_N, dist, rad_ray_freq, max_points_svec_reached)
  ! Linearly interpolates R,z on the ray which corresponds
  ! to the s values given oprintn an equidistant grid.
  ! The cold resonance is detected automatically and the step size is choosen correspondingly.
  ! Once the R,z values are obtained the other ray parameters are aqcuired using sub_local_params.
  ! Every time this routine is called the svec values from i_start to i_end will be filled.
  ! i referst to the svec, while N refers to the ray_segment
  ! This routine also conrolls the step size.
  use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, plasma_params_type, &
                                        ray_element_full_type, max_points_svec, &
                                        spl_type_1d, eps_svec_max_length, rad_diag_ch_mode_ray_freq_type
  use mod_ecfm_refr_utils,        only: binary_search
  use f90_kind
  use mod_ecfm_refr_interpol,     only: make_1d_spline, deallocate_1d_spline, spline_1d_get_roots
  implicit none
  type(plasma_params_type), intent(in)                                       :: plasma_params
  real(rkind), intent(in)                                                    :: omega
  real(rkind), dimension(:), intent(in)                                      :: Y_res! array holding all resonances that require special treatment
  type(rad_diag_ch_mode_ray_freq_svec_type), dimension(:), intent(inout)     :: svec
  type(ray_element_full_type), dimension(max_points_svec), intent(in)        :: ray !the ray
  integer(ikind), intent(out)                                                :: total_LOS_points
  integer(ikind),   intent(in)                                               :: last_N
  real(rkind), dimension(:), intent(in)                                      :: dist ! distance to interpolated, starting point of interpolation
  type(rad_diag_ch_mode_ray_freq_type), intent(inout), optional              :: rad_ray_freq
  logical, intent(out), optional                                             :: max_points_svec_reached
  type(spl_type_1d)                                                          :: spl
  real(rkind), dimension(last_N)                                             :: flush_ray_s, flush_ray_y, roots, s_dense, s_dense_debug
  integer(ikind)                                                             :: i, N_roots, N_s_dense, i_root, roots_processed, grid, &
                                                                                i_next, i_next_root, i_svec, i_last, N_interval
  real(rkind)                                                                :: s_next, a, b, i_interval, cur_dist
  logical                                                                    :: ray_finished
  flush_ray_s = ray(1:last_N)%s
  flush_ray_y = ray(1:last_N)%omega_c / omega
  N_s_dense = 1
  do i = 1,size(Y_res)
  ! Upshift -> resonance at magnetic field slightly lower than at cold resonance
    call make_1d_spline(spl, last_N, flush_ray_s, flush_ray_y - (Y_res(i) + plasma_params%up_shift))
    call spline_1d_get_roots(spl, roots, N_roots)
    do i_root= 1, N_roots
      if(roots(i_root) > 0.d0) then
        s_dense(N_s_dense) = roots(i_root)
        N_s_dense = N_s_dense + 1
      end if
    end do
  ! Downshift -> resonance at magnetic field slightly higher than at cold resonance
    call make_1d_spline(spl, last_N, flush_ray_s, flush_ray_y - (Y_res(i) + plasma_params%down_shift))
    call spline_1d_get_roots(spl, roots, N_roots)
    do i_root= 1, N_roots
      if(roots(i_root) > 0.d0) then
        s_dense(N_s_dense) = roots(i_root)
        N_s_dense = N_s_dense + 1
      end if
    end do
  end do
  if(present(rad_ray_freq)) then
    call make_1d_spline(spl, last_N, flush_ray_s, flush_ray_y - plasma_params%Y_res)
    call spline_1d_get_roots(spl, roots, N_roots)
    if(N_roots == 0) then
      rad_ray_freq%s_res = -1.0
    else
      rad_ray_freq%s_res = minval(roots(1:N_roots))! Convert from ray coordintes to svec coordinates
    end if
  end if
  N_s_dense = N_s_dense - 1 ! Fortran stuff
  s_dense_debug(1:N_s_dense) = s_dense(1:N_s_dense)
  call deallocate_1d_spline(spl)
  grid = 2 ! Grid size of -1 step, usually first step is large that means the previous is small
  !Check if first point already in dense
  do i = 1,size(Y_res)
    if(ray(1)%omega_c / omega < Y_res(i) + plasma_params%down_shift .and. &
       ray(1)%omega_c / omega > Y_res(i) + plasma_params%up_shift) grid = 1 !-> next grid size needs to be large hence intial grid size (-1) grid size is large
  end do
  ray_finished = .false.
  i_next = 1 ! -> for ray
  i_last = 0 ! -> only for diagnostic purposes
  i_svec = 1 ! for svec
  a = 0.d0
  roots_processed = 0
  s_next = -1.d0 ! For initialization
  ! Span the grid
  do while(.not. ray_finished)
    i_last = i_next
    if(roots_processed < N_s_dense) then
    ! Process one interval
      i_next_root = minloc(s_dense(1:N_s_dense), dim=1)
      s_next = s_dense(i_next_root)
      s_dense(i_next_root) = 1.d99 !-> already processed -> Replacement of python pop function
      roots_processed = roots_processed + 1
      if(s_next == 1.d99) then
        print*, "Bug in make_s_grid ran out of s_dense"
        call abort()
      end if
      if(grid == 1) then ! Determines grid size for next step
          grid = 2
      else
          grid = 1
      end if
      do while(s_next <= a .and. roots_processed < N_s_dense)
        ! If the last step caused us to skip a segment we need to find another segment
        ! This only happens in case of the small step size so it does not affect numerical stability negatively
        i_next_root = minloc(s_dense(1:N_s_dense), dim=1)
        s_next = s_dense(i_next_root)
        s_dense(i_next_root) = 1.d99 !-> already processed -> Replacement of python pop function
        roots_processed = roots_processed + 1
        if(s_next == 1.d99) then
          print*, "Bug in make_s_grid ran out of s_dense"
          call abort()
        end if
        ! Still need to swap grids since we need to keep up the sequence of dense and sparse grids
        if(grid == 1) then ! Swap grid size again since we skipped a step
            grid = 2
        else
            grid = 1
        end if
      end do
      cur_dist = dist(grid) ! No more swaps between grid sizes from here on out
      ! If we needed to skip ahead we might have already finished all intervals by now -> check for that here
      if(s_next <= a .and. roots_processed >= N_s_dense) then
        i_next = last_N
        N_interval = floor( (flush_ray_s(i_next) - a)/ cur_dist) + 1
        ! To avoid interpolation errors svec is slightly shorter than the initial ray -> eps_svec_max_length
        cur_dist = (flush_ray_s(i_next) - a - eps_svec_max_length)/ real(N_interval,8) ! slightly smaller than dist(grid)
        ray_finished = .true.
      else if(grid == 1) then
        i_next = binary_search(flush_ray_s(1:last_N), s_next, i_last, last_N)
        if(i_next < i_last .or. i_next > last_N) then
          print*, "Something went wrong when spanning the svec grid"
          print*, "s_1, s, s_2", flush_ray_s(i_last), s_next, flush_ray_s(last_N)
          call abort()
        end if
        N_interval = floor( (flush_ray_s(i_next) - flush_ray_s(i_last))/ cur_dist)
      else ! prefer small over large grid
        i_next = binary_search(flush_ray_s(1:last_N), s_next, i_last, last_N) + 1 ! One extra cell here
        N_interval = floor( (flush_ray_s(i_next) - a)/ cur_dist) + 1 ! Also overlap here
        if(a + N_interval * cur_dist >= flush_ray_s(last_N)) then
            cur_dist = (flush_ray_s(i_next) - a - eps_svec_max_length)/ real(N_interval,8) ! slightly smaller than dist(grid)
            ray_finished = .true.
        end if
      end if
    else
    ! Swap one last time
      if(grid == 1) then ! Swap grid size again since we skipped a step
          grid = 2
      else
          grid = 1
      end if
      cur_dist = dist(grid)
      i_next = last_N
      N_interval = floor( (flush_ray_s(i_next) - a)/ cur_dist) + 1
      ! To avoid interpolation errors svec is slightly shorter than the initial ray -> eps_svec_max_length
      cur_dist = (flush_ray_s(i_next) - a - eps_svec_max_length)/ real(N_interval,8) ! slightly smaller than dist(grid)
      ray_finished = .true.
      !print*, "All roots processed moving to end of ray"
    end if
    do i_interval = 1, N_interval
    ! Process the current interval
      if(i_svec + plasma_params%int_step_cnt > max_points_svec) then
        if(present(max_points_svec_reached)) then
          max_points_svec_reached = .true.
          return
        else
          print*, "max_points_svec has been reached, increase max_point_svec in the input file"
          print*,"===================== DEBUG INFO ========================================"
          print*, "First and last point of the ray in [m]", flush_ray_s(1), flush_ray_s(last_N)
          print*, "First and last point of svec", svec(1)%s, svec(i_svec - 1)%s
          print*, "Current position in the ray and end of current segment: ", a, cur_dist + a
          print*, "Current step size and dist", cur_dist, dist
          print*, "i_svec and N_interval", i_svec
          print*, "Dense points", s_dense_debug(1:N_s_dense)
          open(69, file="s_grid.dat")
          do i = 1, i_svec
            write(69, "(E18.10E3)") svec(i)%s
          end do
          close(69)
          print*, "Created s_grid.dat containing the current s grid"
          call abort()
        end if
      end if
      b = a + cur_dist    
      svec(i_svec :i_svec + plasma_params%int_step_cnt - 1)%s = cur_dist * plasma_params%Int_absz(:) + a
      i_svec = i_svec + plasma_params%int_step_cnt
      a = b      
    end do
  end do
  total_LOS_points = i_svec - 1
  svec(1:total_LOS_points)%s = svec(1:total_LOS_points)%s - svec(1)%s ! svec(1) != due to the way that plasma_params%Int_absz is set up
  end subroutine make_s_grid

  subroutine interpolate_svec(plasma_params, svec, ray, omega, total_LOS_points, N, rad_ray_freq, svec_extra_output)
  use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, plasma_params_type, &
                                        ray_element_full_type, max_points_svec, &
                                        SOL_ne, SOL_Te, spl_type_1d, rad_diag_ch_mode_ray_freq_type, &
                                        rad_diag_ch_mode_ray_freq_svec_extra_output_type
  use mod_ecfm_refr_interpol,     only: make_1d_spline, deallocate_1d_spline, spline_1d
  use constants,                  only: pi, e0, mass_e, eps0, c0
  use mod_ecfm_refr_utils,        only: sub_remap_coords, binary_search
  use f90_kind
  implicit none
  type(plasma_params_type), intent(in)                                       :: plasma_params
  type(rad_diag_ch_mode_ray_freq_svec_type), dimension(:), intent(inout)     :: svec
  type(ray_element_full_type), dimension(max_points_svec), intent(in)        :: ray !temporary ray
  real(rkind), intent(in)                                                    :: omega
  integer(ikind),   intent(in)                                               :: N, total_LOS_points
  type(rad_diag_ch_mode_ray_freq_type), intent(inout), optional              :: rad_ray_freq
  type(rad_diag_ch_mode_ray_freq_svec_extra_output_type), dimension(:), intent(inout), optional :: svec_extra_output
  type(spl_type_1d)                                                          :: spl
  real(rkind), dimension(N)                                                  :: flush_ray_s, flush_ray_y
  real(rkind), dimension(total_LOS_points)                                   :: flush_svec_s, flush_svec_y
  integer(ikind)                                                             :: i, N_plasma, i1, i2
  real(rkind), dimension(3)                                                  :: x_res!, B_res
  ! Use arrays that are continuous in the memory
  ! The ray is expected to already reverted from the ray tracing such that the propagation directions is TOWARDS the antenna
  flush_ray_s = ray(1:N)%s
  flush_svec_s = svec(1:total_LOS_points)%s
!  print*, "N", total_LOS_points, N
  ! First quantities that are present on the entire ray
  do i=1,3
    flush_ray_y = ray(1:N)%x_vec(i)
    call make_1d_spline(spl, N, flush_ray_s, flush_ray_y)
    call spline_1d(spl, flush_svec_s, flush_svec_y)
    svec(1:total_LOS_points)%x_vec(i) = flush_svec_y
    if(present(rad_ray_freq)) then
      if(rad_ray_freq%s_res >= 0.d0) call spline_1d(spl, rad_ray_freq%s_res, x_res(i))
    end if
    flush_ray_y = -ray(1:N)%N_vec(i)
    call make_1d_spline(spl, N, flush_ray_s, flush_ray_y)
    call spline_1d(spl, flush_svec_s, flush_svec_y)
    svec(1:total_LOS_points)%N_vec(i) = flush_svec_y
  end do
  svec(1:total_LOS_points)%R = sqrt(svec(1:total_LOS_points)%x_vec(1)**2 + svec(1:total_LOS_points)%x_vec(2)**2)
  svec(1:total_LOS_points)%z = svec(1:total_LOS_points)%x_vec(3)
  if(present(rad_ray_freq)) then
    if(rad_ray_freq%s_res < 0.d0) then
      rad_ray_freq%R_res = -1.d0
      rad_ray_freq%z_res = -1.d0
    else
      rad_ray_freq%R_res = sqrt(x_res(1)**2 + x_res(2)**2)
      rad_ray_freq%z_res = x_res(3)
    end if
  end if
  N_plasma = size(pack(flush_ray_s, ray(1:N)%rhop > 0))
  flush_ray_s(1:N_plasma) = pack(flush_ray_s, ray(1:N)%rhop > 0)
  i1 = minloc(abs(flush_ray_s(1) - flush_svec_s), dim=1)
  if(flush_svec_s(i1) <  flush_ray_s(1)) i1 = i1 + 1
  i2 = minloc(abs(flush_ray_s(N_plasma) - flush_svec_s), dim=1)
  if(flush_svec_s(i2) >  flush_ray_s(N_plasma)) i2 = i2 - 1
  !Initialize everything for scrape of layer values
  ! Later overwrite with plasma parameters
  do i=1,3
    svec(1:total_LOS_points)%B_vec(i) = 0.d0 ! Could be replaced with vacuum values but seems unneccessary
  end do
  svec(1:total_LOS_points)%rhop = -1.d0
  svec(1:total_LOS_points)%ne = SOL_ne
  svec(1:total_LOS_points)%Te = SOL_Te
  svec(1:total_LOS_points)%freq_2X = 0.d0
  svec(1:total_LOS_points)%theta = 0.d0
  svec(1:total_LOS_points)%cos_theta = 1.d0
  svec(1:total_LOS_points)%sin_theta = 0.d0
  svec(1:total_LOS_points)%N_cold = 1.d0
  svec(1:total_LOS_points)%v_g_perp = 0.d0
  svec(1:total_LOS_points)%plasma = .false.
  if(present(svec_extra_output)) then
    svec_extra_output(1:total_LOS_points)%H = 0.d0
    svec_extra_output(1:total_LOS_points)%N_ray = 1.d0
  end if
  flush_svec_s(1:i2 - i1 + 1) = flush_svec_s(i1:i2)
  svec(i1:i2)%plasma = .true.
  !if(present(rad_ray_freq)) B_res(:) = 0.d0
  do i=1,3
    flush_ray_y(1:N_plasma) = pack(ray(1:N)%B_vec(i), ray(1:N)%rhop > 0)
    call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
    call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
    svec(i1:i2)%B_vec(i) = flush_svec_y(1:i2 - i1 + 1)
!    if(present(rad_ray_freq) .and. rad_ray_freq%s_res >  flush_ray_s(1) .and. &
!                                   rad_ray_freq%s_res <  flush_ray_s(N_plasma)) then
!      call spline_1d(spl, rad_ray_freq%s_res, B_res(i))
!    end if
  end do
  flush_ray_y(1:N_plasma) = pack(ray(1:N)%rhop, ray(1:N)%rhop > 0)
  call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
  call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
  svec(i1:i2)%rhop = flush_svec_y(1:i2 - i1 + 1)
  if(present(rad_ray_freq)) then
    if(rad_ray_freq%s_res < 0.d0 .or. rad_ray_freq%s_res > maxval(flush_ray_s(1:N_plasma))) then
      rad_ray_freq%rhop_res = -1.d0
    else
      call spline_1d(spl, rad_ray_freq%s_res, rad_ray_freq%rhop_res)
    end if
  end if
  flush_ray_y(1:N_plasma) = pack(ray(1:N)%n_e, ray(1:N)%rhop > 0)
  if(plasma_params%prof_log_flag) then
    flush_ray_y(1:N_plasma) = log(flush_ray_y(1:N_plasma) * 1.e-19)
  end if
  call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
  call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
  svec(i1:i2)%ne = flush_svec_y(1:i2 - i1 + 1)
  if(plasma_params%prof_log_flag) then
    svec(i1:i2)%ne = exp(svec(i1:i2)%ne) * 1.e19
  end if
  flush_ray_y(1:N_plasma) = pack(ray(1:N)%T_e, ray(1:N)%rhop > 0)
  if(plasma_params%prof_log_flag) then
    flush_ray_y(1:N_plasma) = log(flush_ray_y(1:N_plasma))
  end if
  call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
  call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
  svec(i1:i2)%Te = flush_svec_y(1:i2 - i1 + 1)
  if(plasma_params%prof_log_flag) then
    svec(i1:i2)%Te = exp(svec(i1:i2)%Te)
  end if
  flush_ray_y(1:N_plasma) = pack(ray(1:N)%N_s, ray(1:N)%rhop > 0)
  call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
  call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
  svec(i1:i2)%N_cold = flush_svec_y(1:i2 - i1 + 1)
  flush_ray_y(1:N_plasma) = pack(ray(1:N)%v_g_perp, ray(1:N)%rhop > 0)
  call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
  call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
  svec(i1:i2)%v_g_perp = flush_svec_y(1:i2 - i1 + 1)
  if(present(svec_extra_output)) then
    flush_ray_y(1:N_plasma) = pack(ray(1:N)%hamil, ray(1:N)%rhop > 0)
    call make_1d_spline(spl, N_plasma, flush_ray_s(1:N_plasma), flush_ray_y(1:N_plasma))
    call spline_1d(spl, flush_svec_s(1:i2 - i1 + 1), flush_svec_y(1:i2 - i1 + 1))
    svec_extra_output(i1:i2)%H = flush_svec_y(1:i2 - i1 + 1)
    svec_extra_output(i1:i2)%N_ray = 0.d0
    do i=1,3
      svec_extra_output(i1:i2)%N_ray = svec_extra_output(i1:i2)%N_ray + svec(i1:i2)%N_vec(i)**2
    end do
    svec_extra_output(i1:i2)%N_ray = sqrt(svec_extra_output(i1:i2)%N_ray)
  end if
  svec(i1:i2)%freq_2X = e0 * sqrt(svec(i1:i2)%B_vec(1)**2 + svec(i1:i2)%B_vec(2)**2 + &
                                  svec(i1:i2)%B_vec(3)**2) / (mass_e * pi)
!  if(present(rad_ray_freq)) then
!    if(rad_ray_freq%s_res >= 0.d0)  print*, "n at resonance", omega / (e0 * sqrt(B_res(1)**2 + B_res(2)**2 + &
!                                           B_res(3)**2) / mass_e)
!  end if
  svec(i1:i2)%cos_theta = 0.d0
  do i=1,3
    svec(i1:i2)%cos_theta = svec(i1:i2)%cos_theta + svec(i1:i2)%N_vec(i) * svec(i1:i2)%B_vec(i)
  end do
  svec(i1:i2)%cos_theta = svec(i1:i2)%cos_theta / sqrt(svec(i1:i2)%B_vec(1)**2 + svec(i1:i2)%B_vec(2)**2 + &
                                               svec(i1:i2)%B_vec(3)**2)
  svec(i1:i2)%cos_theta = svec(i1:i2)%cos_theta / sqrt(svec(i1:i2)%N_vec(1)**2 + svec(i1:i2)%N_vec(2)**2 + &
                                               svec(i1:i2)%N_vec(3)**2)
  svec(i1:i2)%theta = acos(svec(i1:i2)%cos_theta)
  svec(i1:i2)%sin_theta = sin(svec(i1:i2)%theta)
  ! Do this for all points
  svec(1:total_LOS_points)%ibb = (omega / ( 2.d0 * pi))**2 * e0 * &
                svec(1:total_LOS_points)%Te / c0**2
  call deallocate_1d_spline(spl)
  end subroutine interpolate_svec

  subroutine span_svecs(plasma_params)
  ! Creates the svecs for a all diags
  ! Also allocates all vectors related to output_level = .true.
  use mod_ecfm_refr_types,        only: rad, ant, rad_diag_type, plasma_params_type, &
                                        N_ray, N_freq, modes, mode_cnt, output_level, pnts_BPD, &
                                        max_points_svec, Hamil, straight, largest_svec, ray_element_full_type
  use f90_kind
  use constants,                  only: pi,e0, mass_e
  implicit none
  type(plasma_params_type), intent(inout)                       :: plasma_params
  integer(ikind)                                                :: idiag, last_N, ich, ir, &
                                                                   ifreq, imode,  N, N_init, mode
  real(rkind), dimension(2)                                     :: dist
  type(ray_element_full_type), dimension(:), allocatable        :: ray_segment
  real(rkind)                                                   :: omega, temp, X, Y
  real(rkind), dimension(1)                                     :: Y_res_O
  real(rkind), dimension(2)                                     :: Y_res_X
  integer(ikind)                                                :: wall_hits
  logical                                                       :: LOS_end
  if(output_level) then
    if(.not. straight) then
      print*, "Preparing LOS including refraction - this will take a moment"
    else
      print*, "Preparing straight LOSs - this should take only a moment"
    end if
  end if
  ! Only X1, X2, and O1 have strong absorption near the cold resonance
  ! for all other harmonics the large step size should be sufficient
  Y_res_O(1) = 1.d0 ! O1
  Y_res_X(1) = 1.d0 ! X1
  Y_res_X(2) = 0.5d0 ! X2
  ! Far away from strongly absorbing resonance -> large steps
  dist(1) = plasma_params%dist_large * (plasma_params%int_step_cnt)
  ! Close to strongly absorbing resonance -> small steps (usually 10 times smaller)
  dist(2) = plasma_params%dist_small * (plasma_params%int_step_cnt)
  glob_plasma_params = plasma_params
  do idiag = 1, ant%N_diag
#ifdef OMP
    !$omp parallel private(ich, imode, ir, ifreq, &
    !$omp                 N_init, last_N, wall_hits, N, &
    !$omp                 omega, temp, X, Y, &
    !$omp                 ray_segment, mode) default(shared)
#endif
    allocate(ray_segment(max_points_svec), x_loc_vec(3), N_loc_vec(3))
#ifdef OMP
    !$omp do
#endif
    do ich = 1, ant%diag(idiag)%N_ch
      debug_level = plasma_params%debug_level
      do imode = 1, mode_cnt
        if(Hamil == "Dani") then
          mode = rad%diag(idiag)%ch(ich)%mode(imode)%mode
        else
          mode = -rad%diag(idiag)%ch(ich)%mode(imode)%mode
        end if
        omega = ant%diag(idiag)%ch(ich)%f_ECE * 2.d0 * pi
        rad%diag(idiag)%ch(ich)%mode(imode)%s_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%R_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%z_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res = 0.d0
        do ir = 1, N_ray
          ray_segment(:)%h = plasma_params%h
          ifreq = 1 ! Raytrace only central frequency
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res = 0.d0
          if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%use_external_pol_coeff .and. &
             rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%pol_coeff == 0.d0) then
             cycle
          end if
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%max_points_svec_reached = .false.
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes = .true.
          ray_segment(1)%x_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec ! get launching position
          ray_segment(1)%N_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%N_vec ! get launching angles
          if(abs(ray_segment(1)%x_vec(2)) < 1.d-17) then
          ! Rotate everything by 180 degrees to avoid discontinuity at phi = 0.d0
            print*, "Rotating current launch by 180.d0 degrees in phi to avoid discontinuity at phi=0"
            temp = ray_segment(1)%x_vec(1)
            ray_segment(1)%x_vec(1) = ray_segment(1)%x_vec(2)
            ray_segment(1)%x_vec(2) = temp
            temp = ray_segment(1)%N_vec(1)
            ray_segment(1)%N_vec(1) = ray_segment(1)%N_vec(2)
            ray_segment(1)%N_vec(2) = temp
          end if
          call find_first_point_in_plasma(plasma_params, omega, mode, ray_segment, last_N, wall_hits, LOS_end)
          if(LOS_end) wall_hits = 2
          N_init = last_N
          if(debug_level > 0 .and. output_level .and. .not. LOS_end) then
            print*, "First point in plasma",ray_segment(last_N)%R_vec
          end if
          if(last_N  + 1 <= max_points_svec .and. .not. LOS_end) then
            call make_ray_segment(20.d0, plasma_params, omega, mode, ray_segment, last_N, wall_hits, N_init)
          else if(.not. LOS_end) then
            print*,"Ray reached maximum length when searching for first point in vessel"
            print*, "Most likely something is very wrong the launching geometry of the diagnostic"
            print*, "Current diagnostic", ant%diag(idiag)%diag_name
            print*, "position and launch vector in Carthesian coordinates", ray_segment(1)%x_vec, &
              ray_segment(1)%N_vec
            stop "Error when finding first point in plasma in mod_raytrace.f90"
          else
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes = .false.
            if(output_level) print*, "Warning a ray did not pass through the vessel"
          end if
          if(last_N >= max_points_svec) then
            print*, "WARNING a ray did not reach the plasma wall"
            print*, "From here on  output is only for debugging purposes"
            wall_hits = 2
          end if
          if( wall_hits < 2) then
            debug_level = 1
            call find_first_point_in_plasma(plasma_params, omega, mode, ray_segment, last_N, wall_hits, LOS_end)
            N_init = last_N
            call make_ray_segment(20.d0, plasma_params, omega, mode, ray_segment, last_N, wall_hits, N_init)
            print*, "Ray in span_svecs did not end at a wall"
            print*, ray_segment(1)%R_vec(1), ray_segment(1)%R_vec(2), ray_segment(1)%R_vec(3)
            print*, ray_segment(last_N)%R_vec(1), ray_segment(last_N)%R_vec(2), ray_segment(last_N)%R_vec(3)
            print*, "Traveled distance", ray_segment(last_N)%s - ray_segment(N_init)%s
            stop "Error with rays in mod_raytrace.f90"
          end if
          !print*, "plasma", last_N
          if(output_level) then
          ! Copy ray information to the ray_extra_output array
            if(.not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output)) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
            if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s)) then
               deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta)
            end if
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s(last_N), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x(last_N), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y(last_N), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z(last_N), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H(last_N), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray(last_N), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold(last_N), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop(last_N), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta(last_N))
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N = last_N
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s = ray_segment(1:last_N)%s
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x = ray_segment(1:last_N)%x_vec(1)
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y = ray_segment(1:last_N)%x_vec(2)
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z = ray_segment(1:last_N)%x_vec(3)
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H = ray_segment(1:last_N)%Hamil
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray = ray_segment(1:last_N)%N_s
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop = ray_segment(1:last_N)%rhop
            rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta = ray_segment(1:last_N)%theta
            do N=1, last_N
              X = func_X(plasma_params,omega, ray_segment(N)%n_e, ray_segment(N)%T_e)
              Y = func_Y(plasma_params,omega, ray_segment(N)%omega_c / e0 * mass_e, ray_segment(N)%T_e)
              rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold(N) = func_N(X, Y, &
                ray_segment(N)%theta, -rad%diag(idiag)%ch(ich)%mode(imode)%mode)
            end do
            if(.not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)) then
            ! This is important in the case we want some extra output for the last ida optimization
              allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
            end if
          end if
!          print*, "Ray start s", ray_segment(1)%s
!          print*, "Ray start", ray_segment(1)%x_vec
!          print*, "ray launch", ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec
!          print*, "ray extra_output start s", rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s(1)
!          print*, "ray extra_output start", rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x(1), &
!          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y(1), &
!          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z(1)
          ! Reverse the entire ray such that the propagation direction is towards the antenna
          ray_segment(1:last_N) = ray_segment(last_N:1:-1)
          ray_segment(1:last_N)%s = ray_segment(1)%s - ray_segment(1:last_N)%s
          !last_N = last_N - N_init - 1 ! exclude straight line part in vacuum
          if(rad%diag(idiag)%ch(ich)%mode(imode)%mode > 0) then
          ! X-mode -> first and second harmonic problematic
            call make_s_grid(plasma_params, omega, Y_res_X, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, &
                             ray_segment, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, &
                             last_N, dist, rad_ray_freq=rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq))
          else
          ! O-mode -> Only first harmonic has strong absorption
            call make_s_grid(plasma_params, omega, Y_res_O, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, &
                             ray_segment, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, &
                             last_N, dist, rad_ray_freq=rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq))
          end if
          if(output_level) then
            call interpolate_svec(plasma_params, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, ray_segment, omega, &
                                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, last_N, &
                                  rad_ray_freq=rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                                  svec_extra_output=rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)
          else
            call interpolate_svec(plasma_params, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, ray_segment, omega, &
                                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, last_N, &
                                  rad_ray_freq=rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq))
          end if

        end do !ir
        do ir=1, N_ray
          ifreq = 1
          if(.not. (rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%use_external_pol_coeff .and. &
             rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff == 0.d0)) then
            if(output_level) then
            ! Prepare the ray_extra_output arrays for radiation transport quantities (size of svec)
              if(.not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad)) then
                allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(max_points_svec), &
                        rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(max_points_svec))
              end if
            end if
            if(output_level) then
              if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
                print*, "Channel", ich, "ray", ir, "X-mode ", "initialized", "(", &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, "points in total)"
              else
                print*, "Channel", ich, "ray", ir, "O-mode ", "initialized", "(", &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, "points in total)"
              end if
            end if
            if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points == 0) cycle
            !print*, "s_res central", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%s_res
            do ifreq = 2, N_freq ! Copy ray to the other frequencies
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%max_points_svec_reached = .false.
              if(output_level .and. .not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)) then
                ! This is important in the case we want some extra output for the last ida optimization
                allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(:) = &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(:)
              end if
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:) = &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(:)
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points = &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
               ! each frequency and ray has an individual resonance position, here we obtain the average using the weights
            end do ! ifreq
            do ifreq = 1, N_freq
              !if(output_level) print*, "s_res ", ifreq, "-th frequency", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%s_res
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res + &
                ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%s_res
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res + &
                ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%R_res
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res + &
                ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%z_res
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res + &
                ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%rhop_res
            end do
          end if ! Skip exact determination of resonance if polarization filter perfect
          rad%diag(idiag)%ch(ich)%mode(imode)%s_res = rad%diag(idiag)%ch(ich)%mode(imode)%s_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res
          rad%diag(idiag)%ch(ich)%mode(imode)%R_res = rad%diag(idiag)%ch(ich)%mode(imode)%R_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res
          rad%diag(idiag)%ch(ich)%mode(imode)%z_res = rad%diag(idiag)%ch(ich)%mode(imode)%z_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res
          rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res
        end do ! ir
        ! stop "check ray_launch"
      end do ! imode
      if(modes == 1) then ! X-mode
      ! set the resonance position of the channel to the resonance position of the X-mode
      ! if multiple mode analysis is switched on this will be overwritten later
      ! The step is neccessary to quickly perform classical ECE analysis within IDA
        rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%mode(1)%s_res
        rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%mode(1)%R_res
        rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%mode(1)%z_res
        rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
      end if
    end do !ich
#ifdef OMP
    !$omp end do
#endif
    deallocate(ray_segment, x_loc_vec, N_loc_vec)
#ifdef OMP
    !$omp end parallel
#endif
    rad%diag(idiag)%ch(:)%eval_ch = .true. ! start with all channels set to true
    if(output_level) then
      do ich = 1, ant%diag(idiag)%N_ch
        if(.not. allocated(rad%diag(idiag)%ch(ich)%mode_extra_output)) then
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
        end if
        do imode = 1, mode_cnt
          if(.not. allocated(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s)) then
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Te(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD(pnts_BPD))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(pnts_BPD))
            allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary(pnts_BPD))
          end if
        end do
      end do
    end if
    largest_svec = 0
    do ich = 1, ant%diag(idiag)%N_ch
      do imode = 1, mode_cnt
        do ir = 1, N_ray
          do ifreq = 1, N_freq
            ! Do this outside the openMP loop for thread safety
            if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points > largest_svec) &
              largest_svec = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
          end do
        end do
      end do
    end do
  end do ! idiag
  !if(output_level) print*, "largest svec", largest_svec
  !stop "Early ray end?"
  end subroutine span_svecs

  subroutine reinterpolate_svec(svec, max_points_svec_reached, total_LOS_points, omega,  mode, plasma_params, ds1, ds2, svec_extra_output)
  use f90_kind
  use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, &
                                        plasma_params_type, output_level, &
                                        ray_element_full_type, &
                                        rad_diag_ch_mode_ray_freq_svec_extra_output_type
  use f90_kind
  use constants,                  only: pi
  implicit none
  type(rad_diag_ch_mode_ray_freq_svec_type), dimension(:), intent(inout)  :: svec
  logical, intent(inout)                                                  :: max_points_svec_reached
  integer(ikind), intent(inout)                                           :: total_LOS_points
  real(rkind), intent(in)                                                 :: omega
  integer(ikind), intent(in)                                              :: mode
  type(plasma_params_type), intent(in)                                    :: plasma_params
  real(rkind), intent(in)                                                 :: ds1, ds2
  type(rad_diag_ch_mode_ray_freq_svec_extra_output_type), dimension(:), intent(inout), optional  :: svec_extra_output
  type(ray_element_full_type), dimension(total_LOS_points)                :: ray_segment
  integer(ikind)                                                          :: last_N, k
  real(rkind), dimension(2)                                               :: dist, Y_res_X
  real(rkind), dimension(1)                                               :: Y_res_O
  logical                                                                 :: max_points_in_svec_reached
  if(max_points_svec_reached) then
    print*, "reinterpolate_rad_ray_freq should not be called for frequencies which were"
    print*, "identified to have insufficient amount of max_points_svecs to reinterpolate"
    call abort()
  end if
  dist(1) = ds1 * (plasma_params%int_step_cnt)
  dist(2) = ds2 * (plasma_params%int_step_cnt)
  Y_res_O(1) = 1.d0
  Y_res_X(1) = 1.d0
  Y_res_X(2) = 0.5d0
  ray_segment%s = svec(1:total_LOS_points)%s
  do k = 1, 3
    ray_segment%x_vec(k) = svec(1:total_LOS_points)%x_vec(k)
    ! Turn this around to make this a ray that goes from outside the plasma into the plasma
    ! Theta is computed from N and B
    ray_segment%N_vec(k) = -svec(1:total_LOS_points)%N_vec(k)
    ray_segment%B_vec(k) = svec(1:total_LOS_points)%B_vec(k)
  end do
  ray_segment%rhop = svec(1:total_LOS_points)%rhop
  ray_segment%omega_c = svec(1:total_LOS_points)%freq_2X * pi
  ray_segment%T_e = svec(1:total_LOS_points)%Te
  ray_segment%n_e = svec(1:total_LOS_points)%ne
  ray_segment%v_g_perp = svec(1:total_LOS_points)%v_g_perp
  ray_segment%N_s = svec(1:total_LOS_points)%N_cold
  if(output_level) then
    if(.not. present(svec_extra_output)) then
      print*, "Reinterplate_svec has to be called with svec_extra_output present if output_level is true"
      call abort()
    end if
    ray_segment%Hamil= svec_extra_output(1:total_LOS_points)%H
  end if
  max_points_in_svec_reached = .false.
  last_N = total_LOS_points
  if(mode > 0) then
    call make_s_grid(plasma_params, omega, Y_res_X, svec, ray_segment, total_LOS_points, last_N, &
                     dist, max_points_svec_reached=max_points_in_svec_reached)
  else
    call make_s_grid(plasma_params, omega, Y_res_O, svec, ray_segment, total_LOS_points, last_N, &
                     dist, max_points_svec_reached=max_points_in_svec_reached)
  end if
  if(max_points_in_svec_reached) then
    print*, "Tried to reduce step size in radiation transport equation, but insufficient amount of points for svec"
    print*, "Current large/small step size", ds1, " / ",  ds2
    print*, "Largest Te/ne on grid", maxval(svec(1:total_LOS_points)%Te), " / ", &
                                     maxval(svec(1:total_LOS_points)%ne)
    if(output_level) then
      print*, "If Te/ne is reasonable increase max_points_svec and rerun program"
      call abort()
    else
        total_LOS_points = last_N
        svec(1:total_LOS_points)%s = ray_segment%s
      print*, "Flagging this channel as at maximum amount of steps and continuing"
      print*, "------------------------------------------------------------------------------------------"
      print*, "|  Warning: ECFM encountered irresolvable numerical difficulties in radiation transport!  |"
      print*, "|   Since this problem might not affect the final result of IDA, computation continues!   |"
      print*, "------------------------------------------------------------------------------------------"
    end if
    max_points_svec_reached = .true.
    return
  else
    if(present(svec_extra_output)) then
      call interpolate_svec(plasma_params, svec, ray_segment, omega, &
                            total_LOS_points, last_N, &
                            svec_extra_output=svec_extra_output)      
    else
      call interpolate_svec(plasma_params, svec, ray_segment, omega, &
                            total_LOS_points, last_N)
    end if
  end if
  !print*, "s_max after", svec(total_LOS_points)%s
  !print*, "total los points", total_LOS_points
  end subroutine reinterpolate_svec





  subroutine dealloc_rad(rad)
  use mod_ecfm_refr_types,        only: rad_type, ant, ray_init, N_ray, N_freq, mode_cnt, stand_alone, output_level
  use f90_kind
  use constants,                  only: pi,e0, mass_e, c0
  implicit none
  type(rad_type), intent(inout)    :: rad
  integer(ikind)                                                :: idiag, ich, ir, ifreq, imode
  if(stand_alone) then
    print*, "There is no reason to call dealloc_rad in mod_raytrace.f90 in stand_alone mode"
    stop "stand_alone = T in dealloc_rad"
  end if
  do idiag = 1, ant%N_diag
    do ich = 1, ant%diag(idiag)%N_ch
      do imode = 1, mode_cnt
        do ir = 1, N_ray
          do ifreq = 1, N_freq
            deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec)
            if(output_level) then
              if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)) then
                deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)
              end if
            end if
          end do
          deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq)
          if(output_level .and. ray_init) then
            if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s)) then
              deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold)
              deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary, &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary)
            end if
          end if
        end do
        deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray)
        if(output_level) then
          if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output)) then
            deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output)
          end if
          if(ray_init) then
            if(allocated(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s)) then
              deallocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Te, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm)
            end if
          end if
        end if
      end do
      deallocate(rad%diag(idiag)%ch(ich)%mode)
      if(output_level .and. allocated(rad%diag(idiag)%ch(ich)%mode_extra_output)) then
        deallocate(rad%diag(idiag)%ch(ich)%mode_extra_output)
      end if
    end do
    deallocate(rad%diag(idiag)%ch)
  end do
  deallocate(rad%diag)
  end subroutine dealloc_rad

end module mod_ecfm_refr_raytrace
