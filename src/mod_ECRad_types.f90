!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ECRad_types
! Many of the values set here are overwritten during intialization
! values marked with a * are hardcoded here
use f90_kind
use constants,                  only: pi
use magconfig,                only: Scenario_type
! #ifdef OMP
! use omp_lib
! #endif
#ifdef NAG
USE nag_spline_2d             , only: nag_spline_2d_comm_wp => nag_spline_2d_comm_dp
USE nag_spline_1d             , only: nag_spline_1d_comm_wp => nag_spline_1d_comm_dp
#endif
implicit none
integer(ikind), dimension(:), allocatable :: omit_ch


type spl_type_2d
  !integer*4, dimension(3)               :: iopt
  integer*4                             :: nu, nv, nuest, nvest, lwrk, &
                                           kwrk, lwrk_in_spl, kwrk_in_spl
                                            ! spline knot count, length of real and integer work array
  real*8, dimension(:), allocatable     :: tu, tv ! spline knot positions
  real*8, dimension(:,:), allocatable   :: c !spline
  real*8, dimension(:), allocatable     :: wrk ! work array for the spline routine
  integer*4, dimension(:), allocatable  :: iwrk ! integer work array for fitpack
  integer*4                             :: iopt_int = 0! needed for rectangular grid
  real*8                               :: x_start, x_end
  real*8                               :: y_start, y_end
end type spl_type_2d

type spl_type_1d
  !integer*4, dimension(3)               :: iopt, ider
  integer*4                             :: n, nest, lwrk, &
                                           kwrk, k
                                            ! spline knot count, length of real and integer work array
  real*8, dimension(:), allocatable     :: t ! spline knot positions
  real*8, dimension(:), allocatable   :: c !spline
  real*8, dimension(:), allocatable     :: wrk ! work array for the spline routine
  integer*4, dimension(:), allocatable  :: iwrk ! integer work array for fitpack
  integer*4                             :: iopt_int = 0 ! needed for profile interpolation
  real*8                               :: x_start, x_end
end type spl_type_1d
! antenna and LOS parameters (not time-dependent)

type ant_diag_ch_launch_type
  real(rkind), dimension(3)        :: x_vec, N_vec
  real(rkind)                      :: R, phi, z, theta_pol, phi_tor ! TORBEAM launching angles
  real(rkind)                      :: weight     ! weighting of ray with gaussian
end type ant_diag_ch_launch_type


type ant_diag_ch_type
  !real(rkind)                      :: Tunc_ECE   ! SNR & calibration error
  real(rkind)                      :: f_ECE      ! mean frequency of ECE channel [Hz]
  real(rkind)                      :: df_ECE     ! bandwidth of ECE frequencies  [Hz]
  real(rkind)                      :: dist_focus, width ! distance between launch and focus and FWHM of the antenna
  real(rkind)                      :: focus_shift
  !integer(ikind)                  :: N_freq     ! number of frequencies per ECE channels
  !integer(ikind)                  :: N_ray      ! number of rays
  real(rkind),  dimension(:), allocatable :: freq, freq_weight, freq_int_weight ! analyzed frequencies, if band shape and weights for the integration
  type(ant_diag_ch_launch_type),  dimension(:), allocatable :: ray_launch  ! launching angles for multiple rays
end type ant_diag_ch_type

type ant_diag_type
  integer(ikind)                   :: N_ch       ! number of ECE channels
  type(ant_diag_ch_type),               dimension(:), allocatable :: ch   ! ECE channel
  character(8)                     :: expname
  character(3)                     :: diag_name
  integer(ikind)                   :: edition
end type ant_diag_type

type ant_type
  integer(ikind)                   :: N_diag       ! number of ECE channels
  type(ant_diag_type),               dimension(:), allocatable :: diag   ! ECE channel
end type ant_type



type Use_3D_vessel_type
  real(rkind), dimension(:,:), allocatable :: vessel_data_R
  real(rkind), dimension(:,:), allocatable :: vessel_data_z
  real(rkind), dimension(:), allocatable   :: phi
  integer(ikind)                           :: n_phi, n_contour
  real(rkind)                              :: phi_max
end type

! radiation and LOS parameters (time-dependent)
type rad_diag_ch_mode_ray_freq_svec_type ! S. Denk 4. 2013
! This structure saves all the data neccessary for the integration along LOS
! for straight lines this avoids interpolations in the most inner loop of the program
! for bend lines of sight this is mainly used to store data for debug and interpretation purposes
  real(rkind)                                 :: s
  real(rkind)                                 :: R
  real(rkind)                                 :: z
  real(rkind)                                 :: theta
  real(rkind)                                 :: cos_theta
  real(rkind)                                 :: sin_theta
  real(rkind)                                 :: rhop
  real(rkind)                                 :: freq_2X
  real(rkind)                                 :: Te
  real(rkind)                                 :: ne
  real(rkind)                                 :: ibb
  real(rkind)                                 :: N_cold
  real(rkind)                                 :: v_g_perp ! group velocity for Farina's absorption coefficient
  real(rkind), dimension(3)                   :: x_vec, N_vec, B_vec ! for polarization filter
  logical                                     :: plasma ! only true if within vessel
end type rad_diag_ch_mode_ray_freq_svec_type

type rad_diag_ch_mode_ray_freq_svec_spline_type
  type(spl_type_1d)                           :: x, y, z, rhop, Nx, Ny, Nz, Bx,By, Bz, theta, freq_2X, Te, ne, v_g_perp
  real(rkind)                                 :: s_min, s_max, s_center ! end points of the region on which the splines are defined
end type rad_diag_ch_mode_ray_freq_svec_spline_type

type rad_diag_ch_mode_ray_freq_svec_extra_output_type
  real(rkind)                                 :: N_warm
  real(rkind)                                 :: N_cor
  real(rkind)                                 :: Trad
  real(rkind)                                 :: Trad_secondary
  real(rkind)                                 :: em
  real(rkind)                                 :: em_secondary
  real(rkind)                                 :: T
  real(rkind)                                 :: T_secondary
  real(rkind)                                 :: ab
  real(rkind)                                 :: ab_secondary
  real(rkind)                                 :: H
  real(rkind)                                 :: N_ray
end type rad_diag_ch_mode_ray_freq_svec_extra_output_type

type rad_diag_ch_mode_ray_freq_type
  real(rkind)                      :: Trad   ! resonance position [m]
  real(rkind)                      :: Trad_secondary! resonance position
  real(rkind)                      :: tau
  real(rkind)                      :: tau_secondary
  real(rkind)                      :: pol_coeff
  real(rkind)                      :: pol_coeff_secondary
!  real(rkind)                      :: X_refl ! fraction of O-mode converted to X-mode
!  real(rkind)                      :: X_refl_secondary ! fraction of O-mode converted to X-mode
  real(rkind)                      :: s_res ! each ray and frequency have individual resonances
  real(rkind)                      :: R_res ! each ray and frequency have individual resonances
  real(rkind)                      :: z_res ! each ray and frequency have individual resonances
  real(rkind)                      :: rhop_res ! each ray and frequency have individual resonances
                                                  ! For conistent treament of rays the resonances have to be also averaged for one channel
  real(rkind)                      :: rel_s_res ! relativistic resonance position [m]
  real(rkind)                      :: rel_rhop_res
  real(rkind)                      :: rel_R_res
  real(rkind)                      :: rel_z_res
  logical                          :: use_external_pol_coeff = .false.
  logical                          :: max_points_svec_reached  = .false. ! In case it is not possible to refine grid size further
  integer(ikind)                   :: total_LOS_points
  type(rad_diag_ch_mode_ray_freq_svec_type), dimension(:), allocatable  :: svec ! Vector of length total_LOS_points
  type(rad_diag_ch_mode_ray_freq_svec_extra_output_type), dimension(:), allocatable  :: svec_extra_output ! Vector of length total_LOS_points
  type(rad_diag_ch_mode_ray_freq_svec_spline_type)                  :: svec_spline
end type rad_diag_ch_mode_ray_freq_type

type rad_diag_ch_mode_ray_type
  type(rad_diag_ch_mode_ray_freq_type), dimension(:), allocatable :: freq ! Vector of length N_freq
  logical                          :: low_field_side
  real(rkind)                      :: s_axis
  real(rkind)                      :: Trad
  real(rkind)                      :: Trad_secondary
  real(rkind)                      :: tau
  real(rkind)                      :: tau_secondary
  real(rkind)                      :: s_res ! each ray and frequency have individual resonances
  real(rkind)                      :: R_res ! each ray and frequency have individual resonances
  real(rkind)                      :: z_res ! each ray and frequency have individual resonances
  real(rkind)                      :: rhop_res ! each ray and frequency have individual resonances
  real(rkind)                      :: rel_s_res ! relativistic resonance position [m]
  real(rkind)                      :: rel_rhop_res
  real(rkind)                      :: rel_R_res  ! relativistic resonance position [m]
  real(rkind)                      :: rel_z_res  ! relativistic resonance position [m]
  real(rkind)                      :: rel_s_res_secondary ! relativistic resonance position [m]
  real(rkind)                      :: rel_rhop_res_secondary
  real(rkind)                      :: rel_R_res_secondary  ! relativistic resonance position [m]
  real(rkind)                      :: rel_z_res_secondary  ! relativistic resonance position [m]
  logical                          :: contributes ! False for channels that do not cross the vessel
end type rad_diag_ch_mode_ray_type

type rad_diag_ch_mode_ray_extra_output_type
  real(rkind), dimension(:), allocatable :: s
  real(rkind), dimension(:,:), allocatable :: x_vec, B_vec, N_vec
  real(rkind), dimension(:), allocatable :: H, N_ray, N_cold, rhop, theta, Te, ne, v_g_perp
  integer(ikind)                         :: N
  real(rkind), dimension(:), allocatable :: Trad, Trad_secondary, em, em_secondary, ab, ab_secondary, T, T_secondary, BPD, BPD_secondary
end type rad_diag_ch_mode_ray_extra_output_type

type rad_diag_ch_mode_type
  real(rkind)                      :: Trad
  real(rkind)                      :: Trad_secondary
  real(rkind)                      :: tau
  real(rkind)                      :: tau_secondary
  real(rkind)                      :: pol_coeff
  real(rkind)                      :: pol_coeff_secondary
!  real(rkind)                      :: X_refl ! fraction of O-mode converted to X-mode
!  real(rkind)                      :: X_refl_secondary ! fraction of O-mode converted to X-mode
  integer(ikind)                   :: mode ! - 1 O-mode, +1 X-mode
  real(rkind)                      :: s_res
  real(rkind)                      :: R_res ! modes have different propagation
  real(rkind)                      :: z_res
  real(rkind)                      :: rhop_res
  real(rkind)                      :: rel_s_res ! relativistic resonance position [m] -makes sense to do these for each mode individually
  real(rkind)                      :: rel_rhop_res
  real(rkind)                      :: rel_R_res
  real(rkind)                      :: rel_z_res
  real(rkind)                      :: rel_s_res_secondary ! relativistic resonance position [m] -makes sense to do these for each mode individually
  real(rkind)                      :: rel_rhop_res_secondary
  real(rkind)                      :: rel_R_res_secondary
  real(rkind)                      :: rel_z_res_secondary
  type(rad_diag_ch_mode_ray_type),      dimension(:), allocatable :: ray
  type(rad_diag_ch_mode_ray_extra_output_type),      dimension(:), allocatable :: ray_extra_output
  real(rkind)                      :: Trad_mode_frac, Trad_mode_frac_secondary
end type rad_diag_ch_mode_type

type rad_diag_ch_mode_extra_output_type
  real(rkind), dimension(:), allocatable :: s, R,z, Trad, Trad_secondary, &
                                            em, em_secondary, ab, ab_secondary, &
                                            T, T_secondary, N_cold, N_cor, N_warm, Te, ne
  real(rkind), dimension(:), allocatable :: rhop_BPD, BPD, BPD_secondary
end type rad_diag_ch_mode_extra_output_type

type rad_diag_ch_type
  real(rkind)                      :: T_ECE   ! radiation temperature from ECE [eV]
  real(rkind)                      :: Trad    ! total radiation temperature
  real(rkind)                      :: tau     ! optical thickness
  real(rkind)                      :: Trad_secondary    ! total radiation temperature (secondary model)
  real(rkind)                      :: tau_secondary     ! optical thickness (secondary model)
  real(rkind)                      :: rel_s_res ! relativistic resonance position [m]
  real(rkind)                      :: rel_rhop_res
  real(rkind)                      :: rel_R_res
  real(rkind)                      :: rel_z_res
  real(rkind)                      :: rel_s_res_secondary ! relativistic resonance position [m]
  real(rkind)                      :: rel_rhop_res_secondary
  real(rkind)                      :: rel_R_res_secondary
  real(rkind)                      :: rel_z_res_secondary
  real(rkind)                      :: s_res
  real(rkind)                      :: R_res ! each ray and frequency have individual resonances
  real(rkind)                      :: z_res ! each ray and frequency have individual resonances
  real(rkind)                      :: rhop_res! resonance position
  type(rad_diag_ch_mode_type),      dimension(:), allocatable :: mode !  different rays within beam
  type(rad_diag_ch_mode_extra_output_type),      dimension(:), allocatable :: mode_extra_output
  logical                          :: eval_ch
  real(rkind)                      :: pol_coeff_norm, pol_coeff_secondary_norm
end type rad_diag_ch_type


type rad_diag_type
  type(rad_diag_ch_type),          dimension(:), allocatable :: ch
end type rad_diag_type

type rad_type
  type(rad_diag_type),             dimension(:), allocatable :: diag
end type rad_type



type ffp_type
  real(rkind)                        :: u_step, pitch_step
  real(rkind)                        :: rhop_min, rhop_max, u_min = 5.d-2, u_max, pitch_min, pitch_max
  real(rkind)                        :: u_max_data, f_min = 1.e-22
                                        ! u_min is harcoded to avoid numerical difficulties near u = 0
  integer(ikind)                     :: N_rhop, N_pitch, N_u
  real(rkind)                        :: delta_u_par = 1.d-4 !for differentiation
  real(rkind)                        :: delta_u_perp = 1.d-4 !for differentiation
  real(rkind), dimension(:,:,:), allocatable :: f
  real(rkind), dimension(:), allocatable :: rhop, u, pitch
#ifdef NAG
  logical                            :: nag = .true.
#else
  logical                            :: nag = .false.
#endif
  logical                            :: restart_spline
  type(spl_type_1d), dimension(:, :), allocatable :: f_rhop_spl ! spline structure for interpolation along rhop
  integer(ikind)                     :: N_B_min = 100 !&
  real(rkind), dimension(:), allocatable :: rhop_B_min, B_min
  type(spl_type_1d)                  :: B_min_spl ! spline structure for interpolation along rhop
#ifdef NAG
  type(nag_spline_1d_comm_wp)        :: B_min_nag_spl
  type(nag_spline_1d_comm_wp) , dimension(:, :), allocatable   :: f_rhop_nag_spl
#endif
  logical                            :: LUKE
  logical                            :: restrict_u_max = .false. !* for distribution with bad date in the high energy regime
end type ffp_type


type fgene_type
  real(rkind)                        :: vpar_step, mu_step
  real(rkind)                        :: rhop_min, rhop_max, vpar_min, vpar_max, mu_min, mu_max
  real(rkind)                        :: vpar_max_data
                                        ! u_min is harcoded to avoid numerical difficulties near vpar = 0
  integer(ikind)                     :: N_rhop, N_mu, N_vpar
  real(rkind)                        :: delta_u_par = 1.d-4 !for differentiation
  real(rkind)                        :: delta_u_perp = 1.d-4 !for differentiation
  real(rkind)                        :: B0 ! Used for the coordinate remapping from mu to beta_perp
  real(rkind), dimension(:), allocatable :: rhop, vpar,mu
  real(rkind), dimension(:,:,:), allocatable :: g
#ifdef NAG
  type(nag_spline_2d_comm_wp)        :: f0_nag_spl
  logical                            :: nag = .true.
#else
  logical                            :: nag = .false.
#endif
  logical                            :: restart_spline
  type(spl_type_2d)                  :: f0_spl ! f0 for comparison to unperturbed distribution
  type(spl_type_1d), dimension(:, :), allocatable :: g_rhop_spl ! spline structure for interpolation along rhop
  type(spl_type_1d)                  :: Te_perp_spl, Te_par_spl ! For sanity testing the spline interplation
#ifdef NAG
  type(nag_spline_1d_comm_wp) , dimension(:, :), allocatable   :: g_rhop_nag_spl
  type(nag_spline_1d_comm_wp)        :: Te_perp_nag_spl, Te_par_nag_spl ! For sanity testing the spline interplation
#endif
  real(rkind), dimension(:), allocatable :: Te_perp_vec, Te_par_vec
end type fgene_type

type bi_max_type
  type(spl_type_1d)                      :: ne_spl
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: ne_nag_spl
#endif
  real(rkind)                            :: rhop_max, Te_perp, Te_par
end type bi_max_type

type runaway_type
  type(spl_type_1d)                      :: rel_ne_spl
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: rel_ne_nag_spl
#endif
  real(rkind)                            :: Z_eff, E_E_c, rhop_max
end type runaway_type

type drift_m_type
  type(spl_type_1d)                      :: ne_spl
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: ne_nag_spl
#endif
  real(rkind)                            :: rhop_max, Te_perp, Te_par, u_par_drift, u_perp_drift
end type drift_m_type

type multi_slope_type
  type(spl_type_1d)                      :: Te_slope_spl
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: Te_slope_nag_spl
#endif
  real(rkind)                            :: rhop_max, gamma_switch
end type multi_slope_type

type Spitzer_type
  type(spl_type_1d)                      :: j_spl
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: j_nag_spl
#endif
  real(rkind)                            :: rhop_max
end type Spitzer_type

type non_therm_params_type
  real(rkind)                      :: Te_perp, Te_par
  real(rkind)                      :: ne, j
  real(rkind)                      :: cur_rel_ne
  real(rkind)                      :: u_par_drift, u_perp_drift
  real(rkind)                      :: mu_slope, norm, total_norm
  real(rkind)                      :: vd
  real(rkind)                      :: B_min ! Current B_min for Fokker-Planck distributions
end type

type grad_type
  real(rkind)               :: dR, dphi, dz
end type grad_type

type contour_type
  real(rkind), dimension(:), allocatable              :: x, y
end type contour_type

type ripple_params_type
! Parametrisation of the field ripple
  real(rkind)               :: A0, A1, A2
  real(rkind)               :: B0, B1, B2
  real(rkind)               :: C0, C1, C2
  real(rkind)               :: D0, D1, D2
  real(rkind)               :: E0, E1, K1, K2
end type ripple_params_type

type ray_element_full_type
  real(rkind), dimension(3)                         :: x_vec, N_vec, R_vec, N_R_vec, B_vec
  real(rkind)                                       :: s, h
  real(rkind)                                       :: theta
  real(rkind)                                       :: Hamil, N_s, n_e, omega_c, rhop, T_e, v_g_perp ! group velocity
  real(rkind)                                       :: sigma
  integer(ikind)                                    :: first_N_plasma, last_N_plasma
end type ray_element_full_type

type ext_ray_type
  type(ray_element_full_type), dimension(:), allocatable :: ray
  integer(ikind)                                         :: N_steps
end type ext_ray_type

type plasma_params_type
  logical                                           :: on_the_fly_raytracing = .false.
                                                       ! If .true. raytracing and radiation transport are solved simultaneously
                                                       ! if .false. the rays are calculated first during initializiation
#ifdef NAG
  type(nag_spline_2d_comm_wp)                       :: rhop_spline_nag, B_R_spline_nag, B_t_spline_nag, &
                                                       B_z_spline_nag, Te_spline_nag_2D, ne_spline_nag_2D
#endif
  type(spl_type_2d)                                 :: rhop_spline, B_R_spline, B_t_spline, &
                                                       B_z_spline, T_e_spline_2D, n_e_spline_2D
  real(rkind), dimension(:,:), allocatable          :: rhop
  type(spl_type_1d)                                 :: ne_spline, Te_spline
  ! The spline results produces bad derivatives - do not use
  !real(rkind), dimension(:,:,:), allocatable        :: rhop_spline, B_r_spline, B_t_spline, B_z_spline
  !integer(ikind), dimension(:,:), allocatable       :: R_z_ipoint
#ifdef NAG
  type(nag_spline_1d_comm_wp)                       :: ne_spline_nag, Te_spline_nag
#endif
 real(rkind), dimension(:), allocatable             :: IDA_rhop_knots_Te, IDA_T_e, &
                                                       IDA_T_e_dx2, IDA_rhop_knots_ne, IDA_n_e, &
                                                       IDA_n_e_dx2
  type(contour_type)                                :: vessel_poly ! polynome describing the vessel (2D)
  integer(ikind)                                    :: shot, eq_ed
  integer(ikind)                                    :: ida_time_indx
  real(rkind)                                       :: time
  real(rkind), dimension(:), allocatable            :: rhop_vec_ne, n_e_prof, rhop_vec_Te, T_e_prof, R, z, B_t
  real(rkind), dimension(:), allocatable            :: Int_absz, Int_weights ! for radiation transport w. ray tracing
  real(rkind)                                       :: R_min, R_max, z_min, z_max,R_step, z_step, rhop_max = 1.2
                                                       ! rhop_max will be overwritten if a flux matrix is read
  integer(ikind)                                    :: m, n, m_n_e_prof, m_T_e_prof, m_vessel_bd
  character(4)                                      :: eq_exp
  character(3)                                      :: eq_diag
  real(rkind)                                       :: z_lens
  integer(ikind)                                    :: mode = 0! mode = -1 -> X,  mode = +1 -> O
  logical                                           :: w_ripple = .true. ! Input
                                                    ! .true. -> weakly relativistic. false -> cold disperion for ray tracing
  real(rkind)                                       :: Btf0 = -1.d2, R0 = -1.d0 ! Required for the ripple correction
  logical                                           :: precise_interpolation = .false. !*
                                                    ! if true use splines for svec interpolation - SLOW!
                                                    ! if false use linear interpolation (much faster)
  logical                                           :: Te_ne_mat = .false. ! If true an externally given matrix of Te and ne is used instead of the
                                                                              ! Te/ ne profile
  integer(ikind)                                    :: Debug_level = 1 !* Controls the amount of output of raytracing
                                                      ! (0) No output, (1) some output regarding spline interpolation, (2) much output and stop at first interesting point
  real(rkind)                                       :: h = 1.d-3 !* step size for raytracing - WARNING small values will increase the error!
  real(rkind)                                       :: ab_switch = 1.d0 !* If ab * (b - a) > ab_switch switch to the smallest step size in rad_int (default 5)
  real(rkind)                                       :: rhop_scale_te = 1.00, rhop_scale_ne = 1.00 ! rhop will be multiplied with the corresponding scale when evaluating Te and ne
                                                       ! does not apply for analytical data and does not affect the resonance positions
  real(rkind)                                       :: H_last = 1.d-10, trigger_level = 10 !* Used to detect large jump in H
  real(rkind)                                       :: h_min = 5.e-5, h_max = 4.e-3 !* for the adaptive step size we need boundaries
  real(rkind)                                       :: R_ax, z_ax, B_ax ! for HFS, LFS distinction
  real(rkind)                                       :: rhop_max_default = 1.2d0 !*
  real(rkind)                                       :: delta_rhop_exit = 0.005!* Used to trigger an early exit of raytracing
                                                                               ! in case of profiles that stop at the LCFS
  real(rkind)                                       :: X_entry = 0.04 !*
  real(rkind)                                       :: rhop_emit = 1.03 !* Specifies upper limit for the fine grid
                                                       ! This avoids having the small grid for channels with resonances in the SOL
  real(rkind)                                       :: down_shift = 0.005, up_shift = -0.005 ! Defines the region of the dense grid Y_res + down_shift < Y < Y_res + upshift
  real(rkind)                                       :: dist_large = 0.0025d0, dist_small = 0.00025d0
                                                       !dist_large = 0.0002d0, dist_small = 0.00002d0
                                                       ! default is 2 mm large and .2 mm small
                                                       ! For strongly radiating plasmas .2 and .02 mm are recommended
                                                       !* De  fine the large and the small step size for the radiation transport
                                                       ! The amoint of steps defined in the next line has no influence on this
  integer(ikind)                                    :: int_step_cnt = 4 ! Has to be multiple of four since we use Rk4
  real(rkind)                                       :: tau_max = 9 !* If tau > tau_max the radiation transport is deemed finished
                                                       ! tau = 5 corresponds to a reabsorption of more than 99% of all radiation
                                                       ! => good enough for thermal plasmas
  real(rkind)                                       :: angle_threshold = 90.d0 / 180.d0 * pi !2.d0 * pi !90.d0 / 180.d0 * pi !* maximum allowed rotation of wave vector k with respect to the launch
                                                       ! 2.d0 * pi ! no threshhold
                                                       ! propagation with respect to launch (avoids internal reflections)
  real(rkind)                                       :: rp_min
  logical                                           :: No_ne_te = .false. ! True for initialization with straight rays
  logical                                           :: prof_log_flag = .true. ! If True Te and ne interpolated by Exp(Spl) instead of Spl directly
  real(rkind)                                       :: R_shift = 0.d0, z_shift = 0.d0 ! Allows shifting the equilbrium - moves entire flux matrix
  real(rkind)                                       :: theta_pol_cor = 0.0 !-1.0 / 180.0 * pi
  type(Scenario_type)                               :: Scenario ! Configuration info of the 3D Equilibrium
  type(Use_3D_vessel_type)                          :: Use_3D_vessel ! Contains 3D Vessel information
  integer(8), dimension(:), allocatable             :: mconf_addresses
end type plasma_params_type

  type(ant_type)                     :: ant
  logical                            :: use_ext_rays = .false. ! Can be overwritten by the input file, if it is long enough
  character(200)                     :: ext_ray_folder
  type(ext_ray_type), dimension(:), allocatable :: ext_rays
  real(rkind)                        :: Te_min = 1.d0
  type(rad_type)                     :: rad
  type(plasma_params_type)           :: plasma_params
  type(ripple_params_type)           :: ripple_params
  character(5)                       :: ode_integrator = "Rk4" !*
  character(200)                     :: data_folder
  character(200)                     :: ray_out_folder ! folder that saves the result of ray tracing
  logical                            :: non_maxwellian = .false.!.true.
  character(7)                       :: dstf = "numeric"!"relamax"!"maxwell"!
  character(3)                       :: dst_data_folder
  character(5)                       :: Ich_name
  character(2)                       :: dstf_comp, dstf_str
  logical                            :: output_level
  type(ffp_type)                     :: ffp
  type(fgene_type)                   :: fgene
  type(bi_max_type)                  :: bi_max
  type(drift_m_type)                 :: drift_m
  type(multi_slope_type)             :: multi_slope
  type(Spitzer_type)                 :: Spitzer
  type(runaway_type)                 :: runaway
  real(rkind)                        :: ratio_for_third_harmonic = -1.d0! When to include third harmonic
  integer(ikind)                     :: N_max = 3
  integer(ikind)                     :: m_eq = 65,n_eq = 129
  real(rkind)                        :: min_density = 1.e18, min_Te = 100.d0, max_Te=1.e6 ! Used to double check the input
  real(rkind)                        :: SOL_ne = 1.e14, SOL_Te = 20.0d-3 ! should be smaller than min_density and min_Te
  real(rkind)                        :: ignore_Te = 10.0, ignore_ne = 1.e16 ! Radiation transport ignores points below this threshhold
                                                                            ! Should be decently larger than the SOL values
  real(rkind)                        :: ne_max = 1.d21 ! Grid points with densities larger than this will be ignored
  real(rkind)                        :: reflec_X
  real(rkind)                        :: reflec_O
  integer(ikind)                     :: reflec_model = 0 ! # mode one is a ppor approximation of an isotropic background
  real(rkind)                        :: vessel_plasma_ratio
  logical                            :: straight
  logical                            :: Lambda_star = .False. ! Lambda star  = Lambda *  Denominator of dispersion relation = no UH resonacne -> Doesn't seem to work properly!!
  integer(ikind)                     :: pnts_BPD = 2000 !*
  real(rkind)                        :: max_rhop_BPD = 1.05
  logical                            :: use_maximum_for_warm_res = .true.
                                        !* If false uses first moment of BPD
  integer(ikind)                     :: N_freq, N_ray ! Number of frequencies in IF, number of rays
  integer(ikind)                     :: max_points_svec! Maximum allowed points on LOS
  integer(ikind)                     :: largest_svec ! Stores largest svecs -> currently only computed but not output
#ifdef OMP
  integer(ikind)                     :: OMP_thread_count
#endif
  logical                            :: warm_plasma ! Input
  real(rkind)                        :: mode_conv = -1.d0 ! conversion coefficient from O- to X-mode for wall reflections
  integer(ikind)                     :: modes = 0 ! whichs mode(s) are considered
                                                             ! 1 -> X, 2 -> O, 3 -> O, X
  integer(ikind)                     :: mode_cnt ! 1,2
  character(2)                       :: eq_mode ! 1D -> Normal case with 1D profiles, 2D -> 2D profiles, 3D -> 3D equilibrium
  real(rkind)                        :: plasma_vac_boundary = 1.05d0 ! boundary at which the polarization of X and O-mode will be calculated
  character(3), dimension(:), allocatable :: diagnostics ! Separates individual channel bunches into separate diagnostics
  real(rkind)                        :: h_x_glob = 1.d-5 !*
  character(4)                       :: Hamil ="Dani"!"Stix"! !*! !
  logical                            :: UH_stop = .true. !* If true raytracing concludes when approaching the UH resonance (only X-mode affected)
  logical                            :: one_sigma_width = .true. !* Rays are spanned over one beam widths -- useful for benchmarking, but could lead wrong results
  real(rkind)                        :: ECE_fcous_shift = 0.2d0 ! Shifts ECE focus point towards HFS -> better agreement with TORBEAM rays at magnetic axis
  real(rkind)                        :: h_check
  real(rkind)                        :: h_glob = 1.d-4 !*
  real(rkind)                        :: eps_svec_max_length = 1.d-6 !* svec is shorter than the ray by eps_svec_max_length to avoid interpolation errors
  real(rkind)                        :: eps = 1.d-4!* for minimal substitution in S for Hamil=Stix
  real(rkind)                        :: time_smth = 1.d-3 !* smoothing time for Bt_vac at R0 = 1.65
  character(200)                     :: n_e_filename, T_e_filename ! ne and Te file
  character(200)                     :: Vessel_bd_filename! File that provides vessel boundary
  character(200)                     :: ray_launch_file ! File with launch geomertry
  character(200)                     :: working_dir ! Working directory
  character(200)                     :: data_name, data_secondary_name ! Filenames of primary and secondary Trad
  logical                            :: ray_init = .false. !* Controlls whether rays were computed once -> only relevant for IDA usage
  logical                            :: static_grid = .false. ! Prevents the grid for radiation transport to be recomputed
                                                              ! in the case of numerical instability
  logical                            :: double_check_splines = .false. !*
  logical                            :: output_all_ray_data = .true.
  ! The next two options are overwritten if ida is used
  logical                            :: stand_alone = .false. !* Whether the model is used as part of IDA or as a standalone program
  logical                            :: use_ida_spline_Te = .false., use_ida_spline_ne = .false.  !* whether IDA splines are used
  integer(ikind)                     :: eval, not_eval
  real(rkind)                        :: tau_ignore  = 1.d-8!* Threshold that decides whether the absorption of current point is considered
  real(rkind)                        :: tau_thick = 3.d0 !* If tau of two steps is larger than this, it is assumed that Trad = Te is a
                                                         ! good approximation for thermal plasmas. Neccessary to model strongly radiating plasmas.
  real(rkind)                        :: min_level_log_ne = 1.d-15
  integer(ikind)                     :: N_absz = 24, N_absz_large = 128
  logical                            :: new_IO = .false. ! Will be true if the new IO from ECRad_GUI will be used
  logical                            :: use_3D = .false.
end module mod_ECRad_types

