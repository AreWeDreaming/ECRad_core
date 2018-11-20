!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module ida_types
! module ida_global_params
! module opt_types
! module opt_global_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dummy for IDA integration

module ida_types
use f90_kind
implicit none

type ida_ece_type
  character(8)     :: expnam                          ! experiment for ECE measurement
  character(3)     :: diag                            ! diagnostic for ECE measurement
  integer(ikind)   :: edition                         ! edition of ECE measurement
  integer          :: N_ch                            ! # ECE spatial channels
  real(rkind)      :: data_unc_rel_std                ! ECE data uncertainty parameters: standard relative uncertainty (default: 0.07)
  real(rkind)      :: data_unc_offset                 ! ECE data uncertainty parameters: data_unc_offset (eV; default: 25)
  integer          :: rhopol_scal_te_fit              ! 0 / 1 .. const / fit parameter
  real(rkind)      :: rhopol_scal_te                  ! scale of ECE rhopol axis used for temperature profile (1 / <1 / >1 .. no scale / stretching / compression)
  integer          :: rhopol_scal_ne_fit              ! 0 / 1 .. const / fit parameter
  real(rkind)      :: rhopol_scal_ne                  ! scale of ECE rhopol axis used for density profile (1 / <1 / >1 .. no scale / stretching / compression)
  ! Modifikationen von Sylvia Rathgeber
  integer          :: flag_ne_rhop_scal               ! 0: no scaling; 1, val: const scaling with val; 2: variable scaling to Separatrix
  real(rkind)      :: ne_rhop_scal                    ! scale of density rhopol axis relative to ECE (1 / <1 / >1 .. no scale / stretching / compression)
  real(rkind)      :: ne_rhop_scal2                   ! scale of density rhopol axis relative to ECE (1 / <1 / >1 .. no scale / stretching / compression)
  integer(ikind)   :: lfs_mode                        ! use ECE data (0/1/2 .. all / lfs only / hfs only)
  integer(ikind)   :: forward_model                   ! ECE ECE forward model (0 Te=Trad (classic) / 1 band width (bw) / 2 bw+doppler / 3 bw+relativistic / 4 bw+doppler+relativistic)
  integer(ikind)   :: likelihood                      ! ECE likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: likelihood_cauchy_a0            ! ECE likelihood (only for likelihood=2 .. Cauchy a0)
  real(rkind)      :: emrel_theta_beg                 ! interval in theta when calculation ...
  real(rkind)      :: emrel_theta_end                 ! ... of emissivity switched to pure relativistic
  character(20)    :: eedf                            ! thermal electron energy distribution function: maxwell, maxwell_juettner
  character(9)     :: integration_scheme              ! Integration scheme of radiation transport equation (euler or rk4)
  integer(ikind)   :: N_ray                           ! Number of rays for antenna pattern
  integer(ikind)   :: ray_N_aufp                      ! Number of cells on each ray (default = 400)
  integer(ikind)   :: N_freq                          ! Number of sub-frequencies for each channel for antenna pattern (default = 3)
  logical          :: ray_tracing                     ! if 'true' different rays of each LOS vary in whole geometry
                                                      ! if 'false' different rays of each LOS vary only in theta
  logical          :: torbeam                         ! use torbeam code to calculate ECE LOS
  integer(ikind)   :: flag_reflec                     ! 0: single-pass, 1: intensity modified according to Figini 
                                                      ! with fixed reflection coefficient, 2: with paramterized reflection coefficient
  real(rkind)      :: reflec_X_mode                   ! fixed reflection coefficient for X-mode 0...1
  real(rkind)      :: reflec_O_mode                   ! fixed reflection coefficient for O-mode 0...1
  real(rkind)      :: opt_depth_threshold             ! scale of ECE optical depth for masking ECE data (ne*Te > ne_Te_crit = opt_depth_threshold * 3.d21 is used only)
  real(rkind)      :: cutoff_threshold                ! scale of ECE ne_cutoff for masking ECE data (ne < cutoff_threshold * ne_cutoff is used only)
  integer(ikind)   :: ece_1O_flag                     ! contribution from 1st O-mode 0:none, 1: add them to error of 2X, 2: Trad=Trad_2X+Trad_1O
  real(rkind)      :: ece_1O_rhop_min                 ! consider 1O contr outside of rhop_min
  logical          :: ece_3X_flag                     ! ECE contribution from 3st X-mode (default: F)
  real(rkind)      :: data_rhopol_max                 ! ECE data used oonly for rhopol <= rhopol_max (default = 1.1)
  logical          :: use_new_Rz_coord                ! T/F .. calc. / read from CEC   (R,z)-values
  logical          :: rztime_ok                       ! read (R,z)-values from CEC (to be checked if status is ok)
  logical          :: compare_rz_values               ! compare calculated with CEC (R,z)-values (only if CEC_rztime_ok=T)
  character(19)    :: btf_mode                        ! method to calculate Btf: "old_Btf", "BTFABB"
  real(rkind)      :: ch_fm_eval_opt_depth_threshold  ! use ECFM for an ECE channel if opt. depth is smaller than opt_depth_threshold_scale*opt_depth_crit
  real(rkind)      :: ch_fm_eval_rel_diff_Te_Trad     ! use ECFM for an ECE channel if rel. diff. between Te and Trad(ECFM) is larger than value given
  real(rkind)      :: rmc_time_step                   ! RMC temporal binning: time_step (distance between two time points)
  real(rkind)      :: rmc_time_smooth                 ! RMC temporal binning: time_smooth (bin size)
  integer(ikind)   :: ecrad_modes                     ! modes to be used in forward model: 1/2/3 = X/O/X+O  (default = 1)  (Note: for X or O only, the polarization filter is neglected)
  real(rkind)      :: ecrad_O2X_mode_conversion       ! fraction of O-mode converted to X-mode due to wall reflection (used only for ece_ecrad_modes = 3) (default = 0.0 ; not tested comprehensively)
  real(rkind)      :: ecrad_ds_large                  ! upper limit for step size in radiation transport (default = 2.5e-4, 2.5e-5)  (value is reduced automatically if numerical precision is insufficient)
  real(rkind)      :: ecrad_ds_small                  ! upper limit for step size in radiation transport (default = 2.5e-4, 2.5e-5)  (value is reduced automatically if numerical precision is insufficient)
  integer(ikind)   :: ecrad_max_points_svec           ! upper limit of # points on LOS for radiation transport (default = 20000)  (-> main memory of ECrad code!)
  real(rkind)      :: ecrad_ratio_for_third_harmonic  ! 3rd harmonic is considered for omega_c/omega > ece_ecrad_ratio_for_third_harmonic (default = 0.4 !> 1/3!)
  logical          :: ecrad_weak_rel                  ! raytracing using weakly relativistic corrections (default = F ; do not set to T if you are not sure what you are doing!) 
  logical          :: ecrad_Bt_ripple                 ! magnetic field ripple used (default = T)
  logical          :: ecrad_verbose                   ! addt'l output: birthplace distributions, warm resonances (default = F) (future: to be written in shotfile)
  real(rkind)      :: ecrad_R_shift                   ! shifts equilibrium for ECE only by R_shift [m] to the outside (default = 0.0)
  real(rkind)      :: ecrad_z_shift                   ! shifts equilibrium for ECE only by z_shift [m] upwards        (default = 0.0)
end type ida_ece_type

type faraday_set_los_pos_type
  real(rkind) :: r       ! [m] LOS vertical position
  real(rkind) :: z       ! [m] LOS vertical position
end type faraday_set_los_pos_type

type faraday_set_los_type
  real(rkind)  :: rb, re        ! [m] LOS radial   position begin and end 
  real(rkind)  :: zb, ze        ! [m] LOS vertical position begin and end 
  real(rkind)  :: scal          ! [rad m**2 / T] faraday_rotation_measure * wavelength**2
  real(rkind)  :: dx            ! [m] cell size
  character(8) :: nam           ! []    name of LOS
  real(rkind), dimension(2) :: dir  ! (r,z)-direction in the poloidal plane (normalized)
  type(faraday_set_los_pos_type), dimension(:), allocatable :: pos  ! position
end type faraday_set_los_type

type faraday_set_type
  logical        :: use_data
  integer(ikind) :: N_los                          ! # LOS to be evaluated (not all of them are fitted)
  integer(ikind) :: N_pos                          ! # discretization cells along LOS
  integer(ikind) :: N_los_fit                      ! # LOS to be fitted
  type(faraday_set_los_type), dimension(:), allocatable :: los      ! 
  integer(ikind),             dimension(:), allocatable :: los_fit  ! LOS indices to be fitted
end type faraday_set_type

type faraday_res_los_pos_type
  real(rkind) :: rhopol                  ! [] 
  real(rkind) :: ne                      ! [m^-3] electron density
  real(rkind) :: B_pol_parallel_to_los   ! [T]    poloidal magnetic field parallel to SOL
end type faraday_res_los_pos_type

type faraday_res_los_type
  real(rkind)  :: rot                     ! [rad] modelled faraday rotation
  type(faraday_res_los_pos_type), dimension(:), allocatable :: pos  ! position
end type faraday_res_los_type

type faraday_res_type
  type(faraday_res_los_type), dimension(:), allocatable :: los      ! 
  integer(ikind),             dimension(:), allocatable :: los_fit  ! LOS indices to be fitted
end type faraday_res_type

type ida_heb_line_type
  real(rkind)      :: wavelength                   ! [m] wavelength of emission line
  logical          :: fitted                       !     to be fitted
end type ida_heb_line_type

type ida_heb_type
  character(8)     :: expnam                       ! experiment for ECE measurement
  character(3)     :: diag                         ! diagnostic for ECE measurement
  integer(ikind)   :: edition                      ! edition of ECE measurement
  !integer          :: N_ch                         ! # ECE spatial channels
  !integer(ikind)   :: forward_model                ! 
  integer(ikind)   :: likelihood                   ! ECE likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: likelihood_cauchy_a0         ! ECE likelihood (only for likelihood=2 .. Cauchy a0)
  logical          :: mod_data_scale_param_fitted  ! T/F to be fitted or kept constant with value in input file
  real(rkind)      :: mod_data_scale_param         ! data scale parameter
  integer(ikind)   :: N_lines_all                  ! all calculated emission lines
  type(ida_heb_line_type), dimension(:), allocatable :: line
end type ida_heb_type

type ida_time_lib_ch_prof_type
  real(rkind) :: dat_amp   ! emission amplitude
  real(rkind) :: dat_unc   ! emission amplitude uncertainty
  real(rkind) :: dat_mes   ! emission amplitude measured (uncalibrated)
  real(rkind) :: dat_off   ! emission amplitude offset   (uncalibrated)
  real(rkind) :: residue   ! (ida%time(i)%ch(j)%dat_amp_lib(k)-ida%time(i)%ch(j)%dat_model_lib) / ida%time(i)%ch(j)%dat_unc_lib(k)
end type ida_time_lib_ch_prof_type

type ida_time_lib_ch_type
  real(rkind)             :: x     ! channel position
  type(ida_time_lib_ch_prof_type), dimension(:), allocatable :: prof    ! N_profs emission profiles
  real(rkind)             :: rhopol ! rhopol(x)
  real(rkind)             :: Rref   ! R coordinate of reflectometry 
  real(rkind)             :: zref   ! z coordinate of reflectometry 
  real(rkind)             :: dat_model_lib  ! simulated data identical for all N_profs emission profiles
end type ida_time_lib_ch_type

type ida_time_lib_ch_temp_type
  real(rkind)             :: temp     ! temperature       [eV]
  real(rkind)             :: xtemp    ! at position xtemp [cm]
end type ida_time_lib_ch_temp_type

type ida_time_lib_type
  integer(ikind)          :: N_dat_profs           ! # successively measured emission profiles to be used for estimating one density profile
  logical                 :: used                  ! evaluate LIB data for time point
  logical                 :: elmflag               ! flag if ELM activity
  real(rkind),                     dimension(:), allocatable :: vec_tim    ! time of N_profs emission profiles
  real(rkind),                     dimension(:), allocatable :: chi2_profs ! chi2 of N_profs LIB emission profiles
  type(ida_time_lib_ch_type),      dimension(:), allocatable :: ch       ! data
  type(ida_time_lib_ch_temp_type), dimension(:), allocatable :: ch_temp  ! temperature temp at position xtemp
  real(rkind)             :: ln_likelihood_up
  real(rkind)             :: ln_likelihood_mp
  real(rkind)             :: chi2_norm         ! chi2/N_dat of LiBeam data
  real(rkind)             :: fit_scale             ! scale of LiBe data
  real(rkind)             :: N2_z0                 ! N2(z=0)  boundary condition for 2p
end type ida_time_lib_type

type ida_time_dcn_los_prof_type
  real(rkind) :: dat_amp       ! emission amplitude             for N_profs emission profiles
  real(rkind) :: dat_unc       ! emission amplitude uncertainty for N_profs emission profiles
  real(rkind) :: residue       ! (ida%time(i)%dcn%los(j)%prof(k)%dat_amp-ida%time(i)%dcn%los(j)%dat_model) / ida%time(i)%dcn%los(j)%prof(k)%dat_unc
  real(rkind) :: fringe_fract  ! difference of dat_amp and dat_mod as a fraction of a fringe jump
end type ida_time_dcn_los_prof_type

type ida_time_dcn_los_type
  real(rkind)                            :: dat_model     ! simulated data   IDENTICAL for all N_profs emission profiles
  type(ida_time_dcn_los_prof_type), dimension(:), allocatable :: prof    ! N_profs emission profiles
end type ida_time_dcn_los_type

type ida_time_dcn_type
  integer(ikind)          :: N_dat_profs        ! # DCN data sets for joined fit
  integer(ikind)          :: dat_indx_lo        ! # DCN data lower indx for the present time interval
  integer(ikind)          :: dat_indx_up        ! # DCN data upper indx for the present time interval
  integer(ikind)          :: dat_indx_coo_lo    ! # COO data lower indx for the present time interval
  integer(ikind)          :: dat_indx_coo_up    ! # COO data upper indx for the present time interval
  real(rkind)             :: ln_likelihood
  real(rkind)             :: chi2_norm          ! chi2/N_dat of DCN data
  real(rkind),                 dimension(:), allocatable :: vec_tim     ! time of N_profs DCN data
  real(rkind),                 dimension(:), allocatable :: chi2_profs  ! chi2 of N_profs DCN 
  type(ida_time_dcn_los_type), dimension(:), allocatable :: los         ! all for DCN interferometer: line of sight (los)
end type ida_time_dcn_type

type ida_time_ece_ch_prof_type
  real(rkind)                            :: dat_amp       ! corrected data uncertainty
  real(rkind)                            :: dat_unc_corr  ! corrected data uncertainty
  real(rkind)                            :: residue       ! (ida%time(i)%ece%ch(j)%prof(k)%dat_amp-ida%time(i)%ece%ch(j)%dat_model) / ida%time(i)%ece%ch(j)%prof(k)%dat_unc
  logical                                :: dat_use       ! T/F .. use data if ok: density cutoff, opt.depth and data<>0
end type ida_time_ece_ch_prof_type

type ida_time_ece_ch_type
  type(ida_time_ece_ch_prof_type), dimension(:), allocatable :: prof  ! data to be analyzed jointly
  real(rkind)                            :: dat_model     ! simulated data   IDENTICAL for all N_profs emission profiles
  real(rkind)                            :: rhopol        ! mean rhopol
  real(rkind)                            :: rhotor        ! corresponding rho-toroidal
  real(rkind)                            :: R_res_2X      ! R at position of cold resonance of 2X-mode
  real(rkind)                            :: z_res_2X      ! z at position of cold resonance of 2X-mode
  logical                                :: fm_flag       ! T = use forward modelling / F = no fm (Trad=Te)
end type ida_time_ece_ch_type

type ida_time_ece_type
  integer(ikind)                                        :: N_dat_profs       ! # Te-profiles for joined fit
  integer(ikind)                                        :: N_ch              ! # ECE channels
  integer(ikind),             dimension(:), allocatable :: tim_ind           ! time indices of ECE emission profiles
  real(rkind),                dimension(:), allocatable :: chi2_profs        ! chi2 of N_profs ECE
  ! Modifikationen von Sylvia Rathgeber
  real(rkind)                                           :: rhop_iw       !  position of inner wall
  real(rkind)                                           :: rhop_ow       !  position of outer wall
  !real(rkind)                                           :: rp_min            ! minimal ECE rhopol for useful channels to extrapolate within rhopol = (0, rp_min)
  type(ida_time_ece_ch_type), dimension(:), allocatable :: ch    ! channel
  real(rkind)             :: chi2_norm         ! chi2/N_dat of ECE data
end type ida_time_ece_type

type ida_time_ts_set_poly_ch_prof_type
  real(rkind)                            :: dat_amp       ! 
  real(rkind)                            :: dat_unc       ! 
  real(rkind)                            :: rhopol        ! 
  real(rkind)                            :: rhopol_shift  ! 
  real(rkind)                            :: dat_mod       !
  real(rkind)                            :: residue       ! (ida%time(i)%ts%set(l)%poly(j)%prof(k)%dat_amp-ida%time(i)%ts%set(l)%poly(j)%dat_model) / ida%time(i)%ts%set(l)%poly(j)%prof(k)%dat_unc
  logical                                :: dat_use       ! T/F .. use data if ok
end type ida_time_ts_set_poly_ch_prof_type

type ida_time_ts_set_poly_ch_type
  real(rkind)                            :: dat_mod       ! simulated data   IDENTICAL for all N_profs data sets profiles
  type(ida_time_ts_set_poly_ch_prof_type),  dimension(:), allocatable :: prof         ! (N_poly) group delay
end type ida_time_ts_set_poly_ch_type

type ida_time_ts_set_poly_type
  real(rkind)                            :: rhopol        ! 
  real(rkind)                            :: rhopol_shift  !  shifted/scaled rhopol
  type(ida_time_ts_set_poly_ch_type),  dimension(:), allocatable :: ch       ! (ne, Te)
end type ida_time_ts_set_poly_type

type ida_time_ts_set_type
  integer(ikind)                                         :: N_dat_profs  ! # (ne, Te) pairs for joined fit
  integer(ikind),              dimension(:), allocatable :: tim_ind      ! time indices of REF group-delay profiles
  real(rkind),                 dimension(:), allocatable :: vec_tim      ! time of N_profs emission profiles
  real(rkind),                 dimension(:), allocatable :: chi2_profs   ! chi2 of N_profs REF
  type(ida_time_ts_set_poly_type),  dimension(:), allocatable :: poly    ! (N_poly) group delay
end type ida_time_ts_set_type

type ida_time_ts_type
  type(ida_time_ts_set_type),  dimension(:), allocatable :: set       ! (2) set: core or edge
end type ida_time_ts_type

type ida_time_ref_ch_prof_type
  real(rkind)                            :: gd            ! group delay
  real(rkind)                            :: gd_unc        ! group delay uncertainty
  real(rkind)                            :: residue       ! (ida%time(i)%ref%ch(j)%prof(k)%dat_amp-ida%time(i)%ref%ch(j)%dat_model) / ida%time(i)%ref%ch(j)%prof(k)%dat_unc
  logical                                :: dat_use       ! T/F .. use data if ok: density cutoff, opt.depth and data<>0
end type ida_time_ref_ch_prof_type

type ida_time_ref_ch_type
  type(ida_time_ref_ch_prof_type), dimension(:), allocatable :: prof  ! (N_dat_prof) group delay data to be analyzed jointly
  real(rkind)                            :: dat_model     ! simulated data   IDENTICAL for all N_profs group-delay profiles
  real(rkind)                            :: freq          ! probing frequency
end type ida_time_ref_ch_type

type ida_time_ref_los_type
  real(rkind)                                           :: Rreferenz    ! reference major radius from where group delay calculation starts
end type ida_time_ref_los_type

type ida_time_ref_type
  integer(ikind)                                         :: N_dat_profs  ! # group-delay profiles for joined fit
  integer(ikind)                                         :: Nf_lfs       ! # frequency channels
  type(ida_time_ref_los_type), dimension(:), allocatable :: los          ! (N_ref_los)   line of sight
  integer(ikind),              dimension(:), allocatable :: tim_ind      ! time indices of REF group-delay profiles
  real(rkind),                 dimension(:), allocatable :: chi2_profs   ! chi2 of N_profs REF
  type(ida_time_ref_ch_type),  dimension(:), allocatable :: ch_lfs       ! (Nf) group delay
  real(rkind),                 dimension(:), allocatable :: vec_tim      ! time of N_profs emission profiles
  real(rkind)             :: chi2_norm_lfs         ! chi2/N_dat of REF data
end type ida_time_ref_type

type ida_time_bes_type
end type ida_time_bes_type

type ida_time_heb_type
  !integer(ikind)          :: N_dat_profs        ! # HEB data sets for joined fit
  !integer(ikind)          :: dat_indx_lo        ! # HEB data lower indx for the present time interval
  !integer(ikind)          :: dat_indx_up        ! # HEB data upper indx for the present time interval
  real(rkind)             :: ln_likelihood
  real(rkind)             :: chi2_norm          ! chi2/N_dat of HEB data
  !real(rkind),                 dimension(:), allocatable :: vec_tim     ! time of N_profs HEB data
  !real(rkind),                 dimension(:), allocatable :: chi2_profs  ! chi2 of N_profs HEB 
  !type(ida_time_heb_los_type), dimension(:), allocatable :: los         ! all for HEB interferometer: line of sight (los)
end type ida_time_heb_type

type ida_time_ne_type
  real(rkind)             :: x         ! channel position
  real(rkind)             :: amp       ! emission amplitude
  real(rkind)             :: unc       ! emission amplitude uncertainty
  real(rkind)             :: unc2      ! emission amplitude uncertainty
  real(rkind)             :: amp_lo    ! emission amplitude lower "limit"
  real(rkind)             :: amp_up    ! emission amplitude upper "limit"
end type ida_time_ne_type

type ida_time_rp_type
  real(rkind)             :: rhopol    ! rho-poloidal
  real(rkind)             :: rhotor    ! rho-toroidal
  real(rkind)             :: psipol    ! [Vs] poloidal flux
  real(rkind)             :: psitor    ! [Vs] toroidal flux
  real(rkind)             :: r_out     ! R in outer midplane
  real(rkind)             :: z_out     ! z in outer midplane
  real(rkind)             :: drp_dr    ! derivative of rhopol w.r.t. r_out in outer midplane
  real(rkind)             :: ne        ! density 
  real(rkind)             :: ne_unc    ! density uncertainty
  real(rkind)             :: ne_lo     ! density lower "limit"
  real(rkind)             :: ne_up     ! density upper "limit"
  real(rkind)             :: dne_drp   ! derivative of density wrt. rhopol 
  real(rkind)             :: dne_dr    ! derivative of density wrt. r
  real(rkind)             :: Te        ! temperature
  real(rkind)             :: Te_unc    ! temperature uncertainty
  real(rkind)             :: Te_lo     ! temperature lower "limit"
  real(rkind)             :: Te_up     ! temperature upper "limit"
  real(rkind)             :: dTe_drp   ! derivative of temperature wrt. rhopol 
  real(rkind)             :: dTe_dr    ! derivative of temperature wrt. r
  real(rkind)             :: pe        ! pressure
  real(rkind)             :: pe_unc    ! pressure uncertainty
  real(rkind)             :: dpe_drp   ! derivative of pressure wrt. rhopol 
  real(rkind)             :: dpe_dr    ! derivative of pressure wrt. r
  real(rkind)             :: opt_depth ! optical depth
end type ida_time_rp_type

type ida_time_type
  real(rkind)             :: time                  ! time
  real(rkind)             :: time_beg, time_end    ! time interval for analysis [s]
  logical                 :: used                  ! evaluate time point
  real(rkind)             :: deviance              ! chi2 of data fit (?)
  real(rkind)             :: prior_scale_lib, prior_monot_ne, prior_curv_ne
  real(rkind)             :: prior_monot_Te, prior_curv_Te, prior_Te_sep
  real(rkind)             :: prior_Te_sol          ! Te value in SOL  (ecfm)
  real(rkind)             :: prior_Trad_sol          ! Te value in SOL  (ecfm)
  real(rkind)             :: prior_Te_wall         ! Te value at wall (ecfm)
  real(rkind)             :: prior_monot_pe
  real(rkind)             :: prior_time_smooth
  real(rkind)             :: prior_ne_values       ! ne values at arbitrary points
  real(rkind)             :: prior_Te_values       ! Te values at arbitrary points
  ! Modifikationen von Sylvia Rathgeber
  real(rkind)             :: ece_ne_rhop_scal                        ! scale of ECE rhopol axis relative to LiBeam (1 / <1 / >1 .. no scale / stretching / compression) stretching means that ECE separatrix is at larger rhopol values
  real(rkind)             :: ece_reflec                              ! variable reflection coefficient 0...1 (ecfm)
  real(rkind)             :: H2farrot              ! Faraday rotation angle of H-2 interferometry LOS
  real(rkind)             :: te_rp_min             ! minimum rhopol value of all data determining Te -> will be updated with each diagnostic used providing Te data
  type(ida_time_ne_type),      dimension(:), allocatable :: ne       ! density
  type(ida_time_ne_type),      dimension(:), allocatable :: ne_hres  ! density high resolution
  type(ida_time_rp_type),      dimension(:), allocatable :: rp       ! diff. values as a function of rho-poloidal
  real(rkind),                 dimension(:), allocatable :: ne_knot_pos_rhopol  ! rhopol knot positions (only for ne_interpolation_type=2/3/4)
  real(rkind),                 dimension(:), allocatable :: par      ! parameter
  real(rkind),                 dimension(:), allocatable :: par_scal ! scale of parameter
  type(ida_time_lib_type)                                :: lib      ! all for LIB
  type(ida_time_dcn_type)                                :: dcn      ! all for DCN interferometer
  type(ida_time_ece_type)                                :: ece      ! all for ECE
  type(ida_time_ts_type)                                 :: ts      ! all for Thomson scattering
  type(ida_time_ref_type)                                :: ref      ! all for Reflectometry
  type(ida_time_bes_type)                                :: bes      ! all for BES
  type(ida_time_heb_type)                                :: heb      ! all for HEB
  type(faraday_res_type)                                 :: faraday_res         ! results
end type ida_time_type

type ida_time_interval_skipped_type
  real(rkind)             :: tbeg, tend    ! time interval not used [s]
end type ida_time_interval_skipped_type

type ida_prior_ne
  real(rkind)             :: amp           ! density amplitude   [m^-3]
  real(rkind)             :: unc           ! density uncertainty [m^-3]
  real(rkind)             :: rhopol        ! position
end type ida_prior_ne

type ida_prior_Te
  real(rkind)             :: amp           ! temperature amplitude   [eV]
  real(rkind)             :: unc           ! temperature uncertainty [eV]
  real(rkind)             :: rhopol        ! position
end type ida_prior_Te

type randomize_data_type
  logical          :: used                 !     randomize diagnostic data (set T iff data of one time point should be analized multiple times with noise added)
  real(rkind)      :: time                 ! [s] time
end type randomize_data_type

type ida_type
  character(99999) :: str_ida_inp                      ! unformated string with the content of ida.inp 
  integer          :: shot
  integer(ikind)   :: time_indx_begin, time_indx_end, time_indx_step ! time indices to be analyzed and index_step (actual time depends on tbeg, tend, and timstep)
                                                       ! (0, ?) means time points from tbeg / (?, 0) means time points until tend
  logical          :: lib_only                         ! only LIB -> LIN
  logical          :: write_shot_file_flag             ! write in shot file (T/F); file depends on "lib_only"
  character(8)     :: experiment_write                 ! Experiment to be written: RRF, AUGD, AUGE, LIBE, ...
  character(8)     :: program_version                  ! program version
  character(99)    :: fs_tmp                           ! file system for temporary storage: /afs/ipp/u/ (afs) ; /ptmp1/ (on TOK cluster)
  logical          :: use_lib                          ! F/T .. LiBe diagnostic
  logical          :: use_dcn                          ! F/T .. DCN interferometry diagnostic
  logical          :: use_ece                          ! F/T .. ECE diagnostic
  logical          :: use_ts                           ! F/T .. Thomson scattering diagnostic
  logical          :: use_ref                          ! F/T .. Reflectometry diagnostic
  logical          :: use_bes                          ! F/T .. beam emission spectroscopy
  logical          :: use_heb                          ! F/T .. thermal helium beam emission spectroscopy
  logical          :: use_lib_mandatory                ! F/T .. LiBe diagnostic data mandatory
  logical          :: use_dcn_mandatory                ! F/T .. DCN interferometry diagnostic data mandatory
  logical          :: use_ece_mandatory                ! F/T .. ECE diagnostic data mandatory
  logical          :: use_ts_mandatory                 ! F/T .. Thomson scattering diagnostic data mandatory
  logical          :: use_ref_mandatory                ! F/T .. Reflectometry diagnostic data mandatory
  logical          :: use_bes_mandatory                ! F/T .. beam emission spectroscopy data mandatory
  logical          :: use_heb_mandatory                ! F/T .. thermal helium beam emission spectroscopy data mandatory
  logical          :: estimate_ne                      ! F/T .. estimate ne-profile
  logical          :: estimate_Te                      ! F/T .. estimate Te-profile
  character(8)     :: exp_eq                           ! <-  experiment for equilibrium data 
  character(3)     :: diag_eq                          ! <-  diagnostic for equilibrium data
  integer          :: ed_eq                            ! <-> edition for equilibrium data
  integer          :: interpol_eq                      ! <-> interpolation for equilibrium data (1/2 .. next/linear)
  character(8)     :: exp_ip                           ! <-  experiment for current measurement
  character(3)     :: diag_ip                          ! <-  diagnostic for current measurement
  integer          :: ed_ip                            ! -> edition of current measurement
  integer          :: ed_ip_shot                       ! <- edition of current measurement
  integer          :: ed_ip_cal_shot                   ! <- edition of current measurement
  character*8      :: name_ip                          ! <-  signal name for current measurement
    ! Modifikationen von Sylvia Rathgeber
  character(8)     :: elm_exper                        ! experiment of ELM shotfile
  character(3)     :: elm_diag                         ! diagnostic of ELM shotfile
  integer(ikind)   :: elm_edit                         ! edition    of ELM shotfile 
  real(rkind)      :: dtelm                            ! min time to ELM [s]
  integer          :: sig_max_ip                       ! maximum length of Ip data array
  real             :: ip_limit                         ! minimum Ip in suitable discharge
  logical          :: cal_ip                           ! flag to read calibrated signal
  real(rkind)      :: tbeg                             ! <-> first time
  real(rkind)      :: tend                             ! <-> last time
  integer          :: N_prof_rp                        ! number of rho_pol values (-> density values)
  logical          :: ne_rigorous_monot_flag           ! ne = int exp(spline)
  logical          :: prior_scale_flag, prior_ne_monot_flag, prior_ne_curv_flag,  &
                      prior_Te_monot_flag, prior_Te_curv_flag, prior_Te_sep_flag, &
                      prior_time_smooth_flag, prior_pe_monot_flag
  logical          :: prior_Te_sol_flag                ! flag for prior_Te_sol (ecfm)
  logical          :: prior_Trad_sol_flag              ! flag for prior_Trad_sol (ecfm Te<Trad)
  real(rkind)      :: prior_ne_monot_scale             ! scale factor of density monotonicity penalization
  real(rkind)      :: prior_ne_curv_scale              ! scale factor of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp1                ! rhopol left of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp1_decay          ! rhopol exponential decay length of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp1_scal           ! scale factor of decay component of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp2                ! rhopol right of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp2_decay          ! rhopol exponential decay length of density curvature penalization
  real(rkind)      :: prior_ne_curv_rp2_scal           ! scale factor of decay component of density curvature penalization
  real(rkind)      :: prior_Te_monot_scale             ! scale factor of temperature monotonicity penalization
  real(rkind)      :: prior_Te_curv_scale              ! scale factor of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp1                ! rhopol left of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp1_decay          ! rhopol exponential decay length of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp1_scal           ! scale factor of decay component of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp2                ! rhopol right of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp2_decay          ! rhopol exponential decay length of temperature curvature penalization
  real(rkind)      :: prior_Te_curv_rp2_scal           ! scale factor of decay component of temperature curvature penalization
  integer(ikind)   :: prior_Te_sep_pdf                 ! pdf of likelihood (1/2 .. Gauss/Cauchy)
  character(1)     :: prior_Te_sep_constr              ! pdf constr: u/l/b .. upper/lower/both constraints  (default: u)
  real(rkind)      :: prior_Te_sep_mean                ! mean Te value at separatrix [eV]
  real(rkind)      :: prior_Te_sep_std                 ! "standard deviation" of Te value at separatrix [eV]
  integer(ikind)   :: prior_Te_sol_pdf                 ! pdf of likelihood for prior_Te_sol (ecfm) (1/2 .. Gauss/Cauchy)
  integer(ikind)   :: prior_Trad_sol_pdf               ! pdf of likelihood for prior_Trad_sol (ecfm) (1/2 .. Gauss/Cauchy)
  real(rkind)      :: prior_Te_sol_rhop                ! rhop for SOL prior (ecfm)
  real(rkind)      :: prior_Trad_sol_rhop              ! rhop for SOL prior (ecfm)
  real(rkind)      :: prior_Te_sol_mean                ! mean Te value outside separatrix [eV] (ecfm)
  real(rkind)      :: prior_Te_sol_std                 ! "standard deviation" of Te value outside separatrix [eV]
  real(rkind)      :: prior_Trad_sol_std               ! "standard deviation" of Te value outside separatrix [eV]
  logical          :: prior_Te_wall_flag               ! flag for prior_Te_wall (ecfm)
  integer(ikind)   :: prior_Te_wall_pdf                ! pdf of likelihood (1/2 .. Gauss/Cauchy)
  real(rkind)      :: prior_Te_wall_mean               ! mean Te value outside wall [eV]
  real(rkind)      :: prior_Te_wall_std                ! "standard deviation" of Te value outside wall [eV]
  real(rkind)      :: prior_pe_monot_scale             ! scale factor of pressure monotonicity penalization
  integer(ikind)   :: time_base_type                   ! time base for ida structure: 1..LIB structure; 2..continuous
  integer          :: isw_mode                         ! 1/2/3 .. optimization / MCMC / new optimization
  integer          :: optimizer_type                   ! 1/2 .. NAG/MINUIT optimization
  integer          :: ne_interpolation_type            ! 1/2/3/4 .. pointwise/cubic B-spline(NAG)/piecewise cubic Hermite function/cubic spline(NumRecip)
  integer(ikind)   :: ne_N_spline_knots_beam           ! # interpolation knots for ne(Li beam coord.) (only for ne_interpolation_type=2/3/4)
  integer(ikind)   :: ne_N_spline_knots_rhopol         ! # interpolation knots for ne(rhopol) (only for ne_interpolation_type=2/3/4)
  integer(ikind)   :: Te_par_init                      ! initialize Te spline amplitude parameters: 0/1 .. parameter file / ECE radiation temperature (= ECE measurements)
  integer          :: Te_interpolation_type            ! 1 .. cubic spline(NumRecip)
  integer(ikind)   :: Te_N_spline_knots_rhopol         ! # interpolation knots for Te(rhopol)
  real(rkind)      :: Te_sys_unc_abs                   ! Te minimal absolute systematic uncertainty
  real(rkind)      :: Te_sys_unc_rel                   ! Te minimal relative systematic uncertainty
  integer          :: ne_spline_knot_selection_lib_beam_coord  ! ne spline knot position selection on lithium beam coordinate
  integer          :: isw_ne_err_calc                  ! evaluate error on density     (0/1 .. no/chi2)
  integer          :: isw_Te_err_calc                  ! evaluate error on temperature (0/1 .. no/chi2)
  integer          :: N_ne_hres                        ! # high resolution density values
  real(rkind)      :: rhopol_max_spline_knot           ! maximum rhopol at outermost spline knot to be analyzed
  real(rkind)      :: lib_rate_coef_temp_lo            ! lower temperature for LIB rate coefficients [eV]
  real(rkind)      :: lib_rate_coef_temp_up            ! upper temperature for LIB rate coefficients; compared to TEMax [eV]
  logical          :: lib_use_upper_optic              ! use data from upper optic (LIB or LIZ)
  logical          :: lib_use_lower_optic              ! use data from lower optic (LIZ only)
  integer(ikind)   :: lib_likelihood                   ! LIB likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: lib_likelihood_cauchy_a0         ! LIB likelihood (only for dcn_likelihood=2 .. Cauchy a0)
  integer(ikind)   :: dcn_likelihood                   ! DCN likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: dcn_likelihood_cauchy_a0         ! DCN likelihood (only for dcn_likelihood=2 .. Cauchy a0)
  real(rkind)      :: dcn_LOS_z_shift                  ! DCN vertical shift of lines of sight in [m] (default should be 0.00 m)
  real(rkind)      :: btf_corr_fact_ext                ! external Btf correction factor (to shift the ECE channels radially: R_new = R * btf_corr_fact_ext; dR = R*(btf_corr_fact_ext-1))
  integer(ikind)   :: ts_likelihood                    ! TS likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: ts_likelihood_cauchy_a0          ! TS likelihood (only for dcn_likelihood=2 .. Cauchy a0)
  logical          :: ts_core_data_used                ! TS core data set used
  logical          :: ts_core_ne_used                  ! TS core density used
  logical          :: ts_core_te_used                  ! TS core temperature used
  logical          :: ts_edge_data_used                ! TS edge data set used
  logical          :: ts_edge_ne_used                  ! TS edge density used
  logical          :: ts_edge_te_used                  ! TS edge temperature used
  integer          :: ts_rhop_scal_fit                 ! 0 / 1 .. const / fit parameter
  real(rkind)      :: ts_rhop_scal                     ! scale of TS rhopol axis relative to LiBeam (1 / <1 / >1 .. no scale / stretching / compression) stretching means that ECE separatrix is at larger rhopol values
  logical          :: ts_ne_scal_fit                   ! T/F .. fit/keep constant
  real(rkind)      :: ts_ne_scal                       ! scale of TS ne profile
  real(rkind)      :: ts_time_interval_extention       ! extend time interval [s] for TS data to be used: [time_beg - ts_time_interval_extention; time_end + ts_time_interval_extention]
  real(rkind)      :: ts_R_core_shift_common           ! shift R-coordinates of core TS channels [m] (+dR=outward; +dz=upward)
  real(rkind)      :: ts_z_core_shift_common           ! shift z-coordinates of core TS channels [m] (+dR=outward; +dz=upward)
  real(rkind)      :: ts_R_edge_shift_common           ! shift R-coordinates of edge TS channels [m] (+dR=outward; +dz=upward)
  real(rkind)      :: ts_z_edge_shift_common           ! shift z-coordinates of edge TS channels [m] (+dR=outward; +dz=upward)
  integer(ikind)   :: ref_likelihood                   ! REF likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: ref_likelihood_cauchy_a0         ! REF likelihood (only for ref=2 .. Cauchy a0)
  integer(ikind)   :: bes_likelihood                   ! BES likelihood (1/2 .. Gaussian / Cauchy (outlier))
  real(rkind)      :: bes_likelihood_cauchy_a0         ! BES likelihood (only for dcn_likelihood=2 .. Cauchy a0)
  logical          :: bes_mod_data_scale_params_fitted ! T/F to be fitted or kept constant with value in input file
  integer(ikind)   :: bes_mod_data_scale_params_number ! # BES energy component scale parameters
  real(rkind), dimension(:), allocatable :: bes_mod_data_scale_params ! BES modelled data scale factors
  logical          :: bes_zeff_profile_params_fitted   ! BES Zeff profile parameters to be fitted or kept constant with value in input file
  character(99)    :: bes_zeff_profile_params_type     ! BES Zeff profile type "linear" 
  integer(ikind)   :: bes_zeff_profile_params_number   ! # BES Zeff profile parameters
  real(rkind), dimension(:), allocatable :: bes_zeff_profile_params   ! BES Zeff profile parameters
  real(rkind)      :: bes_rhop_scal_ne                 ! scale of BES rhopol axis for electron density relative to LiBeam (1 / <1 / >1 .. no scale / stretching / compression); stretching means that separatrix is at larger rhopol values
  integer          :: mcmc_res_start                   ! res_start.dat: 1/2 .. use as is / calculate from params_all.res
  real(rkind), dimension(:), allocatable :: Te_knot_pos_rhopol   ! Te(rhopol) knot positions
  real(rkind), dimension(:), allocatable :: ne_knot_pos_lib_beam ! knot positions on beam coord. (only for ne_interpolation_type=2/3/4)
  real(rkind), dimension(:), allocatable :: time_eqi      ! EQI/EQH times
  real(rkind), dimension(:), allocatable :: time_eqi_thr  ! EQI/EQH time intervall boundaries
  type(ida_time_interval_skipped_type), dimension(:), allocatable :: time_interval_skipped
  type(ida_time_interval_skipped_type), dimension(:), allocatable :: time_interval_not_written
  type(ida_prior_ne),  dimension(:), allocatable :: prior_ne
  type(ida_prior_Te),  dimension(:), allocatable :: prior_Te
  type(ida_time_type), dimension(:), allocatable :: time
  type(ida_ece_type)        :: ece
  type(ida_heb_type)        :: heb
  type(faraday_set_type)    :: faraday_set
  type(randomize_data_type) :: randomize_data
end type ida_type

end module ida_types

!********************************************************************
!********************************************************************
!********************************************************************

MODULE ida_global_params
use f90_kind
use ida_types,            only: ida_type
implicit none
type(ida_type)         :: ida
character(90)          :: fn_home

end module ida_global_params

!********************************************************************
!********************************************************************
!********************************************************************

module opt_types
! every information needed for optimization
use f90_kind
implicit none

type opt_ece_ch_prof_type
  logical                                :: dat_use              ! T/F .. use data if ok: density cutoff, opt.depth and data<>0
  real(rkind)                            :: Te, Te_unc           ! 
  logical                                :: low_field_side       ! T/F .. lfs/hfs
end type opt_ece_ch_prof_type

type opt_ece_ch_type
  type(opt_ece_ch_prof_type), dimension(:), allocatable :: prof  ! data to be analyzed jointly
  real(rkind)                            :: rhopol               ! mean rhopol
  logical                                :: dat_use              ! T/F .. use data if ok: density cutoff, opt.depth and data<>0
end type opt_ece_ch_type

type opt_ece_type
  integer(ikind)                                   :: N_ch              ! # ECE channels
  integer(ikind)                                   :: N_dat_profs       ! # Te-profiles for joined fit
  real(rkind)                                      :: rp_min            ! minimal ECE rhopol for useful channels to extrapolate within rhopol = (0, rp_min)
  type(opt_ece_ch_type), dimension(:), allocatable :: ch    ! channel
end type opt_ece_type

type opt_type
  !type(opt_lib_type)                                :: lib      ! all for LIB
  !type(opt_dcn_type)                                :: dcn      ! all for DCN interferometer
  type(opt_ece_type)                                :: ece      ! all for ECE
end type opt_type

end module opt_types

!********************************************************************
!********************************************************************
!********************************************************************

MODULE opt_global_params
use f90_kind
use opt_types,            only: opt_type
implicit none
type(opt_type)         :: opt

end module opt_global_params

