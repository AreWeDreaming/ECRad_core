module ece_types
use f90_kind
!use nag_spline_1d,             only: nag_spline_1d_comm_dp
implicit none
logical,                                  public :: static_ne_for_ece    ! T/F .. use predefined density parameters / use current parameters
! par_ece_ne is initialized in subroutine init_ECE_model_for_density.
! It is used for a constant density profile during Te optimization from ECE data
! for reasons of numerical stability
real(rkind), dimension(:),   allocatable, public :: par_ece_ne

real(rkind), dimension(:,:), allocatable :: dat_vec_ece, sig_vec_ece, sig_vec_ece_corr  ! (ece%N_ch*N_dat_profs)
real(rkind), dimension(:),   allocatable :: dat_model_ece        ! (ece%N_ch)
real(rkind), dimension(:,:), allocatable :: ddat_model_ece_dpar  ! (ece%N_ch,size(par))
real(rkind)    :: ln_likelihood_ece
logical        :: flag_ddat_model_ece_dpar
!type(nag_spline_1d_comm_dp) :: map_ece_spline

type ece_time_ch_type
  real(rkind) :: Te                                ! Trad-A
  real(rkind) :: Te_unc                            ! Trad-A uncertainty
  real(rkind) :: rhopol                            ! rhopol at position of cold resonance of 2X-mode
  ! Modifikationen von Sylvia Rathgeber
  real(rkind) :: rhop_iw                           ! rho_pol inner wall
  real(rkind) :: rhop_ow                           ! rho_pol outer wall
  real(rkind) :: rhop_res_1O                       ! rhopol at position of cold resonance of 1O-mode
  real(rkind) :: R_res_2X                          ! R      at position of cold resonance of 2X-mode
  real(rkind) :: z_res_2X                          ! z      at position of cold resonance of 2X-mode
  real(rkind) :: Te_1O                             ! Te at position of cold resonance of 1O-mode
  real(rkind) :: ne_1O                             ! ne at position of cold resonance of 1O-mode
  logical     :: low_field_side                    ! flag if channel is on low-field side wrt magn. axis
  logical     :: dat_use                           ! T/F .. use data if ok
end type ece_time_ch_type

type ece_time_type
  real(rkind)             :: time                  ! time
  type(ece_time_ch_type), dimension(:), allocatable :: ch      ! ch(1:N_ch_useful) for CEC
end type ece_time_type

type ece_timerz_ch_type
  real(rkind) :: ra                                ! R-A
  real(rkind) :: za                                ! z-A
  real(rkind) :: rhopol                            ! corresponding rho-poloidal 
  real(rkind) :: fpf                               ! corresponding Poloidal Flux    [Vs]
  logical     :: X3_contrib                        ! possibly an additional X3 contribution
end type ece_timerz_ch_type

type ece_timerz_type
  real(rkind) :: time                              ! rztime 
  real(rkind) :: R_mag                             ! R magnetic axis
  real(rkind) :: z_mag                             ! z magnetic axis
  real(rkind) :: PF_mag                            ! pol. flux magnetic axis
  real(rkind) :: R_sxp                             ! R separatrix
  real(rkind) :: z_sxp                             ! z separatrix
  real(rkind) :: PF_sxp                            ! pol. flux separatrix
  type(ece_timerz_ch_type), dimension(:), allocatable :: ch      ! ch(1:N_ch_useful) for ece
end type ece_timerz_type

type ece_ch_type
  real(rkind)    :: Trad_rel_unc                   ! relative uncertainty of Trad (SNR & calibration error)
  real(rkind)    :: freq                           ! [Hz] measurement frequency
  real(rkind)    :: dfreq                          ! [Hz] measurement bandwidth frequency
  real(rkind)    :: Btot                           ! [T]  magnetic field
  real(rkind)    :: snr                            ! []   total SNR evaluated considering the average signal to noise ratio of the channel and taking into account also the calibration error (the error in %: 100/SNR)
  integer(ikind) :: ifgroup                        ! [-1, 0-12]   -1 .. undefined
  integer(ikind) :: waveguide                      ! [-1, 0-12]   -1 .. undefined
  integer(ikind) :: availabl                       ! [0/1/-1]  no/yes/undefined data written in Trad-A, R-A, z-A
  real(rkind)    :: ne_cut_off                     ! [m^-3] cut-off density
    ! Modifikationen von Sylvia Rathgeber
  real(rkind)    :: ne_cut_off_1O                  ! [m^-3] cut-off density of 1st O-mode
  logical        :: cutoff_2X_ok
  logical        :: cutoff_1O_ok
end type ece_ch_type

type ece_ifgroup_type
  logical        :: used                           ! T/F .. use ifgroup
  integer(ikind) :: indx                           ! [-1, 0-12]   -1 .. undefined
  integer(ikind) :: N_ch                           ! number of channels
  real(rkind)    :: freq_mean                      ! [Hz] mean measurement frequency
end type ece_ifgroup_type

type ece_type
  !integer          :: shot
  character(8)     :: exp_ece                      ! <-  experiment for ECE measurement
  character(3)     :: diag_ece                     ! <-  diagnostic for ECE measurement
  integer(ikind)   :: ed_ece                       ! <-> edition of ECE measurement
  integer          :: N_ch                         ! # ECE spatial channels
  integer(ikind)   :: harmonic                     ! harmonic number [-1, 0-10]   -1 .. undefined
  real(rkind)      :: z_lens                       ! [m] vertical position of ECE lens
  integer,               dimension(:), allocatable :: ch_indx_masked   !  channel numbers not used (large uncertainty)
  integer(ikind),        dimension(:), allocatable :: ch_indx_output   !  channel numbers with output
  integer(ikind),        dimension(:), allocatable :: waveguid  ! (3)    [-1, 0-12]   -1 .. undefined
  type(ece_ch_type),     dimension(:), allocatable :: ch      ! channel (dim=N_ch_useful)
  type(ece_ifgroup_type),dimension(:), allocatable :: ifgroup ! 1-3
  type(ece_time_type),   dimension(:), allocatable :: time    ! time "time-A" for CEC
  type(ece_timerz_type), dimension(:), allocatable :: timerz  ! time "zrtime" for CEC
end type ece_type

end module ece_types

!********************************************************************
!********************************************************************
!********************************************************************

MODULE ece_global_params
use f90_kind
use ece_types,            only: ece_type
implicit none
type(ece_type)         :: ece

end module ece_global_params

