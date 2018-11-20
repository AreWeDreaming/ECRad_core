module eqi_types
use f90_kind
implicit none

type eqi_time_type
  real(rkind) :: time                              ! time
  real(rkind) :: R_mag                             ! time
  real(rkind) :: z_mag                             ! z magnetic axis
  real(rkind) :: PF_mag                            ! pol. flux magnetic axis
  real(rkind) :: R_sxp                             ! R separatrix
  real(rkind) :: z_sxp                             ! z separatrix
  real(rkind) :: PF_sxp                            ! pol. flux separatrix
  real(rkind) :: Br_mag                            ! radial   magnetic field at magnetic axis
  real(rkind) :: Bz_mag                            ! vertical magnetic field at magnetic axis
  real(rkind) :: Bt_mag                            ! toroidal magnetic field at magnetic axis
  real(rkind) :: Btot_mag                          ! total    magnetic field at magnetic axis
  real(rkind) :: Btot_mag_ece                      ! total    magnetic field at magnetic axis incl. ECE field ripple correction
  real(rkind) :: Btf0_eq                           ! vacuum   magnetic field from kkrzBrzt
  real(rkind) :: Btf0                              ! vacuum   magnetic field from MBI:BTFABB if btf_mode == "BTFABB" else Btf0 = Btf0_eq
end type eqi_time_type

type eqi_type
  integer(ikind)   :: shot
  character(8)     :: exper                        ! <-  experiment for magnetic eqilibrium
  character(3)     :: diag                         ! <-  diagnostic for magnetic eqilibrium
  integer(ikind)   :: edit                         ! <-> edition    of  magnetic eqilibrium
  type(eqi_time_type),   dimension(:), allocatable :: time    ! time
end type eqi_type

end module eqi_types

!********************************************************************
!********************************************************************
!********************************************************************

MODULE eqi_global_params
use f90_kind
use eqi_types,            only: eqi_type
implicit none
type(eqi_type)         :: eqi

end module eqi_global_params

