!****************************************************************
! MODULE mod_density
! MODULE mod_temperature
!****************************************************************
! Reduced version of libe_utils from ida
! these routines do nothing but provide empty functions that act as dummies for the usage within ida



MODULE mod_density
use f90_kind
implicit none

public :: make_density

CONTAINS

!*****************************************************************

subroutine make_density(par, par_scal, ne, x_eval, x_mode, dnedx)
use f90_kind
implicit none
real(rkind), dimension(:),   intent(in)  :: par
real(rkind), dimension(:),   intent(in)  :: par_scal
real(rkind), dimension(:),   intent(out) :: ne
real(rkind), dimension(:),   intent(in),  optional :: x_eval
integer(ikind),              intent(in),  optional :: x_mode   ! 1/2 .. x_LiB/x_rhopol
real(rkind), dimension(:),   intent(out), optional :: dnedx    ! derivative of ne wrt. spline abszissae
ne = 0.d0 ! to supress warnings
stop "sub make_density must not be called in the stand_alone version"
end subroutine make_density

!********************************************************************
!********************************************************************
!********************************************************************
end module mod_density

MODULE mod_temperature
use f90_kind
implicit none

public :: make_temperature

CONTAINS

!*****************************************************************

subroutine make_temperature(par, Te, x_eval, rp_min, dTedx)
use f90_kind
implicit none
real(rkind), dimension(:),   intent(in)  :: par
real(rkind), dimension(:),   intent(out) :: Te
real(rkind), dimension(:),   intent(in),  optional :: x_eval
real(rkind),                 intent(in),  optional :: rp_min   ! for extrapolation in [0, rp_min]
real(rkind), dimension(:),   intent(out), optional :: dTedx    ! derivative of Te wrt. spline abszissae
stop "sub make_temperature must not be called in the stand_alone version"
Te = 0.d0 ! to supress warnings
end subroutine make_temperature

!****************************************************************

END MODULE mod_temperature

