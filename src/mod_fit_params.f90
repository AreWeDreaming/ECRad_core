! dummy module required by mod_ecfm_refr

MODULE mod_fit_params
use f90_kind
implicit none


public ::  map_par_to_ece_reflec
           

CONTAINS

!****************************************************************

subroutine map_par_to_ece_reflec(par, ece_reflec, par_scal)
use f90_kind
implicit none
real(rkind), dimension(:), intent(in)            :: par            ! parameter to be fitted
real(rkind),               intent(out)           :: ece_reflec     ! ECE wall reflection parameter
real(rkind), dimension(:), intent(in),  optional :: par_scal       ! scale of parameter to be fitted
ece_reflec = -1.d0
print*, "map_par_to_ece_reflec contained in mod_fit_params is a dummy routine to be replaced by IDA - do not use in standalone ECRAD"
stop "sub map_par_to_ece_reflec"
end subroutine map_par_to_ece_reflec

subroutine map_par_to_ece_rp_scal_te(par, ece_rp_scal_te, par_scal)
use f90_kind
implicit none
real(rkind), dimension(:), intent(in)            :: par            ! parameter to be fitted
real(rkind),               intent(out)           :: ece_rp_scal_te ! ECE index of rhopol scale of temperature profile
real(rkind), dimension(:), intent(in),  optional :: par_scal       ! scale of parameter to be fitted
ece_rp_scal_te = -1.d0
print*, "map_par_to_ece_rp_scal_te contained in mod_fit_params is a dummy routine to be replaced by IDA - do not use in standalone ECRAD"
stop "sub map_par_to_ece_rp_scal_te"
end subroutine map_par_to_ece_rp_scal_te

!****************************************************************

subroutine map_par_to_ece_rp_scal_ne(par, ece_rp_scal_ne, par_scal)
use f90_kind
implicit none
real(rkind), dimension(:), intent(in)            :: par            ! parameter to be fitted
real(rkind),               intent(out)           :: ece_rp_scal_ne ! ECE index of rhopol scale of density profile
real(rkind), dimension(:), intent(in),  optional :: par_scal       ! scale of parameter to be fitted
ece_rp_scal_ne = -1.d0
print*, "map_par_to_ece_rp_scal_ne contained in mod_fit_params is a dummy routine to be replaced by IDA - do not use in standalone ECRAD"
stop "sub map_par_to_ece_rp_scal_ne"
end subroutine map_par_to_ece_rp_scal_ne

END MODULE mod_fit_params
