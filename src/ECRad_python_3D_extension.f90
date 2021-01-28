!
! ECRad_python_3D_extension.f90
!
!  Created on: Jan 26, 2021
!      Author: denk

module ECRad_python_3D_extension

contains
subroutine initialize_ECRad_3D(N_ch, N_Te_spline_knots, N_ne_spline_knots, &
                               equilibrium_file, equilibrium_type, use_mesh, &
                               use_symmetry, B_ref, s_plus, s_max, &
                               interpolation_acc, fourier_coeff_trunc, &
                               h_mesh, delta_phi_mesh, vessel_filename, &
                               rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ECRad,        only: initialize_ECRad_3D_f2py
implicit none
integer, intent(in)      :: N_ch
integer(4), intent(in) :: N_Te_spline_knots, N_ne_spline_knots
CHARACTER(*), intent(in) :: equilibrium_file
CHARACTER(*), intent(in) :: equilibrium_type
logical, intent(in)      :: use_mesh, use_symmetry
real(8), intent(in)  :: B_ref, s_plus, s_max, h_mesh, delta_phi_mesh, &
                            interpolation_acc, fourier_coeff_trunc
CHARACTER(*), intent(in) :: vessel_filename
real(8), dimension(N_ch), intent(out)  :: rhopol_out
call initialize_ECRad_3D_f2py(N_Te_spline_knots, N_ne_spline_knots, &
                              equilibrium_file, equilibrium_type, use_mesh, &
                              use_symmetry, B_ref, s_plus, s_max, &
                              interpolation_acc, fourier_coeff_trunc, &
                              h_mesh, delta_phi_mesh, vessel_filename, &
                              rhopol_out)
end subroutine initialize_ECRad_3D
end module ECRad_python_3D_extension
