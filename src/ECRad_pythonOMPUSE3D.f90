module ECRad_python

contains

subroutine pre_initialize_ECRad(ecrad_verbose, dstf, ray_tracing, ecrad_Bt_ripple, &
                                rhopol_max_spline_knot, ecrad_weak_rel, &
                                ecrad_ratio_for_third_harmonic, &
                                ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                                ecrad_max_points_svec, N_pts_BPD, &
                                ecrad_O2X_mode_conversion, &
                                rhopol_scal_te, rhopol_scal_ne, &
                                ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &
                                ecrad_z_shift, &
                                ecrad_N_ray, ecrad_N_freq, log_flag, N_vessel, vessel_R, vessel_z, &
                                f, df, R, phi, z, tor, pol, dist_foc, width, pol_coeff)
! Everything that is absolutely static in time is done over here
use mod_ECRad,      only: pre_initialize_ECRad_f2py
implicit none
character(*), intent(in)        :: dstf
real(kind=8), intent(in)      :: rhopol_max_spline_knot, ecrad_ratio_for_third_harmonic, &
                                              reflec_X_mode, reflec_O_mode, ecrad_O2X_mode_conversion, &
                                              rhopol_scal_te, rhopol_scal_ne, &
                                              ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, ecrad_z_shift
integer, intent(in)   :: ecrad_modes, ecrad_max_points_svec, N_pts_BPD, ecrad_N_ray, &
                         ecrad_N_freq,ece_1O_flag
logical, intent(in)           :: ecrad_verbose, ecrad_Bt_ripple, ray_tracing, ecrad_weak_rel, log_flag
integer, intent(in)    :: N_vessel
real(kind=8), dimension(:), intent(in) :: vessel_R, vessel_z ! vessel contour
real(kind=8), dimension(:), intent(in), optional :: f, df, R, phi, z, tor, pol, dist_foc, width, pol_coeff
call pre_initialize_ECRad_f2py(ecrad_verbose, dstf, ray_tracing, ecrad_Bt_ripple, &
                               rhopol_max_spline_knot, ecrad_weak_rel, &
                               ecrad_ratio_for_third_harmonic, &
                               ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                               ecrad_max_points_svec, N_pts_BPD ,& ! (modes = 1 -> pure X-mode, 2 -> pure O-mode, 3 both modes and filter
                               ecrad_O2X_mode_conversion, & ! mode conversion ratio from O-X due to wall reflections
                               ! Scaling of rhop axis for shifting on ne or Te
                               ! Every rhop value obtained in ray tracing will be multiplied by the corresponding scaling value when evaluating Te/ne
                               rhopol_scal_te, rhopol_scal_ne, &
                               ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_z_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_N_ray, ecrad_N_freq, log_flag, N_vessel, vessel_R, vessel_z, &
                               f, df, R, phi, z, tor, pol, dist_foc, width, pol_coeff)
end subroutine pre_initialize_ECRad

subroutine set_ECRad_thread_count(num_threads)
use mod_ECRad, only: set_omp_threads_ECRad_f2py
implicit None
  integer, intent(in) :: num_threads
  call set_omp_threads_ECRad_f2py(num_threads)
end subroutine set_ECRad_thread_count

subroutine reset_ECRad()
use mod_ECRad,      only: clean_up_ECRad
implicit none
  call clean_up_ECRad()
end subroutine reset_ECRad

subroutine initialize_ECRad(N_ch, N_Te_spline_knots, N_ne_spline_knots, &
                            R, z, rhop, Br, Bt, Bz, R_ax, z_ax, rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ECRad,        only: initialize_ECRad_f2py
implicit none
integer, intent(in)                        		     :: N_ch
integer, intent(in)                                :: N_Te_spline_knots, N_ne_spline_knots
real(kind=8), intent(in), optional                 :: R_ax, z_ax
real(kind=8), dimension(:), intent(in)             :: R, z
real(kind=8), dimension(:,:), intent(in)           :: rhop, Br, Bt, Bz
real(kind=8), dimension(N_ch), intent(out)         :: rhopol_out
call initialize_ECRad_f2py(N_Te_spline_knots, N_ne_spline_knots, &
                           R, z, rhop, Br, Bt, Bz, R_ax, z_ax)
end subroutine initialize_ECRad

subroutine initialize_ECRad_2D_profs(N_ch, N_Te_spline_knots, N_ne_spline_knots, &
                            R, z, T_e, n_e, rhop, Br, Bt, Bz, R_ax, z_ax)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ECRad,        only: initialize_ECRad_f2py
implicit none
integer, intent(in)                                :: N_ch
integer, intent(in)                                :: N_Te_spline_knots, N_ne_spline_knots
real(kind=8), intent(in), optional                 :: R_ax, z_ax
real(kind=8), dimension(:), intent(in)             :: R, z
real(kind=8), dimension(:,:), intent(in)           :: rhop, Br, Bt, Bz, T_e, n_e
call initialize_ECRad_f2py(N_Te_spline_knots, N_ne_spline_knots, &
                           R, z, rhop, Br, Bt, Bz, R_ax, z_ax, T_e, n_e)
end subroutine initialize_ECRad_2D_profs

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

subroutine set_ECRad_FFP_dist(rho, u, pitch, f)
  use mod_ECRad,        only: set_ECRad_FFP_dist_f2py
  implicit none
  real(kind=8), dimension(:), intent(in)             :: rho, u, pitch
  real(kind=8), dimension(:,:,:), intent(in)         :: f
  call set_ECRad_FFP_dist_f2py(rho, u, pitch, f)
end subroutine set_ECRad_FFP_dist

subroutine set_ECRad_GENE_dist(rho, vpar, mu, f_0, g)
  use mod_ECRad,        only: set_ECRad_GENE_dist_f2py
  implicit none
  real(kind=8), dimension(:), intent(in)             :: rho, vpar, mu
  real(kind=8), dimension(:,:), intent(in)           :: f_0
  real(kind=8), dimension(:,:,:), intent(in)         :: g
  call set_ECRad_GENE_dist_f2py(rho, vpar, mu, f_0, g)
end subroutine set_ECRad_GENE_dist

subroutine make_rays_ECRad(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, rhop_res)
use mod_ECRad,        only: make_rays_ECRad_f2py
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ECRad_f2py(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
					                T_e=T_e,  rhop_res=rhop_res)
end subroutine make_rays_ECRad

subroutine make_rays_ECRad_2D(N_ch, rhop_res)
use mod_ECRad,        only: make_rays_ECRad_F2py
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
! For 2D profiles Te and ne are already in place
call make_rays_ECRad_f2py(rhop_res=rhop_res)
end subroutine make_rays_ECRad_2D

subroutine make_rays_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                  rhop_res)
use mod_ECRad,        only: make_rays_ECRad_f2py
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ECRad_f2py(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, rhop_res)
end subroutine make_rays_ECRad_spline

subroutine make_dat_model_ECRad(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ECRad,        only: make_dat_model_ece_ECRad_f2py
implicit none
integer, intent(in)             :: N_ch
real(kind=8), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(kind=8),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
logical,     dimension(:), intent(in)  :: ece_fm_flag_ch
real(kind=8), intent(out) :: dat_model_ece(N_ch)
real(kind=8), intent(out)  :: tau(N_ch)
logical,      intent(in)     :: verbose
logical, intent(in)          :: set_grid_dynamic
integer                         :: ich
call make_dat_model_ece_ECRad_f2py(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
								                   T_e=T_e, ne_rhop_scal=ne_rhop_scal, reflec_X_new=reflec_X_new, & ! in
                                   reflec_O_new=reflec_O_new, ece_fm_flag_ch=ece_fm_flag_ch, rp_min=rp_min, &
                                   dat_model_ece=dat_model_ece, tau=tau, set_grid_dynamic=set_grid_dynamic, &
                                   verbose=verbose)
end subroutine make_dat_model_ECRad



subroutine make_TRad_direct()
use mod_ECRad,        only: make_ece_rad_temp, set_for_single_eval
implicit none
integer :: idiag
idiag = 1
! For single call. For usage with fixed Te/ne profiles.
call set_for_single_eval()
call make_ece_rad_temp()
end subroutine make_TRad_direct

subroutine make_dat_model_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ECRad,        only: make_dat_model_ece_ECRad_f2py
implicit none
integer, intent(in)             :: N_ch
real(kind=8), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(kind=8),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
real(kind=8), dimension(:), intent(in)  :: n_e_dx2, T_e_dx2
logical,     dimension(:), intent(in)   :: ece_fm_flag_ch
real(kind=8), intent(out) :: dat_model_ece(N_ch)
real(kind=8), intent(out)  :: tau(N_ch)
logical,      intent(in)     :: verbose
logical, intent(in)          :: set_grid_dynamic
integer                         :: ich
  call make_dat_model_ece_ECRad_f2py(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                ne_rhop_scal, reflec_X_new, & ! in
                                reflec_O_new, ece_fm_flag_ch, rp_min, &
                                dat_model_ece, tau, set_grid_dynamic, verbose)
end subroutine make_dat_model_ECRad_spline


subroutine update_Te_ne_ECRad(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
use mod_ECRad,               only: update_Te_ne
  implicit none
  real(kind=8), dimension(:), intent(in)   :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
  real(kind=8), dimension(:), intent(in), optional   :: n_e_dx2, T_e_dx2
  call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
end subroutine update_Te_ne_ECRad

subroutine get_Trad_resonances_basic(imode, N_ch, Trad, tau, s_res, R_res, z_res, rho_res)
use mod_ECRad,        only: get_Trad_resonances_basic_f2py
implicit None
integer, intent(in)   :: imode, N_ch
! <= 0 -> average over modes
! > 0 corresponding to imode in ECRad
real(kind=8), dimension(N_ch), intent(out) :: Trad, tau, s_res, R_res, z_res, rho_res
  call get_Trad_resonances_basic_f2py(imode, N_ch, Trad, tau, s_res, R_res, z_res, rho_res)
end subroutine get_Trad_resonances_basic

subroutine get_Trad_resonances_extra_output(imode, N_ch, Trad_secondary, tau_secondary, &
                                            rel_s_res, rel_rho_res, &
                                            rel_R_res, rel_z_res, &
                                            rel_s_res_secondary, rel_rho_res_secondary, &
                                            rel_R_res_secondary, rel_z_res_secondary)
use mod_ECRad,        only: get_Trad_resonances_extra_output_f2py
implicit None
integer, intent(in)   :: imode, N_ch
! <= 0 -> average over modes
! > 0 corresponding to imode in ECRad
real(kind=8), dimension(N_ch), intent(out) :: Trad_secondary, tau_secondary, &
                                              rel_s_res, rel_rho_res, &
                                              rel_R_res, rel_z_res, &
                                              rel_s_res_secondary, rel_rho_res_secondary, &
                                              rel_R_res_secondary, rel_z_res_secondary

  call get_Trad_resonances_extra_output_f2py(imode, Trad_secondary, tau_secondary, &
                                             rel_s_res, rel_rho_res, &
                                             rel_R_res, rel_z_res, &
                                             rel_s_res_secondary, rel_rho_res_secondary, &
                                             rel_R_res_secondary, rel_z_res_secondary)
end subroutine get_Trad_resonances_extra_output

subroutine get_BPD(ich, imode, N_BPD, rho, BPD, BPD_second)
use mod_ECRad,        only: get_BPD_f2py
implicit None
integer, intent(in)   :: ich, imode, N_BPD
real(kind=8), dimension(N_BPD), intent(out) :: rho, BPD, BPD_second
  call get_BPD_f2py(ich, imode, rho, BPD, BPD_second)
end subroutine get_BPD

subroutine get_ray_length(ich, imode, ir, N_LOS)
use mod_ECRad,        only: get_ray_length_f2py
implicit None
integer, intent(in)   :: ich, imode, ir
integer, intent(out)  :: N_LOS
integer               :: idiag
  call get_ray_length_f2py(ich, imode, ir, N_LOS)
end subroutine get_ray_length

subroutine get_ray_data(ich, imode, ir, N_LOS, s, x, y, z, Nx, Ny, Nz, &
                        Bx, By, Bz, rho, T_e, n_e, theta, N_cold, H, v_g_perp, &
                        Trad, Trad_secondary, em, em_secondary, ab, ab_secondary, T, &
                        T_secondary, BPD, BPD_secondary)
use mod_ECRad,        only: get_ray_data_f2py
implicit None
integer, intent(in)   :: ich, imode, ir, N_LOS
real(kind=8), dimension(N_LOS), intent(out) :: s, x, y, z, Nx, Ny, Nz, Bx, By, Bz, rho, &
                                               T_e, n_e, theta, N_cold, H, v_g_perp, &
                                               Trad, Trad_secondary, em, em_secondary, &
                                               ab, ab_secondary, T, T_secondary, BPD, BPD_secondary
  call get_ray_data_f2py(ich, imode, ir, s, x, y, z, Nx, Ny, Nz, &
                        Bx, By, Bz, rho, T_e, n_e, theta, N_cold, H, v_g_perp, &
                        Trad, Trad_secondary, em, em_secondary, ab, ab_secondary, T, &
                        T_secondary, BPD, BPD_secondary)
end subroutine get_ray_data

subroutine get_mode_weights(N_ch, imode, pol_coeff, pol_coeff_secondary)
  use mod_ECRad,        only:  get_mode_weights_f2py
  implicit none
    integer, intent(in) :: N_ch, imode
    real(kind=8), dimension(N_ch), intent(out) :: pol_coeff, pol_coeff_secondary
    call get_mode_weights_f2py(N_ch, imode, pol_coeff, pol_coeff_secondary)
end subroutine get_mode_weights

subroutine get_weights(N_ray, N_freq, ich,  ray_weights, freq_weights)
use mod_ECRad,        only: get_weights_f2py
implicit none
  integer, intent(in)   :: N_ray, N_freq, ich
  real(kind=8), dimension(N_ray), intent(out) :: ray_weights
  real(kind=8), dimension(N_freq), intent(out) :: freq_weights
  call get_weights_f2py(ich, ray_weights, freq_weights)
end subroutine get_weights

!subroutine make_BPD_w_res_ch_ECRad(pnts_BPD, idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
!                                   ne_rhop_scal, reflec_X_new, &
!                                   reflec_O_new, rp_min, &
!                                   rhop, BPD, rhop_res_warm)
!use mod_ECRad,        only: make_BPD_w_res_ch
!implicit none
!integer, intent(in)             :: pnts_BPD
!integer, intent(in)             :: idiag, ich
!real(kind=8), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e, n_e_dx2, T_e_dx2
!real(kind=8),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
!real(kind=8),  intent(out) :: rhop(pnts_BPD), BPD(pnts_BPD)
!real(kind=8), intent(out)               :: rhop_res_warm
!real(kind=8), dimension(:), allocatable :: rhop_ECRad_out, BPD_ECRad_out
!call make_BPD_w_res_ch(idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
!                       ne_rhop_scal, reflec_X_new, & ! in
!                       reflec_O_new, rp_min, &
!                       rhop_ECRad_out, BPD_ECRad_out, rhop_res_warm)
!rhop(1:pnts_BPD) = rhop_ECRad_out(1:pnts_BPD)
!BPD(1:pnts_BPD) = BPD_ECRad_out(1:pnts_BPD)
!deallocate(rhop_ECRad_out, BPD_ECRad_out)
!end subroutine make_BPD_w_res_ch_ECRad

end module ECRad_python
