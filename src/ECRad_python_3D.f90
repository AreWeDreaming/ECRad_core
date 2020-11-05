module ECRad_python

contains

subroutine pre_initialize_ECRad(working_dir_in, flag, ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
                               rhopol_max_spline_knot, ecrad_weak_rel, &
                               ecrad_ratio_for_third_harmonic, &
                               ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                               ecrad_max_points_svec, &
                               ecrad_O2X_mode_conversion, &
                               rhopol_scal_te, rhopol_scal_ne, btf_corr_fact_ext, &
                               ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &
                               ecrad_z_shift, &
                               ecrad_N_ray, ecrad_N_freq, log_flag, parallelization_mode, &
                               f, df, R, phi, z, tor, pol, dist_foc, width)
! Everything that is absolutely static in time is done over here
use mod_ecfm_refr,      only: pre_initialize_ecfm
implicit none
character(*), intent(in)      :: working_dir_in, flag
real(kind=8), intent(in)      :: rhopol_max_spline_knot, ecrad_ratio_for_third_harmonic, &
                                              reflec_X_mode, reflec_O_mode, ecrad_O2X_mode_conversion, &
                                              rhopol_scal_te, rhopol_scal_ne, btf_corr_fact_ext, &
                                              ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, ecrad_z_shift
integer, intent(in)   :: ecrad_modes, ecrad_max_points_svec, ecrad_N_ray, &
                                              ecrad_N_freq,ece_1O_flag
logical, intent(in)           :: ecrad_verbose, ecrad_Bt_ripple, ray_tracing, ecrad_weak_rel, log_flag
integer, intent(in), optional  :: parallelization_mode
real(kind=8), dimension(:), intent(in), optional :: f, df, R, phi, z, tor, pol, dist_foc, width
call pre_initialize_ecfm(working_dir_in, flag, ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
                               rhopol_max_spline_knot, ecrad_weak_rel, &
                               ecrad_ratio_for_third_harmonic, &
                               ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                               ecrad_max_points_svec, & ! (modes = 1 -> pure X-mode, 2 -> pure O-mode, 3 both modes and filter
                               ecrad_O2X_mode_conversion, & ! mode conversion ratio from O-X due to wall reflections
                               ! Scaling of rhop axis for shifting on ne or Te
                               ! Every rhop value obtained in ray tracing will be multiplied by the corresponding scaling value when evaluating Te/ne
                               rhopol_scal_te, rhopol_scal_ne, btf_corr_fact_ext, &
                               ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_z_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_N_ray, ecrad_N_freq, log_flag, parallelization_mode, &
                               f, df, R, phi, z, tor, pol, dist_foc, width)
end subroutine pre_initialize_ECRad

subroutine initialize_ECRad(flag, N_ch, N_Te_spline_knots, N_ne_spline_knots, &
                            R, z, rhop, Br, Bt, Bz, R_ax, z_ax, rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ecfm_refr,        only: initialize_ecfm
implicit none
character(*), intent(in)                           :: flag
integer, intent(in)                        		   :: N_ch
integer, intent(in)                                :: N_Te_spline_knots, N_ne_spline_knots
real(kind=8), intent(in), optional                 :: R_ax, z_ax
real(kind=8), dimension(:), intent(in)             :: R, z
real(kind=8), dimension(:,:), intent(in)           :: rhop, Br, Bt, Bz
real(kind=8), dimension(N_ch), intent(out)         :: rhopol_out
call initialize_ecfm(flag, N_Te_spline_knots, N_ne_spline_knots, &
                     R, z, rhop, Br, Bt, Bz, R_ax, z_ax, rhopol_out)
end subroutine initialize_ECRad
#ifdef USE3D
subroutine initialize_ECRad_3D(flag, N_ch, N_Te_spline_knots, N_ne_spline_knots, &
                               equilibrium_file, equilibrium_type, use_mesh, &
                               use_symmetry, B_ref, s_plus, s_max, &
                               interpolation_acc, fourier_coeff_trunc, &
                               h_mesh, delta_phi_mesh, vessel_filename, &
                               rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ecfm_refr,        only: initialize_ecfm_3D
implicit none
character(*), intent(in) :: flag
integer, intent(in)      :: N_ch
integer(4), intent(in) :: N_Te_spline_knots, N_ne_spline_knots
CHARACTER(*), intent(in) :: equilibrium_file
CHARACTER(*), intent(in) :: equilibrium_type
logical, intent(in)      :: use_mesh, use_symmetry
real(8), intent(in)  :: B_ref, s_plus, s_max, h_mesh, delta_phi_mesh, &
                            interpolation_acc, fourier_coeff_trunc
CHARACTER(*), intent(in) :: vessel_filename
real(8), dimension(N_ch), intent(out)  :: rhopol_out
call initialize_ecfm_3D(flag, N_Te_spline_knots, N_ne_spline_knots, &
                        equilibrium_file, equilibrium_type, use_mesh, &
                        use_symmetry, B_ref, s_plus, s_max, &
                        interpolation_acc, fourier_coeff_trunc, &
                        h_mesh, delta_phi_mesh, vessel_filename, &
                        rhopol_out)
end subroutine initialize_ECRad_3D
#endif

subroutine make_rays_ECRad(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, rhop_res)
use mod_ecfm_refr,        only: make_rays_ecfm
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ecfm(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
					T_e=T_e,  rhop_res=rhop_res)
end subroutine make_rays_ECRad

subroutine make_rays_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                          rhop_res)
use mod_ecfm_refr,        only: make_rays_ecfm
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ecfm(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, rhop_res)
end subroutine make_rays_ECRad_spline

subroutine make_dat_model_ECRad(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ecfm_refr,        only: make_dat_model_ece_ecfm_refr
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
call make_dat_model_ece_ecfm_refr(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
								  T_e=T_e, ne_rhop_scal=ne_rhop_scal, reflec_X_new=reflec_X_new, & ! in
                                  reflec_O_new=reflec_O_new, ece_fm_flag_ch=ece_fm_flag_ch, rp_min=rp_min, &
                                  dat_model_ece=dat_model_ece, tau=tau, set_grid_dynamic=set_grid_dynamic, &
                                  verbose=verbose)
end subroutine make_dat_model_ECRad

subroutine make_dat_model_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ecfm_refr,        only: make_dat_model_ece_ecfm_refr
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
call make_dat_model_ece_ecfm_refr(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                  ne_rhop_scal, reflec_X_new, & ! in
                                  reflec_O_new, ece_fm_flag_ch, rp_min, &
                                  dat_model_ece, tau, set_grid_dynamic, verbose)
end subroutine make_dat_model_ECRad_spline

subroutine make_BPD_w_res_ch_ECRad(pnts_BPD, idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                   ne_rhop_scal, reflec_X_new, &
                                   reflec_O_new, rp_min, &
                                   rhop, BPD, rhop_res_warm)
use mod_ecfm_refr,        only: make_BPD_w_res_ch
implicit none
integer, intent(in)             :: pnts_BPD
integer, intent(in)             :: idiag, ich
real(kind=8), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e, n_e_dx2, T_e_dx2
real(kind=8),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
real(kind=8),  intent(out) :: rhop(pnts_BPD), BPD(pnts_BPD)
real(kind=8), intent(out)               :: rhop_res_warm
real(kind=8), dimension(:), allocatable :: rhop_ECRad_out, BPD_ECRad_out
call make_BPD_w_res_ch(idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                       ne_rhop_scal, reflec_X_new, & ! in
                       reflec_O_new, rp_min, &
                       rhop_ECRad_out, BPD_ECRad_out, rhop_res_warm)
rhop(1:pnts_BPD) = rhop_ECRad_out(1:pnts_BPD)
BPD(1:pnts_BPD) = BPD_ECRad_out(1:pnts_BPD)
deallocate(rhop_ECRad_out, BPD_ECRad_out)
end subroutine make_BPD_w_res_ch_ECRad

subroutine update_Te_ne_ECRad(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
use mod_ecfm_refr,               only: update_Te_ne
  implicit none
  real(kind=8), dimension(:), intent(in)   :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
  real(kind=8), dimension(:), intent(in), optional   :: n_e_dx2, T_e_dx2
  call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
end subroutine update_Te_ne_ECRad

end module ECRad_python