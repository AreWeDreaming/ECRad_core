module ECRad_python

contains

subroutine pre_initialize_ECRad(ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
                                rhopol_max_spline_knot, ecrad_weak_rel, &
                                ecrad_ratio_for_third_harmonic, &
                                ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                                ecrad_max_points_svec, &
                                ecrad_O2X_mode_conversion, &
                                rhopol_scal_te, rhopol_scal_ne, &
                                ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &
                                ecrad_z_shift, &
                                ecrad_N_ray, ecrad_N_freq, log_flag, N_vessel, vessel_R, vessel_z, &
                                f, df, R, phi, z, tor, pol, dist_foc, width)
! Everything that is absolutely static in time is done over here
use mod_ECRad,      only: pre_initialize_ECRad_f2py
implicit none
real(kind=8), intent(in)      :: rhopol_max_spline_knot, ecrad_ratio_for_third_harmonic, &
                                              reflec_X_mode, reflec_O_mode, ecrad_O2X_mode_conversion, &
                                              rhopol_scal_te, rhopol_scal_ne, &
                                              ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, ecrad_z_shift
integer, intent(in)   :: ecrad_modes, ecrad_max_points_svec, ecrad_N_ray, &
                                              ecrad_N_freq,ece_1O_flag
logical, intent(in)           :: ecrad_verbose, ecrad_Bt_ripple, ray_tracing, ecrad_weak_rel, log_flag
integer, intent(in)    :: N_vessel
real(kind=8), dimension(:), intent(in) :: vessel_R, vessel_z ! vessel contour
real(kind=8), dimension(:), intent(in), optional :: f, df, R, phi, z, tor, pol, dist_foc, width
call pre_initialize_ECRad_f2py(ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
                               rhopol_max_spline_knot, ecrad_weak_rel, &
                               ecrad_ratio_for_third_harmonic, &
                               ecrad_modes, reflec_X_mode, reflec_O_mode, ece_1O_flag, &
                               ecrad_max_points_svec, & ! (modes = 1 -> pure X-mode, 2 -> pure O-mode, 3 both modes and filter
                               ecrad_O2X_mode_conversion, & ! mode conversion ratio from O-X due to wall reflections
                               ! Scaling of rhop axis for shifting on ne or Te
                               ! Every rhop value obtained in ray tracing will be multiplied by the corresponding scaling value when evaluating Te/ne
                               rhopol_scal_te, rhopol_scal_ne, &
                               ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_z_shift, &    ! Allows shifting the equilbrium - moves entire flux matrix
                               ecrad_N_ray, ecrad_N_freq, log_flag, N_vessel, vessel_R, vessel_z, &
                               f, df, R, phi, z, tor, pol, dist_foc, width)
end subroutine pre_initialize_ECRad

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

subroutine make_rays_ECRad_f2py(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, rhop_res)
use mod_ECRad,        only: make_rays_ECRad
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ECRad(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
					           T_e=T_e,  rhop_res=rhop_res)
end subroutine make_rays_ECRad_f2py

subroutine make_rays_ECRad_2D_f2py(N_ch, rhop_res)
use mod_ECRad,        only: make_rays_ECRad
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
! For 2D profiles Te and ne are already in place
call make_rays_ECRad(rhop_res=rhop_res)
end subroutine make_rays_ECRad_2D_f2py

subroutine make_rays_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                          rhop_res)
use mod_ECRad,        only: make_rays_ECRad
implicit none
integer, intent(in)                        :: N_ch
real(kind=8), dimension(:), intent(in) :: rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2
real(kind=8), dimension(N_ch),  intent(out) :: rhop_res
integer                          :: idiag, ich
call make_rays_ECRad(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, rhop_res)
end subroutine make_rays_ECRad_spline

subroutine make_dat_model_ECRad(N_ch, rhop_knots_ne, n_e, rhop_knots_Te, T_e, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ECRad,        only: make_dat_model_ece_ECRad
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
call make_dat_model_ece_ECRad(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, &
								  T_e=T_e, ne_rhop_scal=ne_rhop_scal, reflec_X_new=reflec_X_new, & ! in
                                  reflec_O_new=reflec_O_new, ece_fm_flag_ch=ece_fm_flag_ch, rp_min=rp_min, &
                                  dat_model_ece=dat_model_ece, tau=tau, set_grid_dynamic=set_grid_dynamic, &
                                  verbose=verbose)
end subroutine make_dat_model_ECRad

subroutine make_TRad_direct()
use mod_ECRad,        only: make_ece_rad_temp
implicit none
! For single call. For usage with fixed Te/ne profiles.
call make_ece_rad_temp()
end subroutine make_TRad_direct

subroutine make_dat_model_ECRad_spline(N_ch, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                        ne_rhop_scal, reflec_X_new, &
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ECRad,        only: make_dat_model_ece_ECRad
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
  call make_dat_model_ece_ECRad(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                ne_rhop_scal, reflec_X_new, & ! in
                                reflec_O_new, ece_fm_flag_ch, rp_min, &
                                dat_model_ece, tau, set_grid_dynamic, verbose)
end subroutine make_dat_model_ECRad_spline

subroutine get_ray_length(ich, imode, ir, N_LOS)
use mod_ECRad_types,        only: rad
implicit None
integer, intent(in)   :: ich, imode, ir
integer, intent(out)   :: N_LOS
integer              :: idiag
  idiag = 1
  N_LOS = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
end subroutine get_ray_length

subroutine get_ray_data(ich, imode, ir, s, x, y, z, Nx, Ny, Nz, &
                        Bx, By, Bz, rho, T_e, n_e, theta, N_cold, H, v_g_perp, &
                        Trad, Trad_secondary, em, em_secondary, ab, ab_secondary, T, &
                        T_secondary, BPD, BPD_secondary)
use mod_ECRad_types,        only: rad
implicit None
integer, intent(in)   :: ich, imode, ir
real(kind=8), dimension(:), intent(out) :: s, x, y, z, Nx, Ny, Nz, Bx, By, Bz, rho, &
                                           T_e, n_e, theta, N_cold, H, v_g_perp, &
                                           Trad, Trad_secondary, em, em_secondary, &
                                           ab, ab_secondary, T, T_secondary, BPD, BPD_secondary
integer            :: idiag, N_points_LOS
  idiag = 1
  N_points_LOS = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
  s(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%s
  x(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%x_vec(1)
  y(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%x_vec(2)
  z(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%x_vec(3)
  Nx(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%N_vec(1)
  Ny(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%N_vec(1)
  Nz(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%N_vec(1)
  Bx(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%B_vec(1)
  By(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%B_vec(2)
  Bz(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%B_vec(3)
  rho(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%rhop
  T_e(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%Te
  n_e(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%ne
  theta(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%theta
  N_cold(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%N_cold
  H(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H(1:N_points_LOS)
  v_g_perp(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(1:N_points_LOS)%v_g_perp
  Trad(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(1:N_points_LOS)
  Trad_secondary(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary(1:N_points_LOS)
  em(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(1:N_points_LOS)
  em_secondary(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(1:N_points_LOS)
  ab(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(1:N_points_LOS)
  ab_secondary(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(1:N_points_LOS)
  T(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(1:N_points_LOS)
  T_secondary(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(1:N_points_LOS)
  BPD(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(1:N_points_LOS)
  BPD_secondary(1:N_points_LOS) = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(1:N_points_LOS)
end subroutine get_ray_data

subroutine get_BPD(ich, imode, rho, BPD, BPD_second)
use mod_ECRad_types,        only: rad
implicit None
integer, intent(in)   :: ich, imode
real(kind=8), dimension(:), intent(out) :: rho, BPD, BPD_second
integer :: idiag
  idiag = 1
  rho = rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD
  BPD = rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD
  BPD_second = rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary
end subroutine get_BPD


subroutine get_Trad_resonances_basic(imode, Trad, tau, s_res, R_res, z_res, rho_res)
use mod_ECRad_types,        only: rad, ant
implicit None
integer, intent(in)   :: imode
! <= 0 -> average over modes
! > 0 corresponding to imode in ECRad
real(kind=8), dimension(:), intent(out) :: Trad, tau, s_res, R_res, z_res, rho_res
integer :: idiag, ich
  if(imode <= 0) then
    Trad = rad%diag(idiag)%ch(:)%Trad
    tau = rad%diag(idiag)%ch(:)%tau
    s_res = rad%diag(idiag)%ch(:)%s_res
    R_res = rad%diag(idiag)%ch(:)%R_res
    z_res = rad%diag(idiag)%ch(:)%z_res
    rho_res = rad%diag(idiag)%ch(:)%rhop_res
  else
    do ich = 1, ant%diag(idiag)%N_ch
      Trad(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%Trad
      tau(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%tau
      s_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%s_res
      R_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%R_res
      z_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%z_res
      rho_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res
    end do
  end if

end subroutine get_Trad_resonances_basic

subroutine get_Trad_resonances_extra_output(imode, Trad_secondary, tau_secondary, &
                                            rel_s_res, rel_rho_res, rel_R_res, rel_z_res, &
                                            rel_s_res_secondary, rel_rho_res_secondary, &
                                            rel_R_res_secondary, rel_z_res_secondary)
use mod_ECRad_types,        only: rad, ant
implicit None
integer, intent(in)   :: imode
! <= 0 -> average over modes
! > 0 corresponding to imode in ECRad
real(kind=8), dimension(:), intent(out) :: Trad_secondary, tau_secondary, &
                                           rel_s_res, rel_rho_res, rel_R_res, rel_z_res, &
                                           rel_s_res_secondary, rel_rho_res_secondary, &
                                           rel_R_res_secondary, rel_z_res_secondary
integer(kind=8) :: idiag, ich
  if(imode <= 0) then
    Trad_secondary = rad%diag(idiag)%ch(:)%Trad_secondary
    tau_secondary = rad%diag(idiag)%ch(:)%tau_secondary
    rel_s_res = rad%diag(idiag)%ch(:)%rel_s_res
    rel_rho_res = rad%diag(idiag)%ch(:)%rel_rhop_res
    rel_R_res = rad%diag(idiag)%ch(:)%rel_R_res
    rel_z_res = rad%diag(idiag)%ch(:)%rel_z_res
    rel_s_res_secondary = rad%diag(idiag)%ch(:)%rel_s_res_secondary
    rel_rho_res_secondary = rad%diag(idiag)%ch(:)%rel_rhop_res_secondary
    rel_R_res_secondary = rad%diag(idiag)%ch(:)%rel_R_res_secondary
    rel_z_res_secondary = rad%diag(idiag)%ch(:)%rel_z_res_secondary
  else
    do ich = 1, ant%diag(idiag)%N_ch
      Trad_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary
      tau_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary
      rel_s_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res
      rel_rho_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res
      rel_R_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res
      rel_z_res(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res
      rel_s_res_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res_secondary
      rel_rho_res_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res_secondary
      rel_R_res_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res_secondary
      rel_z_res_secondary(ich) = rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res_secondary
    end do
  end if

end subroutine get_Trad_resonances_extra_output



subroutine make_BPD_w_res_ch_ECRad(pnts_BPD, idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                   ne_rhop_scal, reflec_X_new, &
                                   reflec_O_new, rp_min, &
                                   rhop, BPD, rhop_res_warm)
use mod_ECRad,        only: make_BPD_w_res_ch
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
use mod_ECRad,               only: update_Te_ne
  implicit none
  real(kind=8), dimension(:), intent(in)   :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
  real(kind=8), dimension(:), intent(in), optional   :: n_e_dx2, T_e_dx2
  call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
end subroutine update_Te_ne_ECRad

end module ECRad_python
