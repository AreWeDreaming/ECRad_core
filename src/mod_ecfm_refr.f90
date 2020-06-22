
module mod_ecfm_refr

use f90_kind
#ifdef OMP
  use omp_lib
#endif

public :: initialize_stand_alone, &
          initialize_ecfm, &
          make_rays_ecfm, &
          make_dat_model_ece_ecfm_refr, &
          pre_initialize_ecfm, &
          make_ece_rad_temp

private :: save_data_to_ASCII

contains



subroutine initialize_stand_alone(working_dir, flag)
! Simulates the structure used in IDA
! While in the stand alone make_ece_rad_temp is only called once it is called many times in IDA
! Hence, to keep the structure similiar all initizalization is performed here
use mod_ecfm_refr_types,        only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, data_name, &
                                      dstf_comp, plasma_params, warm_plasma, data_secondary_name, &
                                      data_folder, Ich_name, dstf_comp, straight, stand_alone, ffp, N_absz, &
                                      N_absz_large, new_IO, ray_init, use_ext_rays
use mod_ecfm_refr_utils,      only: read_input_file, prepare_ECE_diag_new_IO, &
                                    prepare_ECE_diag_old_IO, init_non_therm, &
                                    read_external_rays
use mod_ecfm_refr_raytrace_initialize,    only: init_raytrace
use mod_ecfm_refr_raytrace,               only: span_svecs
use mod_ecfm_refr_abs_Al,         only: abs_Al_init,abs_Al_clean_up
use constants,                    only: pi
implicit none
character(150), intent(in)    :: working_dir
character(*),  intent(in)      :: flag
  if(.not. stand_alone) stop "This function may only be called in stand alone"
  if(trim(flag) == "init") then
    call read_input_file(working_dir)
    if(.not. new_IO) then
      call prepare_ECE_diag_old_IO()
    else
      call prepare_ECE_diag_new_IO()
    end if
    if(use_ext_rays) call read_external_rays()
    dstf_comp = "Th"
    if(trim(dstf) == "Th") then
      dstf = "Th"
      dstf_comp = "DF"! Second absorption coefficient according to D. Farina's paper
      data_name = "TRadM_therm.dat"
      data_secondary_name = "TRadM_Farina.dat"
      Ich_name = "IchTh"
    else if(trim(dstf) == "Re") then
      dstf = "numeric"
      dst_data_folder = "fRe"
      data_name = "TRadM_RELAX.dat"
      Ich_name = "IchRe"
      data_secondary_name = "TRadM_therm.dat"
      ffp%LUKE = .false.
    else if(trim(dstf) == "Lu") then
      dstf = "numeric"
      dst_data_folder = "fLu"
      data_name = "TRadM_LUKE.dat"
      Ich_name = "IchLu"
      data_secondary_name = "TRadM_therm.dat"
      ffp%LUKE = .true.
    else if(trim(dstf) == "Ge") then
      dstf = "gene"
      dst_data_folder = "fGe"
      data_name = "TRadM_GENE.dat"
      Ich_name = "IchGe"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "GB") then
      dstf = "gcomp"
      dst_data_folder = "fGB"
      data_name = "TRadM_GENE.dat"
      Ich_name = "IchGB"
      data_secondary_name = "TRadM_GComp.dat"
    else if(trim(dstf) == "SH") then
      dstf = "Spitzer"
      data_name = "TRadM_SpitH.dat"
      Ich_name = "IchSH"
      data_secondary_name = "TRadM_therm.dat"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "BJ") then
      dstf = "Bi_MaxJ"
      data_name = "TRadM_BiMnJ.dat"
      Ich_name = "IchBJ"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "BM") then
      dstf = "Bi_Maxw"
      data_name = "TRadM_BiMax.dat"
      Ich_name = "IchBM"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "DM") then
      dstf = "drift_m"
      data_name = "TRadM_Drift.dat"
      Ich_name = "IchDM"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "MS") then
      dstf = "multi_s"
      data_name = "TRadM_MultS.dat"
      Ich_name = "IchMS"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "RA") then
      dstf = "runaway"
      data_name = "TRadM_RunAw.dat"
      Ich_name = "IchRA"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "Ta") then
      dstf = "Tanalyt"
      data_name = "TRadM_analy.dat"
      Ich_name = "IchTa"
      data_secondary_name = "TRadM_therm.dat"
    else
      print*, "Invalid flag for dstf: ", dstf
      stop "Input Error"
    end if
    if(output_level) print*,"Chosen distribution is: ", dstf
    if(output_level) flush(6)
    if(.not. new_IO) then
      data_folder = trim(working_dir) //  "ecfm_data" // "/"
    else
      data_folder = trim(working_dir) //  "ECRad_data" // "/"
    end if
    if(output_level) print*, "Data will arrive in: ", data_folder
    if(dstf == "numeric" .or. trim(dstf) == "gene" .or. trim(dstf) == "gcomp") then
      call abs_Al_init(N_absz_large) ! Initializes the weights and abszissae for the gaussian quadrature
    else
      call abs_Al_init(N_absz) ! Initializes the weights and abszissae for the gaussian quadrature
    end if
    if(output_level) then
      print*, "Integrated ray tracing enabled"
      print*,"Rays are preinitialized"
    end if
    ray_out_folder = trim(data_folder) // "ray/"
    call init_raytrace(plasma_params)
    call init_non_therm() ! Reads input data for non-thermal distributions
    if(output_level) then
      print*, "----- Options for raytracing ------"
      print*, "Force straight lines of sight: ", straight
      print*, "Include ripple: ", plasma_params%w_ripple
      print*, "Use weakly relativistic cut off correction for ray tracing: ", warm_plasma
    end if
    call span_svecs(plasma_params)
    ray_init = .true.
  else if(trim(flag) == "clean") then
    !TODO: Implement clean up routine
    call abs_Al_clean_up()
  else
    print*, "Inappropriate flag in initialize in mod_ecfm_ECE_rad"
    stop "Interal error"
  end if
end subroutine initialize_stand_alone

subroutine pre_initialize_ecfm(working_dir_in, flag, ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
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
! Everything that is absolutely static in time is done over here
use mod_ecfm_refr_types,      only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, &
                                      dstf_comp, plasma_params, N_ray, N_freq, ray_init, working_dir, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, straight, stand_alone
use mod_ecfm_refr_utils,      only: parse_ecfm_settings_from_ida, &
                                    prepare_ECE_diag_IDA, dealloc_ant
use mod_ecfm_refr_abs_Al,     only: abs_Al_init,abs_Al_clean_up
use mod_ecfm_refr_raytrace,   only: dealloc_rad
use mod_ecfm_refr_raytrace_initialize, only: dealloc_raytrace
implicit none
character(*), intent(in)      :: working_dir_in, flag
real(rkind), intent(in)       :: rhopol_max_spline_knot, ecrad_ratio_for_third_harmonic, &
                                              reflec_X_mode, reflec_O_mode, ecrad_O2X_mode_conversion, &
                                              rhopol_scal_te, rhopol_scal_ne, btf_corr_fact_ext, &
                                              ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, ecrad_z_shift
integer(ikind), intent(in)    :: ecrad_modes, ecrad_max_points_svec, ecrad_N_ray, &
                                              ecrad_N_freq,ece_1O_flag
logical, intent(in)           :: ecrad_verbose, ecrad_Bt_ripple, ray_tracing, ecrad_weak_rel, log_flag
integer(ikind), intent(in), optional  :: parallelization_mode
real(rkind), dimension(:), intent(in), optional :: f, df, R, phi, z, tor, pol, dist_foc, width
integer(ikind)                :: idiag
if(trim(flag) == "init" .or. trim(flag) == "load") then
! A few things need to be done both times
  working_dir = trim(working_dir_in)  // "/"
  if(present(parallelization_mode)) then
    call parse_ecfm_settings_from_ida(plasma_params, &
                                        ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
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
                                        ecrad_N_ray, ecrad_N_freq, log_flag, parallelization_mode)
  else
    call parse_ecfm_settings_from_ida(plasma_params, &
                                        ecrad_verbose, ray_tracing, ecrad_Bt_ripple, &
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
                                        ecrad_N_ray, ecrad_N_freq, log_flag, 0)
  end if
  ! Hard coded since this is anyways just for testing purposes
  data_folder = trim(working_dir) //  "ecfm_data" // "/"
  ray_out_folder = trim(working_dir) // "ecfm_data/" // "ray/"
  data_folder = trim(working_dir) // "ecfm_data" // "/"
  if(output_level) print*, "ECFM data will arrive in: ", data_folder
  if(dstf == "Th") call abs_Al_init(32) ! Initializes the weights and abszissae for the gaussian quadrature
  if(trim(flag) == "init") then
    if(.not. present(f)) then
      print*, "pre_initialize_ecfm must be called with ece_strut present if flag == init"
      call abort()
    end if
    call prepare_ECE_diag_IDA(working_dir=working_dir, f=f, df=df, R=R, &
    					  phi=phi, z=z, tor=tor, pol=pol, dist_foc=dist_foc, &
    					  width=width)
    ! parses diag info then creates two input files
  else
    call prepare_ECE_diag_IDA(working_dir=working_dir)
    ! load the files the routine 4 lines above creates
  end if
else if(trim(flag) == "clean") then
  call abs_Al_clean_up()
  call dealloc_rad(rad)
  call dealloc_ant(ant)
  call dealloc_raytrace(plasma_params)
  ray_init = .false.
else
  print*, "Inappropriate flag in initialize in mod_ecfm_ECE_rad"
  stop "Internal error"
end if
end subroutine pre_initialize_ecfm


subroutine clean_up_ecfm()
use mod_ecfm_refr_types,      only: plasma_params,rad, ant, ray_init
use mod_ecfm_refr_utils,      only: dealloc_ant
use mod_ecfm_refr_abs_Al,     only: abs_Al_clean_up
use mod_ecfm_refr_raytrace,   only: dealloc_rad
use mod_ecfm_refr_raytrace_initialize, only: dealloc_raytrace
implicit None
  call abs_Al_clean_up()
  call dealloc_rad(rad)
  call dealloc_ant(ant)
  call dealloc_raytrace(plasma_params)
  ray_init = .false.
end subroutine clean_up_ecfm

subroutine initialize_ecfm(flag, N_Te_spline_knots, N_ne_spline_knots, &
                           R, z, rhop, Br, Bt, Bz, R_ax, z_ax, rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ecfm_refr_types,        only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, &
                                      dstf_comp, plasma_params, N_ray, N_freq, ray_init, warm_plasma, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, straight, stand_alone, &
                                      use_ida_spline_Te, use_ida_spline_ne
use mod_ecfm_refr_raytrace_initialize, only: init_raytrace, dealloc_raytrace
implicit none
character(*), intent(in)                          :: flag
integer(ikind), intent(in), optional              :: N_Te_spline_knots, N_ne_spline_knots
real(rkind), intent(in), optional                 :: R_ax, z_ax
real(rkind), dimension(:), intent(in), optional   :: R, z
real(rkind), dimension(:,:), intent(in), optional :: rhop, Br, Bt, Bz
real(rkind), dimension(:), intent(out), optional  :: rhopol_out
integer(ikind)                       :: idiag
logical                                           :: old_straight
  if(trim(flag) == "init") then
    call init_raytrace(plasma_params, R, z, rhop, Br, Bt, Bz, R_ax, z_ax)
    allocate(plasma_params%IDA_rhop_knots_ne(N_ne_spline_knots), plasma_params%IDA_n_e(N_ne_spline_knots), &
             plasma_params%IDA_n_e_dx2(N_ne_spline_knots), plasma_params%IDA_rhop_knots_Te(N_Te_spline_knots), &
             plasma_params%IDA_T_e(N_Te_spline_knots), plasma_params%IDA_T_e_dx2(N_Te_spline_knots))
    if(present(rhopol_out)) then
      ! straight los => Neither Te nor ne matters
      plasma_params%No_ne_te = .true.
      old_straight = straight
      straight = .True.
      call make_rays_ecfm(rhopol_out) ! initial resonances without ray tracing
      plasma_params%No_ne_te = .false.
      straight = old_straight
    end if
    ! Two modes one for shot file load one for loading from files
    if(output_level) then
      print*, "----- Options for raytracing ------"
      print*, "Force straight lines of sight: ", straight
      print*, "Include ripple: ", plasma_params%w_ripple
      print*, "Use weakly relativistic cut off correction for ray tracing: ", warm_plasma
    end if
  else if(trim(flag) == "clean") then
    !TODO: Implement clean up routine
    call dealloc_raytrace(plasma_params)
    ray_init = .false.
    deallocate(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
               plasma_params%IDA_n_e_dx2, plasma_params%IDA_rhop_knots_Te, &
               plasma_params%IDA_T_e, plasma_params%IDA_T_e_dx2)
  else
    print*, "Inappropriate flag in initalize in mod_ecfm_ECE_rad"
    stop "Interal error"
  end if
end subroutine initialize_ecfm

#ifdef USE_3D
subroutine initialize_ecfm_3D(flag, N_Te_spline_knots, N_ne_spline_knots, &
                              equilibrium_file, equilibrium_type, use_mesh, &
                              use_symmetry, B_ref, s_plus, s_max, &
                              interpolation_acc, fourier_coeff_trunc, &
                              h_mesh, delta_phi_mesh, vessel_filename, &
                              rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ecfm_refr_types,        only: plasma_params, use_3D, Vessel_bd_filename, &
									  straight, output_level, ray_init, warm_plasma
use mod_ecfm_refr_raytrace_initialize, only: init_raytrace, dealloc_raytrace
implicit none
character(*), intent(in) :: flag
integer(ikind), intent(in), optional :: N_Te_spline_knots, N_ne_spline_knots
CHARACTER(*), intent(in) :: equilibrium_file
CHARACTER(*), intent(in) :: equilibrium_type
logical, intent(in)      :: use_mesh, use_symmetry
real(rkind), intent(in)  :: B_ref, s_plus, s_max, h_mesh, delta_phi_mesh, &
                            interpolation_acc, fourier_coeff_trunc
CHARACTER(*), intent(in) :: vessel_filename                            
real(rkind), dimension(:), intent(out), optional  :: rhopol_out
integer(ikind)                                    :: idiag
logical                                           :: old_straight
  if(trim(flag) == "init") then
    use_3D = .true.
    plasma_params%Scenario%name_config = equilibrium_file
    plasma_params%Scenario%format_config = equilibrium_type
    plasma_params%Scenario%useMesh = use_mesh
    plasma_params%Scenario%useSymm = use_symmetry
    plasma_params%Scenario%B_ref = B_ref
    plasma_params%Scenario%splus = s_plus
    plasma_params%Scenario%smax = s_max
    plasma_params%Scenario%accbooz = interpolation_acc
    plasma_params%Scenario%tolharm = fourier_coeff_trunc
    plasma_params%Scenario%hgrid = h_mesh
    plasma_params%Scenario%dphic = delta_phi_mesh! Degrees
    vessel_bd_filename = vessel_filename
    call init_raytrace(plasma_params)
    allocate(plasma_params%IDA_rhop_knots_ne(N_ne_spline_knots), plasma_params%IDA_n_e(N_ne_spline_knots), &
             plasma_params%IDA_n_e_dx2(N_ne_spline_knots), plasma_params%IDA_rhop_knots_Te(N_Te_spline_knots), &
             plasma_params%IDA_T_e(N_Te_spline_knots), plasma_params%IDA_T_e_dx2(N_Te_spline_knots))
    if(present(rhopol_out)) then
      ! straight los => Neither Te nor ne matters
      plasma_params%No_ne_te = .true.
      old_straight = straight
      straight = .True.
      call make_rays_ecfm(rhopol_out) ! initial resonances without ray tracing
      plasma_params%No_ne_te = .false.
      straight = old_straight
    end if
    ! Two modes one for shot file load one for loading from files
    if(output_level) then
      print*, "----- Options for raytracing ------"
      print*, "Force straight lines of sight: ", straight
      print*, "Include ripple: ", plasma_params%w_ripple
      print*, "Use weakly relativistic cut off correction for ray tracing: ", warm_plasma
    end if
  else if(trim(flag) == "clean") then
    !TODO: Implement clean up routine
    call dealloc_raytrace(plasma_params)
    ray_init = .false.
    deallocate(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
               plasma_params%IDA_n_e_dx2, plasma_params%IDA_rhop_knots_Te, &
               plasma_params%IDA_T_e, plasma_params%IDA_T_e_dx2)
  else
    print*, "Inappropriate flag in initalize in mod_ecfm_ECE_rad"
    stop "Interal error"
  end if
end subroutine initialize_ecfm_3D
#endif

subroutine make_rays_ecfm(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                          rhop_res)
! Simulates the structure used in IDA
! While in the stand alone make_ece_rad_temp is only called once it is called many times in IDA
! Hence, to keep the structure similiar all initizalization is performed here
use mod_ecfm_refr_types,        only: plasma_params, ant, rad, ray_init, straight, modes
use mod_ecfm_refr_raytrace,               only: span_svecs
implicit none
real(rkind), dimension(:), intent(in), optional   :: rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2
real(rkind), dimension(:),  intent(out), optional :: rhop_res
integer(ikind)                          :: idiag, ich
! Updates the values of the splines: Te, ne
if(.not. plasma_params%No_ne_te ) then
  call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
end if
call span_svecs(plasma_params)
! for second harmonic (n -1)/n =(2 -1)/2 = 0.5d0
ray_init = .true.
if( plasma_params%No_ne_te) then
  ! For the first run of the ecfm we need resonance positions that work best for classical ECE analysis
  ! Hence, we use the resonances of the X-mode and the combination of O and X mode is overwritten.
  if(modes == 3) then
    do idiag = 1, ant%N_diag
      do ich = 1, ant%diag(idiag)%N_ch
        rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%mode(1)%s_res
        rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%mode(1)%R_res
        rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%mode(1)%z_res
        rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
      end do
    end do
  end if
end if
if(present(rhop_res)) then
  do idiag = 1, ant%N_diag
    rhop_res = rad%diag(idiag)%ch(:)%rhop_res
  end do
end if
end subroutine make_rays_ecfm

subroutine make_dat_model_ece_ecfm_refr(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                        ne_rhop_scal, reflec_X_new, & ! in
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, tau, set_grid_dynamic, verbose)
use mod_ecfm_refr_types,        only: reflec_X, reflec_O, plasma_params, rad, ant, ray_init, static_grid
use mod_ecfm_refr_utils,        only: retrieve_T_e
implicit none
real(rkind), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
real(rkind),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
real(rkind), dimension(:), intent(in), optional  :: n_e_dx2, T_e_dx2
logical,     dimension(:), intent(in)  :: ece_fm_flag_ch
real(rkind), dimension(:), intent(out) :: dat_model_ece
real(rkind), dimension(:), intent(out), optional  :: tau
logical,      intent(in), optional     :: verbose
logical, intent(in), optional          :: set_grid_dynamic
integer(ikind)                         :: ich
if(.not. ray_init) then
  print*, "Something wrong with the sequencing!!"
  print*, "make_dat_model_ece_ecfm_refr was called before make_rays_ecfm"
  print*, "Critical error in mod_ecfm_refr.f90 - check stack trace!"
  call abort
end if
if(present(set_grid_dynamic)) then
  if(set_grid_dynamic) then
    static_grid = .false.
  else
    static_grid = .true.
  end if
else
  static_grid = .true.
end if
reflec_X = reflec_X_new
reflec_O = reflec_O_new
rad%diag(1)%ch(:)%eval_ch = ece_fm_flag_ch
plasma_params%rp_min = rp_min
plasma_params%rhop_scale_ne =  ne_rhop_scal
if(.not. present(T_e_dx2) .and. .not. present(n_e_dx2)) then
! Use univariate spline for both
   call update_svecs(rad, rhop_knots_ne=rhop_knots_ne, n_e=n_e, &
                     rhop_knots_Te=rhop_knots_Te, T_e=T_e)
else if(.not. present(T_e_dx2)) then
! Use IDA spline for ne but univariate spline for Te
  call update_svecs(rad, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e)
else if(.not. present(n_e_dx2)) then
! Use IDA spline for Te but univariate spline  for ne -> important ne in units of 1.e19 m^-3
   call update_svecs(rad, rhop_knots_ne, n_e, rhop_knots_Te, T_e, T_e_dx2)
else
  call update_svecs(rad, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
end if
! Perform the classical analysis using the resonance of the X-mode
if(any(ece_fm_flag_ch .eqv. .false.)) then
  call retrieve_T_e(plasma_params, abs(rad%diag(1)%ch(:)%rhop_res), dat_model_ece)
  if(any(dat_model_ece /= dat_model_ece)) then
    print*, "Nan in ECE forward model with classical analysis"
    print*, "Used rho_pol"
    print*, rad%diag(1)%ch(:)%rhop_res
    print*, dat_model_ece
    call abort
  end if
end if
if(any(ece_fm_flag_ch)) then
  call make_ece_rad_temp()
  if(any(rad%diag(1)%ch(:)%Trad /= rad%diag(1)%ch(:)%Trad)) then
    print*, "Nan in ECE forward model with forward modeled Trad"
    print*, rad%diag(1)%ch(:)%Trad
    call save_data_to_ASCII()
    call abort
  end if
end if
where (ece_fm_flag_ch) dat_model_ece = rad%diag(1)%ch(:)%Trad
where(rad%diag(1)%ch(:)%rhop_res < 0.d0) dat_model_ece = 0.d0 ! cut off
if(present(tau)) then
  if(size(tau) /= size(dat_model_ece)) then
    print*, "If provided tau must have the same shape as dat_model_ece"
    print*, "Input error in make_dat_model_ece_ecfm_refr"
    call abort()
  end if
  tau(:) = -1.d0
  where (ece_fm_flag_ch) tau = rad%diag(1)%ch(:)%tau
end if
if(present(verbose)) then
  if(verbose) call save_data_to_ASCII()
end if
!do ich =1, ant%diag(1)%N_ch
!  if(ece_fm_flag_ch(ich)) print*, rad%diag(1)%ch(ich)%Trad
!end do
end subroutine make_dat_model_ece_ecfm_refr

subroutine make_BPD_w_res_ch(idiag, ich, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                        ne_rhop_scal, reflec_X_new, & ! in
                                        reflec_O_new, rp_min, &
                                        rhop, BPD, rhop_res_warm)
use mod_ecfm_refr_types,        only: rad, ant, mode_cnt, pnts_BPD
use mod_ecfm_refr_utils,        only: retrieve_T_e
implicit none
integer(ikind), intent(in)             :: idiag, ich
real(rkind), dimension(:), intent(in)  :: rhop_knots_ne, n_e, rhop_knots_Te, T_e, n_e_dx2, T_e_dx2
real(rkind),               intent(in)  :: ne_rhop_scal, reflec_X_new, reflec_O_new, rp_min
real(rkind), dimension(:), allocatable, intent(out) :: rhop, BPD
real(rkind), intent(out)               :: rhop_res_warm
logical, dimension(ant%diag(idiag)%N_ch)  :: ece_fm_flag_ch
real(rkind), dimension(ant%diag(idiag)%N_ch) :: dat_model_ece_dummy
integer(ikind)                               :: imode
 ece_fm_flag_ch(:) = .false.
 ece_fm_flag_ch(ich) = .true.
 call make_dat_model_ece_ecfm_refr(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2, &
                                   ne_rhop_scal, reflec_X_new, & ! in
                                   reflec_O_new, ece_fm_flag_ch, rp_min, &
                                   dat_model_ece_dummy, set_grid_dynamic = .true.)
  call make_BPD_and_warm_res(idiag, ich)
  allocate(rhop(pnts_BPD))
  allocate(BPD(pnts_BPD))
  rhop = rad%diag(idiag)%ch(ich)%mode_extra_output(1)%rhop_BPD
  if(mode_cnt == 2) then
    BPD(:) = 0.d0
    do imode = 1, mode_cnt
      BPD(:) = BPD(:) + rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(:) * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
      deallocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD, &
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD)
    end do
  else
    BPD = rad%diag(idiag)%ch(ich)%mode_extra_output(1)%BPD
  end if
  deallocate(rad%diag(idiag)%ch(ich)%mode_extra_output)
  rhop_res_warm = rad%diag(idiag)%ch(ich)%rel_rhop_res
end subroutine make_BPD_w_res_ch

subroutine update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
use mod_ecfm_refr_types,               only: plasma_params, use_ida_spline_Te, use_ida_spline_ne, SOL_ne, SOL_Te
use mod_ecfm_refr_interpol,            only: make_1d_spline
use f90_kind
  implicit none
  real(rkind), dimension(:), intent(in)   :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
  real(rkind), dimension(:), intent(in), optional   :: n_e_dx2, T_e_dx2
  real(rkind), dimension(size(n_e))         :: ne_dummy
  real(rkind), dimension(size(T_e))         :: Te_dummy
  if(present(n_e_dx2)) then
  ! In the case that IDA uses e^(spline) the knot values and their second derivatives are expected to be logarithmic
    plasma_params%IDA_rhop_knots_ne = rhop_knots_ne
    plasma_params%IDA_n_e = n_e
    plasma_params%IDA_n_e_dx2 = n_e_dx2
    use_ida_spline_ne = .true.
  else
    if(plasma_params%prof_log_flag) then
    ! Remove any zeros
      ne_dummy(:) = n_e(:)
      where(ne_dummy < SOL_ne) ne_dummy = SOL_ne
      call make_1d_spline( plasma_params%ne_spline, int(size(rhop_knots_ne),4), rhop_knots_ne, log(ne_dummy * 1.d-19))
    else
      call make_1d_spline( plasma_params%ne_spline, int(size(rhop_knots_ne),4), rhop_knots_ne, n_e)
    end if
    plasma_params%ne_spline%iopt_int = 1
    use_ida_spline_ne = .false.
  end if
  if(present(T_e_dx2)) then
  ! In the case that IDA uses e^(spline) the knot values and their second derivatives are expected to be logarithmic
    plasma_params%IDA_rhop_knots_Te = rhop_knots_Te
    plasma_params%IDA_T_e = T_e
    plasma_params%IDA_T_e_dx2 = T_e_dx2
    use_ida_spline_Te = .true.
  else
    if(plasma_params%prof_log_flag) then
      ! Remove any zeros
      Te_dummy(:) = T_e(:)
      where(Te_dummy < SOL_Te) Te_dummy = SOL_Te
      call make_1d_spline( plasma_params%Te_spline, int(size(rhop_knots_Te),4), rhop_knots_Te, log(Te_dummy))
    else
      call make_1d_spline( plasma_params%Te_spline, int(size(rhop_knots_Te),4), rhop_knots_Te, T_e)
    end if
     plasma_params%Te_spline%iopt_int = 1
     use_ida_spline_Te = .false.
  end if
  plasma_params%rhop_max = min(maxval(rhop_knots_ne), maxval(rhop_knots_Te))
end subroutine update_Te_ne

subroutine update_svecs(rad, rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
  ! This routine has to be called every time before Trad is calculated.
  use mod_ecfm_refr_types,               only: rad_type, ant, plasma_params, &
                                               N_ray, N_freq, mode_cnt, stand_alone, &
                                               use_ida_spline_Te, output_level, max_points_svec
  use f90_kind
  use constants,                         only: pi,e0, mass_e, c0
  use mod_ecfm_refr_utils,               only: retrieve_n_e, retrieve_T_e
  implicit none
  real(rkind), dimension(:), intent(in)   :: rhop_knots_ne, n_e, rhop_knots_Te, T_e
  real(rkind), dimension(:), intent(in), optional   :: n_e_dx2, T_e_dx2
  type(rad_type), intent(inout)    :: rad
  integer(ikind)                                  :: idiag, last_N, grid_size, ich, ir, &
                                                     ifreq, imode, i
  real(rkind), dimension(max_points_svec)         :: x_temp, y_temp
!  if(stand_alone) then
!    print*, "There is no reason to call update_svecs in mod_raytrace.f90 in stand_alone mode"
!    stop "stand_alone = T in update_svecs"
!  end if
  if(plasma_params%Te_ne_mat) then
    print*, "The routine update_svecs does not support Te/ne matrices"
    stop "Te_ne_mat = T in update_svecs"
  end if
  ! If no second derivatives for Te and ne provided -> univariate spline
  if(.not. present(T_e_dx2) .and. .not. present(n_e_dx2)) then
     call update_Te_ne(rhop_knots_ne=rhop_knots_ne, n_e=n_e, rhop_knots_Te=rhop_knots_Te, T_e=T_e)
  else if(.not. present(T_e_dx2)) then
    call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e)
  else if(.not. present(n_e_dx2)) then
     call update_Te_ne(rhop_knots_ne, n_e, rhop_knots_Te, T_e, T_e_dx2)
  else
    call update_Te_ne(rhop_knots_ne, n_e, n_e_dx2, rhop_knots_Te, T_e, T_e_dx2)
  end if
  do idiag = 1, ant%N_diag
    do ich = 1, ant%diag(idiag)%N_ch
      if(.not. rad%diag(idiag)%ch(ich)%eval_ch) cycle ! no need to parse if we skip the channel
      do imode = 1, mode_cnt
        do ir = 1, N_ray
          do ifreq = 1, N_freq
            ! The two temporaray arrays are necessary to assure that rhop, ne and Te are well aligned in the memory
            x_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points) = &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)%rhop
            call retrieve_n_e(plasma_params, x_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points), &
                              y_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points))
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)%ne = &
              y_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)
            call retrieve_T_e(plasma_params, x_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points), &
                              y_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points))
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)%Te = &
              y_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)%Ibb = &
                y_temp(:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points) * &
                ant%diag(idiag)%ch(ich)%freq(ifreq)**2 * e0 / c0**2
          end do
        end do
      end do
    end do
  end do
end subroutine update_svecs

subroutine make_ece_rad_temp()
use mod_ecfm_refr_types,        only: dstf, reflec_X, reflec_O, mode_cnt, N_ray, N_freq, plasma_params, &
                                      rad, ant, output_level, max_points_svec, mode_conv, reflec_model, &
                                      vessel_plasma_ratio, stand_alone
use mod_ecfm_refr_rad_transp,   only: calculate_Trad
use constants,                  only: e0, c0, pi
use mod_ecfm_refr_utils,        only: binary_search, bin_ray_BPD_to_common_rhop, make_warm_res_mode, bin_freq_to_ray
use mod_ecfm_refr_raytrace,     only: reinterpolate_svec
use mod_ecfm_refr_interpol,     only: spline_1d
implicit none

real(rkind)     :: ds_small, ds_large
integer(ikind)  :: idiag, ich, imode, ifreq, ir, error
real(rkind), dimension(max_points_svec) :: tau_array, tau_secondary_array
real(rkind) :: I0_X, I0_O, I_X, I_O, T_X, T_O
character(120)  :: out_str
do idiag = 1, ant%N_diag
#ifdef OMP
  !$omp parallel private(ich, imode, ir, ifreq, error, ds_small, ds_large)  &
  !$omp          firstprivate(out_str, tau_array, tau_secondary_array) default(shared)
  !$omp do schedule(static)
#endif
  do ich = 1, ant%diag(idiag)%N_ch
    rad%diag(idiag)%ch(ich)%Trad    = 0.d0
    rad%diag(idiag)%ch(ich)%Trad_secondary   = 0.d0
    rad%diag(idiag)%ch(ich)%tau    = 0.d0
    rad%diag(idiag)%ch(ich)%tau_secondary   = 0.d0
    rad%diag(idiag)%ch(ich)%Trad = 0.d0
    if(.not. rad%diag(idiag)%ch(ich)%eval_ch) cycle
    rad%diag(idiag)%ch(ich)%s_res = 0.d0
    rad%diag(idiag)%ch(ich)%R_res = 0.d0
    rad%diag(idiag)%ch(ich)%z_res = 0.d0
    rad%diag(idiag)%ch(ich)%rhop_res = 0.d0
    do imode = 1, mode_cnt
      rad%diag(idiag)%ch(ich)%mode(imode)%Trad = 0.d0
      if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary = 0.d0
      rad%diag(idiag)%ch(ich)%mode(imode)%tau = 0.d0
      if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary = 0.d0
      ! No need to calculate what is gonna end up to be zero
      if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff .and. &
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff == 0.d0) then
         rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff = 0.d0
         if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary = 0.d0
         cycle
      else if(mode_cnt == 2 .and. .not. rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff) then
        rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary = 0.d0
      else if(mode_cnt == 2 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff) then
        rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff
        rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff_secondary
      end if
      do ir = 1, N_ray
        if(.not. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes) cycle
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad = 0.d0
        if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad_secondary = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau = 0.d0
        if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau_secondary = 0.d0
        do ifreq = 1, N_freq
          if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%use_external_pol_coeff .and. &
             rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff == 0.d0) then
             cycle
          end if
          error = 1
          ds_large = plasma_params%dist_large
          ds_small = plasma_params%dist_small
          do while(error > 0)
            ! initialization
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad =  0.d0
            if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary =  0.d0
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau =  0.d0
            if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary  =  0.d0
            if(plasma_params%on_the_fly_raytracing) then
              stop "On the fly ray tracing currently not supported"
            else if(output_level .or. (dstf /= "relamax" .and. dstf /= "Hu")) then
              call calculate_Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                     ant%diag(idiag)%ch(ich)%freq(ifreq), &
                     ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%mode, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary, &
                     tau_array, tau_secondary_array, &
                     error)
              if(error < 0) then
                print*, "An error occured while solving radiation transport"
                call calculate_Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                    ant%diag(idiag)%ch(ich)%freq(ifreq), &
                    ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                    rad%diag(idiag)%ch(ich)%mode(imode)%mode, &
                    rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary, &
                    rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary, &
                    tau_array, tau_secondary_array, &
                    error, debug = .true.)
                call abort
              end if
            else if((ifreq /= 1 .or. N_freq == 1) .and. (ir /= 1 .or. N_ray == 1)) then
             call calculate_Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                   ant%diag(idiag)%ch(ich)%freq(ifreq), &
                   ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%mode, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary, &
                   tau_array, tau_secondary_array, &
                   error)
              if(error < 0) then
                call calculate_Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                   ant%diag(idiag)%ch(ich)%freq(ifreq), &
                   ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%mode, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary, &
                   tau_array, tau_secondary_array, &
                   error)
                call abort
              end if
            else
            ! Don't calculate anything for the purely diagnostic central ray and central frequency if output_level=false
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = 0.d0
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau = 0.d0
              error = 0
            end if
            if(error == 1) then
              ds_large = ds_large / 2.d0
              ds_small = ds_small / 2.d0
              if(output_level) then
                print*, "Too low resolution on the LOS"
                print*, "Reinterpolating svec on higher resolution grid"
                print*, "New grid size (large/small) [mm]", ds_large * 1.e3, ds_small * 1.e3
                call reinterpolate_svec(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%max_points_svec_reached, &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, &
                                      ant%diag(idiag)%ch(ich)%freq(ifreq) * 2.d0 * pi, rad%diag(idiag)%ch(ich)%mode(imode)%mode, plasma_params, &
                                      ds_large, ds_small, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)
              else
                 call reinterpolate_svec(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec, &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%max_points_svec_reached, &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, &
                                      ant%diag(idiag)%ch(ich)%freq(ifreq) * 2.d0 * pi, rad%diag(idiag)%ch(ich)%mode(imode)%mode, plasma_params, &
                                      ds_large, ds_small)
              end if
            end if
          end do
          if(mode_cnt == 2 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes .and. .not. &
             rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff) rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff = rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff + &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff * ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight
          if(mode_cnt == 2 .and. output_level .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes .and.  .not. &
             rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff) rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary + &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff_secondary * ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight
          if(output_level) then
          ! Consider reflections ONLY if tau < tau_max  .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary < plasma_params%tau_max
          ! This only works if there is no mode conversion !!
            if(reflec_model == 0) then
              if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1  .and. reflec_O /= 0.d0) then
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary * (1.d0 /  &
                      (1.d0 - reflec_O * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary)))
              else if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1  .and. reflec_X /= 0.d0) then
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary * ( 1.d0 / &
                                        (1.d0 - reflec_X * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary)))
              end if
            else if(reflec_model == 1) then
           ! Should not be used!
              if(mode_cnt == 2 .and. mode_conv > 0.d0 .and. imode == mode_cnt) then
                ! In case of mode conversion we can only hande reflections once Trad of both modes is known.
                ! First just the terms without mode conversion
                ! Careful this assumes imode  = 1 -> X-mode, imode = 2 -> O-mode
                ! X-mode
                I_X = rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad_secondary
                I_O = rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad_secondary
                T_X = exp(-rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%tau_secondary)
                T_O = exp(-rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%tau_secondary)
                I0_X =  -((I_O*mode_conv + I_X*(1.d0 + vessel_plasma_ratio - reflec_O - T_O)) / &
                         (mode_conv**2 - (-1.d0 - vessel_plasma_ratio + reflec_O + T_O)*(-1.d0 - vessel_plasma_ratio + reflec_X + T_X)))
                I0_O = -((I_O + vessel_plasma_ratio*I_O + I_X*mode_conv - I_O*reflec_X - I_O*T_X) / &
                        (-1.d0 - 2*vessel_plasma_ratio - vessel_plasma_ratio**2 + reflec_O + vessel_plasma_ratio*reflec_O + &
                        mode_conv**2 + reflec_X + vessel_plasma_ratio*reflec_X - reflec_O*reflec_X + &
                        T_O + vessel_plasma_ratio*T_O - reflec_X*T_O + T_X + vessel_plasma_ratio*T_X - reflec_O*T_X - T_O*T_X))
                rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad_secondary = &
                  rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad_secondary + I0_X * T_X
                rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad_secondary = &
                  rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad_secondary + I0_O * T_O
              else if(mode_cnt == 1 .or. mode_conv <= 0.d0) then ! Only go here if not both modes considered
              ! No mode conversion -> both modes independent of each other
                if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1 .and. reflec_O > 0) then ! O-mode
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary = &
                    rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary * &
                    (1.d0 + exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary) / &
                    (1.d0 + vessel_plasma_ratio * (1.d0 - reflec_O) - exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary)))
                else if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1 .and. reflec_X > 0) then ! X-mode
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary = &
                    rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary * &
                    (1.d0 + exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary) / &
                    (1.d0 + vessel_plasma_ratio * (1.d0 - reflec_X) - exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary)))
                end if
              end if
            endif
          end if
          if(reflec_model == 0) then
          ! Mode specific wall reflections
            if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1  .and. reflec_O /= 0.d0 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad /= 0.d0 ) then
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * (1.d0 /  &
                    (1.d0 - reflec_O * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
            else if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1  .and. reflec_X /= 0.d0 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad /= 0.d0 ) then
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * ( 1.d0 / &
                                      (1.d0 - reflec_X * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
            end if
          else if(reflec_model == 1) then
          ! Should not be used
            if(mode_cnt == 2 .and. mode_conv > 0.d0 .and. imode == mode_cnt) then
              ! In case of mode conversion we can only hande reflections once Trad of both modes is known.
              ! First just the terms without mode conversion
              ! Careful this assumes imode  = 1 -> X-mode, imode = 2 -> O-mode
              ! X-mode
              I_X = rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad
              I_O = rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad
              T_X = exp(-rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%tau)
              T_O = exp(-rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%tau)
              I0_X =  -((I_O*mode_conv + I_X*(1.d0 + vessel_plasma_ratio - reflec_O - T_O)) / &
                       (mode_conv**2 - (-1.d0 - vessel_plasma_ratio + reflec_O + T_O)*(-1.d0 - vessel_plasma_ratio + reflec_X + T_X)))
              I0_O = -((I_O + vessel_plasma_ratio*I_O + I_X*mode_conv - I_O*reflec_X - I_O*T_X) / &
                      (-1.d0 - 2*vessel_plasma_ratio - vessel_plasma_ratio**2 + reflec_O + vessel_plasma_ratio*reflec_O + &
                      mode_conv**2 + reflec_X + vessel_plasma_ratio*reflec_X - reflec_O*reflec_X + &
                      T_O + vessel_plasma_ratio*T_O - reflec_X*T_O + T_X + vessel_plasma_ratio*T_X - reflec_O*T_X - T_O*T_X))
              rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad = &
                rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ifreq)%Trad + I0_X * T_X
              rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad = &
                rad%diag(idiag)%ch(ich)%mode(2)%ray(ir)%freq(ifreq)%Trad + I0_O * T_O
            else if(mode_cnt == 1 .or. mode_conv <= 0.d0) then ! Only go here if not both modes considered
            ! No mode conversion -> both modes independent of each other
              if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1 .and. reflec_O > 0) then ! O-mode
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * &
                  (1.d0 + exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau) / &
                  (1.d0 + vessel_plasma_ratio - reflec_O - exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
              else if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1 .and. reflec_X > 0) then ! X-mode
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * &
                  (1.d0 + exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)  / &
                  (1.d0 + vessel_plasma_ratio - reflec_X - exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
              end if
            end if
          end if
!          print*, "Frequency ", ifreq, "of ", N_freq, " finished - rho", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%rhop_res
!          print*, "Ray rho", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res
!          print*, "ifreq Trad", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad
!          print*, "R z, freq, ray", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%R_res, &
!            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%z_res, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res, &
!            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res, &
!            ant%diag(idiag)%ch(ich)%freq_weight(ifreq), " s_res ray", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)% s_res
!           print*, "freq weight", ant%diag(idiag)%ch(ich)%freq_weight(ifreq)
!           print*, "Trad freq", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad
!           print*, "Trad ray in code", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad
        end do !ifreq
!        do ifreq = 1, N_freq
!          if(mode_cnt == 2 .and. rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes) &
!            rad%diag(idiag)%ch(ich)%mode(imode)%X_refl = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%X_refl + &
!            rad%diag(idiag)%ch(ich)%mode(imode)%X_refl * ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * &
!            ant%diag(idiag)%ch(ich)%freq_weight(ifreq)
!          if(output_level) then
!            if(mode_cnt == 2 .and. rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1 .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes  ) then
!              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%X_refl_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%X_refl_secondary + &
!              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%X_refl_secondary * ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * &
!              ant%diag(idiag)%ch(ich)%freq_weight(ifreq)
!            end if
!          end if
!        end do! ifreq
!        print*, "Weight", ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight
!        print*, "Trad ray", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad
!       print*, "average Trad per ray and total Trad", sum(rad%diag(idiag)%ch(ich)%mode(imode)%ray(:)%Trad) / &
!         (max(1,size(rad%diag(idiag)%ch(ich)%mode(imode)%ray(:)%Trad)))
!      print*, "weights", ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, "freq_weights", ant%diag(idiag)%ch(ich)%freq_weight(:)
      end do !ir
    end do !imode
    ! Before the modes are summed up the polarization coefficients need to be renormalized due to possible numerical issues:
!    rad%diag(idiag)%ch(ich)%pol_coeff_norm = 0.d0
!    if(output_level) rad%diag(idiag)%ch(ich)%pol_coeff_secondary_norm = 0.d0
!    do imode = 1, mode_cnt
!      rad%diag(idiag)%ch(ich)%pol_coeff_norm  = rad%diag(idiag)%ch(ich)%pol_coeff_norm + rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff
!      if(output_level) rad%diag(idiag)%ch(ich)%pol_coeff_secondary_norm  = rad%diag(idiag)%ch(ich)%pol_coeff_secondary_norm + &
!                       rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary
!    end do
    !if(rad%diag(idiag)%ch(ich)%pol_coeff_norm == 0.d0) then
!    print*, "pol coeff norm:", rad%diag(idiag)%ch(ich)%pol_coeff_norm
    do imode = 1, mode_cnt
      if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff .and. &
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff == 0.d0) cycle
      do ir = 1, N_ray
        if(.not. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes) cycle
        do ifreq = 1, N_freq
        ! This Trad is to indicate the contribution of the individual modes
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad + &
            ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad
          if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad_secondary + &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau + &
            ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau
          if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau_secondary + &
            ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary
          ! This Trad is the combination of both modes
          rad%diag(idiag)%ch(ich)%Trad = rad%diag(idiag)%ch(ich)%Trad + &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff * ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * &
                ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad
          if(output_level) rad%diag(idiag)%ch(ich)%Trad_secondary = rad%diag(idiag)%ch(ich)%Trad_secondary + &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%pol_coeff_secondary * ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary
        end do ! ifreq
        rad%diag(idiag)%ch(ich)%mode(imode)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%Trad + &
          ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad
        if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary + &
          ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%Trad_secondary
        rad%diag(idiag)%ch(ich)%mode(imode)%tau = rad%diag(idiag)%ch(ich)%mode(imode)%tau + &
          ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau
        if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary + &
          ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%tau_secondary
      end do !N_ray
    end do !imode
    if(mode_cnt == 2) then
      if(mode_conv > 0.d0 .and. reflec_model == 0) then
        ! O -> X mode conversion
        ! Use average polarization coefficients to include mode conversion
        ! Gain of X-mode
          rad%diag(idiag)%ch(ich)%Trad = rad%diag(idiag)%ch(ich)%Trad + &
            rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff * mode_conv * rad%diag(idiag)%ch(ich)%mode(2)%Trad
          if(output_level) rad%diag(idiag)%ch(ich)%Trad_secondary = rad%diag(idiag)%ch(ich)%Trad_secondary + &
                           rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff_secondary * mode_conv * rad%diag(idiag)%ch(ich)%mode(2)%Trad_secondary
        ! Loss of O-mode
          rad%diag(idiag)%ch(ich)%Trad = rad%diag(idiag)%ch(ich)%Trad - &
            rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff * mode_conv * rad%diag(idiag)%ch(ich)%mode(2)%Trad
          if(output_level) rad%diag(idiag)%ch(ich)%Trad_secondary = rad%diag(idiag)%ch(ich)%Trad_secondary - &
                           rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff_secondary * mode_conv * rad%diag(idiag)%ch(ich)%mode(2)%Trad_secondary
      end if
      if(output_level) then
        !if(rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff + rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff > 1.d0) then
        !  print*, "WARNING!! Something went wrong with mode fraction in the primary model for channel:", ich
        print*, "X-mode fraction primary model", rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff
        print*, "X-mode Trad including reflections primary model", rad%diag(idiag)%ch(ich)%mode(1)%Trad
        print*, "O-mode fraction primary model", rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff
        print*, "O-mode Trad including reflections primary model", rad%diag(idiag)%ch(ich)%mode(2)%Trad
        !end if
        !if(rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff_secondary + rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff_secondary > 1.d0) then
        !  print*, "WARNING!! Something went wrong with mode fraction in the secondary model for channel:", ich
          print*, "X-mode fraction secondary model", rad%diag(idiag)%ch(ich)%mode(1)%pol_coeff_secondary
          print*, "X-mode Trad including reflections secondary model", rad%diag(idiag)%ch(ich)%mode(1)%Trad_secondary
          print*, "O-mode Trad including reflections secondary model", rad%diag(idiag)%ch(ich)%mode(2)%Trad_secondary
          print*, "O-mode fraction secondary model", rad%diag(idiag)%ch(ich)%mode(2)%pol_coeff_secondary
        !end if
      end if
    end if
    do imode = 1, mode_cnt
      if(mode_cnt == 2) then
        if(rad%diag(idiag)%ch(ich)%Trad > 0.d0) then
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac = rad%diag(idiag)%ch(ich)%mode(imode)%Trad * &
                           rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff / &
                           rad%diag(idiag)%ch(ich)%Trad
        else
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac = 0.d0 ! Because 0.d0 Trad means zero contribution of either mode
        end if
        if(rad%diag(idiag)%ch(ich)%Trad_secondary > 0.d0) then
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary * &
                           rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary / &
                           rad%diag(idiag)%ch(ich)%Trad_secondary
        else
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary = 0.d0 ! Because 0.d0 Trad means zero contribution of either mode
        end if
        ! Do not add optical depths directly -  rather add transmittances
        rad%diag(idiag)%ch(ich)%tau = rad%diag(idiag)%ch(ich)%tau + &
                                      exp(-rad%diag(idiag)%ch(ich)%mode(imode)%tau) * &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
        if(output_level) rad%diag(idiag)%ch(ich)%tau_secondary = rad%diag(idiag)%ch(ich)%tau_secondary + &
                                                                 exp(-rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary) * &
                                                                 rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary
        if(rad%diag(idiag)%ch(ich)%mode(imode)%s_res /= 0.d0) then
          rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%s_res + &
            rad%diag(idiag)%ch(ich)%mode(imode)%s_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%R_res + &
            rad%diag(idiag)%ch(ich)%mode(imode)%R_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%z_res + &
            rad%diag(idiag)%ch(ich)%mode(imode)%z_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%rhop_res + &
            rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
        end if
      else
        rad%diag(idiag)%ch(ich)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%Trad
        if(output_level) rad%diag(idiag)%ch(ich)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary
        rad%diag(idiag)%ch(ich)%tau = rad%diag(idiag)%ch(ich)%mode(imode)%tau
        if(output_level) rad%diag(idiag)%ch(ich)%tau_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary
        rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%mode(imode)%s_res
        rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%mode(imode)%R_res
        rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%mode(imode)%z_res
        rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res
      end if
    end do !imode
    if(mode_cnt == 2) then
    ! Transform back from transmittances to optical depths
      if(output_level) then
        if(rad%diag(idiag)%ch(ich)%Trad == 0.d0) then
          rad%diag(idiag)%ch(ich)%tau_secondary = 0.d0
        else
          rad%diag(idiag)%ch(ich)%tau_secondary = -log(rad%diag(idiag)%ch(ich)%tau_secondary)
        end if
      end if
      if(rad%diag(idiag)%ch(ich)%Trad == 0.d0) then
        rad%diag(idiag)%ch(ich)%tau = 0.d0
      else
        rad%diag(idiag)%ch(ich)%tau = -log(rad%diag(idiag)%ch(ich)%tau)
      end if
    end if
    if(rad%diag(idiag)%ch(ich)%s_res == 0.d0) then
    ! if both X mode and O mode Trad are zero use primary mode cold resonance (either X or O)
      rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%mode(1)%s_res
      rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%mode(1)%R_res
      rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%mode(1)%z_res
      rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
    end if
    if(output_level) then
      rad%diag(idiag)%ch(ich)%rel_s_res = 0.d0
      rad%diag(idiag)%ch(ich)%rel_R_res = 0.d0
      rad%diag(idiag)%ch(ich)%rel_z_res = 0.d0
      rad%diag(idiag)%ch(ich)%rel_rhop_res = 0.d0
      do imode = 1, mode_cnt
        rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res = 0.d0
      ! For interpretation and debug we definately want to see the contribution of each mode individually
      ! and not an inseperable superposition
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(:)%s ! unnessecceray
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(:)%R ! unnessecceray
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(:)%z ! unnessecceray
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(:)%N_cold
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec_extra_output(:)%N_cor
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec_extra_output(:)%N_warm
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD = 0.d0
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary = 0.d0
        if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff .and. &
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff == 0.d0) cycle
        do ir = 1, N_ray
          if(.not. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%contributes) cycle
          call bin_freq_to_ray(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir), &
                               ant%diag(idiag)%ch(ich)%freq_weight, &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:), &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:)%total_LOS_points)
          ! Carry out the sum over the frequencies for all rays
        end do ! ir
        call bin_ray_BPD_to_common_rhop(plasma_params, &
                                        rad%diag(idiag)%ch(ich)%mode(imode), &
                                        ant%diag(idiag)%ch(ich)%f_ECE, &
                                        ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                                        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD,  &
                                        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD, &
                                        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary)
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%Trad
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%em
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%T
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%ab
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%Trad_secondary
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%em_secondary
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%T_secondary
        rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(1)%ab_secondary
        call make_warm_res_mode(rad%diag(idiag)%ch(ich)%mode(imode), ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                                ant%diag(idiag)%ch(ich)%f_ECE, .true.)
        if(mode_cnt == 2) then
        ! warm resonances are always unique for each mode
          rad%diag(idiag)%ch(ich)%rel_s_res = rad%diag(idiag)%ch(ich)%rel_s_res + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%rel_R_res = rad%diag(idiag)%ch(ich)%rel_R_res + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
            rad%diag(idiag)%ch(ich)%rel_z_res = rad%diag(idiag)%ch(ich)%rel_z_res + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%rel_rhop_res = rad%diag(idiag)%ch(ich)%rel_rhop_res + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
          rad%diag(idiag)%ch(ich)%rel_s_res_secondary = rad%diag(idiag)%ch(ich)%rel_s_res_secondary + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res_secondary * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary
          rad%diag(idiag)%ch(ich)%rel_R_res_secondary = rad%diag(idiag)%ch(ich)%rel_R_res_secondary + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res_secondary * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary
          rad%diag(idiag)%ch(ich)%rel_z_res_secondary = rad%diag(idiag)%ch(ich)%rel_z_res_secondary + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res_secondary * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary
          rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary = rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary + &
              rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res_secondary * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac_secondary
          ! Cold resonances are unique for each mode
        else
          rad%diag(idiag)%ch(ich)%rel_s_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res
          rad%diag(idiag)%ch(ich)%rel_R_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res
          rad%diag(idiag)%ch(ich)%rel_z_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res
          rad%diag(idiag)%ch(ich)%rel_rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res
          rad%diag(idiag)%ch(ich)%rel_s_res_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res_secondary
          rad%diag(idiag)%ch(ich)%rel_R_res_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res_secondary
          rad%diag(idiag)%ch(ich)%rel_z_res_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res_secondary
          rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res_secondary
        end if
      end do !imode
      if(rad%diag(idiag)%ch(ich)%rel_s_res == 0.d0 .and. output_level) then
        ! if both X mode and O mode Trad are zero use first mode cold resonance (can be either X or O)
        rad%diag(idiag)%ch(ich)%rel_s_res = rad%diag(idiag)%ch(ich)%mode(1)%s_res
        rad%diag(idiag)%ch(ich)%rel_R_res = rad%diag(idiag)%ch(ich)%mode(1)%R_res
        rad%diag(idiag)%ch(ich)%rel_z_res = rad%diag(idiag)%ch(ich)%mode(1)%z_res
        rad%diag(idiag)%ch(ich)%rel_rhop_res = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
      end if
      if(rad%diag(idiag)%ch(ich)%rel_s_res_secondary == 0.d0 .and. output_level) then
        ! if both X mode and O mode Trad are zero use first mode cold resonance (can be either X or O)
        rad%diag(idiag)%ch(ich)%rel_s_res_secondary = rad%diag(idiag)%ch(ich)%mode(1)%s_res
        rad%diag(idiag)%ch(ich)%rel_R_res_secondary = rad%diag(idiag)%ch(ich)%mode(1)%R_res
        rad%diag(idiag)%ch(ich)%rel_z_res_secondary = rad%diag(idiag)%ch(ich)%mode(1)%z_res
        rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
      end if
      if(rad%diag(idiag)%ch(ich)%s_res == 0.d0 .and. output_level) then
        ! If no cold resonance on LOS
        print*, "WARNING: Channel:",ich, " has no cold resonance on the LOS!"
        print*, "WARNING: Using warm resonances instead"
        rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%rel_s_res
        rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%rel_R_res
        rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%rel_z_res
        rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%rel_rhop_res
      end if
    end if ! output_level
    if(output_level) then
      write(out_str, ("(A6F7.3A8F7.3A20F7.3A4)")) "rhop: ", rad%diag(idiag)%ch(ich)%rel_rhop_res,", Trad: ", rad%diag(idiag)%ch(ich)%Trad / 1.d3, &
                                    " keV, Trad (scnd):  ", rad%diag(idiag)%ch(ich)%Trad_secondary / 1.d3, " keV"
      print*,"Channel ", ich, "complete ", trim(out_str)
    end if
  end do !ich
#ifdef OMP
  !$omp end do
  !$omp end parallel
#endif
end do ! N_diag
if( stand_alone .and. output_level) call save_data_to_ASCII() !output_level .and.
!if(stand_alone .and. .not. output_level) call make_BPD_and_warm_res(1, 20)
end subroutine make_ece_rad_temp

subroutine make_BPD_and_warm_res(idiag, ich) ! to be used within IDA, calculates all birthplace distribution for one diagnostic
use mod_ecfm_refr_types,        only: plasma_params, rad, ant, mode_cnt, N_ray, &
                                      N_freq, pnts_BPD, max_points_svec
use mod_ecfm_refr_rad_transp,   only: calculate_Trad, get_em_T_fast
use constants,                  only: e0, c0, pi
use mod_ecfm_refr_utils,        only: binary_search, bin_ray_BPD_to_common_rhop, make_warm_res_mode, bin_freq_to_ray
use mod_ecfm_refr_raytrace,     only: reinterpolate_svec
implicit none
integer(ikind), intent(in) :: idiag, ich
integer(ikind)  :: imode, ifreq, ir
  allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
  do imode = 1, mode_cnt
    allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD(pnts_BPD))
    allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(pnts_BPD))
    allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
    do ir = 1, N_ray
      do ifreq = 1, N_freq
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points))
        call get_em_T_fast(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                           ant%diag(idiag)%ch(ich)%freq(ifreq), &
                           ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                           rad%diag(idiag)%ch(ich)%mode(imode)%mode)
      end do
      allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(max_points_svec), &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(max_points_svec), &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(max_points_svec), &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(max_points_svec), &
               rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(max_points_svec))
    end do
    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD = 0.d0
    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary = 0.d0
    if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff .and. &
       rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff == 0.d0) cycle
    do ir = 1, N_ray
      call bin_freq_to_ray(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir), &
                           ant%diag(idiag)%ch(ich)%freq_weight, &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:), &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:)%total_LOS_points)
      ! Carry out the sum over the frequencies for all rays
      do ifreq = 1, N_freq
        deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)
      end do
    end do ! ir
    call bin_ray_BPD_to_common_rhop(plasma_params, &
                                    rad%diag(idiag)%ch(ich)%mode(imode), &
                                    ant%diag(idiag)%ch(ich)%f_ECE, &
                                    ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                                    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD,  &
                                    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD)
    call make_warm_res_mode(rad%diag(idiag)%ch(ich)%mode(imode), ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                            ant%diag(idiag)%ch(ich)%f_ECE, .false.)
    do ir = 1, N_ray
      deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad, &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em, &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab, &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T, &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD)
    end do
    deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output)
    if(mode_cnt == 2) then
    ! warm resonances are always unique for each mode
      rad%diag(idiag)%ch(ich)%rel_s_res = rad%diag(idiag)%ch(ich)%rel_s_res + &
          rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
      rad%diag(idiag)%ch(ich)%rel_R_res = rad%diag(idiag)%ch(ich)%rel_R_res + &
          rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
        rad%diag(idiag)%ch(ich)%rel_z_res = rad%diag(idiag)%ch(ich)%rel_z_res + &
          rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
      rad%diag(idiag)%ch(ich)%rel_rhop_res = rad%diag(idiag)%ch(ich)%rel_rhop_res + &
          rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res * rad%diag(idiag)%ch(ich)%mode(imode)%Trad_mode_frac
      ! Cold resonances are unique for each mode
    else
      rad%diag(idiag)%ch(ich)%rel_s_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_s_res
      rad%diag(idiag)%ch(ich)%rel_R_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_R_res
      rad%diag(idiag)%ch(ich)%rel_z_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_z_res
      rad%diag(idiag)%ch(ich)%rel_rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%rel_rhop_res
    end if
  end do !imode
end subroutine make_BPD_and_warm_res

subroutine save_data_to_ASCII()
use mod_ecfm_refr_types,        only: mode_cnt, N_ray, N_freq, data_name, output_level, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, ray_out_folder, data_secondary_name, &
                                      output_all_ray_data, new_IO
use mod_ecfm_refr_utils,        only: export_all_ece_data
use constants,                  only: c0, e0
implicit none
Character(200)               :: filename, ich_filename, Och_filename
character(120)               :: cur_filename
character(20)                :: ich_str, ray_str, fmt_string
integer(ikind)               :: idiag, ich, imode, ifreq, ir, i, &
                                ich_tot, N_ray_output
! usually only called if output_level true
! for diagnostic purposes can also be called with output_level = false
ich_tot = 1
if(output_level) then
  filename = trim(data_folder) // data_secondary_name
  open(67, file=filename)
  filename = trim(data_folder) // "sres_rel.dat"
  open(68, file=filename)
end if
filename = trim(data_folder) // data_name
open(66, file=filename)
filename = trim(data_folder) // "sres.dat"
open(69, file=filename)
! write(filename, "(A64,A7)") data_folder, "tau.dat"
! open(67, file=filename)
! open(67, file="TRadM_s.dat")
do idiag = 1, ant%N_diag
  do ich = 1,ant%diag(idiag)%N_ch
    write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
    rad%diag(idiag)%ch(ich)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau
    if(output_level) write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
                            rad%diag(idiag)%ch(ich)%TRad_secondary/ 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau_secondary
    if(output_level) write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") &
                           rad%diag(idiag)%ch(ich)%rel_s_res, " ", rad%diag(idiag)%ch(ich)%rel_R_res, " ",&
                           rad%diag(idiag)%ch(ich)%rel_z_res, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res, " ", &
                           rad%diag(idiag)%ch(ich)%rel_s_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_R_res_secondary, " ",&
                           rad%diag(idiag)%ch(ich)%rel_z_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary
    write(69,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%s_res, " ", rad%diag(idiag)%ch(ich)%R_res, " ",&
                          rad%diag(idiag)%ch(ich)%z_res, " ", rad%diag(idiag)%ch(ich)%rhop_res
  end do
end do
close(66)
if(output_level) then
  close(67)
  close(68)
end if
close(69)
if(mode_cnt > 1 .and. output_level) then
  do imode = 1, 2
    do idiag = 1, ant%N_diag
      do ich = 1,ant%diag(idiag)%N_ch
        if(rad%diag(idiag)%ch(ich)%mode(imode)%mode > 0 .and. idiag == 1 .and. ich == 1) then
          filename = trim(data_folder) // "X_" // data_name
          open(66, file=filename)
          filename = trim(data_folder) // "X_" // data_secondary_name
          open(67, file=filename)
        else if(idiag == 1 .and. ich == 1) then
          filename = trim(data_folder) // "O_" // data_name
          open(66, file=filename)
          filename = trim(data_folder) // "O_" // data_secondary_name
          open(67, file=filename)
        end if
        write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res," ",&
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad / 1000.0d0, " ", &
          rad%diag(idiag)%ch(ich)%mode(imode)%tau, " ", rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff
        write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res," ",&
          rad%diag(idiag)%ch(ich)%mode(imode)%Trad_secondary / 1000.0d0 , " ", &
          rad%diag(idiag)%ch(ich)%mode(imode)%tau_secondary, " ", rad%diag(idiag)%ch(ich)%mode(imode)%pol_coeff_secondary
      end do ! ich
    end do ! idiag
    close(66)
    close(67)
  end do !imode
end if
if(output_level) then
  filename = trim(data_folder) // "ray_weights.dat"
  open(83, file=filename)
  filename = trim(data_folder) // "freq_weights.dat"
  open(84, file=filename)
  do idiag = 1, ant%N_diag
    do ich = 1, ant%diag(idiag)%N_ch
      write(fmt_string, "(I3)") N_ray
      fmt_string = "(" // trim(fmt_string) // "(E18.10E3))"
!      print*,fmt_string
      write(83, fmt_string) ((ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight), ir =1,N_ray)
      write(fmt_string, "(I3)") N_freq
      fmt_string = "(" // trim(fmt_string) // "(E18.10E3))"
      write(84, fmt_string) ((ant%diag(idiag)%ch(ich)%freq_weight(ifreq)), ifreq =1,N_freq)
!      print*,fmt_string
      do imode = 1, mode_cnt
        if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1 .and. output_level) then
          write(ich_str, "(I3.3)") ich_tot
          ich_filename= trim(data_folder) // Ich_name // "/Irhopch" // trim(ich_str) // ".dat"
          open(66, file=trim(ich_filename))
          ich_filename = trim(data_folder) // Ich_name // "/Trhopch" //  trim(ich_str) // ".dat"
          open(67, file=trim(ich_filename))
          ich_filename = trim(data_folder) // Ich_name // "/BPDX" //  trim(ich_str) // ".dat"
          open(70, file=trim(ich_filename))
          if(dstf_comp == "DF") then
            cur_filename = trim(data_folder) // Ich_name // "/Nch" // trim(ich_str) // "_X.dat"
            open(75, file=cur_filename)
          end if
        else if(output_level) then
          write(ich_str, "(I3.3A2)") ich_tot
          Och_filename = trim(data_folder) // Ich_name // "/IrhoOch" // trim(ich_str) // ".dat"
          open(66, file=trim(Och_filename))
          Och_filename = trim(data_folder) // Ich_name // "/TrhoOch" //  trim(ich_str) // ".dat"
          open(67, file=trim(Och_filename))
          Och_filename = trim(data_folder) // Ich_name // "/BPDO" //  trim(ich_str) // ".dat"
          open(70, file=trim(Och_filename))
          if(dstf_comp == "DF") then
            cur_filename = trim(data_folder) // Ich_name // "/Nch" // trim(ich_str) // "_O.dat"
            open(75, file=cur_filename)
          end if
        end if
        if(output_all_ray_data .eqv. .true.) then
          N_ray_output = N_ray
        else
          N_ray_output = 1
        end if
        if(output_level) then
          do ir = 1, N_ray_output
            write(ray_str,  "(I3.3)") ir
            write(ich_str,  "(I3.3)") ich_tot
            if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
              cur_filename =  trim(ray_out_folder) // "Ray" // trim(ray_str) // "ch"// trim(ich_str) // "_X.dat"
              ich_filename =   trim(data_folder) // Ich_name // "/BPD_ray" // trim(ray_str) // "ch"// trim(ich_str) // "_X.dat"
            else
              write(ray_str,  "(I3.3)") ir
              write(ich_str,  "(I3.3)") ich_tot
              cur_filename =  trim(ray_out_folder) // "Ray" // trim(ray_str) // "ch" // trim(ich_str) // "_O.dat"
              ich_filename =   trim(data_folder) // Ich_name // "/BPD_ray" // trim(ray_str) // "ch"// trim(ich_str) // "_O.dat"
            end if
            open(74, file=trim(cur_filename))
            open(98, file=trim(ich_filename))
            do i = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N
              write(74,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold(i), " ", &
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop(i), " ",&
                                          rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta(i)
            end do
            close(74)
            do i = 1,rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
            if(new_IO) then
              write(98,"(26E19.10E3)") &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%s, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(1), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(2), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(3), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%rhop, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(i), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%ne, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%Te, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(i)%H, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(i)%N_ray, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_cold, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%theta, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(1), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(2), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(3), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(1), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(2), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(3), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%v_g_perp
              else
                write(98,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%s, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(1), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(2), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%x_vec(3), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%rhop, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(i), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(i), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%Te, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(i)%H, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(i)%N_ray, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_cold, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%theta, " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(1), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(2), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%N_vec(3), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(1), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(2), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(3), " ", &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%v_g_perp
              end if
            end do
            close(98)
          end do !ir
        end if
        do i = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points
           write(66,"(E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E18.10E3)") &
                                       rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(i)%s, " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary(i)," ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em(i) - &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab(i) - &
                                       rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary(i), " ", &
                                       rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(i)%N_cold
          write(67,"(E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3)") &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(i)%s, " ",&
                                      rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T(i)," ", &
                                      rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary(i), " ", &
                                      rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(i)%Te
          if(dstf_comp == "TB".or. dstf_comp == "TO") write(75,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%svec(i)%rhop, " ", rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold(i), " ",&
                rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor(i)," ", rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm(i)
        end do ! i
        do i = 1, size(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD)
          write(70,"(E14.6E3,A1,E14.6E3,A1,E14.6E3)") &
              rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD(i), " ", &
              rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(i), " ", &
              rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary(i)
        end do
        close(66)
        close(67)
        close(70)
        close(75)
      end do ! imode
      ich_tot = ich_tot + 1
    end do ! ich
    close(83)
    close(84)
  end do ! idiag
end if
call export_all_ECE_data()
end subroutine save_data_to_ASCII
!*******************************************************************************

end module mod_ecfm_refr
