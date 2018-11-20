
module mod_ecfm_refr

use f90_kind
#ifdef OMP
  use omp_lib
#endif
public :: simulate_ida, &
          initialize_stand_alone, &
          initialize_ecfm, &
          make_rays_ecfm, &
          make_dat_model_ece_ecfm_refr, &
          make_dat_model_ece_ecfm_refr_from_Te, &
          make_ece_rad_temp

private :: save_data_to_ASCII

contains

subroutine simulate_ida(working_dir)
use mod_ecfm_refr_types,        only: n_e_filename, T_e_filename, &
                                      ant, plasma_params, rad, output_level, &
                                      Ich_name, data_folder, ray_out_folder, &
                                      N_ray, N_freq, use_ida_spline_ne, use_ida_spline_Te, &
                                      data_name, data_secondary_name
use mod_ecfm_refr_utils,        only: parse_ecfm_settings_from_ida, prepare_ece_strut, retrieve_n_e
use ece_types,                  only: ece_type
use ida_global_params,          only: ida
#ifdef NAG
USE nag_spline_1d,              only: nag_spline_1d_interp
#endif
use mod_ecfm_refr_interpol,     only: make_1d_spline, deallocate_1d_spline
implicit none
character(200), intent(in)    :: working_dir
type(ece_type)                        :: ece_strut
real(rkind), dimension(:), allocatable :: rhop, T_e, n_e, dat_model_ece, par, par_ne, par_scal, ece_rhop, ne_test, rho_pol
real(rkind)                            :: reflec_X, reflec_O, rp_min
logical,     dimension(:), allocatable :: ECE_fm_flag_ch
integer(ikind)                          :: i, j, m, n, itime, idiag, ich
character(1)                            :: sep
character(250)                          :: filename
  ida%shot = 32080
  allocate(ida%time(1))
  itime = 1
  ida%time(itime)%time = 2.57d0
  ida%time(itime)%time_beg =  ida%time(itime)%time - 5.d-4! only used for magnetic on axis
  ida%time(itime)%time_end =  ida%time(itime)%time + 5.d-4
  ida%ece%reflec_X_mode = 0.9
  ida%ece%reflec_O_mode = 0.95
  ida%ece%rhopol_scal_te = 1.0
  ida%ece%rhopol_scal_ne = 1.0
  ida%ece%rhopol_scal_te_fit = 0
  ida%ece%N_ray = 1
  ida%ece%N_freq = 1
  ida%ece%expnam = "AUGD"
  ida%ece%diag = "CEC"
  ida%ece%edition = 0
  ida%diag_eq = "EQH"
  ida%exp_eq =  "AUGD"
  ida%ed_eq = 0
  ida%btf_corr_fact_ext = 1.005d0
  ida%ece%btf_mode = "BTFABB"
  reflec_X = 0.92d0
  reflec_O = 0.95d0
  ! eq diag !
  n_e_filename = trim(working_dir) // "ecfm_data/" // "ne_file.dat"
  T_e_filename = trim(working_dir) // "ecfm_data/" // "Te_file.dat"
  data_name = "TRadM_therm.dat"
  data_secondary_name = "TRadM_TBeam.dat"
  Ich_name = "IchTB"
  data_folder = trim(working_dir) //  "ecfm_data" // "/"
  ray_out_folder = trim(working_dir) // "ecfm_data/" // "ray/"
  ! Simulate loaded ece_strut
  call prepare_ece_strut(ida, ece_strut)
  ! Make a topfile to communicate equilibrium to ecfm
  open(66, file = n_e_filename)
  read(66, "(I7.7)") m
  allocate(rhop(m), n_e(m), ne_test(m))
  do i = 1, m
    read(66,"(E19.12E2A1E19.12E2)") rhop(i), sep, n_e(i)
  end do
  close(66)
  open(66, file = T_e_filename)
  read(66, "(I7.7)") n
  if(m /= n) stop "Te and ne have to be given on the same rhop axis"
  allocate(T_e(m))
  do i = 1, m
    read(66,"(E19.12E2A1E19.12E2)") rhop(i), sep, T_e(i)
  end do
  close(66)
  ! do the spline interpolation here, since we want to simulate an IDA run, where Te and ne are provided by
  ! make_temperature and make_density  - this is only for tests!
  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
#ifdef NAG
  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
#endif
  ida%rhopol_max_spline_knot = maxval(rhop)
  rp_min = maxval(rhop)
  ! initialize_ecfm needs a rework to use the spline structure instead of rhop, T_e, n_e
  ! ideally both par and par_scale and rhop, T_e, n_e can be used
  ! output_level =  true collects the information about ray, los, birthplace distribution
  output_level = .false.
  call pre_initialize_ecfm(working_dir, ida, "init", ece_strut)
  call pre_initialize_ecfm(working_dir, ida, "clean")
  call pre_initialize_ecfm(working_dir, ida, "load")
  allocate(rho_pol(size(ece_strut%ch)))
  call initialize_ecfm(ida, itime, "init", rho_pol)
  print*, "Initialized rho_pol",rho_pol
  allocate(dat_model_ece(size(ece_strut%ch)), ece_fm_flag_ch(size(ece_strut%ch)), ece_rhop(size(ece_strut%ch)))
  ece_fm_flag_ch(:) = .false.
  use_ida_spline_ne = .false.
  use_ida_spline_Te = .false.
  allocate(par(30), par_ne(30), par_scal(30)) ! we need them allocated
  ! First two are hardcoded to be ne and Te scale
  par(1) = plasma_params%rhop_scale_ne
  par(2) = plasma_params%rhop_scale_te
  par_ne = par
  par_scal(:) = 1.d0
!  call retrieve_n_e(plasma_params, rhop(1:m-10), ne_test(1:m-10))
!  print*, "rhop", rhop(1:m-10)
!  print*, "ne", n_e(1:m-10)
!  print*, "ne _interp", ne_test(1:m-10)
!  print*, "Starting optimization"
!  print*, ida%ece%rhopol_scal_te
!  print*, ida%ece%rhopol_scal_ne
  do i =1, 3
    do j =1, 5
      call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
      call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
#ifdef NAG
      call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
      call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
#endif
      call make_dat_model_ece_ecfm_refr(par, par_ne, par_scal, reflec_X, reflec_O, & ! in
                               ece_fm_flag_ch, rp_min, dat_model_ece)
      print*, j, "-th optimization step of", i, "-th optimization"
      T_e(:) = T_e(:) * 1.01d0
      n_e(:) = n_e(:) * 1.01d0
    end do
    call make_rays_ecfm(par, par_scal, ece_rhop)
    ece_fm_flag_ch(:) = .true.
    ! update_plasma_params and span_svecs needs to be called to reinitialize the LOS
    ! another routine that should work with both rhop, T_e, n_e and par, par_scal
  end do
  print*,"Calculating chi**2"
!  call initialize_ecfm( ida, itime, "clean")
!  call pre_initialize_ecfm(working_dir, ida, "clean")
!  call make_rays_ecfm(par, par_scal, ece_rhop)
  do i = 1, 3
    call make_dat_model_ece_ecfm_refr_from_Te(par_ne, par_scal, rhop, T_e, reflec_X, reflec_O, & ! in
                                              ece_fm_flag_ch, rp_min, dat_model_ece)
    use_ida_spline_Te = .false. ! is reset everytime
    T_e(:) = T_e(:) * 1.01d0
  end do
  output_level = .true.
  call initialize_ecfm(ida, itime, "clean")
  call initialize_ecfm(ida, itime, "init")
  ! Test if everything is really deallocated by first clean and then init again
  call initialize_ecfm(ida, itime, "clean")
  call initialize_ecfm(ida, itime, "init")
  use_ida_spline_ne = .false.
  use_ida_spline_Te = .false.
  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
#ifdef NAG
  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
#endif
  call make_rays_ecfm(par, par_scal, ece_rhop)
  call initialize_ecfm(ida, itime, "clean")
  call initialize_ecfm(ida, itime, "init")
  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
#ifdef NAG
  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
#endif
  use_ida_spline_ne = .false.
  use_ida_spline_Te = .false.
  call make_rays_ecfm(par, par_scal, ece_rhop)
  ! Debug information in case of missing shot file
  print*, ant%N_diag, " diagnostics to model"
  ! multiple diagnostics possible, but currently only one for IDA
  print*, ant%diag(1)%N_ch, " channels for first diagnostics"
  print*, N_freq, " frequencies considered"
  print*, N_ray, " rays considered"
  ! makes the radiation temperature
  ! call make_ece_rad_temp()
  ! This has to be done in EVERY optimastion step
  ! Again this routine should also work with par, par_scale
  ! call update_svecs(rad, rhop, T_e, n_e)
  ! Make rad_temp again
  ! call make_ece_rad_temp()
  ! update_plasma_params and span_svecs needs to be called to reinitialize the LOS
  ! another routine that should work with both rhop, T_e, n_e and par, par_scal
  ! call update_plasma_params(plasma_params, rhop, T_e, n_e)
  ! call span_svecs(plasma_params)
  ! Make rad temp again this time with output_level = True
  ! Very slow - should never be used within optimization
  call make_dat_model_ece_ecfm_refr(par, par_ne, par_scal, reflec_X, reflec_O, & ! in
                               ece_fm_flag_ch, rp_min, dat_model_ece)
  ! Deallocate everything
  filename = trim(data_folder) // data_name
  open(66, file=filename)
  filename = trim(data_folder) // data_secondary_name
  open(67, file=filename)
  filename = trim(data_folder) // "sres_rel.dat"
  open(68, file=filename)
  filename = trim(data_folder) // "sres.dat"
  open(69, file=filename)
  ! write(filename, "(A64,A7)") data_folder, "tau.dat"
  ! open(67, file=filename)
  ! open(67, file="TRadM_s.dat")
  do idiag = 1, ant%N_diag
    do ich = 1,ant%diag(idiag)%N_ch
      write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
      rad%diag(idiag)%ch(ich)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau
      write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
        rad%diag(idiag)%ch(ich)%TRad_secondary/ 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau_secondary
      write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") &
                            rad%diag(idiag)%ch(ich)%rel_s_res, " ", rad%diag(idiag)%ch(ich)%rel_R_res, " ",&
                            rad%diag(idiag)%ch(ich)%rel_z_res, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res, " ", &
                            rad%diag(idiag)%ch(ich)%rel_s_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_R_res_secondary, " ",&
                            rad%diag(idiag)%ch(ich)%rel_z_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary
      write(69,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%s_res, " ", rad%diag(idiag)%ch(ich)%R_res, " ",&
                            rad%diag(idiag)%ch(ich)%z_res, " ", rad%diag(idiag)%ch(ich)%rhop_res
    end do
  end do
  close(66)
  close(67)
  close(68)
  close(69)
  call initialize_ecfm(ida, itime, "clean")
  call pre_initialize_ecfm(working_dir, ida, "clean")
  deallocate(rhop, n_e, T_e, par, par_ne, par_scal, ne_test)
end subroutine simulate_ida

subroutine initialize_stand_alone(working_dir, flag)
! Simulates the structure used in IDA
! While in the stand alone make_ece_rad_temp is only called once it is called many times in IDA
! Hence, to keep the structure similiar all initizalization is performed here
use mod_ecfm_refr_types,        only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, data_name, modes, &
                                      dstf_comp, plasma_params, OERT, N_ray, N_freq, warm_plasma, data_secondary_name, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, straight, stand_alone, ffp, N_absz, &
                                      N_absz_large, reflec_model
use mod_ecfm_refr_utils,      only: read_input_file, prepare_ECE_diag, &
                                    import_all_ece_data, make_ecfm_LOS_grid, &
                                    init_non_therm, read_wall_Trad
use mod_ecfm_refr_raytrace_initialize,    only: init_raytrace
use mod_ecfm_refr_raytrace,               only: span_svecs, create_svec_splines, find_cold_resonance
use mod_ecfm_refr_em_Hu,                  only: radiation_gauss_init,radiation_gauss_clean_up
use mod_ecfm_refr_abs_Al,         only: abs_Al_init,abs_Al_clean_up
use constants,                    only: pi
implicit none
character(150), intent(in)    :: working_dir
character(*),  intent(in)      :: flag
integer(ikind)                :: idiag, ich
  if(.not. stand_alone) stop "This function may only be called in stand alone"
  if(trim(flag) == "init") then
    call read_input_file(working_dir)
    call prepare_ECE_diag(working_dir=working_dir)
    dstf_comp = "Th"
    if(trim(dstf) == "Th") then
      dstf = "relamax"
      dstf_comp = "Hu"
      data_name = "TRadM_therm.dat"
      data_secondary_name = "TRadM_thrms.dat"
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
    else if(trim(dstf) == "Mx") then
      dstf = "relamax"
      data_name = "TRadM_therm.dat"
      dstf_comp = "Mx"
      data_secondary_name = "TRadM_Maxwl.dat"
      Ich_name = "IchMx"
    else if(trim(dstf) == "TB") then
      dstf = "relamax"
      dstf_comp = "TB"
      data_name = "TRadM_therm.dat"
      data_secondary_name = "TRadM_TBeam.dat"
      Ich_name = "IchTB"
    else if(trim(dstf) == "TO") then
      dstf = "relamax"
      dstf_comp = "TO"
      data_name = "TRadM_thrms.dat"
      data_secondary_name = "TRadM_TBold.dat"
      Ich_name = "IchTO"
    else if(trim(dstf) == "Al") then
      dstf = "relamax"
      dstf_comp = "Al"
      data_name = "TRadM_Maxwl.dat"
      data_secondary_name = "TRadM_oNBes.dat"
      Ich_name = "IchAl"
    else if(trim(dstf) == "OM") then
      dstf = "Hu_nbes"
      dstf_comp = "Mx"
      data_name = "TRadM_noBes.dat"
      data_secondary_name = "TRadM_olMax.dat"
      Ich_name = "IchOM"
    else if(trim(dstf) == "OB") then
      dstf = "Hu_bess"
      dstf_comp = "NB"
      data_name = "TRadM_olBes.dat"
      data_secondary_name = "TRadM_noBes.dat"
      Ich_name = "IchOB"
    else if(trim(dstf) == "O1") then
      dstf = "relamax"
      dstf_comp = "O1"
      data_name = "TRadM_therm.dat"
      data_secondary_name = "TRadM_O1mod.dat"
      Ich_name = "IchO1"
    else if(trim(dstf) == "SH") then
      dstf = "Spitzer"
      data_name = "TRadM_SpitH.dat"
      Ich_name = "IchSH"
      data_secondary_name = "TRadM_therm.dat"
    else if(trim(dstf) == "Pd") then
      dstf = "numeric"
      dst_data_folder = "fPd"
      data_name = "TRadM_Pdstl.dat"
      Ich_name = "IchPd"
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
    if(output_level) print*,"Calculating Radiation profile for ", plasma_params%shot," at t = ",plasma_params%time," s"
    if(output_level) print*,"Chosen distribution is: ", dstf
    if(output_level) flush(6)
    data_folder = trim(working_dir) //  "ecfm_data" // "/"
    if(output_level) print*, "Data will arrive in: ", data_folder
    if(dstf == "numeric" .or. trim(dstf) == "gene" .or. trim(dstf) == "gcomp") then
      call abs_Al_init(N_absz_large) ! Initializes the weights and abszissae for the gaussian quadrature
    else
      if(dstf_comp == "Hu" .or. dstf_comp == "TO" .or. dstf_comp == "Al" .or. &
         trim(dstf) == "Hu_nbes" .or. trim(dstf) == "Hu_bess") call radiation_gauss_init(N_absz)
      call abs_Al_init(N_absz) ! Initializes the weights and abszissae for the gaussian quadrature
    end if
    if(OERT) then
      if(output_level) then
        print*, "Integrated ray tracing enabled"
        print*,"Rays are preinitialized"
      end if
      ray_out_folder = trim(working_dir) // "ecfm_data/" // "ray/"
      call init_raytrace(plasma_params)
      call init_non_therm() ! Reads input data for non-thermal distributions
      if(reflec_model == 2) call read_wall_Trad()
      if(output_level) then
        print*, "----- Options for raytracing ------"
        print*, "Force straight lines of sight: ", straight
        print*, "Include ripple: ", plasma_params%w_ripple
        print*, "Use weakly relativistic cut off correction for ray tracing: ", warm_plasma
      end if
      call span_svecs(plasma_params)
      !call create_svec_splines(plasma_params)
    else
      if(output_level) print*, "Ray trajectories from external s-vectors"
      call import_all_ece_data()
      call make_ecfm_LOS_grid("initialize")
      if(N_ray /= 1 .or. N_freq /= 1 .or. modes /= 1) then
        print*, "If external s-vectors are used there must not be more than 1 ray and 1 frequency per ECE channel"
        print*, "Only X-mode is supported at this moment."
        stop "Incorrect amount of rays and/or frequencies for external s-vectors"
      end if
      do ich = 1, ant%diag(1)%N_ch
      ! Only 1 ray and 1 frequency supported
        call find_cold_resonance(plasma_params, ant%diag(1)%ch(ich)%freq(1) * 2.d0 * pi, rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1), &
                                 1, rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%total_LOS_points)
        rad%diag(1)%ch(ich)%mode(1)%ray(1)%s_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%s_res
        rad%diag(1)%ch(ich)%mode(1)%ray(1)%R_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%R_res
        rad%diag(1)%ch(ich)%mode(1)%ray(1)%z_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%z_res
        rad%diag(1)%ch(ich)%mode(1)%ray(1)%rhop_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%rhop_res
        rad%diag(1)%ch(ich)%mode(1)%s_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%s_res
        rad%diag(1)%ch(ich)%mode(1)%R_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%R_res
        rad%diag(1)%ch(ich)%mode(1)%z_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%z_res
        rad%diag(1)%ch(ich)%mode(1)%rhop_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%rhop_res
        rad%diag(1)%ch(ich)%s_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%s_res
        rad%diag(1)%ch(ich)%R_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%R_res
        rad%diag(1)%ch(ich)%z_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%z_res
        rad%diag(1)%ch(ich)%rhop_res = rad%diag(1)%ch(ich)%mode(1)%ray(1)%freq(1)%rhop_res
      end do
    end if
  else if(trim(flag) == "clean") then
    !TODO: Implement clean up routine
    if(dstf_comp == "Hu" .or. dstf_comp == "Al") call radiation_gauss_clean_up()
    call abs_Al_clean_up()
  else
    print*, "Inappropriate flag in initialize in mod_ecfm_ECE_rad"
    stop "Interal error"
  end if
end subroutine initialize_stand_alone

subroutine pre_initialize_ecfm(working_dir_in, ida, flag, ece_strut, parallelization_mode)
! Everything that is absolutely static in time is done over here

use ece_types,                  only: ece_type
use mod_ecfm_refr_types,        only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, &
                                      dstf_comp, plasma_params, OERT, N_ray, N_freq, ray_init, working_dir, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, straight, stand_alone
use mod_ecfm_refr_utils,      only: parse_ecfm_settings_from_ida, load_ECE_diag_data, &
                                    prepare_ECE_diag, dealloc_ant, &
                                    import_all_ece_data, make_ecfm_LOS_grid
use mod_ecfm_refr_em_Hu,                  only: radiation_gauss_init,radiation_gauss_clean_up
use mod_ecfm_refr_abs_Al,                only: abs_Al_init,abs_Al_clean_up
use mod_ecfm_refr_raytrace,               only: dealloc_rad
use ida_types,                  only:ida_type
implicit none
character(*), intent(in)      :: working_dir_in, flag
type(ida_type), intent(in)    :: ida
type(ece_type), intent(in), optional  :: ece_strut
integer(ikind), intent(in), optional  :: parallelization_mode
integer(ikind)                :: idiag
!call abort("ECFM disabled")
if(trim(flag) == "init" .or. trim(flag) == "load") then
! A few things need to be done both times
  working_dir = trim(working_dir_in)  // "/"
  if(present(parallelization_mode)) then
    call parse_ecfm_settings_from_ida(ida, plasma_params, parallelization_mode)
  else
    call parse_ecfm_settings_from_ida(ida, plasma_params)
  end if
  ! Hard coded since this is anyways just for testing purposes
  data_folder = trim(working_dir) //  "ecfm_data" // "/"
  ray_out_folder = trim(working_dir) // "ecfm_data/" // "ray/"
  data_folder = trim(working_dir) // "ecfm_data" // "/"
  if(output_level) print*, "ECFM data will arrive in: ", data_folder
  if((dstf == "Hu" .or. dstf_comp == "TO" .or. dstf_comp == "Al" .or. &
     dstf_comp == "Th" .or. trim(dstf) == "Hu_nbes" .or. trim(dstf) == "Hu_bess") .and. trim(dstf) /= "gene") call radiation_gauss_init(32)
  if(dstf_comp /= "TO" .and. dstf /= "Hu" .and. trim(dstf) /= "Hu_nbes" .and. trim(dstf) /= "Hu_bess") call abs_Al_init(32) ! Initializes the weights and abszissae for the gaussian quadrature
  if(trim(flag) == "init") then
    if(.not. present(ece_strut)) then
      print*, "pre_initialize_ecfm must be called with ece_strut present if flag == init"
      call abort()
    end if
    call prepare_ECE_diag(working_dir=working_dir, ida=ida, ece_strut=ece_strut)
    ! loads everything from the ida structure and then creates an input files
  else
    call load_ECE_diag_data(ida, ant, rad)
    ! load the file the routine 4 lines above create
  end if
else if(trim(flag) == "clean") then
  !TODO: Implement clean up routine
  if(dstf == "Hu") then
    call radiation_gauss_clean_up()
  else
    call abs_Al_clean_up()
  end if
  call dealloc_rad(rad)
  call dealloc_ant(ant)
  ray_init = .false.
else
  print*, "Inappropriate flag in initalize in mod_ecfm_ECE_rad"
  stop "Interal error"
end if
end subroutine pre_initialize_ecfm

subroutine initialize_ecfm(ida, itime, flag, rhopol_out)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use ece_types,                  only: ece_type
use mod_ecfm_refr_types,        only: dstf, dst_data_folder, Ich_name, ray_out_folder, output_level, &
                                      dstf_comp, plasma_params, OERT, N_ray, N_freq, ray_init, warm_plasma, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, straight, stand_alone
use mod_ecfm_refr_raytrace_initialize,    only: init_raytrace, dealloc_raytrace
use ida_types,                  only:ida_type
implicit none
character(*), intent(in)             :: flag
type(ida_type), intent(in)           :: ida
integer(ikind), intent(in)           :: itime
real(rkind), dimension(:), intent(out), optional :: rhopol_out
integer(ikind)                       :: idiag
real(rkind), dimension(1)            :: par, par_scal
  if(trim(flag) == "init") then
    plasma_params%ida_time_indx = itime
    plasma_params%time = ida%time(itime)%time
    plasma_params%time_beg = ida%time(itime)%time_beg ! only used for magnetic on axis
    plasma_params%time_end = ida%time(itime)%time_end
    call init_raytrace(plasma_params, flag = "shotfile")
    plasma_params%rhop_max = ida%rhopol_max_spline_knot
    if(present(rhopol_out)) then
      ! straight los => Neither Te nor ne matters
      plasma_params%No_ne_te = .true.
      call make_rays_ecfm(par, par_scal, rhopol_out, new_straight=.true.) ! new resonances with raytracing
      plasma_params%No_ne_te = .false.
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
  else
    print*, "Inappropriate flag in initalize in mod_ecfm_ECE_rad"
    stop "Interal error"
  end if
end subroutine initialize_ecfm

subroutine make_rays_ecfm(par, par_scal, rhop, new_straight)
! Simulates the structure used in IDA
! While in the stand alone make_ece_rad_temp is only called once it is called many times in IDA
! Hence, to keep the structure similiar all initizalization is performed here
use mod_ecfm_refr_types,        only: plasma_params, ant, rad, ray_init, straight, modes
use mod_ecfm_refr_raytrace_initialize,    only: update_plasma_params
use mod_ecfm_refr_raytrace,               only: span_svecs
use mod_ecfm_refr_utils,                  only: make_ece_rhopol_scal_te, make_ece_rhopol_scal_ne
use ece_types,                            only: ece_type
use ida_global_params,                    only: ida
implicit none
real(rkind), dimension(:), intent(in), optional    :: par, par_scal
real(rkind), dimension(:),  intent(out), optional :: rhop
logical, intent(in), optional           :: new_straight
integer(ikind)                          :: idiag, ich
! Generates new LOS based on the new Te/ne profile
if(present(new_straight)) then
  if(new_straight) then
    straight = .True.
    plasma_params%No_ne_te = .True.
  end if
end if
! Updates the values of the splines: par and par scal
if(.not. plasma_params%No_ne_te ) then
  call make_ece_rhopol_scal_te(plasma_params, par, par_scal)
  call make_ece_rhopol_scal_ne(plasma_params, par, par_scal)
  call update_plasma_params(plasma_params, par, par, par_scal)! here par and par_ne are identical
end if
call span_svecs(plasma_params)
! for second harmonic (n -1)/n =(2 -1)/2 = 0.5d0
ray_init = .true.
if(present(new_straight)) then
  if(new_straight) then
  ! For the first run of the ecfm we need resonance positions that work best for classical ECE analysis
  ! Hence, we use the resonances of the X-mode.
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
  straight = .not. ida%ece%ray_tracing
  plasma_params%No_ne_te = .False.
end if
if(present(rhop)) then
  rhop = rad%diag(1)%ch(:)%rhop_res
end if
!print*, "resonances from raytracing", rhop
end subroutine make_rays_ecfm

subroutine make_dat_model_ece_ecfm_refr(par, par_ne, par_scal, reflec_X_new, & ! in
                                        reflec_O_new, ece_fm_flag_ch, rp_min, &
                                        dat_model_ece, set_grid_dynamic)
use mod_ecfm_refr_types,        only: reflec_X, reflec_O, plasma_params, rad, ant, ray_init, static_grid
use mod_ecfm_refr_utils,               only: retrieve_T_e, make_ece_rhopol_scal_te, make_ece_rhopol_scal_ne
use mod_ecfm_refr_raytrace_initialize,    only: update_plasma_params
implicit none
real(rkind), dimension(:), intent(in)  :: par, par_ne, par_scal
real(rkind),               intent(in)  :: reflec_X_new, reflec_O_new, rp_min
logical,     dimension(:), intent(in)  :: ece_fm_flag_ch
real(rkind), dimension(:), intent(out) :: dat_model_ece
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
!    print*, "Grid is dynamic for this evaluation"
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
call update_svecs(rad, par= par, par_ne = par_ne, par_scal = par_scal)
! Perform the classical analysis using the resonance of the X-mode
if(any(ece_fm_flag_ch == .false.)) then
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
    print*, "Nan in ECE forward model with forward modelled Trad"
    call abort
  end if
end if
where (rad%diag(1)%ch(:)%Trad > 1.d-6 .and. ece_fm_flag_ch) dat_model_ece = rad%diag(1)%ch(:)%Trad
where(rad%diag(1)%ch(:)%rhop_res < 0.d0) dat_model_ece = 0.d0 ! cut off
!do ich =1, ant%diag(1)%N_ch
!  if(ece_fm_flag_ch(ich)) print*, rad%diag(1)%ch(ich)%Trad
!end do
end subroutine make_dat_model_ece_ecfm_refr

subroutine make_dat_model_ece_ecfm_refr_from_Te(par_ne, par_scal, rhop, Te, reflec_X_new, & ! in
                                                reflec_O_new, ece_fm_flag_ch, rp_min, dat_model_ece, set_grid_dynamic)
use mod_ecfm_refr_types,        only: reflec_X, reflec_O, plasma_params, rad, ray_init, ant, static_grid, use_ida_spline_Te
use mod_ecfm_refr_utils,               only: retrieve_T_e, make_ece_rhopol_scal_ne
use mod_ecfm_refr_raytrace_initialize,    only: update_plasma_params
implicit none
real(rkind), dimension(:), intent(in)  :: par_ne, par_scal, rhop, Te
real(rkind),               intent(in)  :: reflec_X_new, reflec_O_new, rp_min
logical,     dimension(:), intent(in)  :: ece_fm_flag_ch
real(rkind), dimension(:), intent(out) :: dat_model_ece
logical, intent(in), optional          :: set_grid_dynamic
real(rkind)                            :: rhop_max_old
integer(ikind)                         :: ich, i
if(.not. ray_init) then
  print*, "Something wrong with the sequencing!!"
  print*, "make_dat_model_ece_ecfm_refr was called before make_rays_ecfm"
  print*, "Critical error in mod_ecfm_refr.f90 - check stack trace!"
  call abort
end if
if(present(set_grid_dynamic)) then
  if(set_grid_dynamic) then
    static_grid = .false.
    print*, "Grid is dynamic for this evaluation"
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
! Possible bug in ida -> rhop_max changes when calculating errors
rhop_max_old = plasma_params%rhop_max
plasma_params%rhop_max = maxval(rhop)
dat_model_ece(:) = 25.d-3 ! room temperature
call update_svecs(rad, par_ne = par_ne, par_scal =  par_scal, rhop = rhop, Te = Te)
if(.not. all(ece_fm_flag_ch)) then
  use_ida_spline_Te = .false.
  do ich = 1, ant%diag(1)%N_ch
    if(rad%diag(1)%ch(ich)%rhop_res < plasma_params%rhop_max) call retrieve_T_e(plasma_params, &
      abs(rad%diag(1)%ch(ich)%rhop_res), dat_model_ece(ich))
  end do
  use_ida_spline_Te = .true.
end if
if(any(ece_fm_flag_ch)) call make_ece_rad_temp()
where (rad%diag(1)%ch(:)%Trad > 1.d-6 .and. ece_fm_flag_ch) dat_model_ece = rad%diag(1)%ch(:)%Trad
! Change back if we work again with a par instead of Te profile
plasma_params%rhop_max = rhop_max_old
!write(96, *) (dat_model_ece(ich), ich = 1, ant%diag(1)%N_ch)
end subroutine make_dat_model_ece_ecfm_refr_from_Te

subroutine update_svecs(rad, par, par_ne, par_scal, rhop, Te)
  ! This routine has to be called every time before Trad is calculated.
  ! It updates Te and ne using par, par_scal
  ! TODO: replace vectors with parameter set
  use mod_ecfm_refr_types,               only: rad_type, ant, plasma_params, &
                                               N_ray, N_freq, mode_cnt, stand_alone, &
                                               use_ida_spline_Te, output_level, max_points_svec
  use f90_kind
  use constants,                         only: pi,e0, mass_e, c0
  use mod_ecfm_refr_raytrace_initialize, only: update_plasma_params
  use mod_ecfm_refr_utils,               only: retrieve_n_e, retrieve_T_e, &
                                               make_ece_rhopol_scal_Te, make_ece_rhopol_scal_ne
#ifdef NAG
  USE nag_spline_1d,                     only: nag_spline_1d_interp
#endif
  use mod_ecfm_refr_interpol,            only: make_1d_spline
  implicit none
  type(rad_type), intent(inout)    :: rad
  real(rkind), dimension(:), intent(in), optional :: par, par_ne, par_scal, rhop, Te
  integer(ikind)                                  :: idiag, last_N, grid_size, ich, ir, &
                                                     ifreq, imode, i
  real(rkind), dimension(max_points_svec)         :: x_temp, y_temp
  if(stand_alone) then
    print*, "There is no reason to call update_svecs in mod_raytrace.f90 in stand_alone mode"
    stop "stand_alone = T in update_svecs"
  end if
  if(plasma_params%Te_ne_mat) then
    print*, "The routine update_svecs does not support Te/ne matrices"
    stop "Te_ne_mat = T in update_svecs"
  end if
  if(present(par) .and. present(par_ne) .and. present(par_scal)) then
    call update_plasma_params(plasma_params, par, par_ne, par_scal)
    call make_ece_rhopol_scal_Te(plasma_params, par, par_scal)
    call make_ece_rhopol_scal_ne(plasma_params, par, par_scal)
  else if(present(par_ne) .and. present(par_scal) .and. present(rhop) .and. present(Te)) then
    if(.not. allocated(plasma_params%par_ne)) then
      if(size(par_ne) == 0) then
        print*, "The vector par_ne containing static fit parameters for ne in ida was not allocated when update_plasma_params was called"
        stop "Critical error in update_plasma_params in mod_ecfm_refr_raytrace_initialize.f90"
      end if
      if(size(par_scal) == 0) then
        print*, "The vector par_scal containing the scaling of the fit parameters in ida was not allocated when update_plasma_params was called"
        stop "Critical error in update_plasma_params in mod_ecfm_refr_raytrace_initialize.f90"
      end if
      allocate(plasma_params%par_ne(size(par_ne)), plasma_params%par_scal(size(par_scal)))
    end if
    call make_ece_rhopol_scal_ne(plasma_params, par_ne, par_scal)
    plasma_params%par_ne = par_ne
    plasma_params%par_scal = par_scal
#ifdef NAG
    if(output_level) call nag_spline_1d_interp(rhop, Te, plasma_params%Te_spline_nag)
#endif
    call make_1d_spline( plasma_params%Te_spline, int(size(Te),4), rhop, Te, plasma_params%Te_spline%iopt_int)
    plasma_params%Te_spline%iopt_int = 1
    plasma_params%rhop_max = min(maxval(rhop), plasma_params%rhop_max)
    use_ida_spline_Te = .false.
  else
    print*, "Incorrect usage of update_svecs in mod_ecfm_refr.f90"
    print*, "Update svecs must be called with par_ne, par_scal always present and either par or rhop and Te present"
    call abort("Incorrect usage of update_svecs")
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
  if(present(par_ne) .and. present(par_scal) .and. present(rhop) .and. present(Te)) use_ida_spline_Te = .True.
end subroutine update_svecs

subroutine make_ece_rad_temp()
use mod_ecfm_refr_types,        only: dstf, reflec_X, reflec_O, OERT, mode_cnt, modes, N_ray, N_freq, plasma_params, use_maximum_for_warm_res, &
                                      rad, ant, data_folder, output_level, Ich_name, dstf_comp, max_points_svec, mode_conv, reflec_equ, reflec_model, &
                                      vessel_plasma_ratio
use mod_ecfm_refr_rad_transp,   only: calculate_Trad, calculate_Trad_LSODE
use constants,                  only: e0, c0, pi
use mod_ecfm_refr_utils,        only: binary_search, bin_ray_BPD_to_common_rhop, make_warm_res_mode, bin_freq_to_ray
use mod_ecfm_refr_raytrace,     only: reinterpolate_svec
use mod_ecfm_refr_interpol,     only: spline_1d
implicit none

real(rkind)     :: freq_2X, Trad_1O, rhop_final, sTrad, Trad_comp, max_eff_em
real(rkind)     :: ab, ab_secondary,  Trad, Trad_secondary, pol_coeff_X, pol_coeff_X_secondary, X_refl, X_refl_secondary
real(rkind), dimension(max_points_svec) :: em, em_secondary, T, T_secondary
real(rkind)     :: em_freq, ab_freq, em_secondary_freq, ab_secondary_freq, T_freq, T_secondary_freq, Trad_freq, Trad_secondary_freq, diff
real(rkind)     :: s_freq, R_freq, z_freq, rhop_freq, s, R, z, rhop, Trad_sum, max_BDOP, &
                   dummy1, dummy2, dummy3, dummy4 ! for relativsitic resonance positions
real(rkind), dimension(3) :: normal, ray_vec, x_vec
real(rkind)     :: s_prime, lin_int_res, denominator, rho_cur, max_rho, &
                   ds_small, ds_large
integer(ikind)  :: idiag, ich, imode, ifreq, ir, error
real(rkind), dimension(:,:), allocatable :: x_vec_array
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
            if(reflec_model /= 2) then
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad =  0.d0
              if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary =  0.d0
            else
              if(ant%diag(idiag)%ch(ich)%freq(ifreq) < reflec_equ%f_min .or. &
                 ant%diag(idiag)%ch(ich)%freq(ifreq) > reflec_equ%f_max) then
                 print*, "Wall plasma equilibrium wall reflection model selected (reflec_model = 1)"
                 print*, "ECE frequency not in range of f_min < f_ECE < f_max"
                 print*, "Extent wall_Trad.dat with the appropriate frequencies"
                 call abort()
              end if
              if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
#ifdef NAG
                call spline_1d(reflec_equ%X_Trad_equ_spl, &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%freq, &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, &
                               nag_spline = reflec_equ%X_Trad_equ_spl_nag)
#else
                call spline_1d(reflec_equ%X_Trad_equ_spl, &
                               ant%diag(idiag)%ch(ich)%freq(ifreq), &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad)
#endif
              else
#ifdef NAG
                call spline_1d(reflec_equ%O_Trad_equ_spl, &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%freq, &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad, &
                               nag_spline = reflec_equ%O_Trad_equ_spl_nag)
#else
                call spline_1d(reflec_equ%O_Trad_equ_spl, &
                               ant%diag(idiag)%ch(ich)%freq(ifreq), &
                               rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad)
#endif
              end if
              if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad_secondary =  &
                                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad
            end if
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau =  0.d0
            if(output_level) rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau_secondary  =  0.d0
            if(OERT .and. plasma_params%on_the_fly_raytracing) then
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
          ! Mode conversion
            if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == -1  .and. reflec_O /= 0.d0) then
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * (1.d0 /  &
                    (1.d0 - reflec_O * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
            else if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1  .and. reflec_X /= 0.d0) then
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad * ( 1.d0 / &
                                      (1.d0 - reflec_X * exp(-rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%tau)))
            end if
          else if(reflec_model == 1) then
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
        call make_warm_res_mode(plasma_params, rad%diag(idiag)%ch(ich)%mode(imode), ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                                ant%diag(idiag)%ch(ich)%f_ECE)
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
if(output_level) call save_data_to_ASCII()
end subroutine make_ece_rad_temp

subroutine make_BPD(idiag) ! to be used within IDA, calculates all birthplace distribution for one diagnostic
use mod_ecfm_refr_types,        only: dstf, reflec_X, reflec_O, OERT, mode_cnt, modes, N_ray, N_freq, plasma_params, use_maximum_for_warm_res, &
                                      rad, ant, data_folder, output_level, Ich_name, dstf_comp, max_points_svec, mode_conv
use mod_ecfm_refr_rad_transp,   only: calculate_Trad, calculate_Trad_LSODE, get_em_T_fast
use constants,                  only: e0, c0, pi
use mod_ecfm_refr_utils,        only: binary_search, bin_ray_BPD_to_common_rhop, make_warm_res_mode, bin_freq_to_ray
use mod_ecfm_refr_raytrace,     only: reinterpolate_svec
implicit none

real(rkind)     :: freq_2X, Trad_1O, rhop_final, sTrad, Trad_comp, max_eff_em
real(rkind)     :: ab, ab_secondary,  Trad, Trad_secondary, pol_coeff_X, pol_coeff_X_secondary, X_refl, X_refl_secondary
real(rkind), dimension(max_points_svec) :: em, em_secondary, T, T_secondary
real(rkind)     :: em_freq, ab_freq, em_secondary_freq, ab_secondary_freq, T_freq, T_secondary_freq, Trad_freq, Trad_secondary_freq, diff
real(rkind)     :: s_freq, R_freq, z_freq, rhop_freq, s, R, z, rhop, Trad_sum, max_BDOP, &
                   dummy1, dummy2, dummy3, dummy4 ! for relativsitic resonance positions
real(rkind), dimension(3) :: normal, ray_vec, x_vec
real(rkind)     :: s_prime, lin_int_res, denominator, rho_cur, max_rho, &
                   ds_small, ds_large
integer(ikind)  :: idiag, ich, imode, ifreq, ir, error
real(rkind), dimension(:,:), allocatable :: x_vec_array
real(rkind), dimension(max_points_svec) :: tau_array, tau_secondary_array
do ich = 1, ant%diag(idiag)%N_ch
  do imode = 1, mode_cnt
    do ir = 1, N_ray
      do ifreq = 1, N_freq
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points))
        call get_em_T_fast(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), &
                           ant%diag(idiag)%ch(ich)%freq(ifreq), &
                           ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec, &
                           rad%diag(idiag)%ch(ich)%mode(imode)%mode, &
                           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%Trad)
      end do
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
                                    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD, &
                                    rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary)
    call make_warm_res_mode(plasma_params, rad%diag(idiag)%ch(ich)%mode(imode), ant%diag(idiag)%ch(ich)%ray_launch(:)%weight, &
                            ant%diag(idiag)%ch(ich)%f_ECE)
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
end do !ich
end subroutine make_BPD

subroutine make_ece_flag_reflec(ece_flag_reflec, ida_ece_reflec, par, par_scal, ece_reflec)
! evaluates: ece_reflec
! Copied 1 : 1 from mod_ece_ecfm.f90
use f90_kind
use mod_fit_params,            only: map_par_to_ece_reflec
implicit none
integer(ikind),            intent(in)  :: ece_flag_reflec
real(rkind),               intent(in)  :: ida_ece_reflec
real(rkind), dimension(:), intent(in)  :: par, par_scal
real(rkind),               intent(out) :: ece_reflec

if (ece_flag_reflec == 0) then
  ece_reflec = 0.d0
elseif (ece_flag_reflec == 1) then
  ece_reflec = ida_ece_reflec
elseif (ece_flag_reflec == 2) then
  call map_par_to_ece_reflec(par, ece_reflec, par_scal)
else
  print*, 'ece_flag_reflec not defined properly'
endif ! (ece_flag_reflec)

!stop "sub make_ece_flag_reflec"
end subroutine make_ece_flag_reflec

subroutine save_data_to_ASCII()
use mod_ecfm_refr_types,        only: dstf, OERT, mode_cnt, N_ray, N_freq, plasma_params, data_name, &
                                      rad, ant, data_folder, Ich_name, dstf_comp, ray_out_folder, data_secondary_name, &
                                      output_all_ray_data
use mod_ecfm_refr_utils,        only: export_all_ece_data
use constants,                  only: c0, e0
implicit none
Character(200)               :: filename, ich_filename, Och_filename
character(120)               :: cur_filename
character(20)                :: ich_str, ray_str
integer(ikind)               :: idiag, ich, imode, ifreq, ir, i, ia, ia_start, &
                                ia_end, ip, ic, i_closest, N_ch_prior, ich_tot, N_ray_output
real(rkind)                  :: BPD, BPD_second
ich_tot = 1
filename = trim(data_folder) // data_name
open(66, file=filename)
filename = trim(data_folder) // data_secondary_name
open(67, file=filename)
filename = trim(data_folder) // "sres_rel.dat"
open(68, file=filename)
filename = trim(data_folder) // "sres.dat"
open(69, file=filename)
! write(filename, "(A64,A7)") data_folder, "tau.dat"
! open(67, file=filename)
! open(67, file="TRadM_s.dat")
do idiag = 1, ant%N_diag
  do ich = 1,ant%diag(idiag)%N_ch
    write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
    rad%diag(idiag)%ch(ich)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau
    write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
      rad%diag(idiag)%ch(ich)%TRad_secondary/ 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau_secondary
    write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") &
                          rad%diag(idiag)%ch(ich)%rel_s_res, " ", rad%diag(idiag)%ch(ich)%rel_R_res, " ",&
                          rad%diag(idiag)%ch(ich)%rel_z_res, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res, " ", &
                          rad%diag(idiag)%ch(ich)%rel_s_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_R_res_secondary, " ",&
                          rad%diag(idiag)%ch(ich)%rel_z_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary
    write(69,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%s_res, " ", rad%diag(idiag)%ch(ich)%R_res, " ",&
                          rad%diag(idiag)%ch(ich)%z_res, " ", rad%diag(idiag)%ch(ich)%rhop_res
  end do
end do
close(66)
close(67)
close(68)
close(69)
if(mode_cnt > 1) then
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
do idiag = 1, ant%N_diag
  do ich = 1, ant%diag(idiag)%N_ch
    do imode = 1, mode_cnt
      if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
        write(ich_str, "(I3.3)") ich_tot
        ich_filename= trim(data_folder) // Ich_name // "/Irhopch" // trim(ich_str) // ".dat"
        open(66, file=trim(ich_filename))
        ich_filename = trim(data_folder) // Ich_name // "/Trhopch" //  trim(ich_str) // ".dat"
        open(67, file=trim(ich_filename))
        ich_filename = trim(data_folder) // Ich_name // "/BPDX" //  trim(ich_str) // ".dat"
        open(70, file=trim(ich_filename))
        if(dstf_comp == "TB".or. dstf_comp == "TO") then
          cur_filename = trim(data_folder) // Ich_name // "/Nch" // trim(ich_str) // "_X.dat"
          open(75, file=cur_filename)
        end if
      else
        write(ich_str, "(I3.3A2)") ich_tot
        Och_filename = trim(data_folder) // Ich_name // "/IrhoOch" // trim(ich_str) // ".dat"
        open(66, file=trim(Och_filename))
        Och_filename = trim(data_folder) // Ich_name // "/TrhoOch" //  trim(ich_str) // ".dat"
        open(67, file=trim(Och_filename))
        Och_filename = trim(data_folder) // Ich_name // "/BPDO" //  trim(ich_str) // ".dat"
        open(70, file=trim(Och_filename))
        if(dstf_comp == "TB".or. dstf_comp == "TO") then
          cur_filename = trim(data_folder) // Ich_name // "/Nch" // trim(ich_str) // "_O.dat"
          open(75, file=cur_filename)
        end if
      end if
      if(output_all_ray_data == .true.) then
        N_ray_output = N_ray
      else
        N_ray_output = 1
      end if
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
        if(OERT) then
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
          do i = 1,rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
            write(98,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
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
                       rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec(i)%B_vec(3)
          end do
          close(74)
          close(98)
        end if
      end do !ir
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
end do ! idiag
call export_all_ECE_data()
end subroutine save_data_to_ASCII
!*******************************************************************************

end module mod_ecfm_refr
