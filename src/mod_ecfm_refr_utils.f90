!subroutine export_all_ece_data

!******************************************************************************
!******************************************************************************
!******************************************************************************

module mod_ecfm_refr_utils

use f90_kind

implicit none
public :: read_input_file, &
          prepare_ECE_diag_old_IO, &
          prepare_ECE_diag_new_IO, &
#ifdef IDA
          prepare_ECE_diag_IDA, &
          parse_ecfm_settings_from_ida, &
#endif
          dealloc_ant, &
          func_in_poly, &
          distance_to_poly, &
          binary_search, &
          BrentRoots, &
          sub_remap_coords, &
          func_calc_phi, &
          retrieve_T_e_mat_single, &
          retrieve_T_e_mat_vector, &
          retrieve_n_e_mat_single, &
          retrieve_n_e_mat_vector, &
          bin_ray_BPD_to_common_rhop, &
          bin_freq_to_ray, &
          make_warm_res_mode, &
          init_non_therm, &
          ffp_clean_up, &
          fgene_clean_up, &
          export_all_ece_data ! Writes all the data needed for the forward modelling into files

private :: make_launch, &
           func_is_left, &
           RootBracketed, &
           Minimum, &
           retrieve_T_e_single, &
           retrieve_T_e_vector, &
           retrieve_n_e_single, &
           retrieve_n_e_vector, &
           read_bi_max_data, &
           read_runaway_data, &
           read_drift_m_data, &
           read_multi_slope_data, &
           read_ffp_data, &
           read_fgene_data, &
           read_Gene_Bi_max_data, &
           read_Spitzer_data

    interface retrieve_T_e
      module procedure retrieve_T_e_single, retrieve_T_e_vector
    end interface retrieve_T_e

    interface retrieve_n_e
      module procedure retrieve_n_e_single, retrieve_n_e_vector
    end interface retrieve_n_e

#ifdef IDA
    interface
      subroutine make_temperature(par, Te, x_eval, rp_min, dTedx)
      use f90_kind
      implicit none
      real(rkind), dimension(:),   intent(in)  :: par
      real(rkind), dimension(:),   intent(out) :: Te
      real(rkind), dimension(:),   intent(in),  optional :: x_eval
      real(rkind),                 intent(in),  optional :: rp_min   ! for extrapolation in [0, rp_min]
      real(rkind), dimension(:),   intent(out), optional :: dTedx
      end subroutine make_temperature
    end interface

    interface
      subroutine make_density(par, par_scal, ne, x_eval, x_mode, dnedx)
      use f90_kind
      implicit none
      real(rkind), dimension(:),   intent(in)  :: par
      real(rkind), dimension(:),   intent(in)  :: par_scal
      real(rkind), dimension(:),   intent(out) :: ne
      real(rkind), dimension(:),   intent(in),  optional :: x_eval
      integer(ikind),              intent(in),  optional :: x_mode   ! 1/2 .. x_LiB/x_rhopol
      real(rkind), dimension(:),   intent(out), optional :: dnedx
      end subroutine make_density
    end interface

    interface
      subroutine map_par_to_ece_rp_scal_te(par, ece_rp_scal_te, par_scal)
      use f90_kind
      implicit none
      real(rkind), dimension(:), intent(in)            :: par            ! parameter to be fitted
      real(rkind),               intent(out)           :: ece_rp_scal_te ! ECE index of rhopol scale of temperature profile
      real(rkind), dimension(:), intent(in),  optional :: par_scal
      end subroutine map_par_to_ece_rp_scal_te
    end interface

    interface
      subroutine map_par_to_ece_rp_scal_ne(par, ece_rp_scal_ne, par_scal)
      use f90_kind
      implicit none
      real(rkind), dimension(:), intent(in)            :: par            ! parameter to be fitted
      real(rkind),               intent(out)           :: ece_rp_scal_ne ! ECE index of rhopol scale of density profile
      real(rkind), dimension(:), intent(in),  optional :: par_scal
      end subroutine map_par_to_ece_rp_scal_ne
    end interface
#endif

contains

subroutine read_input_file( working_dir_in)
! Reads the input file for stand-alone ECFM model
use mod_ecfm_refr_types,        only : ant, rad, plasma_params, dstf, diagnostics, output_level, &
                                       n_e_filename, T_e_filename, vessel_bd_filename, ray_launch_file, &
                                       ratio_for_third_harmonic, straight, reflec_X, reflec_O, modes, mode_cnt, working_dir, &
                                       mode_conv, N_freq, N_ray, warm_plasma, max_points_svec, &
                                       reflec_model, vessel_plasma_ratio, new_IO, use_3D
implicit none
character(*), intent(in)          :: working_dir_in
character(50)                     :: diag_str
integer(ikind)                    :: len_diag_str
character(200)                    :: input_filename
integer(ikind)                    :: IOstatus
! Temporary structure of input file - something more sophisticated later...
! shot number
! time point
! forward modelling mode dstf
! diagstr - for now just CEC or RMD permitted
! output_level (T/F)
! equilibrium experiment
! equilibrium diagnostic
! equilibrium edition
! straight
! max harmonic
! wall reflection coefficient
! which mode to consider  - 1 -> X-mode , 2 -> O-mode , 3 -> X- and O-mode
! if O-mode how much is expected to be reflected to X-mode
! number of frequencies
! number of rays - not tested and properly not correctly implemented -> 1 highly recommended
  working_dir = trim(working_dir_in)
  input_filename = trim(working_dir) // "ECRad_data/" // "ECRad.inp"
  inquire(file=input_filename, EXIST=new_IO)
  if(.not. new_IO) input_filename = trim(working_dir) // "ecfm_data/" // "ECFM.inp"
  open(66,file = trim(input_filename))
  read(66, "(A2)") dstf
  ! Only single diagnostics, ECRad should not distinguish between diagnostics
  if(.not. new_IO) then
    read(66, "(A12)") diag_str
    len_diag_str = len_trim(diag_str)
    !print*, diag_str
    if(mod(len_diag_str,3) /= 0) then
      print*, "Input Error"
      print*, "Chosen diagnostics are ", trim(diag_str)
      print*, "But this input string must have a length that is a multiple of 3"
      print*,"INPUT ERROR"
      call abort
    else
      allocate(diagnostics(len_diag_str / 3))
      ant%N_diag = 0
      if(index(trim(diag_str),"ECE") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "ECE"
      end if
      if(index(trim(diag_str),"CTC") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "CTC"
      end if
      if(index(trim(diag_str),"CTA") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "CTA"
      end if
      if(index(trim(diag_str),"IEC") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "IEC"
      end if
      if(index(trim(diag_str),"ECN") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "ECN"
      end if
      if(index(trim(diag_str),"ECO") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "ECO"
      end if
      if(index(trim(diag_str),"EXT") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "EXT"
      end if
      if(index(trim(diag_str),"VCE") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "VCE"
      end if
      if(index(trim(diag_str),"UCE") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "UCE"
      end if
      if(index(trim(diag_str),"LCE") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "LCE"
      end if
      if(index(trim(diag_str),"CCE") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "CCE"
      end if
      if(index(trim(diag_str),"REF") >= 1) then
        ant%N_diag = ant%N_diag + 1
        diagnostics(ant%N_diag) = "REF"
      end if
    end if
    allocate(ant%diag(ant%N_diag), rad%diag(ant%N_diag))
  else
    ant%N_diag = 1
    allocate(ant%diag(ant%N_diag), rad%diag(ant%N_diag), diagnostics(ant%N_diag))
    diagnostics(1) = "EXT"
  end if
  read(66,"(L1)") output_level
  if(.not. new_IO) then
    read(66,"(A4)") plasma_params%eq_exp
    plasma_params%eq_exp = trim(plasma_params%eq_exp)
    read(66,"(A3)") plasma_params%eq_diag
    plasma_params%eq_diag = trim(plasma_params%eq_diag)
    read(66,"(I2)") plasma_params%eq_ed
  end if
  read(66,"(L1)") straight
  read(66,"(L1)") plasma_params%w_ripple
  read(66,"(L1)") warm_plasma
  read(66,"(E19.12E2)") ratio_for_third_harmonic
  read(66,"(I1)") reflec_model! Obsolete has to be removed
  read(66,"(E19.12E2)") reflec_X
  read(66,"(E19.12E2)") reflec_O
  if(.not. new_IO) then
    read(66,"(E19.12E2)") vessel_plasma_ratio
  end if
  if(.not. new_IO) then
    read(66,"(E19.12E2)") plasma_params%btf_corr_fact_ext
  else
    plasma_params%btf_corr_fact_ext = 1.d0 ! Not used anyways
  end if
  read(66,"(I1)") modes
  mode_cnt = 1
  if(modes == 3) mode_cnt = 2
  read(66,"(E19.12E2)") mode_conv
  if(.not. new_IO) then
    read(66,"(E19.12E2)") plasma_params%rhop_scale_Te ! Deprecated -> remove
    read(66,"(E19.12E2)") plasma_params%rhop_scale_ne ! Deprecated -> remove
  else
    plasma_params%rhop_scale_Te = 1.0
    plasma_params%rhop_scale_ne = 1.0
  end if
  read(66,"(I4)"), N_freq ! N_freq > 1 -> consider bandwith
  if( int(real(N_freq + 1,8)/ 2.d0) /= real(N_freq + 1,8)/ 2.d0) then
    print*, "Please chose an odd N_freq in the input file"
    call abort
  end if
  read(66,"(I4)"), N_ray ! N_ray > 1 -> consider VOS instead of LOS
  if( int(sqrt(real(N_ray - 1,8))) /= sqrt(real(N_ray - 1,8))) then
    print*, "Please chose an N_ray = N**2 + 1 with N an integer in the input file"
    call abort
  end if
  read(66,"(E19.12E2)") plasma_params%dist_large
  read(66,"(E19.12E2)") plasma_params%dist_small
  read(66,"(E19.12E2)") plasma_params%R_shift ! Obsolete
  read(66,"(E19.12E2)") plasma_params%z_shift ! Obsolete
  read(66,"(I10)") max_points_svec
#ifdef USE_3D
    read(66,"(L1)",IOSTAT=IOstatus) use_3D
    if(IOstatus /= 0) use_3D = .false.
    if(use_3D .and. output_level) then
        print*, "Using 3-dimensional equilibrium"
    end if
#endif
  close(66)
  if(new_IO) then
    n_e_filename = trim(working_dir) // "ECRad_data/" // "ne_file.dat"
    T_e_filename = trim(working_dir) // "ECRad_data/" // "Te_file.dat"
    Vessel_bd_filename = trim(working_dir) // "ECRad_data/" // "vessel_bd.txt"
    ray_launch_file = trim(working_dir) // "ECRad_data/" // "ray_launch.dat"
  else
    n_e_filename = trim(working_dir) // "ecfm_data/" // "ne_file.dat"
    T_e_filename = trim(working_dir) // "ecfm_data/" // "Te_file.dat"
    Vessel_bd_filename = trim(working_dir) // "ecfm_data/" // "vessel_bd.txt"
    ray_launch_file = trim(working_dir) // "ecfm_data/" // "ray_launch.dat"
  end if
  if(output_level) then
    print*,"Chosen ECRad mode ", dstf
    if (.not. new_IO) print*, "Found ", ant%N_diag," diagnostics to model"
    if(straight) then
      print*, "NO RAYTRACING"
    else
      print*, "RAYTRACING ENABLED"
    end if
  end if

end subroutine read_input_file

#ifdef IDA
subroutine parse_ecfm_settings_from_ida(plasma_params, &
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
use mod_ecfm_refr_types,        only : plasma_params_type, ant, rad, output_level, &
                                       N_freq, N_ray, diagnostics, dstf,&
                                       straight, ratio_for_third_harmonic, reflec_X, reflec_O, warm_plasma, &
                                       modes, mode_cnt, mode_conv, vessel_bd_filename, &
                                       stand_alone, dstf_comp, use_ida_spline_Te, use_ida_spline_ne, &
                                       max_points_svec, reflec_model, data_name, data_secondary_name
implicit none
type(plasma_params_type), intent(inout)    :: plasma_params
real(rkind), intent(in)                    :: rhopol_max_spline_knot, ecrad_ratio_for_third_harmonic, &
                                              reflec_X_mode, reflec_O_mode, ecrad_O2X_mode_conversion, &
                                              rhopol_scal_te, rhopol_scal_ne, btf_corr_fact_ext, &
                                              ecrad_ds_large, ecrad_ds_small, ecrad_R_shift, ecrad_z_shift
integer(ikind), intent(in)                 :: ecrad_modes, ecrad_max_points_svec, ecrad_N_ray, &
                                              ecrad_N_freq,ece_1O_flag
logical, intent(in)                        :: ecrad_verbose, ecrad_Bt_ripple, ray_tracing, ecrad_weak_rel, log_flag
integer(ikind), intent(in), optional       :: parallelization_mode
integer(ikind)                             :: istat
  allocate(diagnostics(1))
  diagnostics(1)= "IDA"
  ant%N_diag = 1
  allocate(ant%diag(ant%N_diag), rad%diag(ant%N_diag))
  output_level = ecrad_verbose
  dstf = "Th"
  dstf_comp = "DF"
  plasma_params%eq_diag = "IDA"
  data_name = "TRadM_therm.dat"
  data_secondary_name = "TRadM_TBeam.dat"
#ifdef OMP
  if(present(parallelization_mode)) then
  ! If multiple IDA runs are done simulatenously we want no parallezation of the model
    if(parallelization_mode > 1) then
      call omp_set_num_threads(1)
    end if
  end if
#endif
  straight = .not. ray_tracing! ida%ece%ray_tracing
  ! Overwrite the settings for stand_alone usage with the settings for library usage/ ida usage
  use_ida_spline_Te = .true.
  use_ida_spline_ne = .true.
  stand_alone = .false.
  plasma_params%w_ripple = ecrad_Bt_ripple !ida%ece%ecrad_Bt_ripple
  plasma_params%rhop_max = rhopol_max_spline_knot !ida%rhopol_max_spline_knot
  plasma_params%prof_log_flag = log_flag
  warm_plasma            = ecrad_weak_rel !ida%ece%ecrad_weak_rel
  ! cut off correction incompatible with cold dispersion used in alabajar model
  ratio_for_third_harmonic = ecrad_ratio_for_third_harmonic !ida%ece%ecrad_ratio_for_third_harmonic ! (omega_c / omega = 0.4 => n = 2.5)
  reflec_X = reflec_X_mode !ida%ece%reflec_X_mode
  reflec_O = reflec_O_mode! ida%ece%reflec_O_mode
  reflec_model = 0 ! Infinite reflection model
  max_points_svec = ecrad_max_points_svec ! ida%ece%ecrad_max_points_svec
  modes =  ecrad_modes ! ida%ece%ecrad_modes ! (modes = 1 -> pure X-mode, 2 -> pure O-mode, 3 both modes and filter
  mode_cnt = 1
  if(modes == 3) mode_cnt = 2
  mode_conv = ecrad_O2X_mode_conversion ! ida%ece%ecrad_O2X_mode_conversion ! mode conversion ratio from O-X due to wall reflections
  ! Scaling of rhop axis for shifting on ne or Te
  ! Every rhop value obtained in ray tracing will be multiplied by the corresponding scaling value when evaluating Te/ne
  plasma_params%rhop_scale_Te = rhopol_scal_te ! ida%ece%rhopol_scal_te
  plasma_params%rhop_scale_ne = rhopol_scal_ne ! ida%ece%rhopol_scal_ne
  plasma_params%btf_corr_fact_ext = btf_corr_fact_ext! ida%btf_corr_fact_ext
  plasma_params%dist_large        = ecrad_ds_large ! ida%ece%ecrad_ds_large
  plasma_params%dist_small        = ecrad_ds_small ! ida%ece%ecrad_ds_small
  plasma_params%R_shift           = ecrad_R_shift! ida%ece%ecrad_R_shift    ! Allows shifting the equilbrium - moves entire flux matrix
  plasma_params%z_shift           = ecrad_z_shift! ida%ece%ecrad_z_shift    ! Allows shifting the equilbrium - moves entire flux matrix
  N_freq = ecrad_N_freq ! ida%ece%N_freq N_freq > 1 -> consider bandwith
  if( int(real(N_freq + 1,8)/ 2.d0) /= real(N_freq + 1,8)/ 2.d0) then
    print*, "Please chose an odd N_freq in the input file"
    print*, "N_freq must be odd"
    call abort
  end if
  N_ray = ecrad_N_ray ! ida%ece% N_ray > 1 -> consider VOS instead of LOS
  if( int(sqrt(real(N_ray - 1,8))) /= sqrt(real(N_ray - 1,8))) then
    print*, "Please chose an N_ray = N**2 + 1 with N an integer in the input file"
    print*, "sqrt(N_ray - 1) must be an integer"
    call abort
  end if
  ! Files used to intersect LOS with vessel wall
  vessel_bd_filename = "vessel_bd.txt"
end subroutine parse_ecfm_settings_from_ida
#endif

subroutine prepare_ECE_diag_old_IO()
! Reads launch information for several diagnostics (old file format) and allocates everything
  use mod_ecfm_refr_types,            only: ant, rad, &
                                          max_points_svec, mode_cnt, modes, N_ray, N_freq, &
                                          diagnostics, ray_launch_file, &
                                          output_level, stand_alone, &
                                          working_dir
  use constants,                      only: pi
  integer(ikind)                      :: idiag, ich, imode, ir, ifreq
  character(200)                      :: cur_filename
  character(1)                        :: sep
  open(66,file = trim(ray_launch_file))
  write(66,"(A122)") "f [GHz]        df [GHz]      x [cm]        y [cm]        z [cm]        tor [deg]     pol [deg]     dist foc[cm]  width [cm]"
  do idiag = 1, ant%N_diag ! N_diag should always be one ECRad does not distinguish between different ECE diagnostics
    ant%diag(idiag)%diag_name = diagnostics(idiag)
    cur_filename = trim(working_dir) // "ecfm_data/" // trim(ant%diag(idiag)%diag_name) // "_launch.dat"
    if(output_level .and. stand_alone) print*, "Reading ", ant%diag(idiag)%diag_name, " from external file: ", cur_filename
    open(77, file=trim(cur_filename))
    read(77, "(I5.5)") ant%diag(idiag)%N_ch
    if(trim(ant%diag(idiag)%diag_name) == "CTA" .or. &
       trim(ant%diag(idiag)%diag_name) == "CTC" .or. &
       trim(ant%diag(idiag)%diag_name) == "IEC" .or. &
       trim(ant%diag(idiag)%diag_name) == "EXT") then
      open(78, file = trim(working_dir) // "ecfm_data/" // trim(ant%diag(idiag)%diag_name) // "_pol_coeff.dat")
    end if
    allocate(ant%diag(idiag)%ch(ant%diag(idiag)%N_ch))
    allocate(rad%diag(idiag)%ch(ant%diag(idiag)%N_ch))
    do ich = 1, ant%diag(idiag)%N_ch
      allocate(ant%diag(idiag)%ch(ich)%freq(N_freq))
      allocate(ant%diag(idiag)%ch(ich)%freq_weight(N_freq))
      allocate(rad%diag(idiag)%ch(ich)%mode(mode_cnt))
      if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
      do imode = 1, mode_cnt
        if((imode == 2 .and. modes == 3) .or. &
            modes == 2) then
          rad%diag(idiag)%ch(ich)%mode(imode)%mode = -1 ! O-mode
        else
          rad%diag(idiag)%ch(ich)%mode(imode)%mode = +1 ! X-mode
        end if
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(N_ray))
        if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
        do ir = 1, N_ray
          allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
          do ifreq = 1, N_freq
              allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
              if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
          end do
        end do
      end do
    end do !Ich
    allocate(ant%diag(idiag)%ch(ich)%ray_launch(N_ray))
    do ich = 1, ant%diag(idiag)%N_ch
    ! Launching is mode independent - read this only once
      read(77, "(E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2)")  &
        ant%diag(idiag)%ch(ich)%f_ECE, sep, &
        ant%diag(idiag)%ch(ich)%df_ECE, sep, &
        ant%diag(idiag)%ch(ich)%ray_launch(1)%R, sep, &
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi, sep, &
        ant%diag(idiag)%ch(ich)%ray_launch(1)%z, sep, &
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor, sep, &
        ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol, sep, &
        ant%diag(idiag)%ch(ich)%width, sep, & ! 1/e^2 of the beam
        ant%diag(idiag)%ch(ich)%dist_focus !distance between launch and focus
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi / 180.d0 * pi
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor / 180.d0 * pi
        ! Use TORBEAM convention for input
        ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol = (ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol) / 180.d0 * pi
        rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff = -1.d0 ! Non-steerable ECE
        if(trim(ant%diag(idiag)%diag_name) == "CTA" .or. &
         trim(ant%diag(idiag)%diag_name) == "CTC" .or. &
         trim(ant%diag(idiag)%diag_name) == "IEC" .or. &
       trim(ant%diag(idiag)%diag_name) == "EXT") read(78, "(E17.10E2)") rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff
    end do
    call make_launch(idiag)
    do ich = 1, ant%diag(idiag)%N_ch
      do ir = 1, N_ray
        write(66,"(E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2A1E13.6E2)") &
              ant%diag(idiag)%ch(ich)%f_ECE, " ", &
              ant%diag(idiag)%ch(ich)%df_ECE, " ", &
              1.d2 * ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec(1), " ", &
              1.d2 * ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec(2)," ", &
              1.d2 * ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec(3), " ", &
              ant%diag(idiag)%ch(ich)%ray_launch(ir)%phi_tor / pi * 1.8d2, " ", &
              ant%diag(idiag)%ch(ich)%ray_launch(ir)%theta_pol / pi * 1.8d2, " ", &
              1.d2 * ant%diag(idiag)%ch(ich)%dist_focus, " ", &
              1.d2 * ant%diag(idiag)%ch(ich)%width
      end do !N_ray
!      do ifreq = 1, N_freq
!        write(77,"(E13.6E2A1E13.6E2)") &
!            ant%diag(idiag)%ch(ich)%freq(ifreq), " ", &
!            ant%diag(idiag)%ch(ich)%freq_weight(ifreq)
!      end do
    end do !N_ch
  end do !diag
  close(66)
end subroutine prepare_ECE_diag_old_IO

subroutine prepare_ECE_diag_new_IO()
! Reads f, df from shot files and also prepares the launching positions and angles
  use mod_ecfm_refr_types,            only: ant, rad, &
                                         	  max_points_svec, mode_cnt, modes, N_ray, N_freq, &
                                         	  ray_launch_file, output_level
  use constants,                      only: pi
  integer(ikind)                             :: idiag, ich, imode, ir, ifreq
  character(1)                               :: sep
  idiag = 1 ! Only one diag -> currently still diag array for historical reasons, should be removed eventually
  open(77, file=trim(ray_launch_file))
  ! Only one diag with new IO
  read(77, "(I5.5)") ant%diag(1)%N_ch
  allocate(rad%diag(idiag)%ch(ant%diag(idiag)%N_ch))
  allocate(ant%diag(idiag)%ch(ant%diag(1)%N_ch))
  do ich = 1, ant%diag(idiag)%N_ch
    allocate(ant%diag(idiag)%ch(ich)%freq(N_freq))
    allocate(ant%diag(idiag)%ch(ich)%freq_weight(N_freq))
    allocate(rad%diag(idiag)%ch(ich)%mode(mode_cnt))
    if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
    do imode = 1, mode_cnt
      if((imode == 2 .and. modes == 3) .or. &
          modes == 2) then
        rad%diag(idiag)%ch(ich)%mode(imode)%mode = -1 ! O-mode
      else
        rad%diag(idiag)%ch(ich)%mode(imode)%mode = +1 ! X-mode
      end if
      allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(N_ray))
      if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
      do ir = 1, N_ray
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
        do ifreq = 1, N_freq
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
            if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
        end do
      end do
    end do
    allocate(ant%diag(idiag)%ch(ich)%ray_launch(N_ray))
  end do
  do ich = 1, ant%diag(idiag)%N_ch
    read(77, "(E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2)")  &
              ant%diag(idiag)%ch(ich)%f_ECE, sep, &
              ant%diag(idiag)%ch(ich)%df_ECE, sep, &
              ant%diag(idiag)%ch(ich)%ray_launch(1)%R, sep, &
              ant%diag(idiag)%ch(ich)%ray_launch(1)%phi, sep, &
              ant%diag(idiag)%ch(ich)%ray_launch(1)%z, sep, &
              ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor, sep, &
              ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol, sep, &
              ant%diag(idiag)%ch(ich)%width, sep, & ! 1/e^2 of the beam
              ant%diag(idiag)%ch(ich)%dist_focus, sep, & !distance between launch and focus
              rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff
      ant%diag(idiag)%ch(ich)%ray_launch(1)%phi = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi / 180.d0 * pi
      ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor / 180.d0 * pi
      ! Use TORBEAM convention for input
      ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol = (ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol) / 180.d0 * pi
  end do ! ich
  close(77)
  call make_launch()
end subroutine prepare_ECE_diag_new_IO

!*******************************************************************************

#ifdef IDA
subroutine prepare_ECE_diag_IDA(working_dir, f, df, R, phi, z, tor, pol, dist_foc, width)
! Reads f, df from shot files and also prepares the launching positions and angles
  use mod_ecfm_refr_types,            only: plasma_params, ant, rad, &
                                            max_points_svec, mode_cnt, modes, N_ray, N_freq, &
                                            diagnostics, output_level, stand_alone, one_sigma_width
  use constants,                      only: pi
  character(*), intent(in)                   :: working_dir
  real(rkind), dimension(:), optional        :: f, df, R, phi, z, tor, pol, dist_foc, width ! angles in deg.
  integer(ikind)                             :: idiag, ich, imode, ir, ifreq
  character(1)                               :: sep
  character(200)                             :: ray_launch_file
  idiag = 1 ! Only one diag -> currently still diag array for historical reasons, should be removed eventually
  ray_launch_file = trim(working_dir) // "ecfm_data/" // "ray_launch.dat"
  open(77,file = trim(ray_launch_file))
  if(present(f)) then
    ant%diag(idiag)%N_ch = size(f)
    write(77, "(i5.5)")  ant%diag(idiag)%N_ch
  else
    read(77, "(I5.5)") ant%diag(1)%N_ch
  end if
  mode_cnt = 1
  if(modes == 3) mode_cnt = 2
  ! Allocate everything
  allocate(rad%diag(idiag)%ch(ant%diag(idiag)%N_ch))
  allocate(ant%diag(idiag)%ch(ant%diag(1)%N_ch))
  do ich = 1, ant%diag(idiag)%N_ch
    allocate(ant%diag(idiag)%ch(ich)%freq(N_freq))
    allocate(ant%diag(idiag)%ch(ich)%freq_weight(N_freq))
    allocate(rad%diag(idiag)%ch(ich)%mode(mode_cnt))
    if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
    do imode = 1, mode_cnt
      if((imode == 2 .and. modes == 3) .or. &
          modes == 2) then
        rad%diag(idiag)%ch(ich)%mode(imode)%mode = -1 ! O-mode
      else
        rad%diag(idiag)%ch(ich)%mode(imode)%mode = +1 ! X-mode
      end if
      allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(N_ray))
      if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
      do ir = 1, N_ray
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
        do ifreq = 1, N_freq
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
            if(output_level) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
        end do
      end do
    end do
    allocate(ant%diag(idiag)%ch(ich)%ray_launch(N_ray))
  end do
  ant%diag(idiag)%diag_name = diagnostics(idiag)
  if(present(f)) then
    do ich = 1, ant%diag(idiag)%N_ch
      ant%diag(idiag)%ch(ich)%f_ECE = f(ich)
      ant%diag(idiag)%ch(ich)%df_ECE = df(ich)
      ant%diag(idiag)%ch(ich)%ray_launch(1)%R = R(ich)
      ant%diag(idiag)%ch(ich)%ray_launch(1)%phi = phi(ich) ! angles in IDA already in rad
      ant%diag(idiag)%ch(ich)%ray_launch(1)%z = z(ich)
      ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor = tor(ich) ! angles in IDA already in rad
      ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol = pol(ich) ! angles in IDA already in rad
      ant%diag(idiag)%ch(ich)%width = width(ich)
      ant%diag(idiag)%ch(ich)%dist_focus = dist_foc(ich)
      rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff = -1.d0 ! For now hardcoded to include polarizer
      write(77, "(E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2)")  &
                ant%diag(idiag)%ch(ich)%f_ECE, " ", &
                ant%diag(idiag)%ch(ich)%df_ECE, " ", &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%R, " ", &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%phi * 180.d0 / pi, " ", &! deg.
                ant%diag(idiag)%ch(ich)%ray_launch(1)%z, " ", &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor * 180.d0 / pi, " ", &! deg.
                ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol * 180.d0 / pi, " ", &! deg.
                ant%diag(idiag)%ch(ich)%width, " ", & ! 1/e^2 of the beam
                ant%diag(idiag)%ch(ich)%dist_focus, " ", & !distance between launch and focus
                rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff
    end do
  else
    do ich = 1, ant%diag(idiag)%N_ch
      read(77, "(E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2A1E17.10E2)")  &
                ant%diag(idiag)%ch(ich)%f_ECE, sep, &
                ant%diag(idiag)%ch(ich)%df_ECE, sep, &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%R, sep, &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%phi, sep, &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%z, sep, &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor, sep, &
                ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol, sep, &
                ant%diag(idiag)%ch(ich)%width, sep, & ! 1/e^2 of the beam
                ant%diag(idiag)%ch(ich)%dist_focus, sep, & !distance between launch and focus
                rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi / 180.d0 * pi
        ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor = ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor / 180.d0 * pi
        ! Use TORBEAM convention for input
        ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol = (ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol) / 180.d0 * pi
    end do ! ich
  end if
  close(77)
  call make_launch()
end subroutine prepare_ECE_diag_IDA
#endif

!*******************************************************************************
subroutine make_launch(idiag_in)
! Reads f, df from shot files and also prepares the launching positions and angles
  use mod_ecfm_refr_types,            only: ant, rad, &
                                            mode_cnt, modes, N_ray, N_freq, one_sigma_width
  use constants,                      only: pi, c0
#ifdef NAG
  use nag_quad_util,                  only: nag_quad_gs_wt_absc
  USE nag_error_handling
#endif
  use quadrature,                     only: cgqf, cdgqf
  integer(ikind), optional                   :: idiag_in
  integer(ikind)                             :: idiag, ich, imode, ir
  integer(ikind)                             :: i, j, N_x
  real(rkind)                                :: N_abs, temp, temp1, temp2
  real(rkind), dimension(3)                  :: x_0, N_0, x1_orth, x2_orth, x_temp, R_temp, x_focus
#ifdef NAG
  real(rkind), dimension(:), allocatable     :: x, dx, x_check, dx_check, freq_check, freq_weight_check
#else
  real(rkind), dimension(:), allocatable     :: x, dx
#endif
  real(rkind)                                :: w,norm, w0, w_focus, dist_waist
#ifdef NAG
  type(nag_error)                            :: error
#endif
  if(present(idiag_in)) then
    idiag = idiag_in
  else
    idiag = 1
  end if
  if (N_ray > 1) then
    N_x = int(sqrt(real(N_ray,8) - 1.d0))
#ifdef NAG
    allocate(x(N_x), dx(N_x), x_check(N_x), dx_check(N_x))
#else
    allocate(x(N_x), dx(N_x))
#endif
  end if
#ifdef NAG
  if(N_freq > 1) allocate(freq_check(N_freq - 1), freq_weight_check(N_freq - 1))
#endif
  ! Calculate k_vector from the launching angles, distance to focues and beam waist and distribute information accross all rays
  do ich = 1, ant%diag(idiag)%N_ch
    do imode = 1,mode_cnt
      ! Distribute external polarization info
      rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff = .false.
      if(imode == 1 .and. modes /= 2) then
        if(rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff < 0.d0) cycle ! Calculate pol coeff internally
        rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%use_external_pol_coeff = .true. ! Otherwise use the externally provided one
      else if(imode == 2 .and. modes == 3) then
        if(rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff < 0.d0) cycle  ! Calculate pol coeff internally
        rad%diag(idiag)%ch(ich)%mode(2)%ray(1)%freq(1)%pol_coeff = 1.d0 - &
          rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff
        rad%diag(idiag)%ch(ich)%mode(2)%ray(1)%freq(1)%use_external_pol_coeff = .true. ! Otherwise use the externally provided one
      else if(modes == 2) then
        if(rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff < 0.d0) cycle  ! Calculate pol coeff internally
        rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff = 1.d0 - & ! Swotch from X-mode to O-mode
          rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%pol_coeff
        rad%diag(idiag)%ch(ich)%mode(1)%ray(1)%freq(1)%use_external_pol_coeff = .true. ! Otherwise use the externally provided one
      end if
      do ir = 1, N_ray
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:)%pol_coeff = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:)%pol_coeff_secondary = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%pol_coeff
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(:)%use_external_pol_coeff = rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%use_external_pol_coeff
      end do
    end do
    ant%diag(idiag)%ch(ich)%ray_launch(1)%weight = 1.d0
    ! Do the coordinate transformation for the central ray
    ! Convert position from cylindrical to spherical coordinates
    ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(1) = ant%diag(idiag)%ch(ich)%ray_launch(1)%R * cos(ant%diag(idiag)%ch(ich)%ray_launch(1)%phi)
    ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(2) = ant%diag(idiag)%ch(ich)%ray_launch(1)%R * sin(ant%diag(idiag)%ch(ich)%ray_launch(1)%phi)
    ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(3) = ant%diag(idiag)%ch(ich)%ray_launch(1)%z
    ! Launching angle depends on position - also set N_vector to position vector
    ! Also normalize it and revert its sign so it points inwards
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(:) = -ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(:) / &
        sqrt(ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(1)**2 + ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(2)**2 + &
             ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec(3)**2)
    ! Now convert to spherical polar coordinates for the rotations
    temp1 = atan2(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2), ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(1))
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(1) = 1.d0
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2) = temp1
    ! Poloidal angle is aligned with R direction
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3) = pi/2.d0
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2) = ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2) + &
                                                      ant%diag(idiag)%ch(ich)%ray_launch(1)%phi_tor
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3) = ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3) + &
                                                      ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol
    ! Back to carthesian coordinates
    temp =  cos(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2)) * sin(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3))
    temp1 =  sin(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2)) * sin(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3))
    temp2 =  cos(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3))
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(1) = temp
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2) = temp1
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3) = temp2
    N_abs = sqrt(ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(1)**2 + &
           ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(2)**2 + &
           ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec(3)**2)
    ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec = ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec / N_abs
    if(N_freq > 1) then
      if(ant%diag(idiag)%ch(ich)%df_ECE < 1.d-10) then
        print*, "Zero or negative bandwidth in calculation with finite bandwidth -- N_freq =", N_freq
        print*, "Problematic channel", ich
        print*, "Check input and rerun the program!"
        call abort()
      end if
#ifdef NAG
      CALL nag_set_error(error, halt_level=4)
      call nag_quad_gs_wt_absc( 0, ant%diag(idiag)%ch(ich)%f_ECE - ant%diag(idiag)%ch(ich)%df_ECE / 2.d0, &
                           ant%diag(idiag)%ch(ich)%f_ECE + ant%diag(idiag)%ch(ich)%df_ECE / 2.d0, &
                           freq_weight_check, &
                           freq_check, error=error)
      if (error%level >= 1) print*, error%msg
      if (error%level >= 1) then
        print*, "Nag Screwed up when setting weights for frequency discretization"
        call abort
      end if
#endif
      call cgqf(int(N_freq - 1,kind=4), int(1,kind=4), 0.d0, 0.d0, &
                ant%diag(idiag)%ch(ich)%f_ECE - ant%diag(idiag)%ch(ich)%df_ECE / 2.d0, &
                ant%diag(idiag)%ch(ich)%f_ECE + ant%diag(idiag)%ch(ich)%df_ECE / 2.d0, &
                ant%diag(idiag)%ch(ich)%freq(2:N_freq), ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq))
      ant%diag(idiag)%ch(ich)%freq(1) = ant%diag(idiag)%ch(ich)%f_ECE
      ant%diag(idiag)%ch(ich)%freq_weight(1) = 0.d0 ! central frequency not considered in the averaging
      ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq) = ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq) / &
                                                      sum(ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq)) ! normalization to one
#ifdef NAG
      freq_weight_check(:) = freq_weight_check(:) / sum(freq_weight_check(:)) ! normalization to one
      if(sum((ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq) - freq_weight_check(1:N_freq - 1))**2) > 1.d-2) then
        print*, "Weights deviate by more than 1.d-4"
        do i = 2, N_freq
          print*, ant%diag(idiag)%ch(ich)%freq_weight(i), freq_weight_check(i - 1)
        end do
        call abort()
      end if
      if(sum((ant%diag(idiag)%ch(ich)%freq(2:N_freq) - freq_check(1:N_freq - 1))**2) > 1.d-2) then
        print*, "Abszissae deviate by more than 1.d-4"
        do i = 2, N_freq
          print*, ant%diag(idiag)%ch(ich)%freq(i), freq_check(i - 1)
        end do
        call abort()
      end if
!        print*, "freq", ant%diag(idiag)%ch(ich)%freq(2:N_freq)
!        print*, "freq weight",  ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq)
!        print*, "sum",  sum(ant%diag(idiag)%ch(ich)%freq_weight(2:N_freq))
!        stop "check weights"
      !print*, Int_absz(1)
#endif
    else
      ant%diag(idiag)%ch(ich)%freq(1) = ant%diag(idiag)%ch(ich)%f_ECE
      ant%diag(idiag)%ch(ich)%freq_weight(1) = 1.d0 ! only central frequency
    end if
    if (N_ray > 1) then
      w = ant%diag(idiag)%ch(ich)%width! beam FWHM -> w
      x_0 =  ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec
      N_0 =  ant%diag(idiag)%ch(ich)%ray_launch(1)%N_vec
      ! Works only if the diagnostic can be assumed as perfectly focused
      ! The 1D heterodyne radiometer is not, hence this will be overwritten later
      ! Diagnsotics
      ! call sub_remap_coords(x_focus, R_focus)
      ! print*, "Focus point in cyl. coordinates [m] u. [rad]", R_focus
      !     print*, "R central at launch",ant%diag(idiag)%ch(ich)%ray_launch(1)%R, &
      !       ant%diag(idiag)%ch(ich)%ray_launch(1)%phi * 180.d0 / pi, &
      !       ant%diag(idiag)%ch(ich)%ray_launch(1)%z
      ! The beam needs to be spanned perpendicular to the central ray
      ! First perpendicular vector
      x1_orth(1) = 1.0
      x1_orth(2) = 0.0
      x1_orth(3) = -N_0(1)/N_0(3)
      ! Second perpendicular vector
      x2_orth(1) = - N_0(1) * N_0(2) / N_0(3)
      x2_orth(2) = N_0(3) - N_0(1)
      x2_orth(3) = - N_0(2)
      x1_orth = x1_orth/ sqrt(x1_orth(1)**2 + x1_orth(3)**2) ! need unit vectors here
      x2_orth = x2_orth/ sqrt(x2_orth(1)**2 + x2_orth(2)**2 + x2_orth(3)**2) ! need unit vectors here
      if(one_sigma_width) then
        call cgqf(int(N_x,kind=4), int(1,kind=4), 0.d0, 0.d0, -w, w, x, dx)
        norm = 1.408163862361065 / (pi * w ** 2) ! norm of the beam gaussian
      else
        call cgqf ( int(N_x,kind=4), int(1,kind=4), 0.d0, 0.d0, 0.d0, 2.d0 / (w**2), x, dx )
        norm = 2.d0 / (pi * w ** 2) ! norm of the beam gaussian
      end if
      w0 = (c0*ant%diag(idiag)%ch(ich)%dist_focus*w) / &
           Sqrt(c0**2*ant%diag(idiag)%ch(ich)%dist_focus**2 + ant%diag(idiag)%ch(ich)%f_ECE**2*Pi**2*w**4)
      dist_waist = (ant%diag(idiag)%ch(ich)%f_ECE**2*Pi**2*ant%diag(idiag)%ch(ich)%dist_focus*w**4) / &
                   (c0**2*ant%diag(idiag)%ch(ich)%dist_focus**2 + ant%diag(idiag)%ch(ich)%f_ECE**2*Pi**2*w**4)
      ! Compute gaussian beam width at focus point
      w_focus = w0*Sqrt(1.d0 + (c0**2*((ant%diag(idiag)%ch(ich)%dist_focus + ant%diag(idiag)%ch(ich)%focus_shift) - dist_waist)**2) / &
                            (ant%diag(idiag)%ch(ich)%f_ECE**2*Pi**2*w0**4))
#ifdef NAG
      call nag_quad_gs_wt_absc( 4, 0.d0, 2.d0 /( w**2), dx_check, x_check, c= 0.d0)
      if(sum((dx - dx_check)**2) > 1.d-5) then
        print*, "Weights deviate by more than 1.d-10"
        print*, "Norm of the integral over the weights"
        print*, "GPL lib", sum(dx)
        print*, "NAG lib", sum(dx_check)
        print*, " x GPL lib", x
        print*, " x NAG lib", x_check
        do i = 1, N_x
          print*, dx(i), dx_check(i)
        end do
        call abort()
      end if
      if(sum((x - x_check)**2) > 1.d-5) then
        print*, "Abszissae deviate by more than 1.d-10"
        do i = 1, N_x
          print*, x(i), x_check(i)
        end do
        call abort()
      end if
#endif
      do i = 1, N_x
        do j = 1, N_x
          x_temp(:) = (x1_orth(:) * x(i) + x2_orth(:) * x(j)) + x_0(:)
          ! Mimic the finite beam width due to diffraction by accounting for vacuum diffraction
          x_focus(:) = ant%diag(idiag)%ch(ich)%ray_launch(1)%x_vec + &
                       (ant%diag(idiag)%ch(ich)%dist_focus + ant%diag(idiag)%ch(ich)%focus_shift) * N_0
          x_focus(:) = x_focus(:) + (x1_orth(:) * x(i) + x2_orth(:) * x(j)) * w_focus/w
          call sub_remap_coords(x_temp, R_temp)
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%R = R_temp(1)
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%phi = R_temp(2)
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%z = R_temp(3)
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%weight = dx(i) * dx(j)
          if(one_sigma_width) ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%weight =  &
                                 ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%weight * exp(-(x(i)**2 + x(j)**2)/w**2)
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%x_vec = x_temp
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec = x_focus - &
            ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%x_vec
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec =  &
            ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec / &
            sqrt(ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(1)**2 + &
                 ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(2)**2 + &
                 ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(3)**2)
          if(ant%diag(idiag)%ch(ich)%dist_focus < 0) ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec = -1.0 * &
                                                     ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%phi_tor = &
            atan(ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(2) / &
                 ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(1))
          ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%theta_pol = ant%diag(idiag)%ch(ich)%ray_launch(1)%theta_pol - &
            acos(ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%N_vec(3))
          !print*, "Weight of ray", 1 + (i - 1) * N_x + j ," at launch", ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%weight * &
          !        norm
          !sum_of_weights = sum_of_weights + ant%diag(idiag)%ch(ich)%ray_launch(1 + (i - 1) * N_x + j )%weight
        end do
      end do
      ant%diag(idiag)%ch(ich)%ray_launch(:)%weight = ant%diag(idiag)%ch(ich)%ray_launch(:)%weight * &
                                                     norm ! normalization
      !sum_of_weights = sum_of_weights * norm
      ! Central ray not part of the integration - only for mapping
      ant%diag(idiag)%ch(ich)%ray_launch(1)%weight = 0.d0
      ! Renormalization
      ! Necessary, because the ray weights do not sum up to one in case of low ray count (truncation error)
      ant%diag(idiag)%ch(ich)%ray_launch(:)%weight = ant%diag(idiag)%ch(ich)%ray_launch(:)%weight / sum(ant%diag(idiag)%ch(ich)%ray_launch(:)%weight)
      !print*, "Sum of weights from ray discretisation:", sum_of_weights
      !print*, "Adjusting weights so that sum of weights gives 1"
      !ant%diag(idiag)%ch(ich)%ray_launch(:)%weight = ant%diag(idiag)%ch(ich)%ray_launch(:)%weight /  sum_of_weights
      !print*, "weights", ant%diag(idiag)%ch(ich)%ray_launch(:)%weight
    else
      ant%diag(idiag)%ch(ich)%ray_launch(1)%weight = 1.d0
    end if ! N_ray > 1
  end do ! ich
#ifdef NAG
  if(N_freq > 1) deallocate(freq_check, freq_weight_check)
  if (N_ray > 1) deallocate(x, dx, x_check, dx_check)
#endif

end subroutine make_launch
!*******************************************************************************

subroutine dealloc_ant(ant)
use mod_ecfm_refr_types,        only :  ant_type, diagnostics
implicit none
  type(ant_type), intent(inout)        :: ant
  integer(ikind)                       :: idiag, ich
  do idiag = 1, ant%N_diag
      do ich = 1, ant%diag(idiag)%N_ch
        deallocate(ant%diag(idiag)%ch(ich)%freq)
        deallocate(ant%diag(idiag)%ch(ich)%freq_weight)
        deallocate(ant%diag(idiag)%ch(ich)%ray_launch)
      end do
      deallocate(ant%diag(idiag)%ch)
  end do
  deallocate(ant%diag)
  deallocate(diagnostics)
  end subroutine dealloc_ant

  function func_is_left(x1, y1, x2, y2, x3, y3)
  use f90_kind
  implicit none
  real(rkind), intent(in) :: x1, y1, x2, y2, x3, y3
  real(rkind)                  :: func_is_left
  func_is_left = (x2 - x1) * (y3 - y1) - &
            (x3 - x1) * (y2 - y1)
  return
  end function func_is_left

  function func_in_poly(x_poly, y_poly, x, y)
  ! Check whether point inside or outside polynome using the winding number test
  ! Algorithm taken from http://geomalgorithms.com/a03-_inclusion.html
  use f90_kind
  implicit none
  real(rkind), dimension(:), intent(in)      :: x_poly, y_poly
  real(rkind), intent(in)                    :: x, y
  logical                                    :: func_in_poly
  real(rkind)                                :: point_x, point_y
  integer(ikind)                             :: wn, i
  point_x = x
  point_y = y
  wn = 0
  func_in_poly = .True.
  do i = 1, size(x_poly) - 1
    if(y_poly(i) <= y) then
        if(y_poly(i + 1) > y) then
            if (func_is_left(x_poly(i), y_poly(i), x_poly(i + 1), y_poly(i + 1), &
                             point_x, point_y) > 0) then
                wn = wn + 1
            end if
        end if
    else
        if(y_poly(i + 1) <= y) then
            if (func_is_left(x_poly(i), y_poly(i), x_poly(i + 1), y_poly(i + 1), &
                             point_x, point_y) < 0) then
                wn = wn - 1
            end if
        end if
    end if
  end do
  if(wn == 0) func_in_poly = .False.
  !print*, "Point",x, y, "is in poly: ", func_in_poly
  return
  end function func_in_poly

  function distance_to_poly(x_poly, y_poly, x, y)
  ! Computes the distance of a point to the curve spanned by poly
  ! Uses interpolation -> quite expensive
  use f90_kind
  use mod_ecfm_refr_interpol, only : make_1d_spline, spline_1d, spline_1d_get_roots, deallocate_1d_spline
  use mod_ecfm_refr_types, only : spl_type_1d
  implicit none
  real(rkind), dimension(:), intent(in)      :: x_poly, y_poly
  real(rkind), intent(in)                    :: x, y
  real(rkind)                                :: distance_to_poly
  type(spl_type_1d)                          :: spl, d_spl
  real(rkind), dimension(size(x_poly))       :: s, f
  real(rkind), dimension(500)                :: s_high_res, dummy, df, dist
  real(rkind), dimension(1000)               :: roots
  integer(ikind)                             :: i, root_cnt
  do i = 1, size(x_poly)
    s(i) = real(i - 1) / real(size(x_poly) - 1)
  end do
  do i = 1, size(s_high_res)
    s_high_res(i) = real(i - 1) / real(size(s_high_res) - 1)
  end do
  f(:) = Sqrt((x_poly - x)**2 + (y_poly - y)**2)
  call make_1d_spline(spl, size(x_poly), s, f)
  call spline_1d(spl, s, dummy, df)
  call make_1d_spline(d_spl, size(s_high_res), s_high_res, df)
  root_cnt = size(roots)
  call spline_1d_get_roots(d_spl, roots, root_cnt)
  call spline_1d(spl, roots(1:root_cnt), dist(1:root_cnt))
  distance_to_poly = minval(dist(1:root_cnt))
  call deallocate_1d_spline(spl)
  call deallocate_1d_spline(d_spl)
  return
  end function distance_to_poly

  function binary_search(array, element, N1, N2, debug)
  ! Simple binary search
  ! returns the integer so that array(i) < element
  ! array must be already sorted!
  ! N1, N2 limits the array search
  ! they are not estimates, but must be accurate
  use f90_kind
  implicit None
  real(rkind), dimension(:), intent(in) :: array
  real(rkind), intent(in)               :: element
  integer(ikind), intent(in)            :: N1, N2
  logical, intent(in), optional         :: debug
  integer(ikind)                        :: binary_search
  integer(ikind)                        :: i_split, i1, i2
  logical                               :: extra_output
  if(present(debug)) then
    extra_output = debug
  else
    extra_output = .false.
  end if
   i1 = N1
   i2 = N2
   if(array(N1) > element .or. array(N2) < element) then
    binary_search = -1
    return
   end if
    do while(.true.)
      if(extra_output) print*, "N1, N2, element, i1, i2, y1,y2 ", N1, N2, element, i1, i2, array(i1), array(i2)
      if(i2 - i1 < 2) then
        if(array(i2) == element) then
          binary_search = i2
          if(extra_output) print*, "result", i2
          if(i2 == N2) then
            if(extra_output) print*, "Value, interval", element, array(binary_search - 1), array(binary_search)
          else
            if(extra_output) print*, "Value, interval", element, array(binary_search), array(binary_search + 1)
          end if
        else if(array(i2) > element .and. array(i1) < element) then
          binary_search = i1
          if(extra_output) print*,"result", i1
          if(extra_output) print*, "Value, interval", element, array(binary_search), array(binary_search + 1)
        else if(array(i1) == element) then
          binary_search = i1
          if(extra_output)  print*,"result", i1
          if(extra_output) print*, "Value, interval", element, array(binary_search), array(binary_search + 1)
        else
          print*, "interval input and found", N1, N2, i1, i2
          print*, "Value, interval", element, array(i1), array(i2)
          binary_search = -1
          print*,"result", -1
        end if
        return
      end if
      i_split = int((real(i2,8) - real(i1,8)) / 2.d0) + 1
      if(extra_output) print*, "i_split", i_split
      if(array(i_split  + i1 - 1) == element) then
        binary_search = i_split + i1 - 1
        if(extra_output) print*,"result", i_split + i1 - 1
        if(extra_output) print*, "Value, interval", element, array(binary_search), array(binary_search + 1)
        return
      end if
      if(array(i_split + i1 - 1) > element) then
        i2 = i_split + i1 - 1
      else
        i1 = i_split + i1 - 1
      end if
    end do
  end function binary_search

  ! TRUE if x1*x2 negative
  function RootBracketed(x1,x2)
    use f90_kind
    implicit none
    real(rkind) :: x1,x2
    integer(ikind) :: resultat, RootBracketed
    if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then
      resultat = 0
    else
      resultat = 1
    endif
    RootBracketed = resultat
  end function RootBracketed

  ! returns the minimum of two real numbers
  function Minimum(x1,x2)
    use f90_kind
    implicit none
    real(rkind) :: x1,x2,resultat, Minimum
    if (x1 < x2) then
      resultat = x1
    else
      resultat = x2
    endif
    Minimum = resultat
  end function Minimum


  !*****************************************************
  !*              Brent Method Function                *
  !* ------------------------------------------------- *
  !* The purpose is to find a real root of a real      *
  !* function f(x) using Brent method.                 *
  !*                                                   *
  !* INPUTS:  x1,x2     : interval of root             *
  !*          Tolerance : desired accuracy for root    *
  !*          maxIter   : maximum number of iterations *
  !*                                                   *
  !* OUTPUTS: The function returns the root value      *
  !*          ValueAtRoot : value of f(root)           *
  !*          niter    : number of done iterations     *
  !*          error    : =0, all OK                    *
  !*                   : =1, no root found in interval *
  !*                   : =2, no more iterations !      *
  !*****************************************************
  function BrentRoots( x1, x2, Tolerance,  &
                       maxIterations, func, &
                       valueAtRoot,        &
                       niter, error )
  implicit none
    real(rkind), parameter :: FPP = 1.d-11, nearzero = 1.d-20
    real(rkind) , intent(in) :: x1,x2,Tolerance
    real(rkind), intent(out) :: valueAtRoot
    integer(ikind), intent(in) :: maxIterations
    INTERFACE
      FUNCTION func(p)
      USE f90_kind
      IMPLICIT NONE
      REAL(rkind), INTENT(IN) :: p
      REAL(rkind) :: func
      END FUNCTION func
    END INTERFACE
    integer(ikind), intent(out) :: error, niter
    real(rkind)                 :: BrentRoots
    real(rkind) :: resultat, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm
    integer(ikind) ::  i, done
    i = 0; done = 0;   error = 0
    AA = x1;  BB = x2;  FA = func(AA); FB = func(BB)
    if(FA /= FA .or. FB /= FB) then
      print*, "NaN in function BrentRoots"
      print*, "x1, x2, f(x1), f(x2)", x1, x2, FA, FB
      call abort()
    end if
    if (RootBracketed(FA,FB).eq.0) then
      error = 1
    else
      CC = BB
      FC = FB
      do while (done.eq.0.and.i < maxIterations)
        if (RootBracketed(FC,FB).eq.0) then
          CC = AA; FC = FA; DD = BB - AA; EE = DD
        endif
        if (dabs(FC) < dabs(FB)) then
          AA = BB; BB = CC; CC = AA
          FA = FB; FB = FC; FC = FA
        endif
        Tol1 = 2.0 * FPP * dabs(BB) + 0.5 * Tolerance
        xm = 0.5 * (CC-BB)
        if ((dabs(xm) <= Tol1).or.(dabs(FA) < nearzero)) then
          ! A root has been found
          resultat = BB;
          done = 1
          valueAtRoot = func(resultat)
        else
          if ((dabs(EE) >= Tol1).and.(dabs(FA) > dabs(FB))) then
            SS = FB/ FA;
            if (dabs(AA - CC) < nearzero) then
              PP = 2.0 * xm * SS;
              QQ = 1.0 - SS;
            else
              QQ = FA/FC;
              RR = FB /FC;
              PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
              QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
            endif
            if (PP > nearzero) QQ = -QQ;
            PP = dabs(PP);
            if ((2.0 * PP) < Minimum(3.0*xm *QQ-dabs(Tol1 * QQ), dabs(EE * QQ))) then
              EE = DD;  DD = PP/QQ;
            else
              DD = xm;   EE = DD;
            endif
          else
            DD = xm;
            EE = DD;
          endif
          AA = BB;
          FA = FB;
          if (dabs(DD) > Tol1) then
            BB = BB + DD;
          else
            if (xm > 0) then
        BB = BB + dabs(Tol1)
            else
        BB = BB - dabs(Tol1)
            endif
          endif
          FB = func(BB)
          i=i+1
        endif
    end do
      if (i >= maxIterations) error = 2
    endif
    niter = i
    BrentRoots = resultat
  end function BrentRoots! BrentRoots()

  subroutine sub_remap_coords(x_vec, R_vec)
  ! Simple remapping routine from LOS coords (carthesian) to tokamak coords (cylindrical)
    USE f90_kind
    USE constants, only: pi
    real(rkind), dimension(:),  intent(in)   :: x_vec
    real(rkind), dimension(:),  intent(out)  :: R_vec
    R_vec(1) = sqrt(x_vec(1)**2 + x_vec(2)**2)
    ! Use high precision implementation of atan2
    ! this is also differentiable, the downside is that we have to use an if clause
    R_vec(2) = func_calc_phi(x_vec(1),x_vec(2))
    R_vec(3) = x_vec(3)
  end subroutine sub_remap_coords

  function func_calc_phi(x, y)
  use f90_kind
  use constants, only: pi
  implicit none
    real(rkind), intent(in)    :: x, y
    real(rkind)                :: func_calc_phi
    if(y /= 0.d0) then
      func_calc_phi = 2.d0 * atan((Sqrt(x**2 + y**2) - x) / y)
    else if( x > 0.d0) then
      func_calc_phi = 0.d0
    else if(x < 0.d0) then
      func_calc_phi = pi
    else
      print*, "encountered atan(0/0)"
      call abort
    end if
  end function func_calc_phi

  subroutine retrieve_T_e_single(plasma_params, rhop, T_e, grad_T_e)
    use f90_kind
#ifdef IDA
    use mod_ecfm_refr_types,        only: plasma_params_type, use_ida_spline_Te, output_level, SOL_Te
#else
    use mod_ecfm_refr_types,        only: plasma_params_type, SOL_Te, output_level
#endif
    use mod_ecfm_refr_interpol,     only: spline_1d
    implicit none
    type(plasma_params_type), intent(in)    :: plasma_params
    real(rkind),               intent(in)   :: rhop
    real(rkind),               intent(out)  :: T_e
    real(rkind),               intent(out), optional :: grad_T_e
    real(rkind)                              :: scaled_rhop
    scaled_rhop = rhop * plasma_params%rhop_scale_te
    if(present(grad_T_e)) grad_T_e = 0.d0
    T_e = SOL_Te
    if(scaled_rhop > plasma_params%rhop_max) return
#ifdef IDA
    if(use_ida_spline_Te) then
      if(present(grad_T_e)) then
        call spline_1d(plasma_params%IDA_rhop_knots_Te, plasma_params%IDA_T_e, &
                       plasma_params%IDA_T_e_dx2, scaled_rhop, T_e, grad_T_e)
      else
        call spline_1d(plasma_params%IDA_rhop_knots_Te, plasma_params%IDA_T_e, &
                       plasma_params%IDA_T_e_dx2, scaled_rhop, T_e)
      end if
    else
#endif
      if(present(grad_T_e)) then
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%Te_spline, scaled_rhop, &
                         T_e, grad_T_e, plasma_params%Te_spline_nag)
#else
          call spline_1d(plasma_params%Te_spline, scaled_rhop, &
                         T_e, grad_T_e)
#endif
        else
          call spline_1d(plasma_params%Te_spline, scaled_rhop, &
   				               T_e, grad_T_e)
        end if
      else
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%Te_spline, scaled_rhop, &
    			               T_e, nag_spline = plasma_params%Te_spline_nag)
#else
          call spline_1d(plasma_params%Te_spline, scaled_rhop, &
         			           T_e)
#endif
        else
          call spline_1d(plasma_params%Te_spline, scaled_rhop, T_e)
        end if
      end if
#ifdef IDA
    end if
#endif
    if(plasma_params%prof_log_flag) then
      T_e = exp(T_e)
      if(present(grad_T_e)) then
        grad_T_e = grad_T_e * T_e
      end if
    end if
  end subroutine retrieve_T_e_single

  subroutine retrieve_T_e_vector(plasma_params, rhop, T_e, grad_T_e)
    use f90_kind
#ifdef IDA
    use mod_ecfm_refr_types,        only: plasma_params_type, use_ida_spline_Te, output_level, SOL_Te
#else
    use mod_ecfm_refr_types,        only: plasma_params_type, SOL_Te, output_level
#endif
    use mod_ecfm_refr_interpol,     only: spline_1d
    implicit none
    type(plasma_params_type), intent(in)    :: plasma_params
    real(rkind), dimension(:), intent(in)   :: rhop
    real(rkind), dimension(:), intent(out)  :: T_e
    real(rkind), dimension(:), intent(out), optional :: grad_T_e
    real(rkind), dimension(size(rhop))            :: rhop_aux
    if(present(grad_T_e)) grad_T_e = 0.d0
    T_e = 0.d0
    rhop_aux = rhop * plasma_params%rhop_scale_te
    where(rhop_aux > plasma_params%rhop_max) rhop_aux = plasma_params%rhop_max
    where(rhop_aux < 0.d0) rhop_aux = plasma_params%rhop_max
#ifdef IDA
    if(use_ida_spline_Te) then
      if(present(grad_T_e)) then
        call spline_1d(plasma_params%IDA_rhop_knots_Te, plasma_params%IDA_T_e, &
                       plasma_params%IDA_T_e_dx2, rhop_aux, T_e, grad_T_e)
      else
        call spline_1d(plasma_params%IDA_rhop_knots_Te, plasma_params%IDA_T_e, &
                       plasma_params%IDA_T_e_dx2, rhop_aux, T_e)
      end if
    else
#endif
      if(present(grad_T_e)) then
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%Te_spline, rhop_aux, &
                         T_e, grad_T_e, plasma_params%Te_spline_nag)
#else
          call spline_1d(plasma_params%Te_spline, rhop_aux, &
                         T_e, grad_T_e)
#endif
        else
          call spline_1d(plasma_params%Te_spline, rhop_aux, &
                         T_e, grad_T_e)
        end if
      else
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%Te_spline, rhop_aux, &
                         T_e, nag_spline = plasma_params%Te_spline_nag)
#else
          call spline_1d(plasma_params%Te_spline, rhop_aux, &
                         T_e)
#endif
        else
          call spline_1d(plasma_params%Te_spline, rhop_aux, T_e)
        end if
      end if
#ifdef IDA
    end if
#endif
    if(plasma_params%prof_log_flag) then
      T_e = exp(T_e)
      if(present(grad_T_e)) then
        grad_T_e = grad_T_e * T_e
      end if
    end if
    where(rhop * plasma_params%rhop_scale_te > plasma_params%rhop_max) T_e = SOL_Te
    where(rhop < 0.d0) T_e = SOL_Te
  end subroutine retrieve_T_e_vector

  subroutine retrieve_T_e_mat_single(plasma_params, x_vec, T_e, spatial_grad_Te)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type
    use constants,                  only: pi
    use mod_ecfm_refr_interpol,     only: rect_spline
    implicit none
    type(plasma_params_type), intent(in)     :: plasma_params
    real(rkind), dimension(:), intent(in)    :: x_vec
    real(rkind),              intent(out)   :: T_e
    real(rkind), dimension(:),intent(out), optional   :: spatial_grad_Te
    real(rkind), dimension(3)               :: R_vec
    real(rkind)                             :: dTe_dR, dTe_dz
    R_vec(1) = sqrt(x_vec(1)**2 + x_vec(2)**2)
    if(x_vec(2) /= 0.d0) then
      R_vec(2) = 2.d0 * atan((Sqrt(x_vec(1)**2 + x_vec(2)**2) - x_vec(1)) / x_vec(2))
    else if( x_vec(1) > 0.d0) then
      R_vec(2) = 0.d0
    else if(x_vec(1) < 0.d0) then
      R_vec(2) = pi
    else
      print*,"encountered atan(0/0)"
      call abort
    end if
    R_vec(3) = x_vec(3)
    if(present(spatial_grad_Te)) then
#ifdef NAG
      call rect_spline(plasma_params%T_e_spline_2D, R_vec(1), R_vec(3), T_e, dTe_dR, dTe_dz, nag_spline = plasma_params%Te_spline_nag_2D)
#else
      call rect_spline(plasma_params%T_e_spline_2D, R_vec(1), R_vec(3), T_e, dTe_dR, dTe_dz)
#endif
      T_e = exp(T_e)
      spatial_grad_Te(1) = x_vec(1) / sqrt(x_vec(1)**2 + x_vec(2)**2) *  dTe_dR * T_e
      spatial_grad_Te(2) = x_vec(2) / sqrt(x_vec(1)**2 + x_vec(2)**2) * dTe_dR * T_e
      spatial_grad_Te(3) = dTe_dz * T_e
    else
#ifdef NAG
      call rect_spline(plasma_params%T_e_spline_2D, R_vec(1), R_vec(3), T_e, nag_spline = plasma_params%Te_spline_nag_2D)
#else
      call rect_spline(plasma_params%T_e_spline_2D, R_vec(1), R_vec(3), T_e)
#endif
      T_e = exp(T_e)
    end if
  end subroutine retrieve_T_e_mat_single

  subroutine retrieve_T_e_mat_vector(plasma_params, x_vec, T_e, spatial_grad_Te)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:,:), intent(in)   :: x_vec
    real(rkind), dimension(:),  intent(out)   :: T_e
    real(rkind), dimension(:,:), intent(out), optional   :: spatial_grad_Te
    integer(ikind)  :: i
    do i = 1, size(x_vec, dim=1)
      if(present(spatial_grad_Te)) then
        call retrieve_T_e_mat_single(plasma_params, x_vec(i,:), T_e(i), spatial_grad_Te(i,:))
      else
        call retrieve_T_e_mat_single(plasma_params, x_vec(i,:), T_e(i))
      end if
    end do
  end subroutine retrieve_T_e_mat_vector

  subroutine retrieve_n_e_single(plasma_params, rhop, n_e, grad_n_e)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type, output_level, use_ida_spline_ne, SOL_ne
    use mod_ecfm_refr_interpol,     only: spline_1d
    implicit none
    type(plasma_params_type), intent(in)    :: plasma_params
    real(rkind),               intent(in)   :: rhop
    real(rkind),               intent(out)  :: n_e
    real(rkind),               intent(out), optional :: grad_n_e
    real(rkind)                              :: scaled_rhop
    scaled_rhop = rhop * plasma_params%rhop_scale_te
    if(present(grad_n_e)) grad_n_e = 0.d0
    n_e = SOL_ne
    if(scaled_rhop > plasma_params%rhop_max) return
#ifdef IDA
    if(use_ida_spline_ne) then
      if(present(grad_n_e)) then
        call spline_1d(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
                       plasma_params%IDA_n_e_dx2, scaled_rhop, n_e, grad_n_e)
      else
        call spline_1d(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
                       plasma_params%IDA_n_e_dx2, scaled_rhop, n_e)
      end if
    else
#endif
      if(present(grad_n_e)) then
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%ne_spline, scaled_rhop, &
                         n_e, grad_n_e, plasma_params%ne_spline_nag)
#else
          call spline_1d(plasma_params%ne_spline, scaled_rhop, &
                         n_e, grad_n_e)
#endif
        else
          call spline_1d(plasma_params%ne_spline, scaled_rhop, &
                         n_e, grad_n_e)
        end if
      else
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%ne_spline, scaled_rhop, &
                         n_e, nag_spline = plasma_params%ne_spline_nag)
#else
          call spline_1d(plasma_params%ne_spline, scaled_rhop, &
                         n_e)
#endif
        else
          call spline_1d(plasma_params%ne_spline, scaled_rhop, n_e)
        end if
      end if
#ifdef IDA
    end if
#endif
    if(plasma_params%prof_log_flag) then
      n_e = exp(n_e)
      if(use_ida_spline_ne) then
        n_e = n_e * 1.e6
      else
        n_e = n_e * 1.e19
      end if
      if(present(grad_n_e)) then
        grad_n_e = grad_n_e * n_e
      end if
    end if
  end subroutine retrieve_n_e_single

  subroutine retrieve_n_e_vector(plasma_params, rhop, n_e, grad_n_e)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type, output_level, use_ida_spline_ne, SOL_ne
    use mod_ecfm_refr_interpol,     only: spline_1d
    implicit none
    type(plasma_params_type), intent(in)    :: plasma_params
    real(rkind), dimension(:), intent(in)   :: rhop
    real(rkind), dimension(:), intent(out)  :: n_e
    real(rkind), dimension(:), intent(out), optional :: grad_n_e
    real(rkind), dimension(size(rhop))            :: rhop_aux
    if(present(grad_n_e)) grad_n_e = 0.d0
    n_e = SOL_ne
    rhop_aux = rhop * plasma_params%rhop_scale_ne
    where(rhop_aux > plasma_params%rhop_max) rhop_aux = plasma_params%rhop_max
    where(rhop_aux < 0.d0) rhop_aux = plasma_params%rhop_max
#ifdef IDA
    if(use_ida_spline_ne) then
      if(present(grad_n_e)) then
        call spline_1d(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
                       plasma_params%IDA_n_e_dx2, rhop_aux, n_e, grad_n_e)
      else
        call spline_1d(plasma_params%IDA_rhop_knots_ne, plasma_params%IDA_n_e, &
                       plasma_params%IDA_n_e_dx2, rhop_aux, n_e)
      end if
    else
#endif
      if(present(grad_n_e)) then
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%ne_spline, rhop_aux, &
                         n_e, grad_n_e, plasma_params%ne_spline_nag)
#else
          call spline_1d(plasma_params%ne_spline, rhop_aux, &
                         n_e, grad_n_e)
#endif
        else
          call spline_1d(plasma_params%ne_spline, rhop_aux, &
                         n_e, grad_n_e)
        end if
      else
        if(output_level) then
#ifdef NAG
          call spline_1d(plasma_params%ne_spline, rhop_aux, &
                         n_e, nag_spline = plasma_params%ne_spline_nag)
#else
          call spline_1d(plasma_params%ne_spline, rhop_aux, &
                         n_e)
#endif
        else
          call spline_1d(plasma_params%ne_spline, rhop_aux, n_e)
        end if
      end if
#ifdef IDA
    end if
#endif
    if(plasma_params%prof_log_flag) then
      n_e = exp(n_e)
      if(use_ida_spline_ne) then
        n_e = n_e * 1.e6
      else
        n_e = n_e * 1.e19
      end if
      if(present(grad_n_e)) then
        grad_n_e = grad_n_e * n_e
      end if
    end if
    where(rhop * plasma_params%rhop_scale_ne > plasma_params%rhop_max) n_e = SOL_ne
    where(rhop < 0.d0) n_e = SOL_ne
  end subroutine retrieve_n_e_vector

  subroutine retrieve_n_e_mat_single(plasma_params, x_vec, n_e, spatial_grad_ne)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type
    use mod_ecfm_refr_interpol, only: rect_spline
    use constants,                  only: pi
    implicit none
    type(plasma_params_type), intent(in)     :: plasma_params
    real(rkind), dimension(:), intent(in)    :: x_vec
    real(rkind),              intent(out)   :: n_e
    real(rkind), dimension(:),intent(out), optional   :: spatial_grad_ne
    real(rkind), dimension(3)               :: R_vec
    real(rkind)                             :: dne_dR, dne_dz
    R_vec(1) = sqrt(x_vec(1)**2 + x_vec(2)**2)
    if(x_vec(2) /= 0.d0) then
      R_vec(2) = 2.d0 * atan((Sqrt(x_vec(1)**2 + x_vec(2)**2) - x_vec(1)) / x_vec(2))
    else if( x_vec(1) > 0.d0) then
      R_vec(2) = 0.d0
    else if(x_vec(1) < 0.d0) then
      R_vec(2) = pi
    else
      print*, "encountered atan(0/0)"
      call abort
    end if
    R_vec(3) = x_vec(3)
    if(present(spatial_grad_ne)) then
#ifdef NAG
      call rect_spline(plasma_params%n_e_spline_2D, R_vec(1), R_vec(3), n_e, dne_dR, dne_dz, nag_spline = plasma_params%ne_spline_nag_2D)
#else
      call rect_spline(plasma_params%n_e_spline_2D, R_vec(1), R_vec(3), n_e, dne_dR, dne_dz)
#endif
      n_e = exp(n_e) * 1.e19
      spatial_grad_ne(1) = x_vec(1) / sqrt(x_vec(1)**2 + x_vec(2)**2) *  dne_dR * n_e
      spatial_grad_ne(2) = x_vec(2) / sqrt(x_vec(1)**2 + x_vec(2)**2) * dne_dR * n_e
      spatial_grad_ne(3) = dne_dz  * n_e
    else
#ifdef NAG
      call rect_spline(plasma_params%n_e_spline_2D, R_vec(1), R_vec(3), n_e, nag_spline = plasma_params%ne_spline_nag_2D)
#else
      call rect_spline(plasma_params%n_e_spline_2D, R_vec(1), R_vec(3), n_e)
#endif
      n_e = exp(n_e) * 1.e19
    end if
  end subroutine retrieve_n_e_mat_single

  subroutine retrieve_n_e_mat_vector(plasma_params, x_vec, n_e, spatial_grad_ne)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type
    implicit none
    type(plasma_params_type), intent(in)     :: plasma_params
    real(rkind), dimension(:,:), intent(in)   :: x_vec
    real(rkind), dimension(:),  intent(out)   :: n_e
    real(rkind), dimension(:,:), intent(out), optional   :: spatial_grad_ne
    integer(ikind)  :: i
    do i = 1, size(x_vec, dim=1)
      if(present(spatial_grad_ne)) then
        call retrieve_n_e_mat_single(plasma_params, x_vec(i,:), n_e(i), spatial_grad_ne(i,:))
      else
        call retrieve_n_e_mat_single(plasma_params, x_vec(i,:), n_e(i))
      end if
    end do
  end subroutine retrieve_n_e_mat_vector

subroutine bin_ray_BPD_to_common_rhop(plasma_params, rad_mode, center_freq, weights, rhop, BPD, BPD_secondary)
use f90_kind
use mod_ecfm_refr_types,       only: plasma_params_type, rad_diag_ch_mode_type, N_ray, spl_type_1d, max_points_svec, max_rhop_BPD
use mod_ecfm_refr_interpol,   only: make_1d_spline,  spline_1d, spline_1d, spline_1d_get_roots, deallocate_1d_spline, spline_1d_integrate
use constants,                 only: pi, mass_e, e0, c0
implicit none
type(plasma_params_type), intent(in)  :: plasma_params
type(rad_diag_ch_mode_type), intent(inout) :: rad_mode
real(rkind), dimension(:),  intent(in)     :: weights ! ray weights
real(rkind), intent(in)                    :: center_freq
real(rkind), dimension(:), intent(inout)   :: rhop, BPD
real(rkind), dimension(:), intent(inout), optional   :: BPD_secondary
integer(ikind)                             :: ir, i, i_start, i_end, root_cnt, i_am_root, &
                                              i_seg_end, i_seg_start, useful_N
real(rkind)                                :: max_s_seg, min_s_seg, rhop_min_seg, rhop_max_seg, &
                                              BPD_norm, BPD_secondary_norm
type(spl_type_1d)                          :: rhop_spl, BPD_spl, BPD_secondary_spl
real(rkind), dimension(100)                :: roots, BPD_step, BPD_secondary_step
real(rkind), dimension(max_points_svec) :: s_arr, rhop_arr, BPD_arr, BPD_second_arr
real(rkind), dimension(10) :: s_aux_arr, rhop_aux_arr! used to increase points to suffice for root finding
logical                                    :: make_secondary_BPD
  make_secondary_BPD = .false.
  if(present(BPD_secondary)) make_secondary_BPD = .true.
  do i = 1, size(rhop)
    rhop(i) = 2.d0 * real(i - 1, 8) * max_rhop_BPD / &
              real(size(rhop) - 1, 8) - max_rhop_BPD
  end do
  BPD(:) = 0.d0
  if(make_secondary_BPD) BPD_secondary(:) = 0.d0
  BPD_step(:) = 0.d0
  if(make_secondary_BPD) BPD_secondary_step(:) = 0.d0
  do ir = 1, N_ray
    if(rad_mode%ray(ir)%Trad == 0.d0) then
      if(.not. make_secondary_BPD) then
        !print*, "Skipped creating BPD for a channel with zero Trad"
        cycle
      else if(rad_mode%ray(ir)%Trad_secondary == 0.d0) then
        !print*, "Skipped creating BPD for a channel with zero Trad"
        cycle
      end if
    end if
    i_start = 1
    i_end = rad_mode%ray(ir)%freq(1)%total_LOS_points
    do while((rad_mode%ray(ir)%freq(1)%svec(i_start)%rhop > plasma_params%rhop_max .or. &
              rad_mode%ray(ir)%freq(1)%svec(i_start)%rhop == -1.d0) .and. &
              i_start < i_end)
      i_start = i_start + 1
    end do
    do while((rad_mode%ray(ir)%freq(1)%svec(i_end)%rhop > plasma_params%rhop_max .or. &
              rad_mode%ray(ir)%freq(1)%svec(i_end)%rhop == -1.d0) .and. &
              i_start < i_end)
      i_end = i_end - 1
    end do
    if(i_end > max_points_svec) then
      print*, "length of svec",  max_points_svec
      print*, "i_end>total_LOS_points", i_end
      call abort()
    end if
    if(i_start >= i_end) then
      print*, "Critical error when calculating BPD"
      print*, "Could not find a single valid rhop along LOS"
      stop "mod_ecfm_refr_utils.f90 - invalid LOS or equilibrium"
    end if
    useful_N = i_end - i_start + 1 ! amount of points inside the flux matrix and profiles
    s_arr(1:useful_N) = rad_mode%ray(ir)%freq(1)%svec(i_start:i_end)%s
    if(any(s_arr(1:useful_N) /= s_arr(1:useful_N))) then
      print*, "Something very wrong with central LOS for making mode BPD"
      call abort()
    end if
    BPD_arr(1:useful_N) = rad_mode%ray_extra_output(ir)%em(i_start:i_end) * &
                                             c0**2 / (center_freq**2 * e0) * &
                                             rad_mode%ray_extra_output(ir)%T(i_start:i_end)
    if(make_secondary_BPD) then
      BPD_second_arr(1:useful_N) = rad_mode%ray_extra_output(ir)%em_secondary(i_start:i_end) * &
                                                c0**2 / (center_freq**2 * e0) * &
                                                rad_mode%ray_extra_output(ir)%T_secondary(i_start:i_end)
    end if
    rhop_arr(1:useful_N) = rad_mode%ray(ir)%freq(1)%svec(i_start:i_end)%rhop
    call make_1d_spline(BPD_spl, int(useful_N, 4), &
                          s_arr(1:useful_N), &
                          BPD_arr(1:useful_N), k=1, iopt=0) ! We want linear splines here to avoid ringing
    if(make_secondary_BPD) then
      call make_1d_spline(BPD_secondary_spl, int(useful_N, 4), &
                            s_arr(1:useful_N), &
                            BPD_second_arr(1:useful_N), k = 1, iopt=0) ! We want linear splines here to avoid ringing
    end if
    where(rad_mode%ray(ir)%freq(1)%svec(i_start:i_end)%freq_2X *  pi * mass_e / e0 > abs(plasma_params%B_ax)) &
        rhop_arr(1:useful_N) = -rhop_arr(1:useful_N)
    ! To get avoid interpolation errors we have to avoid the discontuity from positive to negative rho_pol
    ! Hence HFS and LFS side have to be treated separately
    i_seg_start = i_start
    i_seg_end = i_seg_start + 1
    do while(i_seg_end < i_end)
      do while((rad_mode%ray(ir)%freq(1)%svec(i_seg_start)%freq_2X *  pi * mass_e / e0 -  abs(plasma_params%B_ax)) * &
               (rad_mode%ray(ir)%freq(1)%svec(i_seg_end)%freq_2X *  pi * mass_e / e0 -  abs(plasma_params%B_ax)) > 0.d0 .and. &
               i_seg_end < i_end)
        i_seg_end = i_seg_end + 1
      end do
      i_seg_end = i_seg_end - 1
      rhop_min_seg = minval(rhop_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))
      rhop_max_seg = maxval(rhop_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))
      !print*, "i_seg/i_end", i_seg_end, "/", i_end
      !print*, "segment length", i_seg_end - i_seg_start
      if(sum(exp(BPD_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))) / rad_mode%ray(ir)%Trad < 1.d-8) then
         if(.not. make_secondary_BPD) then
           i_seg_start = i_seg_end + 1
           i_seg_end = i_seg_start + 1
           cycle
         else if(sum(exp(BPD_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1)) + &
             exp(BPD_second_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))) &
              / max(rad_mode%ray(ir)%Trad, rad_mode%ray(ir)%Trad_secondary) < 1.d-8) then
        !print*, "section skipped - BPD sum", sum(exp(BPD_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1)) + &
        !                                         exp(BPD_second_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1)))
           i_seg_start = i_seg_end + 1
           i_seg_end = i_seg_start + 1
           cycle
         end if
      end if
      if(i_seg_start + 1 >= i_seg_end) then
        if(.not. make_secondary_BPD .and. sum(exp(BPD_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))) / rad_mode%ray(ir)%Trad < 1.d-6) then
           print*, "Two points on the HFS were skipped for the birthplace distribution"
           print*, "It might be necessary to rerun the forward model with higher resolution in s if BDP is important"
        else if(sum(exp(BPD_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1)) + &
                exp(BPD_second_arr(i_seg_start - i_start + 1:i_seg_end - i_start + 1))) &
                / max(rad_mode%ray(ir)%Trad, rad_mode%ray(ir)%Trad_secondary) < 1.d-6) then
           print*, "Two points on the HFS were skipped for the birthplace distribution"
           print*, "It might be necessary to rerun the forward model with higher resolution in s if BDP is important"
           i_seg_start = i_seg_end + 1
           i_seg_end = i_seg_start + 1
        end if
        i_seg_start = i_seg_end + 1
        i_seg_end = i_seg_start + 1
        cycle
      end if
  !    print*, "B_ax", plasma_params%B_ax
  !    print*, "B_tot", rad_mode%ray(ir)%freq(1)%svec(i_start:i_end)%freq_2X *  pi * mass_e / e0
      ! Continouus arrays for rhop -> speed up in loop over rhop
      !print*, rhop_arr(1:i_end - i_start)
      ! Continouus array for BPD -> suppresses runtime warnings
      !print*, "s_arr", s_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1)
      if(i_seg_end - i_seg_start + 1 >= 8) then
        rhop_spl%iopt_int = 0
        do i = 1, size(rhop)
          if(.not. (rhop_min_seg < rhop(i) .and. rhop(i) < rhop_max_seg)) then
            !print*, "seg boundaries", rhop(i), rhop_min_seg, rhop_max_seg
            cycle
          end if
          call make_1d_spline(rhop_spl, int(i_seg_end - i_seg_start + 1, 4), &
                              s_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1), &
                              rhop_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1) - rhop(i), rhop_spl%iopt_int)
          rhop_spl%iopt_int = 1
          !print*, "Now making roots with", i_start, 1 + i_seg_start - i_start, i_seg_end - i_start + 1, i_end
          call spline_1d_get_roots(rhop_spl, roots, root_cnt)
          !print*, "rhop, amount of crossings", rhop(i), root_cnt
          if(root_cnt > 0) then
            call spline_1d(BPD_spl, roots(1:root_cnt), BPD_step(1:root_cnt))
            if(make_secondary_BPD) call spline_1d(BPD_secondary_spl, roots(1:root_cnt), BPD_secondary_step(1:root_cnt))
            do i_am_root = 1, root_cnt
              if(rad_mode%ray(ir)%Trad > 0.d0) BPD(i) = BPD(i) + weights(ir) * BPD_step(i_am_root) / rad_mode%ray(ir)%Trad
              if(make_secondary_BPD) then
                if(rad_mode%ray(ir)%Trad_secondary > 0.d0) BPD_secondary(i) = &
                  BPD_secondary(i) + weights(ir) * BPD_secondary_step(i_am_root) / rad_mode%ray(ir)%Trad_secondary
              end if
            end do
          end if
        end do
      else
      ! Not enough points for the roots method -> create auxiliary array with sufficient points first and find the roots in this one
        min_s_seg = minval(s_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1))
        max_s_seg =  maxval(s_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1))
        if(max_s_seg - min_s_seg == 0.d0) then
          i_seg_start = i_seg_end + 1
          i_seg_end = i_seg_start + 1
          cycle
        end if
        do i = 1, size(s_aux_arr)
          s_aux_arr(i) = min_s_seg + real(i - 1,8) / real(size(s_aux_arr) -1,8) * (max_s_seg - min_s_seg)
        end do
        !print*, "s_aux_arr", s_aux_arr
        call make_1d_spline(rhop_spl, int(i_seg_end - i_seg_start + 1, 4), &
                          s_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1), &
                          rhop_arr(1 + i_seg_start - i_start :i_seg_end - i_start + 1))
        call spline_1d(rhop_spl, s_aux_arr, rhop_aux_arr)
        call deallocate_1d_spline(rhop_spl)
        rhop_spl%iopt_int = 0
        if(any(rhop_aux_arr /= rhop_aux_arr)) then
          print*, "Nan in rhop for BPD binning"
          call abort()
        end if
        do i = 1, size(rhop)
          if(.not. (rhop_min_seg < rhop(i) .and. rhop(i) < rhop_max_seg)) then
            !print*, "seg boundaries", rhop(i), rhop_min_seg, rhop_max_seg
            cycle
          end if
          call make_1d_spline(rhop_spl, int(size(s_aux_arr), 4), &
                              s_aux_arr, rhop_aux_arr - rhop(i), rhop_spl%iopt_int)
          rhop_spl%iopt_int = 1
          call spline_1d_get_roots(rhop_spl, roots, root_cnt)
          !print*, "rhop, amount of crossings", rhop(i), root_cnt
          if(root_cnt > 0) then
            call spline_1d(BPD_spl, roots(1:root_cnt), BPD_step(1:root_cnt))
            if(make_secondary_BPD)  call spline_1d(BPD_secondary_spl, roots(1:root_cnt), BPD_secondary_step(1:root_cnt))
            do i_am_root = 1, root_cnt
              if(rad_mode%ray(ir)%Trad > 0.d0) BPD(i) = BPD(i) + weights(ir) * BPD_step(i_am_root) / rad_mode%ray(ir)%Trad
              if(make_secondary_BPD) then
                if(rad_mode%ray(ir)%Trad_secondary > 0.d0) BPD_secondary(i) = &
                  BPD_secondary(i) + weights(ir) * BPD_secondary_step(i_am_root) / rad_mode%ray(ir)%Trad_secondary
              end if
            end do
          end if
        end do
      end if
      call deallocate_1d_spline(rhop_spl)
      i_seg_start = i_seg_end + 1
      i_seg_end = i_seg_start + 1
    end do ! rhops
  !print*, "Ray ", ir, " finished"
  call deallocate_1d_spline(BPD_spl)
  if(make_secondary_BPD) call deallocate_1d_spline(BPD_secondary_spl)
end do ! ir
call deallocate_1d_spline(BPD_spl)
if(make_secondary_BPD) call deallocate_1d_spline(BPD_secondary_spl)
! Normalization for rhop coordinate system
call make_1d_spline(BPD_spl, int(size(rhop), 4), &
                    rhop, BPD)
if(make_secondary_BPD) call make_1d_spline(BPD_secondary_spl, int(size(rhop), 4), &
                    rhop, BPD_secondary)
call spline_1d_integrate(BPD_spl, -1.d0, 1.d0, BPD_norm)
if(make_secondary_BPD) call spline_1d_integrate(BPD_secondary_spl, -1.d0, 1.d0, BPD_secondary_norm)
if(BPD_norm > 0) BPD = BPD / BPD_norm
if(make_secondary_BPD)  then
  if(BPD_secondary_norm > 0) BPD_secondary = BPD_secondary / BPD_secondary_norm
end if
call deallocate_1d_spline(BPD_spl)
if(make_secondary_BPD)  call deallocate_1d_spline(BPD_secondary_spl)
end subroutine bin_ray_BPD_to_common_rhop

subroutine bin_freq_to_ray(ray_extra_output, freq_weight, rad_freq, total_LOS_points)
! This routine averages the radiation transport results from the individual sub-frequencies
! of the bandwidth discretication and stores them in the ray_extra_output structure.
! The routine assumes that the R(s) and z(s) are independent of the subfrequency from the bandwidth
! discretication. The step size of the individual svec grids can be arbitrary.
use f90_kind
use mod_ecfm_refr_types,       only: rad_diag_ch_mode_ray_extra_output_type, &
                                     rad_diag_ch_mode_ray_freq_type, &
                                     N_freq, spl_type_1d, &
                                     output_level
use mod_ecfm_refr_interpol,   only: make_1d_spline, spline_1d, deallocate_1d_spline, spline_1d_integrate
use constants,                only: c0, e0
implicit none
type(rad_diag_ch_mode_ray_extra_output_type), intent(inout)  :: ray_extra_output
real(rkind), dimension(:), intent(in)      :: freq_weight
type(rad_diag_ch_mode_ray_freq_type), dimension(:),  intent(in) :: rad_freq
integer(ikind), dimension(:), intent(in)   :: total_LOS_points
integer(ikind)                             :: ifreq
type(spl_type_1d)                          :: spl
real(rkind), dimension(:), allocatable     :: s_arr, val, s_freq, quant
real(rkind)                                :: BPD_norm, BPD_secondary_norm
  allocate(s_arr(total_LOS_points(1)), quant(maxval(total_LOS_points(:))))
  s_arr(:) = rad_freq(1)%svec(1:total_LOS_points(1))%s
  if(N_freq == 1) then
    ray_extra_output%Trad(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%Trad
    ray_extra_output%em(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%em
    ray_extra_output%T(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%T
    ray_extra_output%ab(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%ab
    if(output_level) then
      ray_extra_output%Trad_secondary(1:total_LOS_points(1))= rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%Trad_secondary
      ray_extra_output%em_secondary(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%em_secondary
      ray_extra_output%T_secondary(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%T_secondary
      ray_extra_output%ab_secondary(1:total_LOS_points(1)) = rad_freq(1)%svec_extra_output(1:total_LOS_points(1))%ab_secondary
    end if
  else
    ray_extra_output%Trad(:) = 0.d0
    ray_extra_output%em(:) = 0.d0
    ray_extra_output%T(:) = 0.d0
    ray_extra_output%ab(:) = 0.d0
    if(output_level) then
      ray_extra_output%Trad_secondary(:) = 0.d0
      ray_extra_output%em_secondary(:) = 0.d0
      ray_extra_output%T_secondary(:) = 0.d0
      ray_extra_output%ab_secondary(:) = 0.d0
    end if
    allocate(val(total_LOS_points(1)), s_freq(maxval(total_LOS_points(:))))
    !print*, "s_arr", s_arr(1:total_LOS_points(1))
    do ifreq = 2, N_freq
      if(any(rad_freq(ifreq)%svec(1:total_LOS_points(ifreq))%s /= &
             rad_freq(ifreq)%svec(1:total_LOS_points(ifreq))%s)) then
        print*, "Something very wrong with freq LOS for making ray BPD"
        call abort()
      end if
      s_freq(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec(1:total_LOS_points(ifreq))%s
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%Trad
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%Trad(1:total_LOS_points(1)) = ray_extra_output%Trad(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%em
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%em(1:total_LOS_points(1)) = ray_extra_output%em(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%T
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%T(1:total_LOS_points(1)) = ray_extra_output%T(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%ab
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%ab(1:total_LOS_points(1)) = ray_extra_output%ab(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      if(.not. output_level) cycle
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%Trad_secondary
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%Trad_secondary(1:total_LOS_points(1)) = ray_extra_output%Trad_secondary(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%em_secondary
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%em_secondary(1:total_LOS_points(1)) = ray_extra_output%em_secondary(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%T_secondary
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%T_secondary(1:total_LOS_points(1)) = ray_extra_output%T_secondary(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
      quant(1:total_LOS_points(ifreq)) = rad_freq(ifreq)%svec_extra_output%ab_secondary
      call make_1d_spline(spl, int(total_LOS_points(ifreq), 4), &
                            s_freq(1:total_LOS_points(ifreq)), &
                            quant(1:total_LOS_points(ifreq)), k=1)
      call spline_1d(spl, s_arr, val)
      ray_extra_output%ab_secondary(1:total_LOS_points(1)) = ray_extra_output%ab_secondary(1:total_LOS_points(1)) + freq_weight(ifreq) * val
      call deallocate_1d_spline(spl)
    end do ! ifreq =1, N_freq
    deallocate(val, s_freq)
  end if ! N_freq == 1
  ray_extra_output%BPD(1:total_LOS_points(1)) = ray_extra_output%T(1:total_LOS_points(1)) * ray_extra_output%em(1:total_LOS_points(1)) * &
                                                c0**2 / (100.d9**2 * e0) ! brings integral closer to 1
  if(output_level) then
    ray_extra_output%BPD_secondary(1:total_LOS_points(1)) = ray_extra_output%T_secondary(1:total_LOS_points(1)) * ray_extra_output%em_secondary(1:total_LOS_points(1)) * &
                                                            c0**2 / (100.d9**2 * e0) ! brings integral closer to 1
  end if
  quant(1:total_LOS_points(1)) = ray_extra_output%BPD(1:total_LOS_points(1))
  call make_1d_spline(spl, int(total_LOS_points(1), 4), &
                            s_arr(1:total_LOS_points(1)), &
                            quant(1:total_LOS_points(1)))
  call spline_1d_integrate(spl, s_arr(1), s_arr(total_LOS_points(1)), BPD_norm)
  call deallocate_1d_spline(spl)
  if(BPD_norm > 0.d0) ray_extra_output%BPD(1:total_LOS_points(1)) = ray_extra_output%BPD(1:total_LOS_points(1)) / BPD_norm
  if(output_level) then
    quant(1:total_LOS_points(1)) = ray_extra_output%BPD_secondary(1:total_LOS_points(1))
    call make_1d_spline(spl, int(total_LOS_points(1), 4), &
                            s_arr(1:total_LOS_points(1)), &
                            quant(1:total_LOS_points(1)))
    call spline_1d_integrate(spl, s_arr(1), s_arr(total_LOS_points(1)), BPD_secondary_norm)
    if(BPD_secondary_norm > 0.d0) ray_extra_output%BPD_secondary(1:total_LOS_points(1)) = ray_extra_output%BPD_secondary(1:total_LOS_points(1)) / BPD_secondary_norm
  end if
  deallocate(s_arr, quant)
end subroutine bin_freq_to_ray

subroutine make_warm_res_mode(rad_mode, weights, f_ECE, make_secondary)
use f90_kind
use mod_ecfm_refr_types,       only: rad_diag_ch_mode_type, N_ray, use_maximum_for_warm_res
use constants,                 only: c0, e0
implicit none
type(rad_diag_ch_mode_type), intent(inout) :: rad_mode
real(rkind), dimension(:), intent(in)      :: weights
logical, intent(in)                        :: make_secondary
real(rkind), intent(in)                    :: f_ECE ! ray weights, freq_weights
integer(ikind)                             :: ir, i, i_max_BDOP, i_max_BDOP_secondary
real(rkind)                                :: sTrad, Trad, sTrad_secondary, Trad_secondary, &
                                              BPD_cur, BPD_cur_secondary, ds
  sTrad = 0.d0
  Trad = 0.d0
  i_max_BDOP = 1
  rad_mode%rel_s_res = 0.d0
  rad_mode%rel_R_res = 0.d0
  rad_mode%rel_z_res = 0.d0
  rad_mode%rel_rhop_res = 0.d0
  if(make_secondary) then
    sTrad_secondary = 0.d0
    Trad_secondary= 0.d0
    i_max_BDOP_secondary = 1
    rad_mode%rel_s_res_secondary = 0.d0
    rad_mode%rel_R_res_secondary = 0.d0
    rad_mode%rel_z_res_secondary = 0.d0
    rad_mode%rel_rhop_res_secondary = 0.d0
  end if
  do ir = 1, N_ray
    ! two steps - first find warm res for each ray and then sum the rays
    rad_mode%ray(ir)%rel_s_res = 0.d0
    rad_mode%ray(ir)%rel_R_res = 0.d0
    rad_mode%ray(ir)%rel_z_res = 0.d0
    rad_mode%ray(ir)%rel_rhop_res = 0.d0
    if(make_secondary) then
      rad_mode%ray(ir)%rel_s_res_secondary = 0.d0
      rad_mode%ray(ir)%rel_R_res_secondary = 0.d0
      rad_mode%ray(ir)%rel_z_res_secondary = 0.d0
      rad_mode%ray(ir)%rel_rhop_res_secondary = 0.d0
    end if
    ! Two options here, either point of maximum emission, or mean of the distribution
    ! the mean has the disadvantage to be horribly wrong for a non-unimodal BPD
    ! the maximum can give one a false sence of  for a non-unimodal BPD

    if(use_maximum_for_warm_res) then
      i_max_BDOP = maxloc(rad_mode%ray_extra_output(ir)%em(1:rad_mode%ray(ir)%freq(1)%total_LOS_points) * &
                          rad_mode%ray_extra_output(ir)%T(1:rad_mode%ray(ir)%freq(1)%total_LOS_points), dim=1)
      if(make_secondary)  i_max_BDOP_secondary = maxloc(rad_mode%ray_extra_output(ir)%em_secondary(1:rad_mode%ray(ir)%freq(1)%total_LOS_points) * &
                             rad_mode%ray_extra_output(ir)%T_secondary(1:rad_mode%ray(ir)%freq(1)%total_LOS_points), dim=1)
    else
      do i = 1, rad_mode%ray(ir)%freq(1)%total_LOS_points
        BPD_cur = rad_mode%ray_extra_output(ir)%em(i) * &
              rad_mode%ray_extra_output(ir)%T(i) /  &
              rad_mode%ray(ir)%Trad * c0**2 / (f_ECE**2 * e0)
        if(make_secondary)  BPD_cur_secondary = rad_mode%ray_extra_output(ir)%em_secondary(i) * &
                        rad_mode%ray_extra_output(ir)%T_secondary(i) /  &
                        rad_mode%ray(ir)%Trad_secondary * c0**2 / (f_ECE * e0)
      ! Simple Trapezoid integration
        if( i == 1) then
          ds = (rad_mode%ray(ir)%freq(1)%svec(i + 1)%s  - &
                rad_mode%ray(ir)%freq(1)%svec(i)%s) * 0.5d0
          sTrad = sTrad + BPD_cur *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          if(make_secondary)  sTrad_secondary = sTrad_secondary + BPD_cur_secondary *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          Trad = Trad + BPD_cur * ds
          if(make_secondary)  Trad_secondary = Trad_secondary + BPD_cur_secondary * ds
            ! normalize with the Trad from the integrated Trad, since Trad from the differential equation might be slightly different (numerical errors)
        else if(i == rad_mode%ray(ir)%freq(1)%total_LOS_points) then
          ds = (rad_mode%ray(ir)%freq(1)%svec(i)%s  - &
                rad_mode%ray(ir)%freq(1)%svec(i - 1)%s) * 0.5d0
          sTrad = sTrad + BPD_cur *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          if(make_secondary)  sTrad_secondary = sTrad_secondary + BPD_cur_secondary *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          Trad = Trad + BPD_cur * ds
          if(make_secondary)  Trad_secondary = Trad_secondary + BPD_cur_secondary * ds
        else
          ds = rad_mode%ray(ir)%freq(1)%svec(i + 1)%s  - &
                rad_mode%ray(ir)%freq(1)%svec(i)%s
          sTrad = sTrad + BPD_cur *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          if(make_secondary)  sTrad_secondary = sTrad_secondary + BPD_cur_secondary *  &
                  rad_mode%ray(ir)%freq(1)%svec(i)%s * ds
          Trad = Trad + BPD_cur * ds
          if(make_secondary) Trad_secondary = Trad_secondary + BPD_cur_secondary * ds
        end if
      end do !i
      i_max_BDOP  = minloc(abs(sTrad / Trad - rad_mode%ray(ir)%freq(1)%svec(1:rad_mode%ray(ir)%freq(1)%total_LOS_points)%s), dim = 1)
      if(make_secondary)  i_max_BDOP_secondary  = minloc(abs(sTrad_secondary / Trad_secondary - rad_mode%ray(ir)%freq(1)%svec(1:rad_mode%ray(ir)%freq(1)%total_LOS_points)%s), dim = 1)
    end if
    !print*, "rhop warm res ir", rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP)%rhop, ir
    rad_mode%ray(ir)%rel_s_res = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP)%s
    rad_mode%ray(ir)%rel_R_res = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP)%R
    rad_mode%ray(ir)%rel_z_res = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP)%z
    rad_mode%ray(ir)%rel_rhop_res = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP)%rhop
    rad_mode%rel_s_res = rad_mode%rel_s_res + rad_mode%ray(ir)%rel_s_res * weights(ir)
    rad_mode%rel_R_res = rad_mode%rel_R_res + rad_mode%ray(ir)%rel_R_res * weights(ir)
    rad_mode%rel_z_res = rad_mode%rel_z_res + rad_mode%ray(ir)%rel_z_res * weights(ir)
    rad_mode%rel_rhop_res = rad_mode%rel_rhop_res + rad_mode%ray(ir)%rel_rhop_res * weights(ir)
    if(make_secondary)  then
      rad_mode%ray(ir)%rel_s_res_secondary = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP_secondary)%s
      rad_mode%ray(ir)%rel_R_res_secondary = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP_secondary)%R
      rad_mode%ray(ir)%rel_z_res_secondary = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP_secondary)%z
      rad_mode%ray(ir)%rel_rhop_res_secondary = rad_mode%ray(ir)%freq(1)%svec(i_max_BDOP_secondary)%rhop
      rad_mode%rel_s_res_secondary = rad_mode%rel_s_res_secondary + rad_mode%ray(ir)%rel_s_res_secondary * weights(ir)
      rad_mode%rel_R_res_secondary = rad_mode%rel_R_res_secondary + rad_mode%ray(ir)%rel_R_res_secondary * weights(ir)
      rad_mode%rel_z_res_secondary = rad_mode%rel_z_res_secondary + rad_mode%ray(ir)%rel_z_res_secondary * weights(ir)
      rad_mode%rel_rhop_res_secondary = rad_mode%rel_rhop_res_secondary + rad_mode%ray(ir)%rel_rhop_res_secondary * weights(ir)
    end if
  end do !ir
end subroutine make_warm_res_mode

!*****************************************************************


subroutine init_non_therm()
use mod_ecfm_refr_types,               only: dstf
implicit None
  if(trim(dstf) == "numeric") call read_ffp_data()
  if(trim(dstf) == "gene") call read_fgene_data()
  if(trim(dstf) == "Bi_Maxw" .or. dstf == "Bi_MaxJ") call read_bi_max_data()
  if(trim(dstf) == "runaway") call read_runaway_data()
  if(trim(dstf) == "drift_m") call read_drift_m_data()
  if(trim(dstf) == "Spitzer") call read_Spitzer_data()
  if(trim(dstf) == "multi_s") call read_multi_slope_data()
  if(trim(dstf) == "gcomp") then
    call read_fgene_data()
    call read_Gene_Bi_max_data()
  end if
end subroutine init_non_therm

subroutine read_bi_max_data()
use mod_ecfm_refr_types, only: bi_max, data_folder, min_level_log_ne
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop, ne_vec
Character(200)  :: cur_filename
  cur_filename =  trim(data_folder) // "bi_max.dat"
  open(67, file = trim(cur_filename))
  read(67,"(2E15.8E2)") bi_max%Te_par, bi_max%Te_perp
  read(67,"(I3.3)") N_rhop
  allocate(rhop(N_rhop), ne_vec(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(2E15.8E2)") rhop(irhop), ne_vec(irhop)
  end do
  close(67)
  bi_max%rhop_max = maxval(rhop)
  where(ne_vec < min_level_log_ne) ne_vec = min_level_log_ne
  ! Interpolate the logarithm for positivity
  call make_1d_spline(bi_max%ne_spl, int(N_rhop,4), rhop, log(ne_vec))
#ifdef NAG
  call nag_spline_1d_interp(rhop, log(ne_vec), bi_max%ne_nag_spl)
#endif
  deallocate(rhop, ne_vec)
end subroutine read_bi_max_data

subroutine read_runaway_data()
use mod_ecfm_refr_types, only: runaway, data_folder, min_level_log_ne
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop, rel_ne
Character(200)  :: cur_filename
  cur_filename = trim(data_folder) // "runaway.dat"
  open(67, file = trim(cur_filename))
  read(67,"(2E15.8E2)")  runaway%Z_eff,runaway%E_E_c
  read(67,"(I3.3)")  N_rhop
  allocate(rhop(N_rhop), rel_ne(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(2E15.8E2)") rhop(irhop), rel_ne(irhop)
  end do
  close(67)
  runaway%rhop_max = maxval(rhop)
  where(rel_ne < min_level_log_ne) rel_ne = min_level_log_ne
  ! Interpolate the logarithm for positivity
  call make_1d_spline(runaway%rel_ne_spl, int(N_rhop,4), rhop, log(rel_ne))
#ifdef NAG
  call nag_spline_1d_interp(rhop, log(rel_ne), runaway%rel_ne_nag_spl)
#endif
  deallocate(rhop, rel_ne)
end subroutine read_runaway_data

subroutine read_drift_m_data()
use mod_ecfm_refr_types, only: drift_m, data_folder, min_level_log_ne
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop, ne_vec
Character(200)  :: cur_filename
  cur_filename = trim(data_folder) // "drift_m.dat"
  open(67, file = trim(cur_filename))
  read(67,"(4E15.8E2)") drift_m%Te_par,drift_m%Te_perp, &
                        drift_m%u_par_drift, drift_m%u_perp_drift
  read(67,"(I3.3)")  N_rhop
  allocate(rhop(N_rhop), ne_vec(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(2E15.8E2)") rhop(irhop), ne_vec(irhop)
  end do
  close(67)
  where(ne_vec < min_level_log_ne) ne_vec = min_level_log_ne
  ! Interpolate the logarithm for positivity
  call make_1d_spline(drift_m%ne_spl, int(N_rhop,4), rhop, log(ne_vec))
#ifdef NAG
  call nag_spline_1d_interp(rhop, log(ne_vec), drift_m%ne_nag_spl)
#endif
  if(drift_m%u_perp_drift /= 0.0 ) stop "No perpendicular drifts allowed"
  drift_m%rhop_max = maxval(rhop)
  deallocate(rhop, ne_vec)
end subroutine read_drift_m_data

subroutine read_multi_slope_data()
use mod_ecfm_refr_types, only: multi_slope, data_folder
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop, Te_slope_vec
Character(200)  :: cur_filename
  cur_filename = trim(data_folder) // "multi_s.dat"
  open(67, file = trim(cur_filename))
  read(67,"(E15.8E2)")  multi_slope%gamma_switch
  read(67,"(I3.3)")  N_rhop
  allocate(rhop(N_rhop), Te_slope_vec(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(2E15.8E2)") rhop(irhop), Te_slope_vec(irhop)
  end do
  close(67)
  call make_1d_spline(multi_slope%Te_slope_spl, int(N_rhop,4), rhop, Te_slope_vec)
#ifdef NAG
  call nag_spline_1d_interp(rhop, Te_slope_vec, multi_slope%Te_slope_nag_spl)
#endif
  multi_slope%rhop_max = maxval(rhop)
  deallocate(rhop, Te_slope_vec)
end subroutine read_multi_slope_data

subroutine read_ffp_data()
  use mod_ecfm_refr_types, only: ffp, data_folder, dst_data_folder
  use mod_ecfm_refr_fp_dist_utils, only: setup_f_rhop_splines, make_rhop_Bmin
  use constants, only: pi
  implicit none
  integer(ikind)  :: irhop, iu, ipitch
  Character(200)  :: cur_filename,filename
  Character(12)  :: irhop_str
  Character(16)  :: format
  filename = trim(data_folder) // dst_data_folder // "/u.dat"
  open(66, file = trim(filename))
  read(66,"(I5.5)") ffp%N_u
  allocate(ffp%u(ffp%N_u))
  do iu = 1, ffp%N_u
    read(66,"(E19.12E2)") ffp%u(iu)
  end do
  close(66)
  ffp%u_max_data =  ffp%u(ffp%N_u)
  ffp%u_max = ffp%u_max_data
  ffp%u_min = max(ffp%u(1), ffp%u_min)
  filename = trim(data_folder) // dst_data_folder // "/pitch.dat"
  open(66, file = trim(filename))
  read(66,"(I5.5)") ffp%N_pitch
  allocate(ffp%pitch(ffp%N_pitch))
  do ipitch = 1, ffp%N_pitch
    read(66,"(E19.12E2)") ffp%pitch(ipitch)
  end do
  close(66)
  ffp%pitch_min =  ffp%pitch(1)
  ffp%pitch_max =  ffp%pitch(ffp%N_pitch)
  filename = trim(data_folder) // dst_data_folder // "/frhop.dat"
  open(67, file = trim(filename))
  read(67,"(I5.5)") ffp%N_rhop
  allocate(ffp%rhop(ffp%N_rhop))
  allocate(ffp%f(ffp%N_rhop, ffp%N_u, ffp%N_pitch)) !
  print*,"Loading drift kinetic distribution function"
  write(format,"(A1I3.3A12)") "(",ffp%N_pitch,"(E15.8E2A1))"
  do irhop = 1, ffp%N_rhop
    read(67,"(E19.12E2)") ffp%rhop(irhop)
    write(irhop_str, "(I3.3)") irhop - 1
    cur_filename = trim(data_folder) // dst_data_folder // "/fu" // trim(irhop_str) // ".dat"
    open(66, file= trim(cur_filename))
    do iu = 1, ffp%N_u
      read(66, *) (ffp%f(irhop, iu, ipitch), ipitch = 1, ffp%N_pitch) ! ipitch = 7, ffp%N_pitch - 6
    end do
    close(66)
  end do
  close(67)
  ffp%rhop_min = ffp%rhop(1)
  ffp%rhop_max = ffp%rhop(ffp%N_rhop) - 1.d-2 ! set boundary in a way that we can fall back to thermal
  call setup_f_rhop_splines(ffp)
  call make_rhop_Bmin()
end subroutine read_ffp_data

subroutine read_fgene_data()
  use mod_ecfm_refr_types, only: fgene, data_folder, dst_data_folder, N_absz_large
  use mod_ecfm_refr_interpol, only: make_rect_spline
  use mod_ecfm_refr_gene_dist_utils, only: setup_fgene_rhop_splines
#ifdef NAG
  USE nag_spline_2d,                only: nag_spline_2d_interp
#endif
  use constants, only: pi
  implicit none
  integer(ikind)  :: irhop, ivpar, imu
  real(rkind), dimension(:,:), allocatable :: f_0
  Character(200)  :: cur_filename,filename
  Character(12)  :: irhop_str
  Character(16)  :: format
  filename = trim(data_folder) // dst_data_folder // "/vpar.dat"
  open(66, file = trim(filename))
  read(66,"(I5.5)") fgene%N_vpar
  fgene%N_vpar = fgene%N_vpar
  allocate(fgene%vpar(fgene%N_vpar))
  do ivpar = 1, fgene%N_vpar
    read(66,"(E19.12E2)") fgene%vpar(ivpar)
  end do
  close(66)
  fgene%vpar_max =  maxval(fgene%vpar)
  fgene%vpar_min = minval(fgene%vpar(:))
  filename = trim(data_folder) // dst_data_folder // "/mu.dat"
  open(66, file = trim(filename))
  read(66,"(I5.5)") fgene%N_mu
  allocate(fgene%mu(fgene%N_mu))
  do imu = 1, fgene%N_mu
    read(66,"(E19.12E2)") fgene%mu(imu)
  end do
  close(66)
  fgene%mu_min =  fgene%mu(1)
  fgene%mu_max =  fgene%mu(fgene%N_mu)
  filename = trim(data_folder) // dst_data_folder // "/grhop.dat"
  open(67, file = trim(filename))
  read(67,"(I5.5)") fgene%N_rhop
  allocate(fgene%rhop(fgene%N_rhop))
  allocate(fgene%g(fgene%N_rhop, fgene%N_vpar, fgene%N_mu)) !
  print*,"Loading GENE distribution function"
  write(format,"(A1I3.3A12)") "(",fgene%N_mu,"(E15.8E2A1))"
  do irhop = 1,fgene%N_rhop
    write(irhop_str, "(I3.3)") irhop - 1
    cur_filename = trim(data_folder) // dst_data_folder // "/gvpar" // trim(irhop_str) // ".dat"
    open(66, file= trim(cur_filename))
    read(67,"(E19.12E2)") fgene%rhop(irhop)
    do ivpar = 1, fgene%N_vpar
      read(66, *) (fgene%g(irhop, ivpar, imu), imu = 1, fgene%N_mu) ! imu = 7, ffp%N_mu - 6
    end do
    close(66)
  end do
  close(67)
   fgene%rhop_min =  fgene%rhop(1)
   fgene%rhop_max = fgene%rhop(fgene%N_rhop) ! set boundary in a way that we can fall back to thermal
  cur_filename = trim(data_folder) // dst_data_folder // "/f0" // ".dat"
  open(66, file= trim(cur_filename))
  ! load the background into the g_inter matrix and make the spline for it
  allocate(f_0(fgene%N_vpar, fgene%N_mu))
  do ivpar = 1, fgene%N_vpar
    read(66, *) (f_0(ivpar, imu), imu = 1, fgene%N_mu) ! imu = 7, ffp%N_mu - 6
  end do
  close(66)
  call make_rect_spline(fgene%f0_spl, int(size(fgene%vpar),4), int(size(fgene%mu),4), &
                        fgene%vpar, fgene%mu, f_0, m_max = N_absz_large)
#ifdef NAG
  call nag_spline_2d_interp(fgene%vpar, fgene%mu, f_0, fgene%f0_nag_spl)
#endif
  deallocate(f_0)
  filename = trim(data_folder) // dst_data_folder // "/B0.dat"
  open(66, file = trim(filename))
  read(66,"(E19.12E2)") fgene%B0
  close(66)
  call setup_fgene_rhop_splines(fgene)
end subroutine read_fgene_data

subroutine read_Gene_Bi_max_data()
use mod_ecfm_refr_types, only: fgene, data_folder
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop
Character(200)  :: cur_filename
Character(1)    :: blanc
  cur_filename =  trim(data_folder) // "Te_perp.dat"
  open(67, file = trim(cur_filename))
  read(67,"(I5.5)") N_rhop
  allocate(rhop(N_rhop), fgene%Te_perp_vec(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(E19.12E2A1E19.12E2)") rhop(irhop), blanc, fgene%Te_perp_vec(irhop)
  end do
  close(67)
  if(abs(minval(rhop) - fgene%rhop_min) > 1.e-5 .or. &
    abs(maxval(rhop) - fgene%rhop_max) > 1.e-5 ) then
    print*, "rhop range for Te_perp not sufficiently similar to rhop range from grid based BiMaxwellian"
    print*, abs(minval(rhop) - fgene%rhop_min), abs(maxval(rhop) - fgene%rhop_max)
    call abort()
  end if
  call make_1d_spline(fgene%Te_perp_spl, int(N_rhop,4), rhop, fgene%Te_perp_vec, k = 1)
#ifdef NAG
  call nag_spline_1d_interp(rhop, fgene%Te_perp_vec, fgene%Te_perp_nag_spl)
#endif
  cur_filename =  trim(data_folder) // "Te_par.dat"
  open(67, file = trim(cur_filename))
  read(67,"(I5.5)") N_rhop
  deallocate(rhop)
  allocate(rhop(N_rhop), fgene%Te_par_vec(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(E19.12E2A1E19.12E2)") rhop(irhop), blanc, fgene%Te_par_vec(irhop)
  end do
  close(67)
  if(abs(minval(rhop) - fgene%rhop_min) > 1.e-5 .or. &
    abs(maxval(rhop) - fgene%rhop_max) > 1.e-5 ) then
    print*, "rhop range for Te_par not sufficiently similar to rhop range from grid based BiMaxwellian"
    print*, abs(minval(rhop) - fgene%rhop_min), abs(maxval(rhop) - fgene%rhop_max)
    call abort()
  end if
  ! Interpolate the logarithm for positivity
  call make_1d_spline(fgene%Te_par_spl, int(N_rhop,4), rhop, fgene%Te_par_vec, k = 1)
#ifdef NAG
  call nag_spline_1d_interp(rhop, fgene%Te_par_vec, fgene%Te_par_nag_spl)
#endif
  deallocate(rhop)
end subroutine read_Gene_Bi_max_data

subroutine read_Spitzer_data()
use mod_ecfm_refr_types, only: data_folder, Spitzer
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
integer(ikind)  :: irhop, N_rhop
real(rkind), dimension(:), allocatable :: rhop, j
Character(200)  :: cur_filename
Character(1)    :: blanc
  cur_filename =  trim(data_folder) // "j.dat"
  open(67, file = trim(cur_filename))
  read(67,"(I5.5)") N_rhop
  allocate(rhop(N_rhop), j(N_rhop))
  do irhop = 1, N_rhop
    read(67,"(E19.12E2A1E19.12E2)") rhop(irhop), blanc, j
  end do
  close(67)
  call make_1d_spline(Spitzer%j_spl, int(N_rhop,4), rhop, j)
#ifdef NAG
  call nag_spline_1d_interp(rhop, j, Spitzer%j_nag_spl)
#endif
  deallocate(rhop,j)
end subroutine read_Spitzer_data

subroutine ffp_clean_up()
  use mod_ecfm_refr_types, only: ffp
  implicit none
  deallocate(ffp%f)
  deallocate(ffp%u)
  deallocate(ffp%pitch)
  deallocate(ffp%rhop)
  deallocate(ffp%B_min)
  deallocate(ffp%rhop_B_min)
!  deallocate(ffp%f_inter)
end subroutine ffp_clean_up

subroutine fgene_clean_up()
  use mod_ecfm_refr_types, only: fgene
  implicit none
  deallocate(fgene%g)
  deallocate(fgene%vpar)
  deallocate(fgene%mu)
  deallocate(fgene%rhop)
end subroutine fgene_clean_up

subroutine export_all_ece_data()
USE mod_ecfm_refr_types, only: ant, rad, data_folder, dstf, N_ray, N_freq, mode_cnt, plasma_params, ffp
implicit none
Character(200)                         :: cur_filename, filename
character(12)                          :: ich_str
Integer(ikind)                         :: idiag, ich, ir, ifreq, iint, imode, N_ch_prior
character(1)   :: sep
sep = " "
filename =  trim(data_folder) // "cnt.dat"
open(67,file=trim(filename))
!filename = trim(data_folder) // "rhopres.dat"
!open(74,file=filename)
!filename = trim(data_folder) // "sres.dat"
!open(76,file=trim(filename))
filename = trim(data_folder) // "f_ECE.dat"
open(79,file=trim(filename))
filename = trim(data_folder) // "diag.dat"
open(84,file=filename)
filename = trim(data_folder) // "parms.dat"
open(80,file=trim(filename))
N_ch_prior = 0
ifreq = 1
ir = 1
N_ray = 1
N_freq = 1
do idiag = 1, ant%N_diag
  do ich = 1, ant%diag(idiag)%N_ch
    write(67,"(I5)") rad%diag(idiag)%ch(ich)%mode(1)%ray(ir)%freq(ir)%total_LOS_points
    write(84,"(A3)") ant%diag(idiag)%diag_name
    do imode = 1, mode_cnt
      ! export only central ray and central freq
      write(ich_str, "(I3.3)") ich + N_ch_prior ! ((ich - 1) * N_freq + ifreq)
      if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
        cur_filename = trim(data_folder) // "chdata" // trim(ich_str) // ".dat"
      else
        cur_filename = trim(data_folder) // "chOdata" // trim(ich_str) // ".dat"
      end if
      open(66, file=cur_filename)
      do iint = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
        write(66,"(9(E18.10E3A1))") rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%s," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%R," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%z," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%rhop," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%ne," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%Te," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%theta," ",&
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%v_g_perp, " ", &
           rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%freq_2X," "
      end do !iint
      close(66)
    end do ! imode
!    write(74,"(E18.10E3)") rad%diag(idiag)%ch(ich)%rhop_res
!    write(76,"(E18.10E3)") rad%diag(idiag)%ch(ich)%s_res
    write(79,"(E18.10E3)") ant%diag(idiag)%ch(ich)%f_ECE
  end do
  N_ch_prior =  N_ch_prior + ant%diag(idiag)%N_ch * N_freq
end do
write(80,"(E18.10E3)") plasma_params%B_ax
!close(74)
!close(76)
close(79)
close(80)
if(dstf == "numeric") then
  filename = trim(data_folder) // "B_min.dat"
  open(90, file=filename)
  do iint = 1, size(ffp%rhop_B_min)
    write(90, fmt="(E18.10E3A1E18.10E3)") ffp%rhop_B_min(iint), " ", ffp%B_min(iint)
  end do
  close(90)
end if
end subroutine export_all_ece_data

end module mod_ecfm_refr_utils
