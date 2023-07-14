module mod_ECRad_IMAS

contains

subroutine pre_initialize_ECRad_IMAS(codeparam_string, wall, &
                                     error_flag, error_message)
  use mod_ECRad, only: pre_initialize_ECRad_f2py
  use xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, &
                        type_xml2eg_document, xml2eg_free_doc, set_verbose
  use ids_schemas, only: ids_wall, ids_is_valid
  implicit none
  ! Input/Output
  character(len=132), pointer, intent(in) :: codeparam_string(:)
  type(ids_wall):: wall
  integer, intent(out) :: error_flag
  character(len=:), pointer, intent(out) :: error_message
  ! Internal
  type(type_xml2eg_document) :: doc
  real(kind=8)     :: ratio_for_3rd_harm, tau_ignore, &
                      reflec_X, reflec_O, mode_conv, &
                      large_ds, small_ds
  character(2)     :: dstf
  integer   :: considered_modes, max_points_svec, N_pts_BPD, N_ray, &
               N_freq, N_max, N_vessel
  logical    :: extra_output, ripple, ray_tracing, weak_rel
  real(kind=8), dimension(:), allocatable :: vessel_R, vessel_z
  ! Parse the "codeparam_string". This means that the data is put into a document "doc"
  call xml2eg_parse_memory(codeparam_string,doc)
  call set_verbose(.TRUE.) ! Only needed if you want to see what's going on in the parsing
   
  call xml2eg_get(doc,'extra_output', extra_output)
  call xml2eg_get(doc,'dstf', dstf)
  call xml2eg_get(doc,'ray_tracing', ray_tracing)
  call xml2eg_get(doc,'ripple', ripple)
  call xml2eg_get(doc,'weak_rel', weak_rel)
  call xml2eg_get(doc,'ratio_for_3rd_harm', ratio_for_3rd_harm)
  call xml2eg_get(doc,'N_max', N_max)
  call xml2eg_get(doc,'tau_ignore', tau_ignore)
  call xml2eg_get(doc,'considered_modes', considered_modes)
  call xml2eg_get(doc,'reflec_X', reflec_X)
  call xml2eg_get(doc,'reflec_O', reflec_O)
  call xml2eg_get(doc,'max_points_svec', max_points_svec)
  call xml2eg_get(doc,'N_pts_BPD', N_pts_BPD)
  call xml2eg_get(doc,'mode_conv', mode_conv)
  call xml2eg_get(doc,'large_ds', large_ds)
  call xml2eg_get(doc,'small_ds', small_ds)
  call xml2eg_get(doc,'N_ray', N_ray)
  call xml2eg_get(doc,'N_freq', N_freq)
  call xml2eg_free_doc(doc)
  if (ids_is_valid(wall%description_2d(1)%limiter%unit(1)%outline%R)) then
    N_vessel = size(wall%description_2d(1)%limiter%unit(1)%outline%R)
    allocate(vessel_R(N_vessel), vessel_z(N_vessel))
    vessel_R = wall%description_2d(1)%limiter%unit(1)%outline%R
    vessel_z = wall%description_2d(1)%limiter%unit(1)%outline%z
    error_flag = 0
  else
    error_flag = -1
    allocate(character(50):: error_message)
    error_message = 'Error in pre_initialize_ECRad_IMAS: input IDS not valid'
    return
  end if
  call pre_initialize_ECRad_f2py(extra_output, dstf, ray_tracing, ripple, &
                                 1.2d0, weak_rel, &
                                 ratio_for_3rd_harm, N_max, tau_ignore, &
                                 considered_modes, reflec_X, reflec_O, 0, &
                                 max_points_svec, N_pts_BPD ,&
                                 mode_conv, &
                                 1.d0, 1.d0, &
                                 large_ds, small_ds, 0.d0, &
                                 0.d0, N_ray, N_freq, .true., N_vessel, vessel_R, vessel_z)
  deallocate(vessel_R, vessel_z)
end subroutine pre_initialize_ECRad_IMAS

subroutine set_ece_ECRad_IMAS(ece, itime, error_flag, error_message)
  use mod_ECRad, only: prepare_ECE_diag_f2py
  use ids_schemas, only: ids_equilibrium, ids_ece, ids_is_valid
  implicit None
  ! Input/Output
  type(ids_ece):: ece
  integer, intent(in) :: itime
  integer, intent(out) :: error_flag
  character(len=:), pointer, intent(out) :: error_message
  ! Internal
  real(kind=8), dimension(:), allocatable :: f, df, R, phi, z, phi_tor, theta_pol, &
                                             dist_focus, width, pol_coeff, x1_vec, x2_vec
  integer :: i, N_ch
  ! CHECK IF INPUT IDS IS VALID
  if (ids_is_valid(ece%line_of_sight%first_point%r)) then
    N_ch = size(ece%channel)
    allocate(f(N_ch), df(N_ch), R(N_ch), phi(N_ch), z(N_ch), &
    phi_tor(N_ch), theta_pol(N_ch), dist_focus(N_ch), width(N_ch), pol_coeff(N_ch), &
              x1_vec(2), x2_vec(2))
    do i = 1, N_ch
      f(i) = ece%channel(i)%frequency%data(itime)
      df(i) =  ece%channel(i)%if_bandwidth
      R(i) = ece%line_of_sight%first_point%r
      phi(i) = ece%line_of_sight%first_point%phi
      z(i) = ece%line_of_sight%first_point%z
      x1_vec(1) = ece%line_of_sight%first_point%r * cos(ece%line_of_sight%first_point%phi)
      x1_vec(2) = ece%line_of_sight%first_point%r * sin(ece%line_of_sight%first_point%phi)
      x2_vec(1) = ece%line_of_sight%second_point%r * cos(ece%line_of_sight%first_point%phi)
      x2_vec(2) = ece%line_of_sight%second_point%r * sin(ece%line_of_sight%first_point%phi)
      theta_pol(i) = atan2((ece%line_of_sight%second_point%z - ece%line_of_sight%first_point%z), &
                           -(ece%line_of_sight%second_point%r - ece%line_of_sight%first_point%r))
      phi_tor(i) = -acos((-x1_vec(1) * (x2_vec(1) - x1_vec(1)) - x1_vec(2) * (x2_vec(2) - x1_vec(2))) / &
                         (R(i) * sqrt(sum((x2_vec - x1_vec)**2))))
      width(i) = 0.1d0
                !  (ece%channel(i)%beam%spot%size%data(1,itime) + &
                !   ece%channel(i)%beam%spot%size%data(1,itime)) / 2.d0 ! Average the ellipse to a circle
      dist_focus(i) = 1.d0
                      ! -(ece%channel(i)%beam%phase%curvature%data(1, itime) + &
                      !   ece%channel(i)%beam%phase%curvature%data(1, itime)) / 2.d0  ! Average the ellipse to a circle
      pol_coeff(i) = -1
      error_flag = 0
    end do
    call prepare_ECE_diag_f2py(f, df, R, phi, z, phi_tor, theta_pol, dist_focus, width, pol_coeff)
    deallocate(f, df, R, phi, z, phi_tor, theta_pol, dist_focus, width, pol_coeff, x1_vec, x2_vec)
  else
      ! ERROR IF THE CODE DOES NOT COMPLETE TO THE END
      error_flag = -1
      allocate(character(50):: error_message)
      error_message = 'Error in set_ece_ECRad_IMAS: input IDS not valid'

  endif

end subroutine

#ifdef OMP
subroutine set_ECRad_thread_count(num_threads)
use mod_ECRad, only: set_omp_threads_ECRad_f2py
implicit None
  integer, intent(in) :: num_threads
  call set_omp_threads_ECRad_f2py(num_threads)
end subroutine set_ECRad_thread_count
#endif

subroutine reset_ECRad()
use mod_ECRad,      only: clean_up_ECRad
implicit none
  call clean_up_ECRad()
end subroutine reset_ECRad

subroutine initialize_ECRad_IMAS(equilibrium, itime, error_flag, error_message)
! Hence, to keep the structure similiar all initizalization is performed here
! Initializations that depend on time are done here
use mod_ECRad,  only: initialize_ECRad_f2py
use ids_schemas, only: ids_equilibrium, ids_ece, ids_is_valid
  ! Input/Output
implicit none
type(ids_equilibrium):: equilibrium
integer(kind=4), intent(in)                       :: itime
integer, intent(out) :: error_flag
character(len=:), pointer, intent(out) :: error_message
real(kind=8)                            :: R_ax, z_ax, psi_ax, psi_sep
real(kind=8), dimension(:), allocatable :: R, z
real(kind=8), dimension(:,:), allocatable :: rhop, Br, Bt, Bz
integer(kind=4) :: i, m, n
  if (ids_is_valid(equilibrium%ids_properties%homogeneous_time) .and. &
      size(equilibrium%time)>0) then
    m = size(equilibrium%time_slice(itime)%profiles_2d(1)%grid%dim1)
    n = size(equilibrium%time_slice(itime)%profiles_2d(1)%grid%dim2)
    allocate(R(m), z(n), rhop(m,n), Br(m,n), Bt(m,n), Bz(m,n))
    R(:) = equilibrium%time_slice(itime)%profiles_2d(1)%grid%dim1(:)
    z(:) = equilibrium%time_slice(itime)%profiles_2d(1)%grid%dim2(:)
    psi_ax = equilibrium%time_slice(itime)%global_quantities%psi_axis
    psi_sep = equilibrium%time_slice(itime)%global_quantities%psi_boundary
    R_ax = equilibrium%time_slice(itime)%global_quantities%magnetic_axis%r
    z_ax = equilibrium%time_slice(itime)%global_quantities%magnetic_axis%z
    do i = 1,n
      rhop(:,i) = equilibrium%time_slice(itime)%profiles_2d(1)%psi(:,i)
      rhop(:,i) = sqrt((rhop(:,i) - psi_ax) / (psi_sep - psi_ax))
      Br(:,i) = equilibrium%time_slice(itime)%profiles_2d(1)%b_field_r(:,i)
      Bt(:,i) = equilibrium%time_slice(itime)%profiles_2d(1)%b_field_tor(:,i)
      Bz(:,i) = equilibrium%time_slice(itime)%profiles_2d(1)%b_field_z(:,i)
    end do
    call initialize_ECRad_f2py(0, 0, R, z, rhop, Br, Bt, Bz, R_ax, z_ax)
    deallocate(R, z, rhop, Br, Bt, Bz)
    error_flag = 0
  else
    ! ERROR IF THE CODE DOES NOT COMPLETE TO THE END
    error_flag = -1
    allocate(character(50):: error_message)
    error_message = 'Error in set_ece_ECRad_IMAS: input IDS not valid'
  endif               
end subroutine initialize_ECRad_IMAS

subroutine get_rho_pol(N, psi, rho_pol)
  implicit None
  integer(kind=4), intent(in) ::  N
  real(kind=8), dimension(:), intent(in) :: psi
  real(kind=8), dimension(:), intent(out) :: rho_pol
  rho_pol = sqrt((psi - psi(1))/ &
                 (psi(N) - psi(1)))
end subroutine

subroutine make_rays_ECRad_IMAS(core_profiles, itime, ece)
use ids_schemas, only: ids_core_profiles, ids_ece
use mod_ECRad,        only: make_rays_ECRad_f2py
implicit none
type(ids_core_profiles), intent(in) :: core_profiles
integer(kind=4), intent(in) :: itime
type(ids_ece), intent(inout) :: ece
integer                      :: N_psi
real(kind=8), dimension(:), allocatable :: rho_pol
N_psi = size(core_profiles%profiles_1d(itime)%grid%psi)
allocate(rho_pol(N_psi))
if(size(core_profiles%profiles_1d(itime)%grid%rho_pol_norm) == 0) then
  call get_rho_pol(N_psi, core_profiles%profiles_1d(itime)%grid%psi, rho_pol)
else
  rho_pol = core_profiles%profiles_1d(itime)%grid%rho_pol_norm
end if
call make_rays_ECRad_f2py(rhop_knots_ne=rho_pol, &
                          n_e=core_profiles%profiles_1d(itime)%electrons%density, &
                          rhop_knots_Te=rho_pol, &
                          T_e=core_profiles%profiles_1d(itime)%electrons%temperature,
                          rhop_res=ece%channel(:)%position%rho_tor_norm)
deallocate(rho_pol)
end subroutine make_rays_ECRad_IMAS


subroutine make_dat_model_ECRad(core_profiles, itime, ece)
use ids_schemas, only: ids_core_profiles, ids_ece
use mod_ECRad,        only: get_N_ch, make_dat_model_ece_ECRad_f2py
implicit none
type(ids_core_profiles), intent(in) :: core_profiles
type(ids_ece), intent(inout) :: ece
integer(kind=4), intent(in) :: itime
integer                      :: N_ch, ich, N_psi
real(kind=8), dimension(:), allocatable :: dat_model_ece, rho_pol
logical, dimension(:), allocatable :: ece_fm_flag_ch
N_ch = get_N_ch()
allocate(dat_model_ece(N_ch), ece_fm_flag_ch(N_ch))
ece_fm_flag_ch(:) = .true.
N_psi = size(core_profiles%profiles_1d(itime)%grid%psi)
allocate(rho_pol(N_psi))
if(size(core_profiles%profiles_1d(itime)%grid%rho_pol_norm) == 0) then
  call get_rho_pol(N_psi, core_profiles%profiles_1d(itime)%grid%psi, rho_pol)
else
  rho_pol = core_profiles%profiles_1d(itime)%grid%rho_pol_norm
end if
call make_dat_model_ece_ECRad_f2py(rhop_knots_ne=rho_pol, &
                                   n_e=core_profiles%profiles_1d(itime)%electrons%density, &
                                   rhop_knots_Te=rho_pol, &
                                   T_e=core_profiles%profiles_1d(itime)%electrons%temperature, &
                                   ne_rhop_scal=1.d0, reflec_X_new=-1.d0, & ! in
                                   reflec_O_new=-1.d0, ece_fm_flag_ch=ece_fm_flag_ch, &
                                   rp_min=minval(core_profiles%profiles_1d(itime)%grid%rho_pol_norm), &
                                   dat_model_ece=dat_model_ece, set_grid_dynamic=.false., &
                                   verbose=.false.)
if(.not. associated(ece%channel)) then
  allocate(ece%channel(N_ch))
end if
ece%ids_properties%homogeneous_time=1
allocate(ece%time(1))
ece%time(1) = core_profiles%time(itime)
do ich = 1, N_ch
  if(.not. associated(ece%channel(ich)%T_e%data)) then
    allocate(ece%channel(ich)%T_e%data(1))
  end if
  ece%channel(ich)%T_e%data(1) = dat_model_ece(ich)
end do
deallocate(dat_model_ece, ece_fm_flag_ch, rho_pol)
end subroutine make_dat_model_ECRad
end module mod_ECRad_IMAS


