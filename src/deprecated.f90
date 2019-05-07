! All the deprecated source code ends up here.
! There individual code chunks are sorted according to their original module
!!!! mod_ecfm_refr_raytrayce_initialize !!!!

#ifdef AUG
  subroutine setup_eq_from_shotfile(plasma_params, R, z, rhop, Br, Bt, Bz, R0_out, Btf0)
    use f90_kind
    use mod_ecfm_refr_types,        only: output_level, m_eq, n_eq, plasma_params_type
    use aug_db_routines,            only: read_mbi, aug_kkEQpfx, aug_kkrzBrzt
    use eqi_mod,                    only: Bt_btf_corrected
#ifdef OMP
    use omp_lib
#endif
    implicit none
    type(plasma_params_type), intent(inout)     :: plasma_params
    real(rkind), dimension(:), allocatable, intent(out) :: R, z
    real(rkind), dimension(:,:), allocatable, intent(out) ::  Br, Bt, Bz
    real(rkind), dimension(:,:), allocatable, intent(out) ::  rhop
    real(rkind)                            , intent(out)  :: R0_out, Btf0
    real(rkind), dimension(:), allocatable      :: z_mat_helper, R_test, z_test,rhop_ne, &
                                                   rhop_Te, mbi_time, BtfABB
    real(4),     dimension(:), allocatable      :: R_temp, z_temp
    real(4), dimension(:,:), allocatable        :: pfm_temp
    real(rkind), dimension(1)                   :: Rv, zv, Br_vac, Bt_vac, Bz_vac
    integer(ikind)                              :: i, error, debugunit, j, k, irhop_extra
    real(4)                                     :: temp_time
    logical                                     :: debug
    character(18)                               :: btf_fmt_str
    integer(4)                                  :: s_error, m_temp, n_temp, mdim
    real(rkind)                                 :: Btf0_eq, R_mag, z_mag, &
                                                   R_sxp, z_sxp, R0
    R0 = 1.65d0
    R0_out = R0
    s_error = 0
    error     = 0
    debug     = .false.
    debugunit = 999
    temp_time = real(plasma_params%time, 4)
    m_temp = int(m_eq,4)
    n_temp = int(n_eq,4)
    mdim = m_temp + 1
    allocate(R_temp(mdim),z_temp(n_temp))
    allocate(pfm_temp(mdim,n_temp))
    if(output_level) print*, "Requested equilibrium edition", plasma_params%eq_ed
    call kkeqpfm(s_error, plasma_params%eq_exp, plasma_params%eq_diag, & ! in
                    int(plasma_params%shot, 4), plasma_params%eq_ed,& ! inout
                    temp_time, &
                    mdim, m_temp, &
                    n_temp, R_temp,  z_temp, pfm_temp)
    if(output_level) print*, "Obtained edition", plasma_params%eq_ed
    !print*, "kkeqpfm finished"
    !print*, "allocating R/z"
    plasma_params%time = real(temp_time, 8)
    if(output_level) print*,"Time of equilibrium is", plasma_params%time
    allocate(R(m_temp),z(n_temp))
    R = real(R_temp(1:m_temp),8)
    z = real(z_temp(1:n_temp),8)
    !print*, "deallocating R_temp/z_temp"
    deallocate(R_temp, z_temp)
    if(s_error /= 0) then
      print*,"equilibrium EXP and diag ", plasma_params%eq_exp, plasma_params%eq_diag
      print*,"shot and time ", int(plasma_params%shot, 4), temp_time
      print*, "Error when calling kkeqpfm: ", s_error
      stop "kk error in mod_raytrace_initialize"
    end if
    !print*, "allocating z_mat_helper"
    allocate(z_mat_helper(n_temp))
    !print*, "allocating rhop"
    allocate(rhop(m_temp,n_temp))
    do i = 1, m_temp
      rhop(i,:) = real(pfm_temp(i,:),8)
    end do
    !print*, "deallocating pfm_temp"
    !deallocate(pfm_temp)
    !print*, "aug_kkEQpfx starting"
    call aug_kkEQpfx(error,                       & ! out
                   debug, debugunit, plasma_params%shot,       & ! in
                   plasma_params%eq_exp, plasma_params%eq_diag, & ! in
                   plasma_params%eq_ed,                        & ! inout
                   temp_time,                                  & ! inout
                   R_mag  = plasma_params%R_ax,                & ! out [m]  radial   coordinates of magnetic axis
                   z_mag  = plasma_params%z_ax,                & ! out [m]  vertical coordinates of magnetic axis
                   pf_mag = plasma_params%pf_mag,               & ! out [Vs] magnetic flux at magnetic axis
                   R_sxp  = plasma_params%R_sep,                & ! out [m]  radial   coordinates of X-point
                   z_sxp  = plasma_params%z_sep,                & ! out [m]  vertical coordinates of X-point
                   pf_sxp = plasma_params%pf_sxp)                 ! out [Vs] magnetic flux at X-point
    !print*, "aug_kkEQpfx finished"
    rhop = sqrt((rhop - plasma_params%pf_mag)/(plasma_params%pf_sxp - plasma_params%pf_mag))
    !print*, "allocating Br"
    allocate(Br(m_temp, n_temp), &
              Bt(m_temp, n_temp), Bz(m_temp, n_temp))

    ! vacuum toroidal magnetic field to be used
    !------------------------------------------
    ! Note that there is a slight distinction to the mod_eqi routine here, since
    ! we do not employ a median filter
    rv(1) = 2.4
    zv(1) = 0.0
    call aug_kkrzBrzt(error,                                   & ! out
                    debug, debugunit, plasma_params%shot,        & ! in
                    plasma_params%eq_exp, plasma_params%eq_diag, & ! in
                    plasma_params%eq_ed,                         & ! inout
                    plasma_params%time,                          & ! inout
                    r  = rv,                                      & ! in
                    z  = zv,                           & ! in
                    Br = Br_vac,                               & ! out
                    Bz = Bz_vac,                                & ! out
                    Bt = Bt_vac,                                & ! out
                    err_stop = .true.)
    Btf0_eq = Bt_vac(1) * rv(1) / R0
    !print*, Btf0_eq
    if (plasma_params%btf_mode == "BTFABB") then
      call read_mbi(plasma_params%shot, "AUGD", plasma_params%time_beg, time_end=plasma_params%time_end, &
                    time=mbi_time, BTFABB=BtfABB)         ! out
      if(size(BtfABB) < 1) stop "Failed to get Btf or BtfABB at R0 in mod_ecfm_refr_raytrace_initialize"
      Btf0 = sum(BtfABB) / size(BtfABB)
      plasma_params%Btf0 = Btf0
    else if (plasma_params%btf_mode == "old_Btf") then
      Btf0 = Btf0_eq
      plasma_params%Btf0 = Btf0_eq
    else
      print*, plasma_params%btf_mode
      stop "sub Blos_along_beam_path: btf_mode not defined properly!"
    endif
    if (abs(Btf0_eq) > 20.d0 .or. abs(Btf0_eq) < 1.d-3) then
      write(6,'(a,f8.4,a,e12.4)')"Btf0_eq(time =", plasma_params%time, " s) = ", Btf0_eq
      stop 'subroutine Blos_along_beam_path: |Btf0_eq| > 20 T or |Btf0_eq| < 1 mT ???'
    endif
    if (abs(Btf0) > 20.d0 .or. abs(Btf0) < 1.d-3) then
      write(6,'(a,f8.4,a,e12.4)')"Btf0(time =", plasma_params%time, " s) = ", Btf0
      stop 'subroutine Blos_along_beam_path: |Btf0| > 20 T or |Btf0| < 1 mT ???'
    endif
    do j = 1, n_temp
      z_mat_helper(:) = z(j)
      call aug_kkrzBrzt(error,                                   & ! out
                    debug, debugunit, plasma_params%shot,        & ! in
                    plasma_params%eq_exp, plasma_params%eq_diag, & ! in
                    plasma_params%eq_ed,                         & ! inout
                    plasma_params%time,                          & ! in
                    r  = R,                                      & ! in
                    z  = z_mat_helper,                           & ! in
                    Br = Br(:,j),                               & ! out
                    Bz = Bz(:,j),                                & ! out
                    Bt = Bt(:,j),                                & ! out
                    err_stop = .true.)
      call Bt_btf_corrected(plasma_params%shot, plasma_params%btf_mode, plasma_params%btf_corr_fact_ext, &
                                               Btf0_eq, Btf0, R0, R, Bt(:,j))
      if (any(abs(Bt(:,j)) > 20.d0)) then
        do i = 1, size(Bt)
          write(6,'(i3,e12.4)')i, Bt(i,j)
        enddo
        stop 'subroutine setup_eq_from_shotfile: |Bt| > 20 T ???'
      endif
    end do
    deallocate(z_mat_helper)
  end subroutine setup_eq_from_shotfile
#else
  subroutine setup_eq_from_shotfile(plasma_params, R, z, rhop, Br, Bt, Bz, R0_out, Btf0)
      use f90_kind
      use mod_ecfm_refr_types,        only: output_level, m_eq, n_eq, plasma_params_type
      use eqi_mod,                    only: Bt_btf_corrected
      implicit none
      type(plasma_params_type), intent(inout)     :: plasma_params
      real(rkind), dimension(:), allocatable, intent(inout) :: R, z
      real(rkind), dimension(:,:), allocatable, intent(inout) ::  Br, Bt, Bz
      real(rkind), dimension(:,:), allocatable, intent(inout) ::  rhop
      real(rkind)                            , intent(out)  :: R0_out, Btf0
      print*, "Subroutine setup_eq_from_shotfile is AUG specific and must not be used on non-AUG devices"
      call abort()
  end subroutine setup_eq_from_shotfile
#endif

!!! mod_ecfm_refr_raytrace
subroutine prepare_svec_segment_eq_dist_grid(plasma_params, omega, Y_res, svec, ray_segment, total_LOS_points, last_N, &
                                               dist, max_points_svec_reached)
  ! Linearly interpolates R,z on the ray which corresponds
  ! to the s values given on an equidistant grid.
  ! The cold resonance is detected automatically and the step size is choosen correspondingly.
  ! Once the R,z values are obtained the other ray parameters are aqcuired using sub_local_params.
  ! Every time this routine is called the svec values from i_start to i_end will be filled.
  ! i referst to the svec, while N refers to the ray_segment
  ! This routine also conrolls the step size.
  use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, plasma_params_type, &
                                        ray_element_full_type, output_level, max_points_svec
  USE interpolation_routines,     only: linear_interpolation
  use constants,                  only: pi, e0, mass_e, eps0, c0
  use mod_ecfm_refr_utils,        only: sub_remap_coords, binary_search
  use f90_kind
  implicit none
  type(plasma_params_type), intent(in)                                       :: plasma_params
  real(rkind), intent(in)                                                    :: omega
  real(rkind), dimension(:), intent(in)                                      :: Y_res! array holding all resonances that require special treatment
  type(rad_diag_ch_mode_ray_freq_svec_type), dimension(:), intent(inout)     :: svec
  type(ray_element_full_type), dimension(max_points_svec), intent(in)        :: ray_segment !the ray
  integer(ikind), intent(out)                                                :: total_LOS_points
  integer(ikind),   intent(in)                                               :: last_N
  real(rkind), dimension(:), intent(in)                                      :: dist ! distance to interpolated, starting point of interpolation
  logical, intent(out), optional                                             :: max_points_svec_reached
  real(rkind)                                                                :: a, b
  integer(ikind)                                                             :: grid_size
  logical                                                                    :: finished_ray
  real(rkind)                                                                :: n_e, Y_last, omega_c, N_abs_1, N_abs_2 ! not safed but important for interpolation
  integer(ikind)                                                             :: i, N, j, k, ifail, dist_last_N, last_good_N, last_grid_size
  real(rkind), dimension(3)                                                  :: x_vec, N_vec, R_vec, B_vec
  logical                                                                    :: dist_ok, extra_output, close_to_resonance, large_skips_resonance
  integer(ikind)                                                             :: ires, step_lower, step_upper, first_N, last_first_N
  first_N = 1
  last_first_N  = 1
  finished_ray = .false.
  a = 0.d0
  i = 1
  grid_size = 1 ! Large
  extra_output = .false.
  large_skips_resonance = .false.
  if(extra_output) print*, "Interpolating things"
  do while(.not. finished_ray)
    if(i + plasma_params%int_step_cnt - 1 >=  max_points_svec) then
      if(present(max_points_svec_reached)) then
        max_points_svec_reached = .true.
        return
      else
        print*, "Failed to span grid for radiation transport due to an insufficient amount of points"
        print*, "Increase max points on LOS or decrease resolution (small and large step size) and rerun"
        call abort()
      end if
    end if
    close_to_resonance = .false.
    large_skips_resonance = .false.
    ! iterate over the grid size to find the optimal one
    b = a + dist(grid_size)
    do while(.true.)
      if(b >= ray_segment(last_N - 1)%s) then ! to assure we have at least two points in the next step - 2
        dist_last_N = last_N
        b = ray_segment(last_N)%s
        finished_ray = .true.
      else
        dist_last_N = binary_search(ray_segment(:)%s, b, first_N, last_N)
        if(ray_segment(dist_last_N)%s < b) dist_last_N = dist_last_N  + 1
        !! the binary search returns the smaller point, since dist_last_N gives the upper limit we want the next larger s value
        ! If a ray point coincides with b the increment is, however, not neccessary
        if(dist_last_N <= 1 .or. dist_last_N > last_N) then ! .or.
          dist_last_N = binary_search(ray_segment(:)%s, a + dist(grid_size), first_N, last_N, .true.)
          print*, "Something wrong with the ray when looking for first point"
          print*, "s value required", a + dist(grid_size)
          print*, "ray s values", ray_segment(1:last_N)%s
          print*, "first_N, last_N, dist_last_N, last_first_N: ", first_N, last_N, dist_last_N, last_first_N
          stop "fatal error in prepare_svec_segment svec segment in mod_raytrace.f90"
        else if(ray_segment(dist_last_N)%s < b .or. a < ray_segment(first_N)%s) then
          dist_last_N = binary_search(ray_segment(:)%s, a + dist(grid_size), first_N, last_N, .true.)
          print*, "b, larger ray value ", b, ray_segment(dist_last_N + 1)%s
          print*, "a, first ray point ", a, ray_segment(first_N)%s
          print*, "first_N, last_N, dist_last_N, last_first_N: ", first_N, last_N, dist_last_N, last_first_N
          stop "fatal error in prepare_svec_segment svec segment in mod_raytrace.f90"
        end if
      end if
      close_to_resonance = .false.
      do ires = 1, size(Y_res)
        if(any(ray_segment(first_N:dist_last_N)%omega_c / omega > Y_res(ires) - 0.01d0) .and. & ! near cold resonance is doppler-shifted regime
           any(ray_segment(first_N:dist_last_N)%omega_c / omega < Y_res(ires) + 0.02d0) .and. & ! Hence, more up-shift than down-shift
           any(ray_segment(first_N:dist_last_N)%rhop < plasma_params%rhop_emit)) then ! no fine grid in SOL
           close_to_resonance = .true.
           exit
         end if
      end do
      if(grid_size == 2 .and. (close_to_resonance .or. large_skips_resonance)) then
        exit
      else if(grid_size == 1 .and. (close_to_resonance .or. large_skips_resonance)) then
        large_skips_resonance = .True.
        grid_size = 2
        cycle
      else if(grid_size == 1 .and. .not. close_to_resonance .and. .not. large_skips_resonance) then
         exit
      else if (grid_size == 2 .and. .not. close_to_resonance .and. .not. large_skips_resonance) then
        grid_size = 1
        cycle
      else
        print*, "Looks like I missed a possibility"
        print*, "grid_size", grid_size
        print*, "large_skips_resonance", large_skips_resonance
        print*, "close_to_resonance", close_to_resonance
        call abort()
      end if
    end do
    if(extra_output) print*,"first ray point, a", ray_segment(1)%s, a
    svec(i:i + plasma_params%int_step_cnt - 1)%s = (b - a) * plasma_params%Int_absz(:) + a
    !print*, "s, a,b , absz",svec(i_start)%s, a, b, plasma_params%Int_absz(1)
    N = first_N ! move to first_N - 1 instead of first_N because we need overlap here to cover the entire interval from 0.0 - 1.0
    ! For Gauss the overlap is not neccessary, because we only need 0.01 - 0.99
    svec(:)%s = svec(:)%s - svec(1)%s ! the LOS coorindate is shifted
                                      ! i.e. point svec(i)%s corresponds to
                                      ! svec(i+ 1)%s
                                      ! This line fixes this and svec(1)%s  = 0.d0
    last_first_N = first_N
    first_N = dist_last_N - 1 ! Usually many interpolations are performed for one ray. Hence, we want to remember our last position in the ray.
    !print*, "first_N, last_first_N", first_N, last_first_N
    !print*, "b, first_N", b, ray_segment(first_N)%s
    a = b
    i = i + plasma_params%int_step_cnt
  end do! finished ray
  total_LOS_points = i - 1
  end subroutine prepare_svec_segment_eq_dist_grid

subroutine preprare_svec_spline(ray_segment, svec_spline, rad_ray_freq, freq)
  use mod_ecfm_refr_types,        only: ray_element_full_type, rad_diag_ch_mode_ray_freq_svec_spline_type, &
                                        rad_diag_ch_mode_ray_freq_type, spl_type_1d, output_level
  use f90_kind
  use constants,                  only: pi,e0, mass_e
  use mod_ecfm_refr_interpol,   only: make_1d_spline,  spline_1d, spline_1d_get_roots, deallocate_1d_spline
  implicit none
  type(ray_element_full_type), dimension(:), intent(in)         :: ray_segment
  type(rad_diag_ch_mode_ray_freq_svec_spline_type), intent(out) :: svec_spline
  type(rad_diag_ch_mode_ray_freq_type), intent(inout)           :: rad_ray_freq
  real(rkind), intent(in)                                       :: freq
  integer(ikind)                                                :: N_ray_segment, root_cnt, N_min, N_max
  real(rkind), dimension(size(ray_segment))                     :: s, temp_array
  type(spl_type_1d)                                             :: res_spline
  real(rkind), dimension(100)                                   :: roots
  real(rkind)                                                   :: x_res, y_res
  call deallocate_1d_spline(svec_spline%x)
  call deallocate_1d_spline(svec_spline%y)
  call deallocate_1d_spline(svec_spline%z)
  call deallocate_1d_spline(svec_spline%rhop)
  call deallocate_1d_spline(svec_spline%Nx)
  call deallocate_1d_spline(svec_spline%Ny)
  call deallocate_1d_spline(svec_spline%Nz)
  call deallocate_1d_spline(svec_spline%Bx)
  call deallocate_1d_spline(svec_spline%By)
  call deallocate_1d_spline(svec_spline%Bz)
  call deallocate_1d_spline(svec_spline%theta)
  call deallocate_1d_spline(svec_spline%freq_2X)
  call deallocate_1d_spline(svec_spline%Te)
  call deallocate_1d_spline(svec_spline%ne)
  call deallocate_1d_spline(svec_spline%v_g_perp)
  N_min = 1
  do while(ray_segment(N_min)%rhop == -1.d0)
    N_min = N_min + 1
  end do
  N_max = size(ray_segment)
  do while(ray_segment(N_max)%rhop == -1.d0)
    N_max = N_max - 1
  end do
  if(N_max <= N_min) then
    print*, "LOS outside the domain where Te, ne and the equilibrium is defined"
    print*, "Check launch geometry"
    call abort()
  end if
  N_ray_segment = N_max - N_min + 1
  svec_spline%s_min = ray_segment(N_min)%s
  svec_spline%s_max = ray_segment(N_max)%s
  s(:N_ray_segment) = ray_segment(N_min:N_max)%s
  if(any(ray_segment(N_min:N_max)%rhop == -1.d0)) then
    print*, "For the splines there must not be any los point outside the domain where Te, ne and the equilibrium is defined"
    print*, "s", ray_segment(:)%s
    print*, "rhop", ray_segment(:)%rhop
    call abort()
  end if
  temp_array = ray_segment(N_min:N_max)%x_vec(1)! supppresses array temporary warning
  call make_1d_spline(svec_spline%x, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%x_vec(2)
  call make_1d_spline(svec_spline%y, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%x_vec(3)
  call make_1d_spline(svec_spline%z, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%rhop
  call make_1d_spline(svec_spline%rhop, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%N_vec(1)
  call make_1d_spline(svec_spline%Nx, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%N_vec(2)
  call make_1d_spline(svec_spline%Ny, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%N_vec(3)
  call make_1d_spline(svec_spline%Nz, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%B_vec(1)
  call make_1d_spline(svec_spline%Bx, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%B_vec(2)
  call make_1d_spline(svec_spline%By, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%B_vec(3)
  call make_1d_spline(svec_spline%Bz, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%n_e
  call make_1d_spline(svec_spline%ne, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%T_e
  call make_1d_spline(svec_spline%Te, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = pi - ray_segment(N_min:N_max)%theta
  call make_1d_spline(svec_spline%theta, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%omega_c / pi
  call make_1d_spline(svec_spline%freq_2X, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  temp_array = ray_segment(N_min:N_max)%omega_c / pi - freq
  call make_1d_spline(res_spline, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  call spline_1d_get_roots(res_spline, roots, root_cnt)
  call deallocate_1d_spline(res_spline)
  if(root_cnt > 0) then
    rad_ray_freq%s_res = roots(root_cnt)
    call spline_1d(svec_spline%x, rad_ray_freq%s_res, x_res)
    call spline_1d(svec_spline%y, rad_ray_freq%s_res, y_res)
    call spline_1d(svec_spline%z, rad_ray_freq%s_res, rad_ray_freq%z_res)
    call spline_1d(svec_spline%rhop, rad_ray_freq%s_res, rad_ray_freq%rhop_res)
    rad_ray_freq%R_res = sqrt(x_res**2 + y_res**2)
  else
    if(output_level) print*, "No cold resonance found -> channel in cut_off"
    rad_ray_freq%s_res = -1.d0
  end if
  temp_array = ray_segment(N_min:N_max)%v_g_perp
  call make_1d_spline(svec_spline%v_g_perp, int(N_ray_segment, 4), s(:N_ray_segment), temp_array(:N_ray_segment))
  svec_spline%s_center = ray_segment(minloc(ray_segment(N_min:N_max)%rhop, dim=1))%s
  end subroutine preprare_svec_spline

  subroutine create_svec_splines(plasma_params)
  ! Creates splines of svec for all diags
  ! Also allocates all vectors related to output_level = .true.
  use mod_ecfm_refr_types,        only: rad, ant, rad_diag_type, plasma_params_type, ray_out_folder, &
                                        N_ray, N_freq, modes, mode_cnt, output_level, pnts_BPD, &
                                        max_points_svec, Hamil, straight, largest_svec, ray_element_full_type
  use f90_kind
  use constants,                  only: pi,e0, mass_e
  implicit none
  type(plasma_params_type), intent(inout)                       :: plasma_params
  type(ray_element_full_type), dimension(max_points_svec)       :: ray_segment
  integer(ikind)                                                :: idiag, last_N, grid_size, ich, ir, &
                                                                   ifreq, imode, i, N, N_init, ich_tot, &
                                                                   mode
  real(rkind), dimension(2)                                     :: dist
  logical                                                       :: finished_ray
  real(rkind)                                                   :: a, b, omega, temp, X, Y, N_cold
  character(200)                                                :: cur_filename
  character(12)                                                 :: ich_tot_str
  integer(ikind)                                                :: wall_hits
  logical                                                       :: LOS_end
  if(output_level) then
    if(.not. straight) then
      print*, "Preparing LOS including refraction - this will take a moment"
    else
      print*, "Preparing straight LOSs - this should take only a moment"
    end if
  end if
  wall_hits = 0
  largest_svec = 0
  ich_tot = 1
  do idiag = 1, ant%N_diag
    do ich = 1, ant%diag(idiag)%N_ch
      if(output_level) then
        if(.not. allocated(rad%diag(idiag)%ch(ich)%mode_extra_output)) then
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(mode_cnt))
        end if
      end if
      do imode = 1, mode_cnt
        if(Hamil == "Dani") then
          mode = rad%diag(idiag)%ch(ich)%mode(imode)%mode
        else
          mode = -rad%diag(idiag)%ch(ich)%mode(imode)%mode
        end if
        omega = ant%diag(idiag)%ch(ich)%f_ECE * 2.d0 * pi
        rad%diag(idiag)%ch(ich)%mode(imode)%s_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%R_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%z_res = 0.d0
        rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res = 0.d0
        do ir = 1, N_ray
          ifreq = 1
          grid_size = 1
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res = 0.d0
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res = 0.d0
          ray_segment(1)%x_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec ! get launching position
          ray_segment(1)%N_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%N_vec ! get launching angles
          call find_first_point_in_plasma(plasma_params, omega, mode, ray_segment, last_N, wall_hits, LOS_end)
          if(LOS_end) wall_hits = 2
          N_init = last_N
          !print*, "vessel", last_N
          !ray_segment(1) = ray_segment(last_N - 1) ! redo a bit of the ray to make sure the ray segment also covers the point a
          if(debug_level > 0 .and. output_level) then
            print*, "First point in plasma",ray_segment(last_N)%R_vec
          end if
          if(last_N  + 1 <= max_points_svec .and. .not. LOS_end) then
            call make_ray_segment(20.d0, plasma_params, omega, mode, ray_segment, last_N, wall_hits, N_init)
          else if(.not. LOS_end) then
            print*,"Ray reached maximum length when searching for first point in vessel"
            print*, "Most likely something is very wrong the launching geometry of the diagnostic"
            print*, "Current diagnostic", ant%diag(idiag)%diag_name
            print*, "position and launch vector in Carthesian coordinates", ray_segment(1)%x_vec, &
              ray_segment(1)%N_vec
            stop "Error when finding first point in plasma in mod_raytrace.f90"
          end if
          if(last_N >= max_points_svec) then
            print*, "WARNING a ray did not reach the plasma wall"
            print*, "From here on  output is only for debugging purposes"
            wall_hits = 2
          end if
          if( wall_hits < 2) then
            debug_level = 1
            call find_first_point_in_plasma(plasma_params, omega, mode, ray_segment, last_N, wall_hits, LOS_end)
            N_init = last_N
            call make_ray_segment(20.d0, plasma_params, omega, mode, ray_segment, last_N, wall_hits, N_init)
            print*, "Ray in span_svecs did not end at a wall"
            print*, ray_segment(1)%R_vec(1), ray_segment(1)%R_vec(2), ray_segment(1)%R_vec(3)
            print*, ray_segment(last_N)%R_vec(1), ray_segment(last_N)%R_vec(2), ray_segment(last_N)%R_vec(3)
            print*, "Travelled distance", ray_segment(last_N)%s - ray_segment(N_init)%s
            stop "Error with rays in mod_raytrace.f90"
          end if
          !print*, "plasma", last_N
          if(output_level) then
            if(.not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output)) allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(N_ray))
            if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s)) then
            ! Necessary for IDA with multiple calls to span_svecs
               deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop, &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta)
               deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary, &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary, &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary, &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary, &
                  rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD, rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary)
            end if
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%s(max_points_svec), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%x(max_points_svec), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%y(max_points_svec), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%z(max_points_svec), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%H(max_points_svec), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_ray(max_points_svec), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold(max_points_svec), rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%rhop(max_points_svec), &
                   rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%theta(max_points_svec))
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(max_points_svec), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(max_points_svec))
            do N=1, last_N
              X = func_X(plasma_params,omega, ray_segment(N)%n_e, ray_segment(N)%T_e)
              Y = func_Y(plasma_params,omega, ray_segment(N)%omega_c / e0 * mass_e, ray_segment(N)%T_e)
              rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%N_cold(N) = func_N(X, Y, &
                ray_segment(N)%theta, -rad%diag(idiag)%ch(ich)%mode(imode)%mode)
            end do
            if(.not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)) then
            ! This is important in the case we want some extra output for the last ida optimization
              allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
            end if
          end if
          ray_segment(1:last_N) = ray_segment(last_N:1:-1)
          ray_segment(1:last_N)%s = ray_segment(1)%s - ray_segment(1:last_N)%s
          call preprare_svec_spline(ray_segment(N_init:last_N), rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_spline, &
                                    rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), omega / (2.d0 * pi))
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points = max_points_svec
          if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points > largest_svec) &
            largest_svec = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
          if(output_level) then
            if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad)) then
            ! Necessary for IDA with multiple calls to span_svecs
              deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD, &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary)

            end if
            allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                     rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%BPD_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points))
          end if
          if(output_level) then
            if(rad%diag(idiag)%ch(ich)%mode(imode)%mode == 1) then
              print*, "Channel", ich, "ray", ir, "X-mode ", "initialized", "(", &
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, "points in total)"
            else
              print*, "Channel", ich, "ray", ir, "O-mode ", "initialized", "(", &
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points, "points in total)"
            end if
          end if
          ! Central frequency
          if(N_freq == 1) then
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res = &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%s_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res = &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%R_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res = &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%z_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res = &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%rhop_res
          end if
          do ifreq = 2, N_freq ! Copy ray to the other frequencies
            if(output_level .and. .not. allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output)) then
              allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(:) = &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_extra_output(:)
            end if
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_spline = &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%svec_spline
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points = &
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points
            call find_cold_resonance(plasma_params, &
                ant%diag(idiag)%ch(ich)%freq(ifreq) * 2.d0 * pi, &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq), 1, &
                rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points)
            if(output_level) print*, "s_res ", ifreq, "-th frequency", rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%s_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res + &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%s_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res + &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%R_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res + &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%z_res
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res + &
              ant%diag(idiag)%ch(ich)%freq_weight(ifreq) * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%rhop_res
          end do ! ifreq
          rad%diag(idiag)%ch(ich)%mode(imode)%s_res = rad%diag(idiag)%ch(ich)%mode(imode)%s_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_res
          rad%diag(idiag)%ch(ich)%mode(imode)%R_res = rad%diag(idiag)%ch(ich)%mode(imode)%R_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%R_res
          rad%diag(idiag)%ch(ich)%mode(imode)%z_res = rad%diag(idiag)%ch(ich)%mode(imode)%z_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%z_res
          rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res = rad%diag(idiag)%ch(ich)%mode(imode)%rhop_res + &
            ant%diag(idiag)%ch(ich)%ray_launch(ir)%weight * rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%rhop_res
        end do ! ir
        ! stop "check ray_launch"
        if(output_level) then
          if(allocated(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s)) then
            deallocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Te, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD, &
            rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary)
          end if
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Te(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm(max_points_svec))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD(pnts_BPD))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(pnts_BPD))
          allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary(pnts_BPD))
        end if
      end do ! imode
      if(modes == 1) then ! X-mode
      ! set the resonance position of the channel to the resonance position of the X-mode
      ! if multiple mode analysis is switched on this will be overwritten later
      ! The step is neccessary to quickly perform classical ECE analysis within IDA
        rad%diag(idiag)%ch(ich)%s_res = rad%diag(idiag)%ch(ich)%mode(1)%s_res
        rad%diag(idiag)%ch(ich)%R_res = rad%diag(idiag)%ch(ich)%mode(1)%R_res
        rad%diag(idiag)%ch(ich)%z_res = rad%diag(idiag)%ch(ich)%mode(1)%z_res
        rad%diag(idiag)%ch(ich)%rhop_res = rad%diag(idiag)%ch(ich)%mode(1)%rhop_res
      end if
      ich_tot = ich_tot + 1
    end do !ich
    rad%diag(idiag)%ch(:)%eval_ch = .true. ! start with all channels set to true
  end do ! idiag
  !stop "Early ray end?"
  end subroutine create_svec_splines


!!!!!! mod_ecfm_refr_utils
subroutine make_CEC_diag_from_shotfile(diag, rad_diag, wg, z_lens)
! calculate weight for each LOS and theta for each aufp along LOS:
!    ant%ch(ich)%ray(ir)%aufp(:)%theta, ant%ch(ich)%ray(ir)%weight
use mod_ecfm_refr_types,        only: ant_diag_type, rad_diag_type, plasma_params, modes, mode_cnt, &
                                      N_ray, N_freq, CEC_exp, CEC_ed, output_level, max_points_svec, &
                                      output_level, stand_alone
use constants,                  only: pi
use aug_db_routines,            only : read_ece, read_rmd
implicit none
type(ant_diag_type) , intent(inout)         :: diag
type(rad_diag_type) , intent(inout)         :: rad_diag
integer(ikind), dimension(:), allocatable, intent(out)   :: wg
real(rkind),  intent(out)                   :: z_lens
!real(rkind)    :: d                                ! distance of antenna to perpendicular
!                                                   ! (through center of lense)
!real(rkind)    :: l1                               ! distance from torus center to intersection
!                                                   ! of LOS with perpendicular
!real(rkind)    :: l2                               ! distance from intersection of LOS with
!                                                   ! perpendicular to lense (= focal length)
integer(ikind)                    :: ich, imode, ir, ifreq, N_ch, useful_ch, useless_cnt, cur_ifgroup, wg_index,N_R_z
integer(ikind), dimension(:), allocatable :: available, wg_temp, ifgroup
real(rkind), dimension(:), allocatable :: f, df
character(4)   :: flag
! calculate angle between LOS and Btor
diag%expname = trim(CEC_exp)
diag%edition = CEC_ed
if(diag%diag_name == "CEC") then
  call read_ece(shot         = plasma_params%shot,           & ! in
              expnam       = diag%expname,     & ! in
              ed           = diag%edition,     & ! in
              ed_out       = diag%edition,       & ! out
              f            = f,              & ! inout
              df           = df,             & ! inout
              N_ch         = N_ch,        & ! out
              zlens        = z_lens,         & ! out
              ifgroup      = ifgroup,        & ! out
              waveguid     = wg_temp,       & ! inout
              availabl     = available) ! inout
else
  call read_RMD(shot       = plasma_params%shot,           & ! in
              expnam       = diag%expname,     & ! in
              ed           = diag%edition,     & ! in
              ed_out       = diag%edition,       & ! out
              f            = f,              & ! inout
              df           = df,             & ! inout
              N_ch         = N_ch,        & ! out
              zlens        = z_lens,         & ! out
              ifgroup      = ifgroup,        & ! out
              waveguid     = wg_temp,       & ! inout
              availabl     = available) ! inout
end if
useful_ch = 0
do ich = 1, N_ch
  if(available(ich) > 0) useful_ch = useful_ch + 1
end do
diag%N_ch = useful_ch
if(output_level .and. stand_alone) print*, "Found", useful_ch, " useful channels for ", diag%diag_name
allocate(diag%ch(useful_ch), wg(useful_ch))
allocate(rad_diag%ch(useful_ch))
useless_cnt = 0
wg_index = 1
cur_ifgroup = ifgroup(1)
do ich = 1, diag%N_ch
  allocate(diag%ch(ich)%freq(N_freq))
  allocate(diag%ch(ich)%freq_weight(N_freq))
  allocate(rad_diag%ch(ich)%mode(mode_cnt))
  if(output_level) allocate(rad_diag%ch(ich)%mode_extra_output(mode_cnt))
  do imode = 1, mode_cnt
    if((imode == 2 .and. modes == 3) .or. &
        modes == 2) then
      rad_diag%ch(ich)%mode(imode)%mode = -1 ! O-mode
    else
      rad_diag%ch(ich)%mode(imode)%mode = +1 ! X-mode
    end if
    allocate(rad_diag%ch(ich)%mode(imode)%ray(N_ray))
    if(output_level) allocate(rad_diag%ch(ich)%mode(imode)%ray_extra_output(N_ray))
    do ir = 1, N_ray
      allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
      do ifreq = 1, N_freq
          allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
          if(output_level) allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
      end do
    end do
  end do
  allocate(diag%ch(ich)%ray_launch(N_ray))
end do
do ich = 1, useful_ch
  do while(available(ich + useless_cnt) == 0)
    useless_cnt = useless_cnt + 1
  end do
  diag%ch(ich)%f_ECE = f(ich + useless_cnt)
  diag%ch(ich)%df_ECE = df(ich + useless_cnt)
  if(diag%ch(ich)%f_ECE < 1.d9) then
    print*,"Warning an ECE channel has a frequency below 1 GHz"
  end if
  if(cur_ifgroup /= ifgroup(ich + useless_cnt)) then
    wg_index = wg_index + 1
    cur_ifgroup = ifgroup(ich + useless_cnt)
  end if
  wg(ich) = wg_temp(wg_index)
end do !ich
deallocate(wg_temp)
end subroutine make_CEC_diag_from_shotfile

#ifdef AUG
subroutine make_j
  USE mod_ecfm_refr_types, only: ant, plasma_params, rad, data_folder, Spitzer
  use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
  USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
  implicit none
  integer(ikind)                :: ich, ifreq, iint,i, i_last
  real(kind=4), dimension(400) :: rhop_r4, jdotB_r4, j_r4
  real(rkind), dimension(:), allocatable   ::d_rhop, d_j, rhop, j
  real(rkind)                   :: R_mag, z_mag, pf_mag, &
                                   R_sxp, z_sxp, pf_sxp, max_cur_d
  character(4)                  :: exp
  character(3)                  :: diag
  character(100)                :: filename
  integer(kind=4)   :: error, ed, nr
  error = 0
  ed = 0
  exp = "AUGD"
  diag = "EQH"
  nr = 400
  max_cur_d = 1.e9
  call kkeqjpar(error, exp, diag, int(plasma_params%shot,4), ed, real(plasma_params%time,4), &
               11, nr, rhop_r4, jdotB_r4, j_r4)
  print*," Found ", nr, " current points"
  allocate(d_rhop(int(nr,kind=8)),d_j(int(nr,kind=8)))
  d_rhop(:) = real(rhop(1:nr),8)
  d_j(:) = real(j_r4(1:nr),8)
  do i = 1, nr
    d_rhop(i) = sqrt((d_rhop(i) - plasma_params%pf_mag)/(plasma_params%pf_sxp - plasma_params%pf_mag))
  end do
  allocate(rhop(int(nr,kind=8)), j(int(nr,kind=8)))
  do i = 1, nr
    rhop(i) = d_rhop(nr - i + 1)
    if(abs(d_j(nr - i + 1)) < max_cur_d  .and. .not. d_j(nr - i + 1) /= d_j(nr - i + 1)) then
      j(i) = d_j(nr - i + 1)
    else
      j(i) = 0.d0
    end if
  end do
  call make_1d_spline(Spitzer%j_spl, nr, rhop, j)
#ifdef NAG
  call nag_spline_1d_interp(rhop, j, Spitzer%j_nag_spl)
#endif
  deallocate(rhop, j)
end subroutine make_j
#else
subroutine make_j()
implicit none
print*, "Make make_j uses AUG specific routines - TCV replacements required"
call abort()
end subroutine make_j
#endif
#ifdef IDA
subroutine make_CEC_diag_from_shotfile(diag, rad_diag, wg, z_lens)
! calculate weight for each LOS and theta for each aufp along LOS:
!    ant%ch(ich)%ray(ir)%aufp(:)%theta, ant%ch(ich)%ray(ir)%weight
use mod_ecfm_refr_types,        only: ant_diag_type, rad_diag_type, plasma_params, modes, mode_cnt, &
                                      N_ray, N_freq, CEC_exp, CEC_ed, output_level, max_points_svec, &
                                      output_level, stand_alone
use constants,                  only: pi
use aug_db_routines,            only : read_ece, read_rmd
implicit none
type(ant_diag_type) , intent(inout)         :: diag
type(rad_diag_type) , intent(inout)         :: rad_diag
integer(ikind), dimension(:), allocatable, intent(out)   :: wg
real(rkind),  intent(out)                   :: z_lens
!real(rkind)    :: d                                ! distance of antenna to perpendicular
!                                                   ! (through center of lense)
!real(rkind)    :: l1                               ! distance from torus center to intersection
!                                                   ! of LOS with perpendicular
!real(rkind)    :: l2                               ! distance from intersection of LOS with
!                                                   ! perpendicular to lense (= focal length)
integer(ikind)                    :: ich, imode, ir, ifreq, N_ch, useful_ch, useless_cnt, cur_ifgroup, wg_index,N_R_z
integer(ikind), dimension(:), allocatable :: available, wg_temp, ifgroup
real(rkind), dimension(:), allocatable :: f, df
character(4)   :: flag
! calculate angle between LOS and Btor
diag%expname = trim(CEC_exp)
diag%edition = CEC_ed
if(diag%diag_name == "CEC") then
  call read_ece(shot         = plasma_params%shot,           & ! in
              expnam       = diag%expname,     & ! in
              ed           = diag%edition,     & ! in
              ed_out       = diag%edition,       & ! out
              f            = f,              & ! inout
              df           = df,             & ! inout
              N_ch         = N_ch,        & ! out
              zlens        = z_lens,         & ! out
              ifgroup      = ifgroup,        & ! out
              waveguid     = wg_temp,       & ! inout
              availabl     = available) ! inout
else
  call read_RMD(shot       = plasma_params%shot,           & ! in
              expnam       = diag%expname,     & ! in
              ed           = diag%edition,     & ! in
              ed_out       = diag%edition,       & ! out
              f            = f,              & ! inout
              df           = df,             & ! inout
              N_ch         = N_ch,        & ! out
              zlens        = z_lens,         & ! out
              ifgroup      = ifgroup,        & ! out
              waveguid     = wg_temp,       & ! inout
              availabl     = available) ! inout
end if
useful_ch = 0
do ich = 1, N_ch
  if(available(ich) > 0) useful_ch = useful_ch + 1
end do
diag%N_ch = useful_ch
if(output_level .and. stand_alone) print*, "Found", useful_ch, " useful channels for ", diag%diag_name
allocate(diag%ch(useful_ch), wg(useful_ch))
allocate(rad_diag%ch(useful_ch))
useless_cnt = 0
wg_index = 1
cur_ifgroup = ifgroup(1)
do ich = 1, diag%N_ch
  allocate(diag%ch(ich)%freq(N_freq))
  allocate(diag%ch(ich)%freq_weight(N_freq))
  allocate(rad_diag%ch(ich)%mode(mode_cnt))
  if(output_level) allocate(rad_diag%ch(ich)%mode_extra_output(mode_cnt))
  do imode = 1, mode_cnt
    if((imode == 2 .and. modes == 3) .or. &
        modes == 2) then
      rad_diag%ch(ich)%mode(imode)%mode = -1 ! O-mode
    else
      rad_diag%ch(ich)%mode(imode)%mode = +1 ! X-mode
    end if
    allocate(rad_diag%ch(ich)%mode(imode)%ray(N_ray))
    if(output_level) allocate(rad_diag%ch(ich)%mode(imode)%ray_extra_output(N_ray))
    do ir = 1, N_ray
      allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
      do ifreq = 1, N_freq
          allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
          if(output_level) allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec_extra_output(max_points_svec))
      end do
    end do
  end do
  allocate(diag%ch(ich)%ray_launch(N_ray))
end do
do ich = 1, useful_ch
  do while(available(ich + useless_cnt) == 0)
    useless_cnt = useless_cnt + 1
  end do
  diag%ch(ich)%f_ECE = f(ich + useless_cnt)
  diag%ch(ich)%df_ECE = df(ich + useless_cnt)
  if(diag%ch(ich)%f_ECE < 1.d9) then
    print*,"Warning an ECE channel has a frequency below 1 GHz"
  end if
  if(cur_ifgroup /= ifgroup(ich + useless_cnt)) then
    wg_index = wg_index + 1
    cur_ifgroup = ifgroup(ich + useless_cnt)
  end if
  wg(ich) = wg_temp(wg_index)
end do !ich
deallocate(wg_temp)
end subroutine make_CEC_diag_from_shotfile

subroutine make_CEC_diag_from_ida(diag, rad_diag, ida, ece_strut, wg, z_lens)
use mod_ecfm_refr_types,        only : rad_diag_type, ant_diag_type, output_level, working_dir, &
                                       modes, mode_cnt, N_ray, N_freq, max_points_svec, stand_alone
use ece_types, only : ece_type
use ida_types, only : ida_type
use constants,                  only : pi
implicit none
type(ant_diag_type), intent(inout)                       :: diag
type(rad_diag_type), intent(inout)                       :: rad_diag
type(ida_type), intent(in)                               :: ida
type(ece_type), intent(in)                               :: ece_strut
integer(ikind), dimension(:), allocatable, intent(out)   :: wg
real(rkind),  intent(out)                                :: z_lens
character(200)                                           :: cur_filename
integer(ikind)                                           :: ich, imode, ir, ifreq
real(rkind)                                              :: N_abs, temp, temp1, temp2, phi_radial
  mode_cnt = 1
  if(modes == 3) mode_cnt = 2
  diag%expname = ece_strut%exp_ece
  diag%diag_name = ece_strut%diag_ece
  diag%edition = ece_strut%ed_ece
  diag%N_ch = ece_strut%N_ch
  allocate(diag%ch(diag%N_ch))
  allocate(rad_diag%ch(diag%N_ch))
  do ich = 1, diag%N_ch
    allocate(diag%ch(ich)%freq(N_freq))
    allocate(diag%ch(ich)%freq_weight(N_freq))
    allocate(rad_diag%ch(ich)%mode(mode_cnt))
    do imode = 1, mode_cnt
      if((imode == 2 .and. modes == 3) .or. &
          modes == 2) then
        rad_diag%ch(ich)%mode(imode)%mode = -1 ! O-mode
      else
        rad_diag%ch(ich)%mode(imode)%mode = +1 ! X-mode
      end if
      allocate(rad_diag%ch(ich)%mode(imode)%ray(N_ray))
      do ir = 1, N_ray
        allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(N_freq))
        do ifreq = 1, N_freq
          allocate(rad_diag%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(max_points_svec))
        end do
      end do
    end do
    allocate(diag%ch(ich)%ray_launch(N_ray))
    diag%ch(ich)%f_ECE = ece_strut%ch(ich)%freq
    diag%ch(ich)%df_ECE = ece_strut%ch(ich)%dfreq
    if(diag%ch(ich)%f_ECE < 1.d9) then
      print*,"Warning an ECE channel has a frequency below 1 GHz"
    end if
  end do
  allocate(wg(diag%N_ch))
  wg = ece_strut%ch(:)%waveguide
  z_lens = ece_strut%z_lens
end subroutine make_CEC_diag_from_ida



# else

subroutine make_CEC_diag_from_shotfile(diag, rad_diag, wg, z_lens)
! calculate weight for each LOS and theta for each aufp along LOS:
!    ant%ch(ich)%ray(ir)%aufp(:)%theta, ant%ch(ich)%ray(ir)%weight
use mod_ecfm_refr_types,        only: ant_diag_type, rad_diag_type, plasma_params, N_ray, N_freq
implicit none
type(ant_diag_type) , intent(inout)      :: diag
type(rad_diag_type) , intent(inout)      :: rad_diag
! To suppress warnings on out variables not being set inout instead of out
integer(ikind), dimension(:),  intent(inout)  :: wg
real(rkind), intent(inout)               :: z_lens
  print*, "The subroutine make_CEC_diag is AUG specific and must not be called on non-AUG systems"
  call abort()
end subroutine make_CEC_diag_from_shotfile

subroutine make_CEC_diag_from_ida(diag, rad_diag, ida, ece_strut, wg, z_lens)
use mod_ecfm_refr_types,        only : plasma_params_type, ant_diag_type, rad_diag_type
use ece_types, only : ece_type
use ida_types, only : ida_type
use constants,                  only : pi
implicit none
type(ida_type), intent(in)               :: ida
type(ece_type), intent(in)               :: ece_strut
type(ant_diag_type) , intent(inout)      :: diag
type(rad_diag_type) , intent(inout)      :: rad_diag
! To suppress warnings on out variables not being set inout instead of out
integer(ikind), dimension(:),  intent(inout)  :: wg
real(rkind) , intent(inout)              :: z_lens
  print*, "The subroutine prepare_ECE_diag_data_from_ida is IDA specific and must not be called on non-AUG systems"
  call abort()
end subroutine make_CEC_diag_from_ida

subroutine make_CEC_diag_geometry(diag, wg, z_lens)
use mod_ecfm_refr_types,        only : plasma_params_type, ant_diag_type
use ece_types, only : ece_type
use ida_types, only : ida_type
use constants,                  only : pi
implicit none
type(ant_diag_type), intent(in)           :: diag
integer(ikind), dimension(:),  intent(in)  :: wg
real(rkind), intent(in)               :: z_lens
  print*, "The subroutine make_CEC_diag_geometry is AUG specific and must not be called on non-AUG systems"
  call abort()
end subroutine make_CEC_diag_geometry
#endif

!mod_ecfm_refr_rad_transp
  use mod_ecfm_refr_types,         only: rad_diag_ch_mode_ray_freq_svec_type, &
                                         rad_diag_ch_mode_ray_freq_svec_extra_output_type
 real(rkind), dimension(116) :: work_lsode
!    real(rkind), dimension(112) :: work_lsode
  integer(ikind), dimension(20) :: iwork_lsode
  integer(ikind)                :: j_LSODE, idiag_LSODE, ich_LSODE, imode_LSODE, iray_LSODE, ifreq_LSODE
  real(rkind)                   :: glob_ab, glob_ab_secondary
  type(rad_diag_ch_mode_ray_freq_svec_type) :: glob_svec
  type(rad_diag_ch_mode_ray_freq_svec_extra_output_type) :: glob_svec_extra_output
subroutine evaluate_em_ab_LSODE(svec, svec_extra_output, s, ds, omega, mode, model, eval_pol_coeff, x_launch, em, ab, pol_coeff)
! Wrapper function for em, ab for the integration of the radiation tranposrt equation with DLSODE
! Since a lot of if clasus are required here this routines cleans up the code below
use mod_ecfm_refr_types,          only: rad_diag_ch_mode_ray_freq_svec_type, rad_diag_ch_mode_ray_freq_svec_extra_output_type, &
                                      output_level, data_folder, Ich_name, &
                                      dstf, dstf_comp, rad, ffp, mode_cnt
use mod_ecfm_refr_em_Hu,                  only: calculate_em, simple_in_cutoff
use constants,                    only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,         only: abs_Albajar, abs_Albajar_fast, abs_Al_tor_abs, func_N_cold, func_rel_N
implicit none
type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout) :: svec
type(rad_diag_ch_mode_ray_freq_svec_extra_output_type), intent(out) :: svec_extra_output
real(rkind),                intent(in)    :: s, omega, ds
integer(ikind),             intent(in)    :: mode
character(*),               intent(in)    :: model
logical,                    intent(in)    :: eval_pol_coeff
real(rkind), dimension(:), intent(in)     :: x_launch
real(rkind),                intent(out)   :: em, ab
real(rkind),                intent(out)   :: pol_coeff
real(rkind)                               :: ab_secondary, dummy
!call benchmark_abs_and_N()
!stop "benchmarking"
em = 0.d0
ab = 0.d0
if(eval_pol_coeff) then
  pol_coeff = 0.d0
end if
! Interpolated svec
if(svec%Te > 1.d0 .and. svec%ne > 1.d15 .and. svec%ne < 1.d21) then
   ! Negative values are possible in the chi^2 calcuation since we use the nag spline in that instance, which may lead to ringing
  svec%N_cold = func_N_cold(omega, svec, mode)! note the convention for mode
  if(trim(model) == "primary") then
    if(output_level .or. (dstf /= "relamax" .and. dstf /= "Hu")) then
      if(output_level) svec_extra_output%N_cor = func_rel_N(omega, svec, mode)
      if(dstf_comp /= "TO") then
       if(eval_pol_coeff .and. mode_cnt == 2) then
          ! Take last point within the separatrix
          call abs_Albajar(svec, omega, mode, ds, ab, em, &
                           pol_coeff = pol_coeff, x_launch = x_launch)
       else
          call abs_Albajar(svec, omega, mode, ds, ab, em)
          !call abs_Albajar_fast(svec, omega, mode, ds2, ab)
          !em = ab * svec%Ibb
        end if
      else
        call calculate_em(svec, omega, em, dummy, ab)                          ! out    [W m^-3 sr^-1 Hz^-1]
      end if
    else
      if(dstf == "relamax") then
        call abs_Albajar_fast(svec, omega, mode, ds, ab)
        em = ab * svec%Ibb
      else if(dstf == "Hu") then
        call calculate_em(svec, omega, em, dummy, ab)                          ! out    [W m^-3 sr^-1 Hz^-1]
        ab =  em / svec%Ibb
      else
        print*, "dstf flag not set correctly - allowed is relamax or Hu"
        print*, "dstf flag is", dstf
        stop "mod_ecfm_rad_int calculate_Trad"
      end if
    end if
  else if(trim(model) == "secondary") then
    if ((dstf_comp == "Al" .or. dstf_comp == "Th") .and. dstf == "relamax") then
      if(mode /= -1) then
        call calculate_em(svec, omega, em, dummy, ab)
      end if
    else if( dstf_comp == "TB" .or. dstf_comp == "TO") then
      if(svec_extra_output%N_warm <= 0 .or. &
        svec_extra_output%N_warm /= svec_extra_output%N_warm) then
        svec_extra_output%N_warm = svec%N_cold ! do not reuse last, but start with cold
      end if
      if(eval_pol_coeff .and. mode_cnt == 2) then
        ! last point within the separatrix
        ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm, &
                    pol_coeff_secondary = pol_coeff, x_launch = x_launch)
      else
        ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm)
      end if
      call abs_Albajar(svec, omega, mode, ds, ab, em, &
                           pol_coeff = pol_coeff, x_launch = x_launch)
      ab = ab_secondary
      em = ab * svec%Ibb
    else if( dstf_comp == "O1") then
      svec_extra_output%N_warm = svec%N_cold
      ab = abs_Al_tor_abs(svec, omega, mode, svec_extra_output%N_warm)
      em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            svec%Te / c0**2
    else
      em = ab * (omega / ( 2.d0 * pi))**2 * e0 * &
            svec%Te / c0**2
    end if
!    if(ab_secondary > 1.d0) then
!      print*, em_secondary , ab_secondary
!    end if
    svec_extra_output%em_secondary = em
    svec_extra_output%ab_secondary = ab
  else
    print*,"Got model", model
    print*,"Expected either 'primary' or 'secondary'"
    stop "Bad arguments in evaluate_em_ab_single_LSODE in mod_ecfm_refr_rad_stranp.f90"
  end if
else
  return
end if
end subroutine evaluate_em_ab_LSODE

subroutine rad_transp(m, s, Trad, dIds)
! RHS of radiation transport equation to be integrated by DLSODE
  use f90_kind
  USE mod_ecfm_refr_types, only : ant, rad, output_level
  use mod_ecfm_refr_interpol,       only: spline_1d
  use constants,                    only: pi, e0, mass_e, eps0, c0
  implicit none
  integer          , intent(in)           :: m
  real(rkind),  intent(in)                :: s
  real(rkind),  dimension(m), intent(in)  :: Trad
  real(rkind),  dimension(m), intent(inout) :: dIds
  real(rkind)                             :: ds, omega, ab, em
  logical                                 :: eval_pol_coeff

  if(s > rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_max) then
    print*, "s in rad_transp >  s_max",s, rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_max
    call abort()
  end if
  glob_svec%s = s
  omega = ant%diag(idiag_LSODE)%ch(ich_LSODE)%freq(ifreq_LSODE) * 2.0 * pi
  if(j_LSODE > 1) then
    ds = s - rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec(j_LSODE - 1)%s
  else
    ds = s - rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%s_min
  end if
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%x, s, &
                  glob_svec%x_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%y, s, &
                 glob_svec%x_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%z, s, &
                 glob_svec%x_vec(3))
  glob_svec%R = sqrt(glob_svec%x_vec(1)**2 + glob_svec%x_vec(2)**2)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%rhop, s, &
                 glob_svec%rhop)
  if(glob_svec%rhop > 0.98 .and. glob_svec%rhop < 1.d0 .and. &
        .not. rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%use_external_pol_coeff) then
     eval_pol_coeff = .true.
  end if
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%theta, s, &
                 glob_svec%theta)
   glob_svec%cos_theta = cos(glob_svec%theta)
   glob_svec%sin_theta = sin(glob_svec%theta)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%freq_2X, s, &
                 glob_svec%freq_2X)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Te, s, &
                 glob_svec%Te)
  glob_svec%Ibb = glob_svec%Te * (omega / (2.d0 * pi))**2 * e0 / c0**2
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%ne, s, &
                 glob_svec%ne)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%v_g_perp, s, &
                 glob_svec%v_g_perp)
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Nx, s, &
                 glob_svec%N_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Ny, s,  &
                 glob_svec%N_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Nz, s,  &
                 glob_svec%N_vec(3))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Bx, s, &
                 glob_svec%B_vec(1))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%By, s,  &
                 glob_svec%B_vec(2))
  call spline_1d(rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%svec_spline%Bz, s,  &
                 glob_svec%B_vec(3))
  call evaluate_em_ab_LSODE( glob_svec, glob_svec_extra_output, &
                             s, ds, omega, rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%mode, "primary", eval_pol_coeff, &
                             ant%diag(idiag_LSODE)%ch(ich_LSODE)%ray_launch(iray_LSODE)%x_vec, em, ab, &
                             rad%diag(idiag_LSODE)%ch(ich_LSODE)%mode(imode_LSODE)%ray(iray_LSODE)%freq(ifreq_LSODE)%pol_coeff)
  glob_ab = ab
  dIds(1) = em  * c0**2 / (omega**2 * e0) * 4.d0 * pi**2 - ab * Trad(1)
end subroutine rad_transp

subroutine Jac(nq,path, y_vec,ml,mu,dfdy,npowpd)
! Empty routine - gradients evaluated numerically by DLSODE
  use f90_kind
  implicit none
  integer,   intent(in)   :: nq,ml,mu,npowpd
  real(rkind), intent(in) :: path
  real(rkind), intent(inout) :: y_vec(:), dfdy(:,:)
   ! empty subroutine for the ODE-solver - if implicit scheme is used use finite differences
  return
end subroutine Jac

subroutine calculate_Trad_LSODE(rad_ray, idiag, ich, imode, ir, ifreq, mode, Trad, Trad_secondary, tau, tau_secondary, tau_array, tau_secondary_array, error, debug)
! In its current state this method does not work. Neither DLSODE mode 4 nor mode 5 produce reasonable results. Even if the resonance is marked as a
! difficult point, the predictor, corrector method of DLSODE is not capable of solving the radiation transport equation more efficiently than the RK4 method below.
! DO NOT USE THIS ROUTINE!!!!
use mod_ecfm_refr_types,        only: ant, rad_diag_ch_mode_ray_type, output_level, data_folder, Ich_name, straight, max_points_svec, &
                                      tau_array, tau_secondary_array, dstf, dstf_comp, rad, ffp, mode_cnt, tau_thick
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,           only: func_N_cold, func_rel_N
use mod_ecfm_refr_interpol,       only: spline_1d
implicit none
type(rad_diag_ch_mode_ray_type), intent(inout) :: rad_ray
integer(ikind),             intent(in)    :: idiag, ich, imode, ir, ifreq
integer(ikind),             intent(in)    :: mode
real(rkind),                intent(out)   :: Trad, Trad_secondary
real(rkind),                intent(out)   :: tau, tau_secondary
real(rkind), dimension(:), intent(inout)  :: tau_array, tau_secondary_array
integer(ikind),             intent(out)   :: error
logical,    intent(in), optional          :: debug
real(rkind)    :: s, ds
real(rkind), dimension(:), allocatable    :: Trad_vec, dIds
integer(ikind)                            :: nTrad, istate, n_s_important, i_important
character(120)               :: cur_filename
character(20)                :: ich_str
real(rkind)    :: h, h_large, h_small, freq_2X_step, freq, s_in, s_switch
real(rkind), dimension(20) :: s_important
logical         :: first, eval_pol_coeff, debug_internal
error = 0
idiag_LSODE = idiag; ich_LSODE = ich; imode_LSODE = imode; iray_LSODE = ir; ifreq_LSODE = ifreq
! initialization
Trad     =  0.d0
Trad_secondary =  0.d0
tau      =  0.d0
tau_secondary  =  0.d0
debug_internal = .false.
rad_ray%freq(ifreq)%svec(:)%rhop = -1.d0 ! to mark all entries in svec that are not used for the radiation transport
if(present(debug)) debug_internal = debug
if(.not. rad_ray%freq(ifreq)%use_external_pol_coeff) then
  rad_ray%freq(ifreq)%pol_coeff = 1.0
  rad_ray%freq(ifreq)%pol_coeff_secondary = 1.0
end if
first = .false.
tau_array(:) = 0.d0
if(output_level) then
  nTrad = 2
  tau_secondary_array(:) = 0.d0
else
  nTrad = 1
end if
allocate(Trad_vec(nTrad), dIds(nTrad))
Trad_vec(:) = 0.d0
freq =  ant%diag(idiag)%ch(ich)%freq(ifreq)
!-----------------------------------------------------------------
j_LSODE = 1
s = rad_ray%freq(ifreq)%svec_spline%s_min
if(dstf == "Re") ffp%restart_spline = .True.
if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE)%N_warm = 0.d0 ! force cold plasma eval
call rad_transp(nTrad, s, Trad_vec, dIds)
rad_ray%freq(ifreq)%svec(j_LSODE) = glob_svec
if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE) = glob_svec_extra_output
h_large = 5.d-2
h_small = 1.d-2
istate = 1
n_s_important = 3
s_important(1) = rad_ray%freq(ifreq)%svec_spline%s_center
s_important(2) = rad_ray%freq(ifreq_LSODE)%s_res
s_important(3) = rad_ray%freq(ifreq)%svec_spline%s_max
!print*, "s start", s_important(1:n_s_important -1)
!print*, "s end", s_important(2:n_s_important)
do while(any(s_important(1:n_s_important -1) > s_important(2:n_s_important)))
  do i_important = 1, n_s_important - 1
    if(s_important(i_important) > s_important(i_important + 1)) then
      s_switch = s_important(i_important)
      s_important(i_important) = s_important(i_important + 1)
      s_important(i_important + 1) = s_switch
    end if
  end do
end do
work_lsode(1) = s_important(1)
i_important = 1
do while(s < rad_ray%freq(ifreq_LSODE)%svec_spline%s_max)
  if(j_LSODE ==  max_points_svec) then
    print*, "Need more maximum points in svec"
    print*, "Please increase max_points_svec and recompile"
    call abort()
  end if
  if(abs(rad_ray%freq(ifreq)%svec(j_LSODE)%freq_2X / freq - 1.d0) <  0.02) then
    h = h_small
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
  else
    h = h_large
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
    call spline_1d(rad_ray%freq(ifreq)%svec_spline%freq_2X, s, freq_2X_step)
    if(abs(freq_2X_step / freq - 1.d0) <  0.02) h = h_small
    if(s + h > rad_ray%freq(ifreq)%svec_spline%s_max) h = rad_ray%freq(ifreq)%svec_spline%s_max - s
  end if
  if(h == 0.d0) then
    print*, "h in rad transp last step was zero", s
    call abort
  end if
  s_in = s
  !print*, "s + h, s_important(i_important), i_important", s + h, s_important(i_important), i_important
  do while(s + h >= s_important(i_important) .and. i_important < n_s_important)
  ! Need do loop here in case s < s(center) < s(res) < s + h
    i_important = i_important + 1
    work_lsode(1) = s_important(i_important)
    !print*, "New s_crit, s_max", s_important(i_important), rad_ray%freq(ifreq)%svec_spline%s_max
  end do
  call dlsode(rad_transp, (/nTrad/) , Trad_vec, s, s + h,       &
                             1, (/1d-7/), (/1.d-7/), 5, istate, 1, &
                             work_lsode, size(work_lsode), iwork_lsode, size(iwork_lsode), Jac,10)
  if (istate /= 2) then
    print*, "h" ,h, "istate", istate, "Trad", Trad_vec
    print*, "LSODE failed to solve the radiation transport equation"
    call abort()
  end if
  ds = s - s_in
  !print*, "s, s_in, ds", s, s_in, ds
  if(glob_svec%s - rad_ray%freq(ifreq)%svec(j_LSODE)%s <= 0) cycle
  j_LSODE = j_LSODE + 1
  if(j_LSODE > 1) then
    tau_array(j_LSODE) = tau_array(j_LSODE) + glob_ab * ds
  else
    tau_array(j_LSODE) = glob_ab  * ds
  end if
  if(output_level) then
    if(j_LSODE > 1) then
      tau_secondary_array(j_LSODE) = tau_secondary_array(j_LSODE - 1)  + glob_ab_secondary  * ds
    else
      tau_secondary_array(j_LSODE) = glob_ab_secondary  * ds
    end if
  end if
  rad_ray%freq(ifreq)%svec(j_LSODE) = glob_svec
  if(output_level) rad_ray%freq(ifreq)%svec_extra_output(j_LSODE) = glob_svec_extra_output
  !print*, "j_LSODE, s", j_LSODE, rad_ray%freq(ifreq)%svec(j_LSODE)%s
  if(Trad_vec(1) /= Trad_vec(1) .or. Trad_vec(1) > max(glob_svec%Te * 10.0, 1.d6) &
                                            .or. tau_array(j_LSODE) < 0.d0) then
    ! Trad >> 1 MeV => much larer than saturation value of most ECEs
    print*, "Trad", Trad_vec(1)
    print*, "I_bb", glob_svec%Te / c0**2 * &
                    (freq**2 * e0)
    print*, "ab", glob_ab
    print*, "tau", tau_array(j_LSODE)
    print*, "optical depth of ds", glob_ab * ds
    print*, "length of ds", ds
    print*, "Te", glob_svec%Te
    print*, "ne", glob_svec%ne
    print*, "theta", glob_svec%theta / pi * 180
    print*, "freq_2X", glob_svec%freq_2X
    print*, "N_abs",  glob_svec%N_cold
    call abort
  end if
  if(output_level) then
    if(Trad_vec(2) /= Trad_vec(2) .or. Trad_vec(2) > max(glob_svec%Te * 10.0, 1.d6) &
                                            .or. tau_secondary_array(j_LSODE) < 0.d0) then
        ! Trad >> 1 MeV => much larer than saturation value of most ECEs
        print*, "Trad secondary", Trad_vec(2)
        print*, "I_bb", glob_svec%Te / c0**2 * &
                        (freq**2 * e0)
        print*, "ab secondary", glob_ab_secondary
        print*, "tau", tau_secondary_array(j_LSODE)
        print*, "optical depth of ds", glob_ab_secondary * ds
        print*, "length of ds", ds
        print*, "Te", glob_svec%Te
        print*, "ne", glob_svec%ne
        print*, "theta", glob_svec%theta / pi * 180
        print*, "freq_2X", glob_svec%freq_2X
        print*, "N_abs",  glob_svec%N_cold
        call abort
      end if
  end if
end do
if(j_LSODE == 1) then
  print*, "LSODE finished after just one step"
  print*, "s", rad_ray%freq(ifreq)%svec_spline%s_min, s, rad_ray%freq(ifreq)%svec_spline%s_max
  call abort()
end if
j_LSODE = j_LSODE - 1
rad_ray%freq(ifreq)%total_LOS_points = j_LSODE
tau = tau_array(j_LSODE)
Trad = Trad_vec(1)
if((ant%diag(idiag)%diag_name == "CTC" .or. ant%diag(idiag)%diag_name == "CTA" .or. &
   ant%diag(idiag)%diag_name == "IEC") .and. .not.  &
        rad_ray%freq(ifreq)%use_external_pol_coeff) then
  if(mode == - 1) then
    rad_ray%freq(ifreq)%pol_coeff = 0.d0
    rad_ray%freq(ifreq)%pol_coeff_secondary = 0.d0
  else
    rad_ray%freq(ifreq)%pol_coeff = 1.d0
    rad_ray%freq(ifreq)%pol_coeff_secondary = 1.d0
  end if
end if
if(output_level) then
  tau_secondary = tau_secondary_array(j_LSODE)
  Trad_secondary = Trad_vec(2)
  tau_secondary = tau_secondary_array(rad_ray%freq(ifreq)%total_LOS_points)
  rad_ray%freq(ifreq)%svec_extra_output(1:rad_ray%freq(ifreq)%total_LOS_points)%T = &
    exp(-(tau - tau_array(1:rad_ray%freq(ifreq)%total_LOS_points)))
  rad_ray%freq(ifreq)%svec_extra_output(1:rad_ray%freq(ifreq)%total_LOS_points)%T_secondary = &
    exp(-(tau_secondary - tau_secondary_array(1:rad_ray%freq(ifreq)%total_LOS_points)))
end if
end subroutine calculate_Trad_LSODE


! entirety of rad_transp_int

module mod_ecfm_refr_rad_transp_int
! THIS MODULE CONTAINS FAILED ATTEMPTS TO SOLVE THE RADIATION TRANSPORT EQUATION EFFICIENTLY.
! ALL ROUTINES DO NOT WORK. THIS MODULE IS CURRENTLY NOT COMPILED WITH THE CURRENT MAKE FILE.
! THIS FILE EXISTS SOLELY FOR RECORD KEEPING PURPOSES. DO NOT USE!
  use f90_kind

  implicit none

  !*******************************************************************************
subroutine bounded_quadratic_interpolation(x0, x1, x2, y0, y1, y2, x, y)
! Quadratic interpolation for monotonically decreasing functions
! y is never smaller than y2
use f90_kind
implicit none
real(rkind) , intent(in)                :: x0, x1, x2, y0, y1, y2, x
real(rkind), intent(out)                :: y
real(rkind)                             :: l0, l1 ,l2 ! lagrance coefficients
! Source wikipedia for n = 2
l0 = (x - x1) * (x - x2) / ((x0 -x1) * (x0 - x2))
l1 = (x - x0) * (x - x2) / ((x1 -x0) * (x1 - x2))
l2 = (x - x0) * (x - x1) / ((x2 -x0) * (x2 - x1))
y = y0 * l0 + y1 * l1 + y2 * l2
y = 0.5d0 * (1.d0 + sign(1.d0, y - y2)) * y + 0.5d0 * (1 + sign( 1.d0 , y2 - y)) * y2 ! should be identical to if(y < y2) y = y2
 ! negative values become zero
!if( y > 1.d0) y = 1.d0
end subroutine bounded_quadratic_interpolation

subroutine integrate_Trad(rad_ray, idiag, ich, imode, ir,  ifreq, mode, Trad, Trad_secondary, tau, tau_secondary)
! WARNING - this routine is known to inaccurate results (numerical instability) DO NOT USE!
! The main problem lies in the fact that the emission is peaked and we cannot be sure where the peak lies without (expensive) calculations
! While the transmittance can be calculated easily and with high resolution, the second integral becomes very expensive, since a small step size
! is neccessary to catch all peaks in the emissivity
! Solves the integral represntation of radiation transport equation using Gaussian quadrature.
! If the integral formulation is used computational power can be saved as only regions with significant,
! observed radiation are considered. I.e. regions, where the emitted radiation is fully reabsorbed are not considered.
! The Integral is solved on many subsegments. These subsegments given by the ray obtained with ray tracing.
use mod_ecfm_refr_types,        only: ant, rad_diag_ch_mode_ray_type, output_level, data_folder, Ich_name, ray_out_folder,&
                                      Ich_name, tau_array,tau_secondary_array, dstf, dstf_comp, plasma_params
use mod_ecfm_refr_em_Hu,                  only: calculate_em, simple_in_cutoff
use constants,                  only: pi, e0, mass_e, eps0, c0
use mod_ecfm_refr_abs_Al,         only: abs_Albajar, abs_Al_tor_abs, func_N_cold, func_rel_N
use mod_ecfm_refr_utils,    only: interpol_LOS
use mod_raytrace,               only: make_ray_segment, find_first_point_in_plasma, prepare_svec_segment_gauss
use interpolation_routines,     only: linear_interpolation
implicit none
type(rad_diag_ch_mode_ray_type), intent(inout) :: rad_ray
integer(ikind),             intent(in)    :: idiag, ich, imode, ir, ifreq
integer(ikind),             intent(in)    :: mode
real(rkind),                intent(out)   :: Trad, Trad_secondary
real(rkind),                intent(out)   :: tau, tau_secondary
real(rkind)                               :: a,b, distance, tau_last, tau_secondary_last, tau_int
character(120)               :: cur_filename
real(rkind), dimension(plasma_params%int_step_cnt) :: em, ab, em_secondary, ab_secondary, tau_step!, dN_abs_2_ds, N_abs_2
real(rkind), dimension(4)                          :: int_dist
integer(ikind)                                     :: grid_size, total_step_cnt
integer(ikind), dimension(4)                       :: grid_size_profiler
real(rkind)    ::  omega, Y_closest_res, dummy, ab_last, N_cor, N_cold, N_gray, tau_secondary_int,T_int, T_secondary_int
integer(ikind) :: i, j, k, last_N, ray_segment_cnt, last_j, int_mode_trans, N
logical         :: transport, wall_hit,  extend_ray, extra_output, cycle_int
Trad     =  0.d0
Trad_secondary =  0.d0
tau      =  0.d0
tau_secondary  =  0.d0
distance = 1.d-1 ! Controlls length of ray segments
extra_output = .false.
int_mode_trans = 2 ! use trapezoid for integration of the transmissivity
j = 1
a = 0.d0 ! s at launching position
omega = ant%diag(idiag)%ch(ich)%freq(ifreq) * 2.d0 * pi
plasma_params%ray_segment(1)%x_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%x_vec ! get launching position
plasma_params%ray_segment(1)%N_vec = ant%diag(idiag)%ch(ich)%ray_launch(ir)%N_vec ! get launching angles
ray_segment_cnt = 0
if(extra_output) print*, "first point before plasma search", plasma_params%ray_segment(1)%s
if(output_level) then
    write(cur_filename, "(A53A5I3.3A4)") ray_out_folder,"Raych", ich,".dat"
    open(74, file=cur_filename)
    write(cur_filename, "(A49A5A4I3.3A4)") data_folder,Ich_name,"/Nch", ich,".dat"
  else
    write(cur_filename, "(A44A5A4I3.3A4)") data_folder,Ich_name,"/Nch", ich,".dat"
  end if
!  !write(cur_filename_3X, "(A64A5A10I3.3A4)") data_folder,Ich_name,"/Irhopch3X", ich,".dat"
  open(75, file=cur_filename)
call find_first_point_in_plasma(plasma_params, omega, plasma_params%ray_segment, last_N) ! find first point in plasma (rhop defined)
if(output_level) then
  do N=1, last_N
    write(74,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
      plasma_params%ray_segment(N)%s, " ", plasma_params%ray_segment(N)%R_vec(1), " ", plasma_params%ray_segment(N)%R_vec(2), " ",&
      plasma_params%ray_segment(N)%R_vec(3)," ", plasma_params%ray_segment(N)%Hamil
  end do
end if
! very coarse grid
! -- debug -- !
  if(extra_output) print*, "last point after plasma search", plasma_params%ray_segment(last_N)%s
  if(extra_output) print*, "last_N", last_N
  ray_segment_cnt = ray_segment_cnt + 1
  if(extra_output) print*,"ray segments:", ray_segment_cnt
! Make the svec also for regions of the wall - neccessary for cold res. positions inside the wall
grid_size = 4
rad_ray%freq(ifreq)%s_res = -1.d0
if(last_N > 1) then
  call prepare_svec_segment_gauss(plasma_params, omega, rad_ray, ifreq, plasma_params%ray_segment, last_N, 1, &
    plasma_params%int_step_cnt, (/2.d0, 2.d0, 2.d0, 2.d0/), grid_size, 0.d0, extend_ray, a,b, .false.)
  ! dist = 2 m since we want to cover the entire ray in one go
  if(extra_output) print*, "a, b, first, last point of ray", a, b, plasma_params%ray_segment(1)%s, plasma_params%ray_segment(last_N)%s
  ! a,b define the edges of the current integral
  ! note that rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s is not the same as a and b,
  ! because the outer abszissae of the gaussian quadrature are not eaxctly -1, 1
  j = plasma_params%int_step_cnt + 1
  ! move to first point of next segment
else
  stop "the LOS of all ECE diagnostics has to start outside of the plasma vessel"
end if
!print*, a, b, plasma_params%ray_segment(1)%s,plasma_params%ray_segment(last_N)%s
transport = .true.
!print*, "First point vessel", plasma_params%ray_segment(last_N)%x_vec
if(.not. rad_ray%freq(ifreq)%use_external_pol_coeff) then
  rad_ray%freq(ifreq)%pol_coeff = 1.0
  rad_ray%freq(ifreq)%pol_coeff_secondary = 1.0
end if
extend_ray = .true. ! if true the radiation transport along the entire ray has been considered and the next ray segment can be started
plasma_params%h = 5.d-4 !  step width for the ray tracer - WARNING: SMALL values of h cause large numerical inaccuracy (many steps with small error -> large error)
rad_ray%freq(ifreq)%svec(1:plasma_params%int_step_cnt + 1)%T = 1.d0 ! Transmittance = 1 if no absorption
if(output_level)  rad_ray%freq(ifreq)%svec(1:plasma_params%int_step_cnt + 1)%T_secondary = 1.d0  ! Transmittance = 1 if no absorption
int_dist(4) = 1.d-1 ! very large step size if far away from any resonance (kinetic or cold)
int_dist(3) = 1.d-2 ! moderate stepsize if ab_last is small and we are near a purely kinetic resonance
int_dist(2) = 2.5d-3 ! smaller step size if near a cold resonance
int_dist(1) = 1.0d-3 ! smallest step size only if (b-a) * ab > plasma_params%ab_switch
grid_size_profiler(:) = 0 ! counts the steps taken for each step size
total_step_cnt = 0 ! counts the integrations
ab_last = 0.d0 ! redo current step if step size is judged as too small
tau_int = 0.d0
tau_secondary_int = 0.d0
j = plasma_params%int_step_cnt + 1
do while(transport)
  a = b ! the last end point becomes the starting point
  if(extend_ray) then ! make next segment
    plasma_params%ray_segment(1) = plasma_params%ray_segment(last_N)
    ! set last point as first point for next ray  - seamless svec
    if(extra_output) print*, "first point before", plasma_params%ray_segment(1)%s
    if(b /= plasma_params%ray_segment(1)%s) then
      print*, "how did we get here"
      print*, "b", b
      print*, "b - int_dist",b - int_dist
      print*, plasma_params%ray_segment(1)%s
      stop "extend_ray true although svec not done with this segment yet"
    end if
    plasma_params%ray_segment(1) = plasma_params%ray_segment(last_N)
    call make_ray_segment(distance, plasma_params, omega, plasma_params%ray_segment, last_N, wall_hit) ! prepare the next segment
    if(output_level) then
      do N=1, last_N
        write(74,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
          plasma_params%ray_segment(N)%s, " ", plasma_params%ray_segment(N)%R_vec(1), " ", plasma_params%ray_segment(N)%R_vec(2), " ",&
          plasma_params%ray_segment(N)%R_vec(3)," ", plasma_params%ray_segment(N)%Hamil
      end do
    end if
    if(extra_output) print*, "first point after", plasma_params%ray_segment(1)%s
    if(extra_output) print*, "last point after", plasma_params%ray_segment(last_N)%s
    ray_segment_cnt = ray_segment_cnt + 1
    if(extra_output) print*,"ray segments:", ray_segment_cnt
    !if(plasma_params%ray_segment(last_N)%s - plasma_params%ray_segment(1)%s > distance * 3) then
    !  print*, "Overshoot ray by: ", plasma_params%ray_segment(last_N)%s - plasma_params%ray_segment(1)%s, distance
    !end if
    if(last_N == 1 .and. wall_hit) then
      if(extra_output) print*, "ran into wall - finishing up"
      exit
    end if
    if(last_N == 1) then
      print*, "ray consists of just one dot"
      stop "this should not happen"
    end if
  end if
  !print*, "after ray",  plasma_params%ray_segment(1)%s
  !print*, a, int_dist, a + int_dist
  !print*,"distance covered, target distance:", plasma_params%ray_segment(last_N)%s - plasma_params%ray_segment(1)%s, distance
  call prepare_svec_segment_gauss(plasma_params, omega, rad_ray, ifreq, plasma_params%ray_segment, &
        last_N, j, j + plasma_params%int_step_cnt - 1, int_dist, grid_size, ab_last, extend_ray, a, b, cycle_int)
  if(extra_output) then
    print*, "j at start of segment", j
    print*, "grid_size", grid_size
    print*, "s first, s_last", rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s
    print*, "Y first, Y_last", rad_ray%freq(ifreq)%svec(j)%freq_2X/omega * Pi, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%freq_2X/omega * Pi
  end if
  ! interpolate everything we need for the next integration step - note that there are many integration steps on one ray segment
  ! print*, "first, last rhop ", rad_ray%freq(ifreq)%svec(j)%rhop, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%rhop
  if(a == b) then
    print*, "ray did not cover any distance or problem in interpolation"
    print*, "ray points", plasma_params%ray_segment(1)%s, plasma_params%ray_segment(last_N)%s
    print*, "svec points", a, b
    stop "step size equal 0"
  end if
  last_j = j ! fall back j for the case we have to redo
  cycle_int = .false. ! if this segment has to be repated,  the large loop has to be repated
  ! this done with the cycle_int variable
  !
  do k = 1, plasma_params%int_step_cnt
    if(dstf == "relamax".or. dstf == "numeric" .or. dstf == "gene" .and. dstf_comp /= "TO") then ! .and. &
      if(j > 1) then
        if(rad_ray%freq(ifreq)%svec(j - 1)%rhop <= 1.01d0 .and. rad_ray%freq(ifreq)%svec(j)%rhop >= 0.98d0 .and. &
            rad_ray%freq(ifreq)%pol_coeff == 1.0 .and. .not.&
            rad_ray%freq(ifreq)%use_external_pol_coeff) then
            if(trim(dstf_comp) == "maxwell") then
              call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k), &
                               rad_ray%freq(ifreq)%pol_coeff, c_abs_secondary = ab_secondary(k))!
            else if(trim(dstf_comp) == "gene") then
              call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k), &
                               rad_ray%freq(ifreq)%pol_coeff, c_abs_secondary = ab_secondary(k), j_secondary = em_secondary(k))!
            else
              call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k), &
                               rad_ray%freq(ifreq)%pol_coeff)!
            end if
        else
          if(trim(dstf_comp) == "maxwell") then
            call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k), &
                             c_abs_secondary = ab_secondary(k))!
          else if(trim(dstf_comp) == "gene") then
            call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k), &
                             c_abs_secondary = ab_secondary(k), j_secondary = em_secondary(k))!
          else
            call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k))!
          end if
        end if
      else
        call abs_Albajar(rad_ray%freq(ifreq)%svec(j), omega, mode, warm_plasma, ab(k), em(k))
      end if
      if(dstf_comp == "maxwell") em_secondary(k) = ab_secondary(k) * rad_ray%freq(ifreq)%svec(j)%Ibb
    else
      call calculate_em(rad_ray%freq(ifreq)%svec(j), omega, em(k), em_secondary(k), ab(k), "relamax" )                          ! out    [W m^-3 sr^-1 Hz^-1]
    end if
    rad_ray%freq(ifreq)%svec(j)%ab = ab(k)
    ab_last = ab(k)
    if(ab_last * (b - a) > plasma_params%ab_switch .and. grid_size > 1) then
      j = last_j
      !print*, "j reset - redo step", j, last_j
      ! This step needs to be redone with a smaller step size
      !print*, "cycled"
      cycle_int = .true.
      extend_ray = .false.
      b = a
      exit
      ! if the optical depth of the last slab is very large and the step size is not already the smallest it is
      ! better to redo this step on a very fine grid
    end if
    if(output_level) then
      if(output_level) then
        N_cor = func_rel_N(omega, rad_ray%freq(ifreq)%svec(j), mode)!rad_ray%freq(ifreq)%svec(j)%N_cold
        N_cold = func_N_cold(omega, rad_ray%freq(ifreq)%svec(j), mode)
        N_gray = N_cor
        rad_ray%freq(ifreq)%svec(j)%N_cold = N_cor ! for Hutch
      end if
      if (dstf_comp == "Al" .or. dstf_comp == "Th") then
        if(mode == -1) then
          em_secondary(k) = 0.d0
          ab_secondary(k) = 0.d0
        else
          call calculate_em(rad_ray%freq(ifreq)%svec(j), omega, em_secondary(k), dummy, ab_secondary(k), dstf_comp )
        end if
      else if( dstf_comp == "TB" .or. dstf_comp == "TO") then
        if(rad_ray%freq(ifreq)%svec(j - 1)%rhop <= 1.d0 .and. rad_ray%freq(ifreq)%svec(j)%rhop >= 1.d0 .and. &
          rad_ray%freq(ifreq)%pol_coeff_secondary == 1.0 .and. .not. &
            rad_ray%freq(ifreq)%use_external_pol_coeff) then
          ab_secondary(k) = abs_Al_tor_abs(rad_ray%freq(ifreq)%svec(j), omega, mode, N_gray, &
                      pol_coeff_secondary = rad_ray%freq(ifreq)%pol_coeff_secondary)
        else
          ab_secondary(k) = abs_Al_tor_abs(rad_ray%freq(ifreq)%svec(j), omega, mode, N_gray)
        end if
        em_secondary(k) = ab_secondary(k) * rad_ray%freq(ifreq)%svec(j)%Ibb
      else if( dstf_comp == "O1") then
        ab_secondary(k) = abs_Al_tor_abs(rad_ray%freq(ifreq)%svec(j), omega, mode, N_gray)
        em_secondary(k) = ab_secondary(k) * rad_ray%freq(ifreq)%svec(j)%Ibb
      else
        em_secondary(k) = ab_secondary(k) * rad_ray%freq(ifreq)%svec(j)%Ibb
      end if
      rad_ray%freq(ifreq)%svec(j)%em_secondary = em_secondary(k)
      rad_ray%freq(ifreq)%svec(j)%ab_secondary = ab_secondary(k)
    end if
    rad_ray%freq(ifreq)%svec(j)%em = em(k)
    if(output_level) write(75,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
        rad_ray%freq(ifreq)%svec(j)%rhop, " ", rad_ray%freq(ifreq)%svec(j)%N_cold, " ", N_cold, " ",&
         N_cor," ", N_gray
    j = j + 1
    !print*, "j incremented"
  end do
  if(cycle_int) cycle
  ! For the integration of indefinite Integral tau(s) we cannot use gaussian quadrature as we not only need tau(s) after one segment,
  ! but also in between segments.
  ! To have enough points to give a smooth integration we have reuse the points of the previous segment. Then we have at least int_step_cnt + 1 points for the integration.
  ! This integration cannot be done using gaussian quadrature, because the abszissae are not at the appropriate spots.
  ! Therefore, to calculate tau(s) for segment pieces the trapezoid rule is used.
  tau_last = tau
  if(output_level) then
    tau_secondary_last = tau_secondary
    tau_secondary_int = tau_secondary
  end if
  rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%T = exp(-tau) !odd point number e.g. 41
  if(output_level) rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%T_secondary = exp(-tau_secondary)
  rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%tau = tau !odd point number e.g. 41
  if(output_level) rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%tau_s = tau_secondary
  do k = 1,plasma_params%int_step_cnt
    tau = tau + plasma_params%Int_weights(k) * ab(k) * (b - a) / 2.d0
    if(output_level) tau_secondary = tau_secondary + plasma_params%Int_weights(k) * ab_secondary(k) * (b - a) / 2.d0
  end do
  rad_ray%freq(ifreq)%svec(j - 1)%T = exp(-tau) ! even point number e.g. 50
  if(output_level) rad_ray%freq(ifreq)%svec(j - 1)%T_secondary = exp(-tau_secondary)
  rad_ray%freq(ifreq)%svec(j - 1)%tau = tau ! even point number e.g. 50
  if(output_level) rad_ray%freq(ifreq)%svec(j - 1)%tau_s = tau_secondary
  ! Do the integration for tau and tau_secondary for the entire segment using gaussian quadrature
!  j = j - plasma_params%int_step_cnt !move to the beginning of the integration e.g. 221 -> 201
!  !print*, "j decremented"
!  tau_int = tau_last !this allows the comparison of the gaussian quadrature with the rectangle rule
!  do k = 1,plasma_params%int_step_cnt
!  ! Interpolate the intermediate values of tau(s) between tau(a) and tau(b) using linear interpolation
!  ! since tau has small curvature linear interpolation works well as long as the step size is appropriate
!    if(int_mode_trans == 2) then
!      call linear_interpolation(rad_ray%freq(ifreq)%svec(j)%s, &
!        rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
!        tau_last, tau, rad_ray%freq(ifreq)%svec(j + k - 1)%s, tau_int)
!      if(tau_int < 0.d0) then
!       print*, "tau small zero", rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
!         tau_last, tau, rad_ray%freq(ifreq)%svec(j + k - 1)%s
!        stop " negative absorption"
!      end if
!    else if(int_mode_trans == 1) then
!      if(k < plasma_params%int_step_cnt) then
!        tau_int = tau_int + 0.5d0 * (rad_ray%freq(ifreq)%svec(j + 1)%s - rad_ray%freq(ifreq)%svec(j)%s) * (ab(k) + ab(k + 1))
! !     else
        tau_int = tau_int + 0.5d0 * (rad_ray%freq(ifreq)%svec(j)%s - rad_ray%freq(ifreq)%svec(j - 1)%s) * (ab(k - 1) + ab(k))
      end if
    else
      print*, "int_mode_trans", int_mode_trans
      stop "invalid int_mode_trans either 1,2"
    end if
    Trad = Trad + plasma_params%Int_weights(k) * em(k) * Exp(-tau_int) * (b - a) / 2.d0
    ! Do the gaussian quadrature for the intensity
    if(Trad < 0.d0 .or. Trad /= Trad) then
      print*, "Trad smaller 0 or NaN", ab(k), em(k),  Exp(-tau_int), plasma_params%Int_weights(k), (b - a)
      stop "radiation transport failed, negative intensities"
    end if
    if(output_level) then
      rad_ray%freq(ifreq)%svec(j + k - 1)%T = Exp(-tau_int)
      if(int_mode_trans == 2) then
        call linear_interpolation(rad_ray%freq(ifreq)%svec(j)%s, &
          rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
          tau_secondary_last, tau_secondary, rad_ray%freq(ifreq)%svec(j + k - 1)%s, tau_int)
        if(tau_int < 0.d0) then
          print*, "tau small zero", rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
            tau_last, tau, rad_ray%freq(ifreq)%svec(j + k - 1)%s
          stop " negative absorption"
        end if
      else if(int_mode_trans == 1) then
        if(k < plasma_params%int_step_cnt) then
          tau_int = tau_int + 0.5d0 * (rad_ray%freq(ifreq)%svec(j + 1)%s - rad_ray%freq(ifreq)%svec(j)%s) * (ab(k) + ab(k + 1))
        else
          tau_int = tau_int + 0.5d0 * (rad_ray%freq(ifreq)%svec(j)%s - rad_ray%freq(ifreq)%svec(j - 1)%s) * (ab(k - 1) + ab(k))
        end if
      else
        print*, "int_mode_trans", int_mode_trans
        stop "invalid int_mode_trans either 1,2"
      end if
      rad_ray%freq(ifreq)%svec(j + k - 1)%T_secondary = Exp(-tau_int)
      Trad_secondary = Trad_secondary + plasma_params%Int_weights(k) * em(k) * Exp(-tau_int) * (b - a) / 2.d0
      if(Trad_secondary < 0.d0 .or. Trad_secondary /= Trad_secondary) then
        print*, "Trad_secondary smaller 0 or NaN", ab_secondary(k), em_secondary(k),  Exp(-tau_int), plasma_params%Int_weights(k), (b - a)
        stop "radiation transport failed, negative intensities"
      end if
      !print*, "j in Trad integral", j + k - 1
    end if
  end do
  j = j  + plasma_params%int_step_cnt! move to the first point for the next integration e.g. 121 instead of 120 if int_step_cnt  = 20
  !print*, "j decremented"
  !print*, "j after Trad integral", j
  if(output_level) then
    rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt: j - 1)%Trad = Trad ! set Trad constant for the integration
    rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt: j - 1)%Trad_s = Trad_secondary  ! set Trad constant for the integration
  end if
  if(tau > plasma_params%tau_max .or. wall_hit) transport = .false.
  ! Stop radiation transport if optical thick or raytracing encountered a wall
  !print*, j, Trad
 ! if(extra_output) print*, "Extend ray? ", extend_ray
  grid_size_profiler(grid_size) = grid_size_profiler(grid_size) + 1
  total_step_cnt = total_step_cnt + 1
end do
 Calculate intensity (second integration)
 We can omit the first plasma_params%int_step_cnt points as they lie in the vessel wall
rad_ray%freq(ifreq)%total_LOS_points = j - 1
do j = plasma_params%int_step_cnt + 1, rad_ray%freq(ifreq)%total_LOS_points - plasma_params%int_step_cnt , plasma_params%int_step_cnt
  do k = 1, plasma_params%int_step_cnt
    call bounded_quadratic_interpolation(rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%s, &
         rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
         rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%tau, &
         rad_ray%freq(ifreq)%svec(j)%tau, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%tau, &
         rad_ray%freq(ifreq)%svec(j + k - 1)%s, tau_int)
    Trad = Trad  + plasma_params%Int_weights(k) * exp(-tau_int) * &
      rad_ray%freq(ifreq)%svec(j + k - 1)%em * &
      (rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s - rad_ray%freq(ifreq)%svec(j)%s) / 2.d0
    if(k /= 1 .and. k /= plasma_params%int_step_cnt) rad_ray%freq(ifreq)%svec(j + k - 1)%T = exp(-tau_int)
    if(output_level) then
      call bounded_quadratic_interpolation(rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%s, &
         rad_ray%freq(ifreq)%svec(j)%s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s, &
         rad_ray%freq(ifreq)%svec(j - plasma_params%int_step_cnt)%tau_s, &
         rad_ray%freq(ifreq)%svec(j)%tau_s, rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%tau_s, &
         rad_ray%freq(ifreq)%svec(j + k - 1)%s, tau_int)
      Trad_secondary = Trad_secondary  + plasma_params%Int_weights(k) * exp(-tau_int) * &
        rad_ray%freq(ifreq)%svec(j + k - 1)%em_secondary * &
        (rad_ray%freq(ifreq)%svec(j + plasma_params%int_step_cnt - 1)%s - rad_ray%freq(ifreq)%svec(j)%s) / 2.d0
      if(k /= 1 .and. k /= plasma_params%int_step_cnt)rad_ray%freq(ifreq)%svec(j + k - 1)%T_secondary = exp(-tau_int)
    end if
  end do
  if(output_level) then
    rad_ray%freq(ifreq)%svec(j:j + plasma_params%int_step_cnt - 1)%Trad = Trad
    rad_ray%freq(ifreq)%svec(j:j + plasma_params%int_step_cnt - 1)%Trad_s = Trad_secondary
  end if
end do
print*, "Super fine steps", grid_size_profiler(1), "/", total_step_cnt
print*, "Resonant steps", grid_size_profiler(2), "/", total_step_cnt
print*, "Kinetic steps", grid_size_profiler(3), "/", total_step_cnt
print*, "Large steps", grid_size_profiler(4), "/", total_step_cnt
if(rad_ray%freq(ifreq)%s_res /= -1.d0) then
  ! if resonance position has already been found we are done at this point
  rad_ray%freq(ifreq)%total_LOS_points = j - 1 ! here we want the last point and not the first point of the next step
  if(extra_output) stop "ray_done"
  if(output_level) then
    close(75)
    close(74)
  end if
  return
end if
! for Doppler-shift dominated channels we need to propagate a little further to find the cold resonance position
if(.not. extend_ray) then
! finish of the remaining ray if there is any
  a = b
  call prepare_svec_segment_gauss(plasma_params, omega, rad_ray, ifreq, plasma_params%ray_segment, last_N, j, &
    j + plasma_params%int_step_cnt - 1, (/2.d0, 2.d0, 2.d0, 2.d0/), grid_size, 0.d0, extend_ray, a, b, .false.)
  j = j + plasma_params%int_step_cnt
  print*, "j incremented"
end if
do while(rad_ray%freq(ifreq)%s_res == -1.d0 .and. wall_hit == .false.) ! no need to look for resonances in the vessel wall
! go along the ray to find the resonance position
  a = b
  if(extend_ray) then
    plasma_params%ray_segment(1) = plasma_params%ray_segment(last_N) ! copy last part to ray to start of next
    call make_ray_segment(distance, plasma_params, omega, plasma_params%ray_segment, last_N, wall_hit) ! prepare the next segment
    if(output_level) then
      do N=1, last_N
        write(74,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
          plasma_params%ray_segment(N)%s, " ", plasma_params%ray_segment(N)%R_vec(1), " ", plasma_params%ray_segment(N)%R_vec(2), " ",&
          plasma_params%ray_segment(N)%R_vec(3)," ", plasma_params%ray_segment(N)%Hamil
      end do
    end if
    ray_segment_cnt = ray_segment_cnt + 1
    if(extra_output) print*,"ray segments (for resonance position search):", ray_segment_cnt
  end if
  !print*,"distance covered, target distance:", plasma_params%ray_segment(last_N)%s - plasma_params%ray_segment(1)%s, distance
  if(extra_output) print*, "interpolating things for the resonance"
  call prepare_svec_segment_gauss(plasma_params, omega, rad_ray, ifreq, plasma_params%ray_segment, last_N, j, &
    j + plasma_params%int_step_cnt - 1, (/2.d0, 2.d0, 2.d0, 2.d0/), grid_size, 0.d0, extend_ray, a, b, .false.)
  ! interpolate the quantities on the ray - this is a general routine and an optizimed routine that only looks for the resonance would be more optimal
  j = j + plasma_params%int_step_cnt
  !print*, "j incremented"
end do
rad_ray%freq(ifreq)%total_LOS_points = j- 1

! If the channel is in cut-off we arrive here. The cold-resonance is not reached by the curved ray
! and the best we can do is assign the measurement to the point on the ray, which is closest to the cold resonance
if(rad_ray%freq(ifreq)%s_res == -1.d0) then
  Y_closest_res = 1.d0
  rad_ray%freq(ifreq)%total_LOS_points = j- 1
  do j = 1, rad_ray%freq(ifreq)%total_LOS_points
    if(abs(rad_ray%freq(ifreq)%svec(j)%freq_2X / omega * Pi - plasma_params%Y_res) < Y_closest_res) then
      rad_ray%freq(ifreq)%s_res = rad_ray%freq(ifreq)%svec(j)%s
      rad_ray%freq(ifreq)%R_res = rad_ray%freq(ifreq)%svec(j)%R
      rad_ray%freq(ifreq)%z_res = rad_ray%freq(ifreq)%svec(j)%z
      rad_ray%freq(ifreq)%rhop_res = rad_ray%freq(ifreq)%svec(j)%rhop
      Y_closest_res = abs(rad_ray%freq(ifreq)%svec(j)%freq_2X / omega * Pi - plasma_params%Y_res)
    end if
  end do
  print*, "Warning channel in cut-off"
  if(extra_output) stop "ray_done"
end if
if(output_level) then
  close(75)
  close(74)
end if
end subroutine integrate_Trad

!subroutine calculate_Trad_segmented(rad_ray_freq, omega, mode, outer_steps, tau_max, Trad, tau, error, debug)
! IDEA DOES NOT WORK - DO NOT USE!
! Solves the radiation transport equation
! For numerical efficiency the radiation transp. equation is split into multiple sections, with the first one being closest to the antenna.
! This allows to quit the radiation transport early if tau_max has been reached.
! Since this method is designed for much speed it is only designed for output_level = .false.
!use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_type, output_level, tau_thick
!use constants,                  only: pi, e0, mass_e, eps0, c0
!use mod_ecfm_refr_abs_Al,       only: abs_Albajar_fast, func_N_cold
!implicit none
!type(rad_diag_ch_mode_ray_freq_type), intent(inout) :: rad_ray_freq
!real(rkind),                intent(in)    :: omega, tau_max
!integer(ikind),             intent(in)    :: mode, outer_steps
!real(rkind),                intent(out)   :: Trad
!real(rkind),                intent(out)   :: tau
!integer(ikind),             intent(out)   :: error
!logical,    intent(in), optional          :: debug
!real(rkind)    :: ds, ds2
!real(rkind), dimension(3) :: em, ab, em_secondary, ab_secondary
!real(rkind), dimension(outer_steps) :: Trad_segment, tau_segment
!logical                             :: cut_off, debug_internal, eval_pol_coeff
!real(rkind)    :: k1, k2, k3, k4, em_diff, dummy, N_cor, max_ab_segment, max_em_segment, &
!                  ab_diff,trans_coeff,trans_coeff_secondary, tau_stop
!integer(ikind) :: i, j, k, outer_pos, i_outer_start, i_outer_end, last_outer_pos
!error = 0
! initialization
!Trad     =  0.d0
!tau      =  0.d0
!Trad_segment(:) = 0.d0
!tau_segment(:) = 0.d0
!cut_off = .false.
!debug_internal = .false.
!if(present(debug)) debug_internal = debug
!outer_pos = 1
!tau_stop = 9.d0 ! a slab with this optical depth absorbs strongly enough to forget all radiation emitted further inside the plasma than the slab
!j = rad_ray_freq%total_LOS_points - rad_ray_freq%total_LOS_points / outer_steps * outer_pos + 1! start close to the antenna
! start with the first segment which is closest to the antenna
!k = 1
!ds2  = abs(rad_ray_freq%svec(j + 1)%s - rad_ray_freq%svec(j)%s)   ! ds for all but last point
!ab(:) = 0.d0
!em(:) = 0.d0
! print*, "ab start, tau start", ab(k), tau_segment(outer_pos)
! loop over all segments starting with the segment closest to the antenna and getting further and further inside of the plasma
!do while(outer_pos <= outer_steps)
!  ! loop over all indices
!  !----------------------
!  i_outer_start = rad_ray_freq%total_LOS_points - rad_ray_freq%total_LOS_points / outer_steps * outer_pos + 1
!  i_outer_end = rad_ray_freq%total_LOS_points - rad_ray_freq%total_LOS_points / outer_steps * (outer_pos - 1)
!  max_ab_segment = 0.d0
!  max_em_segment = 0.d0
!  !print*, "ratio of segment", rad_ray_freq%svec(i_outer_start)%R, rad_ray_freq%svec(i_outer_end)%R
!  if(.not. any(rad_ray_freq%svec(i_outer_start : i_outer_end )%plasma)) then
!    Trad_segment(outer_pos)  = 0.d0
!    tau_segment(outer_pos) = 0.d0
!    outer_pos = outer_pos + 1
!    cycle
!  end if
!  ! We cannot reuse em(1) from the last segment, since we just had a jump
!  k = 1
!  j = i_outer_start
!  if(.not. debug_internal) then
!    call abs_Albajar_fast(rad_ray_freq%svec(j), omega, mode, ds2, ab(k))
!    em(k) = ab(k) * rad_ray_freq%svec(j)%Ibb
!  else
!    if(j > 1) then
!      eval_pol_coeff = rad_ray_freq%svec(j)%rhop <= 1.d0 .and. rad_ray_freq%svec(j - 1)%rhop >= 0.98d0 .and. .not.  &
!        rad_ray_freq%use_external_pol_coeff
!    else
!      eval_pol_coeff = .False.
!    end if
!    if(j > 1)  then
!      rad_ray_freq%svec_extra_output(j)%N_warm = rad_ray_freq%svec_extra_output(j - 1)%N_warm
!    else
!      rad_ray_freq%svec_extra_output(j)%N_warm = rad_ray_freq%svec(j)%N_cold
!    end if
!    call evaluate_em_ab_single(rad_ray_freq, j, omega, mode, ds2, eval_pol_coeff, &
!                     em(k), ab(k), em_secondary(k), ab_secondary(k), rad_ray_freq%pol_coeff, &
!                     rad_ray_freq%pol_coeff_secondary)
!  end if
!  if(j == 1) tau_segment(outer_pos) = 0.5d0 * ab(k) * ds2
!  do i = i_outer_start, i_outer_end-2, 2           ! integration over every second point on LOS
!    ! step size
!    !----------
!    if(.not. all(rad_ray_freq%svec(i + 1:i + 2)%plasma)) then
!      em(:) = 0.d0
!      ab(:) = 0.d0
!      !print*, "Skipped step -  outside of plasma"
!      cycle
!    end if
!    ds  = abs(rad_ray_freq%svec(i+2)%s - rad_ray_freq%svec(i)%s)
!    ds2 = ds / 2.d0 ! ds for all but last point
!    !------------------------------------------------
!    do k = 2, 3             ! 3 Runge-Kutta evaluations; the first comes from previous iteration
!      j = i + k - 1         ! corresponding indices on svec
!      if(.not. debug_internal) then
!        call abs_Albajar_fast(rad_ray_freq%svec(j), omega, mode, ds2, ab(k))
!        em(k) = ab(k) * rad_ray_freq%svec(j)%Ibb
!      else
!        eval_pol_coeff = rad_ray_freq%svec(j)%rhop <= 1.d0 .and. rad_ray_freq%svec(j - 1)%rhop >= 0.98d0 .and. .not.  &
!          rad_ray_freq%use_external_pol_coeff
!        if(output_level) rad_ray_freq%svec_extra_output(j)%N_warm = rad_ray_freq%svec_extra_output(j - 1)%N_warm
!        call evaluate_em_ab_single(rad_ray_freq, j, omega, mode, ds2, eval_pol_coeff, &
!                         em(k), ab(k), em_secondary(k), ab_secondary(k), rad_ray_freq%pol_coeff, &
!                         rad_ray_freq%pol_coeff_secondary)
!      end if
!      if(max_ab_segment < ab(k)) max_ab_segment = ab(k)
!      if(max_em_segment < em(k)) max_em_segment = em(k)
!    enddo !k = 1, 3
!    em(2:3) = em(2:3) * c0**2 / (ant%diag(idiag)%ch(ich)%f_ECE**2 * e0)
!    k1   = em(1) - ab(1) * Trad_segment(outer_pos) !dN_abs_2_ds(2) / N_abs_2(2)
!    k2   = em(2) - ab(2) * (Trad_segment(outer_pos) + k1*ds2)
!    k3   = em(2) - ab(2) * (Trad_segment(outer_pos) + k2*ds2)
!    k4   = em(3) - ab(3) * (Trad_segment(outer_pos) + k3*ds)
!    Trad_segment(outer_pos) = Trad_segment(outer_pos) + ds/6.d0 * (k1 + 2.d0*(k2 + k3) + k4)
!    if(Trad_segment(outer_pos) < 0.d0) then
!      Trad_segment(outer_pos) = rad_ray_freq%svec(j)%Ibb
!      if(debug_internal) print*, "Warning - negative Trad encountered!"
!    end if
!    if(rad_ray_freq%svec(j)%N_cold <= 0) then  ! Cut_pff
!      Trad_segment(outer_pos) = 0.d0
!      tau_segment(outer_pos) = 0.d0
!      cut_off = .true.
!    end if
!    tau_segment(outer_pos) = tau_segment(outer_pos)  + (ab(1) + ab(2)) * ds2
!    if(debug_internal) then
!      j = i
!      if(j + 2 > rad_ray_freq%total_LOS_points) j = rad_ray_freq%total_LOS_points - 2
!      print*, "-----------------", j, "-th step ---------------"
!      print*, "Trad", Trad_segment(outer_pos)
!      print*, "tau", tau_segment(outer_pos)
!      print*, "em",em
!      print*, "ab",ab
!      print*, "Te", rad_ray_freq%svec(j)%Te, rad_ray_freq%svec(j+1)%Te, rad_ray_freq%svec(j+2)%Te
!      print*, "ne", rad_ray_freq%svec(j)%ne, rad_ray_freq%svec(j+1)%ne, rad_ray_freq%svec(j+2)%ne
!    end if
!    if((ab(1) + ab(2)) * ds > tau_thick) then
!    ! For strongly radiating plasmas we get into trouble with numerical stability.
!    ! If, however, the optical depth of the last plasma slab is large we can fall back on just using
!    ! the black body intensity - for tau = 3 the error relative error is about 5 %
!      Trad_segment(outer_pos) = rad_ray_freq%svec(j)%Ibb
!    else if(Trad_segment(outer_pos) > 1.d-5 .and. (ab(1) + ab(2)) * ds > tau_thick * 0.5d0) then
!      print*, "Warning Trad was diverging, although optical depth of last step below tau_thick", (ab(1) + ab(2)) * ds* 0.5d0
!      Trad_segment(outer_pos) = rad_ray_freq%svec(j)%Ibb
!    else if(Trad_segment(outer_pos) > 1.d-5) then
!       Trad_segment(outer_pos) = 1.d0  ! certainly triggers error below
!    end if
!    if(debug_internal) rad_ray_freq%svec_extra_output(i: i + 2)%Trad = Trad  + Trad_segment(outer_pos) * exp(-tau)
!    if(Trad_segment(outer_pos) /= Trad_segment(outer_pos) .or. Trad_segment(outer_pos) > 1.d-4 .or. &
!      tau_segment(outer_pos) < 0.d0) then !/= Trad
!      ! Trad >> 1 MeV => much larer than saturation value of most ECEs
!      j = i
!      if(j + 2 > rad_ray_freq%total_LOS_points) j = rad_ray_freq%total_LOS_points - 2
!      print*, "Trad", Trad_segment(outer_pos)
!      print*, "I_bb", rad_ray_freq%svec(j)%Te / c0**2 * (omega**2 * e0) * 4.d0 * pi**2
!      print*, "em",em
!      print*, "ab",ab
!      print*, "optical depth of ds", ab * ds2
!      print*, "optical depth of segment", tau_segment(outer_pos)
!      print*, "Te", rad_ray_freq%svec(j)%Te, rad_ray_freq%svec(j+1)%Te, rad_ray_freq%svec(j+2)%Te
!      print*, "ne", rad_ray_freq%svec(j)%ne, rad_ray_freq%svec(j+1)%ne, rad_ray_freq%svec(j+2)%ne
!      print*, "theta", rad_ray_freq%svec(j)%theta / pi * 180, rad_ray_freq%svec(j+1)%theta / pi * 180, rad_ray_freq%svec(j+2)%theta / pi * 180
!      print*, "freq_2X", rad_ray_freq%svec(j)%freq_2X, rad_ray_freq%svec(j+1)%freq_2X, rad_ray_freq%svec(j+2)%freq_2X
!      print*, "N_abs",  rad_ray_freq%svec(j)%N_cold, rad_ray_freq%svec(j+1)%N_cold, rad_ray_freq%svec(j+2)%N_cold
!      print*, "j", j, j + 2
!      error = -1
!      return
!    end if
!    ! emissivity and absorption to be used in next iteration
!    !-------------------------------------------------------
!    em(1) = em(3)
!    ab(1) = ab(3)
!  enddo !i = 1, rad_ray_freq%total_LOS_points           ! integration over all points on LOS
!  if(outer_pos == rad_ray_freq%total_LOS_points) then
!    tau_segment(outer_pos) = tau_segment(outer_pos) + 0.5d0 * ab(3) * ds2
!  else
!    ! If we are in cutoff or if this was the last segment on the LOS due to large tau
!    ! we make a mathematical mistake here by omitting the 0.5d0 prefactor.
!    ! However, if we abort the radiation transport due to sufficiently large tau or due to
!    ! cutoff this mistake has no consequennces.
!    tau_segment(outer_pos) = tau_segment(outer_pos) + ab(3) * ds2
!  end if
!  if(debug_internal) rad_ray_freq%svec_extra_output(i_outer_start : i_outer_end )%T = exp(-tau)
!  if(outer_pos == 1) then
!    Trad = Trad_segment(outer_pos)
!    tau = tau_segment(outer_pos)
!  else
!    Trad = Trad + Trad_segment(outer_pos) * exp(-tau) ! cummulative optical depth of all PREVIOUS segments
!    tau = tau + tau_segment(outer_pos)
!  end if
!  if(debug_internal) then
!    rad_ray_freq%svec_extra_output(i_outer_end)%em = em(3)
!    rad_ray_freq%svec_extra_output(i_outer_end)%ab = ab(3)
!  end if
!  print*,"segment No.", outer_pos,  "tau", tau_segment(outer_pos), "Trad [eV]", Trad_segment(outer_pos)
!  print*,"total: ", "tau", tau, "Trad [eV]", Trad
!  print*, "Max ab segment, max em segment", max_ab_segment, max_em_segment
!  if(tau > tau_max .or. cut_off) exit ! leave the loop early because maximum optical depth achieved or cut off
!  outer_pos = outer_pos + 1
!end do ! while not all segments iterated
!end subroutine calculate_Trad_segmented


! Outdated and needs fixing before compileable
!  subroutine calculate_Irad_euler(rad_ray, ich, ifreq, ir, Irad, tau, em_max, R_max)
!    use mod_ecfm_refr_types,        only: ant, rad_diag_ch_mode_ray_type, dstf_comp, output_level, data_folder, Ich_name
!    use mod_ecfm_refr_em_Hu,                  only: calculate_em, calculate_N
!    use constants,                  only: pi, e0, mass_e, eps0, c0
!    use mod_ecfm_refr_abs_Al,         only: abs_Albajar, abs_Al_tor_abs,abs_Al_all_N
!    implicit none
!    type(rad_diag_ch_mode_ray_type), intent(inout) :: rad_ray
!    integer(ikind),             intent(in)    :: ich, ifreq, ir
!    real(rkind),                intent(out)   :: Irad, tau, em_max, R_max
!    character(120)               :: cur_filename, cur_filename_3X
!    real(rkind)    :: em, ab, ds,em_old, abs_old, omega, omega_c, omega_p, N, N_2, N_3, N_4, Irad_tot, tau_tot
!    real(rkind)    :: cputime1, cputime2
!    integer(ikind) :: i
!    omega = ant%diag(idiag)%ch(ich)%freq(ifreq)%freq * 2 * Pi
!    Irad    = 0.d0
!    tau     = 0.d0
!    em_max  = 0.d0
!    R_max   = 0.d0
!    if(output_level .and. ifreq == 1) then
!        end if
!        open(66, file=cur_filename)
!        !write(cur_filename_3X, "(A64A5A11I3.3A4)") data_folder,Ich_name,"/Irhopchs3X", ich,".dat"
!        !open(67, file=cur_filename_3X)
!      end if
!      !call cpu_time(cputime1)
!      do i = 1, rad_ray%freq(ifreq)%total_LOS_points
!        !if(rad%diag(idiag)%freq(ifreq)%svec(i)%in_dense_interval) then
!        !  ds = rad%diag(idiag)%ds2
!        !else
!        !  ds = rad%diag(idiag)%ds1
!        !endif
!        if (i == rad_ray%freq(ifreq)%total_LOS_points ) then
!          ds = abs(rad_ray%freq(ifreq)%svec(i)%s - rad_ray%freq(ifreq)%svec(i-1)%s)
!        else
!          ds = abs(rad_ray%freq(ifreq)%svec(i+1)%s - rad_ray%freq(ifreq)%svec(i)%s)
!        endif
!
!        !write(99,'(i5,l2,3e14.6)')i,rad%diag(idiag)%freq(ifreq)%svec(i)%in_dense_interval, ds, rad%diag(idiag)%freq(ifreq)%svec(i)%s,  rad%diag(idiag)%freq(ifreq)%svec(i)%rhop
!
!        !call interpol_LOS(rad%diag(idiag)%rhop, ich, ir,                  & ! in
!        !                  s_in  = rad%diag(idiag)%freq(ifreq)%svec(i)%s,  & ! in
!        !                  ia_in = rad%diag(idiag)%freq(ifreq)%svec(i)%ia, & ! in
!        !                  ne    = rad%diag(idiag)%freq(ifreq)%svec(i)%ne, & ! out
!        !                  Te    = rad%diag(idiag)%freq(ifreq)%svec(i)%Te, & ! out
!        !                  Ibb   = rad%diag(idiag)%freq(ifreq)%svec(i)%Ibb)  ! out
!        !  ! Uses a predefined ia to avoid search for the s
!        !if (rad_ray%freq(ifreq)%svec(i)%ne   >= ant%diag(idiag)%ch(ich)%freq(ifreq)%cutoff_density_2X .or. &
!        ! rad_ray%freq(ifreq)%svec(i)%rhop >  rad_ray%rhop_ow) then
!          if(dstf_comp /= "TB" .and. dstf_comp /= "Al" .and. dstf_comp /= "O1") then
!            call calculate_em(rad_ray%freq(ifreq)%svec(i),ich,ifreq,i, em, em_old, ab )
!            em = em_old
!            ab = em / rad_ray%freq(ifreq)%svec(i)%Ibb
!          else if( dstf_comp == "Al") then
!            call abs_Albajar(rad_ray%freq(ifreq)%svec(i), omega, ab, 2)
!            em = ab * rad_ray%freq(ifreq)%svec(i)%Ibb
!          else if( dstf_comp == "TB") then
!            ab = abs_Al_tor_abs(rad_ray%freq(ifreq)%svec(i), omega, -1)
!            em = ab * rad_ray%freq(ifreq)%svec(i)%Ibb
!          else
!            ab = abs_Al_tor_abs(rad_ray%freq(ifreq)%svec(i), omega, 1)
!            em = ab * rad_ray%freq(ifreq)%svec(i)%Ibb
!          end if
!          omega_c        = rad_ray%freq(ifreq)%svec(i)%freq_2X * Pi
!          omega_p = sqrt( (rad_ray%freq(ifreq)%svec(i)%ne * e0**2.d0)/(eps0 * mass_e))
!        if(dstf_comp == "O1") then
!            call abs_Al_all_N(omega, omega_c, omega_p, rad_ray%freq(ifreq)%svec(i)%theta, &
!            rad_ray%freq(ifreq)%svec(i)%sin_theta, rad_ray%freq(ifreq)%svec(i)%cos_theta,-1,1, N,N_2,N_3)
!          else
!            call abs_Al_all_N(omega, omega_c, omega_p, rad_ray%freq(ifreq)%svec(i)%theta, &
!            rad_ray%freq(ifreq)%svec(i)%sin_theta, rad_ray%freq(ifreq)%svec(i)%cos_theta,1,2, N,N_2,N_3)
!          end if
!
!        !if (rad_ray%freq(ifreq)%svec(i)%ne >= ant%diag(idiag)%ch(ich)%freq(ifreq)%cutoff_density_2X) Irad = 0.d0
!        if (N == 0.d0) then
!          Irad = 0.d0
!          tau = 0.d0
!        end if
!
!        !rad_ray%freq(ifreq)%svec(i)%em   = em
!        !rad_ray%freq(ifreq)%svec(i)%Irad = Irad
!        !rad_ray%freq(ifreq)%svec(i)%tau  = tau
!
!        if (i == 1 .or. i == rad_ray%freq(ifreq)%total_LOS_points ) then
!          Irad = Irad + (em - ab * Irad) * 0.5d0 * ds
!          tau  = tau  + ab * 0.5d0 * ds
!        else
!          Irad = Irad + (em - ab * Irad) * ds
!          tau  = tau  + ab * ds
!        endif
!        !print*,s,ds2,R,rhop,rad%diag(idiag)%freq(ifreq)%svec(i)%ne,rad%diag(idiag)%freq(ifreq)%svec(i)%Te,rad%diag(idiag)%freq(ifreq)%svec(i)%Ibb,em,ab,Irad
!        !write(98,'(5i4,4e12.4)')i, ich, ifreq, ir, rad%diag(idiag)%freq(ifreq)%svec(i)%ia, rad%diag(idiag)%freq(ifreq)%svec(i)%s, rad%diag(idiag)%freq(ifreq)%svec(i)%ne, rad%diag(idiag)%freq(ifreq)%svec(i)%Te, Ibb
!
!        if (em_max < em) then
!          em_max = em
!          R_max  = rad_ray%freq(ifreq)%svec(i)%R
!        endif
!        write(66,"(E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3,A1,E14.6E3)") &
!          rad_ray%freq(ifreq)%svec(i)%s, " ",&
!          Irad, " ", em, " ",ab, " ", tau, " ", N, " ",N_2," ",N_3
!        !write(98,'(7e12.4)')rad%diag(idiag)%freq(ifreq)%svec(i)%ne,rad%diag(idiag)%freq(ifreq)%svec(i)%Te,em,ab,rad%diag(idiag)%freq(ifreq)%svec(i)%Ibb,Irad,tau
!      enddo
!      !call cpu_time(cputime2)
!      !print*,'time needed for calculate_Irad_euler',cputime2-cputime1,'s'
!    !Irad = Irad * ant%diag(idiag)%ch(ich)%freq(ifreq)%ff
!    !tau  = tau* ant%diag(idiag)%ch(ich)%freq(ifreq)%ff
!    if(output_level .and. ifreq == 1) then
!      close(66)
!      !close(67)
!    end if
!  !stop 'subroutine calculate_Irad_simple_int'
!  end subroutine calculate_Irad_euler

end module mod_ecfm_refr_rad_transp_int
