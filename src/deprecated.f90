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
subroutine sub_spatial_grad_rhop_ana(plasma_params, x_vec, rhop, spatial_grad_rhop)
    USE f90_kind
    USE mod_ecfm_refr_types ,     only: plasma_params_type
    USE ripple3d,                 only: grad_type
    use mod_ecfm_refr_utils,      only: sub_remap_coords
    !USE nag_spline_2d             , only: nag_spline_2d_eval, &
    !                                      nag_error, nag_set_error
    implicit none
    type(plasma_params_type)                   :: plasma_params
    real(rkind), dimension(:)  , intent(in)    :: x_vec
    real(rkind)                , intent(out)   :: rhop
    real(rkind), dimension(3)  , intent(out)   :: spatial_grad_rhop
    real(rkind), dimension(3)                  :: R_vec
    call sub_remap_coords(x_vec, R_vec)
    spatial_grad_rhop(1) = func_dR_dx(x_vec(1),x_vec(2)) * (3.d0*(-1.5d0 + R_vec(1)))/(2.*Sqrt((-1.5d0 + R_vec(1))**2 + R_vec(3)**2/10.d0))
    spatial_grad_rhop(2) = func_dR_dy(x_vec(1),x_vec(2)) * (3.d0*(-1.5d0 + R_vec(1)))/(2.*Sqrt((-1.5d0 + R_vec(1))**2 + R_vec(3)**2/10.d0))
    spatial_grad_rhop(3) = (3.d0*R_vec(3))/(20.d0*Sqrt((-1.5d0 + R_vec(1))**2 + R_vec(3)**2/10.d0))
    rhop =  (3.d0*Sqrt((-1.5d0 + R_vec(1))**2 + R_vec(3)**2/10.d0))/2.d0
  end subroutine sub_spatial_grad_rhop_ana

  subroutine sub_grad_n_e_ana(plasma_params, rhop, n_e,  grad_n_e)
  ! dne / drhop
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type,  h_x_glob
    implicit none
    type(plasma_params_type)   , intent(in)  :: plasma_params
    real(rkind)                , intent(in)  :: rhop
    real(rkind)                , intent(out) :: n_e, grad_n_e
    integer(ikind)                           :: i
    real(rkind)                              :: rhop_debug
    n_e = (3.d19*(1.5625d0 - rhop**2)**2)/exp(rhop**2)
    grad_n_e = (-12.d19*rhop*(1.5625d0 - rhop**2))/exp(rhop**2) - &
               (6.d19*rhop*(1.5625d0- rhop**2)**2)/exp(rhop**2)
!    if(debug_level > 1) then
!      print*, grad_n_e
!      rhop_debug = rhop + h_x_glob
!      print*, ((3.d19*(1.5625d0 - rhop_debug**2)**2)/exp(rhop_debug**2) - (3.d19*(1.5625d0 - rhop**2)**2)/exp(rhop**2))/ h_x_glob
!      stop "ne"
!    end if
  end subroutine sub_grad_n_e_ana

  subroutine sub_grad_T_e_ana(plasma_params, rhop, T_e,  grad_T_e)
  ! dne / drhop
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type,  h_x_glob
    !use interpolation_routines, only : linear_interpolation
    implicit none
    type(plasma_params_type)   , intent(in)  :: plasma_params
    real(rkind)                , intent(in)  :: rhop
    real(rkind)                , intent(out) :: T_e, grad_T_e
    integer(ikind)                           :: i
    real(rkind)                              :: rhop_debug
    T_e = (3.d3*(1.5625d0 - rhop**2)**2)/exp(rhop**2)
    grad_T_e = (-12.d3*rhop*(1.5625d0 - rhop**2))/exp(rhop**2) - &
               (6.d3*rhop*(1.5625d0- rhop**2)**2)/exp(rhop**2)
!    if(debug_level > 1) then
!      print*, grad_n_e
!      rhop_debug = rhop + h_x_glob
!      print*, ((3.d19*(1.5625d0 - rhop_debug**2)**2)/exp(rhop_debug**2) - (3.d19*(1.5625d0 - rhop**2)**2)/exp(rhop**2))/ h_x_glob
!      stop "ne"
!    end if
  end subroutine sub_grad_T_e_ana

 function func_n_e_ana(plasma_params, x_vec)
  ! dne / drhop
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    !USE mod_ecfm_refr_types , only h_rhop ! assuming that taking the derivative
    !                                        without higher order interpolation is sufficient
    implicit none
    type(plasma_params_type)   , intent(in)  :: plasma_params
    real(rkind), dimension(:)  , intent(in)  :: x_vec
    real(rkind)                              :: func_n_e_ana
    real(rkind)                              :: rhop
    rhop = func_rhop_ana(plasma_params, x_vec)
    func_n_e_ana = (3.d19*(1.5625d0 - rhop**2)**2)/exp(rhop**2)
  end function func_n_e_ana

  function func_T_e_ana(plasma_params, x_vec)
  ! dne / drhop
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    !USE mod_ecfm_refr_types , only h_rhop ! assuming that taking the derivative
    !                                        without higher order interpolation is sufficient
    implicit none
    type(plasma_params_type)   , intent(in)  :: plasma_params
    real(rkind), dimension(:)  , intent(in)  :: x_vec
    real(rkind)                              :: func_T_e_ana
    real(rkind)                              :: rhop
    rhop = func_rhop_ana(plasma_params, x_vec)
    func_T_e_ana = (3.d3*(1.5625d0 - rhop**2)**2)/exp(rhop**2) ! same profile shape as ne with 8 keV core Te
  end function func_T_e_ana

  function func_B_abs_ana(plasma_params, x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind)                               :: func_B_abs_ana
    real(rkind), dimension(3)                 :: R_vec, B_x_vec
    B_x_vec(1) = func_B_x_ana(x_vec)
    B_x_vec(2) = func_B_y_ana(x_vec)
    B_x_vec(3) = 0.d0
    func_B_abs_ana = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
  end function func_B_abs_ana

    subroutine sub_N_par_ana(plasma_params, x_vec, N_vec, N_par, N_abs)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type
    Use ripple3d,             only : grad_type
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type), intent(in)      :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec, N_vec
    real(rkind), intent(out)                  :: N_par, N_abs
    real(rkind), dimension(3)                 :: R_vec, B_x_vec, B_R_vec
    real(rkind)                               :: cos_phi, sin_phi, scal_prod, B_abs, ddx, ddy
    integer(ikind)                            :: i
    call sub_remap_coords(x_vec, R_vec)
    cos_phi = cos(R_vec(2))
    sin_phi = sin(R_vec(2))
    N_abs = sqrt(N_vec(1)**2 + N_vec(2)**2 + N_vec(3)**2)
    call sub_B_r_ana( R_vec, B_R_vec)
    B_x_vec(1) = B_R_vec(1) * cos_phi - sin_phi * B_R_vec(2)
    B_x_vec(2) = B_R_vec(1) * sin_phi + cos_phi * B_R_vec(2)
    B_x_vec(3) = B_R_vec(3)
    B_abs = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
    scal_prod = 0.d0
    do i = 1, 3
      scal_prod = scal_prod + N_vec(i) * B_x_vec(i)
    end do
    N_par = scal_prod / B_abs
   end subroutine sub_N_par_ana

subroutine sub_B_r_ana(R_vec, B_R_vec)
    use f90_kind
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(:)  , intent(out)   :: B_R_vec
    B_R_vec(:) = 0.d0
    B_R_vec(2) = 4.25 / R_vec(1)
  end subroutine sub_B_r_ana

  function  func_B_x_ana(x_vec)
    use f90_kind
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3) :: R_vec, B_R_vec
    real(rkind)                :: func_B_x_ana
    call sub_remap_coords(x_vec, R_vec)
    call sub_B_r_ana(R_vec, B_R_vec)
    func_B_x_ana =  B_R_vec(1) * cos(R_vec(2)) - B_R_vec(2) * sin(R_vec(2))
  end function func_B_x_ana

  function  func_B_y_ana(x_vec)
    use f90_kind
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3) :: R_vec, B_R_vec
    real(rkind)                :: func_B_y_ana
    call sub_remap_coords(x_vec, R_vec)
    call sub_B_r_ana(R_vec, B_R_vec)
    func_B_y_ana =  B_R_vec(1) * sin(R_vec(2)) + B_R_vec(2) * cos(R_vec(2))
  end function func_B_y_ana

 function  func_B_x_ana_R(R_vec)
    use f90_kind
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3) :: B_R_vec
    real(rkind)               :: func_B_x_ana_R
    call sub_B_r_ana(R_vec, B_R_vec)
    func_B_x_ana_R =  B_R_vec(1) * cos(R_vec(2)) - B_R_vec(2) * sin(R_vec(2))
  end function func_B_x_ana_R

  function  func_B_y_ana_R(R_vec)
    use f90_kind
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: R_vec
    real(rkind), dimension(3)  :: B_R_vec
    real(rkind)                :: func_B_y_ana_R
    call sub_B_r_ana(R_vec, B_R_vec)
    func_B_y_ana_R =  B_R_vec(1) * sin(R_vec(2)) + B_R_vec(2) * cos(R_vec(2))
  end function func_B_y_ana_R

  subroutine sub_grad_N_par_ana(plasma_params, x_vec, N_vec, N_abs, B_abs, spatial_grad_B_abs, &
                           N_par, N_grad_N_par, spatial_grad_N_par) !spatial_grad_B_x, spatial_grad_B_y, spatial_grad_B_z,
  ! Gradient of vec(B) along LOS coordinates x ripple not (yet) included => dB_vec/dphi = 0
  ! Also calculates N_par, d N_par(theta)/dN_i, d N_par/dx_i, d|B|dx_i since all information is readily available
  ! (Calculating these quantities here safes interpolations)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type, h_x_glob
    USE ripple3d,                 only: grad_type
    use mod_ecfm_refr_utils,  only: sub_remap_coords
    implicit none
    type(plasma_params_type)   , intent(in)   :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec, N_vec
    real(rkind)                , intent(in)   :: N_abs
    real(rkind), dimension(:), intent(out)    :: spatial_grad_B_abs
    real(rkind)                , intent(out)  :: B_abs, N_par
    real(rkind), dimension(:)  , intent(out)  :: N_grad_N_par, spatial_grad_N_par
    real(rkind), dimension(3,3)               :: spatial_grad_B
    real(rkind)                               :: cos_phi_tok, sin_phi_tok, alpha, scal_prod, h_x
    type(grad_type)                           :: dB_r_inter, dB_t_inter, dB_z_inter
    real(rkind), dimension(3)                 :: R_vec, B_R_vec, B_x_vec!, N_vec_norm, B_vec_norm
    real(rkind), dimension(4,3)               :: aux_x, aux_B_R, aux_R
    real(rkind), dimension(3)                 :: dB_x_dR, dB_y_dR, dB_z_dR
    integer(ikind)                            :: i, j, k
    h_x = h_x_glob
    call sub_remap_coords(x_vec, R_vec)
    call sub_B_r_ana(R_vec, B_R_vec)
    cos_phi_tok = cos(R_vec(2))
    sin_phi_tok = sin(R_vec(2))
    B_x_vec(1) = B_R_vec(1) * cos_phi_tok - sin_phi_tok * B_R_vec(2)
    B_x_vec(2) = B_R_vec(1) * sin_phi_tok + cos_phi_tok * B_R_vec(2)
    B_x_vec(3) = B_R_vec(3)
    dB_r_inter%dR = 0.d0
    dB_r_inter%dz = 0.d0
    dB_t_inter%dR = -4.25 * 1.0/R_vec(1)**2
    dB_t_inter%dz = 0.d0
    dB_z_inter%dR = 0.d0
    dB_z_inter%dz = 0.d0
    ! Next dB_x_j/d_x_i --> j (first index) direction of magnetic field
    !                   --> i (second index) direction of the derivative
    ! To do this analyitaclly the multidimensional chain rule is used
    ! First dvec(B_x)dR and dvec(B_x)dz are computed
    ! Since no ripple dvec(B)/dphi= 0, hence we only need to sum over dvec(B)/dR and dvec(B)/dz
    ! dB_x/dR
    dB_x_dR(1) = dB_r_inter%dR * cos_phi_tok - dB_t_inter%dR * sin_phi_tok
    ! dB_x/dphi
    dB_x_dR(2) = -B_R_vec(1) * sin_phi_tok - B_R_vec(2) * cos_phi_tok
    ! dB_x/dz
    dB_x_dR(3) = dB_r_inter%dz * cos_phi_tok - dB_t_inter%dz * sin_phi_tok
    ! dB_y/dR - no ripple yet, hence this term is simpler
    dB_y_dR(1) = dB_r_inter%dR * sin_phi_tok + dB_t_inter%dR * cos_phi_tok
    ! dB_y/dphi - no ripple yet, hence this term is simpler
    dB_y_dR(2) = B_R_vec(1) * cos_phi_tok - B_R_vec(2) * sin_phi_tok
    ! dB_y/dz
    dB_y_dR(3) = dB_r_inter%dz * sin_phi_tok + dB_t_inter%dz * cos_phi_tok
    ! dB_z/dR
    dB_z_dR(1) = dB_z_inter%dR
    ! dB_z/dphi
    dB_z_dR(2) = 0.d0
    ! dB_z/dz
    dB_z_dR(3) = dB_z_inter%dz
    ! Now with the chain rule compute the entries of dB_i/dx_j
    ! dBx/dx third term is zero since dz/dx = 0
    spatial_grad_B(1,1) = dB_x_dR(1) * func_dR_dx(x_vec(1), x_vec(2)) + dB_x_dR(2) * func_dphi_dx(x_vec(1), x_vec(2))
    ! dBx/dy third term is zero since dz/dx = 0
    spatial_grad_B(1,2) = dB_x_dR(1) * func_dR_dy(x_vec(1), x_vec(2)) + dB_x_dR(2) * func_dphi_dy(x_vec(1), x_vec(2))
    ! dBx/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(1,3) = dB_x_dR(3)
    ! dBy/dx third term is zero since dz/dy = 0
    spatial_grad_B(2,1) = dB_y_dR(1) * func_dR_dx(x_vec(1), x_vec(2)) + dB_y_dR(2) * func_dphi_dx(x_vec(1), x_vec(2))
    ! dBy/dy third term is zero since dz/dy = 0
    spatial_grad_B(2,2) = dB_y_dR(1) * func_dR_dy(x_vec(1), x_vec(2)) + dB_y_dR(2) * func_dphi_dy(x_vec(1), x_vec(2))
    ! dBy/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(2,3) = dB_y_dR(3)
    ! dBz/dx third term is zero since dz/dx = 0
    spatial_grad_B(3,1) = dB_z_dR(1) * func_dR_dx(x_vec(1), x_vec(2))
    ! dBz/dy third term is zero since dz/dy = 0
    spatial_grad_B(3,2) = dB_z_dR(1) * func_dR_dy(x_vec(1), x_vec(2))
    ! dBz/dz first and second term is zero since dR/dz = dphi/dz = 0
    spatial_grad_B(3,3) = dB_z_dR(3)
!    print* , "x_vec", x_vec
!    print* , "R_vec", R_vec
    if(debug_level > 2) then
      do i = 1,3
         do j = 1 , 4
            aux_x(j,:) = x_vec
            if(j < 3) aux_x(j,i) = x_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
            if(j >= 3) aux_x(j,i) = x_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
            aux_R(j,:) = R_vec
            if(j < 3) aux_R(j,i) = R_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
            if(j >= 3) aux_R(j,i) = R_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
         end do
         print*,"R",R_vec
         print*,"x ana",spatial_grad_B(1,i)
         print*, "x num",(- func_B_x_ana(aux_x(1,:)) + &
                                    8.d0 *  func_B_x_ana(aux_x(2,:)) - &
                                8.d0 *  func_B_x_ana(aux_x(3,:))  &
                                +   func_B_x_ana(aux_x(4,:))) / (12.d0 *h_x)
         print*,"y ana",spatial_grad_B(2,i)
         print*, "y num",(- func_B_y_ana(aux_x(1,:)) + &
                                    8.d0 *  func_B_y_ana(aux_x(2,:)) - &
                                8.d0 *  func_B_y_ana(aux_x(3,:))  &
                                +   func_B_y_ana(aux_x(4,:))) / (12.d0 *h_x)
        print*,"z ana",spatial_grad_B(3,i)
        print*,"should be zero"
        print*,"dBx/dR ana",dB_x_dR(i)
        print*,"dBx/dR num",(- func_B_x_ana_R(aux_R(1,:)) + &
                                    8.d0 *  func_B_x_ana_R(aux_R(2,:)) - &
                                8.d0 *  func_B_x_ana_R(aux_R(3,:))  &
                                +   func_B_x_ana_R(aux_R(4,:))) / (12.d0 *h_x)
        print*,"dBy/dR ana",dB_y_dR(i)
        print*,"dBy/dR num",(- func_B_y_ana_R(aux_R(1,:)) + &
                                    8.d0 *  func_B_y_ana_R(aux_R(2,:)) - &
                                8.d0 *  func_B_y_ana_R(aux_R(3,:))  &
                                +   func_B_y_ana_R(aux_R(4,:))) / (12.d0 *h_x)
        print*,"dBz/dR ana - should be zero", dB_z_dR(i)
      end do
    end if
!    stop "correct?"
    !---gradient of total magnetic field
    B_abs = sqrt(B_x_vec(1)**2 + B_x_vec(2)**2 + B_x_vec(3)**2)
    ! to obtain d|B|/dx it is split into d|B|/dB_x_j and dB_x_i/dx_j
    ! -> d|B|/dx_i = B_x_j / |B| * dB_j/dx_i
    spatial_grad_B_abs(:) = 0.d0
    do i = 1, 3
        do j = 1,3
          spatial_grad_B_abs(i) = spatial_grad_B_abs(i) + B_x_vec(j) * spatial_grad_B(j,i)
        end do
    end do
    spatial_grad_B_abs = spatial_grad_B_abs / B_abs
    !B_vec_norm = B_x_vec / B_abs
    !N_vec_norm = N_vec / N_abs
    scal_prod = 0.d0
    do i =1, 3
      scal_prod = scal_prod + B_x_vec(i) * N_vec(i)
    end do
    N_par = scal_prod / B_abs
    !------------------------------ Now that all mangetic field gradients are known compute dcos_theta/dx_i-!
    !------------------------------ and dcos_theta/dN_i ----------------------------------------------------!
    ! N_grad_N_par = d/dN Sum_i (N_i * B_i)/ (B_abs) = vec(B)/B(abs)
    !
    N_grad_N_par =  B_x_vec(:) / B_abs
    !----------------------------- Finally do the same thing for the spatial gardiend of N_par -------------------------!
    !----------------------------- Here we resolve dN_par/dx into dN_par/dB_j * dB_j/dx_i  (multidimensional chain rule)!
    spatial_grad_N_par(:) = 0.d0
    do i = 1,3
      do j =1,3
        do k = 1,3
        spatial_grad_N_par(i)  = spatial_grad_N_par(i) + (delta(j,k) * (N_vec(j) * B_abs**2 - B_x_vec(j)**2 * N_vec(j)) - &
                                 (1.d0 - delta(j,k)) * (N_vec(k) * B_x_vec(j) * B_x_vec(k))) * spatial_grad_B(j,i)
        end do
      end do
    end do
    spatial_grad_N_par(:) = spatial_grad_N_par / (B_abs)**3
  end subroutine sub_grad_N_par_ana

  function func_rhop_ana(plasma_params, x_vec)
    USE f90_kind
    USE mod_ecfm_refr_types , only : plasma_params_type, h_x_glob
    USE ripple3d,                 only: grad_type
    use mod_ecfm_refr_utils,      only: sub_remap_coords
    type(plasma_params_type)                  :: plasma_params
    real(rkind), dimension(:)  , intent(in)   :: x_vec
    real(rkind), dimension(3)                 :: R_vec
    real(rkind)                               :: func_rhop_ana
    real(rkind)                               :: ddx, ddy, dummy
    call sub_remap_coords(x_vec, R_vec)
    func_rhop_ana = (3.d0*Sqrt((-1.5d0 + R_vec(1))**2 + R_vec(3)**2/10.d0))/2.d0
  end function func_rhop_ana

    subroutine sub_spatial_grad_N_par_ana(plasma_params, x_vec, N_vec, N_abs, N_par, spatial_grad_N_par, spatial_grad_B_abs, B_abs, N_grad_N_par)
  ! dN_par / dx
      USE f90_kind
      USE mod_ecfm_refr_types, only : plasma_params_type
      implicit none
      type(plasma_params_type), intent(in)        :: plasma_params
      real(rkind), dimension(:),      intent(in)  :: x_vec, N_vec
      real(rkind),                    intent(in)  :: N_abs
      real(rkind),                    intent(out) :: N_par, B_abs
      real(rkind), dimension(:),      intent(out) :: spatial_grad_N_par
      real(rkind), dimension(:),      intent(out) :: spatial_grad_B_abs
      real(rkind), dimension(:),      intent(out) :: N_grad_N_par
      real(rkind), dimension(3)                   :: R_vec
      call sub_grad_N_par_ana(plasma_params, x_vec, N_vec, N_abs, B_abs, spatial_grad_B_abs, &
                             N_par, N_grad_N_par, spatial_grad_N_par)
   end subroutine sub_spatial_grad_N_par_ana

   subroutine sub_spatial_grad_Y_ana(plasma_params, omega, x_vec, B_abs, spatial_grad_B_abs, Y, spatial_grad_Y_ana)
  ! dY^2 / dx
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma
    use constants,            only: pi, e0, mass_e, eps0, c0
    ! corresponds to flux coordinates
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind),                 intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind),                 intent(in)   :: B_abs
    real(rkind), dimension(:),   intent(in)   :: spatial_grad_B_abs
    real(rkind),                 intent(out)  :: Y
    real(rkind), dimension(:),   intent(out)  :: spatial_grad_Y_ana
    integer(ikind)                            :: i
    real(rkind), dimension(3)                 :: B2_grad, spatial_grad_rhop ! grad_x(B^2)
    real(rkind)                               :: rhop, T_e, grad_T_e
    if(warm_plasma) then
      call sub_spatial_grad_rhop_ana(plasma_params, x_vec, rhop, spatial_grad_rhop)
      call sub_grad_T_e_ana(plasma_params, rhop, T_e, grad_T_e)
      spatial_grad_rhop = spatial_grad_rhop * plasma_params%rhop_scale_Te
      spatial_grad_Y_ana(:) =   (e0*((0.2*c0**2*mass_e + e0*T_e) * spatial_grad_B_abs(:) - &
        0.5*e0*B_abs*spatial_grad_rhop*grad_T_e)) / &
        (mass_e*omega*(0.2*c0**2*mass_e + e0*T_e) * &
        Sqrt(1. + (5.*e0*T_e)/(c0**2*mass_e)))
    else
      spatial_grad_Y_ana(:) = spatial_grad_B_abs(:) * e0 / mass_e
      T_e = 0.d0
    end if
    Y = func_Y(plasma_params, omega, B_abs, T_e)
  end subroutine sub_spatial_grad_Y_ana

    subroutine sub_spatial_grad_X_ana(plasma_params, omega, x_vec, X, spatial_grad_X_ana, rhop_out)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, warm_plasma!, h_x_glob
    USE ripple3d,                 only: grad_type
    use constants,                 only : pi, e0, mass_e, eps0, c0
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind)              ,   intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec
    real(rkind)              ,   intent(out)  :: X
    real(rkind), dimension(:),   intent(out)  :: spatial_grad_X_ana
    real(rkind)              ,   intent(out)  :: rhop_out
    real(rkind)                               :: rhop, n_e, grad_n_e, T_e, grad_T_e
    real(rkind), dimension(3)                 :: spatial_grad_rhop_ne, spatial_grad_rhop_Te
    real(rkind), dimension(3)                 :: R_vec
    call sub_spatial_grad_rhop_ana(plasma_params, x_vec, rhop, spatial_grad_rhop_ne)
    call sub_grad_n_e_ana(plasma_params, rhop, n_e, grad_n_e)
    spatial_grad_rhop_Te = spatial_grad_rhop_ne * plasma_params%rhop_scale_Te
    spatial_grad_rhop_ne = spatial_grad_rhop_ne * plasma_params%rhop_scale_ne
    if(warm_plasma) then
      spatial_grad_rhop_Te = spatial_grad_rhop_ne * plasma_params%rhop_scale_Te
      call sub_grad_T_e_ana(plasma_params, rhop, T_e, grad_T_e)
      spatial_grad_X_ana = (e0**2*(2.d0*(c0**2 * mass_e + 5.d0 * e0 * T_e) * &
                       grad_n_e * spatial_grad_rhop_ne - &
                       5.d0 * e0 * n_e * grad_T_e * spatial_grad_rhop_Te)) / &
                       (2.d0*eps0*mass_e * omega**2 * (c0**2*mass_e + 5.d0*e0*T_e) * &
                       Sqrt(1.d0 + (5.d0 * e0*T_e)/(c0**2*mass_e)))
    else
      T_e = 0.d0
      spatial_grad_X_ana = spatial_grad_rhop_ne(:) * grad_n_e * e0**2.d0/(eps0 * mass_e * omega**2)
    end if
    X  =  func_X(plasma_params, omega, n_e, T_e)
    rhop_out = rhop
  end subroutine sub_spatial_grad_X_ana

   subroutine sub_single_step_explicit_RK4(plasma_params, omega, x_vec, N_vec, h, x_vec_out, N_vec_out, B_vec_out, theta_out, Hamil_out, N_s_out, n_e_out, omega_c_out, rhop_out)
    USE f90_kind
    USE mod_ecfm_refr_types, only : plasma_params_type, Hamil
    USE constants,                 only : eps0, mass_e, e0, c0
    implicit None
    type(plasma_params_type),    intent(in)   :: plasma_params
    real(rkind),                 intent(in)   :: omega
    real(rkind), dimension(:),   intent(in)   :: x_vec, N_vec
    real(rkind), dimension(:),   intent(out)  :: B_vec_out
    real(rkind),                 intent(in)   :: h
    real(rkind), dimension(:),   intent(out)  :: x_vec_out, N_vec_out
    real(rkind),                 intent(out)  :: theta_out, Hamil_out, N_s_out, n_e_out, omega_c_out, rhop_out
    real(rkind), dimension(0:4,6)             :: k
    real(rkind), dimension(4)                 :: h_arr
    real(rkind), dimension(6)                 :: y_vec
    real(rkind)                               :: N_abs, k_abs, omega_c, omega_p2, T_e_out
    real(rkind)                               :: X, Y, T_e,  A, B, C, N_par
    integer(ikind)                            :: i
    stop "RUNGE KUTTA 4 is insufficient for the ray tracing problem"
    N_abs = 0.d0
    do i = 1,3
      N_abs = N_abs + N_vec(i)**2
    end do
    N_abs = sqrt(N_abs)
    call sub_local_params(glob_plasma_params, glob_omega, x_vec, N_vec, B_vec_out, N_abs, n_e_out, omega_c_out, T_e_out, theta_out, rhop_out)
    X = func_X(glob_plasma_params,glob_omega, n_e_out,T_e_out)
    Y = func_Y(glob_plasma_params,glob_omega, omega_c * mass_e / e0, T_e_out)
    N_par = cos(theta_out) * N_abs
    h_arr(1) = 0.d0
    h_arr(2) = 0.5d0 * h
    h_arr(3) = 0.5d0 * h
    h_arr(4) = h
    y_vec(1:3) = x_vec
    y_vec(4:6) = N_vec
    k(0,:) = 0.d0
    if(Hamil == "Dani") then
      Hamil_out = func_Lambda_star(N_abs, X, Y, N_par, glob_mode)
      N_s_out = N_abs
      do i = 1, 4
        call f_Lambda(6, 0.d0, y_vec + h_arr(i) * k(i - 1,:), k(i,:))
      end do
    else
      N_s_out = N_abs
      A = func_A(X, Y)
      B = func_B(X, Y, N_par)
      C = func_C(X, Y, N_par)
      Hamil_out = func_H(sqrt(N_abs**2 - N_par**2),A, B, C, glob_mode)
      do i = 1, 4
        call f_H(6, 0.d0, y_vec + h_arr(i) * k(i - 1,:), k(i,:))
      end do
    end if
    x_vec_out = x_vec + h / 6.d0 * (k(1,1:3)  + 2.d0 * k(2,1:3) + 2.d0 * k(3,1:3) + k(4,1:3))
    N_vec_out = N_vec + h / 6.d0 * (k(1,4:6)  + 2.d0 * k(2,4:6) + 2.d0 * k(3,4:6) + k(4,4:6))
  end subroutine sub_single_step_explicit_RK4

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

!*******************************************************************************
subroutine import_all_ece_data()
USE mod_ecfm_refr_types, only: ant, rad, data_folder, dstf, N_freq, plasma_params, ffp
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
  USE nag_spline_1d,                only: nag_spline_1d_interp
#endif
!use constants, only: mass_e, e0
implicit none
integer(ikind)  :: idiag
Integer(ikind)  :: ich,ifreq,iint, imode
real(rkind)     :: dfreq
Character(200)  :: cur_filename, filename
character(12)   :: ich_str
character(1)    :: sep
idiag = 1
imode = 1
filename = trim(data_folder) // "cnt.dat"
open(67,file=trim(filename))
!filename = trim(data_folder) // "rhopres.dat"
!open(74,file=trim(filename))
!filename = trim(data_folder) // "sres.dat"
!open(76,file=trim(filename))
filename = trim(data_folder) // "f_ECE.dat"
open(79,file=trim(filename))
filename = trim(data_folder) // "diag.dat"
open(84,file=trim(filename))
print*, "Importing data for a total of",  ant%diag(idiag)%N_ch, "ECE channels"
read(84,"(A3)") ant%diag(idiag)%diag_name
filename = trim(data_folder) // "parms.dat"
open(80,file=trim(filename))
do ich = 1, ant%diag(idiag)%N_ch
  read(67,"(I5)") rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points
  do ifreq = 1, N_freq
    if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec)) &
       deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec)
    allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
    write(ich_str, "(I3.3)") ich
    cur_filename = trim(data_folder) // "chdata" // trim(ich_str) // ".dat"
    open(66, file=cur_filename)
    do iint = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points
      read(66,"(9(E18.10E3A1))") rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%s,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%R,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%z,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%rhop,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%ne,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%Te,sep,&
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%theta,sep, &
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%v_g_perp, sep, &
         rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(ifreq)%svec(iint)%freq_2x,sep!, &
    end do
    close(66)
  end do
  close(89)
  !read(74,"(E18.10E3)") rad%diag(idiag)%ch(ich)%rhop_res
  !read(76,"(E18.10E3)") rad%diag(idiag)%ch(ich)%s_res
  read(79,"(E18.10E3)") ant%diag(1)%ch(ich)%f_ECE
  ant%diag(1)%ch(ich)%freq(1) = ant%diag(1)%ch(ich)%f_ECE
end do
read(80,"(E18.10E3)") plasma_params%B_ax
close(67)
!close(74)
!close(76)
close(79)
close(80)
close(84)
if(dstf == "numeric") then
  filename = trim(data_folder) // "B_min.dat"
  open(90, file=filename)
  do iint = 1, size(ffp%rhop_B_min)
    read(90, fmt="(E18.10E3A1E18.10E3)") ffp%rhop_B_min(iint), sep, ffp%B_min(iint)
  end do
  close(90)
  call make_1d_spline(ffp%B_min_spl, int(size(ffp%rhop_B_min), kind=4), ffp%rhop_B_min, ffp%B_min)
#ifdef NAG
  call nag_spline_1d_interp(ffp%rhop_B_min, ffp%B_min, ffp%B_min_nag_spl)
#endif
end if
end subroutine import_all_ece_data

subroutine make_ecfm_LOS_grid(flag, sparse_step, dense_step) ! S. Denk 4. 2013
! it = current index in the time vector
! flag = "init" or "terminate " for either allocation or deallocation
! This subroutine initializes the grid for the integration along los
! The parameters that remain constant within the optimization are saved in the
! structure rad_time_ch_ray_svec_type declared in the global params section.
! This structure is memory intensive, but can replace most calls to interpol_LOS.
! S.Denk March 2013
use mod_ecfm_refr_types,                    only: ant,rad,  non_maxwellian, N_ray, N_freq, modes, output_level, &
                                                  dstf, bi_max, data_folder, pnts_BPD, &
                                                  largest_svec
use constants,                              only: e0, c0
!use ecfm_non_therm_abs,                     only: abs_non_therm_init,abs_non_therm_clean_up
implicit none
character(10), intent(in)                   :: flag
integer(ikind), intent(in), optional        :: sparse_step, dense_step
integer(ikind)                              :: idiag
integer(ikind)                              :: ich, imode,ir, iint, ifreq
real(rkind)                                 :: dense_region_l,&
                                               dense_region_u
real(rkind)                                 :: ds1, ds2,s_ow, s_iw
idiag = 1
imode = 1
if(ant%N_diag > 1) then
  stop "Import of svec not supported for multiple diagnostic"
end if
if(modes /= 1) then
  print*, "Loading rays is not supported for O-mode"
  print*, "Input Error in initialize_LOS.f90"
  call abort
end if
if(flag == "initialize") then
  largest_svec = 0
  do ich = 1, ant%diag(idiag)%N_ch
    do ir = 1, N_ray
      do ifreq = 1, N_freq
        ! FIXME: Inconsistent with forward model by Rathgeber
        if(largest_svec < rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points) then
          largest_svec = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
        end if
        do iint = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%Ibb = &
          ant%diag(idiag)%ch(ich)%freq(ifreq)**2 * e0 * &
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%Te / c0**2
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%sin_theta =  &
              sin(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%theta)
          rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%cos_theta = &
              cos(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%theta)
          if(dstf == "numeric" .and. iint < rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points &
            .and. iint > 1) then
            if(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%rhop < &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint + 1)%rhop &
            .and. rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint - 1)%rhop > &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%rhop) &
            rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%s_axis = rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(iint)%s
          end if
        enddo !iint = 1, rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points
        rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(1:rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%total_LOS_points)%plasma = .true.
        where(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:)%rhop == -1.d0) \
              rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec(:)%plasma = .false.
      enddo ! ifreq = 1, ant%diag(idiag)%ch(ich)%N_freq
      if(output_level) then
        allocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%Trad_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%em_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%ab_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points), &
                 rad%diag(idiag)%ch(ich)%mode(imode)%ray_extra_output(ir)%T_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(1)%total_LOS_points))
      end if
    enddo  ! ir = 1, ant%diag(idiag)%ch(ich)%N_ray
    if(output_level) then
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%s(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%R(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%z(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Trad_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%em_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%ab_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%T_secondary(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%Te(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cold(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_cor(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%N_warm(rad%diag(idiag)%ch(ich)%mode(imode)%ray(1)%freq(1)%total_LOS_points))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%rhop_BPD(pnts_BPD))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD(pnts_BPD))
      allocate(rad%diag(idiag)%ch(ich)%mode_extra_output(imode)%BPD_secondary(pnts_BPD))
    end if
  enddo ! ich = 1, ant%diag(idiag)%N_ch
  rad%diag(idiag)%ch(:)%eval_ch = .true. ! stand_alone => evaluate all channels
  !call abs_non_therm_init("f_nor")
else if(flag == "terminate ") then
  do idiag = 1, ant%N_diag
    do ich = 1, ant%diag(idiag)%N_ch
      do ir = 1, N_ray
        do ifreq = 1, N_freq
        ! FIXME: lots of stuff still allocated after this - clean up!
          if(allocated(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec)) then
            deallocate(rad%diag(idiag)%ch(ich)%mode(imode)%ray(ir)%freq(ifreq)%svec)
          end if
        end do ! ifreq = 1, ant%diag(idiag)%ch(ich)%N_freq
      end do ! ir = 1, ant%diag(idiag)%ch(ich)%N_ray
    end do !ich = 1, ant%diag(idiag)%N_ch
  end do
  if(dstf == "numeric") call ffp_clean_up()
  if(trim(dstf) == "gene") call fgene_clean_up()
else
  print*, "Flag has to be init or terminate "
  print*, "It was set to: ", flag
  print*, "Wrong usage of make_ecfm_LOS_grid"
  call abort
end if
end subroutine make_ecfm_LOS_grid

subroutine read_wall_Trad()
use mod_ecfm_refr_types,       only: reflec_equ, data_folder, modes
use mod_ecfm_refr_interpol,    only: make_1d_spline
#ifdef NAG
USE nag_spline_1d,             only: nag_spline_1d_interp
#endif
Character(200)  :: cur_filename
Character(1)    :: blanc
integer(ikind)              :: i
  if(modes == 1 .or. modes == 3) then
    cur_filename = trim(data_folder) // "X_reflec_Trad.dat"
    open(67, file = trim(cur_filename))
    read(67,"(I5.5)") reflec_equ%N_f
    allocate(reflec_equ%f(reflec_equ%N_f), reflec_equ%X_Trad_equ(reflec_equ%N_f))
    do i = 1, reflec_equ%N_f
      read(67,"(E19.12E2A1E19.12E2)") reflec_equ%f(i), blanc, reflec_equ% X_Trad_equ(i)
    end do
    close(67)
    call make_1d_spline(reflec_equ%X_Trad_equ_spl, int(reflec_equ%N_f,4), reflec_equ%f, reflec_equ%X_Trad_equ)
#ifdef NAG
    call nag_spline_1d_interp(reflec_equ%f, reflec_equ%X_Trad_equ, reflec_equ%X_Trad_equ_spl_nag)
#endif
  end if
  if(modes == 2 .or. modes == 3) then
    cur_filename = trim(data_folder) // "O_reflec_Trad.dat"
    open(67, file = trim(cur_filename))
    read(67,"(I5.5)") reflec_equ%N_f
    if(.not. allocated(reflec_equ%f)) allocate(reflec_equ%f(reflec_equ%N_f))
    allocate(reflec_equ%O_Trad_equ(reflec_equ%N_f))
    do i = 1, reflec_equ%N_f
      read(67,"(E19.12E2A1E19.12E2)") reflec_equ%f(i), blanc, reflec_equ%O_Trad_equ(i)
    end do
    close(67)
    call make_1d_spline(reflec_equ%O_Trad_equ_spl, int(reflec_equ%N_f,4), reflec_equ%f, reflec_equ%O_Trad_equ)
#ifdef NAG
    call nag_spline_1d_interp(reflec_equ%f, reflec_equ%O_Trad_equ, reflec_equ%O_Trad_equ_spl_nag)
#endif
  end if
    reflec_equ%f_min = minval(reflec_equ%f)
    reflec_equ%f_max = maxval(reflec_equ%f)
end subroutine read_wall_Trad


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
! mod_ecf_refr_types

type reflec_equ_type
! Parameters for wall plasma equilibrium reflection model - reflec_model = 1
real(rkind), dimension(:), allocatable  :: f, X_Trad_equ, O_Trad_equ
type(spl_type_1d)                       :: X_Trad_equ_spl, O_Trad_equ_spl
real(rkind)                             :: f_min, f_max
integer(ikind)                          :: N_f
#ifdef NAG
  type(nag_spline_1d_comm_wp)            :: X_Trad_equ_spl_nag, O_Trad_equ_spl_nag
#endif
end type

!
!mod_ecfm_refr_rad_transp
subroutine benchmark_abs_and_N()
! produces profiles of both the absorption coefficient and the refractive index
use mod_ecfm_refr_types,         only: rad_diag_ch_mode_ray_freq_svec_type, straight
use f90_kind
use mod_ecfm_refr_abs_Al,             only: abs_Albajar, abs_Al_tor_abs, func_N_cold, func_rel_N, get_upper_limit_tau
use mod_ecfm_refr_em_Hu,                    only: calculate_em
use constants,                    only: pi, e0, mass_e, eps0, c0
#ifdef OMP
use omp_lib
#endif
implicit none
type(rad_diag_ch_mode_ray_freq_svec_type) :: svec
real(rkind)                               :: omega
integer(ikind)                            :: mode
real(rkind)                               :: c_abs_Alb, c_abs_warm_disp, c_abs_Hutch
real(rkind)                               :: N_cold, N_cor, N_gray
real(rkind), dimension(600)               :: freq_2X_prof
integer(ikind)                            :: i, m
real(rkind)                               :: ds, dummy_1, dummy_2, abs_crude_approx!, X, Y,
logical                                   :: warm_plasma
!Te_prof(1) = 200.d0
!Te_prof(200) = 25.d3
!do i  = 2, 199
!  Te_prof(i) = Te_prof(1) + (Te_prof(200) - Te_prof(1)) / 200.d0 * i
!end do
!ne_prof(1) = 0.5d19
!ne_prof(200) = 2.d20
!do i  = 2, 199
!  ne_prof(i) = ne_prof(1) + (ne_prof(200) - ne_prof(1)) / 200.d0 * i
!end do
straight = .true.
m = size(freq_2X_prof)
freq_2X_prof(1) = 65.d9
freq_2X_prof(m) = 220.d9
do i  = 2, m -2
  freq_2X_prof(i) = freq_2X_prof(1) + (freq_2X_prof(m) - freq_2X_prof(1)) / real(m,8) * i
end do
mode = 1
svec%theta = 85.d0/180.d0 * Pi
svec%cos_theta = cos(svec%theta)
svec%sin_theta = sin(svec%theta)
!svec%ne = 1.d19
svec%Te = 8.e3
omega = 140.d9 * 2.d0 * Pi
!X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
!Y = freq_2X_prof(1) * Pi / omega
warm_plasma = .true.
ds = 1.d0
svec%v_g_perp = svec%sin_theta
!print*, "X, Y", X, Y
svec%ne = 10.d19
open(81, file = "v_high_ne_abs_prof.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
svec%ne = 8.d19
open(81, file = "high_ne_abs_prof.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
open(81,  file = "low_ne_abs_prof.dat")
svec%ne = 6.d19
do i = 1, m
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  svec%freq_2X = freq_2X_prof(i)
  !X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  !print*, "X, Y", X, Y
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
omega = 110.d9 * 2.d0 * Pi
svec%ne = 6.d19
open(81, file = "v_high_ne_abs_prof_low_f.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
svec%ne = 4.d19
open(81, file = "high_ne_abs_prof_low_f.dat")
do i = 1, m
  svec%freq_2X = freq_2X_prof(i)
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cor
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp, " ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
open(81,  file = "low_ne_abs_prof_low_f.dat")
svec%ne = 2.d19
do i = 1, m
  !omega = (100.d9 + 50.d9 / 200.d0 * i) * 2.d0 * pi
  svec%freq_2X = freq_2X_prof(i)
  !X =  svec%ne * e0**2 / (eps0 * mass_e) / omega**2
  !Y = svec%freq_2X * Pi / omega
  N_cor = func_rel_N(omega, svec, mode)
  svec%N_cold = N_cor
  N_cold = func_N_cold(omega, svec, mode)
  N_gray = N_cold
  !print*, "X, Y", X, Y
  call abs_Albajar(svec, omega, mode, ds, c_abs_Alb, dummy_1)
  call calculate_em(svec, omega, dummy_1, dummy_2, c_abs_Hutch)
  c_abs_warm_disp = abs_Al_tor_abs(svec, omega, mode, N_gray)
  abs_crude_approx =  get_upper_limit_tau(svec, omega, 1.d0, 0.5d0)
  write(81,"(7(E18.10E3,A1),E18.10E3)") &
          omega / (svec%freq_2X * Pi), " ", c_abs_Alb, " ", c_abs_Hutch, " ",&
          c_abs_warm_disp," ", abs_crude_approx, " ", N_cold, " ", N_cor, " ", N_gray
end do
close(81)
end subroutine benchmark_abs_and_N

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

! mod_ecfm_refr.f90

!#ifdef IDA
!subroutine simulate_ida(working_dir)
!use mod_ecfm_refr_types,        only: n_e_filename, T_e_filename, &
!                                      ant, plasma_params, rad, output_level, &
!                                      Ich_name, data_folder, ray_out_folder, &
!                                      N_ray, N_freq, use_ida_spline_ne, use_ida_spline_Te, &
!                                      data_name, data_secondary_name
!use mod_ecfm_refr_utils,        only: parse_ecfm_settings_from_ida, retrieve_n_e
!#ifdef NAG
!USE nag_spline_1d,              only: nag_spline_1d_interp
!#endif
!use mod_ecfm_refr_interpol,     only: make_1d_spline, deallocate_1d_spline
!implicit none
!character(200), intent(in)    :: working_dir
!type(ece_type)                        :: ece_strut
!real(rkind), dimension(:), allocatable :: rhop, T_e, n_e, dat_model_ece, par, par_ne, par_scal, ece_rhop, ne_test, rho_pol
!real(rkind)                            :: reflec_X, reflec_O, rp_min
!logical,     dimension(:), allocatable :: ECE_fm_flag_ch
!integer(ikind)                          :: i, j, m, n, itime, idiag, ich
!character(1)                            :: sep
!character(250)                          :: filename
!  ida%shot = 32080
!  allocate(ida%time(1))
!  itime = 1
!  ida%time(itime)%time = 2.57d0
!  ida%time(itime)%time_beg =  ida%time(itime)%time - 5.d-4! only used for magnetic on axis
!  ida%time(itime)%time_end =  ida%time(itime)%time + 5.d-4
!  ida%ece%reflec_X_mode = 0.9
!  ida%ece%reflec_O_mode = 0.95
!  ida%ece%rhopol_scal_te = 1.0
!  ida%ece%rhopol_scal_ne = 1.0
!  ida%ece%rhopol_scal_te_fit = 0
!  ida%ece%N_ray = 1
!  ida%ece%N_freq = 1
!  ida%ece%expnam = "AUGD"
!  ida%ece%diag = "CEC"
!  ida%ece%edition = 0
!  ida%diag_eq = "EQH"
!  ida%exp_eq =  "AUGD"
!  ida%ed_eq = 0
!  ida%btf_corr_fact_ext = 1.005d0
!  reflec_X = 0.92d0
!  reflec_O = 0.95d0
!  ! eq diag !
!  n_e_filename = trim(working_dir) // "ecfm_data/" // "ne_file.dat"
!  T_e_filename = trim(working_dir) // "ecfm_data/" // "Te_file.dat"
!  data_name = "TRadM_therm.dat"
!  data_secondary_name = "TRadM_TBeam.dat"
!  Ich_name = "IchTB"
!  data_folder = trim(working_dir) //  "ecfm_data" // "/"
!  ray_out_folder = trim(working_dir) // "ecfm_data/" // "ray/"
!  ! Simulate loaded ece_strut
!  call prepare_ece_strut(ida, ece_strut)
!  ! Make a topfile to communicate equilibrium to ecfm
!  open(66, file = n_e_filename)
!  read(66, "(I7.7)") m
!  allocate(rhop(m), n_e(m), ne_test(m))
!  do i = 1, m
!    read(66,"(E19.12E2A1E19.12E2)") rhop(i), sep, n_e(i)
!  end do
!  close(66)
!  open(66, file = T_e_filename)
!  read(66, "(I7.7)") n
!  if(m /= n) stop "Te and ne have to be given on the same rhop axis"
!  allocate(T_e(m))
!  do i = 1, m
!    read(66,"(E19.12E2A1E19.12E2)") rhop(i), sep, T_e(i)
!  end do
!  close(66)
!  ! do the spline interpolation here, since we want to simulate an IDA run, where Te and ne are provided by
!  ! make_temperature and make_density  - this is only for tests!
!  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
!  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
!#ifdef NAG
!  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
!  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
!#endif
!  ida%rhopol_max_spline_knot = maxval(rhop)
!  rp_min = maxval(rhop)
!  ! initialize_ecfm needs a rework to use the spline structure instead of rhop, T_e, n_e
!  ! ideally both par and par_scale and rhop, T_e, n_e can be used
!  ! output_level =  true collects the information about ray, los, birthplace distribution
!  output_level = .false.
!  call pre_initialize_ecfm(working_dir, ida, "init", ece_strut)
!  call pre_initialize_ecfm(working_dir, ida, "clean")
!  call pre_initialize_ecfm(working_dir, ida, "load")
!  allocate(rho_pol(size(ece_strut%ch)))
!  call initialize_ecfm(ida, itime, "init", rho_pol)
!  print*, "Initialized rho_pol",rho_pol
!  allocate(dat_model_ece(size(ece_strut%ch)), ece_fm_flag_ch(size(ece_strut%ch)), ece_rhop(size(ece_strut%ch)))
!  ece_fm_flag_ch(:) = .false.
!  use_ida_spline_ne = .false.
!  use_ida_spline_Te = .false.
!  allocate(par(30), par_ne(30), par_scal(30)) ! we need them allocated
!  ! First two are hardcoded to be ne and Te scale
!  par(1) = plasma_params%rhop_scale_ne
!  par(2) = plasma_params%rhop_scale_te
!  par_ne = par
!  par_scal(:) = 1.d0
!!  call retrieve_n_e(plasma_params, rhop(1:m-10), ne_test(1:m-10))
!!  print*, "rhop", rhop(1:m-10)
!!  print*, "ne", n_e(1:m-10)
!!  print*, "ne _interp", ne_test(1:m-10)
!!  print*, "Starting optimization"
!!  print*, ida%ece%rhopol_scal_te
!!  print*, ida%ece%rhopol_scal_ne
!  do i =1, 3
!    do j =1, 5
!      call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
!      call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
!#ifdef NAG
!      call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
!      call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
!#endif
!      call make_dat_model_ece_ecfm_refr(par, par_ne, par_scal, reflec_X, reflec_O, & ! in
!                               ece_fm_flag_ch, rp_min, dat_model_ece)
!      print*, j, "-th optimization step of", i, "-th optimization"
!      T_e(:) = T_e(:) * 1.01d0
!      n_e(:) = n_e(:) * 1.01d0
!    end do
!    call make_rays_ecfm(par, par_scal, ece_rhop)
!    ece_fm_flag_ch(:) = .true.
!    ! update_plasma_params and span_svecs needs to be called to reinitialize the LOS
!    ! another routine that should work with both rhop, T_e, n_e and par, par_scal
!  end do
!  print*,"Calculating chi**2"
!!  call initialize_ecfm( ida, itime, "clean")
!!  call pre_initialize_ecfm(working_dir, ida, "clean")
!!  call make_rays_ecfm(par, par_scal, ece_rhop)
!  do i = 1, 3
!    call make_dat_model_ece_ecfm_refr_from_Te(par_ne, par_scal, rhop, T_e, reflec_X, reflec_O, & ! in
!                                              ece_fm_flag_ch, rp_min, dat_model_ece)
!    use_ida_spline_Te = .false. ! is reset everytime
!    T_e(:) = T_e(:) * 1.01d0
!  end do
!  output_level = .true.
!  call initialize_ecfm(ida, itime, "clean")
!  call initialize_ecfm(ida, itime, "init")
!  ! Test if everything is really deallocated by first clean and then init again
!  call initialize_ecfm(ida, itime, "clean")
!  call initialize_ecfm(ida, itime, "init")
!  use_ida_spline_ne = .false.
!  use_ida_spline_Te = .false.
!  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
!  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
!#ifdef NAG
!  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
!  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
!#endif
!  call make_rays_ecfm(par, par_scal, ece_rhop)
!  call initialize_ecfm(ida, itime, "clean")
!  call initialize_ecfm(ida, itime, "init")
!  call make_1d_spline( plasma_params%ne_spline, int(size(rhop),4), rhop, n_e)
!  call make_1d_spline( plasma_params%Te_spline, int(size(rhop),4), rhop, T_e)
!#ifdef NAG
!  call nag_spline_1d_interp(rhop, n_e, plasma_params%ne_spline_nag)
!  call nag_spline_1d_interp(rhop, T_e, plasma_params%Te_spline_nag)
!#endif
!  use_ida_spline_ne = .false.
!  use_ida_spline_Te = .false.
!  call make_rays_ecfm(par, par_scal, ece_rhop)
!  ! Debug information in case of missing shot file
!  print*, ant%N_diag, " diagnostics to model"
!  ! multiple diagnostics possible, but currently only one for IDA
!  print*, ant%diag(1)%N_ch, " channels for first diagnostics"
!  print*, N_freq, " frequencies considered"
!  print*, N_ray, " rays considered"
!  ! makes the radiation temperature
!  ! call make_ece_rad_temp()
!  ! This has to be done in EVERY optimastion step
!  ! Again this routine should also work with par, par_scale
!  ! call update_svecs(rad, rhop, T_e, n_e)
!  ! Make rad_temp again
!  ! call make_ece_rad_temp()
!  ! update_plasma_params and span_svecs needs to be called to reinitialize the LOS
!  ! another routine that should work with both rhop, T_e, n_e and par, par_scal
!  ! call update_plasma_params(plasma_params, rhop, T_e, n_e)
!  ! call span_svecs(plasma_params)
!  ! Make rad temp again this time with output_level = True
!  ! Very slow - should never be used within optimization
!  call make_dat_model_ece_ecfm_refr(par, par_ne, par_scal, reflec_X, reflec_O, & ! in
!                               ece_fm_flag_ch, rp_min, dat_model_ece)
!  ! Deallocate everything
!  filename = trim(data_folder) // data_name
!  open(66, file=filename)
!  filename = trim(data_folder) // data_secondary_name
!  open(67, file=filename)
!  filename = trim(data_folder) // "sres_rel.dat"
!  open(68, file=filename)
!  filename = trim(data_folder) // "sres.dat"
!  open(69, file=filename)
!  ! write(filename, "(A64,A7)") data_folder, "tau.dat"
!  ! open(67, file=filename)
!  ! open(67, file="TRadM_s.dat")
!  do idiag = 1, ant%N_diag
!    do ich = 1,ant%diag(idiag)%N_ch
!      write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
!      rad%diag(idiag)%ch(ich)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau
!      write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
!        rad%diag(idiag)%ch(ich)%TRad_secondary/ 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau_secondary
!      write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") &
!                            rad%diag(idiag)%ch(ich)%rel_s_res, " ", rad%diag(idiag)%ch(ich)%rel_R_res, " ",&
!                            rad%diag(idiag)%ch(ich)%rel_z_res, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res, " ", &
!                            rad%diag(idiag)%ch(ich)%rel_s_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_R_res_secondary, " ",&
!                            rad%diag(idiag)%ch(ich)%rel_z_res_secondary, " ", rad%diag(idiag)%ch(ich)%rel_rhop_res_secondary
!      write(69,"(E16.8E3,A1,E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%s_res, " ", rad%diag(idiag)%ch(ich)%R_res, " ",&
!                            rad%diag(idiag)%ch(ich)%z_res, " ", rad%diag(idiag)%ch(ich)%rhop_res
!    end do
!  end do
!  close(66)
!  close(67)
!  close(68)
!  close(69)
!  call initialize_ecfm(ida, itime, "clean")
!  call pre_initialize_ecfm(working_dir, ida, "clean")
!  deallocate(rhop, n_e, T_e, par, par_ne, par_scal, ne_test)
!end subroutine simulate_ida
!#endif

! module radiation
!        subroutine calculate_abs
!        subroutine calculate_em_tot
!        subroutine calculate_mode_frac
!        subroutine calculate_sf
!        subroutine calculate_sf_rrf
!        subroutine calculate_int_beta_nume
!        subroutine radiation_df_over_f


!******************************************************************************
!******************************************************************************
!******************************************************************************
! Emissivity and absorption according to hutchinson with and without some tweaks.
! Since albajar is generally better, just a bit slower -> Deprecation
module mod_ecfm_refr_em_Hu

  use f90_kind

  implicit none
  public :: calculate_em, &
            radiation_gauss_init, &
            radiation_gauss_clean_up!, &
            !calculate_N
    real(rkind), dimension(:), allocatable :: Int_weights
    real(rkind), dimension(:), allocatable :: Int_absz
    real(rkind), dimension(:), allocatable :: Int_weights_many
    real(rkind), dimension(:), allocatable :: Int_absz_many
    integer(ikind)                         :: total_cnt, shortcut_cnt
  private :: BesselJ, &
             BesselJ_fast_1, &
             BesselJ_fast_2, &
             BesselJ_fast_3, & ! different syntax due to use of spline data
             calculate_em_tot,        & !
             calculate_mode_frac,     & !
             calculate_sf_rrf,        & ! analytical calculation from rrf
             calculate_int_u_gauss, &
             calculate_int_u_gauss_many, &
             calculate_int_u_gauss_arbitrary, &
             bessel_term, &
             Int_weights, &
             Int_absz, &
             Int_weights_many, &
             Int_absz_many, &
             total_cnt, &
             shortcut_cnt

contains

!  subroutine calculate_N(svec, omega, m, N, mode, old_model)
!    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
!    use constants,                  only: pi, e0 , eps0, mass_e
!    implicit None
!    type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout)  :: svec
!    real(rkind), intent(in)                           :: omega
!    integer(ikind), intent(in)                        :: m, mode
!    real(rkind), intent(out)                          :: N
!    logical, intent(in), optional                     :: old_model
!    real(rkind)                   :: omega_bar_2, alpha, beta,omega_p_c_2, omega_c, omega_p
!    real(rkind)                   :: f,rho
!    real(rkind)                   :: omega_LC, omega_RC
!    logical                       :: old_model_flag
!    old_model_flag = .False.
!    if(present(old_model)) old_model_flag = old_model
!    omega_c =  svec%freq_2X * PI
!    omega_p     = sqrt((svec%ne * e0**2.d0)/(eps0 * mass_e))
!    omega_bar_2 = (omega / omega_c)**2
!    alpha = omega_p**2 / omega**2
!    beta = omega_c**2 / omega**2
!    if(int(abs(mode)) /= 1) then
!      print*, "Only -1 (X-mode) and 1 (O-Mode) allowed for variable mode in"
!      print*, "the routine abs_non_therm_N"
!      stop "abs_non_therm_N in mod_ecfm_abs_Albajar.f90"
!    end if
!    omega_RC = sqrt(0.25 * omega_c**2 + omega_p**2)
!    omega_LC = omega_RC - 0.5 * omega_c
!    omega_RC = omega_RC + 0.5 * omega_c
!    omega_p_c_2 = omega_p / omega_c
!    if(alpha > 1.0) then
!      N = 0.0
!    else if(mode == 1) then
!      if((omega**2 - omega_LC**2)*(omega**2 - omega_RC**2) < 0.0) then
!        N = 0.0
!        return
!      end if
!    end if
!    rho =  svec%sin_theta**4 + 4.d0/real(m,8)**2 *(real(m,8)**2 - omega_p_c_2)**2 * svec%cos_theta**2
!    if(rho < 0.d0) then
!      print*, "root of negative number avoided: variable rho"
!      stop "calculate_N in mod_ecfm_radiation.f90"
!    end if
!    rho = sqrt(rho)
!    f = (2.d0 * (reaL(m,8)**2 - omega_p_c_2))/ (2.d0 * (real(m,8)**2 - omega_p_c_2) - (svec%sin_theta**2 + real(mode,8) * rho))
!    N = 1.d0 - omega_p_c_2/real(m,8)**2 * f
!    if(N < 0.d0) then
!      N = 0.0
!      return
!    end if
!    N = sqrt(N)
!    if(old_model_flag) N =  simple_in_cutoff(N,omega / (2 * Pi), svec%ne)
!  end subroutine calculate_N

  function simple_in_cutoff(N_abs,freq, ne)
  use mod_ecfm_refr_types,        only: old_cutoff
  implicit none
  real(rkind), intent(in)       :: N_abs, freq, ne
  logical                       :: simple_in_cutoff
  if(old_cutoff) then
    simple_in_cutoff = (0.620221303d-2 * freq**2 < ne)
    !print*, "cutoff", .620221303d-2 * freq**2, simple_in_cutoff
  else
    simple_in_cutoff = (N_abs <= 0.d0)
  end if
  end function simple_in_cutoff

  subroutine calculate_em(svec,omega, em, em_secondary, abs)
    ! S. Denk 4. 2013
    use mod_ecfm_refr_types,         only: ant, rad, non_maxwellian, &
                                           rad_diag_ch_mode_ray_freq_svec_type, &
                                           bi_max,drift_m, Spitzer, multi_slope, output_level, old_cutoff, &
                                           spl_type_2d, non_therm_params_type, dstf, dstf_comp
    use constants,                   only: pi, e0, mass_e, eps0, c0
    !use mod_ecfm_refr_fp_dist_utils, only: interpolate_f_rhop, export_Te_rhop
    !use ecfm_non_therm_abs,         only: abs_Albajar, abs_non_therm_tor_abs
    use mod_ecfm_radiation_dist,             only: radiation_dist_f_norm, prepare_dist
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(inout)  :: svec
    real(rkind), intent(in)                           :: omega
    real(rkind),    intent(out)                       :: em, em_secondary, abs                  ! [W m^-3 sr^-1 Hz^-1]
    character(30)                                     :: secondary_flag
    real(rkind)                                       :: freq, em_tot, em_tot_secondary, eta_2X, eta_2O, sf, sf_abs , &
                                                         sf_secondary,const, ne_frac, f_rel,f_rel_sq, mu, &
                                                         int_u, int_u_abs, f_norm,a, N, temp_norm , &
                                                         em_omega, abs_omega, em_secondary_omega, int_u_second
    integer(ikind)                                    ::  irhop_TDiff, irhop, m, l_omega
    type(spl_type_2d) :: f_spl
    type(non_therm_params_type)        :: dist_params
    em = 0.d0
    abs = 0.d0
    em_secondary = 0.d0
    if(svec%Te < 1.d0) return ! very low Te => absorption can be ignored
    if(svec%ne < 1.e17) return
    if(svec%N_cold <= 0.d0 .or. svec%N_cold > 1.d0) return
    freq = omega / (2.d0 * Pi)
    f_rel = freq /  svec%freq_2X
    if(dstf /= "relamax" .and. dstf /=  "Hu_nbes" .and.  dstf /= "Hu_bess") call prepare_dist(svec,Int_weights_many, Int_absz_many, f_spl, dist_params)
    do m = 2, 2
      sf = 0.d0
      sf_secondary = 0.d0
      sf_abs = 0.d0
      eta_2X = 0.d0
      eta_2O = 0.d0
      em_tot = 0.d0
      em_tot_secondary = 0.d0
      f_rel_sq = f_rel**2.d0
      const = (mass_e * c0**2.d0) / (2.d0 * e0 * svec%Te)
      mu = 2 * const
! Compare em/abs
      if(svec%cos_theta == 0.0 .or. svec%sin_theta == 0.0) then
          print*, "svec",svec
          print*, f_rel, svec%theta, svec%sin_theta, svec%cos_theta
          print*, "Sub_calculate_em"
          stop "Calculation not possible for theta = 90 deg. or theta == 0 deg."
      end if
      sf_secondary = 0.d0
      if(output_level) then
        if(trim(dstf_comp) == "Mx") then !"relamax"!! "maxwell"
          if(m == 2) then
            call calculate_sf_rrf(svec%Te, mu / 2.d0, svec%theta,       &
                  svec%cos_theta, svec%sin_theta, svec%freq_2X,f_rel, f_rel_sq, int_u_second)
!            print*, "doing the fischer", int_u_second
            call calculate_em_tot_secondary(svec%ne, &
              svec%cos_theta, &
              svec%sin_theta, &
              svec%freq_2X, const,m , em_tot_secondary)
            sf_secondary = int_u_second / (const**2 )
           else
            sf_secondary = 0.d0
          end if
        else if( trim(dstf_comp) == "NB") then
          call calculate_int_u_gauss_no_bessel(svec,mu,f_rel,f_rel_sq, int_u_second, int_u_abs,m,svec%N_cold)!
          call calculate_em_tot_secondary(svec%ne, &!
                    svec%cos_theta, &
                    svec%sin_theta, &
                    svec%freq_2X, const,m , em_tot_secondary)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf_secondary = int_u_second * pi/ freq * f_norm
        else if( trim(dstf_comp) == "Al") then
          call calculate_em_tot(svec%ne, &
            svec%cos_theta, &
            svec%sin_theta, &
            svec%freq_2X, const, m, em_tot_secondary)
          call calculate_int_u_gauss_refr(svec, mu,f_rel,f_rel_sq, int_u_second, int_u_abs,m,svec%N_cold)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf_secondary = int_u_second * pi/ freq * f_norm
        else
          print*, "Unknown dstf_comp:", dstf_comp
          stop "Critical ERROR in calcualte_em"
        end if
        sf = 0.0
        sf_abs = 0.0
        if( dstf == "Hu_nbes") then
          call calculate_int_u_gauss_no_bessel(svec,mu,f_rel,f_rel_sq, int_u, int_u_abs,m,svec%N_cold)!
!          print*, "doing no bessels whatsoever", int_u
          call calculate_em_tot_secondary(svec%ne, &!
            svec%cos_theta, &
            svec%sin_theta, &
            svec%freq_2X, const,m , em_tot)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf = int_u * pi/ freq * f_norm
          sf_abs =  - int_u_abs * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
        else if(dstf == "Hu_bess") then
          call calculate_int_u_gauss(svec, mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = a * (sqrt(mu / (2 * pi))**3)
          f_norm = sqrt(mu / (2 * pi))**3
          sf = int_u * pi/ freq * f_norm
          sf_abs =  - int_u_abs * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
          call calculate_em_tot(svec%ne, & !
                                svec%cos_theta, &
                                svec%sin_theta, &
                                svec%freq_2X, const,m , em_tot)
        else
          call calculate_int_u_gauss_arbitrary(svec, dist_params,mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)
          call calculate_em_tot(svec%ne, &
                                svec%cos_theta, &
                                svec%sin_theta, &
                                svec%freq_2X, const,m , em_tot)
          a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
          f_norm = radiation_dist_f_norm(svec%Te, dstf)
          sf = int_u * pi/ freq * f_norm
          sf_abs =  int_u * mu * pi / (freq**3 * mass_e)
          sf_abs = sf_abs * f_norm
        end if
      else
        call calculate_int_u_gauss(svec,mu,f_rel,f_rel_sq , int_u, int_u_abs,m,svec%N_cold)!_arbitrary
        a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
        f_norm = a * (sqrt(mu / (2 * pi))**3)
        sf = int_u * pi/ freq * f_norm
        call calculate_em_tot(svec%ne, &
                              svec%cos_theta, &
                              svec%sin_theta, &
                              svec%freq_2X, const,m , &
                              em_tot_secondary)
      end if
      call calculate_mode_frac(svec%cos_theta, svec%sin_theta, 0.5d0 / f_rel , eta_2X, eta_2O)
      em = em + eta_2X * sf * em_tot
      abs = abs + eta_2X * sf_abs * em_tot
      em_secondary = em_secondary + em_tot_secondary * eta_2X * sf_secondary
    end do
    if( em /= em  .or. em < 0.d0 .or. abs /= abs .or. abs < 0.d0) then !c_abs < 0.d0 .or.
        print*, "rhop", svec%rhop
        print*, "Te",svec%Te
        print*, "ne",svec%ne
        print*, "freq", freq
        print*, "freq 2",  svec%freq_2X
        print*, "eta_2X",eta_2X
        print*, "em_tot", em_tot
        print*, "int_u_second", int_u_second
        print*, "int_u", int_u
        print*, "sf", sf
        print*, "sf_abs", sf_abs
        print*, "em_tot_secondary",  em_tot_secondary
        print*, "sf_secondary", sf_secondary
        print*, "em", em
        print*, "abs", abs
        print*, "em_secondary", em_secondary
        stop "Nan in calculate_em"
      end if
  end subroutine calculate_em

  !*******************************************************************************

  subroutine calculate_em_tot(ne, cos_theta, sin_theta, freq_2X, const, m, em_tot)

    use constants,                  only: pi, e0, mass_e, eps0, c0

    implicit none

    real(rkind), intent(in)  :: ne, cos_theta, sin_theta, freq_2X, const
    integer(ikind), intent(in) :: m
    real(rkind), intent(out) :: em_tot
    real(rkind) ::  freq_m
    freq_m = freq_2X / 2 * m
    em_tot = (e0 * freq_m)**2.d0 * &
      ne / (eps0 * c0 )
  end subroutine calculate_em_tot

  subroutine calculate_em_tot_secondary(ne, cos_theta, sin_theta, freq_2X, const, m, em_tot)

    use constants,                  only: pi, e0, mass_e, eps0, c0

    implicit none

    real(rkind), intent(in)  :: ne, cos_theta, sin_theta, freq_2X, const
    integer(ikind), intent(in) :: m ! order of the harmonic
    real(rkind), intent(out) :: em_tot
    real(rkind) ::  freq_m
    integer(ikind) :: i, fac
    fac = 1
    do i = 1,m - 1
      fac = fac * i
    end do
    freq_m = freq_2X / 2 * m
    em_tot = (m**(2*(m-1))/fac**2) / (2.d0**(2 * m))* sin_theta**(2*(m-1)) * (cos_theta**2.d0  + 1.d0)
    em_tot = em_tot * (e0 * freq_m)**2.d0 * &
      ne / (eps0 * c0 )
  end subroutine calculate_em_tot_secondary

  !*******************************************************************************

  subroutine calculate_mode_frac(cos_theta, sin_theta, Y, eta_2X, eta_2O)

    implicit none

    real(rkind), intent(in)  :: cos_theta, sin_theta, Y
    real(rkind), intent(out) :: eta_2X, eta_2O
    real(rkind)  :: a1, a2, a3

    a1 = sin_theta**4.d0 / (4.d0 / Y)
    a2 = cos_theta**2.d0
    a3 = (a1 + a2) / ( (a2+1.d0) * sqrt(a2+a1 * Y) )

    eta_2X = 0.5d0 + a3
    eta_2O = 0.5d0 - a3

  end subroutine calculate_mode_frac

  !*******************************************************************************


  !*******************************************************************************

  subroutine calculate_sf_rrf(Te, const, theta, cos_theta, &
    sin_theta, freq_2X, f_rel, f_rel_sq, sf)

    use mod_ecfm_refr_types,        only: ant
    use constants,                  only: pi, e0, mass_e, c0, sPi

    implicit none

    real(rkind),    intent(in)  :: Te, const
    real(rkind),    intent(in)  :: theta, cos_theta, sin_theta
    real(rkind),    intent(in)  :: freq_2X, f_rel, f_rel_sq
    real(rkind),    intent(out) :: sf

    real(rkind)     :: eps, eta, xh1, xh2, prefact, t1, t2, sin_theta_2, num_eps, freq_new
    real(rkind)     :: beta_par_1, beta_par_2, alpha1, alpha2
    real(rkind)     :: e1, e2, sa1, sa2, ca1, ca2, F1, F2

    real(rkind)     :: sf2
    sin_theta_2 = sin_theta**2
    eps = 1.d0 / (const * f_rel_sq)                    ! a small quantity
    eta = 1.d0 + f_rel_sq * cos_theta**2               ! 1 + (freq/freq_2X * cos(theta))**2
    if(1.d0 - f_rel_sq * sin_theta_2 < 0.d0) then
      sf = 0.d0
      return
    end if
    t1 = f_rel_sq * cos_theta                          ! (freq/freq_2X)**2 * cos(theta)
    t2 = sqrt(1.d0 - f_rel_sq * sin_theta_2)           ! sqrt(1-(freq/freq_2X * sin(theta))**2)
    beta_par_1 = (t1 - t2) / eta                       ! lower integration limit
    beta_par_2 = (t1 + t2) / eta
    num_eps = abs(beta_par_2 - beta_par_1) * 1.d-9
    beta_par_1 = beta_par_1 + num_eps
    beta_par_2 = beta_par_2 - num_eps
    sa1 = 1.d0 - beta_par_1 * cos_theta                ! uses beta_par_1
    sa2 = 1.d0 - beta_par_2 * cos_theta                ! uses beta_par_2
    alpha1 = sa1**2                                    ! uses beta_par_1
    alpha2 = sa2**2                                    ! uses beta_par_2
    e1  = exp(-const*(1.d0 - f_rel_sq * alpha1))
    e2  = exp(-const*(1.d0 - f_rel_sq * alpha2))
    ca1 = sqrt(const*f_rel_sq*alpha1)
    ca2 = sqrt(const*f_rel_sq*alpha2)
    xh1 = eta**2 * eps
    xh2 = (3.d0*eta - 2.d0*sin_theta_2/eps) * sqrt(eps)

    call dawson_integral(ca1, F1)
    call dawson_integral(ca2, F2)
    freq_new = f_rel * freq_2X
    prefact = const**1.5d0 / (sPi * freq_new * f_rel_sq * cos_theta**5)

    sf = prefact * ( e1 * (xh1 - sin_theta_2 / sa1 - xh2 * F1) &
      -e2 * (xh1 - sin_theta_2 / sa2 - xh2 * F2) )
    if(sf < 0.d0) sf = 0.d0
  end subroutine calculate_sf_rrf

  !*******************************************************************************

  subroutine dawson_integral(x, daw)
    ! Dawson integral: daw = exp(-x**2) int_0^x exp(u**2) du

    use nr_special_funcs,           only: dawson_s

    implicit none

    real(rkind), intent(in)  :: x
    real(rkind), intent(out) :: daw   ! value of Dawson integral
    real(rkind)  :: y

    if (x >= 4.8d0) then
      y   = 1.d0/(2.d0*x**2)
      daw = (1.d0 + y*(1.d0 + 3.d0*y*(1.d0 + 5.d0*y*(1.d0 + 7.d0*y*(1.d0 + 9.d0*y*(1.d0 + 11.d0*y)))))) / (2.d0*x)
    else
      daw = dawson_s(x)
    endif

  !stop 'subroutine dawson_integral'
  end subroutine dawson_integral

  !*******************************************************************************

subroutine radiation_gauss_init(N_step_sf_int)
#ifdef NAG
    use nag_quad_util,              only: nag_quad_gs_wt_absc
    USE nag_error_handling
#endif
    use quadrature,                 only: cdgqf
    implicit none
    integer(ikind), intent(in)     :: N_step_sf_int
#ifdef NAG
    type(nag_error)                :: error
#endif
    integer(ikind)                 :: i
    real(rkind)                    :: h
    real(rkind), dimension(:), allocatable :: Int_weights_check, Int_absz_check
    allocate( Int_weights(N_step_sf_int),Int_absz(N_step_sf_int))
    allocate( Int_weights_many(64),Int_absz_many(64))
    call cdgqf( int(N_step_sf_int,kind=4), int(1,kind=4), 0.d0, 0.d0, Int_absz, Int_weights)
    call cdgqf( int(64,kind=4), int(1,kind=4), 0.d0, 0.d0, Int_absz_many, Int_weights_many)
#ifdef NAG
    allocate( Int_weights_check(N_step_sf_int), Int_absz_check(N_step_sf_int))
    call nag_quad_gs_wt_absc( 0, -1.d0, 1.d0, Int_weights_check, Int_absz_check)
    if(sum((Int_weights - Int_weights_check)**2) > 1.d-5) then
      print*, "Weights deviate by more than 1.d-10"
      do i = 1, N_step_sf_int
        print*, Int_weights(i), Int_weights_check(i)
      end do
      call abort()
    end if
    if(sum((Int_absz - Int_absz_check)**2) > 1.d-5) then
      print*, "Abszissae deviate by more than 1.d-10"
      do i = 1, N_step_sf_int
        print*, Int_absz(i), Int_absz_check(i)
      end do
      call abort()
    end if
    deallocate(Int_weights_check, Int_absz_check)
    !call nag_quad_gs_wt_absc( 0, -1.d0, 1.d0, Int_weights, Int_abszz)
    if (error%level >= 1) print*, error%msg
    if (error%level >= 1) stop "Failed to setup "
#endif
    total_cnt = 0
    shortcut_cnt = 0
  end subroutine radiation_gauss_init

  subroutine radiation_gauss_clean_up()
    use mod_ecfm_refr_types,        only: output_level
  implicit none

    if(allocated(Int_weights)) deallocate(Int_weights)
    if(allocated(Int_absz)) deallocate(Int_absz)
    if(allocated(Int_weights_many)) deallocate(Int_weights_many)
    if(allocated(Int_absz_many)) deallocate(Int_absz_many)
    if(output_level) print*, shortcut_cnt, " of ", total_cnt, " were calculated quickly"
  end subroutine radiation_gauss_clean_up




   subroutine calculate_int_u_gauss(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((real(m,8)/ (f_rel *2.0 ))**2  < 1.e0 - svec%cos_theta**2) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2
    cos_theta_2 = (cos_theta)**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta  + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 *  f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_JN(m, zeta(k)/ sin_theta))**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
      (Bessel_JN(m - 1, zeta(k)) - Bessel_JN(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * exp( -mu / 2.d0 * ((u_par(k) / gamma(k))**2 + u_perp(k)**2 / gamma(k)**2))* Int_weights(k)
      !* exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights(k)
      !if(stop) print*, u_par(k), int_u_add(k)
      !if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then

        !print*,"trub",
        !k_int = k
       ! stop = .true.
        !goto 1
        !return

      !end if
    enddo
    !if(stop) stop "done"
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
    end do
    int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss

  subroutine calculate_int_u_gauss_no_bessel(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2! * N**2
    cos_theta_2 = (cos_theta)**2! * N**2
    if(m_omega_bar_2 - sin_theta_2 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta  + m_omega_bar !* N
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      int_u_add(k) = u_perp(k)**(2*m) / gamma(k)**(2*m - 2)
      int_u_add(k) = int_u_add(k) * du_dalpha * exp( -mu / 2.d0 * ((u_par(k) / gamma(k))**2 + u_perp(k)**2 / gamma(k)**2)) * Int_weights(k)
      !exp( mu* (1.d0 - gamma(k))) * Int_weights(k)
      !if(stop) print*, u_par(k), int_u_add(k)
!      if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then
!
!        !print*,"trub",
!        !k_int = k
!        stop = .true.
!        !goto 1
!        !return
!
!      end if
    enddo
    !if(stop) stop "done"
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
    end do
    int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_no_bessel


  subroutine calculate_int_u_gauss_arbitrary(svec, dist_params, mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, k_int, non_therm_params_type
    use mod_ecfm_radiation_dist,    only: radiation_dist_f_u, radiation_dist_Rf,radiation_dist_f_norm
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    type(non_therm_params_type), intent(in) :: dist_params
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2, det!,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights_many))     :: u_par,u_perp, gamma, zeta, int_u_add, int_u_abs_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, N_par, N_perp
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    N_par = svec%N_cold * cos_theta !
    N_perp = svec%N_cold * sin_theta !
    m_omega_bar =  real(m,8) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    if(m_omega_bar_2  + N_par**2  - 1.e0 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    !sin_theta_2 = (sin_theta)**2!  * N**2
    !cos_theta_2 = (cos_theta)**2! * N**2
    !det = sqrt(m_omega_bar_2 - sin_theta_2)
    !a = (m_omega_bar*cos_theta - det) / sin_theta_2
    !b = (m_omega_bar*cos_theta + det) / sin_theta_2
!    if(m_omega_bar**2  + N_par**2  - 1.e0 < 0.d0) then
!      Int_u = 0.d0
!      int_u_abs = 0.d0
!      return
!    end if
    a =  (m_omega_bar * N_par - sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    b = (m_omega_bar * N_par + sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    !1 continue
    do k = 1, size(Int_weights_many)
      u_par(k) = Int_absz_many(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * N_par + m_omega_bar
      !gamma(k) = u_par(k) * cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)

      ! Hutchinson
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_Jn(m, zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
      (Bessel_JN(m - 1, zeta(k)) - Bessel_Jn(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * gamma(k)**2 * Int_weights_many(k)
      !if(stop)  print*,  u_par(k),  radiation_dist_Rf(u_par(k), u_perp(k), gamma(k), m_omega_bar, cos_theta, mu, svec)
      !if(u_perp(k) /= u_perp(k) .or. int_u_add(k) /= int_u_add(k)) then
      !  print*, m_omega_bar**2  + N_par**2  - 1.e0
      !stop "nan"
      !end if
      !if(abs(u_par(k)) < 0.1 .and. u_perp(k) < 0.1 .and. svec%Te > 3000 .and. .not. stop) then

        !print*,"trub",
        !k_int = k
        !stop = .true.
        !goto 1
        !return
      !end if
      !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
        !print*, "u_par",u_par(k)
        !print*, "u_perp", u_perp(k)
        !print*, "f", radiation_dist_f_u(u_par(k), u_perp(k), gamma(k), mu, svec)
        !print*, "int_u",int_u_add(k)
      !  print*,u_par(k),radiation_dist_f_u(u_par(k), u_perp(k), gamma(k), mu, svec)
      !end if
    enddo
    int_u_abs_add = int_u_add * radiation_dist_Rf(u_par, u_perp, gamma, m_omega_bar, cos_theta, mu, svec, dist_params)
    int_u_add = int_u_add * radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params)
    !if(stop) stop "done"
    !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
    !  stop "Sense?"
    !end if
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights_many)
      Int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs  +  int_u_abs_add(k)
    end do
    !int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_arbitrary

  subroutine calculate_int_u_gauss_refr(svec, mu, f_rel, f_rel_sq, int_u, int_u_abs, m,N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, k_int, non_therm_params_type
    use mod_ecfm_radiation_dist,    only: radiation_dist_f_u, radiation_dist_Rf,radiation_dist_f_norm
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind)     :: m_omega_bar_2, det!,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add, int_u_abs_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, N_par, N_perp
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta = svec%sin_theta
    N_par = svec%N_cold * cos_theta !
    N_perp = svec%N_cold * sin_theta !
    m_omega_bar =  real(m,8) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    if(m_omega_bar_2  + N_par**2  - 1.e0 < 0.d0) then
      int_u = 0.d0
      int_u_abs = 0.d0
      return
    end if
    !sin_theta_2 = (sin_theta)**2!  * N**2
    !cos_theta_2 = (cos_theta)**2! * N**2
    !det = sqrt(m_omega_bar_2 - sin_theta_2)
    !a = (m_omega_bar*cos_theta - det) / sin_theta_2
    !b = (m_omega_bar*cos_theta + det) / sin_theta_2
!    if(m_omega_bar**2  + N_par**2  - 1.e0 < 0.d0) then
!      Int_u = 0.d0
!      int_u_abs = 0.d0
!      return
!    end if
    a =  (m_omega_bar * N_par - sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    b = (m_omega_bar * N_par + sqrt(m_omega_bar**2  + N_par**2  - 1.e0)) / ( 1.e0 - N_par**2)
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    !!$OMP PARALLEL DO &
    !!$OMP SHARED(mu, cos_theta, m_omega_bar,sin_theta,b,a) PRIVATE(k,u_par, u_perp, zeta,gamma, int_u_add2) &
    !!$OMP REDUCTION(+:int_u)
    !1 continue
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * N_par + m_omega_bar
      !gamma(k) = u_par(k) * cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)

      ! Hutchinson
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* Bessel_Jn(m, zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + ( u_perp(k)/gamma(k))**2  * &
                     (Bessel_JN(m - 1, zeta(k)) - Bessel_Jn(m + 1,zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha * gamma(k)**2 * Int_weights(k)
    enddo
    int_u_abs_add = -int_u_add * mu * exp( -mu / 2.d0 * ((u_par(:) / gamma(:))**2 + u_perp(:)**2 / gamma(:)**2))
    int_u_add = int_u_add  * exp( -mu / 2.d0 * ((u_par(:) / gamma(:))**2 + u_perp(:)**2 / gamma(:)**2))
    !if(stop) stop "done"
    !if(svec%rhop > 0.6 .and. svec%R > 1.7) then
    !  stop "Sense?"
    !end if
    !!$OMP END PARALLEL DO
    do k = 1, size(Int_weights)
      Int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs  +  int_u_abs_add(k)
    end do
    !int_u_abs = -int_u * mu
  end subroutine calculate_int_u_gauss_refr


  function full_relativistics(cos_theta, sin_theta, freq_2X, freq)
   ! Returns true if the resonance curve is in the relativistic regime of velocity space
   ! the relativistic regime is defined here as beta >= 1 %
   ! Therefore gamma may not be larger than approximately 1.005
   use mod_ecfm_refr_types,             only: rad_diag_ch_mode_ray_freq_svec_type
   implicit none
   real(rkind),                 intent(in) :: cos_theta, sin_theta, freq_2X, freq
   logical                                 :: full_relativistics
   real(rkind)                             :: u_par_max
   real(rkind)                             :: gamma_max
   full_relativistics = .true.
   u_par_max = ((freq_2X/freq)*cos_theta + sqrt((freq_2X/freq)**2 - sin_theta**2)) / sin_theta**2
   gamma_max = (freq_2X/freq) + u_par_max * cos_theta
   total_cnt = total_cnt + 1
   if(u_par_max / gamma_max < 1.d-1) full_relativistics = .false.
   !full_relativistics = .true.
  end function full_relativistics

  subroutine calculate_int_u_gauss_opt(cos_theta, sin_theta,mu, f_rel, f_rel_sq, int_u) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none
    real(rkind),    intent(in)  :: cos_theta,sin_theta,mu,f_rel, f_rel_sq
    real(rkind),    intent(out) :: int_u
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind), dimension(size(Int_weights))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: eps, besselj1,besselj2, besselj3
    real(rkind)     :: a, b, du_dalpha, int_u_add2, ab_mean
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_omega_bar =  1.0 / f_rel
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2
    cos_theta_2 = (cos_theta)**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    ab_mean = (b + a) / 2.d0
    int_u       = 0.d0
    do k = 1, size(Int_weights)
      u_par(k) = Int_absz(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k) * cos_theta + m_omega_bar
     if(gamma(k)**2 - u_par(k)**2 - 1.d0 < 0) then
      print*, k
      print*, u_par(k)
      print*, gamma(k)
      stop "imaginary value"
     end if
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta - u_par(k)/gamma(k))* BesselJ_fast_2(zeta(k))/ sin_theta)**2
      int_u_add(k) = int_u_add(k)  + &
        ( u_perp(k)/gamma(k))**2  * (BesselJ_fast_1(zeta(k)) - BesselJ_fast_3(zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha *  exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights(k)
    enddo
    do k = 1, size(Int_weights)
      int_u = int_u + int_u_add(k)
    end do
  end subroutine calculate_int_u_gauss_opt

  function BesselJ_fast_1(x)
  ! Evaluates the Besselfunction of the first kind for n = 1
  ! Fast implementation of BesselJ_fast_1
  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_1
    ! BesselJ = 1/ (n !)
    BesselJ_fast_1 = 1.0 - x**2/8.d0 ! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_1 = BesselJ_fast_1 + x**4/192.d0 ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_1 = BesselJ_fast_1 - x**6/(9216.d0) ! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_1 = BesselJ_fast_1 * x / 2.0
  end function BesselJ_fast_1

  function BesselJ_fast_2(x)
  ! Evaluates the Besselfunction of the first kind for n = 1

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_2
    ! BesselJ = 1/ (n !)                :: n_fac
    BesselJ_fast_2 = 0.5 - x**2/24.d0 ! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_2 = BesselJ_fast_2 + x**4/768.d0 ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_2 = BesselJ_fast_2 - x**6/(46080.d0) ! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_2 = BesselJ_fast_2 * x**2 / 4.0
  end function BesselJ_fast_2

  function BesselJ_fast_3(x)
  ! Evaluates the Besselfunction of the first kind for n = 1

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    real(rkind)                    :: BesselJ_fast_3
    ! BesselJ = 1/ (n !)
    BesselJ_fast_3 = 1/6.d0 - x**2/(96.d0)! x^2/ ((n + 1)! * 2^2)
    BesselJ_fast_3 = BesselJ_fast_3 + x**4/(3840.d0) ! x^4/((n + 2)! * 2^5)
    BesselJ_fast_3 = BesselJ_fast_3 - x**6/(276480.d0)! x^6/((n + 3) ! * 2^7)
    BesselJ_fast_3 = BesselJ_fast_3 * x**3 / (8.0)
  end function BesselJ_fast_3

  function BesselJ(n,x)
  ! Evaluates the Besselfunction of the first kind
  ! For n = 1, the relativ numerical error for x < 1.5 is below 10^-4
  ! For values in between 1.5 < x < 2.5 the numerical error is below 1 %
  ! For values 2.5 > x > 3.0 the numerical error is below 10 %
  ! For values x > 3 this routine should not be used since the truncation error becomes very large
  ! The larger n, the greater the range where this function gives an almost exact result
  ! e.g. for  n = 2 the numerical for x = 3.0 is still below 1 %

  ! Optimization remark: This funciton is intended to be used in a vectorized loop.
    implicit none
    real(rkind),       intent(in)  :: x
    integer(ikind),    intent(in)  :: n
    real(rkind)                    :: besselj
    integer(ikind)                 :: i
    real(rkind)                    :: n_fac
    n_fac = 1.0
    do i= 1,n
      n_fac = real(i) * n_fac
    end do
    BesselJ = 1.0/n_fac
    n_fac = n_fac * real(n + 1)
    BesselJ = BesselJ - x**2/(n_fac * 2.d0**(2))
    n_fac = n_fac * real(n + 2)
    BesselJ = BesselJ + x**4/(n_fac * 2.d0**(5))
    n_fac = n_fac * real(n + 3)
    BesselJ = BesselJ - x**6/(n_fac * 2.d0**(7)  * 3.d0)
    BesselJ = BesselJ * x**n / (2.0)**n
  end function BesselJ



  subroutine calculate_int_u_gauss_many(svec,mu, f_rel, f_rel_sq, int_u, int_u_abs, m, N) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    real(rkind),    intent(in)  :: mu,f_rel, f_rel_sq, N
    real(rkind),    intent(out) :: int_u,int_u_abs
    integer(ikind), intent(in)  :: m
    ! for numerial calculation of integral
    real(rkind)     :: m_omega_bar,cos_theta,sin_theta
    real(rkind), dimension(size(Int_weights_many))     :: u_par,u_perp, gamma, zeta, int_u_add
    real(rkind)     :: m_omega_bar_2,cos_theta_2,sin_theta_2, det
    real(rkind)     :: eps
    real(rkind)     :: a, b, du_dalpha
    integer(ikind)  :: k

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of u integral using gaussian quadrature with boundaries -1.0,1.0 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cos_theta = svec%cos_theta
    sin_theta =svec%sin_theta
    m_omega_bar =  real(m) / (2.0 * f_rel)
    m_omega_bar_2 = (m_omega_bar)**2
    sin_theta_2 = (sin_theta)**2! * N**2
    cos_theta_2 = (cos_theta)**2! * N**2
    det = sqrt(m_omega_bar_2 - sin_theta_2)
    a = (m_omega_bar*cos_theta - det) / sin_theta_2
    b = (m_omega_bar*cos_theta + det) / sin_theta_2
    eps = (b - a) * 1.d-9 !small quantity to assure we stay in the range of the resonance
    ! If eps is equal to zero, there are cases were floating point precision is too low
    ! and the square root of a very small, but negative, number is calculated
    ! At the boundaries u_perp is very small and therefore the single particle emissivity
    ! is neglidibly small.
    ! There should be no noticeable error through this slight tightening of the integral!
    a = a + eps
    b = b - eps
    du_dalpha = (b - a)/2.d0
    int_u       = 0.d0
    int_u_abs   = 0.d0
    do k = 1, size(Int_weights_many)
      u_par(k) = Int_absz_many(k) * du_dalpha + (b + a) / 2.d0
      gamma(k) = u_par(k)  *  cos_theta + m_omega_bar
      u_perp(k) = sqrt(gamma(k)**2 - u_par(k)**2 - 1.d0)
      zeta(k) = 2 * f_rel * u_perp(k) * sin_theta
      int_u_add(k) = ((cos_theta  - u_par(k)/gamma(k))* BesselJ_fast_2(zeta(k))/ sin_theta )**2
      int_u_add(k) = int_u_add(k)  + &
        ( u_perp(k)/gamma(k))**2  * (BesselJ_fast_1(zeta(k)) - BesselJ_fast_3(zeta(k)))**2 / 4.d0
      int_u_add(k) = int_u_add(k) * du_dalpha *  exp( mu* (1.d0 - gamma(k))) * (gamma(k)**2) * Int_weights_many(k)
    end do
    do k = 1, size(Int_weights_many)
      int_u = int_u + int_u_add(k)
      int_u_abs = int_u_abs - int_u_add(k) * mu
    enddo
    if(int_u < 0 .or. int_u /= int_u) then
      print*, "Negative Shape function"
      print*, Int_absz_many(1) * du_dalpha + (b + a) / 2.d0
      print*,(Int_absz_many(1) * du_dalpha + (b + a) / 2.d0)*cos_theta + m_omega_bar
      stop "Bug"
    end if
  end subroutine calculate_int_u_gauss_many

function bessel_term(u_par, u_perp, gamma,zeta, m, cos_theta, sin_theta)
implicit none

real(rkind),    intent(in)  :: u_par, u_perp, gamma,zeta, cos_theta, sin_theta
integer(ikind), intent(in)  :: m
real(rkind)                 :: bessel_term
 bessel_term = ((cos_theta - u_par/gamma)/ (sin_theta))**2 * Bessel_JN(m, zeta)**2
 bessel_term = bessel_term +  ( u_perp/gamma)**2  * &
   (Bessel_JN(m - 1, zeta) - Bessel_JN(m + 1,zeta))**2 / 4.d0

end function bessel_term
end module mod_ecfm_refr_em_Hu

! module abs albajar
  function func_dj_du_par(u_par)
  ! Not used anymore, gives df/dupar to find maximum point of emission
  implicit none
  real(rkind), intent(in) :: u_par
  real(rkind)             :: func_dj_du_par
  real(rkind)             :: m_omega_bar, cos_theta, mu
    m_omega_bar = m_omega_glob
    cos_theta = cos_theta_glob
    mu = mu_glob
    func_dj_du_par = (-1 + (-1 + cos_theta ** 2) * u_par ** 2 + 2 * cos_theta * u_par * m_omega_bar + m_omega_bar ** 2) * &
        (cos_theta ** 4 * u_par ** 3 * mu + 4 * u_par * m_omega_bar + cos_theta ** 3 * u_par ** 2 * (-2 + 3 * mu * m_omega_bar) - &
         cos_theta * (2 + mu * m_omega_bar + 2 * m_omega_bar ** 2 - mu * m_omega_bar ** 3 + u_par ** 2 * (-2 + mu * m_omega_bar)) - &
         cos_theta ** 2 * u_par * (4 * m_omega_bar + mu * (1 + u_par ** 2 - 3 * m_omega_bar ** 2)))
  return
  end function func_dj_du_par

 subroutine abs_Al_pol_fact_perp(svec, t, X, Y, omega_bar, m_0, N_abs, e, mode, m, pol_fact)
  ! According to formula (2a and 2c) in [1]
  ! The 1D ECE is slightly oblique and the Inline and imagining systems are very oblique,
  ! hence the approximation of N_par = 0 is not appropriate
    use mod_ecfm_refr_types,        only: ant, rad_diag_ch_mode_ray_freq_svec_type
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)   :: svec
    real(rkind), dimension(:), intent(in)              :: t
    real(rkind), intent(in)                            ::  X, Y, omega_bar, m_0, N_abs
    complex(r8), dimension(:), intent(in)              :: e
    integer(ikind), intent(in)                         :: mode, m
    real(rkind), dimension(:), intent(out)             :: pol_fact
    real(rkind)                   :: x_m
    real(rkind)                   :: N_eff, Axz_sq, Re_Axz_ey, Re_Axz_ez, Re_ey_ez, N_gray, ey_sq, ez_sq, abs_c
    real(rkind), dimension(3,3)  :: E_mat
    complex(r8), dimension(3)  :: pol_vect
    complex(r8)                :: Axz
    real(rkind), dimension(size(t))  :: bessel_arg, bessel_n_l, bessel_n_2, bessel_n, bessel_n_u, abs_Al_bessel_sqr_deriv
    logical                      :: cold
    cold = .true. ! to see if Gray the polarzation vector increases accuracy -> cold = .false.
    x_m =  N_abs * omega_bar * sqrt((real(m,8)/ m_0)**2 - 1.d0)
    N_eff = 0.d0
    if(cold) then
      pol_vect = e
!      call get_E_factors(X, Y, N_abs, N_perp, N_par, e,  E_mat )
!      ! E_mat (1,2), (2,1) , (3,2) and (2,3) have additional factor of i
!      Axz_sq = E_mat(1,1)!  + N_eff**2 * E_mat(3,3) + 2.d0 * N_eff * E_mat(1,3) ! mixed terms do not vanish!
!      Re_Axz_ey = -E_mat(1,2)! - N_eff * E_mat(3,2) ! this includes the additional i
!      print*, "E_mat(1,2), N_eff * E_mat(3,2)", E_mat(1,2), N_eff * E_mat(3,2)
!      Re_Axz_ez = 0.d0! E_mat(1,3) + N_eff * E_mat(3,3)
!      Re_ey_ez = 0.d0! -E_mat(3,2)
!      ey_sq = E_mat(2,2)
!      ez_sq = 0.d0!E_mat(3,3)
    else
      N_gray = N_abs
      abs_c =  abs_Al_tor_abs(svec, svec%freq_2X * Pi / Y, mode, N_gray, pol_vec =  pol_vect)
    end if
    Axz = e(1)
!    print*, "Re_Axz^2 cold", Axz_sq
    Axz_sq = abs(Axz)**2 ! mixed terms do not vanish!
!    print*, "Re_Axz^2 warm", Axz_sq
!    print*, "Re_Axz_ey cold", Re_Axz_ey
    Re_Axz_ey = real(cmplx(0.d0, 1.d0) * Axz * conjg(pol_vect(2)))
!    print*, "Re_Axz_ey warm", Re_Axz_ey
!    print*, "Re_Axz_ez cold", Re_Axz_ez
    Re_Axz_ez = 0.d0
!    print*, "Re_Axz_ez warm", Re_Axz_ez
!    print*, "Re_ey_ez warm", Re_ey_ez
    Re_ey_ez = 0.d0
    ey_sq = abs(pol_vect(2))**2
    ez_sq = 0.d0
!    print*, "pol_vect cold", sqrt(E_mat(1,1)), sqrt(E_mat(2,2)), sqrt(E_mat(3,3))
!    print*, "normalization  cold", sqrt(E_mat(1,1) + E_mat(2,2) + E_mat(3,3))
!    print*, "pol_vect warm", pol_vect
    bessel_arg = x_m * sqrt(1.d0 - t(:)**2)
    bessel_n_l = BesselJ(m - 1 , bessel_arg)
    bessel_n = BesselJ(m , bessel_arg)
    bessel_n_2 = BesselJ(m , bessel_arg)**2
    bessel_n_u = BesselJ(m + 1 , bessel_arg)
    abs_Al_bessel_sqr_deriv = bessel_arg / x_m * bessel_n * ( bessel_n_l - bessel_n_u)
    pol_fact(:) = ( Axz_sq  + ey_sq) * bessel_n_2
    pol_fact(:) = pol_fact(:) + Re_Axz_ey * x_m / real(m,8) * abs_Al_bessel_sqr_deriv
    pol_fact(:) = pol_fact(:) - (bessel_arg / real(m,8))**2 * &
      ey_sq * bessel_n_l * bessel_n_u
    pol_fact(:) = pol_fact(:) + (x_m / (real(m,8)))**2 * &
      ez_sq * t(:)**2 * bessel_n_2
    pol_fact(:) = pol_fact(:) + x_m / (real(m,8)) * &
      2.d0 * Re_Axz_ez * t(:) * bessel_n_2
    pol_fact(:) = pol_fact(:) + x_m / (real(m,8)) * &
      Re_ey_ez * t(:) * x_m / real(m,8) * abs_Al_bessel_sqr_deriv !
    pol_fact(:) = pol_fact(:)  * (real(m,8)  / (N_abs * omega_bar))**2
  end subroutine abs_Al_pol_fact_perp


! Deprecated - this routine encounters difficulties near perpendicular propagation
  subroutine get_E_factors(X, Y, N_abs, cos_theta, sin_theta, e, E_mat )
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
!    ! Computeses the matrix:
!    ! (ex ex* ex ey* ex ez*)
!    ! (ey ex* ey ey* ey ez*)
!    ! (ez ex* ez ey* ez ez*)
!    ! From Bornatici Review 1983 p. 1198 eq. 3.1.61 and 3.1.62
!    ! No support for first harmonic nor omega_p > omega_c
!    ! ex := A c ey ! c contains all complex factors
!    ! ez := B c' ey ! c' contains all complex factors
!    ! ey ey^* = 1 / (N sqrt(a_sq + b_sq))
    real(rkind), intent(in) :: X, Y, N_abs,  sin_theta, cos_theta
    complex(r8), dimension(:),  intent(in) :: e
    real(rkind), dimension(:,:), intent(out) :: E_mat
    real(rkind)         :: A, B, ey_sq, re_norm,rho, f
    integer(ikind)      :: i
    rho =  Y**2 * sin_theta**4 + 4.d0 * (1.d0 - X)**2 * cos_theta**2
    rho = sqrt(rho)
    f =  (2.d0 * (1.d0 - X)) / (2.d0 * (1.d0 - X) - Y**2 * sin_theta**2 - real(1,8) *  Y* rho)
    A = 1.d0 / Y * (1.d0 - (1 - Y**2) * f)
    B = (N_abs**2 * sin_theta * cos_theta) / ( 1.d0 - X -  N_abs**2 * sin_theta**2) * A
    ey_sq = abs(e(2))**2
!    ! 13,6 Introduced normalization -> |(A 1 B) e_y| = 1
!    !ey_sq = ey_sq / re_norm
    !This what the E_mat looks in the complex realm
    !E_mat(1,1) = CMPLX(A**2,0.d0)
    !E_mat(1,2) = CMPLX(0.d0,A)
    !E_mat(1,3) = CMPLX(-A * B, 0.d0)
    !E_mat(2,1) = CMPLX(0.d0,-A)
    !E_mat(2,2) = CMPLX(1.d0, 0.d0)
    !E_mat(2,3) = CMPLX(0.d0,B)
    !E_mat(3,1) = CMPLX(-A * B, 0.d0)
    !E_mat(3,2) = CMPLX(0.0d, -B)
    !E_mat(3,3) = CMPLX(B**2, 0.d0)
     E_mat(1,1) = A**2 ! real
     E_mat(1,2) = A ! imaginary (* i)
     E_mat(1,3) = -A * B ! real
     E_mat(2,1) = -A !imaginary (* i)
     E_mat(2,2) = 1.d0 ! real
     E_mat(2,3) = B !imaginary (* i)
     E_mat(3,1) = -A * B ! real
     E_mat(3,2) = -B !imaginary (* i)
     E_mat(3,3) = B**2 ! real
     E_mat(:,:) = E_mat(:,:) * ey_sq
    ! Square of norm enters this matrix
    ! Although almost normalized - normalize it here
    !re_norm = abs(E_mat(1,1)) + abs(E_mat(2,2)) +abs(E_mat(3,3))
    !E_mat(:,:) = E_mat(:,:) /  re_norm
    !print*, "re_norm", re_norm
    if(Any(abs(E_mat(:,:)) > 1.d99) .or. Any(E_mat(:,:) /= E_mat(:,:))) then
      print*, A, Y
      stop "E_mat"
    end if
  !if(Y > 0.498 .and. Y < 0.502 .and. sqrt(X) > 0.3) then
  !  print*,sqrt(E_mat(1,1)),sqrt(E_mat(2,2)),sqrt(E_mat(3,3)),sqrt(E_mat(1,1) + E_mat(2,2) + E_mat(3,3))
  !end if
  end subroutine get_E_factors


!  subroutine get_E_factors_cmplx(X, Y, N_abs, N_perp, N_par, f, a_sq, b_sq, E_mat )
!    use constants,                  only: pi, e0, mass_e, eps0, c0
!    implicit none
!    ! Computeses the matrix:
!    ! (ex ex* ex ey* ex ez*)
!    ! (ey ex* ey ey* ey ez*)
!    ! (ez ex* ez ey* ez ez*)
!    ! From Bornatici Review 1983 p. 1198 eq. 3.1.61 and 3.1.62
!    ! No support for first harmonic nor omega_p > omega_c
!    ! ex := A c ey ! c contains all complex factors
!    ! ez := B c' ey ! c' contains all complex factors
!    ! ey ey^* = 1 / (N sqrt(a_sq + b_sq))
!    real(rkind), intent(in) :: X, Y, N_abs, N_perp, N_par, f, a_sq, b_sq
!    complex(r8), dimension(:,:), intent(out) :: E_mat
!    real(rkind)         :: A, B, ey_sq, re_norm
!    integer(ikind)      :: i
!    A = 1.d0 / Y * (1.d0 - (1 - Y**2) * f)
!    B = (N_par * N_perp) / ( 1.d0 - X -  N_perp**2) * A
!    ey_sq = 1.d0 / (N_abs *  sqrt(a_sq + b_sq))
!    E_mat(1,1) = cmplx(A**2, 0.d0)
!    E_mat(1,2) = cmplx(0.d0, -A)
!    E_mat(1,3) = cmplx(- A * B, 0.d0)
!    E_mat(2,1) = cmplx(0.d0, A)
!    E_mat(2,2) = cmplx(1.d0, 0.d0)
!    E_mat(2,3) = cmplx(0.d0, -B)
!    E_mat(3,1) = cmplx(- A * B, 0.d0)
!    E_mat(3,2) = cmplx(0.d0, B)
!    E_mat(3,3) = cmplx(B**2, 0.d0)
!    !E_mat(3,:) = 0.d0
!    !E_mat(:,3) = 0.d0 ! Quasiperpendicular X-mode
!    !re_norm = 1.d0 !/ (sqrt(ey_sq) * sqrt(1.d0 + E_mat(1,1) + E_mat(3,3)))
!    E_mat(:,:) = E_mat(:,:) * ey_sq
!    !do i = 1, 3
!    !  E_mat(:,i) = E_mat(:,i) * ey_sq * re_norm
     ! if(Any(abs(E_mat(:,i)) > 1.d99)) then
     !   print*, A, Y, f
     !   print*, a_sq, b_sq, N_abs, ey_sq
     !   stop "E_mat"
     ! end if
!    !end do
!  !if(Y > 0.498 .and. Y < 0.502 .and. sqrt(X) > 0.3) then
!  !  print*,sqrt(E_mat(1,1)),sqrt(E_mat(2,2)),sqrt(E_mat(3,3)),sqrt(E_mat(1,1) + E_mat(2,2) + E_mat(3,3))
!  !end if
!  end subroutine get_E_factors_cmplx

!  subroutine get_pol_vec_cmplx(X, Y, N_abs, N_perp, N_par, f, a_sq, b_sq, pol_vec )
!    use constants,                  only: pi, e0, mass_e, eps0, c0
!    implicit none
!    ! Computeses the polarization vector, where as e_y is defined to be purely imaginary:
!    ! From Bornatici Review 1983 p. 1198 eq. 3.1.61 and 3.1.62
!    ! ex := A c ey ! c contains all complex factors
!    ! ez := B c' ey ! c' contains all complex factors
!    ! ey ey^* = 1 / (N sqrt(a_sq + b_sq))
!    real(rkind), intent(in) :: X, Y, N_abs, N_perp, N_par, f, a_sq, b_sq
!    complex(r8), dimension(:), intent(out) :: pol_vec
!    real(rkind)         :: A, B, ey_sq, re_norm
!    integer(ikind)      :: i
!    A = 1.d0 / Y * (1.d0 - (1 - Y**2) * f)
!    B = (N_par * N_perp) / ( 1.d0 - X -  N_perp**2) * A
!    ey_sq = 1.d0 / (N_abs *  sqrt(a_sq + b_sq))
!    pol_vec(1) = cmplx(-A, 0.d0)
!    pol_vec(2) = cmplx(0.d0,1.d0) ! define E_y as complex
!    pol_vec(3) = cmplx(B, 0.d0)
!    pol_vec(:) = pol_vec(:) * sqrt(ey_sq)
!    pol_vec(:) = pol_vec(:) / sqrt(abs(pol_vec(1))**2 + abs(pol_vec(2))**2 + abs(pol_vec(3))**2)
!  end subroutine get_pol_vec_cmplx

!function get_filter_transmittance_wrong(omega, X, Y, sin_theta, cos_theta, mode, x_vec, N_vec, B_vec, x_launch, pol_vec_ext)
!  ! THIS ROUTINE IS MOST LIKELY INCORRECT - DO NOT USE!
!  ! The wave pass through the filter perpendicularly. However, because of the optic setup the wave vector at separatrix are not
!  ! the same as the ones in the plasma. Hence, the polarization vector has to be rotated so that the components of the electric field lie
!  ! in the same plane as the polarizer.
!    use constants,                  only: pi, e0, mass_e, eps0, c0
!    use mod_ecfm_refr_utils,        only: sub_remap_coords
!    use mod_ecfm_refr_types,        only: output_level
!    implicit none
!    real(rkind), intent(in)               :: omega, X, Y, sin_theta, cos_theta
!    integer(ikind), intent(in)            :: mode
!    real(rkind), dimension(:), intent(in) :: x_vec, N_vec, B_vec, x_launch
!    complex(r8), dimension(:), intent(in), optional :: pol_vec_ext
!    real(rkind)                 :: get_filter_transmittance_wrong, N_abs_ray, B_abs, cos_phi, B_t, B_pol, &
!                                   sin_phi, N_abs, f, a_sq, b_sq, N_par, N_perp, temp, phi_filter, sigma
!    real(rkind), dimension(3,3) :: wave_rot_mat, pol_perp_rot_mat, wave_rot_mat_transp
!    real(rkind), dimension(3)   :: N_filter, norm_vec_N_rot_plane, pol_vec_real, &
!                                   pol_vector_lab, pol_vec_perp, R_vec, N_vec_perp_test, &
!                                   e_x, e_y, e_z, N_vec_norm ! unit vectors of the coordinate system where the polarization vector is defined
!    complex(r8), dimension(3)   :: pol_vec
!    real(rkind), dimension(2)   :: Jones_vector, filtered_Jones_vector
!    real(rkind), dimension(2,2) :: filter_mat
!    integer(ikind)              :: i, j
!    logical                     :: debug
!    debug = .true.
!    !debug = .false.
!    if(present(pol_vec_ext)) then
!      pol_vec = pol_vec_ext
!    else
!      call abs_Al_N_with_pol_coeff( omega, X, Y, sin_theta, cos_theta, mode, N_abs, f, a_sq, b_sq)
!      if(N_abs == 0.d0) then
!        get_filter_transmittance = 0.d0 ! cut off - polarization vector undefined
!        return
!      end if
!      N_perp = N_abs * sin_theta
!      N_par = N_abs * cos_theta
!      call get_pol_vec_cmplx(X, Y, N_abs, N_perp, N_par, f, a_sq, b_sq, pol_vec )
!    end if
!    pol_vec_real(1) = real(pol_vec(1))
!    pol_vec_real(2) = aimag(pol_vec(2))
!    pol_vec_real(3) = real(pol_vec(3))
!    if(present(pol_vec_ext)) then
!      if(any(pol_vec_real /= pol_vec_real) .and. output_level) then
!        print*, "Something wrong with the external pol_coeff"
!        print*, pol_vec
!      end if
!    else
!      if(any(pol_vec_real /= pol_vec_real)) then
!        print*, "X", X
!        print*, "Y", Y
!        print*, "N_abs", N_abs
!        print*, "N_perp", N_perp
!        print*, "f", f
!        print*, "omega", omega
!        call abort()
!      end if
!    end if
!    if(debug) print*, "pol_vec", pol_vec_real(:)
!    ! First we rotate the polarization vector into the carthesian coordinate system of the machine
!    ! Now calculate rotation from Cartesian to reference frame of the wave and the polarization vector.
!    ! The polarization coefficients are given in a way so that
!    ! vec(N) = (vec(e_x) * sin(theta) + vec(e_z) * cos(theta)) * N_abs
!    ! vec(B) = B e_z
!    ! Hence, the unit vectors are given by:
!    ! vec(e_x) = vec(N) / N_abs - cos(theta) * vec(B) / B_abs
!    ! vec(e_y) = vec(e_x) x vec(e_z)
!    ! vec(e_z) = vec(B) / B_abs
!    ! The rotation is then given according to https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotations_and_motions
!    N_abs_ray = sqrt(sum(N_vec**2))
!    N_vec_norm = N_vec / N_abs_ray
!    B_abs = sqrt(sum(B_vec**2))
!    if(debug) then
!      print*, "X", X
!      print*, "Y", Y
!      print*, "omega", omega
!      print*, "x_vec", x_vec
!      print*, "N_vec", N_vec
!      print*, "B_vec", B_vec
!      print*, "theta", acos(cos_theta) / pi * 180.d0
!    end if
!    if(mode > 0 .and. debug) print*, "X-mode"
!    if(mode < 0 .and. debug) print*, "O-mode"
!    ! N_vec already points towards the antenna when its is copied into svec
!    e_x = N_vec_norm - cos_theta  * B_vec / B_abs
!    e_x = e_x / sqrt(sum(e_x**2))
!    e_z = B_vec/B_abs
!    ! e_y = e_z x e_x
!    e_y(1) = e_z(2) * e_x(3) - e_z(3) * e_x(2)
!    e_y(2) = e_z(3) * e_x(1) - e_z(1) * e_x(3)
!    e_y(3) = e_z(1) * e_x(2) - e_z(2) * e_x(1)
!    e_y(:) = e_y(:) / sqrt(sum(e_y**2)) ! not necessary because e_x and e_z perpendicular
!    if(debug) then
!      print*, "e_x", e_x
!      print*, "e_y", e_y
!      print*, "e_z", e_z
!      print*, "e_x . e_y", sum(e_x * e_y)
!      print*, "e_x . e_z", sum(e_x * e_y)
!      print*, "e_y . e_z", sum(e_y * e_z)
!      print*, "e_x.N_vec", sum(e_x * N_vec_norm)
!      print*, "e_y.N_vec", sum(e_y * N_vec_norm)
!      print*, "e_z.N_vec", sum(e_z * N_vec_norm)
!    end if
!    pol_vector_lab(:) = pol_vec_real(1) * e_x + &
!                        pol_vec_real(2) * e_y + &
!                        pol_vec_real(3) * e_z
!    ! Now remove portion that points along N_vec
!    pol_vector_lab(:) = pol_vector_lab(:) - N_vec_norm(:) * sum(N_vec_norm(:) * pol_vector_lab(:))
!    pol_vector_lab(:) = pol_vector_lab(:) / sqrt(sum(pol_vector_lab(:)**2))
!    if(debug) print*, "Polvec in laboratory frame", pol_vector_lab(:)
!    if(debug) print*, "Dot product N_vec pol vec in lab frame", sum(N_vec * pol_vector_lab)
!    ! Next the rotation of k by the quasi-optical system:
!    ! Normalized vector perpendicular to the filter
!    N_filter(:) = x_launch(:)
!    ! We do not want a z component here
!    N_filter(3) = 0.d0
!    N_filter(:) = N_filter(:) / sqrt(sum(N_filter**2))
!    if(debug) print*, "N_vec norm", N_vec_norm
!    if(debug) print*, "N_filter", N_filter
!    ! Rotation matrix around angle sigma with axis N_filter x N_vec
!    sigma = -acos(sum(N_vec_norm * N_filter))
!    if(debug) print*, "Sigma [deg.]", sigma / pi * 180.d0
!    norm_vec_N_rot_plane(:) = 0.d0
!    ! Rotate the polarization vector in the plane spanned by N_filter and N_vec_norm by the angle sigma
!    norm_vec_N_rot_plane(1) = N_filter(2) * N_vec_norm(3) - N_filter(3) * N_vec_norm(2)
!    norm_vec_N_rot_plane(2) = N_filter(3) * N_vec_norm(1) - N_filter(1) * N_vec_norm(3)
!    norm_vec_N_rot_plane(3) = N_filter(1) * N_vec_norm(2) - N_filter(2) * N_vec_norm(1)
!    norm_vec_N_rot_plane(:) = norm_vec_N_rot_plane / sqrt(sum(norm_vec_N_rot_plane**2))
!    if(debug) print*, "Axis of rotation", norm_vec_N_rot_plane
!    ! First index selects column second index row
!    pol_perp_rot_mat(1,1) = norm_vec_N_rot_plane(1)**2 + (norm_vec_N_rot_plane(2)**2 + &
!                            norm_vec_N_rot_plane(3)**2)*Cos(sigma)
!    pol_perp_rot_mat(2,1) = norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(2) - &
!                            norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(2)*Cos(sigma) - &
!                            norm_vec_N_rot_plane(3)*Sin(sigma)
!    pol_perp_rot_mat(3,1) = norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(3) - &
!                            norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(3)*Cos(sigma) + &
!                            norm_vec_N_rot_plane(2)*Sin(sigma)
!    pol_perp_rot_mat(1,2) = norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(2) - &
!                            norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(2)*Cos(sigma) + &
!                            norm_vec_N_rot_plane(3)*Sin(sigma)
!    pol_perp_rot_mat(2,2) = norm_vec_N_rot_plane(2)**2*(1 - Cos(sigma)) + Cos(sigma)
!    pol_perp_rot_mat(3,2) = norm_vec_N_rot_plane(2)*norm_vec_N_rot_plane(3) - &
!                            norm_vec_N_rot_plane(2)*norm_vec_N_rot_plane(3)*Cos(sigma) - &
!                            norm_vec_N_rot_plane(1)*Sin(sigma)
!    pol_perp_rot_mat(1,3) = norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(3) - &
!                            norm_vec_N_rot_plane(1)*norm_vec_N_rot_plane(3)*Cos(sigma) - &
!                            norm_vec_N_rot_plane(2)*Sin(sigma)
!    pol_perp_rot_mat(2,3) = norm_vec_N_rot_plane(2)*norm_vec_N_rot_plane(3) - &
!                            norm_vec_N_rot_plane(2)*norm_vec_N_rot_plane(3)*Cos(sigma) + &
!                            norm_vec_N_rot_plane(1)*Sin(sigma)
!    pol_perp_rot_mat(3,3) = norm_vec_N_rot_plane(3)**2*(1 - Cos(sigma)) + Cos(sigma)
!    if(debug) print*, "Polarizer rotation matrix"
!    do i = 1,3
!      pol_vec_perp(i) = sum(pol_perp_rot_mat(:,i) * pol_vector_lab(:))
!      if(debug) print*, pol_perp_rot_mat(:,i)
!    end do
!    do i = 1,3
!      N_vec_perp_test(i) = sum(pol_perp_rot_mat(:,i) * N_vec_norm(:))
!    end do
!    if(debug) print*, "pol vec perp filter", pol_vec_perp
!    if(debug) print*, "dot product rotated N_vec and N_filter vec - should be one", sum(N_vec_perp_test * N_filter)
!    if(debug) print*, "dot product rotated polarization vector and N_filter vec - should be 0", sum(pol_vec_perp * N_filter)
!    call sub_remap_coords(x_vec, R_vec)
!    cos_phi = cos(R_vec(2))
!    sin_phi = sin(R_vec(2))
!    Jones_vector(1) = sin_phi * pol_vec_perp(1) + cos_phi *  pol_vec_perp(2)
!    Jones_vector(2) = pol_vec_perp(3)
!    if(debug) print*, "Jones_vector:" , Jones_vector
!    filter_mat(:,:) = 0.d0
!    filter_mat(2,2) = 1.d0
!    do i = 1,2
!      filtered_Jones_vector(i) = sum(filter_mat(i,:) * Jones_vector(:))
!    end do
!    if(debug) print*, "filtered Jones_vector:" , filtered_Jones_vector
!    get_filter_transmittance = sum(filtered_Jones_vector)**2
!    if(debug) print*, "Transmittance", get_filter_transmittance
!    if(get_filter_transmittance /= get_filter_transmittance .or. get_filter_transmittance < 0.0 .or. &
!       get_filter_transmittance > 1.d0) then
!      print*, "pol_vec in original coordinates", pol_vec
!      print*, "pol_vec in carth. coordinates", pol_vec
!      print*, "phi", R_vec(2)
!      print*, "cos_theta", cos_theta
!      print*, "X", X
!      print*, "Y", Y
!      print*, "N_abs", N_abs
!      print*, "Transmittance", get_filter_transmittance
!      stop "Bad value in calculation of polarization filter rejection"
!    end if
!  end function get_filter_transmittance_wrong!

! mod mod_ecfm_refr_dist


function radiation_dist_f_beta(beta, beta_perp, Te, omega_bar, N_par, const)
    use constants,                  only: pi, e0
    use mod_ecfm_refr_types,        only: dstf
    implicit none
    real(rkind), intent(in)            :: beta, beta_perp,Te, omega_bar, N_par, const
    real(rkind)                        :: radiation_dist_f_beta
    real(rkind)                        :: mu, gamma, beta_total
    !real(rkind)     :: du_dbeta, du_dbeta_perp, int_beta_add, beta_perp_u, jacobian
    if(dstf == "maxwell") then
      radiation_dist_f_beta =  Exp(-const * (beta**2 + beta_perp**2))
    else if(dstf == "relamax" .or. dstf == "numeric") then
      !beta_total = sqrt(beta**2 + beta_perp**2)
      gamma = 1 / sqrt(1 - beta**2 - beta_perp**2)
      mu = 2 * const
      !radiation_f =  beta_perp* exp( mu * (1 - 1 / (sqrt(1 - beta**2 - beta_perp**2)) )) * gamma**0
      radiation_dist_f_beta =  exp( mu* (1 - gamma))*gamma**5
      !if(abs(radiation_f) > 1.d-2) then
      !  print*,mu, beta, f_rel, N_par
      !  stop "value"
      !end if
    !else if (dstf == "numeric") then
    !  stop "not implemented"
    end if
  end function radiation_dist_f_beta

function radiation_dist_f_norm_beta(const,Te)
    use constants,                  only:  pi, e0, mass_e, c0
    use mod_ecfm_refr_types,        only: dstf
    implicit none
    real(rkind), intent(in)        :: const, Te
    real(rkind)                    :: radiation_dist_f_norm_beta
    real(rkind)                    :: a, mu
    radiation_dist_f_norm_beta = -1.d0
    if(dstf == "maxwell") then
      radiation_dist_f_norm_beta = sqrt(const/(pi*c0**2))**3
    else if(dstf == "relamax" .or. dstf == "numeric") then
      mu = const*2
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      radiation_dist_f_norm_beta = a * (sqrt(mu / (2.0d0 * pi))**3) ! *c0**2c0)
    else
      stop "the distribution defined in globals does not exist"
    end if
  end function radiation_dist_f_norm_beta


