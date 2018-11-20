! All the deprecated source code ends up here.
! There individual code chunks are sorted according to their original module
!!!! mod_ecfm_refr.raytrayce_initialize !!!!
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
