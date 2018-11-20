module mod_ecfm_refr_raytrace_initialize
    use f90_kind
    use mod_ecfm_refr_types,       only: spl_type_2d
    implicit none
    ! Two interfaces to allow interpolation of Te and ne using both single rhop values and vectors of rhop values
    public :: init_raytrace
    type(spl_type_2d)                       :: Psi_spline
    contains

  subroutine get_psi_ax(R, z, i_ax, j_ax, psi_sep, psi_ax, R_ax, z_ax)
  ! Interpolates
    use f90_kind
    use mod_ecfm_refr_types,       only: output_level
    use nr_mod, only                        : powell
    implicit none
    real(rkind), dimension(:), intent(in)  :: R, z
    integer(ikind)           , intent(in)  :: i_ax, j_ax
    real(rkind),               intent(in)  :: psi_sep
    real(rkind),               intent(out) :: psi_ax
    real(rkind),               intent(out) :: R_ax
    real(rkind),               intent(out) :: z_ax
    real(rkind), dimension(2)              :: p, p_h
    real(rkind), dimension(2,2)            :: xi
    real(rkind)                            :: h, ftol
    integer(ikind)                         :: iter
    p(1) = R(i_ax)
    p(2) = z(j_ax)
    h = 1.d-5
    p_h  = p
    p_h(1) = p_h(1) + h
    xi(:,1) = (eval_Psi(p_h) - eval_Psi(p)) / h
    p_h  = p
    p_h(2) = p_h(2) + h
    xi(:,2) = (eval_Psi(p_h) - eval_Psi(p)) / h
    ftol = 1.d-6
    !print*, "Before opt",  p(1), p(2)
    call powell(p,xi,ftol,iter, psi_ax, eval_Psi)
    R_ax = p(1)
    z_ax = p(2)
    if(output_level) print*, "Magnetic axis position and rhop on axis",  p(1), p(2), sqrt((eval_Psi(p) - psi_ax)/(psi_sep - psi_ax))
  end subroutine get_psi_ax

  function eval_Psi(R)
    use f90_kind
    use mod_ecfm_refr_interpol,    only: rect_spline
    implicit none
    real(rkind), dimension(:), intent(in)   :: R
    real(rkind)                             :: eval_Psi
    call rect_spline(Psi_spline, R(1), R(2), eval_Psi)
  end function eval_Psi

  subroutine read_topfile(plasma_params, R, z, rhop, Br, Bt, Bz, R0, Btf0, itime)
      use f90_kind
      use mod_ecfm_refr_types,       only: plasma_params_type, data_folder
      use mod_ecfm_refr_interpol,    only: make_rect_spline
#ifdef OMP
      use omp_lib
#endif
      implicit none
      type(plasma_params_type) , intent(inout) :: plasma_params
      real(rkind), dimension(:), allocatable, intent(out) :: R, z
      real(rkind), dimension(:,:), allocatable, intent(out) ::  Br, Bt, Bz
      real(rkind), dimension(:,:), allocatable, intent(out) ::  rhop
      real(rkind)                            , intent(out)  :: R0, Btf0
      integer(ikind), intent(in), optional                  :: itime
      character(250) :: filename
      integer(ikind) :: i, j, k, i_ax, j_ax
      real(rkind) :: radin, radout, sgn
      character(70) :: line
      if(present(itime)) then
        write(filename, fmt = "(A7I5.5)") "topfile",  itime
        filename = trim(data_folder) // trim(filename)
      else
        filename =  trim(data_folder) // "topfile"
      end if
      open(70,status='unknown',file=filename)
      read(70,'(A)') line
      read(70,*) plasma_params%m, plasma_params%n
      !print*, "Dimensions of external data: ", plasma_params%m, plasma_params%n
      allocate(R(plasma_params%m), z(plasma_params%n), Br(plasma_params%m,plasma_params%n), &
               Bt(plasma_params%m,plasma_params%n), Bz(plasma_params%m,plasma_params%n), &
               rhop(plasma_params%m, plasma_params%n))
      read(70,'(A)') line
      read(70,*) radin,radout, plasma_params%pf_sxp
! ... Grid: X-coordinates
      read(70,'(A)') line  ;  read(70,fmt=*)(R(i),i=1,plasma_params%m)

! ... Grid: Z-coordinates
      read(70,'(A)') line  ;  read(70,fmt=*)(z(j),j=1,plasma_params%n)

! ... Magnetic field: B_X, B_Y (toroidal), B_z (vertical)
      read(70,'(A)') line  ;  read(70,fmt=*)((Br(i,j),i=1,plasma_params%m),j=1,plasma_params%n)
      read(70,'(A)') line  ;  read(70,fmt=*)((Bt(i,j),i=1,plasma_params%m),j=1,plasma_params%n)
      read(70,'(A)') line  ;  read(70,fmt=*)((Bz(i,j),i=1,plasma_params%m),j=1,plasma_params%n)

!     Poloidal flux: psi
      read(70,'(A)') line  ;  read(70,fmt=*)((rhop(i,j),i=1,plasma_params%m),j=1,plasma_params%n)

      close(70)
      if (rhop(2, plasma_params%n/2) > rhop(plasma_params%m/2, plasma_params%n/2)) then
         sgn = 1.e0_rkind
      else
         sgn = -1.e0_rkind
      endif
      plasma_params%pf_mag = 1.d99
      i_ax = -1
      do j = 1, plasma_params%n
        rhop(:,j) = rhop(:,j) * sgn
        if(plasma_params%pf_mag > minval(rhop(:,j))) then
           plasma_params%pf_mag = minval(rhop(:,j))
          i_ax = minloc(rhop(:,j), dim = 1)
          j_ax = j
        end if
      end do
      if(i_ax == -1) then
        print*, "Could not find magnetic axis in read_topfile in mod_raytrace_initialize.f90"
        print*, "Check if psi matrix dimensions correct: R, z", plasma_params%m, plasma_params%n
        stop "Input error - check topifle"
      end if
      R0 = R(i_ax)
      Btf0 = Bt(i_ax, j_ax)
      plasma_params%Btf0 = Btf0
      plasma_params%R0 = R0
      !print*, "Retrieving magn axis"
      call make_rect_spline(Psi_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), R, z, rhop)
      call get_psi_ax(R, z, i_ax, j_ax, plasma_params%pf_sxp, plasma_params%pf_mag, plasma_params%R_ax, plasma_params%z_ax)
      rhop = sqrt((rhop - plasma_params%pf_mag)/(plasma_params%pf_sxp - plasma_params%pf_mag))
      if(any(rhop /= rhop)) then
        print*, "critical error when reading topfile, NAN in rho_pol found"
        print*, rhop
        stop "sub read topfile in mod_raytrace_intialize"
      end if
  end subroutine read_topfile

  subroutine read_Te_ne_matrix(plasma_params, R, z, Te, ne)
      use f90_kind
      use mod_ecfm_refr_types,       only: plasma_params_type, data_folder, min_density, max_Te, min_Te
      implicit none
      type(plasma_params_type) , intent(inout) :: plasma_params
      real(rkind), dimension(:), allocatable, intent(in) :: R, z
      real(rkind), dimension(:,:), allocatable, intent(out) ::  Te, ne
      real(rkind), dimension(size(R))                       :: R_check
      real(rkind), dimension(size(z))                       :: z_check
      character(250) :: filename
      integer(ikind) :: m_check, n_check, i,j,k
      character(70) :: line
      filename =  trim(data_folder) // "Te_ne_matfile"
      open(71,status='unknown',file=filename)
      read(71,'(A)') line
      read(71,*) m_check, n_check
      if(m_check /= plasma_params%m .or. n_check /= plasma_params%n) then
        print*, "Te and ne grid input is enabled"
        print*, "The Te and ne grid exact same dimension as the magnetic field and flux matrix"
        print*, "dimensions of flux matrix", plasma_params%m, plasma_params%n
        print*, "dimensions of Te/ne matrix", m_check, n_check
        stop "Input Error! Check input and run again"
      end if
      allocate(Te(plasma_params%m,plasma_params%n), &
               ne(plasma_params%m,plasma_params%n))
! ... Grid: X-coordinates
      read(71,'(A)') line  ;  read(71,fmt=*)(R_check(i),i=1,plasma_params%m)
      if(.not. all(R == R_check)) then
        print*, "Te and ne grid input is enabled"
        print*, "R grid of Te/ne is not identical to R grid of flux matrix"
        stop "Input Error! Check input and run again"
      end if
! ... Grid: Z-coordinates
      read(71,'(A)') line  ;  read(71,fmt=*)(z_check(j),j=1,plasma_params%n)
      if(.not. all(z ==z_check)) then
        print*, "Te and ne grid input is enabled"
        print*, "z grid of Te/ne is not identical to z grid of flux matrix"
        stop
      end if
! ... Te
      read(71,'(A)') line  ;  read(71,fmt=*)((Te(i,j),i=1,plasma_params%m),j=1,plasma_params%n)
      if(any(Te > max_Te) .or. any(Te < 0) .or. all(Te < min_Te)) then
        print*, "Te and ne grid input is enabled"
        if(any(Te > max_Te)) then
          print*, "Encoutered Te >", max_Te
          print*, "Largest Te", maxval(Te)
        else if(any(Te < 0.d0)) then
          print*, "Encoutered  negative values for Te"
          print*, "Negative Te", minval(Te)
        else if(all(Te < min_Te)) then
          print*, "All Te smaller than the minium Te: ", min_Te
        end if
        stop "Input Error! Check input and run again"
      end if
! ... ne
      read(71,'(A)') line  ;  read(71,fmt=*)((ne(i,j),i=1,plasma_params%m),j=1,plasma_params%n)
      if( all(ne < min_density)) then
        print*, "Te and ne grid input is enabled"
        print*, "All density points are smaller than:", min_density
        stop "Input Error! Check input and run again"
      end if
  end subroutine read_Te_ne_matrix
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
  subroutine make_topfile(R, z, rhop, Br, Bt, Bz, working_dir, itime)
    use f90_kind
    implicit none
    real(rkind), dimension(:), allocatable, intent(in) :: R, z
    real(rkind), dimension(:,:), allocatable, intent(in) ::  Br, Bt, Bz
    real(rkind), dimension(:,:), allocatable, intent(in) ::  rhop
    character(*), intent(in)                             :: working_dir
    integer(ikind)                                       ::  itime
    character(250) :: filename
    integer(ikind) :: i,j,k
    real(rkind) :: radin,radout,psi_sep,sgn, psi_ax
    character(70) :: line
    character(12) :: format_str
    write(filename, fmt = "(A7I5.5)") "topfile",  itime
    filename = trim(working_dir) // trim(filename)
    open(70,status='unknown',file=filename)
    write(70,'(A)') "Number of radial and vertical grid points"
    write(70,*) size(R), size(z)
    write(70,'(A)') "Inside and Outside radius and psi_sep"
    write(70,*) R(1), R(size(R)), 1.0 ! since we use psi = rhop**2
! ... Grid: X-coordinates
    write(70,'(A)') "Radial grid coordinates"
    write(70,fmt=*)(R(i),i=1,size(R))

! ... Grid: Z-coordinates
    write(70,'(A)') "Vertical grid coordinates"
    write(70,fmt=*)(z(j),j=1,size(z))

! ... Magnetic field: B_X, B_Y (toroidal), B_z (vertical)
    write(70,'(A)') "B_r values"  ;  write(70,fmt=*)((Br(i,j),i=1,size(R)),j=1,size(z))
    write(70,'(A)') "B_t values"  ;  write(70,fmt=*)((Bt(i,j),i=1,size(R)),j=1,size(z))
    write(70,'(A)') "B_z values"  ;  write(70,fmt=*)((Bz(i,j),i=1,size(R)),j=1,size(z))

!     Poloidal flux: psi
    write(70,'(A)') "psi values"  ;  write(70,fmt=*)((rhop(i,j)**2,i=1,size(R)),j=1,size(z))
    close(70)
  end subroutine make_topfile

!  subroutine extend_profile(rhop_vec, profile, m, step, rhop_max, limiting_value, rhop_vec_ext, profile_ext)
!  ! Extends profile linearly so that it is close to limiting_value at rhop_max
!  use f90_kind
!  implicit none
!  real(rkind), dimension(:),  intent(in) :: rhop_vec, profile
!  real(rkind), dimension(:),  allocatable, intent(out) :: rhop_vec_ext, profile_ext
!  integer(ikind), intent(inout)                         :: m
!  real(rkind), intent(in)                               :: rhop_max, step, limiting_value
!  integer(ikind)                                        :: i, irhop_extra
!    if(rhop_max < 1.05) then
!      print*, "Largest rho poloidal found is smaller than 1.05 - aborting"
!      stop "Please check flux matrix (shot file or external)"
!    end if
!    irhop_extra = (rhop_max - rhop_vec(m)) / (step) + 1
!    if(irhop_extra > 0) then
!      allocate(rhop_vec_ext(m + irhop_extra), profile_ext(m + irhop_extra))
!      rhop_vec_ext(1: m) = rhop_vec
!      profile_ext(1: m) = profile
!      do i = m, m + irhop_extra
!        rhop_vec_ext(i) = rhop_vec_ext(i - 1) + step
!        profile_ext(i) = profile_ext(m) * &
!           real(m + irhop_extra - i, 8) / &
!           real(m + irhop_extra, 8) + limiting_value
!      end do
!    else
!      irhop_extra = 0
!      allocate(rhop_vec_ext(m), profile_ext(m))
!      rhop_vec_ext = rhop_vec
!      profile_ext = profile_ext
!    end if
!    m = m + irhop_extra
!  end subroutine extend_profile

  subroutine update_plasma_params(plasma_params, par, par_ne, par_scal)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type
#ifdef NAG
    USE nag_spline_1d,              only: nag_spline_1d_interp
#endif
    implicit none
    type(plasma_params_type), intent(inout) :: plasma_params
    real(rkind), dimension(:), intent(in)   :: par, par_ne, par_scal
    if(.not. allocated(plasma_params%par)) then
      if(size(par) == 0) then
        print*, "The vector par containing the fit parameters in ida was not allocated when update_plasma_params was called"
        stop "Critical error in update_plasma_params in mod_ecfm_refr_raytrace_initialize.f90"
      end if
      if(size(par_scal) == 0) then
        print*, "The vector par_scal containing the scaling of the fit parameters in ida was not allocated when update_plasma_params was called"
        stop "Critical error in update_plasma_params in mod_ecfm_refr_raytrace_initialize.f90"
      end if
      if(size(par_ne) == 0) then
        print*, "The vector par_ne containing static fit parameters for ne in ida was not allocated when update_plasma_params was called"
        stop "Critical error in update_plasma_params in mod_ecfm_refr_raytrace_initialize.f90"
      end if
      allocate(plasma_params%par(size(par)), plasma_params%par_ne(size(par_ne)), plasma_params%par_scal(size(par_scal)))
    end if
    plasma_params%par = par
    plasma_params%par_ne = par_ne
    plasma_params%par_scal = par_scal
  end subroutine update_plasma_params

  subroutine setup_plasma_params(plasma_params,eq_from_shotfile)
    use f90_kind
    use mod_ecfm_refr_types,        only: plasma_params_type, h_x_glob, h_check, output_level, double_check_splines, &
                                          stand_alone, output_level, vessel_bd_filename, max_points_svec
    use mod_ecfm_refr_interpol,     only: make_rect_spline, make_1d_spline, rect_spline
    use constants,                  only: pi
    use quadrature,                 only: cdgqf
#ifdef NAG
    USE nag_spline_2d,              only: nag_spline_2d_interp, &
    									                    nag_error, nag_set_error, nag_spline_2d_eval
    USE nag_spline_1d,              only: nag_spline_1d_interp
    USE nag_error_handling
#endif
    use ripple3d,                   only: init_ripple!, validate_ripple_grad
    implicit none
    type(plasma_params_type), intent(inout)     :: plasma_params
    logical, intent(in)                         :: eq_from_shotfile
    real(rkind), dimension(:,:), allocatable    :: B_r, B_t, B_z, T_e, n_e
    integer(ikind), dimension(:), allocatable   :: R_index_lower, z_index_lower, R_index_upper, z_index_upper
    real(rkind), dimension(:), allocatable      :: rhop_Te, rhop_ne, ne_ext, Te_ext
    real(rkind)                                 :: R_mag, z_mag, pf_mag, R_sxp, z_sxp, pf_sxp, rhop_dx_dummy, &
                                                   rhop_dy_dummy, rhop_dxx_dummy, rhop_dxy_dummy, rhop_dyy_dummy, &
                                                   R0, Btf0, rhop_X, B_last
    integer(ikind)                              :: i, j
    character(20) :: format_str
    character(1)  :: sep
    if(.not. eq_from_shotfile) then
      call read_topfile(plasma_params, plasma_params%R, plasma_params%z, plasma_params%rhop, B_r, B_t, B_z, R0, Btf0)
      if(plasma_params%Te_ne_mat) then
        call read_Te_ne_matrix(plasma_params, plasma_params%R, plasma_params%z, T_e, n_e)
        plasma_params%rhop_max = plasma_params%rhop_entry !
      end if
    else
      call setup_eq_from_shotfile(plasma_params, &
                                  plasma_params%R, plasma_params%z, plasma_params%rhop, B_r, B_t, B_z, R0, Btf0)
    end if
    plasma_params%m = size(plasma_params%R)
    plasma_params%n = size(plasma_params%z)
!    do j = 1, plasma_params%n
!      print*, B_t(:,j)
!    end do
!    stop "B_t weird"
    h_check = h_x_glob * (1.d0 + 1.e-5)
    plasma_params%R_min = plasma_params%R(1)
    plasma_params%R_max = plasma_params%R(plasma_params%m)
    plasma_params%z_min = plasma_params%z(1)
    plasma_params%z_max = plasma_params%z(plasma_params%n)
    plasma_params%R_step = (plasma_params%R_max - plasma_params%R_min) / real(plasma_params%m - 1,8)
    plasma_params%z_step = (plasma_params%z_max - plasma_params%z_min) / real(plasma_params%n - 1,8)
    if(.not. plasma_params%Te_ne_mat .and. stand_alone) then
#ifdef NAG
      if(double_check_splines .and. output_level) then
        call nag_spline_1d_interp(plasma_params%rhop_vec_ne, plasma_params%n_e_prof, plasma_params%ne_spline_nag)
        call nag_spline_1d_interp(plasma_params%rhop_vec_Te, plasma_params%T_e_prof, plasma_params%Te_spline_nag)
      end if
#endif
      call make_1d_spline( plasma_params%ne_spline, int(size(plasma_params%rhop_vec_ne),4), plasma_params%rhop_vec_ne, plasma_params%n_e_prof)
      call make_1d_spline( plasma_params%Te_spline, int(size(plasma_params%rhop_vec_Te),4), plasma_params%rhop_vec_Te, plasma_params%T_e_prof)
    end if
    call init_ripple(R0, Btf0)
    call make_rect_spline(plasma_params%rhop_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, plasma_params%rhop)
    call make_rect_spline(plasma_params%B_r_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_r)
    call make_rect_spline(plasma_params%B_t_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_t)
    call make_rect_spline(plasma_params%B_z_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_z)
    if(plasma_params%Te_ne_mat) then
      call make_rect_spline(plasma_params%T_e_spline_2D, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, log(T_e))
      ! Warning: ne is scaled down to avoid numerical problems!
      call make_rect_spline(plasma_params%n_e_spline_2D, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, log(n_e * 1.e-19))
#ifdef NAG
      if(double_check_splines .and. output_level) then
        call nag_spline_2d_interp(plasma_params%R, plasma_params%z, log(T_e), plasma_params%Te_spline_nag_2D)
        call nag_spline_2d_interp(plasma_params%R, plasma_params%z, log(n_e * 1.e-19), plasma_params%ne_spline_nag_2D)
      end if
#endif
    end if
#ifdef NAG
    if(output_level .and. double_check_splines) then
      call nag_spline_2d_interp(plasma_params%R, plasma_params%z, plasma_params%rhop, plasma_params%rhop_spline_nag)
      if(eq_from_shotfile) then
        call nag_spline_2d_eval(plasma_params%rhop_spline_nag, plasma_params%R_sep, plasma_params%z_sep, rhop_X)
        print*,"Position of X point and rho pol value from spline", plasma_params%R_sep, plasma_params%z_sep, rhop_X
      end if
      call nag_spline_2d_interp(plasma_params%R, plasma_params%z, B_r, plasma_params%B_r_spline_nag)
      call nag_spline_2d_interp(plasma_params%R, plasma_params%z, B_t, plasma_params%B_t_spline_nag)!_cor
      call nag_spline_2d_interp(plasma_params%R, plasma_params%z, B_z, plasma_params%B_z_spline_nag)
    end if
#endif
    plasma_params%B_ax = 0.d0
    call rect_spline(plasma_params%B_r_spline,plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    call rect_spline(plasma_params%B_t_spline,plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    call rect_spline(plasma_params%B_z_spline,plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    plasma_params%B_ax = Sqrt(plasma_params%B_ax)
    allocate(R_index_lower(plasma_params%n), z_index_lower(plasma_params%m), R_index_upper(plasma_params%n), z_index_upper(plasma_params%m))
    R_index_lower(:) = 1
    z_index_lower(:) = 1
    R_index_upper(:) = plasma_params%m
    z_index_upper(:) = plasma_params%n
    plasma_params%int_step_cnt = plasma_params%rad_trans_sections * plasma_params%rad_transp_solver_order
    allocate( plasma_params%Int_weights(plasma_params%int_step_cnt),plasma_params%Int_absz(plasma_params%int_step_cnt))
    if(plasma_params%on_the_fly_raytracing) then
      call cdgqf( int(plasma_params%int_step_cnt,kind=4), int(1,kind=4), 0.d0, 0.d0, plasma_params%Int_weights, plasma_params%Int_absz)
    else
      do i = 1, plasma_params%int_step_cnt
        plasma_params%Int_absz(i) = real(i ,8) / real(plasma_params%int_step_cnt ,8) ! spans the points in a way so that there is no overlap
                                                                                     ! Int_absz(1) > 0.d0
      end do
    end if
!    open(80,file = "Te_mat.dat")
!    write(format_str, fmt = "(A1,I4,A13)"), "(", plasma_params%m, "(E13.6E2,A1))"
!    !print*, ormat_str
!    do j = 1,plasma_params%n
!      write(80, fmt = format_str),(T_e(i,j), " ", i = 1, plasma_params%m)
!    end do
!    close(80)
    deallocate(B_r, B_z, B_t, R_index_lower,z_index_lower,R_index_upper,z_index_upper)!, B_t
    if(plasma_params%Te_ne_mat) deallocate(T_e, n_e)
    ! Load the polygon describing the vessel wall
    open(66, file = vessel_bd_filename)
    read(66, "(I7.7)") plasma_params%m_vessel_bd
    allocate(plasma_params%vessel_poly(plasma_params%m_vessel_bd))
    do i = 1, plasma_params%m_vessel_bd
      read(66,"(E19.12E2A1E19.12E2)") plasma_params%vessel_poly(i)%x, sep,  plasma_params%vessel_poly(i)%y
    end do
    close(66)

!    ! Checks if interpolation values are available for the forward derivatives
  end subroutine setup_plasma_params

  subroutine read_input_lists(plasma_params)
  ! Read the input data n_e, T_e and wall file
    use mod_ecfm_refr_types,       only: plasma_params_type, n_e_filename, T_e_filename, &
                                          vessel_bd_filename, min_density, min_Te, max_Te
    use f90_kind
    use constants, only               : pi
    implicit none
    type(plasma_params_type), intent(inout)           :: plasma_params
    integer(ikind)                                    :: i, j
    character(1)                                      :: sep
    logical                                           :: new_wg
    real(rkind)                                       :: rhop_max_Te, rhop_max_ne, last_rhop
    if(.not. plasma_params%Te_ne_mat) then
      open(66, file = n_e_filename)
      read(66, "(I7.7)") plasma_params%m_n_e_prof
      allocate(plasma_params%rhop_vec_ne(plasma_params%m_n_e_prof), &
               plasma_params%n_e_prof(plasma_params%m_n_e_prof))
      rhop_max_ne = -1.d0
      do i = 1, plasma_params%m_n_e_prof
        read(66,"(E19.12E2A1E19.12E2)") plasma_params%rhop_vec_ne(i), sep, plasma_params%n_e_prof(i)
        if(i > 1) then
          if(last_rhop > plasma_params%rhop_vec_ne(i)) then
            print*, "Error while reading the density profile from file"
            print*, "The rho poloidal has to be strictly increasing!"
            print*, "Current line", i + 1
            print*, "Rho pol from previous line and rho_pol from last line", last_rhop, plasma_params%rhop_vec_ne(i)
            call abort()
          end if
        end if
        last_rhop = plasma_params%rhop_vec_ne(i)
      end do
      if( all(plasma_params%n_e_prof < min_density)) then
        print*, "Error while reading the density profile from file"
        print*, "Not a single point in the density profile is larger than",  min_density
        print*, "The density profile must be given in m^-3"
        print*, "Plasmas with densities below",  min_density, " are not supported at the moment"
        call abort()
      end if
      close(66)
      open(66, file = T_e_filename)
      read(66, "(I7.7)") plasma_params%m_T_e_prof
      allocate(plasma_params%rhop_vec_Te(plasma_params%m_T_e_prof), &
               plasma_params%T_e_prof(plasma_params%m_T_e_prof))
      do i = 1, plasma_params%m_T_e_prof
        read(66,"(E19.12E2A1E19.12E2)") plasma_params%rhop_vec_Te(i), sep, plasma_params%T_e_prof(i)
        if(i > 1) then
          if(last_rhop > plasma_params%rhop_vec_Te(i)) then
            print*, "Error while reading the temperature profile from file"
            print*, "The rho poloidal has to be strictly increasing!"
            print*, "Current line", i + 1
            print*, "Rho pol from previous line and rho_pol from last line", last_rhop, plasma_params%rhop_vec_Te(i)
            call abort()
          end if
        end if
        last_rhop = plasma_params%rhop_vec_Te(i)
      end do
      if( all(plasma_params%T_e_prof < min_Te)) then
        print*, "Error while reading the temperature profile from file"
        print*, "Not a single point in the temperature profile is larger than", min_Te
        print*, "The temperature profile must be given in eV"
        print*, "Plasmas with electron temperatures below",  min_Te, " are not supported at the moment"
        call abort()
      end if
      if( any(plasma_params%T_e_prof > max_Te)) then
        print*, "Error while reading the temperature profile from file"
        print*, "A point of temperature profile is larger than", max_Te, " < ", maxval(plasma_params%T_e_prof)
        print*, "The temperature profile must be given in eV"
        print*, "Plasmas with electron temperatures higher than",  max_Te, " are not supported at the moment"
        call abort()
      end if
      rhop_max_Te = maxval(plasma_params%rhop_vec_Te, dim = 1)
      rhop_max_ne = maxval(plasma_params%rhop_vec_ne, dim = 1)
      plasma_params%rhop_max = min(rhop_max_te, rhop_max_ne)
      close(66)
      if(plasma_params%m_T_e_prof < 40 .or. plasma_params%m_n_e_prof < 40) then
        print*, "points for density and Te profile:", plasma_params%m_T_e_prof
        stop "both density and temperature profile have to have MORE than 40 points"
      end if
    end if
  end subroutine read_input_lists

  subroutine init_raytrace(plasma_params, flag)
    use f90_kind
    use mod_ecfm_refr_types,       only: plasma_params_type, stand_alone
    Use constants                , only : pi
    implicit none
    type(plasma_params_type), intent(inout)                     :: plasma_params
    character(*), intent(in), optional :: flag
    logical                                         :: eq_from_shotfile
    if(present(flag)) then
      if(flag == "shotfile") then
        eq_from_shotfile = .true.
      else
        eq_from_shotfile = .false.
      end if
    else
      if(plasma_params%eq_exp == "Ext" .or. plasma_params%eq_exp == "EXT" .or.&
          plasma_params%eq_exp == "ext") then
        eq_from_shotfile = .false.
        if(plasma_params%eq_diag == "E2D") then
          plasma_params%Te_ne_mat = .true.
        end if
      else
        eq_from_shotfile = .true.
      end if
    end if
    if(stand_alone) then
      if(.not. plasma_params%Te_ne_mat) call read_input_lists(plasma_params)
      call setup_plasma_params(plasma_params, eq_from_shotfile)
    else
      call setup_plasma_params(plasma_params, eq_from_shotfile)
    end if
  end subroutine init_raytrace

  subroutine dealloc_raytrace(plasma_params)
  use f90_kind
  use mod_ecfm_refr_types,       only: plasma_params_type, stand_alone, output_level, &
                                       double_check_splines
#ifdef NAG
  USE nag_lib_support,           only : nag_deallocate
#endif
  use mod_ecfm_refr_interpol,       only: deallocate_rect_spline, deallocate_1d_spline
  implicit none
  type(plasma_params_type), intent(inout)                      :: plasma_params
  integer(ikind)                                             :: i
    if(stand_alone) then
      deallocate(plasma_params%n_e_prof, plasma_params%rhop_vec_ne)
      deallocate(plasma_params%T_e_prof, plasma_params%rhop_vec_Te)
      call deallocate_1d_spline(plasma_params%Te_spline)
      call deallocate_1d_spline(plasma_params%ne_spline)
#ifdef NAG
      if(output_level .and. double_check_splines) then
        call nag_deallocate(plasma_params%ne_spline_nag)
        call nag_deallocate(plasma_params%Te_spline_nag)
      end if
#endif
    else
      if(allocated(plasma_params%par)) deallocate(plasma_params%par)
      if(allocated(plasma_params%par_ne)) deallocate(plasma_params%par_ne)
      if(allocated(plasma_params%par_scal)) deallocate(plasma_params%par_scal)
    end if
    deallocate(plasma_params%Int_absz, plasma_params%Int_weights)
    deallocate(plasma_params%vessel_poly)
    call deallocate_rect_spline(plasma_params%rhop_spline)
    call deallocate_rect_spline(plasma_params%B_R_spline)
    call deallocate_rect_spline(plasma_params%B_t_spline)
    call deallocate_rect_spline(plasma_params%B_z_spline)
    if(plasma_params%Te_ne_mat) then
      call deallocate_rect_spline(plasma_params%T_e_spline_2D)
      call deallocate_rect_spline(plasma_params%n_e_spline_2d)
    end if
#ifdef NAG
    if(output_level .and. double_check_splines) then
      call nag_deallocate(plasma_params%rhop_spline_nag)
      call nag_deallocate(plasma_params%B_R_spline_nag)
      call nag_deallocate(plasma_params%B_t_spline_nag)
      call nag_deallocate(plasma_params%B_z_spline_nag)
    end if
    if(plasma_params%Te_ne_mat .and. output_level .and. double_check_splines) then
      call nag_deallocate(plasma_params%ne_spline_nag_2D)
      call nag_deallocate(plasma_params%Te_spline_nag_2D)
    end if
#endif
  end subroutine dealloc_raytrace

end module mod_ecfm_refr_raytrace_initialize
