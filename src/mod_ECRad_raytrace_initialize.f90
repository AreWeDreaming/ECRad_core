module mod_ECRad_raytrace_initialize
    use f90_kind
    use mod_ECRad_types,       only: spl_type_2d
    implicit none
    ! Two interfaces to allow interpolation of Te and ne using both single rhop values and vectors of rhop values
    type(spl_type_2d)                       :: Psi_spline

    public :: init_raytrace, dealloc_raytrace

    private :: get_psi_ax, &
               eval_Psi, &
               read_topfile, &
               read_Te_ne_matrix, &
               make_topfile, &
               setup_plasma_params, &
               read_3D_vessel_file, &
               read_input_lists



    contains

  subroutine get_psi_ax(R, z, i_ax, j_ax, R_ax, z_ax)
  ! Interpolates
    use f90_kind
    use mod_ECRad_types,       only: output_level
    use nr_mod,                    only: powell ! Minimizer -> Psi_ax has to be a minimum
    use mod_ECRad_interpol,    only: rect_spline
    implicit none
    real(rkind), dimension(:), intent(in)  :: R, z
    integer(ikind)           , intent(in)  :: i_ax, j_ax
    real(rkind),               intent(out) :: R_ax
    real(rkind),               intent(out) :: z_ax
    real(rkind), dimension(2)              :: p
    real(rkind), dimension(2,2)            :: xi
    real(rkind)                            :: ftol, psi_out, dx, dy
    integer(ikind)                         :: iter
    p(1) = R(i_ax)
    p(2) = z(j_ax)
    call rect_spline(Psi_spline, p(1), p(2), psi_out, dfdx=dx, dfdy=dy)
    xi(:,1) = dx
    xi(:,2) = dy
    ftol = 1.d-6
    call powell(p,xi,ftol,iter, psi_out, eval_Psi)
    R_ax = p(1)
    z_ax = p(2)
    !if(output_level) print*, "Magnetic axis position and rhop on axis",  p(1), p(2), sqrt((eval_Psi(p) - psi_ax)/(psi_sep - psi_ax))
    if(output_level) print*, "Magnetic axis position and rhop on axis",  p(1), p(2), sqrt(eval_Psi(p))
  end subroutine get_psi_ax

  function eval_Psi(R)
    use f90_kind
    use mod_ECRad_interpol,    only: rect_spline
    implicit none
    real(rkind), dimension(:), intent(in)   :: R
    real(rkind)                             :: eval_Psi
    call rect_spline(Psi_spline, R(1), R(2), eval_Psi)
  end function eval_Psi

  subroutine read_topfile(plasma_params, R, z, rhop, Br, Bt, Bz, itime)
      use f90_kind
      use mod_ECRad_types,       only: plasma_params_type, data_folder, output_level
      use mod_ECRad_interpol,    only: make_rect_spline, deallocate_rect_spline
      implicit none
      type(plasma_params_type) , intent(inout) :: plasma_params
      real(rkind), dimension(:), allocatable, intent(out) :: R, z
      real(rkind), dimension(:,:), allocatable, intent(out) ::  Br, Bt, Bz
      real(rkind), dimension(:,:), allocatable, intent(out) ::  rhop
      integer(ikind), intent(in), optional                  :: itime
      character(250) :: filename
      integer(ikind) :: i, j, i_ax, j_ax
      real(rkind) :: radin, radout, pf_sxp, pf_mag
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
      if(output_level) print*, "Dimensions of external data: ", plasma_params%m, plasma_params%n
      allocate(R(plasma_params%m), z(plasma_params%n), Br(plasma_params%m,plasma_params%n), &
               Bt(plasma_params%m,plasma_params%n), Bz(plasma_params%m,plasma_params%n), &
               rhop(plasma_params%m, plasma_params%n))
      read(70,'(A)') line
      read(70,*) radin,radout, pf_sxp
      if(abs(pf_sxp -1.e0) > 1.e-6) then
        print*, "The Psi provided to ECRad MUST be normalized such that rhop = Sqrt(PSI)"
        print*, "Hence, psi at the separatrix must be 1"
        print*, "Psi at the separatrix provided is", pf_sxp
        call abort()
      end if
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
!      if(pf_sxp < rhop(plasma_params%m/2, plasma_params%n/2)) then
!         sgn = -1.e0_rkind
!         pf_sxp = pf_sxp * sgn
!      else
!         sgn = 1.e0_rkind
!      endif
      pf_mag = 1.d99
      i_ax = -1
      do j = 1, plasma_params%n
        if(pf_mag > minval(rhop(:,j))) then
          pf_mag = minval(rhop(:,j))
          i_ax = minloc(rhop(:,j), dim = 1)
          j_ax = j
        end if
      end do
      if(i_ax == -1) then
        print*, "Could not find magnetic axis in read_topfile in mod_raytrace_initialize.f90"
        print*, "Check if psi matrix dimensions correct: R, z", plasma_params%m, plasma_params%n
        stop "Input error - check topfile"
      end if
      if(output_level) print*, "Retrieving magn axis"
     ! ECRad wants normalized psi for which rhop = sqrt(Psi)!!!
      call make_rect_spline(Psi_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), R, z, rhop)
      call get_psi_ax(R, z, i_ax, j_ax, plasma_params%R_ax, plasma_params%z_ax)
      rhop = sqrt(rhop) ! - pf_mag)/(pf_sxp - pf_mag)
      if(any(rhop /= rhop) .or. minval(rhop) < 0.d0) then
        print*, "critical error when reading topfile, NAN or negative value in normalized Psi found"
        print*, rhop
        stop "sub read topfile in mod_raytrace_intialize"
      end if
      call deallocate_rect_spline(Psi_spline)
  end subroutine read_topfile

  subroutine read_Te_ne_matrix(plasma_params, R, z, Te, ne)
      use f90_kind
      use mod_ECRad_types,       only: plasma_params_type, data_folder, min_density, max_Te, min_Te
      implicit none
      type(plasma_params_type) , intent(inout) :: plasma_params
      real(rkind), dimension(:), allocatable, intent(in) :: R, z
      real(rkind), dimension(:,:), allocatable, intent(out) ::  Te, ne
      real(rkind), dimension(size(R))                       :: R_check
      real(rkind), dimension(size(z))                       :: z_check
      character(250) :: filename
      integer(ikind) :: m_check, n_check, i,j
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
  
  subroutine read_3D_vessel_file(plasma_params)
  use f90_kind
    use mod_ECRad_types,       only: Vessel_bd_filename, plasma_params_type
    implicit none
    type(plasma_params_type), intent(inout) :: plasma_params
    integer(ikind) :: i, j
    character(200) :: dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8
    open(87, file=Vessel_bd_filename)
    read(87,*) dummy1
    read(87,*) plasma_params%Use_3D_vessel%n_phi, plasma_params%Use_3D_vessel%n_contour, &
               dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, dummy8
    allocate(plasma_params%Use_3D_vessel%vessel_data_R(plasma_params%Use_3D_vessel%n_contour, &
                                                     plasma_params%Use_3D_vessel%n_phi), &
             plasma_params%Use_3D_vessel%vessel_data_z(plasma_params%Use_3D_vessel%n_contour, &
                                                     plasma_params%Use_3D_vessel%n_phi), &
             plasma_params%Use_3D_vessel%phi(plasma_params%Use_3D_vessel%n_phi))
    do j = 1, plasma_params%Use_3D_vessel%n_phi
      read(87,*) plasma_params%Use_3D_vessel%phi(j)
      do i = 1, plasma_params%Use_3D_vessel%n_contour
        read(87,*) plasma_params%Use_3D_vessel%vessel_data_R(i,j), &
                   plasma_params%Use_3D_vessel%vessel_data_z(i,j), &
                   dummy1
      end do
    end do
    close(87)
    if(plasma_params%Use_3D_vessel%n_phi == 1) then
      plasma_params%Use_3D_vessel%phi_max = 361.d0 !
    else
      plasma_params%Use_3D_vessel%phi_max = maxval(plasma_params%Use_3D_vessel%phi) ! Necessary to use symmetry
    end if
  end subroutine read_3D_vessel_file

  subroutine make_topfile(R, z, rhop, Br, Bt, Bz, itime)
    use f90_kind
    use mod_ECRad_types,       only: data_folder
    implicit none
    real(rkind), dimension(:), allocatable, intent(in) :: R, z
    real(rkind), dimension(:,:), allocatable, intent(in) ::  Br, Bt, Bz
    real(rkind), dimension(:,:), allocatable, intent(in) ::  rhop
    integer(ikind), intent(in), optional                  ::  itime
    character(250) :: filename
    integer(ikind) :: i,j
    if(present(itime)) then
      write(filename, fmt = "(A7I5.5)") "topfile",  itime
      filename = trim(data_folder) // trim(filename)
    else
      filename =  trim(data_folder) // "topfile"
    end if
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

  subroutine setup_plasma_params(plasma_params, R_in, z_in, rhop_in, B_r_in, &
                                 B_t_in, B_z_in, R_ax, z_ax, T_e_mat_in, n_e_mat_in)
    use f90_kind
#ifdef NAG
    use mod_ECRad_types,        only: double_check_splines, plasma_params_type, h_x_glob, h_check, &
                                          stand_alone, vessel_bd_filename, use_3D, working_dir
#else
    use mod_ECRad_types,        only: plasma_params_type, h_x_glob, h_check, &
                                          stand_alone, vessel_bd_filename, use_3D, working_dir
#endif
    use mod_ECRad_interpol,     only: make_rect_spline, make_1d_spline, rect_spline
    use constants,                  only: pi
    use quadrature,                 only: cdgqf
#ifdef NAG
    USE nag_spline_2d,              only: nag_spline_2d_interp, &
    									                    nag_error, nag_set_error, nag_spline_2d_eval
    USE nag_spline_1d,              only: nag_spline_1d_interp
    USE nag_error_handling
#endif
    use ripple3d,                   only: init_ripple!, validate_ripple_grad
#ifdef USE_3D
    use magconfig,                only: MConf_Setup_Config, MConf_Load_MagConfig, mcClone
#endif
#ifdef OMP
    use omp_lib,                    only: omp_get_max_threads
#endif
    implicit none
    type(plasma_params_type), intent(inout)     :: plasma_params
	  real(rkind), dimension(:), intent(in), optional :: R_in, z_in
    real(rkind), dimension(:,:), intent(in), optional :: rhop_in, B_r_in, B_t_in, B_z_in, T_e_mat_in, n_e_mat_in
    real(rkind), intent(in), optional                 :: R_ax, z_ax
    real(rkind), dimension(:,:), allocatable    :: B_r, B_t, B_z, T_e_mat, n_e_mat
    integer(ikind), dimension(:), allocatable   :: R_index_lower, z_index_lower, R_index_upper, z_index_upper
    real(rkind)                                 :: B_last, B_vac_R0
    integer(ikind)                              :: i
    integer(ikind)                              :: loaded
    character(1)  :: sep
    character(200)                              :: filename
    integer(ikind)                              :: num_threads
    plasma_params%int_step_cnt = plasma_params%rad_trans_sections * plasma_params%rad_transp_solver_order
    if(plasma_params%int_step_cnt == 0) then
        print*, "For some reason plasma_params%int_step_cnt is zero!"
        print*, "This makes no sense since it would mean that you split the radiation transport into N * 0 steps"
        print*, "Please check setup_plasma_params in mod_ECRad_raytrace_initialize"
        print*, plasma_params%rad_trans_sections, plasma_params%rad_transp_solver_order
        call abort()
    end if
    allocate( plasma_params%Int_weights(plasma_params%int_step_cnt),plasma_params%Int_absz(plasma_params%int_step_cnt))
    if(plasma_params%on_the_fly_raytracing) then
      call cdgqf( int(plasma_params%int_step_cnt,kind=4), int(1,kind=4), 0.d0, 0.d0, plasma_params%Int_weights, plasma_params%Int_absz)
    else
      do i = 1, plasma_params%int_step_cnt
        plasma_params%Int_absz(i) = real(i ,8) / real(plasma_params%int_step_cnt, 8) ! spans the points in a way so that there is no overlap
                                                                                     ! Int_absz(1) > 0.d0
      end do
    end if
    if(.not. plasma_params%Te_ne_mat .and. stand_alone) then
#ifdef NAG
      if(double_check_splines .and. output_level) then
        call nag_spline_1d_interp(plasma_params%rhop_vec_ne, plasma_params%n_e_prof, plasma_params%ne_spline_nag)
        call nag_spline_1d_interp(plasma_params%rhop_vec_Te, plasma_params%T_e_prof, plasma_params%Te_spline_nag)
      end if
#endif
      if(plasma_params%prof_log_flag) then
        call make_1d_spline( plasma_params%ne_spline, int(size(plasma_params%rhop_vec_ne),4), plasma_params%rhop_vec_ne, log(plasma_params%n_e_prof * 1.d-19))
        call make_1d_spline( plasma_params%Te_spline, int(size(plasma_params%rhop_vec_Te),4), plasma_params%rhop_vec_Te, log(plasma_params%T_e_prof))
      else
        call make_1d_spline( plasma_params%ne_spline, int(size(plasma_params%rhop_vec_ne),4), plasma_params%rhop_vec_ne, plasma_params%n_e_prof)
        call make_1d_spline( plasma_params%Te_spline, int(size(plasma_params%rhop_vec_Te),4), plasma_params%rhop_vec_Te, plasma_params%T_e_prof)
      end if
    end if
#ifdef USE_3D
    if(use_3D) then
        plasma_params%Scenario%work_dir = working_dir
        call read_3D_vessel_file(plasma_params)
        call MConf_Setup_Config(plasma_params%Scenario)
        call MConf_Load_MagConfig(loaded, plasma_params%Scenario%mConfAddr)
        if(loaded /= 1) then
            print*, "Failed to load equilibrium"
            print*, "Please check that", plasma_params%Scenario%name_config, " exists!"
            print*, "MConf error code", loaded
            call abort()
        end if
#ifdef OMP
      !$omp parallel
      num_threads = omp_get_max_threads()
      !$omp end parallel
#else
      num_threads = 1
#endif
      allocate(plasma_params%mconf_addresses(num_threads))
      plasma_params%mconf_addresses(1) = plasma_params%Scenario%mConfAddr
#ifdef OMP
      do i = 2, num_threads
        plasma_params%mconf_addresses(i) = mcClone(plasma_params%mconf_addresses(1))
      end do
#endif
      return
    end if
#endif
    ! Equilibrium/Magnetic field splines are set up in the following
    if(.not. present(R_in)) then
      call read_topfile(plasma_params, plasma_params%R, plasma_params%z, plasma_params%rhop, B_r, B_t, B_z)
      if(plasma_params%Te_ne_mat) then
        call read_Te_ne_matrix(plasma_params, plasma_params%R, plasma_params%z, T_e_mat, n_e_mat)
      end if
    else
      plasma_params%m = size(R_in)
      plasma_params%n = size(z_in)
      allocate(plasma_params%R(plasma_params%m), &
               plasma_params%z(plasma_params%n), &
               plasma_params%rhop(plasma_params%m,plasma_params%n), &
               B_r(plasma_params%m,plasma_params%n), &
               B_t(plasma_params%m,plasma_params%n), &
               B_z(plasma_params%m,plasma_params%n))
      if(present(T_e_mat_in)) then
        allocate(T_e_mat(plasma_params%m,plasma_params%n), &
                 n_e_mat(plasma_params%m,plasma_params%n))
      end if
      plasma_params%R=R_in
      plasma_params%z=z_in
      do i =1, plasma_params%m
        plasma_params%rhop(i,:)=rhop_in(i,:)
        B_r(i,:) = B_r_in(i,:)
        B_t(i,:) = B_t_in(i,:)
        B_z(i,:) = B_z_in(i,:)
        if(present(T_e_mat_in)) then
          T_e_mat(i,:) = T_e_mat_in(i,:)
          n_e_mat(i,:) = n_e_mat_in(i,:)
        end if
      end do
      plasma_params%R_ax = R_ax
      plasma_params%z_ax = z_ax
    end if
    ! No information on where the profiles start to be good go with the default value
    if(plasma_params%Te_ne_mat) plasma_params%rhop_max = plasma_params%rhop_entry !
    h_check = h_x_glob * (1.d0 + 1.e-5)
    plasma_params%R_min = plasma_params%R(1)
    plasma_params%R_max = plasma_params%R(plasma_params%m)
    plasma_params%z_min = plasma_params%z(1)
    plasma_params%z_max = plasma_params%z(plasma_params%n)
    plasma_params%R_step = (plasma_params%R_max - plasma_params%R_min) / real(plasma_params%m - 1,8)
    plasma_params%z_step = (plasma_params%z_max - plasma_params%z_min) / real(plasma_params%n - 1,8)
    call make_rect_spline(plasma_params%rhop_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, plasma_params%rhop)
    call make_rect_spline(plasma_params%B_r_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_r)
    call make_rect_spline(plasma_params%B_t_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_t)
    call make_rect_spline(plasma_params%B_z_spline, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, B_z)
    if(plasma_params%Te_ne_mat) then
      call make_rect_spline(plasma_params%T_e_spline_2D, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, log(T_e_mat))
      ! Warning: ne is scaled down to avoid numerical problems!
      call make_rect_spline(plasma_params%n_e_spline_2D, int(plasma_params%m, 4), int(plasma_params%n, 4), plasma_params%R, plasma_params%z, log(n_e_mat * 1.e-19))
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
    call rect_spline(plasma_params%B_r_spline, plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    call rect_spline(plasma_params%B_t_spline,plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    call rect_spline(plasma_params%B_z_spline,plasma_params%R_ax,plasma_params%z_ax,B_last)
    plasma_params%B_ax = plasma_params%B_ax + B_last**2
    plasma_params%B_ax = Sqrt(plasma_params%B_ax)
    B_vac_R0 = B_t(plasma_params%m, int(plasma_params%n/2)) ! Outer most point of Bt -> pure vacuum
    B_vac_R0 = B_vac_R0 * plasma_params%R(plasma_params%m) / plasma_params%R_ax
    if(plasma_params%w_ripple) call init_ripple(plasma_params%R_ax, B_vac_R0)
    allocate(R_index_lower(plasma_params%n), z_index_lower(plasma_params%m), R_index_upper(plasma_params%n), z_index_upper(plasma_params%m))
    R_index_lower(:) = 1
    z_index_lower(:) = 1
    R_index_upper(:) = plasma_params%m
    z_index_upper(:) = plasma_params%n
    plasma_params%int_step_cnt = plasma_params%rad_trans_sections * plasma_params%rad_transp_solver_order
    deallocate(B_r, B_z, B_t, R_index_lower,z_index_lower,R_index_upper,z_index_upper)!, B_t
    if(plasma_params%Te_ne_mat) deallocate(T_e_mat, n_e_mat)
    if(stand_alone) then
      ! Load the polygon describing the vessel wall
      open(66, file = vessel_bd_filename)
      read(66, "(I7.7)") plasma_params%m_vessel_bd
      allocate(plasma_params%vessel_poly%x(plasma_params%m_vessel_bd), plasma_params%vessel_poly%y(plasma_params%m_vessel_bd))
      do i = 1, plasma_params%m_vessel_bd
        read(66,"(E19.12E2A1E19.12E2)") plasma_params%vessel_poly%x(i), sep,  plasma_params%vessel_poly%y(i)
      end do
      close(66)
     end if
  end subroutine setup_plasma_params

  subroutine read_input_lists(plasma_params)
  ! Read the input data n_e, T_e and wall file
    use mod_ECRad_types,       only: plasma_params_type, n_e_filename, T_e_filename, &
                                     min_density, min_Te, max_Te, output_level
    use f90_kind
    use constants, only               : pi
    implicit none
    type(plasma_params_type), intent(inout)           :: plasma_params
    integer(ikind)                                    :: i
    character(1)                                      :: sep
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
      if(output_level) print*, "Largest rhopol in profiles", plasma_params%rhop_max
      if(plasma_params%rhop_exit > plasma_params%rhop_max-plasma_params%delta_rhop_exit) then
        plasma_params%rhop_exit = plasma_params%rhop_max-plasma_params%delta_rhop_exit
      end if
      close(66)
      if(plasma_params%m_T_e_prof < 40 .or. plasma_params%m_n_e_prof < 40) then
        print*, "points for density and Te profile:", plasma_params%m_T_e_prof
        stop "both density and temperature profile have to have MORE than 40 points"
      end if
    end if
  end subroutine read_input_lists

  subroutine init_raytrace(plasma_params, R, z, rhop, Br, Bt, Bz, R_ax, z_ax, T_e_mat, n_e_mat)
    use f90_kind
    use mod_ECRad_types,       only: plasma_params_type, stand_alone, use_3D, eq_mode
    Use constants                , only : pi
    implicit none
    type(plasma_params_type), intent(inout)           :: plasma_params
    real(rkind), dimension(:), intent(in), optional   :: R, z
    real(rkind), dimension(:,:), intent(in), optional :: rhop, Br, Bt, Bz, T_e_mat, n_e_mat
    real(rkind), intent(in), optional                 :: R_ax, z_ax
    if(eq_mode == "2D") then
      plasma_params%Te_ne_mat = .true.
    end if
    if(stand_alone) then
      if(.not. plasma_params%Te_ne_mat) call read_input_lists(plasma_params)
      call setup_plasma_params(plasma_params)
    else
      if(present(T_e_mat)) then
        call setup_plasma_params(plasma_params, R, z, rhop, Br, Bt, Bz, R_ax, z_ax, T_e_mat, n_e_mat)
      else
        call setup_plasma_params(plasma_params, R, z, rhop, Br, Bt, Bz, R_ax, z_ax)
      end if
    end if
  end subroutine init_raytrace

  subroutine dealloc_raytrace(plasma_params)
  use f90_kind
#ifdef NAG
  use mod_ECRad_types,       only: plasma_params_type, stand_alone, double_check_splines, &
                                   Use_3D, ray_init
  USE nag_lib_support,           only : nag_deallocate
#else
  use mod_ECRad_types,       only: plasma_params_type, stand_alone, Use_3D, ray_init
#endif
  use mod_ECRad_interpol,       only: deallocate_rect_spline, deallocate_1d_spline
  implicit none
  type(plasma_params_type), intent(inout)                    :: plasma_params
  ! Need to do this first since we can have this allocated but not anything else
  if(.not. use_3D) then
    if(allocated(plasma_params%vessel_poly%x)) deallocate(plasma_params%vessel_poly%x, plasma_params%vessel_poly%y)
  end if
    if(.not. (allocated(plasma_params%R) .or. allocated(plasma_params%Use_3D_vessel%vessel_data_R))) return
    if(stand_alone) then
      deallocate(plasma_params%n_e_prof, plasma_params%rhop_vec_ne)
      deallocate(plasma_params%T_e_prof, plasma_params%rhop_vec_Te)
#ifdef NAG
      if(output_level .and. double_check_splines) then
        call nag_deallocate(plasma_params%ne_spline_nag)
        call nag_deallocate(plasma_params%Te_spline_nag)
      end if
#endif
    else
      call deallocate_1d_spline(plasma_params%Te_spline)
      call deallocate_1d_spline(plasma_params%ne_spline)
    end if
    deallocate(plasma_params%Int_absz, plasma_params%Int_weights)
    deallocate(plasma_params%R, plasma_params%z, plasma_params%rhop)
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
  if(Use_3D) deallocate(plasma_params%Use_3D_vessel%vessel_data_R, &
                        plasma_params%Use_3D_vessel%vessel_data_z, &
                        plasma_params%Use_3D_vessel%phi)
  ray_init = .false.
  end subroutine dealloc_raytrace

end module mod_ECRad_raytrace_initialize
