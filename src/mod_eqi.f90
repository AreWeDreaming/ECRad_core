MODULE eqi_mod
use f90_kind
implicit none
#ifdef AUG
public  :: btf_corr_factor,                          &
           Bt_btf_corrected
#endif
CONTAINS

!*****************************************************************
#ifdef AUG
subroutine load_eqi_for_ECRad(time, time_beg, time_end, btf_mode, btf_corr_fact_ext, m_eq, n_eq, &
                              R, z, rhop, Br, Bt, Bz, R0_out, Btf0, &
                              R_mag, z_mag, pf_mag, R_sxp, z_sxp, pf_sxp)
    use aug_db_routines,            only: read_mbi, aug_kkEQpfx, aug_kkrzBrzt
    use eqi_global_params,          only: eqi
    implicit none
    real(rkind), intent(in)                             :: time
    real(rkind), intent(inout)                          :: time_beg, time_end ! Needs to inout, because of lack of in out statements in libdd
    character(8), intent(in)                            :: btf_mode
    real(rkind), intent(in)                             :: btf_corr_fact_ext
    integer(ikind), intent(in)                          :: m_eq, n_eq
    real(rkind), dimension(:), allocatable, intent(out) :: R, z
    real(rkind), dimension(:,:), allocatable, intent(out) ::  Br, Bt, Bz
    real(rkind), dimension(:,:), allocatable, intent(out) ::  rhop
    real(rkind), intent(out)                              :: R0_out, Btf0, R_mag, z_mag, &
                                                             pf_mag, R_sxp, z_sxp, pf_sxp
    real(rkind), dimension(:), allocatable      :: z_mat_helper, R_test, z_test,rhop_ne, &
                                                   rhop_Te, mbi_time, BtfABB
    real(4),     dimension(:), allocatable      :: R_temp, z_temp
    real(4), dimension(:,:), allocatable        :: pfm_temp
    real(rkind), dimension(1)                   :: Rv, zv, Br_vac, Bt_vac, Bz_vac
    real(rkind)                                 :: R0, Btf0_eq
    integer(ikind)                              :: i, error, debugunit, j, k, irhop_extra
    real(4)                                     :: temp_time
    logical                                     :: debug
    character(18)                               :: btf_fmt_str
    integer(4)                                  :: s_error, m_temp, n_temp, mdim
    R0 = 1.65d0
    R0_out = R0
    s_error = 0
    error     = 0
    debug     = .false.
    debugunit = 999
    temp_time = real(time, 4)
    m_temp = int(m_eq,4)
    n_temp = int(n_eq,4)
    mdim = m_temp + 1
    allocate(R_temp(mdim),z_temp(n_temp))
    allocate(pfm_temp(mdim,n_temp))
    call kkeqpfm(s_error, eqi%exper, eqi%diag,  & ! in
                    int(eqi%shot, 4), eqi%edit, & ! inout
                    temp_time, &
                    mdim, m_temp, &
                    n_temp, R_temp,  &
                    z_temp, pfm_temp)
    allocate(R(m_temp),z(n_temp))
    R = real(R_temp(1:m_temp),8)
    z = real(z_temp(1:n_temp),8)
    deallocate(R_temp, z_temp)
    if(s_error /= 0) then
      print*,"equilibrium EXP and diag ", eqi%exper, eqi%diag
      print*,"shot and time ", int(eqi%shot, 4), temp_time
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
                   debug, debugunit, eqi%shot,       & ! in
                   eqi%exper, eqi%diag, & ! in
                   eqi%edit,                        & ! inout
                   temp_time,                                  & ! inout
                   R_mag  = R_mag,                & ! out [m]  radial   coordinates of magnetic axis
                   z_mag  = z_mag,                & ! out [m]  vertical coordinates of magnetic axis
                   pf_mag = pf_mag,               & ! out [Vs] magnetic flux at magnetic axis
                   R_sxp  = R_sxp,                & ! out [m]  radial   coordinates of X-point
                   z_sxp  = z_sxp,                & ! out [m]  vertical coordinates of X-point
                   pf_sxp = pf_sxp)                 ! out [Vs] magnetic flux at X-point
    !print*, "aug_kkEQpfx finished"
    rhop = sqrt((rhop - pf_mag)/(pf_sxp - pf_mag))
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
                    debug, debugunit, eqi%shot,        & ! in
                    eqi%exper, eqi%diag,                         &  ! in
                    eqi%edit,                                       & ! inout
                    time,                                        & ! inout
                    r  = rv,                                      & ! in
                    z  = zv,                           & ! in
                    Br = Br_vac,                               & ! out
                    Bz = Bz_vac,                                & ! out
                    Bt = Bt_vac,                                & ! out
                    err_stop = .true.)
    Btf0_eq = Bt_vac(1) * rv(1) / R0
    !print*, Btf0_eq
    if (btf_mode == "BTFABB") then
      call read_mbi(eqi%shot, "AUGD", time_beg, time_end=time_end, &
                    time=mbi_time, BTFABB=BtfABB)         ! out
      if(size(BtfABB) < 1) stop "Failed to get Btf or BtfABB at R0 in mod_ecfm_refr_raytrace_initialize"
      Btf0 = sum(BtfABB) / size(BtfABB)
    else if (btf_mode == "old_Btf") then
      Btf0 = Btf0_eq
    else
      print*, btf_mode
      stop "sub Blos_along_beam_path: btf_mode not defined properly!"
    endif
    if (abs(Btf0_eq) > 20.d0 .or. abs(Btf0_eq) < 1.d-3) then
      write(6,'(a,f8.4,a,e12.4)')"Btf0_eq(time =", time, " s) = ", Btf0_eq
      stop 'subroutine Blos_along_beam_path: |Btf0_eq| > 20 T or |Btf0_eq| < 1 mT ???'
    endif
    if (abs(Btf0) > 20.d0 .or. abs(Btf0) < 1.d-3) then
      write(6,'(a,f8.4,a,e12.4)')"Btf0(time =", time, " s) = ", Btf0
      stop 'subroutine Blos_along_beam_path: |Btf0| > 20 T or |Btf0| < 1 mT ???'
    endif
    do j = 1, n_temp
      z_mat_helper(:) = z(j)
      call aug_kkrzBrzt(error,                                   & ! out
                    debug, debugunit, eqi%shot,        & ! in
                    eqi%exper, eqi%diag,                         &  ! in
                    eqi%edit,                                       & ! inout
                    time,                          & ! in
                    r  = R,                                      & ! in
                    z  = z_mat_helper,                           & ! in
                    Br = Br(:,j),                               & ! out
                    Bz = Bz(:,j),                                & ! out
                    Bt = Bt(:,j),                                & ! out
                    err_stop = .true.)
      call Bt_btf_corrected(eqi%shot, btf_mode, btf_corr_fact_ext, &
                                               Btf0_eq, Btf0, R0, R, Bt(:,j))
      if (any(abs(Bt(:,j)) > 20.d0)) then
        do i = 1, size(Bt)
          write(6,'(i3,e12.4)')i, Bt(i,j)
        enddo
        stop 'subroutine setup_eq_from_shotfile: |Bt| > 20 T ???'
      endif
    end do
    deallocate(z_mat_helper)
end subroutine load_eqi_for_ECRad

function btf_corr_factor(shot, btf_mode, btf_corr_fact_ext)
! ripple correction of toroidal magnetic field (-> fppadjb2rz.c from Wolfgang Suttrop)
! -> /afs/ipp/u/eced/CEC/libece/btot_los.c
use f90_kind
implicit none
integer(ikind), intent(in)           :: shot               ! shot number
character(*),   intent(in)           :: btf_mode           ! method to calculate Btf: "old_Btf", "BTFABB"
real(rkind),    intent(in), optional :: btf_corr_fact_ext  ! user defined correction factor
real(rkind) :: btf_corr_factor

if (trim(btf_mode) == "old_Btf") then
  if ( shot <= 8645 ) then      ! no correction for divertor I
    btf_corr_factor = 1.0d0
  else if ( shot < 9219 ) then  ! correction for divertor II
    btf_corr_factor = 0.98d0    !  MAI BTF factor 0.4260 but MAI BTF factor 0.4495 used in FPP
  else if ( shot < 10392 ) then
    btf_corr_factor = 0.98d0    !   new MAI BTF factor 0.4495
  else if ( shot < 28524 ) then
    btf_corr_factor = 1.00      ! new MAI BTF factor 0.44051
    btf_corr_factor = btf_corr_factor/0.99d0  ! normalization factor 0.99 of ripple correction at R=2.2m
  else
    ! this factor turned out to be wrong with the measurements of the currents in the toroidal field coils
    !btf_corr_factor = 1.015     ! Changed Wed. 16 Jan.2013 (1.5 % error). Bt calibration pulse made 30th Nov.2012
    !btf_corr_factor = 1.0175     ! Changed Wed. 16 Jan.2013 (1.5 % error). Bt calibration pulse made 30th Nov.2012
    !                            ! Change calibration factor for all the pulses of this campaign from Sept 2012
    btf_corr_factor = 1.d0/0.99d0      ! Changed Wed. 15 Oct.2013 (1.0 % error).
  endif

else if (trim(btf_mode) == "BTFABB") then
  btf_corr_factor = 1.d0        ! no correction on Btf necessary

else
  stop "fun btf_corr_factor: btf_mode not properly defined!"
endif

! user defined correction factor
!-------------------------------
if (present(btf_corr_fact_ext)) btf_corr_factor = btf_corr_factor * btf_corr_fact_ext    

!stop "fun btf_corr_factor"
end function btf_corr_factor

subroutine Bt_btf_corrected(shot, btf_mode, btf_corr_fact_ext, &
                                               btf0_eq, btf0, R0, R, Bt)
use f90_kind
implicit none
integer(ikind),            intent(in)    :: shot    ! shot number
character(*),              intent(in)    :: btf_mode  ! method to calculate Btf: "old_Btf", "BTFABB"
real(rkind),               intent(in)    :: btf_corr_fact_ext  ! [] user defined correction factor
real(rkind),               intent(in)    :: Btf0_eq ! [T] vacuum toroidal magnetic field at R0 from equilibrium
real(rkind),               intent(in)    :: Btf0    ! [T] vacuum toroidal magnetic field at R0 to be used
real(rkind),               intent(in)    :: R0      ! [T] vacuum toroidal magnetic field at R0
real(rkind), dimension(:), intent(in)    :: R       ! [m] radial   coordinate  along beam path
real(rkind), dimension(:), intent(inout) :: Bt      ! [T] toroidal magnetic field along beam path
real(rkind)    :: Btok, Btok_eq, Bdia, Badj
real(rkind)    :: xh,xh_no_ripple
integer(ikind) :: i, N

!open(66,file="ece_bdia.dat",position='append')
N = size(R)
do i = 1, N
  ! correction factor of toroidal magnetic field
  !---------------------------------------------
  xh = btf_corr_factor(shot              = shot,            & ! in
                       btf_mode          = btf_mode,        & ! in
                       btf_corr_fact_ext = btf_corr_fact_ext) ! in
  !print*, "BT CORRECTION FACTOR", xh
  ! adjust only 1/R part of toroidal field with btf_corr_factor and toroidal ripple
  !Btok = Btf0 * R0 / r(i)    ! vacuum toroidal field
  !Bdia = Bt(i) - Btok
  !Badj = (Btok * xh) + Bdia

  Btok_eq = Btf0_eq * R0 / R(i)    ! vacuum toroidal field from equilibrium
  if (btf_mode == "BTFABB") then
    Btok  = Btf0    * R0 / R(i)    ! vacuum toroidal field to be used
  else if (btf_mode == "old_Btf") then
    Btok  = Btok_eq                ! vacuum toroidal field from equilibrium
  else
    stop "sub Bt_btf_and_field_ripple_corrected: btf_mode not defined properly!"
  endif
  Bdia = Bt(i) - Btok_eq           ! subtract vacuum toroidal field from equilibrium to obtain diamagnetic field
  Bt(i) = (Btok * xh) + Bdia ! add corrected vacuum toroidal field to be used
  !write(66,'(4e14.6)')r(i), Btok, Bt(i), Bdia
  !write(66,'(3e14.6)')r(i), Br(i), Bz(i)
enddo
!close(66)

!stop "sub Btot_btf_and_field_ripple_corrected"
end subroutine Bt_btf_corrected

!****************************************************************
    

#endif
END MODULE eqi_mod
