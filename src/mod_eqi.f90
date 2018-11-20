MODULE eqi_mod
use f90_kind
implicit none

public  :: btf_corr_factor,                          &
           Bt_btf_corrected

CONTAINS

!*****************************************************************


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
    


END MODULE eqi_mod
