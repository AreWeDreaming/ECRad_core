
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
