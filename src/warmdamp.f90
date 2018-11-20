module farina_abs
use f90_kind
implicit none
public :: warmdamp
contains

  subroutine warmdamp(op, oc, Nr, theta, te, m, imod, imNprwsq, pol_coeff, pol_vect)

    ! ---------------------------------------------------------------------
    ! Wrapper for D. Farina' warmdisp routine.
    ! ---------------------------------------------------------------------

    use ecdisp , only : set_extv, warmdisp, larmornumber
    use constants, only: mass_e, c0, e0
    ! ... implicit none statement ...
    implicit none

    ! ... paramenters ...
    integer, parameter :: imax = 3 ! max harmonics number
    integer, parameter :: fast = 3
    integer, parameter :: r8 = selected_real_kind(15,300)

    ! ... input variables ...
    integer, intent(in)  ::   &
         imod, &                   ! = -1 for O-mode, = 1 for X-mode
         m
    real(r8), intent(in) ::   &
         op,                  & ! = omega_p^2 / omega^2 
         oc,                  & ! = omega_c^2 / omega^2
         theta,               & ! angle between N and the magnetic field
         te                     ! electron temperature in keV
    real(r8), intent(inout) ::  Nr! modulus of the refractive index N
    
    ! ... output variables ...
    real(r8), intent(out) ::  &
         imNprwsq                 ! Im(N_perp^2)
    real(r8), intent(out),optional :: pol_coeff
    complex(r8), intent(out), dimension(:), optional :: pol_vect
    ! ... local variables ...
    integer  ::               &
         err,                 & ! error flag
         nharm,               & ! hamrnic number
         lrm                    ! effective harmonic number
    real(r8) ::               &
         sqoc,                & ! = omega_c / omega
         mu,                  & ! = 511.e0_r8/te
         Nll, Npr               ! parallel and perpendicular N
    complex(r8) ::            &
         Nprw,                & ! Npr obtained from the disp. rel.
         ex, ey, ez             ! polarization unit vector
         
! ... f2py directives ...
!f2py integer intent(in) :: imod    
!f2py real*8 intent(in) :: op, oc, Nr, theta, te
!f2py real*8 intent(out) :: imNprw    

    ! =====================================================================
    ! Executable statements 

    ! ... initialize common variables of the module ecdisp ...
    call set_extv

    ! ... definition of parameters ...
    sqoc = sqrt(oc)
    Nll  = Nr * cos(theta)
    Npr  = sqrt(Nr**2 - Nll**2)
    mu   = mass_e * c0**2/(te * e0)

    ! ... find the harmonic number ...
    nharm = larmornumber(sqoc, Nll, mu)
    !if(nharm == 1) then
    !   imNprw = 0.
    !   return
    !end if

    lrm = min(imax, nharm)
    
    ! ... estimate the warm-plasma dispersion function ... 
    call warmdisp(op, sqoc, Nll, mu, imod, fast, lrm, Npr,                &
         Nprw, ex, ey, ez, err)
    
    ! ... extract imaginary part of the refractive index ...
    imNprwsq = aimag(Nprw**2)
    Nr =  real(Nprw) / sin(theta) ! we only have N_perp here !
    if(present( pol_coeff))  pol_coeff = (abs(ex)**2 + abs(ey)**2)/(abs(ex)**2 + abs(ey)**2 + abs(ez)**2)
    if(present(pol_vect)) then
      pol_vect(1) = ex
      pol_vect(2) = ey
      pol_vect(3) = ez
    end if
    ! ... do not allow negative values, put 0 instead ...
    if (imNprwsq < 0.) then
       imNprwsq = 0.
       Nr = 0.d0
    end if

    return
  end subroutine warmdamp

  end module farina_abs
