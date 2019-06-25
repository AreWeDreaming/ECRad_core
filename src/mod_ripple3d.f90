module ripple3d
use f90_kind

type ripple_params_type
! Parametrisation of the field ripple
! Parametrisation and parameters for analytic rippple correction provided by R. Dux
  real(rkind)               :: A0, A1, A2
  real(rkind)               :: B0, B1, B2
  real(rkind)               :: C0, C1, C2
  real(rkind)               :: D0, D1, D2
  real(rkind)               :: E0, E1, K0, K1
  real(rkind)               :: Btf0, R0, R_min, R_max
  integer(ikind)            :: periods
end type ripple_params_type
! Uncomment this type for use separate from the electron cyclotron forward model
type grad_type
  real(rkind)               :: dR, dphi, dz
end type grad_type
type(ripple_params_type)           :: ripple_params
public :: init_ripple, get_ripple, get_ripple_w_grad
private :: ripple_params

contains
  subroutine init_ripple(R0, Btf0)
    implicit none
    real(rkind), intent(in)           :: R0, Btf0
    ! Btf0 vacuum toroidal magnetic field at R0
    if(abs(Btf0) > 20.d0 .or. abs(Btf0) < 20.d-3 .or. R0 < 0.d0) then
      print*, "Nonsensical Btf0 or R0 in init_ripple in mod_ripple3d.f90"
      print*, "R0 [m], Btf0 [T]", R0, Btf0
      stop "Input error in init_ripple in mod_ripple3d.f90"
    end if
    ripple_params%R_max = 2.25d0 ! If R >  R_max the ripple will be ignored
    ! Routine not tested for HFS -> TODO
    ripple_params%R_min = 1.7d0 ! If R <  R_min_boundary the ripple will be ignored
    ripple_params%periods = 16 ! Symmetry of the machine
    ripple_params%R0 = R0; ripple_params%Btf0 = Btf0
    ripple_params%A0 = -21.959; ripple_params%A1 = 4.653; ripple_params%A2 = 1.747
    ripple_params%B0 = 7.891; ripple_params%B1 = -1.070; ripple_params%B2 = -0.860
    ripple_params%C0 = -23.315; ripple_params%C1 = 5.024; ripple_params%C2 = 0.699
    ripple_params%D0 = 8.196; ripple_params%D1 = -1.362; ripple_params%D2 = -0.367
    ripple_params%E0 = -19.539; ripple_params%E1 = 3.905;
    ripple_params%K0 = 6.199; ripple_params%K1 = -0.767
    ! Uncomment if you wish to test the gradients
!    call validate_ripple_grad()
!    stop "validated ripple"
  end subroutine init_ripple

  subroutine get_ripple(R_vec, B_ripple)
    USE f90_kind
    use constants, only : pi
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: R_vec ! R, psi, z; psi in range of [0, 2 N pi[
    real(rkind), dimension(:)  , intent(out)  :: B_ripple
    real(rkind) :: Btf0, psi
    if(R_vec(1) < ripple_params%R_min .or. R_vec(1) > ripple_params%R_max) then
      B_ripple(:) = 0.d0
      return
    end if
    Btf0 = ripple_params%Btf0 / R_vec(1) * ripple_params%R0
    ! 16 coils at ASDEX Upgrade, hence we have also 16 periods in the ripple function when going once around the torus
    psi = R_vec(2) * 16.d0
    B_ripple(1) = Btf0*exp(ripple_params%C0 + ripple_params%C1*R_vec(3)**2 + ripple_params%C2*R_vec(3)**4 + &
      R_vec(1)*(ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4))*Sin(psi)
    B_ripple(2) = Btf0*exp(ripple_params%A0 + ripple_params%A1*R_vec(3)**2 + ripple_params%A2*R_vec(3)**4 + &
                    R_vec(1)*(ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4))*Cos(psi)
    B_ripple(3) = Btf0*exp(ripple_params%E0 + ripple_params%E1*R_vec(3)**2 + R_vec(1)*(ripple_params%K0 + &
                  ripple_params%K1*R_vec(3)**2))*R_vec(3)*Sin(psi)
  end subroutine get_ripple

  subroutine get_ripple_w_grad(R_vec, B_ripple, grad_B_r_ripple, grad_B_t_ripple, grad_B_z_ripple)
    USE f90_kind
    ! Comment the following line for standalone usage
    ! use fmece_global_params,        only: grad_type
    use constants, only : pi
    implicit none
    real(rkind), dimension(:)  , intent(in)   :: R_vec ! R, psi, z; psi in range of [0, 2 N pi[
    real(rkind), dimension(:)  , intent(out)  :: B_ripple
    type(grad_type), intent(out)              :: grad_B_r_ripple, grad_B_t_ripple, grad_B_z_ripple
    real(rkind)                               :: Btf0, psi
    if(R_vec(1) < ripple_params%R_min .or. R_vec(1) > ripple_params%R_max) then
      B_ripple(:) = 0.d0
      grad_B_r_ripple%dR = 0.d0
      grad_B_r_ripple%dphi = 0.d0
      grad_B_r_ripple%dz = 0.d0
      grad_B_t_ripple%dR = 0.d0
      grad_B_t_ripple%dphi = 0.d0
      grad_B_t_ripple%dz = 0.d0
      grad_B_z_ripple%dR = 0.d0
      grad_B_z_ripple%dphi = 0.d0
      grad_B_z_ripple%dz = 0.d0
      return
    end if
    Btf0 = ripple_params%Btf0 / R_vec(1) * ripple_params%R0
    ! 16 coils at ASDEX Upgrade, hence we have also 16 periods in the ripple function when going once around the torus
    psi = R_vec(2) * 16.d0
    B_ripple(1) = Btf0*exp(ripple_params%C0 + ripple_params%C1*R_vec(3)**2 + ripple_params%C2*R_vec(3)**4 + &
      R_vec(1)*(ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4))*Sin(psi)
    B_ripple(2) = Btf0*exp(ripple_params%A0 + ripple_params%A1*R_vec(3)**2 + ripple_params%A2*R_vec(3)**4 + &
                    R_vec(1)*(ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4))*Cos(psi)
    B_ripple(3) = Btf0*exp(ripple_params%E0 + ripple_params%E1*R_vec(3)**2 + R_vec(1)*(ripple_params%K0 + &
                  ripple_params%K1*R_vec(3)**2))*R_vec(3)*Sin(psi)
    grad_B_r_ripple%dR =   Btf0*exp(ripple_params%C0 + ripple_params%C1*R_vec(3)**2 + ripple_params%C2*R_vec(3)**4 + &
                           R_vec(1)*(ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4)) * &
                           (ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4)*Sin(psi)
    grad_B_r_ripple%dR = grad_B_r_ripple%dR - B_ripple(1) / R_vec(1)
    grad_B_r_ripple%dphi = Btf0*exp(ripple_params%C0 + ripple_params%C1*R_vec(3)**2 + ripple_params%C2*R_vec(3)**4 + &
                           R_vec(1)*(ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4))*Cos(psi)
    grad_B_r_ripple%dphi = grad_B_r_ripple%dphi* 16.d0
    grad_B_r_ripple%dz = Btf0*exp(ripple_params%C0 + ripple_params%C1*R_vec(3)**2 + ripple_params%C2*R_vec(3)**4 + &
                         R_vec(1)*(ripple_params%D0 + ripple_params%D1*R_vec(3)**2 + ripple_params%D2*R_vec(3)**4))* &
                         (2*ripple_params%C1*R_vec(3) + 4*ripple_params%C2*R_vec(3)**3 + R_vec(1)*(2*ripple_params%D1*R_vec(3) + &
                         4*ripple_params%D2*R_vec(3)**3))*Sin(psi)
    grad_B_t_ripple%dR = Btf0*exp(ripple_params%A0 + ripple_params%A1*R_vec(3)**2 + ripple_params%A2*R_vec(3)**4 + &
                         R_vec(1)*(ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4))* &
                         (ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4)*Cos(psi)
    grad_B_t_ripple%dR = grad_B_t_ripple%dR - B_ripple(2) / R_vec(1)
    grad_B_t_ripple%dphi =-(Btf0*exp(ripple_params%A0 + ripple_params%A1*R_vec(3)**2 + ripple_params%A2*R_vec(3)**4 + &
                          R_vec(1)*(ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4))*Sin(psi))
    grad_B_t_ripple%dphi = grad_B_t_ripple%dphi* 16.d0
    grad_B_t_ripple%dz =  Btf0*exp(ripple_params%A0 + ripple_params%A1*R_vec(3)**2 + ripple_params%A2*R_vec(3)**4 + &
                          R_vec(1)*(ripple_params%B0 + ripple_params%B1*R_vec(3)**2 + ripple_params%B2*R_vec(3)**4))* &
                          (2*ripple_params%A1*R_vec(3) + 4*ripple_params%A2*R_vec(3)**3 + R_vec(1)*(2*ripple_params%B1*R_vec(3) + &
                          4*ripple_params%B2*R_vec(3)**3))*Cos(psi)
    grad_B_z_ripple%dR = Btf0*exp(ripple_params%E0 + ripple_params%E1*R_vec(3)**2 + &
                         R_vec(1)*(ripple_params%K0 + ripple_params%K1*R_vec(3)**2))*R_vec(3)* &
                         (ripple_params%K0 + ripple_params%K1*R_vec(3)**2)*Sin(psi)
    grad_B_z_ripple%dR = grad_B_z_ripple%dR - B_ripple(3) / R_vec(1)
    grad_B_z_ripple%dphi = Btf0*exp(ripple_params%E0 + ripple_params%E1*R_vec(3)**2 + R_vec(1)*(ripple_params%K0 + &
                           ripple_params%K1*R_vec(3)**2))*R_vec(3)*Cos(psi)
    grad_B_z_ripple%dphi = grad_B_z_ripple%dphi* 16.d0
    grad_B_z_ripple%dz = Btf0*exp(ripple_params%E0 + ripple_params%E1 * R_vec(3)**2 + R_vec(1)*(ripple_params%K0 + &
                         ripple_params%K1*R_vec(3)**2))*Sin(psi) + &
                         Btf0*exp(ripple_params%E0 + ripple_params%E1 * R_vec(3)**2 + R_vec(1)*(ripple_params%K0 + &
                         ripple_params%K1*R_vec(3)**2))*R_vec(3)*(2*ripple_params%E1*R_vec(3) + &
                         2*ripple_params%K1*R_vec(1)*R_vec(3))* Sin(psi)
  end subroutine get_ripple_w_grad
!
! To validate the gradient of the ripple the following routine calculates the finite differences
! Uncomment if you wish to test the gradients
  subroutine validate_ripple_grad()
  use f90_kind
  !use fmece_global_params,        only: grad_type
  use constants, only: pi
  implicit none
  real(rkind), dimension(3)                 :: R_vec, B_ripple
  real(rkind), dimension(4,3)               :: B_w_ripple_aux
  real(rkind), dimension(4,3)               :: aux_R
  type(grad_type)                           :: grad_B_r_ripple, grad_B_t_ripple, grad_B_z_ripple
  real(rkind)                               :: R0, Btf0, h_x
  integer(ikind)                            :: i, j
  Btf0 = - 2.5d0
  R0 = 1.67d0
  R_vec(1) = 2.130
  R_vec(2) = (16.d0 - 8.5d0) * 22.5 / 180.0 * pi
  R_vec(3) = 0.035d0
  h_x = 1.d-4
  call get_ripple_w_grad(R_vec, B_ripple, grad_B_r_ripple, grad_B_t_ripple, grad_B_z_ripple)
    do i =1,3
      do j = 1 , 4
        aux_R(j,:) = R_vec(:)
        if(j < 3) aux_R(j,i) = R_vec(i) + (3 - j) * h_x !aux_x(1) = x + 2*h, aux_2(2) = x + h
        if(j >= 3) aux_R(j,i) = R_vec(i) + (2 - j) * h_x !aux_x(3) = x - h, aux_2(4) = x - 2*h
        call get_ripple(aux_R(j,:), B_w_ripple_aux(j,:))
      end do
      if(i == 1) then
        print*, "dBr_rip/dR ana", grad_B_r_ripple%dR
      else if(i == 2) then
        print*, "dBr_rip/dphi ana", grad_B_r_ripple%dphi
      else
        print*, "dBr_rip/dz ana", grad_B_r_ripple%dz
      end if
      print*, "dBr_rip/dR num fwd. :", (B_w_ripple_aux(2,1) - B_ripple(1))/h_x
      print*, "dBr_rip/dR num ctr. :" , (-B_w_ripple_aux(1,1) + 8.d0 *  B_w_ripple_aux(2,1) - &
                          8.d0 *  B_w_ripple_aux(3,1) +  B_w_ripple_aux(4,1)) / (12.d0 *h_x)
      if(i == 1) then
        print*, "dBt_rip/dR ana", grad_B_t_ripple%dR
      else if(i == 2) then
        print*, "dBt_rip/dphi ana", grad_B_t_ripple%dphi
      else
        print*, "dBt_rip/dz ana", grad_B_t_ripple%dz
      end if
      print*, "dBt_rip/dR num fwd. :", (B_w_ripple_aux(2,2) - B_ripple(2))/h_x
      print*, "dBt_rip/dR num ctr. :" , (-B_w_ripple_aux(1,2) + 8.d0 *  B_w_ripple_aux(2,2) - &
                          8.d0 *  B_w_ripple_aux(3,2) +  B_w_ripple_aux(4,2)) / (12.d0 *h_x)
      if(i == 1) then
        print*, "dBz_rip/dR ana", grad_B_z_ripple%dR
      else if(i == 2) then
        print*, "dBz_rip/dphi ana", grad_B_z_ripple%dphi
      else
        print*, "dBz_rip/dz ana", grad_B_z_ripple%dz
      end if
      print*, "dBz_rip/dR num fwd. :", (B_w_ripple_aux(2,3) - B_ripple(3))/h_x
      print*, "dBz_rip/dR num ctr. :" , (-B_w_ripple_aux(1,3) + 8.d0 *  B_w_ripple_aux(2,3) - &
                          8.d0 *  B_w_ripple_aux(3,3) +  B_w_ripple_aux(4,3)) / (12.d0 *h_x)
      !h_x = 1.d-4
    end do
!    stop "ok?"
  end subroutine validate_ripple_grad

  end module ripple3d
