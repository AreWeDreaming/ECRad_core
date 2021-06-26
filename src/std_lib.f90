! Light version of the std_lib used for IDA and similar programs at ASDEX Upgrade
! Used to compile ECRad on non IPP infrastructure


!     Last change:  Rainer Fischer   19 Feb 98    2:24 pm
!#######################################################################
!  f90_kind.f
!
!  Module defining the KIND numbers for the NAGWare f90 Compiler.
!
!#######################################################################


module f90_kind
!=======================================================

    intrinsic kind,selected_int_kind,selected_real_kind  ! we use these,
    private   kind,selected_int_kind,selected_real_kind  ! but do not force

!-----------------------------------------------------------------------
! Indicator that the KIND= is not available for this compiler/host
!-----------------------------------------------------------------------

    integer, parameter :: not_available = -1

!-----------------------------------------------------------------------
! Real and Complex numbers
!-----------------------------------------------------------------------

!   Single precision
    integer, parameter :: single  = kind(1.0)
!   Double precision
    integer, parameter :: double  = kind(1.0d0)
!   Quadruple precision
    integer, parameter :: quad    = selected_real_kind(p=30)

!-----------------------------------------------------------------------
! Integers numbers
!-----------------------------------------------------------------------

!   Single byte integer
    integer, parameter :: int_def = kind(1)
!   Single byte integer
    integer, parameter :: int8    = selected_int_kind(2)
!   Two byte integer
    integer, parameter :: int16   = selected_int_kind(4)
!   Four byte integer
    integer, parameter :: int32   = selected_int_kind(9)
!   Eight byte integer
    integer, parameter :: int64   = selected_int_kind(18)

!-----------------------------------------------------------------------
! Logical values
!-----------------------------------------------------------------------

!   Single byte logical
    integer, parameter :: byte    = 1
!   Two byte logical
    integer, parameter :: twobyte = 2
!   Four byte logical
    integer, parameter :: word    = kind(.TRUE.)

!-----------------------------------------------------------------------
! Character type
!-----------------------------------------------------------------------

!   Normal single byte character (ASCII sequence)
    integer, parameter :: ascii   = kind('x')

!-----------------------------------------------------------------------
! Working precision
!-----------------------------------------------------------------------

    integer, parameter :: rkind   = double
    integer, parameter :: ikind   = int_def
    integer, parameter :: lkind   = word
    integer, parameter :: ckind   = ascii
integer, parameter :: I4B     = SELECTED_INT_KIND(9)
integer, parameter :: I2B     = SELECTED_INT_KIND(4)
integer, parameter :: I1B     = SELECTED_INT_KIND(2)
integer, parameter :: IDK     = KIND(1)              ! default kind
integer, parameter :: SP      = KIND(1.0)
integer, parameter :: DP      = KIND(1.0D0)
integer, parameter :: SPC     = KIND((1.0,1.0))
integer, parameter :: DPC     = KIND((1.0D0,1.0D0))
integer, parameter :: LGT     = KIND(.true.)

TYPE sprs2_sp
  INTEGER(I4B) :: n,len
  REAL(SP), DIMENSION(:), POINTER :: val
  INTEGER(I4B), DIMENSION(:), POINTER :: irow
  INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_sp
TYPE sprs2_dp
  INTEGER(I4B) :: n,len
  REAL(DP), DIMENSION(:), POINTER :: val
  INTEGER(I4B), DIMENSION(:), POINTER :: irow
  INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_dp

end module f90_kind

!     Last change:  Rainer Fischer  13 Nov 98    9:40
MODULE constants
!=======================================================
  USE f90_kind
  IMPLICIT NONE
!integer, parameter :: SP      = KIND(1.0)       ! already in f90_kind
!integer, parameter :: DP      = KIND(1.0D0)     ! already in f90_kind
  REAL(SP), PARAMETER :: zero_sp = 0.0_sp
  REAL(SP), PARAMETER :: one_sp  = 1.0_sp
  REAL(SP), PARAMETER :: two_sp  = 2.0_sp
  REAL(DP), PARAMETER :: zero    = 0.0_dp
  REAL(DP), PARAMETER :: one     = 1.0_dp
  REAL(DP), PARAMETER :: two     = 2.0_dp

  REAL(rkind), PARAMETER :: Pi    = 3.141592653589793238462643_rkind
  REAL(rkind), PARAMETER :: s2Pi  = 2.506628274_rkind
  REAL(rkind), PARAMETER :: sPi   = sqrt(Pi)
  REAL(rkind), PARAMETER :: twoPi = 6.28318530717958647692528676655900576839_rkind
  REAL(rkind), PARAMETER :: RtwoPi = 1.d0/twoPi
  REAL(rkind), PARAMETER :: RPi   = 1.d0/Pi
  REAL(rkind), PARAMETER :: ln_2  = 0.6931471805599453_rkind
  REAL(rkind), PARAMETER :: ln_10 = 2.302585092994046_rkind
  REAL(SP), PARAMETER :: PI_sp    = 3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2_sp  = 1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI_sp = 6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2_sp = 1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER_sp = 0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D     = 3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D   = 1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D  = 6.283185307179586476925286766559005768394_dp

  REAL(rkind), PARAMETER :: Ry                 = 13.6056981e+00_rkind        ! [eV]
  REAL(rkind), PARAMETER :: Rydberg_c          = 1.0973731534e+07_rkind      ! [m^-1]
  REAL(rkind), PARAMETER :: mass_e             = 9.1093897e-31_rkind         ! [kg]
  REAL(rkind), PARAMETER :: mass_p             = 1.6726231e-27_rkind         ! [kg]
  REAL(rkind), PARAMETER :: mass_n             = 1.6749286e-27_rkind         ! [kg]
  REAL(rkind), PARAMETER :: mass_u             = 1.6605402e-27_rkind         ! [kg]
  REAL(rkind), PARAMETER :: e0                 = 1.60217733e-19_rkind        ! [C]
  REAL(rkind), PARAMETER :: k_B                = 1.380658e-23_rkind          ! [J/K]
  REAL(rkind), PARAMETER :: c0                 = 2.99792458e+08_rkind        ! [m/s]
  REAL(rkind), PARAMETER :: eps0               = 8.854187817e-12_rkind       ! [F/m]
  REAL(rkind), PARAMETER :: mu0                = 1.25663706144e-6_rkind      ! [H/m]
  REAL(rkind), PARAMETER :: h_Planck           = 6.6260755e-34_rkind         ! [Js]
  REAL(rkind), PARAMETER :: h_Planck_bar       = 1.05457266e-34_rkind        ! [Js]
  REAL(rkind), PARAMETER :: norm_fall_acc      = 9.80665_rkind               ! [m/s^2]
  REAL(rkind), PARAMETER :: gravitation_c      = 6.67259e-11_rkind           ! [m^3/(kg*s^2)]
  REAL(rkind), PARAMETER :: molar_gas_c        = 8.314510e+00_rkind          ! [J/(mol*K)]
  REAL(rkind), PARAMETER :: avogadro_c         = 6.0221367e+23_rkind         ! [1/mol]
  REAL(rkind), PARAMETER :: loschmidt_c        = 2.686754e+25_rkind          ! [1/m^3]
  REAL(rkind), PARAMETER :: stefan_boltzmann_c = 5.67051e-08_rkind           ! [W/(m^2*K^4)]
  REAL(rkind), PARAMETER :: faraday_c          = 9.6485309e+04_rkind         ! [C/mol]
  REAL(rkind), PARAMETER :: fine_structure_c   = 7.29735308e-03_rkind        ! [1]
  REAL(rkind), PARAMETER :: magn_flux_quant    = 2.06783461e-15_rkind        ! [Wb]
  REAL(rkind), PARAMETER :: bohr_magneton      = 9.2740154e-24_rkind         ! [J/T]
  REAL(rkind), PARAMETER :: lambda_compton     = 2.42631058e-12_rkind        ! [m]
  REAL(rkind), PARAMETER :: classical_e_radius = 2.81794092e-15_rkind        ! [m]

  REAL(rkind), PARAMETER :: a_Bohr             = 0.529177249e-10_rkind       ! [m]
  REAL(rkind), PARAMETER :: hc                 = h_Planck*c0*1.e9_rkind/e0   ! [eV*nm]

  real(rkind), parameter :: refractive_indx_const = e0**2 / (4.d0*Pi**2*eps0*mass_e)     ! [m^3/s^2]

END MODULE constants

module nr_utils
! only assert_eq,nrerror, arth, geop needed
!--------------
use f90_kind
implicit none
INTEGER(ikind), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8
  INTEGER(ikind), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
INTEGER(ikind), PARAMETER :: NPAR_CUMSUM=16
INTEGER(ikind), PARAMETER :: NPAR_POLY=8


private

public array_copy,        &
       unit_matrix,       &
       upper_triangle,    &
       get_diag,          &
       imaxloc,           &
       swap,              &
       reallocate,        &
       assert,            &
       assert_eq,         &
       nrerror,           &
       arth,              &
       geop,              &
       outerdiff,         &
       outerand,          &
       outerprod,         &
       tridag,            &
       locate,            &
       cumsum,            &
       iminloc,           &
       diagadd,           &
       poly,              &
       diagmult

  INTERFACE array_copy
    MODULE PROCEDURE array_copy_r, array_copy_i
  END INTERFACE

  INTERFACE imaxloc
    MODULE PROCEDURE imaxloc_r, imaxloc_i
  END INTERFACE

  INTERFACE swap
    MODULE PROCEDURE swap_i, swap_r, swap_rv, swap_cv
  END INTERFACE

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_rv,reallocate_rm,&
      reallocate_iv,reallocate_im,reallocate_hv
  END INTERFACE

  INTERFACE assert
    MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE

  INTERFACE assert_eq
    MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE

  INTERFACE arth
    MODULE PROCEDURE arth_r, arth_i
  END INTERFACE

  INTERFACE geop
   MODULE PROCEDURE geop_r, geop_i, geop_c, geop_dv
  END INTERFACE

  INTERFACE outerdiff
    MODULE PROCEDURE outerdiff_r,outerdiff_i
  END INTERFACE

  INTERFACE cumsum
    MODULE PROCEDURE cumsum_r,cumsum_i
  END INTERFACE

  INTERFACE diagadd
    MODULE PROCEDURE diagadd_rv,diagadd_r
  END INTERFACE

  INTERFACE diagmult
    MODULE PROCEDURE diagmult_rv,diagmult_r
  END INTERFACE

  INTERFACE poly
    MODULE PROCEDURE poly_dd,poly_ddv,poly_msk_ddv
    !MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
    !  poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
  END INTERFACE

contains

  SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: src
  REAL(rkind), DIMENSION(:), INTENT(OUT) :: dest
  INTEGER(ikind), INTENT(OUT) :: n_copied, n_not_copied
  n_copied=min(size(src),size(dest))
  n_not_copied=size(src)-n_copied
  dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_r
!BL
  SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
  INTEGER(ikind), DIMENSION(:), INTENT(IN) :: src
  INTEGER(ikind), DIMENSION(:), INTENT(OUT) :: dest
  INTEGER(ikind), INTENT(OUT) :: n_copied, n_not_copied
  n_copied=min(size(src),size(dest))
  n_not_copied=size(src)-n_copied
  dest(1:n_copied)=src(1:n_copied)
  END SUBROUTINE array_copy_i
!BL
  SUBROUTINE unit_matrix(mat)
  REAL(rkind), DIMENSION(:,:), INTENT(OUT) :: mat
  INTEGER(ikind) :: i,n
  n=min(size(mat,1),size(mat,2))
  mat(:,:)=0.0_rkind
  do i=1,n
    mat(i,i)=1.0_rkind
  end do
  END SUBROUTINE unit_matrix
!BL
  FUNCTION upper_triangle(j,k,extra)
  INTEGER(ikind), INTENT(IN) :: j,k
  INTEGER(ikind), OPTIONAL, INTENT(IN) :: extra
  LOGICAL(lkind), DIMENSION(j,k) :: upper_triangle
  INTEGER(ikind) :: n
  n=0
  if (present(extra)) n=extra
  upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle
!BL
  FUNCTION get_diag(mat)
  REAL(rkind), DIMENSION(:,:), INTENT(IN) :: mat
  REAL(rkind), DIMENSION(size(mat,1)) :: get_diag
  INTEGER(ikind) :: j
  j=assert_eq2(size(mat,1),size(mat,2),'get_diag')
  do j=1,size(mat,1)
    get_diag(j)=mat(j,j)
  end do
  END FUNCTION get_diag
!BL
  FUNCTION reallocate_rv(p,n)
  REAL(rkind), DIMENSION(:), POINTER :: p, reallocate_rv
  INTEGER(ikind), INTENT(IN) :: n
  INTEGER(ikind) :: nold,ierr
  allocate(reallocate_rv(n),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_rv: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p)
  reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
  deallocate(p)
  END FUNCTION reallocate_rv
!BL
  FUNCTION reallocate_iv(p,n)
  INTEGER(ikind), DIMENSION(:), POINTER :: p, reallocate_iv
  INTEGER(ikind), INTENT(IN) :: n
  INTEGER(ikind) :: nold,ierr
  allocate(reallocate_iv(n),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_iv: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p)
  reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
  deallocate(p)
  END FUNCTION reallocate_iv
!BL
  FUNCTION reallocate_hv(p,n)
  CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
  INTEGER(ikind), INTENT(IN) :: n
  INTEGER(ikind) :: nold,ierr
  allocate(reallocate_hv(n),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_hv: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p)
  reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
  deallocate(p)
  END FUNCTION reallocate_hv
!BL
  FUNCTION reallocate_rm(p,n,m)
  REAL(rkind), DIMENSION(:,:), POINTER :: p, reallocate_rm
  INTEGER(ikind), INTENT(IN) :: n,m
  INTEGER(ikind) :: nold,mold,ierr
  allocate(reallocate_rm(n,m),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_rm: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p,1)
  mold=size(p,2)
  reallocate_rm(1:min(nold,n),1:min(mold,m))=&
    p(1:min(nold,n),1:min(mold,m))
  deallocate(p)
  END FUNCTION reallocate_rm
!BL
  FUNCTION reallocate_im(p,n,m)
  INTEGER(ikind), DIMENSION(:,:), POINTER :: p, reallocate_im
  INTEGER(ikind), INTENT(IN) :: n,m
  INTEGER(ikind) :: nold,mold,ierr
  allocate(reallocate_im(n,m),stat=ierr)
  if (ierr /= 0) call &
    nrerror('reallocate_im: problem in attempt to allocate memory')
  if (.not. associated(p)) RETURN
  nold=size(p,1)
  mold=size(p,2)
  reallocate_im(1:min(nold,n),1:min(mold,m))=&
    p(1:min(nold,n),1:min(mold,m))
  deallocate(p)
  END FUNCTION reallocate_im
!BL
  FUNCTION imaxloc_r(arr)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(ikind) :: imaxloc_r
  INTEGER(ikind), DIMENSION(1) :: imax
  imax=maxloc(arr(:))
  imaxloc_r=imax(1)
  END FUNCTION imaxloc_r
!BL
  FUNCTION imaxloc_i(iarr)
  INTEGER(ikind), DIMENSION(:), INTENT(IN) :: iarr
  INTEGER(ikind), DIMENSION(1) :: imax
  INTEGER(ikind) :: imaxloc_i
  imax=maxloc(iarr(:))
  imaxloc_i=imax(1)
  END FUNCTION imaxloc_i
!BL
  SUBROUTINE swap_i(a,b)
  INTEGER(ikind), INTENT(INOUT) :: a,b
  INTEGER(ikind) :: dum
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap_i
!BL
  SUBROUTINE swap_r(a,b)
  REAL(rkind), INTENT(INOUT) :: a,b
  REAL(rkind) :: dum
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap_r
!BL
  SUBROUTINE swap_rv(a,b)
  REAL(rkind), DIMENSION(:), INTENT(INOUT) :: a,b
  REAL(rkind), DIMENSION(SIZE(a)) :: dum
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap_rv
!BL
  SUBROUTINE swap_cv(a,b)
  complex(rkind), DIMENSION(:), INTENT(INOUT) :: a,b
  complex(rkind), DIMENSION(SIZE(a)) :: dum
  dum=a
  a=b
  b=dum
  END SUBROUTINE swap_cv
!BL
  SUBROUTINE assert1(n1,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1
  if (.not. n1) then
    write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
    STOP 'program terminated by assert1'
  end if
  END SUBROUTINE assert1
!BL
  SUBROUTINE assert2(n1,n2,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2
  if (.not. (n1 .and. n2)) then
    write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
    STOP 'program terminated by assert2'
  end if
  END SUBROUTINE assert2
!BL
  SUBROUTINE assert3(n1,n2,n3,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2,n3
  if (.not. (n1 .and. n2 .and. n3)) then
    write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
    STOP 'program terminated by assert3'
  end if
  END SUBROUTINE assert3
!BL
  SUBROUTINE assert4(n1,n2,n3,n4,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, INTENT(IN) :: n1,n2,n3,n4
  if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
    write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
    STOP 'program terminated by assert4'
  end if
  END SUBROUTINE assert4
!BL
  SUBROUTINE assert_v(n,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  LOGICAL, DIMENSION(:), INTENT(IN) :: n
  if (.not. all(n)) then
    write (*,*) 'nrerror: an assertion failed with this tag:', &
      string
    STOP 'program terminated by assert_v'
  end if
  END SUBROUTINE assert_v
!BL
  FUNCTION assert_eq2(n1,n2,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, INTENT(IN) :: n1,n2
  INTEGER :: assert_eq2
  if (n1 == n2) then
    assert_eq2=n1
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      string
    STOP 'program terminated by assert_eq2'
  end if
  END FUNCTION assert_eq2
!BL
  FUNCTION assert_eq3(n1,n2,n3,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, INTENT(IN) :: n1,n2,n3
  INTEGER :: assert_eq3
  if (n1 == n2 .and. n2 == n3) then
    assert_eq3=n1
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      string
    STOP 'program terminated by assert_eq3'
  end if
  END FUNCTION assert_eq3
!BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, INTENT(IN) :: n1,n2,n3,n4
  INTEGER :: assert_eq4
  if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
    assert_eq4=n1
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      string
    STOP 'program terminated by assert_eq4'
  end if
  END FUNCTION assert_eq4
!BL
  FUNCTION assert_eqn(nn,string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  INTEGER, DIMENSION(:), INTENT(IN) :: nn
  INTEGER :: assert_eqn
  if (all(nn(2:) == nn(1))) then
    assert_eqn=nn(1)
  else
    write (*,*) 'nrerror: an assert_eq failed with this tag:', &
      string
    STOP 'program terminated by assert_eqn'
  end if
  END FUNCTION assert_eqn
!BL
  SUBROUTINE nrerror(string)
  CHARACTER(LEN=*), INTENT(IN) :: string
  write (*,*) 'nrerror: ',string
  STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
!BL
  FUNCTION arth_r(first,increment,n)
  REAL(rkind), INTENT(IN) :: first,increment
  INTEGER(ikind), INTENT(IN) :: n
  REAL(rkind), DIMENSION(n) :: arth_r
  INTEGER(ikind) :: k,k2
  REAL(rkind) :: temp
  if (n > 0) arth_r(1)=first
  if (n <= NPAR_ARTH) then
    do k=2,n
      arth_r(k)=arth_r(k-1)+increment
    end do
  else
    do k=2,NPAR2_ARTH
      arth_r(k)=arth_r(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
      if (k >= n) exit
      k2=k+k
      arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
      temp=temp+temp
      k=k2
    end do
  end if
  END FUNCTION arth_r
!BL
  FUNCTION arth_i(first,increment,n)
  INTEGER(ikind), INTENT(IN) :: first,increment,n
  INTEGER(ikind), DIMENSION(n) :: arth_i
  INTEGER(ikind) :: k,k2,temp
  if (n > 0) arth_i(1)=first
  if (n <= NPAR_ARTH) then
    do k=2,n
      arth_i(k)=arth_i(k-1)+increment
    end do
  else
    do k=2,NPAR2_ARTH
      arth_i(k)=arth_i(k-1)+increment
    end do
    temp=increment*NPAR2_ARTH
    k=NPAR2_ARTH
    do
      if (k >= n) exit
      k2=k+k
      arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
      temp=temp+temp
      k=k2
    end do
  end if
  END FUNCTION arth_i
!BL
  FUNCTION geop_r(first,factor,n)
  REAL(rkind), INTENT(IN) :: first,factor
  INTEGER(ikind), INTENT(IN) :: n
  REAL(rkind), DIMENSION(n) :: geop_r
  INTEGER(ikind) :: k,k2
  REAL(rkind) :: temp
  if (n > 0) geop_r(1)=first
  if (n <= NPAR_GEOP) then
    do k=2,n
      geop_r(k)=geop_r(k-1)*factor
    end do
  else
    do k=2,NPAR2_GEOP
      geop_r(k)=geop_r(k-1)*factor
    end do
    temp=factor**NPAR2_GEOP
    k=NPAR2_GEOP
    do
      if (k >= n) exit
      k2=k+k
      geop_r(k+1:min(k2,n))=temp*geop_r(1:min(k,n-k))
      temp=temp*temp
      k=k2
    end do
  end if
  END FUNCTION geop_r
!BL
  FUNCTION geop_i(first,factor,n)
  INTEGER(ikind), INTENT(IN) :: first,factor,n
  INTEGER(ikind), DIMENSION(n) :: geop_i
  INTEGER(ikind) :: k,k2,temp
  if (n > 0) geop_i(1)=first
  if (n <= NPAR_GEOP) then
    do k=2,n
      geop_i(k)=geop_i(k-1)*factor
    end do
  else
    do k=2,NPAR2_GEOP
      geop_i(k)=geop_i(k-1)*factor
    end do
    temp=factor**NPAR2_GEOP
    k=NPAR2_GEOP
    do
      if (k >= n) exit
      k2=k+k
      geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
      temp=temp*temp
      k=k2
    end do
  end if
  END FUNCTION geop_i
!BL
  FUNCTION geop_c(first,factor,n)
  COMPLEX(rkind), INTENT(IN) :: first,factor
  INTEGER(ikind), INTENT(IN) :: n
  COMPLEX(rkind), DIMENSION(n) :: geop_c
  INTEGER(ikind) :: k,k2
  COMPLEX(rkind) :: temp
  if (n > 0) geop_c(1)=first
  if (n <= NPAR_GEOP) then
    do k=2,n
      geop_c(k)=geop_c(k-1)*factor
    end do
  else
    do k=2,NPAR2_GEOP
      geop_c(k)=geop_c(k-1)*factor
    end do
    temp=factor**NPAR2_GEOP
    k=NPAR2_GEOP
    do
      if (k >= n) exit
      k2=k+k
      geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
      temp=temp*temp
      k=k2
    end do
  end if
  END FUNCTION geop_c
!BL
  FUNCTION geop_dv(first,factor,n)
  REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
  INTEGER(ikind), INTENT(IN) :: n
  REAL(DP), DIMENSION(size(first),n) :: geop_dv
  INTEGER(ikind) :: k,k2
  REAL(DP), DIMENSION(size(first)) :: temp
  if (n > 0) geop_dv(:,1)=first(:)
  if (n <= NPAR_GEOP) then
    do k=2,n
      geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
    end do
  else
    do k=2,NPAR2_GEOP
      geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
    end do
    temp=factor**NPAR2_GEOP
    k=NPAR2_GEOP
    do
      if (k >= n) exit
      k2=k+k
      geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
        spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
      temp=temp*temp
      k=k2
    end do
  end if
  END FUNCTION geop_dv
!BL
  FUNCTION outerdiff_r(a,b)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: a,b
  REAL(rkind), DIMENSION(size(a),size(b)) :: outerdiff_r
  outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
    spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
!BL
  FUNCTION outerdiff_i(a,b)
  INTEGER(ikind), DIMENSION(:), INTENT(IN) :: a,b
  INTEGER(ikind), DIMENSION(size(a),size(b)) :: outerdiff_i
  outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
    spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i
!BL
  FUNCTION outerand(a,b)
  LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
  LOGICAL, DIMENSION(size(a),size(b)) :: outerand
  outerand = spread(a,dim=2,ncopies=size(b)) .and. &
    spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand
!BL
  FUNCTION outerprod(a,b)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: a,b
  REAL(rkind), DIMENSION(size(a),size(b)) :: outerprod
  outerprod = spread(a,dim=2,ncopies=size(b)) * &
    spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod
!BL

  SUBROUTINE tridag(a,b,c,r,u)
    ! tridag_ser
  IMPLICIT NONE
  REAL(rkind), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(rkind), DIMENSION(:), INTENT(OUT) :: u
  REAL(rkind), DIMENSION(size(b)) :: gam
  INTEGER(ikind) :: n,j
  REAL(rkind) :: bet
  n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag')
  bet=b(1)
  if (bet == 0.0) call nrerror('tridag: Error at code stage 1')
  u(1)=r(1)/bet
  do j=2,n
    gam(j)=c(j-1)/bet
    bet=b(j)-a(j-1)*gam(j)
    if (bet == 0.0) &
      call nrerror('tridag: Error at code stage 2')
    u(j)=(r(j)-a(j-1)*u(j-1))/bet
  end do
  do j=n-1,1,-1
    u(j)=u(j)-gam(j+1)*u(j+1)
  end do
  END SUBROUTINE tridag

  FUNCTION locate(xx,x)
  IMPLICIT NONE
  REAL(rkind), DIMENSION(:), INTENT(IN) :: xx
  REAL(rkind), INTENT(IN) :: x
  INTEGER(ikind) :: locate
  INTEGER(ikind) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
    if (ju-jl <= 1) exit
    jm=(ju+jl)/2
    if (ascnd .eqv. (x >= xx(jm))) then
      jl=jm
    else
      ju=jm
    end if
  end do
  if (x == xx(1)) then
    locate=1
  else if (x == xx(n)) then
    locate=n-1
  else
    locate=jl
  end if
  END FUNCTION locate

  RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
  IMPLICIT NONE
  REAL(rkind), DIMENSION(:), INTENT(IN) :: arr
  REAL(rkind), OPTIONAL, INTENT(IN) :: seed
  REAL(rkind), DIMENSION(size(arr)) :: ans
  INTEGER(ikind) :: n,j
  REAL(rkind) :: sd
  n=size(arr)
  if (n == 0_ikind) RETURN
  sd=0.0_rkind
  if (present(seed)) sd=seed
  ans(1)=arr(1)+sd
  if (n < NPAR_CUMSUM) then
    do j=2,n
      ans(j)=ans(j-1)+arr(j)
    end do
  else
    ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
    ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
  end if
  END FUNCTION cumsum_r

  RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
  IMPLICIT NONE
  INTEGER(ikind), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(ikind), OPTIONAL, INTENT(IN) :: seed
  INTEGER(ikind), DIMENSION(size(arr)) :: ans
  INTEGER(ikind) :: n,j,sd
  n=size(arr)
  if (n == 0_ikind) RETURN
  sd=0_ikind
  if (present(seed)) sd=seed
  ans(1)=arr(1)+sd
  if (n < NPAR_CUMSUM) then
    do j=2,n
      ans(j)=ans(j-1)+arr(j)
    end do
  else
    ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
    ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
  end if
  END FUNCTION cumsum_i

  FUNCTION iminloc(arr)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(ikind), DIMENSION(1) :: imin
  INTEGER(ikind) :: iminloc
  imin=minloc(arr(:))
  iminloc=imin(1)
  END FUNCTION iminloc

  SUBROUTINE diagadd_rv(mat,diag)
  REAL(rkind), DIMENSION(:,:), INTENT(INOUT) :: mat
  REAL(rkind), DIMENSION(:), INTENT(IN) :: diag
  INTEGER(ikind) :: j,n
  n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
  do j=1,n
    mat(j,j)=mat(j,j)+diag(j)
  end do
  END SUBROUTINE diagadd_rv

  SUBROUTINE diagadd_r(mat,diag)
  REAL(rkind), DIMENSION(:,:), INTENT(INOUT) :: mat
  REAL(rkind), INTENT(IN) :: diag
  INTEGER(ikind) :: j,n
  n = min(size(mat,1),size(mat,2))
  do j=1,n
    mat(j,j)=mat(j,j)+diag
  end do
  END SUBROUTINE diagadd_r


  SUBROUTINE diagmult_rv(mat,diag)
  REAL(rkind), DIMENSION(:,:), INTENT(INOUT) :: mat
  REAL(rkind), DIMENSION(:), INTENT(IN) :: diag
  INTEGER(ikind) :: j,n
  n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
  do j=1,n
    mat(j,j)=mat(j,j)*diag(j)
  end do
  END SUBROUTINE diagmult_rv

  SUBROUTINE diagmult_r(mat,diag)
  REAL(rkind), DIMENSION(:,:), INTENT(INOUT) :: mat
  REAL(rkind), INTENT(IN) :: diag
  INTEGER(ikind) :: j,n
  n = min(size(mat,1),size(mat,2))
  do j=1,n
    mat(j,j)=mat(j,j)*diag
  end do
  END SUBROUTINE diagmult_r



  FUNCTION poly_dd(x,coeffs)
  REAL(rkind), INTENT(IN) :: x
  REAL(rkind), DIMENSION(:), INTENT(IN) :: coeffs
  REAL(rkind) :: poly_dd
  REAL(rkind) :: pow
  REAL(rkind), DIMENSION(:), ALLOCATABLE :: vec
  INTEGER(ikind) :: i,n,nn
  n=size(coeffs)
  if (n <= 0) then
    poly_dd=0.0_rkind
  else if (n < NPAR_POLY) then
    poly_dd=coeffs(n)
    do i=n-1,1,-1
      poly_dd=x*poly_dd+coeffs(i)
    end do
  else
    allocate(vec(n+1))
    pow=x
    vec(1:n)=coeffs
    do
      vec(n+1)=0.0_rkind
      nn=ishft(n+1,-1)
      vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
      if (nn == 1) exit
      pow=pow*pow
      n=nn
    end do
    poly_dd=vec(1)
    deallocate(vec)
  end if
  END FUNCTION poly_dd

  FUNCTION poly_ddv(x,coeffs)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: coeffs,x
  REAL(rkind), DIMENSION(size(x)) :: poly_ddv
  INTEGER(ikind) :: i,n,m
  m=size(coeffs)
  n=size(x)
  if (m <= 0) then
    poly_ddv=0.0_rkind
  else if (m < n .or. m < NPAR_POLY) then
    poly_ddv=coeffs(m)
    do i=m-1,1,-1
      poly_ddv=x*poly_ddv+coeffs(i)
    end do
  else
    do i=1,n
      poly_ddv(i)=poly_dd(x(i),coeffs)
    end do
  end if
  END FUNCTION poly_ddv

  FUNCTION poly_msk_ddv(x,coeffs,mask)
  REAL(rkind), DIMENSION(:), INTENT(IN) :: coeffs,x
  LOGICAL, DIMENSION(:), INTENT(IN) :: mask
  REAL(rkind), DIMENSION(size(x)) :: poly_msk_ddv
  poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  END FUNCTION poly_msk_ddv

end module nr_utils

module nr_special_funcs
!--------------
use f90_kind
implicit none

private

public :: dawson_s

contains

!------------------------------------------------------------------------------


  FUNCTION dawson_s(x)
  !USE nrtype; USE nrutil, ONLY : arth,geop
  use nr_utils, only: arth, geop
  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: x
  REAL(rkind) :: dawson_s
  INTEGER(ikind), PARAMETER :: NMAX=7   ! original 6
  REAL(rkind), PARAMETER :: H=0.4_rkind,A1=2.0_rkind/3.0_rkind,A2=0.4_rkind,&
    A3=2.0_rkind/7.0_rkind
  INTEGER(ikind) :: i,n0
  REAL(rkind) :: ec,x2,xp,xx
  REAL(rkind), DIMENSION(NMAX) :: d1,d2,e1
  REAL(rkind), DIMENSION(NMAX), SAVE :: c=(/ (0.0_rkind,i=1,NMAX) /)
  if (c(1) == 0.0) c(1:NMAX)=exp(-(arth(1,2,NMAX)*H)**2)
  if (abs(x) < 0.2_rkind) then
    x2=x**2
    dawson_s=x*(1.0_rkind-A1*x2*(1.0_rkind-A2*x2*(1.0_rkind-A3*x2)))
  else
    xx=abs(x)
    n0=2*nint(0.5_rkind*xx/H)
    xp=xx-real(n0,rkind)*H
    ec=exp(2.0_rkind*xp*H)
    d1=arth(n0+1,2,NMAX)
    d2=arth(n0-1,-2,NMAX)
    e1=geop(ec,ec**2,NMAX)
    dawson_s=0.5641895835477563_rkind*sign(exp(-xp**2),x)*&
      sum(c*(e1/d1+1.0_rkind/(d2*e1)))
  end if
  END FUNCTION dawson_s

end module nr_special_funcs

MODULE f1dim_mod
  USE f90_kind
  INTEGER(ikind) :: ncom
  REAL(rkind), DIMENSION(:), POINTER :: pcom,xicom
#ifdef OMP
  !$OMP THREADPRIVATE(ncom,pcom,xicom)
#endif
CONTAINS
!BL
  FUNCTION f1dim(x,func)
  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: x
  REAL(rkind) :: f1dim
  INTERFACE
    FUNCTION func(x)
    USE f90_kind
    REAL(rkind), DIMENSION(:), INTENT(IN) :: x
    REAL(rkind) :: func
    END FUNCTION func
  END INTERFACE
  REAL(rkind), DIMENSION(:), ALLOCATABLE :: xt
  allocate(xt(ncom))
  xt(:)=pcom(:)+x*xicom(:)
  f1dim=func(xt)
  deallocate(xt)
  END FUNCTION f1dim
END MODULE f1dim_mod

Module nr_mod
public powell
private linmin, mnbrak2, brent2

contains


FUNCTION brent2(ax,bx,cx,func,tol,xmin,func2)
  USE f90_kind
  USE nr_utils, ONLY : nrerror
  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: ax,bx,cx,tol
  REAL(rkind), INTENT(OUT) :: xmin
  REAL(rkind) :: brent2
  INTERFACE
    FUNCTION func2(p)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), DIMENSION(:), INTENT(IN) :: p
    REAL(rkind) :: func2
    END FUNCTION func2
  END INTERFACE

  INTERFACE
    FUNCTION func(x,func2)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), INTENT(IN) :: x
    INTERFACE
      FUNCTION func2(p)
      USE f90_kind
      IMPLICIT NONE
      REAL(rkind), DIMENSION(:), INTENT(IN) :: p
      REAL(rkind) :: func2
      END FUNCTION func2
    END INTERFACE
    REAL(rkind) :: func
    END FUNCTION func
  END INTERFACE

  INTEGER(ikind), PARAMETER :: ITMAX=100
  REAL(rkind), PARAMETER :: CGOLD=0.3819660d0,ZEPS=1.0d-3*epsilon(ax)
  INTEGER(ikind) :: iter
  REAL(rkind) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
  a=min(ax,cx)
  b=max(ax,cx)
  v=bx
  w=v
  x=v
  e=0.0
  fx=func(x,func2)
  fv=fx
  fw=fx
  do iter=1,ITMAX
    xm=0.5d0*(a+b)
    tol1=tol*abs(x)+ZEPS
    tol2=2.0d0*tol1
    if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
      xmin=x
      brent2=fx
      RETURN
    end if
    if (abs(e) > tol1) then
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0d0*(q-r)
      if (q > 0.0) p=-p
      q=abs(q)
      etemp=e
      e=d
      if (abs(p) >= abs(0.5d0*q*etemp) .or. &
        p <= q*(a-x) .or. p >= q*(b-x)) then
        e=merge(a-x,b-x, x >= xm )
        d=CGOLD*e
      else
        d=p/q
        u=x+d
        if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
      end if
    else
      e=merge(a-x,b-x, x >= xm )
      d=CGOLD*e
    end if
    u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
    fu=func(u,func2)
    if (fu <= fx) then
      if (u >= x) then
        a=x
      else
        b=x
      end if
      call shft(v,w,x,u)
      call shft(fv,fw,fx,fu)
    else
      if (u < x) then
        a=u
      else
        b=u
      end if
      if (fu <= fw .or. w == x) then
        v=w
        fv=fw
        w=u
        fw=fu
      else if (fu <= fv .or. v == x .or. v == w) then
        v=u
        fv=fu
      end if
    end if
  end do
  call nrerror('brent2: exceed maximum iterations')
  CONTAINS
!BL
  SUBROUTINE shft(a,b,c,d)
  REAL(rkind), INTENT(OUT) :: a
  REAL(rkind), INTENT(INOUT) :: b,c
  REAL(rkind), INTENT(IN) :: d
  a=b
  b=c
  c=d
  END SUBROUTINE shft
  END FUNCTION brent2

  SUBROUTINE mnbrak2(ax,bx,cx,fa,fb,fc,func,func2)
  USE f90_kind
  USE nr_utils, ONLY : swap
  IMPLICIT NONE
  REAL(rkind), INTENT(INOUT) :: ax,bx
  REAL(rkind), INTENT(OUT) :: cx,fa,fb,fc
  INTERFACE
    FUNCTION func2(p)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), DIMENSION(:), INTENT(IN) :: p
    REAL(rkind) :: func2
    END FUNCTION func2
  END INTERFACE

  INTERFACE
    FUNCTION func(x,func2)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), INTENT(IN) :: x
    INTERFACE
      FUNCTION func2(p)
      USE f90_kind
      IMPLICIT NONE
      REAL(rkind), DIMENSION(:), INTENT(IN) :: p
      REAL(rkind) :: func2
      END FUNCTION func2
    END INTERFACE
    REAL(rkind) :: func
    END FUNCTION func
  END INTERFACE
  REAL(rkind), PARAMETER :: GOLD=1.618034_rkind,GLIMIT=100.0_rkind,TINY=1.0e-20_rkind
  REAL(rkind) :: fu,q,r,u,ulim
  fa=func(ax,func2)
  fb=func(bx,func2)
  if (fb > fa) then
    call swap(ax,bx)
    call swap(fa,fb)
  end if
  cx=bx+GOLD*(bx-ax)
  fc=func(cx,func2)
  do
    if (fb < fc) RETURN
    r=(bx-ax)*(fb-fc)
    q=(bx-cx)*(fb-fa)
    u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_rkind*sign(max(abs(q-r),TINY),q-r))
    ulim=bx+GLIMIT*(cx-bx)
    if ((bx-u)*(u-cx) > 0.0) then
      fu=func(u,func2)
      if (fu < fc) then
        ax=bx
        fa=fb
        bx=u
        fb=fu
        RETURN
      else if (fu > fb) then
        cx=u
        fc=fu
        RETURN
      end if
      u=cx+GOLD*(cx-bx)
      fu=func(u,func2)
    else if ((cx-u)*(u-ulim) > 0.0) then
      fu=func(u,func2)
      if (fu < fc) then
        bx=cx
        cx=u
        u=cx+GOLD*(cx-bx)
        call shft(fb,fc,fu,func(u,func2))
      end if
    else if ((u-ulim)*(ulim-cx) >= 0.0) then
      u=ulim
      fu=func(u,func2)
    else
      u=cx+GOLD*(cx-bx)
      fu=func(u,func2)
    end if
    call shft(ax,bx,cx,u)
    call shft(fa,fb,fc,fu)
  end do
  CONTAINS
!BL
  SUBROUTINE shft(a,b,c,d)
  REAL(rkind), INTENT(OUT) :: a
  REAL(rkind), INTENT(INOUT) :: b,c
  REAL(rkind), INTENT(IN) :: d
  a=b
  b=c
  c=d
  END SUBROUTINE shft
  END SUBROUTINE mnbrak2

SUBROUTINE linmin(p,xi,fret,func)
  USE f90_kind
  USE nr_utils, ONLY : assert_eq
! USE nr, ONLY : mnbrak,brent
  USE f1dim_mod
  IMPLICIT NONE
  REAL(rkind), INTENT(OUT) :: fret
  REAL(rkind), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
  INTERFACE
    FUNCTION func(p)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), DIMENSION(:), INTENT(IN) :: p
    REAL(rkind) :: func
    END FUNCTION func
  END INTERFACE
  REAL(rkind), PARAMETER :: TOL=1.0d-4
  REAL(rkind) :: ax,bx,fa,fb,fx,xmin,xx
  ncom=assert_eq(size(p),size(xi),'linmin')
  pcom=>p
  xicom=>xi
  ax=0.0
  xx=1.0
  call mnbrak2(ax,xx,bx,fa,fx,fb,f1dim,func)
  fret=brent2(ax,xx,bx,f1dim,TOL,xmin,func)
  if (xmin == 0.) xmin = 1.d-30   ! rrf: some small value for having a new search direction,
                                  !      otherwise the new search direction would be (0,..,0)
                                  !      (06.02.01)
  xi=xmin*xi
  p=p+xi
  END SUBROUTINE linmin

SUBROUTINE powell(p,xi,ftol,iter,fret,func)
  use f90_kind
  USE nr_utils, ONLY : assert_eq,nrerror
! USE nr, ONLY : linmin
  IMPLICIT NONE
  REAL(rkind), DIMENSION(:), INTENT(INOUT) :: p
  REAL(rkind), DIMENSION(:,:), INTENT(INOUT) :: xi
  INTEGER(ikind), INTENT(OUT) :: iter
  REAL(rkind), INTENT(IN) :: ftol
  REAL(rkind), INTENT(OUT) :: fret
  INTERFACE
    FUNCTION func(p)
    USE f90_kind
    IMPLICIT NONE
    REAL(rkind), DIMENSION(:), INTENT(IN) :: p
    REAL(rkind) :: func
    END FUNCTION func
  END INTERFACE
  INTEGER(ikind), PARAMETER :: ITMAX=200
  REAL(rkind), PARAMETER :: TINY=1.0d-25
  INTEGER(ikind) :: i,ibig,n
  REAL(rkind) :: del,fp,fptt,t
  REAL(rkind), DIMENSION(size(p)) :: pt,ptt,xit
  n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
  fret=func(p)
  pt(:)=p(:)
  iter=0
  do
    iter=iter+1
    fp=fret
    ibig=0
    del=0.0
    do i=1,n
      xit(:)=xi(:,i)
      fptt=fret
      call linmin(p,xit,fret,func)
      if (fptt-fret > del) then
        del=fptt-fret
        ibig=i
      end if
    end do
    if (2.0d0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN
    if (iter == ITMAX) call &
      nrerror('powell exceeding maximum iterations')
    ptt(:)=2.0_rkind*p(:)-pt(:)
    xit(:)=p(:)-pt(:)
    pt(:)=p(:)
    fptt=func(ptt)
    if (fptt >= fp) cycle
    t=2.0d0*(fp-2.0d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
    if (t >= 0.0) cycle
    call linmin(p,xit,fret,func)
    xi(:,ibig)=xi(:,n)
    xi(:,n)=xit(:)
  end do
  END SUBROUTINE powell
end module nr_mod

MODULE interpolation_routines
USE f90_kind
IMPLICIT NONE

PRIVATE
PUBLIC :: linear_interpolation

CONTAINS

SUBROUTINE linear_interpolation(xa, xe, ya, ye, x, y)
!====================================================
! linear interpolation
use f90_kind
implicit none
real(rkind), intent(in)  :: xa, xe, ya, ye, x
real(rkind), intent(out) :: y

y = ya + (ye - ya) * (x - xa) / (xe - xa)


END SUBROUTINE linear_interpolation

end MODULE interpolation_routines
