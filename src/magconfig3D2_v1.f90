!########################################################################

 MODULE magconfig3D

!########################################################################
!-------
 IMPLICIT NONE
!----------------- PRIVATE variables ------------------------------------
 CHARACTER(Len=120),PRIVATE :: format_conf   ! format of the magnetic config.
 CHARACTER(Len=120),PRIVATE :: name_conf     ! filename of the magnetic config.
!
 LOGICAL,           PRIVATE :: useMesh,useSymm
!
 INTEGER(8),   PRIVATE :: mconf=0            ! address of mconf object,
!                                            ! see cpp2for.cpp
 REAL(8),   PRIVATE :: hgrid, dR, dZ, dphi ! [m], [rad]
 REAL(8),   PRIVATE :: accbooz, tolharm
 real(8),  PRIVATE  :: comp_eps = 1.d-12
!
 INTEGER,    PRIVATE :: Nsym
 REAL(8), PRIVATE :: SSplus,SSmax
 LOGICAL,           PRIVATE :: trace_scrout=.true.

TYPE Scenario_type
 CHARACTER(Len=120) :: work_dir
 CHARACTER(Len=120) :: name_config
 CHARACTER(Len=120) :: format_config
 logical :: useMesh, useSymm

 real(8) :: splus, smax, hgrid, dphic, accbooz, tolharm
 integer(8) :: mag_conf_addr
END TYPE Scenario_type
#ifdef USE_3D

!########################################################################
!
 INTERFACE
   INTEGER(4) FUNCTION mcload(mconf,name) 
     INTEGER(8)   :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcLoad
   !
   INTEGER(4) FUNCTION mcLoadVMEC(mconf,name) 
     INTEGER(8)   :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcLoadVMEC
   !
   INTEGER(4) FUNCTION mcLoadLHD(mconf,name) 
     INTEGER(8)   :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcLoadLHD
   !
   INTEGER(4) FUNCTION mcLoadEFIT(mconf,name) 
     INTEGER(8)   :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcLoadEFIT
   !
   INTEGER(4) FUNCTION mcIsOK(mconf) 
     INTEGER(8) :: mconf
   END FUNCTION mcIsOK
   !
   INTEGER(4) FUNCTION mcIsMeshOK(mconf) 
     INTEGER(8) :: mconf
   END FUNCTION mcIsMeshOK
   !
   INTEGER(4) FUNCTION mcXYZisINSIDE(mconf,XYZ) 
     INTEGER(8) :: mconf
     REAL(8)    :: XYZ(3)
   END FUNCTION mcXYZisINSIDE
   !
   INTEGER(4) FUNCTION mcCYLisINSIDE(mconf,CYL) 
     INTEGER(8) :: mconf
     REAL(8)    :: CYL(3)
   END FUNCTION mcCYLisINSIDE
   !
   INTEGER(4) FUNCTION mcCreateMesh(mconf,dR,dZ,dphi) 
     INTEGER(8) :: mconf
     REAL(8)    :: dR,dZ,dphi
   END FUNCTION mcCreateMesh
   !
   INTEGER(4) FUNCTION mcCreateMeshUsingSymmetry(mconf,dR,dZ,dphi) 
     INTEGER(8) :: mconf
     REAL(8)    :: dR,dZ,dphi
   END FUNCTION mcCreateMeshUsingSymmetry
   !
   INTEGER(4) FUNCTION mcWrite(mconf,name) 
     INTEGER(8) :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcWrite
   !
   INTEGER(4) FUNCTION mcWriteMesh(mconf,name) 
     INTEGER(8) :: mconf
     CHARACTER(*) :: name
   END FUNCTION mcWriteMesh
   !
   INTEGER(4) FUNCTION mcNperiods(mconf) 
     INTEGER(8) :: mconf
   END FUNCTION mcNperiods
   !
   REAL(8) FUNCTION mcrminor(mconf) 
     INTEGER(8) :: mconf
   END FUNCTION mcrminor
   !
   REAL(8) FUNCTION mcR0(mconf) 
     INTEGER(8) :: mconf
   END FUNCTION mcR0
   !
   REAL(8) FUNCTION mcTorFlux(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcTorFlux
   !
   REAL(8) FUNCTION mcVprime(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcVprime
   !
   REAL(8) FUNCTION mcVolume(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcVolume
   !
   REAL(8) FUNCTION mcAverageCrosssectionArea(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcAverageCrosssectionArea
   !   !
   REAL(8) FUNCTION mcIota(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcIota
   !
   REAL(8) FUNCTION mcIotaPrime(mconf,S) 
     INTEGER(8) :: mconf
     REAL(8)    :: S
   END FUNCTION mcIotaPrime
   !
   REAL(8) FUNCTION mcRaxis(mconf,phi) 
     INTEGER(8) :: mconf
     REAL(8)    :: phi
   END FUNCTION mcRaxis
   !
   SUBROUTINE mcgetGradSxyz(mconf,XYZ,gradS) 
     INTEGER(8) :: mconf
     REAL(8)    :: XYZ(3),gradS(3)
   END SUBROUTINE mcgetGradSxyz
   !
   SUBROUTINE mcM3DgetGradSxyz(mconf,XYZ,gradS) 
     INTEGER(8) :: mconf
     REAL(8)    :: XYZ(3),gradS(3)
   END SUBROUTINE mcM3DgetGradSxyz
   !
   SUBROUTINE mcM3DgetSBxyz(mconf,xyz,S,B) 
     INTEGER(8) :: mconf
     REAL(8)    :: xyz(3),S,B(3)
   END SUBROUTINE mcM3DgetSBxyz
   !
   SUBROUTINE mcgetSBxyz(mconf,xyz,S,B) 
     INTEGER(8) :: mconf
     REAL(8)    :: xyz(3),S,B(3)
   END SUBROUTINE mcgetSBxyz
   !
   SUBROUTINE mcGetGradB(mconf,magCoord,gradB) 
     INTEGER(8) :: mconf
     REAL(8)    :: magCoord(3),gradB(3)
   END SUBROUTINE mcGetGradB
   !
   SUBROUTINE MCgetAllxyz(mconf,xyz,S,B,gradS,gradBi)
     INTEGER(8) :: mconf
     REAL(8)    :: xyz(3),S,B(3),gradS(3)
     REAL(8)    :: gradBi(9)
   END SUBROUTINE MCgetAllxyz
   !
   SUBROUTINE MCM3DgetAllxyz(mconf,xyz,S,B,gradS,gradBi)
     INTEGER(8) :: mconf
     REAL(8)    :: xyz(3),S,B(3),gradS(3)
     REAL(8)    :: gradBi(9)
   END SUBROUTINE MCM3DgetAllxyz
   !
   SUBROUTINE mcM3DSetSmax(mConf,ssMax) 
     INTEGER(8) :: mConf
     REAL(8)    :: ssMax
   END SUBROUTINE mcM3DSetSmax
   !
   REAL(8) FUNCTION mcM3DGetSmax(mConf) 
     INTEGER(8) :: mConf
   END FUNCTION mcM3DGetSmax
   !
   REAL(8) FUNCTION mcTorFlux2PolFlux(mconf,Stor) 
     INTEGER(8) :: mconf
     REAL(8)    :: Stor
   END FUNCTION mcTorFlux2PolFlux
   !
   SUBROUTINE mcTruncate(mConf,r) 
     INTEGER(8) :: mConf
     REAL(8)    :: r
   END SUBROUTINE mcTruncate
   !
   SUBROUTINE mcSetAccuracy(mConf,accuracy) 
     INTEGER(8) :: mConf
     REAL(8)    :: accuracy
   END SUBROUTINE mcSetAccuracy
   !
   SUBROUTINE mcSetIpol(mConf,Ipol) 
     INTEGER(8) :: mConf
     REAL(8)    :: Ipol  ! in Ampers
   END SUBROUTINE mcSetIpol
   !
   SUBROUTINE mcScaleB(mConf,factor) 
     INTEGER(8) :: mConf
     REAL(8)    :: factor
   END SUBROUTINE mcScaleB
   !
   SUBROUTINE mcSetB00(mConf,B00) 
     INTEGER(8) :: mConf
     REAL(8)    :: B00
   END SUBROUTINE mcSetB00
   !
   SUBROUTINE mcSetminB0(mConf,B00) 
     INTEGER(8) :: mConf
     REAL(8)    :: B00
   END SUBROUTINE mcSetminB0
   !
   SUBROUTINE mcSetB0(mConf,B0, angle) 
     INTEGER(8) :: mConf
     REAL(8)    :: B0, angle
   END SUBROUTINE mcSetB0
   !
   REAL(8) FUNCTION mcGetB00(mConf) 
     INTEGER(8) :: mConf
   END FUNCTION mcGetB00
   !
   REAL(8) FUNCTION mcBmn(mConf,m,n,S) 
     INTEGER(8) :: mConf
     INTEGER(4) :: m,n
     REAL(8)    :: S
   END FUNCTION mcBmn
   !
   SUBROUTINE mcGetBcyl(mConf,boozer,B) 
     INTEGER(8) :: mConf
     REAL(8)    :: boozer(3)  ! input: Boozer coordinates (s,theta,phi)
     REAL(8)    :: B          ! output Bcyl
   END SUBROUTINE mcGetBcyl
   !
   SUBROUTINE mcMag2Cyl(mConf,magCoord,RpZ) 
     INTEGER(8) :: mConf
     REAL(8)    :: RpZ(3)       ! input:  cyl. coordinate
     REAL(8)    :: magCoord(3)  ! output: magnetic coordinates (s,theta,phi)
   END SUBROUTINE mcMag2Cyl
   !
   SUBROUTINE mcCyl2Mag(mConf,RpZ,magCoord,Bfield) 
     INTEGER(8) :: mConf
     REAL(8)    :: RpZ(3)       ! input:  cyl. coordinate
     REAL(8)    :: magCoord(3)  ! output: magnetic coordinates (s,theta,phi)
     REAL(8)    :: Bfield(3)    ! output: Bfield=(Br,Bfi,Bz) in cyl. coordinates
   END SUBROUTINE mcCyl2Mag
   !
   !
   REAL(8) FUNCTION mcB(mConf,S,polAngle,torAngle) 
     INTEGER(8) :: mConf
     REAL(8)    :: S,polAngle,torAngle  ! input: Boozer coordinates
   END FUNCTION mcB
   !
 END INTERFACE
!########################################################################

 CONTAINS

!########################################################################

 SUBROUTINE MConf3D_Setup_Config(Scen)
!========================================================================
 IMPLICIT NONE
 TYPE(Scenario_type) :: Scen
!========================================================================
! CALL Setup_CompServ(Scen%work_dir)
!---
 name_conf   = Scen%name_config
 format_conf = Scen%format_config
 useMesh     = Scen%useMesh
 useSymm     = Scen%useSymm
!
 ssplus      = Scen%splus
 ssmax       = Scen%smax
 hgrid       = Scen%hgrid;    dR = hgrid;   dZ = hgrid
!
 dphi        = Scen%dphic
!
 accbooz     = Scen%accbooz
 tolharm     = Scen%tolharm
!
!========================================================================
 END SUBROUTINE MConf3D_Setup_Config

!########################################################################
!########################################################################
!########################################################################

 subroutine MConf3D_load_magconf_file(loaded, mcAddr)
!========================================================================
  INTEGER(8), intent(out) :: mcAddr
  INTEGER, intent(out)    :: loaded
  CHARACTER(Len=4) :: fmt
!========================================================================
  fmt = format_conf(1:4)
  mcAddr = mconf
  SELECT CASE (fmt)
  CASE('VMEC')
     loaded = mcLoadVMEC(mcAddr,name_conf)
   CASE('EFIT')
     loaded = mcLoadEFIT(mcAddr,name_conf)
   CASE('LHD ')
     loaded = mcLoadLHD(mcAddr,name_conf)
   CASE default
     loaded = mcLoad(mcAddr,name_conf)
  END SELECT
!========================================================================
 END subroutine MConf3D_load_magconf_file

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_Load_MagConfig(loaded, mConfAddr)
!========================================================================
 IMPLICIT NONE
 integer, intent(out)   :: loaded
 INTEGER(8),OPTIONAL, INTENT(inout) :: mConfAddr  ! address of C3dMesh-object 
 INTEGER :: i
!========================================================================
 IF(.NOT.PRESENT(mConfAddr)) THEN 
   call MConf3D_load_magconf_file(loaded, mconf)
 ELSEIF(mConfAddr==0) THEN 
   call MConf3D_load_magconf_file(loaded, mconf)
 ELSE
   mconf = mConfAddr
 ENDIF
 IF(PRESENT(mConfAddr)) mConfAddr = mconf

 CALL mcTruncate   (mconf,tolharm)
 CALL mcSetAccuracy(mconf,accbooz)
!---
 IF (trace_scrout) WRITE(*,*) 'Magnetic configuration is loaded.'
 IF (useMesh) THEN
   IF(mcIsMeshOK(mconf) /= 1) THEN ! test whether mconf contains mesh
     IF (trace_scrout) WRITE(*,*) '   generating mesh...'
     CALL mcM3DsetSmax(mconf,ssplus)
     IF(useSymm) THEN
       i = mcCreateMeshUsingSymmetry(mconf,dR,dZ,dphi)
     ELSE
       i = mcCreateMesh             (mconf,dR,dZ,dphi)
     ENDIF
     IF(i /= 1) THEN 
       useMesh = .FALSE.
       IF(trace_scrout) WRITE(*,*) 'Error creating mesh.'
     ELSE
       IF(trace_scrout) WRITE(*,*) 'Mesh is generated.'
     ENDIF  
   ENDIF
 ENDIF
!========================================================================
 END SUBROUTINE MConf3D_Load_MagConfig



!########################################################################

 FUNCTION isInsideLCMS(coord,XYZ) RESULT(isInside)
!========================================================================
!  function returns 'true' if the point lies inside LCMS
!========================================================================
 IMPLICIT NONE
 CHARACTER(Len=2), INTENT(in) :: coord
 REAL(8),        INTENT(in) :: XYZ(3)
 REAL(8) :: XYZ1(3)
 LOGICAL    :: isInside
!========================================================================
   XYZ1 = XYZ
   isInside=.FALSE.
   IF (coord == 'cy') THEN
     isInside = MCCYLISINSIDE(mconf,XYZ1) /= 0
   ELSEIF (coord == 'ca') THEN
     isInside = MCXYZISINSIDE(mconf,XYZ1) /= 0
   ELSE
     CALL abort('magconfig3D.f90, isInsideLCMS: called with wrong coordinates system.')
   ENDIF
!========================================================================
 END FUNCTION isInsideLCMS

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_SBgradSgradBi(XYZ, SS,Bi,gradS,gradBi)
!========================================================================
! Retuns S, B(:), gradS, grad(B(:)) in cartesian coordinates
! gradBi(9) = ( dBx/dx, dBy/dx, dBz/dx,
!               dBx/dy, dBy/dy, dBz/dy,
!               dBx/dz, dBy/dz, dBz/dz )
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in)  :: XYZ(3)
 REAL(8), INTENT(out) :: SS,Bi(1:3),gradS(1:3),gradBi(1:3,1:3)
 REAL(8) :: SS1,XYZ1(1:3),Bi1(1:3),gradS1(1:3),gradBi1(1:9)
!========================================================================
 XYZ1(:) = XYZ(:)
!----------
  IF (useMesh) THEN
   CALL MCM3DgetAllxyz(mConf,XYZ1, SS1,Bi1,gradS1,gradBi1)  ! use 3d-mesh
  ELSE
   CALL MCgetAllxyz(mConf,XYZ1, SS1,Bi1,gradS1,gradBi1) 
  ENDIF
!----------
 SS = SS1
 Bi(:) = Bi1(:)
 gradS(:) = gradS1(:)
!---
 gradBi(:,1) = (/gradBi1(1),gradBi1(1+3),gradBi1(1+6)/)   ! <- grad(Bx)
 gradBi(:,2) = (/gradBi1(2),gradBi1(2+3),gradBi1(2+6)/)   ! <- grad(By)
 gradBi(:,3) = (/gradBi1(3),gradBi1(3+3),gradBi1(3+6)/)   ! <- grad(Bz)
!========================================================================
 END SUBROUTINE MConf3D_SBgradSgradBi

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_gradS(XYZ) RESULT(grads)
!========================================================================
! function retuns grad(s) in cartesian coordinates
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: XYZ(3)
 REAL(8)  :: grads(3)
 REAL(8) :: XYZ1(3),grads1(3),Bi(3),b(3),SS,modB,parB
!========================================================================
 XYZ1(:) = XYZ(:)
!----------
  IF (useMesh) THEN
   CALL mcM3DgetGradSxyz(mconf,XYZ1,grads1)   ! use 3d-mesh
  ELSE
   CALL mcgetGradSxyz   (mconf,XYZ1,grads1)   ! use Newton's method
  ENDIF
  grads(:) = grads1(:)
!========================================================================
! correction to be sure that grad(SS) is perpendicular B
!========================================================================
 CALL MConf3D_FluxLabel_Bfield(XYZ, SS,Bi)
 modB = sqrt(max(comp_eps,dot_product(Bi,Bi)))
 b(:) = Bi(:)/modB
 parB = dot_product(grads,b)
!---
 grads(:) = grads(:)-parB*b(:)
!========================================================================
 END FUNCTION MConf3D_gradS

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_FluxLabel_Bfield(XYZ, S,B)
!========================================================================
! Function retuns the flux surface label S and the 
! magnetic field vector B for the given point XYZ
! Input:
!  XYZ  - array of size 3 containing cartesian coordinate of the point
! Output:
!  S - normalized toroidal flux, 
!       S > 1 if the point XYZ lies outside the last closed surface 
!  B - array with the components of magnetic field in Tesla, optional parameter
!========================================================================
 IMPLICIT NONE
 REAL(8),           INTENT(in)  :: XYZ(3)
 REAL(8),           INTENT(out) :: S
 REAL(8), OPTIONAL, INTENT(out) :: B(3)
 REAL(8)  :: RpZ(3),RpZ0(3),XYZ0(3),Bt(3),x,r,q,cs
 REAL(8) :: XYZ1(3),S1,B1(3)
 XYZ1(:) = XYZ(:)
!----------
 IF (useMesh) THEN
   CALL mcM3DgetSBxyz(mconf,XYZ1, S1,B1)   ! use 3d-mesh
 ELSE
   CALL mcgetSBxyz   (mconf,XYZ1, S1,B1)   ! use Newton's method
 ENDIF
!----------
 S = S1
 IF(present(B)) B(:) = B1(:)
!========================================================================
 END SUBROUTINE MConf3D_FluxLabel_Bfield

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_Cyl2Magn(RpZ, magCoord,B)
!========================================================================
! Transforms the cylindrical coordinates of the given point in the
! magnetic coordinates. Additionally, the magnetic field vector B(:) is 
! calculated.
! Input:
!  RpZ(1:3) - contains cylindr. coord. of the point, (R,phi,Z)
! Output:
!  magCoord(1:3) - contains magn. coord. of the point (s,the,phi)
!       S = magCoord(1) > 1 if the point lies outside of the volume
!  B(1:3) - cylindrical components of magn.field in Tesla, optional parameter
!========================================================================
 IMPLICIT NONE
 REAL(8),           INTENT(in)  :: RpZ(3)
 REAL(8),           INTENT(out) :: magCoord(3)
 REAL(8), OPTIONAL, INTENT(out) :: B(3)
 REAL(8) :: RpZ1(3),magCoord1(3),B1(3)
!========================================================================
! standard model (interpolation)
!========================================================================
 RpZ1(:) = RpZ(:)
 CALL mcCyl2Mag(mconf,RpZ1, magCoord1,B1)
!----------
 magCoord(:) = magCoord1(:)
 IF(present(B)) B(:) = B1(:)
!========================================================================
 END SUBROUTINE MConf3D_Cyl2Magn

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_Magn2Cyl(magCoord, RpZ)
!========================================================================
! Input:  magCoord(1:3) - Booser coord. of point, (s,the,phi)
! Output: RpZ(1:3)      - cylindr. coord. of point, (R,phi,Z)
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in)  :: magCoord(3)
 REAL(8), INTENT(out) :: RpZ(3)
 REAL(8) :: RpZ1(3),magCoord1(3)
!========================================================================
! standard model (interpolation)
!========================================================================
 magCoord1(:) = magCoord(:)
 CALL mcMag2Cyl(mconf,magCoord1, RpZ1)
!----------
 RpZ(:) = RpZ1(:)
!========================================================================
 END SUBROUTINE MConf3D_Magn2Cyl

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf3D_B_Jacobian(S,polAngle,torAngle,B,jac)
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S,polAngle,torAngle
 REAL(8), INTENT(out) :: B,jac
 REAL(8) :: S1,polAngle1,torAngle1,B1,jac1
!========================================================================
 S1        = S
 polAngle1 = polAngle
 torAngle1 = torAngle
 CALL mcBandJacobian(mconf,S1,polAngle1,torAngle1,B1,jac1)
 B   = B1
 jac = jac1
!========================================================================
 END SUBROUTINE MConf3D_B_Jacobian

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_Bmod(S,polAngle,torAngle) RESULT(Bmod)
!========================================================================
! Returns the absolute value of B
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S,polAngle,torAngle
 REAL(8)  :: Bmod
 REAL(8) :: S1,polAngle1,torAngle1
!========================================================================
 S1        = S
 polAngle1 = polAngle
 torAngle1 = torAngle
 Bmod = mcB(mconf,S1,polAngle1,torAngle1)
!========================================================================
 END FUNCTION MConf3D_Bmod

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_tor2pol(S_tor) RESULT(S_pol)
!========================================================================
! Transforms the toroidal flux in poloidal one
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S_tor
 REAL(8)  :: S_pol
 REAL(8) :: S1
!========================================================================
 S1    = S_tor
 S_pol = mcTorflux2Polflux(mconf,S1)
!========================================================================
 END FUNCTION MConf3D_tor2pol

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_dVdSt(S) RESULT(dVdSt)
!========================================================================
! Calculates |dV(s)/dStor| 
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S
 REAL(8)  :: dVdSt
 REAL(8) :: S1
!========================================================================
 S1    = S
 dVdSt = ABS(mcVprime(mconf,S1))
!========================================================================
 END FUNCTION MConf3D_dVdSt

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_iota(S) RESULT(iota)
!========================================================================
! Calculates the iota(s)
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S
 REAL(8)  :: iota
 REAL(8) :: S1
!========================================================================
 S1   = S
 iota = mcIota(mconf,S1)
!========================================================================
 END FUNCTION MConf3D_iota

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_iotaprime(S) RESULT(iotaprime)
!========================================================================
! Calculates the iota(s)
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S
 REAL(8)  :: iotaprime
 REAL(8) :: S1
!========================================================================
 S1   = S
 iotaprime = mcIotaPrime(mconf,S1)
!========================================================================
 END FUNCTION MConf3D_iotaprime

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_torflux(S) RESULT(psi)
!========================================================================
! Function returns toroidal flux for label S
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S
 REAL(8)  :: psi
 REAL(8) :: S1
!========================================================================
 S1  = S
 psi = mcTorFlux(mconf,S1)
!========================================================================
 END FUNCTION MConf3D_torflux

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_volume(S) RESULT(vol)
!========================================================================
! Function returns ABS(Volume) inside flux surface S
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: S
 REAL(8)  :: vol
 REAL(8) :: S1
!========================================================================
 S1  = S
 vol = ABS(mcVolume(mconf,S1))
!========================================================================
 END FUNCTION MConf3D_volume

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf3D_R0major() RESULT(R0)
!========================================================================
! Returns major radius
!========================================================================
 IMPLICIT NONE
 REAL(8) :: R0
!========================================================================
 R0 = mcR0(mconf)
!========================================================================
 END FUNCTION MConf3D_R0major

 FUNCTION MConf3D_Raxis(phi) RESULT(Raxi)
!========================================================================
! For the given toroidal angle calculates the radial position of axis.
!========================================================================
 IMPLICIT NONE
 REAL(8), INTENT(in) :: phi
 REAL(8)             :: Raxi
!========================================================================
 Raxi = mcRaxis(mconf,phi)
!========================================================================
 END FUNCTION MConf3D_Raxis

!########################################################################
!########################################################################
!########################################################################


!########################################################################
#endif
 END MODULE magconfig3D

!########################################################################
