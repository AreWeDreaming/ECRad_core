!########################################################################

 MODULE magconfig

!########################################################################
!-------
 IMPLICIT NONE
!----------------- PRIVATE variables ------------------------------------
 CHARACTER(Len=120),PRIVATE :: format_conf   ! format of the magnetic config.
 CHARACTER(Len=120),PRIVATE :: name_conf     ! filename of the magnetic config.
!
 LOGICAL,           PRIVATE :: useMesh,useSymm
!
 REAL(8),   PRIVATE :: hgrid, dR, dZ, dphi ! [m], [rad]
 REAL(8),   PRIVATE :: accbooz, tolharm
 real(8),  PRIVATE  :: comp_eps = 1.d-12
 real(8),  PRIVATE  :: B_ref
!
 INTEGER,    PRIVATE :: Nsym
 REAL(8), PRIVATE :: SSplus,SSmax
 LOGICAL,           PRIVATE :: trace_scrout=.true.

TYPE Scenario_type
 CHARACTER(Len=120) :: work_dir
 CHARACTER(Len=120) :: name_config
 CHARACTER(Len=120) :: format_config
 logical :: useMesh, useSymm

 real(8) :: B_ref, splus, smax, hgrid, dphic, accbooz, tolharm
 integer(8) :: mConfAddr = 0
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

   INTEGER(8) FUNCTION mcClone(mconf)
    INTEGER(8) :: mconf
   END FUNCTION mcClone
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
   REAL(8) FUNCTION mcRaxis(mconf,phi) 
     INTEGER(8) :: mconf
     REAL(8)    :: phi
   END FUNCTION mcRaxis
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
   REAL(8) FUNCTION mcGetB0(mconf,phi) ! rad
     INTEGER(8) :: mconf
     REAL(8)    :: phi
   END FUNCTION mcGetB0

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

 SUBROUTINE MConf_Setup_Config(Scen)
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
 B_ref       = Scen%B_ref
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
 END SUBROUTINE MConf_Setup_Config

!########################################################################
!########################################################################
!########################################################################

 subroutine MConf_load_magconf_file(loaded, mcAddr)
!========================================================================
  INTEGER(8), intent(out) :: mcAddr
  INTEGER, intent(out)    :: loaded
  CHARACTER(Len=4) :: fmt
!========================================================================
  fmt = format_conf(1:4)
  loaded = mcLoad(mcAddr,name_conf)
!========================================================================
 END subroutine MConf_load_magconf_file

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf_Load_MagConfig(loaded, mConfAddr)
!========================================================================
 IMPLICIT NONE
 integer, intent(out)      :: loaded
 INTEGER(8), INTENT(inout) :: mConfAddr  ! address of C3dMesh-object
 INTEGER :: i
!========================================================================
 IF(mConfAddr==0) THEN
   call MConf_load_magconf_file(loaded, mConfAddr)
 ENDIF
 CALL mcTruncate   (mConfAddr,tolharm)
 CALL mcSetAccuracy(mConfAddr,accbooz)
 CALL mcScaleB(mConfAddr,B_ref)
!---
 IF (useMesh) THEN
   IF(mcIsMeshOK(mConfAddr) /= 1) THEN ! test whether mConfAddr contains mesh
     IF (trace_scrout) WRITE(*,*) '   generating mesh...'
     CALL mcM3DsetSmax(mConfAddr,ssplus)
     IF(useSymm) THEN
       i = mcCreateMeshUsingSymmetry(mConfAddr,dR,dZ,dphi)
     ELSE
       i = mcCreateMesh             (mConfAddr,dR,dZ,dphi)
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
 END SUBROUTINE MConf_Load_MagConfig



!########################################################################

 FUNCTION MConf_isInsideLCMS(mConfAddr, coord, XYZ) RESULT(isInside)
!========================================================================
!  function returns 'true' if the point lies inside LCMS
!========================================================================
 IMPLICIT NONE
 Integer(8), Intent(in)        :: mConfAddr
 CHARACTER(Len=2), INTENT(in) :: coord
 REAL(8),        INTENT(in) :: XYZ(3)
 LOGICAL    :: isInside
!========================================================================
   isInside=.FALSE.
   IF (coord == 'cy') THEN
     isInside = MCCYLISINSIDE(mConfAddr,XYZ) /= 0
   ELSEIF (coord == 'ca') THEN
     isInside = MCXYZISINSIDE(mConfAddr,XYZ) /= 0
   ELSE
     CALL abort()
   ENDIF
!========================================================================
 END FUNCTION MConf_isInsideLCMS

!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf_SBgradSgradBi(mConfAddr, XYZ, SS,Bi,gradS,gradBi)
!========================================================================
! Retuns S, B(:), gradS, grad(B(:)) in cartesian coordinates
! gradBi(9) = ( dBx/dx, dBy/dx, dBz/dx,
!               dBx/dy, dBy/dy, dBz/dy,
!               dBx/dz, dBy/dz, dBz/dz )
!========================================================================
 IMPLICIT NONE
 Integer(8), Intent(in)        :: mConfAddr
 REAL(8), INTENT(in)  :: XYZ(3)
 REAL(8), INTENT(out) :: SS,Bi(1:3),gradS(1:3),gradBi(1:3,1:3)
 REAL(8) :: SS1,XYZ1(1:3),Bi1(1:3),gradS1(1:3),gradBi1(1:9)
!========================================================================
 XYZ1(:) = XYZ(:)
!----------
  IF (useMesh) THEN
   CALL MCM3DgetAllxyz(mConfAddr,XYZ1, SS1,Bi1,gradS1,gradBi1)  ! use 3d-mesh
  ELSE
   CALL MCgetAllxyz(mConfAddr,XYZ1, SS1,Bi1,gradS1,gradBi1)
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
 END SUBROUTINE MConf_SBgradSgradBi


!########################################################################
!########################################################################
!########################################################################

 SUBROUTINE MConf_FluxLabel_Bfield(mConfAddr, XYZ, S,B)
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
 Integer(8), Intent(in)          :: mConfAddr
 REAL(8),           INTENT(in)  :: XYZ(3)
 REAL(8),           INTENT(out) :: S
 REAL(8), OPTIONAL, INTENT(out) :: B(3)
 REAL(8)  :: RpZ(3),RpZ0(3),XYZ0(3),Bt(3),x,r,q,cs
 REAL(8) :: XYZ1(3),S1,B1(3)
 XYZ1(:) = XYZ(:)
!----------
 IF (useMesh) THEN
   CALL mcM3DgetSBxyz(mConfAddr,XYZ1, S1,B1)   ! use 3d-mesh
 ELSE
   CALL mcgetSBxyz   (mConfAddr,XYZ1, S1,B1)   ! use Newton's method
 ENDIF
!----------
 S = S1
 IF(present(B)) B(:) = B1(:)
!========================================================================
 END SUBROUTINE MConf_FluxLabel_Bfield

!########################################################################
!########################################################################
!########################################################################

 FUNCTION MConf_Raxis(mConfAddr, phi) RESULT(Raxi)
!========================================================================
! For the given toroidal angle calculates the radial position of axis.
!========================================================================
 IMPLICIT NONE
 Integer(8), Intent(in)        :: mConfAddr
 REAL(8), INTENT(in) :: phi
 REAL(8)             :: Raxi
!========================================================================
 Raxi = mcRaxis(mConfAddr,phi)
!========================================================================
 END FUNCTION MConf_Raxis

!########################################################################
!########################################################################
!########################################################################


!########################################################################
#endif
 END MODULE magconfig

!########################################################################
