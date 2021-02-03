module mod_ECRad_abs_Fa
! All credit for these routines goes to Daniela Farina and Lorenzo Figini
! Based on Farina, Daniela. "Relativistic dispersion relation of electron cyclotron waves." Fusion science and technology 53.1 (2008): 130-138.
#ifdef OMP
  use omp_lib
#endif
  implicit none


  integer, parameter :: r8=selected_real_kind(15,300)
  real(r8), parameter :: pi=3.1415926535897932384626433832795_r8
  real(r8), parameter :: pi_sqrt=1.7724538509055160272981674833411_r8
  complex(r8), parameter :: ui=(0._r8,1._r8)
  integer, parameter :: ntv=501
  real(r8), parameter :: tmax=5.0_r8, dt=2.0_r8*tmax/(ntv-1)
  real(r8), save :: ttv(ntv),extdtv(ntv)

contains

  subroutine set_extv
    implicit none
    integer :: i
    do i = 1, ntv
      ttv(i) = -tmax+(i-1)*dt
      extdtv(i) = exp(-ttv(i)*ttv(i))*dt
    end do
  end subroutine set_extv

  function expei(x)
! derived from calcei for usee with int=3 only
! result ->  y = expei(x)
!                               
! computes the exponential integral exp(-x)*ei(x)  for real arguments x
! where
!
!          /  integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
! ei(x) = {
!          \ -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
!
!  and where the first integral is a principal value integral.
!
    implicit none
    real(r8) :: expei
! arguments
    real(r8), intent(in) :: x
! local variables
    integer :: i, ip1, im1
    real(r8) :: frac, sump, sumq, t, w, xmx0, y, ysq
    real(r8), dimension(10) :: px, qx
    real(r8), parameter :: zero=0.0_r8,  p037=0.037_r8, half=0.5_r8,   &
        one=1.0_r8, two=2.0_r8, three=3.0_r8, four=4.0_r8,   six=6.0_r8,    &
        twelve=12.0_r8,         two4=24.0_r8, x01=381.5_r8,  x11=1024.0_r8, &
        x02=-5.1182968633365538008e-5_r8,     x0=0.37250741078136663466_r8
! machine-dependent constants
    real(r8), parameter :: xinf=1.79e+308_r8
! coefficients  for -1.0 <= x < 0.0
    real(r8), dimension(7), parameter ::                         &
       a=(/1.1669552669734461083368d2, 2.1500672908092918123209d3,    &
           1.5924175980637303639884d4, 8.9904972007457256553251d4,    &
           1.5026059476436982420737d5,-1.4815102102575750838086d5,    &
           5.0196785185439843791020d0/)
    real(r8), dimension(6), parameter ::                         &
       b=(/4.0205465640027706061433d1, 7.5043163907103936624165d2,    &
           8.1258035174768735759855d3, 5.2440529172056355429883d4,    &
           1.8434070063353677359298d5, 2.5666493484897117319268d5/)
! coefficients for -4.0 <= x < -1.0
    real(r8), dimension(9), parameter ::                         &
       c=(/3.828573121022477169108d-1, 1.107326627786831743809d+1,    &
           7.246689782858597021199d+1, 1.700632978311516129328d+2,    &
           1.698106763764238382705d+2, 7.633628843705946890896d+1,    &
           1.487967702840464066613d+1, 9.999989642347613068437d-1,    &
           1.737331760720576030932d-8/),                              &
       d=(/8.258160008564488034698d-2, 4.344836335509282083360d+0,    &
           4.662179610356861756812d+1, 1.775728186717289799677d+2,    &
           2.953136335677908517423d+2, 2.342573504717625153053d+2,    &
           9.021658450529372642314d+1, 1.587964570758947927903d+1,    &
           1.000000000000000000000d+0/)
! coefficients for x < -4.0
    real(r8), dimension(10), parameter ::                        &
       e=(/1.3276881505637444622987d+2,3.5846198743996904308695d+4,   &
           1.7283375773777593926828d+5,2.6181454937205639647381d+5,   &
           1.7503273087497081314708d+5,5.9346841538837119172356d+4,   &
           1.0816852399095915622498d+4,1.0611777263550331766871d03,   &
           5.2199632588522572481039d+1,9.9999999999999999087819d-1/), &
       f=(/3.9147856245556345627078d+4,2.5989762083608489777411d+5,   &
           5.5903756210022864003380d+5,5.4616842050691155735758d+5,   &
           2.7858134710520842139357d+5,7.9231787945279043698718d+4,   &
           1.2842808586627297365998d+4,1.1635769915320848035459d+3,   &
           5.4199632588522559414924d+1,1.0d0/)
! coefficients for rational approximation to ln(x/a), |1-x/a| < .1
    real(r8), dimension(4), parameter ::                         &
       plg=(/-2.4562334077563243311d+01,2.3642701335621505212d+02,    &
             -5.4989956895857911039d+02,3.5687548468071500413d+02/),  &
       qlg=(/-3.5553900764052419184d+01,1.9400230218539473193d+02,    &
             -3.3442903192607538956d+02,1.7843774234035750207d+02/)
! coefficients for  0.0 < x < 6.0,
! ratio of chebyshev polynomials
    real(r8), dimension(10), parameter ::                         &
       p=(/-1.2963702602474830028590d01,-1.2831220659262000678155d03,  &
           -1.4287072500197005777376d04,-1.4299841572091610380064d06,  &
           -3.1398660864247265862050d05,-3.5377809694431133484800d08,  &
            3.1984354235237738511048d08,-2.5301823984599019348858d10,  &
            1.2177698136199594677580d10,-2.0829040666802497120940d11/),&
       q=(/ 7.6886718750000000000000d01,-5.5648470543369082846819d03,  &
            1.9418469440759880361415d05,-4.2648434812177161405483d06,  &
            6.4698830956576428587653d07,-7.0108568774215954065376d08,  &
            5.4229617984472955011862d09,-2.8986272696554495342658d10,  &
            9.8900934262481749439886d10,-8.9673749185755048616855d10/)
! j-fraction coefficients for 6.0 <= x < 12.0
    real(r8), dimension(10), parameter ::                        &
       r=(/-2.645677793077147237806d00,-2.378372882815725244124d00,   &
           -2.421106956980653511550d01, 1.052976392459015155422d01,   &
            1.945603779539281810439d01,-3.015761863840593359165d01,   &
            1.120011024227297451523d01,-3.988850730390541057912d00,   &
            9.565134591978630774217d00, 9.981193787537396413219d-1/)
    real(r8), dimension(9), parameter ::                         &
       s=(/ 1.598517957704779356479d-4, 4.644185932583286942650d00,   &
            3.697412299772985940785d02,-8.791401054875438925029d00,   &
            7.608194509086645763123d02, 2.852397548119248700147d01,   &
            4.731097187816050252967d02,-2.369210235636181001661d02,   &
            1.249884822712447891440d00/)
! j-fraction coefficients for 12.0 <= x < 24.0
    real(r8), dimension(10), parameter ::                        &
       p1=(/-1.647721172463463140042d00,-1.860092121726437582253d01,  &
            -1.000641913989284829961d01,-2.105740799548040450394d01,  &
            -9.134835699998742552432d-1,-3.323612579343962284333d01,  &
             2.495487730402059440626d01, 2.652575818452799819855d01,  &
            -1.845086232391278674524d00, 9.999933106160568739091d-1/)
    real(r8), dimension(9), parameter ::                         &
       q1=(/ 9.792403599217290296840d01, 6.403800405352415551324d01,  &
             5.994932325667407355255d01, 2.538819315630708031713d02,  &
             4.429413178337928401161d01, 1.192832423968601006985d03,  &
             1.991004470817742470726d02,-1.093556195391091143924d01,  &
             1.001533852045342697818d00/)
! j-fraction coefficients for  x >= 24.0
    real(r8), dimension(10), parameter ::                        &
       p2=(/ 1.75338801265465972390d02,-2.23127670777632409550d02,    &
            -1.81949664929868906455d01,-2.79798528624305389340d01,    &
            -7.63147701620253630855d00,-1.52856623636929636839d01,    &
            -7.06810977895029358836d00,-5.00006640413131002475d00,    &
            -3.00000000320981265753d00, 1.00000000000000485503d00/)
    real(r8), dimension(9), parameter ::                         &
       q2=(/ 3.97845977167414720840d04, 3.97277109100414518365d00,    &
             1.37790390235747998793d02, 1.17179220502086455287d02,    &
             7.04831847180424675988d01,-1.20187763547154743238d01,    &
            -7.99243595776339741065d00,-2.99999894040324959612d00,    &
             1.99999999999048104167d00/)

    if (x == zero) then
      expei = -xinf
    else if (x < zero) then
!
! calculate ei for negative argument or for e1.
!
      y = abs(x)
      if (y <= one) then
        sump = a(7) * y + a(1)
        sumq = y + b(1)
        do i = 2, 6
          sump = sump * y + a(i)
          sumq = sumq * y + b(i)
        end do
        expei = (log(y) - sump / sumq ) * exp(y)
      else if (y <= four) then
        w = one / y
        sump = c(1)
        sumq = d(1)
        do i = 2, 9
          sump = sump * w + c(i)
          sumq = sumq * w + d(i)
        end do
        expei = - sump / sumq
      else
        w = one / y
        sump = e(1)
        sumq = f(1)
        do i = 2, 10
          sump = sump * w + e(i)
          sumq = sumq * w + f(i)
        end do
        expei = w * (w * sump / sumq - one)
      end if
    else if (x < six) then
!
! to improve conditioning, rational approximations are expressed
! in terms of chebyshev polynomials for 0 <= x < 6, and in
! continued fraction form for larger x.
!
      t = (x + x) / three - two
      px(1) = zero
      qx(1) = zero
      px(2) = p(1)
      qx(2) = q(1)
      do i = 2, 9
        ip1 = i+1
        im1 = i-1
        px(ip1) = t * px(i) - px(im1) + p(i)
        qx(ip1) = t * qx(i) - qx(im1) + q(i)
      end do
      sump = half * t * px(10) - px(9) + p(10)
      sumq = half * t * qx(10) - qx(9) + q(10)
      frac = sump / sumq
      xmx0 = (x - x01/x11) - x02
      if (abs(xmx0) >= p037) then
        expei = exp(-x) * ( log(x/x0) + xmx0 * frac )
      else
!
! special approximation to  ln(x/x0)  for x close to x0
!
        y = xmx0 / (x + x0)
        ysq = y*y
        sump = plg(1)
        sumq = ysq + qlg(1)
        do i = 2, 4
          sump = sump*ysq + plg(i)
          sumq = sumq*ysq + qlg(i)
        end do
        expei = exp(-x) * (sump / (sumq*(x+x0)) + frac) * xmx0
      end if
    else if (x < twelve) then
      frac = zero
      do i = 1, 9
        frac = s(i) / (r(i) + x + frac)
      end do
      expei = (r(10) + frac) / x
    else if (x <= two4) then
      frac = zero
      do i = 1, 9
        frac = q1(i) / (p1(i) + x + frac)
      end do
      expei = (p1(10) + frac) / x
    else
      y = one / x
      frac = zero
      do i = 1, 9
        frac = q2(i) / (p2(i) + x + frac)
      end do
      frac = p2(10) + frac
      expei = y + y * y * frac
    end if
  end function expei

  function fact(k)
! factorial function
    implicit none
    real(r8) :: fact
! arguments
    integer, intent(in) :: k
! local variables
    real(r8), dimension(0:15), parameter :: precomp = (/ &
        1._r8, 1._r8, 2._r8, 6._r8, 24._r8, 120._r8, 720._r8, 5040._r8, &
        40320._r8, 362880._r8, 3628800._r8, 39916800._r8, 479001600._r8, &
        6227020800._r8, 87178291200._r8, 1307674368000._r8 /)
    integer :: i

    if(k<0) then
      fact=0.0_r8
!    else
!      fact=1.0_r8
!      do i=2,k
!        fact=fact*i
!      end do
!    end if
    else if(k<16) then
      fact=precomp(k)
    else
      fact=precomp(15)
      do i=16,k
        fact=fact*i
      end do
    end if
  end function fact

  function gammln(x)
! Returns the value ln[Gamma(x)] for x > 0.
    implicit none
    real(r8) :: gammln
! arguments
    real(r8), intent(in) :: x
! parameters
    real(r8), parameter :: stp = 2.5066282746310005_r8
    real(r8), dimension(6), parameter :: cof = (/ 76.18009172947146_r8, &
        -86.50532032941677_r8, 24.01409824083091_r8, -1.231739572450155_r8,  &
        0.1208650973866179e-2_r8, -0.5395239384953e-5_r8 /)
! local variables
    integer :: j
    real(r8) :: ser, tmp, y

    y = x
    tmp = x + 5.5_r8
    tmp = (x + 0.5_r8)*log(tmp) - tmp
    ser = 1.000000000190015_r8
    do j = 1, 6
      y = y + 1._r8
      ser = ser + cof(j)/y
    end do
    gammln = tmp + log(stp*ser/x)
  end function gammln

  function ssbi(zz,n,l)
    implicit none
    real(r8), intent(in) :: zz
    integer, intent(in) :: n, l
    real(r8), dimension(l+2) :: ssbi
    integer :: k, m, mm
    real(r8) :: z2q, c0, c1, sbi
    real(r8), parameter :: eps = 1.e-10_r8

    z2q = 0.25_r8*zz**2
    do m = n, l+2
      c0 = 1.0_r8/exp(gammln(m + 1.5_r8))
!      c0 = 1.0_r8/gammah(m + 1.5_r8)
      sbi = c0
      do k = 1, 50
        c1 = c0*z2q/((m + k) + 0.5_r8)/k
        sbi = sbi + c1
        if(c1/sbi < eps) exit
        c0 = c1
      end do
      mm=m-n+1
      ssbi(mm)=sbi
    end do
  end function ssbi

!  function zetac(xi,yi)!,zr,zi,iflag)
!! PLASMA DISPERSION FUNCTION Z of complex argument
!! Z(z) = i sqrt(pi) w(z)
!! Function w(z) from:
!! algorithm 680, collected algorithms from acm.
!! this work published in transactions on mathematical software,
!! vol. 16, no. 1, pp. 47.
!!
!! given a complex number z = (xi,yi), this subroutine computes
!! the value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
!! where erfc is the complex complementary error-function and i
!! means sqrt(-1).
!! the accuracy of the algorithm for z in the 1st and 2nd quadrant
!! is 14 significant digits; in the 3rd and 4th it is 13 significant
!! digits outside a circular region with radius 0.126 around a zero
!! of the function.
!! all real variables in the program are double precision.
!!
!!
!! the code contains a few compiler-dependent parameters :
!! rmaxreal = the maximum value of rmaxreal equals the root of
!!            rmax = the largest number which can still be
!!            implemented on the computer in double precision
!!            floating-point arithmetic
!! rmaxexp  = ln(rmax) - ln(2)
!! rmaxgoni = the largest possible argument of a double precision
!!            goniometric function (dcos, dsin, ...)
!! the reason why these parameters are needed as they are defined will
!! be explained in the code by means of comments
!!
!!
!! parameter list
!! xi     = real      part of z
!! yi     = imaginary part of z
!! u      = real      part of w(z)
!! v      = imaginary part of w(z)
!! iflag  = an error flag indicating whether overflow will
!!          occur or not; type integer;
!!          the values of this variable have the following
!!          meaning :
!!          iflag=0 : no error condition
!!          iflag=1 : overflow will occur, the routine
!!                    becomes inactive
!! xi, yi       are the input-parameters
!! u, v, iflag  are the output-parameters
!!
!! furthermore the parameter factor equals 2/sqrt(pi)
!!
!! the routine is not underflow-protected but any variable can be
!! put to 0 upon underflow;
!!
!! reference - gpm poppe, cmj wijers; more efficient computation of
!! the complex error-function, acm trans. math. software.
!!
!  implicit none
!  complex(r8) :: zetac
!! arguments
!  real(r8), intent(in)  :: xi, yi
!!  real(r8), intent(out) :: zr, zi
!!  integer      , intent(out) :: iflag
!! local variables
!  integer       :: n, j, i, kapn, nu, np1
!  real(r8) :: u, v, x, y, xabs, yabs, qrho, xabsq, xquad, yquad,   &
!                   xsum, ysum, xaux, u1, v1, daux, u2, v2, h, h2,       &
!                   qlambda, rx, ry, sx, sy, tx, ty, c, w1
!! parameters
!  real(r8), parameter :: factor = 1.12837916709551257388_r8, &
!                              rpi    = 2.0_r8/factor!,             &
!!                              rmaxreal = 0.5e+154_r8,             &
!!                              rmaxexp = 708.503061461606_r8,      &
!!                              rmaxgoni = 3.53711887601422e+15_r8
!
!!    iflag = 0
!    xabs = abs(xi)
!    yabs = abs(yi)
!    x    = xabs/6.3_r8
!    y    = yabs/4.4_r8
!!
!! the following if-statement protects
!! qrho = (x**2 + y**2) against overflow
!!
!!    if ((xabs>rmaxreal).or.(yabs>rmaxreal)) then
!!      iflag = 1
!!      return
!!    end if
!    qrho = x**2 + y**2
!    xabsq = xabs**2
!    xquad = xabsq - yabs**2
!    yquad = 2*xabs*yabs
!    if (qrho<0.085264_r8) then
!!
!! if (qrho<0.085264_r8) then the faddeeva-function is evaluated
!! using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
!! n is the minimum number of terms needed to obtain the required
!! accuracy
!!
!      qrho  = (1-0.85_r8*y)*sqrt(qrho)
!      n     = nint(6 + 72*qrho)
!      j     = 2*n+1
!      xsum  = 1.0_r8/j
!      ysum  = 0.0_r8
!      do i=n, 1, -1
!        j    = j - 2
!        xaux = (xsum*xquad - ysum*yquad)/i
!        ysum = (xsum*yquad + ysum*xquad)/i
!        xsum = xaux + 1.0_r8/j
!      end do
!      u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0_r8
!      v1   =  factor*(xsum*xabs - ysum*yabs)
!      daux =  exp(-xquad)
!      u2   =  daux*cos(yquad)
!      v2   = -daux*sin(yquad)
!      u    = u1*u2 - v1*v2
!      v    = u1*v2 + v1*u2
!    else
!!
!! if (qrho>1.0) then w(z) is evaluated using the laplace
!! continued fraction
!! nu is the minimum number of terms needed to obtain the required
!! accuracy
!!
!! if ((qrho>0.085264d0).and.(qrho<1.0)) then w(z) is evaluated
!! by a truncated taylor expansion, where the laplace continued fraction
!! is used to calculate the derivatives of w(z)
!! kapn is the minimum number of terms in the taylor expansion needed
!! to obtain the required accuracy
!! nu is the minimum number of terms of the continued fraction needed
!! to calculate the derivatives with the required accuracy
!!
!      if (qrho>1.0_r8) then
!        h    = 0.0_r8
!        kapn = 0
!        qrho = sqrt(qrho)
!        nu   = int(3 + (1442/(26*qrho+77)))
!      else
!        qrho = (1-y)*sqrt(1-qrho)
!        h    = 1.88_r8*qrho
!        h2   = 2*h
!        kapn = nint(7  + 34*qrho)
!        nu   = nint(16 + 26*qrho)
!      endif
!      if (h>0.0_r8) qlambda = h2**kapn
!      rx = 0.0_r8
!      ry = 0.0_r8
!      sx = 0.0_r8
!      sy = 0.0_r8
!      do n=nu, 0, -1
!        np1 = n + 1
!        tx  = yabs + h + np1*rx
!        ty  = xabs - np1*ry
!        c   = 0.5_r8/(tx**2 + ty**2)
!        rx  = c*tx
!        ry  = c*ty
!        if ((h>0.0_r8).and.(n<=kapn)) then
!          tx = qlambda + sx
!          sx = rx*tx - ry*sy
!          sy = ry*tx + rx*sy
!          qlambda = qlambda/h2
!        endif
!      end do
!      if (h==0.0_r8) then
!        u = factor*rx
!        v = factor*ry
!      else
!        u = factor*sx
!        v = factor*sy
!      end if
!      if (yabs==0.0_r8) u = exp(-xabs**2)
!    end if
!!
!! evaluation of w(z) in the other quadrants
!!
!    if (yi<0.0_r8) then
!      if (qrho<0.085264_r8) then
!        u2    = 2*u2
!        v2    = 2*v2
!      else
!        xquad =  -xquad
!!
!! the following if-statement protects 2*exp(-z**2)
!! against overflow
!!
!!        if ((yquad>rmaxgoni).or.(xquad>rmaxexp)) then
!!          iflag=1
!!          return
!!        end if
!        w1 =  2.0_r8*exp(xquad)
!        u2  =  w1*cos(yquad)
!        v2  = -w1*sin(yquad)
!      end if
!      u = u2 - u
!      v = v2 - v
!      if (xi>0.0_r8) v = -v
!    else
!      if (xi<0.0_r8) v = -v
!    end if
!!    zr = -v*rpi
!!    zi =  u*rpi
!    zetac = cmplx(-v*rpi,u*rpi,kind=r8)
!  end function zetac
!
!  subroutine fsup(yg,anpl,amu,cefp,cefm,lrm)
!    implicit none
!    real(r8), intent(in) :: yg,anpl,amu
!    integer, intent(in) :: lrm
!    complex(r8), intent(out) :: cefp(0:lrm,0:2),cefm(0:lrm,0:2)
!    real(r8), parameter :: soglia = 0.7_r8
!    integer :: l,is,isa,ir
!    real(r8) :: anpl2hm1,alpha,phi2,phim,psi,apsi
!    real(r8) :: xp,yp,xm,ym,x0,y0
!    complex(r8) :: czp,czm,cf12,cf32,cphi,cz0,cdz0,cf0,cf1,cf2
!  
!    cefp=(0._r8,0._r8)
!    cefm=(0._r8,0._r8)
!    anpl2hm1=anpl**2/2.0_r8-1.0_r8
!    psi=sqrt(0.5_r8*amu)*anpl
!    apsi=abs(psi)
!    do is=-lrm,lrm
!      alpha=anpl2hm1+is*yg
!      phi2=amu*alpha
!      phim=sqrt(abs(phi2))
!      if (alpha>=0) then
!        xp=psi-phim
!        yp=0.0_r8
!        xm=-psi-phim
!        ym=0.0_r8
!        x0=-phim
!        y0=0.0_r8
!      else
!        xp=psi
!        yp=phim
!        xm=-psi
!        ym=phim
!        x0=0.0_r8
!        y0=phim
!      end if
!      czp=zetac(xp,yp)
!      czm=zetac(xm,ym)
!      if (alpha>0) then
!        cf12=-(czp+czm)/(2.0_r8*phim)
!      else if (alpha<0) then
!        cf12=-ui*(czp+czm)/(2.0_r8*phim)
!      else
!        cf12=(0.0_r8,0.0_r8)
!      end if
!      if(apsi>soglia) then
!        cf32=-(czp-czm)/(2.0_r8*psi)
!      else
!        cphi=phim
!        if(alpha<0) cphi=-ui*phim
!        cz0=zetac(x0,y0)
!        cdz0=2.0_r8*(1.0_r8-cphi*cz0)
!        cf32=cdz0
!      end if
!      cf0=cf12
!      cf1=cf32
!      if (is==0) then
!        cefp(0,0)=cf32
!        cefm(0,0)=cf32
!      end if
!      isa=abs(is)
!      do l=1,isa+2
!        if(apsi>soglia) then
!          cf2=(1.0_r8+phi2*cf0-(l-0.5_r8)*cf1)/psi**2
!        else
!          cf2=(1.0_r8+phi2*cf1)/(l+0.5_r8)
!        end if
!        ir=l-isa
!        if(ir>=0) then
!          cefp(isa,ir)=cefp(isa,ir)+cf2
!          if(is>0) then
!            cefm(isa,ir)=cefm(isa,ir)+cf2
!          else
!            cefm(isa,ir)=cefm(isa,ir)-cf2
!          end if
!        end if
!        cf0=cf1
!        cf1=cf2
!      end do
!    end do
!  end subroutine fsup
!  
!  subroutine dieltens_maxw_wr(xg,yg,anpl,amu,e330,epsl,lrm)
!! weakly relativistic dielectric tensor computation
!! Krivenski and Orefice, JPP 30,125 (1983)
!    real(r8), intent(in) :: xg,yg,anpl,amu
!    integer, intent(in) :: lrm
!    complex(r8), intent(out) :: e330, epsl(3,3,lrm)
!    integer :: l,lm,is,k
!    real(r8) :: anpl2,fcl,asl,bsl
!    complex(r8) :: ca11,ca12,ca13,ca22,ca23,ca33
!    complex(r8) :: cq0p,cq0m,cq1p,cq1m,cq2p
!    complex(r8) :: cefp(0:lrm,0:2),cefm(0:lrm,0:2)
!
!    anpl2=anpl**2
!    call fsup(yg,anpl,amu,cefp,cefm,lrm)
!    do l=1,lrm
!      lm=l-1
!      fcl=0.5_r8**l*((1.0_r8/yg)**2/amu)**lm*fact(2*l)/fact(l)
!      ca11=(0._r8,0._r8)
!      ca12=(0._r8,0._r8)
!      ca13=(0._r8,0._r8)
!      ca22=(0._r8,0._r8)
!      ca23=(0._r8,0._r8)
!      ca33=(0._r8,0._r8)
!      do is=0,l
!        k=l-is
!        asl=real((-1)**k,kind=r8)/(fact(is+l)*fact(l-is))
!        bsl=asl*(is**2+real(2*k*lm*(l+is),kind=r8)/(2*l-1))
!        cq0p=amu*cefp(is,0)
!        cq0m=amu*cefm(is,0)
!        cq1p=amu*anpl*(cefp(is,0)-cefp(is,1))
!        cq1m=amu*anpl*(cefm(is,0)-cefm(is,1))
!        cq2p=cefp(is,1)+amu*anpl2*(cefp(is,2)+cefp(is,0)-2.0_r8*cefp(is,1))
!        ca11=ca11+is**2*asl*cq0p
!        ca12=ca12+is*l*asl*cq0m
!        ca22=ca22+bsl*cq0p
!        ca13=ca13+is*asl*cq1m/yg
!        ca23=ca23+l*asl*cq1p/yg
!        ca33=ca33+asl*cq2p/yg**2
!      end do
!      epsl(1,1,l) =  - xg*ca11*fcl
!      epsl(1,2,l) =  + ui*xg*ca12*fcl
!      epsl(2,2,l) =  - xg*ca22*fcl
!      epsl(1,3,l) =  - xg*ca13*fcl
!      epsl(2,3,l) =  - ui*xg*ca23*fcl
!      epsl(3,3,l) =  - xg*ca33*fcl
!    end do
!    cq2p=cefp(0,1)+amu*anpl2*(cefp(0,2)+cefp(0,0)-2.0_r8*cefp(0,1))
!    e330=1.0_r8-xg*amu*cq2p
!    epsl(1,1,1) = 1._r8 + epsl(1,1,1)
!    epsl(2,2,1) = 1._r8 + epsl(2,2,1)
!    do l=1,lrm
!      epsl(2,1,l) = - epsl(1,2,l)
!      epsl(3,1,l) =   epsl(1,3,l)
!      epsl(3,2,l) = - epsl(2,3,l)
!    end do
!  end subroutine dieltens_maxw_wr
  
  subroutine hermitian(yg,anpl,amu,rr,lrm,iwarm)
    implicit none
    real(r8), intent(in) :: yg,amu,anpl
    integer, intent(in) :: lrm,iwarm
    real(r8), intent(out) :: rr(-lrm:lrm,0:2,0:lrm)
    integer :: n,k,m,llm,i,n1,nn
    real(r8) :: t,x
    real(r8) :: cmxw,cr
    real(r8) :: upl,upl2,gx,exdx,gr,s,ffe,zm,fe0m
    real(r8) :: sy1,sy2,sy3,sy4,bth2,bth4,bth6,bth8,anpl2,anpl4
    real(r8) :: rxt,bth,amu2,amu4,amu6,amu8,zm2,zm3

    rr(:,:,:) = 0._r8
    cmxw = 1._r8 + 15._r8/(8._r8*amu) + 105._r8/(128._r8*amu**2)
    cr = -amu*amu/(pi_sqrt*cmxw)
    llm=min(3,lrm)
    bth2=2.0_r8/amu
    bth=sqrt(bth2)
    amu2=amu*amu
    amu4=amu2*amu2
    amu6=amu4*amu2
    amu8=amu4*amu4
    if(iwarm>2) then
      n1=-llm
    else
      n1=1
    end if
    do i = 1, ntv
      t = ttv(i)
      rxt = sqrt(1.0_r8 + t*t/(2.0_r8*amu))
      x = t*rxt
      upl2 = bth2*x**2
      upl = bth*x
      gx = 1.0_r8 + t*t/amu
      exdx = cr*extdtv(i)*gx/rxt
      do n = n1, llm
        nn = abs(n)
        gr = anpl*upl + n*yg
        zm = -amu*(gx-gr)
        s = amu*(gx+gr)
        zm2 = zm**2
        zm3 = zm2*zm
        fe0m = expei(zm)
        do m = nn, llm
          if(m==0) then
            rr(0,2,0) = rr(0,2,0) - exdx*fe0m*upl2
          else 
            if (m==1) then
              ffe=(1._r8+s*(1._r8-zm*fe0m))/amu2
            else if (m==2) then
              ffe=(6._r8-2._r8*zm+4._r8*s+s*s*(1._r8+zm-zm2*fe0m))/amu4
            else if (m==3) then
              ffe=(18._r8*s*(s+4._r8-zm)+6._r8*(20._r8-8._r8*zm+zm2) &
                  +s**3*(2._r8+zm+zm2-zm3*fe0m))/amu6
            else
              ffe=(96._r8*s*(30._r8+s**2-10._r8*zm+zm**2) &
                  +24._r8*(210._r8-6._r8*s**2*(-5._r8+zm) &
                  -90._r8*zm+15._r8*zm**2-zm**3) &
                  +s**4*(6._r8+2._r8*zm+zm**2+zm**3-fe0m*zm**4))/amu8
            end if
            rr(n,0,m) = rr(n,0,m) + exdx*ffe
            rr(n,1,m) = rr(n,1,m) + exdx*ffe*upl
            rr(n,2,m) = rr(n,2,m) + exdx*ffe*upl2
          end if
        end do
      end do
    end do
    if(iwarm>2) return
    sy1 = 1._r8 + yg
    sy2 = 1._r8 + yg*2._r8
    sy3 = 1._r8 + yg*3._r8
    sy4 = 1._r8 + yg*4._r8
    anpl2=anpl**2
    anpl4=anpl2**2
    bth4=bth2*bth2
    bth6=bth4*bth2
    bth8=bth4*bth4
    rr(0,2,0) = -(1.0d0+bth2*(-1.25d0+1.5d0*anpl2) &
                +bth4*(1.71875d0-6.0d0*anpl2+3.75d0*anpl2*anpl2) &
                +bth6*3.0d0*(-65.0d0+456.0d0*anpl2-660.d0*anpl4 &
                + 280.0d0*anpl2*anpl4)/64.0d0 &
                + bth6*bth2*15.0d0*(252.853d3-2850.816d3*anpl2 &
                +6942.720d3*anpl4-6422.528d3*anpl4*anpl2 &
                +2064.384d3*anpl4*anpl4)/524.288d3)
    rr(0,1,1) = -anpl*bth2*(1.0d0+bth2*(-2.25d0+1.5d0*anpl2) &
                +bth4*9.375d-2*(6.1d1-9.6d1*anpl2+4.d1*anpl4 &
                +bth2*(-184.5d0+4.92d2*anpl2-4.5d2*anpl4 &
                +1.4d2*anpl2*anpl4)))
    rr(0,2,1) = -bth2*(1.0d0+bth2*(-0.5d0+1.5d0*anpl2) &
                +0.375d0*bth4*(3.d0-15.d0*anpl2+10.d0*anpl4) + &
                3.d0*bth6*(-61.d0+471.d0*anpl2-680*anpl4+280.d0*anpl2*anpl4)/64.d0)      
    rr(-1,0,1) = -2.0d0/sy1*(1.0d0+bth2/sy1*(-1.25d0+0.5d0*anpl2/sy1) &
                 +bth4/sy1*(-0.46875d0+(2.1875d0+0.625d0*anpl2)/sy1 &
                 -2.625d0*anpl2/sy1**2+0.75d0*anpl4/sy1**3) + bth6/sy1* &
                 (0.234375d0+(1.640625d0+0.234375d0*anpl2)/sy1 + &
                 (-4.921875d0-4.921875d0*anpl2)/sy1**2 + &
                 2.25d0*anpl2*(5.25d0 + anpl2)/sy1**3 - 8.4375d0*anpl4/sy1**4 + &  
                 1.875d0*anpl2*anpl4/sy1**5)+bth6*bth2/sy1*(0.019826889038*sy1 &
                 -0.06591796875d0+(-0.7177734375d0 - 0.1171875d0*anpl2)/sy1+ &
                 (-5.537109375d0 - 2.4609375d0*anpl2)/sy1**2 + &
                 (13.53515625d0 + 29.53125d0*anpl2 + 2.8125d0*anpl4)/sy1**3 + &
                 (-54.140625d0*anpl2 - 32.6953125d0*anpl4)/sy1**4 + &
                 (69.609375d0*anpl4 + 9.84375d0*anpl2*anpl4)/sy1**5 - &
                 36.09375d0*anpl2*anpl4/sy1**6 + 6.5625d0*anpl4**2/sy1**7))       
    rr(-1,1,1) = -anpl*bth2/sy1**2*(1.0d0+bth2*(1.25d0-3.5d0/sy1+ &
                 1.5d0*anpl2/sy1**2)+bth4*9.375d-2*((5.0d0-7.d1/sy1+ &
                 (126.d0+48.d0*anpl2)/sy1**2-144.d0*anpl2/sy1**3+4.d1*anpl4/sy1**4)+ &
                 bth2*(-2.5d0-3.5d1/sy1+(3.15d2+6.d1*anpl2)/sy1**2+ &
                 (-4.62d2-5.58d2*anpl2)/sy1**3+(9.9d2*anpl2+2.1d2*anpl4)/sy1**4 &
                 -6.6d2*anpl4/sy1**5+1.4d2*anpl4*anpl2/sy1**6)))   
    rr(-1,2,1) = -bth2/sy1*(1.0d0+bth2*(1.25d0-1.75d0/sy1+1.5d0*anpl2/sy1**2)+ &
                 bth4*3.d0/32.d0*(5.d0-35.d0/sy1+(42.d0+48.d0*anpl2)/sy1**2- &
                 108.d0*anpl2/sy1**3+40.0d0*anpl4/sy1**4 + &
                 0.5d0*bth2*(-5.d0-35.d0/sy1+(21.d1+12.d1*anpl2)/sy1**2- &
                 (231.d0+837.d0*anpl2)/sy1**3+12.d0*anpl2*(99.d0+35.d0*anpl2)/sy1**4 - &
                 1100.d0*anpl4/sy1**5 + 280.d0*anpl2*anpl4/sy1**6)))

    if(llm==1) return

    rr(0,0,2) = -4.0d0*bth2*(1.0d0+bth2*(-0.5d0+0.5d0*anpl2)+ &
                bth4*(1.125d0-1.875d0*anpl2+0.75d0*anpl4)+bth6* &
                3.0d0*(-61.d0+157.d0*anpl2-136.d0*anpl4+40.d0*anpl2*anpl4)/64.d0)   
    rr(0,1,2) = -2.0d0*anpl*bth4*(1.0d0+bth2*(-1.5d0+1.5d0*anpl2)+ &
                bth4*(39.d0-69.d0*anpl2+30.0d0*anpl4)/8.0d0)       
    rr(0,2,2) = -2.0d0*bth4*(1.0d0+bth2*(0.75d0+1.5d0*anpl2)+ bth4* &
                (13.d0-48.d0*anpl2 +40.0d0*anpl4)*3.0d0/32.0d0)
    rr(-1,0,2) = -4.0d0*bth2/sy1*(1.0d0+bth2* &
                 (1.25d0-1.75d0/sy1+0.5d0*anpl2/sy1**2)+bth4* &
                 (0.46875d0-3.28125d0/sy1+(3.9375d0+1.5d0*anpl2)/sy1**2 &
                 -3.375d0*anpl2/sy1**3+0.75d0*anpl4/sy1**4)+bth4*bth2*3.0d0/64.d0* &
                 (-5.0d0-35.0d0/sy1+(210.d0+40.d0*anpl2)/sy1**2-3.0d0*(77.d0+93.d0*anpl2)/sy1**3+ &
                 (396.0*anpl2+84.d0*anpl4)/sy1**4-220.d0*anpl4/sy1**5+40.d0*anpl4*anpl2/sy1**6))
    rr(-1,1,2) = -2.0d0*bth4*anpl/sy1**2*(1.0d0+bth2* &
                 (3.0d0-4.5d0/sy1+1.5d0*anpl2/sy1**2) + bth4* &
                 (20.d0-93.d0/sy1+(99.d0+42.d0*anpl2)/sy1**2- &
                 88.d0*anpl2/sy1**3+20.d0*anpl4/sy1**4)*3.0d0/16.0d0)          
    rr(-1,2,2) = -2.0d0*bth4/sy1*(1.0d0+bth2* &
                 (3.0d0-2.25d0/sy1+1.5d0*anpl2/sy1**2)+ bth4* &
                 (40.d0*anpl4/sy1**4-132.0d0*anpl2/sy1**3+ &
                 (66.d0+84.d0*anpl2)/sy1**2-93.d0/sy1+40.0d0)*3.0d0/32.0d0)
    rr(-2,0,2) = -4.0d0*bth2/sy2*(1.0d0+bth2* &
                 (1.25d0-1.75d0/sy2+0.5d0*anpl2/sy2**2)+bth4* &
                 (0.46875d0-3.28125d0/sy2+(3.9375d0+1.5d0*anpl2)/sy2**2 &
                 -3.375d0*anpl2/sy2**3+0.75d0*anpl4/sy2**4)+bth4*bth2*3.0d0/64.d0* &
                 (-5.0d0-35.0d0/sy2+(210.d0+40.d0*anpl2)/sy2**2-3.0d0*(77.d0+93.d0*anpl2)/sy2**3+ &
                 (396.0*anpl2+84.d0*anpl4)/sy2**4-220.d0*anpl4/sy2**5+40.d0*anpl4*anpl2/sy2**6))  
    rr(-2,1,2) = -2.0d0*bth4*anpl/sy2**2*(1.0d0+bth2* &
                 (3.0d0-4.5d0/sy2+1.5d0*anpl2/sy2**2) + bth4* &
                 (20.d0-93.d0/sy2+(99.d0+42.d0*anpl2)/sy2**2- &
                 88.d0*anpl2/sy2**3+20.d0*anpl4/sy2**4)*3.0d0/16.0d0)   
    rr(-2,2,2) = -2.0d0*bth4/sy2*(1.0d0+bth2* &
                 (3.0d0-2.25d0/sy2+1.5d0*anpl2/sy2**2) + bth4* &
                 (40.d0*anpl4/sy2**4-132.0d0*anpl2/sy2**3+ &
                 (66.d0+84.d0*anpl2)/sy2**2-93.d0/sy2+40.0d0)*3.0d0/32.0d0)   
    
    if(llm==2) return
    
    rr(0,0,3) = -12.0d0*bth4*(1.0d0+bth2*(0.75d0+0.5d0*anpl2)+bth4* &
                (1.21875d0-1.5d0*anpl2+0.75d0*anpl2*anpl2)) 
    rr(0,1,3) = -6.0d0*anpl*bth6*(1+bth2*(-0.25d0+1.5d0*anpl2))
    rr(0,2,3) = -6.0d0*bth6*(1.0d0+bth2*(2.5d0+1.5d0*anpl2))
    rr(-1,0,3) = -12.0d0*bth4/sy1* &
                 (1.0d0+bth2*(3.0d0-2.25d0/sy1+0.5d0*anpl2/sy1**2)+ &
                 bth4*(3.75d0-8.71875d0/sy1+(6.1875d0+2.625d0*anpl2) &
              /sy1**2-4.125d0*anpl2/sy1**3+0.75*anpl2*anpl2/sy1**4))
    rr(-1,1,3) = -6.0d0*anpl*bth6/sy1**2* &
                 (1.0d0+bth2*(5.25d0-5.5d0/sy1+1.5d0*anpl2/sy1**2))
    rr(-1,2,3) = -6.0d0*bth6/sy1* &
                 (1.0d0+bth2*(5.25d0-2.75d0/sy1+1.5d0*anpl2/sy1**2))
    rr(-2,0,3) = -12.0d0*bth4/sy2* &
                  (1.0d0+bth2*(3.0d0-2.25d0/sy2+0.5d0*anpl2/sy2**2)+ &
                  bth4*(3.75d0-8.71875d0/sy2+(6.1875d0+2.625d0*anpl2) &
              /sy2**2-4.125d0*anpl2/sy2**3+0.75*anpl2*anpl2/sy2**4))
    rr(-2,1,3) = -6.0d0*anpl*bth6/sy2**2* &
                  (1.0d0+bth2*(5.25d0-5.5d0/sy2+1.5d0*anpl2/sy2**2))
    rr(-2,2,3) = -6.0d0*bth6/sy2* &
                  (1.0d0+bth2*(5.25d0-2.75d0/sy2+1.5d0*anpl2/sy2**2))
    rr(-3,0,3) = -12.0d0*bth4/sy3* &
                  (1.0d0+bth2*(3.0d0-2.25d0/sy3+0.5d0*anpl2/sy3**2)+ &
                  bth4*(3.75d0-8.71875d0/sy3+(6.1875d0+2.625d0*anpl2) &
              /sy3**2-4.125d0*anpl2/sy3**3+0.75*anpl2*anpl2/sy3**4))
    rr(-3,1,3) = -6.0d0*anpl*bth6/sy3**2* &
                  (1.0d0+bth2*(5.25d0-5.5d0/sy3+1.5d0*anpl2/sy3**2))
    rr(-3,2,3) = -6.0d0*bth6/sy3* &
                  (1.0d0+bth2*(5.25d0-2.75d0/sy3+1.5d0*anpl2/sy3**2))

    if(llm==3) return

    rr(0,0,4) = -4.8d1*bth6*(1.0d0+bth2*(2.5d0+0.5d0*anpl2))
    rr(0,1,4) = -2.4d1*anpl*bth8
    rr(0,2,4) = -2.4d1*bth8
    rr(-1,0,4) = -4.8d1*bth6*(1.0d0+bth2*(5.25d0-2.75d0/sy1+0.5d0*anpl2/sy1**2))
    rr(-1,1,4) = -2.4d1*anpl*bth8/sy1**2
    rr(-1,2,4) = -2.4d1*bth8/sy1
    rr(-2,0,4) = -4.8d1*bth6*(1.0d0+bth2*(5.25d0-2.75d0/sy2+0.5d0*anpl2/sy2**2))
    rr(-2,1,4) = -2.4d1*anpl*bth8/sy2**2
    rr(-2,2,4) = -2.4d1*bth8/sy2
    rr(-3,0,4) = -4.8d1*bth6*(1.0d0+bth2*(5.25d0-2.75d0/sy3+0.5d0*anpl2/sy3**2))
    rr(-3,1,4) = -2.4d1*anpl*bth8/sy3**2
    rr(-3,2,4) = -2.4d1*bth8/sy3
    rr(-4,0,4) = -4.8d1*bth6*(1.0d0+bth2*(5.25d0-2.75d0/sy4+0.5d0*anpl2/sy4**2))
    rr(-4,1,4) = -2.4d1*anpl*bth8/sy4**2
    rr(-4,2,4) = -2.4d1*bth8/sy4
  end subroutine hermitian

  subroutine antihermitian(yg,anpl,amu,ri,lrm)
    implicit none
    real(r8), intent(in) :: yg,anpl,amu
    integer, intent(in) :: lrm
    real(r8), intent(out) :: ri(lrm,0:2,lrm)
    integer :: n, k, m, mm
    real(r8), dimension(lrm+2) :: fsbi
    real(r8) :: cmu,dnl,cmxw,ci,ygn,rdu2,rdu
    real(r8) :: du,ub,aa,up,um,gp,gm,xp,xm,eep,eem,ee,cm,cim
    real(r8) :: fi0p0,fi1p0,fi2p0,fi0m0,fi1m0,fi2m0
    real(r8) :: fi0p1,fi1p1,fi2p1,fi0m1,fi1m1,fi2m1,fi0m,fi1m,fi2m

    ri(:,:,:) = 0._r8
    dnl = 1.0_r8 - anpl**2
    cmu = anpl*amu
    cmxw = 1._r8 + 15._r8/(8._r8*amu) + 105._r8/(128._r8*amu**2)
    ci = sqrt(2._r8*pi*amu)*amu**2/cmxw
    do n=1,lrm
      ygn=n*yg
      rdu2=ygn**2-dnl
      if(rdu2>0.0_r8) then
        rdu=sqrt(rdu2)
        du=rdu/dnl
        ub=anpl*ygn/dnl
        aa=amu*anpl*du
        ! Change by S. Denk - it is possible to get values for gamma that are smaller than 1 - this line fixes that
        !if(ygn + anpl* (ub+du) < 1.d0 .or. (ygn + anpl* ub-du) < 1.d0) cycle! gamma < 1 possible for N_par large
        ! => resonance condition not full filled
        if (abs(aa)>5.0_r8) then !if (abs(aa)>1.0_r8) then
          up=ub+du
          um=ub-du
          gp=anpl*up+ygn
          gm=anpl*um+ygn
          xp=up+1.0_r8/cmu
          xm=um+1.0_r8/cmu
          eem=exp(-amu*(gm-1.0_r8))
          eep=exp(-amu*(gp-1.0_r8))
          fi0p0=-1.0_r8/cmu
          fi1p0=-xp/cmu
          fi2p0=-(1.0_r8/cmu**2+xp*xp)/cmu
          fi0m0=-1.0_r8/cmu
          fi1m0=-xm/cmu
          fi2m0=-(1.0_r8/cmu**2+xm*xm)/cmu
          do m=1,lrm
            fi0p1=-2.0_r8*m*(fi1p0-ub*fi0p0)/cmu
            fi0m1=-2.0_r8*m*(fi1m0-ub*fi0m0)/cmu
            fi1p1=-((1.0_r8+2*m)*fi2p0-2.0_r8*(m+1)*ub*fi1p0    &
                   +up*um*fi0p0)/cmu
            fi1m1=-((1.0_r8+2*m)*fi2m0-2.0_r8*(m+1)*ub*fi1m0    &
                   +up*um*fi0m0)/cmu
            fi2p1=(2.0_r8*(1+m)*fi1p1-2.0_r8*m*                 &
                 (ub*fi2p0-up*um*fi1p0))/cmu
            fi2m1=(2.0_r8*(1+m)*fi1m1-2.0_r8*m*                 &
                 (ub*fi2m0-up*um*fi1m0))/cmu
            if(m>=n) then
              ri(n,0,m)=0.5_r8*ci*dnl**m*(fi0p1*eep-fi0m1*eem)
              ri(n,1,m)=0.5_r8*ci*dnl**m*(fi1p1*eep-fi1m1*eem)
              ri(n,2,m)=0.5_r8*ci*dnl**m*(fi2p1*eep-fi2m1*eem)
            end if
            fi0p0=fi0p1
            fi1p0=fi1p1
            fi2p0=fi2p1
            fi0m0=fi0m1
            fi1m0=fi1m1
            fi2m0=fi2m1
          end do
        else
          ! Change by S. Denk - it is possible to get values for gamma that are smaller than 1 - this line fixes that
          !if((ygn- +anpl*ub) < 1.d0) cycle! gamma < 1 possible for N_par large
          ! => resonance condition not full filled
          ee=exp(-amu*(ygn-1.0_r8+anpl*ub))
          fsbi(1:lrm+2)=ssbi(aa,n,lrm)
          do m=n,lrm
            cm=pi_sqrt*fact(m)*du**(2*m+1)
            cim=0.5_r8*ci*dnl**m
            mm=m-n+1
            fi0m=cm*fsbi(mm)
            fi1m=-0.5_r8*aa*cm*fsbi(mm+1)
            fi2m=0.5_r8*cm*(fsbi(mm+1)+0.5_r8*aa*aa*fsbi(mm+2))
            ri(n,0,m)=cim*ee*fi0m
            ri(n,1,m)=cim*ee*(du*fi1m+ub*fi0m)
            ri(n,2,m)=cim*ee*(du*du*fi2m+2.0_r8*du*ub*fi1m+ub*ub*fi0m)
          end do
        end if
      end if
    end do
  end subroutine antihermitian

  subroutine dieltens_maxw_fr(xg,yg,anpl,amu,e330,epsl,lrm,iwarm)
! dielectric tensor elements using fully relativistic expressions
! Expansion to lrm-th order in Larmor radius
    implicit none
    real(r8), intent(in) :: xg,yg,anpl,amu
    integer, intent(in) :: lrm,iwarm
    complex(r8), intent(out) :: e330, epsl(3,3,lrm)
    integer :: i,j,l,lm,is,k
    complex(r8) :: ca11,ca12,ca22,ca13,ca23,ca33
    complex(r8) :: cq0p,cq0m,cq1p,cq1m,cq2p
    real(r8) :: rr(-lrm:lrm,0:2,0:lrm),ri(lrm,0:2,lrm)
    real(r8) :: fal,asl,bsl

    epsl(:,:,:) = (0._r8,0._r8)
    call hermitian(yg,anpl,amu,rr,lrm,iwarm)
    call antihermitian(yg,anpl,amu,ri,lrm)
    do l=1,lrm
      lm=l-1
      fal=-0.25_r8**l*fact(2*l)/(fact(l)**2*yg**(2*lm))
      ca11=(0._r8,0._r8)
      ca12=(0._r8,0._r8)
      ca13=(0._r8,0._r8)
      ca22=(0._r8,0._r8)
      ca23=(0._r8,0._r8)
      ca33=(0._r8,0._r8)
      do is=0,l
        k=l-is
        asl=real((-1)**k,kind=r8)/(fact(is+l)*fact(l-is))
        bsl=asl*(is*is+real(2*k*lm*(l+is),kind=r8)/(2*l-1))
        if(is>0) then
          cq0p=rr(is,0,l)+rr(-is,0,l)+ui*ri(is,0,l)
          cq0m=rr(is,0,l)-rr(-is,0,l)+ui*ri(is,0,l)
          cq1p=rr(is,1,l)+rr(-is,1,l)+ui*ri(is,1,l)
          cq1m=rr(is,1,l)-rr(-is,1,l)+ui*ri(is,1,l)
          cq2p=rr(is,2,l)+rr(-is,2,l)+ui*ri(is,2,l)
        else
          cq0p=rr(is,0,l)
          cq0m=rr(is,0,l)
          cq1p=rr(is,1,l)
          cq1m=rr(is,1,l)
          cq2p=rr(is,2,l)
        end if
        ca11=ca11+is**2*asl*cq0p
        ca12=ca12+is*l*asl*cq0m
        ca22=ca22+bsl*cq0p
        ca13=ca13+is*asl*cq1m/yg
        ca23=ca23+l*asl*cq1p/yg
        ca33=ca33+asl*cq2p/yg**2
      end do
      epsl(1,1,l) =  - xg*ca11*fal
      epsl(1,2,l) =  + ui*xg*ca12*fal
      epsl(2,2,l) =  - xg*ca22*fal
      epsl(1,3,l) =  - xg*ca13*fal
      epsl(2,3,l) =  - ui*xg*ca23*fal
      epsl(3,3,l) =  - xg*ca33*fal
    end do
    cq2p=rr(0,2,0)
    e330=1.0_r8+xg*cq2p
    epsl(1,1,1) = 1._r8 + epsl(1,1,1)
    epsl(2,2,1) = 1._r8 + epsl(2,2,1)
    do l=1,lrm
      epsl(2,1,l) = - epsl(1,2,l)
      epsl(3,1,l) =   epsl(1,3,l)
      epsl(3,2,l) = - epsl(2,3,l)
    end do
  end subroutine dieltens_maxw_fr

  subroutine warmdisp(xg,yg,anpl,amu,sox,iwarm,lrm,anprc,anpr,ex,ey,ez,ierr)
    implicit none
    real(r8), intent(in) :: xg,yg,anpl,amu,anprc
    integer, intent(in) :: sox,iwarm,lrm
    complex(r8), intent(out) :: anpr
    complex(r8), intent(out) :: ex,ey,ez
    integer, intent(out) :: ierr
    integer, parameter :: imx=100
    complex(r8) :: cc0,cc2,cc4,rr,anpr2a,anpra,anpr2,den
    complex(r8), dimension(3,3,lrm) :: epsl
    complex(r8), dimension(3,3) :: sepsl
    complex(r8) :: e330
    complex(r8) :: e11,e22,e12,e13,e23!,e33,e21,e31,e32
    complex(r8) :: a13,a31,a23,a32,a33
    real(r8) :: errnpr,s,anpl2,ex2,ey2,ez2,enx2
    integer :: i,j,k,ilrm

    ierr=0
    errnpr=1.0_r8
    anpr2a=anprc**2
    anpl2=anpl*anpl

!    if (iwarm==1) then
!      call dieltens_maxw_wr(xg,yg,anpl,amu,e330,epsl,lrm)
!    else
      call dieltens_maxw_fr(xg,yg,anpl,amu,e330,epsl,lrm,iwarm)
!    end if

    do i=1,imx
      do j=1,3
        do k=1,3
          sepsl(k,j)=cmplx(0.0_r8,0.0_r8,kind=r8)
          do ilrm=1,lrm
            sepsl(k,j)=sepsl(k,j)+epsl(k,j,ilrm)*anpr2a**(ilrm-1)
          end do
        end do
      end do
      anpra=sqrt(anpr2a)

      e11=sepsl(1,1)
      e22=sepsl(2,2)
      e12=sepsl(1,2)
      a33=sepsl(3,3)
      a13=sepsl(1,3)
      a23=sepsl(2,3)
      a31=a13
      a32=-a23
!     e33=e330+anpr2a*a33
      e13=anpra*a13
      e23=anpra*a23
!     e21=-e12
!     e31=e13
!     e32=-e23

!omaj - change: reduce the threshold (the first line is the original one)
!      if(i>2 .and. errnpr<1.0e-3_r8) exit
      if(i>2 .and. errnpr<1.0e-4_r8) exit
!omaj - end change

      cc4=(e11-anpl2)*(1.0_r8-a33)+(a13+anpl)*(a31+anpl)
      cc2=-e12*e12*(1.0_r8-a33) &
          -a32*e12*(a13+anpl)+a23*e12*(a31+anpl) &
          -(a23*a32+e330+(e22-anpl2)*(1.0_r8-a33))*(e11-anpl2) &
          -(a13+anpl)*(a31+anpl)*(e22-anpl2)
      cc0=e330*((e11-anpl2)*(e22-anpl2)+e12*e12)
      rr=cc2*cc2-4.0_r8*cc0*cc4
      if(yg>1.0_r8) then
        s=sox
        if(aimag(rr)<=0.0_r8) s=-s
      else
        s=-sox
        if(real(rr)<=0._r8 .and. aimag(rr)>=0._r8) s=-s
      end if
      anpr2=(-cc2+s*sqrt(rr))/(2.0_r8*cc4)
      errnpr=abs(1.0_r8-abs(anpr2)/abs(anpr2a))
      anpr2a=anpr2
    end do

    if(real(anpr2)<0.0_r8 .and. aimag(anpr2)<0.0_r8) then
!      print*,'  X, Y nperp2=',xg,yg,anpr2,'   nperp2 < 0'
      anpr2=0.0_r8
      ierr=99
    end if
    if(i>imx) then
!      print*,'    i>imx ',yg,errnpr,i
      ierr=100
    end if

    anpr=sqrt(anpr2)

    ex=(0.0_r8,0.0_r8)
    ey=(0.0_r8,0.0_r8)
    ez=(0.0_r8,0.0_r8)
    if (abs(anpl)>1.0e-6_r8) then
      den=e12*e23-(e13+anpr*anpl)*(e22-anpr2-anpl2)
      ey=-(e12*(e13+anpr*anpl)+(e11-anpl2)*e23)/den
      ez=(e12*e12+(e22-anpr2-anpl2)*(e11-anpl2))/den
      ez2=abs(ez)**2
      ey2=abs(ey)**2
      enx2=1.0_r8/(1.0_r8+ey2+ez2)
      ex=cmplx(sqrt(enx2),0.0_r8,kind=r8)
      ez2=ez2*enx2
      ey2=ey2*enx2
      ez=ez*ex
      ey=ey*ex
    else
      if(sox<0) then
        ez=(1.0_r8,0.0_r8)
        ez2=abs(ez)**2
      else
        ex2=1.0_r8/(1.0_r8+abs(-e11/e12)**2)
        ex=sqrt(ex2)
        ey=-ex*e11/e12
        ey2=abs(ey)**2
        ez2=0.0_r8
      end if
    end if
  end subroutine warmdisp

  subroutine colddisp(xg,yg,npl,sox,nprf)
!
! solution cold dispersion relation
!
    implicit none
!
    real(r8), intent(in) :: xg           ! X=omegap^2/omega^2
    real(r8), intent(in) :: yg           ! Y=omegac/omega
    real(r8), intent(in) :: npl          ! N parallel to B
    integer, intent(in) :: sox           ! sox=-1 ==> O mode, sox=1 ==> X mode
    real(r8), intent(out) :: nprf        ! N perpendicular to B (cold)
!
    real(r8) :: yg2,npl2,nprf2,n_f2,dnl,del  
!
    npl2=npl**2
    dnl=1.0_r8-npl2
    yg2=yg**2
    del=sqrt(dnl**2+4.0_r8*npl2*(1.0_r8-xg)/yg2)
    n_f2=1.0_r8-xg-xg*yg2*(1.0_r8+npl2+sox*del)/(1.0_r8-xg-yg2)/2.0_r8
!
    nprf2=n_f2-npl2                                                                                                        
    nprf=sqrt(nprf2)
!
  end subroutine colddisp

  function larmornumber(yg,npl,mu) result(nharm)
!
! computation of local harmonic number
!
    implicit none
!  
    integer :: nharm                    ! n=number of the armonic
    real(r8), parameter :: expcr=15._r8 ! maximum value for mu*(gamma-1) above
                                        ! which the distribution function is
                                        ! considered to be 0
!    real(r8), parameter :: eps=1.d-8    ! small number to have a correct rounding 
!                                        ! when ygnc/yg is an integer
    real(r8), intent(in) :: yg          ! Y=omegac/omega
    real(r8), intent(in) :: npl         ! parallel N
    real(r8), intent(in) :: mu          ! me*c^2/Te
    real(r8) :: ygn                     ! n*Y
    real(r8) :: gg                      ! gamma=sqrt(1+u^2)
    real(r8) :: dnl,rdu2
    real(r8) :: argexp
    integer :: lrm,larm, imax
!
    dnl=1.0_r8-npl**2
!      
    imax = 1
    nharm=int(1.0_r8/yg)
    if(nharm*yg < 1.0_r8) nharm=nharm+1
!    ygnc=sqrt(dnl)             ! critical value of Y*n to have resonance
!    nharm=int(ygnc/yg-eps)+1 ! first resonant harmonic
!    nharm=int(1.0_r8/yg-eps)+1 ! first harmonic with u1,u2 of opposite sign
!  
    do
      ygn=nharm*yg   
      rdu2=ygn**2-dnl
      gg=(ygn-sqrt(npl**2*rdu2))/dnl
      argexp=mu*(gg-1.0_r8)
      if(argexp.gt.expcr) exit
      nharm=nharm+1
      imax = imax + 1
      if(imax > 100) then
        nharm = int(yg)
        exit
      end if
    enddo
!    if((yg*(nharm-1))**2 >= dnl) nharm=nharm-1
!
  end function larmornumber
! Interface written by Omar Maj
  subroutine warmdamp(op, oc, Nr, theta, te, m, imod, N_perp_cmplx, pol_vec)

    ! ---------------------------------------------------------------------
    ! Wrapper for D. Farina' warmdisp routine.
    ! ---------------------------------------------------------------------

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
    real(r8), intent(in) ::  Nr! modulus of the refractive index N

    ! ... output variables ...
    complex(r8), intent(out) ::  &
        N_perp_cmplx                 ! complex refractive index
                                     ! N_perp_cmplx obtained from the disp. rel.
    complex(r8), intent(out), dimension(:), optional :: pol_vec
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
         ex, ey, ez             ! polarization unit vector

! ... f2py directives ...
!f2py integer intent(in) :: imod
!f2py real*8 intent(in) :: op, oc, Nr, theta, te
!f2py real*8 intent(out) :: imNprw

    ! =====================================================================
    ! Executable statements

    ! ... initialize common variables of the module ecdisp ...
    !call set_extv

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
         N_perp_cmplx, ex, ey, ez, err)

    ! ... extract imaginary part of the refractive index ...
    if(present(pol_vec)) then
      pol_vec(1) = ex
      pol_vec(2) = ey
      pol_vec(3) = ez
    end if
    ! ... do not allow negative values, put 0 instead ...
    if (aimag(N_perp_cmplx) < 0.) then
       N_perp_cmplx = cmplx(0., 0.)
       if(present(pol_vec)) pol_vec(:) = cmplx(0.d0, 0.d0)
    end if

    return
  end subroutine warmdamp


end module mod_ECRad_abs_Fa
