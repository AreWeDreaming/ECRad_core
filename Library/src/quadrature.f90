module quadrature
  public :: cdgqf, cgqf
  private :: gamma_func, class_matrix, imtqlx, parchk, scqf, sgqf

  contains

  subroutine cdgqf ( nt, rule, alpha, beta, t, wts )

  !*****************************************************************************80
  !
  !! CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
  !
  !  Discussion:
  !
  !    This routine computes all the knots and weights of a Gauss quadrature
  !    formula with a classical weight function with default values for A and B,
  !    and only simple knots.
  !
  !    There are no moments checks and no printing is done.
  !
  !    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 January 2010
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NT, the number of knots.
  !
  !    Input, integer ( kind = 4 ) RULE, the rule.
  !    1, Legendre,             (a,b)       1.0
  !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
  !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
  !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
  !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
  !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
  !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
  !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
  !
  !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
  !
  !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
  !
  !    Output, real ( kind = 8 ) T(NT), the knots.
  !
  !    Output, real ( kind = 8 ) WTS(NT), the weights.
  !
    implicit none

    integer ( kind = 4 ) nt

    real ( kind = 8 ) aj(nt)
    real ( kind = 8 ) alpha
    real ( kind = 8 ) beta
    real ( kind = 8 ) bj(nt)
    integer ( kind = 4 ) rule
    real ( kind = 8 ) t(nt)
    real ( kind = 8 ) wts(nt)
    real ( kind = 8 ) zemu

    call parchk ( rule, 2 * nt, alpha, beta )
  !
  !  Get the Jacobi matrix and zero-th moment.
  !
    call class_matrix ( rule, nt, alpha, beta, aj, bj, zemu )
  !
  !  Compute the knots and weights.
  !
    call sgqf ( nt, aj, bj, zemu, t, wts )

    return
  end subroutine cdgqf

  subroutine cgqf ( nt, rule, alpha, beta, a, b, t, wts )

  !*****************************************************************************80
  !
  !! CGQF computes knots and weights of a Gauss quadrature formula.
  !
  !  Discussion:
  !
  !    The user may specify the interval (A,B).
  !
  !    Only simple knots are produced.
  !
  !    Use routine EIQFS to evaluate this quadrature formula.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    16 February 2010
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NT, the number of knots.
  !
  !    Input, integer ( kind = 4 ) rule, the rule.
  !    1, Legendre,             (a,b)       1.0
  !    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
  !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
  !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
  !    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
  !    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
  !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
  !    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
  !    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
  !
  !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
  !
  !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
  !
  !    Input, real ( kind = 8 ) A, B, the interval endpoints, or
  !    other parameters.
  !
  !    Output, real ( kind = 8 ) T(NT), the knots.
  !
  !    Output, real ( kind = 8 ) WTS(NT), the weights.
  !
    implicit none

    integer ( kind = 4 ) nt

    real ( kind = 8 ) a
    real ( kind = 8 ) alpha
    real ( kind = 8 ) b
    real ( kind = 8 ) beta
    integer ( kind = 4 ) i
    integer ( kind = 4 ) rule
    integer ( kind = 4 ), allocatable :: mlt(:)
    integer ( kind = 4 ), allocatable :: ndx(:)
    real ( kind = 8 ) t(nt)
    real ( kind = 8 ) wts(nt)
  !
  !  Compute the Gauss quadrature formula for default values of A and B.
  !
    call cdgqf ( nt, rule, alpha, beta, t, wts )
  !
  !  Prepare to scale the quadrature formula to other weight function with
  !  valid A and B.
  !
    allocate ( mlt(1:nt) )

    mlt(1:nt) = 1

    allocate ( ndx(1:nt) )

    do i = 1, nt
      ndx(i) = i
    end do

    call scqf ( nt, t, mlt, wts, nt, ndx, wts, t, rule, alpha, beta, a, b )

    deallocate ( mlt )
    deallocate ( ndx )

    return
  end subroutine cgqf

  function gamma_func ( x )
  !*****************************************************************************80
  !
  !! GAMMA_FUNC evaluates Gamma(X) for a real argument.
  !
  !  Discussion:
  !
  !    This routine calculates the gamma function for a real argument X.
  !
  !    Computation is based on an algorithm outlined in reference 1.
  !    The program uses rational functions that approximate the gamma
  !    function to at least 20 significant decimal digits.  Coefficients
  !    for the approximation over the interval (1,2) are unpublished.
  !    Those for the approximation for 12 <= X are from reference 2.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    11 February 2008
  !
  !  Author:
  !
  !    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    William Cody,
  !    An Overview of Software Development for Special Functions,
  !    in Numerical Analysis Dundee, 1975,
  !    edited by GA Watson,
  !    Lecture Notes in Mathematics 506,
  !    Springer, 1976.
  !
  !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
  !    Charles Mesztenyi, John Rice, Henry Thatcher,
  !    Christoph Witzgall,
  !    Computer Approximations,
  !    Wiley, 1968,
  !    LC: QA297.C64.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) X, the argument of the function.
  !
  !    Output, real ( kind = 8 ) gamma_func, the value of the function.
  !
    implicit none
    real ( kind = 8 ), dimension ( 7 ) :: c = (/ &
     -1.910444077728D-03, &
      8.4171387781295D-04, &
     -5.952379913043012D-04, &
      7.93650793500350248D-04, &
     -2.777777777777681622553D-03, &
      8.333333333333333331554247D-02, &
      5.7083835261D-03 /)
    real ( kind = 8 ), parameter :: eps = 2.22D-16
    real ( kind = 8 ) fact
    integer ( kind = 4 ) i
    integer ( kind = 4 ) n
    real ( kind = 8 ), dimension ( 8 ) :: p = (/ &
      -1.71618513886549492533811D+00, &
       2.47656508055759199108314D+01, &
      -3.79804256470945635097577D+02, &
       6.29331155312818442661052D+02, &
       8.66966202790413211295064D+02, &
      -3.14512729688483675254357D+04, &
      -3.61444134186911729807069D+04, &
       6.64561438202405440627855D+04 /)
    logical parity
    real ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
    real ( kind = 8 ), dimension ( 8 ) :: q = (/ &
      -3.08402300119738975254353D+01, &
       3.15350626979604161529144D+02, &
      -1.01515636749021914166146D+03, &
      -3.10777167157231109440444D+03, &
       2.25381184209801510330112D+04, &
       4.75584627752788110767815D+03, &
      -1.34659959864969306392456D+05, &
      -1.15132259675553483497211D+05 /)
    real ( kind = 8 ) gamma_func
    real ( kind = 8 ) res
    real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
    real ( kind = 8 ) sum
    real ( kind = 8 ) x
    real ( kind = 8 ), parameter :: xbig = 171.624D+00
    real ( kind = 8 ) xden
    real ( kind = 8 ), parameter :: xinf = 1.0D+30
    real ( kind = 8 ), parameter :: xminin = 2.23D-308
    real ( kind = 8 ) xnum
    real ( kind = 8 ) y
    real ( kind = 8 ) y1
    real ( kind = 8 ) ysq
    real ( kind = 8 ) z

    parity = .false.
    fact = 1.0D+00
    n = 0
    y = x
  !
  !  Argument is negative.
  !
    if ( y <= 0.0D+00 ) then

      y = - x
      y1 = aint ( y )
      res = y - y1

      if ( res /= 0.0D+00 ) then

        if ( y1 /= aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
          parity = .true.
        end if

        fact = - pi / sin ( pi * res )
        y = y + 1.0D+00

      else

        res = xinf
        gamma_func = res
        return

      end if

    end if
  !
  !  Argument is positive.
  !
    if ( y < eps ) then
  !
  !  Argument < EPS.
  !
      if ( xminin <= y ) then
        res = 1.0D+00 / y
      else
        res = xinf
        gamma_func = res
        return
      end if

    else if ( y < 12.0D+00 ) then

      y1 = y
  !
  !  0.0 < argument < 1.0.
  !
      if ( y < 1.0D+00 ) then

        z = y
        y = y + 1.0D+00
  !
  !  1.0 < argument < 12.0.
  !  Reduce argument if necessary.
  !
      else

        n = int ( y ) - 1
        y = y - real ( n, kind = 8 )
        z = y - 1.0D+00

      end if
  !
  !  Evaluate approximation for 1.0 < argument < 2.0.
  !
      xnum = 0.0D+00
      xden = 1.0D+00
      do i = 1, 8
        xnum = ( xnum + p(i) ) * z
        xden = xden * z + q(i)
      end do

      res = xnum / xden + 1.0D+00
  !
  !  Adjust result for case  0.0 < argument < 1.0.
  !
      if ( y1 < y ) then

        res = res / y1
  !
  !  Adjust result for case 2.0 < argument < 12.0.
  !
      else if ( y < y1 ) then

        do i = 1, n
          res = res * y
          y = y + 1.0D+00
        end do

      end if

    else
  !
  !  Evaluate for 12.0 <= argument.
  !
      if ( y <= xbig ) then

        ysq = y * y
        sum = c(7)
        do i = 1, 6
          sum = sum / ysq + c(i)
        end do
        sum = sum / y - y + sqrtpi
        sum = sum + ( y - 0.5D+00 ) * log ( y )
        res = exp ( sum )

      else

        res = xinf
        gamma_func = res
        return

      end if

    end if
  !
  !  Final adjustments and return.
  !
    if ( parity ) then
      res = - res
    end if

    if ( fact /= 1.0D+00 ) then
      res = fact / res
    end if

    gamma_func = res

    return
  end function gamma_func

  subroutine class_matrix ( rule, m, alpha, beta, aj, bj, zemu )

  !*****************************************************************************80
  !
  !! CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
  !
  !  Discussion:
  !
  !    This routine computes the diagonal AJ and sub-diagonal BJ
  !    elements of the order M tridiagonal symmetric Jacobi matrix
  !    associated with the polynomials orthogonal with respect to
  !    the weight function specified by rule.
  !
  !    For weight functions 1-7, M elements are defined in BJ even
  !    though only M-1 are needed.  For weight function 8, BJ(M) is
  !    set to zero.
  !
  !    The zero-th moment of the weight function is returned in ZEMU.
  !
  !    Thanks to Larz White for pointing out that "class_matrix" had somehow
  !    been damaged to 'atrix', 02 February 2013.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 December 2009
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) rule, the rule.
  !    1, Legendre,             (a,b)       1.0
  !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
  !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
  !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
  !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
  !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
  !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
  !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
  !
  !    Input, integer ( kind = 4 ) M, the order of the Jacobi matrix.
  !
  !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
  !
  !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
  !
  !    Output, real ( kind = 8 ) AJ(M), BJ(M), the diagonal and subdiagonal
  !    of the Jacobi matrix.
  !
  !    Output, real ( kind = 8 ) ZEMU, the zero-th moment.
  !
    implicit none

    integer ( kind = 4 ) m

    real ( kind = 8 ) a2b2
    real ( kind = 8 ) ab
    real ( kind = 8 ) aba
    real ( kind = 8 ) abi
    real ( kind = 8 ) abj
    real ( kind = 8 ) abti
    real ( kind = 8 ) aj(m)
    real ( kind = 8 ) alpha
    real ( kind = 8 ) apone
    real ( kind = 8 ) beta
    real ( kind = 8 ) bj(m)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) rule
    real ( kind = 8 ), parameter :: pi = 3.14159265358979323846264338327950D+00
    real ( kind = 8 ) temp
    real ( kind = 8 ) temp2
    real ( kind = 8 ) zemu

    temp = epsilon ( temp )

    call parchk ( rule, 2 * m - 1, alpha, beta )

    temp2 = gamma_func ( 0.5D+00 )

    if ( 500.0D+00 * temp < abs ( temp2 * temp2 - pi ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLASS_MATRIX - Fatal error!'
      write ( *, '(a)' ) '  Gamma function does not match machine parameters.'
      stop
    end if

    if ( rule == 1 ) then

      ab = 0.0D+00

      zemu = 2.0D+00 / ( ab + 1.0D+00 )

      aj(1:m) = 0.0D+00

      do i = 1, m
        abi = i + ab * mod ( i, 2 )
        abj = 2 * i + ab
        bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 2 ) then

      zemu = pi

      aj(1:m) = 0.0D+00

      bj(1) =  sqrt ( 0.5D+00 )
      bj(2:m) = 0.5D+00

    else if ( rule == 3 ) then

      ab = alpha * 2.0D+00
      zemu = 2.0D+00**( ab + 1.0D+00 ) * ( gamma_func ( alpha + 1.0D+00 ) )**2 &
        / gamma_func ( ab + 2.0D+00 )

      aj(1:m) = 0.0D+00
      bj(1) = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
      do i = 2, m
        bj(i) = i * ( i + ab ) &
          / ( 4.0D+00 * ( i + alpha ) * ( i + alpha ) - 1.0D+00 )
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 4 ) then

      ab = alpha + beta
      abi = 2.0D+00 + ab
      zemu = 2.0D+00**( ab + 1.0D+00 ) * gamma_func ( alpha + 1.0D+00 ) &
        * gamma_func ( beta + 1.0D+00 ) / gamma_func ( abi )
      aj(1) = ( beta - alpha ) / abi
      bj(1) = 4.0D+00 * ( 1.0 + alpha ) * ( 1.0D+00 + beta ) &
        / ( ( abi + 1.0D+00 ) * abi * abi )
      a2b2 = beta * beta - alpha * alpha

      do i = 2, m
        abi = 2.0D+00 * i + ab
        aj(i) = a2b2 / ( ( abi - 2.0D+00 ) * abi )
        abi = abi * abi
        bj(i) = 4.0D+00 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) &
          / ( ( abi - 1.0D+00 ) * abi )
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 5 ) then

      zemu = gamma_func ( alpha + 1.0D+00 )

      do i = 1, m
        aj(i) = 2.0D+00 * i - 1.0D+00 + alpha
        bj(i) = i * ( i + alpha )
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 6 ) then

      zemu = gamma_func ( ( alpha + 1.0D+00 ) / 2.0D+00 )

      aj(1:m) = 0.0D+00

      do i = 1, m
        bj(i) = ( i + alpha * mod ( i, 2 ) ) / 2.0D+00
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 7 ) then

      ab = alpha
      zemu = 2.0D+00 / ( ab + 1.0D+00 )

      aj(1:m) = 0.0D+00

      do i = 1, m
        abi = i + ab * mod ( i, 2 )
        abj = 2 * i + ab
        bj(i) = abi * abi / ( abj * abj - 1.0D+00 )
      end do
      bj(1:m) =  sqrt ( bj(1:m) )

    else if ( rule == 8 ) then

      ab = alpha + beta
      zemu = gamma_func ( alpha + 1.0D+00 ) * gamma_func ( - ( ab + 1.0D+00 ) ) &
        / gamma_func ( - beta )
      apone = alpha + 1.0D+00
      aba = ab * apone
      aj(1) = - apone / ( ab + 2.0D+00 )
      bj(1) = - aj(1) * ( beta + 1.0D+00 ) / ( ab + 2.0D+00 ) / ( ab + 3.0D+00 )
      do i = 2, m
        abti = ab + 2.0D+00 * i
        aj(i) = aba + 2.0D+00 * ( ab + i ) * ( i - 1 )
        aj(i) = - aj(i) / abti / ( abti - 2.0D+00 )
      end do

      do i = 2, m - 1
        abti = ab + 2.0D+00 * i
        bj(i) = i * ( alpha + i ) / ( abti - 1.0D+00 ) * ( beta + i ) &
          / ( abti * abti ) * ( ab + i ) / ( abti + 1.0D+00 )
      end do

      bj(m) = 0.0D+00
      bj(1:m) =  sqrt ( bj(1:m) )

    end if

    return
  end subroutine class_matrix

  subroutine imtqlx ( n, d, e, z )

  !*****************************************************************************80
  !
  !! IMTQLX diagonalizes a symmetric tridiagonal matrix.
  !
  !  Discussion:
  !
  !    This routine is a slightly modified version of the EISPACK routine to
  !    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
  !
  !    The authors thank the authors of EISPACK for permission to use this
  !    routine.
  !
  !    It has been modified to produce the product Q' * Z, where Z is an input
  !    vector and Q is the orthogonal matrix diagonalizing the input matrix.
  !    The changes consist (essentially) of applying the orthogonal
  !    transformations directly to Z as they are generated.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 December 2009
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !    Roger Martin, James Wilkinson,
  !    The Implicit QL Algorithm,
  !    Numerische Mathematik,
  !    Volume 12, Number 5, December 1968, pages 377-383.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
  !    On output, the information in D has been overwritten.
  !
  !    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the
  !    matrix, in entries E(1) through E(N-1).  On output, the information in
  !    E has been overwritten.
  !
  !    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
  !    the value of Q' * Z, where Q is the matrix that diagonalizes the
  !    input symmetric tridiagonal matrix.
  !
    implicit none

    integer ( kind = 4 ) n

    real ( kind = 8 ) b
    real ( kind = 8 ) c
    real ( kind = 8 ) d(n)
    real ( kind = 8 ) e(n)
    real ( kind = 8 ) f
    real ( kind = 8 ) g
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ii
    integer ( kind = 4 ), parameter :: itn = 30
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real ( kind = 8 ) p
    real ( kind = 8 ) prec
    real ( kind = 8 ) r
    real ( kind = 8 ) s
    real ( kind = 8 ) z(n)

    prec = epsilon ( prec )

    if ( n == 1 ) then
      return
    end if

    e(n) = 0.0D+00

    do l = 1, n

      j = 0

      do

        do m = l, n

          if ( m == n ) then
            exit
          end if

          if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
            exit
          end if

        end do

        p = d(l)

        if ( m == l ) then
          exit
        end if

        if ( itn <= j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IMTQLX - Fatal error!'
          write ( *, '(a)' ) '  Iteration limit exceeded.'
          write ( *, '(a,i8)' ) '  J = ', j
          write ( *, '(a,i8)' ) '  L = ', l
          write ( *, '(a,i8)' ) '  M = ', m
          write ( *, '(a,i8)' ) '  N = ', n
          stop
        end if

        j = j + 1
        g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
        r =  sqrt ( g * g + 1.0D+00 )
        g = d(m) - p + e(l) / ( g + sign ( r, g ) )
        s = 1.0D+00
        c = 1.0D+00
        p = 0.0D+00
        mml = m - l

        do ii = 1, mml

          i = m - ii
          f = s * e(i)
          b = c * e(i)

          if ( abs ( g ) <= abs ( f ) ) then
            c = g / f
            r =  sqrt ( c * c + 1.0D+00 )
            e(i+1) = f * r
            s = 1.0D+00 / r
            c = c * s
          else
            s = f / g
            r =  sqrt ( s * s + 1.0D+00 )
            e(i+1) = g * r
            c = 1.0D+00 / r
            s = s * c
          end if

          g = d(i+1) - p
          r = ( d(i) - g ) * s + 2.0D+00 * c * b
          p = s * r
          d(i+1) = g + p
          g = c * r - b
          f = z(i+1)
          z(i+1) = s * z(i) + c * f
          z(i) = c * z(i) - s * f

        end do

        d(l) = d(l) - p
        e(l) = g
        e(m) = 0.0D+00

      end do

    end do
  !
  !  Sorting.
  !
    do ii = 2, n

      i = ii - 1
      k = i
      p = d(i)

      do j = ii, n
        if ( d(j) < p ) then
          k = j
          p = d(j)
        end if
      end do

      if ( k /= i ) then
        d(k) = d(i)
        d(i) = p
        p = z(i)
        z(i) = z(k)
        z(k) = p
      end if

    end do

    return
  end subroutine imtqlx

  subroutine parchk ( rule, m, alpha, beta )

  !*****************************************************************************80
  !
  !! PARCHK checks parameters ALPHA and BETA for classical weight functions.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 December 2009
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) rule, the rule.
  !    1, Legendre,             (a,b)       1.0
  !    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
  !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
  !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
  !    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
  !    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
  !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
  !    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
  !
  !    Input, integer ( kind = 4 ) M, the order of the highest moment to
  !    be calculated.  This value is only needed when KIND = 8.
  !
  !    Input, real ( kind = 8 ) ALPHA, BETA, the parameters, if required
  !    by the value of KIND.
  !
    implicit none

    real ( kind = 8 ) alpha
    real ( kind = 8 ) beta
    integer ( kind = 4 ) rule
    integer ( kind = 4 ) m
    real ( kind = 8 ) tmp

    if ( rule <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  rule <= 0.'
      stop
    end if
  !
  !  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
  !
    if ( 3 <= rule .and. alpha <= -1.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  3 <= rule and ALPHA <= -1.'
      stop
    end if
  !
  !  Check BETA for Jacobi.
  !
    if ( rule == 4 .and. beta <= -1.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARCHK - Fatal error!'
      write ( *, '(a)' ) '  rule == 4 and BETA <= -1.0.'
      stop
    end if
  !
  !  Check ALPHA and BETA for rational.
  !
    if ( rule == 8 ) then
      tmp = alpha + beta + m + 1.0D+00
      if ( 0.0D+00 <= tmp .or. tmp <= beta ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARCHK - Fatal error!'
        write ( *, '(a)' ) '  rule == 8 but condition on ALPHA and BETA fails.'
        stop
      end if
    end if

    return
  end subroutine parchk

  subroutine scqf ( nt, t, mlt, wts, nwts, ndx, swts, st, rule, alpha, beta, a, &
    b )

  !*****************************************************************************80
  !
  !! SCQF scales a quadrature formula to a nonstandard interval.
  !
  !  Discussion:
  !
  !    The arrays WTS and SWTS may coincide.
  !
  !    The arrays T and ST may coincide.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    27 December 2009
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NT, the number of knots.
  !
  !    Input, real ( kind = 8 ) T(NT), the original knots.
  !
  !    Input, integer ( kind = 4 ) MLT(NT), the multiplicity of the knots.
  !
  !    Input, real ( kind = 8 ) WTS(NWTS), the weights.
  !
  !    Input, integer ( kind = 4 ) NWTS, the number of weights.
  !
  !    Input, integer ( kind = 4 ) NDX(NT), used to index the array WTS.
  !    For more details see the comments in CAWIQ.
  !
  !    Output, real ( kind = 8 ) SWTS(NWTS), the scaled weights.
  !
  !    Output, real ( kind = 8 ) ST(NT), the scaled knots.
  !
  !    Input, integer ( kind = 4 ) rule, the rule.
  !    1, Legendre,             (a,b)       1.0
  !    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
  !    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
  !    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
  !    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
  !    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
  !    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
  !    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
  !    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
  !
  !    Input, real ( kind = 8 ) ALPHA, the value of Alpha, if needed.
  !
  !    Input, real ( kind = 8 ) BETA, the value of Beta, if needed.
  !
  !    Input, real ( kind = 8 ) A, B, the interval endpoints.
  !
    implicit none

    integer ( kind = 4 ) nt
    integer ( kind = 4 ) nwts

    real ( kind = 8 ) a
    real ( kind = 8 ) al
    real ( kind = 8 ) alpha
    real ( kind = 8 ) b
    real ( kind = 8 ) be
    real ( kind = 8 ) beta
    integer ( kind = 4 ) i
    integer ( kind = 4 ) k
    integer ( kind = 4 ) rule
    integer ( kind = 4 ) l
    integer ( kind = 4 ) mlt(nt)
    integer ( kind = 4 ) ndx(nt)
    real ( kind = 8 ) p
    real ( kind = 8 ) shft
    real ( kind = 8 ) slp
    real ( kind = 8 ) st(nt)
    real ( kind = 8 ) swts(nwts)
    real ( kind = 8 ) t(nt)
    real ( kind = 8 ) temp
    real ( kind = 8 ) tmp
    real ( kind = 8 ) wts(nwts)

    temp = epsilon ( temp )

    call parchk ( rule, 1, alpha, beta )

    if ( rule == 1 ) then

      al = 0.0D+00
      be = 0.0D+00

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    else if ( rule == 2 ) then

      al = -0.5D+00
      be = -0.5D+00

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    else if ( rule == 3 ) then

      al = alpha
      be = alpha

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    else if ( rule == 4 ) then

      al = alpha
      be = beta

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    else if ( rule == 5 ) then

      if ( b <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  B <= 0'
        stop
      end if

      shft = a
      slp = 1.0D+00 / b
      al = alpha
      be = 0.0D+00

    else if ( rule == 6 ) then

      if ( b <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  B <= 0.'
        stop
      end if

      shft = a
      slp = 1.0D+00 / sqrt ( b )
      al = alpha
      be = 0.0D+00

    else if ( rule == 7 ) then

      al = alpha
      be = 0.0D+00

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    else if ( rule == 8 ) then

      if ( a + b <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  A + B <= 0.'
        stop
      end if

      shft = a
      slp = a + b
      al = alpha
      be = beta

    else if ( rule == 9 ) then

      al = 0.5D+00
      be = 0.5D+00

      if ( abs ( b - a ) <= temp ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SCQF - Fatal error!'
        write ( *, '(a)' ) '  |B - A| too small.'
        stop
      end if

      shft = ( a + b ) / 2.0D+00
      slp = ( b - a ) / 2.0D+00

    end if

    p = slp**( al + be + 1.0D+00 )

    do k = 1, nt

      st(k) = shft + slp * t(k)
      l = abs ( ndx(k) )

      if ( l /= 0 ) then
        tmp = p
        do i = l, l + mlt(k) - 1
          swts(i) = wts(i) * tmp
          tmp = tmp * slp
        end do
      end if

    end do

    return
  end subroutine scqf

  subroutine sgqf ( nt, aj, bj, zemu, t, wts )

  !*****************************************************************************80
  !
  !! SGQF computes knots and weights of a Gauss Quadrature formula.
  !
  !  Discussion:
  !
  !    This routine computes all the knots and weights of a Gauss quadrature
  !    formula with simple knots from the Jacobi matrix and the zero-th
  !    moment of the weight function, using the Golub-Welsch technique.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 January 2010
  !
  !  Author:
  !
  !    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    Sylvan Elhay, Jaroslav Kautsky,
  !    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
  !    Interpolatory Quadrature,
  !    ACM Transactions on Mathematical Software,
  !    Volume 13, Number 4, December 1987, pages 399-415.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) NT, the number of knots.
  !
  !    Input, real ( kind = 8 ) AJ(NT), the diagonal of the Jacobi matrix.
  !
  !    Input/output, real ( kind = 8 ) BJ(NT), the subdiagonal of the Jacobi
  !    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
  !
  !    Input, real ( kind = 8 ) ZEMU, the zero-th moment of the weight function.
  !
  !    Output, real ( kind = 8 ) T(NT), the knots.
  !
  !    Output, real ( kind = 8 ) WTS(NT), the weights.
  !
    implicit none

    integer ( kind = 4 ) nt

    real ( kind = 8 ) aj(nt)
    real ( kind = 8 ) bj(nt)
    integer ( kind = 4 ) i
    real ( kind = 8 ) t(nt)
    real ( kind = 8 ) wts(nt)
    real ( kind = 8 ) zemu
  !
  !  Exit if the zero-th moment is not positive.
  !
    if ( zemu <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGQF - Fatal error!'
      write ( *, '(a)' ) '  ZEMU <= 0.'
      stop
    end if
  !
  !  Set up vectors for IMTQLX.
  !
    t(1:nt) = aj(1:nt)

    wts(1) = sqrt ( zemu )
    wts(2:nt) = 0.0D+00
  !
  !  Diagonalize the Jacobi matrix.
  !
    call imtqlx ( nt, t, bj, wts )

    wts(1:nt) = wts(1:nt)**2

    return
  end subroutine sgqf
end module quadrature
