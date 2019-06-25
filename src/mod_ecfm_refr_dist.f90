module mod_ecfm_radiation_dist
  use f90_kind
  implicit none
  logical                            :: fall_back_thermal
  !$OMP THREADPRIVATE(fall_back_thermal)
  public :: prepare_dist, &
            make_f_and_Rf_along_line, &
            radiation_dist_f_u, &
            radiation_dist_f_norm, &
            radiation_dist_Rf
  private :: fall_back_thermal, &
             make_norm_multi_slope, &
             radiation_dist_bi_maxJ, &
             radiation_dist_bi_maxw, &
             radiation_dist_Spitzer, &
             radiation_dist_bi_maxJ_Rf, &
             radiation_dist_bi_Maxw_Rf, &
             radiation_dist_Spitzer_Rf
  contains

subroutine prepare_dist(svec, Int_absz_many, Int_weights_many, f_spl, dist_params)
    use mod_ecfm_refr_types,            only: dstf, rad_diag_ch_mode_ray_freq_svec_type, min_level_log_ne, &
                                             bi_max, drift_m, Spitzer, multi_slope, runaway, ffp, fgene, Spitzer, &
                                             spl_type_2d, non_therm_params_type
    use constants,                      only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_fp_dist_utils,    only: make_B_min_and_f_inter
    use mod_ecfm_refr_gene_dist_utils,  only: make_g_inter
    use mod_ecfm_refr_interpol,         only: spline_1d
    implicit none
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)  :: svec
    real(rkind), dimension(:), intent(in) :: Int_absz_many, Int_weights_many
    type(spl_type_2d), intent(inout) :: f_spl
    type(non_therm_params_type), intent(out)               :: dist_params
    real(rkind)                                       :: mu
    fall_back_thermal = .true.
    if(dstf == "Th") then
      print*, "Called prepare_dist for a thermal distribution function"
      call abort()
    end if
    if (dstf == "numeric") then
      if(svec%rhop >= ffp%rhop_max) return
    else if (trim(dstf) == "gene" .or. trim(dstf) == "gcomp" ) then
      if(svec%rhop >= fgene%rhop_max .or. svec%rhop <= fgene%rhop_min) return
    else if (dstf == "Bi_MaxJ" .or. dstf == "Bi_Maxw") then
      if(svec%rhop >= bi_max%rhop_max) return
    else if (dstf == "drift_m") then
      if(svec%rhop >= drift_m%rhop_max) return
    else if(dstf == "multi_s") then
      if(svec%rhop >= multi_slope%rhop_max) return
    else if(dstf == "runaway") then
      if(svec%rhop >= runaway%rhop_max) return
    else if(dstf == "Spitzer") then
      if(svec%rhop > Spitzer%rhop_max) return
    end if
    fall_back_thermal = .false.
    mu = (mass_e * c0**2.d0) / (e0 * svec%Te)
    dist_params%B_min = 0.d0
    if(dstf == "numeric") then
      call make_B_min_and_f_inter(svec, f_spl, dist_params%B_min)
      return
    else if(trim(dstf) == "gene" .or. trim(dstf) == "gcomp") then
      call make_g_inter(svec, f_spl)
      if(trim(dstf) == "gcomp") then
#ifdef NAG
        call spline_1d(fgene%Te_perp_spl, svec%rhop, dist_params%Te_perp, nag_spline = fgene%Te_perp_nag_spl)
        call spline_1d(fgene%Te_par_spl, svec%rhop, dist_params%Te_par, nag_spline = fgene%Te_perp_nag_spl)
#else
        call spline_1d(fgene%Te_perp_spl, svec%rhop, dist_params%Te_perp)
        call spline_1d(fgene%Te_par_spl, svec%rhop, dist_params%Te_par)
#endif
      end if
      return
    else if(dstf == "Bi_Maxw" .or. dstf == "Bi_MaxJ") then
      if(svec%rhop < Bi_max%rhop_max) then
#ifdef NAG
        call spline_1d(bi_max%ne_spl, svec%rhop, dist_params%ne, nag_spline = bi_max%ne_nag_spl)
#else
        call spline_1d(bi_max%ne_spl, svec%rhop, dist_params%ne)
#endif
      dist_params%ne = exp(dist_params%ne)
      if(dist_params%ne < min_level_log_ne) dist_params%ne = 0.d0
      else
        dist_params%ne = 0.d0
      end if
      dist_params%Te_perp = bi_max%Te_perp
      dist_params%Te_par = bi_max%Te_par
    else if(dstf == "drift_m") then
      if(svec%rhop < drift_m%rhop_max) then
#ifdef NAG
        call spline_1d(drift_m%ne_spl, svec%rhop, dist_params%ne, nag_spline = drift_m%ne_nag_spl)
#else
        call spline_1d(drift_m%ne_spl, svec%rhop, dist_params%ne)
#endif
        dist_params%ne = exp(dist_params%ne)
        if(dist_params%ne < min_level_log_ne) dist_params%ne = 0.d0
      else
        dist_params%ne = 0.d0
      end if
      dist_params%Te_perp = drift_m%Te_perp
      dist_params%Te_par = drift_m%Te_par
      dist_params%u_perp_drift = drift_m%u_perp_drift
      dist_params%u_par_drift = drift_m%u_par_drift
    else if(dstf == "runaway") then
      if(svec%rhop < runaway%rhop_max) then
#ifdef NAG
        call spline_1d(runaway%rel_ne_spl, svec%rhop, dist_params%cur_rel_ne, nag_spline = runaway%rel_ne_nag_spl)
#else
        call spline_1d(runaway%rel_ne_spl, svec%rhop, dist_params%cur_rel_ne)
#endif
        dist_params%cur_rel_ne = exp(dist_params%cur_rel_ne)
        if(dist_params%cur_rel_ne < min_level_log_ne) dist_params%cur_rel_ne = 0.d0
      else
        dist_params%cur_rel_ne = 0.d0
      end if
    else if(dstf == "Spitzer") then
      if(svec%rhop < Spitzer%rhop_max)  then
#ifdef NAG
        call spline_1d(Spitzer%j_spl, svec%rhop, dist_params%vd, nag_spline = Spitzer%j_nag_spl)
#else
        call spline_1d(Spitzer%j_spl, svec%rhop, dist_params%vd)
#endif
      else
        dist_params%vd = 0.d0
      end if
      dist_params%vd = dist_params%vd  /(svec%ne* e0)
    else if(dstf == "multi_s") then
      if(svec%rhop > multi_slope%rhop_max) then
        dist_params%mu_slope = mu
        dist_params%norm = 1.d0
        dist_params%total_norm = 1.d0 ! Normalized with thermal distribution norm later
      else
#ifdef NAG
        call spline_1d(multi_slope%Te_slope_spl, svec%rhop, dist_params%mu_slope, nag_spline = multi_slope%Te_slope_nag_spl)
#else
        call spline_1d(multi_slope%Te_slope_spl, svec%rhop, dist_params%mu_slope)
#endif
        dist_params%mu_slope = (mass_e * c0**2.d0) / ( e0 * dist_params%mu_slope)
        ! Logarithm to avoid underflow
        dist_params%norm = mu * (1.d0 - multi_slope%gamma_switch) - &
                           dist_params%mu_slope * (1.d0 - multi_slope%gamma_switch)
        dist_params%norm = exp(dist_params%norm)
        call make_norm_multi_slope(svec,mu, Int_absz_many, Int_weights_many, dist_params)
      end if
    end if
end subroutine prepare_dist

subroutine make_norm_multi_slope(svec, mu, Int_absz_many, Int_weights_many, dist_params) !
    use mod_ecfm_refr_types,        only: rad_diag_ch_mode_ray_freq_svec_type, &
                                          non_therm_params_type
    use constants,                  only: pi
    implicit none

    type(rad_diag_ch_mode_ray_freq_svec_type),    intent(in)  :: svec
    type(non_therm_params_type), intent(inout)                :: dist_params
    real(rkind),    intent(in)  :: mu
    ! for numerial calculation of integral
    real(rkind), dimension(:), intent(in) :: Int_absz_many, Int_weights_many
    real(rkind), dimension(size(Int_weights_many))     :: u_par,u_perp, gamma, int_u_add
    real(rkind)     :: a, b, du_dalpha,a_2, b_2, du_dalpha_2, int_u_perp, gamma_max, u_max, Integral
    integer(ikind)  :: k,l

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! numerial calculation of the norm of the two slope distribution using gaussian quadrature with boundaries -u_max, u_max !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dist_params%total_norm = radiation_dist_f_norm(svec%Te, "relamax") ! Causes the integral to be already almost normalized
    int_u_perp = 0.0
    int_u_add(:) = 0
    gamma_max = (10.0 / max(mu, dist_params%mu_slope)) + 1.0
    u_max = sqrt(gamma_max**2 - 1.0)
    a =  -u_max
    b =  u_max
    du_dalpha = (b - a)/2.d0
    a_2 =  0
    b_2 =  u_max
    du_dalpha_2 = (b_2 - a_2)/2.d0
    Integral = 0.d0
    do k = 1, size(Int_weights_many)
      u_par(:) = Int_absz_many(k) * du_dalpha + (b + a) / 2.d0
      u_perp(:) = Int_absz_many(:) * du_dalpha_2 + (b_2 + a_2) / 2.d0
      gamma(:) = sqrt(1 + u_par(:)**2 + u_perp(:)**2)
      int_u_add(:) = radiation_dist_f_u(u_par, u_perp(:), gamma(:), mu, svec, dist_params) * 2 * Pi &
          * u_perp(:) * Int_weights_many(:)* du_dalpha_2
      do l = 1 ,size(Int_weights_many)
        int_u_perp = int_u_perp + int_u_add(l)
      end do
      Integral  = Integral + int_u_perp * Int_weights_many(k) * du_dalpha
      int_u_perp = 0.0
    enddo
    dist_params%total_norm  = 1.d0 / Integral
    !
    if(dist_params%total_norm /= dist_params%total_norm) then
      print*,dist_params%total_norm
      stop "nan in norm"
    else if (dist_params%total_norm < 0.0) then
      print*,dist_params%total_norm
      stop "- in norm"
    end if
  end subroutine make_norm_multi_slope

  subroutine make_f_and_Rf_along_line(u_par, u_perp, gamma, m_omega_bar, N_par, mu, svec, f_spl, dist_params, dstf, f, Rf)
    use constants,                   only: pi, e0,c0, mass_e
    use mod_ecfm_refr_types,         only: spl_type_2d, rad_diag_ch_mode_ray_freq_svec_type, non_therm_params_type, multi_slope
    use mod_ecfm_refr_fp_dist_utils, only: make_f_and_f_grad_along_line
    use mod_ecfm_refr_gene_dist_utils, only: make_gene_f_and_gene_f_grad_along_line, &
                                             make_gene_f0_and_gene_f0_grad_along_line
    implicit none
    real(rkind), dimension(:),  intent(in)            :: u_par, u_perp, gamma
    real(rkind),                intent(in)            :: m_omega_bar, N_par, mu
    type(spl_type_2d), intent(in)                     :: f_spl
    type(non_therm_params_type), intent(in)           :: dist_params
    character(*),  intent(in)                         :: dstf
    real(rkind), dimension(:),  intent(out)           :: f, Rf
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in) :: svec
    real(rkind), dimension(size(u_par)) :: df_du_par, df_du_perp
    real(rkind)                         :: a
    if(dstf == "maxwell") then
      f =  Exp(-mu/2.0 * (u_par**2 + u_perp**2))
      Rf = - mu * f
    else if(trim(dstf) == "Th" .or. fall_back_thermal) then
      f =  exp( mu* (1 - gamma))
      Rf = - mu * f
    else if (dstf == "numeric") then
      call make_f_and_f_grad_along_line(u_par, u_perp, svec, f_spl, dist_params%B_min, f, df_du_par, df_du_perp)!
      Rf = (m_omega_bar * df_du_perp / u_perp  + &
        N_par *  df_du_par)
    else if( trim(dstf) == "gene" .or. trim(dstf) == "genef0" .or. trim(dstf) == "gcomp") then
      ! With GENE background
      if(trim(dstf) == "gene".or. trim(dstf) == "gcomp") then
        call make_gene_f_and_gene_f_grad_along_line(u_par, u_perp, svec, f_spl, f, df_du_par, df_du_perp)
      else
        call make_gene_f0_and_gene_f0_grad_along_line(u_par, u_perp, svec, f, df_du_par, df_du_perp)
      end if
      Rf = (m_omega_bar * df_du_perp / u_perp + &
            N_par *  df_du_par)
    else if( trim(dstf) == "gcomp0") then
      f = radiation_dist_bi_maxJ(u_par, u_perp, gamma, dist_params%Te_perp, dist_params%Te_par)
      Rf = radiation_dist_bi_maxJ_Rf(u_par, u_perp, gamma, m_omega_bar, N_par, dist_params%Te_par, dist_params%Te_perp)
    else if( trim(dstf) == "Bi_MaxJ") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      f = (1.0 - dist_params%ne) * a * sqrt(mu / (2.d0 * pi))**3 * &
                                        exp( mu* (1.d0 - gamma)) ! Thermal bulk
      f = radiation_dist_bi_maxJ(u_par, u_perp, gamma, dist_params%Te_perp, dist_params%Te_par)
      Rf = -(1.0 - dist_params%ne) * a * (sqrt(mu / (2.d0 * pi))**3) * &
                                         exp( mu* (1.d0 - gamma)) * mu
      Rf = Rf + radiation_dist_bi_maxJ_Rf(u_par, u_perp, gamma, m_omega_bar, &
                                          N_par, dist_params%Te_perp, dist_params%Te_par)
    else if( trim(dstf) == "Bi_Maxw") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      f = (1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
                                       exp( mu* (1 - gamma)) ! Thermal bulk
      f = f + dist_params%ne * radiation_dist_bi_maxw(u_par, u_perp, gamma, dist_params%Te_par, dist_params%Te_perp) ! BiMaxwellian Tail
      Rf = -(1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
                                         exp( mu* (1 - gamma)) * mu
      Rf = Rf + dist_params%ne * radiation_dist_bi_maxw_Rf(u_par, u_perp,gamma,m_omega_bar, &
                                                           N_par, dist_params%Te_par, &
                                                           dist_params%Te_perp)
    else if(trim(dstf) == "multi_s") then
      f = exp( mu* (1.d0 - gamma))
      where(gamma > multi_slope%gamma_switch) f = exp(dist_params%mu_slope* (1.d0 - gamma)) * dist_params%norm
      f = f * dist_params%total_norm
      Rf = - f * mu
      where(gamma > multi_slope%gamma_switch) &
        Rf = - f * dist_params%mu_slope
    else
      print*,"Selected distribution is ", dstf
      stop "Unknown distribution flag in make_f_and_Rf_along_line"
    end if
    if(sum(Rf) > 0  .and. sum(Rf) < 1.d-5) then
      f(:) = 0.d0
      Rf(:) = 0.d0
    end if
    if(any(f  /=  f) .or.  any(Rf /= Rf)) then
      print*,  "rhop", svec%rhop
      print*,  "u_perp", u_perp
      print*,  "u_par", u_par
      print*,  "Rf",  Rf
      print*,  "Rf",  f
      if(trim(dstf) == "multi_s") then
        print*, dist_params%total_norm
        print*, dist_params%mu_slope
        print*, dist_params%norm
        print*, multi_slope%gamma_switch
      end if
      call abort()
    end if
end subroutine make_f_and_Rf_along_line

function radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params)
    use constants,                   only: pi, e0, mass_e, c0
    use mod_ecfm_refr_types,         only: dstf, non_therm_params_type, &
                                          rad_diag_ch_mode_ray_freq_svec_type, multi_slope, runaway
    implicit none
    real(rkind), dimension(:),  intent(in)            :: u_par, u_perp, gamma
    real(rkind), intent(in)                           :: mu
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in)  :: svec
    type(non_therm_params_type), intent(in)                :: dist_params
    real(rkind), dimension(size(u_par))               :: radiation_dist_f_u
    real(rkind)                        :: a, lnLambda, alpha, cZ! additional gamma to evaluate f_thermal at the boundaries
    if(dstf == "maxwell") then
      radiation_dist_f_u =  Exp(-mu/2.0 * (u_par**2 + u_perp**2))
    else if(dstf == "Th".or. fall_back_thermal) then
      radiation_dist_f_u =  exp( mu* (1 - gamma))
    else if (dstf == "numeric") then
      print*, "Numeric distribution function will be only evaluated along the entire resonance line"
      stop "bad usage of radiation_dist_f_u"
    else if (dstf == "Bi_MaxJ") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      !print*, "therm part", (1.0 - bi_max%ne) * a * (sqrt(mu / (2 * pi))**3) * &
       ! exp( mu* (1 - gamma))
      !print*, "non-therm part", bi_max%ne * radiation_dist_bi_maxJ(u_par, u_perp, gamma)
      radiation_dist_f_u = (1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
      exp( mu* (1 - gamma)) + dist_params%ne * radiation_dist_bi_maxJ(u_par, u_perp, gamma, &
                                                                      dist_params%Te_perp, &
                                                                      dist_params%Te_par)
    else if (dstf == "Bi_Maxw") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      radiation_dist_f_u = (1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
      exp( mu* (1 - gamma)) + dist_params%ne * radiation_dist_bi_maxw(u_par, u_perp, gamma, &
                                                                      dist_params%Te_perp, &
                                                                      dist_params%Te_par)
#ifdef NAG
    else if (dstf == "drift_m") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      radiation_dist_f_u = a * (sqrt(mu / (2 * pi))**3) * (1.0 - dist_params%ne) * &
        exp( mu* (1 - gamma))  + dist_params%ne * radiation_dist_drift_m(u_par, u_perp, gamma, &
                                                    dist_params%Te_par, &
                                                    dist_params%Te_perp, &
                                                    dist_params%u_par_drift, &
                                                    dist_params%u_perp_drift)
#endif
    else if (dstf == "Spitzer") then
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      radiation_dist_f_u = a * (sqrt(mu / (2 * pi))**3) * &
      exp( mu* (1 - gamma)) * (1.d0 + radiation_dist_Spitzer(u_par, mu, gamma, u_perp, dist_params%vd))
    else if(dstf == "multi_s") then
      radiation_dist_f_u = exp( mu* (1 - gamma))
      where(gamma > multi_slope%gamma_switch) &
                      radiation_dist_f_u = exp(dist_params%mu_slope* (1 - gamma)) * dist_params%norm
      radiation_dist_f_u = radiation_dist_f_u * dist_params%total_norm
    else if(dstf == "runaway") then
      lnLambda = 14.9-0.5*log(svec%ne / 1d20) + log(svec%Te)
      alpha = (runaway%E_E_c - 1.d0) / (runaway%Z_eff+1.d0)
      cZ = sqrt(3.d0*(runaway%Z_eff + 5.d0)/pi)*lnLambda
      radiation_dist_f_u = alpha * dist_params%cur_rel_ne/ (2.d0 * pi * cZ * u_par) * &
        exp( -u_par/(cZ * lnLambda) - 0.5*alpha * u_perp**2 / u_par )
      where(u_par < 0) radiation_dist_f_u = 0.d0
      where(radiation_dist_f_u < 0) radiation_dist_f_u = 0.d0
    else
      print*,"Selected distribution is ", dstf
      stop "Unknown distribution flag in radiation_dist_f_u"
    end if
    if(any(radiation_dist_f_u  < 0.0)  .and. any(( abs(radiation_dist_f_u) > 1.d-10 )) .or. &
      any(radiation_dist_f_u /= radiation_dist_f_u) .or.  any(radiation_dist_f_u  > 1.d6)) then
      print*,  svec%rhop, u_perp, u_par, radiation_dist_f_u
      if (dstf == "Bi_MaxJ") then
        print*, "ne bi maxJ", dist_params%ne
        print*, "f bi maxJ", radiation_dist_bi_maxJ(u_par, u_perp, gamma, &
                                                    dist_params%Te_par, &
                                                    dist_params%Te_perp)
      else if (dstf == "Bi_Maxw") then
        print*, "ne bi max", dist_params%ne
        print*, "f bi maxJ", radiation_dist_bi_maxw(u_par, u_perp, gamma, &
                                                    dist_params%Te_par, &
                                                    dist_params%Te_perp)
#ifdef NAG
      else if (dstf == "drift_m") then
        print*, "ne drift", dist_params%ne
        print*, "f drift", radiation_dist_drift_m(u_par, u_perp, gamma, &
                                                    dist_params%Te_par, &
                                                    dist_params%Te_perp, &
                                                    dist_params%u_par_drift, &
                                                    dist_params%u_perp_drift)
#endif
      else if (dstf == "Spitzer") then
       print*, "f Spitz", radiation_dist_Spitzer(u_par, mu, gamma, u_perp, dist_params%vd)
      else if(dstf == "multi_s") then
        print*, "f mult s", radiation_dist_f_u
      else if(dstf == "runaway") then
        print*, "runaway ne", dist_params%cur_rel_ne
        print*, "f runaway", radiation_dist_f_u
      end if
      print*, "WARNING - Negative distribution in radiation_dist_f_u"
      stop "too large negative distribution or NAN"
    end if
  end function radiation_dist_f_u

  function radiation_dist_f_norm(Te, dstf)
    use constants,                  only:  pi, e0, mass_e, c0
    implicit none
    real(rkind), intent(in)        :: Te
    character(*),  intent(in)      :: dstf
    real(rkind)                    :: radiation_dist_f_norm
    real(rkind)                    :: a, mu
    radiation_dist_f_norm = -1.d0
    if(dstf == "Th" .or. fall_back_thermal .or. dstf == "multi_s") then
      mu = mass_e *c0**2 / (e0 * Te)
      a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
      radiation_dist_f_norm = a * (sqrt(mu / (2 * pi))**3) ! *c0**2c0)
    else if (dstf == "numeric" .or. dstf == "Bi_Maxw" .or. dstf == "drift_m" .or. &
        dstf == "Spitzer" .or. dstf == "runaway" .or. dstf == "Bi_MaxJ"  .or. dstf == "gene" .or. &
        trim(dstf) == "gcomp") then
      radiation_dist_f_norm = 1.d0
      return
    else
      print*,"Selected distribution is ", dstf
      stop "the distribution specified is not implemented"
    end if
  end function radiation_dist_f_norm


  function radiation_dist_Rf(u_par, u_perp, gamma, m_omega_bar, N_par, mu, svec, dist_params)
    use constants,                   only: pi, e0,c0, mass_e
    use mod_ecfm_refr_types,         only: dstf, non_therm_params_type, &
                                          rad_diag_ch_mode_ray_freq_svec_type, multi_slope, runaway
    implicit none
    real(rkind), dimension(:),  intent(in)  :: u_par, u_perp, gamma
    real(rkind), intent(in)            :: m_omega_bar, N_par, mu
    type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
    type(non_therm_params_type), intent(in)                :: dist_params
    real(rkind), dimension(size(u_par)) :: radiation_dist_Rf
    real(rkind)                        :: lnLambda, alpha, cZ, a
    a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
    if (dstf =="Th" .or. fall_back_thermal) then
       radiation_dist_Rf = -radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params) * mu
    else if(dstf == "Bi_Maxw") then
       radiation_dist_Rf = -(1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
                                                exp( mu* (1 - gamma)) * mu + &
                           dist_params%ne * radiation_dist_bi_maxJ_Rf(u_par, u_perp, gamma, m_omega_bar, N_par, &
                                                                      dist_params%Te_par, dist_params%Te_perp)
    else if(dstf == "Bi_MaxJ") then
       radiation_dist_Rf = -(1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
                                                exp( mu* (1 - gamma)) * mu + &
                              dist_params%ne * radiation_dist_bi_maxw_Rf(u_par, u_perp, gamma, m_omega_bar, N_par, &
                                                                        dist_params%Te_par, dist_params%Te_perp)
#ifdef NAG
    else if(dstf == "drift_m") then
       radiation_dist_Rf =-(1.0 - dist_params%ne) * a * (sqrt(mu / (2 * pi))**3) * &
                                                exp( mu* (1 - gamma)) * mu + &
                              dist_params%ne * radiation_dist_drift_m_Rf(u_par, u_perp, gamma,m_omega_bar, N_par, &
                                                                         dist_params%Te_par, dist_params%Te_perp, &
                                                                         dist_params%u_par_drift, dist_params%u_perp_drift)
#endif
    else if(dstf == "multi_s") then
      radiation_dist_Rf = - radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params) * mu
      where(gamma > multi_slope%gamma_switch) &
        radiation_dist_Rf = - radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params) * dist_params%mu_slope
    else if(dstf == "runaway") then
      lnLambda = 14.9-0.5*log(svec%ne / 1d20) + log(svec%Te)
      alpha = (runaway%E_E_c - 1.e0) / (runaway%Z_eff+1)
      cZ = sqrt(3.d0*(runaway%Z_eff + 5.d0)/pi)*lnLambda
      radiation_dist_Rf = (alpha*exp(-(u_par/(cZ*lnLambda)) - (0.5*alpha*u_perp**2)/u_par)*dist_params%cur_rel_ne * &
        (-0.5*alpha*cZ*lnLambda*u_par*m_omega_bar + &
        (-0.5*cZ*lnLambda*u_par - 0.5*u_par**2 + &
        0.25*alpha*cZ*lnLambda*u_perp**2)*N_par)) / &
        (cZ**2*lnLambda*pi*u_par**3)
      where(radiation_dist_Rf > 0 .or. u_par < 0) radiation_dist_Rf = 0.d0
      if(any(radiation_dist_Rf /= radiation_dist_Rf)) then
        print*, u_perp, u_par,-(u_par/(cZ*lnLambda)) - (0.5*alpha*u_perp**2)/u_par
        stop "R f NAN"
      end if
    else if(dstf == "Spitzer") then
       radiation_dist_Rf = radiation_dist_f_u(u_par, u_perp, gamma, mu, svec, dist_params) * ( -mu + &
                            radiation_dist_Spitzer_Rf(u_par, mu, gamma, u_perp, dist_params%vd) * N_par)
    else if (dstf == "numeric") then
      print*, "Numeric distribution function will be only evaluated along the entire resonance line"
      stop "bad usage of radiation_dist_f_u"
    else
      print*,"Selected distribution is ", dstf
      stop "the distribution defined in globals does not exist"
    end if
    if(any(radiation_dist_Rf /= radiation_dist_Rf)) then
      print*,  svec%rhop, u_perp, u_par,radiation_dist_Rf
      print*, "WARNING - NAN in R[f]"
        radiation_dist_Rf = 0.0
    end if
  end function radiation_dist_Rf

function radiation_dist_bi_maxJ(u_par, u_perp, gamma, Te_par, Te_perp)
!   The relativistic Bi-Mawxwellian distribution has been taken from:
!"Numerical computation of Thomson scattering spectra for non-Maxwellian or anisotropic electrion distribution functions"
! By I. Pastor, J. Guasp, R.F., et al.
! Published in Nuclear Fusion 52 (2012)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_types,        only: bi_max
    implicit none
    real(rkind), dimension(:),  intent(in)            :: u_par, u_perp, gamma
    real(rkind), intent(in), optional :: Te_perp, Te_par
    real(rkind), dimension(size(u_par))      :: radiation_dist_bi_maxJ
    real(rkind), dimension(size(u_par))      :: gamma_bi
    real(rkind)                        :: r,s, T0, mu, a
    if(present(Te_perp) .and. present(Te_par)) then
      bi_max%Te_perp = Te_perp
      bi_max%Te_par = Te_par
    end if
!    print*, "Te_perp, Te_par", bi_max%Te_perp, bi_max%Te_par
    T0 = bi_max%Te_par**(1.0d0/3.0d0) * bi_max%Te_perp**(2.0d0/3.0d0)
    mu =(mass_e * c0**2.d0) / (e0 * T0)
    r = T0/bi_max%Te_par
    s = T0/bi_max%Te_perp
    gamma_bi = sqrt(1.0d0 +r*u_par**2 + s*u_perp**2)
    a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))
    radiation_dist_bi_maxJ =a * (sqrt(mu / (2 * pi))**3)* Exp(mu*(1.d0 - gamma_bi))
end function radiation_dist_bi_maxJ

function radiation_dist_bi_maxw(u_par, u_perp, gamma, Te_par, Te_perp)
!   The relativistic Bi-Mawxwellian distribution has been taken from:
!"Numerical computation of Thomson scattering spectra for non-Maxwellian or anisotropic electrion distribution functions"
! By I. Pastor, J. Guasp, R.F., et al.
! Published in Nuclear Fusion 52 (2012)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    real(rkind), dimension(:),  intent(in)            :: u_par, u_perp, gamma
    real(rkind), intent(in)                           :: Te_perp, Te_par
    real(rkind), dimension(size(u_par)) :: radiation_dist_bi_maxw
    real(rkind)                        ::mu_perp, mu_par
    mu_perp =(mass_e * c0**2.d0) / (e0 * Te_perp)
    mu_par =(mass_e * c0**2.d0) / (e0 * Te_par)
    radiation_dist_bi_maxw =(mu_perp / (2.d0 * pi)) *sqrt(mu_par/ (2.d0 * pi))* Exp(-mu_perp/2.d0 * u_perp **2 - mu_par/2.d0 * u_par**2)
end function radiation_dist_bi_maxw

#ifdef NAG
function radiation_dist_drift_m(u_par, u_perp, gamma, Te_par, Te_perp, u_par_drift, u_perp_drift)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_types,        only: drift_m
    implicit none
    real(rkind), dimension(:),  intent(in)            :: u_par, u_perp, gamma
    real(rkind), dimension(size(u_par))         :: radiation_dist_drift_m
    real(rkind), intent(in)                     :: Te_perp, Te_par, u_par_drift, u_perp_drift
    real(rkind)                        :: r,s, T0, mu, gamma_drift_m
    real(rkind)                        :: a
    radiation_dist_drift_m =  radiation_dist_drifting_maxwellian(u_par, u_perp, &
                                Te_par, Te_perp, u_par_drift, u_perp_drift)
end function radiation_dist_drift_m
#endif

function radiation_dist_Spitzer(u, mu, gamma, u_perp, vd) ! u = u_par here
! This function computes the non-thermal Spitzer-Deformation Term
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    real(rkind), dimension(:),  intent(in) :: u, u_perp, gamma
    real(rkind), intent(in)            :: mu
    real(rkind), dimension(size(u)) :: radiation_dist_Spitzer
    real(rkind), dimension(size(u)) :: v
    real(rkind), intent(in)         :: vd
    real(rkind)                        :: u_th, R, gamma_th
    u_th =sqrt(2.d0 / mu + (1.d0/mu**2)) !note here relativistic computation of thermal velocity
    gamma_th = sqrt(1.d0 + u_th**2)
    R = vd/c0 / u_th* gamma_th
    v = (u/gamma)/(u_th/gamma_th)
    radiation_dist_Spitzer =  v/abs(v) * R * &
      0.7619*(0.09476232*v**4 -0.08852586*v**3 +   1.32003051*v**2 -0.19511956*v)
    where(abs(radiation_dist_Spitzer) > 1.0) radiation_dist_Spitzer = 0.d0

end function radiation_dist_Spitzer
#ifdef NAG
function radiation_dist_drifting_maxwellian(u_par, u_perp, Te_par, Te_perp, u_par_drift, u_perp_drift)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_types,        only: drift_m
    use nag_err_fun,                only: nag_erf
    implicit none
    real(rkind), dimension(:),  intent(in) :: u_par, u_perp
    real(rkind), intent(in)            :: Te_perp,Te_par
    real(rkind), intent(in)            :: u_par_drift, u_perp_drift
    real(rkind), dimension(size(u_par)) :: radiation_dist_drifting_maxwellian
    real(rkind)                        :: r,s, T0, mu,cnst
    real(rkind)                        :: a
    real(rkind), dimension(size(u_par)) ::  gamma_drift_m
    !cnst = (mass_e * c0**2.d0) / (e0)
    !mu = cnst /( Te_perp**(2.0/3.0)*Te_par**(1.0/3.0))
    T0 = Te_par**(1.0d0/3.0d0) * Te_perp**(2.0d0/3.0d0)
    mu =(mass_e * c0**2.d0) / (e0 * T0)
    r = T0/Te_par
    s = T0/Te_perp
    gamma_drift_m = sqrt(1.0d0 +r*(u_par - u_par_drift)**2 + s*u_perp**2)

    a = 1.0d0/(1.0d0 + 105.0d0/(128.0d0 * mu**2) + 15.0d0/(8.0d0 * mu))

    radiation_dist_drifting_maxwellian =a * (sqrt(mu / (2 * pi))**3)* Exp(mu*(1- gamma_drift_m))
    !u_th = sqrt(2*abs(bi_max%Te_par)/510998.910 + (bi_max%Te_par/510998.910)**2)
    !a * (sqrt(mu / (2 * pi))**3)* &
    !  Exp(-cnst*((u_par- u_par_drift)**2/ Te_par + (u_perp- u_perp_drift)**2/ Te_perp))
        !R_spitz = bi_max%j/(u_th/sqrt(1.0 + u_th**2)*299792458)
    !v_par_norm = abs(u_par)/sqrt(1 + u_par**2) / u_th * sqrt(1 + u_th**2);
    !a =  1.0 /  (exp(-u_par_drift**2 *mu) + sqrt(pi* u_par_drift* &
    !    (1 + nag_erf(u_par_drift* sqrt(mu)))))
end function radiation_dist_drifting_maxwellian
#endif

function radiation_dist_bi_maxJ_Rf(u_par, u_perp,gamma,n_omega_bar, N_par, Te_par, Te_perp)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_types,        only: bi_max
    implicit none
    real(rkind), dimension(:),  intent(in) :: u_par, u_perp, gamma
    real(rkind), intent(in)            :: n_omega_bar, N_par
    real(rkind),  intent(in), optional :: Te_perp, Te_par
    real(rkind), dimension(size(u_par)) :: radiation_dist_bi_maxJ_Rf
    real(rkind)                        :: r,s, T0, mu
    real(rkind), dimension(size(u_par)) :: gamma_bi_max
    if(present(Te_perp) .and. present(Te_par)) then
      bi_max%Te_perp = Te_perp
      bi_max%Te_par = Te_par
    end if
!    print*, "Te_perp, Te_par Rf", bi_max%Te_perp, bi_max%Te_par
    T0 = bi_max%Te_par**(1.d0/3.d0) * bi_max%Te_perp**(2.d0/3.d0)
    mu = (mass_e * c0**2.d0) / (e0 *T0)
    r = T0/bi_max%Te_par
    s = T0/bi_max%Te_perp
    gamma_bi_max = sqrt(1.d0 +r*u_par**2 + s*u_perp**2)
    radiation_dist_bi_maxJ_Rf = - radiation_dist_bi_maxJ(u_par,u_perp,gamma, Te_par, Te_perp) * &
                                  mu*(n_omega_bar * Te_par + N_par * u_par * Te_perp) / &
                                  gamma_bi_max * (T0/(Te_par * Te_perp))
end function radiation_dist_bi_maxJ_Rf

function radiation_dist_bi_maxw_Rf(u_par, u_perp,gamma,n_omega_bar, N_par, Te_par, Te_perp)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    real(rkind), dimension(:),  intent(in) :: u_par, u_perp, gamma
    real(rkind), intent(in)            :: n_omega_bar, N_par, Te_perp, Te_par
    real(rkind), dimension(size(u_par))   :: radiation_dist_bi_maxw_Rf
    real(rkind)                        :: mu_perp, mu_par
    mu_perp = (mass_e * c0**2.d0) / (e0 * Te_perp)
    mu_par = (mass_e * c0**2.d0) / (e0 * Te_par)
    radiation_dist_bi_maxw_Rf = -radiation_dist_bi_maxw(u_par,u_perp,gamma, Te_par, Te_perp) * &
                               (n_omega_bar * mu_perp + N_par * u_par * mu_par)
end function radiation_dist_bi_maxw_Rf

#ifdef NAG
function radiation_dist_drift_m_Rf(u_par, u_perp,gamma, n_omega_bar, &
                                  N_par, Te_par, Te_perp, u_par_drift, &
                                  u_perp_drift)
    use constants,                  only: pi, e0, mass_e, eps0, c0
    use mod_ecfm_refr_types,        only: drift_m
    implicit none
    real(rkind), dimension(:),  intent(in) :: u_par, u_perp, gamma
    real(rkind), intent(in)            :: n_omega_bar, N_par
    real(rkind), dimension(size(u_par)) :: radiation_dist_drift_m_Rf
    real(rkind), intent(in)            :: Te_perp,Te_par
    real(rkind), intent(in)            :: u_par_drift, u_perp_drift
    real(rkind)                        :: mu, T0, r,s
    real(rkind), dimension(size(u_par)) :: gamma_drift_m
    T0 = Te_par**(1.d0/3.d0) * Te_perp**(2.d0/3.d0)
    mu = (mass_e * c0**2.d0) / (e0 *T0)
    r = T0/Te_par
    s = T0/Te_perp
    gamma_drift_m= sqrt(1 +r*(u_par- u_par_drift) **2 + s*u_perp**2)
    radiation_dist_drift_m_Rf = -radiation_dist_drifting_maxwellian(u_par,u_perp,&
                                Te_par,Te_perp,u_par_drift, &
                                u_perp_drift)  * &
                                mu*(n_omega_bar * Te_par + N_par * (u_par - &
                                u_par_drift) * Te_perp) / &
                                gamma_drift_m * (T0/(Te_par * Te_perp))
    ! -radiation_dist_drifting_maxwellian(u_par,u_perp,&
    !                            drift_m%Te_par_supra,drift_m%Te_perp_supra,drift_m%u_par_supra_drift, &
    !                            drift_m%u_perp_supra_drift) * mu_supra * &
    !                            (2 * n_omega_bar * (u_perp + drift_m%u_perp_supra_drift) / u_perp &
    !                            + N_par * (2 * u_par - drift_m%u_par_supra_drift)) + &
    !                            radiation_dist_drifting_maxwellian(u_par,u_perp, &
    !                            drift_m%Te_par_ECRH,drift_m%Te_perp_ECRH,drift_m%u_par_ECRH_drift, &
    !                            drift_m%u_perp_ECRH_drift) * mu_ECRH* &
    !                            (2 * n_omega_bar * (u_perp + drift_m%u_perp_ECRH_drift) / u_perp &
    !                            + N_par * (2 * u_par - drift_m%u_par_ECRH_drift))
end function radiation_dist_drift_m_Rf
#endif NAg
function radiation_dist_Spitzer_Rf(u, mu, gamma, u_perp, vd) ! u = u_par here
! This function computes the non-thermal Spitzer-Deformation Term
    use constants,                  only: pi, e0, mass_e, eps0, c0
    implicit none
    real(rkind), dimension(:),  intent(in) :: u, u_perp, gamma
    real(rkind), intent(in)            :: mu
    real(rkind), dimension(size(u)) :: radiation_dist_Spitzer_Rf
    real(rkind), intent(in)            :: vd
    real(rkind)                        :: u_th
    real(rkind), dimension(size(u)) :: R, v, gamma_th, chain_rule
    u_th =sqrt(2.d0 / mu + (1.d0/mu**2)) !note here relativistic computation of thermal velocity
    gamma_th = sqrt(1.d0 + u_th**2 + u_perp**2)
    R = vd/c0 / u_th** gamma_th
    v = (u/gamma)/(u_th/gamma_th)
    chain_rule = (1 + u_perp**2) / (gamma *(u_th/gamma_th)) ! from d/du v
    radiation_dist_Spitzer_Rf = mu * 2.d0 * u / gamma * radiation_dist_Spitzer(u, mu, gamma, u_perp, vd)
    ! From f_Sh *d/du_{||} f_MJ first term in product rule
    radiation_dist_Spitzer_Rf = radiation_dist_Spitzer_Rf +  v/abs(v) * R * &
        chain_rule * 0.7619 * (4.d0 * 0.09476232*v**3 -3.d0 * 0.08852586*v**3 + &
          2.d0 * 1.32003051*v - 0.19511956)
end function radiation_dist_Spitzer_Rf
end module mod_ecfm_radiation_dist
