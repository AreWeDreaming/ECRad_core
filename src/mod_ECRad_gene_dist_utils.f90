module mod_ECRad_gene_dist_utils
  use f90_kind
  implicit none

  public :: setup_fgene_rhop_splines, &
            make_g_inter, &
            make_g_and_g_grad_along_line, &
            make_gene_f_and_gene_f_grad_along_line, &
            make_gene_f0_and_gene_f0_grad_along_line

  private :: v_par_mu_to_cyl

contains

subroutine setup_fgene_rhop_splines(fgene)
  use mod_ECRad_types,          only: fgene_type
  use mod_ECRad_interpol,       only: make_1d_spline
#ifdef NAG
  USE nag_spline_1d,                only: nag_spline_1d_interp
#endif
  implicit none
  type(fgene_type), intent(inout)   :: fgene
  integer(ikind)                  :: ivpar, imu
  ! N_vpar x N_mu splines are set up along rhop
  allocate(fgene%g_rhop_spl(fgene%N_vpar, fgene%N_mu))
#ifdef NAG
  allocate(fgene%g_rhop_nag_spl(fgene%N_vpar, fgene%N_mu))
#endif
  do ivpar = 1, fgene%N_vpar
    do imu = 1, fgene%N_mu
    ! No nag spline here, because wwe do not want compare linear and cubic splines
      call make_1d_spline(fgene%g_rhop_spl(ivpar, imu), fgene%N_rhop, fgene%rhop, fgene%g(:, ivpar, imu), k = 1)
    end do
  end do
end subroutine setup_fgene_rhop_splines

subroutine make_g_inter(svec, g_spl)
  use mod_ECRad_types,          only: fgene, rad_diag_ch_mode_ray_freq_svec_type, N_absz_large
  use mod_ECRad_interpol,       only: make_rect_spline, deallocate_rect_spline, rect_spline, rect_spline_vec, spline_1d, spl_type_2d
#ifdef NAG
  USE nag_spline_2d,                only: nag_spline_2d_interp
#endif
  use constants,                    only: pi, mass_e, e0, c0
  implicit none
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  type(spl_type_2d), intent(inout) :: g_spl
  integer(ikind)                  :: ivpar, imu
  real(rkind)                     :: rhop_best
  real(rkind), dimension(fgene%N_vpar, fgene%N_mu) :: g_inter
  rhop_best = svec%rhop
  if(rhop_best > fgene%rhop_max) return ! automatic switch to thermal distributions when the distribution function is evaluated
  if(rhop_best < fgene%rhop_min) rhop_best = fgene%rhop_min ! if rhop grid sensible this should be fine
!  irhop = minloc(abs(fgene%rhop-rhop_best), dim=1)
  do ivpar = 1, fgene%N_vpar
!    fgene%g_inter(ivpar,:) = fgene%g(irhop,ivpar,:)
    do imu = 1, fgene%N_mu
    ! No nag spline here, because wwe do not want compare linear and cubib splines
      call spline_1d(fgene%g_rhop_spl(ivpar, imu), rhop_best, g_inter(ivpar, imu))
    end do
  end do
  g_spl%iopt_int = 0 ! restart
  call make_rect_spline(g_spl, int(size(fgene%vpar),4), int(size(fgene%mu),4), &
       fgene%vpar, fgene%mu, g_inter, iopt = 0, m_max = N_absz_large)
!#ifdef NAG
!  call nag_spline_2d_interp(fgene%vpar, fgene%mu, g_inter, fgene%g_nag_spl)
!#endif
end subroutine make_g_inter

subroutine v_par_mu_to_cyl(u_par, u_perp, svec, vpar, mu)
  use mod_ECRad_types,          only: fgene,rad_diag_ch_mode_ray_freq_svec_type, fgene
  use constants,                    only: mass_e, pi, e0, c0
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec ! Not used here, but could be used in future to compute actual B_0
  real(rkind),  dimension(size(u_par)), intent(out)         :: vpar, mu
  real(rkind), dimension(size(u_par))                       :: gam
  gam = sqrt(1.d0 + u_par**2 + u_perp**2)
  vpar = u_par / gam
  mu = (u_perp/gam)**2 / (2.d0 * fgene%B0)
end subroutine v_par_mu_to_cyl

subroutine make_g_and_g_grad_along_line(u_par, u_perp, svec, g_spl, g, dg_du_par, dg_du_perp, debug)
  use mod_ECRad_types,          only: fgene, rad_diag_ch_mode_ray_freq_svec_type
  use mod_ECRad_interpol,       only: spl_type_2d
  use constants,                    only: mass_e, pi, e0, c0
  use mod_ECRad_interpol,       only: rect_spline_vec
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  type(spl_type_2d), intent(in) :: g_spl
  real(rkind), dimension(:), intent(out) :: g, dg_du_par, dg_du_perp
  logical, intent(in), optional          :: debug
  logical                                :: dbg
  integer(ikind)                  :: i
  real(rkind), dimension(size(u_par)) :: temp_vpar, temp_mu, vpar, mu, dg_dvpar, dg_dmu,gam, &
                                         dvpar_du_par, dvpar_du_perp, &
                                         dmu_du_par, dmu_du_perp, vpar_step, mu_step
  gam = Sqrt(1.d0 + u_perp**2 + u_par**2)
  dbg = .true.
  if(present(debug)) dbg = debug
  call v_par_mu_to_cyl(u_par, u_perp, svec, vpar, mu)
  temp_vpar = vpar
  temp_mu = mu
  ! Avoid interpolation errors - these values will be replaced by zeros later
  where(temp_vpar <= fgene%vpar_min) temp_vpar = fgene%vpar_min + 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_vpar >= fgene%vpar_max) temp_vpar = fgene%vpar_max - 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_mu <= fgene%mu_min) temp_mu = fgene%mu_min + 1.d-5 * (fgene%mu_max - fgene%mu_min)
  where(temp_mu >= fgene%mu_max) temp_mu = fgene%mu_max - 1.d-5 * (fgene%mu_max - fgene%mu_min)
!#ifdef NAG
!  call rect_spline_vec(fgene%g_spl, temp_vpar, temp_mu, g, dg_dvpar, dg_dmu, fgene%g_nag_spl)
!#else
  call rect_spline_vec(g_spl, temp_vpar, temp_mu, g, dg_dvpar, dg_dmu)
!#endif
  dvpar_du_par = (1.d0 + u_perp**2) / gam**3
  dvpar_du_perp = -(u_par * u_perp) / gam**3
  dmu_du_par = -(u_par * u_perp**2) / (fgene%B0 * gam**4)
  dmu_du_perp = (1.d0 + u_par**2) * u_perp / (fgene%B0 * gam**4)
  if(dbg) then
    open(93, file="finite_difg_par.txt")
    open(94, file="finite_difg_perp.txt")
    call v_par_mu_to_cyl(u_par + 1.d-4, u_perp, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*,"dpar", u_par(i), dmu_du_par(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_par(i), (vpar_step - vpar(i)) / 1.d-4
      write(93,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_par(i), " ",&
              dmu_du_par(i),  " ",&
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_par(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4

    end do
    call v_par_mu_to_cyl(u_par, u_perp + 1.d-4, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*, "dperp",u_perp(i), dmu_du_perp(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_perp(i), (vpar_step - vpar(i)) / 1.d-4
      write(94,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_perp(i), " ",&
              dmu_du_perp(i), " ", &
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_perp(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4
   end do
  close(93)
  close(94)
  end if
  dg_du_par = (dg_dvpar * dvpar_du_par + dmu_du_par * dg_dmu) * g
  dg_du_perp = (dg_dvpar * dvpar_du_perp + dmu_du_perp * dg_dmu) * g
  do i = 1, size(u_par)
  ! Set all values that are not in the defined region of the distribution pertubation to 0
    if(vpar(i) >= fgene%vpar_max .or. vpar(i) <= fgene%vpar_min .or. mu(i) <= fgene%mu_min .or. mu(i) >= fgene%mu_max) then
    ! No pertubation of distribution
      g(i) = 0.d0
      dg_du_par(i) = 0.d0
      dg_du_perp(i) = 0.d0
    end if
  end do
end subroutine make_g_and_g_grad_along_line

subroutine make_gene_f_and_gene_f_grad_along_line(u_par, u_perp, svec, g_spl, f, df_du_par, df_du_perp, debug)
  use mod_ECRad_types,          only: fgene, rad_diag_ch_mode_ray_freq_svec_type
  use mod_ECRad_interpol,       only: spl_type_2d
  use constants,                    only: mass_e, pi, e0, c0
  use mod_ECRad_interpol,       only: rect_spline_vec
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  type(spl_type_2d), intent(in) :: g_spl
  real(rkind), dimension(:), intent(out) :: f, df_du_par, df_du_perp
  logical, intent(in), optional          :: debug
  logical                                :: dbg
  integer(ikind)                  :: i
  real(rkind), dimension(size(u_par)) :: temp_vpar, temp_mu, vpar, mu, df_dvpar, df_dmu,gam, &
                                         dvpar_du_par, dvpar_du_perp, &
                                         dmu_du_par, dmu_du_perp, vpar_step, mu_step
  gam = Sqrt(1.d0 + u_perp**2 + u_par**2)
  dbg = .false.
  if(present(debug)) dbg = debug
  call v_par_mu_to_cyl(u_par, u_perp, svec, vpar, mu)
  temp_vpar(:) = vpar(:)
  temp_mu(:) = mu(:)
  ! Avoid interpolation errors - these values will be replaced by zeros later
  where(temp_vpar <= fgene%vpar_min) temp_vpar = fgene%vpar_min + 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_vpar >= fgene%vpar_max) temp_vpar = fgene%vpar_max - 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_mu <= fgene%mu_min) temp_mu = fgene%mu_min + 1.d-5 * (fgene%mu_max - fgene%mu_min)
  where(temp_mu >= fgene%mu_max) temp_mu = fgene%mu_max - 1.d-5 * (fgene%mu_max - fgene%mu_min)
!#ifdef NAG
!  call rect_spline_vec(fgene%g_spl, temp_vpar, temp_mu, f, df_dvpar, df_dmu, fgene%g_nag_spl)
!#else
  call rect_spline_vec(g_spl, temp_vpar, temp_mu, f, df_dvpar, df_dmu)
!#endif
  f = exp(f)
  dvpar_du_par = (1.d0 + u_perp**2) / gam**3
  dvpar_du_perp = -(u_par * u_perp) / gam**3
  dmu_du_par = -(u_par * u_perp**2) / (fgene%B0 * gam**4)
  dmu_du_perp = (1.d0 + u_par**2) * u_perp / (fgene%B0 * gam**4)
  df_du_par = (df_dvpar * dvpar_du_par + dmu_du_par * df_dmu) * f
  df_du_perp = (df_dvpar * dvpar_du_perp + dmu_du_perp * df_dmu) * f
  do i = 1, size(u_par)
  ! Set all values that are not in the defined region of the distribution pertubation to 0
    if(vpar(i) >= fgene%vpar_max .or. vpar(i) <= fgene%vpar_min .or. mu(i) <= fgene%mu_min .or. mu(i) >= fgene%mu_max) then
    ! No pertubation of distribution
      f(i) = 0.d0
      df_du_par(i) = 0.d0
      df_du_perp(i) = 0.d0
    end if
  end do
  if(dbg) then
    open(93, file="finite_difg_par.txt")
    open(94, file="finite_difg_perp.txt")
    call v_par_mu_to_cyl(u_par + 1.d-4, u_perp, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*,"dpar", u_par(i), dmu_du_par(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_par(i), (vpar_step - vpar(i)) / 1.d-4
      write(93,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_par(i), " ",&
              dmu_du_par(i),  " ",&
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_par(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4
    end do
    call v_par_mu_to_cyl(u_par, u_perp + 1.d-4, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*, "dperp",u_perp(i), dmu_du_perp(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_perp(i), (vpar_step - vpar(i)) / 1.d-4
      write(94,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_perp(i), " ",&
              dmu_du_perp(i), " ", &
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_perp(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4
   end do
  close(93)
  close(94)
  end if
end subroutine make_gene_f_and_gene_f_grad_along_line

subroutine make_gene_f0_and_gene_f0_grad_along_line(u_par, u_perp, svec, f0, df0_du_par, df0_du_perp, debug)
  use mod_ECRad_types,          only: fgene, rad_diag_ch_mode_ray_freq_svec_type
  use constants,                    only: mass_e, pi, e0, c0
  use mod_ECRad_interpol,       only: rect_spline_vec
  implicit none
  real(rkind), dimension(:), intent(in)         :: u_par, u_perp!
  type(rad_diag_ch_mode_ray_freq_svec_type), intent(in):: svec
  real(rkind), dimension(:), intent(out) :: f0, df0_du_par, df0_du_perp
  logical, intent(in), optional          :: debug
  logical                                :: dbg
  integer(ikind)                  :: i
  real(rkind), dimension(size(u_par)) :: temp_vpar, temp_mu, vpar, mu, df_dvpar, df_dmu,gam, &
                                         dvpar_du_par, dvpar_du_perp, &
                                         dmu_du_par, dmu_du_perp, vpar_step, mu_step
  gam = Sqrt(1.d0 + u_perp**2 + u_par**2)
  dbg = .false.
  if(present(debug)) dbg = debug
  call v_par_mu_to_cyl(u_par, u_perp, svec, vpar, mu)
  temp_vpar(:) = vpar(:)
  temp_mu(:) = mu(:)
  ! Avoid interpolation errors - these values will be replaced by zeros later
  where(temp_vpar <= fgene%vpar_min) temp_vpar = fgene%vpar_min + 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_vpar >= fgene%vpar_max) temp_vpar = fgene%vpar_max - 1.d-5 * (fgene%vpar_max - fgene%vpar_min )
  where(temp_mu <= fgene%mu_min) temp_mu = fgene%mu_min + 1.d-5 * (fgene%mu_max - fgene%mu_min)
  where(temp_mu >= fgene%mu_max) temp_mu = fgene%mu_max - 1.d-5 * (fgene%mu_max - fgene%mu_min)
#ifdef NAG
  call rect_spline_vec(fgene%f0_spl, temp_vpar, temp_mu, f0, df_dvpar, df_dmu, fgene%f0_nag_spl)
#else
  call rect_spline_vec(fgene%f0_spl, temp_vpar, temp_mu, f0, df_dvpar, df_dmu)
#endif
  f0 = exp(f0)
  dvpar_du_par = (1.d0 + u_perp**2) / gam**3
  dvpar_du_perp = -(u_par * u_perp) / gam**3
  dmu_du_par = -(u_par * u_perp**2) / (fgene%B0 * gam**4)
  dmu_du_perp = (1.d0 + u_par**2) * u_perp / (fgene%B0 * gam**4)
  df0_du_par = (df_dvpar * dvpar_du_par + dmu_du_par * df_dmu) * f0
  df0_du_perp = (df_dvpar * dvpar_du_perp + dmu_du_perp * df_dmu) * f0
  do i = 1, size(u_par)
  ! Set all values that are not in the defined region of the distribution pertubation to 0
    if(vpar(i) >= fgene%vpar_max .or. vpar(i) <= fgene%vpar_min .or. mu(i) <= fgene%mu_min .or. mu(i) >= fgene%mu_max) then
    ! No pertubation of distribution
      f0(i) = 0.d0
      df0_du_par(i) = 0.d0
      df0_du_perp(i) = 0.d0
    end if
  end do
  if(dbg) then
    open(93, file="finite_difg_par.txt")
    open(94, file="finite_difg_perp.txt")
    call v_par_mu_to_cyl(u_par + 1.d-4, u_perp, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*,"dpar", u_par(i), dmu_du_par(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_par(i), (vpar_step - vpar(i)) / 1.d-4
      write(93,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_par(i), " ",&
              dmu_du_par(i),  " ",&
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_par(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4

    end do
    call v_par_mu_to_cyl(u_par, u_perp + 1.d-4, svec, vpar_step, mu_step)
    do i = 1, size(u_par)
      !print*, "dperp",u_perp(i), dmu_du_perp(i), (mu_step - mu(i)) / 1.d-4, dvpar_du_perp(i), (vpar_step - vpar(i)) / 1.d-4
      write(94,"(E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3,A1,E18.10E3)") &
              u_perp(i), " ",&
              dmu_du_perp(i), " ", &
              (mu_step(i) - mu(i)) / 1.d-4, " ", &
              dvpar_du_perp(i), " ", &
              (vpar_step(i) - vpar(i)) / 1.d-4
   end do
  close(93)
  close(94)
  end if
end subroutine make_gene_f0_and_gene_f0_grad_along_line

end module mod_ECRad_gene_dist_utils
