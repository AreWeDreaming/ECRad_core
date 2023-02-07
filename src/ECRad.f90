program ECRad
    use f90_kind
    use mod_ECRad_types,        only: rad, data_folder, output_level, &
                                          ant, stand_alone, eval, not_eval, &
                                          data_name, modes
    use mod_ECRad, only : initialize_stand_alone, make_ece_rad_temp!, simulate_ida
    use mod_ECRad_utils, only : export_all_ece_data
! #ifdef OMP
!     use omp_lib
! #endif
    implicit none
    Character(len=200) :: filename
    Character(len=200) :: working_dir
    integer(ikind) :: idiag, ich!, ifreq, ir
    call getarg(1, working_dir)
    ! To test the ECRad integration into ida
    ! For stand alone usage comment out the following line
    stand_alone = .True.
    !stand_alone = .False.
    eval = 0
    not_eval = 0
!    if(stand_alone) then
    !stop "Sanity check"
      call initialize_stand_alone(working_dir, "init")
      call make_ece_rad_temp()
      !write(filename, "(A64,A7)") data_folder, "tau.dat"
      !open(67, file=filename)
      !open(67, file="TRadM_s.dat")
      if(.not. output_level) then
        filename = trim(data_folder) // trim(data_name)
        open(66, file=filename)
        filename = trim(data_folder) // "sres.dat"
        open(69, file=filename)
        filename = trim(data_folder) // "diag.dat"
        open(70,file=filename)
        if(modes == 1 .or. modes == 3) then
          filename = trim(data_folder) // "X_" // trim(data_name)
          open(67, file=filename)
        end if
        if(modes == 2 .or. modes == 3) then
          filename = trim(data_folder) // "O_" // trim(data_name)
          open(68, file=filename)
        end if
        do idiag = 1, ant%N_diag
          do ich = 1,ant%diag(idiag)%N_ch
            write(66,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%rhop_res," ",&
            rad%diag(idiag)%ch(ich)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%tau
            write(69,"(E16.8E3A1E16.8E3A1E16.8E3A1E16.8E3)") rad%diag(idiag)%ch(ich)%s_res, " ", rad%diag(idiag)%ch(ich)%R_res, " ",&
                                    rad%diag(idiag)%ch(ich)%z_res, " ", rad%diag(idiag)%ch(ich)%rhop_res
            write(70,"(A3)") ant%diag(idiag)%diag_name
            if(modes == 1) then
               write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(1)%rhop_res," ",&
                          rad%diag(idiag)%ch(ich)%mode(1)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%mode(1)%tau
            else if(modes == 2) then
              write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(1)%rhop_res," ",&
                rad%diag(idiag)%ch(ich)%mode(1)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%mode(1)%tau
            else
              write(67,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(1)%rhop_res," ",&
                rad%diag(idiag)%ch(ich)%mode(1)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%mode(1)%tau
              write(68,"(E16.8E3,A1,E16.8E3,A1,E16.8E3)") rad%diag(idiag)%ch(ich)%mode(2)%rhop_res," ",&
                rad%diag(idiag)%ch(ich)%mode(2)%TRad / 1000.0d0, " ", rad%diag(idiag)%ch(ich)%mode(2)%tau
             end if
            !write(67,"(E12.4E3,A1,E12.4E3)") rad%diag(idiag)%ch(ich)%ray(1)%freq(1)%rhop_res_2X," ",&
            !rad%diag(idiag)%ch(ich)%tau
            !write(67,"(E12.4E3,A1,E12.4E3)") rad%diag(idiag)%ch(ich)%ray(1)%freq(2)%s_res_2X," ",&
            !TRad(ich) / 1000.0d0
          !print*, rad%diag(idiag)%ch(ich)%ray(1)%freq(2)%rhop_res_2X, TRad(ich) / 1000.0d0
          end do
        end do
        close(66)
        if(modes == 1 .or. modes == 3) close(67)
        if(modes == 2 .or. modes == 3) close(68)
        close(69)
        close(70)
      end if
!    else
!      call simulate_ida(working_dir)
!    end if
!    if(output_level .and. eval + not_eval > 0)  print*, "Saved", not_eval, " of ", not_eval + eval, "evaluations: ", real(not_eval,8) / real( eval + not_eval,8) * 100.d0, "%"
!    if(dstf == "numeric" .or. dstf_comp == "O1") then
!      dstf = "relamax"
!      output_level = .False.
!      call make_ece_rad_temp(1,reflec,TRad)
!      write(filename, "(A64,A15)") data_folder, "TRadM_therm.dat"
!      open(66, file=filename)
!      !write(filename, "(A64,A7)") data_folder, "tau_therm.dat"
!      !open(67, file=filename)
!      !open(67, file="TRadM_s.dat")
!      do ich = 1,ant%N_ch
!        write(66,"(E12.4E3,A1,E12.4E3)") rad%diag(idiag)%ch(ich)%ray(1)%freq(1)%rhop_res_2X," ",&
!        TRad(ich) / 1000.0d0
!        !write(67,"(E12.4E3,A1,E12.4E3)") rad%diag(idiag)%ch(ich)%ray(1)%freq(1)%rhop_res_2X," ",&
!        !rad%diag(idiag)%ch(ich)%tau
!        !write(67,"(E12.4E3,A1,E12.4E3)") rad%diag(idiag)%ch(ich)%ray(1)%freq(2)%s_res_2X," ",&
!        !TRad(ich) / 1000.0d0
!      !print*, rad%diag(idiag)%ch(ich)%ray(1)%freq(2)%rhop_res_2X, TRad(ich) / 1000.0d0
!      end do
!      close(66)
!      !close(67)
!    end if
    if(output_level) print*,"ECRad successful"
end program ECRad
