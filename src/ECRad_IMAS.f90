program ECRad_IMAS
    use f90_kind
#ifdef OMP
  use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
                             set_ECRad_thread_count, initialize_ECRad_IMAS, &
                             reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#else
  use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
                             initialize_ECRad_IMAS, &
                             reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#endif

    use ids_schemas, only: ids_equilibrium, ids_core_profiles, ids_ece, ids_wall
    use ids_routines, only: imas_open_env,imas_create_env,imas_close,ids_get,ids_put,ids_deallocate
    use f90_file_reader, only: file2buffer
    implicit none
   
    integer:: idx,error_flag
    type(ids_wall):: wall
    type(ids_equilibrium):: equilibrium
    type(ids_core_profiles):: core_profiles
    type(ids_ece):: ece_in, ece_out
    character(len=200):: user,machine
    character(len=:), pointer:: error_message
    character(len=132), pointer :: buffer(:) => NULL()
    integer :: io_unit = 1
    
    call file2buffer("../src/code_params.xml", io_unit, buffer)

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'iter'

    ! OPEN INPUT DATAFILE FROM OFFICIAL IMAS SCENARIO DATABASE
    write(*,*) '=> Read input IDSs'
    call imas_open_env('ids',116000,2,idx,'public','iter','3')
    call ids_get(idx,'wall', wall)
    call imas_close(idx)
    write(*,*) 'Finished reading input IDSs'

    call imas_open_env('ids',104103,13,idx,'public','iter','3')
    call ids_get(idx,'equilibrium', equilibrium)
    call imas_close(idx)
    write(*,*) 'Finished reading input IDSs'

    call imas_open_env('ids',104103,13,idx,'public','iter','3')
    call ids_get(idx,'core_profiles', core_profiles)
    call imas_close(idx)
    write(*,*) 'Finished reading input IDSs'

    call imas_open_env('ids',150601,1,idx,'public','iter','3')
    call ids_get(idx,'ece', ece_in)
    call imas_close(idx)
    write(*,*) 'Finished reading input IDSs'

   
    ! EXECUTE PHYSICS CODE
    call pre_initialize_ECRad_IMAS(buffer, wall, error_flag, error_message)
    call set_ece_ECRad_IMAS(ece_in, 0, error_flag, error_message)
    call initialize_ECRad_IMAS(equilibrium, 0, error_flag,  &
                               error_message)
    call make_rays_ECRad_IMAS(core_profiles, 0)
    if(error_flag.eq.0) then
      
       ! EXPORT RESULTS TO LOCAL DATABASE
       write(*,*) '=> Export output IDSs to local database'
       call imas_create_env('ids',104103,1,0,0,idx,trim(user),trim(machine),'3')
       call ids_put(idx,'ece', ece_out)
       call make_dat_model_ECRad(core_profiles, 0, ece_out)
       call imas_close(idx)
       write(*,*) 'Done exporting.'
       write(*,*) ' '
       write(*,*) 'End of standalone'
    else
       write(*,*) error_message
       write(*,*) '=> Program stopped.'
   
    endif
   
    ! RELEASING THE MEMORY ALLOCATED FOR EACH IDS
    call ids_deallocate(wall)
    call ids_deallocate(core_profiles)
    call ids_deallocate(ece_in)
    call ids_deallocate(ece_out)
end program ECRad_IMAS
