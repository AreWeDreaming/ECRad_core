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
  use mod_ECRad_actor_IMAS
  use ids_schemas, only: ids_equilibrium, ids_core_profiles, &
       ids_ece, ids_wall, ids_parameters_input
  use ids_routines, only: imas_open_env,imas_create_env,imas_close, &
       ids_get,ids_get_slice,ids_put,ids_deallocate
  use f90_file_reader, only: file2buffer
  use mod_codeparam_standalone_IMAS
  implicit none

  type(ids_wall):: wall
  type(ids_equilibrium):: equilibrium
  type(ids_core_profiles):: core_profiles
  type(ids_ece):: ece_in, ece_out
  type(ids_parameters_input):: codeparam_standalone,codeparam_ecrad
  type(type_codeparam_standalone):: codeparam_standalone_data  
  integer:: io_unit = 1, idx
  integer:: shot_scenario,run_scenario,shot_wall,run_wall,shot_ece,run_ece,run_out
  character(len=200):: user_scenario,db_scenario,user_wall
  character(len=200):: db_wall,user_ece,db_ece,local_db,local_user
  double precision:: time_slice
  character(len=:), pointer:: output_message
  integer:: output_flag

  ! READ STANDALONE XML INPUT FILE, DEFINED BY SPECIFIC XSD FILE
  call file2buffer('input/standalone.xml',io_unit, codeparam_standalone%parameters_value)
  call parse_standalone_codeparam(codeparam_standalone%parameters_value,codeparam_standalone_data)
  shot_scenario = codeparam_standalone_data%shot_scenario
  run_scenario  = codeparam_standalone_data%run_scenario
  user_scenario = codeparam_standalone_data%user_scenario
  db_scenario   = codeparam_standalone_data%db_scenario
  shot_wall     = codeparam_standalone_data%shot_wall
  run_wall      = codeparam_standalone_data%run_wall
  user_wall     = codeparam_standalone_data%user_wall
  db_wall       = codeparam_standalone_data%db_wall
  shot_ece      = codeparam_standalone_data%shot_ece
  run_ece       = codeparam_standalone_data%run_ece
  user_ece      = codeparam_standalone_data%user_ece
  db_ece        = codeparam_standalone_data%db_ece
  time_slice    = codeparam_standalone_data%time_slice
  run_out       = codeparam_standalone_data%run_out
  local_db      = codeparam_standalone_data%local_db

  write(*,*) '------------------------------------'
  write(*,*) 'Parameters read from input xml file:'
  write(*,*) '------------------------------------'
  write(*,*) '- SCENARIO'
  write(*,'(a21,i7)')  '     shot_scenario = ',codeparam_standalone_data%shot_scenario
  write(*,'(a21,i7)')  '     run_scenario  = ',codeparam_standalone_data%run_scenario
  write(*,'(a21,a30)') '     user_scenario = ',codeparam_standalone_data%user_scenario
  write(*,'(a21,a30)') '     db_scenario   = ',codeparam_standalone_data%db_scenario
  write(*,*) '- WALL MACHINE DESCRIPTION'
  write(*,'(a21,i7)')  '     shot_wall     = ',codeparam_standalone_data%shot_wall
  write(*,'(a21,i7)')  '     run_wall      = ',codeparam_standalone_data%run_wall
  write(*,'(a21,a30)') '     user_wall     = ',codeparam_standalone_data%user_wall
  write(*,'(a21,a30)') '     db_wall       = ',codeparam_standalone_data%db_wall
  write(*,*) '- ECE MACHINE DESCRIPTION'
  write(*,'(a21,i7)')  '     shot_ece      = ',codeparam_standalone_data%shot_ece
  write(*,'(a21,i7)')  '     run_ece       = ',codeparam_standalone_data%run_ece
  write(*,'(a21,a30)') '     user_ece      = ',codeparam_standalone_data%user_ece
  write(*,'(a21,a30)') '     db_ece        = ',codeparam_standalone_data%db_ece
  write(*,*) '- OTHER PARAMETERS'
  write(*,'(a21,f7.3)')'     time_slice    = ',codeparam_standalone_data%time_slice
  write(*,'(a21,i7)')  '     run_out       = ',codeparam_standalone_data%run_out
  write(*,'(a21,a30)') '     local_db      = ',codeparam_standalone_data%local_db
  write(*,*) '------------------------------------'

  ! USERNAME DEFINED BY ENVIRONMENT VARIABLE USERNAME
  call getenv('USER',local_user)

  ! OPEN INPUT DATAFILE FROM OFFICIAL IMAS SCENARIO DATABASE
  write(*,*) 'Read input IDSs'

  write(*,*) '  --> equilibrium IDS'
  call imas_open_env('ids',shot_scenario,run_scenario,idx,user_scenario,db_scenario,'3')
  call ids_get_slice(idx,'equilibrium', equilibrium,time_slice,1,output_flag)
  write(*,*) '  --> core_profiles IDS'
  call ids_get_slice(idx,'core_profiles', core_profiles,time_slice,1,output_flag)
  call imas_close(idx)

  write(*,*) '  --> wall IDS'
  call imas_open_env('ids',shot_wall,run_wall,idx,user_wall,db_wall,'3')
  call ids_get(idx,'wall', wall)
  call imas_close(idx)

  write(*,*) '  --> ece IDS'
  call imas_open_env('ids',shot_ece,run_ece,idx,user_ece,db_ece,'3')
  call ids_get(idx,'ece', ece_in)
  call imas_close(idx)

  write(*,*) 'Finished reading input IDSs'
  write(*,*) '------------------------------------'

  ! PREPARE AND RUN ECRAD
  call file2buffer("input/code_params.xml", io_unit, codeparam_ecrad%parameters_value)
  call ECRad_actor_IMAS(equilibrium,core_profiles,wall,ece_in,ece_out, &
       codeparam_ecrad,output_flag,output_message)
  
  !call file2buffer("input/code_params.xml", io_unit, buffer)
  !call pre_initialize_ECRad_IMAS(buffer, wall, error_flag, error_message)
  !call set_ece_ECRad_IMAS(ece_in, 1, error_flag, error_message)
  !call initialize_ECRad_IMAS(equilibrium, 1, error_flag, error_message)
  !call make_rays_ECRad_IMAS(core_profiles, equilibrium, 1)
  !if(error_flag.eq.0) call make_dat_model_ECRad(core_profiles,1, ece_out)
     
  ! EXPORT RESULTS TO LOCAL DATABASE
  if(output_flag.eq.0) then
     write(*,*) '=> Export ece IDS to local database'
     call imas_create_env('ids',shot_scenario,run_scenario,0,0,idx,local_user,local_db,'3')
     call ids_put(idx,'ece', ece_out)
     call imas_close(idx)
     write(*,*) 'Done exporting.'
     write(*,*) ' '
     write(*,*) 'End of standalone'
  else
     write(*,*) output_message
     write(*,*) '=> Program stopped.'

  endif

  ! RELEASING THE MEMORY ALLOCATED FOR EACH IDS
  call ids_deallocate(wall)
  call ids_deallocate(core_profiles)
  call ids_deallocate(ece_in)
  call ids_deallocate(ece_out)

end program ECRad_IMAS
