module mod_ECRad_actor_IMAS

contains

  subroutine ECRad_actor_IMAS(equilibrium,core_profiles,wall,ece_in,ece_out, &
       codeparam_ecrad,output_flag,output_message)

    use ids_schemas, only: ids_equilibrium, ids_core_profiles, &
         ids_ece, ids_wall, ids_parameters_input
    use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
         initialize_ECRad_IMAS, make_rays_ECRad_IMAS, make_dat_model_ECRad

    implicit none

    character(len=:), pointer:: output_message
    integer:: output_flag, io_unit = 1
    type(ids_wall):: wall
    type(ids_equilibrium):: equilibrium
    type(ids_core_profiles):: core_profiles
    type(ids_ece):: ece_in, ece_out
    type(ids_parameters_input):: codeparam_ecrad

    call pre_initialize_ECRad_IMAS(codeparam_ecrad%parameters_value, wall, output_flag, output_message)
    call set_ece_ECRad_IMAS(ece_in, 1, output_flag, output_message)
    call initialize_ECRad_IMAS(equilibrium, 1, output_flag, output_message)
    call make_rays_ECRad_IMAS(core_profiles, 1)
    if(output_flag.eq.0) call make_dat_model_ECRad(core_profiles,1, ece_out)

  end subroutine ECRad_actor_IMAS

end module mod_ECRad_actor_IMAS


