module mod_codeparam_standalone_IMAS

  implicit none

  type type_codeparam_standalone
     integer:: shot_scenario,run_scenario
     integer:: shot_wall,run_wall
     integer:: shot_ece,run_ece
     integer:: run_out
     character(len=200):: user_scenario,db_scenario
     character(len=200):: user_wall,db_wall
     character(len=200):: user_ece,db_ece
     character(len=200):: local_db
     double precision:: time_slice
  end type type_codeparam_standalone

contains

  subroutine parse_standalone_codeparam(codeparam,codeparam_data)

    use xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, &
         type_xml2eg_document, xml2eg_free_doc, set_verbose

    ! Input/Output
    character(len=132), pointer :: codeparam(:)
    type(type_codeparam_standalone), intent(out) :: codeparam_data

    ! Internal
    type(type_xml2eg_document) :: doc

    ! Parse the "codeparam" string. This means that the data is put into a document "doc"
    call xml2eg_parse_memory(codeparam,doc)
    call set_verbose(.TRUE.) ! Only needed if you want to see what's going on in the parsing
    
    call xml2eg_get(doc,'shot_scenario',codeparam_data%shot_scenario)
    call xml2eg_get(doc,'run_scenario',codeparam_data%run_scenario)
    call xml2eg_get(doc,'user_scenario',codeparam_data%user_scenario)
    call xml2eg_get(doc,'db_scenario',codeparam_data%db_scenario)

    call xml2eg_get(doc,'shot_wall',codeparam_data%shot_wall)
    call xml2eg_get(doc,'run_wall',codeparam_data%run_wall)
    call xml2eg_get(doc,'user_wall',codeparam_data%user_wall)
    call xml2eg_get(doc,'db_wall',codeparam_data%db_wall)

    call xml2eg_get(doc,'shot_ece',codeparam_data%shot_ece)
    call xml2eg_get(doc,'run_ece',codeparam_data%run_ece)
    call xml2eg_get(doc,'user_ece',codeparam_data%user_ece)
    call xml2eg_get(doc,'db_ece',codeparam_data%db_ece)
    
    call xml2eg_get(doc,'run_out',codeparam_data%run_out)
    call xml2eg_get(doc,'local_db',codeparam_data%local_db)
    call xml2eg_get(doc,'time_slice',codeparam_data%time_slice)

    ! Make sure to clean up after you!!
    call xml2eg_free_doc(doc)

    return
  end subroutine parse_standalone_codeparam

end module mod_codeparam_standalone_IMAS
