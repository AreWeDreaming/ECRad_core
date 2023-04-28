module mod_MUSCLE3_helper

   public :: get_ids, &
             send_error, &
             send_message, &
             get_task

   private :: get_wall_ids, &
              get_equilibrium_ids, &
              get_core_profiles_ids, &
              get_ece_ids
   
   interface get_ids
      module procedure get_wall_ids, get_equilibrium_ids, &
                       get_core_profiles_ids, get_ece_ids
   end interface get_ids
   

   contains



subroutine send_error(instance, error_message)
   use ymmsl
   use libmuscle
   implicit None
   type(LIBMUSCLE_Instance), intent(inout) :: instance
   character(*), intent(in) :: error_message
   type(LIBMUSCLE_Message) :: msg
   type(LIBMUSCLE_Data) :: data_intent_out
   real(kind=8)  :: real_time
   data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
   call LIBMUSCLE_Data_set_item(data_intent_out, 1, error_message)
   call send_message(instance, data_intent_out)
   call LIBMUSCLE_Message_free(msg)
   call LIBMUSCLE_Data_free(data_intent_out)
end subroutine send_error

subroutine send_message(instance, send_data)
   use ymmsl
   use libmuscle
   implicit None
   type(LIBMUSCLE_Instance), intent(inout) :: instance
   type(LIBMUSCLE_Data), intent(inout) :: send_data
   type(LIBMUSCLE_Message) :: msg
   real(kind=8)  :: real_time
   call cpu_time(real_time)
   msg = LIBMUSCLE_Message_create(real_time, send_data)
   call LIBMUSCLE_Instance_send(instance, 'ECRad_report', msg)
   call LIBMUSCLE_Message_free(msg)
   call LIBMUSCLE_Data_free(send_data)
end subroutine send_message

function get_task(instance, data_intent_in)
   use ymmsl
   use libmuscle
   implicit None
   type(LIBMUSCLE_Instance), intent(inout) :: instance
   type(LIBMUSCLE_DataConstRef), intent(in)   :: data_intent_in
   type(LIBMUSCLE_DataConstRef)               :: arg_intent_in
   character(len=:), allocatable  :: get_task
   arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 1)
   if (.not. LIBMUSCLE_DataConstRef_is_a_character(arg_intent_in)) then
      call send_error(instance, "First argument must always be a string")
      get_task = "Error"
      return
   end if
   get_task = LIBMUSCLE_DataConstRef_as_string(arg_intent_in)
end function get_task

subroutine get_wall_ids(data_to_decode, wall)
   use ymmsl
   use libmuscle
   use ids_schemas, only: ids_wall
   use ids_routines, only: ids_deallocate
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(in) :: data_to_decode
   type(ids_wall), intent(inout) :: wall
   character(len=1), dimension(:), allocatable :: serialized_ids
   call ids_deallocate(wall)
   allocate(serialized_ids(LIBMUSCLE_DataConstRef_size(data_to_decode)))
   call LIBMUSCLE_DataConstRef_as_byte_array(data_to_decode, serialized_ids)
   call LIBMUSCLE_DataConstRef_free(data_to_decode)
   call ids_deserialize(serialized_ids, wall)
   deallocate(serialized_ids)
end subroutine get_wall_ids

subroutine get_equilibrium_ids(data_to_decode, equilibrium)
   use ymmsl
   use libmuscle
   use ids_schemas, only: ids_equilibrium
   use ids_routines, only: ids_deallocate
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(in) :: data_to_decode
   type(ids_equilibrium), intent(inout) :: equilibrium
   character(len=1), dimension(:), allocatable :: serialized_ids
   call ids_deallocate(equilibrium)
   allocate(serialized_ids(LIBMUSCLE_DataConstRef_size(data_to_decode)))
   call LIBMUSCLE_DataConstRef_as_byte_array(data_to_decode, serialized_ids)
   call LIBMUSCLE_DataConstRef_free(data_to_decode)
   call ids_deserialize(serialized_ids, equilibrium)
   deallocate(serialized_ids)
end subroutine get_equilibrium_ids

subroutine get_core_profiles_ids(data_to_decode, core_profiles)
   use ymmsl
   use libmuscle
   use ids_schemas, only: ids_core_profiles
   use ids_routines, only: ids_deallocate
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(in) :: data_to_decode
   type(ids_core_profiles), intent(inout) :: core_profiles
   character(len=1), dimension(:), allocatable :: serialized_ids
   call ids_deallocate(core_profiles)
   allocate(serialized_ids(LIBMUSCLE_DataConstRef_size(data_to_decode)))
   call LIBMUSCLE_DataConstRef_as_byte_array(data_to_decode, serialized_ids)
   call LIBMUSCLE_DataConstRef_free(data_to_decode)
   call ids_deserialize(serialized_ids, core_profiles)
   deallocate(serialized_ids)
end subroutine get_core_profiles_ids

subroutine get_ece_ids(data_to_decode, ece)
   use ymmsl
   use libmuscle
   use ids_schemas, only: ids_ece
   use ids_routines, only: ids_deallocate
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(in) :: data_to_decode
   type(ids_ece), intent(inout) :: ece
   character(len=1), dimension(:), allocatable :: serialized_ids
   call ids_deallocate(ece)
   allocate(serialized_ids(LIBMUSCLE_DataConstRef_size(data_to_decode)))
   call LIBMUSCLE_DataConstRef_as_byte_array(data_to_decode, serialized_ids)
   call LIBMUSCLE_DataConstRef_free(data_to_decode)
   call ids_deserialize(serialized_ids, ece)
   deallocate(serialized_ids)
end subroutine get_ece_ids

end module mod_MUSCLE3_helper


program ECRad_MUSCLE3
#ifdef OMP
  use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
       set_ECRad_thread_count, initialize_ECRad_IMAS, &
       reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#else
  use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
       initialize_ECRad_IMAS, &
       reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#endif
  use ids_schemas, only: ids_equilibrium, ids_core_profiles, &
                         ids_ece, ids_wall, ids_parameters_input
  use ids_routines, only: ids_deallocate
  use f90_file_reader, only: file2buffer
  use mod_MUSCLE3_helper, only: get_ids, send_error, &
                                send_message, get_task
  use ymmsl
  use libmuscle
implicit None

! ---------------------------------------
! PURPOSE: SIMPLE PSEUDO PHYSICS CODE 
! (IT READS AN EQUILIBRIUM IDS,
! MODIFIES ONE VARIABLE IN THE IDS, 
! AND WRITES THE NEW EQUILIBRIUM IDS)
! USING MUSCLE MANAGER
! ---------------------------------------

type(ids_wall):: wall
type(ids_equilibrium):: equilibrium
type(ids_core_profiles):: core_profiles
type(ids_ece):: ece_in, ece_out
type(ids_parameters_input):: codeparam
character(len=:), allocatable :: task
type(LIBMUSCLE_PortsDescription) :: ports
type(LIBMUSCLE_Instance) :: instance
type(LIBMUSCLE_Message) :: msg_in, msg_out
type(LIBMUSCLE_DataConstRef) :: data_intent_in, arg_intent_in
type(LIBMUSCLE_Data) :: data_intent_out
character(len=1), dimension(:), allocatable :: serialized_ids
logical                                     :: init_success, time_point_set
integer:: error_flag,iounit, itime_equilibrium, itime_core_profiles, error_state_out
character(len=:), pointer:: error_message
character(len=200):: xml_path
character(len=132), pointer :: codeparam_string
! Create ports for connection to the outside world
ports = LIBMUSCLE_PortsDescription_create()

! We need initial state, and I/O for intermediate states
! Initial state (wall, equilibrium, ece hardware)
call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_F_INIT,'ECRad_init')
! Task, i.e. switch time points, do raytracing, or radiationt ransport
call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_S, 'ECRad_task')
! Result of the task, i.e. sucess, or success + output ece IDS
call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_O_I, 'ECRad_report')
! Set up the ports
instance = LIBMUSCLE_Instance_create(ports)
! Deallocate port description
call LIBMUSCLE_PortsDescription_free(ports)

! Set up code parameters here. This is managed by IDA itself
xml_path = LIBMUSCLE_Instance_get_setting_as_character(instance, 'xml_path')
call file2buffer(xml_path, iounit, codeparam_string)
! Unless we critically fail to initialize we might be able to recover
init_success = .false.
time_point_set = .false.
do while (LIBMUSCLE_Instance_reuse_instance(instance))
   if(.not. init_success) then
      msg_in = LIBMUSCLE_Instance_receive(instance, 'ECRad_init')
      data_intent_in = LIBMUSCLE_Message_get_data(msg_in)
      if(LIBMUSCLE_DataConstRef_size(data_intent_in) /= 4) then
         call send_error(instance, "Wrong amount of arguments for INIT")
         cycle
      end if
      task = get_task(instance, data_intent_in)
      if(task == "Error") cycle
      if(task /= "INIT") then
         call send_error(instance, "First task must be INIT not " // trim(task))
         cycle
      end if
      ! This has to be the wall ids
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 2)
      call get_ids(arg_intent_in, wall)
      ! This has to be the equilibrium ids
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 3)
      call get_ids(arg_intent_in, equilibrium)
      ! This has to be the ece ids 
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 4)
      call get_ids(arg_intent_in, ece_in)
      ! Finally we have all the information to init ECRad and move beyond the initialization phase
      call pre_initialize_ECRad_IMAS(codeparam_string, wall, &
                                     error_flag, error_message)
      if(error_flag /= 0) then
         call send_error(instance, error_message)
         cycle
      end if
      call set_ece_ECRad_IMAS(ece_in, 1, error_flag, error_message)
      if(error_flag /= 0) then
         call send_error(instance, error_message)
      else
         data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
         call LIBMUSCLE_Data_set_item(data_intent_out, 1, "Init success")
         call send_message(instance, data_intent_out)
         call LIBMUSCLE_Data_free(data_intent_out)
         init_success = .true.
      end if
      cycle
   end if
   msg_in = LIBMUSCLE_Instance_receive(instance, 'ECRad_task')
   data_intent_in = LIBMUSCLE_Message_get_data(msg_in)
   task = get_task(instance, data_intent_in)
   if(task == "Timepoint") then
      call reset_ECRad()
      call set_ece_ECRad_IMAS(ece_in, 1, error_flag, error_message)
      if(LIBMUSCLE_DataConstRef_size(data_intent_in) /= 5) then
         call send_error(instance, "Wrong amount of arguments for Timepoint")
            cycle
      end if
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 2)
      itime_equilibrium = int(LIBMUSCLE_as_int8(arg_intent_in),4)
      call initialize_ECRad_IMAS(equilibrium, itime, output_flag, output_message)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 3)
      if(associated(core_profiles)) call ids_deallocate(core_profiles)
      call get_ids(arg_intent_in, core_profiles)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 4)
      itime_core_profiles = int(LIBMUSCLE_as_int8(arg_intent_in),4)
      call make_rays_ECRad_IMAS(core_profiles, itime)
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      call LIBMUSCLE_Data_set_item(data_intent_out, 1, "Timepoint success")
      call send_message(instance, data_intent_out)
      call LIBMUSCLE_Data_free(data_intent_out)
   else if(.not. time_point_set) then
      call send_error(instance, "Need to set time point first before further execution")
         cycle
   else if(task == "Run") then
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, 2)
      call get_ids(arg_intent_in, core_profiles)
      call make_dat_model_ECRad(core_profiles, itime_core_profiles, ece_out)
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      call LIBMUSCLE_Data_set_item(data_intent_out, 1, "Run success")
      call ids_serialize(ece_out, serialized_ids)
      call LIBMUSCLE_Data_set_item(data_intent_out, 2, serialized_ids)
      call send_message(instance, data_intent_out)
      call LIBMUSCLE_Data_free(data_intent_out)
   else
      call send_error(instance, "Follow-up tasks must be either 'Timepoint' or 'Run' not " // trim(task))
   end if
enddo
call ids_deallocate(wall)
call ids_deallocate(equilibrium)
call ids_deallocate(core_profiles)
call ids_deallocate(ece_in)
call ids_deallocate(ece_out)
end program ECRad_MUSCLE3