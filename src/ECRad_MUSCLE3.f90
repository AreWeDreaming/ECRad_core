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
   type(LIBMUSCLE_Data) :: data_intent_out
   real(kind=8)  :: real_time
   data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
   write(96,*), "Encountered an Error: " // error_message
   flush(96)
   call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), error_message)
   call send_message(instance, data_intent_out)
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
   arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(1, LIBMUSCLE_size))
   if (.not. LIBMUSCLE_DataConstRef_is_a_character(arg_intent_in)) then
      call send_error(instance, "First argument must always be a string")
      get_task = "Error"
      return
   end if
   get_task = LIBMUSCLE_DataConstRef_as_character(arg_intent_in)
   write(96,*), "Received task " // get_task
   flush(96)
end function get_task

subroutine get_wall_ids(data_to_decode, wall)
   use ymmsl
   use libmuscle
   use ids_schemas, only: ids_wall
   use ids_routines, only: ids_deallocate, ids_deserialize
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(inout) :: data_to_decode
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
   use ids_routines, only: ids_deallocate, ids_deserialize
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(inout) :: data_to_decode
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
   use ids_routines, only: ids_deallocate, ids_deserialize
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(inout) :: data_to_decode
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
   use ids_routines, only: ids_deallocate, ids_deserialize
   implicit None
   type(LIBMUSCLE_DataConstRef), intent(inout) :: data_to_decode
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
      allocate_ece_ids, set_ECRad_thread_count, initialize_ECRad_IMAS, &
      reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#else
  use mod_ECRad_IMAS, only : pre_initialize_ECRad_IMAS, set_ece_ECRad_IMAS, &
      allocate_ece_ids, initialize_ECRad_IMAS, &
      reset_ECRad, make_rays_ECRad_IMAS, make_dat_model_ECRad
#endif
  use ids_schemas, only: ids_equilibrium, ids_core_profiles, &
                         ids_ece, ids_wall, ids_parameters_input
  use ids_routines, only: ids_deallocate, ids_serialize
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
type(LIBMUSCLE_Data) :: data_intent_out, data_segment_intent_out
type(ids_parameters_input):: codeparam_ecrad
character(len=1), dimension(:), allocatable :: serialized_ids
logical                                     :: init_success, time_point_set, verbose=.false.
integer:: io_unit = 1, error_flag,iounit, itime_equilibrium, itime_core_profiles, error_state_out
character(len=:), pointer:: error_message
character(len=200):: xml_path
if(verbose) then
   open(96, file = "ECRad_worker.log")
   write(96,*), "Starting ECRad"
   flush(96)
   ! Create ports for connection to the outside world
   write(96,*), "Setting up ports"
   flush(96)
end if
ports = LIBMUSCLE_PortsDescription_create()

! We need initial state, and I/O for intermediate states
! Task, i.e. switch time points, do raytracing, or radiationt ransport
call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_F_INIT, 'ECRad_task')
! Result of the task, i.e. sucess, or success + output ece IDS
call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_O_F, 'ECRad_report')
if(verbose) write(96,*), "Setting up instance"
if(verbose) flush(96)
! Set up the ports
instance = LIBMUSCLE_Instance_create(ports)
! Deallocate port description
call LIBMUSCLE_PortsDescription_free(ports)

! Set up code parameters here. This is managed by IDA itself
! Unless we critically fail to initialize we might be able to recover
init_success = .false.
time_point_set = .false.
if(verbose) write(96,*), "Starting loop"
if(verbose) flush(96)
do while (LIBMUSCLE_Instance_reuse_instance(instance))
   if(verbose) write(96,*), "Iterating loop"
   if(verbose) flush(96)
   if(.not. init_success) then
      if(verbose) write(96,*), "Getting XML path data"
      if(verbose) flush(96)
      xml_path = LIBMUSCLE_Instance_get_setting_as_character(instance, 'xml_path')
      if(verbose) write(96,*), "Loading xml at" // xml_path
      if(verbose) flush(96)
      call file2buffer(xml_path, io_unit, codeparam_ecrad%parameters_value)
      if(verbose) write(96,*), "Waiting for init message"
      if(verbose) flush(96)
      msg_in = LIBMUSCLE_Instance_receive(instance, 'ECRad_task')
      if(verbose) write(96,*), "Received task"
      if(verbose) flush(96)
      data_intent_in = LIBMUSCLE_Message_get_data(msg_in)
      if(LIBMUSCLE_DataConstRef_size(data_intent_in) /= 4) then
         call send_error(instance, "Wrong amount of arguments for INIT")
         cycle
      end if
      task = get_task(instance, data_intent_in)
      if(trim(task) == "Error") cycle
      if(trim(task) /= "INIT") then
         call send_error(instance, "First task must be INIT not " // trim(task))
         cycle
      end if
      if(verbose) write(96,*), "Got call for INIT"
      if(verbose) flush(96)
      ! This has to be the wall ids
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(2, LIBMUSCLE_size))
      call get_ids(arg_intent_in, wall)
      ! This has to be the equilibrium ids
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(3, LIBMUSCLE_size))
      call get_ids(arg_intent_in, equilibrium)
      ! This has to be the ece ids 
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(4, LIBMUSCLE_size))
      call get_ids(arg_intent_in, ece_in)
      ! Finally we have all the information to init ECRad and move beyond the initialization phase
      call pre_initialize_ECRad_IMAS(codeparam_ecrad%parameters_value, wall, &
                                     error_flag, error_message)
      if(error_flag /= 0) then
         call send_error(instance, error_message)
         cycle
      end if
      call set_ece_ECRad_IMAS(ece_in, 1, error_flag, error_message)
      call allocate_ece_ids(size(ece_in%channel), ece_out)
      if(error_flag /= 0) then
         call send_error(instance, error_message)
      else
         data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
         if(verbose) write(96,*), "Finished INIT successfully"
         if(verbose) flush(96)
         call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), "Init success")
         call send_message(instance, data_intent_out)
         init_success = .true.
      end if
      cycle
   end if
   if(verbose) write(96,*), "Waiting for new task after INIT"
   if(verbose) flush(96)
   msg_in = LIBMUSCLE_Instance_receive(instance, 'ECRad_task')
   data_intent_in = LIBMUSCLE_Message_get_data(msg_in)
   task = get_task(instance, data_intent_in)
   if(trim(task) == "Timepoint") then
      if(verbose) write(96,*), "Starting work on timepoint"
      if(verbose) flush(96)
      if(LIBMUSCLE_DataConstRef_size(data_intent_in) /= 4) then
         call send_error(instance, "Wrong amount of arguments for Timepoint")
            cycle
      end if
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(2, LIBMUSCLE_size))
      itime_equilibrium = int(LIBMUSCLE_DataConstRef_as_int(arg_intent_in),4)
      call initialize_ECRad_IMAS(equilibrium, itime_equilibrium, error_flag, error_message, reinitialize = time_point_set)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(3, LIBMUSCLE_size))
      call ids_deallocate(core_profiles)
      call get_ids(arg_intent_in, core_profiles)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(4, LIBMUSCLE_size))
      itime_core_profiles = int(LIBMUSCLE_DataConstRef_as_int(arg_intent_in),4)
      call make_rays_ECRad_IMAS(core_profiles, itime_core_profiles, ece_out)
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      if(verbose) write(96,*), "Finished work on timepoint"
      if(verbose) flush(96)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), "Timepoint success")
      call ids_serialize(ece_out, serialized_ids)
      data_segment_intent_out = LIBMUSCLE_Data_create_byte_array(serialized_ids)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(2, LIBMUSCLE_size), data_segment_intent_out)
      call send_message(instance, data_intent_out)
      deallocate(serialized_ids)
      time_point_set = .true.
   else if(.not. time_point_set) then
      call send_error(instance, "Need to set time point first before further execution")
         cycle
   else if(trim(task) == "Retrace") then
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(2, LIBMUSCLE_size))
      call ids_deallocate(core_profiles)
      call get_ids(arg_intent_in, core_profiles)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(3, LIBMUSCLE_size))
      itime_core_profiles = int(LIBMUSCLE_DataConstRef_as_int(arg_intent_in),4)
      call make_rays_ECRad_IMAS(core_profiles, itime_core_profiles, ece_out)
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      if(verbose) write(96,*), "Finished work on Retrace"
      if(verbose) flush(96)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), "Retrace success")
      call ids_serialize(ece_out, serialized_ids)
      data_segment_intent_out = LIBMUSCLE_Data_create_byte_array(serialized_ids)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(2, LIBMUSCLE_size), data_segment_intent_out)
      call send_message(instance, data_intent_out)
      deallocate(serialized_ids)
   else if(trim(task) == "Run") then
      if(verbose) write(96,*), "Started work on run"
      if(verbose) flush(96)
      arg_intent_in = LIBMUSCLE_DataConstRef_get_item(data_intent_in, int(2, LIBMUSCLE_size))
      call get_ids(arg_intent_in, core_profiles)
      call make_dat_model_ECRad(core_profiles, itime_core_profiles, ece_out)
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), "Run success")
      call ids_serialize(ece_out, serialized_ids)
      data_segment_intent_out = LIBMUSCLE_Data_create_byte_array(serialized_ids)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(2, LIBMUSCLE_size), data_segment_intent_out)
      if(verbose) write(96,*), "Finished work on run"
      if(verbose) flush(96)
      call send_message(instance, data_intent_out)
      deallocate(serialized_ids)
   else if(trim(task) == "Exit") then
      data_intent_out = LIBMUSCLE_Data_create_nils(2_LIBMUSCLE_size)
      call LIBMUSCLE_Data_set_item(data_intent_out, int(1, LIBMUSCLE_size), "Exiting")
      if(verbose) write(96,*), "Exiting as requested"
      if(verbose) flush(96)
      call send_message(instance, data_intent_out)
      exit
   else
      call send_error(instance, "Follow-up tasks must be either 'Timepoint' or 'Run' not " // trim(task))
   end if
enddo
call ids_deallocate(wall)
call ids_deallocate(equilibrium)
call ids_deallocate(core_profiles)
call ids_deallocate(ece_in)
call ids_deallocate(ece_out)
call reset_ECRad()
if(verbose) close(96)
end program ECRad_MUSCLE3