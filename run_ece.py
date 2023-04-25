#!/usr/bin/env python

# -------------------------------
# PYTHON SCRIPT TO EXECUTE ECRAD
# -------------------------------

# NEEDED MODULES
import os,imas,sys,copy,yaml,argparse
from scipy import interpolate
from time_compute import time_compute
import numpy as np

# MANAGEMENT OF INPUT ARGUMENTS
parser = argparse.ArgumentParser(description=\
                   '---- Run ECRad according to parameters from input configuration folder')
parser.add_argument('-c','--config_folder',help='Configuration folder',required=False,type=str)

# CONFIGURATION FOLDER
args = vars(parser.parse_args())
config_folder = args["config_folder"]
if config_folder is None:
    config_folder = 'input'

scenario = config_folder+'/scenario.yaml'
codeparam_ece = config_folder+'/code_params.xml'
codeparam_chease = config_folder+'/input_chease.xml'

# ----------------------------------------------------------------------
  
# INPUT/OUTPUT CONFIGURATION
file = open(scenario, 'r')    

# INPUT/OUTPUT CONFIGURATION
cc = yaml.load(file,Loader=yaml.CLoader)
file.close()
if cc['output_user_or_path'] == 'default':
    cc['output_user_or_path'] = os.getenv('USER')

# TO SUPPRESS THE ACCESS LAYER WARNINGS
# (THE LAZY WAY, OTHERWISE SEE WARNOFF BELOW)
if cc['warnoff'] == 1:
    import logging
    log = logging.getLogger()
    log.setLevel(logging.ERROR)

# DISPLAY SIMULATION INFORMATION
print('---------------------------------')
print('shot_scenario   = ',cc['shot_scenario'])
print('run_scenario    = ',cc['run_scenario'])
print('user_scenario   = ',cc['user_scenario'])
print('db_scenario     = ',cc['db_scenario'])
print('----')
print('shot_wall       = ',cc['shot_wall'])
print('run_wall        = ',cc['run_wall'])
print('user_wall       = ',cc['user_wall'])
print('db_wall         = ',cc['db_wall'])
print('----')
print('shot_ece        = ',cc['shot_ece'])
print('run_ece         = ',cc['run_ece'])
print('user_ece        = ',cc['user_ece'])
print('db_ece          = ',cc['db_ece'])
print('----')
print('output_database = ',cc['output_database'])
print('run_out         = ',cc['run_out'])
print('time_slice      = ',cc['time_slice'])
print('ntimes          = ',cc['ntimes'])
print('dt_required     = ',cc['dt_required'],' s')
print('run_chease      = ',cc['run_chease'])
print('run_ecrad       = ',cc['run_ecrad'])
print('debug_chease    = ',cc['debug_chease'])
print('debug_ecrad     = ',cc['debug_ecrad'])
print('---------------------------------')

# FUNCTION TO FIND INDEX OF NEAREST TIME SLICE IN TIME ARRAY
def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx],idx

# OPEN INPUT DATAFILE TO GET DATA FROM INPUT DATABASE
print('=> Open scenario datafile')
scenario = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,cc['db_scenario'],cc['shot_scenario'],cc['run_scenario'],cc['user_scenario'])
scenario.open()

# READ FULL TIME VECTOR OF EQUILIBRIUM IDS TO GET THE TIME BASE
time_array,it = time_compute(scenario,cc['time_slice'],cc['ntimes'],cc['dt_required'])

# ECE MACHINE DESCRIPTION DATA
print('   ---> Get ece')
md_ece = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,cc['db_ece'],cc['shot_ece'],cc['run_ece'],cc['user_ece'])
md_ece.open()
input_ece = md_ece.get('ece')
md_ece.close()

# WALL MACHINE DESCRIPTION DATA
print('   ---> Get wall')
md_wall = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,cc['db_wall'],cc['shot_wall'],cc['run_wall'],cc['user_wall'])
md_wall.open()
input_wall = md_wall.get('wall')
md_wall.close()

# IF LOCAL DATABASE DOES NOT EXIST: CREATE IT
local_database = os.getenv("HOME") + "/public/imasdb/" + cc['output_database'] + "/3/0"
if os.path.isdir(local_database) == False:
    print("-- Create local database " + cc['local_database'])
    os.makedirs(cc['local_database'])

# CREATE OUTPUT DATAFILE
print('=> Create output datafile')
output = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND,\
                      cc['output_database'],cc['shot_scenario'],cc['run_out'],cc['output_user_or_path'])
output.create()

# INITIALISATIONS
if cc['run_chease'] ==1:
    from chease.actor import chease as chease_actor
    chease = chease_actor()
    code_parameters = chease.get_code_parameters()
    runtime_settings = chease.get_runtime_settings()
    if cc['debug_chease'] == 1:
        from chease.common.runtime_settings import DebugMode
        runtime_settings.debug_mode = DebugMode.STANDALONE
    code_parameters.parameters_path = codeparam_chease
    chease.initialize(code_parameters=code_parameters,runtime_settings=runtime_settings)

if cc['run_ecrad'] == 1:
    from ecrad.actor import ecrad as ecrad_actor
    ecrad = ecrad_actor()
    code_parameters = ecrad.get_code_parameters()
    runtime_settings = ecrad.get_runtime_settings()
    if cc['debug_ecrad'] == 1:
        from ecrad.common.runtime_settings import DebugMode
        runtime_settings.debug_mode = DebugMode.STANDALONE
    code_parameters.parameters_path = codeparam_ece
    ecrad.initialize(code_parameters=code_parameters,runtime_settings=runtime_settings)

# TIME LOOP TO EXECUTE THE ECRAD CODE
for itime in range(it,it+cc['ntimes']):

  if len(time_array)>1:
    print('Time = %5.2f' % time_array[itime],'s, itime = ',itime,'/',it+cc['ntimes']-1)
    time = time_array[itime]
  else:
    time = time_array[0]

  # READ IDSS FROM INPUT SCENARIO
  print('=> Read input IDSs')
  print('   ---> equilibrium')
  input_equilibrium = scenario.get_slice('equilibrium',time_array[itime],1)
  print('   ---> core_profiles')
  input_core_profiles = scenario.get_slice('core_profiles',time_array[itime],1)

  # FILL RHO_POL_NORM IF IT DOES NOT EXIST
  if len(input_core_profiles.profiles_1d[0].grid.rho_pol_norm) == 0:
      npsi = len(input_core_profiles.profiles_1d[0].grid.psi)
      psiaxis = input_core_profiles.profiles_1d[0].grid.psi[0]
      if imas.imasdef.isFieldValid(input_core_profiles.profiles_1d[0].grid.psi_boundary):
          psiedge = input_core_profiles.profiles_1d[0].grid.psi_boundary
      elif imas.imasdef.isFieldValid(input_equilibrium.time_slice[0].global_quantities.psi_boundary):
          psiedge = input_equilibrium.time_slice[0].global_quantities.psi_boundary
      else:
          print('rho_pol_norm missing and impossible to re-calculate')
          print('=> Execusion aborted.')
          exit()
      print('rho_pol_norm missing but re-calculated from core_profiles.profiles_1d[0].grid.psi')
      input_core_profiles.profiles_1d[0].grid.rho_pol_norm.resize(npsi)
      input_core_profiles.profiles_1d[0].grid.rho_pol_norm = \
          (input_core_profiles.profiles_1d[0].grid.psi-psiaxis) \
          /(psiedge-psiaxis)

  # EXECUTE CHEASE
  if cc['run_chease'] == 1:
      print('=> Execute CHEASE')
      output_equilibrium = chease(input_equilibrium)
      input_equilibrium = copy.deepcopy(output_equilibrium)

  # EXECUTE ECRAD
  if cc['run_ecrad'] == 1:
      print('=> Execute ECRAD')
      output_ece = ecrad(input_equilibrium,input_core_profiles,input_wall,input_ece)

  # WRITE THE OUTPUT TIME SLICE
  output.put_slice(input_equilibrium)
  output.put_slice(input_core_profiles)
  output.put_slice(input_wall)
  if cc['run_ecrad'] == 1:
      output.put_slice(output_ece)
  else:
      output.put_slice(input_ece)

  if cc['run_ecrad'] == 1:
      print('------------------------')
      print('Output time = %5.2f s' % (output_ece.time[0]))
      print('------------------------')

# CLOSE FILES
if cc['run_ecrad'] == 1:
    ecrad.finalize()
scenario.close()
output.close()
print('End of Python script.')




