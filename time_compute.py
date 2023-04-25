import numpy as np

# FUNCTION TO FIND INDEX OF NEAREST TIME SLICE IN TIME ARRAY
def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx],idx

def time_compute(input,time_slice,ntimes,dt_required):

    # READ FULL TIME VECTOR OF EQUILIBRIUM IDS TO GET THE TIME BASE
    print('=> Open input datafile and to read time array from equilibrium IDS')
    time_array = input.partial_get(ids_name='equilibrium',data_path='time')
    if time_slice == 0:
        ntimes = len(time_array)
        it = 0
        print('-----------------------------------------------------')
        print('Use the whole time array of the input equilibrium IDS')
        print('- t_init = ',time_array[0],'s')
        print('- ntimes = ',ntimes)
        print('-----------------------------------------------------')
    elif time_slice < 0:
        it = 0
        time_slice = time_array[0]
        if ntimes > len(time_array):
          time_array = np.linspace(time_array[0],time_array[0]+ntimes*dt_required,ntimes)
          print('-------------------------------------------------------------------------------------')
          print('Use the first time slice of the input equilibrium IDS with ntimes defined by the user')
          print('... However ntimes is greater than the input time array')
          print('--> Artificial time array generated with dt defined by the user:')
          print('    - t_init = ',time_array[it],'s')
          print('    - ntimes = ',ntimes)
          print('    - dt_required = ',dt_required)
          print('-------------------------------------------------------------------------------------')
        else:
            print('--------------------------------------------------------------------------------------')
            print('Use the first time slice of the input equilibrium IDS with ntimes defined by the user:')
            print('- t_init = ',time_array[0],'s')
            print('- ntimes = ',ntimes)
            print('--------------------------------------------------------------------------------------')
    else:
      # FIND INDEX OF NEAREST TIME SLICE IN TIME ARRAY
      [tc,it] = find_nearest(time_array,time_slice)
      if ntimes > len(time_array) or it+ntimes>len(time_array):
          time_array = np.linspace(time_array[it],time_array[it]+ntimes*dt_required,ntimes)
          it = 0
          print('----------------------------------------------------------------')
          print('Required ntimes greater than the input time array')
          print('or required t_init too large for the required ntimes')
          print('--> Artificial time array generated with dt defined by the user:')
          print('    - Required t_init = ',time_slice,'s')
          print('    - Actual t_init = ',time_array[it],'s')
          print('    - ntimes = ',ntimes)
          print('    - dt_required = ',dt_required)
          print('----------------------------------------------------------------')
      else:
          print('--------------------------------------------------------------------')
          print('Use the specificed first time slice with ntimes defined by the user:')
          print('- Required t_init = ',time_slice,'s')
          print('- Actual t_init = ',time_array[it],'s')
          print('- ntimes = ',ntimes)
          print('--------------------------------------------------------------------')

    return time_array,it
