Initialise event list and determine stride boundaries
=====================================================

First of all, we create a table similar than the one made by the LIGO Scientific Collaboration (LSC) where all the information will be stored. Such table is commonly know as ``lsctables``. A pre-defined LSC table can be constructed using ``New`` function from the `glue.ligolw.lsctables <https://github.com/ligo-cbc/pycbc-glue/blob/master/glue/ligolw/lsctables.py>`_ module. We use the ``SnglBurstTable`` function for the type of data to be stored and define all the columns we wish to record. ::

  # Create event list for single burst table
  event_list = lsctables.New(lsctables.SnglBurstTable,
                             ['start_time','start_time_ns','peak_time','peak_time_ns',
                              'duration','bandwidth','central_freq','chisq_dof',
                              'confidence','snr','amplitude','channel','ifo',
                              'process_id','event_id','search','stop_time','stop_time_ns'])

We also need to determine the indexes of both starting and ending times for the first segment to analyse, respectively ``t_idx_min`` and ``t_idx_max``. The default values are considered to be 0 for the starting index and the segment length in sample unit for the ending time index. Also, if the user defines a different starting time than the one from the loaded data, the offset index in sample unit is determined and added the both starting and ending time indexes. ::

  # Determine boundaries of stride in time domain
  t_idx_min, t_idx_max = 0, seg_len
  # Check if user requested starting time is defined
  if args.analysis_start_time is not None:
      # Define the time difference in seconds between data and user requested starting times
      t_idx_off = args.analysis_start_time - ts_data.start_time
      # Calculate the index of the user requested starting point in the data
      t_idx_off = int(t_idx_off * args.sample_rate)
  else:
      # Define index of the starting point as first value in data
      t_idx_off = 0
  # Initialise minimum index values as offset starting index
  t_idx_min += t_idx_off
  # Initialise maximum index values as offset starting index
  t_idx_max += t_idx_off

Finally, the index for the ending time after all the segments have been analysed can be estimated for the user-defined parameter or is defined as the length of the time series data ``ts_data``. ::

  # Check if user requested end time is defined
  if args.analysis_end_time is not None:
      # Define the time difference between data and user requested ending times
      t_idx_max_off = args.analysis_end_time - ts_data.start_time
      # Calculate the index of the user requested starting point in the data
      t_idx_max_off = int(t_idx_max_off * args.sample_rate)
  else:
      # Define index of the ending point as the length of data array
      t_idx_max_off = len(ts_data)
