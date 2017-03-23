.. _analysingblocks:

Define analysing blocks
-----------------------

The first thing we do is to calculate the time series for the segment that is covered (``tmp_ts_data``) and redefined the metadata, especially the time of the first sample in seconds which is defined by the ``epoch`` argument and is different for every segment. After plotting the time series for that segment, the data are then converted into frequency series (``fs_data``) using the `to_frequencyseries <http://ligo-cbc.github.io/pycbc/latest/html/pycbc.types.html#pycbc.types.timeseries.TimeSeries.to_frequencyseries>`_ module from the ``pycbc.types.timeseries.TimeSeries`` library. Finally, the frequency data are then whitened. ::

  # Loop over each data within the user requested time period
  while t_idx_max <= t_idx_max_off:
      # Define starting and ending time of the segment in seconds
      start_time = ts_data.start_time + t_idx_min/float(args.sample_rate)
      end_time = ts_data.start_time + t_idx_max/float(args.sample_rate)
      print tprint(t0,t1),"Analyzing block %i to %i (%.2f percent)"%(start_time,end_time,100*float(t_idx_max)/float(idx_max_off))
      # Model a withen time series for the block
      tmp_ts_data = types.TimeSeries(ts_data[t_idx_min:t_idx_max]*window, 1.0/args.sample_rate,epoch=start_time)
      # Save time series in segment repository
      segfolder = 'segments/%i-%i'%(start_time,end_time)
      os.system('mkdir -p '+segfolder)
      plot_ts(tmp_ts_data,fname='%s/ts.png'%(segfolder))
      # Convert times series to frequency series
      fs_data = tmp_ts_data.to_frequencyseries()
      print tprint(t0,t1),"Frequency series data has variance: %s" % fs_data.data.std()**2
      # Whitening (FIXME: Whiten the filters, not the data)
      fs_data.data /= numpy.sqrt(fd_psd) / numpy.sqrt(2 * fd_psd.delta_f)
      print tprint(t0,t1),"Whitened frequency series data has variance: %s" % fs_data.data.std()**2
