Create time-frequency map
=========================

We initialise a 2D zero array for a time-frequency map (``tf_map``) which will be computed for each frequency-domain filter associated to each PSD segment and where the filtered time-series for each frequency channels will be stored. The number of rows corresponds to the total number of frequency channels which is defined by the ``nchans`` variable. The number of columns corresponds to the segment length in samples (i.e. the number of samples covering one segment) which is defined by the ``seg_len`` variable. ::

      # Initialise 2D zero array for time-frequency map
      tf_map = numpy.zeros((nchans, seg_len), dtype=numpy.complex128)

We also initialise a zero vector for a temporary filter bank (``tmp_filter_bank``) that will store, for a given channel, the filter's values from the original filter bank (``filter_bank``) for that channel only. The length of the temporary filter bank is equal to the length of the PSD frequency series (``fd_psd``). ::

      # Initialise 1D zero array 
      tmp_filter_bank = numpy.zeros(len(fd_psd), dtype=numpy.complex128)

We then loop over all the frequency channels. While in the loop, we first re-initialise the temporary filter bank with zero values everywhere along the frequency series. We then determine the first and last frequency of each channel and re-define the values of the filter in that frequency range based on the values from the original channel's filter from the original filter bank. ::

      # Loop over all the channels
      print tprint(t0,t1),"Filtering all %d channels..." % nchans
      for i in range(nchans):
          # Reset filter bank series
          tmp_filter_bank *= 0.0
          # Index of starting frequency
          f1 = int(filter_bank[i].f0/fd_psd.delta_f)
          # Index of ending frequency
          f2 = int((filter_bank[i].f0 + 2*band)/fd_psd.delta_f)+1
          # (FIXME: Why is there a factor of 2 here?)
          tmp_filter_bank[f1:f2] = filter_bank[i].data.data * 2

We then extract the frequency series from the filter bank for that channel, which will be used as a template waveform to filter the actual data from the channel. ::

          # Define the template to filter the frequency series with
          template = types.FrequencySeries(tmp_filter_bank, delta_f=fd_psd.delta_f, copy=False)

Finally, we use the `matched_filter_core <http://ligo-cbc.github.io/pycbc/latest/html/pycbc.filter.html>`_ module from the ``pycbc.filter.matchedfilter`` library to filter the frequency series from the channel. This will return both a time series containing the complex signal-to-noise matched filtered against the data, and a frequency series containing the correlation vector. ::

          # Create filtered series
          filtered_series = filter.matched_filter_core(template,fs_data,h_norm=None,psd=None,
                                                       low_frequency_cutoff=filter_bank[i].f0,
                                                       high_frequency_cutoff=filter_bank[i].f0+2*band)

The `matched filter <http://en.wikipedia.org/wiki/Matched_filter>`_ is the optimal linear filter for maximizing the signal to noise ratio (SNR) in the presence of additive stochastic noise. The filtered time series is stored in the time-frequency map and can be used to produce a spectrogram of the segment of data being analysed. ::

          # Include filtered series in the map
          tf_map[i,:] = filtered_series[0].numpy()

The time-frequency map is a 2D array with a length that corresponds to the number of channels and a width equal to the number of sample present in one segment of data, i.e. segment's length in seconds times the the sampling rate. The map can finally be plotted with a :math:`\Delta t` corresponding to the sampling period of the original dataset (i.e. inverse of the original sampling rate), and :math:`\Delta f` is equal to the bandwidth of one channel. ::

      plot_spectrogram(numpy.abs(tf_map).T,tmp_ts_data.delta_t,fd_psd.delta_f,ts_data.sample_rate,start_time,end_time,fname='%s/tf.png'%(segfolder))
