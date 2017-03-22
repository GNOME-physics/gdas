Excess-Power algorithm
======================

General overview
----------------

The **Excess Power method** is known as the *optimal detection strategy* to search for burst signals for which only the duration and frequency band are known, which is basically the case for GNOME and its search of Axion-Like Particles (ALP). This method was developed and introduced by `Anderson et al. (200) <https://arxiv.org/pdf/gr-qc/0008066v1.pdf>`_ and has been extensively used in the detection of burst sources of gravitational radiation. A more technical documentation was written by `Brady et al. (2007) <http://www.lsc-group.phys.uwm.edu/~siemens/power.pdf>`_ describing how the algorithm used by the LIGO collaboration works and how the theory is translated into code.

We present below a step-by-step procedure followed during the Excess Power search analysis. For a better representation of what is happening, the figure at the end shows how the data is being split and analysed to search for multiple signals of different bandwidth and duration in the time-frequency plane.

- :ref:`Time domain segmentation and PSD estimate <psdestimate>`

   We first estimate the instrument's noise Power Spectral Density (PSD) by splitting the time-series data into multiple overlapping segments. A periodogram for each segment is calculated separately and then averaged, which will reduce the variance of the individual power measurements. The result is a frequency series where samples are separated in frequency space by :math:`\Delta f` equal to the inverse of a segment’s length and with a high end frequency limit equal to the Nyquist limit. The final power spectrum will help reveal the existence, or the absence, of repetitive patterns and correlation structures in a signal process.

- :ref:`Comb of frequency channels <filterbank>`

   We then split the PSD frequency series into multiple channels. For each channel, a frequency domain filter is created with a :math:`\Delta f` determined by the PSD and a total extent in Fourier space that is twice the stated bandwidth of a channel. The result is a list of each channel filter's frequency series.

- :ref:`Creating analysing blocks <analysingblocks>`

   The Excess Power method can lead to moderately-large computational requirements, and it has been found that the computational efficiency of this implementation can be improved upon by considering blocks of data that are much longer than the longest signal time duration. The entire time series is therefore split into separate blocks. We use the length of the segments used for PSD estimate to define the duration of each block. For each block, the time series is c0Aonverted into frequency series which is then filtered by the filter bank throughout all the channels. A time-frequency map is finally created which stores all the filtered frequency series from each channel.

- :ref:`Creating tiles with different bandwidth <tilebandwidth>`

   We can now construct tiles with different bandwidth by summing multiple channels together.

- :ref:`Exploring tiles with different duration <tileduration>`

   For each given tile's bandwidth, one can investigate different tile's duration. This can be done by exploring different number of degrees of freedom, :math:`d`, which can be calculated as follows: :math:`d=2BT` where :math:`B` and :math:`T` are respectively the bandwidth and duration of the tile. Section 2.2.5 of `Brady et al. <http://www.lsc-group.phys.uwm.edu/~siemens/power.pdf>`_ gives a great description of how to interpret the number of degrees of freedom. Therefore, by changing the :math:`d`, one can explore multiple tile's duration for different bandwidth.

- :ref:`Define triggering signal <triggerfinding>`

   The energy of each tile in the time-frequency space is calculated and compare to a user-defined threshold value. After defining a tile false alarm probability threshold in Gaussian noise and using the number of degrees of freedom for each tile, one can define a energy threshold value above which a burst trigger can be identified by comparing the energy threshold with the tile's energy in the time-frequency map. A tile energy time frequency map plot similar to Figure 5 in `Pustelny et al. (2013) <http://budker.berkeley.edu/~vincent/gnome/documentation/2013_pustelny.pdf>`_ can then be made which plots the outlying tile energies present in the data.

.. figure:: ./img/overview.png

	    Overview of the Excess Power method and difference between segments, channels, tiles and blocks.
  
.. _psdestimate:

Estimate Power Spectral Density (PSD)
-------------------------------------

The instrument's noise Power Spectral Density (PSD) will be used to whiten the data and help reveal the existence, or the absence, of repetitive patterns and correlation structures in the signal process. It will also determine the total bandwidth spanned by each of the filters that will subsequently be created. The first thing to do before calculating the PSD is to ensure that the time series data is converted into an array of floating values. ::

  # Convert time series as array of float
  data = ts_data.astype(numpy.float64)

The PSD is calculated by splitting up the signal into overlapping segments and scan through each segment to calculate individual periodogram. The periodograms from each segment are then averaged, reducing the variance of the individual power measurements. In order to proceed, we need to define the average method, ``avg_method``, that will be used to measure the PSD from the data. This can be specified with the ``--psd-estimation`` option. ::

  # Average method to measure PSD from the data
  avg_method = args.psd_estimation

One also needs to specify the length of each segment, ``seg_len``, as well as the separation between 2 consecutive segments, ``seg_stride``. Both parameters can be defined in second units with the ``--psd-segment-length`` and ``--psd-segment-stride`` arguments respectively and can then be converted into sample unit. ::

  # The segment length for PSD estimation in samples
  seg_len = int(args.psd_segment_length * args.sample_rate)
  # The separation between consecutive segments in samples
  seg_stride = int(args.psd_segment_stride * args.sample_rate)

We then use the `Welch's method <https://en.wikipedia.org/wiki/Welch%27s_method>`_ to perform the power spectral density estimate using the `welch <http://ligo-cbc.github.io/pycbc/latest/html/_modules/pycbc/psd/estimate.html#welch>`_ module from the ``pycbc.psd`` library. What this will do is to compute the discrete Fourier transform for each PSD segment to produce invidual periodograms, and then compute the squared magnitude of the result. The individual periodograms are then averaged using the user-defined average method, ``avg_method``, and return the frequency series, ``fd_psd``, which will store the power measurement for each frequency bin. ::

  # Lifted from the psd.from_cli module
  fd_psd = psd.welch(data,avg_method=avg_method,seg_len=seg_len,seg_stride=seg_stride)
  # Plot the power spectral density
  plot_spectrum(fd_psd)
  # We need this for the SWIG functions
  lal_psd = fd_psd.lal()

One can display the power measurements, frequency array and frequency between consecutive samples, :math:`\Delta f` in Hertz, by printing the following variables: ::

  print 'Display power measurements of the first 10 frequency bins'
  print fd_psd[:10]
  print 'Display central frequency of the first 10 bins'
  print fd_psd.sample_frequencies[:10]
  print 'Display the frequency separation between bins'
  print fd_psd.delta_f

:math:`\Delta f` corresponds to the inverse of a segment's length which is the smallest frequency (i.e. highest period) of detectable signals in each segment. The frequency range spans from 0 to the Nyquist frequency, i.e. half de the sampling rate.

Checking filtering settings
---------------------------

The first thing to check is that the frequency of the high-pass filter (if defined) is below the minimum frequency of the filter bank. Indeed, a high-pass filter will only let pass frequency that are higher than the cutoff frequency (here defined by the ``strain_high_pass`` argument). If the high pass frequency is greater from the minimum frequency in the filter bank, the signal with frequencies lower than the cutoff frequency will get attenuated. ::

  if args.min_frequency < args.strain_high_pass:
      print >>sys.stderr, "Warning: strain high pass frequency %f is greater than the tile minimum frequency %f --- this is likely to cause strange output below the bandpass frequency" % (args.strain_high_pass, args.min_frequency)

In case the maximum frequency in the filter bank is not defined, we set it to be equal to the Nyquist frequency, i.e. half the sampling rate, which makes sense as a larger signal will not be able to get easily identifiable. ::

  if args.max_frequency is None:
      args.max_frequency = args.sample_rate / 2.0

If the bandwidth of the finest filter (``--tile-bandwidth`` argument, see section :ref:`construct_args <construct_args>` or the number of frequency channels (=--channels= argument) is not defined but the total spectral band is (``data_band``), one can then determined all the filter settings as follows: ::


  if args.tile_bandwidth is None and args.channels is None:
      # Exit program with error message
      exit("Either --tile-bandwidth or --channels must be specified to set up time-frequency plane")
  else:
      # Define as assert statement that tile maximum frequency larger than its minimum frequency
      assert args.max_frequency >= args.min_frequency
      # Define spectral band of data
      data_band = args.max_frequency - args.min_frequency
      # Check if tile bandwidth or channel is defined
      if args.tile_bandwidth is not None:
          # Define number of possible filter bands
          nchans = args.channels = int(data_band / args.tile_bandwidth)  - 1
      elif args.channels is not None:
          # Define filter bandwidth 
          band = args.tile_bandwidth = data_band / (args.channels + 1)
      assert args.channels > 1

The minimum frequency to be explored can be user-defined by using the ``--min-frequency`` option. ::

  # Lowest frequency of the first filter
  flow = args.min_frequency

Whitening window and spectral correlation
-----------------------------------------

This part determines how much data on either side of the tukey window is to be discarded. Nominally, this means that one will lose ``window_fraction`` * ``args.psd_segment_length`` to corruption from the window, i.e. this is simply discarded. This is tuned to give an integer offset when used with ``args.psd_segment_length`` equal to 8, smaller windows will have fractions of integers, but larger powers of two will still preseve this (probably not a big deal in the end). ::

  window_fraction = 0

The two point spectral correlation is then done with the :ref:`calculate_spectral_correlation <calculate_spectral_correlation>` function which will return both the Tukey window applied to the original time series data and the actual two-point spectral correlation function for the whitened frequency series from the applied whitening window. ::

  # Do two point spectral correlation
  window, spec_corr = calculate_spectral_correlation(seg_len,'tukey',window_fraction=window_fraction)
  window = window.data.data
  window_sigma_sq = numpy.mean(window**2)
  # Pre scale the window by its root mean squared -- see eqn 11 of EP document
  #window /= numpy.sqrt(window_sigma_sq)

.. _filterbank:
  
Computing the filter bank
-------------------------

The filter bank will create band-pass filters for each channel in the PSD frequency domain. The :ref:`create_filter_bank <create_filter_bank>` function will san the bandwidth from the central frequency of the first channel (i.e. flow+band/2) to final frequency of the last channel (i.e. band*nchans) in a increment equal to the frequency band. The filter's total extent in Fourier space is actually twice the stated bandwidth (FWHM). ::

  # Define filters
  filter_bank, fdb = create_filter_bank(fd_psd.delta_f, flow+band/2, band, nchans, fd_psd, spec_corr)

This function will returns 2 arrays: the ``filter_bank`` array which is a list of `COMPLEX16FrequencySeries <http://software.ligo.org/docs/lalsuite/lal/struct_c_o_m_p_l_e_x16_frequency_series.html>`_ arrays corresponding to each channel's filter, and the =fdb= array which provides the time-series from each filter. The length of each array is equal to the total number of channel (i.e. =nchans=). The filter's data, :math:`\Delta f` value, and first and last frequencies of any channel's filter can be displayed as followed: ::

  # Print data of first channel's filter
  print filter_bank[0].data.data
  # Print frequency separation between 2 values in the first channel's filter
  print filter_bank[0].deltaF
  # Print first frequency of the first channel's filter
  print filter_bank[0].f0
  # Print last frequency of the first channel's filter (equal to twice the channel's bandwidth)
  print filter_bank[0].f0+(len(filter_bank[0].data.data)-1)*filter_bank[0].deltaF

Further in the analysis, the following filters will used:
1. ``white_filter_ip``: Whitened filter inner products computed with themselves.
2. ``unwhite_filter_ip``: Unwhitened filter inner products computed with themselves.
3. ``white_ss_ip``: Whitened filter inner products computed between input adjacent filters.
4. ``unwhite_ss_ip``: Unwhitened filter inner products computed between input adjacent filters.

::
     
   # This is necessary to compute the mu^2 normalizations
   white_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, None)
   unwhite_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, lal_psd)
   # These two are needed for the unwhitened mean square sum (hrss)
   white_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, None)
   unwhite_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, lal_psd)
     
Normalization of virtual channel
--------------------------------

The virtual channels will be used during the excesspower analysis to explore different frequency ranges around each PSD segments and look for possible triggers. Each channel is renormalized using the :ref:`compute_channel_renomalization <compute_channel_renomalization>` internal function. ::

  # Initialise dictionary
  mu_sq_dict = {}
  # nc_sum additional channel adds
  for nc_sum in range(0, int(math.log(nchans, 2))):
      min_band = (len(filter_bank[0].data.data)-1) * filter_bank[0].deltaF / 2
      print tprint(t0,t1),"Calculation for %d %d Hz channels" % (nc_sum+1, min_band)
      nc_sum = 2**nc_sum - 1
      mu_sq_dict[nc_sum] = compute_channel_renomalization(nc_sum, filter_bank, spec_corr, nchans)

Initialise event list and determine stride boundaries
-----------------------------------------------------

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

Create time-frequency map for each block
----------------------------------------

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

.. _tilebandwidth:

Constructing tiles of different bandwidth
-----------------------------------------

First and foremost, we define a clipping region in the data to be used to remove window corruption, this is non-zero if the ``window_fraction`` variable is set to a non-zero value. ::

      print tprint(t0,t1),"Beginning tile construction..."
      # Clip the boundaries to remove window corruption
      clip_samples = int(args.psd_segment_length * window_fraction * args.sample_rate / 2)

In order to perform a multi-resolution search, tiles of many different bandwidths and durations will be scanned. We first need to setup a loop such that the maximum number of additional channel is equal to the base 2 logarithm of the total number of channels. The number of narrow band channels to be summed (``nc_sum``) would therefore be equal to 2 to the power of the current quantity of additional channels. ::

      for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]: # nc_sum additional channel adds
          nc_sum = 2**nc_sum - 1
          print tprint(t0,t1,t2),"Summing %d narrow band channels..." % (nc_sum+1)

The undersampling rate for this tile can be calculated using the channel frequency band and the number of narrow band channels to be summed such that the bandwidth of the tile is equal to ``band * (nc_sum + 1)``. ::

          us_rate = int(round(1.0 / (2 * band*(nc_sum+1) * ts_data.delta_t)))
          print >>sys.stderr, "Undersampling rate for this level: %f" % (args.sample_rate/us_rate)

"Virtual" wide bandwidth channels are constructed by summing the samples from multiple channels, and correcting for the overlap between adjacent channel filters. We then define the normalised channel at the current level and create a time frequency map for this tile using the :ref:`make_indp_tiles <make_indp_tiles>` internal function. In other word, we are constructing multiple sub-tiles for which we can determined the respective energy in the given frequency band. ::

          mu_sq = mu_sq_dict[nc_sum]
          sys.stderr.write("\t...calculating tiles...")
          if clip_samples > 0:
              tiles = make_indp_tiles(tf_map[:,clip_samples:-clip_samples:us_rate], nc_sum, mu_sq)
          else:
              tiles = make_indp_tiles(tf_map[:,::us_rate], nc_sum, mu_sq)
          sys.stderr.write(" TF-plane is %dx%s samples... " % tiles.shape)
          print >>sys.stderr, " done"
          print "Tile energy mean: %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))

.. _tileduration:

Explore multiple tile durations
-------------------------------

Now that we create a tile with a specific bandwidth, we can start exploring different durations for the tile. We will start checking if the user manually defined a value for the longest duration tile to compute, which can be done using the ``--max-duration`` argument. If not, the value will be set to 32. ::

          if args.max_duration is not None:
              max_dof = 2 * args.max_duration * (band * (nc_sum+1))
          else:
              max_dof = 32
          assert max_dof >= 2

Since we produce (initially) tiles with 1 degree of freedom, the duration goes as one over twice the bandwidth. ::

          print "\t\t...getting longer durations..."
          #for j in [2**l for l in xrange(1, int(math.log(max_dof, 2))+1)]:
          for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
              sys.stderr.write("\t\tSumming DOF = %d ..." % (2*j))
              #tlen = tiles.shape[1] - j + 1
              tlen = tiles.shape[1] - 2*j + 1 + 1
              if tlen <= 0:
                  print >>sys.stderr, " ...not enough samples."
                  continue
              dof_tiles = numpy.zeros((tiles.shape[0], tlen))
              #:sum_filter = numpy.ones(j)
              # FIXME: This is the correct filter for 50% overlap
              sum_filter = numpy.array([1,0] * (j-1) + [1])
              #sum_filter = numpy.array([1,0] * int(math.log(j, 2)-1) + [1])
              for f in range(tiles.shape[0]):
                  # Sum and drop correlate tiles
                  # FIXME: don't drop correlated tiles
                  #output = numpy.convolve(tiles[f,:], sum_filter, 'valid')
                  dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
              print >>sys.stderr, " done"
              print "Summed tile energy mean: %f, var %f" % (numpy.mean(dof_tiles), numpy.var(dof_tiles))
              level_tdiff = time.time() - tdiff
              print >>sys.stderr, "Done with this resolution, total %f" % level_tdiff

Finally, the bandwidth and duration of the tile can be defined as followed: ::

              # Current bandwidth of the time-frequency map tiles
              current_band = band * (nc_sum + 1)
              # How much each "step" is in the frequency domain -- almost
              # assuredly the fundamental bandwidth
              df = current_band
              # How much each "step" is in the time domain -- under sampling rate
              # FIXME: THis won't work if the sample rate isn't a power of 2
              dt = 1.0 / 2 / (2 * current_band) * 2
              full_band = 250
              dt = current_band / full_band * ts_data.sample_rate
              dt = 1.0/dt
              # Duration is fixed by the NDOF and bandwidth
              duration = j / 2.0 / current_band    

.. _triggerfinding:

Trigger finding
---------------

In order to find any trigger in the data, we first need to set a false alarm probability threshold in Gaussian noise above which signal will be distinguished from the noise. Such threshold can be determined by using the /inverse survival function/ method from the `scipy.stats.chi2 <https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.stats.chi2.html>`_ package. ::

              threshold = scipy.stats.chi2.isf(args.tile_fap, j)
              print "Threshold for this level: %f" % threshold
              #if numpy.any(dof_tiles > threshold):
                  #plot_spectrogram(dof_tiles.T)
                  #import pdb; pdb.set_trace()

Once the threshold is set, one can then run the :ref:`trigger_list_from_map <trigger_list_from_map>` function to quickly find the trigger signal from the ``dof_tiles`` array that ::

              # Since we clip the data, the start time needs to be adjusted accordingly
              window_offset_epoch = fs_data.epoch + args.psd_segment_length * window_fraction / 2
              trigger_list_from_map(dof_tiles, event_list, threshold, window_offset_epoch, filter_bank[0].f0 + band/2, duration, current_band, df, dt, None)
              for event in event_list[::-1]:
                  if event.amplitude != None:
                      continue
                  etime_min_idx = float(event.get_start()) - float(fs_data.epoch)
                  etime_min_idx = int(etime_min_idx / tmp_ts_data.delta_t)
                  etime_max_idx = float(event.get_start()) - float(fs_data.epoch) + event.duration
                  etime_max_idx = int(etime_max_idx / tmp_ts_data.delta_t)
                  # (band / 2) to account for sin^2 wings from finest filters
                  flow_idx = int((event.central_freq - event.bandwidth / 2 - (band / 2) - flow) / band)
                  fhigh_idx = int((event.central_freq + event.bandwidth / 2 + (band / 2) - flow) / band)
                  # TODO: Check that the undersampling rate is always commensurate
                  # with the indexing: that is to say that
                  # mod(etime_min_idx, us_rate) == 0 always
                  z_j_b = tf_map[flow_idx:fhigh_idx,etime_min_idx:etime_max_idx:us_rate]
                  # FIXME: Deal with negative hrss^2 -- e.g. remove the event
                  try:
                      event.amplitude = measure_hrss(z_j_b, unwhite_filter_ip[flow_idx:fhigh_idx], unwhite_ss_ip[flow_idx:fhigh_idx-1], white_ss_ip[flow_idx:fhigh_idx-1], fd_psd.delta_f, tmp_ts_data.delta_t, len(filter_bank[0].data.data), event.chisq_dof)
                  except ValueError:
                      event.amplitude = 0
      
              print "Total number of events: %d" % len(event_list)

Switch to new block
-------------------

The following will move the frequency band to the next segment: ::
    
      tdiff = time.time() - tdiff
      print "Done with this block: total %f" % tdiff
      
      t_idx_min += int(seg_len * (1 - window_fraction))
      t_idx_max += int(seg_len * (1 - window_fraction))

Extracting GPS time range
-------------------------

We use the `LIGOTimeGPS <http://software.ligo.org/docs/lalsuite/lal/group___l_a_l_datatypes.html#ss_LIGOTimeGPS>`_ structure from the =glue.lal= package to /store the starting and ending time in the dataset to nanosecond precision and synchronized to the Global Positioning System time reference/. Once both times are defined, the range of value is stored in a semi-open interval using the `segment <http://software.ligo.org/docs/glue/glue.__segments.segment-class.html>`_ module from the =glue.segments= package. ::

  # Starting epoch relative to GPS starting epoch
  start_time = LIGOTimeGPS(args.analysis_start_time or args.gps_start_time)
  # Ending epoch relative to GPS ending epoch
  end_time = LIGOTimeGPS(args.analysis_end_time or args.gps_end_time)
  # Represent the range of values in the semi-open interval 
  inseg = segment(start_time,end_time)

Prepare output file for given time range
----------------------------------------

::

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  
  ifo = args.channel_name.split(":")[0]
  proc_row = register_to_xmldoc(xmldoc, __program__, args.__dict__, ifos=[ifo],version=glue.git_version.id, cvs_repository=glue.git_version.branch, cvs_entry_time=glue.git_version.date)
  
  # Figure out the data we actually analyzed
  outseg = determine_output_segment(inseg, args.psd_segment_length, args.sample_rate, window_fraction)
  
  ss = append_search_summary(xmldoc, proc_row, ifos=(station,), inseg=inseg, outseg=outseg)
  
  for sb in event_list:
      sb.process_id = proc_row.process_id
      sb.search = proc_row.program
      #sb.ifo, sb.channel = args.channel_name.split(":")
      sb.ifo, sb.channel = station, setname
  
  xmldoc.childNodes[0].appendChild(event_list)
  fname = make_filename(station, inseg)
  
  utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"), verbose=True)

Plot trigger results
--------------------

::

  events = SnglBurstTable.read(fname+'.gz')
  #del events[10000:]
  plot = events.plot('time', 'central_freq', "duration", "bandwidth", color="snr")
  #plot = events.plot('time', 'central_freq', color='snr')
  #plot.set_yscale("log")
  plot.set_ylim(1e-0, 250)
  t0 = 1153742417
  plot.set_xlim(t0 + 0*60, t0 + 1*60)
  #plot.set_xlim(t0 + 28, t0 + 32)
  pyplot.axvline(t0 + 30, color='r')
  cb = plot.add_colorbar(cmap='viridis')
  plot.savefig("triggers.png")

Module Access
=============
   
Extract Magnetic Field Data
---------------------------

Extract magnetic field data from HDF5 files.

.. currentmodule:: gdas.retrieve

.. autosummary::
   :toctree: generated/

   magfield
   file_to_segment
   construct_utc_from_metadata
   generate_timeseries
   create_activity_list
   retrieve_data_timeseries
   retrieve_channel_data
   
Plotting routines
-----------------

Methods to produce time-frequency plots and others

.. currentmodule:: gdas.plots

.. autosummary::
   :toctree: generated/

   plot_activity
   plot_time_series
   plot_asd
   plot_whitening
   plot_ts
   plot_spectrum
   plot_spectrogram
   plot_spectrogram_from_ts
   plot_triggers

Utilities
---------

Independent routines to do various other things

.. currentmodule:: gdas.utils

.. autosummary::
   :toctree: generated/

   create_sound


.. _file_to_segment:

.. ** Extract segment information
.. 
.. The starting and ending UTC times for a specific HDF5 file are determined by using the =Date=, =t0= and =t1= attributes from the metadata. The [[construct_utc_from_metadata][=construct_utc_from_metadata=]] function is then used to calculate the UTC time. Finally, the [[http://software.ligo.org/docs/glue/glue.__segments.segment-class.html][=segment=]] module from the =glue.segments= library is used to represent the range of times in a semi-open interval.
.. 
.. #+BEGIN_SRC python
.. def file_to_segment(hfile,segname):
..     # Extract all atributes from the data
..     attrs = hfile[segname].attrs
..     # Define each attribute
..     dstr, t0, t1 = attrs["Date"], attrs["t0"], attrs["t1"]
..     # Construct GPS starting time from data
..     start_utc = construct_utc_from_metadata(dstr, t0)
..     # Construct GPS starting time from data
..     end_utc = construct_utc_from_metadata(dstr, t1)
..     # Represent the range of times in the semi-open interval
..     return segment(start_utc,end_utc)
.. #+END_SRC
.. 
.. ** Constructing UTC from metadata
.. <<construct_utc_from_metadata>>
.. 
.. #+BEGIN_SRC python
.. def construct_utc_from_metadata(datestr, t0str):
..     instr = "%d-%d-%02dT" % tuple(map(int, datestr.split('/')))
..     instr += t0str
..     t = Time(instr, format='isot', scale='utc')
..     return t.gps
.. #+END_SRC
.. 
.. ** Generate time series
.. <<generate_timeseries>>
.. 
.. #+BEGIN_SRC python
.. def generate_timeseries(data_list, setname="MagneticFields"):
..     full_data = TimeSeriesList()
..     for seg in sorted(data_list):
..         hfile = h5py.File(data_list[seg], "r")
..         full_data.append(retrieve_data_timeseries(hfile, "MagneticFields"))
..         hfile.close()
..     return full_data
.. #+END_SRC
.. 
.. ** Retrieve data time series
.. <<retrieve_data_timeseries>>
.. 
.. #+BEGIN_SRC python
.. def retrieve_data_timeseries(hfile, setname):
..     dset = hfile[setname]
..     sample_rate = dset.attrs["SamplingRate(Hz)"]
..     gps_epoch = construct_utc_from_metadata(dset.attrs["Date"], dset.attrs["t0"])
..     data = retrieve_channel_data(hfile, setname)
..     ts_data = TimeSeries(data, sample_rate=sample_rate, epoch=gps_epoch)
..     return ts_data
.. #+END_SRC
.. 
.. ** Retrieve channel data
.. <<retrieve_channel_data>>
.. 
.. #+BEGIN_SRC python
.. def retrieve_channel_data(hfile, setname):
..     return hfile[setname][:]
.. #+END_SRC
.. 
.. .. _calculate_spectral_correlation:
.. 
.. ** Two point spectral correlation
.. 
.. For our data, we apply a Tukey window whose flat bit corresponds to =window_fraction= (in percentage) of the segment length (in samples) used for PSD estimation (i.e. =fft_window_len=). This can be done by using the [[http://software.ligo.org/docs/lalsuite/lal/_window_8c_source.html#l00597][=CreateTukeyREAL8Window=]] module from the =lal= library. 
.. 
.. #+BEGIN_SRC python
.. def calculate_spectral_correlation(fft_window_len, wtype='hann', window_fraction=None):
..     if wtype == 'hann':
..         window = lal.CreateHannREAL8Window(fft_window_len)
..     elif wtype == 'tukey':
..         window = lal.CreateTukeyREAL8Window(fft_window_len, window_fraction)
..     else:
..         raise ValueError("Can't handle window type %s" % wtype)
.. #+END_SRC
.. 
.. Once the window is built, a new frequency plan is created which will help performing a [[http://fourier.eng.hmc.edu/e101/lectures/fourier_transform_d/node1.html][forward transform]] on the data. This is done with the [[http://software.ligo.org/docs/lalsuite/lal/group___real_f_f_t__h.html#gac4413752db2d19cbe48742e922670af4][=CreateForwardREAL8FFTPlan=]] module which takes as argument the total number of points in the real data and the measurement level for plan creation (here 1 stands for measuring the best plan).
.. 
.. #+BEGIN_SRC python
..     fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)
.. #+END_SRC
.. 
.. We can finally compute and return the two-point spectral correlation function for the whitened frequency series (=fft_plan=) from the window applied to the original time series using the [[http://software.ligo.org/docs/lalsuite/lal/group___time_freq_f_f_t__h.html#ga2bd5c4258eff57cc80103d2ed489e076][=REAL8WindowTwoPointSpectralCorrelation=]] module.
.. 
.. #+BEGIN_SRC python
..     return window, lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)
.. #+END_SRC
.. 
.. ** Create filter bank
.. <<create_filter_bank>>
.. 
.. The construction of a filter bank is fairly simple. For each channel, a frequency domain channel filter function will be created using the [[http://software.ligo.org/docs/lalsuite/lalburst/group___e_p_search__h.html#ga899990cbd45111ba907772650c265ec9][=CreateExcessPowerFilter=]] module from the =lalburst= package. Each channel filter is divided by the square root of the PSD frequency series prior to normalization, which has the effect of de-emphasizing frequency bins with high noise content, and is called "over whitening". The data and metadata are finally stored in the =filter_fseries= and =filter_bank= arrays respectively. Finally, we store on a final array, called =np_filters= the all time-series generated from each filter so that we can plot them afterwards
.. 
.. #+BEGIN_SRC python
.. def create_filter_bank(delta_f, flow, band, nchan, psd, spec_corr):
..     lal_psd = psd.lal()
..     lal_filters, np_filters = [],[]
..     for i in range(nchan):
..         lal_filter = lalburst.CreateExcessPowerFilter(flow + i*band, band, lal_psd, spec_corr)
..         np_filters.append(Spectrum.from_lal(lal_filter))
..         lal_filters.append(lal_filter)
..     return filter_fseries, lal_filters, np_filters
.. #+END_SRC
.. 
.. ** Compute filter inner products with themselves
.. <<compute_filter_ips_self>>
.. #+BEGIN_SRC python
.. def compute_filter_ips_self(lal_filters, spec_corr, psd=None):
..     """
..     Compute a set of inner products of input filters with themselves. If psd
..     argument is given, the unwhitened filter inner products will be returned.
..     """
..     return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, psd) for f in lal_filters])
.. #+END_SRC
.. 
.. ** Compute filter inner products with adjecant filters
.. <<compute_filter_ips_adjacent>>
.. 
.. #+BEGIN_SRC python
.. def compute_filter_ips_adjacent(lal_filters, spec_corr, psd=None):
..     """
..     Compute a set of filter inner products between input adjacent filters.
..     If psd argument is given, the unwhitened filter inner products will be
..     returned. The returned array index is the inner product between the
..     lal_filter of the same index, and its (array) adjacent filter --- assumed
..     to be the frequency adjacent filter.
..     """
..     return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f1, f2, spec_corr, psd) for f1, f2 in zip(lal_filters[:-1], lal_filters[1:])])
.. #+END_SRC
.. 
.. .. _compute_channel_renomalization:
.. 
.. Compute channel renormalization
.. -------------------------------
.. 
.. Compute the renormalization for the base filters up to a given bandwidth.
.. 
.. #+BEGIN_SRC python
.. def compute_channel_renomalization(nc_sum, lal_filters, spec_corr, nchans, verbose=True):
..     mu_sq = (nc_sum+1)*numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None) for f in lal_filters])
..     # Uncomment to get all possible frequency renormalizations
..     #for n in xrange(nc_sum, nchans): # channel position index
..     for n in xrange(nc_sum, nchans, nc_sum+1): # channel position index
..         for k in xrange(0, nc_sum): # channel sum index
..             # FIXME: We've precomputed this, so use it instead
..             mu_sq[n] += 2*lalburst.ExcessPowerFilterInnerProduct(lal_filters[n-k], lal_filters[n-1-k], spec_corr, None)
..     #print mu_sq[nc_sum::nc_sum+1]
..     return mu_sq
.. #+END_SRC
.. 
.. ** Measure root-sum-square strain (hrss)
.. <<measure_hrss>>
.. 
.. #+BEGIN_SRC python
.. def measure_hrss(z_j_b, uw_ss_ii, uw_ss_ij, w_ss_ij, delta_f, delta_t, filter_len, dof):
..     """
..     Approximation of unwhitened sum of squares signal energy in a given EP tile.
..     See T1200125 for equation number reference.
..     z_j_b      - time frequency map block which the constructed tile covers
..     uw_ss_ii   - unwhitened filter inner products
..     uw_ss_ij   - unwhitened adjacent filter inner products
..     w_ss_ij    - whitened adjacent filter inner products
..     delta_f    - frequency binning of EP filters
..     delta_t    - native time resolution of the time frequency map
..     filter_len - number of samples in a fitler
..     dof        - degrees of freedom in the tile (twice the time-frequency area)
..     """
..     s_j_b_avg = uw_ss_ii * delta_f / 2
..     # unwhitened sum of squares of wide virtual filter
..     s_j_nb_avg = uw_ss_ii.sum() / 2 + uw_ss_ij.sum()
..     s_j_nb_avg *= delta_f
..     s_j_nb_denom = s_j_b_avg.sum() + 2 * 2 / filter_len * \
..         numpy.sum(numpy.sqrt(s_j_b_avg[:-1] * s_j_b_avg[1:]) * w_ss_ij)
..     # eqn. 62
..     uw_ups_ratio = s_j_nb_avg / s_j_nb_denom
..     # eqn. 63 -- approximation of unwhitened signal energy time series
..     # FIXME: The sum in this equation is over nothing, but indexed by frequency
..     # I'll make that assumption here too.
..     s_j_nb = numpy.sum(z_j_b.T * numpy.sqrt(s_j_b_avg), axis=0)
..     s_j_nb *= numpy.sqrt(uw_ups_ratio / filter_len * 2)
..     # eqn. 64 -- approximate unwhitened signal energy minus noise contribution
..     # FIXME: correct axis of summation?
..     return math.sqrt(numpy.sum(numpy.absolute(s_j_nb)**2) * delta_t - s_j_nb_avg * dof * delta_t)
.. #+END_SRC
.. 
.. ** Unwhitened inner products filtering
.. <<uw_sum_sq>>
.. 
.. #+BEGIN_SRC python
.. # < s^2_j(f_1, b) > = 1 / 2 / N * \delta_t EPIP{\Theta, \Theta; P}
.. def uw_sum_sq(filter1, filter2, spec_corr, psd):
..     return lalburst.ExcessPowerFilterInnerProduct(filter1, filter2, spec_corr, psd)
.. #+END_SRC
.. 
.. ** Unwhitened sum of squares signal
.. <<measure_hrss_slowly>>
.. 
.. #+BEGIN_SRC python
.. def measure_hrss_slowly(z_j_b, lal_filters, spec_corr, psd, delta_t, dof):
..     """
..     Approximation of unwhitened sum of squares signal energy in a given EP tile.
..     See T1200125 for equation number reference. NOTE: This function is deprecated
..     in favor of measure_hrss, since it requires recomputation of many inner products,
..     making it particularly slow.
..     """
..     # FIXME: Make sure you sum in time correctly
..     # Number of finest bands in given tile
..     nb = len(z_j_b)
..     # eqn. 56 -- unwhitened mean square of filter with itself
..     uw_ss_ii = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i], spec_corr, psd) for i in range(nb)])
..     s_j_b_avg = uw_ss_ii * lal_filters[0].deltaF / 2
..     # eqn. 57 -- unwhitened mean square of filter with adjacent filter
..     uw_ss_ij = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i+1], spec_corr, psd) for i in range(nb-1)])
..     # unwhitened sum of squares of wide virtual filter
..     s_j_nb_avg = uw_ss_ii.sum() / 2 + uw_ss_ij.sum()
..     s_j_nb_avg *= lal_filters[0].deltaF
.. 
..     # eqn. 61
..     w_ss_ij = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i+1], spec_corr, None) for i in range(nb-1)])
..     s_j_nb_denom = s_j_b_avg.sum() + 2 * 2 / len(lal_filters[0].data.data) * \
..         (numpy.sqrt(s_j_b_avg[:-1] * s_j_b_avg[1:]) * w_ss_ij).sum()
.. 
..     # eqn. 62
..     uw_ups_ratio = s_j_nb_avg / s_j_nb_denom
.. 
..     # eqn. 63 -- approximation of unwhitened signal energy time series
..     # FIXME: The sum in this equation is over nothing, but indexed by frequency
..     # I'll make that assumption here too.
..     s_j_nb = numpy.sum(z_j_b.T * numpy.sqrt(s_j_b_avg), axis=0)
..     s_j_nb *= numpy.sqrt(uw_ups_ratio / len(lal_filters[0].data.data) * 2)
..     # eqn. 64 -- approximate unwhitened signal energy minus noise contribution
..     # FIXME: correct axis of summation?
..     return math.sqrt((numpy.absolute(s_j_nb)**2).sum() * delta_t - s_j_nb_avg * dof * delta_t)
.. #+END_SRC
.. 
.. ** Measure root-mean square strain poorly
.. <<measure_hrss_poorly>>
.. 
.. #+BEGIN_SRC python
.. def measure_hrss_poorly(tile_energy, sub_psd):
..     return math.sqrt(tile_energy / numpy.average(1.0 / sub_psd) / 2)
.. #+END_SRC
.. 
.. ** List triggers from map
.. <<trigger_list_from_map>>
.. 
.. #+BEGIN_SRC python
.. def trigger_list_from_map(tfmap, event_list, threshold, start_time, start_freq, duration, band, df, dt, psd=None):
.. 
..     # FIXME: If we don't convert this the calculation takes forever --- but we should convert it once and handle deltaF better later
..     if psd is not None:
..         npy_psd = psd.numpy()
.. 
..     start_time = LIGOTimeGPS(float(start_time))
..     ndof = 2 * duration * band
.. 
..     spanf, spant = tfmap.shape[0] * df, tfmap.shape[1] * dt
..     print "Processing %.2fx%.2f time-frequency map." % (spant, spanf)
.. 
..     for i, j in zip(*numpy.where(tfmap > threshold)):
..         event = event_list.RowType()
.. 
..         # The points are summed forward in time and thus a `summed point' is the
..         # sum of the previous N points. If this point is above threshold, it
..         # corresponds to a tile which spans the previous N points. However, th
..         # 0th point (due to the convolution specifier 'valid') is actually
..         # already a duration from the start time. All of this means, the +
..         # duration and the - duration cancels, and the tile 'start' is, by
..         # definition, the start of the time frequency map if j = 0
..         # FIXME: I think this needs a + dt/2 to center the tile properly
..         event.set_start(start_time + float(j * dt))
..         event.set_stop(start_time + float(j * dt) + duration)
..         event.set_peak(event.get_start() + duration / 2)
..         event.central_freq = start_freq + i * df + 0.5 * band
.. 
..         event.duration = duration
..         event.bandwidth = band
..         event.chisq_dof = ndof
.. 
..         event.snr = math.sqrt(tfmap[i,j] / event.chisq_dof - 1)
..         # FIXME: Magic number 0.62 should be determine empircally
..         event.confidence = -lal.LogChisqCCDF(event.snr * 0.62, event.chisq_dof * 0.62)
..         if psd is not None:
..             # NOTE: I think the pycbc PSDs always start at 0 Hz --- check
..             psd_idx_min = int((event.central_freq - event.bandwidth / 2) / psd.delta_f)
..             psd_idx_max = int((event.central_freq + event.bandwidth / 2) / psd.delta_f)
.. 
..             # FIXME: heuristically this works better with E - D -- it's all
..             # going away with the better h_rss calculation soon anyway
..             event.amplitude = measure_hrss_poorly(tfmap[i,j] - event.chisq_dof, npy_psd[psd_idx_min:psd_idx_max])
..         else:
..             event.amplitude = None
.. 
..         event.process_id = None
..         event.event_id = event_list.get_next_id()
..         event_list.append(event)
.. #+END_SRC
.. 
.. ** Determine output segment
.. <<determine_output_segment>>
.. 
.. #+BEGIN_SRC python
.. def determine_output_segment(inseg, dt_stride, sample_rate, window_fraction=0.0):
..     """
..     Given an input data stretch segment inseg, a data block stride dt_stride, the data sample rate, and an optional window_fraction, return the amount of data that can be processed without corruption effects from the window.
.. 
..     If window_fration is set to 0 (default), assume no windowing.
..     """
..     # Amount to overlap successive blocks so as not to lose data
..     window_overlap_samples = window_fraction * sample_rate
..     outseg = inseg.contract(window_fraction * dt_stride / 2)
.. 
..     # With a given dt_stride, we cannot process the remainder of this data
..     remainder = math.fmod(abs(outseg), dt_stride * (1 - window_fraction))
..     # ...so make an accounting of it
..     outseg = segment(outseg[0], outseg[1] - remainder)
..     return outseg
.. #+END_SRC
.. 
.. ** Make tiles
.. <<make_tiles>>
.. 
.. #+BEGIN_SRC python
.. def make_tiles(tf_map, nc_sum, mu_sq):
..     tiles = numpy.zeros(tf_map.shape)
..     sum_filter = numpy.ones(nc_sum+1)
..     # Here's the deal: we're going to keep only the valid output and
..     # it's *always* going to exist in the lowest available indices
..     for t in xrange(tf_map.shape[1]):
..         # Sum and drop correlate tiles
..         # FIXME: don't drop correlated tiles
..         output = numpy.convolve(tf_map[:,t], sum_filter, 'valid')[::nc_sum+1]
..         #output = fftconvolve(tf_map[:,t], sum_filter, 'valid')[::nc_sum+1]
..         tiles[:len(output),t] = numpy.absolute(output) / math.sqrt(2)
..     return tiles[:len(output)]**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)
.. #+END_SRC
.. 
.. ** Create a time frequency map
.. <<make_indp_tiles>>
.. 
.. In this function, we create a time frequency map with resolution similar than =tf_map= but rescale by a factor of =nc_sum= + 1. All tiles will be independent up to overlap from the original tiling. The =mu_sq= is applied to the resulting addition to normalize the outputs to be zero-mean unit-variance Gaussian variables (if the input is Gaussian).
.. 
.. #+BEGIN_SRC python
.. def make_indp_tiles(tf_map, nc_sum, mu_sq):
..     tiles = tf_map.copy()
..     # Here's the deal: we're going to keep only the valid output and
..     # it's *always* going to exist in the lowest available indices
..     stride = nc_sum + 1
..     for i in xrange(tiles.shape[0]/stride):
..         numpy.absolute(tiles[stride*i:stride*(i+1)].sum(axis=0), tiles[stride*(i+1)-1])
..     return tiles[nc_sum::nc_sum+1].real**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)
.. #+END_SRC
.. 
.. ** Create output filename
.. <<make_filename>>
.. 
.. #+BEGIN_SRC python
.. def make_filename(ifo, seg, tag="excesspower", ext="xml.gz"):
..     if isinstance(ifo, str):
..         ifostr = ifo
..     else:
..         ifostr = "".join(ifo)
..     st_rnd, end_rnd = int(math.floor(seg[0])), int(math.ceil(seg[1]))
..     dur = end_rnd - st_rnd
..     return "%s-%s-%d-%d.%s" % (ifostr, tag, st_rnd, dur, ext)
.. #+END_SRC

