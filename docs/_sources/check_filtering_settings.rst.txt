Checking filtering settings
===========================

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
