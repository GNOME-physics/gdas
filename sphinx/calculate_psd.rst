.. _psdestimate:

Calculate Power Spectral Density (PSD)
======================================

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

Code access
-----------

.. currentmodule:: gdas.epower

.. autosummary::
   :toctree: generated/

   calculate_psd



