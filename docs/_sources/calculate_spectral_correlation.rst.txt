Two point spectral correlation
==============================

This part determines how much data on either side of the tukey window is to be discarded. Nominally, this means that one will lose ``window_fraction`` * ``args.psd_segment_length`` to corruption from the window, i.e. this is simply discarded. This is tuned to give an integer offset when used with ``args.psd_segment_length`` equal to 8, smaller windows will have fractions of integers, but larger powers of two will still preseve this (probably not a big deal in the end). ::

  window_fraction = 0

The two point spectral correlation is then done with the :ref:`calculate_spectral_correlation <calculate_spectral_correlation>` function which will return both the Tukey window applied to the original time series data and the actual two-point spectral correlation function for the whitened frequency series from the applied whitening window. ::

  # Do two point spectral correlation
  window, spec_corr = calculate_spectral_correlation(seg_len,'tukey',window_fraction=window_fraction)
  window = window.data.data
  window_sigma_sq = numpy.mean(window**2)
  # Pre scale the window by its root mean squared -- see eqn 11 of EP document
  #window /= numpy.sqrt(window_sigma_sq)
