.. _filterbank:
  
Excess Power - Step 4: Computing the filter bank
================================================

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

Code access
-----------

.. currentmodule:: gdas.epower

.. autosummary::
   :toctree: generated/

   create_filter_bank
