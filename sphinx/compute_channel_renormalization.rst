Normalization of virtual channel
================================

The virtual channels will be used during the excesspower analysis to explore different frequency ranges around each PSD segments and look for possible triggers. Each channel is renormalized using the :ref:`compute_channel_renomalization <compute_channel_renomalization>` internal function. ::

  # Initialise dictionary
  mu_sq_dict = {}
  # nc_sum additional channel adds
  for nc_sum in range(0, int(math.log(nchans, 2))):
      min_band = (len(filter_bank[0].data.data)-1) * filter_bank[0].deltaF / 2
      print tprint(t0,t1),"Calculation for %d %d Hz channels" % (nc_sum+1, min_band)
      nc_sum = 2**nc_sum - 1
      mu_sq_dict[nc_sum] = compute_channel_renomalization(nc_sum, filter_bank, spec_corr, nchans)
