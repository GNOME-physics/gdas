Submodules
==========

.. _file_to_segment:

file_to_segment
---------------

The starting and ending UTC times for a specific HDF5 file are determined by using the =Date=, =t0= and =t1= attributes from the metadata. The [[construct_utc_from_metadata][=construct_utc_from_metadata=]] function is then used to calculate the UTC time. Finally, the [[http://software.ligo.org/docs/glue/glue.__segments.segment-class.html][=segment=]] module from the =glue.segments= library is used to represent the range of times in a semi-open interval. ::

  def file_to_segment(hfile,segname):
      # Extract all atributes from the data
      attrs = hfile[segname].attrs
      # Define each attribute
      dstr, t0, t1 = attrs["Date"], attrs["t0"], attrs["t1"]
      # Construct GPS starting time from data
      start_utc = construct_utc_from_metadata(dstr, t0)
      # Construct GPS starting time from data
      end_utc = construct_utc_from_metadata(dstr, t1)
      # Represent the range of times in the semi-open interval
      return segment(start_utc,end_utc)

.. _construct_utc_from_metadata:
      
construct_utc_from_metadata
---------------------------
::

  def construct_utc_from_metadata(datestr, t0str):
      instr = "%d-%d-%02dT" % tuple(map(int, datestr.split('/')))
      instr += t0str
      t = Time(instr, format='isot', scale='utc')
      return t.gps

.. _generate_timeseries:
      
generate_timeseries
-------------------
::
   

  def generate_timeseries(data_list, setname="MagneticFields"):
      full_data = TimeSeriesList()
      for seg in sorted(data_list):
          hfile = h5py.File(data_list[seg], "r")
          full_data.append(retrieve_data_timeseries(hfile, "MagneticFields"))
          hfile.close()
      return full_data

.. _retrieve_data_timeseries:
      
retrieve_data_timeseries
------------------------
::

  def retrieve_data_timeseries(hfile, setname):
      dset = hfile[setname]
      sample_rate = dset.attrs["SamplingRate(Hz)"]
      gps_epoch = construct_utc_from_metadata(dset.attrs["Date"], dset.attrs["t0"])
      data = retrieve_channel_data(hfile, setname)
      ts_data = TimeSeries(data, sample_rate=sample_rate, epoch=gps_epoch)
      return ts_data

.. _retrieve_channel_data:
      
retrieve_channel_data
---------------------
::
   
  def retrieve_channel_data(hfile, setname):
      return hfile[setname][:]

.. _calculate_spectral_correlation:

calculate_spectral_correlation
------------------------------

For our data, we apply a Tukey window whose flat bit corresponds to ``window_fraction`` (in percentage) of the segment length (in samples) used for PSD estimation (i.e. ``fft_window_len``). This can be done by using the `CreateTukeyREAL8Window <http://software.ligo.org/docs/lalsuite/lal/_window_8c_source.html#l00597>`_ module from the ``lal`` library. ::

  def calculate_spectral_correlation(fft_window_len, wtype='hann', window_fraction=None):
      if wtype == 'hann':
          window = lal.CreateHannREAL8Window(fft_window_len)
      elif wtype == 'tukey':
          window = lal.CreateTukeyREAL8Window(fft_window_len, window_fraction)
      else:
          raise ValueError("Can't handle window type %s" % wtype)
  
Once the window is built, a new frequency plan is created which will help performing a `forward transform <http://fourier.eng.hmc.edu/e101/lectures/fourier_transform_d/node1.html>`_ on the data. This is done with the `CreateForwardREAL8FFTPlan <http://software.ligo.org/docs/lalsuite/lal/group___real_f_f_t__h.html#gac4413752db2d19cbe48742e922670af4>`_ module which takes as argument the total number of points in the real data and the measurement level for plan creation (here 1 stands for measuring the best plan). ::

  fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)

We can finally compute and return the two-point spectral correlation function for the whitened frequency series (``fft_plan``) from the window applied to the original time series using the `REAL8WindowTwoPointSpectralCorrelation <http://software.ligo.org/docs/lalsuite/lal/group___time_freq_f_f_t__h.html#ga2bd5c4258eff57cc80103d2ed489e076>`_ module. ::

  return window, lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)

.. _create_filter_bank:

create_filter_bank
------------------

The construction of a filter bank is fairly simple. For each channel, a frequency domain channel filter function will be created using the [[http://software.ligo.org/docs/lalsuite/lalburst/group___e_p_search__h.html#ga899990cbd45111ba907772650c265ec9][=CreateExcessPowerFilter=]] module from the =lalburst= package. Each channel filter is divided by the square root of the PSD frequency series prior to normalization, which has the effect of de-emphasizing frequency bins with high noise content, and is called "over whitening". The data and metadata are finally stored in the =filter_fseries= and =filter_bank= arrays respectively. Finally, we store on a final array, called =np_filters= the all time-series generated from each filter so that we can plot them afterwards. ::

  def create_filter_bank(delta_f, flow, band, nchan, psd, spec_corr):
      lal_psd = psd.lal()
      lal_filters, np_filters = [],[]
      for i in range(nchan):
          lal_filter = lalburst.CreateExcessPowerFilter(flow + i*band, band, lal_psd, spec_corr)
          np_filters.append(Spectrum.from_lal(lal_filter))
          lal_filters.append(lal_filter)
      return filter_fseries, lal_filters, np_filters

.. _compute_filter_ips_self:

compute_filter_ips_self
-----------------------
::
   
  def compute_filter_ips_self(lal_filters, spec_corr, psd=None):
      """
      Compute a set of inner products of input filters with themselves. If psd
      argument is given, the unwhitened filter inner products will be returned.
      """
      return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, psd) for f in lal_filters])

.. _compute_filter_ips_adjacent:
      
compute_filter_ips_adjacent
---------------------------
::

  def compute_filter_ips_adjacent(lal_filters, spec_corr, psd=None):
      """
      Compute a set of filter inner products between input adjacent filters.
      If psd argument is given, the unwhitened filter inner products will be
      returned. The returned array index is the inner product between the
      lal_filter of the same index, and its (array) adjacent filter --- assumed
      to be the frequency adjacent filter.
      """
      return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f1, f2, spec_corr, psd) for f1, f2 in zip(lal_filters[:-1], lal_filters[1:])])

.. _compute_channel_renomalization:

Compute channel renormalization
-------------------------------

Compute the renormalization for the base filters up to a given bandwidth. ::

  def compute_channel_renomalization(nc_sum, lal_filters, spec_corr, nchans, verbose=True):
      mu_sq = (nc_sum+1)*numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None) for f in lal_filters])
      # Uncomment to get all possible frequency renormalizations
      #for n in xrange(nc_sum, nchans): # channel position index
      for n in xrange(nc_sum, nchans, nc_sum+1): # channel position index
          for k in xrange(0, nc_sum): # channel sum index
              # FIXME: We've precomputed this, so use it instead
              mu_sq[n] += 2*lalburst.ExcessPowerFilterInnerProduct(lal_filters[n-k], lal_filters[n-1-k], spec_corr, None)
      #print mu_sq[nc_sum::nc_sum+1]
      return mu_sq

.. _measure_hrss:
      
measure_hrss
------------
::

  def measure_hrss(z_j_b, uw_ss_ii, uw_ss_ij, w_ss_ij, delta_f, delta_t, filter_len, dof):
      """
      Approximation of unwhitened sum of squares signal energy in a given EP tile.
      See T1200125 for equation number reference.
      z_j_b      - time frequency map block which the constructed tile covers
      uw_ss_ii   - unwhitened filter inner products
      uw_ss_ij   - unwhitened adjacent filter inner products
      w_ss_ij    - whitened adjacent filter inner products
      delta_f    - frequency binning of EP filters
      delta_t    - native time resolution of the time frequency map
      filter_len - number of samples in a fitler
      dof        - degrees of freedom in the tile (twice the time-frequency area)
      """
      s_j_b_avg = uw_ss_ii * delta_f / 2
      # unwhitened sum of squares of wide virtual filter
      s_j_nb_avg = uw_ss_ii.sum() / 2 + uw_ss_ij.sum()
      s_j_nb_avg *= delta_f
      s_j_nb_denom = s_j_b_avg.sum() + 2 * 2 / filter_len * \
          numpy.sum(numpy.sqrt(s_j_b_avg[:-1] * s_j_b_avg[1:]) * w_ss_ij)
      # eqn. 62
      uw_ups_ratio = s_j_nb_avg / s_j_nb_denom
      # eqn. 63 -- approximation of unwhitened signal energy time series
      # FIXME: The sum in this equation is over nothing, but indexed by frequency
      # I'll make that assumption here too.
      s_j_nb = numpy.sum(z_j_b.T * numpy.sqrt(s_j_b_avg), axis=0)
      s_j_nb *= numpy.sqrt(uw_ups_ratio / filter_len * 2)
      # eqn. 64 -- approximate unwhitened signal energy minus noise contribution
      # FIXME: correct axis of summation?
      return math.sqrt(numpy.sum(numpy.absolute(s_j_nb)**2) * delta_t - s_j_nb_avg * dof * delta_t)

.. _uw_sum_sq:
      
uw_sum_sq
---------
::
   
  # < s^2_j(f_1, b) > = 1 / 2 / N * \delta_t EPIP{\Theta, \Theta; P}
  def uw_sum_sq(filter1, filter2, spec_corr, psd):
      return lalburst.ExcessPowerFilterInnerProduct(filter1, filter2, spec_corr, psd)

.. _measure_hrss_slowly:

measure_hrss_slowly
-------------------
::

   def measure_hrss_slowly(z_j_b, lal_filters, spec_corr, psd, delta_t, dof):
      """
      Approximation of unwhitened sum of squares signal energy in a given EP tile.
      See T1200125 for equation number reference. NOTE: This function is deprecated
      in favor of measure_hrss, since it requires recomputation of many inner products,
      making it particularly slow.
      """
      # FIXME: Make sure you sum in time correctly
      # Number of finest bands in given tile
      nb = len(z_j_b)
      # eqn. 56 -- unwhitened mean square of filter with itself
      uw_ss_ii = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i], spec_corr, psd) for i in range(nb)])
      s_j_b_avg = uw_ss_ii * lal_filters[0].deltaF / 2
      # eqn. 57 -- unwhitened mean square of filter with adjacent filter
      uw_ss_ij = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i+1], spec_corr, psd) for i in range(nb-1)])
      # unwhitened sum of squares of wide virtual filter
      s_j_nb_avg = uw_ss_ii.sum() / 2 + uw_ss_ij.sum()
      s_j_nb_avg *= lal_filters[0].deltaF
  
      # eqn. 61
      w_ss_ij = numpy.array([uw_sum_sq(lal_filters[i], lal_filters[i+1], spec_corr, None) for i in range(nb-1)])
      s_j_nb_denom = s_j_b_avg.sum() + 2 * 2 / len(lal_filters[0].data.data) * \
          (numpy.sqrt(s_j_b_avg[:-1] * s_j_b_avg[1:]) * w_ss_ij).sum()
  
      # eqn. 62
      uw_ups_ratio = s_j_nb_avg / s_j_nb_denom
  
      # eqn. 63 -- approximation of unwhitened signal energy time series
      # FIXME: The sum in this equation is over nothing, but indexed by frequency
      # I'll make that assumption here too.
      s_j_nb = numpy.sum(z_j_b.T * numpy.sqrt(s_j_b_avg), axis=0)
      s_j_nb *= numpy.sqrt(uw_ups_ratio / len(lal_filters[0].data.data) * 2)
      # eqn. 64 -- approximate unwhitened signal energy minus noise contribution
      # FIXME: correct axis of summation?
      return math.sqrt((numpy.absolute(s_j_nb)**2).sum() * delta_t - s_j_nb_avg * dof * delta_t)

.. _measure_hrss_poorly:

measure_hrss_poorly
-------------------
::
   
  def measure_hrss_poorly(tile_energy, sub_psd):
      return math.sqrt(tile_energy / numpy.average(1.0 / sub_psd) / 2)

.. _trigger_list_from_map:

trigger_list_from_map
---------------------
::

  def trigger_list_from_map(tfmap, event_list, threshold, start_time, start_freq, duration, band, df, dt, psd=None):
  
      # FIXME: If we don't convert this the calculation takes forever --- but we should convert it once and handle deltaF better later
      if psd is not None:
          npy_psd = psd.numpy()
  
      start_time = LIGOTimeGPS(float(start_time))
      ndof = 2 * duration * band
  
      spanf, spant = tfmap.shape[0] * df, tfmap.shape[1] * dt
      print "Processing %.2fx%.2f time-frequency map." % (spant, spanf)
  
      for i, j in zip(*numpy.where(tfmap > threshold)):
          event = event_list.RowType()
  
          # The points are summed forward in time and thus a `summed point' is the
          # sum of the previous N points. If this point is above threshold, it
          # corresponds to a tile which spans the previous N points. However, th
          # 0th point (due to the convolution specifier 'valid') is actually
          # already a duration from the start time. All of this means, the +
          # duration and the - duration cancels, and the tile 'start' is, by
          # definition, the start of the time frequency map if j = 0
          # FIXME: I think this needs a + dt/2 to center the tile properly
          event.set_start(start_time + float(j * dt))
          event.set_stop(start_time + float(j * dt) + duration)
          event.set_peak(event.get_start() + duration / 2)
          event.central_freq = start_freq + i * df + 0.5 * band
  
          event.duration = duration
          event.bandwidth = band
          event.chisq_dof = ndof
  
          event.snr = math.sqrt(tfmap[i,j] / event.chisq_dof - 1)
          # FIXME: Magic number 0.62 should be determine empircally
          event.confidence = -lal.LogChisqCCDF(event.snr * 0.62, event.chisq_dof * 0.62)
          if psd is not None:
              # NOTE: I think the pycbc PSDs always start at 0 Hz --- check
              psd_idx_min = int((event.central_freq - event.bandwidth / 2) / psd.delta_f)
              psd_idx_max = int((event.central_freq + event.bandwidth / 2) / psd.delta_f)
  
              # FIXME: heuristically this works better with E - D -- it's all
              # going away with the better h_rss calculation soon anyway
              event.amplitude = measure_hrss_poorly(tfmap[i,j] - event.chisq_dof, npy_psd[psd_idx_min:psd_idx_max])
          else:
              event.amplitude = None
  
          event.process_id = None
          event.event_id = event_list.get_next_id()
          event_list.append(event)

.. _determine_output_segment:

determine_output_segment
------------------------
::

  def determine_output_segment(inseg, dt_stride, sample_rate, window_fraction=0.0):
      """
      Given an input data stretch segment inseg, a data block stride dt_stride, the data sample rate, and an optional window_fraction, return the amount of data that can be processed without corruption effects from the window.
  
      If window_fration is set to 0 (default), assume no windowing.
      """
      # Amount to overlap successive blocks so as not to lose data
      window_overlap_samples = window_fraction * sample_rate
      outseg = inseg.contract(window_fraction * dt_stride / 2)
  
      # With a given dt_stride, we cannot process the remainder of this data
      remainder = math.fmod(abs(outseg), dt_stride * (1 - window_fraction))
      # ...so make an accounting of it
      outseg = segment(outseg[0], outseg[1] - remainder)
      return outseg

.. _make_tiles:

make_tiles
----------
::

  def make_tiles(tf_map, nc_sum, mu_sq):
      tiles = numpy.zeros(tf_map.shape)
      sum_filter = numpy.ones(nc_sum+1)
      # Here's the deal: we're going to keep only the valid output and
      # it's *always* going to exist in the lowest available indices
      for t in xrange(tf_map.shape[1]):
          # Sum and drop correlate tiles
          # FIXME: don't drop correlated tiles
          output = numpy.convolve(tf_map[:,t], sum_filter, 'valid')[::nc_sum+1]
          #output = fftconvolve(tf_map[:,t], sum_filter, 'valid')[::nc_sum+1]
          tiles[:len(output),t] = numpy.absolute(output) / math.sqrt(2)
      return tiles[:len(output)]**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)

.. _make_indp_tiles:

make_indp_tiles
---------------

In this function, we create a time frequency map with resolution similar than =tf_map= but rescale by a factor of =nc_sum= + 1. All tiles will be independent up to overlap from the original tiling. The =mu_sq= is applied to the resulting addition to normalize the outputs to be zero-mean unit-variance Gaussian variables (if the input is Gaussian). ::

  def make_indp_tiles(tf_map, nc_sum, mu_sq):
      tiles = tf_map.copy()
      # Here's the deal: we're going to keep only the valid output and
      # it's *always* going to exist in the lowest available indices
      stride = nc_sum + 1
      for i in xrange(tiles.shape[0]/stride):
          numpy.absolute(tiles[stride*i:stride*(i+1)].sum(axis=0), tiles[stride*(i+1)-1])
      return tiles[nc_sum::nc_sum+1].real**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)

.. _make_filename:

make_filename
-------------
::
   
  def make_filename(ifo, seg, tag="excesspower", ext="xml.gz"):
      if isinstance(ifo, str):
          ifostr = ifo
      else:
          ifostr = "".join(ifo)
      st_rnd, end_rnd = int(math.floor(seg[0])), int(math.ceil(seg[1]))
      dur = end_rnd - st_rnd
      return "%s-%s-%d-%d.%s" % (ifostr, tag, st_rnd, dur, ext)
