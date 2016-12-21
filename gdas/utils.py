"""Filtering methods"""

def check_filtering_settings(sample_rate,min_frequency,max_frequency,channels=None,tile_bandwidth=None):
    """
    Check filtering settings and define the total number of channels
    and bandwidth to use for filter bank.

    Parameters
    ----------
    sample_rate : float
      Sampling rate in Hz of the data retrieved from the metadata
    min_frequency : float
      Lowest frequency of the filter bank
    max_frequency : float
      Highest frequency of the filter bank
    channels : int
      Number of frequency channels to use
    tile_bandwidth : float
      Bandwidth of the finest filters

    Return
    ------
    nchans, band, flow : int, float, float
      Number of channels, filter bandwidth, initial frequency offset
    """
    # Check if tile maximum frequency is not defined
    if max_frequency is None or max_frequency>sample_rate/2.:
        # Set the tile maximum frequency equal to the Nyquist frequency (i.e. half the sampling rate)
        max_frequency = sample_rate / 2.0
    # Check whether or not tile bandwidth and channel are defined
    if tile_bandwidth is None and channels is None:
        # Exit program with error message
        exit("Either --tile-bandwidth or --channels must be specified to set up time-frequency plane")
    else:
        # Define as assert statement that tile maximum frequency larger than its minimum frequency
        assert max_frequency >= min_frequency
        # Define spectral band of data
        data_band = max_frequency - min_frequency
        # Check if tile bandwidth or channel is defined
        if tile_bandwidth is not None:
            # Define number of possible filter bands
            band = tile_bandwidth
            nchans = channels = int(data_band / tile_bandwidth) - 1
        elif channels is not None:
            # Define filter bandwidth 
            band = tile_bandwidth = data_band / channels
            nchans = channels - 1
        assert channels > 1
    # Lowest frequency of the filter bank
    flow = min_frequency
    return nchans,band,flow

def calculate_spectral_correlation(fft_window_len, wtype='hann', window_fraction=None):
    """
    Calculate the two point spectral correlation introduced by windowing
    the data before transforming to the frequency domain -- valid choices
    are 'hann' and 'tukey'. The window_fraction parameter only has meaning
    for wtype='tukey'.
    """
    if wtype == 'hann':
        window = lal.CreateHannREAL8Window(fft_window_len)
    elif wtype == 'tukey':
        window = lal.CreateTukeyREAL8Window(fft_window_len, window_fraction)
    else:
        raise ValueError("Can't handle window type %s" % wtype)
    fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)
    return window, lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)

def create_filter_bank(delta_f, flow, band, nchan, psd, spec_corr,fmin,fmax):
    lal_psd = psd.lal()
    lal_filters, np_filters = [], []
    for i in range(nchan):
        lal_filter = lalburst.CreateExcessPowerFilter(flow + i*band, band, lal_psd, spec_corr)
        np_filters.append(Spectrum.from_lal(lal_filter))
        lal_filters.append(lal_filter)
    return lal_filters, np_filters

def compute_filter_ips_self(lal_filters, spec_corr, psd=None):
    """
    Compute a set of inner products of input filters with themselves. If psd
    argument is given, the unwhitened filter inner products will be returned.
    """
    return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, psd) for f in lal_filters])

def compute_filter_ips_adjacent(lal_filters, spec_corr, psd=None):
    """
    Compute a set of filter inner products between input adjacent filters.
    If psd argument is given, the unwhitened filter inner products will be
    returned. The returned array index is the inner product between the
    lal_filter of the same index, and its (array) adjacent filter --- assumed
    to be the frequency adjacent filter.
    """
    return numpy.array([lalburst.ExcessPowerFilterInnerProduct(f1, f2, spec_corr, psd) for f1, f2 in zip(lal_filters[:-1], lal_filters[1:])])

def compute_channel_renormalization(nc_sum, lal_filters, spec_corr, nchans, verbose=True):
    """
    Compute the renormalization for the base filters up to a given bandwidth.
    """
    mu_sq = (nc_sum+1)*numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None) for f in lal_filters])
    # Uncomment to get all possible frequency renormalizations
    #for n in xrange(nc_sum, nchans): # channel position index
    for n in xrange(nc_sum, nchans, nc_sum+1): # channel position index
        for k in xrange(0, nc_sum): # channel sum index
            # FIXME: We've precomputed this, so use it instead
            mu_sq[n] += 2*lalburst.ExcessPowerFilterInnerProduct(lal_filters[n-k], lal_filters[n-1-k], spec_corr, None)
    #print mu_sq[nc_sum::nc_sum+1]
    return mu_sq

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

# < s^2_j(f_1, b) > = 1 / 2 / N * \delta_t EPIP{\Theta, \Theta; P}
def uw_sum_sq(filter1, filter2, spec_corr, psd):
    return lalburst.ExcessPowerFilterInnerProduct(filter1, filter2, spec_corr, psd)

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

def measure_hrss_poorly(tile_energy, sub_psd):
    return math.sqrt(tile_energy / numpy.average(1.0 / sub_psd) / 2)

def trigger_list_from_map(tfmap, event_list, threshold, start_time, start_freq, duration, band, df, dt, psd=None):
    
    # FIXME: If we don't convert this the calculation takes forever ---
    # but we should convert it once and handle deltaF better later
    if psd is not None:
        npy_psd = psd.numpy()
        
    start_time = LIGOTimeGPS(float(start_time))
    ndof = 2 * duration * band
    
    for i, j in zip(*numpy.where(tfmap > threshold)):
        event = event_list.RowType()
        
        # The points are summed forward in time and thus a `summed point' is the
        # sum of the previous N points. If this point is above threshold, it
        # corresponds to a tile which spans the previous N points. However, the
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

#
# Optimization plan: If we keep the summed complex TF plane in known indices,
# we can save ourselves individual sums at wider frequency resolutions.
# Caveats:
#   1.  We have to keep track of where we're storing things
#   2.  We have to do it from the finest resolution (for *all* t0s) and work our way up
#
# In the end, I think this is a Haar wavelet transform. Look into it
#
def make_indp_tiles(tf_map, nc_sum, mu_sq):
    """
    Create a time frequency map with resolution of tf_map binning
    divided by nc_sum + 1. All tiles will be independent up to
    overlap from the original tiling. The mu_sq is applied to the
    resulting addition to normalize the outputs to be zero-mean
    unit-variance Gaussian variables (if the input is Gaussian).
    """
    tiles = tf_map.copy()
    # Here's the deal: we're going to keep only the valid output and
    # it's *always* going to exist in the lowest available indices
    stride = nc_sum + 1
    for i in xrange(tiles.shape[0]/stride):
        numpy.absolute(tiles[stride*i:stride*(i+1)].sum(axis=0), tiles[stride*(i+1)-1])
    return tiles[nc_sum::nc_sum+1].real**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)

def make_filename(ifo, seg, tag="excesspower", ext="xml.gz"):
    if isinstance(ifo, str):
        ifostr = ifo
    else:
        ifostr = "".join(ifo)
    st_rnd, end_rnd = int(math.floor(seg[0])), int(math.ceil(seg[1]))
    dur = end_rnd - st_rnd
    #return "%s-%s-%d-%d.%s" % (ifostr, tag, st_rnd, dur, ext)
    return "%s.%s" % (tag, ext)

def tprint(t0=False,t1=False,t2=False,t3=False):
    if t0: text  = "{:>8} |".format('%.3f'%(time.time()-t0))
    if t1: text += "{:>8} |".format('%.3f'%(time.time()-t1))
    if t2: text += "{:>8} |".format('%.3f'%(time.time()-t2))
    if t3: text += "{:>8} |".format('%.3f'%(time.time()-t3))
    return text
