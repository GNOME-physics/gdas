import time,math,os,scipy,glue
from lal import *
from lalburst import *
from glue.ligolw.utils.search_summary import append_search_summary
from glue.ligolw.utils.process        import register_to_xmldoc
from glue.ligolw   import lsctables,ligolw,utils
from glue.lal      import LIGOTimeGPS
from glue.segments import segment
from scipy.signal  import fftconvolve
from pycbc import psd,types,filter
from utils import *
from plots import *

def excess_power(ts_data,psd_segment_length,psd_segment_stride,psd_estimation,window_fraction,tile_fap,station,nchans=None,band=None,fmin=0,fmax=None,max_duration=None):
    """
    Perform excess-power search analysis on magnetic field data.

    Parameters
    ----------
    ts_data : TimeSeries
      Time Series from magnetic field data
    psd_segment_length : float
      Length of each segment in seconds
    psd_segment_stride : float
      Separation between 2 consecutive segments in seconds
    psd_estimation : string
      Average method
    window_fraction : float
      Withening window fraction
    tile_fap : float
      Tile false alarm probability threshold in Gaussian noise.
    nchans : int
      Total number of channels
    band : float
      Tile bandwidth
    fmin : float
      Lowest frequency of the filter bank.
    fmax : float
      Highest frequency of the filter bank

    Notes
    -----
    This method will produce a bunch of time-frequency plots for every
    tile duration and bandwidth analysed as well as a XML file identifying
    all the triggers found in the selected data within the user-defined
    time range.
    """
    #print strain.insert_strain_option_group.__dict__
    print psd.insert_psd_option_group.__dict__
    sample_rate = ts_data.sample_rate
    nchans,band,flow = check_filtering_settings(sample_rate,nchans,band,fmin,fmax)
    seg_len,fd_psd,lal_psd = calculate_psd(ts_data,sample_rate,psd_segment_length,psd_segment_stride,psd_estimation)
    window, spec_corr = calculate_spectral_correlation(seg_len,'tukey',window_fraction=window_fraction)
    window = window.data.data
    window_sigma_sq = numpy.mean(window**2)
    # Pre scale the window by its root mean squared -- see eqn 11 of EP document
    #window /= numpy.sqrt(window_sigma_sq)
    filter_bank, fdb = create_filter_bank(fd_psd.delta_f,flow+band/2,band,nchans,fd_psd,spec_corr,fmin,fmax)
    # This is necessary to compute the mu^2 normalizations
    #white_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, None)
    #unwhite_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, lal_psd)
    # These two are needed for the unwhitened mean square sum (hrss)
    #white_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, None)
    #unwhite_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, lal_psd)
    tdb = convert_to_time_domain(fdb,sample_rate)
    plot_bank(fdb)
    plot_filters(tdb,flow,band)
    mu_sq_dict = compute_channel_renormalization(filter_bank, spec_corr, nchans)
    event_list = lsctables.New(lsctables.SnglBurstTable,
                               ['start_time','start_time_ns','peak_time','peak_time_ns',
                                'duration','bandwidth','central_freq','chisq_dof',
                                'confidence','snr','amplitude','channel','ifo',
                                'process_id','event_id','search','stop_time','stop_time_ns'])
    t_idx_min, t_idx_max = 0, seg_len
    os.system('mkdir -p segments/time-frequency')
    os.system('mkdir -p segments/time-series')
    while t_idx_max <= len(ts_data):        
        start_time, end_time, tmp_ts_data, fs_data = identify_block(ts_data,fd_psd,window,t_idx_min,t_idx_max)
        tf_map = create_tf_plane(fd_psd,nchans,seg_len,filter_bank,band,fs_data)
        plot_spectrogram(numpy.abs(tf_map).T,tmp_ts_data.delta_t,band,ts_data.sample_rate,start_time,end_time,fname='segments/time-frequency/%i-%i.png'%(start_time,end_time))
        for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]: # nc_sum additional channel adds
            nc_sum = 2**nc_sum - 1
            mu_sq = mu_sq_dict[nc_sum]
            max_dof, tiles, us_rate, dt, df = construct_tiles(nc_sum,mu_sq,band,ts_data,tf_map,psd_segment_length,window_fraction,max_duration)
            t3 = time.time()
            for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
                duration  = j / 2.0 / df
                dof_tiles = create_tile_duration(j,df,duration,tiles)
                plot_spectrogram(dof_tiles.T,dt,df,ts_data.sample_rate,start_time,end_time,fname='segments/%i-%i/tf_%02ichans_%02idof.png'%(start_time,end_time,nc_sum+1,2*j))
                threshold = scipy.stats.chi2.isf(tile_fap, j)
                print "|------ Threshold for this level: %f" % threshold
                spant, spanf = dof_tiles.shape[1] * dt, dof_tiles.shape[0] * df
                print "|------ Processing %.2fx%.2f time-frequency map." % (spant, spanf)
                # Since we clip the data, the start time needs to be adjusted accordingly
                window_offset_epoch = fs_data.epoch + psd_segment_length * window_fraction / 2
                trigger_list_from_map(dof_tiles, event_list, threshold, window_offset_epoch, filter_bank[0].f0 + band/2, duration, df, df, dt, None)
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
                    event.amplitude = 0
                print "|------ Total number of events: %d" % len(event_list)
        t_idx_min += int(seg_len * (1 - window_fraction))
        t_idx_max += int(seg_len * (1 - window_fraction))
    create_xml(ts_data,psd_segment_length,window_fraction,event_list,station)
        
def check_filtering_settings(sample_rate,channels,tile_bandwidth,fmin,fmax):
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
    if fmax is None or fmax>sample_rate/2.:
        # Set the tile maximum frequency equal to the Nyquist frequency (i.e. half the sampling rate)
        fmax = sample_rate / 2.0
    # Check whether or not tile bandwidth and channel are defined
    if tile_bandwidth is None and channels is None:
        # Exit program with error message
        exit("Either --tile-bandwidth or --channels must be specified to set up time-frequency plane")
    else:
        # Define as assert statement that tile maximum frequency larger than its minimum frequency
        assert fmax >= fmin
        # Define spectral band of data
        data_band = fmax - fmin
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
    flow = fmin
    return nchans,band,flow

def calculate_psd(ts_data,sample_rate,psd_segment_length,psd_segment_stride,psd_estimation):
    """
    Estimate Power Spectral Density (PSD)

    Parameters
    ----------
    ts_data : TimeSeries
      Time series of magnetic field data
    sample_rate : float
      Sampling rate of data
    psd_segment_length : float
      Length of data segment in seconds
    psd_segment_stride : float
      Separation between consecutive segments in seconds
    psd_estimation : string
      Average method to measure PSD from the data

    Return
    ------
    seg_len, fd_psd, lal_psd : Segment length in sample unit, PSD results
    in 2 different formats.

    Notes
    -----
    Need to contact Chris Pankow for more information on the 2 formats.
    """
    print "|- Estimating PSD from segments of time %.2f s in length, with %.2f s stride..." % (psd_segment_length, psd_segment_stride)
    # Convert time series as array of float
    data = ts_data.astype(numpy.float64)
    # Average method to measure PSD from the data
    avg_method = psd_estimation
    # The segment length for PSD estimation in samples
    seg_len = int(psd_segment_length * sample_rate)
    # The separation between consecutive segments in samples
    seg_stride = int(psd_segment_stride * sample_rate)
    # Lifted from the psd.from_cli module
    fd_psd = psd.welch(data,avg_method=avg_method,seg_len=seg_len,seg_stride=seg_stride)
    # Plot the power spectral density
    plot_spectrum(fd_psd)
    # We need this for the SWIG functions
    lal_psd = fd_psd.lal()
    return seg_len,fd_psd,lal_psd
    
def calculate_spectral_correlation(fft_window_len,wtype='hann',window_fraction=None):
    """
    Calculate the two point spectral correlation introduced by windowing
    the data before transforming to the frequency domain -- valid choices
    are 'hann' and 'tukey'. The window_fraction parameter only has meaning
    for wtype='tukey'.
    """
    print "|- Whitening window and spectral correlation..."
    if wtype == 'hann':
        window = lal.CreateHannREAL8Window(fft_window_len)
    elif wtype == 'tukey':
        window = lal.CreateTukeyREAL8Window(fft_window_len, window_fraction)
    else:
        raise ValueError("Can't handle window type %s" % wtype)
    fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)
    return window, lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)

def create_filter_bank(delta_f,flow,band,nchan,psd,spec_corr,fmin=0,fmax=None):
    """
    Create filter bank
    
    Parameters
    ----------
    delta_f : float
      Bandwidth of each filter
    flow : float
      Lowest frequency of the filter bank
    band : 
    """
    print "|- Create filter..."
    lal_psd = psd.lal()
    lal_filters, np_filters = [], []
    for i in range(nchan):
        lal_filter = lalburst.CreateExcessPowerFilter(flow + i*band, band, lal_psd, spec_corr)
        np_filters.append(Spectrum.from_lal(lal_filter))
        lal_filters.append(lal_filter)
    return lal_filters, np_filters

def convert_to_time_domain(fdb,sample_rate):
    """
    Convert filter bank from frequency to time domain
    
    Parameters
    ----------
    fdb : list
      List of filters from the filter bank in frequency domain
    sample_rate : float
      Sampling rate of magnetic field data
    
    Return
    ------
    tdb : list
      List of filters from the filter bank in time domain
    """
    print "|- Convert all the frequency domain to the time domain..."
    tdb = []
    for fdt in fdb:
        zero_padded = numpy.zeros(int((fdt.f0 / fdt.df).value) + len(fdt))
        st = int((fdt.f0 / fdt.df).value)
        zero_padded[st:st+len(fdt)] = numpy.real_if_close(fdt.value)
        n_freq = int(sample_rate / 2 / fdt.df.value) * 2
        tdt = numpy.fft.irfft(zero_padded, n_freq) * math.sqrt(sample_rate)
        tdt = numpy.roll(tdt, len(tdt)/2)
        tdt = TimeSeries(tdt, name="", epoch=fdt.epoch, sample_rate=sample_rate)
        tdb.append(tdt)
    return tdb

def identify_block(ts_data,fd_psd,window,t_idx_min,t_idx_max):
    """
    Get frequency series of the current block
    
    Parameters
    ----------
    ts_data : TimeSeries
      Time series of magnetic field data
    fd_psd : 
      Power Spectrum Density
    window : 
    t_idx_min : float
      Index in time series of first data point
    t_idx_max : float
      Index in time series of last data point

    Return
    ------
    start_time : float
      Starting time of the block
    end_time : float
      Ending time of the block
    tmp_ts_data : TimeSeries
      Time series magnetic data of the block
    fs_data : FrequencySeries
      Frequency series magnetic data of the block
    """
    # Define starting and ending time of the segment in seconds
    start_time = ts_data.start_time + t_idx_min/float(ts_data.sample_rate)
    end_time = ts_data.start_time + t_idx_max/float(ts_data.sample_rate)
    print "|-- Analyzing block %i to %i (%.2f percent)"%(start_time,end_time,100*float(t_idx_max)/len(ts_data))
    # Model a withen time series for the block
    tmp_ts_data = types.TimeSeries(ts_data[t_idx_min:t_idx_max]*window,delta_t=1./ts_data.sample_rate,epoch=start_time)
    # Save time series in relevant repository
    segfolder = 'segments/%i-%i'%(start_time,end_time)
    os.system('mkdir -p '+segfolder)
    plot_ts(tmp_ts_data,fname='segments/time-series/%i-%i.png'%(start_time,end_time))
    # Convert times series to frequency series
    fs_data = tmp_ts_data.to_frequencyseries()
    print "|-- Frequency series data has variance: %s" % fs_data.data.std()**2
    # Whitening (FIXME: Whiten the filters, not the data)
    fs_data.data /= numpy.sqrt(fd_psd) / numpy.sqrt(2 * fd_psd.delta_f)
    print "|-- Whitened frequency series data has variance: %s" % fs_data.data.std()**2
    return start_time, end_time, tmp_ts_data, fs_data

def create_tf_plane(fd_psd,nchans,seg_len,filter_bank,band,fs_data):
    """
    Create time-frequency map

    Parameters
    ----------
    fd_psd : array
      Power Spectrum Density
    """
    print "|-- Create time-frequency plane for current block"
    # Return the complex snr, along with its associated normalization of the template,
    # matched filtered against the data
    #filter.matched_filter_core(types.FrequencySeries(tmp_filter_bank,delta_f=fd_psd.delta_f),fs_data,h_norm=1,psd=fd_psd,low_frequency_cutoff=filter_bank[0].f0,high_frequency_cutoff=filter_bank[0].f0+2*band)
    print "|-- Filtering all %d channels..." % nchans
    # Initialise 2D zero array 
    tmp_filter_bank = numpy.zeros(len(fd_psd), dtype=numpy.complex128)
    # Initialise 2D zero array for time-frequency map
    tf_map = numpy.zeros((nchans, seg_len), dtype=numpy.complex128)
    # Loop over all the channels
    for i in range(nchans):
        # Reset filter bank series
        tmp_filter_bank *= 0.0
        # Index of starting frequency
        f1 = int(filter_bank[i].f0/fd_psd.delta_f)
        # Index of ending frequency
        f2 = int((filter_bank[i].f0 + 2*band)/fd_psd.delta_f)+1
        # (FIXME: Why is there a factor of 2 here?)
        tmp_filter_bank[f1:f2] = filter_bank[i].data.data * 2
        # Define the template to filter the frequency series with
        template = types.FrequencySeries(tmp_filter_bank, delta_f=fd_psd.delta_f, copy=False)
        # Create filtered series
        filtered_series = filter.matched_filter_core(template,fs_data,h_norm=None,psd=None,low_frequency_cutoff=filter_bank[i].f0,high_frequency_cutoff=filter_bank[i].f0+2*band)
        # Include filtered series in the map
        tf_map[i,:] = filtered_series[0].numpy()
    return tf_map

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

def compute_channel_renormalization(filter_bank, spec_corr, nchans):
    """
    Compute the renormalization for the base filters up to a given bandwidth.
    """
    mu_sq_dict = {}
    for nc_sum in range(0, int(math.log(nchans, 2))):
        min_band = (len(filter_bank[0].data.data)-1) * filter_bank[0].deltaF / 2
        print "|- Calculation for %d %dHz channels" % (nc_sum+1, min_band)
        nc_sum = 2**nc_sum - 1
        mu_sq = (nc_sum+1)*numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None) for f in filter_bank])
        # Uncomment to get all possible frequency renormalizations
        #for n in xrange(nc_sum, nchans): # channel position index
        for n in xrange(nc_sum, nchans, nc_sum+1): # channel position index
            for k in xrange(0, nc_sum): # channel sum index
                # FIXME: We've precomputed this, so use it instead
                mu_sq[n] += 2*lalburst.ExcessPowerFilterInnerProduct(filter_bank[n-k], filter_bank[n-1-k], spec_corr, None)
        #print mu_sq[nc_sum::nc_sum+1]
        mu_sq_dict[nc_sum] = mu_sq
    return mu_sq_dict

def measure_hrss(z_j_b, uw_ss_ii, uw_ss_ij, w_ss_ij, delta_f, delta_t, filter_len, dof):
    """
    Approximation of unwhitened sum of squares signal energy in a given EP tile.
    See T1200125 for equation number reference.

    Parameters
    ----------
    z_j_b      : time frequency map block which the constructed tile covers
    uw_ss_ii   : unwhitened filter inner products
    uw_ss_ij   : unwhitened adjacent filter inner products
    w_ss_ij    : whitened adjacent filter inner products
    delta_f    : frequency binning of EP filters
    delta_t    : native time resolution of the time frequency map
    filter_len : number of samples in a fitler
    dof        : degrees of freedom in the tile (twice the time-frequency area)
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

def uw_sum_sq(filter1, filter2, spec_corr, psd):
    # < s^2_j(f_1, b) > = 1 / 2 / N * \delta_t EPIP{\Theta, \Theta; P}
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

def make_indp_tiles(tf_map, nc_sum, mu_sq):
    """
    Create a time frequency map with resolution of tf_map binning
    divided by nc_sum + 1. All tiles will be independent up to
    overlap from the original tiling. The mu_sq is applied to the
    resulting addition to normalize the outputs to be zero-mean
    unit-variance Gaussian variables (if the input is Gaussian).
    
    Notes
    -----
    Optimization plan: If we keep the summed complex TF plane in known
    indices, we can save ourselves individual sums at wider frequency
    resolutions.
    Caveats:
      1.  We have to keep track of where we're storing things
      2.  We have to do it from the finest resolution (for *all* t0s) and work our way up
    In the end, I think this is a Haar wavelet transform. Look into it
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

def construct_tiles(nc_sum,mu_sq,band,ts_data,tf_map,psd_segment_length,window_fraction,max_duration):
    """
    Constructing tile and calculate their energy
    """
    # Clip the boundaries to remove window corruption
    clip_samples = int(psd_segment_length * window_fraction * ts_data.sample_rate / 2)
    print "|--- Constructing tile with %d summed channels..." % (nc_sum+1)
    # Current bandwidth of the time-frequency map tiles
    df = band * (nc_sum + 1)
    # How much each "step" is in the time domain -- under sampling rate
    dt = 1.0 / (2 * df)
    us_rate = int(round(dt / ts_data.delta_t))
    print "|--- Undersampling rate for this level: %f" % (ts_data.sample_rate/us_rate)
    print "|--- Calculating tiles..."
    if clip_samples > 0: # because [0:-0] does not give the full array
        tiles = make_indp_tiles(tf_map[:,clip_samples:-clip_samples:us_rate], nc_sum, mu_sq)
    else:
        tiles = make_indp_tiles(tf_map[:,::us_rate], nc_sum, mu_sq)
    print "|--- TF-plane is %dx%s samples" % tiles.shape
    print "|--- Tile energy mean %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))
    if max_duration is not None:
        max_dof = 2 * max_duration * (band * (nc_sum+1))
    else:
        max_dof = 32
    assert max_dof >= 2
    return max_dof, tiles, us_rate, dt, df

def create_tile_duration(j,df,duration,tiles):
    # Duration is fixed by the NDOF and bandwidth
    duration = j / 2.0 / df
    print "|----- Explore signal duration of %f s..." % duration
    print "|----- Summing DOF = %d ..." % (2*j)
    tlen = tiles.shape[1] - 2*j + 1 + 1
    dof_tiles = numpy.zeros((tiles.shape[0], tlen))
    sum_filter = numpy.array([1,0] * (j-1) + [1])
    for f in range(tiles.shape[0]):
        # Sum and drop correlate tiles
        dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
    print "|----- Summed tile energy mean: %f, var %f" % (numpy.mean(dof_tiles), numpy.var(dof_tiles))
    return dof_tiles
    
def create_xml(ts_data,psd_segment_length,window_fraction,event_list,station,setname="MagneticFields"):
    __program__ = 'pyburst_excesspower'
    start_time = LIGOTimeGPS(int(ts_data.start_time))
    end_time = LIGOTimeGPS(int(ts_data.end_time))
    inseg = segment(start_time,end_time)
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    ifo = 'H1'#channel_name.split(":")[0]
    straindict = psd.insert_psd_option_group.__dict__
    proc_row = register_to_xmldoc(xmldoc, __program__,straindict, ifos=[ifo],version=glue.git_version.id, cvs_repository=glue.git_version.branch, cvs_entry_time=glue.git_version.date)
    outseg = determine_output_segment(inseg, psd_segment_length, ts_data.sample_rate, window_fraction)
    ss = append_search_summary(xmldoc, proc_row, ifos=(station,), inseg=inseg, outseg=outseg)
    for sb in event_list:
        sb.process_id = proc_row.process_id
        sb.search = proc_row.program
        sb.ifo, sb.channel = station, setname
    xmldoc.childNodes[0].appendChild(event_list)
    fname = make_filename(station, inseg)
    utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"))

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
