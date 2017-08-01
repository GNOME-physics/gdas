import time,math,os,scipy,lal,lalburst
from glue                             import git_version
from glue.lal                         import LIGOTimeGPS
from glue.ligolw                      import lsctables,ligolw,utils
from glue.ligolw.utils.search_summary import append_search_summary
from glue.ligolw.utils.process        import register_to_xmldoc
from glue.segments                    import segment
from gwpy.spectrum                    import Spectrum
from gwpy.timeseries                  import TimeSeries
from plots                            import *
from pycbc                            import psd,types,filter
from scipy.signal                     import fftconvolve
from utils                            import *

def excess_power(ts_data,                      # Time series from magnetic field data 
                 band=None,                    # Channel bandwidth
                 channel_name='channel-name',  # Channel name
                 fmin=0,                       # Lowest frequency of the filter bank.
                 fmax=None,                    # Highest frequency of the filter bank.
                 impulse=False,                # Impulse response
                 make_plot=True,               # Condition to produce plots
                 max_duration=None,            # Maximum duration of the tile
                 nchans=256,                   # Total number of channels
                 psd_estimation='median-mean', # Average method
                 psd_segment_length=60,        # Length of each segment in seconds
                 psd_segment_stride=30,        # Separation between 2 consecutive segments in seconds
                 station='station-name',       # Station name
                 tile_fap=1e-7,                # Tile false alarm probability threshold in Gaussian noise.
                 verbose=True,                 # Print details
                 window_fraction=0,            # Withening window fraction
                 wtype='tukey'):               # Whitening type, can tukey or hann
    """
    Perform excess-power search analysis on magnetic field data.
    This method will produce a bunch of time-frequency plots for every
    tile duration and bandwidth analysed as well as a XML file identifying
    all the triggers found in the selected data within the user-defined
    time range.

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
      Channel bandwidth
    fmin : float
      Lowest frequency of the filter bank.
    fmax : float
      Highest frequency of the filter bank
    """
    # Determine sampling rate based on extracted time series
    sample_rate = ts_data.sample_rate
    # Check if tile maximum frequency is not defined
    if fmax is None or fmax>sample_rate/2.:
        # Set the tile maximum frequency equal to the Nyquist frequency
        # (i.e. half the sampling rate)
        fmax = sample_rate / 2.0
    # Check whether or not tile bandwidth and channel are defined
    if band is None and nchans is None:
        # Exit program with error message
        exit("Either bandwidth or number of channels must be specified...")
    else:
        # Check if tile maximum frequency larger than its minimum frequency
        assert fmax >= fmin
        # Define spectral band of data
        data_band = fmax - fmin
        # Check whether tile bandwidth or channel is defined
        if band is not None:
            # Define number of possible filter bands
            nchans = int(data_band / band)
        elif nchans is not None:
            # Define filter bandwidth
            band = data_band / nchans
            nchans -= 1
        # Check if number of channels is superior than unity
        assert nchans > 1
    # Print segment information
    if verbose: print '|- Estimating PSD from segments of',
    if verbose: print '%.2f s, with %.2f s stride...'%(psd_segment_length, psd_segment_stride)
    # Convert time series as array of float
    data = ts_data.astype(numpy.float64)
    # Define segment length for PSD estimation in sample unit
    seg_len = int(psd_segment_length * sample_rate)
    # Define separation between consecutive segments in sample unit
    seg_stride = int(psd_segment_stride * sample_rate)
    # Minimum frequency of detectable signal in a segment
    delta_f = 1. / psd_segment_length
    # Calculate PSD length counting the zero frequency element
    fd_len = fmax / delta_f + 1
    # Calculate the overall PSD from individual PSD segments
    if impulse:
        # Produce flat data
        flat_data = numpy.ones(int(fd_len)) * 2. / fd_len
        # Create PSD frequency series
        fd_psd = types.FrequencySeries(flat_data, 1./psd_segment_length, ts_data.start_time)
    else:
        # Create overall PSD using Welch's method
        fd_psd = psd.welch(data,avg_method=psd_estimation,seg_len=seg_len,seg_stride=seg_stride)
    if make_plot:
        # Plot the power spectral density
        plot_spectrum(fd_psd)
    # We need this for the SWIG functions
    lal_psd = fd_psd.lal()
    # Create whitening window
    if verbose: print "|- Whitening window and spectral correlation..."
    if wtype == 'hann': window = lal.CreateHannREAL8Window(seg_len)
    elif wtype == 'tukey': window = lal.CreateTukeyREAL8Window(seg_len, window_fraction)
    else: raise ValueError("Can't handle window type %s" % wtype)
    # Create FFT plan
    fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)
    # Perform two point spectral correlation
    spec_corr = lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)
    # Determine length of individual filters
    filter_length = int(2*band/fd_psd.delta_f)+1
    # Initialise filter bank
    if verbose: print "|- Create bank of %i filters of %i Hz bandwidth..."%(nchans,filter_length)
    # Initialise array to store filter's frequency series and metadata
    lal_filters = []
    # Initialise array to store filter's time series
    fdb = []
    # Loop over the channels
    for i in range(nchans):
        # Define central position of the filter
        freq = fmin + band/2 + i*band
        # Create excess power filter
        lal_filter = lalburst.CreateExcessPowerFilter(freq, band, lal_psd, spec_corr)
        # Testing spectral correlation on filter
        #print lalburst.ExcessPowerFilterInnerProduct(lal_filter, lal_filter, spec_corr, None)
        # Append entire filter structure
        lal_filters.append(lal_filter)
        # Append filter's spectrum
        fdb.append(Spectrum.from_lal(lal_filter))
    if make_plot:
        # Plot filter bank
        plot_bank(fdb)
        # Convert filter bank from frequency to time domain
        if verbose: print "|- Convert all the frequency domain to the time domain..."
        tdb = []
        # Loop for each filter's spectrum
        for fdt in fdb:
            zero_padded = numpy.zeros(int((fdt.f0 / fdt.df).value) + len(fdt))
            st = int((fdt.f0 / fdt.df).value)
            zero_padded[st:st+len(fdt)] = numpy.real_if_close(fdt.value)
            n_freq = int(sample_rate / 2 / fdt.df.value) * 2
            tdt = numpy.fft.irfft(zero_padded, n_freq) * math.sqrt(sample_rate)
            tdt = numpy.roll(tdt, len(tdt)/2)
            tdt = TimeSeries(tdt, name="", epoch=fdt.epoch, sample_rate=sample_rate)
            tdb.append(tdt)
        # Plot time series filter
        plot_filters(tdb,fmin,band)
    # Computer whitened inner products of input filters with themselves
    white_filter_ip = numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None) for f in lal_filters])
    # Computer unwhitened inner products of input filters with themselves
    unwhite_filter_ip = numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, lal_psd) for f in lal_filters])
    # Computer whitened filter inner products between input adjacent filters
    white_ss_ip = numpy.array([lalburst.ExcessPowerFilterInnerProduct(f1, f2, spec_corr, None) for f1, f2 in zip(lal_filters[:-1], lal_filters[1:])])
    # Computer unwhitened filter inner products between input adjacent filters
    unwhite_ss_ip = numpy.array([lalburst.ExcessPowerFilterInnerProduct(f1, f2, spec_corr, lal_psd) for f1, f2 in zip(lal_filters[:-1], lal_filters[1:])])
    # Check filter's bandwidth is equal to user defined channel bandwidth
    min_band = (len(lal_filters[0].data.data)-1) * lal_filters[0].deltaF / 2
    assert min_band==band
    # Create an event list where all the triggers will be stored
    event_list = lsctables.New(lsctables.SnglBurstTable,
                               ['start_time','start_time_ns','peak_time','peak_time_ns',
                                'duration','bandwidth','central_freq','chisq_dof',
                                'confidence','snr','amplitude','channel','ifo',
                                'process_id','event_id','search','stop_time','stop_time_ns'])
    # Create repositories to save TF and time series plots
    os.system('mkdir -p segments/time-frequency')
    os.system('mkdir -p segments/time-series')
    # Define time edges
    t_idx_min, t_idx_max = 0, seg_len
    # Loop over each segment
    while t_idx_max <= len(ts_data):        
        # Define first and last timestamps of the block
        start_time = ts_data.start_time + t_idx_min/float(ts_data.sample_rate)
        end_time   = ts_data.start_time + t_idx_max/float(ts_data.sample_rate)
        if verbose: print "\n|- Analyzing block %i to %i (%.2f percent)"%(start_time,end_time,100*float(t_idx_max)/len(ts_data))
        # Debug for impulse response
        if impulse:
            for i in range(t_idx_min, t_idx_max):
                ts_data[i] = 1000. if i == (t_idx_max + t_idx_min)/2 else 0.
        # Model a withen time series for the block
        tmp_ts_data = types.TimeSeries(ts_data[t_idx_min:t_idx_max]*window.data.data,
                                       delta_t=1./ts_data.sample_rate,epoch=start_time)
        # Save time series in relevant repository
        segfolder = 'segments/%i-%i'%(start_time,end_time)
        os.system('mkdir -p '+segfolder)
        if make_plot:
            # Plot time series
            plot_ts(tmp_ts_data,fname='segments/time-series/%i-%i.png'%(start_time,end_time))
        # Convert times series to frequency series
        fs_data = tmp_ts_data.to_frequencyseries()
        if verbose: print "|- Frequency series data has variance: %s" % fs_data.data.std()**2
        # Whitening (FIXME: Whiten the filters, not the data)
        fs_data.data /= numpy.sqrt(fd_psd) / numpy.sqrt(2 * fd_psd.delta_f)
        if verbose: print "|- Whitened frequency series data has variance: %s" % fs_data.data.std()**2
        if verbose: print "|- Create time-frequency plane for current block"
        # Return the complex snr, along with its associated normalization of the template,
        # matched filtered against the data
        #filter.matched_filter_core(types.FrequencySeries(tmp_filter_bank,delta_f=fd_psd.delta_f),
        #                           fs_data,h_norm=1,psd=fd_psd,low_frequency_cutoff=lal_filters[0].f0,
        #                           high_frequency_cutoff=lal_filters[0].f0+2*band)
        if verbose: print "|- Filtering all %d channels...\n" % nchans,
        # Initialise 2D zero array
        tmp_filter_bank = numpy.zeros(len(fd_psd), dtype=numpy.complex128)
        # Initialise 2D zero array for time-frequency map
        tf_map = numpy.zeros((nchans, seg_len), dtype=numpy.complex128)
        # Loop over all the channels
        for i in range(nchans):
            # Reset filter bank series
            tmp_filter_bank *= 0.0
            # Index of starting frequency
            f1 = int(lal_filters[i].f0/fd_psd.delta_f)
            # Index of last frequency bin
            f2 = int((lal_filters[i].f0 + 2*band)/fd_psd.delta_f)+1
            # (FIXME: Why is there a factor of 2 here?)
            tmp_filter_bank[f1:f2] = lal_filters[i].data.data * 2
            # Define the template to filter the frequency series with
            template = types.FrequencySeries(tmp_filter_bank, delta_f=fd_psd.delta_f, copy=False)
            # Create filtered series
            filtered_series = filter.matched_filter_core(template,fs_data,h_norm=None,psd=None,
                                                         low_frequency_cutoff=lal_filters[i].f0,
                                                         high_frequency_cutoff=lal_filters[i].f0+2*band)
            # Include filtered series in the map
            tf_map[i,:] = filtered_series[0].numpy()
        if make_plot:
            # Plot spectrogram
            plot_spectrogram(numpy.abs(tf_map).T,tmp_ts_data.delta_t,band,ts_data.sample_rate,start_time,
                             end_time,fname='segments/time-frequency/%i-%i.png'%(start_time,end_time))
        # Loop through powers of 2 up to number of channels
        for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]:
            # Calculate total number of summed channels
            nc_sum = 2**nc_sum
            if verbose: print "\n\t|- Contructing tiles containing %d narrow band channels"%nc_sum
            # Compute full bandwidth of virtual channel
            df = band * nc_sum
            # Compute minimal signal's duration in virtual channel
            dt = 1.0 / (2 * df)
            # Compute under sampling rate
            us_rate = int(round(dt / ts_data.delta_t))
            if verbose: print "\t|- Undersampling rate for this level: %f" % (ts_data.sample_rate/us_rate)
            if verbose: print "\t|- Calculating tiles..."
            # Clip the boundaries to remove window corruption
            clip_samples = int(psd_segment_length * window_fraction * ts_data.sample_rate / 2)
            # Undersample narrow band channel's time series 
            # Apply clipping condition because [0:-0] does not give the full array
            tf_map_temp = tf_map[:,clip_samples:-clip_samples:us_rate] \
                          if clip_samples > 0 else tf_map[:,::us_rate]
            # Initialise final tile time-frequency map
            tiles = numpy.zeros(((nchans+1)/nc_sum,tf_map_temp.shape[1]))
            # Loop over tile index
            for i in xrange(len(tiles)):
                # Sum all inner narrow band channels
                ts_tile = numpy.absolute(tf_map_temp[nc_sum*i:nc_sum*(i+1)].sum(axis=0))
                # Define index of last narrow band channel for given tile
                n = (i+1)*nc_sum - 1
                n = n-1 if n==len(lal_filters) else n
                # Computer withened inner products of each input filter with itself
                mu_sq = nc_sum*lalburst.ExcessPowerFilterInnerProduct(lal_filters[n], lal_filters[n], spec_corr, None)
                #kmax = nc_sum-1 if n==len(lal_filters) else nc_sum-2
                # Loop over the inner narrow band channels
                for k in xrange(0, nc_sum-1):
                    # Computer whitened filter inner products between input adjacent filters
                    mu_sq += 2*lalburst.ExcessPowerFilterInnerProduct(lal_filters[n-k],lal_filters[n-1-k],spec_corr,None)
                # Normalise tile's time series
                tiles[i] = ts_tile.real**2 / mu_sq
            if verbose: print "\t|- TF-plane is %dx%s samples" % tiles.shape
            if verbose: print "\t|- Tile energy mean %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))
            # Define maximum number of degrees of freedom and check it larger or equal to 2
            max_dof = 32 if max_duration==None else int(max_duration / dt)
            assert max_dof >= 2
            # Loop through multiple degrees of freedom
            for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
                # Duration is fixed by the NDOF and bandwidth
                duration = j * dt
                if verbose: print "\n\t\t|- Summing DOF = %d ..." % (2*j)
                if verbose: print "\t\t|- Explore signal duration of %f s..." % duration
                # Construct filter
                sum_filter = numpy.array([1,0] * (j-1) + [1])
                if verbose: print sample_rate,dt,df,us_rate,j
                # Calculate length of filtered time series
                tlen = tiles.shape[1] - sum_filter.shape[0] + 1
                # Initialise filtered time series array
                dof_tiles = numpy.zeros((tiles.shape[0], tlen))
                # Loop over tiles
                for f in range(tiles.shape[0]):
                    # Sum and drop correlate tiles
                    dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
                if verbose: print "\t\t|- Summed tile energy mean: %f" % (numpy.mean(dof_tiles))
                if verbose: print "\t\t|- Variance tile energy: %f" % (numpy.var(dof_tiles))
                if make_plot:
                    plot_spectrogram(dof_tiles.T,dt,df,ts_data.sample_rate,start_time,end_time,
                                     fname='segments/%i-%i/tf_%02ichans_%02idof.png'%(start_time,end_time,nc_sum,2*j))
                threshold = scipy.stats.chi2.isf(tile_fap, j)
                if verbose: print "\t\t|- Threshold for this level: %f" % threshold
                spant, spanf = dof_tiles.shape[1] * dt, dof_tiles.shape[0] * df
                if verbose: print "\t\t|- Processing %.2fx%.2f time-frequency map." % (spant, spanf)
                # Since we clip the data, the start time needs to be adjusted accordingly
                window_offset_epoch = fs_data.epoch + psd_segment_length * window_fraction / 2
                window_offset_epoch = LIGOTimeGPS(float(window_offset_epoch))
                for i, j in zip(*numpy.where(dof_tiles > threshold)):
                    event = event_list.RowType()
                    # The points are summed forward in time and thus a `summed point' is the
                    # sum of the previous N points. If this point is above threshold, it
                    # corresponds to a tile which spans the previous N points. However, the
                    # 0th point (due to the convolution specifier 'valid') is actually
                    # already a duration from the start time. All of this means, the +
                    # duration and the - duration cancels, and the tile 'start' is, by
                    # definition, the start of the time frequency map if j = 0
                    # FIXME: I think this needs a + dt/2 to center the tile properly
                    event.set_start(window_offset_epoch + float(j * dt))
                    event.set_stop(window_offset_epoch + float(j * dt) + duration)
                    event.set_peak(event.get_start() + duration / 2)
                    event.central_freq = lal_filters[0].f0 + band/2 + i * df + 0.5 * df
                    event.duration = duration
                    event.bandwidth = df
                    event.chisq_dof = 2 * duration * df
                    event.snr = math.sqrt(dof_tiles[i,j] / event.chisq_dof - 1)
                    # FIXME: Magic number 0.62 should be determine empircally
                    event.confidence = -lal.LogChisqCCDF(event.snr * 0.62, event.chisq_dof * 0.62)
                    event.amplitude = None
                    event.process_id = None
                    event.event_id = event_list.get_next_id()
                    event_list.append(event)
                for event in event_list[::-1]:
                    if event.amplitude != None:
                        continue
                    etime_min_idx = float(event.get_start()) - float(fs_data.epoch)
                    etime_min_idx = int(etime_min_idx / tmp_ts_data.delta_t)
                    etime_max_idx = float(event.get_start()) - float(fs_data.epoch) + event.duration
                    etime_max_idx = int(etime_max_idx / tmp_ts_data.delta_t)
                    # (band / 2) to account for sin^2 wings from finest filters
                    flow_idx = int((event.central_freq - event.bandwidth / 2 - (df / 2) - fmin) / df)
                    fhigh_idx = int((event.central_freq + event.bandwidth / 2 + (df / 2) - fmin) / df)
                    # TODO: Check that the undersampling rate is always commensurate
                    # with the indexing: that is to say that
                    # mod(etime_min_idx, us_rate) == 0 always
                    z_j_b = tf_map[flow_idx:fhigh_idx,etime_min_idx:etime_max_idx:us_rate]
                    # FIXME: Deal with negative hrss^2 -- e.g. remove the event
                    try:
                        event.amplitude = measure_hrss(z_j_b, unwhite_filter_ip[flow_idx:fhigh_idx],
                                                       unwhite_ss_ip[flow_idx:fhigh_idx-1],
                                                       white_ss_ip[flow_idx:fhigh_idx-1],
                                                       fd_psd.delta_f, tmp_ts_data.delta_t,
                                                       len(lal_filters[0].data.data), event.chisq_dof)
                    except ValueError:
                        event.amplitude = 0
                if verbose: print "\t\t|- Total number of events: %d" % len(event_list)
        t_idx_min += int(seg_len * (1 - window_fraction))
        t_idx_max += int(seg_len * (1 - window_fraction))
    setname="MagneticFields"
    __program__ = 'pyburst_excesspower_gnome'
    start_time = LIGOTimeGPS(int(ts_data.start_time))
    end_time = LIGOTimeGPS(int(ts_data.end_time))
    inseg = segment(start_time,end_time)
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    ifo = channel_name.split(":")[0]
    straindict = psd.insert_psd_option_group.__dict__
    proc_row = register_to_xmldoc(xmldoc, __program__,straindict, ifos=[ifo],version=git_version.id,
                                  cvs_repository=git_version.branch, cvs_entry_time=git_version.date)
    dt_stride = psd_segment_length
    sample_rate = ts_data.sample_rate
    # Amount to overlap successive blocks so as not to lose data
    window_overlap_samples = window_fraction * sample_rate
    outseg = inseg.contract(window_fraction * dt_stride / 2)
    # With a given dt_stride, we cannot process the remainder of this data
    remainder = math.fmod(abs(outseg), dt_stride * (1 - window_fraction))
    # ...so make an accounting of it
    outseg = segment(outseg[0], outseg[1] - remainder)
    ss = append_search_summary(xmldoc, proc_row, ifos=(station,), inseg=inseg, outseg=outseg)
    for sb in event_list:
        sb.process_id = proc_row.process_id
        sb.search = proc_row.program
        sb.ifo, sb.channel = station, setname
    xmldoc.childNodes[0].appendChild(event_list)
    ifostr = ifo if isinstance(ifo, str) else "".join(ifo)
    st_rnd, end_rnd = int(math.floor(inseg[0])), int(math.ceil(inseg[1]))
    dur = end_rnd - st_rnd    
    fname = "%s-excesspower-%d-%d.xml.gz" % (ifostr, st_rnd, dur)
    utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"))

def measure_hrss(z_j_b, uw_ss_ii, uw_ss_ij, w_ss_ij, delta_f, delta_t, filter_len, dof):
    """
    Approximation of unwhitened sum of squares signal energy in a given EP tile. See T1200125 for equation number reference.
    z_j_b - time frequency map block which the constructed tile covers
    uw_ss_ii - unwhitened filter inner products
    uw_ss_ij - unwhitened adjacent filter inner products
    w_ss_ij - whitened adjacent filter inner products
    delta_f - frequency binning of EP filters
    delta_t - native time resolution of the time frequency map
    filter_len - number of samples in a fitler
    dof - degrees of freedom in the tile (twice the time-frequency area)
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
    Approximation of unwhitened sum of squares signal energy in a given EP tile. See T1200125 for equation number reference.

    NOTE: This function is deprecated in favor of measure_hrss, since it requires recomputation of many inner products, making it particularly slow.
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

