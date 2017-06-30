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

def excess_power2(ts_data,            # Time series from magnetic field data
                  psd_segment_length, # Length of each segment in seconds
                  psd_segment_stride, # Separation between 2 consecutive segments in seconds
                  psd_estimation,     # Average method
                  window_fraction,    # Withening window fraction
                  tile_fap,           # Tile false alarm probability threshold in Gaussian noise.
                  station,            # Station
                  nchans=None,        # Total number of channels
                  band=None,          # Channel bandwidth
                  fmin=0,             # Lowest frequency of the filter bank.
                  fmax=None,          # Highest frequency of the filter bank.
                  max_duration=None,  # Maximum duration of the tile
                  wtype='tukey'):     # Whitening type, can tukey or hann
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
            nchans = int(data_band / band) - 1
        elif nchans is not None:
            # Define filter bandwidth 
            band = data_band / nchans
            nchans = nchans - 1
        # Check if number of channels is superior than unity
        assert nchans > 1
    # Print segment information
    print '|- Estimating PSD from segments of time',
    print '%.2f s in length, with %.2f s stride...'%(psd_segment_length, psd_segment_stride)
    # Convert time series as array of float
    data = ts_data.astype(numpy.float64)
    # Define segment length for PSD estimation in sample unit
    seg_len = int(psd_segment_length * sample_rate)
    # Define separation between consecutive segments in sample unit
    seg_stride = int(psd_segment_stride * sample_rate)
    # Calculate the overall PSD from individual PSD segments
    fd_psd = psd.welch(data,avg_method=psd_estimation,seg_len=seg_len,seg_stride=seg_stride)
    # We need this for the SWIG functions...
    lal_psd = fd_psd.lal()
    # Plot the power spectral density
    plot_spectrum(fd_psd)
    # Create whitening window
    print "|- Whitening window and spectral correlation..."
    if wtype == 'hann':
        window = lal.CreateHannREAL8Window(seg_len)
    elif wtype == 'tukey':
        window = lal.CreateTukeyREAL8Window(seg_len, window_fraction)
    else:
        raise ValueError("Can't handle window type %s" % wtype)
    # Create FFT plan
    fft_plan = lal.CreateForwardREAL8FFTPlan(len(window.data.data), 1)
    # Perform two point spectral correlation
    spec_corr = lal.REAL8WindowTwoPointSpectralCorrelation(window, fft_plan)
    # Initialise filter bank 
    print "|- Create filter..."
    filter_bank, fdb = [], []
    # Loop for each channels
    for i in range(nchans):
        channel_flow = fmin+band/2 + i*band
        channel_width = band
        # Create excess power filter
        lal_filter = lalburst.CreateExcessPowerFilter(channel_flow, channel_width, lal_psd, spec_corr)
        filter_bank.append(lal_filter)
        fdb.append(Spectrum.from_lal(lal_filter))
    # Calculate the minimum bandwidth
    min_band = (len(filter_bank[0].data.data)-1) * filter_bank[0].deltaF / 2
    # Plot filter bank
    plot_bank(fdb)
    # Convert filter bank from frequency to time domain
    print "|- Convert all the frequency domain to the time domain..."
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
    # Compute the renormalization for the base filters up to a given bandwidth.
    mu_sq_dict = {}
    # Loop through powers of 2 up to number of channels
    for nc_sum in range(0, int(math.log(nchans, 2))):
        nc_sum = 2**nc_sum - 1
        print "|- Calculating renormalization for resolution level containing %d %fHz channels"%(nc_sum+1,min_band)
        mu_sq = (nc_sum+1)*numpy.array([lalburst.ExcessPowerFilterInnerProduct(f, f, spec_corr, None)
                                        for f in filter_bank])
        # Uncomment to get all possible frequency renormalizations
        #for n in xrange(nc_sum, nchans): # channel position index
        for n in xrange(nc_sum, nchans, nc_sum+1): # channel position index
            for k in xrange(0, nc_sum): # channel sum index
                # FIXME: We've precomputed this, so use it instead
                mu_sq[n] += 2*lalburst.ExcessPowerFilterInnerProduct(filter_bank[n-k],
                                                                     filter_bank[n-1-k],
                                                                     spec_corr,None)
        #print mu_sq[nc_sum::nc_sum+1]
        mu_sq_dict[nc_sum] = mu_sq
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
    while t_idx_max <= len(ts_data):        
        # Define starting and ending time of the segment in seconds
        start_time = ts_data.start_time + t_idx_min/float(ts_data.sample_rate)
        end_time   = ts_data.start_time + t_idx_max/float(ts_data.sample_rate)
        print "\n|-- Analyzing block %i to %i (%.2f percent)"%(start_time,end_time,100*float(t_idx_max)/len(ts_data))
        # Model a withen time series for the block
        tmp_ts_data = types.TimeSeries(ts_data[t_idx_min:t_idx_max]*window.data.data,
                                       delta_t=1./ts_data.sample_rate,epoch=start_time)
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
        print "|-- Create time-frequency plane for current block"
        # Return the complex snr, along with its associated normalization of the template,
        # matched filtered against the data
        #filter.matched_filter_core(types.FrequencySeries(tmp_filter_bank,delta_f=fd_psd.delta_f),
        #                           fs_data,h_norm=1,psd=fd_psd,low_frequency_cutoff=filter_bank[0].f0,
        #                           high_frequency_cutoff=filter_bank[0].f0+2*band)
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
            filtered_series = filter.matched_filter_core(template,fs_data,h_norm=None,psd=None,
                                                         low_frequency_cutoff=filter_bank[i].f0,
                                                         high_frequency_cutoff=filter_bank[i].f0+2*band)
            # Include filtered series in the map
            tf_map[i,:] = filtered_series[0].numpy()
        # Plot spectrogram
        plot_spectrogram(numpy.abs(tf_map).T,tmp_ts_data.delta_t,band,ts_data.sample_rate,start_time,
                         end_time,fname='segments/time-frequency/%i-%i.png'%(start_time,end_time))
        # Loop through all summed channels
        for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]:
            nc_sum = 2**nc_sum - 1
            mu_sq = mu_sq_dict[nc_sum]
            # Clip the boundaries to remove window corruption
            clip_samples = int(psd_segment_length * window_fraction * ts_data.sample_rate / 2)
            # Constructing tile and calculate their energy
            print "\n|--- Constructing tile with %d summed channels..." % (nc_sum+1)
            # Current bandwidth of the time-frequency map tiles
            df = band * (nc_sum + 1)
            dt = 1.0 / (2 * df)
            # How much each "step" is in the time domain -- under sampling rate
            us_rate = int(round(dt / ts_data.delta_t))
            print "|--- Undersampling rate for this level: %f" % (ts_data.sample_rate/us_rate)
            print "|--- Calculating tiles..."
            # Making independent tiles
            # because [0:-0] does not give the full array
            tf_map_temp = tf_map[:,clip_samples:-clip_samples:us_rate] \
                          if clip_samples > 0 else tf_map[:,::us_rate]
            tiles = tf_map_temp.copy()
            # Here's the deal: we're going to keep only the valid output and
            # it's *always* going to exist in the lowest available indices
            stride = nc_sum + 1
            for i in xrange(tiles.shape[0]/stride):
                numpy.absolute(tiles[stride*i:stride*(i+1)].sum(axis=0), tiles[stride*(i+1)-1])
            tiles = tiles[nc_sum::nc_sum+1].real**2 / mu_sq[nc_sum::nc_sum+1].reshape(-1, 1)
            print "|--- TF-plane is %dx%s samples" % tiles.shape
            print "|--- Tile energy mean %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))
            # Define maximum number of degrees of freedom and check it larger or equal to 2
            max_dof = 32 if max_duration==None else 2 * max_duration * df
            assert max_dof >= 2
            # Loop through multiple degrees of freedom
            for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
                # Duration is fixed by the NDOF and bandwidth
                duration = j * dt
                print "\n|----- Explore signal duration of %f s..." % duration
                print "|----- Summing DOF = %d ..." % (2*j)
                tlen = tiles.shape[1] - 2*j + 1 + 1
                dof_tiles = numpy.zeros((tiles.shape[0], tlen))
                sum_filter = numpy.array([1,0] * (j-1) + [1])
                for f in range(tiles.shape[0]):
                    # Sum and drop correlate tiles
                    dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
                print "|----- Summed tile energy mean: %f, var %f" % (numpy.mean(dof_tiles), numpy.var(dof_tiles))
                plot_spectrogram(dof_tiles.T,dt,df,ts_data.sample_rate,start_time,end_time,
                                 fname='segments/%i-%i/tf_%02ichans_%02idof.png'%(start_time,end_time,nc_sum+1,2*j))
                threshold = scipy.stats.chi2.isf(tile_fap, j)
                print "|------ Threshold for this level: %f" % threshold
                spant, spanf = dof_tiles.shape[1] * dt, dof_tiles.shape[0] * df
                print "|------ Processing %.2fx%.2f time-frequency map." % (spant, spanf)
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
                    event.central_freq = filter_bank[0].f0 + band/2 + i * df + 0.5 * df
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
                    event.amplitude = 0
                print "|------ Total number of events: %d" % len(event_list)
        t_idx_min += int(seg_len * (1 - window_fraction))
        t_idx_max += int(seg_len * (1 - window_fraction))
    setname="MagneticFields"
    __program__ = 'pyburst_excesspower'
    start_time = LIGOTimeGPS(int(ts_data.start_time))
    end_time = LIGOTimeGPS(int(ts_data.end_time))
    inseg = segment(start_time,end_time)
    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    ifo = 'H1'#channel_name.split(":")[0]
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
    fname = 'excesspower.xml.gz'
    utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"))
