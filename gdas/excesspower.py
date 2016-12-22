from utils import *
from plots import *

class excesspower(object):

    def __init__(self):
    
        #======================================================================
        # Estimate Power Spectral Density (PSD)
        #======================================================================
        
        print tprint(t0),"Estimating PSD from segments of time %.2f s in length, with %.2f s stride..."  % (args.psd_segment_length, args.psd_segment_stride)
        # Convert time series as array of float
        data = ts_data.astype(numpy.float64)
        # Average method to measure PSD from the data
        avg_method = args.psd_estimation
        # The segment length for PSD estimation in samples
        seg_len = int(args.psd_segment_length * args.sample_rate)
        # The separation between consecutive segments in samples
        seg_stride = int(args.psd_segment_stride * args.sample_rate)
        # Lifted from the psd.from_cli module
        fd_psd = psd.welch(data,avg_method=avg_method,seg_len=seg_len,seg_stride=seg_stride)
        # Plot the power spectral density
        plot_spectrum(fd_psd)
        # We need this for the SWIG functions
        lal_psd = fd_psd.lal()
    
        #======================================================================
        # Whitening window and spectral correlation
        #======================================================================
        
        print tprint(t0),"Determine how much data on either side of the tukey window to discard..."
        # NOTE: Nominally, this means that one will lose window_fraction * args.psd_segment_length to
        # corruption from the window: this is simply discarded
        # NOTE 2: This is tuned to give an integer offset when used with args.psd_segment_length = 8,
        # smaller windows will have fractions of integers, but larger powers of two will
        # still preseve this --- probably not a big deal in the end
        # Do two point spectral correlation
        window, spec_corr = calculate_spectral_correlation(seg_len,'tukey',window_fraction=args.window_fraction)
        window = window.data.data
        window_sigma_sq = numpy.mean(window**2)
        # Pre scale the window by its root mean squared -- see eqn 11 of EP document
        #window /= numpy.sqrt(window_sigma_sq)
        
        #======================================================================
        # Filter generation
        #======================================================================
        
        print tprint(t0),"Create filter..."
        # Define filters
        filter_bank, fdb = create_filter_bank(fd_psd.delta_f, flow+band/2, band, nchans, fd_psd, spec_corr,args.min_frequency,args.max_frequency)
        # This is necessary to compute the mu^2 normalizations
        #white_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, None)
        #unwhite_filter_ip = compute_filter_ips_self(filter_bank, spec_corr, lal_psd)
        # These two are needed for the unwhitened mean square sum (hrss)
        #white_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, None)
        #unwhite_ss_ip = compute_filter_ips_adjacent(filter_bank, spec_corr, lal_psd)
        # Print data of first channel's filter
        print filter_bank[0].data.data
        # Print frequency separation between 2 values in the first channel's filter
        print filter_bank[0].deltaF
        # Print first frequency of the first channel's filter
        print filter_bank[0].f0
        # Print last frequency of the first channel's filter
        print filter_bank[0].f0+len(filter_bank[0].data.data)*filter_bank[0].deltaF
        
        #======================================================================
        # Convert all the frequency domain to the time domain and plot them
        #======================================================================
        
        print tprint(t0),"Convert all the frequency domain to the time domain..."
        tdb = []
        for fdt in fdb:
            zero_padded = numpy.zeros(int((fdt.f0 / fdt.df).value) + len(fdt))
            st = int((fdt.f0 / fdt.df).value)
            zero_padded[st:st+len(fdt)] = numpy.real_if_close(fdt.value)
            n_freq = int(args.sample_rate / 2 / fdt.df.value) * 2
            tdt = numpy.fft.irfft(zero_padded, n_freq) * math.sqrt(args.sample_rate)
            tdt = numpy.roll(tdt, len(tdt)/2)
            tdt = TimeSeries(tdt, name="", epoch=fdt.epoch, sample_rate=args.sample_rate)
            tdb.append(tdt)
            
        pyplot.figure()
        for i, fdt in enumerate(fdb[:5]):
            pyplot.plot(fdt.frequencies, fdt, 'k-')
        pyplot.grid()
        pyplot.xlabel("frequency [Hz]")
        pyplot.savefig('bank.png')
        pyplot.close()
        
        pyplot.figure()
        pyplot.subplots_adjust(left=0.2,right=0.95,bottom=0.15,top=0.95,hspace=0,wspace=0)
        for i, tdt in enumerate(tdb[:8:3]):
            ax = pyplot.subplot(3, 1, i+1)
            ax.plot(tdt.times.value - 2., numpy.real_if_close(tdt.value), 'k-')
            c_f = flow+band/2 + 3 * (band*i) + 2.
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("%d Hz" % c_f)
            ax.set_xlim(25.0, 31.0)
            ax.set_ylim([-max(tdt.value), max(tdt.value)])
            if i!=2: pyplot.setp(ax.get_xticklabels(), visible=False)
        pyplot.savefig('filters.png')
        pyplot.close()        
        
        #======================================================================
        # Compute normalization for virtual channel
        #======================================================================
        
        print tprint(t0),"Compute normalization for virtual channel..."
        # Initialise dictionary
        mu_sq_dict = {}
        t1 = time.time()
        # nc_sum additional channel adds
        for nc_sum in range(0, int(math.log(nchans, 2))):
            min_band = (len(filter_bank[0].data.data)-1) * filter_bank[0].deltaF / 2
            print tprint(t0,t1),"Calculation for %d %dHz channels" % (nc_sum+1, min_band)
            nc_sum = 2**nc_sum - 1
            mu_sq_dict[nc_sum] = compute_channel_renormalization(nc_sum, filter_bank, spec_corr, nchans)
        
        #======================================================================
        # Setting up event list and boundaries...
        #======================================================================
        
        print tprint(t0),"Setting up event list and boundaries..."
        # Create event list for single burst table
        event_list = lsctables.New(lsctables.SnglBurstTable,
                                   ['start_time','start_time_ns','peak_time','peak_time_ns',
                                    'duration','bandwidth','central_freq','chisq_dof',
                                    'confidence','snr','amplitude','channel','ifo',
                                    'process_id','event_id','search','stop_time','stop_time_ns'])
        
        # Determine boundaries of stride in time domain
        t_idx_min, t_idx_max = 0, seg_len
        # Check if user requested starting time is defined
        if args.analysis_start_time is not None:
            # Define the time difference in seconds between data and user requested starting times
            t_idx_off = args.analysis_start_time - ts_data.start_time
            # Calculate the index of the user requested starting point in the data
            t_idx_off = int(t_idx_off * args.sample_rate)
        else:
            # Define index of the starting point as first value in data
            t_idx_off = 0
        # Initialise minimum index values as offset starting index
        t_idx_min += t_idx_off
        # Initialise maximum index values as offset starting index
        t_idx_max += t_idx_off
        
        # Check if user requested end time is defined
        if args.analysis_end_time is not None:
            # Define the time difference between data and user requested ending times
            t_idx_max_off = args.analysis_end_time - ts_data.start_time
            # Calculate the index of the user requested starting point in the data
            t_idx_max_off = int(t_idx_max_off * args.sample_rate)
        else:
            # Define index of the ending point as the length of data array
            t_idx_max_off = len(ts_data)
        
        #======================================================================
        # Loop over every segments
        #======================================================================
    
    
        print tprint(t0),"Loop over each segments..."
        os.system('mkdir -p segments/time-frequency')
        os.system('mkdir -p segments/time-series')
        # Loop over each data within the user requested time period
        t1 = time.time()
        while t_idx_max <= t_idx_max_off:
            
            #======================================================================
            # Get frequency series of the current block
            #======================================================================
            
            # Define starting and ending time of the segment in seconds
            start_time = ts_data.start_time + t_idx_min/float(args.sample_rate)
            end_time = ts_data.start_time + t_idx_max/float(args.sample_rate)
            print tprint(t0,t1),"Analyzing block %i to %i (%.2f percent)"%(start_time,end_time,100*float(t_idx_max)/float(t_idx_max_off))
            # Model a withen time series for the block
            tmp_ts_data = types.TimeSeries(ts_data[t_idx_min:t_idx_max]*window,delta_t=1./args.sample_rate,epoch=start_time)
            # Save time series in relevant repository
            segfolder = 'segments/%i-%i'%(start_time,end_time)
            os.system('mkdir -p '+segfolder)
            plot_ts(tmp_ts_data,fname='segments/time-series/%i-%i.png'%(start_time,end_time))
            # Convert times series to frequency series
            fs_data = tmp_ts_data.to_frequencyseries()
            print tprint(t0,t1),"Frequency series data has variance: %s" % fs_data.data.std()**2
            # Whitening (FIXME: Whiten the filters, not the data)
            fs_data.data /= numpy.sqrt(fd_psd) / numpy.sqrt(2 * fd_psd.delta_f)
            print tprint(t0,t1),"Whitened frequency series data has variance: %s" % fs_data.data.std()**2
            
            #======================================================================
            # Create Time-Frequency plane for current block
            #======================================================================
            
            # Return the complex snr, along with its associated normalization of the template,
            # matched filtered against the data
            #filter.matched_filter_core(types.FrequencySeries(tmp_filter_bank,delta_f=fd_psd.delta_f),fs_data,h_norm=1,psd=fd_psd,low_frequency_cutoff=filter_bank[0].f0,high_frequency_cutoff=filter_bank[0].f0+2*band)
            print tprint(t0,t1),"Filtering all %d channels..." % nchans
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
            plot_spectrogram(numpy.abs(tf_map).T,tmp_ts_data.delta_t,band,ts_data.sample_rate,start_time,end_time,fname='segments/time-frequency/%i-%i.png'%(start_time,end_time))
    
            #======================================================================
            # Constructing tiles with different bandwidth
            #======================================================================    
            
            t2 = time.time()
            # Clip the boundaries to remove window corruption
            clip_samples = int(args.psd_segment_length * args.window_fraction * args.sample_rate / 2)
            for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]: # nc_sum additional channel adds
                
                #======================================================================
                # Define frequency width and undersample rate
                #======================================================================
                
                nc_sum = 2**nc_sum - 1
                print tprint(t0,t1),"Constructing tile with %d summed channels..." % (nc_sum+1)
                mu_sq = mu_sq_dict[nc_sum]
                # Current bandwidth of the time-frequency map tiles
                df = band * (nc_sum + 1)
                # How much each "step" is in the time domain -- under sampling rate
                dt = 1.0 / (2 * df)
                us_rate = int(round(dt / ts_data.delta_t))
                print tprint(t0,t1,t2),"Undersampling rate for this level: %f" % (args.sample_rate/us_rate)
                
                #======================================================================
                # Calculate tile energy
                #======================================================================
                
                print tprint(t0,t1,t2),"Calculating tiles..."
                if clip_samples > 0: # because [0:-0] does not give the full array
                    tiles = make_indp_tiles(tf_map[:,clip_samples:-clip_samples:us_rate], nc_sum, mu_sq)
                else:
                    tiles = make_indp_tiles(tf_map[:,::us_rate], nc_sum, mu_sq)
                print tprint(t0,t1,t2),"TF-plane is %dx%s samples" % tiles.shape
                print tprint(t0,t1,t2),"Tile energy mean %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))
                if args.max_duration is not None:
                    max_dof = 2 * args.max_duration * (band * (nc_sum+1))
                else:
                    max_dof = 32
                assert max_dof >= 2
                
                #======================================================================
                # Plot
                #======================================================================
                
                # bins = numpy.linspace(0, 40, 100)
                # cnt = numpy.zeros(bins.shape[0]-1)
                # for i, tdf in enumerate(tdb[:nchans]):
                #     us_rate = int(1.0 / (2 * band*nc_sum * ts_data.dt.value))
                #     pyplot.figure(0, figsize=(10, 10))
                #     pyplot.subplot(nchans, 1, i+1)
                #     white = tmp_ts_data.whiten(64, 32, asd=numpy.sqrt(cdata_psd_tmp), window='boxcar') * args.sample_rate/4
                #     snr_1dof = numpy.convolve(tdf, white, "valid")
                #     # Undersample the data
                #     snr_1dof = snr_1dof[::us_rate]**2
                #     # Sum semi-adjacent samples to get 2 DOF tiles
                #     snr_2dof = numpy.convolve(snr_1dof, numpy.array([1, 0, 1, 0]))
                #     t = TimeSeries(snr_2dof, epoch=white.epoch, sample_rate=int(1.0/(us_rate * tmp_ts_data.dt.value)))
                #     pyplot.plot(t.times + len(tdf)/2 * tdf.dt, snr_2dof, 'k-')
                #     pyplot.axvline(random_time)
                #     tmp, _ = numpy.histogram(snr_2dof, bins=bins)
                #     cnt += tmp
                # plot_spectrogram(dof_tiles.T,fname='%s/tf_%ichans_%02idof.png'%(segfolder,nc_sum+1,2*j))
                # plot.savefig("%s/bands.png"%(segfolder))
                
                #======================================================================
                # Exploring different durations
                #======================================================================
                
                t3 = time.time()
                for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
                    
                    # Duration is fixed by the NDOF and bandwidth
                    duration = j / 2.0 / df
                    print tprint(t0,t1,t2),"Explore signal duration of %f s..." % duration
                    
                    #======================================================================
                    # Creating tile array with different duration
                    #======================================================================
                    
                    print tprint(t0,t1,t2,t3),"Summing DOF = %d ..." % (2*j)
                    tlen = tiles.shape[1] - 2*j + 1 + 1
                    dof_tiles = numpy.zeros((tiles.shape[0], tlen))
                    sum_filter = numpy.array([1,0] * (j-1) + [1])
                    for f in range(tiles.shape[0]):
                        # Sum and drop correlate tiles
                        dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
                    print tprint(t0,t1,t2,t3),"Summed tile energy mean: %f, var %f" % (numpy.mean(dof_tiles), numpy.var(dof_tiles))
                    
                    #======================================================================
                    # Trigger finding
                    #======================================================================
                    
                    threshold = scipy.stats.chi2.isf(args.tile_fap, j)
                    print tprint(t0,t1,t2,t3),"Threshold for this level: %f" % threshold
                    plot_spectrogram(dof_tiles.T,dt,df,ts_data.sample_rate,start_time,end_time,fname='segments/%i-%i/tf_%02ichans_%02idof.png'%(start_time,end_time,nc_sum+1,2*j))
                    spant, spanf = dof_tiles.shape[1] * dt, dof_tiles.shape[0] * df
                    print tprint(t0,t1,t2,t3),"Processing %.2fx%.2f time-frequency map." % (spant, spanf)
                    # Since we clip the data, the start time needs to be adjusted accordingly
                    window_offset_epoch = fs_data.epoch + args.psd_segment_length * args.window_fraction / 2
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
                    print tprint(t0,t1,t2,t3),"Total number of events: %d" % len(event_list)
                    
            t_idx_min += int(seg_len * (1 - args.window_fraction))
            t_idx_max += int(seg_len * (1 - args.window_fraction))
    
        #======================================================================
        # Finalize analysis
        #======================================================================
        
        print tprint(t0),"Create output file..."
        __program__ = 'pyburst_excesspower'
        # Represent the range of values in the semi-open interval
        start_time = LIGOTimeGPS(args.analysis_start_time or args.gps_start_time or int(ts_data.start_time))
        end_time = LIGOTimeGPS(args.analysis_end_time or args.gps_end_time or int(ts_data.end_time))
        inseg = segment(start_time,end_time)
        xmldoc = ligolw.Document()
        xmldoc.appendChild(ligolw.LIGO_LW())
        ifo = args.channel_name.split(":")[0]
        proc_row = register_to_xmldoc(xmldoc, __program__, args.__dict__, ifos=[ifo],version=glue.git_version.id, cvs_repository=glue.git_version.branch, cvs_entry_time=glue.git_version.date)
        # Figure out the data we actually analyzed
        outseg = determine_output_segment(inseg, args.psd_segment_length, args.sample_rate, args.window_fraction)
        ss = append_search_summary(xmldoc, proc_row, ifos=(station,), inseg=inseg, outseg=outseg)
        for sb in event_list:
            sb.process_id = proc_row.process_id
            sb.search = proc_row.program
            sb.ifo, sb.channel = station, setname
        xmldoc.childNodes[0].appendChild(event_list)
        fname = make_filename(station, inseg)
        utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"), verbose=True)
        
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
