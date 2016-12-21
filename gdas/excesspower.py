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
    
    
