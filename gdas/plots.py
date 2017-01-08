import matplotlib,numpy
from astropy.units import Quantity
from matplotlib    import pyplot

def plot_activity(full_seglist):
    """
    Plot full activity period for station.

    Parameters
    ----------
    full_seglist : dictionary
      Continuous list of available data in the selected time period
    """
    from gwpy.plotter import SegmentPlot
    # Import gwpy tools
    plot = SegmentPlot()
    # Initialize plotting figure
    ax = plot.gca()
    # Plot all segment in figure
    ax.plot(full_seglist)
    # Save figure
    pyplot.savefig("activity.png",dpi=500)

def plot_time_series(station,ts_list,seglist=None,hp=None):
    """
    Generate a plot of the whole data time series
    """
    from gwpy.plotter  import TimeSeriesPlot
    from gwpy.segments import SegmentList
    plot = TimeSeriesPlot()
    ax = plot.gca()
    # Loop over all the time series 
    for ts in ts_list:
        # Plot time series for each segment
        ax.plot(ts, color='black')
    # Display title
    ax.set_title('$\mathrm{'+station+'}$')
    # Plot fake signal
    if hp!=None:
        ax.plot(hp, color='red')
    # Plot activity segments
    if seglist!=None:
        plot.add_state_segments(SegmentList(seglist[station].active),plotargs={'label':'data present','facecolor': 'g','edgecolor': 'k'})
    # Fix exceeded cell block limit error
    matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000
    # Save figure
    plot.savefig('time_series.png',dpi=500)

def plot_asd(station,ts_list):
    """
    Plot Amplitude Spectral Density. AGG complexity starts to complain
    with large numbers of points. And we somehow invoke precision issues
    that need to be ameliorated.
    """
    from gwpy.plotter import SpectrumPlot
    if station!='fake':
        for d in ts_list:
            d.x0 = Quantity(int(d.x0.value * 500), d.xunit)
            d.dx = Quantity(1, d.xunit)
        ts_list.coalesce()
        for d in ts_list:
            d.x0 = Quantity(d.x0.value / 500, d.xunit)
            d.dx = Quantity(0.002, d.xunit)
    # Initialize plotting functionality
    plot = SpectrumPlot()
    # Loop over all the time series
    for d in ts_list:
        # Generate 8 seconds per FFT with 4 second (50%) overlap
        spectrum = d.asd(8, 4)
        # Create plotting axis
        ax = plot.gca()
        # Plot square root of the spectrum
        ax.plot(numpy.sqrt(spectrum))
    # Set x axis to log scale
    ax.set_xscale('log')
    # Set y axis to log scale
    ax.set_yscale('log')
    # Set x axis limits
    ax.set_xlim(1e-1, 500)
    # Save figure
    plot.savefig("asd.png",dpi=500)
    
def plot_whitening(station,ts_list,seglist=None):
    """
    Generate a spectrogram plot and normalized spectrogram
    norm: \sqrt{S(f,t)} / \sqrt{\overbar{S(f)}}
    """
    from gwpy.plotter import SpectrogramPlot
    stride,fftlength,overlap = 20,6,3
    plot = SpectrogramPlot()
    ax = plot.gca()
    white_plot = SpectrogramPlot()
    wax = white_plot.gca()
    for ts in ts_list:
        if (len(ts) * ts.dt).value < stride:
            continue
        spec = ts.spectrogram(stride, fftlength=fftlength, overlap=overlap)
        ax.plot(spec)
        wspec = spec.ratio('median')
        wax.plot(wspec, vmin=0.1, vmax=100)
    ax.set_title('$\mathrm{'+station+'}$')
    ax.set_ylim(0.1, ts.sample_rate.value/2.)
    ax.set_yscale('log')
    wax.set_title('$\mathrm{'+station+'}$')
    wax.set_ylim(0.1, ts.sample_rate.value/2.)
    wax.set_yscale('log')
    plot.add_colorbar(label='Amplitude')
    white_plot.add_colorbar(label='Amplitude')
    if seglist!=None:
        plot.add_state_segments(SegmentList(seglist[station].active),plotargs={'label':'data present','facecolor':'g','edgecolor':'k'})
        white_plot.add_state_segments(SegmentList(seglist[station].active),plotargs={'label':'data present','facecolor':'g','edgecolor':'k'})
    plot.savefig("spectrogram.png",dpi=500)
    white_plot.savefig("whitened.png",dpi=500)

def plot_bank(fdb):
    pyplot.figure()
    for i, fdt in enumerate(fdb[:5]):
        pyplot.plot(fdt.frequencies, fdt, 'k-')
    pyplot.grid()
    pyplot.xlabel("frequency [Hz]")
    pyplot.savefig('bank.png')
    pyplot.close()

def plot_filters(tdb,flow,band):
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
    
def plot_ts(ts, fname="ts.png"):
    from gwpy.plotter    import TimeSeriesPlot
    from gwpy.timeseries import TimeSeries
    plot = TimeSeriesPlot()
    ax = plot.gca()
    ax.plot(TimeSeries(ts, sample_rate=1.0/ts.delta_t, epoch=ts.start_time))
    ax.set_xlim(ts.start_time,ts.end_time)
    pyplot.savefig(fname)
    pyplot.close()

def plot_spectrum(fd_psd):
    from gwpy.plotter  import SpectrumPlot
    from gwpy.spectrum import Spectrum
    plot = SpectrumPlot()
    ax = plot.gca()
    ax.plot(Spectrum(fd_psd, df=fd_psd.delta_f))
    #pyplot.ylim(1e-10, 1e-3)
    pyplot.xlim(0.1, 500)
    pyplot.loglog()
    pyplot.savefig("psd.png")
    pyplot.close()

def plot_spectrogram(spec,dt,df,sample_rate,start_time,end_time,fname="specgram.png"):
    from gwpy.plotter     import SpectrogramPlot
    from gwpy.spectrogram import Spectrogram
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec,dt=dt,df=df,epoch=start_time), cmap='viridis')
    plot.add_colorbar(label='Amplitude')
    pyplot.xlim(start_time,end_time)
    pyplot.ylim(0,sample_rate/2.)
    pyplot.savefig(fname)
    pyplot.close()

def plot_spectrogram_from_ts(ts):
    from gwpy.plotter     import SpectrogramPlot
    from gwpy.spectrogram import Spectrogram
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec))
    #pyplot.ylim(1e-9, 1e-2)
    #pyplot.xlim(0.1, 500)
    #pyplot.loglog()
    pyplot.savefig("specgram.png")
    pyplot.close()

def plot_triggers():
    from gwpy.table.lsctables import SnglBurstTable
    events = SnglBurstTable.read('excesspower.xml.gz')
    #plot = events.plot('time', 'central_freq', "duration", "bandwidth", color='snr')
    plot = events.plot('time','central_freq',color='snr',edgecolor='none')
    #plot.set_xlim(time_start,time_end)
    #plot.set_ylim(band, sample_rate/2.)
    plot.set_ylabel('Frequency [Hz]')
    plot.set_yscale('log')
    #plot.set_title('GNOME '+station+' station event triggers')
    plot.add_colorbar(cmap='copper_r',label='Tile Energy')
    pyplot.savefig("triggers.png",dpi=400)

def plot_tiles():    
    from gwpy.timeseries import TimeSeries
    bins = numpy.linspace(0, 40, 100)
    cnt = numpy.zeros(bins.shape[0]-1)
    for i, tdf in enumerate(tdb[:nchans]):
        us_rate = int(1.0 / (2 * band*nc_sum * ts_data.dt.value))
        pyplot.figure(0, figsize=(10, 10))
        pyplot.subplot(nchans, 1, i+1)
        white = tmp_ts_data.whiten(64, 32, asd=numpy.sqrt(cdata_psd_tmp), window='boxcar') * sample_rate/4
        snr_1dof = numpy.convolve(tdf, white, "valid")
        # Undersample the data
        snr_1dof = snr_1dof[::us_rate]**2
        # Sum semi-adjacent samples to get 2 DOF tiles
        snr_2dof = numpy.convolve(snr_1dof, numpy.array([1, 0, 1, 0]))
        t = TimeSeries(snr_2dof, epoch=white.epoch, sample_rate=int(1.0/(us_rate * tmp_ts_data.dt.value)))
        pyplot.plot(t.times + len(tdf)/2 * tdf.dt, snr_2dof, 'k-')
        pyplot.axvline(random_time)
        tmp, _ = numpy.histogram(snr_2dof, bins=bins)
        cnt += tmp
    plot_spectrogram(dof_tiles.T,fname='%s/tf_%ichans_%02idof.png'%(segfolder,nc_sum+1,2*j))
    plot.savefig("%s/bands.png"%(segfolder))
