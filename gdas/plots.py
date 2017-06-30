import matplotlib,numpy
matplotlib.use('Agg')
from astropy.units        import Quantity
from matplotlib           import pyplot
from gwpy.plotter         import SegmentPlot,TimeSeriesPlot
from gwpy.plotter         import SpectrumPlot,SpectrogramPlot
from gwpy.segments        import SegmentList
from gwpy.spectrum        import Spectrum
from gwpy.spectrogram     import Spectrogram
from gwpy.table.lsctables import SnglBurstTable
from gwpy.timeseries      import TimeSeries
from pylab                import *
from scipy                import signal

def plot_activity(full_seglist):
    """
    Plot full activity period for station.

    Parameters
    ----------
    full_seglist : dictionary
      Continuous list of available data in the selected time period
    """
    # Import gwpy tools
    plot = SegmentPlot()
    # Initialize plotting figure
    ax = plot.gca()
    # Plot all segment in figure
    ax.plot(full_seglist)
    # Save figure
    pyplot.savefig("activity.png",dpi=300)

def plot_time_series(station,ts_list,start_time,end_time,
                     seglist=None,hp=None):
    """
    Generate a plot of the whole data time series
    """
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
        plot.add_state_segments(SegmentList(seglist[station].active),
                                plotargs={'label':'data present',
                                          'facecolor': 'g',
                                          'edgecolor': 'k'})
    # Set limits
    plot.axes[0].set_epoch(start_time)
    plot.axes[1].set_epoch(start_time)
    ax.set_xlim(start_time,end_time)
    # Fix exceeded cell block limit error
    matplotlib.pyplot.rcParams['agg.path.chunksize'] = 20000
    # Save figure
    plot.savefig('time_series.png',dpi=300)

def plot_asd(station,ts_list):
    """
    Plot Amplitude Spectral Density. AGG complexity starts to complain
    with large numbers of points. And we somehow invoke precision issues
    that need to be ameliorated.
    """
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
    ax.set_xlabel('Frequency [Hz]')
    # Set y axis to log scale
    ax.set_yscale('log')
    ax.set_ylabel('Amplitude [pT]')
    # Set x axis limits
    ax.set_xlim(1e-1, 500)
    import matplotlib.ticker as ticker
    x = ax.get_xticklabels()
    def myticks(x,pos):
        if x == 0: return "$0$"
        exponent = int(np.log10(x))
        coeff = x/10**exponent
        if coeff==1:
            return r"$10^{{ {:2d} }}$".format(exponent)
        else:
            return r"${:2.0f} \times 10^{{ {:2d} }}$".format(coeff,
                                                             exponent)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(myticks))
    # Save figure
    plot.savefig("asd.png",dpi=300)
    
def plot_whitening(station,ts_list,start_time,end_time,seglist=None):
    """
    Generate a spectrogram plot and normalized spectrogram
    norm: \sqrt{S(f,t)} / \sqrt{\overbar{S(f)}}
    """
    stride,fftlength,overlap = 20,6,3
    plot = SpectrogramPlot()
    ax = plot.gca()
    white_plot = SpectrogramPlot()
    wax = white_plot.gca()
    for ts in ts_list:
        if (len(ts) * ts.dt).value < stride:
            continue
        spec = ts.spectrogram(stride, fftlength=fftlength, overlap=overlap)
        ax.plot(spec, cmap='jet',norm=matplotlib.colors.LogNorm())
        wspec = spec.ratio('median')
        wax.plot(wspec, vmin=0.1, vmax=100,cmap='jet',
                 norm=matplotlib.colors.LogNorm())
    ax.set_title('$\mathrm{'+station+'}$')
    ax.set_ylim(0.1, ts.sample_rate.value/2.)
    ax.set_yscale('log')
    wax.set_title('$\mathrm{'+station+'}$')
    wax.set_ylim(0.1, ts.sample_rate.value/2.)
    wax.set_yscale('log')
    plot.add_colorbar(label='Amplitude')
    white_plot.add_colorbar(label='Amplitude')
    if seglist!=None:
        plot.add_state_segments(SegmentList(seglist[station].active),
                                plotargs={'label':'data present',
                                          'facecolor':'g',
                                          'edgecolor':'k'})
        white_plot.add_state_segments(SegmentList(seglist[station].active),
                                      plotargs={'label':'data present',
                                                'facecolor':'g',
                                                'edgecolor':'k'})
    # Set limits
    plot.axes[0].set_epoch(start_time)
    plot.axes[2].set_epoch(start_time)
    #plot.axes[1].set_epoch(start_time)
    white_plot.axes[0].set_epoch(start_time)
    white_plot.axes[2].set_epoch(start_time)
    ax.set_xlim(start_time,end_time)
    wax.set_xlim(start_time,end_time)
    # Save figures
    plot.savefig("spectrogram.png",dpi=300)
    white_plot.savefig("whitened.png",dpi=300)

def plot_bank(fdb):
    pyplot.figure()
    for i, fdt in enumerate(fdb[:5]):
        pyplot.plot(fdt.frequencies, fdt, 'k-')
    pyplot.grid()
    #xmin = fdb[0].frequencies[0].value
    #xmax = fdb[-1].frequencies[-1].value
    #pyplot.xlim([xmin,xmax])
    pyplot.xlabel("frequency [Hz]")
    pyplot.savefig('bank.png',dpi=300)
    pyplot.close()

def plot_filters(tdb,fmin,band):
    pyplot.figure()
    pyplot.subplots_adjust(left=0.2,right=0.95,
                           bottom=0.15,top=0.95,
                           hspace=0,wspace=0)
    for i, tdt in enumerate(tdb[:8:3]):
        ax = pyplot.subplot(3, 1, i+1)
        ax.plot(tdt.times.value - 2., numpy.real_if_close(tdt.value), 'k-')
        c_f = fmin + band/2 + 3 * (band*i) + 2.
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("%d Hz" % c_f)
        ax.set_xlim(25.0, 31.0)
        ax.set_ylim([-max(tdt.value), max(tdt.value)])
        if i!=2: pyplot.setp(ax.get_xticklabels(), visible=False)
    pyplot.savefig('filters.png',dpi=300)
    pyplot.close()
    
def plot_ts(ts, fname="ts.png"):
    plot = TimeSeriesPlot()
    ax = plot.gca()
    ax.plot(TimeSeries(ts, sample_rate=1.0/ts.delta_t, epoch=ts.start_time))
    ax.set_xlim(ts.start_time,ts.end_time)
    pyplot.savefig(fname,dpi=300)
    pyplot.close()

def plot_spectrum(fd_psd):
    plot = SpectrumPlot()
    ax = plot.gca()
    ax.plot(Spectrum(fd_psd, df=fd_psd.delta_f))
    #pyplot.ylim(1e-10, 1e-3)
    pyplot.xlim(0.1, 500)
    pyplot.loglog()
    pyplot.savefig("psd.png",dpi=300)
    pyplot.close()

def plot_spectrogram(spec,dt,df,sample_rate,start_time,end_time,
                     fname="specgram.png"):
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec,dt=dt,df=df,epoch=start_time), cmap='viridis')
    plot.add_colorbar(label='Amplitude')
    pyplot.xlim(start_time,end_time)
    pyplot.ylim(0,sample_rate/2.)
    pyplot.savefig(fname)#,dpi=300)
    pyplot.close()

def plot_spectrogram_from_ts(ts):
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec))
    #pyplot.ylim(1e-9, 1e-2)
    #pyplot.xlim(0.1, 500)
    #pyplot.loglog()
    pyplot.savefig("specgram.png",dpi=300)
    pyplot.close()

def plot_triggers(filename='excesspower.xml.gz'):
    events = SnglBurstTable.read(filename)
    #plot = events.plot('time', 'central_freq', "duration", "bandwidth", color='snr')
    plot = events.plot('time','central_freq',color='snr',edgecolor='none')
    #plot.set_xlim(time_start,time_end)
    #plot.set_ylim(band, sample_rate/2.)
    plot.set_ylabel('Frequency [Hz]')
    plot.set_yscale('log')
    #plot.set_title('GNOME '+station+' station event triggers')
    plot.add_colorbar(cmap='copper_r',label='Tile Energy')
    pyplot.savefig("triggers.png",dpi=300)

def plot_tiles():    
    bins = numpy.linspace(0, 40, 100)
    cnt = numpy.zeros(bins.shape[0]-1)
    for i, tdf in enumerate(tdb[:nchans]):
        us_rate = int(1.0 / (2 * band*nc_sum * ts_data.dt.value))
        pyplot.figure(0, figsize=(10, 10))
        pyplot.subplot(nchans, 1, i+1)
        white = tmp_ts_data.whiten(64,32,
                                   asd=numpy.sqrt(cdata_psd_tmp),
                                   window='boxcar') * sample_rate/4
        snr_1dof = numpy.convolve(tdf, white, "valid")
        # Undersample the data
        snr_1dof = snr_1dof[::us_rate]**2
        # Sum semi-adjacent samples to get 2 DOF tiles
        snr_2dof = numpy.convolve(snr_1dof, numpy.array([1, 0, 1, 0]))
        t = TimeSeries(snr_2dof,epoch=white.epoch,
                       sample_rate=int(1.0/(us_rate * tmp_ts_data.dt.value)))
        pyplot.plot(t.times + len(tdf)/2 * tdf.dt, snr_2dof, 'k-')
        pyplot.axvline(random_time)
        tmp, _ = numpy.histogram(snr_2dof, bins=bins)
        cnt += tmp
    plot_spectrogram(dof_tiles.T,
                     fname='%s/tf_%ichans_%02idof.png'%(segfolder,
                                                        nc_sum+1,
                                                        2*j))
    plot.savefig("%s/bands.png"%(segfolder),dpi=300)

def wavelet(ts_data):
    import mlpy
    z = numpy.array([float(i) for i in ts_data])
    t = numpy.array([float(i) for i in ts_data.sample_times])
    # Decimate magnetic field data to 1 sample/second
    rate = [5,10,10] if ts_data.sample_rate==500 else [8,8,8]
    for i in rate:
        z = signal.decimate(z,i,zero_phase=True)
    # Extract time every 500 sample
    t = [t[n*ts_data.sample_rate] for n in range(len(t)/ts_data.sample_rate)]
    # Convert every timing points to scale (hr,min,sec) units
    s = 60.
    t = [(t[i]-t[0])/s for i in range(len(t))]
    # Do wavelet analysis
    omega0 = 6
    fct    = "morlet"
    scales = mlpy.wavelet.autoscales(N=len(z),dt=1,dj=0.05,wf=fct,p=omega0)
    spec   = mlpy.wavelet.cwt(z,dt=1,scales=scales,wf=fct,p=omega0)
    freq   = (omega0 + numpy.sqrt(2.0 + omega0 ** 2)) / \
             (4 * numpy.pi * scales[1:]) * 1000
    idxs   = numpy.where(numpy.logical_or(freq<0.1,1000<freq))[0]
    spec   = numpy.delete(spec,idxs,0)
    freq   = numpy.delete(freq,idxs,0)
    # Initialise axis
    fig = figure(figsize=(12,8))
    plt.subplots_adjust(left=0.1,right=1,bottom=0.1,
                        top=0.94,hspace=0,wspace=0)
    ax1 = fig.add_axes([0.10,0.75,0.70,0.20])
    ax2 = fig.add_axes([0.10,0.10,0.70,0.60], sharex=ax1)
    ax3 = fig.add_axes([0.83,0.10,0.03,0.60])
    # Plot time series
    ax1.plot(t,abs(z)-numpy.average(abs(z)),'k')
    ax1.set_ylabel('Magnetic Fields [uT]')
    # Set up axis range for spectrogram
    twin_ax = ax2.twinx()
    twin_ax.set_yscale('log')
    twin_ax.set_xlim(t[0], t[-1])
    twin_ax.set_ylim(freq[-1], freq[0])
    twin_ax.tick_params(which='both', labelleft=True,
                        left=True, labelright=False)
    # Plot spectrogram
    img = ax2.imshow(numpy.abs(spec)**2,extent=[t[0],t[-1],freq[-1],freq[0]],
                     aspect='auto',interpolation='nearest',
                     cmap=cm.jet,norm=mpl.colors.LogNorm()) # cm.cubehelix
    ax2.tick_params(which='both', labelleft=False, left=False)
    ax2.set_xlabel('Time [mins]')
    ax2.set_ylabel('Frequency [mHz]',labelpad=50)
    fig.colorbar(img, cax=ax3)
    plt.savefig('wavelet.png',dpi=300)
