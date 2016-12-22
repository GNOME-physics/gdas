import matplotlib,numpy
from astropy.units   import Quantity
from matplotlib      import pyplot
from gwpy.plotter    import SegmentPlot,TimeSeriesPlot,SpectrumPlot,SpectrogramPlot
from gwpy.segments   import SegmentList
from gwpy.timeseries import TimeSeries

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
    pyplot.savefig("activity.png",dpi=500)

def plot_time_series(station,ts_list,seglist=None,hp=None):
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
    
def plot_whitening(station,ts_list,sample_rate,seglist=None):
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
        ax.plot(spec)
        wspec = spec.ratio('median')
        wax.plot(wspec, vmin=0.1, vmax=100)
    ax.set_title('$\mathrm{'+station+'}$')
    ax.set_ylim(0.1, sample_rate/2.)
    ax.set_yscale('log')
    wax.set_title('$\mathrm{'+station+'}$')
    wax.set_ylim(0.1, sample_rate/2.)
    wax.set_yscale('log')
    plot.add_colorbar(label='Amplitude')
    white_plot.add_colorbar(label='Amplitude')
    if seglist!=None:
        plot.add_state_segments(SegmentList(seglist[station].active),plotargs={'label':'data present','facecolor':'g','edgecolor':'k'})
        white_plot.add_state_segments(SegmentList(seglist[station].active),plotargs={'label':'data present','facecolor':'g','edgecolor':'k'})
    plot.savefig("spectrogram.png",dpi=500)
    white_plot.savefig("whitened.png",dpi=500)
    
def plot_ts(ts, fname="ts.png"):
    plot = TimeSeriesPlot()
    ax = plot.gca()
    ax.plot(TimeSeries(ts, sample_rate=1.0/ts.delta_t, epoch=ts.start_time))
    ax.set_xlim(ts.start_time,ts.end_time)
    pyplot.savefig(fname)
    pyplot.close()

def plot_spectrum(fd_psd):
    plot = SpectrumPlot()
    ax = plot.gca()
    ax.plot(Spectrum(fd_psd, df=fd_psd.delta_f))
    #pyplot.ylim(1e-10, 1e-3)
    pyplot.xlim(0.1, 500)
    pyplot.loglog()
    pyplot.savefig("psd.png")
    pyplot.close()

def plot_spectrogram(spec,dt,df,sample_rate,start_time,end_time,fname="specgram.png"):
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec,dt=dt,df=df,epoch=start_time), cmap='viridis')
    plot.add_colorbar(label='Amplitude')
    pyplot.xlim(start_time,end_time)
    pyplot.ylim(0,sample_rate/2.)
    pyplot.savefig(fname)
    pyplot.close()

def plot_spectrogram_from_ts(ts):
    plot = SpectrogramPlot()
    ax = plot.gca()
    ax.plot(Spectrogram(spec))
    #pyplot.ylim(1e-9, 1e-2)
    #pyplot.xlim(0.1, 500)
    #pyplot.loglog()
    pyplot.savefig("specgram.png")
    pyplot.close()

def plot_triggers(sample_rate):
    events = SnglBurstTable.read('excesspower.xml.gz')
    #plot = events.plot('time', 'central_freq', "duration", "bandwidth", color='snr')
    plot = events.plot('time','central_freq',color='snr',edgecolor='none')
    #plot.set_xlim(time_start,time_end)
    plot.set_ylim(band, sample_rate/2.)
    plot.set_ylabel('Frequency [Hz]')
    plot.set_yscale('log')
    plot.set_title('GNOME '+station+' station event triggers')
    plot.add_colorbar(cmap='copper_r',label='Tile Energy')
    pyplot.savefig("triggers.png",dpi=400)

    
