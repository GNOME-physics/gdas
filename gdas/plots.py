def plot_spectrum(fd_psd):
    plot = SpectrumPlot()
    ax = plot.gca()
    ax.plot(Spectrum(fd_psd, df=fd_psd.delta_f))
    #pyplot.ylim(1e-10, 1e-3)
    pyplot.xlim(0.1, 500)
    pyplot.loglog()
    pyplot.savefig("psd.png")
    pyplot.close()

def plot_ts(ts, fname="ts.png"):
    plot = TimeSeriesPlot()
    ax = plot.gca()
    ax.plot(TimeSeries(ts, sample_rate=1.0/ts.delta_t, epoch=ts.start_time))
    ax.set_xlim(ts.start_time,ts.end_time)
    pyplot.savefig(fname)
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

def plot_triggers():
    
    print tprint(t0),"Plot trigger results..."
    events = SnglBurstTable.read('excesspower.xml.gz')
    #plot = events.plot('time', 'central_freq', "duration", "bandwidth", color='snr')
    plot = events.plot('time','central_freq',color='snr',edgecolor='none')
    #plot.set_xlim(time_start,time_end)
    plot.set_ylim(band, args.sample_rate/2.)
    plot.set_ylabel('Frequency [Hz]')
    plot.set_yscale('log')
    plot.set_title('GNOME '+station+' station event triggers')
    plot.add_colorbar(cmap='copper_r',label='Tile Energy')
    pyplot.savefig("triggers.png",dpi=400)
