import numpy
from matplotlib      import pyplot
from pylab           import *
from datetime        import datetime
from lalsimulation   import SimBurstGaussian,SimBurstSineGaussian
from gwpy.timeseries import TimeSeries
from pycbc           import types

def create_sound(ts):
    """
    Create sound based on the data

    Parameters
    ----------
    ts : TimeSeries
      Time-series magnetic field data
    """
    wout = wave.open("pure_tone.wav", "w")
    wout.setnchannels(1) # mono
    wout.setsampwidth(4) # 32 bit audio
    wout.setframerate(1000)
    wout.writeframes(ts[:])
    wout.close()

def burst_inject(ts_data,loc=0.5,duration=0.1,hrss=0.0275,amp=100.,plot=True,sine=False):
    '''
    Inject burst signal in time series data

    Parameters
    ----------
    ts_data : TimeSeries
      Time series magnetic field data
    duration : float
      Duration of the burst
    hrss : float
      hrss
    
    '''
    # Define sampling rate
    sample_rate = float(ts_data.sample_rate)
    # Define time period from sampling rate
    delta_t = 1.0 / sample_rate
    # Define start time of the time series
    t0 = float(ts_data.start_time)
    # Define final time of the time series
    t1 = t0 + len(ts_data) / sample_rate
    if sine:
        # First method to create sine gaussian burst
        f_0 = 18
        filter_band = 4
        q = math.sqrt(2)*f_0/filter_band * 2
        # Create sine gaussian burst
        hp, hx = SimBurstSineGaussian(q * 2, f_0, hrss, 1, 0, delta_t)
    else:
        # Create  gaussian burst
        hp, hx = SimBurstGaussian(duration, hrss, delta_t)
    hp = TimeSeries.from_lal(hp)
    hx = TimeSeries.from_lal(hx)
    # We rescale the amplitude to hide or expose it in the data a bit better
    hp *= amp
    if plot:
        # Plot fake burst signal
        pyplot.plot(hp.times, hp, 'k-')
        #pyplot.xlim([-0.5, 0.5])
        #pyplot.ylim([-0.1, 0.1]);
        pyplot.xlabel('Time (s)')
        pyplot.ylabel('Magnitude')
        pyplot.savefig('fakesignal.png')
        pyplot.close()
    # Define burst epoch
    hp.epoch = int(t0+loc*(t1-t0))
    # Define burst first timestamp
    st = int((hp.epoch.value-t0)*sample_rate-len(hp)/2)
    # Define burst final timestamp
    en = st+len(hp)
    # Convert time series into gwpy.timeseries.TimeSeries format
    ts_data = TimeSeries(ts_data,sample_rate=sample_rate,epoch=t0)
    # Include burst in data
    ts_data[st:en] += hp
    # Convert back to pycbc.types.TimeSeries format
    ts_data = types.TimeSeries(ts_data.value,delta_t=1.0/sample_rate,epoch=t0)
    # Return time series with burst included
    return ts_data,hp
