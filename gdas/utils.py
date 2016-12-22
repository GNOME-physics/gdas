"""Filtering methods"""

def check_filtering_settings(sample_rate,min_frequency,max_frequency,channels=None,tile_bandwidth=None):
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
    if max_frequency is None or max_frequency>sample_rate/2.:
        # Set the tile maximum frequency equal to the Nyquist frequency (i.e. half the sampling rate)
        max_frequency = sample_rate / 2.0
    # Check whether or not tile bandwidth and channel are defined
    if tile_bandwidth is None and channels is None:
        # Exit program with error message
        exit("Either --tile-bandwidth or --channels must be specified to set up time-frequency plane")
    else:
        # Define as assert statement that tile maximum frequency larger than its minimum frequency
        assert max_frequency >= min_frequency
        # Define spectral band of data
        data_band = max_frequency - min_frequency
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
    flow = min_frequency
    return nchans,band,flow

def create_sounds(ts):
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
    
