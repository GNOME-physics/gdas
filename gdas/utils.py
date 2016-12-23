"""Other routines"""

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
