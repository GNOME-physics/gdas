import os,glob,h5py,astropy,numpy,scipy
from astropy.time    import Time
from datetime        import datetime,timedelta
from glue.segments   import segment,segmentlist
from gwpy.segments   import DataQualityDict,DataQualityFlag
from gwpy.timeseries import TimeSeries,TimeSeriesList
from pycbc           import types

def impulse_data(epoch=1153742417.0,sample_rate=512,psd_segment_length=60):
    """
    Create fake time series data. The flux data is generated using a
    random Gaussian distribution.

    Parameters
    ----------
    sample_rate : int
      Sampling rate of fake data
    psd_segment_length : int
      Length of each segment in seconds
    """
    ts_data = numpy.zeros(sample_rate * psd_segment_length)
    ts_data = types.TimeSeries(ts_data, delta_t=1.0/sample_rate, epoch=epoch)
    return ts_data

def fake_data(epoch=1153742417.0,sample_rate=512,psd_segment_length=60,nsegs=16):
    """
    Create fake time series data. The flux data is generated using a
    random Gaussian distribution.

    Parameters
    ----------
    sample_rate : int
      Sampling rate of fake data
    psd_segment_length : int
      Length of each segment in seconds
    nsegs : int
      Number of segments present in time series
    """
    ts_data = numpy.random.normal(0,1,sample_rate*psd_segment_length*nsegs)
    ts_data = types.TimeSeries(ts_data,delta_t=1.0/sample_rate,epoch=epoch)
    return ts_data

def get_data(station,starttime,endtime,rep='/GNOMEDrive/gnome/serverdata/',resample=None):
    """
    Glob all files withing user-defined period and extract data.
    
    Parameters
    ----------
    station : str
      Name of the station to be analysed
    t0 : int
      GPS timestamp of the first required magnetic field data
    t1 : int
      GPS timestamp of the last required magnetic field data
    
    Return
    ------
    ts_data, ts_list, activity : TimeSeries, dictionary, list
      Time series data for selected time period, list of time series
      for each segment, sampling rate of the retrieved data
    """
    # Define data attribute to be extracted from HDF5 files
    setname   = "MagneticFields"
    dstr      = ['%Y','%m','%d','%H','%M','%S','%f']
    dsplit    = '-'.join(dstr[:starttime.count('-')+1])
    start     = datetime.strptime(starttime,dsplit)
    dsplit    = '-'.join(dstr[:endtime.count('-')+1])
    end       = datetime.strptime(endtime,dsplit)
    dataset   = []
    for date in numpy.arange(start,end,timedelta(minutes=1)):
        date = date.astype(datetime)
        path1 = rep+station+'/'+date.strftime("%Y/%m/%d/")
        path2 = station+'_'+date.strftime("%Y%m%d_%H*.hdf5")
        fullpath = os.path.join(path1,path2)
        dataset += glob.glob(fullpath)
    if len(dataset)==0:
        print "ERROR: No data files were found..."
        quit()
    file_order,data_order = {},{}
    for fname in dataset:
        hfile = h5py.File(fname, "r")
        # Extract all atributes from the data
        attrs = hfile[setname].attrs
        # Define each attribute
        dstr, t0, t1 = attrs["Date"], attrs["t0"], attrs["t1"]
        # Construct GPS starting time from data
        start_utc = construct_utc_from_metadata(dstr, t0)
        # Construct GPS ending time from data
        end_utc = construct_utc_from_metadata(dstr, t1)
        # Represent the range of times in the semi-open interval
        segfile = segment(start_utc,end_utc)
        file_order[segfile] = fname
        data_order[segfile] = hfile
    # Create list of time series from every segment
    ts_list = TimeSeriesList()
    for seg in sorted(file_order):
        hfile = h5py.File(file_order[seg], "r")
        dset = hfile[setname]
        sample_rate = dset.attrs["SamplingRate(Hz)"]
        gps_epoch = construct_utc_from_metadata(dset.attrs["Date"], dset.attrs["t0"])
        data = hfile[setname][:]
        ts_data = TimeSeries(data, sample_rate=sample_rate, epoch=gps_epoch)
        ts_list.append(ts_data)
        hfile.close()
    # Generate an ASCII representation of the GPS timestamped segments of time covered by the input data
    seglist = segmentlist(data_order.keys())
    # Sort the segment list
    seglist.sort()
    # Initialise dictionary for segment information
    activity = DataQualityDict()
    # Save time span for each segment in ASCII file
    with open("segments.txt", "w") as fout:
        for seg in seglist:
            print >>fout, "%10.9f %10.9f" % seg
    # FIXME: Active should be masked from the sanity channel
    activity[station] = DataQualityFlag(station,active=seglist.coalesce(),known=seglist.coalesce())
    # Generate an ASCII representation of the GPS timestamped segments of time covered by the input data
    seglist = segmentlist(data_order.keys())
    # Sort the segment list
    seglist.sort()
    # Retrieve channel data for all the segments
    full_data = numpy.hstack([data_order[seg][setname][:] for seg in seglist])
    new_sample_rate = float(sample_rate) if resample==None else float(resample)
    new_data_length = len(full_data)*new_sample_rate/float(sample_rate)
    full_data = scipy.signal.resample(full_data,int(new_data_length))
    # Models a time series consisting of uniformly sampled scalar values
    ts_data = types.TimeSeries(full_data,delta_t=1./new_sample_rate,epoch=seglist[0][0])
    for v in data_order.values():
        v.close()
    return ts_data,ts_list,activity

def time_convert(starttime,endtime):
    dstr      = ['%Y','%m','%d','%H','%M','%S','%f']
    dsplit    = '-'.join(dstr[:starttime.count('-')+1])
    start     = datetime.strptime(starttime,dsplit)
    starttime = construct_utc_from_metadata(start.strftime("%Y/%m/%d"),
                                            start.strftime("%H:%M:%S.%f"))
    dsplit    = '-'.join(dstr[:endtime.count('-')+1])
    end       = datetime.strptime(endtime,dsplit)
    endtime   = construct_utc_from_metadata(end.strftime("%Y/%m/%d"),
                                            end.strftime("%H:%M:%S.%f"))
    return starttime,endtime

def construct_utc_from_metadata(datestr, t0str):
    """
    .. _construct_utc_from_metadata:

    Constructing UTC timestamp from metadata

    Parameters
    ----------
    datestr : str
      Date of the extracted data
    t0str : str
      GPS time
    """
    instr = "%d-%d-%02dT" % tuple(map(int, datestr.split('/')))
    instr += t0str
    t = astropy.time.Time(instr, format='isot', scale='utc')
    return t.gps
