""" Retrieving magnetic field data. """

import os,glob,h5py,astropy,numpy,astropy,scipy
from astropy.time    import Time
from datetime        import datetime,timedelta
from glue.segments   import segment,segmentlist
from gwpy.segments   import DataQualityDict,DataQualityFlag
from gwpy.timeseries import TimeSeries,TimeSeriesList
from pycbc           import types

def fake_data(sample_rate,psd_segment_length):
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
    epoch   = 1153742447.0 - 10
    ts_data = numpy.random.normal(0,1,sample_rate*psd_segment_length*16)
    ts_data = types.TimeSeries(ts_data,delta_t=1.0/sample_rate,epoch=epoch)
    return ts_data

def get_data(station,starttime,endtime,activity=False,
            rep='/GNOMEDrive/gnome/serverdata/',resample=None):
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
    setname   = "MagneticFields"
    dstr      = ['%Y','%m','%d','%H','%M']
    dsplit    = '-'.join(dstr[:starttime.count('-')+1])
    start     = datetime.strptime(starttime,dsplit)
    starttime = construct_utc_from_metadata(start.strftime("%Y/%m/%d"),
                                            start.strftime("%H:%M:%S.%d"))
    dsplit    = '-'.join(dstr[:endtime.count('-')+1])
    end       = datetime.strptime(endtime,dsplit)
    endtime   = construct_utc_from_metadata(end.strftime("%Y/%m/%d"),
                                            end.strftime("%H:%M:%S.%d"))
    dataset   = []
    for date in numpy.arange(start,end,timedelta(minutes=1)):
        date = date.astype(datetime)
        path1 = rep+station+'/'+date.strftime("%Y/%m/%d/")
        path2 = station+'_'+date.strftime("%Y%m%d_%H%M*.hdf5")
        fullpath = os.path.join(path1,path2)
        dataset += glob.glob(fullpath)
    if len(dataset)==0:
        print "ERROR: No data files were found..."
        quit()
    file_order,data_order = {},{}
    for fname in dataset:
        hfile = h5py.File(fname, "r")
        segfile = file_to_segment(hfile,setname)
        file_order[segfile] = fname
        data_order[segfile] = hfile
    # Extract sample rate from metadata of last read data file
    sample_rate = hfile[setname].attrs["SamplingRate(Hz)"]
    # Estimate full segment activity list
    activity = create_activity_list(station,data_order)
    # Generate an ASCII representation of the GPS timestamped
    # segments of time covered by the input data
    seglist = segmentlist(data_order.keys())
    # Sort the segment list
    seglist.sort()
    # Create list of time series from every segment
    ts_list = generate_timeseries(file_order,setname)
    # Retrieve channel data for all the segments
    full_data = numpy.hstack([retrieve_channel_data(data_order[seg],setname)
                              for seg in seglist])
    new_sample_rate = sample_rate if resample==None else resample
    new_data_length = len(full_data)/float(sample_rate)*new_sample_rate
    full_data = scipy.signal.resample(full_data,int(new_data_length))
    # Models a time series consisting of uniformly sampled scalar values
    ts_data = types.TimeSeries(full_data,delta_t=1./new_sample_rate,
                               epoch=seglist[0][0])
    for v in data_order.values():
        v.close()
    return ts_data,ts_list,activity,int(starttime),int(endtime)

def file_to_segment(hfile,segname):
    """
    .. _file_to_segment:

    Define length of data segment. The starting and ending UTC times
    for a specific HDF5 file are determined by using the ``Date``,
    ``t0`` and ``t1`` attributes from the metadata. The
    :ref:`construct_utc_from_metadata <construct_utc_from_metadata>`
    function is then used to calculate the UTC time. Finally, the
    `segment <http://software.ligo.org/docs/glue/glue.__segments.segment-class.html>`_
    module from the ``glue.segments`` library is used to represent
    the range of times in a semi-open interval.
    
    Parameters
    ----------
    hfile : HDF5 file object
      HDF5 data file preloaded with the h5py package
    segname : str
      Attribute name of the metadata to extract.
    """
    # Extract all atributes from the data
    attrs = hfile[segname].attrs
    # Define each attribute
    dstr, t0, t1 = attrs["Date"], attrs["t0"], attrs["t1"]
    # Construct GPS starting time from data
    start_utc = construct_utc_from_metadata(dstr, t0)
    # Construct GPS ending time from data
    end_utc = construct_utc_from_metadata(dstr, t1)
    # Represent the range of times in the semi-open interval
    return segment(start_utc,end_utc)

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

def generate_timeseries(data_list, setname="MagneticFields"):
    """
    Generate time series using list of HDF5 data file paths

    Parameters
    ----------
    data_list : dictionary
      Dictionary that stored the path to the HDF5 data file
      for each segment of data available.
    
    Returns
    -------
    full_data : Array of segment's name
    """
    full_data = TimeSeriesList()
    for seg in sorted(data_list):
        hfile = h5py.File(data_list[seg], "r")
        full_data.append(retrieve_data_timeseries(hfile, "MagneticFields"))
        hfile.close()
    return full_data

def create_activity_list(station,data_order):
    """
    Create consecutive list of available data segment.

    Parameters
    ----------
    station : string
      Name of the station
    data_order : dictionary
      List of all the HDF5 data file for each segment

    Return
    ------
    full_seglist : dictionary
      Ordered list of segment
    """
    # Generate an ASCII representation of the GPS timestamped segments of time covered by the input data
    seglist = segmentlist(data_order.keys())
    # Sort the segment list
    seglist.sort()
    # Initialise dictionary for segment information
    full_seglist = DataQualityDict()
    # Save time span for each segment in ASCII file
    with open("segments.txt", "w") as fout:
        for seg in seglist:
            print >>fout, "%10.9f %10.9f" % seg
    # FIXME: Active should be masked from the sanity channel
    full_seglist[station] = DataQualityFlag(station,active=seglist.coalesce(),known=seglist.coalesce())
    return full_seglist

def retrieve_data_timeseries(hfile, setname):
    """
    Retrieve data time series from HDF5 data file

    Parameters
    ----------
    hfile : h5py file object
      Metadata from the HDF5 data file
    setname : string
      Attribute of the channel to retrieve data from
    """
    dset = hfile[setname]
    sample_rate = dset.attrs["SamplingRate(Hz)"]
    gps_epoch = construct_utc_from_metadata(dset.attrs["Date"], dset.attrs["t0"])
    data = retrieve_channel_data(hfile, setname)
    ts_data = TimeSeries(data, sample_rate=sample_rate, epoch=gps_epoch)
    return ts_data

def retrieve_channel_data(hfile, setname):
    """
    Retrieve the data from specific channel

    Parameters
    ----------
    hfile : h5py file object
      Metadata from the HDF5 data file
    setname : string
      Attribute of the channel to retrieve data from

    Return
    ------
    data : array
      Data from setname channel
    """
    return hfile[setname][:]
