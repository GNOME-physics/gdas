""" Retrieving magnetic field data.""" 

import os,glob,h5py,astropy,numpy,astropy
from astropy.time import Time
from datetime import datetime,timedelta
from glue.segments import segment,segmentlist
from gwpy.timeseries import TimeSeries,TimeSeriesList
from pycbc import types

def convertdate(start,end):
    """
    Convert start and end dates to GPS time

    Parameters
    ----------
    start : str
      Starting date in format YYYY-MM-DD[-HH-MM]
    end : str
      Ending date in format YYYY-MM-DD[-HH-MM]
    """
    dstr   = ['%Y','%m','%d','%H','%M']
    date1  = '-'.join(dstr[:start.count('-')+1])
    date1  = datetime.strptime(start,date1)
    start  = int((date1-datetime(1980,1,6)).total_seconds())
    date2  = '-'.join(dstr[:end.count('-')+1])
    date2  = datetime.strptime(end,date2)
    end    = int((date2-datetime(1980,1,6)).total_seconds())
    return start,end

def magfield(station,starttime,endtime):
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
    """    
    setname = "MagneticFields"
    start = datetime(1980,1,6)+timedelta(seconds=starttime)
    end   = datetime(1980,1,6)+timedelta(seconds=endtime)
    dataset = []
    for date in numpy.arange(start,end,timedelta(minutes=1)):
        date = date.astype(datetime)
        year,month,day,hour,minute = date.year,date.month,date.day,date.hour,date.minute
        path1 = '/Users/vincent/ASTRO/data/GNOMEDrive/gnome/serverdata/'
        path2 = "%s/%s/%02i/%02i/"%(station,year,month,day)
        path3 = "%s_%s%02i%02i_%02i%02i*.hdf5"%(station,year,month,day,hour,minute)
        fullpath = os.path.join(path1,path2,path3)
        dataset += glob.glob(fullpath)
    file_order,data_order = {},{}
    for fname in dataset:
        hfile = h5py.File(fname, "r")
        segfile = file_to_segment(hfile,setname)
        file_order[segfile] = fname
        data_order[segfile] = hfile
    # Extract sample rate from metadata of last read data file
    sample_rate = hfile[setname].attrs["SamplingRate(Hz)"]
    # Generate an ASCII representation of the GPS timestamped segments
    # of time covered by the input data
    seglist = segmentlist(data_order.keys())
    # Sort the segment list
    seglist.sort()
    # Generate time series for the ensemble of data
    data_list = generate_timeseries(file_order,setname)
    # Retrieve channel data for all the segments
    full_data = numpy.hstack([retrieve_channel_data(data_order[seg],setname) for seg in seglist])
    # Models a time series consisting of uniformly sampled scalar values
    ts_data = types.TimeSeries(full_data,delta_t=1/sample_rate,epoch=seglist[0][0])
    for v in data_order.values():
        v.close()        
    return ts_data

def file_to_segment(hfile,segname):
    """
    Define length of data segment
    
    Parameters
    ----------
    hfile : HDF5 file object
      HDF5 data file preloaded with the h5py package
    segname : str
      Attribute name of the metadata to extract.
    """
    attrs = hfile[segname].attrs
    dstr, t0, t1 = attrs["Date"], attrs["t0"], attrs["t1"]
    start_utc = construct_utc_from_metadata(dstr, t0)
    end_utc = construct_utc_from_metadata(dstr, t1)
    return segment(start_utc,end_utc)

def construct_utc_from_metadata(datestr, t0str):
    """
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
    full_data : 
    """
    full_data = TimeSeriesList()
    for seg in sorted(data_list):
        hfile = h5py.File(data_list[seg], "r")
        full_data.append(retrieve_data_timeseries(hfile, "MagneticFields"))
        hfile.close()
    return full_data

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
    """
    return hfile[setname][:]
