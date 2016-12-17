""" Retrieving magnetic field data.""" 

import os,glob,h5py,astropy,numpy,astropy
from astropy.time import Time
from datetime import datetime,timedelta
from glue.segments import segment

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
    
    # Get start date in GPS time
    date1 = start.split('-')
    year1, month1, day1 = [int(i) for i in date1[:3]]
    hour1 = int(date1[3]) if len(date1)>3 else 0 
    min1  = int(date1[4]) if len(date1)>4 else 0
    date1 = datetime(year1,month1,day1,hour1,min1)
    start = int((date1-datetime(1980,1,6)).total_seconds())
    
    # Get end date in GPS time
    date2 = end.split('-')
    year2, month2, day2 = [int(i) for i in date2[:3]]
    hour2 = int(date2[3]) if len(date2)>3 else 0
    min2  = int(date2[4]) if len(date2)>4 else 0
    date2 = datetime(year2,month2,day2,hour2,min2)
    end   = int((date2-datetime(1980,1,6)).total_seconds())

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
    data_order = {}
    for fname in dataset:
        hfile = h5py.File(fname, "r")
        segfile = file_to_segment(hfile,setname)
        data_order[segfile] = hfile
        
    return data_order

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
