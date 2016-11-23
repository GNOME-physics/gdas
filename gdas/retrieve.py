""" Retrieving magnetic field data.""" 
from __future__ import division, print_function, unicode_literals

import os,glob,h5py
from astropy.time import Time

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
    end = datetime(1980,1,6)+timedelta(seconds=endtime)
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
    sample_rate = hfile[setname].attrs["SamplingRate(Hz)"]

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
    t = Time(instr, format='isot', scale='utc')
    return t.gps
