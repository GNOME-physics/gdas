Working Example
===============

In this section we present an example on how to use the analysis tools available on gdas to create your own data analysis script. Either on a Jupyter notebook or on a Python script, the first thing to do is to import the ``gdas`` package that contain all the modules present in the GNOME software. That can be done easily by doing the following::

  import gdas

In order to retrieve a specific chunk of data to be analyzed for a particular station, the name of the station along with the start and end dates should be specified::
  
  station    = 'fribourg01'
  start_time = '2016-11-03-04'
  end_time   = '2016-11-03-04-2'

where the start and end times should always have at least the year, month and day specified, and with the values separated by a dash symbol. Hour and minute can also be specified.

If you are not working on the server and the data are located in a different repository than ``/GNOMEDrive/gnome/serverdata/``, a custom path can be defined. For instance::
  
  datapath = '/Users/vincent/data/GNOMEDrive/gnome/serverdata/'

The magnetic field data can then be retrieve as follows::
  
  ts_data,ts_list,activity = gdas.magfield(station,start_time,end_time,rep=datapath)

The ``gdas.magfield`` method will return 3 arrays of data that can then be used to produce different plots::
  
  gdas.plot_activity(activity)
  gdas.plot_time_series(station,ts_list,seglist=activity)
  gdas.plot_asd(station,ts_list)
  gdas.plot_whitening(station,ts_list,activity)

This is a script to do Excess Power analysis::
  
  psd_segment_length = 60
  psd_segment_stride = 30
  psd_estimation     = 'median-mean'
  window_fraction    = 0
  tile_fap           = 1e-5
  channels           = 250
  
  gdas.excess_power(ts_data,psd_segment_length,psd_segment_stride,psd_estimation,window_fraction,tile_fap,station,nchans=channels)
  gdas.plot_triggers()
