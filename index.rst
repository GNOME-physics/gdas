GNOME Data Analysis Software
============================

.. toctree::
   :maxdepth: 2

   index.rst

Introduction
============

This package contains functions useful for magnetic field signal processing, with a focus on Excess Power search analysis and application on the data for the GNOME collaboration, see `Pustelny et al. (2013) <https://arxiv.org/abs/1303.5524>`_. This documentation details all the available functions and tasks available through this software. Here are some example tasks that can (or will soon to) be handled:

* Plot usual time series and spectrogram of magnetic field data.
* Perform excess power analysis and plot detected triggers in time-frequency map.
* Create artificial data for testing data analysis.
* Inject fake signal of different bandwidth and durations.
* Cross-correlation of continuous sine wave signals.
* Perform Allan Standard deviation.

.. raw:: html

 <a href="https://github.com/GNOME-physics/gdas"><img style="position: absolute; top: 0; right: 0; border: 0;" src="https://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" alt="Fork me on GitHub"></a>

Installation
============

The program requires the following general packages to run: `Numpy <http://numpy.scipy.org/>`_, `Matplotlib <http://matplotlib.sourceforge.net/>`_, `Scipy <http://www.scipy.org/>`_ and `Astropy <http://www.astropy.org/>`_. The following LIGO-related packages are also required for full functionality: `Gwpy <https://gwpy.github.io/>`_, `PyCBC <https://github.com/ligo-cbc/pycbc>`_, `Glue <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_, `LAL <http://software.ligo.org/docs/lalsuite/lal/index.html>`_, `LALburst <http://software.ligo.org/docs/lalsuite/lalburst/index.html>`_ and `LALsimulation <http://software.ligo.org/docs/lalsuite/lalsimulation/index.html>`_.

While most of the packages can be installed automatically using `pip <http://www.pip-installer.org/en/latest/index.html>`_, some LIGO packages (Glue, LAL, LALburst and LALsimulation) must be installed separately beforehand as they contain several C routines that need specific compilation. However, these packages are already included in a bigger package called `LALsuite <https://wiki.ligo.org/DASWG/LALSuite>`_ which can be installed fairly easily on Debian (Linux) and Mac OS machines.

LALsuite tools
--------------

Some useful pages on how to download and install the LIGO software can be found `here <https://wiki.ligo.org/DASWG/HowToDocs>`_.

MacPorts (Mac)
~~~~~~~~~~~~~~

For Mac users, the installation is pretty easy, detailed information can be found on `this page <https://wiki.ligo.org/DASWG/MacPorts>`_. You need to have `MacPorts <https://www.macports.org/install.php>`_ installed. The following commands should suffice to install the LALsuite package on your machine::

  sudo port install lscsoft-deps
  sudo port install glue
  sudo port install lalapps

The first command will install all the dependencies needed for the LIGO software to be installed. The following 2 commands will install the actual packages.

apt-get (Debian)
~~~~~~~~~~~~~~~~

Since the LIGO software is not a default package in the apt package manager system on Debian machine, additional steps will be needed. The first step is to add the following links to the source list located at ``/etc/apt/sources.list``::

  deb [arch=amd64] http://software.ligo.org/lscsoft/debian jessie contrib
  deb-src [arch=amd64] http://software.ligo.org/lscsoft/debian jessie contrib

Note that the ``[arch=amd64]`` is needed to fix the architecture problem in case it tries to install i386 version on 64-bit Debian. Once the sources have been added, you must first install all the dependencies as follows::

  apt-get install build-essential automake autoconf libtool devscripts

The LIGO software can finally be installed using the following command::

  apt-get install lscsoft-all
  
Main Program
------------

The best way to install the GNOME software along with the rest of the dependencies is by using `pip`::

   pip install gdas

(You may need to put a ``sudo`` in front of this). For this to work
you need to have `pip
<http://www.pip-installer.org/en/latest/index.html>`_ installed. This
method allows for easy uninstallation.

You can also simply download the tarball from the PyPI website, unpack it and then do::

   python setup.py install

The latest stable package can be downloaded from PyPI: https://pypi.python.org/pypi/gdas.
The development version can be downloaded from `here <https://github.com/GNOME-physics/gdas>`_.

Multi-user Server
=================

A GNOME JupyterHub, or multi-user server has been created to allow each member to access the entire available dataset. Member who do not have access to the server but wish to access it should send a request to Dr. Sam Afach. Member who are not part of the GNOME collaboration will not be granted access to the dataset but are free to use our software on their own data.

The server can be accessed in two ways, either by acceding the `server's webpage <https://budker.uni-mainz.de:8000/hub/login>`_, or from your terminal through SSH::

  ssh -X username@budker.uni-mainz.de -p 8022

While SSH is very handy for people using UNIX-like operating systems, this can become more complicated for those working on Windows machines. Fortunately, access to a terminal is also possible through the webpage, which means directly from your internet browser! This can be done by clicking on the New tab after login and select Terminal:

.. figure:: img/jupyter1.png
	    :width: 70%
	    :align: center

You can then use the terminal window to access files and create new Python scripts for your analysis.

.. figure:: img/jupyter2.png
	    :width: 70%
	    :align: center

Getting Started
===============

Either on your own computer or on the server, on a Jupyter notebook or on a Python script, the first thing to do is to import the ``gdas`` package that contain all the modules present in the GNOME software. That can be done easily by doing the following::

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


Functionalities
===============
   
Extract Magnetic Field Data
---------------------------

Extract magnetic field data from HDF5 files.

.. currentmodule:: gdas.retrieve

.. autosummary::
   :toctree: generated/

   magfield
   file_to_segment
   construct_utc_from_metadata
   generate_timeseries
   create_activity_list
   retrieve_data_timeseries
   retrieve_channel_data
   
Plotting routines
-----------------

Methods to produce time-frequency plots and others

.. currentmodule:: gdas.plots

.. autosummary::
   :toctree: generated/

   plot_activity
   plot_time_series
   plot_asd
   plot_whitening
   plot_ts
   plot_spectrum
   plot_spectrogram
   plot_spectrogram_from_ts
   plot_triggers

Excess Power Search Analysis
----------------------------

Main class to do excess-power search analysis

.. currentmodule:: gdas.epower

.. autosummary::
   :toctree: generated/

   excess_power
   check_filtering_settings
   calculate_psd
   calculate_spectral_correlation
   create_filter_bank
   convert_to_time_domain
   identify_block
   create_tf_plane
   compute_filter_ips_self
   compute_filter_ips_adjacent
   compute_channel_renormalization
   measure_hrss
   measure_hrss_slowly
   measure_hrss_poorly
   trigger_list_from_map
   determine_output_segment
   make_tiles
   make_indp_tiles
   make_filename
   construct_tiles
   create_tile_duration
   create_xml
   
Utilities
---------

Independent routines to do various other things

.. currentmodule:: gdas.utils

.. autosummary::
   :toctree: generated/

   create_sound
