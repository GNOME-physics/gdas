GNOME Data Analysis Software
============================

This package contains functions useful for magnetic field signal processing, with a focus on Excess Power search analysis and application on the data for the GNOME collaboration, see `Pustelny et al. (2013) <https://arxiv.org/abs/1303.5524>`_. This documentation details all the available functions and tasks available through this software. Here are some example tasks that can be handled:

* Plot usual time series and spectrogram of magnetic field data.
* Perform excess power analysis and plot detected triggers in time-frequency map.
* Create artificial data for testing data analysis.
* Inject fake signal of different bandwidth and durations.
* Cross-correlation of continuous sine wave signals.
* Perform Allan Standard deviation.

It requires the following general packages to run: `Numpy <http://numpy.scipy.org/>`_, `Matplotlib <http://matplotlib.sourceforge.net/>`_, `Scipy <http://www.scipy.org/>`_ and `Astropy <http://www.astropy.org/>`_. The following LIGO-related packages are also required for full functionality: `Glue <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_, `Gwpy <https://gwpy.github.io/>`_, `PyCBC <https://github.com/ligo-cbc/pycbc>`_, `lal <http://software.ligo.org/docs/lalsuite/lal/index.html>`_, `lalburst <http://software.ligo.org/docs/lalsuite/lalburst/index.html>`_ and `LALsimulation <http://software.ligo.org/docs/lalsuite/lalsimulation/index.html>`_

The best way to install the package is by using `pip`::

   pip install gdas

(You may need to put a ``sudo`` in front of this). For this to work
you need to have `pip
<http://www.pip-installer.org/en/latest/index.html>`_ installed. This
method allows for easy uninstallation.

You can also simply download the tarball from the PyPI website, unpack
it and then do::

   python setup.py install

The latest stable package can be downloaded from PyPI: https://pypi.python.org/pypi/gdas.
The development version can be downloaded from `here <https://github.com/GNOME-physics/gdas>`_.

.. raw:: html

 <a href="https://github.com/GNOME-physics/gdas"><img style="position: absolute; top: 0; right: 0; border: 0;" src="https://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" alt="Fork me on GitHub"></a>
  
.. toctree::
   :maxdepth: 2

   index.rst

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
   create_xml
   
Utilities
---------

Independent routines to do various other things

.. currentmodule:: utils

.. autosummary::
   :toctree: generated/

   create_sound
