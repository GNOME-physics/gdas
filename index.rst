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

The software is fairly easy to use, you can either use the main gdas script present in the `scripts folder <https://github.com/GNOME-physics/gdas/scripts>`_ or you can follow the steps from one of the notebooks available `here <https://github.com/GNOME-physics/gdas/notebooks>`_ to write your own script. The notebooks can also be used as a guide on how the analysis is done. Finally, a much more thorough description on how the theory behind the Excess Power search analysis is being translated into code is provided at the following link: https://gnome-physics.github.io/

.. raw:: html

 <a href="https://github.com/GNOME-physics/gdas"><img style="position: absolute; top: 0; right: 0; border: 0;" src="https://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" alt="Fork me on GitHub"></a>
  
.. toctree::
   :maxdepth: 2

   index.rst

arguments
---------

Construct argument array based on user-defined parameters.

.. currentmodule:: gdas.arguments

.. autosummary::
   :toctree: generated/

   construct_args

retrieve
--------

Extract magnetic field data from HDF5 files.

.. currentmodule:: gdas.retrieve

.. autosummary::
   :toctree: generated/

   convertdate
   magfield
   file_to_segment
   construct_utc_from_metadata
   generate_timeseries
   retrieve_data_timeseries
   retrieve_channel_data
   
