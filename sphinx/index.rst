GNOME Data Analysis Software
============================

This package contains functions useful for magnetic field signal processing, with a focus on Excess Power search analysis and application on the data for the GNOME collaboration, see `Pustelny et al. (2013) <https://arxiv.org/abs/1303.5524>`_. This documentation details all the available functions and tasks available through this software. Here are some example tasks that can (or will soon to) be handled:

* Plot usual time series and spectrogram of magnetic field data.
* Perform excess power analysis and plot detected triggers in time-frequency map.
* Create artificial data for testing data analysis.
* Inject fake signal of different bandwidth and durations.
* Cross-correlation of continuous sine wave signals.
* Perform Allan Standard deviation.

Should you have any questions or suggested corrections to be made, do not hesitate to `contact me <mailto:vincentdumont11@gmail.com>`_.

.. raw:: html

   <a href="https://github.com/GNOME-physics/gdas"><img style="position: absolute; top: 0; right: 0; border: 0;" src="https://s3.amazonaws.com/github/ribbons/forkme_right_orange_ff7600.png" alt="Fork me on GitHub"></a>

Getting Started
---------------

.. toctree::
   :maxdepth: 1

   installation
   server
   example

Excess Power Search Analysis
----------------------------

.. toctree::
   :maxdepth: 1

   excess_power
   excess_power_submodules
