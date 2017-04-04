Excess Power Overview
=====================

The **Excess Power method** is known as the *optimal detection strategy* to search for burst signals for which only the duration and frequency band are known, which is basically the case for GNOME and its search of Axion-Like Particles (ALP). This method was developed and introduced by `Anderson et al. (200) <https://arxiv.org/pdf/gr-qc/0008066v1.pdf>`_ and has been extensively used in the detection of burst sources of gravitational radiation. A more technical documentation was written by `Brady et al. (2007) <http://www.lsc-group.phys.uwm.edu/~siemens/power.pdf>`_ describing how the algorithm used by the LIGO collaboration works and how the theory is translated into code.

We present below a step-by-step procedure followed during the Excess Power search analysis. For a better representation of what is happening, the figure at the end shows how the data is being split and analysed to search for multiple signals of different bandwidth and duration in the time-frequency plane.

- :ref:`Time domain segmentation and PSD estimate <psdestimate>`

   We first estimate the instrument's noise Power Spectral Density (PSD) by splitting the time-series data into multiple overlapping segments. A periodogram for each segment is calculated separately and then averaged, which will reduce the variance of the individual power measurements. The result is a frequency series where samples are separated in frequency space by :math:`\Delta f` equal to the inverse of a segment’s length and with a high end frequency limit equal to the Nyquist limit. The final power spectrum will help reveal the existence, or the absence, of repetitive patterns and correlation structures in a signal process.

- :ref:`Comb of frequency channels <filterbank>`

   We then split the PSD frequency series into multiple channels. For each channel, a frequency domain filter is created with a :math:`\Delta f` determined by the PSD and a total extent in Fourier space that is twice the stated bandwidth of a channel. The result is a list of each channel filter's frequency series.

- :ref:`Creating analysing blocks <analysingblocks>`

   The Excess Power method can lead to moderately-large computational requirements, and it has been found that the computational efficiency of this implementation can be improved upon by considering blocks of data that are much longer than the longest signal time duration. The entire time series is therefore split into separate blocks. We use the length of the segments used for PSD estimate to define the duration of each block. For each block, the time series is c0Aonverted into frequency series which is then filtered by the filter bank throughout all the channels. A time-frequency map is finally created which stores all the filtered frequency series from each channel.

- :ref:`Creating tiles with different bandwidth <tilebandwidth>`

   We can now construct tiles with different bandwidth by summing multiple channels together.

- :ref:`Exploring tiles with different duration <tileduration>`

   For each given tile's bandwidth, one can investigate different tile's duration. This can be done by exploring different number of degrees of freedom, :math:`d`, which can be calculated as follows: :math:`d=2BT` where :math:`B` and :math:`T` are respectively the bandwidth and duration of the tile. Section 2.2.5 of `Brady et al. <http://www.lsc-group.phys.uwm.edu/~siemens/power.pdf>`_ gives a great description of how to interpret the number of degrees of freedom. Therefore, by changing the :math:`d`, one can explore multiple tile's duration for different bandwidth.

- :ref:`Define triggering signal <triggerfinding>`

   The energy of each tile in the time-frequency space is calculated and compare to a user-defined threshold value. After defining a tile false alarm probability threshold in Gaussian noise and using the number of degrees of freedom for each tile, one can define a energy threshold value above which a burst trigger can be identified by comparing the energy threshold with the tile's energy in the time-frequency map. A tile energy time frequency map plot similar to Figure 5 in `Pustelny et al. (2013) <http://budker.berkeley.edu/~vincent/gnome/documentation/2013_pustelny.pdf>`_ can then be made which plots the outlying tile energies present in the data.

.. figure:: ./img/overview.png

	    Overview of the Excess Power method and difference between segments, channels, tiles and blocks.
