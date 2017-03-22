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
