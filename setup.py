from distutils.core import setup
import os,gdas.__init__
from glob import glob

def get_data_names(root):
    """
    Return list of all filenames (not directories) under root.
    """
    temp = []
    for dirpath, dirnames, filenames in os.walk(root):
        temp.extend((os.path.join(dirpath, d, '*') for d in dirnames))
    names = []
    for path in temp:
        if any(os.path.isfile(f) for f in glob(path)):
            names.append(path[5:])
    return names

#del os.link
setup(
    name="gdas",
    version=gdas.__version__,
    author="Vincent Dumont",
    author_email="vincentdumont@gmail.com",
    packages=["gdas"],
    requires=["numpy","matplotlib","scipy","astropy","gwpy","pycbc"],
    package_data = {'gdas' : get_data_names('gdas/ligo')},
    include_package_data=True,
    url="https://gnome-physics.github.io/gdas",
    description="GNOME Data Analysis package",
    install_requires=[]
)
