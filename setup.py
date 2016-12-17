from distutils.core import setup
import os,gdas.__init__
#del os.link
setup(
    name="gdas",
    version=gdas.__version__,
    author="Vincent Dumont",
    author_email="vincentdumont@gmail.com",
    packages=["gdas"],
    requires=["numpy","matplolib","scipy","astropy","glue","gwpy","pycbc"],
    include_package_data=True,
    url="https://github.com/GNOME-physics/gdas",
    description="GNOME Data Analysis package",
    install_requires=[]
)
