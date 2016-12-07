""" Retrieve user-defined arguments""" 

import numpy,logging
from argparse import ArgumentParser
from pycbc import strain,psd

def construct_args():
    
    """
    Construct argument array
    """
    
    numpy.random.seed(0)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    argp = ArgumentParser()
    argp.add_argument("--tile-fap", type=float, default=1e-7,help="Tile false alarm probability threshold in Gaussian noise. Default is 1e-7")
    argp.add_argument("--verbose", action="store_true", help="Be verbose")
    argp.add_argument("--mainplot", action="store_true", help="Plot general figures")
    argp.add_argument("--excesspower", action="store_true", help="Do excesspower")
    argp.add_argument("--triggerplot", action="store_true", help="Plot trigger figure")
    argp.add_argument("--min-frequency", type=float, default=0, help="Lowest frequency of the filter bank, default is 0 Hz.")
    argp.add_argument("--max-frequency", type=float, default=None, help="Highest frequency of the filter bank, default is None, meaning use Nyquist.")
    argp.add_argument("--max-duration", type=float, default=None, help="Longest duration tile to compute.")
    argp.add_argument("--tile-bandwidth", type=float, default=None, help="Bandwidth of the finest filters. Default is None, and would be inferred from the data bandwidth and number of channels.")
    argp.add_argument("--channels", type=int, default=None, help="Number of frequency channels to use. Default is None, and would be inferred from the data bandwidth and tile bandwidth.")
    argp.add_argument("--analysis-start-time", type=str, default=None, help="Start analysis from this GPS time instead of the --gps-start-time")
    argp.add_argument("--analysis-end-time", type=str, default=None, help="End analysis at this GPS time instead of the --gps-end-time")
    argp.add_argument("--station", type=str, default=None, help="Name of the station")
    argp.add_argument("--window-fraction", type=float, default=None, help="Withening window fraction")
    strain.insert_strain_option_group(argp)
    psd.insert_psd_option_group(argp)
    args = argp.parse_args()

    return args
