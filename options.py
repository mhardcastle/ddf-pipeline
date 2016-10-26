# Options for killms pipeline code

import ConfigParser
import os

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile:
        lines=len(open(nodefile).readlines())
        return lines
    else:
        import multiprocessing
        return multiprocessing.cpu_count()

def options(filename):

    # option_list format is: section, name, type, default
    # names must be unique -- section names are not used in output dict

    option_list = ( ( 'machine', 'NCPU_DDF', int, getcpus() ),
                    ( 'machine', 'NCPU_killms', int, getcpus() ),
                    ( 'data', 'mslist', str, None ),
                    ( 'data', 'full_mslist', str, None ),
                    ( 'solutions', 'ndir', int, 30 ),
                    ( 'solutions', 'NChanSols', int, 1 ),
                    ( 'solutions', 'dt', float, 1. ),
                    ( 'solutions', 'LambdaKF', float, 0.5 ),
                    ( 'image', 'imsize', int, 20000 ),
                    ( 'image', 'cellsize', float, 1.5 ),
                    ( 'image', 'robust', float, -0.15 ),
                    ( 'image', 'final_robust', float, -0.5 ),
                    ( 'image', 'psf_arcsec', float, None ),     # Force restore with this value if set, otherwise use default
                    ( 'image', 'final_psf_arcsec', float, None ),
                    ( 'masking', 'ga', int, 25 ),
                    ( 'masking', 'phase', int, 20 ),
                    ( 'masking', 'ampphase', int, 10 ),
                    ( 'masking', 'full', int, 5 ),
                    ( 'control', 'quiet', bool, False ),
                    ( 'control', 'logging', str, 'logs' ),
                    ( 'control', 'dryrun', bool, False ),
                    ( 'control', 'restart', bool, True ),
                    ( 'control', 'clearcache', bool, True ),
                    ( 'control', 'bootstrap', bool, False ),
                    ( 'bootstrap', 'use_mpi', bool, False) ,
                    ( 'bootstrap', 'bscell', float, 4.5) ,
                    ( 'bootstrap', 'bsimsize', int, 6000) )

    odict = {}
    config=ConfigParser.SafeConfigParser()
    config.read(filename)
    cased={int: config.getint, float: config.getfloat, bool: config.getboolean, str: config.get}
    for (section, name, otype, default) in option_list:
        try:
            result=cased[otype](section,name)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            result=default
        odict[name]=result

    return odict

if __name__=='__main__':
    
    print options('test.cfg')
