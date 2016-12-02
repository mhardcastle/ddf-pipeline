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

option_list = ( ( 'machine', 'NCPU_DDF', int, getcpus() ),
                ( 'machine', 'NCPU_killms', int, getcpus() ),
                ( 'data', 'mslist', str, None ),
                ( 'data', 'full_mslist', str, None ),
                ( 'data', 'colname', str, 'CORRECTED_DATA' ),
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
                ( 'image', 'low_psf_arcsec', float, None ),
                ( 'image', 'low_robust', float, -0.20 ),
                ( 'image', 'low_cell', float, 4.5 ),
                ( 'image', 'low_imsize', int, None ),
                ( 'image', 'do_decorr', bool, True ),
                ( 'masking', 'ga', int, 25 ),
                ( 'masking', 'phase', int, 20 ),
                ( 'masking', 'ampphase', int, 10 ),
                ( 'masking', 'full', int, 5 ),
                ( 'masking', 'tgss', str, None ),
                ( 'masking', 'tgss_radius', float, 8.0 ), # radius in pix
                ( 'masking', 'tgss_flux', float, 500 ), # peak flux in mJy
                ( 'control', 'quiet', bool, False ),
                ( 'control', 'logging', str, 'logs' ),
                ( 'control', 'dryrun', bool, False ),
                ( 'control', 'restart', bool, True ),
                ( 'control', 'clearcache', bool, True ),
                ( 'control', 'bootstrap', bool, False ),
                ( 'control', 'stagedir', str, None ),
                ( 'bootstrap', 'use_mpi', bool, False) ,
                ( 'bootstrap', 'bscell', float, 4.5) ,
                ( 'bootstrap', 'bsimsize', int, 6000 ) ,
                ( 'bootstrap', 'groups', list, None ), 
                ( 'bootstrap', 'frequencies', list, None ), 
                ( 'bootstrap', 'names', list, None ), 
                ( 'bootstrap', 'radii', list, None ), 
                ( 'bootstrap', 'catalogues', list, None ) )

def options(filename):

    # option_list format is: section, name, type, default
    # names must be unique -- section names are not used in output dict

    odict = {}
    config=ConfigParser.SafeConfigParser()
    config.read(filename)
    cased={int: config.getint, float: config.getfloat, bool: config.getboolean, str: config.get, list: lambda x,y: eval(config.get(x,y))}
    for (section, name, otype, default) in option_list:
        try:
            result=cased[otype](section,name)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            result=default
        odict[name]=result
    if odict['logging']=='None':
        odict['logging']=None
    return odict

def typename(s):
    return str(s).replace("type '","").replace("'","")

def print_options():
    # expected to be called if a config file is not specified. Print a
    # list of options
    sections=set(x[0] for x in option_list)
    for s in sections:
        print '\n[%s]' % s
        for (section, name, otype, default) in option_list:
            if section==s:
                print '%-16s = %-20s (default %s)' % (name, typename(otype), str(default))

if __name__=='__main__':
    
    print options('example.cfg')
    #print_options()
