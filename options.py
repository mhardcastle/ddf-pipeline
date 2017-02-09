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

option_list = ( ( 'machine', 'NCPU_DDF', int, getcpus(),
                  'Number of CPUS to use for DDF'),
                ( 'machine', 'NCPU_killms', int, getcpus(),
                  'Number of CPUS to use for KillMS' ),
                ( 'data', 'mslist', str, None,
                  'Initial measurement set list to use -- must be specified' ),
                ( 'data', 'full_mslist', str, None,
                  'Full-bandwidth measurement set to use for final step, if any' ),
                ( 'data', 'colname', str, 'CORRECTED_DATA', 'MS column to use' ),
                ( 'solutions', 'ndir', int, 60, 'Number of directions' ),
                ( 'solutions', 'NChanSols', int, 1, None ),
                ( 'solutions', 'dt', float, 1., 'Time interval for killMS' ),
                ( 'solutions', 'LambdaKF', float, 0.5, None ),
                ( 'solutions', 'NIterKF', list, [1, 6, 6], 'Kalman filter iterations for killMS for the three self-cal steps' ),
                ( 'solutions', 'normalize', list, ['AbsAnt', 'AbsAnt', 'Abs'], 'How to normalize solutions for the three self-cal steps' ),
                ( 'image', 'imsize', int, 20000, 'Image size' ),
                ( 'image', 'msmf_threshold', float, 10e-3 ),
                ( 'image', 'cellsize', float, 1.5 ),
                ( 'image', 'robust', float, -0.15 ),
                ( 'image', 'final_robust', float, -0.5 ),
                ( 'image', 'psf_arcsec', float, None, 'Force restore with this PSF size if set, otherwise use default' ),
                ( 'image', 'final_psf_arcsec', float, 'Final image restored with this PSF size' ),
                ( 'image', 'low_psf_arcsec', float, None ),
                ( 'image', 'low_robust', float, -0.20 ),
                ( 'image', 'low_cell', float, 4.5 ),
                ( 'image', 'low_imsize', int, None ),
                ( 'image', 'do_decorr', bool, True ),
                ( 'image', 'HMPsize', int, None ),
                ( 'masking', 'thresholds', list, [25,20,10,5],
                  'sigmas to use in (auto)masking for initial clean and 3 self-cals'),
                ( 'masking', 'tgss', str, None, 'Path to TGSS catalogue file' ),
                ( 'masking', 'tgss_radius', float, 8.0, 'TGSS mask radius in pixels' ), 
                ( 'masking', 'tgss_flux', float, 500, 'Use TGSS components with peak flux in catalogue units (mJy) above this value' ),
                ( 'masking', 'tgss_extended', bool, False, 'Make extended regions for non-pointlike TGSS sources' ),
                ( 'masking', 'tgss_pointlike', float, 30, 'TGSS source considered pointlike if below this size in arcsec' ),
                ( 'masking', 'region', str, None, 'ds9 region to merge with mask'),
                ( 'masking', 'extended_size', int, None ),
                ( 'masking', 'extended_rms', float, 3.0 ), 
                ( 'control', 'quiet', bool, False ),
                ( 'control', 'nobar', bool, False ),
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
    for o in option_list:
        if len(o)==4:
            (section, name, otype, default)=o
        else:
            (section, name, otype, default,_)=o
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
        for o in option_list:
            if len(o)==4:
                (section, name, otype, default)=o
                doc=None
            else:
                (section, name, otype, default, doc)=o
            if section==s:
                print '%-16s = %-10s (default %s)' % (name, typename(otype), str(default))
                if doc is not None:
                    print ' '*18,doc

if __name__=='__main__':
    
    #print options('example.cfg')
    print_options()
