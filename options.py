# Options for killms pipeline code

import ConfigParser
import os
import struct
import re

def _get_terminal_size_linux():
    ''' From https://gist.github.com/jtriley/1108174 '''
    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import termios
            cr = struct.unpack('hh',
                               fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
            return cr
        except:
            pass
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            return None
    return int(cr[1]), int(cr[0])

def get_physical_cpus():
    # find total number of physical cores, ignoring hyperthreading
    lines=open('/proc/cpuinfo').readlines()
    cpus=[]
    for l in lines:
        bits=l.split(':')
        if 'physical id' in l:
            phys=int(bits[1])
        if 'core id' in l:
            core=int(bits[1])
        if l=='\n':
            cpus.append((phys,core))
    cpus=set(cpus)
    return len(cpus)

def getcpus():
    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile:
        lines=len(open(nodefile).readlines())
        return lines
    else:
        return get_physical_cpus()

option_list = ( ( 'machine', 'NCPU_DDF', int, getcpus(),
                  'Number of CPUS to use for DDF'),
                ( 'machine', 'NCPU_killms', int, getcpus(),
                  'Number of CPUS to use for KillMS' ),
                ( 'data', 'mslist', str, None,
                  'Initial measurement set list to use -- must be specified' ),
                ( 'data', 'full_mslist', str, None,
                  'Full-bandwidth measurement set to use for final step, if any' ),
                ( 'data', 'colname', str, 'CORRECTED_DATA', 'MS column to use' ),
                ( 'solutions', 'ndir', int, 45, 'Number of directions' ),
                ( 'solutions', 'NChanSols', int, 1, 'NChanSols for killMS' ),
                ( 'solutions', 'dt', float, 1., 'Time interval for killMS (minutes)' ),
                ( 'solutions', 'LambdaKF', float, 0.5, 'Kalman filter lambda for killMS' ),
                ( 'solutions', 'NIterKF', list, [1, 6, 6], 'Kalman filter iterations for killMS for the three self-cal steps' ),
                ( 'solutions', 'PowerSmooth', float, 1.0, 'Underweighting factor for missing baseliness in killMS'),
                ( 'solutions', 'normalize', list, ['BLBased', 'BLBased', 'BLBased'], 'How to normalize solutions for the three self-cal steps' ),
                ( 'solutions', 'uvmin', float, None, 'Minimum baseline length to use in self-calibration (km)' ),
                ( 'solutions', 'auto_uvmin', bool, False, 'Optimize uv distance in self-calibration automatically from model and measurement sets. If uvmin is not None, use this value as a lower bound.' ),
                ( 'solutions', 'wtuv', float, None, 'Factor to apply to fitting weights of data below uvmin. None implies, effectively, zero.'),
                ( 'solutions', 'robust', float, None, 'Briggs robustness to use in killMS. If None, natural weighting is used.'),
                ( 'solutions', 'smoothing', int, None, 'Smoothing interval for amplitudes, in units of dt. Must be odd.'),
                ( 'image', 'imsize', int, 20000, 'Image size in pixels' ),
                ( 'image', 'cellsize', float, 1.5, 'Pixel size in arcsec' ),
                ( 'image', 'robust', float, -0.15, 'Imaging robustness' ),
                ( 'image', 'final_robust', float, -0.5, 'Final imaging robustness' ),
                ( 'image', 'psf_arcsec', float, None, 'Force restore with this PSF size in arcsec if set, otherwise use default' ),
                ( 'image', 'final_psf_arcsec', float, None, 'Final image restored with this PSF size in arcsec' ),
                ( 'image', 'final_psf_minor_arcsec', float, None, 'Final image restored with this PSF minor axis in arcsec' ),
                ( 'image', 'final_psf_pa_deg', float, None, 'Final image restored with PSF with this PA in degrees' ),
                ( 'image', 'low_psf_arcsec', float, None, 'Low-resolution restoring beam in arcsec' ),
                ( 'image', 'low_robust', float, -0.20, 'Low-resolution image robustness' ),
                ( 'image', 'low_cell', float, 4.5, 'Low-resolution image pixel size in arcsec' ),
                ( 'image', 'low_imsize', int, None, 'Low-resolution image size in pixels' ),
                ( 'image', 'do_decorr', bool, True, 'Use DDF\'s decorrelation mode' ),
                ( 'image', 'HMPsize', int, 10, 'Island size to use HMP initialization' ),
                ( 'image', 'uvmin', float, 0.1, 'Minimum baseline length to use in imaging (km)'),
                ( 'image', 'apply_weights', list, [False, True, True, True], 'Use IMAGING_WEIGHT column from killms'),
                ( 'masking', 'thresholds', list, [25,20,10,5],
                  'sigmas to use in (auto)masking for initial clean and 3 self-cals'),
                ( 'masking', 'tgss', str, None, 'Path to TGSS catalogue file' ),
                ( 'masking', 'tgss_radius', float, 8.0, 'TGSS mask radius in pixels' ), 
                ( 'masking', 'tgss_flux', float, 300, 'Use TGSS components with peak flux in catalogue units (mJy) above this value' ),
                ( 'masking', 'tgss_extended', bool, False, 'Make extended regions for non-pointlike TGSS sources' ),
                ( 'masking', 'tgss_pointlike', float, 30, 'TGSS source considered pointlike if below this size in arcsec' ),
                ( 'masking', 'region', str, None, 'ds9 region to merge with mask'),
                ( 'masking', 'extended_size', int, None,
                  'If generating a mask from the bootstrap low-res images, use islands larger than this size in pixels' ),
                ( 'masking', 'extended_rms', float, 3.0,
                  'Threshold value defining an island in the extended mask'),
                ( 'masking', 'rmsfacet', bool, False, 'If True calculate one rms per facet rather than per image when making the extended rms maps' ),  
                ( 'control', 'quiet', bool, False, 'If True, do not log to screen' ),
                ( 'control', 'nobar', bool, False, 'If True, do not print progress bars' ),
                ( 'control', 'logging', str, 'logs', 'Name of directory to save logs to, or \'None\' for no logging' ),
                ( 'control', 'dryrun', bool, False, 'If True, don\'t run anything, just print what would be run' ),
                ( 'control', 'restart', bool, True, 'If True, skip steps that would re-generate existing files' ),
                ( 'control', 'cache_dir', str, None, 'Directory for ddf cache files -- default is working directory'),
                ( 'control', 'clearcache', bool, True, 'If True, clear all DDF cache before running' ),
                ( 'control', 'bootstrap', bool, False, 'If True, do bootstrap' ),
                ( 'control', 'second_selfcal', bool, False, 'If True, do second round of selfcal on full bandwidth' ),
                ( 'control', 'catch_signal', bool, True, 'If True, catch SIGUSR1 as graceful exit signal -- stops when control returns to the pipeline.'),
                ( 'control', 'exitafter', str, None, 'Step to exit after -- dirin, phase, ampphase, fulllow'),
                ( 'bootstrap', 'bscell', float, 4.5, 'Bootstrap image cell size') ,
                ( 'bootstrap', 'bsimsize', int, 6000, 'Bootstrap image size' ) ,
                ( 'bootstrap', 'catalogues', list, None, 'File names of catalogues for doing bootstrap' ),
                ( 'bootstrap', 'groups', list, None, 'Group numbers for catalogues. At least one match must be found in each group. Optional -- if not present each catalogue is in a different group.' ), 
                ( 'bootstrap', 'frequencies', list, None, 'Frequencies for catalogues (Hz)' ), 
                ( 'bootstrap', 'names', list, None, 'Short names for catalogues' ), 
                ( 'bootstrap', 'radii', list, None, 'Crossmatch radii for catalogues (arcsec)' ),
                ( 'offsets', 'method', str, None, 'Offset correction method to use. None -- no correction'),
                ( 'offsets', 'fit', str, 'mcmc', 'Histogram fit method' ),
                ( 'offsets', 'mode', str, 'normal', 'Mode of operation: normal or test' ) )

def options(optlist):

    # option_list format is: section, name, type, default
    # section names are used in the output dict only if names are not unique

    odict = {}
    config=ConfigParser.SafeConfigParser()
    filenames=[]
    cmdlineset=[]
    if isinstance(optlist,str):
        optlist=[optlist]

    for o in optlist:
        if o[:2]=='--':
            optstring=o[2:]
            result=re.match('(\w*)-(\w*)\s*=\s*(.*)',optstring)
            if result is None:
                print 'Cannot parse option',optstring
            else:
                cmdlineset.append(result.groups())
        else:
            filenames.append(o)

    config.read(filenames)
    for c in cmdlineset:
        try:
            config.add_section(c[0])
        except ConfigParser.DuplicateSectionError:
            pass
        config.set(c[0],c[1],c[2])
    cased={int: config.getint, float: config.getfloat, bool: config.getboolean, str: config.get, list: lambda x,y: eval(config.get(x,y))}
    for o in option_list:
        (section, name, otype, default)=o[:4]
        # if this name is duplicated in another section, we need to know
        count=0
        for o2 in option_list:
            if o2[1]==name: count+=1
        # get result
        try:
            result=cased[otype](section,name)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            result=default
        if count>1:
            odict[section+'_'+name]=result
        else:
            odict[name]=result
    if odict['logging']=='None':
        odict['logging']=None
    return odict

def typename(s):
    return str(s).replace("type '","").replace("'","")

def print_options():
    import textwrap
    from auxcodes import bcolors
    # expected to be called if a config file is not specified. Print a
    # list of options
    width,height=_get_terminal_size_linux()
    sections=sorted(set(x[0] for x in option_list))
    klen=max([len(x[1]) for x in option_list])
    tlen=max([len(typename(x[2])) for x in option_list])
    fstring='%-'+str(klen)+'s = %-'+str(tlen)+'s (default %s)'
    indent=' '*(klen+3)
    for s in sections:
        print bcolors.OKBLUE+'\n[%s]' % s+bcolors.ENDC
        for o in option_list:
            if len(o)==4:
                (section, name, otype, default)=o
                doc=None
            elif len(o)==5:
                (section, name, otype, default, doc)=o
            else:
                print 'Oops!',o
                continue
            if section==s:
                
                print bcolors.BOLD+fstring % (name, typename(otype), str(default))+bcolors.ENDC
                if doc is not None:
                    print textwrap.fill(doc,width-1,initial_indent=indent,subsequent_indent=indent)

if __name__=='__main__':
    import sys
    config=sys.argv[1:]
    print options(config)

