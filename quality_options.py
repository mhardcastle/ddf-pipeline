# Options for quality control pipeline code

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

option_list = ( ( 'machine', 'NCPU', int, getcpus() ),
                ( 'image', 'pbimage', str, 'image_full_ampphase1m.smooth.int.restored.fits' ),
                ( 'image', 'nonpbimage', str, 'image_full_ampphase1m.app.restored.fits' ),
                ( 'image', 'catprefix', str, 'image_full_ampphase1m' ),
                ( 'control', 'sfind', bool, True ),
                ( 'control', 'sfind_pixel_fraction', float, 0.6 ),
                ( 'control', 'quiet', bool, False ),
                ( 'control', 'logging', str, 'logs' ),
                ( 'control', 'dryrun', bool, False ),
                ( 'control', 'restart', bool, True ),
                ( 'comparison_cats', 'list', list, None ),
                ( 'comparison_cats', 'filenames', list, None),
                ( 'comparison_cats', 'radii', list, None),
                ( 'comparison_cats', 'fluxfactor', list, None),
                ( 'comparison_cats', 'TGSS', str, None ),
                ( 'comparison_cats', 'TGSS_matchrad', float, 10.0 ),
                ( 'comparison_cats', 'TGSS_match_majkey1', float, 'Maj_1' ),
                ( 'comparison_cats', 'TGSS_match_majkey2', float, 'Maj_2' ),
                ( 'comparison_cats', 'TGSS_filtersize', float, 40.0 ),
                ( 'comparison_cats', 'TGSS_fluxfactor', float, 1000.0 ),
                ( 'comparison_cats', 'FIRST', str, None ),
                ( 'comparison_cats', 'FIRST_matchrad', float, 10.0 ),
                ( 'comparison_cats', 'FIRST_match_majkey1', float, 'Maj' ),
                ( 'comparison_cats', 'FIRST_match_majkey2', float, 'MAJOR' ),
                ( 'comparison_cats', 'FIRST_filtersize', float, 10.0 ),
                ( 'comparison_cats', 'FIRST_fluxfactor', float, 1.0 ) )

def options(filename):

    # option_list format is: section, name, type, default
    # names must be unique -- section names are not used in output dict

    odict = {}
    config=ConfigParser.SafeConfigParser()
    config.read(filename)
    cased={int: config.getint, float: config.getfloat, bool: config.getboolean, str: config.get, list: lambda x,y: eval(config.get(x
,y))}
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
