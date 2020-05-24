from __future__ import print_function
from __future__ import absolute_import
# Options for killms pipeline code

from future import standard_library
standard_library.install_aliases()
#from builtins import str
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os
import re
import sys
from termsize import get_terminal_size_linux

def options(optlist,option_list):

    # option_list format is: section, name, type, default
    # section names are used in the output dict only if names are not unique

    odict = {}
    config=configparser.SafeConfigParser()
    filenames=[]
    cmdlineset=[]
    if isinstance(optlist,str):
        optlist=[optlist]

    for o in optlist:
        if o[:2]=='--':
            optstring=o[2:]
            result=re.match('(\w*)-(\w*)\s*=\s*(.*)',optstring)
            if result is None:
                print('Cannot parse option',optstring)
            else:
                cmdlineset.append(result.groups())
        else:
            filenames.append(o)

    for f in filenames:
        if not os.path.isfile(f):
            print('Config file',f,'does not exist!')
            sys.exit(1)

    config.read(filenames)
    for c in cmdlineset:
        for o in option_list:
            (section, name, otype, default)=o[:4]
            if c[0]==section and c[1]==name:
                break
        else:
            print('Option %s-%s does not exist!' % (c[0],c[1]))
            sys.exit(2)
        try:
            config.add_section(c[0])
        except configparser.DuplicateSectionError:
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
        except (configparser.NoSectionError, configparser.NoOptionError):
            result=default
        # python2 compatibility, force back to str
        try:
            if isinstance(result,unicode):
                result=result.encode("utf-8")
        except NameError: # unicode type doesn't exist, we are in py3
            pass
        if otype is str and result=="None":
            result=None
        if count>1:
            odict[section+'_'+name]=result
        else:
            odict[name]=result
    return odict

def typename(s):
    return str(s).replace("type '","").replace("class ","").replace("'","")

def print_options(option_list):
    import textwrap
    from auxcodes import bcolors
    # expected to be called if a config file is not specified. Print a
    # list of options
    option_list=sorted(option_list,key=lambda x:x[1])
    width,height=get_terminal_size_linux()
    if width is None:
        width=80
    sections=sorted(set(x[0] for x in option_list))
    klen=max([len(x[1]) for x in option_list])
    tlen=max([len(typename(x[2])) for x in option_list])
    fstring='%-'+str(klen)+'s = %-'+str(tlen)+'s (default %s)'
    indent=' '*(klen+3)
    for s in sections:
        print(bcolors.OKBLUE+'\n[%s]' % s+bcolors.ENDC)
        for o in option_list:
            if len(o)==4:
                (section, name, otype, default)=o
                doc=None
            elif len(o)==5:
                (section, name, otype, default, doc)=o
            else:
                print('Oops!',o)
                continue
            if section==s:
                
                print(bcolors.BOLD+fstring % (name, typename(otype), str(default))+bcolors.ENDC)
                if doc is not None:
                    print(textwrap.fill(doc,width-1,initial_indent=indent,subsequent_indent=indent))

if __name__=='__main__':
    import sys
    from parset import option_list
    config=sys.argv[1:]
    print(options(config,option_list))

