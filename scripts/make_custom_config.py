#!/usr/bin/env python

import os
import sys
from surveys_db import SurveysDB

def make_custom_config(name,workdir,do_field,averaged=False):
    if do_field:
        with SurveysDB() as sdb:
            idd=sdb.get_field(name)
            no_wenss=((idd['decl']<32) | (idd['decl']>72))
            no_tgss=(idd['no_tgss']==1)
    else:
        no_wenss=False
        no_tgss=False

    if no_wenss:
        template=os.environ['DDF_DIR']+'/ddf-pipeline/examples/tier1-jul2018-NVSS.cfg'
    else:
        template=os.environ['DDF_DIR']+'/ddf-pipeline/examples/tier1-jul2018.cfg'

    lines=open(template).readlines()
    outfile=open(workdir+'/tier1-config.cfg','w')
    for l in lines:
        if 'colname' in l and averaged:
            outfile.write('colname=DATA\n')
        elif '[control]' in l and no_tgss:
            outfile.write(l+'redo_DI=True\n')
        else:
            outfile.write(l)

if __name__=='__main__':
    averaged=False
    try:
        if sys.argv[2]:
            averaged=True
    except:
        pass
    make_custom_config(sys.argv[1],'.',True,averaged)
    
