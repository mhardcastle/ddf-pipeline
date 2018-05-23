#!/usr/bin/env python

# Remove bootstrap solutions from a list of mss

import sys
import pyrap.tables as pt

def remove_columns(mslist_name,colnames=['SCALED_DATA']):

    mslist=[s.strip() for s in open(mslist_name).readlines()]
    for ms in mslist:
        t = pt.table(ms)
        for colname in colnames:
            print 'Removing',colname,'column in',mslist_name
            try:
                dummy=t.getcoldesc(colname)
            except RuntimeError:
                dummy=None
            t.close()
            if dummy is not None:
                print 'Removing',colname,' from',ms
                t=pt.table(ms,readonly=False)
                t.removecols(colname)
                t.close()
            else:
                print 'Table',ms,'has no',colname,'column'

if __name__=='__main__':
    remove_columns(sys.argv[1])
