# Remove bootstrap solutions from a list of mss

import sys
import pyrap.tables as pt

def remove_columns(mslist_name):

    print 'Removing SCALED_DATA column in',mslist_name
    mslist=[s.strip() for s in open(mslist_name).readlines()]
    for ms in mslist:
        t = pt.table(ms)
        try:
            dummy=t.getcoldesc('SCALED_DATA')
        except RuntimeError:
            dummy=None
        t.close()
        if dummy is not None:
            print 'Removing SCALED_DATA from',ms
            t=pt.table(ms,readonly=False)
            t.removecols('SCALED_DATA')
            t.close()
        else:
            print 'Table',ms,'has no SCALED_DATA column'

if __name__=='__main__':
    remove_columns(sys.argv[1])
