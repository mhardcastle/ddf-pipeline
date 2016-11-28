#!/usr/bin/python

import glob
import os

def average():
    mss=glob.glob('*.pre-cal.ms')
    for m in mss:
        print 'Averaging',m
        outms=m.replace('pre-cal.ms','ave.ms')
        outfile=open('NDPPP.in','w')
        outfile.write('msin=['+m+']\nmsin.datacolumn = CORRECTED_DATA\nmsin.baseline = [CR]S*\nmsout='+outms+'\nmsout.datacolumn = DATA\nsteps = [count,avg]\navg.type = average\navg.freqstep = 2\navg.timestep = 2\n')
        outfile.close()
        os.system('NDPPP NDPPP.in')


if __name__=='__main__':
    average()
