#!/usr/bin/python

import glob
import os

def unpack():
    # unpack all the files in current working directory
    destdir='prefactor/results/'
    files=glob.glob('*.tar.gz')
    for f in files:
        print 'Unpacking',f
        os.system('tar xf '+f)
        os.system('mv '+destdir+'* .')

if __name__=='__main__':
    unpack()
