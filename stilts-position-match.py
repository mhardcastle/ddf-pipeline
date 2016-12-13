#!/opt/local/bin/python2.7
import stilts
import math
import os,sys


cat1 = str(raw_input('catname 1: '))
cat2 = str(raw_input('catname 2: '))
maxdist = float(raw_input('matching distance (arcsec): '))
outname = str(raw_input('outname: '))
           
scat1 = stilts.tread(cat1)
scat2 = stilts.tread(cat2)
catalog = stilts.tskymatch2(in1=scat1, in2=scat2, ra1='RA', dec1='DEC', ra2='RA', dec2='DEC', error=maxdist, join='1and2',find='all')

catalog.write(outname)

