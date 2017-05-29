import os,sys
import numpy as np
import argparse

def sepn(r1,d1,r2,d2):
    """
    Calculate the separation between 2 sources, RA and Dec must be
    given in radians. Returns the separation in radians [TWS]
    """
    cos_sepn=np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
    sepn = np.arccos(cos_sepn)
    return sepn    

def read_pointingfile(pointingfilename):

	pointingdict = {}
	infile = open(pointingfilename,'r')
	for line in infile:
		line = line[:-1]
		line = line.split(',')
		while '' in line:
			line.remove('')
		if '#' in line[0]:
			continue
		pointingname,radeg,decdeg,integrationtime,pointingids = line[0],line[1],line[2],line[3],line[4]
		pointingdict[pointingname] = [pointingids,radeg,decdeg,integrationtime]
	return pointingdict

def find_pointings_to_mosaic(pointingdict):
	seps = np.array([])
	tomospointings = np.array([])
	ra1,dec1 = pointingdict[mospointingname][1],pointingdict[mospointingname][2]
	for j in pointingdict:
		ra2,dec2 = pointingdict[j][1],pointingdict[j][2]
		seps = np.append(seps,(sepn(float(ra1)*deg2rad,float(dec1)*deg2rad,float(ra2)*deg2rad,float(dec2)*deg2rad)*rad2deg))
		tomospointings = np.append(tomospointings,j)
	sortedseps = np.sort(seps)
	# For now include the 6 closest pointings plus keep including others until we have all that are within the distance of 1.1 times the 6th closest pointing. At later stage try to take into account e.g. elongation of beam for low declination observations.
	endmosindex = 7
	for i in range(7,20):
		if (sortedseps[i])/(sortedseps[6]) < 1.1:
			endmosindex = i
	mosaicpointings = tomospointings[np.where(seps < sortedseps[endmosindex])]
	mosseps = seps[np.where(seps < sortedseps[endmosindex])]
	return mosaicpointings,mosseps





parser = argparse.ArgumentParser()
parser.add_argument('pointingfile', type=str, help='LoTSS pointing progress file')
parser.add_argument('mospointingname', type=str, help='Mosiac central pointing name')

args = parser.parse_args()
pointingfilename = args.pointingfile
mospointingname = args.mospointingname

deg2rad = np.pi/180.0
rad2deg = 180.0/np.pi

pointingdict = read_pointingfile(pointingfilename)

mosaicpointings,mosseps = find_pointings_to_mosaic(pointingdict)

printline = ''
printlineseps=''
for i in range(0,len(mosaicpointings)):
	printline += (',%s'%mosaicpointings[i])
	printlineseps += (',%s'%mosseps[i])

print 'Separations: %s'%printlineseps[1:]
print 'Pointings: %s'%printline[1:]
