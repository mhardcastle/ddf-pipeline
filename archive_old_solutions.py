# Tidy up solutions from a previous run

from options import options
import os
import sys

def rename(src,dest):
    if not(os.path.isfile(src)):
        print 'file',src,'does not exist!'
    else:
        print 'renaming',src,'as',dest
        os.rename(src,dest)

o=options(sys.argv[1])
if o['mslist'] is None:
    die('MS list must be specified')

with open(o['mslist'],'r') as f:
    msnames=[l.strip() for l in f.readlines()]

# find version number in case this has been run already

prefix='old'
count=1
while os.path.isfile(msnames[0]+'/'+prefix+'.killms_p1.sols.npz'):
    count+=1
    prefix='old%i' % count

print 'prefix is',prefix

for ms in msnames:
    for s in ['p1','ap1']:
        for file in ['npz','parset']:
            rename(ms+'/killMS.killms_'+s+'.sols.'+file,ms+'/'+prefix+'.killms_'+s+'.sols.'+file)
    
if o['full_mslist'] is not None:
    with open(o['full_mslist'],'r') as f:
        msnames=[l.strip() for l in f.readlines()]
    for ms in msnames:
        s='f_ap1'
        for file in ['npz','parset']:
            rename(ms+'/killMS.killms_'+s+'.sols.'+file,ms+'/'+prefix+'.killms_'+s+'.sols.'+file)
