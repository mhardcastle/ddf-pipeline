from __future__ import print_function
import os
import glob
import numpy as np

# Code from Tim for fixing symbolic links

def get_solutions_timerange(sols):
    print('Reading %s'%sols)
    S=np.load(sols)
    t = np.concatenate([S["Sols"]["t0"],S["Sols"]["t1"]])
    # Fix for when first and last time slots are offset.
    tdiff = np.diff(t, 2)
    if not np.allclose(tdiff, 0):
        offset_start = tdiff[0]
        offset_end = tdiff[-1]
        print('Found offset of %d, %d for first and last time slots.'%(offset_start,offset_end))
        # +1 seems needed to match the time in the solution filename.
        t[0] -= offset_start + 1
        t[-1] -= offset_end
    return np.min(t),np.max(t)

def fixsymlinks(ddsols,workdir='.',stype='smoothed',verbose=True,delete_existing=False):
    #dds3smoothed = glob.glob('SOLSDIR/*/*killMS.DDS3_full_smoothed*npz')
    dds3 = glob.glob(workdir+'/SOLSDIR/*/killMS.' + ddsols + '.sols.npz')
    if verbose:
        print('About to process',len(dds3),'solution lines')
    for i in range(0,len(dds3)):
        symsolname = dds3[i].split('killMS.' + ddsols + '.sols.npz')[0] + 'killMS.'+ddsols+'_'+stype+'.sols.npz' 
        solname = dds3[i]

        start_time,t1 = get_solutions_timerange(solname)
        # Rounding different on different computers which is a pain.
        # find start time generally for any type of filename
        filename = glob.glob(workdir+'/%s_%s*_%s.npz'%(ddsols,int(start_time),stype))[0]
        bits=filename.split('_')
        for b in bits[1:]:
            try:
                start_time=str(float(b))
                break
            except ValueError:
                pass
        else:
            raise RuntimeError('Could not find time')
                
        if os.path.islink(symsolname):
            if verbose: print('Symlink ' + symsolname + ' already exists, recreating')
            os.unlink(symsolname)
            os.symlink(os.path.relpath('../../%s_%s_%s.npz'%(ddsols,start_time,stype)),symsolname)
        else:
            if verbose: print('Symlink ' + symsolname + ' does not yet exist, creating')
            if os.path.isfile(symsolname):
                if verbose: print('Deleting existing real file in this location, since you asked me to!')
                os.unlink(symsolname)
            os.symlink(os.path.relpath('../../%s_%s_%s.npz'%(ddsols,start_time,stype)),symsolname)
            
    return

if __name__=='__main__':
    fixsymlinks('DDS3_full',delete_existing=True)
    
