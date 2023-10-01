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
        print('Symsolname is',symsolname,'and solname is',solname)
        start_time,t1 = get_solutions_timerange(solname)
        # Rounding different on different computers which is a pain.
        # find start time generally for any type of filename
        # this is thoroughly broken, so
        # find all the files of this form and try to find one that matches!
        globst=workdir+'/%s_*_%s.npz'%(ddsols,stype)
        print('Looking for filenames of the form',globst)
        for filename in glob.glob(globst):
            bits=filename.split('_')
            f_start_time=None
            for b in bits[1:]:
                try:
                    f_start_time=float(b)
                    break
                except ValueError:
                    pass
            offset=f_start_time-start_time
            if np.abs(offset)<3600: # adjust to taste
                print('Taking %s to be a match (offset %.2f seconds)' % (filename,offset))
                break
        else:
            raise RuntimeError('Failed to find a match!')
                
        if os.path.islink(symsolname):
            if verbose: print('Symlink ' + symsolname + ' already exists, recreating')
            os.unlink(symsolname)
            os.symlink(os.path.relpath('../../'+filename),symsolname)
        else:
            if verbose: print('Symlink ' + symsolname + ' does not yet exist, creating')
            if os.path.isfile(symsolname):
                if verbose: print('Deleting existing real file in this location, since you asked me to!')
                os.unlink(symsolname)
            os.symlink(os.path.relpath('../../'+filename),symsolname)
            
    return

if __name__=='__main__':
    fixsymlinks('DDS3_full',delete_existing=True)
    
