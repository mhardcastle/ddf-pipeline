import os
import glob
import numpy as np

# Code from Tim for fixing symbolic links for DDS3_

def get_solutions_timerange(sols):
    t = np.load(sols)['BeamTimes']
    return np.min(t),np.max(t)

def fixsymlinks(ddsols,workdir='.',stype='smoothed'):
    #dds3smoothed = glob.glob('SOLSDIR/*/*killMS.DDS3_full_smoothed*npz')
    dds3 = glob.glob(workdir+'/SOLSDIR/*/killMS.' + ddsols + '.sols.npz')
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
            print('Symlink ' + symsolname + ' already exists, recreating')
            os.unlink(symsolname)
            os.symlink(os.path.relpath('../../%s_%s_%s.npz'%(ddsols,start_time,stype)),symsolname)
        else:
            print('Symlink ' + symsolname + ' does not yet exist, creating')
            os.symlink(os.path.relpath('../../%s_%s_%s.npz'%(ddsols,start_time,stype)),symsolname)
            
    return

