import os
from mpi4py import MPI
import mpi4py.futures
import concurrent.futures
import numpy as np
import time
from pathlib import Path
import argparse
import shlex
import glob

#Copied from kMS.py
def read_options():
    OP = argparse.ArgumentParser(description="Parser for kMS.py arguments")
    
    OP.add_argument('--MSName',type=str,help='Input MS to draw [no default]')
    OP.add_argument('--TChunk',help='Time Chunk in hours. ')
    OP.add_argument('--InCol',help='Column to work on. ')
    OP.add_argument('--OutCol',help='Column to write to. ')
    OP.add_argument('--FreePredictColName',type=str,help=' . ')
    OP.add_argument('--FreePredictGainColName',type=str,help=' . ')
    OP.add_argument('--Parallel',type=int,help=' . ')

    OP.add_argument('--SkyModel',help='List of targets [no default]')
    OP.add_argument('--kills',help='Name or number index of sources to kill')
    OP.add_argument('--invert',help='Invert the selected sources to kill')
    OP.add_argument('--Decorrelation',type=str,help=' . ')
    OP.add_argument('--FreeFullSub',type=int,help=' . ')
    OP.add_argument('--SkyModelCol',type=str,help=' . ')

    OP.add_argument('--BaseImageName')
    OP.add_argument('--ImagePredictParset')
    OP.add_argument('--DicoModel')
    OP.add_argument('--OverS')
    OP.add_argument('--wmax')
    OP.add_argument('--MaskImage')
    OP.add_argument('--NodesFile')
    OP.add_argument('--MaxFacetSize')
    OP.add_argument('--MinFacetSize')
    OP.add_argument('--DDFCacheDir')
    OP.add_argument('--RemoveDDFCache')
    OP.add_argument('--FilterNegComp')
    OP.add_argument('--ThSolve',help="If the tessel has an apparant SumFlux bellow ThSolve*MaxSumFlux (Max over tessels), it will be unsolved (J=1)")

    OP.add_argument('--CompressionMode',help='Only Auto implemented. ')
    OP.add_argument('--CompressionDirFile',help='Directions in which to do the compression.')
    OP.add_argument('--MergeStations',help='Merge stations into a single one. Use --MergeStations=[CS] to merge all core stations.')
    
    OP.add_argument('--UVMinMax',help='Baseline length selection in km. For example UVMinMax=0.1,100 selects baseline with length between 100 m and 100 km. ')
    OP.add_argument('--ChanSlice',type=str,help='Channel selection option. ')
    OP.add_argument('--FlagAnts',type=str,help='FlagAntenna patern. ')
    OP.add_argument('--DistMaxToCore',type=float,help='Maximum distance to core in km. ')
    OP.add_argument('--FillFactor',type=float)
    OP.add_argument('--FieldID',type=int)
    OP.add_argument('--DDID',type=int)

    OP.add_argument('--BeamModel',type=str,help='Apply beam model, Can be set to: None/LOFAR. ')
    OP.add_argument('--BeamAt',type=str,help='Where to apply beam model, Can be set to: tessel/facet. ')
    OP.add_argument('--LOFARBeamMode',type=str,help='LOFAR beam mode. "AE" sets the beam model to Array and Element. ')
    OP.add_argument('--DtBeamMin',type=float,help='Estimate the beam every this interval [in minutes]. ')
    OP.add_argument('--CenterNorm',type=str,help='Normalise the beam at the field center. ')
    OP.add_argument('--NChanBeamPerMS',type=int,help='Number of channel in the Beam Jones matrix. ')
    OP.add_argument('--FITSParAngleIncDeg',type=float,help='Estimate the beam every this PA change [in deg]. ')
    OP.add_argument('--FITSFile',type=str,help='FITS beam mode filename template. ')
    OP.add_argument('--FITSLAxis',type=str,help='L axis of FITS beam. ')
    OP.add_argument('--FITSMAxis',type=str,help='L axis of FITS beam. ')
    OP.add_argument('--FITSFeed',type=str,help='FITS feed. xy or rl or None to take from MS. ')
    OP.add_argument('--FITSFeedSwap',type=int,help='Swap the feeds around. ')
    OP.add_argument('--FITSVerbosity',type=int,help='Verbosity of debug messages. ')
    OP.add_argument('--ApplyPJones',type=int,help='derotate visibility data (only when FITS beam is active and also time sampled)')
    OP.add_argument('--FlipVisibilityHands',type=int,help='apply anti-diagonal matrix if FITS beam is enabled effectively swapping X and Y or R and L and their respective hands')
    OP.add_argument('--FeedAngle',type=float,help='offset feed angle to add to parallactic angle')
    OP.add_argument('--FITSFrame', type=str, help=' coordinate frame for FITS beams. Currently, alt-az, equatorial and zenith mounts are supported. #options:altaz|altazgeo|equatorial|zenith . ')

    OP.add_argument('--PreApplySols',type=str,help='Pre-apply killMS solutions in the predict step. Has to be a list. ')
    OP.add_argument('--PreApplyMode',type=str,help='Mode for the pre-applied killMS solutions ("A", "P" and "AP" for Amplitude, Phase and Amplitude+Phase). Has to be a list. ')


    OP.add_argument('--Resolution',type=float,help='Resolution in arcsec. ')
    OP.add_argument('--WeightInCol',type=str,help='Weighting column to take into account to weight the visibilities in the solver. ')
    OP.add_argument('--Weighting',type=str,help='Weighting scheme. ')
    OP.add_argument('--Robust',type=float,help='Briggs Robust parameter. ')
    OP.add_argument('--WeightUVMinMax',help='Baseline length selection in km for full weight. For example WeightUVMinMax=0.1,100 selects baseline with length between 100 m and 100 km. ')
    OP.add_argument('--WTUV',type=float,help='Scaling factor to apply to weights outside range of WeightUVMinMax. ')
    
    OP.add_argument('--DoPlot',type=int,help='Plot the solutions, for debugging. ')
    OP.add_argument('--SubOnly',type=int,help='Subtract selected sources. ')
    OP.add_argument('--DoBar',help=' Draw progressbar. ',default="1")
    OP.add_argument('--NCPU',type=int,help='Number of cores to use.  ')
    OP.add_argument('--NThread',type=int,help='Number of OMP/BLAS/etc. threads to use.  ', default=1)
    OP.add_argument('--UpdateWeights',help='Update imaging weights. ',default="1")
    OP.add_argument('--DebugPdb',type=int,help='Drop into Pdb on error. ')

    OP.add_argument('--ExtSols',type=str,help='External solution file. If set, will not solve.')
    OP.add_argument('--ApplyMode',type=str,help='Subtract selected sources. ')
    OP.add_argument('--ClipMethod',type=str,help='Clip data in the IMAGING_WEIGHT column. Can be set to Resid, DDEResid or ResidAnt . ')
    OP.add_argument('--OutSolsName',type=str,help='If specified will save the estimated solutions in this file. ')
    OP.add_argument('--ApplyToDir',type=int,help='Apply direction averaged gains to residual data in the mentioned direction. If ApplyCal=-1 takes the mean gain over directions. -2 if off. ')
    OP.add_argument('--MergeBeamToAppliedSol',type=int,help='Use the beam in applied solution. ')
    OP.add_argument('--SkipExistingSols',type=int,help='Skipping existing solutions if they exist. ')
    OP.add_argument('--SolsDir',type=str,help='Directory in which to save the solutions. ')
    
    OP.add_argument('--SolverType',help='Name of the solver to use (CohJones/KAFCA)')
    OP.add_argument('--PrecisionDot',help='Dot product Precision (S/D). .',type=str)
    OP.add_argument('--PolMode',help='Polarisation mode (Scalar/IFull).')
    OP.add_argument('--dt',type=float,help='Time interval for a solution [minutes].')
    OP.add_argument('--NChanSols',type=int,help='Number of solutions along frequency axis.')
    
    OP.add_argument('--NIterLM',type=int,help=' Number of iterations for the solve.')
    OP.add_argument('--LambdaLM',type=float,help=' Lambda parameter for CohJones.')
    OP.add_argument('--LambdaTk',type=float,help=' Tikhonov regularisation parameter.')
    

    OP.add_argument('--NIterKF',type=int,help=' Number of iterations for the solve.')
    OP.add_argument('--LambdaKF',type=float,help=' Lambda parameter for KAFCA.')
    OP.add_argument('--InitLM',type=int,help='Initialise Kalman filter with Levenberg Maquardt. ')
    OP.add_argument('--InitLMdt',type=float,help='Time interval in minutes. ')
    OP.add_argument('--CovP',type=float,help='Initial prior Covariance in fraction of the initial gain amplitude. ') 
    OP.add_argument('--CovQ',type=float,help='Intrinsic process Covariance in fraction of the initial gain amplitude. ') 
    OP.add_argument('--PowerSmooth',type=float,help='When an antenna has missing baselines (like when using UVcuts) underweight its Q matrix. ') 
    OP.add_argument('--evPStep',type=int,help='Start calculation evP every evP_Step after that step. ')
    OP.add_argument('--evPStepStart',type=int,help='Calculate (I-KJ) matrix every evP_Step steps. ')
    OP.add_argument('--EvolutionSolFile',type=str,help='Evolution solution file. ')

    #FBO specific
    OP.add_argument('--KMDir',type=str,required=True,help='The prefix for MPI jobs directories')
    OP.add_argument('--MSTESTDir',type=str,required=True,help='Where to find the MS and result from previous commands')
    OP.add_argument('--Container',type=str,required=True,help='Name of the container to run the cmd in')
    OP.add_argument('--CPath',type=str,required=True,help='Name of the container to run the cmd in')
    return OP.parse_args()

def reconstruct_command(args):
    args_dict = vars(args)
    cmd = "kMS.py"
    for key, value in args_dict.items():
        if key not in ["KMDir","MSTESTDir","MSName","Container","CPath"] and value is not None:
            print(f"{key}: {value}")
            cmd += '  --'+key+' '+str(value)
    return cmd

def killms_job(args,base_cmd,ms,dry_run=False):

    
    kms_wdir= args.KMDir+ms+".tmp"
    #os.makedirs(kms_wdir, exist_ok=True)
    worker_sh = 'if [ -d "'+kms_wdir+'" ]; then'+'\n'
    worker_sh += '   rm -rf '+kms_wdir+'\n'
    worker_sh += 'fi'+'\n'
    worker_sh += "mkdir "+kms_wdir+'\n'
    tp = args.CPath
    wdir = args.MSTESTDir
    files1 = [f.name for f in Path(wdir).glob('*.fits')]
    files2 = [f.name for f in Path(wdir).glob('*.DicoModel')]
    files3 = [f.name for f in Path(wdir).glob('*.npy')]
    directories = [entry.name for entry in os.scandir(wdir) if entry.is_dir() and entry.name.endswith('.MS')]
    files = files1 + files2 + files3 + directories
    worker_sh += "CUR_DIR=$(pwd)"+'\n'
    worker_sh += "echo $CUR_DIR"+'\n'
    worker_sh += "cd "+kms_wdir+'\n'

    for f in files:
        target = os.path.join(wdir,f)  # The actual file or directory      
        link_name = os.path.join(tp,kms_wdir,f)
        print("link name",link_name,"target",target)
        try:
            #os.symlink(target, link_name)
            worker_sh += "ln -s $CUR_DIR/"+target+'\n'
        except Exception as e:
            print("Symb Link error",str(e))
            pass
    #allowed_dir = os.environ.get("SINGULARITY_ALLOWED_DIR")
    worker_sh += "/usr/local/src/killMS/killMS/"+base_cmd + ' --MSName '+ms+'\n'
    worker_sh +="cd $CUR_DIR"+'\n'
    shell_name = "kms_worker_"+ms+".sh"
    fw = open(shell_name,"w")
    fw.write(worker_sh)
    fw.close()
    cmd = 'singularity run -B .:'+tp+' ${SINGULARITY_ALLOWED_DIR}/'+args.Container+' bash ./'+shell_name
    print("Running killMS: ",cmd)
    if not dry_run:
        os.system(cmd)
    else:
        print("EXECUTING",cmd)
   
if __name__ == '__main__':
    
    args = read_options()
    base_command = reconstruct_command(args)
    #print(base_command)

    work_dir =  args.MSTESTDir
    print('work_dir=',work_dir)
    ms_dirs = [os.path.basename(d) for d in glob.glob(os.path.join(work_dir, "*.MS")) if os.path.isdir(d)]
    print(ms_dirs)
    #Creation des répertoires de run de kMS
    #for ms in ms_dirs:
    #    killms_job(args,base_command,ms,dry_run=True)


    with mpi4py.futures.MPIPoolExecutor() as mpi_executor:
        size_mpi = MPI.COMM_WORLD.Get_size()
        killmsfs = []
        for i in range(1,len(ms_dirs)):
            killmsfs.append(mpi_executor.submit(killms_job,args,base_command,ms_dirs[i],dry_run=True))
        killms_job(args,base_command,ms_dirs[0],dry_run=True)
        mpi4py.futures.wait(killmsfs)








          
