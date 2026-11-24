#!/usr/bin/env python

# Code to bootstrap the flux density scale using killMS/DDF
# Should run standalone or as part of the pipeline

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
import os,sys
import os.path
import copy
import socket
# When launched as a pipeline subprocess, run() prefixes DDF_FORCE_NOT_USE_MPI=1.
# Scrub the inherited launcher env as well, or this process' MPI_Init would
# rejoin the parent job instead of starting a clean singleton.
if os.environ.get("DDF_FORCE_NOT_USE_MPI") == "1":
    for _k in [k for k in os.environ if k.startswith(("OMPI_", "PMI_", "PMIX_", "HYDRA_", "I_MPI_", "PRTE_", "MV2_", "SLURM_"))]:
        del os.environ[_k]
try:
    import bdsf as bdsm
except ImportError:
    from lofar import bdsm
import pyrap.tables as pt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import shutil
from astropy.io import fits

LOCAL_DEV = os.environ.get("DDF_LOCAL_DEV", "0") == "1"
if not LOCAL_DEV:
    from auxcodes import run,warn,die,MSList
    from pipeline import ddf_image,make_external_mask
    import mpi_manager
else:
    from utils.auxcodes import run,warn,die,MSList
    from scripts.pipeline import ddf_image,make_external_mask
    from utils import mpi_manager

def logfilename(s):
    if o['logging'] is not None:
        return o['logging']+'/'+s 
    else:
        return None

def run_bootstrap_mpi(MPI_Manager,o):
    ListJobs=[]
    if o is not None:
        for Node in MPI_Manager.DicoNode2mslist.keys():
            oc=copy.deepcopy(o)
            oc["mslist"]=MPI_Manager.DicoNode2mslist[Node]
            ListJobs.append((Node,run_bootstrap,(oc,),{}))
        
    mpi_manager.callParallel(ListJobs)

def run_bootstrap(o):

    # guess colname. This is necesssary because skip_di means there is
    # no DI corrected column. Will break if e.g. sub_square and
    # skip_di are used together.
    
    if o['skip_di']:
        colname=o['colname']
    else:
        colname='DATA_DI_CORRECTED'
        
    if o['mslist'] is None:
        die('MS list must be specified')

    if o['logging'] is not None and not os.path.isdir(o['logging']):
        os.mkdir(o['logging'])

    # check the data supplied
    if o['frequencies'] is None or o['catalogues'] is None:
        die('Frequencies and catalogues options must be specified')

    if "DDF_PIPELINE_CATALOGS" not in list(os.environ.keys()):
        warn("You need to define the environment variable DDF_PIPELINE_CATALOGS where your catalogs are located")
        sys.exit(2)

    if o["tgss"] is None:
        die('TGSS catalogue path must be specified (tgss option)')

    o["tgss"]=o["tgss"].replace("$$",os.environ["DDF_PIPELINE_CATALOGS"])
    o["catalogues"]=[l.replace("$$",os.environ["DDF_PIPELINE_CATALOGS"]) for l in o["catalogues"]]
    lCat=o["catalogues"]+[o["tgss"]]
    for fCat in lCat:
        if not os.path.isfile(fCat):
            warn("Catalog %s does not exist"%fCat)
            sys.exit(2)

    cl=len(o['catalogues'])
    if o['names'] is None:
        o['names']=[os.path.basename(x).replace('.fits','') for x in o['catalogues']]
    if o['radii'] is None:
        o['radii']=[10]*cl
    if o['groups'] is None:
        o['groups']=list(range(cl))
    if (len(o['frequencies'])!=cl or len(o['radii'])!=cl or
        len(o['names'])!=cl or len(o['groups'])!=cl):
        die('Names, groups, radii and frequencies entries must be the same length as the catalogue list')


    # If dedicated bootstrap keys are given, use those.
    if o['bootstrap_psf_arcsec'] is not None:
        boot_psf_arcsec  = o['bootstrap_psf_arcsec']
        boot_uvrange     = [o['bootstrap_uvmin'], 2.5*206.0 / boot_psf_arcsec]
        boot_robust      = o['bootstrap_robust']
        boot_cell        = o['bootstrap_cell']
        boot_inner_taper = o['bootstrap_inner_taper']
        boot_outer_taper = o['bootstrap_outer_taper']

        # temp bootstrap debug cleanmode option. Allows for different cleanmode than main process,
        # since it is separate and doesn't share any DicoModel data with the rest.
        if o['boot_cleanmode_override'] != 'None':
            boot_cleanmode = o['boot_cleanmode_override']
        else:
            boot_cleanmode = o['cleanmode']
        
        if o['bootstrap_imsize'] is not None:
            boot_imsize = o['bootstrap_imsize']
        else:
            # fallback to using same extent as main image, but differnt cellsize
            boot_imsize = int(o['imsize'] * o['cellsize'] / boot_cell)

    # If no bootstrap keys found, or not given, default to image_low parameters
    else:
        if o['low_psf_arcsec'] is None: die('bootstrap needs low_psf_arcsec set (e.g. 50.) in the [image] section')

        boot_psf_arcsec  = o['low_psf_arcsec']
        boot_uvrange     = [o['image_uvmin'], 2.5*206.0 / boot_psf_arcsec]
        boot_robust      = o['low_robust']
        boot_cell        = o['low_cell']
        boot_inner_taper = o['low_inner_taper']
        boot_outer_taper = o['low_outer_taper']

        boot_cleanmode = o['cleanmode']
        
        if o['low_imsize'] is not None:
            boot_imsize = o['low_imsize'] # allow over-ride
        else:
            boot_imsize = int(o['imsize'] * o['cellsize'] / boot_cell)
    

    # Clear the shared memory
    run('CleanSHM.py',dryrun=o['dryrun'])

    # We use the individual ms in mslist.
    m=MSList(o['mslist'])
    Uobsid = sorted(set(m.obsids))
    
    for obsid in Uobsid:
        
        warn('Running bootstrap for obsid %s' % obsid)

        freqs=[]
        omslist=[]
        for ms,ob,f in zip(m.mss,m.obsids,m.freqs):
            if ob==obsid:
                omslist.append(ms)
                freqs.append(f)

        if len(freqs)<4:
            die('Not enough frequencies to bootstrap. Check your mslist or MS naming scheme')

        # sort to work in frequency order

        freqs,omslist = (list(x) for x in zip(*sorted(zip(freqs, omslist), key=lambda pair: pair[0])))

        for f,ms in zip(freqs,omslist):
            print(ms,f)

        # generate the sorted input mslist
        with open('temp_mslist.txt','w') as f:
            for line in omslist:
                f.write(line+'\n')

        # Clean in cube mode
        # As for the main pipeline, first make a dirty map
        ddf_image('image_bootstrap_'+obsid+'_init','temp_mslist.txt',
                  cleanmask=None,
                  cleanmode=boot_cleanmode,
                  ddsols='DDS0',
                  applysols=o['apply_sols'][6],
                  majorcycles=0,
                  robust=boot_robust,
                  uvrange=boot_uvrange,
                  beamsize=boot_psf_arcsec,
                  imsize=boot_imsize,
                  cellsize=boot_cell,
                  inner_taper=boot_inner_taper,
                  outer_taper=boot_outer_taper,
                  options=o,
                  colname=colname,
                  automask=True,
                  automask_threshold=15,
                  smooth=True,
                  cubemode=True,
                  conditional_clearcache=True,
                  STEP=0)

        external_mask='bootstrap_external_mask.fits'
        make_external_mask(external_mask,'image_bootstrap_'+obsid+'_init.dirty.fits',use_tgss=True,clobber=False,cellsize='low_cell',options=o)

        # Deep clean with this external mask and automasking
        ddf_image('image_bootstrap_'+obsid,'temp_mslist.txt',
                  cleanmask=external_mask,
                  reuse_psf=True,
                  reuse_dirty=True,
                  cleanmode=boot_cleanmode,
                  ddsols='DDS0',
                  applysols=o['apply_sols'][6],
                  majorcycles=5,
                  robust=boot_robust,
                  uvrange=boot_uvrange,
                  beamsize=boot_psf_arcsec,
                  imsize=boot_imsize,
                  cellsize=boot_cell,
                  inner_taper=boot_inner_taper,
                  outer_taper=boot_outer_taper,
                  options=o,
                  colname=colname,
                  automask=True,
                  automask_threshold=15,
                  smooth=True,
                  cubemode=True,
                  wscms_peakfactor=0.01,
                  wscms_rescue_rms_factor=5.0,
                  conditional_clearcache=False,
                  STEP=0)

        if (os.path.isfile('image_bootstrap_'+obsid+'.cube.int.restored.pybdsm.srl') or
            os.path.isfile('image_bootstrap_'+obsid+'.cube.int.restored.pybdsf.srl')):

            warn('Source list exists, skipping source extraction')
        else:
            warn('Running PyBDSM, please wait...')
            img=bdsm.process_image('image_bootstrap_'+obsid+'.cube.int.restored.fits',thresh_pix=5,rms_map=True,atrous_do=True,atrous_jmax=2,group_by_isl=True,rms_box=(80,20), adaptive_rms_box=True, adaptive_thresh=80, rms_box_bright=(35,7),mean_map='zero',spectralindex_do=True,specind_maxchan=1,debug=True,kappa_clip=3,flagchan_rms=False,flagchan_snr=False,incl_chan=True,spline_rank=1)
            # Write out in ASCII to work round bug in pybdsm
            img.write_catalog(catalog_type='srl',format='ascii',incl_chan='true')
            img.export_image(img_type='rms',img_format='fits')

        resfile = obsid+'crossmatch-results-2.npy'
        if o['bootstrap_method']=='rfs':

            # Radio-Flux-Solver bootstrap. 
            # The cube above was made from the same data, though RFS fits it per slice.
            # RFS can use arbitrary catalogs, and is generally more accurate than the classic method.
            
            import importlib.util
            from astropy.wcs import WCS

            if importlib.util.find_spec('flux_calibrator_cli') is None:
                die('RFS bootstrap needs the Radio-Flux-Solver package importable in this environment')

            if o['rfs_catalog_dir'] is not None and not os.path.isdir(o['rfs_catalog_dir']):
                die(f"RFS catalogue directory {o['rfs_catalog_dir']} does not exist")

            # freqs is the same sorted list the spline below will use
            cubefile = 'image_bootstrap_'+obsid+'.cube.int.restored.fits'
            cubeheader = fits.getheader(cubefile)
            cubedata = fits.getdata(cubefile)

            # ensure data is as expected
            if cubeheader['NAXIS'] != 4:
                die(f"RFS bootstrap expects a 4D cube (freq,stokes,dec,ra), got NAXIS={cubeheader['NAXIS']}")
            if cubeheader['NAXIS4'] != len(freqs):
                die(f"Cube has {cubeheader['NAXIS4']} planes but the mslist gives {len(freqs)} frequencies")

            # a clean 2D sky header for the per-band planes
            planeheader = WCS(cubeheader).celestial.to_header()
            for keyword in ['BUNIT','BMAJ','BMIN','BPA','CLEANBMJ','CLEANBMN','CLEANBPA']:
                if keyword in cubeheader:
                    planeheader[keyword]=cubeheader[keyword]
            chanwidth = abs(cubeheader['CDELT4'])

            def loadcorrection(path):
                # Read a band fit: [scale, error] as a valid (1,2) array, else None
                try:
                    result = np.load(path)
                except Exception:
                    return None
                if result.shape == (1,2) and np.all(np.isfinite(result)) and result[0,0] > 0:
                    return result
                return None

            corrections = []
            refitted = []
            for i,freq in enumerate(freqs):
                band = f'{obsid}band{i+1}'
                planefile = os.path.abspath(band+'.fits')
                outdir = os.path.abspath(band)
                npyfile = os.path.join(outdir,band+'.npy')

                planefreq = cubeheader['CRVAL4'] + (i+1 - cubeheader['CRPIX4']) * cubeheader['CDELT4']
                if abs(planefreq - freq) > chanwidth / 2.0:
                    die(f'Cube plane {i+1} is at {planefreq:g} Hz but the mslist says {freq:g} Hz')

                if not o['dryrun'] and not os.path.isfile(planefile):
                    fits.PrimaryHDU(data=cubedata[i,0,:,:], header=planeheader.copy()).writeto(planefile)

                # check if the files already exist
                if os.path.isfile(npyfile):
                    correction = loadcorrection(npyfile)
                    if correction is not None:
                        warn(f'{band} correction exists, skipping RFS fit')
                        corrections.append(correction)
                        continue
                    warn(f'{band} exists but is not valid, refitting')

                warn(f'Running RFS for {band} ({freq/1e6:.6g} MHz)')
                cmd = (f'{sys.executable} -m flux_calibrator_cli {planefile} '
                     f'--freq {float(freq)!r} --freq-unit Hz '
                     f'--catalogs {o["rfs_preset"]} '
                     f'--anchor-name {band} --output-dir {outdir} '
                     f'--save-correction-factors --save-csv --output-name {band} '
                     f'--seed 0 --minimum-position-error 1.5 '
                     f'--minimum-frequency-spacing {2e6:.1f} '
                     f'--unit-fallback --spatial-filter')

                if o['rfs_catalog_dir'] is not None:
                    cmd += f' --reference-catalog-dir {os.path.abspath(o["rfs_catalog_dir"])}'

                if not o['dryrun']:
                    os.makedirs(outdir, exist_ok=True)

                run(cmd, dryrun=o['dryrun'], log=logfilename('RFS-'+band+'.log'))

                if o['dryrun']:
                    continue

                correction = loadcorrection(npyfile)
                if correction is None:
                    die(f'RFS did not produce any results for {band}')

                corrections.append(correction)
                refitted.append(band)

            resfile = obsid + 'crossmatch-results_rfs.npy'
            if len(corrections) == len(freqs) and not o['dryrun']:
                # rebuild unless nothing was refitted and the saved file is still good
                reassemble = refitted or not os.path.isfile(resfile)
                if not reassemble:
                    saved = np.load(resfile)
                    reassemble = (saved.shape != (len(freqs), 2)
                                or not np.all(np.isfinite(saved))
                                or not np.all(saved[:, 0] > 0))
                if reassemble:
                    np.save(resfile, np.vstack(corrections))
                    warn(f'Assembled {resfile}')

        else:

            if not LOCAL_DEV:
                from make_fitting_product import make_catalogue
                import fitting_factors
                import find_outliers
            else:
                from utils.make_fitting_product import make_catalogue
                from utils import fitting_factors,find_outliers

            # generate the fitting product
            if os.path.isfile(obsid+'crossmatch-1.fits'):
                warn('Crossmatch table exists, skipping crossmatch')
            else:
                t = pt.table(omslist[0]+ '::FIELD', readonly=True, ack=False)
                direction = t[0]['PHASE_DIR']
                ra, dec = direction[0]

                if (ra<0):
                    ra+=2*np.pi
                ra*=180.0/np.pi
                dec*=180.0/np.pi

                cats=list(zip(o['catalogues'],o['names'],o['groups'],o['radii']))
                #make_catalogue('image_bootstrap_'+obsid+'.cube.int.restored.pybdsf.srl',ra,dec,2.5,cats,outnameprefix=obsid)
                filenames=['image_bootstrap_'+obsid+'.cube.int.restored.pybdsm.srl','image_bootstrap_'+obsid+'.cube.int.restored.pybdsf.srl']
                for filename in filenames:
                    if os.path.isfile(filename):
                        break
                else:
                    die('Catalogue file does not exist!')

                make_catalogue(filename,ra,dec,2.5,cats,outnameprefix=obsid)

    
            freqlist=open(obsid+'frequencies.txt','w')
            for n,f in zip(o['names'],o['frequencies']):
                freqlist.write('%f %s_Total_flux %s_E_Total_flux False\n' % (f,n,n))
            for i,f in enumerate(freqs):
                freqlist.write('%f Total_flux_ch%i E_Total_flux_ch%i True\n' % (f,i+1,i+1))
            freqlist.close()

            # Now call the fitting code

            if os.path.isfile(obsid+'crossmatch-results-1.npy'):
                warn('Results 1 exists, skipping first fit')
            else:
                fitting_factors.run_all(1, name=obsid)

            nreject=-1 # avoid error if we fail somewhere
            if os.path.isfile(obsid+'crossmatch-2.fits'):
                warn('Second crossmatch exists, skipping outlier rejection')
            else:
                nreject=find_outliers.run_all(1, name=obsid)
    
            if os.path.isfile(obsid+'crossmatch-results-2.npy'):
                warn('Results 2 exists, skipping second fit')
            else:
              if nreject==0:
                  shutil.copyfile(obsid+'crossmatch-results-1.npy',obsid+'crossmatch-results-2.npy')
            if os.path.isfile(obsid+'crossmatch-results-2.npy'):
                warn('Results 2 exists, skipping first fit')
            else:
                fitting_factors.run_all(2, name=obsid)

        # Now apply corrections

        if o['full_mslist'] is None:
            die('Need big mslist to apply corrections')
        if not(o['dryrun']):
            warn('Applying corrections to MS list')
            scale=np.load(resfile)[:,0]
            # InterpolatedUS gives us linear interpolation between points
            # and extrapolation outside it
            spl = InterpolatedUnivariateSpline(freqs, scale, k=1)
            
            bigmslist=[s.strip() for s in open(o['full_mslist']).readlines()]
            obigmslist = [ms for ms in bigmslist if obsid in ms]
            
            for ms in obigmslist:
                t = pt.table(ms)
                try:
                    dummy=t.getcoldesc('SCALED_DATA')
                except RuntimeError:
                    dummy=None
                t.close()
                if dummy is not None:
                    warn('Table '+ms+' has already been corrected, skipping')
                else:
                    # in this version we need to scale both the original data and the data in colname
                    t = pt.table(ms+'::SPECTRAL_WINDOW', readonly=True, ack=False)
                    frq=t[0]['REF_FREQUENCY']
                    factor=spl(frq)
                    print(frq,factor)
                    t=pt.table(ms,readonly=False)
                    if o['do_wide']:
                        warn('Using DATA_SUB column in bootstrap')
                        desc=t.getcoldesc('DATA_SUB')
                    else:
                        desc=t.getcoldesc(o['colname'])
                    desc['name']='SCALED_DATA'
                    t.addcols(desc)
                    if o['do_wide']:
                        warn('Using DATA_SUB column in bootstrap')
                        d=t.getcol('DATA_SUB')
                    else:
                        d=t.getcol(o['colname'])
                    d*=factor
                    t.putcol('SCALED_DATA',d)
                    try:
                        dummy=t.getcoldesc(colname)
                    except RuntimeError:
                        dummy=None
                    if dummy is not None:
                        desc=t.getcoldesc(colname)
                        newname=colname+'_SCALED'
                        desc['name']=newname
                        t.addcols(desc)
                        d=t.getcol(colname)
                        d*=factor
                        t.putcol(newname,d)

                    t.close()
    if os.path.isfile('image_bootstrap.app.mean.fits'):
        warn('Mean bootstrap image exists, not creating it')
    else:
        warn('Creating mean bootstrap image')
        hdus=[]
        for obsid in Uobsid:
            hdus.append(fits.open('image_bootstrap_'+obsid+'.app.restored.fits'))
        for i in range(1,len(Uobsid)):
            hdus[0][0].data+=hdus[i][0].data
        hdus[0][0].data/=len(Uobsid)
        hdus[0].writeto('image_bootstrap.app.mean.fits')
                
def driver():
    if not LOCAL_DEV:
        from parset import option_list
        from options import options
    else:
        from utils.parset import option_list
        from utils.options import options
    global o
    o=options(sys.argv[1:],option_list)
    SetMS=mpi_manager.MSSet(o['mslist'])
    MainNode=None
    if mpi_manager.RANK==0: MainNode=socket.gethostname()
    if mpi_manager.USE_MPI:
        MainNode=mpi_manager.MPI.COMM_WORLD.bcast(MainNode, root=0)
    MPI_Manager=mpi_manager.mpi_manager(o, SetMS, None, MainNode)
    if MPI_Manager.UseMPI:
        run_bootstrap_mpi(MPI_Manager, o)
    else:
        run_bootstrap(o)

if __name__=='__main__':
    driver()
