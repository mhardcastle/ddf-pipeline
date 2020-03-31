#!/usr/bin/env python

from __future__ import print_function
from builtins import range
import os
from auxcodes import separator,warn
from surveys_db import SurveysDB
from shutil import copy2
from subprocess import call
from time import sleep
from datetime import datetime
import glob

archive='/disks/paradata/shimwell/LoTSS-DR2/archive/'

def do_rsync(s):
    while True:
        retval=call(s,shell=True)
        if retval!=0:
            print('Non-zero return value',retval)
        if retval!=30:
            break
        else:
            print('Retry on timeout')
    return retval
   

def download_file(dir,name):
    target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl:'+archive
    s='rsync -avz --progress --timeout=10 '+target+dir+'/'+name+' '+dir
    return do_rsync(s)

def link(name,id,lroot,anchor,wdir):
    if os.path.isfile(wdir+id+'/'+name):
        return '<a href="%s%s">%s</a>' % (lroot,name,anchor)
    else:
        return '&mdash;'
    
if __name__=='__main__':
    workdir='/data/lofar/DR2'
    skip_construct=False
    while True:
        print(datetime.now())
        print()

        # make status plot
        separator('Making plot')
        os.system('plot_db_projection.py /home/mjh/lofar-surveys/static/Tier1-dbstatus.png')

        # make JSON file for HIPS
        separator('Make JSON file')
        os.chdir('/home/mjh/lofar-surveys/static/lotss_aladin')
        os.system('python survey_status_mysql.py')

        if not skip_construct:
            # sync mosaics directory
            separator('Mosaic sync')
            os.chdir(workdir+'/mosaics')
            for d in ['RA13h_field','RA0h_field']:
                s="rsync --progress --timeout=10 -avz --exclude 'astroblank-*' --exclude '*.out' --exclude '*.py' --exclude 'background' --exclude '*~' --exclude 'old' --exclude '*.sh' --exclude 'low-mosaic-weights.fits' --exclude 'mosaic.fits' --exclude 'wavelet' --exclude 'model' --exclude 'residual' --exclude 'mosaic.pybdsmmask.fits' --exclude 'mosaic-weights.fits' --exclude 'reproject-*.fits' --exclude 'weight-*.fits' --exclude 'low-reproject-*.fits' --exclude 'low-weight-*.fits' --exclude 'low-mosaic.fits' %s@ssh.strw.leidenuniv.nl:/disks/paradata/shimwell/LoTSS-DR2/mosaics/%s/ ." % (os.environ['DDF_PIPELINE_LEIDENUSER'],d)
                do_rsync(s)
            # now go through all archived and completed fields and make sure they're in the DR2 directory

        with SurveysDB() as sdb:
            sdb.cur.execute('select * from fields left join quality on fields.id=quality.id where status="Archived" or status="Complete" order by ra')
            result=sdb.cur.fetchall()

        print('There are',len(result),'complete datasets')

        if not skip_construct:
            separator('Preparing release directory')
            releasefiles=['image_full_low_stokesV.dirty.fits','image_full_vlow_QU.cube.dirty.corr.fits.fz','image_full_low_QU.cube.dirty.corr.fits.fz','image_full_vlow_QU.cube.dirty.fits.fz','image_full_low_QU.cube.dirty.fits.fz','image_full_low_m.int.restored.fits','image_full_low_m.app.restored.fits','image_full_ampphase_di_m.NS.tessel.reg','image_full_ampphase_di_m.NS_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_shift.app.facetRestored.fits','image_full_ampphase_di_m.NS_Band0_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_Band1_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_Band0_shift.app.facetRestored.fits','image_full_ampphase_di_m.NS_Band1_shift.app.facetRestored.fits','image_full_ampphase_di_m.NS_Band2_shift.app.facetRestored.fits','astromap.fits','DynSpec*.tgz']

            os.chdir(workdir+'/fields')
            for r in result:
                id=r['id']
                print('Doing',id)
                if not os.path.isdir(id):
                    warn('Directory %s does not exist, making it' % id)
                    os.mkdir(id)
                tdir=workdir+'/fields/'+id
                if r['clustername']=='Herts' and r['location']!="":
                    location=r['location']
                    resolved_release=[]
                    for f in releasefiles:
                        if '*' in f:
                            resolved_release+=[os.path.basename(g) for g in glob.glob(location+'/'+f)]
                        else:
                            resolved_release.append(f)                       

                    if location:
                        for f in resolved_release:
                            source=location+'/'+f
                            if not os.path.isfile(tdir+'/'+f) or (os.path.isfile(source)  and os.path.getmtime(source)>os.path.getmtime(tdir+'/'+f)):
                                if os.path.isfile(source):
                                    print('Need to copy',f,'to',tdir)
                                    copy2(source,tdir)
                                else:
                                    warn('Source file %s does not exist' % source)
                else:
                    # get from archive if necessary
                    if r['status']!='Archived' or r['archive_version']<4:
                        continue
                    else:
                        for f in releasefiles:
                            if '*' in f:
                                g=glob.glob(tdir+'/'+f)
                                if len(g)==0:
                                    download_file(id,f)
                            else:
                                if not os.path.isfile(tdir+'/'+f):
                                    print('Need to download',id+'/'+f,'from archive')
                                    download_file(id,f)

        separator('Write web page')

        outfile=open('/home/mjh/lofar-surveys/templates/dr2-mosaics.html','w')
        for r in result:
            if r['dr2']==0:
                continue
            id=r['id']
            if os.path.isdir(workdir+'/mosaics/'+id) and os.path.isfile(workdir+'/mosaics/'+id+'/mosaic-blanked.fits'):
                root='downloads/DR2/mosaics/'+id+'/'
                f=root+'mosaic-blanked.fits'
                rms=root+'mosaic.rms.fits'
                resid=root+'mosaic.resid.fits'
                low=root+'low-mosaic-blanked.fits'
                cat=root+'mosaic.cat.fits'
                outfile.write('<tr><td>%s</td><td>%.3f</td><td>%.3f</td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td></tr>\n' % (id,r['ra'],r['decl'],f,rms,resid,low,cat))
        outfile.close()

        outfile=open('/home/mjh/lofar-surveys/templates/dr2-fields.html','w')

        for r in result:
            id=r['id']
            lroot='downloads/DR2/fields/'+id+'/'
            if os.path.isdir(workdir+'/fields/'+id):
                fint=link('image_full_ampphase_di_m.NS_shift.int.facetRestored.fits',id,lroot,'True',workdir+'/fields/')
                fapp=link('image_full_ampphase_di_m.NS_shift.app.facetRestored.fits',id,lroot,'App',workdir+'/fields/')
                lowint=link('image_full_low_m.int.restored.fits',id,lroot,'True',workdir+'/fields/')
                lowapp=link('image_full_low_m.app.restored.fits',id,lroot,'App',workdir+'/fields/')
                band=[]
                for i in range(3):
                    band.append(link('image_full_ampphase_di_m.NS_Band%i_shift.int.facetRestored.fits' % i,id,lroot,'%i' %i, workdir+'/fields/'))
                stokesv=link('image_full_low_stokesV.dirty.fits',id,lroot,'Download',workdir+'/fields/')
                stokesqu=link('image_full_low_QU.cube.dirty.corr.fits.fz',id,lroot,'Low true',workdir+'/fields/')
                stokesquvlow=link('image_full_vlow_QU.cube.dirty.corr.fits.fz',id,lroot,'Vlow true',workdir+'/fields/')
                stokesqu_app=link('image_full_low_QU.cube.dirty.fits.fz',id,lroot,'Low app',workdir+'/fields/')
                stokesquvlow_app=link('image_full_vlow_QU.cube.dirty.fits.fz',id,lroot,'Vlow app',workdir+'/fields/')
                if r['nvss_scale'] is None:
                    scale='&mdash;'
                else:
                    scale="%.3f" % (5.9124/r['nvss_scale'])
                outfile.write('<tr><td>%s</td><td>%.3f</td><td>%.3f</td><td>%s</td><td>%s</td><td>%s, %s</td><td>%s, %s</td><td>%s, %s, %s</td><td>%s</td><td>%s, %s, %s, %s</td></tr>\n' % (id,r['ra'],r['decl'],r['end_date'],scale,fint,fapp,lowint,lowapp,band[0],band[1],band[2],stokesv,stokesqu,stokesquvlow,stokesqu_app,stokesquvlow_app))

        outfile.close()
        separator('Quality pipeline')
        os.system('queue_quality.py')
        separator('Sleeping')

        sleep(7200)
