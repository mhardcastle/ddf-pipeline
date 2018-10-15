#!/usr/bin/python

import os
from auxcodes import separator,warn
from surveys_db import SurveysDB
from shutil import copy2
from subprocess import call
from time import sleep
import datetime

archive='/disks/paradata/shimwell/LoTSS-DR2/archive/'

def download_file(dir,name):
    target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl:'+archive
    s='rsync -avz --progress --timeout=10 '+target+dir+'/'+name+' '+dir
    while True:
        retval=call(s,shell=True)
        if retval!=0:
            print 'Non-zero return value',retval
        if retval!=30:
            break
        else:
            print 'Retry on timeout'
    return retval

def link(name,id,lroot,anchor,wdir):
    if os.path.isfile(wdir+id+'/'+name):
        return '<a href="%s%s">%s</a>' % (lroot,name,anchor)
    else:
        return '&mdash;'
    

workdir='/data/lofar/DR2'
skip_construct=False
while True:
    print(datetime.datetime.now())
    print

    # make status plot
    separator('Making plot')
    os.system('plot_db_projection.py /home/mjh/lofar-surveys/static/Tier1-dbstatus.png')

    if not skip_construct:
        # sync mosaics directory
        separator('Mosaic sync')
        os.chdir(workdir+'/mosaics')
        os.system("rsync -avz --exclude '*.out' --include 'P*' --include 'mosaic-blanked.fits' --include 'mosaic.cat.fits' --include 'mosaic.resid.fits' --include 'mosaic.rms.fits' --include 'low-mosaic-blanked.fits' --exclude '*' 'hardcastle@ssh.strw.leidenuniv.nl:/disks/paradata/shimwell/LoTSS-DR2/mosaics/*' .")

        # now go through all archived and completed fields and make sure they're in the DR2 directory


    with SurveysDB() as sdb:
        sdb.cur.execute('select * from fields where status="Archived" or status="Complete" order by ra')
        result=sdb.cur.fetchall()

    print 'There are',len(result),'complete datasets'

    if not skip_construct:
        separator('Preparing release directory')
        releasefiles=['image_full_low_stokesV.dirty.fits','image_full_vlow_QU.cube.dirty.corr.fits.fz','image_full_low_QU.cube.dirty.corr.fits.fz','image_full_low_m.int.restored.fits','image_full_low_m.app.restored.fits','image_full_ampphase_di_m.NS_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_shift.app.facetRestored.fits','image_full_ampphase_di_m.NS_Band0_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_Band1_shift.int.facetRestored.fits','image_full_ampphase_di_m.NS_Band2_shift.int.facetRestored.fits','astromap.fits']

        os.chdir(workdir+'/fields')
        for r in result:
            id=r['id']
            if not os.path.isdir(id):
                warn('Directory %s does not exist, making it' % id)
                os.mkdir(id)
            tdir=workdir+'/fields/'+id
            if r['clustername']=='Herts':
                location=r['location']
                for f in releasefiles:
                    if not os.path.isfile(tdir+'/'+f):
                        source=location+'/'+f
                        if os.path.isfile(source):
                            print 'Need to copy',f,'to',tdir
                            copy2(source,tdir)
                        else:
                            warn('Source file %s does not exist' % source)
            else:
                # get from archive if necessary
                if r['status']!='Archived':
                    continue
                else:
                    for f in releasefiles:
                        if not os.path.isfile(tdir+'/'+f):
                            print 'Need to download',id+'/'+f,'from archive'
                            download_file(id,f)

    separator('Write web page')

    outfile=open('/home/mjh/lofar-surveys/templates/dr2.html','w')
    outfile.write('''{% extends "layout.html" %}
{% set active_page = "dr2" %}
{% block content %}
<div class="title">
          <h2>DR2: internal data release</h2>
                </div>
</div>
<div id="page" class="container">
<div class="content">
<p>This page gives preliminary internal access to the second internal
  data release (DR2). You can preview the images through a HIPS viewer and download full-quality DR2 mosaics and catalogues or pipeline products for individual fields. Later you will be able to search for images by co-ordinates.</p>
<h3>HIPS previewer</h3>
<ul>
<li><a href="http://lofar.strw.leidenuniv.nl/LoTSS_DR2_high_hips">High-resolution preview</a> (best viewed in the Aladin app rather than a web browser)</li>
<li><a href="http://lofar.strw.leidenuniv.nl/LoTSS_DR2_low_hips">Low-resolution preview</a></li>
</ul>
<h3>Mosaics</h3><br>
<table id="table_id" class="display">
<thead>
<tr><th>Field name</th><th>RA</th><th>Dec</th><th>Full-res mosaic</th><th>Full-res rms map</th><th>Full-res residual</th><th>Low-res mosaic</th><th>Catalogue</th></tr>
</thead>
<tbody>
''')

    for r in result:
        id=r['id']
        if os.path.isdir(workdir+'/mosaics/'+id):
            root='downloads/DR2/mosaics/'+id+'/'
            f=root+'mosaic-blanked.fits'
            rms=root+'mosaic.rms.fits'
            resid=root+'mosaic.resid.fits'
            low=root+'low-mosaic-blanked.fits'
            cat=root+'mosaic.cat.fits'
            outfile.write('<tr><td>%s</td><td>%.3f</td><td>%.3f</td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td><td><a href=\"%s\">Download</a></td></tr>\n' % (id,r['ra'],r['decl'],f,rms,resid,low,cat))

    outfile.write('''</tbody></table>
    <h3>Fields</h3><br>
    <table id="field_id" class="display">
    <thead>
    <tr><th>Field name</th><th>RA</th><th>Dec</th><th>Completed</th><th>Full-res image</th><th>Low-res image</th><th>Band images</th><th>Stokes V</th><th>Stokes QU</th></tr>
    </thead>
    <tbody>
    ''')
    for r in result:
        id=r['id']
        lroot='downloads/DR2/fields/'+id+'/'
        if os.path.isdir(workdir+'/fields/'+id):
            fint=link('image_full_ampphase_di_m.NS_shift.int.facetRestored.fits',id,lroot,'True',workdir+'/fields/')
            fapp=link('image_full_ampphase_di_m.NS_shift.int.facetRestored.fits',id,lroot,'App',workdir+'/fields/')
            lowint=link('image_full_low_m.int.restored.fits',id,lroot,'True',workdir+'/fields/')
            lowapp=link('image_full_low_m.app.restored.fits',id,lroot,'App',workdir+'/fields/')
            band=[]
            for i in range(3):
                band.append(link('image_full_ampphase_di_m.NS_Band%i_shift.int.facetRestored.fits' % i,id,lroot,'%i' %i, workdir+'/fields/'))
            stokesv=link('image_full_low_stokesV.dirty.fits',id,lroot,'Download',workdir+'/fields/')
            stokesqu=link('image_full_low_QU.cube.dirty.corr.fits.fz',id,lroot,'Low',workdir+'/fields/')
            stokesquvlow=link('image_full_vlow_QU.cube.dirty.corr.fits.fz',id,lroot,'Vlow',workdir+'/fields/')

            outfile.write('<tr><td>%s</td><td>%.3f</td><td>%.3f</td><td>%s</td><td>%s, %s</td><td>%s, %s</td><td>%s, %s, %s</td><td>%s</td><td>%s, %s</td></tr>\n' % (id,r['ra'],r['decl'],r['end_date'],fint,fapp,lowint,lowapp,band[0],band[1],band[2],stokesv,stokesqu,stokesquvlow))

    outfile.write('''<script>
$(document).ready( function () { 
$('#field_id').DataTable();                                                 
} );                                                                            
</script>
''')
    outfile.write('''</tbody></table>
</div>
</div>
{% endblock %}

''')
    outfile.close()
    separator('Sleeping')
    
    sleep(7200)
