from builtins import str
import numpy as np
import pyrap.tables as pt
import  os,sys
import glob
from auxcodes import getpos

def redo_dppp_di(o):
        mslist_name = o['full_mslist']
        mslist=[s.strip() for s in open(mslist_name).readlines()]

        name,ra,dec = getpos(mslist[0])
        Radius = 5.0
        SkymodelPath = 'dppp-skymodel.txt'

        if o['redo_DI_type'] == 'GSM':
                os.system("wget -O "+SkymodelPath+ " \'https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord="+str(ra)+","+str(dec)+"&radius="+str(Radius)+"&unit=deg&deconv=y\' ")
        if o['redo_DI_type'] == 'TGSS':
                os.system("wget -O "+SkymodelPath+ " \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv3.cgi?coord="+str(ra)+","+str(dec)+"&radius="+str(Radius)+"&unit=deg&deconv=y\' ")

        # Create skymodel with makesourcedb
        os.system('makesourcedb in=%s out=%s format="<" outtype=blob'%(SkymodelPath,'dppp-skymodel.sourcedb'))

        for inmsfile in mslist:
                calparset = open('dppp-cal.parset','w')
                calparset.write('msin=%s\n'%inmsfile)
                calparset.write('msin.datacolumn                 = DATA\n')
                calparset.write('msin.baseline                   = CS*&; RS*&; CS*&RS*\n')
                calparset.write('msout                          = .\n')
                calparset.write('msout.datacolumn                = CORRECTED_DATA\n')
                calparset.write('steps                           = [filter,gaincal]\n')
                calparset.write('filter.type                              = filter\n')
                calparset.write('filter.blrange                           = [150, 999999]\n')
                calparset.write('gaincal.type                    = gaincal\n')
                calparset.write('gaincal.maxiter                 = 500\n')
                calparset.write('gaincal.caltype                 = phaseonly\n')
                calparset.write('gaincal.nchan                   = 0\n')
                calparset.write('gaincal.solint                  = 1\n')
                calparset.write('gaincal.usebeammodel            = True\n')
                calparset.write('gaincal.usechannelfreq          = True\n')
                calparset.write('gaincal.beammode                = array_factor\n')
                calparset.write('gaincal.sourcedb                = dppp-skymodel.sourcedb\n')
                calparset.write('gaincal.parmdb                  = %s/instrument_new\n'%inmsfile)
                calparset.close()

                os.system('DP3 dppp-cal.parset')

                applyparset = open('apply-cal.parset','w')
                applyparset.write('msin                            = %s\n'%inmsfile)
                applyparset.write('msout                                = .\n')
                applyparset.write('msin.datacolumn                 = DATA\n')
                applyparset.write('msout.datacolumn                = CORRECTED_DATA\n')
                applyparset.write('msout.writefullresflag          = False\n')
                applyparset.write('steps                           = [applycal]\n')
                applyparset.write('applycal.type                   = applycal\n')
                applyparset.write('applycal.correction             = gain\n')
                applyparset.write('applycal.parmdb                 = %s/instrument_new\n'%inmsfile)
                applyparset.close()
        
                os.system('DP3 apply-cal.parset')

if __name__=='__main__':
    from options import options
    from parset import option_list
    o=options(sys.argv[1:],option_list)
    redo_dppp_di(o)
