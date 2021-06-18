#!usr/bin/env python

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import numpy as np
import os, pdb, sys

###################################################
### Changable parameters
################################################### 

do_del    = 1  #delete old image files
 
do_splat  = 1  #splat out rings
do_fit    = 1  #do mvrc fitting and applying
do_slfcl  = 0  #do 1st order phase self-cal

do_image1 = 1  #rPR
do_image2 = 1  #MVRC
do_image3 = 0  #selfcal
do_jmfit  = 1

EXPCLASS    =   "UVDATA"
EXPSEQ      =   1
EXPDISK     =   1


orbits_of_interest = [0,1,2]
number_of_rings    = 3
array              = [1,2,3,4]
refant             = {0:2,1:2,2:2}
middle_sources     = range(number_of_rings)
middle_sources[0]  = 'G0634-2335'
middle_sources[1]  = 'G1336-0829'
middle_sources[2]  = 'G1901-2112'
orbit_sources      = range(number_of_rings)
orbit_sources[0]   = ['J0636-2113',
                      'J0643-2451',
                      'J0620-2515',
                      'J0639-2141',
                      'J0632-2614',
                      'J0629-1959']
orbit_sources[1]   = ['J1354-0206',
                      'J1351-1449',
                      'J1312-0424',
                      'J1406-0848',
                      'J1305-1033',
                      'J1406-0707']
orbit_sources[2]   = ['J1916-1519',
                      'J1848-2718',
                      'J1928-2035',
                      'J1832-2039',
                      'J1916-2708']


###################################################
###################################################

def main():
    #pdb.set_trace()
    AIPS.userno, EXPERIMENT=get_experiment()
    if not 'EXPERIMENT' in globals(): globals().update({'EXPERIMENT': EXPERIMENT})
    indata = AIPSUVData(EXPERIMENT,EXPCLASS,EXPDISK,EXPSEQ)
    if not indata.exists(): sys.exit(str(indata)+ ' does not exist!')
    else: print str(indata)+ ' exists!'
    for i in orbits_of_interest:
        sources = [middle_sources[i]]+orbit_sources[i]
        if do_del==1:
            for k in range(len(sources)):
                _zapimage(sources[k],1)
                _zapimage(sources[k],2)
                _zapimage(sources[k],3)
                _zapimage(sources[k],4)
                _zapimage(sources[k],5)
                _zapimage(sources[k],6)
                _zapimage(sources[k],7)
                _zapimage(sources[k],8)
                _zapimage(sources[k],9)
                _zapimage(sources[k],10)
        if do_splat==1:
            _splat1(sources)
            _splat2(sources)
            _fring_RATE(middle_sources[i],refant[i])
            _clcal(middle_sources[i],1,1)
        else:
            print 'Not running SPLAT'
        if do_fit == 1:
            _fring_PHAS(middle_sources[i],refant[i])
            _tbout(middle_sources[i],'SN',1)
            _tbout(middle_sources[i],'SN',2)
            _suout(middle_sources[i])
            command = './../fit/fit_phase_plane_v6.py ./phases/SN2_'+middle_sources[i][:5]+'.dat'
            os.system(command)
            _tbin(middle_sources[i],'SN',3)
            _snsmo_PHAS(middle_sources[i],3,4)
            _clcal_SELF(middle_sources[i],4,2)
        if do_slfcl == 1:
            _fring_RATE2(middle_sources[i],refant[i])
            _clcal_SELF(middle_sources[i],5,2)
    if do_image1==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],2)
                _zapbeam(orbit_sources[i][j])
            _image(middle_sources[i],middle_sources[i],2)
            _zapbeam(middle_sources[i])
    if do_image2==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],3)
                _zapbeam(orbit_sources[i][j])
            _image(middle_sources[i],middle_sources[i],3)
            _zapbeam(middle_sources[i])
    if do_image3==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],4)
                _zapbeam(orbit_sources[i][j])
            _image(middle_sources[i],middle_sources[i],4)
            _zapbeam(middle_sources[i])
    if do_jmfit==1:
        for i in orbits_of_interest:        
            sources = [middle_sources[i]]+orbit_sources[i]
            for inseq in range(3):
                fitout = './phases/'+middle_sources[i]+'_'+str(inseq+1)+'.jmfit'
                if os.path.exists(fitout): 
                    os.remove(fitout)
                for source in sources:
                    _jmfit(inname=source,inseq=inseq+1,fitout=fitout)

###################################################
### AIPS TASKS
###################################################


def _jmfit(inname='',inseq=1,fitout=''):
    jmfit           =  AIPSTask('jmfit')
    jmfit.default
    image           =  AIPSImage(inname,'ICL001',EXPDISK,inseq)
    jmfit.indata    =  image
    jmfit.fitout    =  fitout
    jmfit.doprint   = -4
    jmfit.blc[1:]  = [216.0, 216.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.trc[1:]  = [296.0, 296.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.niter     =  500
    if image.exists(): jmfit()

def check_sncl(indata,sn,cl):
    if (indata.table_highver('AIPS CL')==cl and
        indata.table_highver('AIPS SN')==sn):
            print('Tables are fine')
    
    if indata.table_highver('AIPS CL')<cl:
        print(indata.table_highver('AIPS CL'),cl)
        raise RuntimeError('Not enough CL tables')

    if indata.table_highver('AIPS CL')>cl:
        while indata.table_highver('AIPS CL')>cl:
            indata.zap_table('AIPS CL', 0)

    if indata.table_highver('AIPS SN')<sn:
        print(indata.table_highver('AIPS SN'),sn)
        raise RuntimeError('Not enough SN tables')

    if indata.table_highver('AIPS SN')>sn:
        while indata.table_highver('AIPS SN')>sn:
            indata.zap_table('AIPS SN', 0)

def _zapimage(source,inseq):
    img=AIPSImage(source,'ICL001',1,inseq)
    if img.exists():
        img.zap()

def _zapbeam(source):
    beam=AIPSImage(source,'IBM001',1,1)
    if beam.exists():
        beam.zap()

def _image(midsource,source,Gin):
    imagr=AIPSTask('imagr')
    imagr.default
    imagr.inname=midsource
    imagr.source[1]=source
    imagr.inclass='RING'
    imagr.inseq=2
    imagr.indisk=1
    imagr.docalib=1
    imagr.gainu=Gin
    imagr.outseq=0
    #imagr.bmaj=1.6e-3
    #imagr.bmin=1.6e-3
    #imagr.nbox=1
    #imagr.clbox[1][1:]=[108,108,148,148]
    imagr.antenna[1:]= array
    imagr.baseline[1:]=array
    imagr.outname=source
    imagr.cell[1:]=[0.5e-4,0.5e-4]
    imagr.imsize[1:]=[512,512]
    imagr.gain=0.3
    imagr.niter=200
    imagr.dotv=-1
    imagr()

def _snsmo_PHAS(midsource,inver,outver):
    snsmo=AIPSTask('snsmo')
    snsmo.default
    data = AIPSUVData(midsource,'RING',EXPDISK,2)
    snsmo.indata=data
    snsmo.bparm[2]=1.0
    snsmo.cparm[2]=1.0
    snsmo.cparm[7]=57.0
    snsmo.smotype='PHAS'
    snsmo.samptype='MWF'
    snsmo.invers=inver
    snsmo.outvers=outver
    snsmo.doblank=-1
    snsmo.dobtween=-1
    snsmo()

def _splat1(sources):
    splat=AIPSTask('splat')
    splat.default
    splat.inname=EXPERIMENT
    splat.inclass=EXPCLASS
    splat.inseq=EXPSEQ
    splat.indisk=EXPDISK
    splat.outseq=1
    splat.outdisk=1
    splat.sources[1:]=sources
    splat.outname=splat.sources[1]
    splat.docalib=1
    splat.gainu=0
    splat.aparm[1]=0
    splat.bchan=0
    splat.echan=0
    splat.outclass='RING'
    data = AIPSUVData(splat.outname,splat.outclass,1,1)
    if data.exists():
        data.zap()
    splat()

def _splat2(sources):
    splat=AIPSTask('splat')
    splat.default
    splat.inname=sources[0]
    splat.inclass='RING'
    splat.inseq=1
    splat.indisk=EXPDISK
    splat.outseq=2
    splat.outdisk=1
    splat.sources[1:]=sources
    splat.outname=splat.sources[1]
    splat.docalib=1
    splat.gainu=0
    splat.aparm[1]=2
    splat.bchan=4
    splat.echan=28
    splat.outclass='RING'
    data = AIPSUVData(splat.outname,splat.outclass,1,2)
    if data.exists():
        data.zap()
    splat()

def _fring_MDEL(midsource,refant):
    fring=AIPSTask('fring')
    fring.default
    data = AIPSUVData(midsource,'RING',1,1)
    fring.indata = data
    fring.calsour[1]=midsource
    fring.snver=1
    fring.solint=-1
    fring.refant=refant
    fring.aparm[1:]=[0,0,0,0,2,0,0,0,0,0]
    fring.dparm[1:]=[0,10,10,0,0,0,0,2,0]
    #fring()

def _fring_RATE(midsource,refant):
    fring=AIPSTask('fring')
    fring.default
    fring.inname=midsource
    fring.calsour[1]=midsource
    fring.inclass='RING'
    fring.inseq=2
    fring.indisk=EXPDISK
    fring.snver=1
    fring.solint=2
    fring.refant=refant
    fring.aparm[1:]=[0,0,0,0,0,0,0,0,0,0]
    fring.dparm[1:]=[2,-1,20,0,0,0,0,0,0]
    fring()

def _fring_RATE2(midsource,refant):
    fring=AIPSTask('fring')
    fring.default
    data                    =  AIPSUVData(midsource,'RING',1,2)
    fring.indata            =  data
    check_sncl(data,4,3)
    fring.calsour[1]        = ''
    fring.docalib           =  1
    fring.gainu             =  2
    fring.snver             =  5
    fring.solint            =  2
    fring.refant            =  refant
    fring.aparm[1:]         = [0,0,0,0,0,0,1,0,0,0]
    fring.dparm[1:]         = [2,-1,20,0,0,0,0,0,0]
    fring()

def _fring_PHAS(midsource,refant):
    fring=AIPSTask('fring')
    fring.default
    fring.inname    =  midsource
    fring.inclass   = 'RING'
    fring.inseq     =  2
    fring.indisk    =  EXPDISK
    fring.snver     =  2
    fring.solint    =  2
    fring.refant    =  refant
    fring.docalib   =  1
    fring.gainu     =  0
    fring.aparm[1:] = [3,0,0,0,0,0,4,0,0,0]
    fring.dparm[1:] = [3,0,0,0,0,0,0,0,0]
    fring()

def _tbout(midsource,table,version):
    tbout=AIPSTask('tbout')
    tbout.default
    tbout.inname  = midsource
    tbout.inclass = 'RING'
    tbout.inseq   = 2
    tbout.indisk  = EXPDISK
    tbout.inext   = table
    tbout.inver   = version
    if not os.path.exists('./phases'): os.mkdir('./phases')
    outpath       = './phases/' + table+str(version)+'_'+midsource[:5]+'.dat'
    tbout.outtext = outpath
    tbout.docrt   = 1000
    if os.path.exists(outpath): os.remove(outpath)
    tbout.go()

def _tbin(midsource,table,ver):
    tbin=AIPSTask('tbin')
    tbin.default
    data = AIPSUVData(midsource,'RING',EXPDISK,2)
    tbin.outdata = data
    check_sncl(data,ver-1,2)
    inpath       = './phases/' + table+str(ver)+'_'+midsource[:5]+'.dat'
    tbin.intext = inpath
    tbin.go()

def _suout(midsource):
    table   = 'SU'
    version =  1 
    prtab   =  AIPSTask('prtab')
    prtab.default
    prtab.inname  = midsource
    prtab.inclass = 'RING'
    prtab.inseq   = 2
    prtab.indisk  = EXPDISK
    prtab.inext   = table
    prtab.inver   = version
    prtab.box[1][1:] = [1,2,11,12]
    if not os.path.exists('./phases'): os.mkdir('./phases')
    outpath       = './phases/' + table+str(version)+'_'+midsource[:5]+'.dat'
    prtab.outprint= outpath
    prtab.docrt   = -1
    if os.path.exists(outpath): os.remove(outpath)
    prtab.go()

def _clcal(midsource,snver,clver):
    clcal=AIPSTask('clcal')
    clcal.default
    clcal.inname=midsource
    clcal.indisk=1
    clcal.interpol='AMBG'
    clcal.inclass='RING'
    clcal.inseq=2
    clcal.snver=snver
    clcal.gainv=clver
    clcal.gainu=clver+1
    clcal()

def _clcal_SELF(midsource,snver,clver):
    clcal=AIPSTask('clcal')
    clcal.default
    clcal.inname=midsource
    clcal.indisk=1
    clcal.inclass='RING'
    clcal.inseq=2
    clcal.interpol='SELF'
    clcal.snver=snver
    clcal.gainv=clver
    clcal.gainu=0
    clcal()

def _calib_AnP(midsource,refant,ring_num,Gin,Sout):
    calib=AIPSTask('calib')
    calib.default
    calib.inname    = midsource
    calib.calsour[1]= '-'+midsource
    calib.inclass   = 'RING'
    calib.inseq     = ring_num+1
    calib.indisk    = EXPDISK
    calib.gainu     = Gin
    calib.docalib   = 1
    calib.snver     = Sout
    calib.refant    = refant
    calib.solint    = 2
    calib.solmode   = 'A&P'
    calib.soltype   = 'L1R'
    calib.aparm[1:] = [0,0,0,0,0,0,4,0,0,0]
    calib()

###################################################
###################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=map(lambda s: s.strip(), g)
    return g

def get_experiment():
    files = os.listdir('./')
    reffile=0
    for i in range(len(files)):
        if "LST" in files[i]:
            reffile = files[i]
            break
    if reffile==0:
        sys.exit('No AIPS output file to read from')
    readfile = get_file(reffile)
    aipsid   = readfile[0].split()[2]
    experiment = readfile[1].split()[2]
    return int(aipsid), experiment

###################################################
###################################################

if __name__=='__main__':
    main()


###################################################
###################################################
