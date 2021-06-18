#!/usr/local/bin/ParselTongue

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import numpy as np
import os, pdb, sys

###################################################
### Changable parameters
################################################### 

do_del   = 1
do_splat = 1
do_fit   = 0

do_image1 = 0
do_image2 = 0
do_image3 = 0
do_image4 = 0
do_jmfit  = 0

EXPCLASS    =  "UVDATA"
EXPSEQ      =   1
EXPDISK     =   1


#orbits_of_interest = [0,4]
orbits_of_interest = [6]
number_of_rings    = 9
array              = [1,2,3,4]
refant             = 1
##################################################################################################################
middle_sources     = range(number_of_rings) ; flux = range(number_of_rings)
middle_sources[0]  = 'G0634-2335' ; flux[0]  = 1.0
middle_sources[1]  = 'G0643-2451' ; flux[1]  = 0.2
middle_sources[2]  = 'G0650-1637' ; flux[2]  = 1.1
middle_sources[3]  = 'G1336-0829' ; flux[3]  = 0.6
middle_sources[4]  = 'G1351-1449' ; flux[4]  = 0.6
middle_sources[5]  = 'G1400-1858' ; flux[5]  = 0.4
middle_sources[6]  = 'G1901-2112' ; flux[6]  = 0.2
middle_sources[7]  = 'G1939-1525' ; flux[7]  = 2.0
middle_sources[8]  = 'G2018-0509' ; flux[8]  = 0.1

##################################################################################################################
orbit_sources      = range(number_of_rings)
orbit_sources[0]   = ['J0636-2113','J0643-2451','J0620-2515','J0639-2141','J0632-2614','J0629-1959']
orbit_sources[1]   = ['J0653-1929','J0648-3044','J0629-1959','J0702-2841','J0620-2515','J0706-2311']
orbit_sources[2]   = ['J0702-1015','J0706-2311','J0621-1402','J0721-1530','J0634-2335']
orbit_sources[3]   = ['J1354-0206','J1351-1449','J1312-0424','J1406-0848','J1305-1033','J1406-0707']
orbit_sources[4]   = ['J1405-1440','J1344-1723','J1349-1132','J1357-1744','J1337-1257']
orbit_sources[5]   = ['J1405-1440','J1409-2315','J1344-1723','J1419-1928','J1342-2051','J1351-1449']
orbit_sources[6]   = ['J1916-1519','J1848-2718','J1928-2035','J1832-2039','J1916-2708']
orbit_sources[7]   = ['J1939-1002','J1949-1957','J1954-1123','J1916-1519','J2000-1748','J1928-2035','J2000-1325']
orbit_sources[8]   = ['J2023-0123','J2025-0735','J2015-0137','J2034-0523','J2008-0418']

##################################################################################################################
##################################################################################################################

def main():
    #pdb.set_trace()
    AIPS.userno, EXPERIMENT=get_experiment()
    if not 'EXPERIMENT' in globals(): globals().update({'EXPERIMENT': EXPERIMENT})
    indata = AIPSUVData(EXPERIMENT,EXPCLASS,EXPDISK,EXPSEQ)
    if not indata.exists(): sys.exit(str(indata)+ ' does not exist!')
    else: print str(indata)+ ' exists!'
    

    for i in orbits_of_interest:
        sources = [middle_sources[i]]+orbit_sources[i]
        ######################################################################################################
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

        ######################################################################################################
        if do_splat==1:
            #_splat_mdel(sources)
            _splat(sources)
            _fring_RATE(inname=middle_sources[i],refant=1,CLin=1,SNout=1) #selfcal SN1
            _clcal(inname=middle_sources[i],SNin=1,CLin=1,CLout=2,interpol='SELF') #apply CL1 + SN1 --> CL2
            _fring_RATE(inname=middle_sources[i],calsour=middle_sources[i],refant=1,CLin=1,SNout=2) #fringefit
            _clcal(inname=middle_sources[i],SNin=2,CLin=1,CLout=3,interpol='AMBG') #apply CL1 + SN2 --> CL3 
            _calib_P(inname=middle_sources[i],refant=1,CLin=3,SNout=3) #phasefit SN3
            _clcal(inname=middle_sources[i],SNin=3,CLin=3,CLout=4,interpol='SELF') #apply to self
        else:
            print 'Not running SPLAT'  
        ###################################################################################################### 
        if do_fit == 1:
            _tbout(middle_sources[i],'SN',3)
            _suout(middle_sources[i])
            command = './../fit/fit_phase_plane_v6.py ./phases/SN3_'+middle_sources[i][:5]+'.dat'
            os.system(command)
            _tbin(middle_sources[i],'SN',4)
            _snsmo_PHAS(inname=middle_sources[i],SNin=4,SNout=5)
            #_tacop(inname=middle_sources[i],inext='CL',inver=3,outver=5)
            _clcal(inname=middle_sources[i],SNin=5,CLin=3,CLout=5,interpol='SELF')


    if do_image1==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],2)
            _image(middle_sources[i],middle_sources[i],2)
    
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

    if do_image4==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],5)
                _zapbeam(orbit_sources[i][j])
            _image(middle_sources[i],middle_sources[i],5)
            _zapbeam(middle_sources[i])

    if do_jmfit==1:
        for i in orbits_of_interest:        
            sources = [middle_sources[i]]+orbit_sources[i]
            for inseq in range(4):
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
    jmfit.blc[1:]  = [236.0, 236.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.trc[1:]  = [276.0, 276.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.niter     =  200
    #jmfit.inp()
    jmfit()


def _tacop(inname='',inext='',inver=0,outver=1):
    tacop           =  AIPSTask('tacop')
    tacop.default
    data            =  AIPSUVData(inname,'RING',EXPDISK,1)
    tacop.indata    =  data
    tacop.inver     =  inver
    tacop.inext     =  inext
    tacop.outver    =  outver
    tacop()

def _splat(sources):
    splat               = AIPSTask('splat')
    splat.default
    splat.inname        = EXPERIMENT
    splat.inclass       = EXPCLASS
    splat.inseq         = EXPSEQ
    splat.indisk        = EXPDISK
    splat.outseq        = 1
    splat.outdisk       = 1
    splat.sources[1:]   = sources
    splat.outname       = splat.sources[1]
    splat.docalib       = 1
    splat.gainu         = 0
    splat.aparm[1]      = 2
    splat.bchan         = 4
    splat.echan         = 28
    splat.outclass      = 'SPLAT'
    data                = AIPSUVData(splat.outname,splat.outclass,1,1)
    if data.exists():
        data.zap()
    splat()
    uvavg               = AIPSTask('uvavg')
    uvavg.default
    uvavg.inname        =  splat.outname
    uvavg.inclass       =  splat.outclass
    uvavg.inseq         =  splat.outseq
    uvavg.indisk        =  splat.outdisk
    uvavg.outseq        =  1
    uvavg.outdisk       =  1
    uvavg.sources[1:]   =  sources
    uvavg.outname       =  splat.sources[1]
    uvavg.outclass      = 'RING'
    uvavg.docalib       =  1
    uvavg.gainu         =  0
    uvavg.outclass      = 'RING'
    uvavg.yinc          =  5
    data2               = AIPSUVData(uvavg.outname,uvavg.outclass,1,1)
    if data2.exists():
        data2.zap()
    uvavg()    
    if data.exists():
        data.zap()

def _calib_P(inname='',calsour='',refant=0,CLin=0,SNout=0):
    calib=AIPSTask('calib')
    calib.default
    calib.inname     =  inname
    calib.calsour[1] =  calsour
    calib.inclass    = 'RING'
    calib.inseq      =  1
    calib.indisk     =  EXPDISK
    calib.gainu      =  CLin
    calib.docalib    =  1
    calib.snver      =  SNout
    calib.refant     =  refant
    calib.solint     =  2
    calib.solmode    = 'P'
    calib.soltype    = 'L1R'
    calib.aparm[1:]  = [0,0,0,0,0,0,4,0,0,0]
    calib()

def _fring_RATE(inname='',calsour='',refant=0,CLin=0,SNout=0):
    fring=AIPSTask('fring')
    fring.default
    fring.inname     =  inname
    fring.calsour[1] =  calsour
    fring.inclass    = 'RING'
    fring.inseq      =  1
    fring.indisk     =  EXPDISK
    fring.snver      =  SNout
    fring.solint     =  2
    fring.docalib    =  1
    fring.gainu      =  CLin
    fring.refant     =  refant
    fring.aparm[1:]  = [0,0,0,0,0,0,0,0,0,0]
    fring.dparm[1:]  = [2,-1,20,0,0,0,0,0,0]
    fring()

def _split(inname,source):
    split               = AIPSTask('split')
    split.default
    split.inname        = inname
    split.inclass       = 'RING'
    split.inseq         = 1
    split.indisk        = 1
    split.outseq        = 99
    split.outdisk       = 1
    split.source[1]     = source
    split.docalib       = 1
    split.gainu         = 0
    data                = AIPSUVData(inname,'SPLIT',1,99)
    if data.exists():
        data.zap()
    split()

def _fittp(inname,path='./difmap/'):
    fittp               = AIPSTask('fittp')
    fittp.default
    data                = AIPSUVData(inname,'SPLIT',1,99)
    fittp.indata        = data
    fittp.dataout       = path+inname+'.SELFCAL.fits'
    if os.path.exists(path+inname+'.SELFCAL.fits'):
        os.remove(path+inname+'.SELFCAL.fits')
    fittp.go()
    data.zap()

def _fring_PHAS(inname,refant,CLin,SNout):
    fring=AIPSTask('fring')
    fring.default
    fring.inname    =  midsource
    fring.inclass   = 'RING'
    fring.inseq     =  1
    fring.indisk    =  EXPDISK
    fring.snver     =  SNout
    fring.solint    =  2
    fring.refant    =  refant
    fring.docalib   =  1
    fring.gainu     =  CLin
    fring.aparm[1:] = [3, 0, 0,0,0,0,4,0,0,0]
    fring.dparm[1:] = [3,-1,10,0,0,0,0,0,0]
    fring()

def _calib_AnP(inname,refant,Gin,Sout,Flux):
    calib=AIPSTask('calib')
    calib.default
    calib.inname     =  inname
    calib.calsour[1] =  inname
    calib.inclass    = 'RING'
    calib.inseq      =  1
    calib.indisk     =  EXPDISK
    calib.gainu      =  CLin
    calib.docalib    =  1
    calib.snver      =  SNout
    calib.refant     =  refant
    calib.solint     =  2
    calib.solmode    = 'A&P'
    calib.smodel[1]  =  Flux
    calib.soltype    = 'L1R'
    calib.aparm[1:]  = [0,0,0,0,0,0,4,0,0,0]
    #calib.inp()
    calib()

def _tbout(inname,table,version):
    tbout=AIPSTask('tbout')
    tbout.default
    tbout.inname  = inname
    tbout.inclass = 'RING'
    tbout.inseq   = 1
    tbout.indisk  = EXPDISK
    tbout.inext   = table
    tbout.inver   = version
    if not os.path.exists('./phases'): os.mkdir('./phases')
    outpath       = './phases/' + table+str(version)+'_'+inname[:5]+'.dat'
    tbout.outtext = outpath
    tbout.docrt   = 1000
    if os.path.exists(outpath): os.remove(outpath)
    tbout.go()

def _tbin(inname,table,ver):
    tbin            = AIPSTask('tbin')
    tbin.default
    data            = AIPSUVData(inname,'RING',EXPDISK,1)
    tbin.outdata    = data
    check_sncl(data,ver-1,4)
    inpath          = './phases/' + table+str(ver)+'_'+inname[:5]+'.dat'
    tbin.intext     = inpath
    tbin.go()

def _suout(inname):
    table   = 'SU'
    version =  1 
    prtab   =  AIPSTask('prtab')
    prtab.default
    prtab.inname  = inname
    prtab.inclass = 'RING'
    prtab.inseq   = 1
    prtab.indisk  = EXPDISK
    prtab.inext   = table
    prtab.inver   = version
    prtab.box[1][1:] = [1,2,11,12]
    if not os.path.exists('./phases'): os.mkdir('./phases')
    outpath       = './phases/' + table+str(version)+'_'+inname[:5]+'.dat'
    prtab.outprint= outpath
    prtab.docrt   = -1
    if os.path.exists(outpath): os.remove(outpath)
    prtab.go()

def _clcal(inname='',SNin=0,CLin=0,CLout=0,interpol=''):
    clcal=AIPSTask('clcal')
    clcal.default
    clcal.inname    = inname
    clcal.indisk    = 1
    clcal.inclass   = 'RING'
    clcal.cutoff    = 15.0
    clcal.inseq     = 1
    clcal.interpol  = interpol
    clcal.snver     = SNin
    clcal.gainv     = CLin
    clcal.gainu     = CLout
    clcal()


def _flag_TIME(midsource):
    os.system('grep -C1 '+midsource+' '+EXPERIMENT+'_C.LST | grep -v 0.000 | grep -v 8.1970 | grep -v 8.213 | grep : | grep 0000 | grep MID >> /tmp/flag')

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

def _image(inname,source,Gin):
    imagr=AIPSTask('imagr')
    imagr.default
    imagr.inname=inname
    imagr.source[1]=source
    imagr.inclass='RING'
    imagr.inseq=1
    imagr.indisk=1
    imagr.docalib=1
    imagr.gainu=Gin
    imagr.outseq=0
    imagr.bmaj=1.6e-3
    imagr.bmin=1.6e-3
    #imagr.nbox=1
    #imagr.clbox[1][1:]=[108,108,148,148]
    imagr.antenna[1:]= array
    imagr.baseline[1:]=array
    imagr.outname=source
    imagr.cell[1:]=[1e-4,1e-4]
    imagr.imsize[1:]=[512,512]
    imagr.gain=0.3
    imagr.niter=200
    imagr.dotv=-1
    imagr()
    _zapbeam(source)

def _snsmo_PHAS(inname='',SNin=0,SNout=0):
    snsmo           =  AIPSTask('snsmo')
    snsmo.default
    data            =  AIPSUVData(inname,'RING',EXPDISK,1)
    snsmo.indata    =  data
    snsmo.bparm[2]  =  1.0
    snsmo.cparm[2]  =  1.0
    snsmo.cparm[7]  =  57.0
    snsmo.smotype   = 'PHAS'
    snsmo.samptype  = 'MWF'
    snsmo.invers    =  SNin
    snsmo.outvers   =  SNout
    snsmo.doblank   = -1
    snsmo.dobtween  = -1
    snsmo()

###################################################
###################################################

def _fring_MDEL(midsource,refant,Gin,Sout):
    fring=AIPSTask('fring')
    fring.default
    fring.inname=midsource
    fring.calsour[1]=midsource
    fring.inclass='RING'
    fring.inseq=1
    fring.indisk=EXPDISK
    fring.snver=Sout
    fring.solint=-1
    fring.docalib=1
    fring.gainu=Gin
    fring.refant=refant
    fring.aparm[1:]=[0,0,0,0,2,0,0,0,0,0]
    fring.dparm[1:]=[0,10,10,0,0,0,0,2,0]
    fring()

def _splat_mdel(sources):
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

###################################################
###################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=map(lambda s: s.strip(), g)
    return g

def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return new_list

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
