#!/usr/local/bin/parseltongue

#  Parseltongue script to interact with AIPS
#  For a LBA schedule with masers demoted 'G'
#  Script will go through VEX file, identify maser sources and output 
#  scalar-time/baseline/pol-averaged spectrum, then find peaks above 20 SNR
#  Once found, script will then output spectral data into a text file
#  Pre-calibration is assumed but I am unsure whether it is necessary as it scalar averages and only
#  uses amplitudes.
################################################### IMPORTS
try:
    from AIPS import AIPS
    from AIPSTask import AIPSTask, AIPSList
    from AIPSData import AIPSUVData, AIPSImage
    AIPS.userno = 534  
    import numpy as np
    #import matplotlib.pyplot as plt
    import os, sys, peakutils
except ImportError:
    print "Are you running the correct PYTHON environment? python2 for python, parseltongue for AIPS"

###################################################
################## MAIN FUNCTION ##################
###################################################

def __main__():
    #rm_file('./prtuv.log')
    exp        = 'B'
    exp2       = 'b'
    path       = '/Users/lucash/V534/V534'+exp+'/SPEC/'
    vex_path   = '/Users/lucash/V534/v534'+exp2+'.vex'
    experiment = "V534"+exp+"_60S"
    expclass   = "UVDATA"
    expdisk    = 1
    expseq     = 1
    uvdata = AIPSUVData(experiment,expclass,expdisk,expseq)
    ant_name   = {1: 'AT', 2: 'CD', 3: 'HB', 4: 'MP', 5: 'PA', 6: 'WA'}
    ant_number = {'AT': 1, 'CD': 2, 'HB': 3, 'MP': 4, 'PA': 5, 'WA': 6}
    antenna    = [1, 2, 3, 4, 5, 6]
    polar      = ['RR', 'LL']
    source_list = make_vector(get_sources(vex_path),1)
    for source in source_list:
        try:
            [heights, locations, identifiers] = calculate_peaks(experiment,expclass,expseq,expdisk,source)
        except IOError:
            continue
        for k in range(len(locations)):
            channel = locations[k]
            get_AIPSUVdata(experiment,expclass,expseq,expdisk,source,channel    ,0-1,'./uvprt.tmp')
            xc     = extract_flux_xcorr('./uvprt.tmp')
            ac     = extract_flux_acorr('./uvprt.tmp')
            get_AIPSUVdata(experiment,expclass,expseq,expdisk,source,channel+150,0-1,'./uuvprt.tmp')
            up_ac  = extract_flux_acorr('./uuvprt.tmp')
            get_AIPSUVdata(experiment,expclass,expseq,expdisk,source,channel-150,0-1,'./duvprt.tmp')
            dn_ac  = extract_flux_acorr('./duvprt.tmp')
            if len(dn_ac) == len(ac):
                for i in range(len(dn_ac)):
                    ac[i][3] = ac[i][3]-0.5*(dn_ac[i][3]+up_ac[i][3])
                vis = xc + ac
            else:
                vis = xc
            print_results(source+'_'+str(identifiers[k])[:7], vis, './prtuv.log')
    
    return

###################################################
################## SUBFUNCTIONS ###################
###################################################

################################################### analysis

def calculate_peaks(inname, inclass, indisk, inseq, source, IF):
    #opens aips, gets scalar-averaged spectrum, finds peaks and then outputs them
    path =  './spec.tmp'
    rm_file(path)
    possm_outprint_scal_av_acor(inname, inclass, inseq, indisk, path, source, IF)
    delete_table(inname, inclass, inseq, indisk, 'PL', -1)
    t = splitt(get_file(path))
    if not t[2][1] == source:
        raise IOError
    T_t     = refine_fits(t,7,None,None)[2:]  # to remove header
    T       = T_t[len(T_t)/2+1:]
    S       = make_float(make_vector(T,5))  #Jy
    VEL     = make_float(make_vector(T,4))  #km/s
    CHANS   =   make_int(make_vector(T,0))
    S_m     = median(S)         #getting baseline and/or noise
    S_rms   = std_approx(S)     #getting rms as 1 sigma
    s       = np.array(S)
    s_norm  = (s-S_m)/S_rms     #SNR only
    indices = peakutils.indexes(s_norm, thres = (20.0/s_norm.max()))
    if len(indices)==0:
        peaks   = []
        chans   = []
        vel     = []
    else:
        peaks   = s[indices].tolist()
        chans   = np.array(CHANS)[indices].tolist()
        vel     = np.array(VEL)[indices].tolist()
    return [peaks, chans, vel]

def get_spectrum(uvdata, source, IF):
    #opens aips, gets scalar-averaged spectrum, finds peaks and then outputs them
    path =  './spec.tmp'
    IF   = 2
    flux = range(len(uvdata.antennas))
    for ant in range(len(uvdata.antennas)):
        rm_file(path)
        possm_outprint_scal_av_acor(uvdata.name, uvdata.klass, 
            uvdata.seq,uvdata.disk,path,source,IF,ant+1)
        delete_table(uvdata.name,uvdata.klass,uvdata.seq,uvdata.disk,'PL',-1)
        t = splitt(get_file(path))
        if not t[2][1] == source:
            raise IOError
        T_t         = refine_fits(t,7,None,None)[2:]  # to remove header
        T           = T_t[len(T_t)/2+1:]
        flux[ant]   = make_float(make_vector(T,5))  #Jy
        VEL         = make_float(make_vector(T,4))  #km/s
        CHANS       =   make_int(make_vector(T,0))

def print_results(source,data,outfile):
    out = open(outfile,'a+')
    print >> out, ''
    print >> out, source,
    for i in range(len(data)):
        for j in range(len(data[i])):
            print >> out, data[i][j],
    return

################################################### BASIC FUNCTIONS

def make_new_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        return
    else:
        os.system('rm -r '+path)
        os.makedirs(path)
        return

def rm_file(path):
    if not os.path.exists(path):
        return
    else:
        os.system('rm -r '+path)
        return

def median(data):
    average=np.percentile(data,50)
    return average

def std_approx(data):
    approx_std = 0.5*(np.percentile(data,84)-np.percentile(data,16))
    return approx_std

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

def make_matrix(vec1, vec2):
    Matrix=range(len(vec1))
    for i in range(len(vec1)):
        Matrix[i]=[vec1[i], vec2[i]]

################################################### SPECIFIC DATA MANIPULATIONS

def refine_fits(old_list, length_cat, symbol_cat, column_retrival_number):
    #to refine an imported list after the use of get_file() 
    #to remove header information or filter the list
    # only works for 2-dim type objects
    #symbol_cat is a single character
    tmp_list=[]
    if not length_cat==None:
        for i in old_list:
            if len(i)==length_cat:
                tmp_list+=[i]
    else:
        tmp_list=old_list
    tmp_list_2=[]
    if not symbol_cat==None:
        for i in tmp_list:
            if not column_retrival_number==None:
                column=i[column_retrival_number-1] #searches specific column for match
            else:
                column=i #searches all columns.
            for row_element in column: #searching rows in columns
                if row_element.count(symbol_cat)>0:
                    tmp_list_2+=[column]
    else:
        tmp_list_2=tmp_list
    return tmp_list_2

def split_block(total_block, key):
    #this arbitrarily splits up a large list into smaller lists via a keyword, 
    #in this case most likely the word 'Header'
    #this will allow the printed files to be separated into different times
    ind=[]
    for i in range(len(total_block)):
        cell=total_block[i]
        if cell.count(key)>0:
            ind+=[i]
    n_blocks=len(ind)
    spl_block=[]
    for i in range(n_blocks-1):
        tmp=[]
        tmp=total_block[ind[i]:ind[i+1]]
        spl_block+=[tmp]
    if n_blocks>1:
        tmp=[]
        tmp=total_block[ind[len(ind)-1]:len(total_block)]
        spl_block+=[tmp]
    else:
        spl_block=[total_block]
    return spl_block

def get_sources(vexfile):
    os.system('cd ~')
    os.system('grep source= '+vexfile+' > sources.txt')
    sources=splitt(get_file('sources.txt'))
    t=make_vector(sources, 0)
    source_list=make_vector(sources,2)
    source_list=clean_list(clean_list(source_list, 'source='), ';')
    t=clean_list(clean_list(t, 'start_day='), ';')
    time=[]
    for i in t:
        time+=[[int(i[5:8]), int(i[9:11]), int(i[12:14]), int(i[15:17])]]
    start_day=min(make_vector(time, 0))
    for i in time:
        i[0]=i[0]-start_day
    timer=[]
    for i in range(len(time)):
        down=time[i]  
        if i==len(time)-1:
            up=[99, 23, 59, 59]
        else:
            up=time[i+1]
        timer+=[down+up]
    source_list3=make_matrix(timer, source_list)
    sl=[s for s in source_list3 if 'G' in s[1]]
    return sl

def remove_repeats(list):
    track_count=[]
    for entry in list:
        if entry in track_count:
            continue
        else:
            track_count+=[entry]
    return track_count

def extract_flux_xcorr(pathtofile):
    #same as extract_flux.pl from Simon
    # but better
    uvfile = get_file(pathtofile)
    [acorr, xcorr] = convert_uvfile(uvfile)

    U   = make_float(make_vector(xcorr,3))
    V   = make_float(make_vector(xcorr,4))
    S   = make_float(make_vector(xcorr,6))
    UV  = sqrt(U,V)
    a   = make_int(make_vector(xcorr, 1))
    b   = make_int(make_vector(xcorr, 2))
    retrn = range(len(a))
    for i in range(len(a)):
        retrn[i] = [a[i], b[i], UV[i], S[i]]
    return retrn

def extract_flux_acorr(pathtofile):
    #same as extract_flux.pl from Simon
    # but better
    uvfile = get_file(pathtofile)
    [acorr, xcorr] = convert_uvfile(uvfile)

    U   = make_float(make_vector(acorr,3))
    V   = make_float(make_vector(acorr,4))
    S   = make_float(make_vector(acorr,6))
    UV  = sqrt(U,V)
    a   = make_int(make_vector(acorr, 1))
    b   = make_int(make_vector(acorr, 2))
    retrn = range(len(a))
    for i in range(len(a)):
        retrn[i] = [a[i], b[i], UV[i], S[i]]
    return retrn

def extract_flux_axcorr(pathtofile):
    #same as extract_flux.pl from Simon
    # but better
    uvfile = get_file(pathtofile)
    [acorr, xcorr] = convert_uvfile(uvfile)
    axcorr = acorr + xcorr
    U   = make_float(make_vector(axcorr,3))
    V   = make_float(make_vector(axcorr,4))
    S   = make_float(make_vector(axcorr,6))
    UV  = sqrt(U,V)
    a   = make_int(make_vector(axcorr, 1))
    b   = make_int(make_vector(axcorr, 2))
    retrn = range(len(a))
    for i in range(len(a)):
        retrn[i] = [a[i], b[i], UV[i], S[i]]
    return retrn

def print_results(source,data,outfile):
    out = open(outfile,'a+')
    print >> out, ''
    print >> out, source,
    for i in range(len(data)):
        for j in range(len(data[i])):
            print >> out, data[i][j],
    return

def convert_uvfile(uvfile):
    #converts uvprt outfiles into U, V, R, L xcorr files.
    #taken from uvplt_printcat.py
    #pretty garbage and could be done better but meh
    #Blame Lucas from 2016
    h=[]
    for i in uvfile:
        if i is not '':
            h+=[i]
    #clean out words
    k=[]
    for i in h:
        if i[0] is '0':
            k+=[i]
        else:
            if i[0] is '1':
                k+=[i]
    #split lines into a list
    l=[]
    for i in k:
        l+=[i.split()]
    #get rid of '-'s
    m=[]
    for i in l:
        element=[]
        for j in i:
            element+=j.split('-')
        m+=[element]
    #get rid of empty cells in list
    for i in m:
        while '' in i:
            i.remove('')
    #to separate out overflowing cells
    p=[]
    for row in m:
        if not len(row)==18:
            new_row=[]
            for element_num in range(len(row)):
                if row[element_num].count('.')==2:
                    element=row[element_num]
                    delim_pos=element.index('.')
                    if element_num<6:
                        brek=delim_pos+3 #if w/amp
                    else:
                        brek=delim_pos+5 #if wt/amp
                    new_row+=[element[:brek],element[brek:]]
                else:
                    new_row+=[row[element_num]]
            p+=[new_row]
        else:
            p+=[row]
    #get rid of autocorrelations
    q=[]
    for row in p:
        if not len(row)==18:
            ''
        else:
            q+=[row]
    xcorr=[]
    acorr=[]
    for i in q:
        if i[1]==i[2]:
            acorr+=[i]
        else:
            xcorr+=[i]
    return [acorr, xcorr]

################################################### BASIC LIST MANIPULATIONS

def printall(list):
    #prints all items in a list (matlab style)
    for i in list:
        print i
    return

def printall2(list):
    #prints all items in a list (matlab style)
    for i in list:
        print i,
    return

def make_vector(m, col):
    x=[]
    for i in m:
        x+=[i[col]]
    return x

def make_str(list):
    for i in range(len(list)):
        list[i]=str(list[i])
    return list

def make_int(list):
    for i in range(len(list)):
        list[i]=int(list[i])
    return list

def make_float(list):
    for i in range(len(list)):
        try:
            list[i]=float(list[i])
        except ValueError:
            try:
                list[i]=list[i-1]
            except ValueError:
                list[i]=list[i+1]
    return list

def diff(vector):
    #takes a list of floats and creates a len(vector)-1 list of the differences between them
    d=[]
    for i in range(len(vector)-1):
        d+=[vector[i]-vector[i+1]]
    return d

def make_matrix(vec1, vec2):
    Matrix=range(len(vec1))
    for i in range(len(vec1)):
        Matrix[i]=[vec1[i], vec2[i]]
    return Matrix

def clean_list(list, key):
    #removes bits of list
    new_list=[]
    for i in list:
        new_list+=[i.strip(key)]
    return new_list

def sqrt(vec1, vec2):
    tmp=[]
    for i in range(len(vec1)):
        tmp+=[(vec1[i]**2+vec2[i]**2)**0.5]
    return tmp

################################################### AIPS TASKS

def possm_outprint_scal_av_xcor(inname, inclass, inseq, indisk, path, source, IF):
    rm_file(path)
    possm              =  AIPSTask('POSSM')
    possm.default
    possm.inname       =  inname
    possm.inclass      =  inclass
    possm.inseq        =  inseq
    possm.indisk       =  indisk
    possm.docalib      =  1
    possm.solint       =  0
    possm.timer        =  [None,0]
    possm.aparm        =  [None, 0-1, 0, 0, 0, 0, 0, 0, 0, 0, 0] #scaler average, 
    possm.bparm        =  [None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    possm.source       =  [None, source]
    possm.antenna      =  [None,0] #ant is a number 0 == 'all'
    possm.stokes       =  'I' #pol is 'RR' or 'LL' or 'I'
    possm.bif          =  IF
    possm.eif          =  IF
    possm.nplots       =  0
    possm.codetype     =  'AMP'
    possm.outtext      =  path
    try:
        possm.go()
    except RuntimeError:
        ''

def possm_outprint_scal_av_acor(inname, inclass, inseq, indisk, path, source, IF, ant):
    rm_file(path)
    possm              =  AIPSTask('POSSM')
    possm.default
    possm.inname       =  inname
    possm.inclass      =  inclass
    possm.inseq        =  inseq
    possm.indisk       =  indisk
    possm.docalib      =  1
    possm.solint       =  0
    possm.timer        =  [None,0]
    possm.aparm        =  [None, 0-1, 0, 0, 0, 0, 0, 1, 0, 0, 0] #scaler avg,auto
    possm.bparm        =  [None,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    possm.source       =  [None, source]
    possm.antenna      =  [None, ant] #ant is a number z0 == 'all'
    possm.stokes       =  'I' #pol is 'RR' or 'LL' or 'I'
    possm.nplots       =  0
    possm.codetype     =  'AMP'
    possm.outtext      =  path
    possm.go()


def get_AIPSUVdata(inname, inclass, inseq, indisk, source, chan, ant, path):
    rm_file(path)
    uvprt              =  AIPSTask('UVPRT')
    uvprt.default
    uvprt.inname       =  inname
    uvprt.inclass      =  inclass
    uvprt.indisk       =  inseq
    uvprt.inseq        =  indisk
    uvprt.sources      =  [None, source]
    uvprt.channel      =  chan
    uvprt.antenna      =  [None, ant]
    #uvprt.baseline     =  uvprt.antenna
    uvprt.stokes       =  'HALF'
    uvprt.docalib      =  1
    uvprt.docrt        = -1
    uvprt.gainu        =  0
    uvprt.outprint     =  path
    uvprt.go()
    return

def delete_table(inname, inclass, inseq, indisk, inext, invers):
    uvdata = AIPSUVData(inname, inclass, inseq, indisk)
    uvdata.zap_table(inext, invers)
    return

################################################### 

__main__()