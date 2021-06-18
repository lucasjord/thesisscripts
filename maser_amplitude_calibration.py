#!usr/local/bin env python

#this script is designed to perform ampliude calibrations
#on maser data when no gain curve is provided. 
#The principle is that all masers have a intrinsic brightness 
#that should be fully captured by all sites simultaneously in 
#THE TOTAL POWER SPECTRA. 
#Equating this over time should be able to identify the changing SEFDs of the telescopes.

#version 3

################################################### IMPORTS
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage

AIPS.userno=534

import scipy.io as sio

apply_calib = 1

import numpy as np
import os
import matplotlib.pyplot as plt
from math import factorial
import copy
import sys
#os.system('sudo sysctl -w kern.tty.ptmx_max=600')

################################################### parameters
exp = 'B'
exp2 = 'b'

base = '/Users/lucash/V534/'

path = base+'V534'+exp+'/SPEC/'
source_path = base+'v534'+exp2+'.vex'

#EXPERIMENT = "V534A_N_60S"
#EXPCLASS = "UVAVG"
EXPERIMENT = "V534"+exp+"_60S"
EXPCLASS = "UVDATA"
EXPSEQ = 1
EXPDISK = 1

ant_name  ={1: 'AT', 2: 'CD', 3: 'HB', 4: 'MP', 5: 'PA', 6: 'WA'}
ant_number={'AT': 1, 'CD': 2, 'HB': 3, 'MP': 4, 'PA': 5, 'WA': 6}
antenna=[1, 2, 3, 4, 5, 6]
polar=['RR', 'LL']

################################################### basic data manipul

def mean(data):
    average=sum(data)/len(data)
    return average

def std(data):
    g=[]
    for i in data:
        g+=[(i-mean(data))**2]
    g2=(sum(g)/(len(data)-1))**0.5
    return g2

def stderr(data):
    d=std(data)/((len(data))**0.5)
    return d

def mean2(data, rng):
    f=[]
    for i in data:
        if i>(mean(data)-rng*std(data)):
            if i<(mean(data)+rng*std(data)):
                f+=[i]
    mf=mean(f)
    return mf

def std2(data, rng):
    f=[]
    for i in data:
        if i>(mean(data)-rng*std(data)):
            if i<(mean(data)+rng*std(data)):
                f+=[i]
    mf=std(f)
    return mf

def stderr2(data, rng):
    f=[]
    for i in data:
        if i>(mean(data)-rng*std(data)):
            if i<(mean(data)+rng*std(data)):
                f+=[i]
    mf=std(f)
    return mf/((len(f))**0.5) 

def w_mean(data, weights):
    #weighted mean
    if not len(data)==len(weights):
        print 'NOT THE SAME LENGTH'
        return
    tmp=[]
    for i in range(len(data)):
        tmp+=[data[i]*weights[i]]
    w_av=sum(tmp)/sum(weights)
    return w_av

def mode(list):
    m=max(set(list), key=list.count)
    if list.count(m)==1:
        raise TypeError('No mode in list')
        return
    else:
        return m

################################################### data manipulations

def get_file(fopen):
    #opens and external file and makes it into a list
    n = path
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

################################################### BASIC LIST MANIPULATIONS

def printall(list):
    #prints all items in a list (matlab style)
    for i in list:
        print i
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

def make_int(list):
    for i in range(len(list)):
        try:
            list[i]=int(list[i])
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
        d+=[-vector[i]+vector[i+1]]
    return d

################################################### vector stuff

def avg_pos_vec(vector):
    #removes an element from the vector and moves the remaing elements to the average position
    #designed to be used for x-axis in derivatives. Need to be float type data
    d=[]
    for i in range(len(vector)-1):
        d+=[(vector[i]+vector[i+1])/2]
    return d    

def deriv(vector1, vector2):
    #takes the derivative of vector 1 wrt vector2
    if not len(vector2)==len(vector1):
        print "\n\n\nWARNING!\n\n\nVECTORS NOT THE SAME LENGTH\n\n\n"
    #we'll choose to 'dump' the last entry
    df=diff(vector1)
    dx=diff(vector2)
    dfdx=[]
    for i in range(len(df)):
         dfdx+=[df[i]/dx[i]]
    new_X=avg_pos_vec(vector2)
    return [new_X, dfdx]

def get_section(vector, up, down):
    #removes section out of vector, assuming the vector is sorted into decreasing order
    tmp=[]
    section=[]
    for i in vector:
        if down<i<up:
            section+=[i]
        #else:
        #    section+=[None]
    return section

def remove_section(vector, up, down):
    #removes section out of vector, assuming the vector is sorted into decreasing order
    tmp=[]
    section=[]
    for i in vector:
        if down<i<up:
            ''
        #    section+=[None]
        else:
            section+=[i]
    return section

################################################### SMOOTHING FUNCTION

def sav_g(y, window_size, order, deriv=0, rate=1):
    y=np.array(y)
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

################################################### LINE FITTING

def get_noise(vector, peak):
   length_scale=len(vector)/10
   peak_index=vector.index(peak)
   try:
       range_for_removal=peak_index-length_scale
       section=vector[range_for_removal-length_scale/10:range_for_removal+length_scale/10]
   except ValueError:
       range_for_removal=peak_index+length_scale
       section=vector[range_for_removal-length_scale/10:range_for_removal+length_scale/10]
   section1=sav_g(section, length_scale/10, 3)
   section2=section-section1
   noise=std(section2)
   return noise

def remove_line(Xvec, Yvec, vel):
    #removes a maserline and fits bandpass
    #masers should be less than 20km/s wide
    up_vel=vel+30
    down_vel=vel-30
    up_line_vel=vel+10
    down_line_vel=vel-10
    cut_bpp_line_x=get_section(Xvec, up_vel, down_vel)
    cut_bpp_line_y=Yvec[Xvec.index(max(cut_bpp_line_x)):Xvec.index(min(cut_bpp_line_x))+1]
    remove_line_x=remove_section(cut_bpp_line_x, up_line_vel, down_line_vel)
    remove_line_y=[]
    for i in remove_line_x:
        remove_line_y+=[cut_bpp_line_y[cut_bpp_line_x.index(i)]]
    return [cut_bpp_line_x, cut_bpp_line_y, remove_line_x, remove_line_y]

def fit_line(Xvec, Yvec, vel):
    [x,y,x1,y1]=remove_line(Xvec, Yvec, vel)
    p=np.polyfit(x1,y1,2)
    fit_bp=[]
    fit_noise=[]
    for i in x:
        fit_bp+=[p[0]*i**2+p[1]*i+p[2]]
    for i in x1:
        fit_noise+=[p[0]*i**2+p[1]*i+p[2]]
    line=list(np.array(y)-np.array(fit_bp))
    peak=max(line)
    velocity=x[line.index(peak)]
    noise=std(list(np.array(y1)-np.array(fit_noise)))
    return [peak, velocity, noise]

###################################################

def plot_sefd(antenna):
    ant=inte[antenna]
    rr=ant['RR']
    ll=ant['LL']
    x1=range(len(rr[0]))
    for i in range(len(rr[0])):
        tmp=rr[0]
        x1[i]=mean(tmp[i])
    x2=range(len(ll[0]))
    for i in range(len(ll[0])):
        tmp=ll[0]
        x2[i]=mean(tmp[i])
    y1=rr[1]
    y2=ll[1]
    print len(y1)
    print len(y2)
    err1=list(np.array(y1)/(np.array(rr[2])))
    err2=list(np.array(y2)/(np.array(ll[2])))
    plt.errorbar(x1, y1, xerr=0, yerr=err1, fmt='r.')
    plt.errorbar(x2, y2, xerr=0, yerr=err2, fmt='b.')
    plt.show()
    return

def get_ant(antenna):
    ant=inte[antenna]
    rr=ant['RR']
    ll=ant['LL']
    x1=rr[0]
    x2=ll[0]
    y1=rr[1]
    y2=ll[1]
    return [x1, y1, x2, y2]

def average_same_entries(time,y):
    T=[]
    Y=[]
    k=range(len(time))
    try:
        for i in k:
            if not i==0:
                if time[i]==time[i-1]:
                    continue #skip already processed entries
            if time[i]==time[i+1]:
                for j in range(10):
                    if not time[i]==time[i+j]: #check for another entry
                        break #break when entries not the same
                T+=[time[i]]
                yentry=mean(y[i:i+j])
                Y+=[yentry]
                #print 'Averaged'
            else:
                T+=[time[i]]
                Y+=[y[i]]
                #print 'not the same'
    except IndexError:
        if i==len(time)-1:
            T+=[time[i]]
            Y+=[y[i]]
            print 'End'
    return [T, Y]  

def convertdhms(time):
    dd=int(time/(24*60*60))
    time=time-dd*24*60*60
    hh=int(time/(60*60))
    time=time-hh*60*60
    mm=int(time/60)
    ss=int(time-mm*60)
    return [dd, hh, mm, ss]

def amplitude_calibration(inname, inclass, inseq, indisk, antenna, pol, timerange, prm, ingain, outgain):
    clcor=AIPSTask('clcor')
    clcor.default
    clcor.inname=inname
    clcor.inclass=inclass
    clcor.inseq=inseq
    clcor.indisk=indisk
    clcor.opcode='GAIN'
    [dd1, hh1, mm1, ss1]=convertdhms(timerange[0])
    [dd2, hh2, mm2, ss2]=convertdhms(timerange[1])
    modifypol={'RR': 'R', 'LL': 'L'}
    clcor.ant=[None, ant]
    clcor.stokes=modifypol[pol]
    clcor.clcorprm=[None, prm]
    clcor.timer=[None, dd1, hh1, mm1, ss1, dd2, hh2, mm2, ss2]
    clcor.gainver=ingain
    clcor.gainuse=outgain
    clcor.go()

def flag(inname, inclass, inseq, indisk, antenna, pol, timerange, reason):
    uvflg=AIPSTask('uvflg')
    uvflg.default
    uvflg.inname=inname
    uvflg.inclass=inclass
    uvflg.inseq=inseq
    uvflg.indisk=indisk
    uvflg.opcode='FLAG'
    #start=convertdhms(timerange[0])
    #stop=convertdhms(timerange[1])
    uvflg.ant=[None, ant]
    uvflg.stokes=pol
    uvflg.outfgver=1
    uvflg.timer=[None]+timerange
    uvflg.reason = reason
    uvflg.go()


###################################################

def sort_data(data):
    #just for specific case of this script
    new_list=[]
    for i in data:
        new_list+=[[i[1], i[2], [float(i[3])*24*60*60+float(i[4])*60*60+float(i[5])*60+float(i[6]),
        float(i[7])*24*60*60+float(i[8])*60*60+float(i[9])*60+float(i[10])], 
        float(i[11]), float(i[11])/float(i[15]), int(i[17])]]
    tmp=copy.copy(new_list)
    for i in range(len(new_list)):
        if not tmp.count(new_list[i])==1:
            tmp.remove(new_list[i])
    return tmp

def get_section(matrix, key):
    section=[s for s in matrix if key in s]
    return section

def clean_list(list, key):
    #removes bits of list
    new_list=[]
    for i in list:
        new_list+=[i.strip(key)]
    return new_list

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

def make_matrix(vec1, vec2):
    Matrix=range(len(vec1))
    for i in range(len(vec1)):
        Matrix[i]=[vec1[i], vec2[i]]
    return Matrix

def convert_s(timerange):
    i=timerange
    new=[float(i[0])*24*60*60+float(i[1])*60*60+float(i[2])*60+float(i[3]), float(i[4])*24*60*60+float(i[5])*60*60+float(i[6])*60+float(i[7])]
    return new

def remove_every_second_line(list1):
    list2 = []
    for i in range(len(list1)/2):
        list2 += [list1[2*i]]
    return list2

################################################### loading in sources

kopen=path+EXPERIMENT+'_true_flux.txt' #path to outfile 1(?)
kkk = open(kopen, 'w')
print >> kkk, 'File to store fluxes for '+EXPERIMENT
kkk.close()

gopen=path+EXPERIMENT+'_heights.txt' #path to infile 1
heights_open=splitt(get_file(gopen))
heights=refine_fits(heights_open, 18, None, None)

source_file = get_sources(source_path)
source_list=make_vector(source_file, 1)

#making dictionaries to pipe to
ATCA    ={'RR':[[], [], []], 'LL':[[], [], []]}
CED     ={'RR':[[], [], []], 'LL':[[], [], []]}
HOB     ={'RR':[[], [], []], 'LL':[[], [], []]}
MOP     ={'RR':[[], [], []], 'LL':[[], [], []]}
PARK    ={'RR':[[], [], []], 'LL':[[], [], []]}
WARK    ={'RR':[[], [], []], 'LL':[[], [], []]}

inte={'AT': ATCA, 'CD': CED, 'HB': HOB, 'MP': MOP, 'PA': PARK, 'WA': WARK}

#source_list=['G348.550-0.9']
source_count=[]

kkk = open(kopen, 'a+')
for ii in range(len(source_file)):
    tmp=source_file[ii]
    timer=tmp[0]
    source=tmp[1]
    if source in source_count:
        continue
    else:
        source_count+=[source]
    section_original =[s for s in heights if source in s]
    

    section=[]
    for line in section_original:
        if float(line[11])>0:
            section+=[line]

    vel         = make_float(make_vector(section, 13))
    vel_median  = np.percentile(vel, 50)

    vel_ref     = mean(make_float(make_vector(section, 14)))
    flux_ref    = mean(make_float(make_vector(section, 12)))
    dx          = mean(make_float(make_vector(section, 16)))
    chans       = make_int(make_vector(section, 17))

    up_ref_vel = vel_median+2*dx
    down_ref_vel = vel_median-2*dx
    good_data=[]
    for i in range(len(vel)):
        if vel[i]<=up_ref_vel and vel[i]>=down_ref_vel:
            good_data+=[section[i]]

    if len(good_data)<=5:
        section_2 = remove_every_second_line(section_original) #remove dodgy method
        section = []
        for line in section_2:
            if float(line[11])>0:
                section+=[line]

        vel         = make_float(make_vector(section, 13))
        vel_median  = np.percentile(vel, 50)

        vel_ref     = mean(make_float(make_vector(section, 14)))
        flux_ref    = mean(make_float(make_vector(section, 12)))
        dx          = mean(make_float(make_vector(section, 16)))
        chans       = make_int(make_vector(section, 17))
        
        good_data=[]
        try:
            chan_mode = mode(chans)
            vel_mode = vel[0]-dx*(chan_mode-chans[0])

            up_ref_chan = chan_mode + 1
            down_ref_chan = chan_mode - 1

            for i in range(len(vel)):
                if chans[i]<=up_ref_chan and chans[i]>=down_ref_chan:
                    good_data+=[section[i]]
        except TypeError:
            ''
        if len(good_data)<=5:
            print 'ERROR: Problem with '+source
            #printall(section)
            print >> kkk, '%-s %10.1f %8.3f %10.3f %5s %4i' % (source, 0.0, 0.0, 0.0, '!!', len(good_data))
            continue
        else:
            ''
    else:
        ''

    reduced_data = sort_data(good_data)
    fluxes  = make_vector(reduced_data, 3)
    channel = int(np.percentile(make_vector(reduced_data, 5), 50)) #!!!!!!!!!!! not sure about this

    #now to choose an appropriate scale for calibration. We will initially use Ceduna or Hobart
    #but if they are not both present (unlikely) we will use meth multibeam and their velocities/fluxes from within
    reforder=[3,2]
    try:
        del scale_tmp
    except NameError:
        ''
    
    key = 'blank'
    for ant_n in reforder:
        try:
            scale_tmp = np.percentile([s[3] for s in reduced_data if ant_name[ant_n] in s], 50)
            key = ant_name[ant_n]
            channel_final = channel
            vel_final = vel_median
            break
        except (ZeroDivisionError, ValueError):
            ''
    try:
        scale = scale_tmp
    except NameError:
        channel_final = int(channel + (vel_median - vel_ref)/dx)
        scale = flux_ref
        vel_final = vel_ref
        key = 'MM'

    scaled_fluxes =list(np.array(fluxes)/scale)
    data=range(len(reduced_data))
    for i in range(len(reduced_data)):
        data[i]=copy.copy(reduced_data[i])
    for i in range(len(scaled_fluxes)):
        g=data[i]
        g[4]=scaled_fluxes[i]
    
    #outprint scaling factor ('true' autocorrelated flux) 
    print >> kkk, '%-s %10.1f %8i %10.1f %5s %4i' % (source, vel_final, channel_final, scale, key, len(good_data))

    #now to collect all the solution together in a big dictionary
    for ant in antenna:
        nm=ant_name[ant]
        var=inte[nm]
        tel_section=get_section(data, nm)
        for pol in polar:
            pol_section=get_section(tel_section, pol)
            var2=var[pol]
            for i in range(len(pol_section)):
                bit=pol_section[i]
                var2[0]+=[bit[2]]   #time
                var2[1]+=[bit[4]]   #scaled
                var2[2]+=[scale]    #scale


######################################################################## applying corrections just determined

if apply_calib==1:
    uvdata=AIPSUVData(EXPERIMENT, EXPCLASS, EXPDISK, EXPSEQ)
    uvdata.zap_table('FG', -1) #delete all FG

    high_cl = uvdata.table_highver('CL')
    r = high_cl - 2
    delete_cl = list(np.array(range(r+1))+2) #don't want to delete CL1
    for ver in delete_cl:
        uvdata.zap_table('CL', ver)

    time_list=make_vector(source_file, 0)
    for ant in antenna:
        [x1, y1, x2, y2]=get_ant(ant_name[ant]) #extract antenna data out of inte

        M1=range(len(x1))
        for i in range(len(x1)):
            M1[i]=[x1[i], y1[i]]
        M2=range(len(x2))
        for i in range(len(x2)):
            M2[i]=[x2[i], y2[i]]
        M1=sorted(M1)
        M2=sorted(M2)
        matrix={'RR': M1, 'LL': M2} #dictionary for polarisation data

        for pol in polar:
            M=matrix[pol]
            time=make_vector(M, 0)
            y=make_vector(M, 1)
            [T,Y]=average_same_entries(time, y) #makes time entries unique

            #converting timer data from sec to dd:hh:mm:ss
            for i in range(len(T)):
                u=T[i]
                u=convertdhms(u[0])+convertdhms(u[1])
                T[i]=u
            total=make_matrix(T,Y)
        
            #defining and making a antenna/polarisation dependent calibration list with all times
            calibration_file=range(len(source_file))
            for i in range(len(source_file)):
                if time_list[i] in T:
                    k=T.index(time_list[i])
                    calibration_file[i]=[T[k], Y[k]]
                else:
                    calibration_file[i]=[time_list[i], 'blank']
        
            # #here I am trying to replace non calibrated data spaces with the previous calibrated data entry
            # for i in range(len(calibration_file)):
            #     entry=calibration_file[i]
            #     if entry[1]=='blank':
            #         new_entry=calibration_file[i-1]
            #         if new_entry[1]=='blank':
            #             new_entry[1]=1
            #         entry[1]=new_entry[1]

            #flagging data that wasn't calibrated
            for i in range(len(calibration_file)):
                entry=calibration_file[i]
                timer=entry[0]
                if entry[1]=='blank':
                    flag(EXPERIMENT, EXPCLASS, EXPDISK, EXPSEQ, ant, pol, timer, 'AMP')

            #now apply to uvdata via parseltongue
            for entry in calibration_file:
                timer=entry[0]
                if entry[1]=='blank':
                    entry[1]=1000000
                amp=1/entry[1]**0.5
                amp=int(amp*10000)/10000.0 #rounding
                if uvdata.table_highver('CL')==1:
                    in_cl=1
                    out_cl=2
                else:
                    in_cl=2
                    out_cl=2
                timer=convert_s(timer) #convert time back to seconds, amplitude_calibration converts it back to dd:hh:mm:ss
                timer[1]=timer[1]-0.1 #make timerange box unique
                amplitude_calibration(EXPERIMENT, EXPCLASS, EXPDISK, EXPSEQ, ant, pol, timer, amp, in_cl, out_cl)

        print '\nCalibrated '+ant_name[ant]+pol+'\n'
else:
    print '\n\n\Finished without calibration\n\n'



