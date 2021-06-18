#!/usr/local/bin/python3

#version 2, removing loops where possible to speed up
#also remove weights and gains

import sys, os, random, time, warnings
import numpy as np
from math import exp, log, pi
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from matplotlib import rc, font_manager
from scipy.special import jv
from difflib import get_close_matches as match

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text',usetex=True)

warnings.filterwarnings('ignore')

######################################################################

def main():
    ##################################################################
    datain  = sys.argv[1] #    './prtuv_ab.log'
    ##################################################################
    start          = time.time()
    dataout1       = './fit2/'+datain[:-4]+'_lsqhalohalo.out'
    dataout2       = './fit2/'+datain[:-4]+'_lsqhalhalo2.out'
    catalogue      = loaduvdata(datain)
    source_list0   = make_vector(catalogue,0)
    data0          = make_vector(catalogue,1)
    mmc            = np.asarray(splitt(get_file('./apjsab06fbt6_mrt.txt'))[33:])
    ref            = np.asarray(splitt(get_file('./v534.su')))
    correct_names0 = closest_match_coords(list(map(lambda x: x.split('_')[0],source_list0)),ref,mmc)

    # combining double observed maser spots using their velocity
    data    = []; correct_names = []; source_list = []
    V       = np.array(list(map(lambda x: x.split('_')[1],source_list0))).astype('float')
    cn      = np.array(correct_names0)
    indx    = np.array(range(len(correct_names0)))
    masers  = remove_repeats(correct_names0)
    for maser in masers:
        M   = (np.matrix(V[indx[cn==maser]]).T*np.matrix(np.ones(shape=V[indx[cn==maser]].shape))).T
        x,y = np.asarray(np.where(abs(M - M.T+0.2*np.identity(V[indx[cn==maser]].shape[0]))<=0.135))
        X   = indx[cn==maser][x[:int(len(x)/2)]]; Y = indx[cn==maser][y[:int(len(y)/2)]]
        for i in indx[cn==maser]:
            if i in X: 
                data          += [[data0[X[X==i][0]][0]+data0[Y[X==i][0]][0],data0[X[X==i][0]][1]+data0[Y[X==i][0]][1]]]
                correct_names += [correct_names0[i]]
                source_list   += [source_list0[i]]
            elif i in Y: pass
            else: 
                data          += [[data0[i][0],data0[i][1]]]
                correct_names += [correct_names0[i]]
                source_list   += [source_list0[i]]
    
    ###################################
    #### create out_file for parms ####
    ###################################
    out_file1=open(dataout1, 'w+')
    print("{:<14s}{:>10s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}".format('Source', 'Velocity',
        'S_h1','S_h2','th_h1','th_h2','chisq','df','err'),file=out_file1)
    out_file1.close(); out_file1=open(dataout1, 'a+')
    print("{:<14s}{:>10s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}".format('','(km/s)','(Jy)',
        '(Jy)','(mas)','(mas)','','','(Jy)'),file=out_file1)
    ###################################
    out_file2=open(dataout2, 'w+')
    print("{:<14s}{:>10s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}".format('Source', 'Velocity',
        'S_h1','S_h2','th_h1','th_h2','chisq','df','err'),file=out_file2)
    out_file2.close(); out_file2=open(dataout2, 'a+')
    print("{:<14s}{:>10s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}{:>12s}".format('','(km/s)','(Jy)',
        '(Jy)','(mas)','(mas)','','','(Jy)'),file=out_file2)
    print('Setup for {:<10.0f} targets in {:10.2f} seconds'.format(len(source_list),time.time()-start))    
    ######################################################################
    ######### least_squares fitting and plotting for 3 models ############
    ######################################################################
    X = np.arange(4e-7,101e3,100)
    n = 0
    for sn in range(len(data)):
        src_data = np.asarray(data[sn])
        D        = src_data[:,src_data[0,:]!=0]
        D2       = src_data[1,src_data[0,:]==0]
        err0     = D2.std(0)
        d_mean   = src_data[1,src_data[0,:]==0].mean(0)
        if d_mean<0.0: continue
        d        = np.append(D,[[0.0],[d_mean]],axis=1)
        d2       = src_data
        ######################################################################
        ######################## fitting and plotting ########################
        ######################################################################
        name, vel = source_list[sn].split('_')
        Name = correct_names[sn]
        res_lsq1 = least_squares(model_core_halo_res1,[1.0,0.3,1.0,0.3],args=(d[0,:],d[1,:]),
            bounds=[[0.3,0.1,0.3,0.1],[np.inf,100,np.inf,300]])
        res_lsq2 = least_squares(model_core_halo_res1,[1.0,0.3,1.0,0.3],args=(d[0,:],d[1,:]),
            bounds=[[0.3,0.1,0.3,0.1],[np.inf,100,np.inf,300]], loss='soft_l1')
        if (res_lsq1.x < [1.0,0.3,1,0.3]).all() or (res_lsq2.x < [1.0,0.3,1,0.3]).all(): continue
        chis2_m10= sum(  ((model_core_halo1(res_lsq1.x,d[0,:])-d[1,:])/err0)**2  )  /  (len(d[0,:])-5) 
        err = err0*(chis2_m10)**0.5
        chis2_m1 = sum(  ((model_core_halo1(res_lsq1.x,d[0,:])-d[1,:])/err)**2  )  /  (len(d[0,:])-5) 
        A_m1 = model_core_halo1(res_lsq1.x,X); A_m2 = model_core_halo1(res_lsq2.x,X)        
        ######################################################################
        S = np.asarray(res_lsq1.x)[np.asarray([0,2])]
        T = np.asarray(res_lsq1.x)[np.asarray([1,3])]
        print("{:<14s}{:>+10.2f}{:>12.1f}{:>12.1f}{:>12.2f}{:>12.2f}{:>12.1f}{:>12.0f}{:>12.2f}".format(Name,
            float(vel),S[T.argmax()],S[T.argmin()],T.max(),T.min(),chis2_m1,len(d[0,:])-5,err),file=out_file1)
        ######################################################################
        S = np.asarray(res_lsq2.x)[np.asarray([0,2])]
        T = np.asarray(res_lsq2.x)[np.asarray([1,3])]
        print("{:<14s}{:>+10.2f}{:>12.1f}{:>12.1f}{:>12.2f}{:>12.2f}{:>12.1f}{:>12.0f}{:>12.2f}".format(Name,
            float(vel),S[T.argmax()],S[T.argmin()],T.max(),T.min(),chis2_m1,len(d[0,:])-5,err),file=out_file2)
        ######################################################################
        plt.close(1); plt.close(2)
        ######################################################################
        fig1 = plt.figure(1,figsize=(5,5))
        ax1  = fig1.add_subplot(111)
        ax1.plot(np.asarray(d[0,:])/1e3,d[1,:],'k.')
        ax1.errorbar(np.asarray(d[0,:])/1e3,d[1,:],yerr=0.5*err,ecolor='k',fmt='none',
            legend=False,alpha=0.25,capsize=2)
        ax1.set_xlabel(R"$uv-$distance $(M\lambda)$")
        ax1.set_ylabel("Correlated Flux Density (Jy)")
        ax1.set_xlim([-1,100])
        ax1.plot(X/1e3,A_m1,'r--')
        ax1.plot(X/1e3,A_m2,'m--')
        ax1.set_ylim([-0.1,1.1*max([max(d[1,:]),max(A_m1),max(A_m2)])])
        ax1.grid(True,alpha=0.5)
        ax1.legend([R'$uv-$data',R'Linear',R'Softl1 Loss']) 
        ax1.hlines(xmin=-10,xmax=50,y=10,color='g',ls='-.',lw=1.0)
        ax1.hlines(xmin=50,xmax=110,y=15,color='g',ls='-.',lw=1.0)
        outname = 'G{:}_{:.2f}'.format(Name,float(vel))
        fig1.savefig('./fit2/{:}.pdf'.format(outname.replace('.',',')))
        ######################################################################

        ######################################################################        
        fig2   = plt.figure(2,figsize=(5,5))
        uplims = np.zeros(d[1,:].shape)
        yerror = 0.5*(err/d[1,:])*1/(2.3025)
        yerror[d[1,:]<=err] = (-np.log10(0.3)+np.log10(err))
        uplims[d[1,:]<=err] = True
        d[1,:][d[1,:]<=err] = err
        ax2 = fig2.add_subplot(111)
        ax2.plot(np.asarray(d[0,:])/1e3,np.log10(d[1,:]),'k.')
        ax2.errorbar(np.asarray(d[0,:])/1e3,np.log10(d[1,:]),yerr=yerror,ecolor='k',fmt='none',
            legend=False,alpha=0.25,capsize=2,uplims=uplims)
        ax2.set_xlabel(R"$uv-$distance $(M\lambda)$")
        ax2.set_ylabel(R"$\log_{10} S_\nu$ (Jy)")
        ax2.set_xlim([-1,100])
        ax2.plot(X/1e3,np.log10(A_m1),'r--')
        #ax2.plot(X/1e3,np.log10(A_m2),'m--')
        ax2.set_yticks(np.log10([0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000]))
        ax2.set_yticklabels(['0.1','0.2','0.5','1','2','5','10','20','50','100','200',
            '500','1000','2000','5000'])
        ax2.set_ylim([np.log10(0.1),np.log10(1.5*max([max(d[1,:]),max(A_m1),max(A_m2)]))])
        ax2.grid(True,alpha=0.5)
        ax2.legend([R'$uv-$data',R'LSQ']) 
        ax2.hlines(xmin=50,xmax=110,y=np.log10(15),color='g',ls='-.',lw=1.0)
        ax2.hlines(xmin=-10,xmax=50,y=np.log10(10),color='g',ls='-.',lw=1.0)
        #ax2.yaxis.tick_right()
        #ax2.yaxis.set_label_position("right")
        #fig2.savefig('./fit2/G{:}_{:.2f}_log10.pdf'.format(Name,float(vel)))
        ######################################################################
        n = n + 2
        ######################################################################
    out_file1.close()
    out_file2.close()
    print('{:<5.0f} fits in {:10.2f} sec'.format(n,time.time()-start))

######################################################################

def model_core_halo1(p, uv_ij):
    factor = 180.0/pi*3600/1*1000/1 #mas
    A_ij = (p[0]*np.exp(-3.56*((1/factor*p[1]*1000)**2*(np.array(uv_ij)**2))) 
                + p[2]*np.exp(-3.56*((1/factor*p[3]*1000)**2*(np.array(uv_ij)**2)))) # Janskys
    return A_ij

def model_core_halo_res1(p, uv_ij, A):
    A_ij = model_core_halo1(p,uv_ij)
    return A_ij-A

######################################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f     = open(fopen, 'r+')
    g     = list(f)
    G     = list(map(lambda s: s.strip(), g))
    return G

def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return new_list

def make_vector(M, col):
    #takes column col from array M (list-list) and turns into 1D vector (list)
    x=[]
    for i in M:
        x+=[i[col]]
    return x

######################################################################

def remove_repeats(list):
    track_count=[]
    for entry in list:
        if entry in track_count:
            continue
        else:
            track_count+=[entry]
    return track_count

def loaduvdata(path):
    #highly specific
    #loads the data from a space-delimited file.
    #assumes format is [source, a1, b1, x1, y1, a2, b2, x2, y2], etc
    #changes to [source, [x], [y]],
    data = splitt(remove_repeats(get_file(path)))
    for row in data:
        if len(row)==0:
            data.remove(row)
    out_data = np.asarray(range(len(data))).tolist()
    for row_number in range(len(data)):
        columns  = len(data[row_number])
        out_data[row_number]=[data[row_number][0],[[],[]]]
        b1=[4*s+1 for s in range(columns//4)] #baseline 1 pointers
        b2=[4*s+2 for s in range(columns//4)] #baseline 2 pointers
        uv=[4*s+3 for s in range(columns//4)] #uv dist pointers
        f =[4*s+4 for s in range(columns//4)] #flux pointers
        for clmm1 in range(columns-1):
            column=clmm1+1 #skip 0 entry
            if column in b1:
                '' #out_data[row_number][1][0]+=[int(data[row_number][column])]
            elif column in b2:
                '' #out_data[row_number][1][1]+=[int(data[row_number][column])]
            elif column in uv:
                out_data[row_number][1][0]   +=[float(data[row_number][column])] #in lambda
            else:
                out_data[row_number][1][1]   +=[float(data[row_number][column])] #in Janskys
    return out_data

def closest_match_coords(targets,ref,mmc):
    sign   = np.array(list(map(lambda x: x.strip(x[1:])+'1.0',mmc[:,5]))).astype(float)
    ra_deg = 15.0*(mmc[:,2].astype('float')+(mmc[:,3].astype('float')+mmc[:,4].astype('float')/60.0)/60.0)
    de_deg = mmc[:,5].astype('float')+sign*((mmc[:,6].astype('float')+mmc[:,7].astype('float')/60.0)/60.0)
    new_names = []
    for name in targets:
        coords = ref[ref[:,0]==name][:,1:].astype('float')
        new_names.append(mmc[(abs(np.array(list(zip(ra_deg, de_deg)))-coords)**2).sum(1).argmin(),1])
    return new_names

######################################################################

if __name__ == '__main__':
    main()

######################################################################
