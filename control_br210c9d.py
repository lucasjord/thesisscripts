##############################################################################
# ParselTongue Script for BeSSeL-Survey calibration                          #
#                                                                            #
# Reads in data, runs indxr, checks antenna tables, writes station_file.inp, #
# control_file.inp and calibrator_file.inp (for continuum data), downloads   #
# TEC-maps and EOPs, and runs TECOR and EOP corrections, runs manual         #
# phase-cal, finge-fit on geo-sources, writes out files for fitting          #
#                                                                            #
# Version 2.3.0 (2016/03/22)                                                 #
#                                                                            #
# Changes:                                                                   #
#                                                                            #
# 1.0.1: Allow shifting of multiple sources                                  #
# 1.0.2: checkatmos added to look for zeros in zenith delay offsets          #
# 1.1.0: Slightly changed structure, full EVN capability                     #
# 2.0.0: Moved main script to definitions file                               #
# 2.0.1: Added DELZN, and imput model for fringe                             #
# 2.1.0: Added plot_tables, new options and experimental self-cal            #
# 2.2.0: Added AIPSVER, niter, get_key_flag, first epoch products            #
# 2.3.0: Added RDBE_check                                                    #
#                                                                            #
##############################################################################

from AIPS import AIPS
import os

def_file    = '/Users/Lucas/bessel/definitions_bessel_current.py'

aipsver     = '31DEC16'

####################
# Input Parameters #
####################

logfile     = 'ParselTongue.log'

file_path   = '/Users/Lucas/bessel/br210c_g037/br210c9d/files/'
eop_path    = '/Users/Lucas/bessel/br210c_g037/br210c9d/eop/'
fit_path    = '/Users/Lucas/bessel/fit/'

AIPS.userno = 379

inter_flag  = 0         # interactive (1) or non-interactive (0) mode
n           = 3         # Number of UV-data files either to download, on disk
                        # or already in AIPS and to be analyzed
defdisk     = 1         # Default AIPS disk to use (can be changed later)

#############################################################################
###                  Do not change or move this part                     ####
[filename, outname, outclass]    = [range(n),range(n),range(n)]
[nfiles, ncount, doconcat]       = [range(n),range(n),range(n)]
[outdisk, flagfile, antabfile]   = [range(n),range(n),range(n)]
for i in range(n):
    [flagfile[i], antabfile[i], outdisk[i]] = ['','',defdisk]
    [nfiles[i], ncount[i], doconcat[i]]     = [0,1,-1]
#############################################################################

antname     = 'VLBA'                 # Antenna order for FITLD    
refant      = 9                      # refant=0 to select refant automatically

calsource   = ''  #'J1855_G037'    # calibrator        '' => automatically
target      = ['']                             # continuum sources '' => automatically
mp_source   = ['']          # fringe finder     '' => automatically
mp_timera   = [0, 0,0,0,0,0,0,0]  # constrain time range for fringe finder?
bandcal     = ['']                   # Bandpass calibrator

delzn_flag  =  0              # Use DELZN (1) or fit_geblocks (0)?
dual_geo    =  0              # Using dual frequency geodetic blocks?

geo_data_nr =  0              # data file with geo data? (<0 for no geoblock)
cont        =  1              # data file with continuum data?
line        =  2              # data file with line data?
                              # if you have only one dataset, use cont = line

channel     = 1002              # channel used for fringe fit
cvelsource  = ['G037.81+0.41']             # line sources '' => calsource
vel         = [19]             # Velocity for cvel

pos_shift   = {'G037.81+0.41'   : [-0.158,0.904]}
                               # Which source should be shifted ''=calsource
                               # e.g. {'G211':[0.2,0.3], 'J0007':[0.1,0.3]}
                               # {'' : [0,0]} => no shift
#################
# Split Options #
#################

smooth       = [0,0,0]          # Smooth during split for line data
split_outcl  = 'SPLIT'          # outclass in SPLIT '' => 'SPLIT'

##################################
# Optional inputs for fringe fit #
##################################

[fr_n, fr_c, fr_d, fr_s] = ['','',1,1]
                                             # Input image to use in FRINGE:
                                             # ['NAME','CLASS',DISK,SEQUENCE]
smodel = [1,0]                 # SMODEL in FRING                                
solint = 0                     # SOLINT in FRING
nmaps  = 1                     # NMAPS in FRING

##############################################
# Optional imput for phased VLA observations #
##############################################

flux = {'' : [0,0,0,0]}         # Flux of cal source (for VLA)

######################
# Imaging Parameters #
######################

cellsize     = 0.0001          # CELLSIZE for IMAGR
imsize       = 1000             # IMSIZE for IMAGR
bchan        = 990            # Start channel for data cube
echan        = 1020            # End channel for data cube
niter        = 200             # niter in imagr
uvwtfn       = 'N'             # UVWTFN in IMAGR
robust       =  0              # ROBUST in IMAGR
beam         = [0,0,0]         # Restoring beam: bmaj, bmin, PA
antennas     = [1,2,3,4,5,6,7,8,9]   # Also used for RPOSSM

######################
#   SAD Parameters   #
######################

min_snr      = 10                            # Fit only Peaks with SNR>min_snr
dyna         = 0.1                           # Don't fit peaks less than 
                                             # dyna*peak in channel
####################
# Self Calibration #
####################

phase_loop = [30./60,20./60.]                # Solints for Phase selfcal
amp_loop   = [1200.,360, 60.,30.]
dofit      = [[0], [0], [0],[0]]


###############
# Input Files #
###############

filename[0] = 'br210c9d_geo.idifits'
outname[0]  = 'Br210c9d-G'
outclass[0] = 'UVDATA'
outdisk[0]  = 1

filename[1] = 'br210c9d_cont.idifits'
outname[1]  = 'Br210c9d-C'
outclass[1] = 'UVDATA'
outdisk[1]  = 1

filename[2] = 'br210c9d_line.idifits'
outname[2]  = 'Br210c9d-L'
outclass[2] = 'UVDATA'
outdisk[2]  = 2

                                    # Optional parameters for each file
# outdisk[2]   = 3                  # AIPS disk for this file (if != defdisk)
# flagfile[2]  = 'flagfile.uvflg'   # flag file for UVFLG
# antabfile[2] = 'antabfile.antab'  # antab file for ANTAB
# nfiles[2]    = 0                  # FITLD parameter NFILES
# ncount[2]    = 6                  # FITLD parameter NCOUNT
# doconcat[2]  = 1                  # FITLD parameter DOCONCAT

#################
# Control Flags #
#################

####################
# Data preparation #
####################
download_flag   = 0        # Download data from archive?
load_flag       = 0        # Load data from disk?
listr_flag      = 0        # Print out LISTR?
get_key_flag    = 0        # Download key-file from archive

#####################
# geoblock analysis #
#####################
RDBE_check      = 0        # Check Geoblock data for RDBE errors?
geo_prep_flag   = 0        # Run TECOR and EOP corrections?
geo_fringe_flag = 0        # Fringe fit the data?
doprt_flag      = 0        # Print out files for fitting?
dofit_flag      = 0        # Start position fitting?
doplot_flag     = 0        # Plot fit_geoblock output?

#####################
# phaseref analysis #
#####################
tasav_flag      = 0        # Run tasav on original tables?
restore_su_flag = 0        # Restore original SU table from TASAV?
restore_fg_flag = 0        # Restore original FG table from TASAV?

pr_prep_flag    = 0        # Run TECOR, EOPs, ATMOS, PANG, and position shift?
apcal_flag      = 0        # Do amplitude calibration?
pr_fringe_flag  = 0        # Do manual phase cal?
do_band_flag    = 0        # Do bandpass calibration?
cvel_flag       = 0        # Run CVEL on line data?
#
possm_flag      = 0        # Run POSSM to select channel?
ma_fringe_flag  = 0        # Fringe fit on maser and calibrate data?
co_fringe_flag  = 0        # Fringe fit on continuum and calibrate data?

snflg_flag      = 0        # Run snflg after fringe?
print_sn_flag   = 0        # Print out last SN table?
plot_tables     =-1        # Plot SN-Tables of which data set (-1 for no plots)

split_flag      = 0       # Split calibrated data?
fittp_flag      = 0       # Write calibrated data to disk?
co_imagr_flag   = 0       # Image continuum (target) sources?
#
ma_imagr_flag   = 0       # Image one channel of line (cvel-)sources?
cube_imagr_flag = 0       # Image data cube of line (cvel-)sources?

rpossm_flag     = 0       # Produce first epoch spectrum?
ma_sad_flag     = 0       # Run SAD on cube?
plot_map        = 1       # Make spot map

phase_cal_flag  = 0        # Phase self-calibration on calibrator (cont)
amp_cal_flag    = 0        # Amplitude self-calibration on calibrator (cont)
apply_selfcal   = 0        # Apply self calibration to target sources?


####################
# Download options #
####################

# Only for archive download use arch_user = 'nopass' for e2earchive

if download_flag==1:
    file        = range(n)
    arch_user   = 'nopass'
    arch_pass   = 'nopass'
    file[0]='VLBA_BR210C9_br210c9_BIN0_SRC0_0_151104T153935.idifits'
    file[1]='VLBA_BR210C9_br210c9_BIN0_SRC0_4_151104T153935.idifits'
    file[2]='VLBA_BR210C9_hiresc9_BIN0_SRC0_3_151104T154136.idifits'

##############################################################################
# Start main script

execfile(r''+def_file)

#
##############################################################################



