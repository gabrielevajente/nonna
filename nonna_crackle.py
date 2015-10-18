#!/usr/bin/python

# Non-stationary Noise Analysis (NonNA) - Gabriele Vajente - version 2015-06-09 
#
# Compute the band-limited RMS of a signal in a given band, and try to make the best 
# possible prediction of its slow time variation using a set of auxiliary channels.
# The script produces two output plots:
#  1) time series of the BLRMS and best fit obtained with the auxiliary channels
#  2) ranking of the channels, in terms of how much they contribute to the reduction of
#     the error in the prediction 
# 
# Most of the configuration parameters are in the following lines (20-53). 

import nds2
import scipy
import scipy.signal
import pylab
from nonna_functions import *

##### CONFIGURATION BITS #################################################################

# NDS connection configuration
nds_server = 'cymac2.caltech.edu'
nds_port = 8088

# target channel and frequency band for BLRMS computation
target = 'X1:KR2-BBOARD_LOCK_SERV_IN1_DQ'
band_freqs = [10., 40.]

# processing and output frequencies
fs = 256	# decimate channel to this frequency before band-passing (this is necessary
			# to avoid the band pass filter to go unstable and return all NaNs)	
outfs = 1	# sampling frequency for the BLRMS production

# list of channels using for the prediction
channels = ['X1:KR2-BBOARD_SENS_X1_OUT_DQ', 'X1:KR2-BBOARD_SENS_X2_OUT_DQ', 
			'X1:KR2-BBOARD_SENS_Y1_OUT_DQ', 'X1:KR2-BBOARD_SENS_Y2_OUT_DQ', 
            #'X1:KR2-BBOARD_SENS_Z1_OUT_DQ', 'X1:KR2-BBOARD_SENS_Z2_OUT_DQ', 
            'X1:KR2-SUSP_SENS_A_OUT_DQ', 'X1:KR2-SUSP_SENS_B_OUT_DQ', 
            'X1:KR2-SUSP_SENS_C_OUT_DQ', 'X1:KR2-SUSP_SENS_D_OUT_DQ', 
            'X1:KR2-SUSP_SENS_E_OUT_DQ', 'X1:KR2-SUSP_SENS_F_OUT_DQ']

# Start time
gps = 1128369557
# Length of data in seconds
dt = 1800

# if you want to remove outliers from the data, you can set the threshold below to 
# represent how much data you want to keep. For example, setting the threshold to 0.9 as
# below means that NonNa will find a value of the BLRMS such that 90% of the data will be 
# below it. The analysis will be performed using only that data.
# Of course, setting the threshold to 1 means that no data section will be applied. 
outlier_threshold = 1

##### READ DATA AND COMPUTE BLRMS ########################################################

# Read the target channel and the auxiliary ones, compute the BLRMS, and downsample 
# everything at the desired output frequency
t,blrms,aux = nonna_get_data(target, channels, gps, dt, band_freqs, outfs, 
							 fs=fs, nds_server=nds_server, nds_port=nds_port)

##### DATA SELECTION #####################################################################

# throw away outliers until we are left with a chosen fraction of the datapoints, set 
# by the threshold outlier_threshold.
idx = nonna_select_data(blrms, outlier_threshold, level='high')

##### LSQ PREDICTION USING ALL CHANNELS ##################################################

# compute the best reconstruction of the BLRMS using the auxiliary signals, with powers
# up to 2.
p, X, cnames = nonna_lsq(blrms, aux, idx=idx, names=channels, order=2)

##### CHANNEL RANKING BY REMOVAL #########################################################

# Same as above, but this time returns also the ranking of the signals in order of 
# descending relevance (how much the inclusion of each channel improves the residual
# error
p, X, cnames, id, de = nonna_lsq_signal_ranking(blrms, aux, idx=idx, names=channels, order=2)

##### PLOT RESULTS #######################################################################

# plot the measured BLRMS and the prediction
pylab.ion()
pylab.figure()
pylab.plot(t, blrms, 'b', t, X*p, 'r')
pylab.legend(('BLRMS', 'Best prediction'))
pylab.grid()

# plot the channel ranking
pylab.figure()
pylab.barh(range(len(de),0,-1), de)
pylab.ylim((0, len(de)+1))
pylab.xlim((0, 1.2*max(de)))
pylab.yticks(scipy.arange(len(de)+0.5,0.5,-1), array(cnames)[id])
pylab.subplots_adjust(left=0.4, right=0.9, top=0.9, bottom=0.1)
pylab.xlabel('Residual error improvement')
pylab.grid()
