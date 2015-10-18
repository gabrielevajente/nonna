from numpy import *
import os.path
import os
import glob

darm_blrms_folder = os.path.expanduser('~/Documents/MATLAB/aLIGO/blrms_o1/LHO')
darm_blrms_name = 'lho_blrms_fft_'
aux_signals_folder = os.path.expanduser('~/Documents/MATLAB/aLIGO/blrms_o1/LHO')
aux_signals_name = 'lho_isignd'
bands = [[10, 20], [20, 30], [30, 40], [43,58], [70, 100], \
		 [140, 160], [200, 220], [240, 260], [305, 370], \
		 [340, 360], [400, 450], [700, 900]]

# load all DARM BLRMS data 
files = glob.glob(darm_blrms_folder + '/' + darm_blrms_name + '*')
b = []
for f in files:
	# get start GPS form file name
	gps0 = int(f.split('/')[-1].split('_')[-1].split('.')[0])
	# load data
	bb = loadtxt(f)
	try:
		# add start GPS time
		bb[:,0] = bb[:,0] + gps0
		# append to existing data
		if len(b) == 0:
			b = bb;
		else:
			b = vstack((b, bb))
	except:
		pass

# time in GPS and days from 2015-09-19 00:00 UTC
t = b[:,0]
T = (t - 1126569617)/86400 
b = b[:,1:]

# load all auxiliary data 
files = glob.glob(darm_blrms_folder + '/' + aux_signals_name + '*')
a = []
for f in files:
	# load data
	aa = loadtxt(f)
	try:
		# append to existing data
		if len(a) == 0:
			a = aa;
		else:
			a = vstack((a, aa))
	except:
		pass

# time in GPS and days from 2015-09-19 00:00 UTC
ta = a[:,0]
Ta = (ta - 1126569617)/86400 
a = a[:,1:]

# interpolate auxiliary signals to the same times as the DARM BLRMS
aux = zeros((len(t), shape(a)[1]))
for i in range(shape(a)[1]):
	aux[:,i] = interp(t, ta, a[:,i])

# y axis limits for good plots
ymax = [5e-21, 5e-22, 1e-19,   1e-20,   3e-19,   1e-18,   2e-18,   2e-18,   1.4e-16,   2e-17, 1e-17,   1e-16]
ymin = [0e-21, 0e-22, 0.5e-19, 0.1e-20, 0.4e-19, 0.2e-18, 0.5e-18, 0.8e-18, 0.8e-16, 0e-17, 0.5e-17, 0.5e-16]

# 