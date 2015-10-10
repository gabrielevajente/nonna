# Auxiliary functions for the Non-stationary Noise Analysis tool
# Gabriele Vajente 2015-06-09

import nds2
import scipy
import scipy.signal
#import pylab
import numpy
from pylal.Fr import frgetvect1d
import subprocess

# Read data, compute the BLRMS and precondition the auxiliary channels ###################
def nonna_get_data(target_channel, aux_channels, gps_start, duration, band_freqs, outfs, 
				   fs = 4096, nds_server = 'nds.ligo-wa.caltech.edu', nds_port = 31200):
	"""
	This function reads data and prepares it for the analysis. It computes the band-limited
	RMS and downsamples the auxiliary channels.
	
	Input arguments:
	target_channel = compute the BLRMS of this channel
	aux_channels   = all the slow channels that will be used to predict the BLRMS time
	                 variation
	gps_start      = start reading data at this time
	duration       = number of seconds of data to read
	band_freqs     = [fmin, fmax] edge frequencies of the band used to compute the BLRMS
	outfs          = output sampling frequency for BLRMS and the other channels
	fs             = first downsample the target channel to this sampling rate, to avoid
					 numerical instabilities of the band-pass filter, which would return 
					 only a bunch of NaNs
	nds_server     = address of the NDS2 server
	nds_port       = port used to contact the server
	
	Output data:
	t     = vector of time values for each sample
	blrms = vector containing the BLRMS (at the rate outfs)
	aux   = matrix of all auxiliary channels, downsampled to the output rate
	
	Note that the first two and last two seconds of data are discarded, to cope with 
	the band-pass and low-pass filter transients.
	"""
	
	##### READ DATA 
	
	# open connection
	conn = nds2.connection(nds_server, nds_port)
	# read all data
	buffers = conn.fetch(gps_start, gps_start+duration, [target_channel]+aux_channels)

	##### COMPUTE BLRMS 

	# decimate target signal
	tg = scipy.signal.decimate(buffers[0].data, int(buffers[0].length/duration/fs))
	# band pass
	b,a = scipy.signal.butter(6, scipy.array(band_freqs)/(fs/2.), btype='bandpass')
	tg = scipy.signal.filtfilt(b, a, tg)
	# square and low pass
	tg = tg**2
	b,a = scipy.signal.butter(4, outfs/(fs/2.), btype='lowpass')
	blrms = scipy.signal.filtfilt(b, a, tg)
	# decimate
	blrms = blrms[::fs/outfs]

	##### DECIMATE THE OTHER CHANNELS 

	aux = scipy.zeros((duration*outfs, len(buffers)-1))

	# loop over all channels
	for i in range(1,len(buffers)):
		# low pass and decimate
		fs_aux = int(buffers[i].length/duration)
		b,a = scipy.signal.butter(4, outfs/(fs_aux/2.), btype='lowpass')
		x = scipy.signal.filtfilt(b, a, buffers[i].data)
		aux[:,i-1] = x[::fs_aux/outfs]

	# get rid of initial and final transients
	aux = aux[2*outfs:-2*outfs]
	blrms = blrms[2*outfs:-2*outfs]
	t = scipy.arange(0,len(blrms))/float(outfs)
	
	# RETURN RESULTS
	return t, blrms, aux

# wrapper around the LIGO function to find where data is, returns a list of files
def find_LIGO_data(observatory, gpsb, gpse):
    o = subprocess.Popen(["/usr/bin/gw_data_find", "-o", observatory[0],
                    		"-t", observatory[0] + "1_R", "-s", str(gpsb), "-e", str(gpse), "-u", "file"],
                        	stdout=subprocess.PIPE).communicate()[0]        
    return o.splitlines()

# read data directly from frame files on disk
def nonna_get_data_from_disk(channel, gps_start, duration, outfs=-1, verbose=False):
	"""
	This function reads data directly from disk, using gw_data_find to locate the gwf
	files.

	Input arguments:
	channel = name of the signal to load
	gps_start      = start reading data at this time
	duration       = number of seconds of data to read
	outfs          = output sampling frequency for BLRMS and the other channels
	"""

	# get the list of files
	files = find_LIGO_data('H1', int(gps_start), int(gps_start+duration))
	# loop over all files
	for i,f in enumerate(files):
		if verbose:
			print 'Reading data from file %s (%d/%d)' % (f, i, len(files))
		# the file name will tell use what times are available
		t = f.split('.')[-2].split('-')[-2:]
		gps_file = int(t[0])
		gps_span = int(t[1])
		# get all the data we can from this file
		gps0 = max(gps_file, int(gps_start))
		gps1 = min(gps_file+gps_span, int(gps_start+duration))
		buffer = frgetvect1d(f[16:], channel, gps0, gps1-gps0, 0)
		if i == 0:
			# at the first read, allocate the vector
			fs = int(1/buffer[3])
			data = numpy.zeros(fs*duration)
		data[(gps0-gps_start)*fs:(gps1-gps_start)*fs] = numpy.array(buffer[0])
	if outfs == -1:
		# return the whole data
		return data
	else:
		# downsample
		if verbose:
			print 'Downsampling'
		infs = (len(data)/duration)
		b,a = scipy.signal.butter(4, outfs/(infs/2.), btype='lowpass')
		data = scipy.signal.lfilter(b, a, data)
		return data[::infs/outfs]

# compute BLRMS of a signal
def nonna_blrms(signal, f1, f2, infs, outfs, remove=10):
	"""
	Compute the band-limited RMS of the input signal.

	Input arguments:
	signal = the input signal
	f1,f2 = corner frequencies of the band
	infs = sampling frequency of signal
	outfs = desired sampling frequency of the output BLRMS
	remove = number of seconds of data to throw away at beginning
	"""
	# check the ratio of BLRMS corner frequency and sampling frequency
	# if too large we have to decimate the signal first
	if infs/f1 > 30:
		print 'Decimating signal'
		decfs = infs/8
                b,a = scipy.signal.butter(4, decfs/(infs/2.), btype='lowpass')
                x = scipy.signal.lfilter(b, a, signal)
                signal = x[::infs/decfs]
		infs = decfs
	
	# define band-pass filter
	b,a = scipy.signal.butter(6, scipy.array([f1, f2])/(infs/2.), btype='bandpass')
	# band pass
	signal = scipy.signal.filtfilt(b, a, signal)
	# square
	signal = signal**2
	# low pass and decimate
	b,a = scipy.signal.butter(4, outfs/(infs/2.), btype='lowpass')
        signal = scipy.signal.filtfilt(b, a, signal)
        # decimate
        signal = signal[::infs/outfs]
	signal = signal[outfs*remove:]
	return signal

# Select data based on outlier removal ###################################################
def nonna_select_data(data, outlier_threshold, level='high'):
	"""
	This function returns a list of indexed after identifying the main outliers. It applies
	a cut on the data to remove exactly a fraction (1-outlier_threshold) of all data points.
	By default the cut is applied only at the higher end of the data values, but the 
	parameter level can be used to change this
	
	Input arguments:
	data              = vector containing all data points
	outlier_threshold = remove outliers until we are left with exactly this fraction of the
	                    original data
	level             = 'high|low|both' determines if the outliers are removed only from the
					    high values end, the low values end of both ends.
					    
	Output:
	idx               = index of selected (good) data
	"""
	
	# histogram all the data values
	n,x = scipy.histogram(data, len(data)/10)
	# compute the cumulative distribution and normalize
	nn = scipy.cumsum(n)
	nn = nn / float(max(nn))
	
	if level=='high':
		# select the value such that a fraction outlier_threshold of the data lies below it
		if outlier_threshold < 1:
			val = x[pylab.find(nn/float(max(nn)) >= outlier_threshold)[0]]
		else:
			val = max(data)
		# use that fraction of data only
		idx = data <= val 
	elif level=='low':
		# select the value such that a fraction outlier_threshold of the data lies above it
		if outlier_threshold < 1:
			val = x[pylab.find(nn/float(max(nn)) <= (1-outlier_threshold))[-1]]
		else:
			val = min(data)
		# use that fraction of data only
		idx = data >= val 
	elif level=='both':		
		# select the value such that a fraction outlier_threshold/2 of the data lies below it
		if outlier_threshold < 1:
			Hval = x[pylab.find(nn/float(max(nn)) >= 1-(1-outlier_threshold)/2)[0]]
		else:
			Hval = max(data)	
		# select the value such that a fraction outlier_threshold/2 of the data lies above it
		if outlier_threshold < 1:
			Lval = x[pylab.find(nn/float(max(nn)) <= (1-outlier_threshold)/2)[-1]]
		else:
			Lval = min(data)
  		# use that fraction of data only
		idx = scipy.logical_and(data >= Lval, data <= Hval) 
	
	return idx
	
# Produce the Least squares prediction of the target signal using the auxiliary signals
# and their powers if requested
def nonna_lsq(target, aux, idx=(), names=(), order=2):
	"""
	This function returns the coefficients of the least square prediction of the target
	signal, using the auxiliary signals and their powers, as specified by the order argument.
	
	Input arguments:
	target = target signal
	aux    = matrix of auxiliary signals
	idx    = boolean vector to select a subset of the data for the LSQ fit
	order  = order of the polynomial of aux signals to be used in the fit, default is 2
	names  = list of the auxiliary signal names
	
	Output:
	p      = list of coefficients
	X      = matrix of the signals used in the reconstruction
	cnames = list of the corresponding signals
	
	Note that the mean will be removed from the auxiliary signals. 
	"""
	# number of auxiliary channels
	naux = scipy.shape(aux[1])
	
	if len(names) == 0:
		# since the user didn't provide signal names, let's build some
		names = map(lambda x: 'S'+str(x), scipy.arange(naux)+1)
		
	if len(idx) == 0:
		# no index means use all
		idx = array(target, dtype=bool)
		idx[:] = True
	
	##### PREPARE CHANNELS FOR LSQ PREDICTION 

	# prepare channels and their squared values
	X = scipy.zeros((scipy.shape(aux)[0], order*scipy.shape(aux)[1]+1))
	cnames = []
	for i in range(scipy.shape(aux)[1]):
		for j in range(order):
			# add the (j+1)th power of the signal after removing the mean
			X[:,order*i+j] = numpy.power((aux[:,i] - scipy.mean(aux[idx,i])), j+1)
			# then remove the mean of the result
			X[:,order*i+j] = X[:,order*i+j] - scipy.mean(X[idx,order*i+j])
			# save the name, including the power
			if j==0:
				cnames.append(names[i])
			else:
				cnames.append(names[i]+'^'+str(j+1))
				
	# add a constant at the end of the list
	X[:,-1] = 1
	cnames.append('1')
	# convert to matrix object for simpler manipulation
	X = scipy.mat(X)
	
	##### best estimate of coefficients to minimize the squared error
	p = scipy.linalg.inv(X[idx,:].T * X[idx,:]) * X[idx,:].T * scipy.mat(target[idx]).T

	# return all the results
	return p, X, cnames
	
# Produce the least square estimation with channel ranking 
def nonna_lsq_signal_ranking(target, aux, idx=(), names=(), order=2):
	"""
	This function returns the coefficients of the least square prediction of the target
	signal, using the auxiliary signals and their powers, as specified by the order argument.
	It also returns a ranking of the signals in terms of their contribution to the 
	reduction of the residual error.
	
	Input arguments:
	target = target signal
	aux    = matrix of auxiliary signals
	idx    = boolean vector to select a subset of the data for the LSQ fit
	order  = order of the polynomial of aux signals to be used in the fit, default is 2
	names  = list of the auxiliary signal names
	
	Output:
	p      = list of coefficients
	X      = matrix of the signals used in the reconstruction
	cnames = list of the corresponding signals
	id     = list of signal indexes, in order of reducing relevance
	de     = list of the residual error reduction provided by including each signal, in
			 the same order as the list above 
	
	Note that the mean will be removed from the auxiliary signals. 
	"""
	
	# first estimation with all channels
	p0, X, cnames = nonna_lsq(target, aux, idx=idx, names=names, order=order)

	# convert B to matrix for convenience and remove the mean (to avoid counting in the
	# constant term in the ranking)
	B = scipy.mat(target - scipy.mean(target[idx]))

	# define the function used to compute the residual error
	def error(p):
		return scipy.mean(scipy.square(B[:,idx].T - X[idx,:]*p))

	# compute the initial error when all channels are used
	e0 = error(p0)
	print '0) initial error %g' % e0

	# init variables to store residuals and indexes at each iteration
	e  = scipy.zeros((scipy.shape(X)[1],))
	id = scipy.zeros((scipy.shape(X)[1],), dtype=int)
	# init all indexes to dummy values at the beginning (no channel removed yet)
	id[:] = -1
	
	# Repeat the estimate of the best fit with all possible reduced set of signals. We'll 
	# remove one at each step
	for i in range(scipy.shape(X)[1]):
		# this is going to be the list of the new residual errors when we removed one
		# additional channel
		newerrors = scipy.zeros((scipy.shape(X)[1],1))
		# loop over all channels and remove one by one
		for j in range(scipy.shape(X)[1]):
			# check if this channel was already removed
			if not any(id == j):
				# remove all the channels that are already in the list, plus the one under
				# consideration
				ind = scipy.setdiff1d(range(scipy.shape(X)[1]), id)
				ind = scipy.setdiff1d(ind, j)
				# start with an empty set of coefficients
				pp = scipy.zeros((scipy.shape(X)[1],1))
				# compute the best estimate of coefficients
				if len(ind) != 0:
					pp[ind] = scipy.linalg.inv(X[idx,:][:,ind].T * X[idx,:][:,ind]) * X[idx,:][:,ind].T * B[:,idx].T
				# and finally compute the new residual errors
				newerrors[j] = error(pp)
			else:
				# we already used this channel, let's make the error infinite so it won't be
				# picked later on
				newerrors[j] = scipy.inf
		
		# Now we have to choose the channel that (when removed) still gives the minimum 
		# residual error
		e[i] = min(newerrors)
		id[i] = scipy.argmin(newerrors)
		# Print out some information
		print '%d) new error %g (removed channel %s)' % (i+1, e[i], cnames[id[i]])
	
	# Final steps, build incremental residual error worsening
	de = scipy.diff(scipy.concatenate((scipy.array([e0]), e[:])))
	
	# sort them out
	ii = scipy.argsort(de)
	de = de[ii[::-1]]
	id = id[ii[::-1]]
	
	# return results
	return p0, X, cnames, id, de

def nonna_readsegmentfile(file):
	f = open(file, 'r')
	L = f.readlines()
	f.close()

	gps1 = []
	gps2 = []
	for line in L:
		x = line.split()
		gps1.append(int(x[1]))
		gps2.append(int(x[2]))
	return gps1, gps2

