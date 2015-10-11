from nonna_functions import *
import numpy
from scipy.signal import *

gps1,gps2 = nonna_readsegmentfile('segments_h1.txt')

# frequency bands for BLRMS
bands = [[10, 20], [15, 25], [20, 30], [30, 40], [43,58], [70, 100], [140, 160], [200, 220], [240, 260], [305, 370], \
	 [340, 360], [400, 450], [700, 900]]
# output period
Tout = 60
# FFT resolution
Tfft = 5
# data samping frequency
fs = 2048

for i in range(len(gps1)):
	print 'Segment %d: %d - %d (%d seconds)' % (i, gps1[i], gps2[i], gps2[i] - gps1[i])
	print 'Read data...'
	d = nonna_get_data_from_disk('H1:CAL-DELTAL_EXTERNAL_DQ', gps1[i], gps2[i] - gps1[i], outfs=fs, verbose=False)
	
	# calibrate removing whitening
	#N,D = zpk2tf(-2*numpy.pi*numpy.array([100, 100, 100, 100, 100]), -2*numpy.pi*numpy.array([1, 1, 1, 1, 1]), 1e-10)	
	#b,a = bilinear(N,D,fs)
	#d = lfilter(b,a,d)

	print 'Compute BLRMS...'
	t,b = nonna_blrms_fft(d, bands, fs, Tout, Tfft)
	print 'Saving...'
	X = numpy.zeros((len(t), len(bands)+1))
	X[:,0] = t
	X[:,1:] = b
	numpy.savetxt('../../blrms_o1/lho_blrms_fft_%d.txt' % (gps1[i]), X)
