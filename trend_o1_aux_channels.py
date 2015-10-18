from nonna_functions import *
import numpy
from scipy.signal import *
import nds2

gps1,gps2 = nonna_readsegmentfile('segments_h1.txt')

# channels
channels = [\
	'H1:ISI-GND_STS_ETMX_X_BLRMS_100M_300M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_10_30.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_1_3.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_300M_1.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_30M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_30M_100M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_X_BLRMS_30_100.mean,m-trend',\
	'H1:ISI-GND_STS_ETMX_X_BLRMS_3_10.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_100M_300M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_10_30.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_1_3.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_300M_1.mean,m-trend',\
	'H1:ISI-GND_STS_ETMX_Y_BLRMS_30M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_30M_100M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_30_100.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Y_BLRMS_3_10.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_100M_300M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_10_30.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_1_3.mean,m-trend',\
	'H1:ISI-GND_STS_ETMX_Z_BLRMS_300M_1.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_30M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_30M_100M.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_30_100.mean,m-trend',\
 	'H1:ISI-GND_STS_ETMX_Z_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_X_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Y_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ETMY_Z_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_X_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Y_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_ITMY_Z_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_X_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Y_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM5_Z_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_X_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Y_BLRMS_3_10.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_100M_300M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_10_30.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_1_3.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_300M_1.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_30M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_30M_100M.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_30_100.mean,m-trend',\
        'H1:ISI-GND_STS_HAM2_Z_BLRMS_3_10.mean,m-trend']

conn = nds2.connection('nds.ligo-wa.caltech.edu')

for i in range(len(gps1)):
	print 'Segment %d of %d: %d - %d (%d seconds)' % (i, len(gps1), gps1[i], gps2[i], gps2[i] - gps1[i])
	print 'Read data...'

	# round GPS times to multiple of 60 seconds
	g1 = gps1[i] - gps1[i]%60 + 60
	g2 = gps2[i] - gps2[i]%60

	if g2<=g1:
		continue

	buf = conn.fetch(g1, g2, channels)
	
	print 'Saving...'
	X = numpy.zeros(((g2-g1)/60, len(channels)+1))
	X[:,0] = arange(g1, g2, 60)
	for j,b in enumerate(buf):
		X[:,j+1] = b.data
	numpy.savetxt('../../blrms_o1/lho_isignd_%d.txt' % (gps1[i]), X)
