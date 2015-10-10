from nonna_functions import *
import numpy
from scipy.signal import *

gps1,gps2 = nonna_readsegmentfile('segments_l1.txt')

for i in range(len(gps1)):
	print 'Segment %d: %d - %d (%d seconds)' % (i, gps1[i], gps2[i], gps2[i] - gps1[i])
	print 'Read data...'
	d = nonna_get_data_from_disk('L1:CAL-DELTAL_EXTERNAL_DQ', gps1[i], gps2[i] - gps1[i], outfs=1024, verbose=False)
	
	# calibrate removing whitening
	N,D = zpk2tf(-2*numpy.pi*numpy.array([100, 100, 100, 100, 100]), -2*numpy.pi*numpy.array([1, 1, 1, 1, 1]), 1e-10)	
	b,a = bilinear(N,D,1024)
	d = lfilter(b,a,d)

	print 'Compute BLRMS...'
	b1 = nonna_blrms(d, 10, 20, 1024, 1)
	b2 = nonna_blrms(d, 20, 30, 1024, 1)
	b9 = nonna_blrms(d, 15, 25, 1024, 1)
	b11 = nonna_blrms(d, 240, 260, 1024, 1)
	b12 = nonna_blrms(d, 340, 360, 1024, 1)
	b3 = nonna_blrms(d, 43, 58, 1024, 1)
	b4 = nonna_blrms(d, 70, 100, 1024, 1)
	b5 = nonna_blrms(d, 140, 160, 1024, 1)
	b6 = nonna_blrms(d, 200, 220, 1024, 1)
	b7 = nonna_blrms(d, 305, 370, 1024, 1)
	b8 = nonna_blrms(d, 400, 450, 1024, 1)
        b10 = nonna_blrms(d, 700, 900, 1024, 1)
	print 'BLRMS 10 20 Hz   = %e' % numpy.mean(b1)
        print 'BLRMS 20 30 Hz   = %e' % numpy.mean(b2)
	print 'BLRMS 15-25 Hz   = %e' % numpy.mean(b9)
        print 'BLRMS 43 58 Hz   = %e' % numpy.mean(b3)
        print 'BLRMS 70 100 Hz  = %e' % numpy.mean(b4)
        print 'BLRMS 140 160 Hz = %e' % numpy.mean(b5)
        print 'BLRMS 200 220 Hz = %e' % numpy.mean(b6)
        print 'BLRMS 305 370 Hz = %e' % numpy.mean(b7)
        print 'BLRMS 240 260 Hz = %e' % numpy.mean(b11)
        print 'BLRMS 340 360 Hz = %e' % numpy.mean(b12)
	print 'BLRMS 400 450 Hz = %e' % numpy.mean(b8)
	print 'BLRMS 700 900 Hz = %e' % numpy.mean(b10)

	print 'Saving result'
	numpy.savetxt('../../blrms_o1/cal_blrms_10_20_%d.txt' % (gps1[i]+10), b1)
	numpy.savetxt('../../blrms_o1/cal_blrms_20_30_%d.txt' % (gps1[i]+10), b2)
	numpy.savetxt('../../blrms_o1/cal_blrms_15_25_%d.txt' % (gps1[i]+10), b9)
	numpy.savetxt('../../blrms_o1/cal_blrms_43_58_%d.txt' % (gps1[i]+10), b3)
	numpy.savetxt('../../blrms_o1/cal_blrms_70_100_%d.txt' % (gps1[i]+10), b4)
	numpy.savetxt('../../blrms_o1/cal_blrms_140_160_%d.txt' % (gps1[i]+10), b5)
        numpy.savetxt('../../blrms_o1/cal_blrms_200_220_%d.txt' % (gps1[i]+10), b6)
        numpy.savetxt('../../blrms_o1/cal_blrms_305_370_%d.txt' % (gps1[i]+10), b7)
        numpy.savetxt('../../blrms_o1/cal_blrms_400_450_%d.txt' % (gps1[i]+10), b8)
	numpy.savetxt('../../blrms_o1/cal_blrms_700_900_%d.txt' % (gps1[i]+10), b10)
        numpy.savetxt('../../blrms_o1/cal_blrms_240_260_%d.txt' % (gps1[i]+10), b11)
        numpy.savetxt('../../blrms_o1/cal_blrms_340_360_%d.txt' % (gps1[i]+10), b12)
