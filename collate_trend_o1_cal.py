from nonna_functions import *
import numpy

gps1,gps2 = nonna_readsegmentfile('segments_h1.txt')

t = []
blrms = [[],[],[],[],[],[],[],[],[],[],[],[]]

for i in range(len(gps1)):
	print 'Segment %d: %d - %d (%d seconds)' % (i, gps1[i], gps2[i], gps2[i] - gps1[i])
	if gps2[i] - gps1[1] > 20:
		t.append(numpy.arange(gps2[i]-gps1[i]-10) + gps1[i]+10)
		blrms[0].append(numpy.loadtxt('../../blrms_o1/cal_blrms_10_20_%d.txt' % (gps1[i]+10)))
		blrms[1].append(numpy.loadtxt('../../blrms_o1/cal_blrms_20_30_%d.txt' % (gps1[i]+10)))
		blrms[2].append(numpy.loadtxt('../../blrms_o1/cal_blrms_15_25_%d.txt' % (gps1[i]+10)))
		blrms[3].append(numpy.loadtxt('../../blrms_o1/cal_blrms_43_58_%d.txt' % (gps1[i]+10)))
		blrms[4].append(numpy.loadtxt('../../blrms_o1/cal_blrms_70_100_%d.txt' % (gps1[i]+10)))
		blrms[5].append(numpy.loadtxt('../../blrms_o1/cal_blrms_140_160_%d.txt' % (gps1[i]+10)))
        	blrms[6].append(numpy.loadtxt('../../blrms_o1/cal_blrms_200_220_%d.txt' % (gps1[i]+10)))
        	blrms[7].append(numpy.loadtxt('../../blrms_o1/cal_blrms_305_370_%d.txt' % (gps1[i]+10)))
       	 	blrms[8].append(numpy.loadtxt('../../blrms_o1/cal_blrms_400_490_%d.txt' % (gps1[i]+10)))
		blrms[9].append(numpy.loadtxt('../../blrms_o1/cal_blrms_700_900_%d.txt' % (gps1[i]+10)))
		blrms[10].append(numpy.loadtxt('../../blrms_o1/cal_blrms_240_260_%d.txt' % (gps1[i]+10)))
		blrms[11].append(numpy.loadtxt('../../blrms_o1/cal_blrms_340_360_%d.txt' % (gps1[i]+10)))

for i in range(len(blrms)):
	blrms[i] = numpy.concatenate(blrms[i])
t = numpy.concatenate(t)

numpy.savetxt('o1_cal_gps.txt', t)
numpy.savetxt('o1_cal_blrms_10_20.txt', blrms[0])
numpy.savetxt('o1_cal_blrms_20_30.txt', blrms[1])
numpy.savetxt('o1_cal_blrms_15_25.txt', blrms[2])
numpy.savetxt('o1_cal_blrms_43_58.txt', blrms[3])
numpy.savetxt('o1_cal_blrms_70_100.txt', blrms[4])
numpy.savetxt('o1_cal_blrms_140_160.txt', blrms[5])
numpy.savetxt('o1_cal_blrms_200_220.txt', blrms[6])
numpy.savetxt('o1_cal_blrms_305_370.txt', blrms[7])
numpy.savetxt('o1_cal_blrms_400_490.txt', blrms[8])
numpy.savetxt('o1_cal_blrms_700_900.txt', blrms[9])
numpy.savetxt('o1_cal_blrms_240_260.txt', blrms[10])
numpy.savetxt('o1_cal_blrms_340_360.txt', blrms[11])
