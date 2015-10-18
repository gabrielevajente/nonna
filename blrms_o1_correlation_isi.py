from numpy import *
from nonna_functions import *
from numpy import *
import os.path

darm_blrms_folder = os.path.expanduser('~/Documents/MATLAB/aLIGO/blrms_o1')
darm_blrms_name = 'lho_blrms_darm.txt'
aux_signals_folder = os.path.expanduser('~/Documents/MATLAB/aLIGO/blrms_o1')
aux_signals_name = 'lho_isi_gnd.txt'
aux_signals_list = 'lho_isi_gnd_channels.txt'
bands = [[10, 20], [20, 30], [30, 40], [43,58], [70, 100], \
		 [140, 160], [200, 220], [240, 260], [305, 370], \
		 [340, 360], [400, 450], [700, 900]]

# load all data
b = loadtxt(darm_blrms_folder + '/' + darm_blrms_name)
a = loadtxt(aux_signals_folder + '/' + aux_signals_name)

t = b[:,0]
b = b[:,1:]
a = a[:,1:]

f = open(aux_signals_folder + '/' + aux_signals_list, 'r')
channels = map(lambda x: x[:-14], f.readlines())
f.close()

##### LSQ PREDICTION USING ALL CHANNELS ##################################################
p = zeros((shape(a)[1]*2+1, shape(b)[1]))
for i in range(shape(b)[1]):
	pp, X, cnames = nonna_lsq(b[:,i], a, names=channels, order=2)
	p[:,i] = ravel(pp)
	
bfit = X * p
	
##### CHANNEL RANKING BY REMOVAL #########################################################

# Same as above, but this time returns also the ranking of the signals in order of 
# descending relevance (how much the inclusion of each channel improves the residual
# error
p = zeros((shape(a)[1]*2+1, shape(b)[1]))
id = zeros((shape(a)[1]*2+1, shape(b)[1]))
de = zeros((shape(a)[1]*2+1, shape(b)[1]))

for i in range(shape(b)[1]):
	pd, X, cnames, iid, dde = nonna_lsq_signal_ranking(b[:,i], a, names=channels, order=2)
	p[:,i] = ravel(pp)
	id[:,i] = ravel(iid)
	de[:,i] = ravel(dde)
	
	pylab.ion()
	pylab.figure()
	pylab.plot(t, b[:,i], '.b', t, X*p, '.r')
	pylab.legend(('BLRMS %d-%d Hz' % (bands[i][0], bands[i][1]), 'Best prediction'))
	pylab.grid()
	savefig('fit_%d_%d.png' % (bands[i][0], bands[i][1]), format='png')

	# plot the channel ranking
	pylab.figure()
	pylab.barh(range(len(de),0,-1), de)
	pylab.xlim((0, 1.2*max(de)))
	pylab.yticks(scipy.arange(len(de)+0.5,0.5,-1), array(cnames)[id])
	pylab.subplots_adjust(left=0.4, right=0.9, top=0.9, bottom=0.1)
	pylab.xlabel('Residual error improvement')
	pylab.grid()
	pylab.ylim((len(de)-20, len(de)+1))
	savefig('rank_%d_%d.png' % (bands[i][0], bands[i][1]), format='png')
