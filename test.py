from nonna_functions import *

gps1,gps2 = nonna_readsegmentfile('segments_h1.txt')
i = 12


d = nonna_get_data_from_disk('H1:CAL-DELTAL_EXTERNAL_DQ', gps1[i], gps2[i] - gps1[i], outfs=1024, verbose=False)
b1 = nonna_blrms(d, 10, 20, 1024, 1)
