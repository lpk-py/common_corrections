'''
Goal: Create dataset1_wri for testing purposes
'''

import numpy as np
from obspy.core.util import locations2degrees
import os
import sys
import utils

#-----------------INPUT--------------------
phase = 'P'
filename = 'dataset1_wri_example'
#------------------------------------------

if os.path.exists(os.path.join('.', filename)):
    sys.exit('ERROR: %s exists' %(filename))

# Example 1:
#source = [[90., 0., 100.]]
#receiver = [[30., 0., 0.]]

# Example tomo:
# used for P crustal correction map
#source = []
#receiver = []
#counter = 1 
#counter_limit = 100
#for i in range(-90, 90, 20):
#    for j in range(-179, 180, 20):
#        rela = 0.; relo = 0.
#        #flag = False
#        source.append([i, j, 10.0])
#        while not 40.0 < locations2degrees(i, j, rela, relo) < 80.0:
#            #print '.',
#            #flag = True
#            rela = np.random.uniform(-89, 89)
#            relo = np.random.uniform(-179, 179)
#        #if flag: print '\n'
#        receiver.append([rela, relo, 0.0])
#        counter += 1
#    print counter,
#    if counter > counter_limit:
#        print "# receiver/source pairs: %s - %s" %(len(receiver), len(source))
#        break

# Example 2:
# used for Pdiff crustal correction map
#source = []
#receiver = []
#for i in range(-89, 90):
#    for j in range(-179, 180):
#        rela = 0.; relo = 0.
#        #flag = False
#        source.append([i, j, 0.0])
#        while not 100.0 < locations2degrees(i, j, rela, relo) < 150.0:
#            #print '.',
#            #flag = True
#            rela = np.random.uniform(-89, 89)
#            relo = np.random.uniform(-179, 179)
#        #if flag: print '\n'
#        receiver.append([rela, relo, 0.0])

# Example 3
#source = []
#receiver = []
#for i in range(-10, -90, -1):
#    receiver.append([i, 10.0, 0.0])
#    source.append([90.0, 10.0, 0.0])
#source = []
#receiver = []
#for i in range(60, -10, -1):
#    receiver.append([i, 10.0, 0.0])
#    source.append([90.0, 10.0, 0.0])

# Example 4
source = []
receiver = []
for i in range(32, 90, 1):
    receiver.append([0.0, i, 0.0])
    source.append([0.0, 0.0, 0.0])
#for i in range(100, 180, 1):
#    receiver.append([0.0, i, 0.0])
#    source.append([0.0, 0.0, 0.0])

# Example 5
#source = []
#receiver = []
##for i in range(30, 90, 1):
##    receiver.append([i, 0.0, 0.0])
##    source.append([0.0, 0.0, 0.0])
#for i in range(80, 0, -1):
#    receiver.append([i, 180.0, 0.0])
#    source.append([0.0, 0.0, 0.0])


# Creating the header (wave definition + filters)
comp_header = utils.dataset_header(phase, filename)
comp_header = comp_header.split('\n')
comp_file = []
for i in range(len(comp_header)-1, -1, -1):
    if not comp_header[i].isspace() and \
            not comp_header[i] == '':
        comp_file.append(comp_header[i])
comp_file.reverse()

# Creating the required source-receiver pairs
srcrcv = utils.source_receiver(source, receiver)
for i in range(len(srcrcv)):
    comp_file.append(srcrcv[i])

fio = open(os.path.join('.', filename), 'w')
for i in range(len(comp_file)):
    fio.writelines(comp_file[i] + '\n')
fio.close()

fio = open(os.path.join('.', 'in.'+filename), 'w')
fio.writelines('IASP91\n')
fio.writelines('1 20\n')
fio.writelines(filename + '\n')
if phase.lower() == 'pdiff':
    fio.writelines('Pdef_Pdiff')
elif phase.lower() == 'p':
    fio.writelines('Pdef_P')
else:
    sys.exit('ERROR: unrecognized phase: %s' %(phase))
fio.close()
