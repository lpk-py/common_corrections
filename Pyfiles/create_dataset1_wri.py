'''
Goal: Create dataset1_wri for testing purposes
'''

import os
import sys
import utils

#-----------------INPUT--------------------
phase = 'Pdiff'
filename = 'dataset1_wri_example'
#------------------------------------------

# Example 1:
#source = [[89., 30., 0.]]
#receiver = [[-29., 30., 0.]]

# Example 2:
source = []
receiver = []
for i in range(-89, 90):
    for j in range(-179, 179):
        source.append([i, j, 0])
        receiver.append([0, 0, 0])

if os.path.exists(os.path.join('.', filename)):
    sys.exit('ERROR: %s exists' %(filename))

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
