"""
Plot the ray which has been traced!
Not completely implemented!

input: output of raydata program
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

filename = os.path.abspath(sys.argv[1])
ICB = 1217.9998
CMB = 3482.0
Rad = 6371.0

# Opening the input file
fio = open(filename, 'r')
fi = fio.readlines()

line_to_start = int(raw_input('Enter the line number: (starting from 1)'))
#print int(fi[line_to_start-1])
number_lines_2_read = int(fi[line_to_start-1].split()[-1])

# r: radius, delta: epicentral distance (radians)
r = []
delta = []
for i in range(line_to_start, line_to_start+number_lines_2_read):
    fi_tmp = fi[i].split()
    if len(fi_tmp) == 7:
        r.append(float(fi_tmp[0]))
        #delta.append(float(fi_tmp[2])*180./np.pi)
        delta.append(float(fi_tmp[2]))

ax = plt.subplot(111, polar=True)
ax.plot(delta, r, color='r', linewidth=3)

# DIRTY! plot ICB and CMB
ax.plot(np.linspace(0, 360, 10000), np.ones(10000)*ICB, lw=1.0, color='k')
ax.plot(np.linspace(0, 360, 10000), np.ones(10000)*CMB, lw=1.0, color='k')

# change the theta direction to clockwise!
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')

ax.set_thetagrids([0, 45, 90, 135, 180, 225, 270, 315], size='x-large', weight='bold')
ax.set_rgrids([ICB, CMB], labels=['ICB', 'CMB'], size='x-large', weight='bold')
ax.set_rmax(Rad)

plt.show()
