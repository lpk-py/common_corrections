"""
Plot the difference between ellipticity corrections (two methods)
"""

import matplotlib.pyplot as plt
import os
import sys

filename = os.path.abspath(sys.argv[1])

fio = open(filename, 'r')
fi = fio.readlines()

diff = []
for i in range(len(fi)):
    diff.append(float(fi[i].split(',')[2]))

plt.hlines(0.1, 0, len(diff), lw=1, linestyle='dashed')
plt.plot(diff, lw=3)
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')
plt.ylabel('abs(getelcor - ellip)', size='x-large', weight='bold')
plt.show()
