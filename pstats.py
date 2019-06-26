#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.genfromtxt("Stats.txt",skip_header=1,delimiter=',')
f = open("Stats.txt")
headers = f.readline().rstrip().split(',')
f.close()

d = {headers[i]: data[:,i] for i in range(data.shape[1])}

plt.rcParams['legend.loc'] = 'best'
plt.rc('font', family='serif')

for el in sys.argv[1:]:
    plt.plot(d['time'],d[el],label=el)

plt.xlabel("Time (s)")

plt.legend()
plt.grid()
plt.show()