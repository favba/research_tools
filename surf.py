#!/usr/bin/env python

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import sys

data = np.loadtxt(sys.argv[1])

Z = data[:,0]
X = data[:,1]
E = data[:,2]

nx = (X == X[0]).sum()
nz = (Z == Z[0]).sum()
s = (nx,nz)

Z = Z.reshape(s)
X = X.reshape(s)
E = E.reshape(s)

fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(X, Z, E, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax = fig.add_subplot(111)
#surf = ax.pcolormesh(X, Z, E, cmap=cm.coolwarm)
#surf = ax.contourf(X, Z, E)
surf = ax.imshow(E)
ax.set_aspect(nz/nx)
fig.colorbar(surf)
plt.show()