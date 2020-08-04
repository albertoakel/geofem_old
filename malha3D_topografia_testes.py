#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 13:14:56 2020

@author: akel
"""
import SimPEG as simpeg
from SimPEG import Mesh, Utils
from discretize.utils import mkvc, refine_tree_xyz
import numpy as np
import matplotlib.pyplot as plt
import time

import os
import psutil

plt.close('all')

dh=50
X=10*1024 #50x1024
Y=10*1024
Z=16*1024


#Dimensao eixo vertical( base 2)
nbcx = 2**int(np.round(np.log(X/dh)/np.log(2.)))
nbcy = 2**int(np.round(np.log(Y/dh)/np.log(2.)))
nbcz = 2**int(np.round(np.log(Z/dh)/np.log(2.)))
# Define base mesh (domain and finest discretization)
hx = [(dh, nbcx)]
hy = [(dh, nbcy)]
hz = [(dh, nbcz)]

xx = np.arange(-6400.0,6400.0,dh)
yy = np.arange(-6400.0,6400.0,dh)
xx, yy = np.meshgrid(xx, yy)

zz = 400*np.sin(2*np.pi*xx/8000) +400*np.sin(2*np.pi*yy/4000) +1 #superficie definida


M=simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')
#
xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz)]
M = refine_tree_xyz(M, xyz, octree_levels=[5,5,5], method='surface', finalize=False)

##plano de camada
M.finalize()

s=np.zeros(M.nC) + 10

k=len(M.gridCC)
print('k-->',k)

X=M.gridCC[:,0];
Y=M.gridCC[:,1];

Z=400*np.sin(2*np.pi*X/8000) +400*np.sin(2*np.pi*Y/4000)
#Z=np.zeros(len(X))

s[(M.gridCC[:,2]  < Z) ] = 100


fig, a0 = plt.subplots()
fig.canvas.set_window_title('mesh 3D')
M.plotGrid(showIt=True)


fig, a1 = plt.subplots()
fig.canvas.set_window_title('Slice Y')
M.plotSlice(np.log10(s),grid=True, normal='y',ax=a1)

fig, a2 = plt.subplots()
fig.canvas.set_window_title('Slice X')
M.plotSlice(np.log10(s),grid=True, normal='x',ax=a2)



fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(xx, yy, zz,cmap=cm.coolwarm,linewidth=0, antialiased=False)

#Customize the z axis,
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

print("\n the mesh has {} cells".format(M.nC))
print(M)
