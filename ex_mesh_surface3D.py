#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 18:20:14 2020

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

##surface method

Msurf=simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')

#Msurf=simpeg.Mesh.TensorMesh([hx, hy, hz], x0='CCC') #tensormesh(caso opte, comentar linhas 54a58)


xx = np.arange(-np.sum(Msurf.hx)/2,np.sum(Msurf.hx)/2,dh)
yy = np.arange(-np.sum(Msurf.hx)/2,np.sum(Msurf.hx)/2,dh)
xx, yy = np.meshgrid(xx, yy)

zz=np.zeros([len(xx),len(xx)])  #superficie 1)

zz2=0.0001*xx**2+0.0001*yy**2-4200 #superficie 2)


#plano de camada
#caso opte por tensor mesh, comentar as linhas(55,56,57,58,59)
xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz)]
Msurf = refine_tree_xyz(Msurf, xyz, octree_levels=[2,2,2], method='surface', finalize=False)
xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz2)]
Msurf = refine_tree_xyz(Msurf, xyz, octree_levels=[3,3,4], method='surface', finalize=False)
Msurf.finalize()
s2=np.zeros(Msurf.nC) + 10


k=len(Msurf.gridCC)
##print('k-->',k)
##
X=Msurf.gridCC[:,0];
Y=Msurf.gridCC[:,1];

Z=np.zeros(len(X)) #superficie 1)

#Z=top_func(X,Y)
Z2=(0.0001*X**2+0.0001*Y**2)-4200 #superficie 2)
s2[(Msurf.gridCC[:,2]  < Z) ] = 10000
#

s2[(Msurf.gridCC[:,2]  < Z) & (Msurf.gridCC[:,2]  >= Z2)]=100
#

fig, a3 = plt.subplots()
fig.canvas.set_window_title('Slice Y - surface Method')
Msurf.plotSlice(np.log10(s2),grid=True, normal='y',ax=a3)

fig, a4 = plt.subplots()
fig.canvas.set_window_title('Slice X - surface Method')
Msurf.plotSlice(np.log10(s2),grid=True, normal='x',ax=a4)

print(Msurf)


