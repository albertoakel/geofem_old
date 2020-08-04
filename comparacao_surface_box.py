#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:56:10 2020

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

#out=load.inputfiles('MT3D_input_ex3.in')
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

##plano_camada--Box function
Mbox = simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')
xp, yp, zp = np.meshgrid( [-np.sum(M.hx)/2, np.sum(M.hx)/2],[-np.sum(M.hy)/2,np.sum(M.hy)/2], [-0-1*M.hz[0],-0+1*M.hz[0]])
xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)] 
Mbox = refine_tree_xyz(Mbox, xyz, octree_levels=[1,1,1], method='box', finalize=False )
Mbox.finalize()

###plano de camada
s=np.zeros(Mbox.nC) + 10
s[(Mbox.gridCC[:,2]  < 0) ] = 100

fig, a1 = plt.subplots()
fig.canvas.set_window_title('Slice Y - Box Method')
Mbox.plotSlice(np.log10(s),grid=True, normal='y',ax=a1)

fig, a2 = plt.subplots()
fig.canvas.set_window_title('Slice X - Box Method')
Mbox.plotSlice(np.log10(s),grid=True, normal='x',ax=a2)

print(Mbox)


##surface method

Msurf=simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')

xx = np.arange(-np.sum(Msurf.hx)/2,np.sum(Msurf.hx)/2,dh) #np.sum(Msurf.hxh2) ou dg --> dá as bordas
yy = np.arange(-np.sum(Msurf.hx)/2,np.sum(Msurf.hx)/2,dh)
xx, yy = np.meshgrid(xx, yy)
zz=np.zeros([len(xx),len(xx)])+0  #função de superficie

xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz)]
Msurf = refine_tree_xyz(Msurf, xyz, octree_levels=[2,2,2], method='surface', finalize=False)



#plano de camada
Msurf.finalize()
s2=np.zeros(Msurf.nC) + 10


k=len(Msurf.gridCC)
##print('k-->',k)
##
X=Msurf.gridCC[:,0];
Y=Msurf.gridCC[:,1];

Z=np.zeros(len(X))

#Z=top_func(X,Y)
s2[(Msurf.gridCC[:,2]  < Z) ] = 100
#

fig, a3 = plt.subplots()
fig.canvas.set_window_title('Slice Y - surface Method')
Msurf.plotSlice(np.log10(s2),grid=True, normal='y',ax=a3)

fig, a4 = plt.subplots()
fig.canvas.set_window_title('Slice X - surface Method')
Msurf.plotSlice(np.log10(s2),grid=True, normal='x',ax=a4)

print(Msurf)