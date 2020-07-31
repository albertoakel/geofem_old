#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 17:43:22 2019

@author: ellen
"""

from SimPEG import Mesh, Utils
from discretize.utils import mkvc, refine_tree_xyz
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


style_list = ['default', 'classic'] + sorted(
        style for style in plt.style.available if style != 'classic')

plt.close('all')

# sphinx_gallery_thumbnail_number = 4

###############################################
# Basic Example
# -------------
#
# Here we demonstrate the basic two step process for creating a 2D tree mesh
# (QuadTree mesh). The region of highest discretization if defined within a
# rectangular box. We use the keyword argument *octree_levels* to define the
# rate of cell width increase outside the box.
#

#size dh x nbc(base2)
#Dimensoes na horizontal
dh = 50   # minimum cell width (base mesh cell width)

#Dimensao eixo vertical( base 2)
nbcx =1024# number of base mesh cells in x
nbcy =1024
nbcz =128
# Define base mesh (domain and finest discretization)
hx = dh*np.ones(nbcx)
hy = dh*np.ones(nbcy)
hz = dh*np.ones(nbcz)


M =  Mesh.TreeMesh([hx,hy]);
M2=Mesh.TensorMesh([hx,hy]);


M.x0 = np.r_[-(nbcx*dh)/2,-15000]


#definir a camada 

#xp, yp = np.meshgrid( [-(nbcx*dh)/2, (nbcx*dh)/2], [-0., -1.]) #layer
#xy = np.c_[mkvc(xp), mkvc(yp)]  # mkvc creates vectors
##
## Discretize to finest cell size within rectangular box
#M = refine_tree_xyz(
#    M, xy, octree_levels=[2,2], method='box', finalize=False
#    )


#definir camada como função
xx = np.linspace(-25600.0,25600.0,10*nbcx)
yy = 400*np.sin((2*xx*np.pi)/10000)
#yy = 160*np.sin((2*xx*np.pi-3000)/100000)



pts = np.c_[mkvc(xx), mkvc(yy)]
M = refine_tree_xyz(M, pts, octree_levels=[2, 2], method='radial', finalize=False
    )


M.finalize()  # Must finalize tree mesh before use

s=np.zeros(M.nC) + 0

k=len(M.gridCC)
xx = np.linspace(-25600.0,25600.0,k)
yy =400*np.sin((2*xx*np.pi)/10000)

for i in range(0,len(yy),1):
#    print('y',[i,yy[i]])
    s[(M.gridCC[:,1]  < yy) ] = 100






fig1 = plt.figure(figsize=(6, 6))
fig1.canvas.set_window_title('Topografia grid 2D')

ax = fig1.add_subplot(111)
M.plotImage(s,grid=True, ax=ax)
ax.set_title('Topografia - Teste octree 2D')










#nl,nc,p
fig2 = plt.figure()
fig2.canvas.set_window_title('Topografia grid 2D - Zoom')

ax1 = fig2.add_subplot(221)
M.plotImage(s, grid=True, ax=ax1)
ax1.set_title('a')
plt.xlim(-20000,-10000)
plt.ylim(-2000,2000)
#
ax2 = fig2.add_subplot(222)
M.plotImage(s, grid=True, ax=ax2)
ax2.set_title('b')
plt.xlim(-10000,0000)
plt.ylim(-2000,2000)
#
ax3 = fig2.add_subplot(223)
M.plotImage(s, grid=True, ax=ax3)
ax3.set_title('c')
plt.xlim(-20000,-10000)
plt.ylim(-1000,1000)

ax4 = fig2.add_subplot(224)
M.plotImage(s, grid=True, ax=ax4)
ax3.set_title('d')
plt.xlim(-5000,0000)
plt.ylim(-600,600)







#
#
#fig, ax = plt.subplots(figsize=(6,6))
#M.plotImage(s,grid=True)
#plt.xlim(-2000,2000)
#plt.ylim(-400,400)



#M.plotGrid(showIt=True)
#ax = plt.gca()
##ax.invert_yaxis()
#plt.show()


print(M)
#print("Aqui!")
#mesh.plotGrid(showIt=True)