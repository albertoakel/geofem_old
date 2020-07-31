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
xx = np.linspace(-25600.0,25600.0,8*nbcx)
yy = 400*np.sin((2*xx*np.pi)/10000)
#yy = 160*np.sin((2*xx*np.pi-3000)/100000)



pts = np.c_[mkvc(xx), mkvc(yy)]
M = refine_tree_xyz(M, pts, octree_levels=[2, 2], method='radial', finalize=False
    )





M.finalize()  # Must finalize tree mesh before use

s=np.zeros(M.nC) + 0

k=len(M.gridCC)
xx = np.linspace(-25600.0,25600.0,k)
yy = 100*np.sin((2*xx*np.pi)/10000)

for i in range(0,len(yy),1):
#    print('y',[i,yy[i]])
    s[(M.gridCC[:,1]  < yy) ] = 100






fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
M.plotImage(s,grid=True, ax=ax)
ax.set_title('teste')

#nl,nc,p
#ax1 = plt.subplot(211)
#M.plotImage(s, grid=True, ax=ax1)
#ax1.set_title('a')
##
#ax2 = plt.subplot(223)
#ax2.margins(2, 2)           # Values >0.0 zoom out
#M.plotImage(s, grid=True, ax=ax2)
#ax2.set_title('b')
#plt.xlim(-1000,1000)
#plt.ylim(-500,500)
##
#ax3 = plt.subplot(224)
#ax3.margins(x=0, y=-0.25)   # Values in (-0.5, 0.0) zooms in to center
#M.plotImage(s, grid=True, ax=ax3)
#ax3.set_title('c')
#plt.xlim(-300,300)
#plt.ylim(-300,300)






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