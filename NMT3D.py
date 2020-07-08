#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:03:20 2020

@author: akel
# Realiza modelagem MT3D
"""


import SimPEG as simpeg
import numpy as np
import matplotlib.pyplot as plt
import sys
from colorama import Fore, Back, Style 
from discretize.utils import mkvc, refine_tree_xyz

def MT3D(input_var):
    
    dx=input_var['dxdydz'][0]
    dy=input_var['dxdydz'][1]
    dz=input_var['dxdydz'][2]
    x_length = input_var['x']    # tamanho do dominio em x
    y_length = input_var['y']    # tamanho do dominio em y
    z_length = input_var['z']    # tamanho do dominio em z
    
#===============Definções de mesh==========================================   
    # Compute number of base mesh cells required in x and y
    nbcx = 2**int(np.round(np.log(x_length/dx)/np.log(2.)))
    nbcy = 2**int(np.round(np.log(y_length/dy)/np.log(2.)))
    nbcz = 2**int(np.round(np.log(z_length/dz)/np.log(2.)))
    
    # Define the base mesh
    hx = [(dx, nbcx)]
    hy = [(dy, nbcy)]
    hz = [(dz, nbcz)]
    
  
    M = simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')
    [xx, yy] = np.meshgrid(M.vectorNx, M.vectorNy)
    [xx, yy,zz] = np.meshgrid([-1000,1000], [-1000,1000],[-2300,-2100])
    zz = 1.*(xx**2 + yy**2)
    
    
    
    pts = np.c_[mkvc(xx), mkvc(yy), mkvc(zz)]
    M = refine_tree_xyz(
        M, pts, octree_levels=[2, 2], method='box', finalize=False
                   )
    
    
    xp, yp, zp = np.meshgrid([-3200., 3200.], [-1000., 1000.], [-1000., -1800.])
    xyz = np.c_[mkvc(xp), mkvc(yp), mkvc(zp)]

    M = refine_tree_xyz(
    M, xyz, octree_levels=[6,6], method='box', finalize=False)
    
    # M.finalize() 
    # M.plotGrid(showIt=True)
    
    def refine(cell):
         if np.sqrt(((np.r_[cell.center]-0.5)**2).sum()) < 0.5:
             return 4
         return 5
    M.refine(refine)
    M.finalize() 
    
   
#==============FIM==========================================   
   
    sig=np.zeros(M.nC) + 1e-12 # define 
    
    #inclusão de camadas, se for o caso    
    if 'layer' in input_var:
        easymodellayer(M,sig,input_var['layer'],input_var['cond'])
        pass
        
    sigBG = sig
  
   #inclusão de estruturas , se for o caso    

    if 'box' in input_var:
        easymodelbox(M,sigBG,input_var['box'])
        pass
    

#    mesh1d = simpeg.Mesh.TensorMesh([M.hz], np.array([M.x0[2]]))
#    sigBG1d = np.zeros(mesh1d.nC) + conds[1]

   
    fig,axes = plt.subplots(num=1,clear=True)
    M.plotSlice(np.log(sig), grid=True, normal='y',ax=axes)
# #    plt.ylim(-6000,2000)
# #    plt.colorbar(sig)
    plt.show()
    
    

    return M,sig


#função para criar os box
def easymodelbox(M,S,B):
    n_box=len(B)/7
    for i in range(0,int(n_box),1):
        x=B[0+int(i*7)]
        Lx=B[1+int(i*7)]
        y=B[2+int(i*7)]
        Ly=B[3+int(i*7)]
        z=B[4+int(i*7)]
        Lz=B[5+int(i*7)]
        aim_cond=B[6+int(i*7)]
        S[(M.gridCC[:,0]  < x+Lx) & (M.gridCC[:,0]  > x) & (M.gridCC[:,1]  < y+Ly) & (M.gridCC[:,1]  > y) & (M.gridCC[:,2]  < z+Lz) & (M.gridCC[:,2]  > z) ]  =  aim_cond
    return

#função para criar as camadas
def easymodellayer(M,S,thi,cond):
    # while True:
    #     if len(cond)!=len(thi):
    #         print()
    #         print(Fore.RED +'++++ERRO NA DIMENSÃO DAS CONDUTIVIDADES X ESPESSURAS++++')
    #         sys.exit() 
    #         break
    #     else:
    n_layers=len(thi)
    S[M.gridCC[:,2] >= 0] =   1e-12  #cond. ar
    c=0
    for i in range(0,n_layers,1):
        c=0+c
        S[(M.gridCC[:,2]  < -c) & (M.gridCC[:,2]  >= -thi[i])]=cond[i]
        c=thi[i]
    #define limite do grid igual a camada anterior    
    S[(M.gridCC[:,2]  < -c) & (M.gridCC[:,2]  >= -M.gridCC[-1,2])]=cond[i]
    return