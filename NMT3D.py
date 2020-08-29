#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 13:03:20 2020

@author: akel
controi modelo e malha MT3D(input,opt='string')
input --> dicionário contendo informações da simulação,modelo geologico e meshgrid.
          veja a função loadmigeo para maiores informações.
          
opt ---> Variavel string que define o tipo de malha usada 'tensor' ou 'octree'
         se não declarado opt o código roda com a malha tensor.
         
         
         

def utilizadas 
         
         easymodelbox    Cria uma estrutura 3D com uma condutividade definida
         easymodellayer  Cria camadas planas 
         layermesh       Realiza mesh nas camadas
         boxmesh         Realiza mesh nas estruturas do tipo box
               
"""

import SimPEG as simpeg
import numpy as np
from discretize.utils import mkvc, refine_tree_xyz

def modelmesh(input_var,**kwargs):
    
    op=kwargs.get('opt') #tensor
    lv=kwargs.get('level') #grau de refinemanto
    
    if lv==None:
        lv=1
        pass

    dx=input_var['dxdydz'][0]
    dy=input_var['dxdydz'][1]
    dz=input_var['dxdydz'][2]
    x_length = input_var['x']    # tamanho do dominio em x
    y_length = input_var['y']    # tamanho do dominio em y
    z_length = input_var['z']    # tamanho do dominio em z
    
#Definções de mesh   
#    # Compute number of base mesh cells required in x and y
    nbcx = 2**int(np.round(np.log(x_length/dx)/np.log(2.)))
    nbcy = 2**int(np.round(np.log(y_length/dy)/np.log(2.)))
    nbcz = 2**int(np.round(np.log(z_length/dz)/np.log(2.)))

    hx = [(dx, nbcx)]
    hy = [(dy, nbcy)]
    hz = [(dz, nbcz)]
    
#    M_tensor=simpeg.Mesh.TensorMesh([hx, hy, hz], x0='CCC')  
##Discretização via treeMesh
    
    if op ==  None :
        M = simpeg.Mesh.TreeMesh([hx, hy, hz], x0='CCC')
        layermesh(M,input_var['layer'],lv,opt='S')
        if 'box' in input_var:
            boxmesh(M,input_var['box'],lv)
            pass
        M.finalize()
        
        
    if op == 'tensor':
        M=simpeg.Mesh.TensorMesh([hx, hy,hz], x0=['C', 'C','C'])
        print(hz)
        pass   
#"Contrução" do modelo (add a condutividade )  
    sig=np.zeros(M.nC) + 1e-12 # define 
    
    #inclusão de camadas, se for o caso    
    if 'layer' in input_var:
        easymodellayer(M,sig,input_var['layer'],input_var['cond'],opt='S')
        
        pass
        
    sigBG = sig
  
   #inclusão de estruturas , se for o caso    
    if 'box' in input_var:
        easymodelbox(M,sigBG,input_var['box'])
        pass
    
#       meshpyvista=simpeg.Mesh.TensorMesh([hx, hy,[(dz, nbcz/2)]
#    ], x0=['C', 'C',-dz*nbcz/2])
    
    
#    if op=='tensor':
#        SS2=sig[(M.gridCC[:,2]>100)]
#        pass
    
    
    ##To work
    #>>add  simulação MT3D
    #>>add  plot com resultados da simulação
    
    return M,sig


"""
Funções(Def) usadas: easymodelbox,easymodellayer(),layermesh,boxmesh
  
"""


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
        S[(M.gridCC[:,0]  < x+Lx) & 
          (M.gridCC[:,0]  > x) &
          (M.gridCC[:,1]  < y+Ly) &
          (M.gridCC[:,1]  > y) &
          (M.gridCC[:,2]  < z+Lz) &
          (M.gridCC[:,2]  > z) ]  =  aim_cond
    return

#função para criar as camadas
def easymodellayer(M,S,camada,cond,**kwargs):
    op=kwargs.get('opt')
    if op=='B':
        print('Model do tipo box')
        S[M.gridCC[:,2] >= 0] = 1e-12  #cond. ar
        c=0
        for i in range(0,len(camada),1):
            c=0+c
            S[(M.gridCC[:,2]  < -c) & (M.gridCC[:,2]  >= -camada[i])]=cond[i]
            c=camada[i]
            
        S[(M.gridCC[:,2]  < -c) ]=cond[i+1]
        pass
    if op=='S':
        print('Model do tipo surface')
#        k=len(M.gridCC)
        X=M.gridCC[:,0];
#        Y=M.gridCC[:,1];
        Z=np.zeros(len(X))
        S[(M.gridCC[:,2]  < Z) ] = 1e-12  #cond. ar
        c=0
        for i in range(0,len(camada),1):
            c=0+c
            S[(M.gridCC[:,2]  < -c) & (M.gridCC[:,2]  >= -camada[i])]=cond[i]
            c=camada[i]
            
        S[(M.gridCC[:,2]  < -c) ]=cond[i+1]
        pass
        
    
    return

def layermesh(M,camada,lv,**kwargs):
    
    op=kwargs.get('opt')
    if op=='B':
        print('Mesh do tipo box')

        #z=0
        xp, yp, zp = np.meshgrid( [-np.sum(M.hx)/2, np.sum(M.hx)/2],
                                   [-np.sum(M.hy)/2,np.sum(M.hy)/2],
                                   [-0-1*M.hz[0],-0+1*M.hz[0]])
    
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)] 
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box',
                finalize=False)
        #add camadas
        for i in range(0,len(camada),1):
            xp, yp, zp = np.meshgrid( [-np.sum(M.hx)/2, np.sum(M.hx)/2],
                                       [-np.sum(M.hy)/2,np.sum(M.hy)/2], 
                                       [-camada[i]-1*M.hz[0],
                                        -camada[i]+1*M.hz[0]])
            xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)]
            M = refine_tree_xyz(M, xyz, octree_levels=[lv,lv,lv], method='box',
                                finalize=False)
        pass
    
    if op=='S':
        print('Mesh do tipo surface')
        #z=0
        #M=simpeg.Mesh.TreeMesh([[(M.hx[0],len(M.hx))], [(M.hy[0],len(M.hy))],[(M.hz[0],len(M.hz))]], x0='CCC') 
        xx = np.arange(-np.sum(M.hx)/2,np.sum(M.hx)/2,M.hx[0]) #np.sum(Msurf.hxh2) ou dg --> dá as bordas
        yy = np.arange(-np.sum(M.hx)/2,np.sum(M.hx)/2,M.hx[0])
        xx, yy = np.meshgrid(xx, yy)
        zz=np.zeros([len(xx),len(xx)])-0
        #função de superficie
        xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz)]
        M = refine_tree_xyz(M, xyz, octree_levels=[lv,lv,lv], method='surface', finalize=False)
        
        #add camadas
        for i in range(0,len(camada),1):
            zz=np.zeros([len(xx),len(xx)])-camada[i]
            xyz = np.c_[mkvc(xx), mkvc(yy),mkvc(zz)]
            M = refine_tree_xyz(M, xyz, octree_levels=[lv,lv,lv], method='surface', finalize=False)
        pass
            
    return

def boxmesh(M,box,lv):
    lv=1
    n_box=len(box)/7
    for i in range(0,int(n_box),1):
        x1=box[0+int(i*7)]
        x2=x1+box[1+int(i*7)]
        y1=box[2+int(i*7)]
        y2=y1+box[3+int(i*7)]
        z1=box[4+int(i*7)]
        z2=z1+box[5+int(i*7)]
    #plano1 XY-ztop
        xp, yp, zp = np.meshgrid( [x1, x2],[y1,y2], [z1-M.hz[0],z1+M.hz[0]])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)]  
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    #plano2 XY-zboton
        xp, yp, zp = np.meshgrid( [x1,x2],[y1,y2], [z2-M.hz[0],z2+M.hz[0]])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)] 
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    #plano3 XZ-yleft
        xp, yp, zp = np.meshgrid( [x1-2*M.hx[0],x1+2*M.hx[0]],[y1,y2], [z2,z1])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)] 
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    #plano4 XZ-yrigth
        xp, yp, zp = np.meshgrid( [x2-2*M.hx[0],x2+2*M.hx[0]],[y1,y2], [z2,z1])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)]  # mkvc creates vectors
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    #plano5 YZ-Xleft
        xp, yp, zp = np.meshgrid( [x1,x2],[y1-2*M.hy[0],y1+2*M.hy[0]], [z2,z1])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)]  # mkvc creates vectors
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    #plano5 YZ-Xrigth
        xp, yp, zp = np.meshgrid( [x1,x2],[y2-2*M.hy[0],y2+2*M.hy[0]], [z2,z1])
        xyz = np.c_[mkvc(xp), mkvc(yp),mkvc(zp)]  # mkvc creates vectors
        M = refine_tree_xyz(
                M, xyz, octree_levels=[lv,lv,lv], method='box', finalize=False
                )
    return


def pyvista_view(input_var):
    import pyvista as pv
    
    dx=input_var['dxdydz'][0]
    dy=input_var['dxdydz'][1]
    dz=input_var['dxdydz'][2]
    x_length = input_var['x']    # tamanho do dominio em x
    y_length = input_var['y']    # tamanho do dominio em y
    z_length = input_var['z']    # tamanho do dominio em z
    
    #Definções de mesh   
   # Compute number of base mesh cells required in x and y
    nbcx = 2**int(np.round(np.log(x_length/dx)/np.log(2.)))
    nbcy = 2**int(np.round(np.log(y_length/dy)/np.log(2.)))
    nbcz = 2**int(np.round(np.log(z_length/dz)/np.log(2.)))

    hx = [(dx, nbcx)]
    hy = [(dy, nbcy)]
    hz = [(dz, nbcz/2)]
    
    M=simpeg.Mesh.TensorMesh([hx, hy,hz], x0=['C', 'C', -dz*nbcz/2])
    
    sig=np.zeros(M.nC) + 1e-12 # define 
    
    #inclusão de camadas, se for o caso    
    if 'layer' in input_var:
        easymodellayer(M,sig,input_var['layer'],input_var['cond'],opt='S')
        
        pass
        
    sigBG = sig
  
   #inclusão de estruturas , se for o caso    
    if 'box' in input_var:
        easymodelbox(M,sigBG,input_var['box'])
        pass
    
    models = {'res': sig}
    dataset = M.toVTK(models)
    p = pv.Plotter(notebook=0)
    p.show_grid(location='outer')
    
    p.add_mesh(dataset.slice('x'), opacity=0.75, name='x-slice')
    p.add_mesh(dataset.slice('y'), opacity=0.75, name='y-slice')
    p.add_mesh(dataset.slice('z'), opacity=0.75, name='z-slice')
    p.add_mesh(dataset.threshold([0.1,0.1]), name='vol')
    p.show() 
    
#    clim = np.abs(np.log10([np.nanmin(models)]))
#    dparams = {'rng': clim, 'cmap': 'viridis', 'show_edges': False}
#    p.add_mesh(dataset.slice('x'), opacity=0.75, name='x-slice',**dparams)
#    p.add_mesh(dataset.slice('y'), opacity=0.75, name='y-slice',**dparams)
#    p.add_mesh(dataset.slice('z'), opacity=0.75, name='z-slice',**dparams)
#    p.add_mesh(dataset.threshold([0.09,0.11]), name='vol',**dparams)
    
  

    
    
   
    