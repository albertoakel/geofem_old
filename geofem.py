#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:47:10 2020

@author: akel
Programa teste com entradas para modelagem MT3D

"""
#import matplotlib.pyplot as plt

import argparse
import NMT3D as model
import loadmigeo as load #contem arquivos de input
#import numpy as np
print('Running geofem test')

#plt.close('all')

##usando o spyder! 
#out=load.inputfiles('MT3D_input.in')


#usando o terminal!
parser = argparse.ArgumentParser()
parser.add_argument("inp", help="input file",type=str)
ARG = parser.parse_args()
out=load.inputfiles(ARG.inp)
#print(out['box'])
#print('numero de box:', len(out['box'])/7)
Me,cond=model.MT3D(out)
