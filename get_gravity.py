#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 11:30:40 2018

@author: Weigang Liang

Generates the Bouguer gravity and gravity gradient maps of the Moon
"""
import pyshtools

import matplotlib.pyplot as plt
import pdb 
import numpy as np
import scipy.io as io
import sys
import time
import copy
from scipy.interpolate import UnivariateSpline
    
#Jeff: the for loop for the Bouguer elevation part starts at line 161
start_time = time.time()
# your code

rho_c = 2550
rho_d = -450
a = 1738000
gm = 4902.80011526323e9
l = 1200#500
nmax = 14
shape = pyshtools.SHCoeffs.from_file('MoonTopo2600p.shape.sh',lmax=l)
omega = 0
f=0
near = 1735000
far = 1742500
dh = far-near
mass = 7.3457892e22


#original
lmax_calc = 500#450#500#l
  
shape0 = copy.deepcopy(shape)
shape_grid0 = shape0.expand(grid='DH2',lmax=1000)
interface0 = shape_grid0.to_array()

temp = shape.to_array()
temp[:,31:,:] = 0
interface2 = pyshtools.expand.MakeGridDH(temp,lmax=lmax_calc,sampling=2)


rho = rho_d+np.zeros(interface2.shape)

size = interface2.shape

b_step = 100


pass_num =  5#low-pass spherical harmonic degree cutoff

added = 1000#elevation above surface for gravity calculation

mmin = np.min(interface2)
mmax = np.max(interface2)


temp0 = shape0.to_array()
interface0 = pyshtools.expand.MakeGridDH(temp0,lmax=lmax_calc,sampling=2)#1000

rho = 600 #Bouguer correction density, in kg/m^3

cilm_inter = pyshtools.expand.SHExpandDH(interface0,lmax_calc=lmax_calc,sampling=2)
clm_inter = pyshtools.SHCoeffs.from_array(cilm_inter,lmax=lmax_calc)
    
bc = pyshtools.SHGravCoeffs.from_shape(clm_inter, rho=rho, gm=gm, lmax=lmax_calc,nmax=9)#3100
bc2 = bc.change_ref(r0=a)

bg = pyshtools.SHGravCoeffs.from_file('gggrx_1200a_bouguer_sha.tab.txt', errors=True, header_units='km')  
mod_clm0 = bg.change_ref(r0=a,gm=gm,lmax=lmax_calc)


mod_clm = mod_clm0 - bc2
mod_cilm = mod_clm.to_array()
mod_cilm[0][0][0] = 1

mod_cilm[:,1:pass_num,:] = 0
temp_clm = mod_cilm[:,:lmax_calc+1,:lmax_calc+1]


potential = np.zeros([interface2.shape[0],interface2.shape[1],np.arange(mmin,mmax,b_step).shape[0]])
ind = 0
new_potential3 = np.zeros(interface2.shape)
bouguer3 = np.zeros(interface2.shape)
bg4 = pyshtools.SHGravCoeffs.from_array(temp_clm,gm,a,omega=0,lmax=lmax_calc)

# Gravity on a surface algorithm
for i in np.arange(mmin,mmax,b_step):
    bg3 = bg4.change_ref(r0=i+added,gm=gm,lmax=lmax_calc)#bg4
    temp = bg3.to_array()

    potential2 = pyshtools.expand.MakeGridDH(temp,lmax=lmax_calc,sampling=2)
    new_potential3[(interface2>=i)*(interface2<=i+b_step)] = potential2[(interface2>=i)*(interface2<=i+b_step)]
    r,t,p,gravb,pot = pyshtools.gravmag.MakeGravGridDH(temp_clm,gm,i+added,lmax=lmax_calc,a=i,f=0,omega=0,normal_gravity=1)
    bouguer3[(interface2>=i)*(interface2<=i+b_step)] = gravb[(interface2>=i)*(interface2<=i+b_step)]
    if ind == 0:
        new_potential3[interface2<=i] = potential2[interface2<=i]
        bouguer3[interface2<=i] = gravb[interface2<=i]
    elif ind == np.arange(mmin,mmax,b_step).shape[0]-1:
        new_potential3[interface2>=i] = potential2[interface2>=i]
        bouguer3[interface2>=i] = gravb[interface2>=i]
    ind+=1

temp_clm1 = pyshtools.expand.SHExpandDH(new_potential3,lmax_calc=lmax_calc,sampling=2)
temp_clm1[:,:5,:] = 0

temp_clm1[0][0][0] = 1

# Accounting for the nearside-farside crustal thickness variations
from_matlab = {}
from_matlab = io.loadmat('rscl.mat')
rscl = np.squeeze(from_matlab['rscl'])

# cosine tapering at the end of the spherical harmonic window to avoid aliasing
L = np.arange(0,lmax_calc+1)
L1 = 400
L2 = lmax_calc

taper0 = (1-np.cos((L2-L)/(L2-L1)*np.pi))/2
taper0[30:110] = rscl[30:110]**0.5
taper = np.tile(taper0,(L2+1,1)).transpose()

taperf = np.zeros((2,L2+1,L2+1))
taperf[0,:,:] = taper
taperf[1,:,:] = taper

temp_clm1[:,30:110,:] = taperf[:,30:110,:]*temp_clm1[:,30:110,:]
temp_clm1[:,L1:,:] = taperf[:,L1:,:]*temp_clm1[:,L1:,:]

#
r,t,p,grav0,pot = pyshtools.gravmag.MakeGravGridDH(temp_clm1,gm,a,lmax=lmax_calc,a=a,f=0,omega=0,normal_gravity=1)


clm0 = pyshtools.SHGravCoeffs.from_array(temp_clm1,gm,a,omega=0,lmax=lmax_calc)
tensor0 = clm0.tensor(a=a,f=0,lmax=lmax_calc)

tensor0.compute_eigh()
eig0 = tensor0.eighh.data


eig22 = np.concatenate((eig0[:,size[1]//4:3*size[1]//4],eig0[:,3*size[1]//4:],eig0[:,:size[1]//4]),axis=1)
grav22 = np.concatenate((grav0[:,size[1]//4:3*size[1]//4],grav0[:,3*size[1]//4:],grav0[:,:size[1]//4]),axis=1)

sys.exit()
#Spherical harmonic localization and calculating the power spectra

tapersn_og, eigenvalues_nog = pyshtools.spectralanalysis.SHReturnTapersMap(maskn0, 40, sampling=2,ntapers=3)#4,3
mtsen_og, sdnog = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersn_og)#temp_clm1

tapersf_og, eigenvalues_fog = pyshtools.spectralanalysis.SHReturnTapersMap(maskf0, 40, sampling=2,ntapers=2)#23,22
mtsef_og, sdfog = pyshtools.spectralanalysis.SHMultiTaperMaskSE(temp_clm1, tapersf_og)


elapsed_time = time.time() - start_time
print(elapsed_time)

