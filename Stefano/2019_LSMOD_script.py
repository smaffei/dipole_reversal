#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 13:50:57 2020

@author: earsmaf
"""


import numpy as np
from scipy import interpolate
from scipy.interpolate import SmoothSphereBivariateSpline as SSBS
import os
import math
import sys
import csv
import glob

sys.path.append('./SpecialFunctions/SphericalHarmonics/')
import SH_library

import subs

import matplotlib.tri as tri

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
from matplotlib import gridspec



# inputs for animation
animation_flag = 0 # 0: no animation; 1: create animation  
yps = 10 # years per second of animation 

# other inputs
tML = 80 # Mono Lake
tL   = 220 # Laschamp
########################################################
# from Brown et al., 2018 the Mono Lake and LAschamp
# excursions are conventionally centered at 34 and 41 ka,
# respectively
#########################################################

r_cmb = 3485.0e3
r_a   = 6371.0e3

# Sulmona colatitute and longitude
colat_SUL = 47.8488
lat_SUL = 90-colat_SUL
lon_SUL   = 13.8229 
# Sidney
lat_Sid = -33.8688
colat_Sid = 90-lat_Sid
lon_Sid = 151.2093
#Quito
lat_Quito = -0.1807
colat_Quito = 90-lat_Quito
lon_Quito = -78.4675
# specific locations
lat_locations = np.array([lat_Quito, lat_SUL, lat_Sid]) # equator, sulmona, southern Hemisphere
theta_locations = -lat_locations*np.pi/180. + np.pi/2
lon_locations = np.array([lon_Quito, lon_SUL, lon_Sid])
phi_locations = lon_locations*np.pi/180.

#input for maps
nlats = 80
nlons = 100




"""
if we use cartopy (bugged!) we need to arrange lat and lon in the following way
lats = np.linspace(-np.pi / 2, np.pi / 2, nlats)
lons = np.linspace(0, 2 * np.pi, nlons)
"""
lats = np.linspace(np.pi/2,-np.pi/2,num=100)
lons = np.linspace(-np.pi,np.pi,num=100)
lons, lats = np.meshgrid(lons, lats)
theta = -lats+np.pi/2

figname = 'LSMOD1'
folder = '../models/LSMOD1/'
figfolder = '../models/LSMOD1/figures/'

line = []
coeffs_B = []
with open(folder+'GaussCoeff.txt', 'rb') as fipgp:
    while True:
        line = fipgp.readline()
        if not line:
            break
        x = [list(map(float, line.split() ))]
        coeffs_B = np.append(coeffs_B,x)
cols = len(x[0])
rows = len(coeffs_B)
coeffs_B = np.reshape(coeffs_B,(rows/cols, cols))

Lmax = -1+np.sqrt(cols) # shouold be cols-1+1
    
#####################################
# coefficients are in nanoTeslas!!!!
#####################################    

n_instants=coeffs_B.shape[0]
times = coeffs_B[:,0]   
age = (-times+1950.)/1000.
# use energy definition from Brown, Korte et al., 2018 (PNAS): eq 1 in the Methods section
ME_cmb_spectra = np.zeros((n_instants,Lmax+1))
ME_cmb_spectra[:,0] = 1.*coeffs_B[:,0]
ME_a_spectra = ME_cmb_spectra.copy()

n_coeffs=coeffs_B.shape[1]-1
Br_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Br_c = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Bt_a = np.zeros((n_instants, theta.shape[0], theta.shape[1])) 
Bp_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))

Br_locations = np.zeros((n_instants, theta_locations.size))
Bt_locations = np.zeros((n_instants, theta_locations.size))
Bp_locations = np.zeros((n_instants, theta_locations.size))



it = 0
l=1
m=0
cs='c'

for ic in range(n_coeffs):

    betalm = coeffs_B[:,ic+1]
    ME_cmb_spectra[:,l] = ME_cmb_spectra[:,l] + (l+1)* (r_a/r_cmb)**(2*l+4) * betalm**2
    ME_a_spectra[:,l] = ME_a_spectra[:,l] +  (l+1)*betalm**2
    
    if cs=='s':
        dthSH = SH_library.DthSchmidtSH(l,m,theta,lons,'s')
        dphSH = SH_library.DphSchmidtSH(l,m,theta,lons,'s')
        SH    = SH_library.SchmidtSH(l,m,theta,lons,'s')
        
        dthSH_locations = SH_library.DthSchmidtSH(l,m,theta_locations,phi_locations,'s')
        dphSH_locations = SH_library.DphSchmidtSH(l,m,theta_locations,phi_locations,'s')
        SH_locations    = SH_library.SchmidtSH(l,m,theta_locations,phi_locations,'s')
        
    else:
        dthSH = SH_library.DthSchmidtSH(l,m,theta,lons,'c')
        dphSH = SH_library.DphSchmidtSH(l,m,theta,lons,'c')
        SH    = SH_library.SchmidtSH(l,m,theta,lons,'c')
        
        dthSH_locations = SH_library.DthSchmidtSH(l,m,theta_locations,phi_locations,'c')
        dphSH_locations = SH_library.DphSchmidtSH(l,m,theta_locations,phi_locations,'c')
        SH_locations    = SH_library.SchmidtSH(l,m,theta_locations,phi_locations,'c')
        
    for it in range(n_instants):
        Br_a[it] = Br_a[it] + (l+1) * (r_a/r_a)**(l+2) * betalm[it] * SH
        Br_c[it] = Br_c[it] + (l+1) * (r_a/r_cmb)**(l+2) * betalm[it] * SH
        Bt_a[it] = Bt_a[it] - (r_a/r_a)**(l+2) * betalm[it] * dthSH
        Bp_a[it] = Bp_a[it] - (r_a/r_a)**(l+2) * betalm[it] * np.divide(dphSH,np.sin(theta))
    
        Br_locations[it] = Br_locations[it] + (l+1) * (r_a/r_a)**(l+2) * betalm[it] * SH_locations
        Bt_locations[it] = Bt_locations[it] - (r_a/r_a)**(l+2) * betalm[it] * dthSH_locations
        Bp_locations[it] = Bp_locations[it] - (r_a/r_a)**(l+2) * betalm[it] * np.divide(dphSH_locations,np.sin(theta_locations))
        
        #print("time = " +str(times[it]) + "; l = " + str(l) + "; m= " + str(m)+ "; cs= " + cs + "; beta = " + str(betalm[it]))
    if m==l:
        if cs=='c':
            cs='s'
        else:
            l=l+1
            m=0
            cs='c'
    else:
        if m==0:            
            m=m+1
            cs='c'
        else:
            if cs=='c':            
                cs='s'
            else:           
                cs='c'
                m=m+1
        
             
F_a = np.sqrt(Br_a**2 + Bt_a**2 + Bp_a**2)
F_locations = np.sqrt(Br_locations**2 + Bt_locations**2 + Bp_locations**2)

g10 = coeffs_B[:,1]
g11 = coeffs_B[:,2]
h11 = coeffs_B[:,3]
m1  = np.sqrt(g10**2 + g11**2 + h11**2)

# geomagnetic dipole
dipole_colatitude = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,m1)))
dipole_lon        =  180/np.pi * np.arctan(np.divide(g11,h11)) 



lats = np.rad2deg(lats)
lons = np.rad2deg(lons)

# find magnetic poles
inclination = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2+Bp_a**2)))
idxN = np.zeros(n_instants)
idxS = np.zeros(n_instants)
idxSAA = np.zeros(n_instants)
N_pole_colatitude = np.zeros(n_instants)
N_pole_lon        = np.zeros(n_instants)
S_pole_colatitude = np.zeros(n_instants)
S_pole_lon        = np.zeros(n_instants)
SAA_colatitude = np.zeros(n_instants)
SAA_lon        = np.zeros(n_instants)

for it in range(n_instants):
    # find poles locations
    idxN[it] = np.nanargmin(np.abs(inclination[it] - np.pi/2))
    idxS[it] = np.nanargmin(np.abs(inclination[it] + np.pi/2))
    N_pole_colatitude[it] = lats.flat[int(idxN[it])]
    N_pole_lon[it]        = lons.flat[int(idxN[it])]
    S_pole_colatitude[it] = lats.flat[int(idxS[it])]
    S_pole_lon[it]        = lons.flat[int(idxS[it])]
    # find minimum intensity (SAA)
    idxSAA[it] = np.nanargmin(F_a[it])
    SAA_colatitude[it] = lats.flat[int(idxSAA[it])]
    SAA_lon[it]        = lons.flat[int(idxSAA[it])]



############################################
# get 2019 CHAOS coefficients for comparison
############################################
CHAOS_2019 = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007.dat' 

coeffs_CHAOS=subs.read_coeffs(CHAOS_2019,1,14)

CHAOS_a_spectra = SH_library.calc_spectra(coeffs_CHAOS,r_a,r_a)
CHAOS_cmb_spectra = SH_library.calc_spectra(coeffs_CHAOS,r_a,r_cmb)




##########################################
# time-series of optimizes quantities
##########################################

# input file
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 10
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 

listing = glob.glob('../models/LSMOD1/*dat')

times_it=np.linspace(1,len(listing),len(listing))
I_opt_SUL = []
I_opt_Sid = []
I_opt_Quito = []
tilt_opt = []
g10_opt = []

# optimise for inclination at SUL
filename_I_opt_SUL = ('dIdt_SUL_LSMOD1_timeseries.txt')
try:
    I_opt_SUL = np.loadtxt(filename_I_opt_SUL)
except:
    print('no incliation file found')
    
if I_opt_SUL == []:
    for it in times_it:
        MODEL = '../models/LSMOD1/'+str(int(it))+'.dat'
        FLAG=5
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        opt_file.readline() # header
        line=opt_file.readline()
        x = [list(map(float, line.split() ))]
        opt_file.close()
        I_opt_SUL=np.append(I_opt_SUL,[it, x[0][0]])
        os.system('rm *DAT')
        
    I_opt_SUL = np.reshape(I_opt_SUL,(len(I_opt_SUL)/2, 2))     
    np.savetxt(filename_I_opt_SUL,I_opt_SUL,fmt='%6d %.10f')


# optimise for inclination at Sidney
filename_I_opt_Sid = ('dIdt_Sid_LSMOD1_timeseries.txt')
try:
    I_opt_Sid = np.loadtxt(filename_I_opt_Sid)
except:
    print('no incliation file found')
    
if I_opt_Sid == []:
    for it in times_it:
        MODEL = '../models/LSMOD1/'+str(int(it))+'.dat'
        FLAG=5
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_Sid,lon_Sid,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        opt_file.readline() # header
        line=opt_file.readline()
        x = [list(map(float, line.split() ))]
        opt_file.close()
        I_opt_Sid=np.append(I_opt_Sid,[it, x[0][0]])
        os.system('rm *DAT')
        
    I_opt_Sid = np.reshape(I_opt_Sid,(len(I_opt_Sid)/2, 2))     
    np.savetxt(filename_I_opt_Sid,I_opt_Sid,fmt='%6d %.10f')
    
    
# optimise for inclination at Quito
filename_I_opt_Quito = ('dIdt_Quito_LSMOD1_timeseries.txt')
try:
    I_opt_Quito = np.loadtxt(filename_I_opt_Quito)
except:
    print('no incliation file found')
    
if I_opt_Quito == []:
    for it in times_it:
        MODEL = '../models/LSMOD1/'+str(int(it))+'.dat'
        FLAG=5
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_Quito,lon_Quito,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        opt_file.readline() # header
        line=opt_file.readline()
        x = [list(map(float, line.split() ))]
        opt_file.close()
        I_opt_Quito=np.append(I_opt_Quito,[it, x[0][0]])
        os.system('rm *DAT')
        
    I_opt_Quito = np.reshape(I_opt_Quito,(len(I_opt_Quito)/2, 2))     
    np.savetxt(filename_I_opt_Quito,I_opt_Quito,fmt='%6d %.10f')
    
    
# optimise for dipole tilt
filename_tilt_opt = ('dtilt_dt_LSMOD1_timeseries.txt')
try:
    tilt_opt = np.loadtxt(filename_tilt_opt)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for it in times_it:
        MODEL = '../models/LSMOD1/'+str(int(it))+'.dat'
        FLAG=0
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        opt_file.readline() # header
        line=opt_file.readline()
        x = [list(map(float, line.split() ))]
        opt_file.close()
        tilt_opt=np.append(tilt_opt,[it, x[0][0]])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/2, 2))     
    np.savetxt(filename_tilt_opt,tilt_opt,fmt='%6d %.10f')
    
# optimise for g10
# IMMAB4 is in milliTesla, not nanoTesla!!!!!!!
filename_g10_opt = ('dg10dt_LSMOD1_timeseries.txt')
try:
    g10_opt = np.loadtxt(filename_g10_opt)
except:
    print('no g10 file found')
    
if g10_opt == []:
    for it in times_it:
        MODEL = '../models/LSMOD1/'+str(int(it))+'.dat'
        FLAG=4
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        opt_file.readline() # header
        line=opt_file.readline()
        x = [list(map(float, line.split() ))]
        opt_file.close()
        g10_opt=np.append(g10_opt,[it, x[0][0]])
        os.system('rm *DAT')
        
    g10_opt = np.reshape(g10_opt,(len(g10_opt)/2, 2))     
    np.savetxt(filename_g10_opt,g10_opt,fmt='%6d %.10f')
    
    
    
    
######################
# plot optimal time series
#####################

# Dipole tilt

fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the excursion times
ax_twin.plot([age[tML], age[tML]],[0.05, 0.2],'--',color='gray')
ax_twin.plot([age[tL], age[tL]],[0.05, 0.2],'--',color='gray')
ax_twin.set_ylim(0.05, 0.2)

ax.set_xlabel('Time / kyr')
ax.set_ylabel(r'Max $|d\theta_d/dt|$ ($deg/yr$)',color='b')

ax.plot(age,tilt_opt[:,1]*180/np.pi,color='b')
ax.set_xlim(age[-1], age[0])
ax_twin.plot(age,g10_opt[:,1]/1000.,color='r')
ax_twin.set_ylabel(r'Max $|dg_1^0/dt|$ ($\mu T/yr$)',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')

plt.title(r'Max dipole tilt and $g_1^0$ rate of change')
plt.savefig(folder+'figures/LSMOD1_opt_dipole.pdf',bbox_inches='tight',pad_inches=0.0)




# Inclination time-series
inclination_locations = np.arctan(np.divide(-Br_locations,np.sqrt(Bt_locations**2+Bp_locations**2)))*180/np.pi

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([age[tML], age[tML]],[0, 12],'--',color='gray')
ax_i.plot([age[tL], age[tL]],[0, 12],'--',color='gray')
ax_i.set_ylim(0, 12)

ax_i.plot(age,I_opt_SUL[:,1]*180/np.pi,color='r',
          label='Sulmona'  )
ax_i.plot(age,I_opt_Quito[:,1]*180/np.pi,color='b',
          label='Quito'  )
ax_i.plot(age,I_opt_Sid[:,1]*180/np.pi,color='g',
          label='Sidney'  )
plt.xlabel('Time / kyr')
plt.ylabel(r'Max $|dI/dt|$ ($deg/yr$)')
ax_i.set_xlim(age[-1], age[0])
ax_i.legend(fontsize=10,loc='upper left')
plt.title('Max inclination rate of change at locations')
plt.show
plt.savefig(folder+'figures/LSMOD1_opt_inclinations.pdf',bbox_inches='tight',pad_inches=0.0)







    
    
############################################
# SENSITIVITY STUDY
# Optimal inclination from 42.85 LSMOD1 model
###########################################


# sensitivity calculations
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 10
MODEL       = '../models/LSMOD1/258.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0
FLAG         = 5
 

# fix LMAX_B_OBS and vary LMAX_U for unrestricted flow
filename = 'dIdt_LSMOD1_42850_LMAX_B_OBS_'+str(LMAX_B_OBS)+'.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')

Iopt_fixed_LB = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for unrestricted flow
LMAX_U      = 25
filename = 'dIdt_LSMOD1_42850_LMAX_U_'+str(LMAX_U)+'.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')
    
Iopt_fixed_LU = np.loadtxt(filename)




# fix LMAX_B_OBS and vary LMAX_U for columnar flow
RESTRICTION  = 3
LMAX_U      = 25
LMAX_B_OBS  = 10
filename = 'dIdt_LSMOD1_42850_LMAX_B_OBS_'+str(LMAX_B_OBS)+'_columnar.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')

Iopt_fixed_LB_columnar = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for columnar flow
LMAX_U      = 25
filename = 'dIdt_LSMOD1_42850_LMAX_U_'+str(LMAX_U)+'_columnar.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')
    
Iopt_fixed_LU_columnar = np.loadtxt(filename)


# plot inclination sensitivity to flow field truncation
fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle('Max $dI/dt$ (in $deg/yr$) at SUL (LSMOD1, 42.85 ka)',fontsize = 20)

#ax_i.set_ylim(-90, 90)

axs[0].plot(Iopt_fixed_LB[:,0],Iopt_fixed_LB[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[0].plot(Iopt_fixed_LB_columnar[:,0],Iopt_fixed_LB_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[0].set_title('Fixed $L_B=13$',fontsize=15)
axs[0].set_xlabel('$L_U$')
#axs[0].set_ylabel('Max $dI/dt$ / (deg/yr)',fontsize=16)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend(fontsize=10,loc='lower right')
axs[0].grid(which='both',alpha=0.5)

# plot inclination sensitivity
axs[1].plot(Iopt_fixed_LU[:,1],Iopt_fixed_LU[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[1].plot(Iopt_fixed_LU_columnar[:,1],Iopt_fixed_LU_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[1].set_title('Fixed $L_U=25$',fontsize=15)
axs[1].set_xlabel('$L_B$')
#axs[1].set_ylabel('Max dI/Dt / (deg/yr)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend(fontsize=10,loc='upper left')
axs[1].grid(which='both',alpha=0.5)

fig.tight_layout(rect=[0,0, 1, 0.94])

plt.show
plt.savefig('dIdt_sensitivity_LSMOD1_42850_SUL.pdf',bbox_inches='tight',pad_inches=0.0)



############################################
# SENSITIVITY STUDY
# Optimal dipole tilt from 42850 LSMOD1 model
###########################################

# sensitivity calculations
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 10
MODEL       = '../models/LSMOD1/258.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0
FLAG         = 0
 

# fix LMAX_B_OBS and vary LMAX_U for unrestricted flow
filename = 'dtiltdt_LSMDO1_42850_LMAX_B_OBS_'+str(LMAX_B_OBS)+'.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')

tilt_opt_fixed_LB = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for unrestricted flow
LMAX_U      = 25
filename = 'dtiltdt_LSMOD1_42850_LMAX_U_'+str(LMAX_U)+'.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')
    
tilt_opt_fixed_LU = np.loadtxt(filename)




# fix LMAX_B_OBS and vary LMAX_U for columnar flow
RESTRICTION  = 3
LMAX_U      = 25
LMAX_B_OBS  = 10
filename = 'dtiltdt_LSMOD1_42850_LMAX_B_OBS_'+str(LMAX_B_OBS)+'_columnar.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')

tilt_opt_fixed_LB_columnar = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for columnar flow
LMAX_U      = 25
filename = 'dtiltdt_LSMOD1_42850_LMAX_U_'+str(LMAX_U)+'_columnar.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if tilt_opt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')
    
tilt_opt_fixed_LU_columnar = np.loadtxt(filename)


# plot tilt sensitivity to flow field truncation
fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle(r'Max $d{\theta}_d/dt$ (in $deg/yr$) (LSMOD1, 42.85 ka)',fontsize = 20)

#ax_i.set_ylim(-90, 90)

axs[0].plot(tilt_opt_fixed_LB[:,0],tilt_opt_fixed_LB[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[0].plot(tilt_opt_fixed_LB_columnar[:,0],tilt_opt_fixed_LB_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[0].set_title('Fixed $L_B=13$',fontsize=15)
axs[0].set_xlabel('$L_U$')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend(fontsize=10,loc='lower right')
axs[0].grid(which='both',alpha=0.5)

# plot tilt sensitivity
axs[1].plot(tilt_opt_fixed_LU[:,1],tilt_opt_fixed_LU[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[1].plot(tilt_opt_fixed_LU_columnar[:,1],tilt_opt_fixed_LU_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[1].set_title('Fixed $L_U=25$',fontsize=15)
axs[1].set_xlabel('$L_B$')
#axs[1].set_ylabel('Max dI/Dt / (deg/yr)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend(fontsize=10,loc='upper left')
axs[1].grid(which='both',alpha=0.5)

fig.tight_layout(rect=[0,0, 1, 0.94])

plt.show
plt.savefig('dtheta_d_dt_sensitivity_LSMOD1_42850.pdf',bbox_inches='tight',pad_inches=0.0)







############################################
# SENSITIVITY STUDY
# Optimal inclination from 49 ka LSMOD1 model
###########################################


# sensitivity calculations
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 10
MODEL       = '../models/LSMOD1/381.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0
FLAG         = 5
 

# fix LMAX_B_OBS and vary LMAX_U for unrestricted flow
filename = 'dIdt_LSMOD1_49000_LMAX_B_OBS_'+str(LMAX_B_OBS)+'.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')

Iopt_fixed_LB = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for unrestricted flow
LMAX_U      = 25
filename = 'dIdt_LSMOD1_49000_LMAX_U_'+str(LMAX_U)+'.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')
    
Iopt_fixed_LU = np.loadtxt(filename)




# fix LMAX_B_OBS and vary LMAX_U for columnar flow
RESTRICTION  = 3
LMAX_U      = 25
LMAX_B_OBS  = 10
filename = 'dIdt_LSMOD1_49000_LMAX_B_OBS_'+str(LMAX_B_OBS)+'_columnar.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')

Iopt_fixed_LB_columnar = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for columnar flow
LMAX_U      = 25
filename = 'dIdt_LSMOD1_49000_LMAX_U_'+str(LMAX_U)+'_columnar.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
Iopt = []
try:
    Iopt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if Iopt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt(filename,Iopt,fmt='%6d %6d %.10f')
    
Iopt_fixed_LU_columnar = np.loadtxt(filename)


# plot inclination sensitivity to flow field truncation
fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle('Max $dI/dt$ (in $deg/yr$) at SUL (LSMOD1, 49 ka)',fontsize = 20)

#ax_i.set_ylim(-90, 90)

axs[0].plot(Iopt_fixed_LB[:,0],Iopt_fixed_LB[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[0].plot(Iopt_fixed_LB_columnar[:,0],Iopt_fixed_LB_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[0].set_title('Fixed $L_B=13$',fontsize=15)
axs[0].set_xlabel('$L_U$')
#axs[0].set_ylabel('Max $dI/dt$ / (deg/yr)',fontsize=16)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend(fontsize=10,loc='lower right')
axs[0].grid(which='both',alpha=0.5)

# plot inclination sensitivity
axs[1].plot(Iopt_fixed_LU[:,1],Iopt_fixed_LU[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[1].plot(Iopt_fixed_LU_columnar[:,1],Iopt_fixed_LU_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[1].set_title('Fixed $L_U=25$',fontsize=15)
axs[1].set_xlabel('$L_B$')
#axs[1].set_ylabel('Max dI/Dt / (deg/yr)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend(fontsize=10,loc='upper left')
axs[1].grid(which='both',alpha=0.5)

fig.tight_layout(rect=[0,0, 1, 0.94])

plt.show
plt.savefig('dIdt_sensitivity_LSMOD1_49000_SUL.pdf',bbox_inches='tight',pad_inches=0.0)



############################################
# SENSITIVITY STUDY
# Optimal dipole tilt from 49 ka LSMOD1 model
###########################################

# sensitivity calculations
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 10
MODEL       = '../models/LSMOD1/381.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0
FLAG         = 0
 

# fix LMAX_B_OBS and vary LMAX_U for unrestricted flow
filename = 'dtiltdt_LSMDO1_49000_LMAX_B_OBS_'+str(LMAX_B_OBS)+'.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')

tilt_opt_fixed_LB = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for unrestricted flow
LMAX_U      = 25
filename = 'dtiltdt_LSMOD1_49000_LMAX_U_'+str(LMAX_U)+'.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')
    
tilt_opt_fixed_LU = np.loadtxt(filename)




# fix LMAX_B_OBS and vary LMAX_U for columnar flow
RESTRICTION  = 3
LMAX_U      = 25
LMAX_B_OBS  = 10
filename = 'dtiltdt_LSMOD1_49000_LMAX_B_OBS_'+str(LMAX_B_OBS)+'_columnar.txt'
LMAX_U_list =np.concatenate((np.linspace(1,25,num=25),np.linspace(30,100,num=8)),axis=0)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no tilt file found')
    
if tilt_opt == []:
    for LMAX_U in LMAX_U_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')

tilt_opt_fixed_LB_columnar = np.loadtxt(filename)

# fix LMAX_U and vary LMAX_B_OBS for columnar flow
LMAX_U      = 25
filename = 'dtiltdt_LSMOD1_49000_LMAX_U_'+str(LMAX_U)+'_columnar.txt'
# define the LMAX_B_OBSs on which to test the code
LMAX_B_OBS_list = np.linspace(1,10,num=10)
tilt_opt = []
try:
    tilt_opt = np.loadtxt(filename)
except:
    print('no incliation file found')
    
if tilt_opt == []:
    for LMAX_B_OBS in LMAX_B_OBS_list:
        subs.write_optimal_flow_input("input_pedagogical_examples",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_pedagogical_examples')
        tilt_opt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        tilt_opt_file.readline() # header
        line=tilt_opt_file.readline()
        x = [list(map(float, line.split() ))]
        tilt_opt_file.close()
        tilt_opt=np.append(tilt_opt,[LMAX_U, LMAX_B_OBS, x[0][0]])
        print(LMAX_U, LMAX_B_OBS, x[0][0])
        os.system('rm *DAT')
        
    tilt_opt = np.reshape(tilt_opt,(len(tilt_opt)/3, 3))     
    np.savetxt(filename,tilt_opt,fmt='%6d %6d %.10f')
    
tilt_opt_fixed_LU_columnar = np.loadtxt(filename)


# plot tilt sensitivity to flow field truncation
fig, axs = plt.subplots(1, 2, figsize=(10,5))
fig.suptitle(r'Max $d{\theta}_d/dt$ (in $deg/yr$) (LSMOD1, 49 ka)',fontsize = 20)

#ax_i.set_ylim(-90, 90)

axs[0].plot(tilt_opt_fixed_LB[:,0],tilt_opt_fixed_LB[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[0].plot(tilt_opt_fixed_LB_columnar[:,0],tilt_opt_fixed_LB_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[0].set_title('Fixed $L_B=13$',fontsize=15)
axs[0].set_xlabel('$L_U$')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend(fontsize=10,loc='lower right')
axs[0].grid(which='both',alpha=0.5)

# plot tilt sensitivity
axs[1].plot(tilt_opt_fixed_LU[:,1],tilt_opt_fixed_LU[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='b',
          label='unrestricted flow')
axs[1].plot(tilt_opt_fixed_LU_columnar[:,1],tilt_opt_fixed_LU_columnar[:,2]* 180.0/np.pi,'--k',marker='s',markerfacecolor='g',
          label='columnar flow')
axs[1].set_title('Fixed $L_U=25$',fontsize=15)
axs[1].set_xlabel('$L_B$')
#axs[1].set_ylabel('Max dI/Dt / (deg/yr)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend(fontsize=10,loc='upper left')
axs[1].grid(which='both',alpha=0.5)

fig.tight_layout(rect=[0,0, 1, 0.94])

plt.show
plt.savefig('dtheta_d_dt_sensitivity_LSMOD1_49000.pdf',bbox_inches='tight',pad_inches=0.0)


    
#######
#######
# plots
#######
#######

######################
# cmb radial field
######################

it = tML

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Br_c[it], cmap='coolwarm')
cs = map1.contour(xs,ys,Br_c[it], 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(age[it], y=1.08)
plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
fname = figfolder +'Br_CMB_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)

# Animation
if animation_flag==1:
    files = []
    
    for it in range(n_instants):
        plt.figure(figsize=(11,6))
        #Mollweide projectionfrom scipy.interpolate import griddata
        map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
        
        # draw coastlines, country boundaries, fill continents.
        map1.drawcoastlines(linewidth=0.50)
        
        # draw lat/lon grid lines every 30 degrees.
        #map1.drawmeridians(np.arange(0,360,30),linewidth=1)
        #map1.drawparallels(np.arange(-90,90,30),linewidth=1)
        bounds = [-800, -600, -400, -200, 0, 200, 400, 600, 800]
    
        xs, ys = map1(lons, lats)
        cf = map1.contourf(xs,ys,Br_c[it], cmap='coolwarm',
                           vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend ='both')
        cs = map1.contour(xs,ys,Br_c[it], 23,
                          colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
        
        plt.colorbar(cf,fraction=0.03, pad=0.04)
        plt.title(age[it], y=1.08)
        plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'Br_CMB_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
#        dt = times[1]-times[0] #assuming a constant timestep
        
    #fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'Br_CMB_%07d.png -vcodec libx264 -y -an ' +folder+ 'Br_CMB.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
    #os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')


######################
# surface radial field
######################
# plot to check things

it = tML

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Br_a[it], cmap='coolwarm')
cs = map1.contour(xs,ys,Br_a[it], 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(age[it], y=1.08)
plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
fname = figfolder +'Br_a_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)

# Animation
if animation_flag==1:
    files = []
    
    for it in range(n_instants):
        plt.figure(figsize=(11,6))
        #Mollweide projectionfrom scipy.interpolate import griddata
        map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
        
        # draw coastlines, country boundaries, fill continents.
        map1.drawcoastlines(linewidth=0.50)
        
        # draw lat/lon grid lines every 30 degrees.
        #map1.drawmeridians(np.arange(0,360,30),linewidth=1)
        #map1.drawparallels(np.arange(-90,90,30),linewidth=1)
        bounds = [-48, -36, -24, -12, 0, 12, 24, 36, 48]
    
        xs, ys = map1(lons, lats)
        cf = map1.contourf(xs,ys,Br_a[it], cmap='coolwarm',
                           vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend ='both')
        cs = map1.contour(xs,ys,Br_a[it], 23,
                          colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
        
        plt.colorbar(cf,fraction=0.03, pad=0.04)
        plt.title(age[it], y=1.08)
        plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'Br_a_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
        #dt = times[1]-times[0] #assuming a constant timestep
        
    #fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'Br_a_%07d.png -vcodec libx264 -y -an ' +folder+ 'Br_a.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
    #os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')




####################
# Surface intensity
###################

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,F_a[it], cmap='jet')
cs = map1.contour(xs,ys,F_a[it], 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(age[it], y=1.08)
plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
fname = figfolder +'F_surface_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)

# animation (can be modified to only plot one, of course)
if animation_flag==1:

    files = []
    
    for it in range(n_instants):
        plt.figure(figsize=(11,6))
        #Mollweide projectionfrom scipy.interpolate import griddata
        map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
        
        # draw coastlines, country boundaries, fill continents.
        map1.drawcoastlines(linewidth=0.50)
        
        # draw lat/lon grid lines every 30 degrees.
        #map1.drawmeridians(np.arange(0,360,30),linewidth=1)
        #map1.drawparallels(np.arange(-90,90,30),linewidth=1)
        bounds = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60]
    
        xs, ys = map1(lons, lats)
        cf = map1.contourf(xs,ys,F_a[it], cmap='jet',
                           vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend = 'both')
        cs = map1.contour(xs,ys,F_a[it], 23,
                          colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
        """
        # find border of SAA
        idxBorder = np.argmin(np.abs(cs.levels -32000))
        cs.collections[idxBorder].set_color('white')
        cs.collections[idxBorder].set_linewidth(3)
        cs.collections[idxBorder].set_alpha(1)
        """
        plt.colorbar(cf,fraction=0.03, pad=0.04)
        plt.title(age[it], y=1.08)
        plt.text(420e5, 9e6, r'$n T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'F_surface_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
    
    
    #dt = times[1]-times[0] #assuming a constant timestep
    #fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'F_surface_%07d.png -vcodec libx264 -y -an ' +folder+ 'F_surface.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
    #os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')


######################
# spectra
#####################

it = tML

L_CHAOS = 10

fig, axs = plt.subplots(1, 2, figsize=(8,5))
fig.suptitle(str(age[it]))

axs[0].set_yscale('log')
axs[0].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_a_spectra[0:L_CHAOS]/1e6,'o',color='grey')
axs[0].plot([1,2,3,4,5,6,7,8,9,10],ME_a_spectra[it,1:]/1e6,'ro')
#axs[0].set_ylim([np.amin(ME_a_spectra[:,1:]), np.amax(ME_a_spectra[:,1:])  ])
#axs[0].set_ylim([1e-2, 1e4 ])
axs[0].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
axs[0].set_title('Surface')
axs[0].set_xlabel('SH degree')
axs[0].set_ylabel('$n T^2$')
axs[0].grid(True)

axs[1].set_yscale('log')
axs[1].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_cmb_spectra[0:L_CHAOS]/1e6,'o',color='grey')
axs[1].plot([1,2,3,4,5,6,7,8,9,10],ME_cmb_spectra[it,1:]/1e6,'ro')
#axs[1].set_ylim([np.amin(ME_cmb_spectra[:,1:]), np.amax(ME_cmb_spectra[:,1:])  ])
#axs[1].set_ylim([1e2, 1e6])
axs[1].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
axs[1].set_title('CMB')
axs[1].set_xlabel('SH degree')
axs[1].set_ylabel('$n T^2$')
axs[1].grid(True)

fig.tight_layout(rect=[0,0, 1, 0.95])
plt.show()


if animation_flag==1:

    files = []
    
    for it in range(155,n_instants):
        fig, axs = plt.subplots(1, 2, figsize=(8,5))
        fig.suptitle(str(age[it]))
        
        axs[0].set_yscale('log')
        axs[0].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_a_spectra[0:L_CHAOS]/1e6,'o',color='grey')
        axs[0].plot([1,2,3,4,5,6,7,8,9,10],ME_a_spectra[it,1:]/1e6,'ro')
        #axs[0].set_ylim([np.amin(ME_a_spectra[:,1:]), np.amax(ME_a_spectra[:,1:])  ])
        #axs[0].set_ylim([1e-2, 1e4 ])
        axs[0].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
        axs[0].set_title('Surface')
        axs[0].set_xlabel('SH degree')
        axs[0].set_ylabel('$n T^2$')
        axs[0].grid(True)
        
        axs[1].set_yscale('log')
        axs[1].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_cmb_spectra[0:L_CHAOS]/1e6,'o',color='grey')
        axs[1].plot([1,2,3,4,5,6,7,8,9,10],ME_cmb_spectra[it,1:]/1e6,'ro')
        #axs[1].set_ylim([np.amin(ME_cmb_spectra[:,1:]), np.amax(ME_cmb_spectra[:,1:])  ])
        #axs[1].set_ylim([1e2, 1e6])
        axs[1].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
        axs[1].set_title('CMB')
        axs[1].set_xlabel('SH degree')
        axs[1].set_ylabel('$n T^2$')
        axs[1].grid(True)
        
        fig.tight_layout(rect=[0,0, 1, 0.95])
    
        # saving the plot
        fname = folder +'figures/spectra_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
    
    
    #dt = times[1]-times[0] #assuming a constant timestep
    #fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'figures/spectra_%07d.png -vcodec libx264 -y -an ' +folder+ 'figures/spectra.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')



# combine multiple plots and animations

if animation_flag==1:

    files = []
    
    for it in range(n_instants):
        fig, axs = plt.subplots(2, 2, figsize=(8,5))
        
        fig.suptitle(str(age[it]))
        
        axs[0,0].set_yscale('log')
        axs[0,0].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_a_spectra[0:L_CHAOS]/1e6,'o',color='grey', label='2019 field')
        axs[0,0].plot([1,2,3,4,5,6,7,8,9,10],ME_a_spectra[it,1:]/1e6,'ro', label='LSMOD1')
        #axs[0].set_ylim([np.amin(ME_a_spectra[:,1:]), np.amax(ME_a_spectra[:,1:])  ])
        #axs[0,0].set_ylim([1e-2, 1e4 ])
        axs[0,0].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
        axs[0,0].set_title('Surface energy spectra')
        axs[0,0].set_xlabel('SH degree')
        axs[0,0].set_ylabel('$n T^2$')
        axs[0,0].grid(True)
        axs[0,0].legend()
        
        axs[1,0].set_yscale('log')
        axs[1,0].plot(np.linspace(1,L_CHAOS,L_CHAOS),CHAOS_cmb_spectra[0:L_CHAOS]/1e6,'o',color='grey', label='2019 field')
        axs[1,0].plot([1,2,3,4,5,6,7,8,9,10],ME_cmb_spectra[it,1:]/1e6,'ro', label='LSMOD1')
        #axs[1].set_ylim([np.amin(ME_cmb_spectra[:,1:]), np.amax(ME_cmb_spectra[:,1:])  ])
        #axs[1,0].set_ylim([1e2, 1e6])
        axs[1,0].set_xticks(np.linspace(1,L_CHAOS,L_CHAOS))
        axs[1,0].set_title('CMB energy spectra')
        axs[1,0].set_xlabel('SH degree')
        axs[1,0].set_ylabel('$n T^2$')
        axs[1,0].grid(True)
        axs[1,0].legend()
        
        plt.subplot(222)
        map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
        
        # draw coastlines, country boundaries, fill continents.
        map1.drawcoastlines(linewidth=0.50)
        bounds = [6, 12, 18, 24, 30, 36, 42, 48, 54, 60]
    
        xs, ys = map1(lons, lats)
        cf = map1.contourf(xs,ys,F_a[it], cmap='jet',
                           vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend = 'both')
        cs = map1.contour(xs,ys,F_a[it], 23,
                          colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
        
        plt.colorbar(cf,fraction=0.03, pad=0.04)
        plt.text(438e5, 10e6, r'$n T$' , fontsize=8)
        plt.title('Surface field intensity')
        
        
        
        plt.subplot(224)
        
        map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

        # draw coastlines, country boundaries, fill continents.
        map1.drawcoastlines(linewidth=0.50)
        
        
        bounds = [-800, -600, -400, -200, 0, 200, 400, 600, 800]
    
        xs, ys = map1(lons, lats)
        cf = map1.contourf(xs,ys,Br_c[it], cmap='coolwarm',
                           vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend ='both')
        cs = map1.contour(xs,ys,Br_c[it], 23,
                          colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
        
        plt.colorbar(cf,fraction=0.03, pad=0.04)
        plt.text(438e5, 9e6, r'$n T$' , fontsize=8)
        plt.title('CMB radial field')
        
        
        
        fig.tight_layout(rect=[0,0, 1, 0.95])

        # saving the plot
        fname = folder +'figures/LSMOD1_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
    
    
    #dt = times[1]-times[0] #assuming a constant timestep
    #fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'figures/LSMOD1_%07d.png -vcodec libx264 -y -an ' +folder+ 'figures/LSMOD1.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -reinit_filter 0')







######################
# time series
#####################

# Dipole tilt

fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the excursions
ax_twin.plot([age[tML], age[tML]],[2, 32],'--',color='gray')
ax_twin.plot([age[tL], age[tL]],[2, 32],'--',color='gray')
ax_twin.set_ylim(2, 32)

ax.set_xlabel('Time / kyr')
ax.set_ylabel('Dipole latitude/$^\circ$',color='b')

ax.plot(age,dipole_colatitude,color='b')
ax.set_xlim(age[-1], age[0])
ax_twin.plot(age,m1/1000.,color='r')
ax_twin.set_ylabel('Dipole intensity/ $\mu T$',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')

plt.title('Dipole latitude and intensity evolution')
plt.savefig(folder+'figures/Dipole_tilt.pdf',bbox_inches='tight',pad_inches=0.0)


# Dipole intensity

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([age[tML], age[tML]],[-35, 0],'--',color='gray')
ax_i.plot([age[tL], age[tL]],[-35, 0],'--',color='gray')
ax_i.set_ylim(-35,0)
plt.xlabel('Time / kyr')
plt.ylabel('$g_1^0$ / $\mu T$')

ax_i.plot(age,g10/1000.,color='b')
ax_i.set_xlim(age[-1], age[0])
ax_i.legend(fontsize=10,loc='lower left')
plt.title('$g_1^0$')
plt.show
plt.savefig(folder+'figures/g10.pdf',bbox_inches='tight',pad_inches=0.0)



# Inclination time-series
inclination_locations = np.arctan(np.divide(-Br_locations,np.sqrt(Bt_locations**2+Bp_locations**2)))*180/np.pi

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([age[tML], age[tML]],[-90, 90],'--',color='gray')
ax_i.plot([age[tL], age[tL]],[-90, 90],'--',color='gray')
ax_i.set_ylim(-90, 90)

ax_i.plot(age,inclination_locations[:,1],color='r',
          label='Sulmona'  )
ax_i.plot(age,inclination_locations[:,0],color='b',
          label='Quito'  )
ax_i.plot(age,inclination_locations[:,2],color='g',
          label='Sidney'  )
plt.xlabel('Time / kyr')
plt.ylabel('Inclination / $^\circ$')
ax_i.set_xlim(age[-1], age[0])
ax_i.legend(fontsize=10,loc='upper left')
plt.title('Inclination at locations')
plt.show
plt.savefig(folder+'figures/inclinations.pdf',bbox_inches='tight',pad_inches=0.0)


