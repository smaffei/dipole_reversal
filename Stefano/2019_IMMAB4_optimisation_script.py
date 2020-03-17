#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:16:01 2020

@author: earsmaf
"""


import numpy as np
from scipy import interpolate
from scipy.interpolate import SmoothSphereBivariateSpline as SSBS
import os
import math
import sys

sys.path.append('./SpecialFunctions/SphericalHarmonics/')
import SH_library

import subs

import matplotlib.tri as tri

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap


######################################################
# define the grid for plotting, interpolating and such
# this is the format accepted by the plotting routines
######################################################

lons_grid = np.linspace(-np.pi,np.pi,num=100)
lats_grid = np.linspace(np.pi/2,-np.pi/2,num=100)

LONS, LATS = np.meshgrid(lons_grid, lats_grid)
THETA = -LATS + np.pi/2

r_a = 6371.0
r_c = 3485.0


##############################
# get initial IMMAB4 model
##############################
MODEL       = '../models/IMMAB4/123.out' # When the reversal starts: 781.8
IMMAB4_init = subs.read_coeffs(MODEL,1,4)
# get age of model
IMMA_file = open(MODEL, 'r')
line=IMMA_file.readline()
t_init = float(line)*1000
IMMA_file.close()


###############################
# OPTIMAL INCLINATION, IMMAB4
# optimise I 
###############################

# input file
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 4
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 
FLAG         = 5


 
folder_test = 'IMMAB4_770.9_max_SUL_inclination'
# prepare the input file
subs.write_optimal_flow_input("input_opt",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_opt')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

# read in the optimal inclination
fincl = open(folder_test+'/OPTIMISED_QUANTITY_DOT.DAT', 'r')
line = fincl.readline()
line = fincl.readline()
fincl.close()

incl_dot = float(line) * 180.0/np.pi

# calculte all SV coefficients
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0

subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
U_file = []
try:
    U_file = np.loadtxt(folder_test+'/OPTIMAL_FLOW.DAT')
except:
    print('no flow file found')
    
if U_file == []:
    os.system('./calc_SV < input_calc_SV')
    os.system('mv *DAT '+folder_test)

REVERSE=1

subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
U_file = []
try:
    U_file = np.loadtxt(folder_test+'/OPTIMAL_FLOW.DAT')
except:
    print('no flow file found')
    
if U_file == []:
    os.system('./calc_SV < input_calc_SV')
    os.system('mv SV.DAT SV_reversed.DAT')
    os.system('mv LM_CMB.DAT LM_CMB_reversed.DAT')
    os.system('mv LM_SURFACE.DAT LM_SURFACE_reversed.DAT')
    os.system('mv *DAT '+folder_test)

################################################################
# OPTIMAL REVERSAL TIMES, IMMAB4
# Optimal inclination , timestepping induction eq.
################################################################

ETA     = 0.7
ALPHA   = 0.6
T0      = 0
TF      = 200
DT      = 1
REVERSE=  0

folder_timestep = 'IMMAB4_770.9_max_SUL_inclination_timestep_dt1'
# prepare the input file
subs.write_timesteping_input("input_timestepping_CMB",LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE)
# run the timestepping fortran code
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    os.system('./timestepping_induction < input_timestepping_CMB')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_incl_evolved = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_incl_ev = coeffs_MF_incl_evolved[:,0]

incl_ev = np.zeros(times_incl_ev.shape)
# calculate inclination time-series
for it in range(coeffs_MF_incl_evolved.shape[0]):
    coeffsB = coeffs_MF_incl_evolved[it,1:]
    beta = SH_library.lin2matCoeffs(coeffsB)
    Br_a, Bt_a, Bp_a=SH_library.calcB(beta,np.array([colat_SUL*np.pi/180.0]),np.array([lon_SUL*np.pi/180.0]),r_a,r_a)
    incl_ev[it] = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2 + Bp_a**2)))*180/np.pi
    

################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal inclination, timestepping/opt 
################################################################

folder_timestep = 'IMMAB4_770.9_max_SUL_inclination_timestep_opt_dt1'

FLAG_U_INIT     = 1
OUT_RATE        = 1
#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_incl_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
coeffs_SV_incl_evolved_opt = subs.read_MF_file(folder_timestep+'/SV_COEFFS.DAT')
times_incl_ev_opt = coeffs_MF_incl_evolved_opt[:,0]

incl_ev_opt = np.zeros(times_incl_ev_opt.shape)
g10_incl_ev_opt = np.zeros(times_incl_ev_opt.shape)
# calculate inclination time-series
for it in range(coeffs_MF_incl_evolved_opt.shape[0]):
    coeffsB = coeffs_MF_incl_evolved_opt[it,1:]
    beta = SH_library.lin2matCoeffs(coeffsB)
    Br_a, Bt_a, Bp_a=SH_library.calcB(beta,np.array([colat_SUL*np.pi/180.0]),np.array([lon_SUL*np.pi/180.0]),r_a,r_a)
    g10_incl_ev_opt[it] = beta[0,2]
    incl_ev_opt[it] = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2 + Bp_a**2)))*180/np.pi

subs.show_flow_streamlines(folder_timestep+'/',0,r_c,r_a)


################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal inclination, timestepping/opt, REVERSED FLOW
################################################################

folder_timestep = 'IMMAB4_770.9_max_SUL_inclination_timestep_opt_dt1_reversed'

FLAG_U_INIT     = 1
OUT_RATE        = 1
REVERSE         = 1
TF              = 5

#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
    
    
################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal inclination, columnar flow, timestepping/opt 
################################################################

folder_timestep = 'IMMAB4_770.9_max_SUL_inclination_columnar_timestep_opt_dt1'

RESTRICTION  = 3
REVERSE      = 0
TF           = 200
#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_incl_columnar_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
coeffs_SV_incl_columnar_evolved_opt = subs.read_MF_file(folder_timestep+'/SV_COEFFS.DAT')
times_incl_columnar_ev_opt = coeffs_MF_incl_columnar_evolved_opt[:,0]

incl_columnar_ev_opt = np.zeros(times_incl_columnar_ev_opt.shape)
g10_incl_columnar_ev_opt = np.zeros(times_incl_columnar_ev_opt.shape)
# calculate inclination time-series
for it in range(coeffs_MF_incl_columnar_evolved_opt.shape[0]):
    coeffsB = coeffs_MF_incl_columnar_evolved_opt[it,1:]
    beta = SH_library.lin2matCoeffs(coeffsB)
    Br_a, Bt_a, Bp_a=SH_library.calcB(beta,np.array([colat_SUL*np.pi/180.0]),np.array([lon_SUL*np.pi/180.0]),r_a,r_a)
    g10_incl_columnar_ev_opt[it] = beta[0,2]
    incl_columnar_ev_opt[it] = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2 + Bp_a**2)))*180/np.pi


################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal inclination, toroidal flow, timestepping/opt 
################################################################

folder_timestep = 'IMMAB4_770.9_max_SUL_inclination_toroidal_timestep_opt_dt1'

RESTRICTION  = 2
#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_incl_toroidal_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
coeffs_SV_incl_toroidal_evolved_opt = subs.read_MF_file(folder_timestep+'/SV_COEFFS.DAT')
times_incl_toroidal_ev_opt = coeffs_MF_incl_toroidal_evolved_opt[:,0]

incl_toroidal_ev_opt = np.zeros(times_incl_toroidal_ev_opt.shape)
g10_incl_toroidal_ev_opt = np.zeros(times_incl_toroidal_ev_opt.shape)
# calculate inclination time-series
for it in range(coeffs_MF_incl_toroidal_evolved_opt.shape[0]):
    coeffsB = coeffs_MF_incl_toroidal_evolved_opt[it,1:]
    beta = SH_library.lin2matCoeffs(coeffsB)
    Br_a, Bt_a, Bp_a=SH_library.calcB(beta,np.array([colat_SUL*np.pi/180.0]),np.array([lon_SUL*np.pi/180.0]),r_a,r_a)
    g10_incl_toroidal_ev_opt[it] = beta[0,2]
    incl_toroidal_ev_opt[it] = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2 + Bp_a**2)))*180/np.pi

    
    
###############################
# OPTIMAL REVERSAL TIMES, IMMAB4
# optimise g10 
###############################

# input file
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 4
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 
FLAG         = 4


 
folder_test = 'IMMAB4_770.9_max_g10'
# prepare the input file
subs.write_optimal_flow_input("input_opt",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
U_file = []
try:
    U_file = np.loadtxt(folder_test+'/OPTIMAL_FLOW.DAT')
except:
    print('no flow file found')
    
if U_file == []:
    os.system('./dipole_tilt_bound < input_opt')
    os.system('mkdir '+folder_test)
    os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

# read in the optimal inclination
fopt = open(folder_test+'/OPTIMISED_QUANTITY_DOT.DAT', 'r')
line = fopt.readline()
line = fopt.readline()
fopt.close()

g10_dot = float(line)

FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=1


################################################################
# OPTIMAL REVERSAL TIMES, IMMAB4
# Optimal g10 , timestepping induction eq.
################################################################

ETA     = 0.7
ALPHA   = 0.6
T0      = 0
TF      = 400
DT      = 1

folder_timestep = 'IMMAB4_770.9_max_g10_timestep_dt1'
# prepare the input file
subs.write_timesteping_input("input_timestepping_CMB",LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE)
# run the timestepping fortran code
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    os.system('./timestepping_induction < input_timestepping_CMB')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
# load MF coefficients
coeffs_MF_g10_evolved = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_g10_ev = coeffs_MF_g10_evolved[:,0]
# load SV coefficients 
coeffs_SV_g10_evolved = subs.read_MF_file(folder_timestep+'/SV_COEFFS.DAT')


################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal g10, timestepping/opt 
################################################################

folder_timestep = 'IMMAB4_770.9_max_g10_timestep_opt_dt1'

FLAG_U_INIT     = 1
OUT_RATE        = 1
#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_g10_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_g10_ev_opt = coeffs_MF_g10_evolved_opt[:,0]




###############################
# OPTIMAL REVERSAL TIMES, IMMAB4
# optimise dipole tilt 
###############################

# input file
colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 4
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 
FLAG         = 0


 
folder_test = 'IMMAB4_770.9_max_dipole_tilt'
# prepare the input file
subs.write_optimal_flow_input("input_opt",colat_SUL,lon_SUL,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
U_file = []
try:
    U_file = np.loadtxt(folder_test+'/OPTIMAL_FLOW.DAT')
except:
    print('no flow file found')
    
if U_file == []:
    os.system('./dipole_tilt_bound < input_opt')
    os.system('mkdir '+folder_test)
    os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

# read in the optimal inclination
fopt = open(folder_test+'/OPTIMISED_QUANTITY_DOT.DAT', 'r')
line = fopt.readline()
line = fopt.readline()
fopt.close()

theta_d_dot = float(line)


FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0


################################################################
# OPTIMAL REVERSAL TIMES, IMMAB4
# Optimal dipole tilt , timestepping induction eq.
################################################################

ETA     = 0.7
ALPHA   = 0.6
T0      = 0
TF      = 400
DT      = 1

folder_timestep = 'IMMAB4_770.9_max_dipole_tilt_timestep_dt1'
# prepare the input file
subs.write_timesteping_input("input_timestepping_CMB",LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE)
# run the timestepping fortran code
os.system('./timestepping_induction < input_timestepping_CMB')
os.system('mkdir '+folder_timestep)
os.system('mv *DAT '+folder_timestep)
# load MF coefficients
coeffs_MF_tilt_evolved = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_tilt_ev = coeffs_MF_tilt_evolved[:,0]

g10   = coeffs_MF_tilt_evolved[:,1]
g11   = coeffs_MF_tilt_evolved[:,2]
h11   = coeffs_MF_tilt_evolved[:,3]
md    = np.sqrt(np.square(g10) + np.square(g11) + np.square(h11))
colatitude_tilt_ev = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,md)))
# load SV coefficients 
coeffs_SV_tilt_evolved = subs.read_MF_file(folder_timestep+'/SV_COEFFS.DAT')

################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal dipole tilt , timestepping/opt 
################################################################

folder_timestep = 'IMMAB4_770.9_max_dipole_tilt_timestep_opt_dt1'

FLAG_U_INIT     = 1
OUT_RATE        = 1
#prepare the input file
subs.write_timesteping_opt_input("input_timestepping_CMB_opt",LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat_SUL,lon_SUL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE)
# this could be a long run, let's skip it if we find the MF file
MF_file = []
try:
    MF_file = np.loadtxt(folder_timestep+'/MF_COEFFS.DAT')
except:
    print('no MF file found')
    
if MF_file == []:
    # run the timestepping fortran code
    os.system('./timestepping_induction_opt < input_timestepping_CMB_opt')
    os.system('mkdir '+folder_timestep)
    os.system('mv *DAT '+folder_timestep)
    
# load MF coefficients
coeffs_MF_tilt_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_tilt_ev_opt = coeffs_MF_tilt_evolved_opt[:,0]

g10   = coeffs_MF_tilt_evolved_opt[:,1]
g11   = coeffs_MF_tilt_evolved_opt[:,2]
h11   = coeffs_MF_tilt_evolved_opt[:,3]
md    = np.sqrt(np.square(g10) + np.square(g11) + np.square(h11))
colatitude_tilt_ev_opt = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,md)))


################
# Plots
################

# define initial coefficients for IMMAB4 model
g10 = IMMAB4_init[0,2]
g11 = IMMAB4_init[1,2]
h11 = IMMAB4_init[1,3]
H1  = np.sqrt(g11**2+h11**2)
md  = np.sqrt(g10**2+g11**2+h11**2)
theta_d = np.arccos(g10/np.sqrt(g10**2 + g11**2 + h11**2))

#############
# Dipole tilt
#############
fig_t,ax_t =  plt.subplots(figsize=(8,5))

colatitude_2 = -(90.0 - 180/np.pi * ( theta_d + (times_tilt_ev-times_tilt_ev[0])*theta_d_dot) )
# from optimal theta_d solutions
ax_t.plot(t_init-times_tilt_ev,colatitude_2,color='b',label='linear extrapolation')
ax_t.plot(t_init-times_tilt_ev,colatitude_tilt_ev,color='m',label='timestepping')
ax_t.plot(t_init-times_tilt_ev_opt,colatitude_tilt_ev_opt,color='r',label='timestepping/optimisation')

# highlight the reversal instant
ax_t.plot(t_init-times_tilt_ev_opt,np.zeros(times_tilt_ev_opt.shape),'--',color='k')
plt.xlabel('Time / yr')
plt.ylabel('Dipole latitude / $^\circ$')
#ax_t.set_xlim([t_init, t_init-TF/1000.])
ax_t.set_ylim(-90, 90)
ax_t.invert_xaxis()

ax_t.legend(fontsize=10)
plt.title('Dipole latitude')
plt.show
plt.savefig('dipole_tilt_decay_IMMAB4_770.pdf',bbox_inches='tight',pad_inches=0.0)

#####
# g10
#####
fig_t,ax_t =  plt.subplots(figsize=(8,5))
# from optimal g10 solutions
ax_t.plot(t_init-times_g10_ev,g10 - ( times_g10_ev - times_g10_ev[0] )*g10_dot,color='b',label='linear extrapolation')
ax_t.plot(t_init-times_g10_ev,coeffs_MF_g10_evolved[:,1],color='m',label='timestepping')
ax_t.plot(t_init-times_g10_ev_opt,coeffs_MF_g10_evolved_opt[:,1],color='r',label='timestepping/optimisation')


# highlight the reversal instant
ax_t.plot(t_init-times_g10_ev_opt,np.zeros(times_g10_ev_opt.shape),'--',color='k')
plt.xlabel('Time / yr')
plt.ylabel('g10 / mT')
#ax_t.set_xlim([t_init, t_init-TF/1000.])
ax_t.invert_xaxis()

ax_t.legend(fontsize=10)
plt.title('$g_1^0$ [mT]')
plt.show
plt.savefig('g10_decay_IMMAB4_770.pdf',bbox_inches='tight',pad_inches=0.0)



#########################
# Inclination time-series
#########################
#inclination_locations = np.arctan(np.divide(-Br_locations,np.sqrt(Bt_locations**2+Bp_locations**2)))*180/np.pi

fig_i,ax_i =  plt.subplots(figsize=(8,5))

Br_a, Bt_a, Bp_a=SH_library.calcB(IMMAB4_init,np.array([colat_SUL*np.pi/180.0]),np.array([lon_SUL*np.pi/180.0]),r_a,r_a)
incl_init = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2 + Bp_a**2)))*180/np.pi


ax_i.set_ylim(-90, 90)

ax_i.plot(t_init-times_incl_ev,incl_init + ( times_incl_ev - times_incl_ev[0] )*incl_dot,color='b',
          label='linear extrapolation')
ax_i.plot(t_init-times_incl_ev,incl_ev,color='m',
          label='timestepping'  )
ax_i.plot(t_init-times_incl_ev_opt,incl_ev_opt,color='r',
          label='timestepping/optimisation'  )
ax_i.plot(t_init-times_incl_columnar_ev_opt,incl_columnar_ev_opt,'--',color='r',
          label='timestepping/optimisation, columnar'  )
ax_i.plot(t_init-times_incl_toroidal_ev_opt,incl_toroidal_ev_opt,'-.',color='r',
          label='timestepping/optimisation, toroidal'  )
plt.xlabel('Time / yr')
plt.ylabel('Inclination / $^\circ$')
# highlight the reversal instant
ax_i.plot(t_init-times_incl_ev_opt,np.zeros(times_incl_ev_opt.shape),'--',color='k')
#ax_i.set_xlim([t_init, t_init-TF/1000.])
ax_i.invert_xaxis()

#ax_i.set_xlim(times[0], times[-1])
ax_i.legend(fontsize=10,loc='bottom right')
plt.title('Inclination at SUL')
plt.show
plt.savefig('inclinations_IMMAB4_770.pdf',bbox_inches='tight',pad_inches=0.0)


# equivalent g10
fig_i,ax_i =  plt.subplots(figsize=(8,5))


ax_i.plot(t_init-times_incl_ev_opt,g10_incl_ev_opt,color='r',
          label='timestepping/optimisation'  )
ax_i.plot(t_init-times_incl_columnar_ev_opt,g10_incl_columnar_ev_opt,'--',color='r',
          label='timestepping/optimisation, columnar'  )
ax_i.plot(t_init-times_incl_toroidal_ev_opt,g10_incl_toroidal_ev_opt,'-.',color='r',
          label='timestepping/optimisation, toroidal'  )
plt.xlabel('Time / yr')
plt.ylabel('g10 / mT')
# highlight the reversal instant
ax_i.plot(t_init-times_incl_ev_opt,np.zeros(times_incl_ev_opt.shape),'--',color='k')
#ax_i.set_xlim([t_init, t_init-TF/1000.])
ax_i.invert_xaxis()

#ax_i.set_xlim(times[0], times[-1])
ax_i.legend(fontsize=10,loc='bottom right')
plt.title('$g_1^0$ [mT]')
plt.show
plt.savefig('g10_inclinations_IMMAB4_770.pdf',bbox_inches='tight',pad_inches=0.0)