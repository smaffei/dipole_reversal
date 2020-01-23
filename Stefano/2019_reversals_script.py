#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:02:12 2019

@author: earsmaf

Unique script to run the calculations of the dipole reversal paper and plot the figures

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
###############################
# PEDAGOGICAL EXAMPLES: 
# optimise g10 from g10 field
###############################

# prepare input file:

colat_Leeds = 36.1992
lon_Leeds   = 358.451    
LMAX_U      = 25
LMAX_B_OBS  = 1
MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_axial_dipole.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 
FLAG         = 4
 
folder_test = 'g10_CHAOS6_2019.0007_max_g10'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")


###############################
# PEDAGOGICAL EXAMPLES: 
# optimise g10 from g11 field
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_g11_dipole.dat' 

folder_test = 'g11_CHAOS6_2019.0007_max_g10'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

###############################
# PEDAGOGICAL EXAMPLES: 
# optimise g11 from g10 field
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_axial_dipole.dat' 
FLAG         = 6

folder_test = 'g10_CHAOS6_2019.0007_max_g11'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

###############################
# PEDAGOGICAL EXAMPLES: 
# optimise g11 from g10 field
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_g11_dipole.dat' 

folder_test = 'g11_CHAOS6_2019.0007_max_g11'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

###############################
# PEDAGOGICAL EXAMPLES: 
# optimise dipole tilt from g10 field
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_g10_small_H.dat' 
FLAG         = 0

folder_test = 'g10_small_H_CHAOS6_2019.0007_max_dipole_tilt'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")


###############################
# PEDAGOGICAL EXAMPLES: 
# optimise inclination in Leeds from g10 field
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_axial_dipole.dat' 
FLAG         = 5

folder_test = 'g10_CHAOS6_2019.0007_max_Leeds_inclination'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")

# plot inclination everywhere else
# calculate SV coefficients first
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)
# calculate B and B_dot at Earth surface
Br_a_g10, Bt_a_g10, Bp_a_g10 = SH_library.calcB(coeffsB,THETA,LONS,r_a,r_a)
g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]

dtBr_a_Leeds, dtBt_a_Leeds, dtBp_a_Leeds = SH_library.calcB(coeffsBdot,THETA,LONS,r_a,r_a)
g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]

F_a_g10 = np.sqrt(Br_a_g10**2 + Bt_a_g10**2 + Bp_a_g10**2)
H_a_g10 = np.sqrt(Bt_a_g10**2 + Bp_a_g10**2) 
# calculate dIdt in rad/yr (supposedly)
dIdt_Leeds = (  Bt_a_g10 * Br_a_g10 * dtBt_a_Leeds  \
              + Bp_a_g10 * Br_a_g10 * dtBp_a_Leeds  \
              - H_a_g10**2 * dtBr_a_Leeds) / ( H_a_g10 * F_a_g10**2 )

subs.show_global_contourf(LONS,LATS,dIdt_Leeds,"viridis",folder_test+"/dIdt_map.pdf")



###############################
# now with purely polidal flow
###############################
RESTRICTION  = 1

folder_test = 'g10_CHAOS6_2019.0007_max_Leeds_inclination_poloidal'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")

# plot inclination everywhere else
# calculate SV coefficients first
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)
# calculate B and B_dot at Earth surface
Br_a_g10, Bt_a_g10, Bp_a_g10 = SH_library.calcB(coeffsB,THETA,LONS,r_a,r_a)
g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]

dtBr_a_Leeds, dtBt_a_Leeds, dtBp_a_Leeds = SH_library.calcB(coeffsBdot,THETA,LONS,r_a,r_a)
g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]

F_a_g10 = np.sqrt(Br_a_g10**2 + Bt_a_g10**2 + Bp_a_g10**2)
H_a_g10 = np.sqrt(Bt_a_g10**2 + Bp_a_g10**2) 
# calculate dIdt in rad/yr (supposedly)
dIdt_Leeds = (  Bt_a_g10 * Br_a_g10 * dtBt_a_Leeds  \
              + Bp_a_g10 * Br_a_g10 * dtBp_a_Leeds  \
              - H_a_g10**2 * dtBr_a_Leeds) / ( H_a_g10 * F_a_g10**2 )

subs.show_global_contourf(LONS,LATS,dIdt_Leeds,"viridis",folder_test+"/dIdt_map.pdf")


###############################
# now with purely toroidal flow
###############################
RESTRICTION  = 2

folder_test = 'g10_CHAOS6_2019.0007_max_Leeds_inclination_toroidal'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")

# plot inclination everywhere else
# calculate SV coefficients first
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)
# calculate B and B_dot at Earth surface
Br_a_g10, Bt_a_g10, Bp_a_g10 = SH_library.calcB(coeffsB,THETA,LONS,r_a,r_a)
g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]

dtBr_a_Leeds, dtBt_a_Leeds, dtBp_a_Leeds = SH_library.calcB(coeffsBdot,THETA,LONS,r_a,r_a)
g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]

F_a_g10 = np.sqrt(Br_a_g10**2 + Bt_a_g10**2 + Bp_a_g10**2)
H_a_g10 = np.sqrt(Bt_a_g10**2 + Bp_a_g10**2) 
# calculate dIdt in rad/yr (supposedly)
dIdt_Leeds = (  Bt_a_g10 * Br_a_g10 * dtBt_a_Leeds  \
              + Bp_a_g10 * Br_a_g10 * dtBp_a_Leeds  \
              - H_a_g10**2 * dtBr_a_Leeds) / ( H_a_g10 * F_a_g10**2 )

subs.show_global_contourf(LONS,LATS,dIdt_Leeds,"viridis",folder_test+"/dIdt_map.pdf")

###############################
# now with purely columnar flow
###############################
RESTRICTION  = 3

folder_test = 'g10_CHAOS6_2019.0007_max_Leeds_inclination_columnar'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")

# plot inclination everywhere else
# calculate SV coefficients first
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)
# calculate B and B_dot at Earth surface
Br_a_g10, Bt_a_g10, Bp_a_g10 = SH_library.calcB(coeffsB,THETA,LONS,r_a,r_a)
g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]

dtBr_a_Leeds, dtBt_a_Leeds, dtBp_a_Leeds = SH_library.calcB(coeffsBdot,THETA,LONS,r_a,r_a)
g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]

F_a_g10 = np.sqrt(Br_a_g10**2 + Bt_a_g10**2 + Bp_a_g10**2)
H_a_g10 = np.sqrt(Bt_a_g10**2 + Bp_a_g10**2) 
# calculate dIdt in rad/yr (supposedly)
dIdt_Leeds = (  Bt_a_g10 * Br_a_g10 * dtBt_a_Leeds  \
              + Bp_a_g10 * Br_a_g10 * dtBp_a_Leeds  \
              - H_a_g10**2 * dtBr_a_Leeds) / ( H_a_g10 * F_a_g10**2 )

subs.show_global_contourf(LONS,LATS,dIdt_Leeds,"viridis",folder_test+"/dIdt_map.pdf")

###############################
# PEDAGOGICAL EXAMPLES: 
# optimise inclination in Leeds from g10 field with a small horizontal components
###############################
RESTRICTION  = 0
MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007_g10_small_H.dat' 

folder_test = 'g10_small_H_CHAOS6_2019.0007_max_Leeds_inclination'
# prepare the input file
subs.write_optimal_flow_input("input_pedagogical_examples",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_pedagogical_examples')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)
# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")
subs.show_flow_spectra(folder_test+"/")

# plot inclination everywhere else
# calculate SV coefficients first
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)
# calculate B and B_dot at Earth surface
Br_a_g10, Bt_a_g10, Bp_a_g10 = SH_library.calcB(coeffsB,THETA,LONS,r_a,r_a)
g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]

dtBr_a_Leeds, dtBt_a_Leeds, dtBp_a_Leeds = SH_library.calcB(coeffsBdot,THETA,LONS,r_a,r_a)
g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]

F_a_g10 = np.sqrt(Br_a_g10**2 + Bt_a_g10**2 + Bp_a_g10**2)
H_a_g10 = np.sqrt(Bt_a_g10**2 + Bp_a_g10**2) 
# calculate dIdt in rad/yr (supposedly)
dIdt_Leeds1 = (  Bt_a_g10 * Br_a_g10 * dtBt_a_Leeds  \
              + Bp_a_g10 * Br_a_g10 * dtBp_a_Leeds  \
              - H_a_g10**2 * dtBr_a_Leeds) / ( H_a_g10 * F_a_g10**2 )

subs.show_global_contourf(LONS,LATS,dIdt_Leeds1,"viridis",folder_test+"/dIdt_map.pdf")

# expected dipole tilt
m1Leeds = np.sqrt(g10**2 + g11**2 + h11**2)
h1Leeds = np.sqrt(g11**2 + h11**2)
exp_tilt_Leeds = -(m1Leeds/h1Leeds) * ( g10dot*(m1Leeds**2) - g10*(g10*g10dot + g11*g11dot + h11*h11dot) ) / m1Leeds**3

###############################
# MAP OF OPTIMISED dI/dT: 
###############################

MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007.dat' 
LMAX_B_OBS  = 13 

Nth = 15

(theta_virt,lons_virt, lats_virt)=subs.equispaced_grid(Nth)
Iopt =[]
x=[]

try:
    Iopt = np.loadtxt('dIdt_Nth'+str(Nth)+'.txt')
except:
    print('no incliation file found')
    
if Iopt == []:
    print('no inclination file found, calculating')
    for it in range(lats_virt.size):
        if lats_virt[it]==0 or lats_virt[it]==180:
            continue
        subs.write_optimal_flow_input("input_inclination_map",lats_virt[it],lons_virt[it],LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
        os.system('./dipole_tilt_bound < input_inclination_map')
        Iopt_file = open('./OPTIMISED_QUANTITY_DOT.DAT', 'r')
        Iopt_file.readline() # header
        line=Iopt_file.readline()
        x = [list(map(float, line.split() ))]
        Iopt_file.close()
        Iopt=np.append(Iopt,[lons_virt[it], lats_virt[it], x[0][0]])
        print(lons_virt[it],lats_virt[it],x[0][0])
        os.system('rm *DAT')
    
    Iopt = np.reshape(Iopt,(len(Iopt)/3, 3))     
    np.savetxt('dIdt_Nth'+str(Nth)+'.txt',Iopt,fmt='%.6f')


# plot the result on a scatter plot:
subs.show_global_scatterplot(Iopt[:, 0],-Iopt[:, 1]+90,Iopt[:,2],"viridis",'dIdt_opt_global_scatter.pdf')

# SH least squares interpolation
# max SH degree is set to Nth*2/3 to avoid overinterpolation

beta, chi2 = SH_library.SHLS(Iopt[:,2],Iopt[:,1]*np.pi/180,Iopt[:,0]*np.pi/180,Nth*2/3)

dIdt_global = SH_library.ForwardSH(beta,THETA,LONS)

subs.show_global_contourf(LONS,LATS,dIdt_global,"viridis","dIdt_opt_global_interp.pdf")


# for comparison, let's see the dipolar component of the background field (CHAOS)

coeffsB_CHAOS_2019    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
Br_a_CHAOS_2019, Bt_a_CHAOS_2019, Bp_a_CHAOS_2019 = SH_library.calcB(coeffsB_CHAOS_2019,THETA,LONS,r_a,r_a)
F_a_CHAOS_2019 = np.sqrt(Br_a_CHAOS_2019**2 +  Bt_a_CHAOS_2019**2 + Bp_a_CHAOS_2019**2)

subs.show_global_contourf(LONS,LATS,Br_a_CHAOS_2019,"coolwarm","models/CHAOS-6/Br_a_CHAOS_2019.pdf")
subs.show_global_contourf(LONS,LATS,F_a_CHAOS_2019,"jet","models/CHAOS-6/F_a_CHAOS_2019.pdf")

# just the dipole 
coeffsB_CHAOS_2019_dipole    = subs.read_coeffs(MODEL,1,1)
Br_a_CHAOS_2019_dipole, Bt_a_CHAOS_2019_dipole, Bp_a_CHAOS_2019_dipole = SH_library.calcB(coeffsB_CHAOS_2019_dipole,THETA,LONS,r_a,r_a)
F_a_CHAOS_2019_dipole = np.sqrt(Br_a_CHAOS_2019_dipole**2 +  Bt_a_CHAOS_2019_dipole**2 + Bp_a_CHAOS_2019_dipole**2)

subs.show_global_contourf(LONS,LATS,Br_a_CHAOS_2019_dipole,"coolwarm","models/CHAOS-6/Br_a_CHAOS_2019_dipole.pdf")
subs.show_global_contourf(LONS,LATS,F_a_CHAOS_2019_dipole,"jet","models/CHAOS-6/F_a_CHAOS_2019_dipole.pdf")

"""
fig = plt.figure(2, figsize=(11,6))

#cf = plt.contourf(LONS*180./np.pi, LATS*180./np.pi, dIdt)
cf = plt.contourf(LONS*180./np.pi, LATS*180./np.pi, dIdt)
plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.show()


fig = plt.figure(figsize=(11,6))
    
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

# draw lat/lon grid lines every 30 degrees.
map1.drawmeridians(np.arange(0,360,30),linewidth=1)
map1.drawparallels(np.arange(-90,90,30),linewidth=1)
    
    
xs, ys = map1(LONS*180./np.pi, LATS*180./np.pi)
#cs = map1.contourf(xs,ys,dIdt)
cs = map1.contourf(xs,ys,dIdt,cmap="coolwarm")
plt.colorbar(cs,fraction=0.03, pad=0.04)
"""

##########################################################
# OPTIMAL REVERSAL TIMES:
# Optimal g10 from 2019 CHAOS model, linear extrapolation
##########################################################

colat_SUL = 47.8488
lon_SUL   = 13.8229 
LMAX_U      = 25
LMAX_B_OBS  = 13
MODEL       = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007.dat' 
TARGET_RMS  = 13 
SCALE_FACTOR = 5e-2
RESTRICTION  = 0
ETA          = 0 
FLAG         = 4
 
folder_test = 'CHAOS6_2019.0007_max_g10'
# prepare the input file
subs.write_optimal_flow_input("input_opt",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_opt')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

# calculte all SV coefficients
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=0
subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)

g10 = coeffsB[0,2]
g10_1 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]
md  = np.sqrt(g10**2+g11**2+h11**2)
H1  = np.sqrt(g11**2+h11**2)
theta_d_1 = np.arccos(g10/md)

g10dot = coeffsBdot[0,2]
g10dot_1 = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]
theta_ddot_1 = -( md/H1 )*( g10dot*md**2 - g10*(g10*g10dot+g11*g11dot+h11*h11dot) ) / md**3

################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal g10 from 2019 CHAOS model, timestepping induction eq.
################################################################

ETA     = 0.7
ALPHA   = 0.6
T0      = 2019.007
TF      = 2219.007
DT      = 1

folder_timestep = 'CHAOS6_2019.0007_max_g10_timestep_dt1'
# prepare the input file
subs.write_timesteping_input("input_timestepping_CMB",LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE)
# run the timestepping fortran code
os.system('./timestepping_induction < input_timestepping_CMB')
os.system('mkdir '+folder_timestep)
os.system('mv *DAT '+folder_timestep)
# load MF coefficients
coeffs_MF_1_evolved = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_1_ev = coeffs_MF_1_evolved[:,0]

################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal g10 from 2019 CHAOS model, timestepping/opt 
################################################################

folder_timestep = 'CHAOS6_2019.0007_max_g10_timestep_opt_dt1'

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
coeffs_MF_1_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_1_ev_opt = coeffs_MF_1_evolved_opt[:,0]

############################################
# OPTIMAL REVERSAL TIMES:
# Optimal dipole tilt from 2019 CHAOS model
###########################################

ETA         = 0
FLAG         = 0
 
folder_test = 'CHAOS6_2019.0007_max_dipole_tilt'
# prepare the input file
subs.write_optimal_flow_input("input_opt",colat_Leeds,lon_Leeds,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG)
# run the fortran code
os.system('./dipole_tilt_bound < input_opt')
os.system('mkdir '+folder_test)
os.system('mv *DAT '+folder_test)

# plot the flow
subs.show_flow(folder_test+"/")
subs.show_flow_global(folder_test+"/")

# calculte all SV coefficients
FLOW = folder_test+'/OPTIMAL_FLOW.DAT'
REVERSE=1

subs.write_calc_SV_input('input_calc_SV',LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE)
os.system('./calc_SV < input_calc_SV')
os.system('mv *DAT '+folder_test)
# load SV coefficients
coeffsB    = subs.read_coeffs(MODEL,1,LMAX_B_OBS)
# these coefficients are in nT/yr
coeffsBdot = subs.read_coeffs(folder_test+'/SV.DAT',0,LMAX_B_OBS+LMAX_U)

g10 = coeffsB[0,2]
g11 = coeffsB[1,2]
h11 = coeffsB[1,3]
md  = np.sqrt(g10**2+g11**2+h11**2)
H1  = np.sqrt(g11**2+h11**2)
theta_d_2 = np.arccos(g10/md)

g10dot = coeffsBdot[0,2]
g11dot = coeffsBdot[1,2]
h11dot = coeffsBdot[1,3]
theta_ddot_2 = -( md/H1 )*( g10dot*md**2 - g10*(g10*g10dot+g11*g11dot+h11*h11dot) ) / md**3

################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal diple tilt from 2019 CHAOS model, timestepping induction eq.
################################################################

ETA     = 0.7
ALPHA   = 0.6
T0      = 2019.007
TF      = 2219.007
DT      = 1

folder_timestep = 'CHAOS6_2019.0007_max_dipolte_tilt_timestep_dt1'
# prepare the input file
subs.write_timesteping_input("input_timestepping_CMB",LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE)
# run the timestepping fortran code
os.system('./timestepping_induction < input_timestepping_CMB')
os.system('mkdir '+folder_timestep)
os.system('mv *DAT '+folder_timestep)

coeffs_MF_2_evolved = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_2_ev = coeffs_MF_2_evolved[:,0]
g10   = coeffs_MF_2_evolved[:,1]
g11   = coeffs_MF_2_evolved[:,2]
h11   = coeffs_MF_2_evolved[:,3]
md    = np.sqrt(np.square(g10) + np.square(g11) + np.square(h11))
colatitude_2_ev = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,md)))
dipole_lon_2_ev =  180/np.pi * np.arctan(np.divide(g11,h11)) 


################################################################
# OPTIMAL REVERSAL TIMES:
# Optimal dipole tilt from 2019 CHAOS model, timestepping/opt 
################################################################

folder_timestep = 'CHAOS6_2019.0007_max_dipole_tilt_timestep_opt_dt1'

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
    
coeffs_MF_2_evolved_opt = subs.read_MF_file(folder_timestep+'/MF_COEFFS.DAT')
times_2_ev_opt = coeffs_MF_2_evolved_opt[:,0]
g10   = coeffs_MF_2_evolved_opt[:,1]
g11   = coeffs_MF_2_evolved_opt[:,2]
h11   = coeffs_MF_2_evolved_opt[:,3]
md    = np.sqrt(np.square(g10) + np.square(g11) + np.square(h11))
colatitude_2_ev_opt = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,md)))
dipole_lon_2_ev_opt =  180/np.pi * np.arctan(np.divide(g11,h11)) 



################
# Plots
################

#############
# Dipole tilt
#############
fig_t,ax_t =  plt.subplots(figsize=(8,5))
colatitude_2 = -(90.0 - 180/np.pi * ( theta_d_2 + (times_2_ev-times_2_ev[0])*theta_ddot_2) )
# from optimal theta_d solutions
ax_t.plot(times_2_ev,colatitude_2,color='b',label='linear extrapolation')
ax_t.plot(times_2_ev,colatitude_2_ev,color='m',label='timestepping')
ax_t.plot(times_2_ev_opt,colatitude_2_ev_opt,color='r',label='timestepping/optimisation')

# highlight the reversal instant
ax_t.plot(times_2_ev_opt,np.zeros(times_2_ev_opt.shape),'--',color='k')

ax_t.set_xlim([2019,2220])
#ax_t.set_ylim([60,95])

ax_t.legend(fontsize=10)
plt.title('Dipole latitude')
plt.show
plt.savefig('dipole_tilt_decay_CHAOS_2019.pdf',bbox_inches='tight',pad_inches=0.0)

#####
# g10
#####
fig_t,ax_t =  plt.subplots(figsize=(8,5))
colatitude_2 = -(90.0 - 180/np.pi * ( theta_d_2 + (times_2_ev-times_2_ev[0])*theta_ddot_2) )
# from optimal g10 solutions
ax_t.plot(times_1_ev,g10_1 + ( times_1_ev - times_1_ev[0] )*g10dot_1,color='b',label='linear extrapolation')
ax_t.plot(times_1_ev,coeffs_MF_1_evolved[:,1],color='m',label='timestepping')
ax_t.plot(times_1_ev_opt,coeffs_MF_1_evolved_opt[:,1],color='r',label='timestepping/optimisation')
# from optimal theta_d solutions
ax_t.plot(times_2_ev,colatitude_2,color='b',label='linear extrapolation')
ax_t.plot(times_2_ev,colatitude_2_ev,color='m',label='timestepping')
ax_t.plot(times_2_ev_opt,colatitude_2_ev_opt,color='r',label='timestepping/optimisation')


# highlight the reversal instant
ax_t.plot(times_1_ev_opt,np.zeros(times_1_ev_opt.shape),'--',color='k')

ax_t.set_xlim([2019,2220])
#ax_t.set_ylim([60,95])

ax_t.legend(fontsize=10)
plt.title('$g_1^0$ [nT]')
plt.show
plt.savefig('g10_decay_CHAOS_2019.pdf',bbox_inches='tight',pad_inches=0.0)