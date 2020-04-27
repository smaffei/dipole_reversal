#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:12:41 2020

@author: Stefano

Collect and analyze paleomagnetic dataset for rapid geomagnetic variations in the past
"""


import numpy as np
from scipy import interpolate
from scipy.interpolate import SmoothSphereBivariateSpline as SSBS
import os
import math
import sys
import csv
import glob
import codecs

sys.path.append('./SpecialFunctions/SphericalHarmonics/')
import SH_library

import subs

import matplotlib.tri as tri

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
from matplotlib import gridspec

from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

import time

# get current year
year_str = time.strftime("%Y")
year_now = int(year_str)

r_cmb = 3485.0e3
r_a   = 6371.0e3

# Sulmona colatitute and longitude
colat_SUL = 47.8488
lat_SUL = 90-colat_SUL
lon_SUL   = 13.8229 
# Sanxing Cave coordinates
lat_C2 = 27.366667
lon_C2 = 107.183333


plt.close('all')
#########################################
########################################
# Sagnotti
########################################
########################################


########################################
# Read in Sagnotti's 2014 data
# from the supplementary material
# Depth has been converted to age by fixing 
# the age of the SUl2-16 and SUL2-22 
# layers from figure 4
########################################
SUL_folder = '../datasets/Sagnotti_2014/'
data_S1_file = SUL_folder + 'Sagnotti_et_al_Table_S1_rock_mag_pmag.csv'

data_S1 = []
depth_S1 =[]
Incl_S1 = []
Decl_S1 = []
with open(data_S1_file) as f1:
    readCSV = csv.reader(f1,delimiter=',')
    for row in readCSV:
        data_S1.append(row)
        # get depth data point
        try:
            depth = float(row[0])
        except:
            # not a string: probably still in the header
            continue
        depth_S1 = np.append(depth_S1,depth)
        # get inclination data point
        try:
            incl = float(row[4])
        except:
            incl = np.nan
        # get the declination data point
        try:
            decl = float(row[5])
        except:
            decl = np.nan
            
        Incl_S1 = np.append(Incl_S1,incl)
        Decl_S1 = np.append(Decl_S1,decl)

f1.close()
    
# convert depth to time
time_SUL2_15 = 773.4 #not really used
time_SUL2_16 = 781.3 # depth: 48.5944 / 51.4056
depth_SUL2_16 = (48.5944 + 51.4056)/2
time_SUL2_22 = 791.9 # depth : 275.0964 / 278.3092
depth_SUL2_22 = (275.0964 + 278.3092)/2


# select only data between these intervals
depth_SUL = depth_S1[26:171]
Incl_SUL  = Incl_S1[26:171]
Decl_SUL  = Decl_S1[26:171]

# deposition rate in cm/kyr
dep_rate  = (depth_SUL2_22 - depth_SUL2_16) / (time_SUL2_22 - time_SUL2_16)
t0        = time_SUL2_16 - depth_SUL2_16/dep_rate
time_SUL  = depth_SUL/dep_rate + t0

depth_SUL2_19 = 175.85 # from Sagnotti 2016
time_SUL2_19  = depth_SUL2_19/dep_rate + t0


########################################
# Read in Sagnotti's 2015 data
# taken from reading in figure 13 with WebPlootDigitizer
########################################
SUL_2015_folder = '../datasets/Sagnotti_2015/'
data_2015_file = SUL_2015_folder + 'Inclination.csv'

data_2015 = []
depth_2015 = []
Incl_2015 = []
with open(data_2015_file) as f1:
    readCSV = csv.reader(f1,delimiter=',')
    for row in readCSV:
        data_2015.append(row)
        # get time stamp of data point
        try:
            time1 = float(row[0])
        except:
            # not a string: probably still in the header
            continue
        depth_2015 = np.append(depth_2015,time1)
        # get inclination data point
        try:
            incl = float(row[1])
        except:
            incl = np.nan
            
        Incl_2015 = np.append(Incl_2015,incl)

f1.close()

# sedimentation rate = 0.28 \pm 0.12 mm / yr = cm/ka
sed_rate_2015 = 28.0
time_2015 = - depth_2015 / sed_rate_2015 + time_SUL2_19


# calculate VGP
VGP_SUL_lat, VGP_SUL_lon = subs.VGP_from_DI(Incl_SUL,Decl_SUL,lat_SUL,lon_SUL)

# characteristic rapid variation in deg/yr
dVGPlat_dt_SUL = (77.2906+72.3984)/(786163-786069)
age_fast_SUL = (786163+786069)/2

# SAgnotti 2016 (same variation but in 13 years)
dVGPlat_dt_SUL_2 = (77.2906+72.3984)/(13)

##################################
# plot Inclination and Declination
##################################
fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the transitional period
#ax_twin.plot([times[twrite], times[twrite]],[0, 40],'--',color='lightgray')
#ax_twin.plot([times[tend], times[tend]],[0, 40],'--',color='lightgray')
#ax_twin.set_ylim(0, 40)

ax.set_xlabel('Age / ka')
ax.set_ylabel('Inclination/$^\circ$',color='b')

ax.plot(time_SUL,Incl_SUL,'o-',color='b')
#ax.set_xlim(times[0], times[-1])

ax_twin.plot(time_SUL,Decl_SUL,'o-',color='r')
ax_twin.set_ylabel('Declination/$^\circ$',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')
plt.show
plt.title('Inclination and Declination at Sulmona')


##################################
# plot VGP, lon and lat
##################################
fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the transitional period
#ax_twin.plot([times[twrite], times[twrite]],[0, 40],'--',color='lightgray')
#ax_twin.plot([times[tend], times[tend]],[0, 40],'--',color='lightgray')
#ax_twin.set_ylim(0, 40)

ax.set_xlabel('Age / ka')
ax.set_ylabel('VGP Latitude/$^\circ$',color='b')

ax.plot(time_SUL,VGP_SUL_lat,'o-',color='b')
#ax.set_xlim(times[0], times[-1])

ax_twin.plot(time_SUL,VGP_SUL_lon,'o-',color='r')
ax_twin.set_ylabel('VGP Longitude/$^\circ$',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')
plt.show
plt.title('VGP position from Sulmona data')



####################
# plot VGP position
####################

lats = np.linspace(np.pi/2,-np.pi/2,num=100)
lons = np.linspace(-np.pi,np.pi,num=100)
lons, lats = np.meshgrid(lons, lats)

lats = np.rad2deg(lats)
lons = np.rad2deg(lons)

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(VGP_SUL_lon, VGP_SUL_lat)

sc = map1.scatter(xs,ys, marker='D',c=time_SUL, cmap='rainbow')

plt.colorbar(sc,fraction=0.03, pad=0.04)
plt.title('Sulmona VGPs', y=1.08)



#########################################
########################################
# IMMAB4
########################################
########################################

# characteristic rapid dipole latitude variation in deg/yr
dtheta_dt_IMMAB4 = np.abs( (77.2790- (-70.662)) / (774000 - 776000) )
age_dthetadt_IMMAB4   = (774000 + 776000) / 2.0 # central age for it (in yr)

# VGP lat at SUL in deg / yr
dVGPlat_dt_SUL_IMMAB4 = np.abs( (58.5598- (-41.962)) / (771200 - 776400) )
age_dVGPlat_dt_SUL_IMMAB4   = (771200 + 776400) / 2.0 # central age for it (in yr)

# optimal solution
dVGPlat_dt_SUL_IMMAB4_opt = 26.786474819337922
age_dVGPlat_dt_SUL_IMMAB4_opt =776400

# solution for optimal dIdt
I_opt_SUL_dVGPlat_dt = 4.0365043002999998
age_I_opt_SUL = 777400

# optimal dtheta_dt
IMMAB4_tilt_opt = 1.8747264771166685
age_IMMAB4_tilt_opt = 774100

#########################################
########################################
# Chou et al., 2018 data (speleothems)
########################################
########################################

Chou_2018_folder = '../datasets/Chou_2018/'
data_C2_file = Chou_2018_folder + 'chou_2018_multidecadal_oscillations_dataset02.csv'

data_C2 = []
depth_C2 =[]
Incl_C2 = []
Decl_C2 = []
Age_C2  = []
with open(data_C2_file) as f1:
#    readCSV = csv.reader(f1,delimiter=',')
    readCSV = csv.reader(codecs.EncodedFile(f1, 'utf8', 'utf_8_sig'), delimiter=',')
    for row in readCSV:
        #data_C2.append(map(float,row))
        data_C2.append(row)
        # get depth data point
        try:
            depth = float(row[0])
        except:
            # not a string: probably still in the header
            continue
        depth_C2 = np.append(depth_C2,depth)
        # get inclination data point
        try:
            incl = float(row[5])
        except:
            incl = np.nan
        # get the declination data point
        try:
            decl = float(row[4])
        except:
            decl = np.nan
        # get the age for  data point
        try:
            age = float(row[3])
        except:
            age = np.nan
            
        Incl_C2 = np.append(Incl_C2,incl)
        Decl_C2 = np.append(Decl_C2,decl)
        Age_C2  = np.append(Age_C2,age)
        
f1.close()

# calculate VGP
# declination path were rotated by 150 counterclockwise (Chou et al., 2018)
# it' all relative, as absolute declination was not possible to determine 
# from the available data (not sure why, but that's what the paper says)
VGP_C2_lat, VGP_C2_lon = subs.VGP_from_DI(Incl_C2,Decl_C2-150,lat_C2,lon_C2)

# characteristic rapid variation in deg/yr
#(from figure S10, difference between two points in the opposite hemispheres.
# I could not calculate proper VGPs for this paper....)
dVGPlat_dt_C2 = (39.5- (-44.3)) / 10
age_fast_C2   = (98226 + 98216) / 2.0 # central age for it (in yr)


##################################
# plot Inclination and Declination
##################################
fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the transitional period
#ax_twin.plot([times[twrite], times[twrite]],[0, 40],'--',color='lightgray')
#ax_twin.plot([times[tend], times[tend]],[0, 40],'--',color='lightgray')
#ax_twin.set_ylim(0, 40)

ax.set_xlabel('Age / ka')
ax.set_ylabel('Inclination/$^\circ$',color='b')

ax.plot(Age_C2/1000,Incl_C2,'o-',color='b')
#ax.set_xlim(times[0], times[-1])
ax_twin.plot(Age_C2/1000,Decl_C2-150,'o-',color='r')
ax_twin.set_ylabel('Declination/$^\circ$',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')
plt.show
plt.title('Inclination and Declination at Sanxing')



####################
# plot VGP position
####################

lats = np.linspace(np.pi/2,-np.pi/2,num=100)
lons = np.linspace(-np.pi,np.pi,num=100)
lons, lats = np.meshgrid(lons, lats)

lats = np.rad2deg(lats)
lons = np.rad2deg(lons)

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(VGP_C2_lon[44:91], VGP_C2_lat[44:91])

sc = map1.scatter(xs,ys, marker='D',c=Age_C2[44:91]/1000, cmap='rainbow')
#plt.gray()

plt.colorbar(sc,fraction=0.03, pad=0.04)
plt.title('Sanxing VGPs', y=1.08)




#########################################
########################################
# Sheep Creek, Bogue et al., 2010
########################################
########################################

lat_SC = 40.73
lon_SC = 243.15

# 1 degree/week (actually 53 degrees in a period between 12 and 18 month, from the paper)
dangle_dt_Bogue = 1*52 # 52 weeks in year
age_Bogue=15580

# indicative (not real, I and D before and after the transitional times)
Incl_SC = [64, 62] 
Decl_SC = [357, 107]


# calculate VGP
VGP_SC_lat, VGP_SC_lon = subs.VGP_from_DI(Incl_SC,Decl_SC,lat_SC,lon_SC)
dVGPlat_dt_SC = np.abs( ( VGP_SC_lat[-1] - VGP_SC_lat[0] ) / 1 ) # this happened in about 1 year



#############################################################
#############################################################
# Nowaczyk 2012: LAschamp excursion from black sea sediments
#############################################################
#############################################################


########################################
# Read in Nowaczyk 2012, figure 9 data
# M72/5-22 stack
########################################

lat_M72_5_22 = 42.225500 
lon_M72_5_22 = 36.492500

BS_folder = '../datasets/Nowaczyk_2012/'
data_M72_5_22_file = BS_folder + 'VGP_lat_stack_fig9.csv'

data_M72_5_22 = []
age_M72_5_22 =[]
VGP_lat_M72_5_22 = []
with open(data_M72_5_22_file) as f1:
    readCSV = csv.reader(f1,delimiter=',')
    for row in readCSV:
        data_M72_5_22.append(row)
        # get depth data point
        try:
            age = float(row[0])
        except:
            # not a string: probably still in the header
            continue
        age_M72_5_22 = np.append(age_M72_5_22,age)
        # get inclination data point
        try:
            VGP_lat = float(row[1])
        except:
            VGP_lat = np.nan            
        VGP_lat_M72_5_22 = np.append(VGP_lat_M72_5_22,VGP_lat)

f1.close()
    


####################################################################
# plot VGP lat time series (should look like fig 9 of the paper)
####################################################################
fig,ax = plt.subplots(figsize=(8,5))

# add vertical lines for the transitional period
#ax_twin.plot([times[twrite], times[twrite]],[0, 40],'--',color='lightgray')
#ax_twin.plot([times[tend], times[tend]],[0, 40],'--',color='lightgray')
#ax_twin.set_ylim(0, 40)

ax.set_xlabel('Age / ka')
ax.set_ylabel('VGP Latitude/$^\circ$',color='k')

ax.plot(age_M72_5_22,VGP_lat_M72_5_22,'o-',color='k')
#ax.set_xlim(times[0], times[-1])

#ax.tick_params('y', colors='k')
plt.show
plt.title('VGP latitude for M72_5_22 samples during Laschamp (Nowaczyk et al., 2012)')

# get fast variation
dVGP_lat_dt_M72_5_22 = np.abs( (VGP_lat_M72_5_22[75] - VGP_lat_M72_5_22[68]) / (1000*age_M72_5_22[75] - 1000*age_M72_5_22[68]) )
age_fast_M72_5_22   = (age_M72_5_22[75] + age_M72_5_22[68]) / 2.0 # central age for it (in yr)


#########################################
########################################
# LSMOD, Laschamp
########################################
########################################

# characteristic rapid dipole latitude variation in deg/yr
dtheta_dt_LSMOD = np.abs( (74.1- 14.6) / (41000 - 41115) )
age_dthetadt_LSMOD   = (41.000 + 41.115) / 2.0 # central age for it (in yr)

# VGP lat at Black Sea in deg / yr
dVGPlat_dt_BS_LSMOD = np.abs( (87.41- (-25.08)) / (40050 - 41200) )
age_dVGPlat_dt_BS_LSMOD   = (40.05 + 41.2) / 2.0 # central age for it (in yr)

# optimal solution
dVGPlat_dt_BS_LSMOD_opt = 9.36626747722275
age_dVGPlat_dt_BS_LSMOD_opt = 41250

# optimal solution if we optimise dIdt
I_opt_M72_5_22_dVGPlat_dt = 3.7299165168999999
age_I_opt_M72_5_22 = 41100

# optimal dtheta_dt
LSMOD_tilt_opt = 2.8668367
age_LSMOD_tilt_opt = 41000

#########################################
########################################
# Baseline: Cals10k
########################################
########################################

line = []
coeffs_B = []
with open('../models/CALS10K1B/cals10k1b_gauss.txt', 'rb') as f:
    while True:
        line = f.readline()
        if not line:
            break
        x = [list(map(float, line.split() ))]
        coeffs_B = np.append(coeffs_B,x)
cols = len(x[0])
rows = len(coeffs_B)
CALS10K1B_Coeffs = np.reshape(coeffs_B,(rows/cols, cols))

time_CALS10K1B = CALS10K1B_Coeffs[:,0]
age_CALS10K1B  = year_now-time_CALS10K1B
g10_CALS10K1B = CALS10K1B_Coeffs[:,1]
g11_CALS10K1B = CALS10K1B_Coeffs[:,2]
h11_CALS10K1B = CALS10K1B_Coeffs[:,3]

# calculate dipole tilt timeseries, in deg
CALS10K1B_theta_d = np.arccos(g10_CALS10K1B/np.sqrt(g10_CALS10K1B**2+g11_CALS10K1B**2+h11_CALS10K1B**2))*180./np.pi
# calculate dipole tilt time derivative (deg/yr) wihth first differences
CALS10K1B_dtheta_d_dt = (CALS10K1B_theta_d[1:] - CALS10K1B_theta_d[0:-1])/ (time_CALS10K1B[1:] -  time_CALS10K1B[0:-1])
# get RMS value
CALS10K1B_dtheta_d_dt_RMS = np.sqrt(np.mean(CALS10K1B_dtheta_d_dt**2))


# plot it
fig,ax = plt.subplots(figsize=(8,5))

ax.set_xlabel('Age / ka')
ax.set_ylabel('$^\circ$',color='k')

ax.plot(age_CALS10K1B,CALS10K1B_theta_d,'-',color='k')
#ax.set_xlim(times[0], times[-1])

#ax.tick_params('y', colors='k')
plt.show
plt.title('CALS10K1b dipole tilt')


###############################################################################

###############################
###############################
# Summary plot
###############################
###############################

'''
# need to add: 

Directional change during a Miocene R‐N geomagnetic polarity reversal recorded by mafic lava flows, Sheep Creek Range, north central Nevada, USA
S. W. Bogue  J. M. G. Glen  N. A. Jarboe
https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017GC007049

Jarboe, 2010

Baseline, maybe from  CALS100k (as a continuous line?)

# use different filling/sybols to indicate the degree of acceptance of different datasets.
'''

# plot parmeters
CS  = 10    # capsize
ELW = 2     # elinewidth
MS  = 40    # s (marker size)
AY = 1.1    # annotate Y position factor
AYL = 0.6   # for annotation below the marker
DX = 30     # annotate X position displacement
#fig,(ax,ax2) = plt.subplots(1,2,sharey=True,figsize=(8,5))
fig,(ax,ax2) = plt.subplots(1,2,sharey=True,figsize=(8,8))



# Chou et al. 2018
# plot temporal span of dataset
#ax.errorbar((np.nanmin(Age_C2) + np.nanmax(Age_C2))/2/1000,dVGPlat_dt_C2,
#            xerr=(np.nanmin(Age_C2) - np.nanmax(Age_C2))/2/1000, 
#            color='b', ecolor='blue', capsize=CS, elinewidth=ELW)

# plot central instant of max variation
ax.scatter(age_fast_C2/1000,dVGPlat_dt_C2,
           marker='^', s=MS, color='k', edgecolors='k')
ax.annotate('Sanxing Cave', (DX+age_fast_C2/1000,dVGPlat_dt_C2*(AYL+0.2) ) )

# Sagnotti et al., 2015
# plot temporal span of dataset
#ax.errorbar((np.nanmin(time_SUL) + np.nanmax(time_SUL))/2,dVGPlat_dt_SUL,
#            xerr=(np.nanmin(time_SUL) - np.nanmax(time_SUL))/2, 
#            color='r', ecolor='r', capsize=CS, elinewidth=ELW)
# plot central instant of max variation
ax.scatter(age_fast_SUL/1000,dVGPlat_dt_SUL,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax.annotate('Sulmona (2014)', (DX+age_fast_SUL/1000,dVGPlat_dt_SUL*AYL ) )

# Sagnotti et al., 2016
# plot temporal span of dataset
#ax.errorbar((np.nanmin(time_SUL) + np.nanmax(time_SUL))/2,dVGPlat_dt_SUL_2,
#            xerr=(np.nanmin(time_SUL) - np.nanmax(time_SUL))/2, 
#            color='orange', ecolor='orange', capsize=CS, elinewidth=ELW)
# plot central instant of max variation
ax.scatter(age_fast_SUL/1000,dVGPlat_dt_SUL_2,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax.annotate('Sulmona (2015)', (DX+age_fast_SUL/1000,dVGPlat_dt_SUL_2*AY ) )

# IMMAB4, VGP at SUlmona
# plot central instant of max variation
ax.scatter(age_dVGPlat_dt_SUL_IMMAB4/1000,dVGPlat_dt_SUL_IMMAB4,
           marker='^', s=MS, color='w', edgecolors='k')
ax.annotate('Sulmona (IMMAB4)', (DX+age_dVGPlat_dt_SUL_IMMAB4/1000,dVGPlat_dt_SUL_IMMAB4*AY ) )

# IMMAB4, dtheta/dt
# plot central instant of max variation
ax.scatter(age_dthetadt_IMMAB4/1000,dtheta_dt_IMMAB4,
           marker='d', s=MS, color='w', edgecolors='k')
ax.annotate('IMMAB4', (DX+age_dVGPlat_dt_SUL_IMMAB4/1000,dtheta_dt_IMMAB4*AY ) )

# Coe, 1995

# plot central instant of max variation
# this is just th e6 degrees per day indicated in the paper: need to translate it to a VGP lat variation
dangle_dt_Coe = 6*356
age_Coe=16200
ax2.scatter(age_Coe,dangle_dt_Coe,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax2.annotate('Steens Mountains', (15750,dangle_dt_Coe*AY ) )

# Bogue et al., 2010. Sheep Creek (thermal history different from Steens mountain)
# problem with this point is that I used "indicative" values for initial and final Incl and Decl
ax2.scatter(age_Bogue,dVGPlat_dt_SC,
           marker='^', s=MS, color='k', edgecolors='k')
ax2.annotate('Sheep Creek', (DX/5+age_Bogue,dVGPlat_dt_SC*AY ) )

# Nowaczyk 2012. Black Sea Laschamp excursion
ax.scatter(age_fast_M72_5_22,dVGP_lat_dt_M72_5_22,
           marker='^', s=MS, color='k', edgecolors='k')
ax.annotate('Black Sea', (DX+age_fast_M72_5_22,dVGP_lat_dt_M72_5_22*AY ) )

# LSMOD. Laschamp excursion
ax.scatter(age_dthetadt_LSMOD,dtheta_dt_LSMOD,
           marker='d', s=MS, color='white', edgecolors='k')
ax.annotate('LSMOD1', (DX+age_dthetadt_LSMOD,dtheta_dt_LSMOD*(AYL) ) )

# LSMOD. Laschamp excursion
ax.scatter(age_dVGPlat_dt_BS_LSMOD,dVGPlat_dt_BS_LSMOD,
           marker='^', s=MS, color='white', edgecolors='k')
ax.annotate('Black Sea \n (LSMOD1)', (DX+age_dVGPlat_dt_BS_LSMOD,dVGPlat_dt_BS_LSMOD*(AYL-0.2) ) )


# CALS10K1B: baseline
ax.plot([0, 50000],[CALS10K1B_dtheta_d_dt_RMS, CALS10K1B_dtheta_d_dt_RMS],'--',color='brown')
ax2.plot([0, 50000],[CALS10K1B_dtheta_d_dt_RMS, CALS10K1B_dtheta_d_dt_RMS],'--',color='brown')
ax2.annotate('CALS10K1B', (16000,CALS10K1B_dtheta_d_dt_RMS*(AY)), color='brown')

# fake points for a legend
l1 = ax.scatter(-1000,-1000,
           marker='d', s=MS, color='w', edgecolors='k',
            label='Dipole tilt')
l2 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='w', edgecolors='k',
            label='VGP latitude')
l3 = ax.scatter(-1000,-1000,
           marker='v', s=MS, color='w', edgecolors='k',
            label='VGP latitude (from $dI/dt$)')

l4 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='k', edgecolors='k')
l5 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='lightgray', edgecolors='k')
           
# create axis spines
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
ax.tick_params(labelright='off')
ax2.yaxis.tick_right()

# draw spines at the break of the x-axis
d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d), (-d,+d), **kwargs)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)

# final touches and legend
ax2.set_xlim(15500, 16500)
ax.set_xlim(0, 1400)

ax.set_xlabel('Age / ka')
#ax.set_ylabel('$d\lambda / dt (^\circ / yr)$')
ax.set_ylabel('$^\circ / yr$')
ax.xaxis.set_label_coords(0.5, 0.05, transform=fig.transFigure)
 
ax.set_ylim(0.001, 10000)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.tick_params('y')
fig.suptitle('Rapid paleomagnetic variations')
#ax.legend(fontsize=9,loc='upper left',title='Rate of change of')
leg1=ax.legend([(l1),(l2, l4, l5)], [l1._label,l2._label], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=9,loc='upper left',title='Rate of change of')

leg2=ax.legend([(l1,l2), (l4), (l5)], ['SH models','data','contested data','optimal solution'], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=9,loc=(0.025, 0.736),title='Source:')
ax.add_artist(leg1)


plt.savefig('Rapid_variations_summary_data.pdf',bbox_inches='tight',pad_inches=0.0)
leg1.remove()
leg2.remove()

# LSMOD optimal solution
ax.scatter(age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt,
           marker='^', s=MS, color='r', edgecolors='k')
ax.annotate('Black Sea (LSMOD1)', (-DX+age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt*(AY+0.1) ) ,color='red')
# VGP solution for optimal dIdt
ax.scatter(age_I_opt_M72_5_22/1000,I_opt_M72_5_22_dVGPlat_dt,
           marker='v', s=MS, color='r', edgecolors='k')
ax.annotate('Black Sea (LSMOD1)', (-DX+age_I_opt_M72_5_22/1000,I_opt_M72_5_22_dVGPlat_dt*(AY+0.1) ) ,color='red')
# optimal dipole tilt rate of change
ax.scatter(age_LSMOD_tilt_opt/1000,LSMOD_tilt_opt,
           marker='d', s=MS, color='r', edgecolors='k')
ax.annotate('LSMOD1', (DX/2+age_LSMOD_tilt_opt/1000,LSMOD_tilt_opt*AYL ) ,color='red' )
#ax.annotate('Optimal \n (LSMOD1)', (DX+age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt*AY ),
#            color='red')

# IMMAB4, optimal
# plot central instant of max variation
ax.scatter(age_dVGPlat_dt_SUL_IMMAB4_opt/1000,dVGPlat_dt_SUL_IMMAB4_opt,
           marker='^', s=MS, color='r', edgecolors='k')
ax.annotate('Sulmona (IMMAB4)', (DX+age_dVGPlat_dt_SUL_IMMAB4_opt/1000,dVGPlat_dt_SUL_IMMAB4_opt*AY ) ,color='red')
# VGP solution for max dIdt
ax.scatter(age_I_opt_SUL/1000,I_opt_SUL_dVGPlat_dt,
           marker='v', s=MS, color='r', edgecolors='k')
ax.annotate('Sulmona (IMMAB4)', (DX+age_I_opt_SUL/1000,I_opt_SUL_dVGPlat_dt*AY )  ,color='red')
# optimal dipole tilt rate of change
ax.scatter(age_IMMAB4_tilt_opt/1000,IMMAB4_tilt_opt,
           marker='d', s=MS, color='r', edgecolors='k')
ax.annotate('IMMAB4', (DX+age_IMMAB4_tilt_opt/1000,IMMAB4_tilt_opt*AY )  ,color='red')

# add optimal solution to legend
l6 = ax.scatter(-1000,-1000,
           marker='d', s=MS, color='r', edgecolors='k')
l7 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='r', edgecolors='k')
l8 = ax.scatter(-1000,-1000,
           marker='v', s=MS, color='r', edgecolors='k')

leg1=ax.legend([(l1, l6),(l2, l4, l5, l7),(l8)], [l1._label,l2._label,l3._label], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=9,loc='upper left',title='Rate of change of')

leg2=ax.legend([(l1,l2), (l4), (l5), (l6,l7,l8)], ['SH models','data','contested data','optimal solution'], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=9,loc=(0.025, 0.67),title='Source:')

ax.add_artist(leg1)

plt.savefig('Rapid_variations_summary.pdf',bbox_inches='tight',pad_inches=0.0)
plt.show







###############################
###############################
# Summary plot, log log version
###############################
###############################

'''
# need to add: 

Directional change during a Miocene R‐N geomagnetic polarity reversal recorded by mafic lava flows, Sheep Creek Range, north central Nevada, USA
S. W. Bogue  J. M. G. Glen  N. A. Jarboe
https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017GC007049

Jarboe, 2010

Baseline, maybe from  CALS100k (as a continuous line?)

# use different filling/sybols to indicate the degree of acceptance of different datasets.
'''

# plot parmeters
CS  = 10    # capsize
ELW = 2     # elinewidth
MS  = 40    # s (marker size)
AY = 1.1    # annotate Y position factor
AYL = 0.6   # for annotation below the marker
DX = 30     # annotate X position displacement
AFS = 12     # annotation font size
LFS = 10    # legend font size
#fig,(ax,ax2) = plt.subplots(1,2,sharey=True,figsize=(8,5))
fig,ax = plt.subplots(1,1,sharey=True,figsize=(8,8))

#plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=AFS)     # fontsize of the axes title
plt.rc('axes', labelsize=AFS)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=AFS)    # fontsize of the tick labels
plt.rc('ytick', labelsize=AFS)    # fontsize of the tick labels
#plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# Chou et al. 2018
# plot temporal span of dataset
#ax.errorbar((np.nanmin(Age_C2) + np.nanmax(Age_C2))/2/1000,dVGPlat_dt_C2,
#            xerr=(np.nanmin(Age_C2) - np.nanmax(Age_C2))/2/1000, 
#            color='b', ecolor='blue', capsize=CS, elinewidth=ELW)

# plot central instant of max variation
ax.scatter(age_fast_C2/1000,dVGPlat_dt_C2,
           marker='^', s=MS, color='k', edgecolors='k')
ax.annotate('Sanxing Cave', (AY*age_fast_C2/1000,dVGPlat_dt_C2*(AYL) ), 
            fontsize=AFS )

# Sagnotti et al., 2015
# plot temporal span of dataset
#ax.errorbar((np.nanmin(time_SUL) + np.nanmax(time_SUL))/2,dVGPlat_dt_SUL,
#            xerr=(np.nanmin(time_SUL) - np.nanmax(time_SUL))/2, 
#            color='r', ecolor='r', capsize=CS, elinewidth=ELW)
# plot central instant of max variation
ax.scatter(age_fast_SUL/1000,dVGPlat_dt_SUL,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax.annotate('Sulmona (2014)', (AY*age_fast_SUL/1000,dVGPlat_dt_SUL*AYL ), 
            fontsize=AFS  )

# Sagnotti et al., 2016
# plot temporal span of dataset
#ax.errorbar((np.nanmin(time_SUL) + np.nanmax(time_SUL))/2,dVGPlat_dt_SUL_2,
#            xerr=(np.nanmin(time_SUL) - np.nanmax(time_SUL))/2, 
#            color='orange', ecolor='orange', capsize=CS, elinewidth=ELW)
# plot central instant of max variation
ax.scatter(age_fast_SUL/1000,dVGPlat_dt_SUL_2,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax.annotate('Sulmona (2015)', (AY*age_fast_SUL/1000,dVGPlat_dt_SUL_2*AY ), 
            fontsize=AFS  )

# IMMAB4, VGP at SUlmona
# plot central instant of max variation
ax.scatter(age_dVGPlat_dt_SUL_IMMAB4/1000,dVGPlat_dt_SUL_IMMAB4,
           marker='^', s=MS, color='w', edgecolors='k')
ax.annotate('Sulmona (IMMAB4)', (AY*age_dVGPlat_dt_SUL_IMMAB4/1000,dVGPlat_dt_SUL_IMMAB4*AY ), 
            fontsize=AFS  )

# IMMAB4, dtheta/dt
# plot central instant of max variation
ax.scatter(age_dthetadt_IMMAB4/1000,dtheta_dt_IMMAB4,
           marker='d', s=MS, color='w', edgecolors='k')
ax.annotate('IMMAB4', (AY*age_dVGPlat_dt_SUL_IMMAB4/1000,dtheta_dt_IMMAB4*AY ), 
            fontsize=AFS  )

# Coe, 1995

# plot central instant of max variation
# this is just th e6 degrees per day indicated in the paper: need to translate it to a VGP lat variation
dangle_dt_Coe = 6*356
age_Coe=16200
ax.scatter(age_Coe,dangle_dt_Coe,
           marker='^', s=MS, color='lightgray', edgecolors='k')
ax.annotate('Steens Mountains', (age_Coe/12,dangle_dt_Coe*AY ), 
            fontsize=AFS  )

# Bogue et al., 2010. Sheep Creek (thermal history different from Steens mountain)
# problem with this point is that I used "indicative" values for initial and final Incl and Decl
ax.scatter(age_Bogue,dVGPlat_dt_SC,
           marker='^', s=MS, color='k', edgecolors='k')
ax.annotate('Sheep Creek', (age_Bogue/6,dVGPlat_dt_SC*AY ), 
            fontsize=AFS  )

# Nowaczyk 2012. Black Sea Laschamp excursion
ax.scatter(age_fast_M72_5_22,dVGP_lat_dt_M72_5_22,
           marker='^', s=MS, color='k', edgecolors='k')
ax.annotate('Black Sea', (AY*age_fast_M72_5_22,dVGP_lat_dt_M72_5_22*AY ), 
            fontsize=AFS  )

# LSMOD. Laschamp excursion
ax.scatter(age_dthetadt_LSMOD,dtheta_dt_LSMOD,
           marker='d', s=MS, color='white', edgecolors='k')
ax.annotate('LSMOD1', (AY*age_dthetadt_LSMOD,dtheta_dt_LSMOD*(AYL) ), 
            fontsize=AFS  )

# LSMOD. Laschamp excursion
ax.scatter(age_dVGPlat_dt_BS_LSMOD,dVGPlat_dt_BS_LSMOD,
           marker='^', s=MS, color='white', edgecolors='k')
ax.annotate('Black Sea \n (LSMOD1)', (AY*age_dVGPlat_dt_BS_LSMOD,dVGPlat_dt_BS_LSMOD*(AYL-0.2) ), 
            fontsize=AFS  )


# CALS10K1B: average
#ax.scatter(np.mean([age_CALS10K1B[-1]/1000., age_CALS10K1B[0]/1000.]),CALS10K1B_dtheta_d_dt_RMS,
 #          marker='d', s=MS, color='white', edgecolors='k')
ax.annotate('CALS10k.1b (rms) ', (0.5*np.mean([age_CALS10K1B[-1]/1000., age_CALS10K1B[0]/1000.]),CALS10K1B_dtheta_d_dt_RMS*(AYL) ),
            fontsize = AFS)
ax.plot([age_CALS10K1B[-1]/1000., age_CALS10K1B[0]/1000.],[CALS10K1B_dtheta_d_dt_RMS, CALS10K1B_dtheta_d_dt_RMS],'k--')
# fake points for a legend
l1 = ax.scatter(-1000,-1000,
           marker='d', s=MS, color='w', edgecolors='k',
            label='Dipole tilt')
l2 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='w', edgecolors='k',
            label='VGP latitude')
l3 = ax.scatter(-1000,-1000,
           marker='v', s=MS, color='w', edgecolors='k',
            label='VGP latitude (from $dI/dt$)')

l4 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='k', edgecolors='k')
l5 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='lightgray', edgecolors='k')
           

ax.set_xlabel('Age / ka')
#ax.set_ylabel('$d\lambda / dt (^\circ / yr)$')
ax.set_ylabel('$^\circ / yr$')
#ax.xaxis.set_label_coords(0.5, 0.05, transform=fig.transFigure)
 
ax.set_ylim(0.001, 10000)
ax.set_xlim(1, 20000)
ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params('y')
fig.suptitle('Rapid paleomagnetic variations')
#ax.legend(fontsize=9,loc='upper left',title='Rate of change of')
leg1=ax.legend([(l1),(l2, l4, l5)], [l1._label,l2._label], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=LFS,loc='upper left',title='Rate of change of')

leg2=ax.legend([(l1,l2), (l4), (l5)], ['SH models','data','contested data','optimal solution'], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=LFS,loc=(0.011, 0.71),title='Source:')
ax.add_artist(leg1)


plt.savefig('Rapid_variations_summary_data_loglog.pdf',bbox_inches='tight',pad_inches=0.0)
leg1.remove()
leg2.remove()

# LSMOD optimal solution
ax.scatter(age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt,
           marker='^', s=MS, color='r', edgecolors='k')
ax.annotate('Black Sea (LSMOD1)', (0.06*age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt*(AY+0.1) ) ,
            color='red',
            fontsize=AFS)
# VGP solution for optimal dIdt
#ax.scatter(age_I_opt_M72_5_22/1000,I_opt_M72_5_22_dVGPlat_dt,
#           marker='v', s=MS, color='r', edgecolors='k')
#ax.annotate('Black Sea (LSMOD1)', (0.06*age_I_opt_M72_5_22/1000,I_opt_M72_5_22_dVGPlat_dt*(AY+0.1) ) ,
#            color='red',
#            fontsize=AFS)
# optimal dipole tilt rate of change
ax.scatter(age_LSMOD_tilt_opt/1000,LSMOD_tilt_opt,
           marker='d', s=MS, color='r', edgecolors='k')
ax.annotate('LSMOD1', (AY*age_LSMOD_tilt_opt/1000,LSMOD_tilt_opt*AYL ) ,
            color='red',
            fontsize=AFS )
#ax.annotate('Optimal \n (LSMOD1)', (AY*age_dVGPlat_dt_BS_LSMOD_opt/1000,dVGPlat_dt_BS_LSMOD_opt*AY ),
#            color='red')

# IMMAB4, optimal
# plot central instant of max variation
ax.scatter(age_dVGPlat_dt_SUL_IMMAB4_opt/1000,dVGPlat_dt_SUL_IMMAB4_opt,
           marker='^', s=MS, color='r', edgecolors='k')
ax.annotate('Sulmona (IMMAB4)', (AY*age_dVGPlat_dt_SUL_IMMAB4_opt/1000,dVGPlat_dt_SUL_IMMAB4_opt*AY ) ,
            color='red',
            fontsize=AFS)
# VGP solution for max dIdt
#ax.scatter(age_I_opt_SUL/1000,I_opt_SUL_dVGPlat_dt,
#           marker='v', s=MS, color='r', edgecolors='k')
#ax.annotate('Sulmona (IMMAB4)', (AY*age_I_opt_SUL/1000,I_opt_SUL_dVGPlat_dt*AY )  ,
#            color='red',
#            fontsize=AFS)
# optimal dipole tilt rate of change
ax.scatter(age_IMMAB4_tilt_opt/1000,IMMAB4_tilt_opt,
           marker='d', s=MS, color='r', edgecolors='k')
ax.annotate('IMMAB4', (AY*age_IMMAB4_tilt_opt/1000,IMMAB4_tilt_opt*AY )  ,
            color='red',
            fontsize=AFS)

# add optimal solution to legend
l6 = ax.scatter(-1000,-1000,
           marker='d', s=MS, color='r', edgecolors='k')
l7 = ax.scatter(-1000,-1000,
           marker='^', s=MS, color='r', edgecolors='k')
l8 = ax.scatter(-1000,-1000,
           marker='v', s=MS, color='r', edgecolors='k')

#leg1=ax.legend([(l1, l6),(l2, l4, l5, l7),(l8)], [l1._label,l2._label,l3._label], 
#               numpoints=1,
#               handler_map={tuple: HandlerTuple(ndivide=None)},
#               fontsize=LFS,loc='upper left',title='Rate of change of')
leg1=ax.legend([(l1, l6),(l2, l4, l5, l7)], [l1._label,l2._label], 
               numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=LFS,loc='upper left',title='Rate of change of')


leg2=ax.legend([(l1,l2), (l4), (l5), (l6,l7,l8)], ['SH models','data','contested data','optimal solution'], numpoints=1,
               handler_map={tuple: HandlerTuple(ndivide=None)},
               fontsize=LFS,
#               loc=(0.011, 0.63),
               loc=(0.011, 0.67),
               title='Source:')

ax.add_artist(leg1)

plt.savefig('Rapid_variations_summary_loglog.pdf',bbox_inches='tight',pad_inches=0.0)
plt.show