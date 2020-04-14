#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 15:42:00 2020

@author: Stefano

Plot some magnetic field quantities, to play around with

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


lats = np.linspace(np.pi/2,-np.pi/2,num=100)
lons = np.linspace(-np.pi,np.pi,num=100)
lons, lats = np.meshgrid(lons, lats)
theta = -lats+np.pi/2

## get the chaos model
CHAOS_2019 = '../models/CHAOS-6/CHAOS-6-x9_core_2019.0007.dat' 

coeffs_CHAOS=subs.read_coeffs(CHAOS_2019,1,1)
coeffs_CHAOS[0,2]=-coeffs_CHAOS[0,2]
Br, Bth, Bphi = SH_library.calcB(coeffs_CHAOS,theta,lons,r_a,r_a)

H = np.sqrt( Bth**2 + Bphi**2 )

Incl = np.arctan(-Br/H)
Decl = np.arctan(-Bphi/Bth)




#######################
# plots
#######################

lats = np.rad2deg(lats)
lons = np.rad2deg(lons)


# plot Br
plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Br, cmap='coolwarm')
cs = map1.contour(xs,ys,Br, 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title('$B_r$', y=1.08)
#plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)




# plot Btheta
plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Bth, cmap='coolwarm')
cs = map1.contour(xs,ys,Bth, 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title('$B_{\theta}$', y=1.08)
#plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)




# polt Inclination
plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,np.rad2deg(Incl), cmap='coolwarm')
cs = map1.contour(xs,ys,np.rad2deg(Incl), 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title('Inclination', y=1.08)
#plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)




# polt Declination
plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,np.rad2deg(Decl), cmap='coolwarm')
cs = map1.contour(xs,ys,np.rad2deg(Decl), 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title('Declination', y=1.08)
#plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)