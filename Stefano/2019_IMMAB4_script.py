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

sys.path.append('./SpecialFunctions/SphericalHarmonics/')
import SH_library

import subs

import matplotlib.tri as tri

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap



# inputs for animation
animation_flag = 0 # 0: no animation; 1: create animation  
yps = 10 # years per second of animation 

# other inputs
twrite = 122 # the snapshot to be stored in file
tend   = 231
##################################################
# In Figure 8 of Leonhardt 2007 the start of the 
# reversal is between t = 781.8 Myr
# and t = 770.9
# I took these values visually from the plot.
##################################################

# Sulmona colatitute and longitude
lat_SUL = 90-47.8488
lon_SUL   = 13.8229 
# Sidney
lat_Sid = -33.8688
lon_Sid = 151.2093
#Quito
lat_Quito = -0.1807
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

figname = 'IMMAB4'
folder = '../models/IMMAB4/'
figfolder = '../models/IMMAB4/figures/'

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
# coefficients are in microTeslas!!!!
#####################################    

n_instants=coeffs_B.shape[0]
times = coeffs_B[:,0]   

n_coeffs=coeffs_B.shape[1]-1

Br_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Br_c = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Bt_a = np.zeros((n_instants, theta.shape[0], theta.shape[1])) 
Bp_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))

Br_locations = np.zeros((n_instants, theta_locations.size))
Bt_locations = np.zeros((n_instants, theta_locations.size))
Bp_locations = np.zeros((n_instants, theta_locations.size))

r_cmb = 3485.0e3
r_e   = 6371.0e3

it = 0
l=1
m=0
cs='c'

for ic in range(n_coeffs):

    betalm = coeffs_B[:,ic+1]
    
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
        Br_a[it] = Br_a[it] + (l+1) * (r_e/r_e)**(l+2) * betalm[it] * SH
        Br_c[it] = Br_c[it] + (l+1) * (r_e/r_cmb)**(l+2) * betalm[it] * SH
        Bt_a[it] = Bt_a[it] - (r_e/r_e)**(l+2) * betalm[it] * dthSH
        Bp_a[it] = Bp_a[it] - (r_e/r_e)**(l+2) * betalm[it] * np.divide(dphSH,np.sin(theta))
    
        Br_locations[it] = Br_locations[it] + (l+1) * (r_e/r_e)**(l+2) * betalm[it] * SH_locations
        Bt_locations[it] = Bt_locations[it] - (r_e/r_e)**(l+2) * betalm[it] * dthSH_locations
        Bp_locations[it] = Bp_locations[it] - (r_e/r_e)**(l+2) * betalm[it] * np.divide(dphSH_locations,np.sin(theta_locations))
        
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


########################################
# Read in Sagnotti's 2014 data
########################################
SUL_folder = '/nfs/see-fs-01_users/earsmaf/geomag/datasets/Sagnotti_2014/'
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

#######
#######
# plots
#######
#######

######################
# cmb radial field
######################

it = tend

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
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
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
        plt.title(times[it], y=1.08)
        plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'Br_CMB_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
        dt = times[1]-times[0] #assuming a constant timestep
        
    fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'Br_CMB_%07d.png -vcodec libx264 -y -an ' +folder+ 'Br_CMB.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
    #os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')


######################
# surface radial field
######################
# plot to check things

it = tend

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
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
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
        plt.title(times[it], y=1.08)
        plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'Br_a_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
        dt = times[1]-times[0] #assuming a constant timestep
        
    fps = int(yps/dt)
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
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
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
        plt.title(times[it], y=1.08)
        plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)
    
        # saving the plot
        fname = folder +'F_surface_%07d.png' % it
        plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
        files.append(fname)
        plt.close()
    
    
    dt = times[1]-times[0] #assuming a constant timestep
    fps = int(yps/dt)
    fps = 20
    os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'F_surface_%07d.png -vcodec libx264 -y -an ' +folder+ 'F_surface.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
    #os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')




######################
# time series
#####################

# Dipole tilt

fig,ax = plt.subplots(figsize=(8,5))
ax_twin = ax.twinx()

# add vertical lines for the transitional period
ax_twin.plot([times[twrite], times[twrite]],[0, 40],'--',color='gray')
ax_twin.plot([times[tend], times[tend]],[0, 40],'--',color='gray')
ax_twin.set_ylim(0, 40)

ax.set_xlabel('Time')
ax.set_ylabel('Dipole latitude/$^\circ$',color='b')

ax.plot(times,dipole_colatitude,color='b')
ax.set_xlim(times[0], times[-1])
ax_twin.plot(times,m1,color='r')
ax_twin.set_ylabel('Dipole intensity/$\mu$T',color='r')
ax_twin.tick_params('y', colors='r')
ax.tick_params('y', colors='b')

plt.title('Dipole latitude and intensity evolution')
plt.savefig(folder+'figures/Dipole_tilt.pdf',bbox_inches='tight',pad_inches=0.0)


# Dipole intensity

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([times[twrite], times[twrite]],[-40, 40],'--',color='gray')
ax_i.plot([times[tend], times[tend]],[-40, 40],'--',color='gray')
ax_i.set_ylim(-40, 40)

ax_i.plot(times,g10,color='b')
ax_i.set_xlim(times[0], times[-1])
ax_i.legend(fontsize=10,loc='lower left')
plt.title('$g_1^0$')
plt.show
plt.savefig(folder+'figures/g10.pdf',bbox_inches='tight',pad_inches=0.0)



# Inclination time-series
inclination_locations = np.arctan(np.divide(-Br_locations,np.sqrt(Bt_locations**2+Bp_locations**2)))*180/np.pi

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([times[twrite], times[twrite]],[-90, 90],'--',color='gray')
ax_i.plot([times[tend], times[tend]],[-90, 90],'--',color='gray')
ax_i.set_ylim(-90, 90)

ax_i.plot(times,inclination_locations[:,1],color='r',
          label='Sulmona'  )
ax_i.plot(times,inclination_locations[:,0],color='b',
          label='Quito'  )
ax_i.plot(times,inclination_locations[:,2],color='g',
          label='Sidney'  )

ax_i.set_xlim(times[0], times[-1])
ax_i.legend(fontsize=10,loc='upper left')
plt.title('Inclination at locations')
plt.show
plt.savefig(folder+'figures/inclinations.pdf',bbox_inches='tight',pad_inches=0.0)



# SUL Inclination time-series. IMMAB4 vs Sagnotti
inclination_locations = np.arctan(np.divide(-Br_locations,np.sqrt(Bt_locations**2+Bp_locations**2)))*180/np.pi

fig_i,ax_i =  plt.subplots(figsize=(8,5))
# add vertical lines for the transitional period
ax_i.plot([times[twrite], times[twrite]],[-90, 90],'--',color='gray')
ax_i.plot([times[tend], times[tend]],[-90, 90],'--',color='gray')
ax_i.set_ylim(-90, 90)

ax_i.plot(times,inclination_locations[:,1],color='r',
          label='IMMAB4'  )

ax_i.plot(time_SUL,Incl_SUL,'--',color='r',
          label='Sagnotti et al., 2014'  )

ax_i.plot(time_SUL-13.95,Incl_SUL,'--',color='c',
          label='Sagnotti et al., 2014 (shifted)'  )

ax_i.set_xlim(times[0], times[-1])
ax_i.legend(fontsize=10,loc='bottom right')
plt.title('Inclination at SUL')
plt.show
plt.savefig(folder+'figures/inclinations_SUL.pdf',bbox_inches='tight',pad_inches=0.0)

