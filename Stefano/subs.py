#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 18:26:15 2019

@author: earsmaf
"""
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.special import lpmv # associate legendre functions, in vector form
import sys
import os
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, './SpecialFunctions/SphericalHarmonics/')
from SH_library import DthSchmidtSH
from SH_library import DphSchmidtSH
from SH_library import SchmidtSH
from SH_library import mat2linCoeffs
from SH_library import calcB
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc
import math
import cartopy.crs as ccrs

c = 3485.0e3

"""
def plot_mollweide_contourf(lon,lat,data,colormap):

    input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
    output_file = folder+'flow_mollweide.pdf'
    
    X,Y,U,V = np.loadtxt(input_file,unpack=True)
    
    # display max flow on input grid
    USQ = []
    for i in range(len(U)):
        USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
    print ('Largest value of |u| on grid is ', max(USQ))
    long = X
    lat = Y
    
    plt.figure(1, figsize=(10,8))
    
    
    #Mollweide projectionfrom scipy.interpolate import griddata
    map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    
    # draw coastlines, country boundaries, fill continents.
    map1.drawcoastlines(linewidth=0.50)
    
    # draw lat/lon grid lines every 30 degrees.
    map1.drawmeridians(np.arange(0,360,30),linewidth=1)
    map1.drawparallels(np.arange(-90,90,30),linewidth=1)
    
    # contourf of the module
    ##############################
    # from interpolation of the vectors in input_file, probably not the best
    
    # define the grid
    grid_x, grid_y = np.mgrid[0:360:200j, 90:-90:400j] # original. Bit weird
    #grid_x, grid_y = np.mgrid[0:360:200j, 0:180:400j]
    
    xg, yg = map1(grid_x, grid_y)
    
    
    XY = np.transpose([X,Y])

    
    
    # define the grid
    lats = np.linspace(np.pi/2,-np.pi/2,num=100)
    lons = np.linspace(-np.pi,np.pi,num=100)
    lats, lons = np.meshgrid(lats,lons)
    lats = np.transpose(lats)
    lons = np.transpose(lons)
    theta = -lats+np.pi/2
    
    Ush = np.zeros(theta.shape) # u_th
    Vsh = np.zeros(theta.shape) # u_ph

        
    xs, ys = map1(lons*180./np.pi, lats*180./np.pi)
    cs = map1.contourf(xs,ys,speed, cmap='Blues')
    plt.colorbar(cs,fraction=0.03, pad=0.04)
    
    
    #####################
    # quiver of the flow
    # make up some data on a regular lat/lon grid.
    uproj,vproj,xx,yy = \
    map1.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)
    
    
    for i in range(0,len(xx)):
        Q = map1.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=600,width=0.002,headwidth=6)
    plt.rcParams['font.size'] = 18
    plt.quiverkey(Q,0.5,-0.2, 20,'20 km/yr',coordinates='axes')
    
    plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)

    plt.close('all')
    
    return;
"""

def write_optimal_flow_input(filename,colat,lon,LMAX_U,LMAX_B_OBS,MODEL,TARGET_RMS,SCALE_FACTOR,RESTRICTION,ETA,FLAG):
    
    f1 = open(filename,"w+")
    f1.write('%(colat)f %(lon)f \n' %{'colat':colat, 'lon':lon} )
    f1.write('%(lmax_u)d\n' %{'lmax_u':LMAX_U} )
    f1.write('%(lmax_b_obs)d\n' %{'lmax_b_obs':LMAX_B_OBS} )
    f1.write('"%(model)s"\n' %{'model':MODEL} )
    f1.write('%(target_rms)d\n' %{'target_rms':TARGET_RMS} )
    f1.write('%(scale_factor)d\n' %{'scale_factor':SCALE_FACTOR} )
    f1.write('%(restriction)d\n' %{'restriction':RESTRICTION} )
    f1.write('%(eta)d\n' %{'eta':ETA} )
    f1.write('%(flag)d\n' %{'flag':FLAG} )
    f1.close()
    
    return;


def write_calc_SV_input(filename,LMAX_U,LMAX_B_OBS,FLOW,MODEL,REVERSE):
    
    f1 = open(filename,"w+")
    f1.write('%(lmax_u)d\n' %{'lmax_u':LMAX_U} )
    f1.write('%(lmax_b_obs)d\n' %{'lmax_b_obs':LMAX_B_OBS} )
    f1.write('"%(flow)s"\n' %{'flow':FLOW} )
    f1.write('"%(model)s"\n' %{'model':MODEL} )
    f1.write('%(reverse)d\n' %{'reverse':REVERSE} )
    f1.close()
    
    return;

def write_timesteping_opt_input(filename,LMAX_U,LMAX_B_OBS,FLOW,MODEL,FLAG_U_INIT,ETA,ALPHA,T0,TF,DT,REVERSE,colat,lon,TARGET_RMS,SCALE_FACTOR,RESTRICTION,FLAG,OUT_RATE):
    """
    to write the input file for timestepping_induction_opt.F90
    """
    f1 = open(filename,"w+")
    f1.write('%(lmax_u)d\n' %{'lmax_u':LMAX_U} )
    f1.write('%(lmax_b_obs)d\n' %{'lmax_b_obs':LMAX_B_OBS} )
    f1.write('"%(flow)s"\n' %{'flow':FLOW} )
    f1.write('"%(model)s"\n' %{'model':MODEL} )
    f1.write('%(flag_u_init)d\n' %{'flag_u_init':FLAG_U_INIT} )    
    f1.write('%(eta)f\n' %{'eta':ETA} )
    f1.write('%(alpha)f\n' %{'alpha':ALPHA} )
    f1.write('%(t0)f\n' %{'t0':T0} )
    f1.write('%(tf)f\n' %{'tf':TF} )
    f1.write('%(dt)f\n' %{'dt':DT} )
    f1.write('%(reverse)d\n' %{'reverse':REVERSE} )
    f1.write('%(colat)f %(lon)f \n' %{'colat':colat, 'lon':lon} )
    f1.write('%(target_rms)f\n' %{'target_rms':TARGET_RMS} )
    f1.write('%(scale_factor)f\n' %{'scale_factor':SCALE_FACTOR} )
    f1.write('%(restriction)d\n' %{'restriction':RESTRICTION} )
    f1.write('%(flag)d\n' %{'flag':FLAG} )
    f1.write('%(out_rate)d\n' %{'out_rate':OUT_RATE} )
    f1.close()
    
    return;

def write_timesteping_input(filename,LMAX_U,LMAX_B_OBS,FLOW,MODEL,ETA,ALPHA,T0,TF,DT,REVERSE):
    """
    to write the input file for timestepping_induction.F90
    """
    f1 = open(filename,"w+")
    f1.write('%(lmax_u)d\n' %{'lmax_u':LMAX_U} )
    f1.write('%(lmax_b_obs)d\n' %{'lmax_b_obs':LMAX_B_OBS} )
    f1.write('"%(flow)s"\n' %{'flow':FLOW} )
    f1.write('"%(model)s"\n' %{'model':MODEL} )
    f1.write('%(eta)f\n' %{'eta':ETA} )
    f1.write('%(alpha)f\n' %{'alpha':ALPHA} )
    f1.write('%(t0)f\n' %{'t0':T0} )
    f1.write('%(tf)f\n' %{'tf':TF} )
    f1.write('%(dt)f\n' %{'dt':DT} )
    f1.write('%(reverse)d\n' %{'reverse':REVERSE} )
    f1.close()
    
    return;

def read_coeffs(filename,hlines,LMAX_U):
    """
    reads in coefficients as stored in the CHAOS model files and ordered as:
        HEADER (hlines number of lines of discarded header)
        1 0 g10 h10
        ...
        l m glm hlm
        ...
    
    """
    input_file = open(filename, 'r')
    for i in range(hlines): 
        input_file.readline() # header
    coeffs = []
    while True:
        line = input_file.readline()
        x = [list(map(float, line.split() ))]
        l = x[0][0]
        m = x[0][1]
        coeffs = np.append(coeffs,x[0][:])
        if m == LMAX_U: break # end of file
    coeffs = np.reshape(coeffs,(len(coeffs)/len(x[0][:]), len(x[0][:])))     
    input_file.close()
    return coeffs;


def read_MF_file(filename_MF):
    coeffs_MF =[]

    with open(filename_MF, 'rb') as fipgp:
        while True:
            line = fipgp.readline()
            if not line:
                break
            x = [list(map(float, line.split() ))]
            coeffs_MF = np.append(coeffs_MF,x)
    cols = len(x[0])
    rows = len(coeffs_MF)
    coeffs_MF = np.reshape(coeffs_MF,(rows/cols, cols))
    
    return coeffs_MF;

def show_flow(folder):
    """
    original script that shows local flows on an orthographic projection centered around the 0,0 location
    """
    input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
    output_file = folder+'flow.pdf'
    
    X,Y,U,V = np.loadtxt(input_file,unpack=True)
    
    # display max flow on input grid
    USQ = []
    for i in range(len(U)):
        USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
    print ('Largest value of |u| on grid is ', max(USQ))
    long = X
    lat = Y
    
    plt.figure( figsize=(10,8))
    
    #Mollweide projection
    map = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    
    # Orthographic projection
    map = Basemap(projection='ortho',lat_0=0,lon_0=0,resolution='l')
    
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.50)
    
    
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30),linewidth=1)
    map.drawparallels(np.arange(-90,90,30),linewidth=1)
    # make up some data on a regular lat/lon grid.
    
    uproj,vproj,xx,yy = \
    map.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)
    
    
    for i in range(0,len(xx)):
        Q = map.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=400,width=0.002,headwidth=6)
    plt.rcParams['font.size'] = 18
    plt.quiverkey(Q,0.5,-0.1, 20,'20 km/yr',coordinates='axes')
    
    #plt.show()
    plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)
    plt.close('all')
    return;

def show_flow_global(folder):
    '''
    show quiver plot superimposed to intensity plot of the velocity fields from the files of the kind
    FLOWS_VECTORS_RANDOM.DAT
    OPTIMAL_FLOW.DAT
    '''
    input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
    output_file = folder+'flow_mollweide.pdf'
    
    X,Y,U,V = np.loadtxt(input_file,unpack=True)
    
    # display max flow on input grid
    USQ = []
    for i in range(len(U)):
        USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
    print ('Largest value of |u| on grid is ', max(USQ))
    long = X
    lat = Y
    
    plt.figure(1, figsize=(10,8))
    
    
    #Mollweide projectionfrom scipy.interpolate import griddata
    map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    
    # draw coastlines, country boundaries, fill continents.
    map1.drawcoastlines(linewidth=0.50)
    
    # draw lat/lon grid lines every 30 degrees.
    map1.drawmeridians(np.arange(0,360,30),linewidth=1)
    map1.drawparallels(np.arange(-90,90,30),linewidth=1)
    
    # contourf of the module
    ##############################
    # from interpolation of the vectors in input_file, probably not the best
    
    # define the grid
    grid_x, grid_y = np.mgrid[0:360:200j, 90:-90:400j] # original. Bit weird
    #grid_x, grid_y = np.mgrid[0:360:200j, 0:180:400j]
    
    xg, yg = map1(grid_x, grid_y)
    
    
    XY = np.transpose([X,Y])
    U_grid = griddata(XY, U, (grid_x, grid_y), method='cubic')
    V_grid = griddata(XY, V, (grid_x, grid_y), method='cubic')
    
    #cf=map1.contourf(xg,yg,np.sqrt(U_grid**2 + V_grid**2), 20, cmap = 'Blues')
    #plt.colorbar(cf,fraction=0.03, pad=0.04)
    
    #################################
    # from the coefficients, probably best
    Ucoeffs = open(folder+'OPTIMAL_FLOW.DAT', 'r')
    header = Ucoeffs.readline()
    lmaxU, lmaxB = [int(s) for s in header.split() if s.isdigit()]
    line = []
    coeffs_U = []
    while True:
        line = Ucoeffs.readline()
        x = [list(map(float, line.split() ))]
        l = x[0][0]
        m = x[0][1]
        coeffs_U = np.append(coeffs_U,x)
        if m == lmaxU: break # end of file
    
    coeffs_U = np.reshape(coeffs_U,(len(coeffs_U)/6, 6)) # 6 elements/row : l,m,s_{lm}^{cos},s_{lm}^{sin},t_{lm}^{cos},t_{lm}^{sin}
    
    Ucoeffs.close()
    
    # define the grid
    lats = np.linspace(np.pi/2,-np.pi/2,num=100)
    lons = np.linspace(-np.pi,np.pi,num=100)
    lats, lons = np.meshgrid(lats,lons)
    lats = np.transpose(lats)
    lons = np.transpose(lons)
    theta = -lats+np.pi/2
    
    Ush = np.zeros(theta.shape) # u_th
    Vsh = np.zeros(theta.shape) # u_ph
    r_cmb = 3485.0e3
    
    for ic in range(len(coeffs_U)):
        l, m, slm_c, slm_s, tlm_c, tlm_s = coeffs_U[ic][0:6]
        dthSH_c = DthSchmidtSH(l,m,theta,lons,'c')
        dthSH_s = DthSchmidtSH(l,m,theta,lons,'s')
        dphSH_c = DphSchmidtSH(l,m,theta,lons,'c')
        dphSH_s = DphSchmidtSH(l,m,theta,lons,'s')
        
        Ush = Ush + tlm_c * np.divide(dphSH_c,np.sin(theta))  \
                  + tlm_s * np.divide(dphSH_s,np.sin(theta))  \
                  + slm_c * dthSH_c  \
                  + slm_s * dthSH_s  
    
        Vsh = Vsh - tlm_c * dthSH_c / r_cmb \
                  - tlm_s * dthSH_s / r_cmb \
                  + slm_c * np.divide(dphSH_c,np.sin(theta))  \
                  + slm_s * np.divide(dphSH_s,np.sin(theta)) 
        
    # adjust the units
    Ush = Ush * r_cmb * 1e-3
    Vsh = Vsh * r_cmb * 1e-3 
    speed = np.sqrt(Ush**2 + Vsh**2)
        
    xs, ys = map1(lons*180./np.pi, lats*180./np.pi)
    cs = map1.contourf(xs,ys,speed, cmap='Blues')
    plt.colorbar(cs,fraction=0.03, pad=0.04)
    
    
    #####################
    # quiver of the flow
    # make up some data on a regular lat/lon grid.
    uproj,vproj,xx,yy = \
    map1.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)
    
    
    for i in range(0,len(xx)):
        Q = map1.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=600,width=0.002,headwidth=6)
    plt.rcParams['font.size'] = 18
    plt.quiverkey(Q,0.5,-0.2, 20,'20 km/yr',coordinates='axes')
    
    plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)

    plt.close('all')
    return;


def show_flow_global(folder):
    '''
    show quiver plot superimposed to intensity plot of the velocity fields from the files of the kind
    FLOWS_VECTORS_RANDOM.DAT
    OPTIMAL_FLOW.DAT
    '''
    input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
    output_file = folder+'flow_mollweide.pdf'
    
    X,Y,U,V = np.loadtxt(input_file,unpack=True)
    
    # display max flow on input grid
    USQ = []
    for i in range(len(U)):
        USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
    print ('Largest value of |u| on grid is ', max(USQ))
    long = X
    lat = Y
    
    plt.figure(1, figsize=(10,8))
    
    
    #Mollweide projectionfrom scipy.interpolate import griddata
    map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    
    # draw coastlines, country boundaries, fill continents.
    map1.drawcoastlines(linewidth=0.50)
    
    # draw lat/lon grid lines every 30 degrees.
    map1.drawmeridians(np.arange(0,360,30),linewidth=1)
    map1.drawparallels(np.arange(-90,90,30),linewidth=1)
    
    # contourf of the module
    ##############################
    # from interpolation of the vectors in input_file, probably not the best
    
    # define the grid
    grid_x, grid_y = np.mgrid[0:360:200j, 90:-90:400j] # original. Bit weird
    #grid_x, grid_y = np.mgrid[0:360:200j, 0:180:400j]
    
    xg, yg = map1(grid_x, grid_y)
    
    
    XY = np.transpose([X,Y])
    U_grid = griddata(XY, U, (grid_x, grid_y), method='cubic')
    V_grid = griddata(XY, V, (grid_x, grid_y), method='cubic')
    
    #cf=map1.contourf(xg,yg,np.sqrt(U_grid**2 + V_grid**2), 20, cmap = 'Blues')
    #plt.colorbar(cf,fraction=0.03, pad=0.04)
    
    #################################
    # from the coefficients, probably best
    Ucoeffs = open(folder+'OPTIMAL_FLOW.DAT', 'r')
    header = Ucoeffs.readline()
    lmaxU, lmaxB = [int(s) for s in header.split() if s.isdigit()]
    line = []
    coeffs_U = []
    while True:
        line = Ucoeffs.readline()
        x = [list(map(float, line.split() ))]
        l = x[0][0]
        m = x[0][1]
        coeffs_U = np.append(coeffs_U,x)
        if m == lmaxU: break # end of file
    
    coeffs_U = np.reshape(coeffs_U,(len(coeffs_U)/6, 6)) # 6 elements/row : l,m,s_{lm}^{cos},s_{lm}^{sin},t_{lm}^{cos},t_{lm}^{sin}
    
    Ucoeffs.close()
    
    # define the grid
    lats = np.linspace(np.pi/2,-np.pi/2,num=100)
    lons = np.linspace(-np.pi,np.pi,num=100)
    lats, lons = np.meshgrid(lats,lons)
    lats = np.transpose(lats)
    lons = np.transpose(lons)
    theta = -lats+np.pi/2
    
    Ush = np.zeros(theta.shape) # u_th
    Vsh = np.zeros(theta.shape) # u_ph
    r_cmb = 3485.0e3
    
    for ic in range(len(coeffs_U)):
        l, m, slm_c, slm_s, tlm_c, tlm_s = coeffs_U[ic][0:6]
        dthSH_c = DthSchmidtSH(l,m,theta,lons,'c')
        dthSH_s = DthSchmidtSH(l,m,theta,lons,'s')
        dphSH_c = DphSchmidtSH(l,m,theta,lons,'c')
        dphSH_s = DphSchmidtSH(l,m,theta,lons,'s')
        
        Ush = Ush + tlm_c * np.divide(dphSH_c,np.sin(theta))  \
                  + tlm_s * np.divide(dphSH_s,np.sin(theta))  \
                  + slm_c * dthSH_c  \
                  + slm_s * dthSH_s  
    
        Vsh = Vsh - tlm_c * dthSH_c / r_cmb \
                  - tlm_s * dthSH_s / r_cmb \
                  + slm_c * np.divide(dphSH_c,np.sin(theta))  \
                  + slm_s * np.divide(dphSH_s,np.sin(theta)) 
        
    # adjust the units
    Ush = Ush * r_cmb * 1e-3
    Vsh = Vsh * r_cmb * 1e-3 
    speed = np.sqrt(Ush**2 + Vsh**2)
        
    xs, ys = map1(lons*180./np.pi, lats*180./np.pi)
    cs = map1.contourf(xs,ys,speed, cmap='Blues')
    plt.colorbar(cs,fraction=0.03, pad=0.04)
    
    
    #####################
    # quiver of the flow
    # make up some data on a regular lat/lon grid.
    uproj,vproj,xx,yy = \
    map1.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)
    
    
    for i in range(0,len(xx)):
        Q = map1.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=600,width=0.002,headwidth=6)
    plt.rcParams['font.size'] = 18
    plt.quiverkey(Q,0.5,-0.2, 20,'20 km/yr',coordinates='axes')
    
    plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)

    plt.close()
    return;
    
def show_flow_Br_global(folder,model,r_cmb,r_a,LMAX_B_OBS,hline_B,time,units,colatS,lonS):
    '''
    show quiver plot superimposed to intensity plot of the velocity fields from the files of the kind
    FLOWS_VECTORS_RANDOM.DAT
    OPTIMAL_FLOW.DAT
    Superimpose do a snapshot of a magnetic field model provided in the file 
    model
    '''
    input_file = folder+'FLOW_VECTORS_RANDOM.DAT'
    output_file = folder+'flow_Br_mollweide.pdf'
    
    
    
    # load magnetic field
    coeffsBmat = read_coeffs(model,hline_B,LMAX_B_OBS)
    
    
    lats = np.linspace(np.pi/2,-np.pi/2,num=100)
    lons = np.linspace(-np.pi,np.pi,num=100)
    lats, lons = np.meshgrid(lats,lons)
    lats = np.transpose(lats)
    lons = np.transpose(lons)
    theta = -lats+np.pi/2
    
    Br, Bt, Bp = calcB(coeffsBmat,theta,lons,r_a,r_cmb)    
    
    
    X,Y,U,V = np.loadtxt(input_file,unpack=True)
    
    # display max flow on input grid
    USQ = []
    for i in range(len(U)):
        USQ.append( ( U[i]**2 + V[i]**2 )**0.5 )
    print ('Largest value of |u| on grid is ', max(USQ))
    long = X
    lat = Y
    
    plt.figure(1, figsize=(10,8))
    
    
    #Mollweide projectionfrom scipy.interpolate import griddata
    map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')
    
    # draw coastlines, country boundaries, fill continents.
    map1.drawcoastlines(linewidth=0.50)
    
    # draw lat/lon grid lines every 30 degrees.
    map1.drawmeridians(np.arange(0,360,30),linewidth=1)
    map1.drawparallels(np.arange(-90,90,30),linewidth=1)
    
    # contourf of the module
    ##############################
    # from interpolation of the vectors in input_file, probably not the best
    
    # define the grid
    #grid_x, grid_y = np.mgrid[0:360:200j, 90:-90:400j] # original. Bit weird
    #grid_x, grid_y = np.mgrid[0:360:200j, 0:180:400j]
    
    #xg, yg = map1(grid_x, grid_y)

    xs, ys = map1(lons*180./np.pi, lats*180./np.pi)
    cs = map1.contourf(xs,ys,Br/1000., cmap='coolwarm')
    plt.colorbar(cs,fraction=0.03, pad=0.04)
    plt.text(42000000, 8500000, units , fontsize=16)

    
    #####################
    # quiver of the flow
    # make up some data on a regular lat/lon grid.
    uproj,vproj,xx,yy = \
    map1.rotate_vector(np.array(V),-np.array(U),np.array(long),np.array(lat),returnxy=True)
    
    
    for i in range(0,len(xx)):
        Q = map1.quiver(xx[i],yy[i],uproj[i],vproj[i],angles='uv',pivot='tail',scale=600,width=0.002,headwidth=6)
    plt.rcParams['font.size'] = 18
    plt.quiverkey(Q,0.5,-0.2, 20,'20 km/yr',coordinates='axes')
    
    #station location
    xS, yS = map1(lonS, colatS)
    if colatS != None:
        map1.scatter(xS,yS,
                    s=80,edgecolor='k',color='w',alpha=1, zorder=3,marker='^')
        
    plt.savefig(output_file, bbox_inches='tight',pad_inches=0.7,dpi=400)

    plt.close()
    return;    
    
    
    
def show_flow_streamlines(folder,it,r_cmb,r_a,units,colatS,lonS):
    '''
    plot flows and CMB field from the data stored in the output of the timestepping/optimisation routines:
        MF_COEFFS.DAT
        U_POL_COEFFS.DAT
        U_TOR_COEFFS.DAT
    INPUT:
        folder: location of the .DAT files
        it : time instant (line of the files) that is to be plotted
    '''
    
    ###################################
    # define the grid
    # input lat-lon resolution
    nlats = 80
    nlons = 150
    lats = np.linspace(-np.pi / 2, np.pi / 2, nlats)
    lons = np.linspace(0, 2 * np.pi, nlons)
    lons, lats = np.meshgrid(lons, lats)
    theta = -lats+np.pi/2
    
    Ush = np.zeros(theta.shape) # u_th
    Vsh = np.zeros(theta.shape) # u_ph
    Br_c = np.zeros(theta.shape) # B_r
    Br_a = np.zeros(theta.shape) # B_r
    Bt_a = np.zeros(theta.shape) # B_r
    Bp_a = np.zeros(theta.shape) # B_r

    # input legends
    legend1 = ''

    #################################
    # read U coefficients
    f_U_pol_coeffs = open(folder+'U_POL_COEFFS.DAT', 'r')
    f_U_tor_coeffs = open(folder+'U_TOR_COEFFS.DAT', 'r')
    lines_pol=f_U_pol_coeffs.readlines()
    lines_tor=f_U_tor_coeffs.readlines()
    f_U_pol_coeffs.close()
    f_U_tor_coeffs.close()
    
    line_pol = lines_pol[it]
    line_tor = lines_tor[it]
    x = [list(map(float, line_pol.split() ))]
    U_pol_coeffs = np.asarray(x[0])
    x = [list(map(float, line_tor.split() ))]
    U_tor_coeffs = np.asarray(x[0])
    LmaxU = int(math.ceil( -1+np.sqrt(U_tor_coeffs.size) ))
    
    ###################################
    # read MF coefficients
    f_MF_coeffs = open(folder+'MF_COEFFS.DAT', 'r')
    lines_MF=f_MF_coeffs.readlines()
    line_MF = lines_MF[it]
    f_MF_coeffs.close()
    
    x = [list(map(float, line_MF.split() ))]
    MF_coeffs = np.asarray(x[0])
    LmaxB = int(math.ceil( -1+np.sqrt(MF_coeffs.size) ))


    # calculate B and U
    l = 1
    m = 0
    cs = 'c'
    KE_pol = np.zeros(LmaxU)
    KE_tor = np.zeros(LmaxU)
    ME_a   = np.zeros(LmaxB)
    ME_cmb = np.zeros(LmaxB)
    
    for ic in range(max([len(U_pol_coeffs[1:-1]), len(MF_coeffs[1:-1])])):    
        dthSH = DthSchmidtSH(l,m,theta,lons,cs)
        dphSH = DphSchmidtSH(l,m,theta,lons,cs)
    
        if ic < len(U_pol_coeffs[1:-1]):
            
            slm = U_pol_coeffs[ic+1]
            tlm = U_tor_coeffs[ic+1]
            
            Ush = Ush + tlm * np.divide(dphSH,np.sin(theta))  \
                      + slm * dthSH  \
        
            Vsh = Vsh - tlm * dthSH / r_cmb \
                      + slm * np.divide(dphSH,np.sin(theta))  \
            # spectra of flow
            # correct for the factor that was multiplied to the coefficients in the input file
            # the coefficients are still normalized to the target u_rms value          
            factor = (r_cmb/1000)
            q_pol = 	slm*factor
            q_tor = 	tlm*factor
    
            # the coefficients should be such that u_rms is the target rms (13 km/yr in Livermore, 2014)  
            norm = 4*np.pi*l*(l+1)/(2*l+1)
            KE_pol[l-1] =  KE_pol[l-1] + norm * (q_pol**2)
            KE_tor[l-1] =  KE_tor[l-1] + norm * (q_tor**2)
           
        if ic < len(MF_coeffs[1:-1]):
            
            glm = MF_coeffs[ic+1]
            SH    = SchmidtSH(l,m,theta,lons,cs)           
            
            Br_c = Br_c + (l+1) * (r_a/r_cmb)**(l+2) *glm * SH                     
            Br_a = Br_a + (l+1) * (r_a/r_a)**(l+2) * glm * SH                    
            Bt_a = Bt_a - (r_a/r_a)**(l+2) * glm * dthSH               
            Bp_a = Bp_a - (r_a/r_a)**(l+2) * glm * np.divide(dphSH,np.sin(theta))   
            
            ME_a[l-1]   =  ME_a[l-1]    + (l+1) * (glm**2)
            ME_cmb[l-1] =  ME_cmb[l-1]  + (l+1) * (r_a/r_cmb)**(2*l+4) * (glm**2)
    
        # updte indices
        if cs == 'c':
            if m==0:
                m=1
            elif m>0:
                cs = 's'
        elif cs == 's':
            if m<l:
                cs='c'
                m = m+1
            elif m == l:
                l=l+1
                cs = 'c'
                m=0
                
                
    F_a = np.sqrt(Br_a**2 + Bt_a**2 + Bp_a**2) 
    Incl = np.rad2deg(np.arctan(-Br_a/np.sqrt(Bt_a**2 + Bp_a**2)))     
    Ush = Ush * r_cmb * 1e-3
    Vsh = Vsh * r_cmb * 1e-3 
    speed = np.sqrt(Ush**2 + Vsh**2)
    
    lats = np.rad2deg(lats)
    lons = np.rad2deg(lons)
    
    #------#
    # plots
    #------#

    # Br-flow map at CMB
    #sinlge plot
    files = []
    
    figBr = plt.figure(figsize=(11, 6))
    #    plt.cla()
    axBr = figBr.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    #axBr.gridlines(linewidth=1, alpha=0.5, linestyle='-')
    cf = axBr.contourf(lons, lats, Br_c/1000, # breaks down with contourlevels >=24
                transform=ccrs.PlateCarree(),
                cmap='coolwarm')
    axBr.contour(lons, lats, Br_c/1000,23, 
                transform=ccrs.PlateCarree(),
                colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
    
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    
    plt.text(24300000, 0, units , fontsize=16)
    lw = 1.9*speed / np.nanmax(speed)
    axBr.streamplot(lons, lats, Vsh, -Ush, 
                transform=ccrs.PlateCarree(),
                linewidth=lw,
                arrowsize=0.5,
                density=6,
                color='k')
    axBr.coastlines()
    axBr.set_global()
    #station location
    if colatS != None:
        axBr.scatter(lonS,colatS,
                   s=80,edgecolor='k',color='w',alpha=1, zorder=3,marker='^',
                   transform=ccrs.PlateCarree())
    
    plt.title(U_pol_coeffs[0], y=1.08)
    fname = folder +'Br_CMB_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()
    
    
    # Surface intensity
    

    files = []
    
    figF = plt.figure(figsize=(11, 6))
    #    plt.cla()
    axF = figF.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    #axF.gridlines(linewidth=1, alpha=0.5, linestyle='-')
    cf = axF.contourf(lons, lats, F_a, # breaks down with contourlevels >=24
                transform=ccrs.PlateCarree(),
                cmap='jet')
    cs=axF.contour(lons, lats, F_a,23, 
                transform=ccrs.PlateCarree(),
                colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    #cf.clim(vmin=-Brmax, vmax=Brmax)
    plt.text(24300000, 0, units , fontsize=16)
    axF.coastlines()
    axF.set_global()
    
    
    plt.title(U_pol_coeffs[0], y=1.08)
    fname = folder +'F_surface_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()

    # flow map at CMB
    #sinlge plot
    files = []
    
    figU = plt.figure(figsize=(11, 6))
    #    plt.cla()
    axU = figU.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    #axU.gridlines(linewidth=1, alpha=0.5, linestyle='-')
    
    cf = axU.contourf(lons, lats, speed, # breaks down with contourlevels >=24
                transform=ccrs.PlateCarree(),
                cmap='Blues')
    axU.contour(lons, lats, speed,23, 
                transform=ccrs.PlateCarree(),
                colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
    
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    
#    plt.text(24300000, 0, r'$\mu T$' , fontsize=16)
    lw = 1.9*speed / np.nanmax(speed)
    axU.streamplot(lons, lats, Vsh, -Ush, 
                transform=ccrs.PlateCarree(),
                linewidth=lw,
                arrowsize=0.5,
                density=6,
                color='k')
    axU.coastlines()
    axU.set_global()
    
    plt.title(U_pol_coeffs[0], y=1.08)
    fname = folder +'U_CMB_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()

    # flow inclination at surface
    #sinlge plot
    files = []
    
    figI = plt.figure(figsize=(11, 6))
    #    plt.cla()
    axI = figI.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    #axU.gridlines(linewidth=1, alpha=0.5, linestyle='-')
    
    cf = axI.contourf(lons, lats, Incl, # breaks down with contourlevels >=24
                transform=ccrs.PlateCarree(),
                cmap='viridis',vmin=-90.0, vmax=90.0)
    axI.contour(lons, lats, Incl,23, 
                transform=ccrs.PlateCarree(),
                colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
    
    #v = np.linspace(-90.0, 2.0, 90.0)
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    plt.text(23000000, 0, 'deg' , fontsize=16)

    axI.coastlines()
    axI.set_global()
    #station location
    if colatS != None:
        axI.scatter(lonS,colatS,
                   s=80,edgecolor='k',color='w',alpha=1, zorder=3,marker='^',
                   transform=ccrs.PlateCarree())
    
    plt.title(U_pol_coeffs[0], y=1.08)
    fname = folder +'Incl_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()

    return

def show_flow_streamlines_inst(folder,flow,model,r_cmb,r_a,LMAX_U,hline_U,LMAX_B_OBS,hline_B,time,units,colatS,lonS):
    '''
    plot flows and CMB field from the data stored in the output of the 
    instantaneous calculation routines:
        OPTIMAL_FLOW.DAT
    and from the field model

    INPUT:
        folder: location of the optimal flow model
        flow:   optimal flow file (OPTIMAL_FLOW.DAT)
        model:  file containing the background magnetic field
    '''
    
    # read magnetix field
    coeffsBmat = read_coeffs(model,hline_B,LMAX_B_OBS)
    
    coeffsBlin = mat2linCoeffs(coeffsBmat,time)
    
    # read flow coefficients
    coeffsUmat = read_coeffs(flow,hline_U,LMAX_U)
    coeffsPOLmat = coeffsUmat[:,0:4]
    coeffsTORmat = np.zeros(coeffsPOLmat.shape)
    coeffsTORmat[:,0:2] = coeffsUmat[:,0:2]
    coeffsTORmat[:,2:4] = coeffsUmat[:,4:6]
    
    coeffsPOLlin = mat2linCoeffs(coeffsPOLmat,time)
    coeffsTORlin = mat2linCoeffs(coeffsTORmat,time)
    
    # save to files that show_flow_streamlines accepts
    np.savetxt(folder+'/MF_COEFFS.DAT', coeffsBlin[None], fmt="%s")
    np.savetxt(folder+'/U_POL_COEFFS.DAT', coeffsPOLlin[None], fmt="%s")
    np.savetxt(folder+'/U_TOR_COEFFS.DAT', coeffsTORlin[None], fmt="%s")
    
    show_flow_streamlines(folder+"/",0,r_cmb,r_a,units,colatS,lonS)
    
    
    return 


def show_flow_spectra(folder):
    # read in coefficients
    line_num = 1
    coeff = []
    with open(folder+'OPTIMAL_FLOW.DAT') as f:
      for line in f:
        if line_num == 1:
          line_num = line_num + 1
        else:
          coeff.append([float(x) for x in line.split()])
          line_num = line_num + 1
          
    lmax = int(coeff[len(coeff)-1][0])
    
    KE_pol = np.zeros(lmax)
    KE_tor = np.zeros(lmax)
    
    for i in range(len(coeff)):
      l = 		int(coeff[i][0])
      m = 		int(coeff[i][1])
      # correct for the factor that was multiplied to the coefficients in the input file
      # the coefficients are still normalized to the target u_rms value
      factor = (c / 1000) 
      q_pol_cos = 	coeff[i][2]*factor
      q_pol_sin = 	coeff[i][3]*factor
      q_tor_cos = 	coeff[i][4]*factor
      q_tor_sin = 	coeff[i][5]*factor
    
      # the coefficients should be such that u_rms is the target rms (13 km/yr in Livermore, 2014)  
      norm = 4*np.pi*l*(l+1)/(2*l+1)
      KE_pol[l-1] =  KE_pol[l-1] + norm * (q_pol_cos**2 + q_pol_sin**2)
      KE_tor[l-1] =  KE_tor[l-1] + norm * (q_tor_cos**2 + q_tor_sin**2)
    
    print 'poloidal kinetic energy: ', sum(KE_pol)
    print 'rms poloidal velocity: ', np.sqrt( sum(KE_pol) / (4*np.pi) ) 
    print 'toroidal kinetic energy: ', sum(KE_tor)
    print 'rms toroidal velocity: ', np.sqrt( sum(KE_tor) / (4*np.pi) ) 
    
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    ax = plt.subplot(111)
    x = range(1,lmax+1)
    ax.plot(x,KE_pol,label='poloidal',color='red')
    ax.plot(x,KE_tor,color='blue',label='toroidal')
    
    # adjust axes
    ax.ticklabel_format(style='sci')
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.xlim(1,lmax)
    #ax.set_yscale('log')
    ax.legend(fontsize=15)
    ax.set_xlabel(r'$l$', fontsize=15)
    ax.set_ylabel(r'$K(l)$',rotation=0,labelpad=20, fontsize=15)
    ax.tick_params(labelsize=15)
    #plt.xticks(np.arange(1, max(x)+1, 1.0))
    
    plt.grid(True)
    #plt.yscale('log')
    plt.savefig(folder+'Flow_spectra.pdf',bbox_inches='tight',pad_inches=0.0)
    plt.close('all')
    return;

def show_global_contourf(LONS,LATS,F,COLORMAP,figname):
    """
    global, mollweide map of the field F
    """
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
    cs = map1.contourf(xs,ys,F,cmap=COLORMAP)
    plt.colorbar(cs,fraction=0.03, pad=0.04)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.0)
    plt.close('all')
    return;

def show_global_scatterplot(LONS,LATS,F,COLORMAP,figname):
    
    fig = plt.figure(1, figsize=(11, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
    
    ax.coastlines()
    ax.set_global()
    
    # geomagnetic poles location
    sc=plt.scatter(LONS,LATS,
               s=20, c=F,edgecolor='k',alpha=1, zorder=3,marker='o',
               transform=ccrs.PlateCarree(), cmap=COLORMAP)
    plt.colorbar(sc)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.0)
    plt.close('all')
    return;


def equispaced_grid(Nth):
# Create equispoaced grid on a sphere
# based on a document by Markus Deserno (How to generate equidistributed points on the surface of a sphere)
    r_a = 6371e3 # radius of the earth

    dth = np.pi/(Nth-1) # distance in lat
    x=[]
    y=[]
    z=[]
    lats=[]
    lons=[]
    for ith in range(Nth):
        rs = r_a * np.sin(ith*dth) # radius of the circle    
        zs = r_a * np.cos(ith*dth) # and its height above equator
        if (ith==0) or (ith==Nth-1):
            nphi = 0
            dphi = 0
            x.append(0) 
            y.append(0)
            z.append(zs)
            lats.append(ith*dth)
            lons.append(0)
        else:
            nphi = math.ceil(2*np.pi*rs/dth/r_a)
            dphi=2*np.pi/nphi; # distance in lon    
            for iphi in range(int(nphi)):
                x.append(rs*np.cos(iphi*dphi))
                y.append(rs*np.sin(iphi*dphi))
                z.append(zs)
                lats.append(ith*dth)
                lons.append(iphi*dphi)
                
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    
    
    theta = -lats+np.pi/2
    
    theta = np.rad2deg(theta)
    lons = np.rad2deg(lons)
    lats = np.rad2deg(lats)
    return theta, lons, lats


def VGP_from_DI(Incl,Decl,lat,lon):
# calculate VGP latitude and longitude from time-series of Inclination, Declination
# INPUTS:
# I: time-series of inclination at location
# D: time-series of declination at location
# lat, lon: coordinate of location in degrees
# OUTPUT:
# VGP_lat, VGP_lon: latitude and longitude of VGP, in degrees
#
# Formulas from Encyclopedia of Geomagnetism and Paleomagnetism

    lat = np.deg2rad(lat)#np.pi*lat/180
    lon = np.deg2rad(lon)#np.pi*lon/180
    Incl = np.deg2rad(Incl)#np.pi*Incl/180
    Decl = np.deg2rad(Decl)#np.pi*Decl/180
    
    VGP_lat= np.zeros(Incl.shape)
    VGP_lat[:] = np.nan 
    VGP_lon= np.zeros(Incl.shape)
    VGP_lon[:] = np.nan 
    
    for i in range(len(Incl)):
        # relative VGP colat
        p = np.arctan(2.0 / np.tan(Incl[i]))  
        if p <0:
            p = p+np.pi
        # VGP latitude
        VGP_lat[i] = np.arcsin( np.sin(lat)*np.cos(p) + np.cos(lat)*np.sin(p)*np.cos(Decl[i]) )
        # major arc distance
        beta = np.arcsin( np.sin(p)*np.sin(Decl[i])/np.cos(VGP_lat[i]) )
        # VGP longitude
        if np.cos(p)>=np.sin(lat)*np.sin(VGP_lat[i]):
            VGP_lon[i] = lon + beta
        else:
            VGP_lon[i] = lon + np.pi - beta
        
    VGP_lon = 180 * VGP_lon/np.pi
    VGP_lat = 180 * VGP_lat/np.pi
    
    
    return VGP_lat, VGP_lon
        

def dVGP_dt_from_DI(Incl,Decl,lat,lon,dIdt,dDdt):
    """
    calculate variations of VGP latitude and longitude (dVGP_lat_dt, dVGP_lon_dt)
    from temporal variations of Inclination and Declination  (dIdt,dDdt) and
    and Inclination and Declination themselves (Incl, Decl)
    measured at the location (lon,lat).
    
    All quantities exept lon and lat have to be arrays
    """
    
    lat = np.deg2rad(lat)#np.pi*lat/180
    lon = np.deg2rad(lon)#np.pi*lon/180
    Incl = np.deg2rad(Incl)#np.pi*Incl/180
    Decl = np.deg2rad(Decl)#np.pi*Decl/180
    dIdt = np.deg2rad(dIdt)
    dDdt = np.deg2rad(dDdt)
    
    VGP_lat= np.zeros(Incl.shape)
    VGP_lat[:] = np.nan 
    VGP_lon= np.zeros(Incl.shape)
    VGP_lon[:] = np.nan 
    dVGP_lat_dt= np.zeros(Incl.shape)
    dVGP_lat_dt[:] = np.nan 
    dVGP_lon_dt= np.zeros(Incl.shape)
    dVGP_lon_dt[:] = np.nan 
    
    p= np.zeros(Incl.shape)
    p[:] = np.nan     
    dpdt= np.zeros(Incl.shape)
    dpdt[:] = np.nan    
    
    beta= np.zeros(Incl.shape)
    beta[:] = np.nan     
    dbetadt= np.zeros(Incl.shape)
    dbetadt[:] = np.nan    
    
    for i in range(len(Incl)):
        # relative VGP colat
        p[i] = np.arctan(2.0 / np.tan(Incl[i]))  
        if p[i] <0:
            p[i] = p[i]+np.pi
        # VGP latitude
        VGP_lat[i] = np.arcsin( np.sin(lat)*np.cos(p[i]) + np.cos(lat)*np.sin(p[i])*np.cos(Decl[i]) )
        # major arc distance
        beta[i] = np.arcsin( np.sin(p[i])*np.sin(Decl[i])/np.cos(VGP_lat[i]) )
        # VGP longitude
        if np.cos(p[i])>=np.sin(lat)*np.sin(VGP_lat[i]):
            VGP_lon[i] = lon + beta[i]
        else:
            VGP_lon[i] = lon + np.pi - beta[i]
            
        # temporal variations
        dpdt[i] = - 2. * dIdt[i] / ( 1 + 3*(np.cos(Incl[i]))**2 )
        dVGP_lat_dt[i] = ( -np.sin(lat)*np.sin(p[i])*dpdt[i] 
                          + np.cos(lat)*np.cos(p[i])*dpdt[i]*np.cos(Decl[i]) 
                          - np.cos(lat)*np.sin(p[i])*np.sin(Decl[i])*dDdt[i] ) / np.abs(np.cos(VGP_lat[i]))
        
        dbetadt[i] = ( np.cos(VGP_lat[i])*np.cos(p[i])*dpdt[i]*np.sin(Decl[i])
                  + np.cos(VGP_lat[i])*np.sin(p[i])*np.cos(Decl[i])*dDdt[i]
                  + np.sin(p[i])*np.sin(Decl[i])*np.sin(VGP_lat[i])*dVGP_lat_dt[i]  ) / ( np.cos(VGP_lat[i])**2 * np.abs(np.cos(beta[i])) )
        if np.cos(p[i])>=np.sin(lat)*np.sin(VGP_lat[i]):
            dVGP_lon_dt[i] = dbetadt[i]
        else:
            dVGP_lon_dt[i] = -dbetadt[i]
            
            
    dVGP_lon_dt = 180 * dVGP_lon_dt/np.pi
    dVGP_lat_dt = 180 * dVGP_lat_dt/np.pi
    """
    # test figure: VGP lat and its time derivatives
    fig,ax = plt.subplots(figsize=(8,5))
    ax_twin = ax.twinx()
    
    ax.set_xlabel('Age / index')
    ax.set_ylabel('VGP latitude/$^\circ$')
    
    ax.plot(beta,color='k', label='beta')
    ax_twin.plot(dbetadt,color='r', label='dbetadt, exact')
    ax_twin.plot(0.01*(beta[1:]-beta[0:-1]),color='b', label='betapdt, first differences')
    ax.legend(fontsize=10,loc='upper left')
    ax_twin.legend(fontsize=10,loc='upper right')
    plt.title('beta and dbetadt')
    plt.show()
    """

    return dVGP_lat_dt, dVGP_lon_dt


def local_RoC_from_B(Br,Bt,Bp,Brdot,Btdot,Bpdot,lat,lon):
    """
    from the geomagnetic field components (Br,Bt,Bp) and time derivatives (Brdot,Btdot,Bpdot)
    at location (lon,lat), calculate Incl,Decl,VGP position and their rate of change.
    All quantities are scalar
    """
    
    H = np.sqrt(Bp**2 + Bt**2)
    F = np.sqrt(Br**2 + Bt**2 + Bp**2)
    
    Incl = np.arctan(np.divide(-Br,np.sqrt(Bt**2+Bp**2)))*180/np.pi
    Decl = math.atan2(-Bt,Bp)*180/np.pi
    VGP_lat, VGP_lon = VGP_from_DI([Incl],[Decl],lat,lon)
    VGP_lon = VGP_lon[0]
    VGP_lat = VGP_lat[0]
    
    dIdt = (Bt*Br*Btdot + Bp*Br*Bpdot - H**2*Brdot) / (H*F**2)
    dDdt = (Bp*Btdot - Bt*Bpdot) / H**2    
    dVGP_lat_dt, dVGP_lon_dt = dVGP_dt_from_DI([Incl],[Decl],lat,lon,[dIdt*180/np.pi],[dDdt*180/np.pi])
    dVGP_lat_dt = dVGP_lat_dt[0]
    dVGP_lon_dt = dVGP_lon_dt[0]

    return Incl, Decl, VGP_lat, VGP_lon, dIdt, dDdt, dVGP_lat_dt, dVGP_lon_dt