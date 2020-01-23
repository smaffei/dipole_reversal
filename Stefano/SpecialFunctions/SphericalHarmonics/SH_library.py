#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 13:38:28 2019

@author: Stefano Maffei
"""

import numpy as np
from scipy.special import lpmv # associate legendre functions, in vector form
import math
#from decimal import *


# define some common variables:
r_a = 6371000
r_c = 3485000    
    
def SchmidtSH(l,m,theta,phi,cs):
    '''
    Calculate Schmidt semi-normalized Spehrical harmonics in the form
    P_{lm}(cos(theta)) * cos(m*phi)    if cs == 'c'
    P_{lm}(cos(theta)) * sin(m*phi)    if cs == 's'
    
    theta and phi should be in meshgrid form with
    0<=theta<=pi is the colatitude
    -180<=phi<=180 is the longitude
    
    l, m = order and degree (scalars)
    '''
    if m == 0: FactorS = 1
    if m > 0: FactorS =np.multiply((-1)**m,np.sqrt(2*np.divide(math.factorial(l-m),float(math.factorial(l+m)))))
#    if m > 0: 
#        FactorS=np.multiply((-1)**m,np.sqrt(2*np.divide(Decimal(math.factorial(l-m)),Decimal(math.factorial(l+m)))))
#        FactorS=float(FactorS)
    if cs == 'c' : Trig = np.cos(m*phi)
    if cs == 's' : Trig = np.sin(m*phi)
    
    SH = FactorS*np.multiply(lpmv(m,l,np.cos(theta)),Trig)
    return SH;


def DthSchmidtSH(l,m,theta,phi,cs):
    '''
    Calculate theta (latitudinal) derivative of Schmidt semi-normalized Spehrical harmonics in the form
    d(P_{lm}(cos(theta)))dtheta * cos(m*phi)    if cs == 'c'
    d(P_{lm}(cos(theta)))dtheta  * sin(m*phi)    if cs == 's'
    
    theta and phi should be in meshgrid form with
    0<=theta<=pi is the colatitude
    -180<=phi<=180 is the longitude
    
    l, m = order and degree (scalars)
    '''
    if m == 0: FactorS = 1
    if m > 0: FactorS =np.multiply((-1)**m,np.sqrt(2*np.divide(math.factorial(l-m),float(math.factorial(l+m)))))
    if cs == 'c' : Trig = np.cos(m*phi)
    if cs == 's' : Trig = np.sin(m*phi)
    
    dPdth = - 0.5* ( (l+m)*(l-m+1) * lpmv(m-1,l,np.cos(theta)) - lpmv(m+1,l,np.cos(theta)) )
    SH = FactorS*np.multiply(dPdth,Trig)
    return SH;


def DphSchmidtSH(l,m,theta,phi,cs):
    '''
    Calculate phi (longitudinal) derivative of Schmidt semi-normalized Spehrical harmonics in the form
    d(P_{lm}(cos(theta)))dtheta * cos(m*phi)    if cs == 'c'
    d(P_{lm}(cos(theta)))dtheta  * sin(m*phi)    if cs == 's'
    
    theta and phi should be in meshgrid form with
    0<=theta<=pi is the colatitude
    -180<=phi<=180 is the longitude
    
    l, m = order and degree (scalars)
    '''
    if m == 0: FactorS = 1
    if m > 0: FactorS =np.multiply((-1)**m,np.sqrt(2*np.divide(math.factorial(l-m),float(math.factorial(l+m)))))
    if cs == 'c' : Trig = -m*np.sin(m*phi)
    if cs == 's' : Trig =  m*np.cos(m*phi)
    
    SH = FactorS*np.multiply(lpmv(m,l,np.cos(theta)),Trig)
    return SH;

def SHLS(d,th,ph,Lmax):
    '''
    least square calculation of spherical harmonics interpolation from scattered data
    INPUT:
        d=the data
        (th,ph)= the spatial coordinates (in colatitude and longitude)
            notation: d[i] =d(th[i],ph[i])
        Lmax = maximum SH degree 
    OUTPUT:
        beta = the set of coefficients ordered as:
            [l,m,glm,hlm]
        chi2 = the residual sum of squares misfits
    '''    

    Nd = len(d)
    Nl = 1+2*Lmax+Lmax**2
    if Nl > Nd: Nl = Nd
    Y = np.zeros((Nd,Nl))

    
    # prepare the coefficient matrix Y, such as:
    # d = Y*beta
    for i in range(Nd): # cycle over the rows: position
        l=0
        m=0
        cs='c'
        for j in range(Nl): # cycle over the columns degree and order of SH
            #print(l,m,cs)       
            Y[i,j] = SchmidtSH(l,m,th[i],ph[i],cs)
            # update indexes
            if l==m:
                if l==0:
                    l=1
                    m=0
                    cs='c'
                elif l>0:
                    if cs=='s':
                        m=0
                        l=l+1
                        cs='c'
                    elif cs=='c':
                        cs='s'
            elif m<l:
                if m==0:
                    m=1
                elif m>0:
                    if cs=='c':
                        cs='s'
                    elif cs=='s':
                        cs='c'
                        m=m+1
                        
    # least squares solution
    YTY = np.matmul(np.transpose(Y),Y)
    b = np.matmul(np.linalg.inv(YTY),np.matmul(np.transpose(Y),d))
    
    # calculate chi2 misfit
    e = np.matmul(Y,b)
    chi2 = 0
    for i in range(len(d)):
        chi2 = chi2 + (d[i]-e[i])**2 / (e[i])**2 
    
    
    # re-arrange coefficients
    Nentry = int( 1 + (l-1) + 0.5*( (l-1)*l ) + m )
    beta = np.zeros((Nentry,4))
    
    l=0
    m=0
    i=0
    for ientry in range(Nentry):
        if m==0:
            beta[ientry,:] = [l,m,b[i],0]
            i = i+1
        else:
            beta[ientry,:] = [l,m,b[i],b[i+1]]
            i = i+2
        # update indexes
        if m==l:
            l=l+1
            m=0
        else:
            m=m+1
        
    return beta, chi2;

def ForwardSH(beta,theta,phi):
    '''
    calculate a field F from its SH expansion beta
    INPUT
    beta = [l,m,glm,hlm]
    theta, phi = meshgrid matrices of the spatial coordinates
    0<theta<pi; 0<phi<2 pi
    OUTPUT
    F = sum(beta_lm*Y_lm(theta,phi))
    '''
    F=np.zeros(theta.shape)
    
    for i in range(beta.shape[0]):
        l = beta[i,0]
        m = beta[i,1]
        glm=beta[i,2]
        hlm=beta[i,3]
        
        F = F + glm * SchmidtSH(l,m,theta,phi,'c') + hlm * SchmidtSH(l,m,theta,phi,'s')
        
    return F;

def calcB(beta,theta,phi,a,r):
    '''
    calculate the geomagnetic field components from its SH expansion beta
    INPUT
    beta = [l,m,glm,hlm]
    theta, phi = meshgrid matrices of the spatial coordinates
    0<theta<pi; 0<phi<2 pi
    OUTPUT
    Br_a, Bt_a, Bp_a = magnetic field components at the Earth's surface
    Br_c, Br_c, Bp_c = magnetic field components at the CMB
    '''

    Br_r = np.zeros(theta.shape)
    Bt_r = np.zeros(theta.shape)
    Bp_r = np.zeros(theta.shape)

    for i in range(beta.shape[0]):
        l = beta[i,0]
        m = beta[i,1]
        glm=beta[i,2]
        hlm=beta[i,3]

        SH_s    = SchmidtSH(l,m,theta,phi,'s')
        SH_c    = SchmidtSH(l,m,theta,phi,'c')
        dthSH_s = DthSchmidtSH(l,m,theta,phi,'s')
        dthSH_c = DthSchmidtSH(l,m,theta,phi,'c')
        dphSH_s = DphSchmidtSH(l,m,theta,phi,'s')
        dphSH_c = DphSchmidtSH(l,m,theta,phi,'c')
        
        Br_r = Br_r + (l+1) *  (a/r)**(l+2) * ( glm * SH_c + hlm * SH_s )
        Bt_r = Bt_r - (a/r)**(l+2) * ( glm * dthSH_c + hlm * dthSH_s )
        Bp_r = Bp_r - (a/r)**(l+2) * ( glm * np.divide(dphSH_c,np.sin(theta)) + hlm * np.divide(dphSH_s,np.sin(theta)) )
         
    return Br_r, Bt_r, Bp_r;