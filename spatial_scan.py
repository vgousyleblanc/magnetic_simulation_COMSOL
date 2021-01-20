# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:52:30 2020

@author: vgousyleblanc
"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

x,y,z,mx,my,mz=np.loadtxt('spatial_scan_3inch.txt',unpack=True,skiprows=9)
#convert=1e4*1e3
#convert_g=12.54

def spatial_scan_2(x,y,z,mx,my,mz,plane):
    x=x/100
    y=y/100
    z=z/100
    
    plane1z=np.where(z==plane)#np.where(z==plane)
    
    xu=np.unique(np.round(x[plane1z],3))
    yu=np.unique(y[plane1z])
    zu=np.unique(z[plane1z])
    
    x2=np.round(x[plane1z],3)
    y2=y[plane1z]
    mx2=mx[plane1z]
    my2=my[plane1z]
    mz2=mz[plane1z]
    
    mg_x=np.zeros((len(xu),len(yu)))
    mg_y=np.zeros((len(xu),len(yu)))
    mg_z=np.zeros((len(xu),len(yu)))
    j=np.where(z==-0.3)
    for i in range(0,len(xu)):
        ind=np.where(x2==xu[i])
        #print(np.shape(ind))
        #print(arxu[i])
        mg_y[i,:]=my2[ind]
        mg_x[i,:]=mx2[ind]
        mg_z[i,:]=mz2[ind]
        #print(np.size(ind))
    rd_bu = plt.get_cmap('RdBu')
    #mg_x=mg_x-mg_xw
    amax_n = max(abs(mg_x.min()), abs(mg_x.max()))
    norm = mpl.colors.Normalize(vmin=-amax_n, vmax=amax_n)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    
    ax.imshow(mg_x, cmap=cmap, norm=norm, origin='lower',extent=[min(xu),max(xu),min(yu),max(yu)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_x$ (z={plane} m ) [mG]")
    fig.colorbar(mapp)
    
    plt.show()
    #mg_y=mg_y-mg_yw
    amax_ny = max(abs(mg_y.min()), abs(mg_y.max()))
    norm = mpl.colors.Normalize(vmin=-amax_ny, vmax=amax_ny)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    ax.imshow(mg_y, cmap=cmap, norm=norm, origin='lower',extent=[min(xu),max(xu),min(yu),max(yu)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_y$ (z={plane}m  ) [mG]")
    fig.colorbar(mapp)
    
    plt.show()
    #mg_z=mg_z-mg_zw
    amax= max(abs(mg_z.min()), abs(mg_z.max()))
    norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    ax.imshow(mg_z, cmap=cmap,norm=norm, origin='lower',extent=[min(xu),max(xu),min(yu),max(yu)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_z$ (z={plane}m) [mG]")
    fig.colorbar(mapp)
    
    plt.show()
    return mg_x,mg_y,mg_z,xu
#%%
mg_x1,mg_y1,mg_z1,xu1=spatial_scan_2(x,y,z,mx,my,mz,-0.2)

mg_x2,mg_y2,mg_z2,xu1=spatial_scan_2(x,y,z,mx,my,mz,-0.3)

def difference(mg_x,mg_y,mg_z):
    
    rd_bu = plt.get_cmap('RdBu')
    amax_ny = max(abs(mg_y.min()), abs(mg_y.max()))
    norm = mpl.colors.Normalize(vmin=-amax_ny, vmax=amax_ny)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    ax.imshow(mg_y, cmap=cmap, norm=norm, origin='lower',extent=[min(xu1),max(xu1),min(xu1),max(xu1)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_y$ [mG]")
    fig.colorbar(mapp)
    
    plt.show()
    #mg_z=mg_z-mg_zw
    amax= max(abs(mg_z.min()), abs(mg_z.max()))
    norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    ax.imshow(mg_z, cmap=cmap,norm=norm, origin='lower',extent=[min(xu1),max(xu1),min(xu1),max(xu1)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_z$[mG]")
    fig.colorbar(mapp)
    
    
    plt.show()
    
    amax_n = max(abs(mg_x.min()), abs(mg_x.max()))
    norm = mpl.colors.Normalize(vmin=-amax_n, vmax=amax_n)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    
    ax.imshow(mg_x, cmap=cmap, norm=norm, origin='lower',extent=[min(xu1),max(xu1),min(xu1),max(xu1)])
    ax.set_xlabel(f"X [m]")
    ax.set_ylabel(f"Y [m]")
    ax.set_title(f"$B_x$ [mG]")
    fig.colorbar(mapp)
    
    plt.show()
    
difference(mg_x1-mg_x2,mg_y1-mg_y2,mg_z1-mg_z2)


