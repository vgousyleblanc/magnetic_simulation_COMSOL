# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:46:20 2020

@author: vgousyleblanc
"""

#Voltage scan !

# This file analyzes the result of a voltage scan
# You can configure what it does using `volt_config.py`
# The input should be a csv file with the following columns:
#
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import matplotlib as mpl
from matplotlib import pyplot as plt

def bilinear(X, m1, m2, b):
    v0, v1 = X
    return (m1*v0) + (m2*v1) + b


def bilinear_fit(v0,v1,mg, dimension: str, coila: int, coilb: int):

    p_up, c_up = curve_fit(bilinear, (v0, v1), mg)  # , sigma=s_up)
    
    e_up = np.sqrt(np.diag(c_up))
   
    sys.stdout.write(f" {GREEN}Done.{RESET} Parameters:\n")
    for (d, p, e) in (('Up', p_up, e_up), ('Down', p_dn, e_dn), ('Side', p_sd, e_sd)):
        print("    %s% 4s%s:  (% 5.3f ± %5.3f) ⋅ V_%i  +  (% 5.3f ± %5.3f) ⋅ V_%i  +  (% 5.3f ± %5.3f)  =  B%s" % (
            BOLD, d, RESET, p[0], e[0], coila, p[1], e[1], coilb, p[2], e[2], dimension
        ))

    sys.stdout.write(f"  {BOLD}Down{RESET} full compensation:")
    sys.stdout.write(" V_%i  =  %5.3f  +  %5.3f ⋅ V_%i\n" % (coilb, (-p_dn[2]/p_dn[1]), (-p_dn[0]/p_dn[1]), coila))


    v0s = tile(linspace(0, upper_lim + 0.5, 1024), 1024)
    v1s = repeat(linspace(0, upper_lim + 0.5, 1024), 1024)

    bup = bilinear((v0s, v1s), *p_up)
    bdn = bilinear((v0s, v1s), *p_dn)
    bsd = bilinear((v0s, v1s), *p_sd)

    err = abs(bdn)

    imin = argmin(err)

    print(f"  Minimum of fit error at {GREEN}({v0s[imin]}V, {v1s[imin]}V): {MAGENTA}{err[imin]} mG{RESET}.")

    bup.resize(1024, 1024)
    bdn.resize(1024, 1024)
    bsd.resize(1024, 1024)

    return bup, bdn, bsd, (-p_dn[0]/p_dn[1]), (-p_dn[2]/p_dn[1])

cmap_diverging = plt.get_cmap('RdBu')
cmap_normal    = plt.get_cmap('viridis')

def optimize(v0,v1,mgz):

    p_avg, c_avg = curve_fit(bilinear, (v0, v1), mgz)  # , sigma=s_up)
    
    e_avg = np.sqrt(np.diag(c_avg))
    print('The error on the parameter is :',e_avg) 
    v0u=np.unique(v0)
    v1u=np.unique(v1)
    print(p_avg)
    
    fitp = bilinear((v0, v1), *p_avg)
    fitm=(-p_avg[1]/p_avg[0])
    fitb=(-p_avg[2]/p_avg[0])
    print('the fitted coils voltage',fitm,fitb) 
    mg_z=np.zeros((len(v0u),len(v1u)))
    #=np.zeros((len(v0u),len(v1u)))
    for i in range(0,len(v0u)):
        ind=np.where(v0==v0u[i])
        #print(np.shape(ind))
        #print(arxu[i])
        #mg_x[i,:]=mx2[ind]
        #mg_y[i,:]=my2[ind]
        mg_z[i,:]=mgz[ind]
        
    #mg_z=mg_z[1:,:]*1000*1000
    #print(mg_z)
    rd_bu = plt.get_cmap('RdBu')
    amax = max(abs(mg_z.min()), abs(mg_z.max()))
    norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
    mapp = mpl.cm.ScalarMappable(norm, rd_bu)
    cmap = rd_bu
    
    
    fig, ax = plt.subplots(figsize=(8,6),dpi=192)
    mapp.set_array([])
    ax.imshow(mg_z, cmap=cmap,origin="lower",norm=norm,extent=[min(v1u),max(v1u),min(v1u),max(v1u)])#extent=[min(v0u[1:]),max(v0u),min(v1u[1:]),max(v1u)])
    #ax.plot(v0,b,'-.')
    #normalize the data
    
    x=[-3, -6]
    y=[-3*fitm+fitb,-6*fitm+fitb]
    ax.plot(x,y,c='black',ls='--',label='Fit')
    ax.set_xlabel(f"Voltage for Coil 1[V]")
    ax.set_ylabel(f"Voltage for Coil 2[V]")
    ax.set_title(f"$B_z$ [mG]")
    
    fig.colorbar(mapp)
    
    plt.show()
    return fitm,fitb

def spatial_scan(x,y,z,mx,my,mz):
    x=x/100
    y=y/100
    z=z/100
    
    
    plane=-0.3
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
#%%
file_name="z_compensation_test.txt"
v0,v1,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=5)
fitm,fitb=optimize(v0,v1,mz)