# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 08:53:15 2020

@author: vgousyleblanc
"""

import numpy as np
import matplotlib.pyplot as plt


from scipy.optimize import curve_fit
#x,y,z,mx,my,mz=np.loadtxt('first_compensation.txt',unpack=True,skiprows=9)
#x=x/100
#y=y/100
#z=z/100

#%%

def differential_fit(x,y,z,mx,my,mz,dim): # define the dimension of the plane 
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
    for i in range(0,len(xu)):
        ind=np.where(x2==xu[i])
        #print(np.shape(ind))
        #print(arxu[i])
        mg_x[i,:]=mx2[ind]
        mg_y[i,:]=my2[ind]
        mg_z[i,:]=mz2[ind]
    if dim=='x':
        mean_low = np.mean(mg_x[0,:])
        mean_hi = np.mean(mg_x[-1,:])
        diff=mean_hi - mean_low
    if dim=='y':
        mean_low = np.mean(mg_y[:,0])
        mean_hi = np.mean(mg_y[:,-1])
        diff=mean_hi - mean_low
    if dim=='z':
        planeu=-0.2
        planed=-0.4
        planezu=np.where(z==planeu)
        planezd=np.where(z==planed)
        mxu=mx[planezu]
        myu=my[planezu]
        mzu=mz[planezu]
        mxd=mx[planezd]
        myd=my[planezd]
        mzd=mz[planezd]
        
        mg_xu=np.zeros((len(xu),len(yu)))
        mg_yu=np.zeros((len(xu),len(yu)))
        mg_zu=np.zeros((len(xu),len(yu)))
        mg_xd=np.zeros((len(xu),len(yu)))
        mg_yd=np.zeros((len(xu),len(yu)))
        mg_zd=np.zeros((len(xu),len(yu)))
        for i in range(1,len(xu)):
            ind=np.where(x2==xu[i])
            #print(np.shape(ind))
            #print(arxu[i])
            mg_xu[i,:]=mxu[ind]
            mg_yu[i,:]=myu[ind]
            mg_zu[i,:]=mzu[ind] 
            mg_zd[i,:]=mzd[ind]    
            
        diff=np.mean(mg_zd)-np.mean(mg_zu)
    
    return diff # we only want to have the differential plo

"""
x,y,z,mx,my,mz=np.loadtxt("spatial_scan_3_x_nG.txt",unpack=True,skiprows=9)
x=x/100
y=y/100
z=z/100
dim="x"
k[0]=differential_fit(x,y,z,mx,my,mz,dim)

x,y,z,mx,my,mz=np.loadtxt("spatial_scan_6_x_nG.txt",unpack=True,skiprows=9)
x=x/100
y=y/100
z=z/100
k[1]=differential_fit(x,y,z,mx,my,mz,dim)

x,y,z,mx,my,mz=np.loadtxt("spatial_scan_9_x_nG.txt",unpack=True,skiprows=9)
x=x/100
y=y/100
z=z/100
k[2]=differential_fit(x,y,z,mx,my,mz,dim)


coil=np.array([3,6,9])
plt.plot(coil,k,'.')
plt.show()
slope=linregress(coil,k)

coil1=-slope.intercept/slope.slope
coil2=fitm*coil1+fitb
plt.plot(coil,k,'.',label='data')
plt.plot(coil,coil*slope.slope+slope.intercept,label=f'Linear regression {dim}  direction')
plt.xlabel(f"coil 5 voltage")
plt.ylabel("differential [mG]")
plt.legend()
plt.show()
"""



