# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 09:05:45 2020

@author: vgousyleblanc
"""

import mph
import numpy as np
from voltage_scan import spatial_scan
from voltage_scan import optimize
from differential_fit import differential_fit
from scipy.stats import linregress
import time
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
client = mph.Client()
model = client.load('magnetic_field_compensation_geo_extended_2.mph')


#model.parameter('z_1_p', '50[cm]')
#model.parameter('z_2_p', '50[cm]')

#model.parameter('x_2_p', '50[cm]')
#model.parameter('x_2_p', '50[cm]')

#model.parameter('y_1_p', '50[cm]')
#model.parameter('y_2_p', '50[cm]')

model.parameter('V_z1', '0[V]')
model.parameter('V_z2', '0[V]')

model.parameter('V_x12', '0[V]')
model.parameter('V_x22', '0[V]')

model.parameter('V_y1', '-2.09[V]')
model.parameter('V_y2', '-2.16[V]')

#%%
model.mesh()
start = time.time()
scan='v'
dim='z'

if scan=='v' and dim=='y':
    print('Starting the voltage scan in the '+str(dim)+'  dimension')
    model.solve('Y_sweep')
    file_name='voltage_scan_y.txt'
    model.export('Y_table',file_name)
    v0,v1,mx,my,mz=np.loadtxt("voltage_scan_y_vf_2.txt",unpack=True,skiprows=5)
    fitm,fitb=optimize(v0,v1,my)
    computation_time = time.time()/60
    print('Time taken to compute the result',computation_time)
if scan=='v' and dim=='x':
    print('Starting the voltage scan in the '+str(dim)+'  dimension')
    #print('Starting the voltage scan')
    model.solve('X_sweep')
    file_name='voltage_scan_x_vf.txt'
    v0,v1,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=5)
    fitm,fitb=optimize(v0,v1,mx)
if scan=='v' and dim=='z':
    model.solve('Z_sweep')
    print('Starting the voltage scan in the '+str(dim)+'  dimension')
    file_name='voltage_scan_z_vf.txt'
    #model.export('Table 1',file_name)
    v0,v1,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=20)
    fitm,fitb=optimize(v0,v1,mz)
#%%    
scan='s'
if scan=='s':   
    
        if dim=='x':
            print('Starting the spatial scan in the '+str(dim))
            v0=np.array([0,1,2])   
            v1=np.array([fitm*v0[0]+fitb,fitm*v0[1]+fitb,fitm*v0[2]+fitb])
            k=np.zeros(len(v0))
            for i in range(0,len(v0)):
                model.parameter('V_x12', str(v0[i])+'[V]')
                model.parameter('V_x22', str(v1[i])+'[V]')  
                model.solve('Study 1')
                file_name='spatial_scan_'+str(dim)+str(v0[i])+'.txt'
                model.export('Full_compensation',file_name)
                x,y,z,mx,my,mz=np.loadtxt('spatial_scan_'+str(dim)+str(v0[i])+'.txt',unpack=True,skiprows=9)
                x=x/100
                y=y/100
                z=z/100
            #dim="x"
                k[i]=differential_fit(x,y,z,mx,my,mz,dim)
            slope=linregress(v0,k)
            coil1=-slope.intercept/slope.slope
            coil2=fitm*coil1+fitb
            print('Tension in coil 1',coil1)
            print('Tension in coil 2',coil2)
            model.parameter('V_x12', str(coil1)+'[V]')
            model.parameter('V_x22', str(coil2)+'[V]')  
            model.solve('Study 1')
            file_name='full_compensation'+str(dim)+'.txt'
            model.export('Full_compensation',file_name)
            x,y,z,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=9)
        
            spatial_scan(x,y,z,mx,my,mz)
        
        if dim=='y':
            print('Starting the spatial scan in the '+str(dim))
            v0=np.array([0,1,2])   
            v1=np.array([fitm*v0[0]+fitb,fitm*v0[1]+fitb,fitm*v0[2]+fitb])
            k=np.zeros(len(v0))
            for i in range(0,len(v0)):
                model.parameter('V_y1', str(v0[i])+'[V]')
                model.parameter('V_y2', str(v1[i])+'[V]')  
                model.solve('Study 1')
                file_name='spatial_scan_'+str(dim)+str(v0[i])+'.txt'
                model.export('Full_compensation',file_name)
                x,y,z,mx,my,mz=np.loadtxt('spatial_scan_'+str(dim)+str(v0[i])+'.txt',unpack=True,skiprows=9)
                x=x/100
                y=y/100
                z=z/100
            #dim="x"
                k[i]=differential_fit(x,y,z,mx,my,mz,dim)
            slope=linregress(v0,k)
            coil1=-slope.intercept/slope.slope
            coil2=fitm*coil1+fitb
            print('Tension in coil 1',coil1)
            print('Tension in coil 2',coil2)
            model.parameter('V_y1', str(coil1)+'[V]')
            model.parameter('V_y2', str(coil2)+'[V]')  
            model.solve('Study 1')
            file_name='full_compensation'+str(dim)+'.txt'
            model.export('Full_compensation',file_name)
            x,y,z,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=9)
            spatial_scan(x,y,z,mx,my,mz)
            end = time.time()
        if dim=='z':
            v0=np.array([0,1,2])   
            v1=np.array([fitm*v0[0]+fitb,fitm*v0[1]+fitb,fitm*v0[2]+fitb])
            k=np.zeros(len(v0))
            for i in range(0,len(v0)):
                model.parameter('V_z1', str(v0[i])+'[V]')
                model.parameter('V_z2', str(v1[i])+'[V]')  
                model.solve('Study 1')
                file_name='spatial_scan_'+str(dim)+str(v0[i])+'.txt'
                model.export('Full_compensation',file_name)
                x,y,z,mx,my,mz=np.loadtxt('spatial_scan_'+str(dim)+str(v0[i])+'.txt',unpack=True,skiprows=9)
            #dim="x"
                k[i]=differential_fit(x,y,z,mx,my,mz,dim)
            slope=linregress(v0,k)
            coil1=-slope.intercept/slope.slope
            coil2=fitm*coil1+fitb
            print('Tension in coil 1',coil1)
            print('Tension in coil 2',coil2)
            
            model.parameter('V_z1', str(coil1)+'[V]')
            model.parameter('V_z2', str(coil2)+'[V]')  
            model.solve('Study 1')
            file_name='full_compensation'+str(dim)+'.txt'
            model.export('Full_compensation',file_name)
            x,y,z,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=9)
          
            spatial_scan(x,y,z,mx,my,mz)
end=time.time()
print('final computation time', start-end)     
    #x,y,z,mx,my,mz=np.loadtxt(file_name,unpack=True,skiprows=20)
    #spatial_scan(x,y,z,mx,my,mz)
    




    

#%%



v0,v1,mx,my,mz=np.loadtxt("voltage_scan_test_dup.txt",unpack=True,skiprows=5)
fitm,fitb=optimize(v0,v1,my)





#%%
client.remove(model)



