# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 09:33:05 2019

@author: mgropp
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import netCDF4 as nc
import scipy
from scipy import interpolate
import math

def roundup(x):
    return int(math.floor(x / 25.0)) * 25

# Ahhere
	 
 # another time
 # Just for testing

def regrid(org_lat, org_press, org_data, new_lat, new_press):
        X, Y = np.meshgrid(org_lat, org_press)
        XI, YI = np.meshgrid(new_lat, new_press)
        new_grid = scipy.interpolate.griddata((X.flatten(),Y.flatten()), org_data.flatten()
                        , (XI,YI))#, method='cubic')
        return new_grid


class Data:
	def __init__(self, var):
		self.var = var
		self.file_path = 'E:/dissertation'
		self.read_data()
		
	def read_data(self):
		for file_name in sorted(os.listdir(self.file_path)):
			if self.var in file_name:
				full_path = self.file_path + '/' + file_name
				self.data = nc.Dataset(full_path)
				
				
			
#u_wind, v_wind = Data('ua').data, Data('va').data

u_wind, a, b, p0, ps, sigma = Data('ua').data.variables['ua'][:], Data('ua').data.variables['a'][:],\
Data('ua').data.variables['b'][:],  Data('ua').data.variables['p0'][:], Data('ua').data.variables['ps'][:],  Data('ua').data.variables['lev'][:]
					
press_coords = np.zeros((40,90,144))

for i in range(90):
	for j in range(144):
			press_coords[:,i,j] = (a * p0 + b * ps[0,i,j]) / 100.


print(press_coords[:-16,0,0])
print(press_coords[:-16,10,20])

new_u_wind = np.zeros((26,90,144))

for i in range(90):
	for j in range(144):
			f = interpolate.interp1d(press_coords[:-16,i,j], u_wind[0,:-16,i,j])
			first = press_coords[0,i,j]
			last = 138.99999857
			new_press_coord = np.linspace(last, first, 25)[::-1]
			new_u_wind[:,i,j] = f(new_press_coord)
			
# Next step is very similar to the find closest height stuff from the tracking algorithm
# ORR interpolate 
