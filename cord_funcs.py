# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 07:23:24 2022

@author: jkoet
"""

import numpy as np

def dist2cord(xy=np.array([]),x=np.array([]),y=np.array([]),corner=np.array([0,0]),f=1):
    #coordinates to meters
        if np.any(xy):
            lon=xy[:,0]/(f*np.cos(corner[1]*np.pi/180)*111111)+corner[0] #meters
            lat=xy[:,1]/(f*111111)+corner[1] #meters
        else:
            lon=x/(f*np.cos(corner[1]*np.pi/180)*111111)+corner[0]#meters
            lat=y/(f*111111)+corner[1] #meters
        return np.vstack((lon,lat)).transpose()

def cord2dist(xy=np.array([]),x=np.array([]),y=np.array([]),corner=np.array([0,0]),f=1):
    #coordinates to meters
        if np.any(xy):
            x=(xy[:,0]-corner[0])*np.cos(corner[1]*np.pi/180)*111111 #meters
            y=(xy[:,1]-corner[1])*111111 #meters
            return np.vstack((x,y)).transpose()*f
        else:
            x=(x-corner[0])*np.cos(corner[1]*np.pi/180)*111111 #meters
            y=(y-corner[1])*111111 #meters
            return x*f,y*f
        