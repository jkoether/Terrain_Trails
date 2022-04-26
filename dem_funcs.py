# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 06:52:41 2022

@author: jkoet
"""
import rasterio as rio
import numpy as np

def getDem(importBounds,res):
    #Load USGS DEM data from 1/3 arcsecond (10m) data files
    elevFiles=[]
    for y in range(importBounds[1],importBounds[3]):
        for x in range(importBounds[0],importBounds[2]):
            if x<0:
                xStr=str(f'W{-x:03d}')
            else:
                xStr=str(f'E{x:03d}')
            if y<0:
                yStr=str(f'S{-y+1:02d}')
            else:
                yStr=str(f'N{y+1:02d}')
            elevFiles.append(yStr + xStr + ".tif")
    fd=[]
    x=[]
    y=[]
    for f in elevFiles:
        if res==10:
            overlap=6
            filename='USGS_13_'+f
        elif res==30:
            overlap=6
            filename='USGS_1_'+f
        else:
            raise NameError  ('resolution must be 10 or 30')
        dataset = rio.open('dem/' + filename)
        z=np.flip(dataset.read(1),axis=0)
        fd.append(z[overlap:-overlap,overlap:-overlap])
        x.append(round(dataset.bounds.left))
        y.append(round(dataset.bounds.top))
    return fd,x,y


class Dem:
    def __init__(self, poly,res,offset,dsf):
        importBounds=np.zeros([4],dtype='int16')
        print('')
        print('Loading Elevation Data...')
        
        #import a buffer around the border, needed for meshing.
        if res==10:
            buffer=10/10800
        else:
            buffer=10/3600
        offset=np.array(offset)
        
        offset[0]=offset[0]/(np.cos(poly[:,1].min()*np.pi/180)*111111)
        offset[1]=offset[1]/111111
        offset=np.flip(offset)#flip for coordinate convention (lat, lon)

        importBounds[0:2]=np.floor(np.min(poly,axis=0)-offset-buffer)
        importBounds[2:4]=np.ceil(np.max(poly,axis=0)-offset+buffer)
        importBounds=importBounds.astype(int)

        fd,x,y=getDem(importBounds,res)



        ux=np.unique(x)
        uy=np.unique(y)

        xIdx=[list(ux).index(i) for i in x]
        yIdx=[list(uy).index(i) for i in y]


        uy=uy-1 #why is it off by 1?

        # n=len(fd[0].heights)
        # m=len(fd[0].heights[0])
        n=fd[0].shape[0]
        m=fd[0].shape[1]


        #shift dem model by offset
        z=np.zeros([(n-1)*len(uy)+1,(m-1)*len(ux)+1],dtype='int16')
        
        
        self.lat=np.linspace(min(uy),max(uy)+1,z.shape[0])+offset[0]
        self.lon=np.linspace(min(ux),max(ux)+1,z.shape[1])+offset[1]
        



        # combine all elevation arrays (from seperate files) into one array.  each new array overlaps by one
        # row since the last row/column is the same as first row/column in the next set.
        for i in range(0,len(fd)):
            z[(n-1)*yIdx[i]:(n-1)*(yIdx[i]+1)+1,(m-1)*xIdx[i]:(m-1)*(xIdx[i]+1)+1]=fd[i]

        #downsample by average 2x2 (or larger) blocks
        if dsf>1:
            print('Downsample factor = ' + str(dsf))
            z=z[0:dsf*int(z.shape[0]/dsf),0:dsf*int(z.shape[1]/dsf)]
            z=z.reshape(z.shape[0],-1,dsf)
            z=np.sum(z,axis=2).T
            z=z.reshape(z.shape[0],-1,dsf)
            z=np.sum(z,axis=2).T/dsf**2

            self.lat=self.lat[0:dsf*int(self.lat.shape[0]/dsf)]
            self.lat=self.lat.reshape(-1,dsf)
            self.lat=np.sum(self.lat,axis=1)/dsf

            self.lon=self.lon[0:dsf*int(self.lon.shape[0]/dsf)]
            self.lon=self.lon.reshape(-1,dsf)
            self.lon=np.sum(self.lon,axis=1)/dsf

        self.z=z

        #Crop to requested area
        lon_range=range(np.argmax(self.lon > np.min(poly[:,0])-buffer)-1,np.argmax(self.lon > np.max(poly[:,0])+buffer)+1)
        # lon_range=(self.lon > np.min(poly[:,1])) & (self.lon < np.max(poly[:,1]))
        lat_range=range(np.argmax(self.lat > np.min(poly[:,1])-buffer)-1,np.argmax(self.lat > np.max(poly[:,1])+buffer)+1)
        # lat_range=(self.lat > np.min(poly[:,0])) & (self.lat < np.max(poly[:,0]))
        self.z=self.z[lat_range,:]
        self.z=self.z[:,lon_range] #crop to range that bounds area of interest.
        self.lat=self.lat[lat_range]
        self.lon=self.lon[lon_range]
        self.corner=[self.lon[0],self.lat[0]]



    def getElev(self,l):
        l=l

        #lon/lat steps for elevation data
        dy=(self.lat[-1]-self.lat[0])/(len(self.lat)-1)
        dx=(self.lon[-1]-self.lon[0])/(len(self.lon)-1)

        #advanced indexing lookup
        z=self.z[np.round((l[:,1]-self.lat[0])/dy).astype('int16'),
                 np.round((l[:,0]-self.lon[0])/dx).astype('int16')]

        return z


    def plotElev(self):
        fig = go.Figure(data=go.Heatmap(
        x=self.lon,
        y=self.lat,
        z=self.z,
        colorscale='Viridis'))
        fig.add_trace(go.Scatter(x=[-80.326680], y=[40.740028],
                              marker = dict(
                                  symbol = 'star',
                                  size = 15)))
        fig.update_layout(
            yaxis = dict(
                scaleanchor = "x",
                scaleratio = 1),
            width=1000, height=1000,
            title='title')
        fig.show()

    def plotElev2(self):
        plt.imshow(np.flip(self.z,axis=0), cmap='viridis')
        plt.scatter(x=[-80.326680], y=[40.740028],marker='*', c='r', s=40)
        plt.colorbar()
        plt.show()
        plt.pause(0.1)