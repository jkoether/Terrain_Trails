# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:50:30 2021

@author: Jeremy Koether
"""


import time as time
import inflator
import os as os
from copy import deepcopy

import overpy
import numpy as np

import plotly.graph_objects as go

from scipy.spatial import Delaunay

import trimesh as tm
import pyclipper as pc
import triangle as tr

import gpxpy as gpxpy

import rasterio as rio
from scipy import ndimage
from subprocess import run



from zipfile import ZipFile

import matplotlib.pyplot as plt




api = overpy.Overpass()


def get10dem(importBounds):
    #Load USGS DEM data from 1/3 arcsecond (10m) data files
    elevFiles=[]
    demPath="D:/GDEM/ASTGTMV003_"
    for y in range(importBounds[0],importBounds[2]):
        for x in range(importBounds[1],importBounds[3]):
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

        filename='USGS_13_'+f
        dataset = rio.open('dem/' + filename)
        fd.append(np.flip(dataset.read(1),axis=0))
        x.append(round(dataset.bounds.left))
        y.append(round(dataset.bounds.top))
    return fd,x,y

def get30dem(importBounds):
    #Load USGS DEM data from 1 arcsecond (30m) data files
    elevFiles=[]
    demPath="D:/GDEM/ASTGTMV003_"
    for y in range(importBounds[0],importBounds[2]):
        for x in range(importBounds[1],importBounds[3]):
            if x<0:
                xStr=str(f'W{-x:03d}')
            else:
                xStr=str(f'E{x:03d}')
            if y<0:
                yStr=str(f'S{-y:02d}')
            else:
                yStr=str(f'N{y:02d}')
            elevFiles.append(demPath + yStr + xStr + ".zip")
    fd=[]
    x=[]
    y=[]
    for f in elevFiles:
        drive, path = os.path.splitdrive(f)
        path, filename = os.path.split(path)
        filename, file_extension = os.path.splitext(filename)
        filename=filename+'_dem.tif'
        if not os.path.isfile('dem/' + filename):
            with ZipFile(f, 'r') as zipObj:

               # Extract all the contents of zip file in different directory
               # zipObj.extractall('temp')
               listOfFileNames = zipObj.namelist()
               for f2 in listOfFileNames:
                   zipObj.extract(f2,path='temp')
        dataset = rio.open('dem/' + filename)
        fd.append(np.flip(dataset.read(1),axis=0))
        x.append(round(dataset.bounds.left))
        y.append(round(dataset.bounds.top))
        return fd,x,y

class Dem:
    def __init__(self, poly,res,offset,dsf):
        importBounds=np.zeros([4],dtype='int16')
        print('')
        print('Loading Elevation Data...')
        offset[0]=offset[0]/(np.cos(poly[:,0].min()*np.pi/180)*111111)
        offset[1]=offset[1]/111111

        importBounds[0:2]=np.floor(np.min(poly,axis=0)-offset)
        importBounds[2:4]=np.ceil(np.max(poly,axis=0)-offset)
        importBounds=importBounds.astype(int)

        if res==10:
            fd,x,y=get10dem(importBounds)
        elif res==30:
            fd,x,y=get30dem(importBounds)
        else:
            raise NameError('resolution must be 10 or 30')

        ux=np.unique(x)
        uy=np.unique(y)

        xIdx=[list(ux).index(i) for i in x]
        yIdx=[list(uy).index(i) for i in y]
        # yIdx=yIdx[::-1]

        uy=uy-1 #why is it off by 1?

        # n=len(fd[0].heights)
        # m=len(fd[0].heights[0])
        n=fd[0].shape[0]
        m=fd[0].shape[1]

        #shift dem model by offset
        z=np.zeros([(n-1)*len(uy)+1,(m-1)*len(ux)+1],dtype='int16')
        self.lat=np.linspace(min(uy),max(uy)+1,z.shape[0])+offset[0]
        self.lon=np.linspace(min(ux),max(ux)+1,z.shape[1])+offset[1]
        # offset[0]=offset[0]/(np.cos(self.lat.min()*np.pi/180)*111111)
        # offset[1]=offset[1]/111111



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
        lon_range=range(np.argmax(self.lon > np.min(poly[:,1]))-1,np.argmax(self.lon > np.max(poly[:,1]))+1)
        # lon_range=(self.lon > np.min(poly[:,1])) & (self.lon < np.max(poly[:,1]))
        lat_range=range(np.argmax(self.lat > np.min(poly[:,0]))-1,np.argmax(self.lat > np.max(poly[:,0]))+1)
        # lat_range=(self.lat > np.min(poly[:,0])) & (self.lat < np.max(poly[:,0]))
        self.z=self.z[lat_range,:]
        self.z=self.z[:,lon_range] #crop to range that bounds area of interest.
        self.lat=self.lat[lat_range]
        self.lon=self.lon[lon_range]




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





def cord2dist(xy=np.array([]),x=np.array([]),y=np.array([]),corner=np.array([0,0]),f=1):
    #coordinates to meters
        if np.any(xy):
            xy[:,0]=(xy[:,0]-corner[0])*np.cos(np.min(xy[:,1])*np.pi/180)*111111 #meters
            xy[:,1]=(xy[:,1]-corner[1])*111111 #meters
            return xy*f
        else:
            x=(x-corner[0])*np.cos(np.min(y)*np.pi/180)*111111 #meters
            y=(y-corner[1])*111111 #meters
            return x*f,y*f


def mergeWays(ways):


    n=0
    for w in ways:
        n=n+len(w.get_nodes())

    print("{:,}".format(n))
    verts =np.full([n,2], np.nan,float) #coords of unique points
    edges = []
    pID=np.zeros(n,dtype='int64') #list of all unique open streetmap point IDs.
    nConn2 = [ [] for _ in range(n)]

    n=0
    #fig3, ax2 = plt.subplots()


    for w in ways:

        #print(n)
        # c.plot(ax)
        l=len(w.get_nodes())
        v = np.zeros((l,2), float)


        v[:,0]=[float(n.lon) for n in w.nodes]
        v[:,1]=[float(n.lat) for n in w.nodes]

        ids=np.array(w._node_ids)

        #intially set up up edges between subsequent points, relative to start of this path
        e=np.arange(0,l-1)
        e=np.vstack((e,e+1))

        while np.unique(ids).shape[0]!=ids.shape[0]:   #look for duplicate points within way (self intersections)

            s=np.sort(ids)
            dup=s[np.argmax(np.diff(s)==0)]
            idx=np.where(ids==dup)[0]

            #idx[1] will be deleted, first occurance will be used.
            e[e==idx[1]]=idx[0]
            e[e>idx[1]]=e[e>idx[1]]-1
            ids=np.delete(ids,idx[1],0)
            v=np.delete(v,idx[1],0)
            l=l-1
        while any(np.isin(ids,pID)): #loop through all duplicates (intersections with other ways)
            if np.isin(105225537,ids):
                stop=1
            j=np.where(np.isin(ids,pID)) #which point duplicates an existing point.
            j=j[0][0] #next first duplicate
            idx=np.where(np.isin(pID,ids[j])) #duplicated existing point
            e[e==j]=idx[0]-n
            e[e>j]=e[e>j]-1
            ids=np.delete(ids,j,0)
            v=np.delete(v,j,0)
            l=l-1

        # print(w._node_ids)
        # print(np.vstack(([i.id for i in w.nodes[0:-1]],[i.id for i in w.nodes[1:]])))
        # ax.plot(nLoc[n:n+l,1],nLoc[n:n+l,0],'c-')
        pID[n:n+l]=ids
        verts[n:n+l,:]=v
        edges.append(np.transpose(e+n))
        n=n+l
        if np.any(np.diff(e,axis=1)==0):
            print("*")
            print(e)
        #print(str(np.argmax(np.isnan(verts[:,0]))-1) + " - " + str(np.vstack(edges).max()))

    #remove unused elements (motorways)
    pID=pID[:n]
    edges=np.vstack(edges)

    verts=verts[:n,:]

    return verts[:n,:],edges

def addLake(ways,corner,scale_factor,clearance,base,height_factor,dem,n):

    #create extruded polygon for lakes, constant thickness.
    verts,edg=mergeWays(ways)

    z=dem.getElev(verts) #perimeter elevation profile
    z=(z-np.min(dem.z))*scale_factor*height_factor+base
    z=np.sort(z)[int(0.25*z.shape[0])] #use the 25th percentile value as water level

    verts=cord2dist(xy=verts,corner=corner,f=scale_factor)
    pco = pc.PyclipperOffset()
    pco.AddPath(verts*1000, 0,pc.ET_CLOSEDPOLYGON)
    #individual stl
    p=np.vstack(pco.Execute(0)).astype('double')/1000

    edg=inflator.seqEdges(p.shape[0])
    face = {
    "vertices": p,
    "segments": edg}

    t = tr.triangulate(face,'p')
    mesh=tm.creation.extrude_triangulation(t['vertices'], t['triangles'],3)

    mesh.export("print_files/waterbody_"+str(n)+".stl") #output individual lakes

    #negatives
    p=np.vstack(pco.Execute(clearance*1000)).astype('double')/1000

    face = {
    "vertices": p,
    "segments": inflator.seqEdges(p.shape[0])}

    t = tr.triangulate(face,'p')
    mesh=tm.creation.extrude_triangulation(t['vertices'], t['triangles'],200)

    # if n>1: #combine all lake negatives to one STL
    #negative for waterways (full height)
    mesh.export("temp/waterbody_N1_"+str(n)+".stl")
    # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o temp\waterbody_N1.stl Waterbody_N1_comb.scad"')

    #negative for terrain (3mm deep)
    translation = [0,0,z-3]  # mean water height-3
    mesh.apply_translation(translation)
    mesh.export("temp/waterbody_N2_"+str(n)+".stl")

    return verts

def getFootpaths(poly,trail_exclude,corner,scale_factor,p_width,s_width,clearance,base):
    print(' ')
    print('Retrieving OSM path nodes...')

    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,0]),np.min(poly[:,1]),np.max(poly[:,0]),np.max(poly[:,1])))
    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            result = api.query("way(" + areaStr + ") [""highway""];    (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)

    inc_ways=[]
    for w in result.ways:
        if w.tags['highway'] in ['path','footway','cycleway'] and not w.id in trail_exclude:
            inc_ways.append(w)
    verts,edg=mergeWays(inc_ways)
    verts=cord2dist(xy=verts,corner=corner,f=scale_factor)
    pMesh(verts,edg,p_width,s_width,clearance,base,'trail')
    return dict([('verts', verts), ('edges', edg)])


def getRoads(poly,roads,corner,scale_factor,p_width,s_width,clearance,base):
    print(' ')
    print('Retrieving OSM road nodes...')

    rd_names=[]
    if any(roads):
        for i in range(len(roads)-1,-1,-1):
            if isinstance(roads[i],str):
                rd_names.append(roads.pop(i))

    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,0]),np.min(poly[:,1]),np.max(poly[:,0]),np.max(poly[:,1])))
    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            result = api.query("way(" + areaStr + ") [""highway""];    (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)
    inc_ways=[]
    for w in result.ways:
        if w.id in roads or ("name" in w.tags and w.tags["name"] in rd_names):
            inc_ways.append(w)
    verts,edg=mergeWays(inc_ways)
    verts=cord2dist(xy=verts,corner=corner,f=scale_factor)
    pMesh(verts,edg,p_width,s_width,clearance,base,'road')
    return dict([('verts', verts), ('edges', edg)])


def getWaterways(poly,waterways,corner,scale_factor,p_width,s_width,clearance,base):
    print(' ')
    print('Retrieving OSM waterway path nodes...')

    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,0]),np.min(poly[:,1]),np.max(poly[:,0]),np.max(poly[:,1])))

    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            result = api.query("way(" + areaStr + ") [""waterway""];    (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)
    inc_ways=[]

    for w in result.ways:
        if  w.id in waterways or ("name" in w.tags and w.tags["name"] in waterways):
            inc_ways.append(w)
    if any(inc_ways):
        verts,edg=mergeWays(inc_ways)
        verts=cord2dist(xy=verts,corner=corner,f=scale_factor)
        pMesh(verts,edg,p_width,s_width,clearance,base,'waterway')
        return dict([('verts', verts), ('edges', edg)])
    else:
        return []

def getWaterbodies(poly,bodies,corner,scale_factor,clearance,base,height_factor,dem):

    print(' ')
    print('Retrieving OSM waterbody nodes...')
    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,0]),np.min(poly[:,1]),np.max(poly[:,0]),np.max(poly[:,1])))
    wb_out=[]

    n=1

    # Lakes can be either relations or individual ways
    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            result = api.query("rel(" + areaStr + ") [""water""];    (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)

    for r in result.relations:
        #find by OSM id or name
        if r.tags['name'] in bodies or r.id in bodies:
            #once we find a relation, extract and merge all the ways
            ways=[]
            for m in r.members:
                if m._type_value == 'way' and m.role == 'outer':
                    num_attempt=0
                    result=[]
                    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
                        try:
                            result = api.query("way(id:"+str(m.ref)+");    (._;>;); out body;")
                        except:
                            num_attempt=num_attempt+1
                            time.sleep(5)
                    ways.append(result.ways[0])
            wb_out.append(addLake(ways,corner,scale_factor,clearance,base,height_factor,dem,n))
            n=n+1
    #now look for lakes made of individual ways
    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            result = api.query("way(" + areaStr + ") [""water""];    (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)
    for w in result.ways:
        if ('name' in w.tags and w.tags['name'] in bodies) or w.id in bodies:
            wb_out.append(addLake([w],corner,scale_factor,clearance,base,height_factor,dem,n))
            n=n+1
    return wb_out



def tMesh(dem,edge_poly,max_Size,height_factor,base_height,water_drop):
    #create terrian mesh that bounds lat/lon polygon.  polygon shape and all path
    #features will be cut out using boolean operations
    x=dem.lon
    y=dem.lat


    corner=np.stack((np.min(x),np.min(y)))

    p=cord2dist(np.flip(deepcopy(edge_poly),1),corner=corner)
    d=np.sqrt(np.sum((p[0,:]-p)**2,axis=1))
    a=np.arctan2(p[:,1]-p[0,1],p[:,0]-p[0,0])
    scale=np.zeros(91)
    for i in range(91):
        #Rotate points in one degee increments, and determine the min scale factor between the two directions
        p=np.vstack((d*np.cos(a+i*np.pi/180),d*np.sin(a+i*np.pi/180))).T

        scale[i]=min([max(max_Size)/(p[:,0].max()-p[:,0].min()),
                      min(max_Size)/(p[:,1].max()-p[:,1].min())]) #scale factor for each angle
    print('print angle: ' +str(np.argmax(scale)))

    scale=scale.max() #find angle that allows the largest print.

    x,y=cord2dist(x=x,y=y,corner=corner)


    #convert meters on map to mm on print.
    x=x*scale
    y=y*scale

    print('DEM resolution: {0:.2f}x{1:.2f} mm'.format(((x.max()-x.min())/(x.shape[0]-1)),(y.max()-y.min())/(y.shape[0]-1)))
    x,y=np.meshgrid(x,y)
    x=x.flatten()
    y=y.flatten()
    z=dem.z.flatten().transpose()*height_factor
    z=(z-z.min())*scale+base_height+1 #+1 provides 1mm extra thickness for min path hieght

    verts=np.vstack((x,y)).transpose()


    tri = Delaunay(verts) #triangulate in 2D

    verts=np.vstack((x,y,z)).transpose() #inflate to 3D


    verts_bot=np.vstack((x,y,np.zeros(z.shape[0], dtype=z.dtype))).transpose() #bottom surface same as top but all z=0.



    faces1= tri.simplices;
    faces1=np.vstack((faces1,faces1+verts.shape[0])) #assign elevation to top surface of mesh


    #perimeter points, CCW starting at 0,0
    s1=np.flatnonzero(verts[:,0]==np.min(x)) #x=0 side
    s2=np.flatnonzero(verts[:,1]==np.max(y)) #y=max side
    s3=np.flip(np.flatnonzero(verts[:,0]==np.max(x)))#x=max side
    s4=np.flip(np.flatnonzero(verts[:,1]==np.min(y))) #y=0 side
    perimeter_vert=np.hstack((s1,s2[1:],s3[1:],s4[1:])) #combine all, first point repeated

    perimeter_faces=np.hstack((
        np.vstack((perimeter_vert[0:-1], perimeter_vert[1:],perimeter_vert[0:-1]+verts.shape[0])),
        np.vstack((perimeter_vert[1:], perimeter_vert[0:-1]+verts.shape[0],perimeter_vert[1:]+verts.shape[0]))))

    faces1=np.vstack((faces1,perimeter_faces.transpose()))
    verts=np.vstack((verts,verts_bot))

    msh = tm.Trimesh(vertices=verts,
                       faces=faces1)
    tm.repair.fix_normals(msh)
    msh.export('temp/terrain.stl')

    #used for waterways
    verts2=verts
    verts2[:,2]=verts2[:,2]-water_drop
    msh = tm.Trimesh(vertices=verts2,
                       faces=faces1)
    tm.repair.fix_normals(msh)
    msh.export('temp/terrain_low1.stl')

   #used to create thicker top section for paths
    verts2=verts
    verts2[:,2]=verts2[:,2]-1
    msh = tm.Trimesh(vertices=verts2,
                       faces=faces1)
    tm.repair.fix_normals(msh)
    msh.export('temp/terrain_low2.stl')


    return scale,corner


def plotPaths(dem,sf,wb,fp,rp,ww,p):


    #Create map figure to show terrain an slected paths.
    fig4, ax = plt.subplots()
    x=dem.lon
    y=dem.lat
    corner=np.stack((np.min(x),np.min(y)))
    x,y=cord2dist(x=x,y=y,corner=corner)
    x=x*sf
    y=y*sf
    grad=ndimage.sobel(dem.z)
    plt.imshow(grad, cmap ='pink',
                 extent =[x.min(), x.max(), y.min(), y.max()],
                    origin ='lower')

    for i in range(fp["edges"].shape[0]):
        ax.plot([fp["verts"][fp["edges"][i,0],0],fp["verts"][fp["edges"][i,1],0]],[fp["verts"][fp["edges"][i,0],1],fp["verts"][fp["edges"][i,1],1]],'g-')
    for i in range(rp["edges"].shape[0]):
        ax.plot([rp["verts"][rp["edges"][i,0],0],rp["verts"][rp["edges"][i,1],0]],[rp["verts"][rp["edges"][i,0],1],rp["verts"][rp["edges"][i,1],1]],'k-')
    if any(ww):
        for i in range(ww["edges"].shape[0]):
            ax.plot([ww["verts"][ww["edges"][i,0],0],ww["verts"][ww["edges"][i,1],0]],[ww["verts"][ww["edges"][i,0],1],ww["verts"][ww["edges"][i,1],1]],'b-')
    for i in range(p.shape[0]):
        ax.plot(np.hstack((p[:,0],p[0,0])),np.hstack((p[:,1],p[0,1])),'r:')

    ax.plot(p[0,0],p[0,1],'kx')

    for w in wb:
        ax.fill(w[:,0],w[:,1],facecolor='teal',alpha=0.75)
    # ax.fill(wb[:,0],wb[:,1],facecolor='teal',alpha=0.75)
    #plt.title(title)
    plt.axis('equal')
    plt.show()
    plt.pause(0.1)


def pMesh(verts,edges,p_width,s_width,clearance,base,fname):
    translation = [0,0,base]  # box offset + plane

    path_mesh=inflator.InflateNetwork(verts,edges,s_width/2)
    path_mesh=tm.creation.extrude_triangulation(path_mesh['vertices'], path_mesh['triangles'], 100)

    path_mesh.apply_translation(translation)
    path_mesh.export('temp/'+fname + '_1.stl')

    path_mesh=inflator.InflateNetwork(verts,edges,(p_width)/2)
    path_mesh=tm.creation.extrude_triangulation(path_mesh['vertices'], path_mesh['triangles'], 100)
    path_mesh.apply_translation(translation)
    path_mesh.export('temp/'+fname + '_2.stl')

    path_mesh=inflator.InflateNetwork(verts,edges,(p_width+clearance)/2)
    path_mesh=tm.creation.extrude_triangulation(path_mesh['vertices'], path_mesh['triangles'], 100)
    tm.repair.fill_holes(path_mesh)
    path_mesh.apply_translation(translation)
    path_mesh.export('temp/'+fname + '_3.stl')


def genPoly(p,offset):

    pco = pc.PyclipperOffset()
    pco.AddPath(p*1000, 0,pc.ET_CLOSEDPOLYGON)

    p=np.vstack(pco.Execute(offset*1000)).astype('double')/1000

    face = {
      "vertices": p,
      "segments": inflator.seqEdges(p.shape[0])}

    t = tr.triangulate(face,'p')

    mesh=tm.creation.extrude_triangulation(t['vertices'], t['triangles'], 200)
    # translation = [0,0,2]  # box offset + plane offset
    # mesh.apply_translation(translation)
    return mesh


def GenerateSTLs(CoordPoly,rd_include=[],trail_exclude=[],waterway_include=[],waterbody=[],
                 path_width=1.2,support_width=0.9,path_clearance=0.1,height_factor=2,base_height=4,
                 edge_width=1.5,max_print_size=[248,198],water_drop=0.5,load_area=[],resolution=30,
                 dem_offset=[0,0],downsample_factor=1,map_only=False):
    """

    INPUTS:

    CoordPoly - polygon to define shape of terrain, input as longitude/lattitude array, or gpx file of points
    rd_include - roads names / ids to inlcude, no roads included by default.
    trail_exclude - footpaths names / ids to inlcude, all included by default.
    waterway_include - water paths to inlcude, not inlcuding polygon water bodies
    path_width - path (road, footpaths and waterway) print top width, 3 total passes works best.
    support_width - width of base (one pass less than top)
    path_clearance - clearance between path prints and cutouts.
    height_factor - exaggeration factor for elevation. 1.0 makes on same scale as horizontal dimensions.
    base_height - minium print thickness
    edge_width - width of terrian border with no path cutouts
    max_print_size - maximum print dimension.  will scale and rotate print to fit this.
    water_drop - water ways printed slightly lower than other paths and terrain
    load_area - overide automatic area selection
    resolution - 30 10 or 30, for 10/30 meter resolution DEM.
    dem_offset - offset DEM relative to OSM data to account for shifts in data.
    downsample_factor, integer factor to reduce resolution of terrain surface, needed for larger models.  1mm resolution is reasonable minimum.
    map_only - stop after generating map to review (no boolean ops), recomend to do this first until all desired features look corrects oi it doesn't get hung up on boolean operations for hours.
    """


    #pull print shape from gpx file if string is input
    if isinstance(CoordPoly,str):
        gpx_file = open(CoordPoly, 'r')

        gpx = gpxpy.parse(gpx_file)
        poly=[]
        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    poly.append([point.latitude, point.longitude])
                    # print('Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation))
        CoordPoly=np.array(poly)
    else:
        CoordPoly=np.array(CoordPoly)

    if not any(load_area): #load area bounds polygon unless larger area is input
        load_area=np.vstack((np.min(CoordPoly,axis=0),np.max(CoordPoly,axis=0)))
    dem=Dem(CoordPoly,resolution,dem_offset,downsample_factor)
    # dem.plotElev2()
    #plotElev(dem)

    #generate terrain, save STL and return scale factor and x/y=0 corner of print in coords
    scale_factor,corner=tMesh(dem,CoordPoly,max_print_size,height_factor,base_height,water_drop)


    Waterbodies=getWaterbodies(CoordPoly,waterbody,corner,scale_factor,path_clearance,base_height,height_factor,dem)
    Footpaths=getFootpaths(CoordPoly,trail_exclude,corner,scale_factor,path_width,support_width,path_clearance,base_height)
    Roads=getRoads(CoordPoly,rd_include,corner,scale_factor,path_width,support_width,path_clearance,base_height)
    Waterways=getWaterways(CoordPoly,waterway_include,corner,scale_factor,path_width,support_width,path_clearance,base_height)

    p=cord2dist(xy=np.flip(CoordPoly),corner=corner,f=scale_factor)

    genPoly(p,-edge_width-path_clearance).export('temp/poly_1.stl') #path boundary.
    genPoly(p,-edge_width).export('temp/poly_2.stl') #path coutout boundary
    genPoly(p,0).export('temp/poly_3.stl') #outer edge


    plotPaths(dem,scale_factor,Waterbodies,Footpaths,Roads,Waterways,p)

    if not map_only:
        #replace with pymadcad?
        #Boolean ops using OpenSCAD, can take hours
        print('Boolean operations may take minutes/hours...')
        print('-terrain')
        os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Terrain.3mf Terrain.scad"')
        print('-waterways')
        os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Waterway.3mf Waterway.scad"')
        # print('-waterbody')
        # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Waterbody.3mf Waterbody.scad"')
        print('-trails')
        os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Trails.3mf Trails.scad"')
        print('-roads')
        os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Roads.3mf Roads.scad"')

    print("COMPLETE")
