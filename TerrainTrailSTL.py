# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:50:30 2021

@author: Jeremy Koether
"""


import time as time
# import inflator
import os as os
from copy import deepcopy
import glob as glob

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
import shapely as shp
import descartes as desc



from zipfile import ZipFile

import matplotlib.pyplot as plt

#TODO
#   some routes can go over a taller peak, these should be cropped.
#   accept two coordinate points for a rectangle
#   Routes could be terminated when they get too far away from the required low elevation.

api = overpy.Overpass()

def getDem(importBounds,res):
    #Load USGS DEM data from 1/3 arcsecond (10m) data files
    elevFiles=[]
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
        if res==10:
            filename='USGS_13_'+f
        elif res==30:
            filename='USGS_1_'+f
        else:
            raise NameError  ('resolution must be 10 or 30')
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


def seqEdges(n):
    edg=np.arange(n-1)
    edg=np.transpose(np.vstack((edg,edg+1)))
    edg=np.vstack((edg,[n-1,0]))
    return edg



def getFootpaths(poly,trail_exclude,trail_gpx,corner,scale_factor,offsets,base):
    print(' ')

    if isinstance(trail_gpx,str):
        print('Loading footpath from gpx file...')

        gpx_file = open(trail_gpx, 'r')
        
        gpx = gpxpy.parse(gpx_file)
        poly=[]
        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    poly.append([point.latitude, point.longitude])
                    # print('Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation))
        verts=np.flip(np.array(poly),axis=1)
        #TODO smoothing Kalman filter
        lines = shp.geometry.MultiLineString([cord2dist(xy=verts,corner=corner,f=scale_factor)])
    else:
        print('Retrieving OSM footpaths...')

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
        coords=[]
        for w in result.ways:
            if w.tags['highway'] in ['path','footway','cycleway'] and not w.id in trail_exclude:
                inc_ways.append(w)
                c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
                c=cord2dist(xy=c,corner=corner,f=scale_factor)
                coords.append(c)

        lines = shp.geometry.MultiLineString(coords)
    # p=offsetLines(lines,offsets)
    # return p
    return lines

def getRoads(poly,roads,corner,scale_factor,offsets,base):
    if len(roads)==0:
        return []
    print(' ')
    print('Retrieving OSM roads...')

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
    coords=[]
    for w in result.ways:
        if w.id in roads or ("name" in w.tags and w.tags["name"] in rd_names):
            inc_ways.append(w)
            c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
            c=cord2dist(xy=c,corner=corner,f=scale_factor)
            coords.append(c)
            lines = shp.geometry.MultiLineString(coords)
    # if any(inc_ways):
    #     return offsetLines(lines,offsets)
    # else:
    #     return []
    if any(inc_ways):
        return lines
    else:
        return []



def getWaterways(poly,waterways,corner,scale_factor,offsets,base,map_only):
    if len(waterways)==0:
        return []
    print(' ')
    print('Retrieving OSM waterways...')

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
    coords=[]
    for w in result.ways:
        if  w.id in waterways or ("name" in w.tags and w.tags["name"] in waterways):
            inc_ways.append(w)
            c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
            c=cord2dist(xy=c,corner=corner,f=scale_factor)
            coords.append(c)
            lines = shp.geometry.MultiLineString(coords)
    # if any(inc_ways):
    #     return offsetLines(lines,offsets)
    # else:
    #     return []
    if any(inc_ways):
        return lines
    else:
        return []
    
def flip(x, y,z):
    """Flips the x and y coordinate values"""
    return y, x, z


def findConnLine(l1,lines):
    endPoint=np.array(l1.xy)[:,-1]
    for l in lines:
        xy=np.array(l.xy)
        if np.all(abs(xy[:,0]-endPoint)<0.1):
            l2=lines.pop(0)
            c1=np.array([tuple(coord) for coord in list(l1.coords)])
            c2=np.array([tuple(coord) for coord in list(l2.coords)])
            l1=shp.geometry.LineString(np.vstack((c1,c2[1:,:])))
            return l1,lines,True
        if np.all(abs(xy[:,-1]-endPoint)<0.1):
            l2=lines.pop(0)
            c1=np.array([tuple(coord) for coord in list(l1.coords)])
            c2=np.array([tuple(coord) for coord in list(l2.coords)])
            l1=shp.geometry.LineString(np.vstack((c1,c2[:-1,:])))
            return l1,lines,True
    return l1,lines,False
        

def getWaterbodies(poly,bodies,corner,scale_factor,clearance,base,height_factor,dem):
    if len(bodies)==0:
        return [],[]
    Thickness=3
    
    print(' ')
    print('Retrieving OSM waterbodies...')
    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,0]),np.min(poly[:,1]),np.max(poly[:,0]),np.max(poly[:,1])))
    wb=[]


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
            wb.append(ways)
            
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
            wb.append(w)
    polys=[]
    polys2=[]
    
    for w in wb:
        if isinstance(w, list):
            lines=[]
            for w2 in w:
                c=np.array([[float(n.lon) for n in w2.nodes],[float(n.lat) for n in w2.nodes]]).transpose()
                z=dem.getElev(c)
                c=cord2dist(xy=c,corner=corner,f=scale_factor)
                c=np.vstack((c.transpose(),z)).transpose()
                lines.append(shp.geometry.LineString(c))
            # c=np.vstack(c)
            # z=dem.getElev(c) #perimeter elevation profile
            # z=(z-np.min(dem.z))*scale_factor*height_factor+base
            # z=np.sort(z)[int(0.25*z.shape[0])] #use the 25th percentile value as water level
            # elev.append(z+base)
            # c=cord2dist(xy=c,corner=corner,f=scale_factor)
            # lines=shp.ops.linemerge(lines)
            if len(lines)==1:
                polys.append(shp.geometry.Polygon(lines[0]))
            else:
                merged_line=lines.pop(0)
                added=True
                while len(lines)>0:
                    merged_line,lines,added=findConnLine(merged_line,lines)
                    if not added:
                        polys.append(shp.geometry.Polygon(merged_line))
                        merged_line=lines.pop(0)
                polys.append(shp.geometry.Polygon(merged_line))

        else:
            c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
            z=dem.getElev(c) #perimeter elevation profile
            
            c=cord2dist(xy=c,corner=corner,f=scale_factor)
            c=np.vstack((c.transpose(),z.transpose())).transpose() #why doesn't hstack work here?!
            polys.append(shp.geometry.Polygon(c))

    elev=[]
    for p in polys:
        c=np.array([tuple(coord) for coord in list(p.exterior.coords)])
        z=c[:,2]
        # z=(z-np.min(dem.z))*scale_factor*height_factor+base
        elev.append(np.sort(z)[int(0.25*z.shape[0])]) #use the 25th percentile value as water level)
    polys=shp.geometry.MultiPolygon(polys)
    
    #add flat spot to dem corresponding to points within the warter bodies.
    #This helps the lake edges, and is needed if any of the edges cut through a lake.
    grid=np.meshgrid(dem.lat,dem.lon)
    x,y=cord2dist(x=grid[1].flatten(),y=grid[0].flatten(),corner=corner,f=scale_factor)
    pts=shp.geometry.MultiPoint(np.array([x,y]).transpose())
    for i in range(len(polys)):
        #find points within each waterbody.
        idx=np.zeros(len(pts),bool)
        for j in range(len(pts)):
            idx[i]=polys[i].contains(pts[j])
        idx=np.reshape(idx,[dem.lon.shape[0],dem.lat.shape[0]]).transpose()
        dem.z[idx]=elev[i] #set points to nominal elevation.
    #find all points within each lake
    return polys,elev,dem

def printScaling(dem,edge_poly,max_Size):
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

    scale=scale.max() #use angle that allows the largest print.  model is not rotated, optimal angle is output and part will need rotated in slicer.
    
    x,y=cord2dist(x=x,y=y,corner=corner)
    print('DEM resolution: {0:.2f}x{1:.2f} mm'.format(((x.max()-x.min())/(x.shape[0]-1)*scale),(y.max()-y.min())/(y.shape[0]-1)*scale))

    return scale,corner

def tMesh(dem,scale,corner,height_factor,base_height,water_drop):
    # #create terrian mesh that bounds lat/lon polygon.  polygon shape and all path
    # #features will be cut out using boolean operations
    x=dem.lon
    y=dem.lat


    # corner=np.stack((np.min(x),np.min(y)))

    # p=cord2dist(np.flip(deepcopy(edge_poly),1),corner=corner)
    # d=np.sqrt(np.sum((p[0,:]-p)**2,axis=1))
    # a=np.arctan2(p[:,1]-p[0,1],p[:,0]-p[0,0])
    # scale=np.zeros(91)
    # for i in range(91):
    #     #Rotate points in one degee increments, and determine the min scale factor between the two directions
    #     p=np.vstack((d*np.cos(a+i*np.pi/180),d*np.sin(a+i*np.pi/180))).T

    #     scale[i]=min([max(max_Size)/(p[:,0].max()-p[:,0].min()),
    #                   min(max_Size)/(p[:,1].max()-p[:,1].min())]) #scale factor for each angle
    # print('print angle: ' +str(np.argmax(scale)))

    # scale=scale.max() #find angle that allows the largest print.

    x,y=cord2dist(x=x,y=y,corner=corner)


    #convert meters on map to mm on print.
    x=x*scale
    y=y*scale

    
    
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


    return scale,corner,msh

def drawPoly(ax,ply,c):
    if ply.geom_type == 'MultiPolygon':
        for p in ply:
            patch = desc.PolygonPatch(p,fc=c,linewidth=0)
            ax.add_patch(patch) 
    else:
        patch = desc.PolygonPatch(ply,fc=c,linewidth=0)
        ax.add_patch(patch) 


def plotPaths(dem,sf,wb,fp,rp,ww,poly,cut):


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

    for i in range(poly.shape[0]):
        ax.plot(np.hstack((poly[:,0],poly[0,0])),np.hstack((poly[:,1],poly[0,1])),'r:')

    ax.plot(poly[0,0],poly[0,1],'kx')
    
    
    if len(rp)>0:
        drawPoly(ax,rp[1],'black')
    if len(fp)>0:
        drawPoly(ax,fp[1],'green')

    if len(wb)>0:
        drawPoly(ax,wb[1],'teal')
    if len(ww)>0:
        drawPoly(ax,ww[1],'blue')
        
    if cut.geom_type == 'MultiPolygon':
        for p in cut:
            patch = desc.PolygonPatch(p,fc='none',linewidth=0.5)
            ax.add_patch(patch) 
    else:
        patch = desc.PolygonPatch(cut,fc='none',linewidth=0.5)
        ax.add_patch(patch) 
    #plt.title(title)
    plt.axis('equal')
    plt.show()
    plt.pause(0.1)


def offsetLines(lines,offsets):
    if not isinstance(offsets, list):
        offsets=[offsets]
    ply=[]
    for o in offsets:
        ply.append(lines.buffer(o))
    return ply
    
def offsetPoly(ply,offsets):
    if not isinstance(offsets, list):
        offsets=[offsets]
    ply2=[]
    for o in offsets:
        if o==0:
            ply2.append(ply)
        else:
            if ply.geom_type == 'MultiPolygon':
                mp=[]
                for j in range(len(ply)):
                    mp.append(ply[j].buffer(o))
                ply2.append(shp.geometry.MultiPolygon(mp))
            else:
                ply2.append(ply.buffer(o))
    return ply2

def genPoly(p,offset):

    pco = pc.PyclipperOffset()
    pco.AddPath(p*1000, 0,pc.ET_CLOSEDPOLYGON)

    p=np.vstack(pco.Execute(offset*1000)).astype('double')/1000

    face = {
      "vertices": p,
      "segments": seqEdges(p.shape[0])}

    t = tr.triangulate(face,'p')

    mesh=tm.creation.extrude_triangulation(t['vertices'], t['triangles'], 200)
    # translation = [0,0,2]  # box offset + plane offset
    # mesh.apply_translation(translation)
    return mesh

def bufferSmooth(ply,offset):
    if not isinstance(offset, list):
        offset=[offset]
    # if ply.geom_type == 'MultiPolygon':
    #     b=[]
    #     for j in range(len(ply)):
    #         b.append=ply[j].buffer(offset[0])
    #         b[j]=ply[j].buffer(-offset[0]-offset[1])
    #         b[j]=ply[j].buffer(offset[1])
    #     ply=
    # else:
    ply=ply.buffer(offset[0])
    if len(offset)==1:
        ply=ply.buffer(-offset[0])
    else:
        ply=ply.buffer(-offset[0]+offset[1])
        ply=ply.buffer(-offset[1])
    ply=ply.simplify(0.05) #, preserve_topology=False
    return ply

def combineCutouts(coutout,waterbodies,base,h,elev):

    msh=meshgen(coutout,base,h,'',elev=0)
    msh.apply_translation([0,0,base])
    # if not isinstance(elev, list):
    #     elev=[elev]*len(ply)
    if len(waterbodies)>0:
        for i in range(len(waterbodies)):
            msh2=tm.creation.extrude_polygon(waterbodies[i])
            msh2.apply_translation([0,0,base+elev[i]-3])
            msh=tm.boolean.union([msh,msh2])
    else:
        msh2=tm.creation.extrude_polygon(waterbodies, h)
        msh2.apply_translation([0,0,base+elev[i]-3])
        msh=tm.boolean.union([msh,msh2])
    msh.export('temp/cutout_1.stl')

            

def meshgen(ply,h,fname,elev=0):
    msh=[]
    if not isinstance(ply, list):
        ply=[ply]
    if not isinstance(h, list):
        h=[h]*len(ply)
    if not isinstance(elev, list):
        elev=[elev]*len(ply)
    if isinstance(ply, list):
        for i in range(len(ply)):
            if ply[i].geom_type == 'MultiPolygon':
                msh=[]
                for j in range(len(ply[i])):
                    msh.append(tm.creation.extrude_polygon(ply[i][j],h[i]))
                msh=tm.boolean.union(msh)
                msh.apply_translation([0,0,elev[i]])
                if len(fname)>0:
                    msh.export(fname + '_' +str(i+1) +'.stl')
            else:
                msh=tm.creation.extrude_polygon(ply[i], h[i])
                msh.apply_translation([0,0,elev[i]])
                if len(fname)>0:
                    msh.export(fname + '_' +str(i+1) +'.stl')
    else:
        print('xxx')
        # if ply.geom_type == 'MultiPolygon':
        #     msh=[]
        #     for j in range(len(ply)):
        #         msh.append(tm.creation.extrude_polygon(ply[j], h))
        #     msh=tm.boolean.union(msh)
        #     msh.apply_translation([0,0,base+elev[0]])
        #     if len(fname)>0:
        #         msh.export(fname + '.stl')
    return msh
    
 

def binaryOps2(Waterbodies,Footpaths,Roads,Waterways,b,pWidth,sWidth,clearance):
    
    # 0-cutout size
    # 1-top size
    # 2-support size

    

    if len(Roads)>0:
        Roads=offsetLines(Roads,(pWidth+clearance)/2)
        Roads[0]=Roads[0].intersection(b[0])
        Roads[0]=bufferSmooth(Roads[0],[.5,-.2])
        # Roads[0]=Roads[0].simplify(pWidth/8)

        Roads.extend(offsetPoly(Roads[0],-clearance/2))# offset down for top
        Roads.extend(offsetPoly(Roads[1],(sWidth-pWidth)/2))# offset down for support
        #nothing to cut out of roads.
    
    #assume there is always a footpath
    Footpaths=offsetLines(Footpaths,(pWidth/2))[0]
    Footpaths=Footpaths.intersection(b[0])

    
    if len(Roads)>0:       
        Footpaths=Footpaths.difference(Roads[1])
        Footpaths=bufferSmooth(Footpaths,[+.4,-.15])
        Footpaths=Footpaths.difference(Roads[1]) #needs repeated to prevent bridging/overlap
        Footpaths=offsetPoly(Footpaths,[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout
        Footpaths[0]=bufferSmooth(Footpaths[0],0.01)#fixed self-intersection.
        Footpaths[2]=bufferSmooth(Footpaths[2],0.01)#fixed self-intersection.
        cutout=Roads[0].union(Footpaths[0])
    else: #no roads
        Footpaths=offsetPoly(Footpaths[0],[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout
        cutout=Footpaths[0]
        

    if len(Waterbodies)>0:
        #already a polygon
        Waterbodies=bufferSmooth(Waterbodies,0.01)#fixed self-intersection.
        Waterbodies=Waterbodies.intersection(b[0])
        
      
        
        Waterbodies=Waterbodies.difference(Roads[1])
        Waterbodies=Waterbodies.difference(Footpaths[1])
        Waterbodies=bufferSmooth(Waterbodies,[+.4,-.15])

        Waterbodies=offsetPoly(Waterbodies,[clearance/2,0])# offset for cutout, no support needed.
        cutout=cutout.union(Waterbodies[0])
        
    if len(Waterways)>0:
        Waterways=offsetLines(Waterways,(pWidth/2))[0]
        Waterways=Waterways.intersection(b[0])
        if len(Roads)>0:    
            Waterways=Waterways.difference(Roads[1])
        Waterways=Waterways.difference(Footpaths[1])
        if len(Waterbodies)>0:
            Waterways=Waterways.difference(Waterbodies[1])
            
        Waterways=bufferSmooth(Waterways,[+.4,-.15])
        if len(Roads)>0:    
            Waterways=Waterways.difference(Roads[1]) #needs repeated to prevent bridging/overlap
        Waterways=Waterways.difference(Footpaths[1])
        if len(Waterbodies)>0:
            Waterways=Waterways.difference(Waterbodies[1])
        
        Waterways=offsetPoly(Waterways,[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout
        Waterways[0]=bufferSmooth(Waterways[0],0.01)#fixed self-intersection.

        
        cutout=cutout.union(Waterways[0])
        cutout=bufferSmooth(cutout,.1)
        Waterways.extend(offsetPoly(Waterways[1],(sWidth-pWidth)/2))# offset down for support

  
    

    return Waterbodies,Footpaths,Roads,Waterways,cutout


def GenerateSTLs(Boundary,rd_include=[],trail_exclude=[],waterway_include=[],waterbody=[],trail_gpx=[],
                 path_width=.7,support_width=0.9,path_clearance=0.1,height_factor=2,base_height=4,top_thickness=1.5,
                 edge_width=1.5,max_print_size=[248,198],water_drop=0.5,load_area=[],resolution=30,
                 dem_offset=[0,0],downsample_factor=1,map_only=False):
    """

    INPUTS:

    Boundary - polygon to define shape of terrain, input as longitude/lattitude array, or gpx file of points
    rd_include - roads names / ids to inlcude, no roads included by default.
    trail_exclude - footpaths names / ids to inlcude, all included by default.
    waterway_include - water paths to inlcude, not inlcuding polygon water bodies
    waterway_include - water areas to include (lakes), current set to print at one level, just below surrounding land.  Will not work well for rivers.
    trail_gpx - gpx file to use for trail path
    path_width - path (road, footpaths and waterway) print top width, 3 total passes works best.
    support_width - width of base (one pass less than top)
    path_clearance - clearance between path prints and cutouts.
    height_factor - exaggeration factor for elevation. 1.0 makes on same scale as horizontal dimensions.
    base_height - minium print thickness
    edge_width - width of terrian border with no path cutouts
    max_print_size - maximum print dimension.  will scale and rotate print to fit this.
    water_drop - water ways printed slightly lower than other paths and terrain
    load_area - overide automatic area selection
    top_thickness - vertical thickness of wider top section of roads, trails and waterways.
    resolution - 30 10 or 30, for 10/30 meter resolution DEM.
    dem_offset - offset DEM relative to OSM data to account for shifts in data.
    downsample_factor, integer factor to reduce resolution of terrain surface, needed for larger models.  1mm resolution is reasonable minimum.
    map_only - stop after generating map to review (no boolean ops), recomend to do this first until all desired features look corrects oi it doesn't get hung up on boolean operations for hours.
    """
    if not os.path.isdir('print_files/'):
        os.mkdir('print_files/')
    if not os.path.isdir('temp/'):
        os.mkdir('temp/')
    # clear temp files
    fileList = glob.glob('temp/*.stl')
    # Iterate over the list of filepaths & remove each file.
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)
            
            
    # clear print files
    fileList = glob.glob('print_files/*.stl')
    # Iterate over the list of filepaths & remove each file.
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)
    
    #pull print shape from gpx file if string is input
    if isinstance(Boundary,str):
        gpx_file = open(Boundary, 'r')

        gpx = gpxpy.parse(gpx_file)
        poly=[]
        for track in gpx.tracks:
            for segment in track.segments:
                for point in segment.points:
                    poly.append([point.latitude, point.longitude])
                    # print('Point at ({0},{1}) -> {2}'.format(point.latitude, point.longitude, point.elevation))
        Boundary=np.array(poly)
    else:
        Boundary=np.array(Boundary)

    if not any(load_area): #load area bounds polygon unless larger area is input
        load_area=np.vstack((np.min(Boundary,axis=0),np.max(Boundary,axis=0)))
        
   
    dem=Dem(Boundary,resolution,dem_offset,downsample_factor)
    # dem.plotElev2()
    #plotElev(dem)

    #find the largest print that fits (rotated) within print size a return scale factor and x/y=0 corner of print in coords'
    scale_factor,corner=printScaling(dem,Boundary,max_print_size,)


    # offsets=[support_width/2,path_width/2,(path_width+path_clearance)/2]
    offsets=[(path_width+path_clearance)/2]
    
    Footpaths=getFootpaths(Boundary,trail_exclude,trail_gpx,corner,scale_factor,offsets,base_height)
    
    Roads=getRoads(Boundary,rd_include,corner,scale_factor,offsets,base_height)
    
    Waterways=getWaterways(Boundary,waterway_include,corner,scale_factor,offsets,base_height,map_only)
    Waterbodies,wb_heights,dem=getWaterbodies(Boundary,waterbody,corner,scale_factor,path_clearance,base_height,height_factor,dem)
    
    
    
    p=cord2dist(xy=np.flip(Boundary),corner=corner,f=scale_factor)
    borders=offsetLines(shp.geometry.Polygon(p),[-edge_width-path_clearance,-edge_width-path_clearance/2,-edge_width,0])
    
    print('2D Boolean Operations')
    Waterbodies,Footpaths,Roads,Waterways,Cutout=binaryOps2(Waterbodies,Footpaths,Roads,Waterways,borders,path_width,support_width,path_clearance)
    if not map_only:
        tMesh(dem,scale_factor,corner,height_factor,base_height,water_drop)
        print('Meshing')
        meshgen(Roads[1:],100,'temp/road',elev=base_height)
        meshgen(Footpaths[1:],100,'temp/trail',elev=base_height)
        meshgen(Waterways[1:],100,'temp/water',elev=base_height)
        if len(Waterbodies)>0:
            meshgen(Waterbodies[1],wb_heights,'print_files/waterbodies')
        meshgen(Cutout,100,'temp/cutout',elev=base_height)
        meshgen(borders[3],100,'temp/border',elev=-1)#shift down 1mm to avoid coplaner faces]
        
        #combineCutouts(Cutout,Waterbodies[0],base_height,100,wb_heights)
    plotPaths(dem,scale_factor,Waterbodies,Footpaths,Roads,Waterways,p,Cutout)



    if not map_only:
        print('3D Boolean Mesh Operations')

        blender_loc='C:\\Program Files\\Blender Foundation\\Blender 3.0\\blender.exe'

        #replace with pymadcad?
        #Boolean ops using OpenSCAD, can take hours
        t1_start = time.perf_counter()
        os.system('"' + blender_loc + '"' +' --background --python genPrints.py -- '+str(top_thickness)+' '+str(water_drop))
        print(str(time.perf_counter()-t1_start))
        # print('-terrain')
        # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Terrain.3mf Terrain.scad"')
        # print('-waterways')
        # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Waterway.3mf Waterway.scad"')
        # # print('-waterbody')
        # # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Waterbody.3mf Waterbody.scad"')
        # print('-trails')
        # t1_start = time.perf_counter()
        # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Trails.3mf Trails.scad"')
        # print("0.1 : " + str(time.perf_counter()-t1_start))
        # print('-roads')
        # os.system('cmd /c "C:\\"Program Files\OpenSCAD\openscad.com" -o print_files\Roads.3mf Roads.scad"')

    print("COMPLETE")
