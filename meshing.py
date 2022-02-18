# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 06:51:38 2022

@author: jkoet
"""


# from scipy.spatial import Delaunay
import numpy as np
# import triangle as tr
import shapely as shp
import trimesh as tm
# import pymeshfix as pmf 
import scipy as scipy
from cord_funcs import *
import matplotlib as mpl
import matplotlib.pyplot as plt


def borderPoints(dem,poly):
    # generate pts along border of polygon along dem grid lines
    pts=[poly] 
    for i in range(poly.shape[0]-1):
        # x crossings
        x=dem.lon[(dem.lon>=min(poly[i:i+2,0])) & (dem.lon<=max(poly[i:i+2,0]))]
        y=poly[i,1]+(x-poly[i,0])*np.diff(poly[i:i+2,1])/np.diff(poly[i:i+2,0])
        pts.append(np.vstack((x,y)).transpose())
        
        # y-crossings
        y=dem.lat[(dem.lat>=min(poly[i:i+2,1])) & (dem.lat<=max(poly[i:i+2,1]))]
        x=poly[i,0]+(y-poly[i,1])*np.diff(poly[i:i+2,0])/np.diff(poly[i:i+2,1])
        pts.append(np.vstack((x,y)).transpose())
            
    pts=np.vstack(pts)
    return pts


def borderPoints2(dem,poly):
    
    # generate pts along border of polygon along dem grid lines
    
    
    pts=[]
    for i in range(poly.shape[0]-1):
        # x crossings
        n=max(sum((dem.lon>=min(poly[i:i+2,0])) & (dem.lon<=max(poly[i:i+2,0]))),
              sum((dem.lat>=min(poly[i:i+2,1])) & (dem.lat<=max(poly[i:i+2,1]))))
        
        x=np.linspace(poly[i,0],poly[i+1,0],n+1)
        y=np.linspace(poly[i,1],poly[i+1,1],n+1)
        
        x=x[1:]
        y=y[1:]


        pts.append(np.vstack((x,y)).transpose())
            
    pts=np.vstack(pts)
    return pts

# def getMeshPermimeter

def triAltitude(t,verts):
    #calc altitude of 3s triangle, not complete.
    
    tv=verts[t,:] #trinagle verts
    
    #distance of each leg in triangle
    d=np.sort(d)
    d=np.sqrt(np.sum(np.diff(tv.take(range(0,4),axis=1,mode='wrap'),axis=1)**2,axis=2)) #lengths of 3 sides of all triangles
    s=np.sum(d,1)/2

    return p
    
    return r

def rotMesh(msh,angle):
    #rotate/scale mesh
    rotmat=np.zeros([4,4])
    rotmat[0,0]=np.cos(angle)
    rotmat[0,1]=-np.sin(angle)
    rotmat[1,0]=np.sin(angle)
    rotmat[1,1]=np.cos(angle)
    rotmat[2,2]=1
    rotmat[3,3]=1
    msh.apply_transform(rotmat)
    return msh
    
def tMesh(dem,Boundary,scale,corner,height_factor,base_height,water_drop):
    # create terrian mesh that bounds lat/lon polygon.
    
    
    
    x=dem.lon
    y=dem.lat

    x,y=cord2dist(x=x,y=y,corner=corner,f=scale)


    gridsize=np.mean([np.diff(x[0:2]),np.diff(y[0:2])])

    
    x,y=np.meshgrid(x,y)
    x=x.flatten()
    y=y.flatten()
    
    
    edge=borderPoints2(dem,dist2cord(xy=np.array(Boundary.exterior.xy).T,corner=dem.corner,f=scale))
    
    pts=shp.geometry.MultiPoint(np.vstack((x,y)).transpose())
    idx=np.zeros(len(pts),bool)
    #poly=shp.geometry.Polygon(cord2dist(Boundary,corner=corner,f=scale))
    poly=Boundary.buffer(-gridsize*.5) #don't include points too close to edge points
    for i in range(len(pts)):
        idx[i]=poly.contains(pts[i])
    x=x[idx]
    y=y[idx]
    idx=np.reshape(idx,[dem.lat.shape[0],dem.lon.shape[0]])
    

    
    z=np.hstack((dem.z[idx],dem.getElev(edge)))
    
    edge=cord2dist(edge,corner=corner,f=scale)
    
    # outer edge of mesh is always a pile of shit triangles, so make a buffer to remove after meshing.
    p=shp.geometry.Polygon(edge)
    edge_buffer = shp.affinity.scale(p,
                              xfact=(edge[:,0].max()+10*gridsize)/edge[:,0].max(),
                              yfact=(edge[:,1].max()+10*gridsize)/edge[:,1].max(),
                              origin=p.centroid)
    edge_buffer=np.array(edge_buffer.exterior.xy).transpose()
    
    z=(z-z.min())*scale*height_factor+base_height+1 #+1 provides 1mm extra thickness for min path height
    
    # tri = mpl.Triangulation(x_test, y_test)
    # ntri = tri.triangles.shape[0]
    verts=np.vstack((x,y)).transpose()
    verts=np.vstack((verts,edge,edge_buffer))
    # verts3=np.vstack((verts.transpose(),z)).transpose()
    # tri=scipy.spatial.Delaunay(verts)
    # msh=tm.creation.extrude_triangulation(tri.points, tri.vertices, 10 )
    # msh.export('dgas.stl')
    
    surf = mpl.tri.Triangulation(verts[:,0],verts[:,1])
    
    # fig1, ax1 = plt.subplots()
    # ax1.set_aspect('equal')
    # ax1.triplot(surf, 'bo-', lw=1)
    # plt.show()
    # plt.pause(0.1)
    
    vert_cutoff=verts.shape[0]-edge_buffer.shape[0]
    
    t=surf.triangles
    
    mask=np.any(t>vert_cutoff,axis=1)
    surf.set_mask(mask)
    msh=tm.creation.extrude_triangulation(verts, surf.get_masked_triangles(), 5 )
    
    #get elevations to map to top surface of mesh
    
    verts=verts[0:vert_cutoff,:]
    verts=dist2cord(verts,corner=corner,f=scale)
    
    
    idx=msh.vertices[:,2]>0
    z=dem.getElev(dist2cord(msh.vertices[idx,0:2],corner=corner,f=scale))
    z=(z-dem.z.min())*scale*height_factor+base_height+1
    msh.vertices[idx,2]=z
    t=surf.get_masked_triangles()

    
    tm.repair.fix_normals(msh)
    return msh

        
def includedPoints(poly,dem,corner,scale_factor):
    grid=np.meshgrid(dem.lon,dem.lat)
    pts=np.vstack((grid[0].flatten(),grid[1].flatten())).transpose()
    pts=cord2dist(xy=pts,corner=corner,f=scale_factor)
    idx=(pts[:,0]>poly.bounds[0]) & (pts[:,0]<poly.bounds[2]) & (pts[:,1]>poly.bounds[1]) & (pts[:,1]<poly.bounds[3])
    pts=shp.geometry.MultiPoint(pts[idx,:])
    c=np.vstack(poly.exterior.xy).transpose();
        
    idx2=np.zeros(len(pts),bool)
    for j in range(len(pts)):
        idx2[j]=poly.contains(pts[j])
    # c=dist2cord(c,corner=corner,f=scale_factor)
    
    idx[idx]=idx2
    
    idx=np.reshape(idx,[dem.lat.shape[0],dem.lon.shape[0]])
    return idx

def bufferSmooth(ply,offset):
    if not isinstance(offset, list):
        offset=[offset]
    ply=ply.buffer(offset[0])
    if len(offset)==1:
        ply=ply.buffer(-offset[0])
    else:
        ply=ply.buffer(-offset[0]-offset[1])
        ply=ply.buffer(+offset[1])
    ply=ply.simplify(0.05) #, preserve_topology=False
    return ply


def meshgen_wb(ply,h,fname,elev):
    msh=[]
    if not isinstance(ply, list):
        ply=[ply]
    msh=[]
    for i in range(len(ply)):
        if ply[i].geom_type == 'MultiPolygon':
            msh.append([])
            for j in range(len(ply[i])):
                if ply[i][j].area>1: #don't use areas less than 0.5 mm^2
                    msh[i].append(tm.creation.extrude_polygon(ply[i][j],h[i]))
                    msh[i][-1].apply_translation([0,0,elev[j]])
            msh[i]=tm.boolean.union(msh[i])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
        else:
            msh.append(tm.creation.extrude_polygon(ply[i], h[i]))
            msh[i].apply_translation([0,0,elev[0]])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
    return msh

                
def meshgen(ply,h,fname,elev=0):

    if not isinstance(ply, list):
        ply=[ply]

    msh=[]
    for i in range(len(ply)):
        if ply[i].geom_type == 'MultiPolygon':
            msh.append([])
            for j in range(len(ply[i])):
                if ply[i][j].area>0.5: #don't use areas less than 0.5 mm^2
                    msh[i].append(tm.creation.extrude_polygon(ply[i][j],h))
                else:
                    stop=1
            msh[i]=tm.boolean.union(msh[i])
            msh[i].apply_translation([0,0,elev])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
        else:
            msh.append(tm.creation.extrude_polygon(ply[i], h))
            msh[i].apply_translation([0,0,elev])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
    return msh

def meshgen2(ply,h,dem,corner,sf,hf,base,fname=[]):
    if not isinstance(ply, list):
        ply=[ply]
    msh=[]
    for i in range(len(ply)):
        if ply[i].geom_type == 'MultiPolygon':
            msh.append([])
            for j in range(len(ply[i])):
                if ply[i][j].area>0.5: #don't use areas less than 0.5 mm^2
                    c=np.vstack(ply[i][j].exterior.xy).transpose();
                    c=dist2cord(c,corner=corner,f=sf)
                    z=(min(dem.getElev(c))-np.min(dem.z))*sf*hf+base+1-2.5
                    msh[i].append(tm.creation.extrude_polygon(ply[i][j],h))
                    msh[i][-1].apply_translation([0,0,z])
            msh[i]=tm.boolean.union(msh[i])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
        else:
            c=np.vstack(ply[i].exterior.xy).transpose();
            c=dist2cord(c,corner=corner,f=sf)
            z=(min(dem.getElev(c))-np.min(dem.z))*sf*hf+base+1-2.5
            msh.append(tm.creation.extrude_polygon(ply[i], h))
            msh[i].apply_translation([0,0,z])
            if len(fname)>0:
                msh[i].export(fname + '_' +str(i+1) +'.stl')
    return msh
    
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
                    buffer_poly=ply[j].buffer(o)
                    if buffer_poly.geom_type == 'MultiPolygon':
                        for p in buffer_poly:
                            mp.append(p)
                    else:
                        mp.append(buffer_poly)
                ply2.append(shp.geometry.MultiPolygon(mp))
            else:
                ply2.append(ply.buffer(o))
    return ply2 

def binaryOps2(Waterbodies,Footpaths,Roads,Waterways,b,pWidth,sWidth,clearance):
    
    # 0-cutout size
    # 1-top size
    # 2-support size

    #this cold probably done with loop and heiarchy

    if len(Roads)>0:
        Roads=offsetLines(Roads,(pWidth)/2)[0]
        Roads=Roads.intersection(b[0])
        Roads=bufferSmooth(Roads,[.5,-.2])
        # Roads[0]=Roads[0].simplify(pWidth/8)

        Roads=offsetPoly(Roads,[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout

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
        Footpaths[1]=bufferSmooth(Footpaths[1],0.01)#fixed self-intersection.
        Footpaths[2]=bufferSmooth(Footpaths[2],0.01)#fixed self-intersection.
        cutout=Roads[0].union(Footpaths[0])
    else: #no roads
        Footpaths=offsetPoly(Footpaths,[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout
        cutout=Footpaths[0]
        

    if len(Waterbodies)>0:
        #already a polygon
        Waterbodies=bufferSmooth(Waterbodies,0.01)#fixed self-intersection.
        Waterbodies=Waterbodies.intersection(b[0])
        
      
        if len(Roads)>0:  
            Waterbodies=Waterbodies.difference(Roads[1])
        Waterbodies=Waterbodies.difference(Footpaths[1])
        Waterbodies=bufferSmooth(Waterbodies,[-.25,+.15])

        Waterbodies=offsetPoly(Waterbodies,[clearance/2,0,-(pWidth-sWidth)/2])# offset for cutout, no support needed.
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
        Waterways=bufferSmooth(Waterways,-0.2) 
        
        Waterways=offsetPoly(Waterways,[clearance/2,0,-(pWidth-sWidth)/2])#add support & cutout
        Waterways[0]=bufferSmooth(Waterways[0],0.01)#fixed self-intersection.

        
        cutout=cutout.union(Waterways[0])
        cutout=bufferSmooth(cutout,.1)
        #Waterways.extend(offsetPoly(Waterways[1],(sWidth-pWidth)/2))# offset down for support


    return Waterbodies,Footpaths,Roads,Waterways,cutout