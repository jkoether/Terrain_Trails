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
import numpy as np
import trimesh as tm
# import pyclipper as pc


import gpxpy as gpxpy


from scipy import ndimage
from subprocess import run
import shapely as shp
import descartes as desc

import matplotlib.pyplot as plt


from dem_funcs import Dem 
from osm_queries import *
from meshing import *
from cord_funcs import *




def printScaling(dem,Boundary,print_size):
    #determine largest size trhat can be printed with a single print
    x=dem.lon
    y=dem.lat


    corner=np.stack((np.min(x),np.min(y)))

    p=cord2dist(Boundary,corner=corner)

    print_angle=np.linspace(0,np.pi/2,181) 
    scale=np.zeros(print_angle.shape[0])
    
    for i in range(print_angle.shape[0]):
        #rotated poly:
        rp=shp.affinity.rotate(shp.geometry.Polygon(p), print_angle[i],use_radians=True)
        scale[i]=min([print_size[0]/(rp.bounds[2]-rp.bounds[0]),print_size[1]/(rp.bounds[3]-rp.bounds[1])])


    print('print angle: {0:.2f} deg'.format(print_angle[np.argmax(scale)]*180/np.pi))

    scale=scale.max() #use angle that allows the largest print.  model is not rotated, optimal angle is output and part will need rotated in slicer.
    
    x,y=cord2dist(x=x,y=y,corner=corner)
    print('DEM resolution: {0:.2f}x{1:.2f} mm'.format(((x.max()-x.min())/(x.shape[0]-1)*scale),(y.max()-y.min())/(y.shape[0]-1)*scale))
    edge_poly=shp.geometry.Polygon(p*scale)
    dem.scale_factor=scale
    return scale,corner,edge_poly

def printScaling_tiled(dem,Boundary,print_size,tiles):
    #determine largest size that can be printed with the given number of tiles (only 2 supported for now.)
    dovetail_gap=-.1 #negative for interference
    dovetail_height=7
    dovetail_spacing=100 #set to zero for just one per edge
    
    cutout=tm.load('dovetail_cutout.stl')
    t=np.eye(4)
    t[2, 2] *= (dovetail_height+2)/10
    cutout.apply_transform(t)
    
    

    x=dem.lon
    y=dem.lat
    
    corner=np.stack((np.min(x),np.min(y)))

    p=cord2dist(Boundary,corner=corner)

    print_angle=np.linspace(0,np.pi,181) 
    scale=np.zeros(print_angle.shape[0])
    corners=np.zeros([2,2,print_angle.shape[0]])
    
    #find the rotation that allows for the largest scale print
    for i in range(print_angle.shape[0]):
        #rotated boundary poly:
        rp=shp.affinity.rotate(shp.geometry.Polygon(p), print_angle[i],use_radians=True,origin=[0,0])
        scale[i]=print_size[1]*2/(rp.bounds[3]-rp.bounds[1])
        #split on line parrallel to x-axis, find scaling assuming y-axis length is limiting 
        mid=(rp.bounds[3]+rp.bounds[1])/2
        split_line=shp.geometry.LineString([[rp.bounds[0]-1,mid],[rp.bounds[2]+1,mid]])
        rp=shp.ops.split(rp, split_line)
        r=np.array([p.bounds for p in rp])
        idx=np.argsort(r[:,1])
        r=r[idx.tolist()]
        max_width=(r[:,2]-r[:,0]).max()
        
        #check if print is to large in x-direction, allowing for offsets between rows, reduce scale if needed.
        if max_width*scale[i]>print_size[0]: 
            scale[i]=print_size[0]/max_width 

        corners[0,0,i]=r[0,0]
        corners[0,1,i]=r[1,1]-print_size[1]/scale[i]
        corners[1,:,i]=r[1,0:2]
    idx=np.argmax(scale)
    print_angle=print_angle[idx]
    
    print('print angle: {0:.2f} deg'.format(print_angle*180/np.pi))
    
    scale=scale[idx] #use angle that allows the largest print.  model is not rotated, optimal angle is output and part will need rotated in slicer.
    corners=corners[:,:,idx]*scale
    

    
    #rotated boundary polygon, used to place inserts
    rp=shp.affinity.rotate(shp.geometry.Polygon(p*scale), print_angle,use_radians=True,origin=[0,0])
    
    
    edge_poly=[]
    pts=np.array([[0, 0],[0,print_size[1]],[print_size[0],print_size[1]],[print_size[0],0]])
    

    cutouts=[]
    
    #layout tile edges and inserts in rotated space.
    for i in range(corners.shape[0]):
        edge_poly.append(shp.geometry.Polygon(corners[i,:]+pts))
        
        #build complete collection of insert cutouts that can be subtracted from each tile to make required cutouts.
        if i>0:
            split_line=shp.geometry.LineString([[rp.bounds[0]-1,corners[i,1]],[rp.bounds[2]+1,corners[i,1]]])
            sp=shp.ops.split(rp, split_line)[0]
            seam=np.vstack(sp.exterior.xy)
            seam=seam[(0,abs(seam[1,:]-corners[i,1])<0.5)] #find points close to seam between tiles
            
            center=(seam.max()+seam.min())/2
            #-2 in Z here to prevent colinear faces during boolean ops.
            cutouts_loc=[[-dovetail_spacing/2+center,corners[i,1],-2],[dovetail_spacing/2+center,corners[i,1],-2]]
            for c in cutouts_loc:
                new_cutout=deepcopy(cutout)
                new_cutout.apply_translation(c)
                cutouts.append(new_cutout)
                
    cutouts=tm.boolean.union(cutouts)
    cutouts=rotMesh(cutouts,-print_angle)

    
    edge_poly=shp.geometry.MultiPolygon(edge_poly)
    edge_poly=shp.affinity.rotate(edge_poly, -print_angle,use_radians=True,origin=[0,0])
    x,y=cord2dist(x=x,y=y,corner=corner)
    print('DEM resolution: {0:.2f}x{1:.2f} mm'.format(((x.max()-x.min())/(x.shape[0]-1)*scale),(y.max()-y.min())/(y.shape[0]-1)*scale))
    dem.scale_factor=scale
    return scale,corner,edge_poly,[cutouts]

def drawPoly(ax,ply,c):
    if ply.geom_type == 'MultiPolygon':
        for p in ply:
            patch = desc.PolygonPatch(p,fc=c,linewidth=0)
            ax.add_patch(patch) 
    else:
        patch = desc.PolygonPatch(ply,fc=c,linewidth=0)
        ax.add_patch(patch) 

def plotPaths(dem,sf,wb,fp,rp,ww,border,cut,cmp,e_poly):


    #Create map figure to show terrain an slected paths.
    fig4, ax = plt.subplots()
    x=dem.lon
    y=dem.lat
    corner=np.stack((np.min(x),np.min(y)))
    x,y=cord2dist(x=x,y=y,corner=corner,f=sf)

    grad=ndimage.sobel(dem.z)
    plt.imshow(grad, cmap ='pink',
                 extent =[x.min(), x.max(), y.min(), y.max()],
                    origin ='lower')

    # for i in range(border.shape[0]):
    #     ax.plot(np.hstack((border[:,0],border[0,0])),np.hstack((border[:,1],border[0,1])),'r:')

    # ax.plot(border[0,0],border[0,1],'kx')
    patch = desc.PolygonPatch(border,fc='none',ec='red',linewidth=0.5)
    ax.add_patch(patch) 
    
    if len(rp)>0:
        drawPoly(ax,rp[1],'black')
    if len(fp)>0:
        drawPoly(ax,fp[1],'green')

    if len(wb)>0:
        drawPoly(ax,wb[1],'teal')
    if len(ww)>0:
        drawPoly(ax,ww[1],'blue')
    if len(cmp)>0:
        if cmp.geom_type == 'MultiPolygon':
            for p in cmp:
                patch = desc.PolygonPatch(p,fc='white',linewidth=0.5)
                ax.add_patch(patch) 
        else:
            patch = desc.PolygonPatch(cmp,fc='white',linewidth=0.5)
            ax.add_patch(patch) 
        
    if cut.geom_type == 'MultiPolygon':
        for p in cut:
            patch = desc.PolygonPatch(p,fc='none',linewidth=0.5)
            ax.add_patch(patch) 
    else:
        patch = desc.PolygonPatch(cut,fc='none',linewidth=0.5)
        ax.add_patch(patch) 
    if e_poly.geom_type == 'MultiPolygon':
        for p in e_poly:
            patch = desc.PolygonPatch(p,fc='none',linewidth=1,linestyle='--')
            ax.add_patch(patch) 
    else:
        patch = desc.PolygonPatch(e_poly,fc='none',linewidth=1,linestyle='--')
        ax.add_patch(patch)      
    #plt.title(title)
    plt.axis('equal')
    plt.show()
    plt.pause(0.1)


def LakeElev(poly,dem):
    c=np.vstack(poly.exterior.xy).transpose();
    c=dist2cord(c,corner=dem.corner,f=dem.scale_factor)
    
    z=dem.getElev(c) #perimeter elevation profile
    z=np.sort(z)[int(0.25*z.shape[0])] #use the 25th percentile value as water level)
    return z

    

def GenerateSTLs(Boundary=[],Rect_Pt=[],Rect_Pt_rot=[],rd_include=[],trail_exclude=[],waterway_include=[],waterbody=[],trail_gpx=[],
                 path_width=.7,support_width=0.9,path_clearance=0.1,height_factor=2,base_height=4,top_thickness=1.5,
                 edge_width=1.5,max_print_size=[248,198],tiles=1,water_drop=0.5,load_area=[],resolution=30,
                 dem_offset=[0,0],downsample_factor=1,map_only=False,compass_loc=[],compass_size=1.0):
    """

    INPUTS:

    Boundary - polygon to define shape of terrain, input as longitude/lattitude array, or gpx file of points
    Rect_Pt - points to bound with rectangular boundary
    Rect_Pt_rot - points to bound with a rotated rectangular boundary
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
    tiles - number of tiles to use for terrain print 
    water_drop - water ways printed slightly lower than other paths and terrain
    load_area - overide automatic area selection
    top_thickness - vertical thickness of wider top section of roads, trails and waterways.
    resolution - 10 or 30, for 10/30 meter resolution DEM.
    dem_offset - offset DEM relative to OSM data to account for shifts in data.
    downsample_factor, integer factor to reduce resolution of terrain surface, needed for larger models.  1mm resolution is reasonable minimum.
    map_only - stop after generating map to review (no boolean ops), recomend to do this first until all desired features look corrects oi it doesn't get hung up on boolean operations for hours.
    compass_loc - XY location for printed compass.
    compass_size - scale factor for size of compass (Letter N might reqiure 0.25 nozzle)
    """
    if not os.path.isdir('print_files/'):
        os.mkdir('print_files/')
    if not os.path.isdir('temp/'):
        os.mkdir('temp/')

    # clear temp files
    fileList = glob.glob('temp/*.stl')
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
    if len(Boundary)==0 and len(Rect_Pt)>1:
        Rect_Pt=np.array(Rect_Pt)
        Rect_Pt=np.vstack((Rect_Pt.min(0),[Rect_Pt[:,0].min(),Rect_Pt[:,1].max()],Rect_Pt.max(0),[Rect_Pt[:,0].max(),Rect_Pt[:,1].min()]))
        Boundary=Rect_Pt
    if len(Boundary)==0 and len(Rect_Pt_rot)>2:
        Rect_Pt_rot=np.flip(Rect_Pt_rot,axis=1)
        c=Rect_Pt_rot.min(axis=0)
        Rect_Pt_rot=shp.geometry.Polygon(cord2dist(xy=Rect_Pt_rot,corner=c))
        Boundary=Rect_Pt_rot.minimum_rotated_rectangle
        Boundary=np.flip(dist2cord(xy=np.array(Boundary.exterior.xy).transpose(),corner=c),axis=1)

    if not any(load_area): #load area bounds polygon unless larger area is input
        load_area=np.vstack((np.min(Boundary,axis=0),np.max(Boundary,axis=0)))
        
    Boundary=np.flip(Boundary,axis=1) #from here on in, we use x,y notation
    dem=Dem(Boundary,resolution,dem_offset,downsample_factor)
    
    # dem.plotElev2()
    #plotElev(dem)

    #find the largest print that fits (rotated) within print size a return scale factor and x/y=0 corner of print in coords'
    if tiles>1:
        scale_factor,corner,edge_poly,cutouts=printScaling_tiled(dem,Boundary,max_print_size,tiles)
    else:
        scale_factor,corner,edge_poly=printScaling(dem,Boundary,max_print_size)
        cutouts=[]

    # offsets=[support_width/2,path_width/2,(path_width+path_clearance)/2]
    offsets=[(path_width+path_clearance)/2]
    
    Footpaths=getFootpaths(Boundary,trail_exclude,trail_gpx,corner,scale_factor,offsets,base_height)
    
    Roads=getRoads(Boundary,rd_include,corner,scale_factor,offsets,base_height)
    
    Waterways=getWaterways(Boundary,waterway_include,corner,scale_factor,offsets,base_height,map_only)
    Waterbodies=getWaterbodies(Boundary,waterbody,corner,scale_factor,path_clearance,base_height,height_factor,dem)
    
    Boundary=shp.geometry.Polygon(cord2dist(xy=Boundary,corner=corner,f=scale_factor))

    borders=offsetPoly(Boundary,[-edge_width-path_clearance,-edge_width-path_clearance/2,-edge_width,0])
    
    print('2D Boolean Operations')
    Waterbodies,Footpaths,Roads,Waterways,Cutout=binaryOps2(Waterbodies,Footpaths,Roads,Waterways,borders,path_width,support_width,path_clearance)
    
    if len(compass_loc)==2: #redundant with other compass section
        print('generating compass')
        cmp=tm.load("compass.stl")
        t=np.eye(4)
        t[:2, :2] *= compass_size
        cmp.apply_transform(t) #scale
        
        cmp.apply_translation([compass_loc[0], compass_loc[1],0])
      
        c_poly = cmp.section(plane_origin=[0,0,-5],plane_normal=[0,0,1])

        # cmp=tm.intersections.mesh_plane(msh, [0,0,1], )
        p=[]
        for e in c_poly.entities:
            p.append(shp.geometry.Polygon(c_poly.vertices[e.points,:]))
        if len(p)>1:
            c_poly=shp.geometry.MultiPolygon(p)
        else:
            c_poly=shp.geometry.Polygon(p[0])
    else:
        c_poly=[]
    
    if not map_only:
        if len(compass_loc)==2:
           #get elevation around edge ofcompass polygon
           edge=[]
           for p in c_poly:
               edge.append(np.array(p.exterior.xy).T)
           edge=dist2cord(np.vstack(edge),corner=corner,f=scale_factor)
           # for i in range(len(pts)):
           #     idx[i]=compass_bounds.contains(pts[i])
           z=dem.getElev(edge)
           z=z-np.min(dem.z)
           z=z*scale_factor*height_factor+base_height+1

           
           #Cut off bottom to just below min elevation
           box=tm.creation.box([200,200,200])
           box.apply_translation([compass_loc[0], compass_loc[1],min(z)-2-100])
           
           cmp.apply_translation([0,0,max(z)])
           cmp.export('temp/c1.stl')
           box.export('temp/b1.stl')
           cmp=tm.boolean.difference([cmp,box])
           cmp.export('print_files/compass.stl')
           
           
           cmp_cut=tm.load("compass_cutout.stl")
           cmp_cut.apply_transform(t)
           cmp_cut.apply_translation([0,0,min(z)+50])
           box=tm.creation.box([200,200,200])
           cmp_cut.apply_translation([compass_loc[0], compass_loc[1],0])
           cmp_cut=[cmp_cut]
           # cmp_cut.export('temp/compass_cut.stl')
        else:
            cmp_cut=[]
        wb_heights2=[]
        if len(Waterbodies)>0:
            grid=np.meshgrid(dem.lat,dem.lon)
            # x=grid[1].flatten()
            # y=grid[0].flatten()
            x,y=cord2dist(x=grid[1].flatten(),y=grid[0].flatten(),corner=corner,f=scale_factor)
            pts=shp.geometry.MultiPoint(np.array([x,y]).transpose())
            
            if Waterbodies[1].geom_type == 'MultiPolygon':
                for i in range(len(Waterbodies[1])):
    
                    z=LakeElev(Waterbodies[1][i],dem)
                    wb_heights2.append((z-np.min(dem.z))*scale_factor*height_factor+1-3+base_height) #+1 to line up terrain
                    idx=includedPoints(Waterbodies[1][i],dem,corner,scale_factor)             
                    dem.z[idx]=z #set points to nominal elevation.     
            else:
                z=LakeElev(Waterbodies[1],dem)
                wb_heights2.append((z-np.min(dem.z))*scale_factor*height_factor+1-3+base_height) #+1 to line up terrain
                idx=includedPoints(Waterbodies[1],dem,corner,scale_factor)             
                dem.z[idx]=z #set points to nominal elevation. 
  
        print('Meshing')
        rd_msh=meshgen2(Roads,100,dem,corner,scale_factor,height_factor,base_height,fname='temp/road')
        fp_msh=meshgen2(Footpaths,100,dem,corner,scale_factor,height_factor,base_height,fname='temp/trail')
        ww_msh=meshgen2(Waterways,100,dem,corner,scale_factor,height_factor,base_height,fname='temp/water')
        wb_msh=meshgen_wb(Waterbodies,[20,1.5,3],'temp/waterbodies',wb_heights2)
        #cut_msh=meshgen(Cutout,100,corner,scale_factor)
        #border_msh=meshgen(borders[3],100,'temp/border',elev=-1)#shift down 1mm to avoid coplaner faces]

        
        if edge_poly.geom_type == 'MultiPolygon':
            terrain=[]
            for e in edge_poly:
                b=e.intersection(Boundary)
                t=tMesh(dem,b,scale_factor,corner,height_factor,base_height,water_drop)
                terrain.append(t)
                # terrain.append(tm.boolean.difference([t,cutouts]))
        else:
            terrain=[tMesh(dem,Boundary,scale_factor,corner,height_factor,base_height,water_drop)]
            
        terrain_all=tMesh(dem,Boundary,scale_factor,corner,height_factor,base_height,water_drop)
        terrain_all.export('temp/terrain.stl')
       
        
        print('Cutting Terrain Model')
        
        cutlist=[]
        i=1
        for m in [rd_msh,fp_msh,ww_msh,wb_msh,cmp_cut,cutouts]:
            
            if len(m)>0:
                cutlist.append(m[0])
                m[0].export('temp/cc'+str(i)+'.stl')
                i=i+1
        for i in range(len(terrain)):
            terrain[i].export('temp/t-'+str(i+1)+'.stl')
            terrain[i]=tm.boolean.difference([terrain[i]]+cutlist)
            terrain[i].export('print_files/terrain-'+str(i+1)+'.stl')

        print('Building Inserts')
        #cut top of top section
        
        fp_msh[1]=tm.boolean.intersection([fp_msh[1],terrain_all])
        
        if len(rd_msh)>0:
            rd_msh[1]=tm.boolean.intersection([rd_msh[1],terrain_all])
        if len(ww_msh)>0:
            ww_msh[1]=tm.boolean.intersection([ww_msh[1],terrain_all])

         #cut top of support section
        terrain_all.apply_translation([0,0,-1])
        
        fp_msh[2]=tm.boolean.intersection([fp_msh[2],terrain_all])
        
        if len(rd_msh)>0:
            rd_msh[2]=tm.boolean.intersection([rd_msh[2],terrain_all])
        if len(ww_msh)>0:
            ww_msh[2]=tm.boolean.intersection([ww_msh[2],terrain_all])

        
        terrain_all.apply_translation([0,0,-1])
             
        # cut off bottom of top section, and join with support section
        fp_msh[1]=tm.boolean.difference([fp_msh[1],terrain_all])
        fp_msh[2]=tm.boolean.union([fp_msh[1],fp_msh[2]])
        fp_msh[2].export('print_files/footpaths.stl')
        
        if len(rd_msh)>0:
            rd_msh[1]=tm.boolean.difference([rd_msh[1],terrain_all])
            rd_msh[2]=tm.boolean.union([rd_msh[1],rd_msh[2]])
            rd_msh[2].export('print_files/roads.stl')

            
        if len(ww_msh)>0:
            ww_msh[1]=tm.boolean.difference([ww_msh[1],terrain_all])
            ww_msh[2]=tm.boolean.union([ww_msh[1],ww_msh[2]])
            ww_msh[2].export('print_files/waterways.stl')

        
        if len(wb_msh)>0:
            wb_msh[1].apply_translation([0,0,1.5])
            wb_msh[1]=tm.boolean.union([wb_msh[1],wb_msh[2]])
            wb_msh[1].export('print_files/waterbodies.stl')


        
        
        
    plotPaths(dem,scale_factor,Waterbodies,Footpaths,Roads,Waterways,Boundary,Cutout,c_poly,edge_poly)
    


    print("COMPLETE")
