# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 16:43:16 2022

@author: jkoet
"""
import trimesh as tm
import shapely as shp
from cord_funcs import *
from copy import deepcopy


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

def yDovetails(p1,p2,rp):
    #locate dovetails between two adjacent polygons, assuming polygons share a horizontal edge with p2 on top
    loc=[]
    if p1.bounds[3]>p2.bounds[1]+10**-6: #adjacent tiles shouldn't overlap, slight buffer for precision error.
        print('tiles overlap (Y)')
        return []
    #lower poly
    p1=p1.intersection(rp)
    if p1.geom_type=='GeometryCollection':
        return[]
    xy=np.array(p1.exterior.xy).transpose()
    xy=xy[xy[:,1]>p1.bounds[3]-.1]
    edge_pts=[xy[:,0].min(),xy[:,0].max()]
    
    #upper poly
    p2=p2.intersection(rp)
    if p2.geom_type=='GeometryCollection':
        return[]
    xy=np.array(p2.exterior.xy).transpose()
    xy=xy[xy[:,1]<p1.bounds[3]+.1]
    edge_pts=np.vstack((edge_pts,[xy[:,0].min(),xy[:,0].max()]))
    
    edge_pts=[np.max(edge_pts[:,0]),np.min(edge_pts[:,1])]
    if abs(np.diff(edge_pts))>150:
        loc.append([edge_pts[0]+np.diff(edge_pts).item()/4,p1.bounds[3],0])
        loc.append([edge_pts[0]+3*np.diff(edge_pts).item()/4,p1.bounds[3],0])
    elif abs(np.diff(edge_pts))>80:
        loc.append([edge_pts[0]+np.diff(edge_pts).item()/2,p1.bounds[3],0])
    return loc


def xDovetails(p1,p2,rp):
    #locate dovetails between two adjacent polygons, assuming polygons share a vertical edge with p2 on the right
    loc=[]
    if p1.bounds[2]>p2.bounds[0]+10**-6: #adjacent tiles shouldn't overlap, slight buffer for precision error.
        print('tiles overlap (X)')
        return []
    #left poly
    p1=p1.intersection(rp)
    if p1.geom_type=='GeometryCollection':
        return[]
    xy=np.array(p1.exterior.xy).transpose()
    xy=xy[xy[:,0]>p1.bounds[2]-.1]
    edge_pts=[xy[:,1].min(),xy[:,1].max()]
    
    #right poly
    p2=p2.intersection(rp)
    if p2.geom_type=='GeometryCollection':
        return[]
    xy=np.array(p2.exterior.xy).transpose()
    xy=xy[xy[:,0]<p1.bounds[2]+.1]
    edge_pts=np.vstack((edge_pts,[xy[:,1].min(),xy[:,1].max()]))
    
    edge_pts=[np.max(edge_pts[:,0]),np.min(edge_pts[:,1])]
    if abs(np.diff(edge_pts))>150:
        loc.append([p1.bounds[2],edge_pts[0]+np.diff(edge_pts).item()/4,0])
        loc.append([p1.bounds[2],edge_pts[0]+3*np.diff(edge_pts).item()/4,0])
    elif abs(np.diff(edge_pts))>80:
        loc.append([p1.bounds[2],edge_pts[0]+np.diff(edge_pts).item()/2,0])
    return loc

def splitPoly(poly,y):
    
    split_line=shp.geometry.LineString([[poly.bounds[0]-1,poly.bounds[1]+y],[poly.bounds[2]+1,poly.bounds[1]+y]])
    split_poly=np.array(shp.ops.split(poly, split_line))
    #find section(s) above line
    idx=np.array([p.bounds[3] for p in split_poly])>poly.bounds[1]+y
    upper=split_poly[idx]
    lower=split_poly[~idx]
    if np.all(idx):
        lower=-1
    elif lower.shape[0]>1: #split section might result in 2 polygons
        lower=shp.geometry.MultiPolygon(lower.tolist())
    else:
        lower=lower[0]
    if not np.any(idx):
        upper=shp.geometry.Polygon([])
    elif upper.shape[0]>1: #split section might result in 2 polygons
        upper=shp.geometry.MultiPolygon(upper.tolist())
    else:
        upper=upper[0]
    return lower,upper

def DovetailInserts(edge_poly,rp,dovetail_height):
    dovetail_gap=-.1 #negative for interference
    dovetail_spacing=100 #set to zero for just one per edge
    
    insert=tm.load('dovetail_insert.stl')
    t=np.eye(4)
    t[2, 2] *= (dovetail_height-0.6)/10
    insert.apply_transform(t)
    insert.export('print_files/dovetail_insert.stl')
    
    cutout=tm.load('dovetail_cutout.stl')
    t=np.eye(4)
    t[2, 2] *= (dovetail_height+2)/10
    cutout.apply_transform(t)
    cutout.apply_translation([0,0,-2]) #shift down 2mm to prevent co-linear edges
    cutouts=[]
    
    cut_loc=[]
    for i in range(len(edge_poly)-1):
        #Y+
        line=shp.geometry.LineString([[edge_poly[i].bounds[0],edge_poly[i].bounds[3]+.1],[edge_poly[i].bounds[2],edge_poly[i].bounds[3]+.1]])
        for j in range(i+1,len(edge_poly)):
            if edge_poly[j].intersects(line):
               cut_loc=cut_loc+yDovetails(edge_poly[i],edge_poly[j],rp)
        #Y-
        line=shp.geometry.LineString([[edge_poly[i].bounds[0],edge_poly[i].bounds[1]-.1],[edge_poly[i].bounds[2],edge_poly[i].bounds[1]-.1]])
        for j in range(i+1,len(edge_poly)):
            if edge_poly[j].intersects(line):
                cut_loc=cut_loc+yDovetails(edge_poly[j],edge_poly[i],rp)

    for c in cut_loc:
        new_cutout=deepcopy(cutout)
        new_cutout.apply_translation(c)
        cutouts.append(new_cutout)
        
    cutout=rotMesh(cutout,np.pi/2) #rotate cutout sideways for all vertical seams.
    cut_loc=[]
    
    for i in range(len(edge_poly)-1):
        #X+
        line=shp.geometry.LineString([[edge_poly[i].bounds[2]+.1,edge_poly[i].bounds[1]+10],[edge_poly[i].bounds[2]+.1,edge_poly[i].bounds[3]-10]])
        for j in range(i+1,len(edge_poly)):
            if edge_poly[j].intersects(line):
               cut_loc=cut_loc+xDovetails(edge_poly[i],edge_poly[j],rp)
        #X-
        line=shp.geometry.LineString([[edge_poly[i].bounds[0]-.1,edge_poly[i].bounds[1]+10],[edge_poly[i].bounds[0]-.1,edge_poly[i].bounds[3]-10]])
        for j in range(i+1,len(edge_poly)):
            if edge_poly[j].intersects(line):
                cut_loc=cut_loc+xDovetails(edge_poly[j],edge_poly[i],rp)

    for c in cut_loc:
        new_cutout=deepcopy(cutout)
        new_cutout.apply_translation(c)
        cutouts.append(new_cutout)
    # cutout=rotMesh(cutout,np.pi):    
    # add X-direction inserts
    
    cutouts=tm.boolean.union(cutouts)
    
    
    return cutouts
def printScaling_tiled(dem,Boundary,print_size,tiles,dovetail_height):
    
    print('Optimizing print size and tile layout.')
    #determine largest size that can be printed with the given number of tiles (only 2 supported for now.)

    

    x=dem.lon
    y=dem.lat
    
    corner=np.stack((np.min(x),np.min(y)))

    bourder_poly=cord2dist(Boundary,corner=corner) #border in meters

    print_angle=np.linspace(0,np.pi-np.pi/180,180) 
   
    
    
    #set up combinations of angles and different layouts. currently all tiles are the same orientation.
    factors=np.arange(int(tiles))+1
    rows=factors[np.array([tiles%f for f in factors])==0]
    rows=np.tile(rows,[print_angle.shape[0],1])
    print_angle=np.tile(print_angle,[rows.shape[1],1]).flatten()
    rows=rows.flatten('F')
    
    scale=np.zeros(print_angle.shape[0])
    corners=np.ones([tiles,2,print_angle.shape[0]])*np.nan
    #find the rotation that allows for the largest scale print

    for i in range(print_angle.shape[0]):
        if i==275:
            stop=0
        col=int(tiles/rows[i])
        
        #rotated boundary poly:
        rp=shp.affinity.rotate(shp.geometry.Polygon(bourder_poly), print_angle[i],use_radians=True,origin=[0,0])
        scale[i]=print_size[1]*rows[i]/(rp.bounds[3]-rp.bounds[1])
        #split on line parrallel to x-axis, find scaling assuming y-axis length is limiting
        scale_acceptable=False
        # row_height=(rp.bounds[3]-rp.bounds[1])/rows[i]
        while not scale_acceptable:
            scale_acceptable=True
            row_height=print_size[1]/scale[i]
            for j in range(rows[i]): #loop through each row
                if j==rows[i]-1: #don't need to split top/last row.
                    poly_row=rp
                else:
                    poly_row,rp=splitPoly(rp,row_height)
                width=(poly_row.bounds[2]-poly_row.bounds[0])
                if scale[i]>print_size[0]*col/width:
                    scale[i]=print_size[0]*col/width 
                    rp=shp.affinity.rotate(shp.geometry.Polygon(bourder_poly), print_angle[i],use_radians=True,origin=[0,0])
                    scale_acceptable=False  #j loop needs re-run if the scale factor needs decreased.
                    break
                #offset print so edge pieces are equal width
                offset=(width-(np.ceil(width/(print_size[0]/scale[i]))*print_size[0]/scale[i]))/2
                corners[j*col:(j+1)*col,0,i]=offset+poly_row.bounds[0]+np.arange(col)*print_size[0]/scale[i]
                corners[j*col:(j+1)*col,1,i]=np.tile(poly_row.bounds[1],[1,col])
                if rp.length==0: #no more rows needed
                    break
                
    idx=np.argmax(scale)
    print_angle=print_angle[idx]
    
    print('print angle: {0:.2f} deg'.format(print_angle*180/np.pi))
    
    scale=scale[idx] #use angle that allows the largest print.  model is not rotated, optimal angle is output and part will need rotated in slicer.
    corners=corners[:,:,idx]*scale
    corners=corners[np.any(~np.isnan(corners),1),:] #drop unused tiles

    rows=rows[idx] 

    #the seams between adjacent tiles need located so we can place the dovetail inserts/
    #rotate the tiles so they are aligned with X/Y axes.
    rp=shp.affinity.rotate(shp.geometry.Polygon(bourder_poly*scale), print_angle,use_radians=True,origin=[0,0])
    
    
    edge_poly=[]

    pts=np.array([[0, 0],[0,print_size[1]],[print_size[0],print_size[1]],[print_size[0],0]])
    
    #set up tile edges
    for i in range(corners.shape[0]):
        edge_poly.append(shp.geometry.Polygon(corners[i,:]+pts))
        if (edge_poly[-1].intersection(rp)).length==0:
            del edge_poly[-1]
    
    if dovetail_height>0:
        cutouts=DovetailInserts(edge_poly,rp,dovetail_height)
        cutouts=rotMesh(cutouts,-print_angle)
        
        
    else:
        cutouts=[]
    
    edge_poly=shp.geometry.MultiPolygon(edge_poly)
    edge_poly=shp.affinity.rotate(edge_poly, -print_angle,use_radians=True,origin=[0,0])
    x,y=cord2dist(x=x,y=y,corner=corner)
    print('DEM resolution: {0:.2f}x{1:.2f} mm'.format(((x.max()-x.min())/(x.shape[0]-1)*scale),(y.max()-y.min())/(y.shape[0]-1)*scale))
    dem.scale_factor=scale
    return scale,corner,edge_poly,[cutouts]