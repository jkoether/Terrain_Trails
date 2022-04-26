import overpy
import numpy as np
import shapely as shp
import time as time
from cord_funcs import *
import gpxpy as gpxpy
api = overpy.Overpass()


        
def getFootpaths(result,trail_exclude,trail_include,trail_gpx,corner,scale_factor,offsets,base):

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
        w_id=np.array([w.id for w in result.ways])
        inc_ways=[]
        coords=[]
        if len(trail_include)>0:
            for w in result.ways:

                if (("highway" in w.tags) and (w.tags['highway'] in ['path','footway','cycleway'])) and ((w.id in trail_include) or ("name" in w.tags and w.tags["name"] in trail_include)):
                    inc_ways.append(w.id)

            for r in result.relations:
                #find by OSM id or name
                if ("name" in r.tags and r.tags['name'] in trail_include) or r.id in trail_include:
                    #once we find a relation, extract and merge all the ways
                    for m in r.members:
                        w=result.ways[np.argmax(np.isin(w_id,m.ref))]
                        # if m._type_value == 'way':
                        if (("highway" in w.tags) and (w.tags['highway'] in ['path','footway','cycleway'])) and w.id not in inc_ways:
                            # outer.append(result.ways[idx])
                            inc_ways.append(w.id)
            
        else:
            for w in result.ways:
                if  ('highway' in w.tags and w.tags['highway'] in ['path','footway','cycleway']) and not (w.id in trail_exclude) and not ("name" in w.tags and w.tags["name"] in trail_exclude):
                    inc_ways.append(w.id)

        for num in inc_ways:
            w=result.ways[np.argmax(np.isin(w_id,num))]
            c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
            c=cord2dist(xy=c,corner=corner,f=scale_factor)
            coords.append(c)
        lines = shp.geometry.MultiLineString(coords)
    # p=offsetLines(lines,offsets)
    # return p
    return lines

def getRoads(result,roads,corner,scale_factor,offsets,base):
    if len(roads)==0:
        return []

    rd_names=[]
    if any(roads):
        for i in range(len(roads)-1,-1,-1):
            if isinstance(roads[i],str):
                rd_names.append(roads.pop(i))

    w_id=np.array([w.id for w in result.ways])
    inc_ways=[]
    coords=[]
    for w in result.ways:
        if ("highway" in w.tags and (w.tags['highway'] not in ['path','footway','cycleway'])) and (w.id in roads or ("name" in w.tags and w.tags["name"] in rd_names)):
            inc_ways.append(w.id)
        if "railway" in w.tags and (w.id in roads or ("name" in w.tags and w.tags["name"] in rd_names)):
            inc_ways.append(w.id)
            
    for r in result.relations:
        #find by OSM id or name
        if ("name" in r.tags and r.tags['name'] in rd_names) or r.id in roads:
            #once we find a relation, extract and merge all the ways
            for m in r.members:
                    w=result.ways[np.argmax(np.isin(w_id,m.ref))]
                    # if m._type_value == 'way':
                    if (("highway" in w.tags) and (w.tags['highway'] not in ['path','footway','cycleway'])) and w.id not in inc_ways:
                        # outer.append(result.ways[idx])
                        inc_ways.append(w.id)
    for num in inc_ways:
        w=result.ways[np.argmax(np.isin(w_id,num))]
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



def getWaterways(result,waterways,corner,scale_factor,offsets,base,map_only):
    if len(waterways)==0:
        return []
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

        

def getWaterbodies(result,bodies,corner,scale_factor,clearance,base,height_factor,dem):
    if len(bodies)==0:
        return []
    Thickness=3

    w_id=np.array([w.id for w in result.ways])
    polys=[]
    for r in result.relations:
        #find by OSM id or name
        if ("name" in r.tags and r.tags['name'] in bodies) or r.id in bodies:
            #once we find a relation, extract and merge all the ways
            outer=[]
            inner=[]
            
            for m in r.members:
                idx=np.argmax(np.isin(w_id,m.ref))
                if m._type_value == 'way' and m.role == 'outer':
                    outer.append(result.ways[idx])
                if m._type_value == 'way' and m.role == 'inner':
                    inner.append(result.ways[idx])
            holes=[]
            p=[]
            for w in inner:
                c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).T
                c=cord2dist(xy=c,corner=corner,f=scale_factor)
                holes.append(c)
            for w in outer:
                p.append(np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).T)
            p=np.vstack(p)
            p=cord2dist(xy=p,corner=corner,f=scale_factor)
            polys.append(shp.geometry.Polygon(p,holes=holes))
       
    #some lakes are single ways
    for w in result.ways:
        #find by OSM id or name
        if w.id in bodies or ("name" in w.tags and w.tags["name"] in bodies):
            #once we find a relation, extract and merge all the ways
            p=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).T
            p=cord2dist(xy=p,corner=corner,f=scale_factor)
            polys.append(shp.geometry.Polygon(p))    


    polys=shp.geometry.MultiPolygon(polys)

    return polys