import overpy
import numpy as np
import shapely as shp
import time as time
from cord_funcs import *

api = overpy.Overpass()


        
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

        areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,1]),np.min(poly[:,0]),np.max(poly[:,1]),np.max(poly[:,0])))
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

    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,1]),np.min(poly[:,0]),np.max(poly[:,1]),np.max(poly[:,0])))
    num_attempt=0
    result=-1
    while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
        try:
            #roads and railways
            #result = api.query("way(" + areaStr + ") [""highway""];   (._;>;); out body;")
            # result = api.query("way(" + areaStr + ") [""railway""];   (._;>;); out body;")
            result = api.query("(way(" + areaStr + ") [""highway""];way(" + areaStr + ") [""railway""];);   (._;>;); out body;")
        except:
            num_attempt=num_attempt+1
            time.sleep(5)
            
    # for i in range(len(result.ways)):
    #     if "name" in result.ways[i].tags and result.ways[i].tags['name']=="P&W Subdivision":
    #         print(i)
    inc_ways=[]
    coords=[]
    for w in result.ways:
        if w.id in roads or ("name" in w.tags and w.tags["name"] in rd_names):
            if ("name" in w.tags and w.tags["name"]=="P&W Subdivision"):
                x=1
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

    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,1]),np.min(poly[:,0]),np.max(poly[:,1]),np.max(poly[:,0])))

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
        return []
    Thickness=3
    
    print(' ')
    print('Retrieving OSM waterbodies...')
    areaStr=str("%f, %f, %f, %f" % (np.min(poly[:,1]),np.min(poly[:,0]),np.max(poly[:,1]),np.max(poly[:,0])))
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

    # for r in result.relations:
    #     #find by OSM id or name
    #     if r.tags['name'] in bodies or r.id in bodies:
    #         #once we find a relation, extract and merge all the ways
    #         ways=[]
    #         for m in r.members:
    #             if m._type_value == 'way' and m.role == 'outer':
    #                 num_attempt=0
    #                 w_res=[]
    #                 while not str(w_res.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
    #                     try:
    #                         w_res = api.query("way(id:"+str(m.ref)+");    (._;>;); out body;")
    #                     except:
    #                         num_attempt=num_attempt+1
    #                         time.sleep(5)
    #                 ways.append(w_res.ways[0])
    #         wb.append(ways)
    w_id=np.array([w.id for w in result.ways])
    polys=[]
    for r in result.relations:
        #find by OSM id or name
        if r.tags['name'] in bodies or r.id in bodies:
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
    #now look for lakes made of individual ways
    # num_attempt=0
    # result=-1
    # while not str(result.__class__)=="<class 'overpy.Result'>" and num_attempt<=5:
    #     try:
    #         result = api.query("way(" + areaStr + ") [""water""];    (._;>;); out body;")
    #     except:
    #         num_attempt=num_attempt+1
    #         time.sleep(5)
    
    # for w in result.ways:
    #     if ('name' in w.tags and w.tags['name'] in bodies) or w.id in bodies:
    #         wb.append(w)
    # polys=[]
    # polys2=[]
    
    # for w in wb:
    #     if isinstance(w, list):
    #         lines=[]
    #         for w2 in w:
    #             c=np.array([[float(n.lon) for n in w2.nodes],[float(n.lat) for n in w2.nodes]]).transpose()
    #             c=cord2dist(xy=c,corner=corner,f=scale_factor)
    #             lines.append(shp.geometry.LineString(c))

    #         if len(lines)==1:
    #             polys.append(shp.geometry.Polygon(lines[0]))
    #         else:
    #             merged_line=lines.pop(0)
    #             added=True
    #             while len(lines)>0:
    #                 merged_line,lines,added=findConnLine(merged_line,lines)
    #                 if not added:
    #                     polys.append(shp.geometry.Polygon(merged_line))
    #                     merged_line=lines.pop(0)
    #             polys.append(shp.geometry.Polygon(merged_line))

    #     else:
    #         c=np.array([[float(n.lon) for n in w.nodes],[float(n.lat) for n in w.nodes]]).transpose()
           
    #         c=cord2dist(xy=c,corner=corner,f=scale_factor)

    #         polys.append(shp.geometry.Polygon(c))

    polys=shp.geometry.MultiPolygon(polys)

    return polys