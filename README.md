# BACKGROUND
This has been an ongoing project of mine for a few months now.  I run trails a lot and original the goal was to create a function to combine trail data from Openstreetmap and geospatial data to create a 3D printable terrain model with the path features as seperate integrated prints. The terrain data is downloaded from USGS (10 or 30m 3DEP data), overpass (overpy) is used for the OSM queries, and OpensSCAD is used for all the boolean operations with a collection of Python functions to do all the processing.  The hardest part was creating a function to inflate (widen) the trail lines to a polygon that could be extruded and meshed.  There are several python libraries that do line/plolygon inflation, but none that really worked for an arbitray network of intersecting and deadend lines so I had start from scratch.  Going into this I knew almost no Python so this has been a learning experience.

The "inflator" function can take almost any complex set of multiple intersecting and/or seperate polylines and inflate to a specified width, handling any resulting overlap, and removing interior polygons that get too small or thin. This is then 2D meshed using a Delaunay triangulation (Triangle library) and extruded to a 3D watertight mesh using pymesh.  The intersection of this mesh and the terrain mesh is used so the top surface mataches the terrain exactly once inserted.  OpenSCAD seems very slow (several hours with 10m DEM) for boolean operations, but it does finish eventually and it does exactly what I want.

Once I had this I realized adding in roads and streams would be easy so the final product takes an input deck and produces a set of ready to print models (not multi-material) for the terrain, with cutouts for trails, roads, streams and lakes and the individual models for the trails, roads, streams and lakes. The clearance has to be tweaked for your setup, but once you get it right everything just clicks together with very few visible gaps.  The paths are all created with a short thicker section at the top so it is easier to get them in with a tight fit.

An input script is needed for a given area to define your printer size, bounding polygon coordinates and the roads, streams and lakes to incude.  The focus was on trails for me so all OSM footpaths are included by default with a list of exluded trails as an input.  I was able to get good results using a 0.25m nozzle for the paths and a path width of around 0.8 mm (3 passes), the terrain is printed with a 0.4mm nozzle.  The terrain is a long print if you maximize the size and use .10 layers, 2-3 days on a Mk3, but the other parts all print very quickly.  I found that if I print with a .10 mm layer height the layer lines are barely noticable, no post processing required.

It's not quite plug and play but you should be able to get working for any area you want with a little work.  The great thing about OSM is that it is user maintained so you can add the trails you need and do any required clean-up on what is there.


This is currently only able to pull from downloaded elevation data.  I have found that 10m data can be easily downloaded from:
https://apps.nationalmap.gov/tnmaccess/products.html 
The current script is set up to read the geotiff files using that search.

Future options I would like to include are:
* Automatic sectioning and dovetailing for terrain prints larger than printbed
* Seprate groups of trails for different colored trails. (can be done with slicer...)
* Additional printed path from GPX file.
* Integrated compass for terrain model.
* Guidance for material changes at specific layers heights to generate contour lines on terrain model.
* Adjust base height each path section based on minimum elevation.  This would reduce priting time and reduce the number of tall skinny path prints that are prone to failing.

The following libraries are used:
* overpy
* numpy trimesh as tm
* import pyclipper
* triangle
* gpxpy
* rasterio

OpenSCAD is also required, exe loacation is assumed to be: "C:\"Program Files\OpenSCAD\openscad.com"


# USAGE:

In general to make a print you need to define a polygon using coordinates for the print area. you can also use a gpx file to input the area, I use caltopo to create these: caltopo.com > +Add > polygon > [draw polygon] > export.

I made this with a focus on trails so it will import all footpaths within the selected area, so you will need to exlude paths by name or OSM id.  Use the openstreetmap.org "query features" tool to get the ids for the paths to exclude.  Roads, waterways and waterbodies have to be explicitly included.

This is currently only able to pull from downloaded elevation data.  I have found that 10m data can be easily downloaded from: https://apps.nationalmap.gov/tnmaccess/products.html.

The 30m data is easier to find but I have found the quality of the 10m data to be much better, and not just because of the resolution.

See the included example for North Park (pittsburgh area): NP_STL.py, NP.gpX, tiff files downloaded by searching bounding box="-80.034,40.577,-79.974,40.6222", NED 1/3 arc-second, GeoTIFF.




# INPUTS:

* CoordPoly - polygon to define shape of terrain, input as longitude/lattitude array, or gpx file of points
* rd_include - roads names / ids to inlcude, no roads included by default.
* trail_exclude - footpath names / ids to inlcude, all included by default.
* waterway_include - water paths to inlcude, not inlcuding polygon water bodies
* waterbody - waterbody polygons (lakes, ponds) rivers will work in threory but the code is not set up to allow for any slope on waterbodies.
* path_width - path (road, footpaths and waterway) print top width, 3 total passes works best.
* support_width - width of base (shoot for one pass less than path_width)
* path_clearance - clearance between path prints (tops) and cutouts.
* height_factor - exaggeration factor for elevation. 1.0 makes on same scale as horizontal dimensions. 2-3 is reasonable factor.
* base_height - minium print thickness
* edge_width - width of terrian border with no path cutouts
* max_print_size - maximum print dimension.  will scale and rotate print to fit this.
* water_drop - water ways printed slightly lower than other paths and terrain
* load_area - overide automatic area selection
* resolution - 10 or 30, for 10/30 meter resolution DEM.
* dem_offset - offset DEM relative to OSM data to account for shifts in data
* downsample_factor - integer factor to reduce resolution of terrain surface, needed for larger models.  1mm resolution is reasonable minimum.
* map_only - stop after generating map to review (no boolean ops), recomend to do this first until all desired features look corrects oi it doesn't get hung up on boolean operations for hours.
