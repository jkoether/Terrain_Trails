# BACKGROUND
This has been an ongoing project of mine for a few months now.  I run trails a lot and original the goal was to create a function to combine trail data from Openstreetmap and geospatial data to create a 3D printable terrain model with the path features as seperate integrated prints. The terrain data is downloaded from USGS (10 or 30m 3DEP data), overpass (overpy) is used for the OSM queries, and shapely and Blender are used for all the boolean operations.  Going into this I knew almost no Python so this has been a learning experience.

An input script is needed for a given area to define your max print size, bounding polygon coordinates and the roads, streams and lakes to incude.  The focus was on trails for me so all OSM footpaths are included by default with a list of exluded trails as an input.  I was able to get good results using a 0.25m nozzle for the paths and a path width of around 0.8 mm (3 passes), the terrain is printed with a 0.4mm nozzle.  The terrain is a long print if you maximize the size and use .10 layers, 2-3 days on a Mk3, but the other parts all print very quickly.  I found that if I print with a .10 mm layer height the layer lines are barely noticable, no post processing required.

It's not quite plug and play but you should be able to get it working for any area you want with a little work.  The great thing about OSM is that it is user maintained so you can add the trails you need and do any required clean-up on what is there.

This is currently only able to pull from downloaded elevation data.  I have found that 10m data can be easily downloaded from:
https://apps.nationalmap.gov/tnmaccess/products.html 
The current script is set up to read the geotiff files using that search.

Future options I would like to include are:
* Automatic sectioning and dovetailing for terrain prints larger than printbed
* Seprate groups of trails for different colored trails. (can be done with slicer...)
* Integrated compass for terrain model.
* Guidance for material changes at specific layers heights to generate contour lines on terrain model.
* Adjust base height each path section based on minimum elevation.  This would reduce priting time and reduce the number of tall skinny path prints that are prone to failing.

The following libraries are used:
* overpy
* Shapely
* numpy trimesh as tm
* triangle
* gpxpy
* rasterio

Blender is also required, exe loacation is assumed to be: "C:\Program Files\Blender Foundation\Blender 3.0\blender.exe


# USAGE:

In general to make a print you need to define a polygon using coordinates for the print area. you can also use a gpx file to input the area, I use caltopo to create these: caltopo.com > +Add > polygon > [draw polygon] > export.

I made this with a focus on trails so it will import all footpaths within the selected area, so you will need to exlude paths by name or OSM id.  Use the openstreetmap.org "query features" tool to get the ids for the paths to exclude.  Roads, waterways and waterbodies have to be explicitly included.  It's up to you how many roads you want to inlcde.

This is currently only able to pull from downloaded elevation data.  I have found that 10m data can be easily downloaded from: https://apps.nationalmap.gov/tnmaccess/products.html.

The 30m data is easier to find but I have found the quality of the 10m (1/3 arcsecond) data to be much better and not just because of the resolution, it is just more accurate and less noisy.

See the included example for North Park (pittsburgh area): NP_STL.py, NP.gpx, tiff files downloaded by searching bounding box="-80.034,40.577,-79.974,40.6222", NED 1/3 arc-second, GeoTIFF.

Final assembly can be a little tricky. You do want to clean off any fuzzies, and maybe cut really big path sections into several pieces.  Don't force them in too hard or they will get stuck partially inserted. A little sanding in the cutouts and on the side of the paths really helps the process.



# INPUTS:

* Boundary - polygon to define shape of terrain, input as longitude/lattitude array, or gpx file of points
* Rect_Pt - points to bound with a rotated rectangular boundary
* rd_include - roads names / ids to inlcude, no roads included by default.
* trail_exclude - footpath names / ids to inlcude, all included by default.
* waterway_include - water paths to inlcude, not inlcuding polygon water bodies
* waterbody - waterbody polygons (lakes, ponds) rivers will work in threory but the code is not set up to allow for any slope on waterbodies.
* path_width - path (road, footpaths and waterway) print top width, 3 total passes works best.
* support_width - width of base (shoot for one pass less than path_width)
* path_clearance - clearance between path prints (tops) and cutouts.
* height_factor - exaggeration factor for elevation. 1.0 makes on same scale as horizontal dimensions. 2-3 is reasonable factor.
* base_height - minium print thickness
* trail_gpx - gpx file to use for trail path
* edge_width - width of terrian border with no path cutouts
* max_print_size - maximum print dimension.  will scale and rotate print to fit this.
* water_drop - water ways printed slightly lower than other paths and terrain
* load_area - overide automatic area selection
* resolution - 10 or 30, for 10/30 meter resolution DEM.
* dem_offset - offset DEM relative to OSM data to account for shifts in data
* downsample_factor - integer factor to reduce resolution of terrain surface, needed for larger models.  1mm resolution is reasonable minimum.
* map_only - stop after generating map to review (no boolean ops), recomend to do this first until all desired features look correct so you can iterate through changes in your input deck quicker
* compass_loc - XY location for printed compass.
* compass_size - scale factor for size of compass 50x50 default

# February 2022 Updates:
Updated to add the following features:
* Base elevation for path cutouts is based on the minimum elevation for that path, this allow higher paths to not require the same depth cutouts which will reduce print times.
* Islands are in waterbodies are preserved.
* An initial tiling method has been implimented, currently only supports 2 tiles. This will rotate 2 prints to maximize the scale and then produce 2 terrain tile files. The tiles will include a dovetail cutout to connect the tiles using an insert.
* A compass can be inlcluded by defining the location and scale. Default scale is a 50x50mm compass.

Fixes:

* Corrected issue with DEM overlap, this was the main cause for the terrain misalignment so the dem_offset input will hopefully not be needed.
* More reliable terrain meshing methodology.
* Code re-organized, blender script removed.

I currently can't get reliable boolean operation using only OpenSCAD or Blender so the code requires both for now.


# January 2022 Updates:
The function has been significantly updated to streamline the path creation and boolean ops and is more flexible.  All of the path lines are now inflated using shapely, with a lot of boolean operation done using shapely polygons rather than the 3D meshes.  The function will now process inputs with no roads, waterways or waterbodies, and footpath can be input as a gpx file in which case no footpaths are read from OSM.  Mesh boolean ops are all performed in Blender now, and with the other changes the function usually runs in a few minutes, rather than hours.  The paths now include some smoothing and the buffer function is used to elimintate sections that are too small or thin to print, although there is still some room for improvement here.
