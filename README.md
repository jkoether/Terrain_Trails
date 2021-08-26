This has been an ongoing project of mine for a few months now.  I run trails a lot and original the goal was to create a function to combine trail data from Openstreetmap and geospatial data to create a 3D printable terrain model with the path features as seperate integrated prints. The terrain data is downloaded from USGS (10 or 30m 3DEP data), overpass (overpy) is used for the OSM queries, and OpensSCAD is used for all the boolean operations with a collection of Python functions to do all the processing.  The hardest part was creating a function to inflate (widen) the trail lines to a polygon that could be extruded and meshed.  There are several python libraries that do line/plolygon inflation, but none that really worked for an arbitray network of intersecting and deadend lines so I had start from scratch.  Going into this I knew almost no Python so this has been a learning experience.

The "inflator" function can take almost any complex set of multiple intersecting and/or seperate polylines and inflate to a specified width, handling any resulting overlap, and removing interior polygons that get too small or thin. This is then 2D meshed using a Delaunay triangulation (Triangle library) and extruded to a 3D watertight mesh using pymesh.  The intersection of this mesh and the terrain mesh is used so the top surface mataches the terrain exactly once inserted.  OpenSCAD seems very slow (several hours with 10m DEM) for boolean operations, but it does finish eventually and it does exactly what I want.

Once I had this I realized adding in roads and streams would be easy so the final product takes an input deck and produces a set of ready to print models (not multi-material) for the terrain, with cutouts for trails, roads, streams and lakes and the individual models for the trails, roads, streams and lakes. The clearance has to be tweaked for your setup, but once you get it right everything just clicks together with very few visible gaps.  The paths are all created with a short thicker section at the top so it is easier to get them in with a tight fit.

An input script is needed for a given area to define your printer size, bounding polygon coordinates and the roads, streams and lakes to incude.  The focus was on trails for me so all OSM footpaths are included by default with a list of exluded trails as an input.  I was able to get good results using a 0.25m nozzle for the paths and a path width of around 0.8 mm (3 passes), the terrain is printed with a 0.4mm nozzle.  The terrain is a long print if you maximize the size and use .10 layers, 2-3 days on a Mk3, but the other parts all print very quickly.  I found that if I print with a .10 mm layer height the layer lines are barely noticable, no post processing required.

It's not quite plug and play but you should be able to get working for any area you want with a little work.  The great thing about OSM is that it is user maintained so you can add the trails you need and do any required clean-up on what is there.


This is currently only able to pull from downloaded elevation data.  I have found that 10m data can be easily downloaded from:
https://apps.nationalmap.gov/tnmaccess/products.html 
The current script is set up to read the geotiff files using that search.

Future options I would like to include are:
*Automatic sectioning and dovetailing for terrain prints larger than printbed
*Seprate groups of trails for different colored trails. (can be done with slicer...)
*Additional printed path from GPX file.
*Integrated compass for terrain model.
*Guidance for material changes at specific layers heights to generate contour lines on terrain model.
*Adjust base height each path section based on minimum elevation.  This would reduce priting time and reduce the number of tall skinny path prints that are prone to failing.

