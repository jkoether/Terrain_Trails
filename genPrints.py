#!BPY

import bpy
import glob
import sys  
import os as os


def terrainShift(ob_name,new_name,offset):

    # get the object
    ob = bpy.context.scene.objects.get(ob_name)

    # copy, name and store new_name
    new_ob = ob.copy()
    new_ob.name = new_name
    new_ob.location.z += offset

    # place it in the collection
    col=ob.users_collection
    col=col[0]
    #bpy.context.scene.collection.children.link(newCol) <-- in case you need to add the collection to the scene
    col.objects.link(new_ob)
    print (col.name) 

def applyBoolean(operation, a, b):
    """For operation use either UNION, DIFFERENCE, or INTERSECT"""
    if isinstance(b,list):
        for bi in b:
            modifier = a.modifiers.new(type="BOOLEAN", name=bi.name+"_mod")
            modifier.object = bi
            modifier.operation = operation
            #modifier.double_threshold=10^-4
            #modifier.solver = 'EXACT'
            bpy.context.view_layer.objects.active = a
            bpy.ops.object.modifier_apply(modifier=bi.name+"_mod")
    else:
        modifier = a.modifiers.new(type="BOOLEAN", name=b.name+"_mod")
        modifier.object = b
        modifier.operation = operation
        #modifier.double_threshold=10^-4
        #modifier.solver = 'EXACT'
        bpy.context.view_layer.objects.active = a
        bpy.ops.object.modifier_apply(modifier=b.name+"_mod")
    # a.modifiers.clear()
    #bpy.ops.object.select_all(action='DESELECT')
    
top_thickness = float(sys.argv[-2])
water_offset = float(sys.argv[-1])


f = open("blender_log.txt", "w")
f.write('***')
f.write(str(top_thickness))
f.write(str(water_offset))
#delete existing meshes    
for o in bpy.context.scene.objects:
    if o.type == 'MESH':
        o.select_set(True)
    else:
        o.select_set(False)
bpy.ops.object.delete()

bpy.ops.import_mesh.stl(filepath='temp/terrain.stl')

bpy.ops.import_mesh.stl(filepath='temp/border_1.stl')
bpy.ops.import_mesh.stl(filepath='temp/cutout_1.stl')

from os.path import exists

if exists('temp/trail_1.stl'):
    bpy.ops.import_mesh.stl(filepath='temp/trail_1.stl')
    bpy.ops.import_mesh.stl(filepath='temp/trail_2.stl')
    includesTrail=True
else:
    includesTrail=False
    
if exists('temp/road_1.stl'):
    bpy.ops.import_mesh.stl(filepath='temp/road_1.stl')
    bpy.ops.import_mesh.stl(filepath='temp/road_2.stl')
    includesRoad=True
else:
    includesRoad=False
    
if exists('temp/water_1.stl'):
    bpy.ops.import_mesh.stl(filepath='temp/water_1.stl')
    bpy.ops.import_mesh.stl(filepath='temp/water_2.stl')
    includesWater=True
else:
    includesWater=False


bpy.ops.object.select_all(action='DESELECT')

# ## Terrain
applyBoolean('INTERSECT',  bpy.data.objects['terrain'],bpy.data.objects['border_1'])
applyBoolean('DIFFERENCE', bpy.data.objects['terrain'], bpy.data.objects['cutout_1'])

# # save STL
bpy.data.objects['terrain'].select_set(True)
bpy.ops.export_mesh.stl(filepath='print_files\Terrain.stl',use_selection=True)


#import fresh copy of terrain model
bpy.ops.object.delete()
bpy.ops.import_mesh.stl(filepath='temp/terrain.stl')
bpy.ops.object.select_all(action='DESELECT')


## form top surface of top sections
if includesTrail:
    applyBoolean('INTERSECT', bpy.data.objects['trail_1'], bpy.data.objects['terrain'])
if includesRoad:
    applyBoolean('INTERSECT', bpy.data.objects['road_1'], bpy.data.objects['terrain'])
if includesWater:
    applyBoolean('INTERSECT', bpy.data.objects['water_1'], bpy.data.objects['terrain'])

bpy.data.objects["terrain"].location.z += -1.0 #shift terrains

# form top surface of support sections
if includesTrail:
    applyBoolean('INTERSECT', bpy.data.objects['trail_2'], bpy.data.objects['terrain'])
if includesRoad:
    applyBoolean('INTERSECT', bpy.data.objects['road_2'], bpy.data.objects['terrain'])
if includesWater:
    applyBoolean('INTERSECT', bpy.data.objects['water_2'], bpy.data.objects['terrain'])

bpy.data.objects["terrain"].location.z += -1.0

if includesTrail:
    applyBoolean('DIFFERENCE', bpy.data.objects['trail_1'], bpy.data.objects['terrain']) # form bottom surface of top section
    applyBoolean('UNION', bpy.data.objects['trail_2'], bpy.data.objects['trail_1']) #combine

    #save to STL
    bpy.data.objects['trail_2'].select_set(True)
    bpy.ops.export_mesh.stl(filepath='print_files\Trails.stl',use_selection=True)
    bpy.data.objects['trail_2'].select_set(False)

if includesRoad:
    applyBoolean('DIFFERENCE', bpy.data.objects['road_1'], bpy.data.objects['terrain'])# form bottom surface of top section
    applyBoolean('UNION', bpy.data.objects['road_2'], bpy.data.objects['road_1']) #combine
    #save to STL
    bpy.data.objects['road_2'].select_set(True)
    bpy.ops.export_mesh.stl(filepath='print_files\Roads.stl',use_selection=True)
    bpy.data.objects['road_2'].select_set(False)

if includesWater:
    applyBoolean('DIFFERENCE', bpy.data.objects['water_1'], bpy.data.objects['terrain'])# form bottom surface of top section
    applyBoolean('UNION', bpy.data.objects['water_2'], bpy.data.objects['water_1']) #combine
    bpy.data.objects['water_2'].select_set(True)
    bpy.ops.export_mesh.stl(filepath='print_files\Waterways.stl',use_selection=True)
    bpy.data.objects['water_2'].select_set(False)



f.close