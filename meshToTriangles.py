import bpy
import math
import numpy as np
import mathutils

bl_info = {
    "name": "Convert mesh into forge objects",
    "author": "Yolomcswag",
    "version": (1, 3),
    "blender": (3, 3, 1),
    "location": "View3D > Side Panel",
    "description": "Adds a 1 button solution to convert a selected mesh into forge objects.",
    "warning": " WARNING: This script requires references to 'Floor Angled Standard A' and 'Primitive Triangle' to function, simply place them in the scene to add a reference.",
}

scaleFactor = 0.025
minimumThickness = 0.05

def forward_up_right_to_euler(forward, up, right):
    forward = np.array(forward)
    up = np.array(up)
    right = np.array(right)

    # Ensure vectors are normalized
    forward /= np.linalg.norm(forward)
    up /= np.linalg.norm(up)
    right /= np.linalg.norm(right)

    # Calculate the rotation matrix
    rotation_matrix = np.column_stack((-forward, right, up))

    # Calculate Euler angles from the rotation matrix
    sy = math.sqrt(rotation_matrix[0, 0]**2 + rotation_matrix[1, 0]**2)
    
    singular = sy < 1e-6

    if not singular:
        x = math.atan2(rotation_matrix[2, 1], rotation_matrix[2, 2])
        y = math.atan2(-rotation_matrix[2, 0], sy)
        z = math.atan2(rotation_matrix[1, 0], rotation_matrix[0, 0])
    else:
        x = math.atan2(-rotation_matrix[1, 2], rotation_matrix[1, 1])
        y = math.atan2(-rotation_matrix[2, 0], sy)
        z = 0

    return x, y, z

def displayVectors(forward, up, right, pivot): #Debug vector visualizer
    # Create a new mesh object
    mesh = bpy.data.meshes.new("newMesh")

    # Create vertices
    vertices = [(0,0,0), up * 1.5, forward, right]

    # Create edges (specify vertex indices for each edge)
    edges = [(0, 1), (0, 2), (0, 3)]

    # Link the mesh to the current scene
    scene = bpy.context.scene
    newMesh = bpy.data.objects.new("vectorVisualizer", mesh)

    # Set the mesh's vertices and edges
    mesh.from_pydata(vertices, edges, [])

    # Update the mesh geometry
    mesh.update()

    # Optionally, set the object's location
    newMesh.location = (pivot)
    bpy.data.collections["Triangles"].objects.link(newMesh)
    
    return 0

def meshToObjects():
    #Create "Triangles" collection if it doesn't exist
    collectionFound = False

    for myCol in bpy.data.collections:
        if myCol.name == "Triangles":
            collectionFound = True
            print ("Collection found in scene")
            break
        
    if not collectionFound:
        bpy.ops.collection.create(name  = "Triangles")
        bpy.context.scene.collection.children.link(bpy.data.collections["Triangles"])

    parent_collection = bpy.data.collections.get("Triangles")
    sub_collection = parent_collection.copy()

    #Reference triangle, to be eventually replaced with a good solution
    rightTriangle = bpy.data.objects['Floor Angled Standard A']
    floorAngledStandardA_ID = rightTriangle.forge_object_id

    isoscelesTriangle = bpy.data.objects['Primitive Triangle']
    primitiveTriangle_ID = isoscelesTriangle.forge_object_id


    objectCount = 0;
    hasNgons = False
    hasMoved = False

    #This is where the fun begins

    # Ensure an object is selected
    if bpy.context.active_object is not None:
        obj = bpy.context.active_object
        initialLocation = obj.matrix_world.translation.copy()
        obj.location = (0,0,0) #My code doesnt work if the mesh is not at 0,0,0 so instead of fixing the code,
        hasMoved = True        #I'm taking the mesh and pushing it somewhere else
        bpy.context.view_layer.update()

        # Check if the selected object is a mesh
        if obj.type == 'MESH':
            # Get the mesh data
            mesh = obj.data

            # Access the faces of the mesh
            faces = mesh.polygons
            
            newCollection = 0
            generateGroups = True
            uiGenerateGroups = bpy.context.scene.subgroupsToggle
            #Generate sub collections for large meshes
            print(uiGenerateGroups)
            
            if len(faces) > 150 and uiGenerateGroups: #First sub collection
                sub_collection = bpy.data.collections.new("Group")
                parent_collection.children.link(sub_collection)
                
            else:
                generateGroups = False
                sub_collection = parent_collection
            
            # Iterate through each face in the mesh
            for face in faces:
                
                isIsosceles = False
                isEquallateral = False
                isRight = False
                isScalene = False
                
                # Access face data
                vertices = face.vertices  # List of vertex indices for the face

                if len(vertices) > 0:
                    #===============================Start of triangle logic======================================
                    vertex_indices = face.vertices
                    
                    if len(vertex_indices) == 3:
                        
                        if newCollection > 148 and generateGroups and uiGenerateGroups:
                            #Adds additional sub collections every ~150 objects
                            newCollection = 0
                            sub_collection = bpy.data.collections.new("Group")
                            parent_collection.children.link(sub_collection)
                        
                        
                       #  Extract vertex positions using their indices
                        vA = obj.matrix_world @ obj.data.vertices[vertices[0]].co
                        vB = obj.matrix_world @ obj.data.vertices[vertices[1]].co
                        vC = obj.matrix_world @ obj.data.vertices[vertices[2]].co
                        

                       #  Create an excessive amount of variables
                        edgeAB = vB - vA
                        edgeAC = vC - vA
                        edgeBC = vC - vB
                        
                        pivotVertex = vA
                        vertexName = "vA"
                        
                        world_matrix = obj.matrix_world
                        normal = face.normal
                        world_normal = (world_matrix @ face.normal).normalized()
                        offset = mathutils.Vector((0.0, 0.0, 0.0))
                        
                        forwardVector = edgeAB.normalized()
                        rightVector = edgeAC.normalized()
                        
                        sides = (edgeAB.length, edgeAC.length, edgeBC.length)
                        edges = (edgeAB, edgeAC, edgeBC)
                        vert = (vC, vB, vA)
                        
                        scaleX = 1
                        scaleY = 1
                        t2scaleX = 1
                            
                        edgeAB = edgeAB.normalized()
                        edgeAC = edgeAC.normalized()
                        edgeBC = edgeBC.normalized()
                        
                        index = 0
                        tempArray = [sides[0], sides[1], sides[2]] 
                        tempVert = [vert[0], vert[1], vert[2]]
                        sidesShifted = sides
                        vertexShifted = vert
                        angle = [0,0,0]
                        offset = vA
                        
                        for side in sides: #Calculate each angle in triangle
                            angle[index] = ((tempArray[1] ** 2) + (tempArray[2] ** 2) - (tempArray[0] ** 2)) / (2 * (tempArray[1] * tempArray[2]))
                            angle[index] = math.acos(angle[index])

                            first = tempArray[0]
                            tempArray[0] = tempArray[1]
                            tempArray[1] = tempArray[2]
                            tempArray[2] = first
                            index = index + 1
                            
                        index = 0    
                        for side in sides: #Shift side angle and vertex arrays until largest side is index 0
                            #I could have sorted it but this was the easiest way I could think of to preserve the relationship of the sides/angles/vertexes by index
                            if sides[index] == max(sides):
                                sidesShifted = tempArray
                                vertexShifted = tempVert
                                break
                            
                            first = tempArray[0]
                            tempArray[0] = tempArray[1]
                            tempArray[1] = tempArray[2]
                            tempArray[2] = first
                            
                            first = tempVert[0]
                            tempVert[0] = tempVert[1]
                            tempVert[1] = tempVert[2]
                            tempVert[2] = first
                            
                            first = angle[0]
                            angle[0] = angle[1]
                            angle[1] = angle[2]
                            angle[2] = first
                            
                            index = index + 1
                        
                        hypotenuse = sidesShifted[0] #Make a hypotenuse variable for readability
                        
                        #==================================Determine triangle type=======================================
                        
                        #Right angle triangle
                        if abs((hypotenuse ** 2) - ((sidesShifted[1] ** 2) + (sidesShifted[2] ** 2))) <= 0.0001:
                            #Using inverse pythagorean theorem to check if right triangle
                            isRight = True

                            pivotVertex = vertexShifted[0]
                            vertexName = "vC"
                            forwardVector = (vertexShifted[0] - vertexShifted[1]).normalized()
                            rightVector = forwardVector.cross(world_normal).normalized()

                            scaleX = sidesShifted[1]
                            scaleY = sidesShifted[2]
                            
                            #Spawn triangle and add forge id
                            new_triangle = bpy.data.objects.new(name="Floor Angled Standard A", object_data=rightTriangle.data.copy())
                            new_triangle.forge_object_id = floorAngledStandardA_ID
                            
                            
                        #Equallateral triangle    
                        elif abs(sides[0] - sides[1]) <= 0.001 and abs(sides[1] - sides[2]) <= 0.0001: #Check if all sides are equal
                            isEquallateral = True
                            
                            offset = vA + forwardVector * (sides[0] / 2)
                            rightVector = forwardVector
                            forwardVector = (vC - offset).normalized()
                            pivotVertex = offset + (forwardVector * (sides[0] / 2)) - (world_normal * .1)
                            scaleX = sides[0]
                            
                            #Spawn triangle and add forge id
                            new_triangle = bpy.data.objects.new(name="Primitive Triangle", object_data=isoscelesTriangle.data.copy())
                            new_triangle.active_material = bpy.data.materials.get("Yellow")
                            new_triangle.forge_object_id = primitiveTriangle_ID
                        
                        
                        #Isosceles triangle    
                        elif abs(sides[0] - sides[1]) <= 0.001 or abs(sides[1] - sides[2]) <= 0.0001: #Check if 2 sides are equal
                            isIsosceles = True
                            index = 0
                            otherSide = sides[0]
                            
                            if abs(sides[0] - sides[1]) <= 0.001: #Find outlier side and save one of the remaining sides for reference
                                index = 2
                                otherSide = sides[0]
                                
                            if abs(sides[1] - sides[2]) <= 0.001:
                                index = 0
                                otherSide = sides[1]
                                
                            if abs(sides[0] - sides[2]) <= 0.001:
                                index = 1
                                otherSide = sides[0]
                                
                            outlier = sides[index]
                            midLine = math.sqrt(abs((otherSide ** 2) - ((outlier / 2) ** 2)))
                            forwardVector = edges[index].normalized()
                            rightVector = forwardVector.cross(world_normal).normalized()
                            
                            offset = rightVector * ((midLine / 2.36602)) #this number is used to correct for 343's terrible choice of object origin
                            pivotVertex = vert[index] + offset - (world_normal * minimumThickness)
                            
                            scaleX = midLine
                            scaleY = outlier
                            
                            #Spawn triangle and add forge id
                            new_triangle = bpy.data.objects.new(name="Primitive Triangle", object_data=isoscelesTriangle.data.copy())
                            new_triangle.forge_object_id = primitiveTriangle_ID
                        
                        
                        #Scalene triangle    
                        else: 
                            isScalene = True
      
                            #Spawn 2 triangles and add forge ids
                            new_triangle = bpy.data.objects.new(name="Floor Angled Standard A", object_data=rightTriangle.data.copy())
                            new_triangle2 = bpy.data.objects.new(name="Floor Angled Standard A", object_data=rightTriangle.data.copy())
                            new_triangle.forge_object_id = floorAngledStandardA_ID
                            new_triangle2.forge_object_id = floorAngledStandardA_ID
                            
                            new_triangle2.active_material = bpy.data.materials.get("Green")
                            
                            forwardVector = (vertexShifted[2] - vertexShifted[1]).normalized()
                            rightVector = forwardVector.cross(world_normal).normalized()
                            bisectLength = sidesShifted[1] * math.sin(angle[2])
                            pivotVertex = vertexShifted[0] - (rightVector * bisectLength)

                            scaleX = (pivotVertex - vertexShifted[1]).length
                            t2scaleX = (pivotVertex - vertexShifted[2]).length
                            scaleY = bisectLength
                            
                            #print("bisectLength: " + str(bisectLength))
                            
                            
                        #Add to object count total
                        if isIsosceles or isEquallateral or isRight:
                            objectCount = objectCount + 1
                            newCollection = newCollection + 1
                        elif isScalene:
                            objectCount = objectCount + 2 #Add 2 to the object count if triangle shape is cringe
                            newCollection = newCollection + 2
                        

                        print(sides)
                        print("isIsosceles = " + str(isIsosceles) + ", isEquallateral = " + str(isEquallateral) + ", isRight = " + str(isRight) + ", Pivot = " + str(vertexName))
                        
                        #=============================== Danger: rotation math ======================================
                        
                        # Define forward, right and up vectors. Some values are inverted to fix object rotation mismatches 
                        
                        if isRight:
                            forward_vector = -forwardVector.normalized()
                            right_vector = -rightVector.normalized()
                            up_vector = world_normal
                        
                        if isEquallateral: 
                            forward_vector = forwardVector.normalized()
                            right_vector = rightVector.normalized()
                            up_vector = world_normal
                            
                        if isIsosceles: 
                            forward_vector = -rightVector.normalized()
                            right_vector = forwardVector.normalized()
                            up_vector = world_normal
                            
                        if isScalene:
                            forward_vector = forwardVector.normalized()
                            right_vector = rightVector.normalized()
                            up_vector = world_normal


                        #Apply rotations   
                        
                        if isRight or isEquallateral or isIsosceles or isScalene: 
                            #print(up_vector)
                            #print(forward_vector)
                            #print(right_vector)
                            
                            #Chat gbt wrote this function, I don't understand it and don't want to
                            new_triangle.rotation_euler = forward_up_right_to_euler(forward_vector, up_vector, right_vector)
                            
                            if isScalene: #Rotate second triangle
                                euler_angles = forward_up_right_to_euler(-forward_vector, -up_vector, right_vector)
                                new_triangle2.rotation_euler = euler_angles
                                
                                #DEBUG - Uncomment to draw edges that match input vectors
                                #displayVectors(-forward_vector, -up_vector, right_vector, pivotVertex + initialLocation)
                                
                            #DEBUG - Uncomment to draw edges that match input vectors
                            #displayVectors(forward_vector, up_vector, right_vector, pivotVertex + initialLocation)
                            
                            if isRight:
                                # Rotate Z axis by 90 degrees to fix bad offset
                                rotation_quaternion = new_triangle.rotation_euler.to_quaternion() @ mathutils.Quaternion((0, 0, 1), (math.pi / 2))
                                new_triangle.rotation_euler = rotation_quaternion.to_euler()
                                
                        #============================================= Postion and Scale ================================================
                    
                        # Place the triangle at pivot vertex and set scale

                        #Initializers, values are unimportant
                        desired_dimension = (scaleX, scaleY, minimumThickness)
                        new_x = 1
                        new_y = 1
                        new_z = minimumThickness
                        
                        #Constrain minimum scale to 0.5 forge units
                        if scaleX < minimumThickness: scaleX = minimumThickness
                        if t2scaleX < minimumThickness: t2scaleX = minimumThickness
                        if scaleY < minimumThickness: scaleY = minimumThickness
                        
                        #Adjust for default object sizes and translate object
                        pivotVertex = pivotVertex + initialLocation #Since the mesh was moved to (0,0,0) this shifts the triangle back to the starting location
                        
                        if isRight or isScalene:
                            new_triangle.location = pivotVertex
                            desired_dimension = (scaleX, scaleY, minimumThickness)
                            new_x = desired_dimension[0] / (16 / scaleFactor)#Default size of Floor Angled Standard A
                            new_y = desired_dimension[1] / (16 / scaleFactor)
                            new_z = desired_dimension[2] / (4 / scaleFactor)
                            
                            if isScalene:
                                new_triangle2.location = pivotVertex - (world_normal * minimumThickness)
                                desired_dimension2 = (t2scaleX, scaleY, minimumThickness)
                                new_x2 = desired_dimension2[0] / (16 / scaleFactor)
                            
                        if isIsosceles or isEquallateral:
                            if isEquallateral:
                                new_triangle.location = pivotVertex
                                desired_dimension = ((scaleX * 0.866025), scaleX, minimumThickness)
                                
                            elif isIsosceles:
                                new_triangle.location = pivotVertex
                                desired_dimension = (scaleX, scaleY, minimumThickness)
                                
                            new_x = desired_dimension[0] / (1.73205 / scaleFactor) #Default size of Primitive Triangle
                            new_y = desired_dimension[1] / (2 / scaleFactor)
                            new_z = desired_dimension[2] / (2 / scaleFactor)
                            
        
                        # Set the object's scale to the new dimensions
                        if isRight or isEquallateral or isIsosceles or isScalene:
                            new_triangle.scale.x = new_x / new_triangle.scale.x
                            new_triangle.scale.y = new_y / new_triangle.scale.y
                            new_triangle.scale.z = new_z / new_triangle.scale.z
                            
                            #Link new triangle to the current collection
                            sub_collection.objects.link(new_triangle)
                            
                            if isScalene:
                                new_triangle2.scale.x = new_x2 / new_triangle2.scale.x
                                new_triangle2.scale.y = new_y / new_triangle2.scale.y
                                new_triangle2.scale.z = new_z / new_triangle2.scale.z
                                #Link second triangle
                                sub_collection.objects.link(new_triangle2)
                                
                    else:
                        hasNgons = True
        else:
            errorMsg = "Selected object is not a mesh."
            print(errorMsg)
            raise Exception(errorMsg)
    else:
        errorMsg = "No object selected."
        print(errorMsg)
        raise Exception(errorMsg)
        
    print("============================================================")

    if hasNgons:
        errorMsg = "Mesh is not triangulated, only triangles are supported"
        print(errorMsg)
        raise Exception(errorMsg)
        
    if hasMoved:
        #Move mesh back to it's initial location
        obj.location = initialLocation
        bpy.context.view_layer.update()
        
    print("Generated " + str(objectCount) + " triangles")
    return 0

#================================ Ui setup and logic ==============================

class MainPanel(bpy.types.Panel):
    bl_label = "Mesh to Forge Objects"
    bl_idname = "VIEW_PT_MainPanel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Mesh to Forge Objects'
    
    def draw(self, context):
        layout = self.layout
        layout.scale_y = 1.5
        
        row = layout.row()
        row = layout.row()
        row.operator("wm.myop" ,icon= 'CUBE', text= "Generate Objects")
        layout.prop(context.scene, "subgroupsToggle")
        
    bpy.types.Scene.subgroupsToggle = bpy.props.BoolProperty(
    name="Generate Subgroups",
    default=True
    )
 
class WM_OT_myOp(bpy.types.Operator):
    """Do the thing"""
    bl_label = "More than 5000 faces detected, continue?"
    bl_idname = "wm.myop"  
    confirmation: bpy.props.BoolProperty(name="Yes, lag my pc")
    
    def execute(self, context):
        if self.confirmation:
            print("Yes")
            meshToObjects()
        else:
            print("No")
            
        self.confirmation = False
        return {'FINISHED'}
    
    def invoke(self, context, event):
        # Get the selected object
        selected_object = bpy.context.active_object
        num_faces = 0

        # Check if the selected object is a mesh
        if selected_object and selected_object.type == 'MESH':
            # Get the number of faces
            num_faces = len(selected_object.data.polygons)
            
        if num_faces > 5000:
            return context.window_manager.invoke_props_dialog(self)
        else:
            meshToObjects()
            return {'FINISHED'}
               
def register():
    bpy.utils.register_class(MainPanel)
    bpy.utils.register_class(WM_OT_myOp)
 
def unregister():
    bpy.utils.unregister_class(MainPanel)
    bpy.utils.unregister_class(WM_OT_myOp)
    del bpy.types.Scene.subgroupsToggle
    
    # This is required in order for the script to run in the text editor    
if __name__ == "__main__":
    register()