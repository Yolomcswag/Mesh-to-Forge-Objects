# Mesh-to-Forge-Objects
An addon for this project: https://github.com/TubbyMcFatDuck/Halo-Infinite-Blender-2-Forge-Printer

A simple addon to convert triangulated meshes into forge triangles. Install as a blender addon and  a side button should appear in blender.
In order for the addon to function it needs references to 'Floor Angled Standard A' and 'Primitive Triangle', to make these simple place both objects in the scene somewhere. (they can have visibility hidden the program just needs a reference) 
The 'Generate Subgroups' option creates a collection for every ~150 triangles so large prints can be easily divided into prefabs.
If all of the objects made with your mesh appear to be off try flipping the normals on your mesh, the normals control which direction is considered 'Up' by the addon so if they are wrong way every triangle will be placed on the wrong side of the faces.
