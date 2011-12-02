#!/usr/bin/env python

import fluidity.diagnostics.annulus_mesh as annulus_mesh
import fluidity.diagnostics.triangletools as triangletools

coordsY = annulus_mesh.SliceCoordsLinear(0.0, 2000e3, 5e3, 70)
coordsX = annulus_mesh.SliceCoordsLinear(0.0, 2890e3, 5e3, 70)
mesh = annulus_mesh.GenerateRectangleMesh(coordsX, coordsY)
triangletools.WriteTriangle(mesh, "square-structured-linear")

