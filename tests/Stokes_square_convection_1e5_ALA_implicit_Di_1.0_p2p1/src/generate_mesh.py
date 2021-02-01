#!/usr/bin/python

import fluidity.diagnostics.annulus_mesh as annulus_mesh
import fluidity.diagnostics.gmshtools as gmshtools

coordsY = annulus_mesh.SliceCoordsLinear(0.0, 1.0, 0.01, 30)
coordsX = annulus_mesh.SliceCoordsLinear(0.0, 1.0, 0.01, 30)
mesh = annulus_mesh.GenerateRectangleMesh(coordsX, coordsY)
gmshtools.WriteMsh(mesh, "square-structured-linear.msh")

