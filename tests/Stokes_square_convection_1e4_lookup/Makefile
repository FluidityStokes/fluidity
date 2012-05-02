MODEL=Mantle_Lookup_Test

default: input

input:  clean
	gmsh -2 -optimize src/Lookup_Mesh.geo

clean:
	rm -f  *.ele *.edge *.node *.vtu *.stat *.msh *.detectors fluidity.* matrixdump* fluidity.log* fluidity.err* nohup.* Parallel*
	rm -R -f $(MODEL)_*
