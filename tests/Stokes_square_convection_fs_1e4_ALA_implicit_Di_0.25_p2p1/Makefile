input: clean
	gmsh -2 src/squareA.geo -o squareA.msh
	gmsh -2 src/squareB.geo -o squareB.msh

clean:
	rm -rf *.msh *.stat *.vtu *.log-0 *.err-0 *checkpoint *.convergence *.log *.node *.edge *.ele *.halo *.pvtu stokes-B_* stokes-A_* fluidity.* *.pyc
