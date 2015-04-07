FLFLAGS = -v2 -l
FLMLMODEL = stokes-1e4-p2p1-ALA-implicit-Di-0.25
FLUIDITY = fluidity
MESH = square-structured-linear

default: input

input: clean
	./src/generate_mesh.py

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).msh
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*_checkpoint.edge $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.node $(FLMLMODEL)_?*_checkpoint.flml $(FLMLMODEL)_?*_?*_checkpoint.vtu $(FLMLMODEL).stat $(FLMLMODEL).detectors fluidity.err-0 fluidity.log-0 matrixdump matrixdump.info
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS)
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS)
