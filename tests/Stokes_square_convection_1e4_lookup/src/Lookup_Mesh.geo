npoints = 100 ;
Point(1) = {0.0, 0.0, 0.0, 100e3};
Extrude {0.0, -1000e3, 0} {
  Point{1}; Layers{npoints};
}
Extrude {2890e3, 0.0, 0.0} {
  Line{1}; Layers{npoints};
}
Physical Line(6) = {1};
Physical Line(7) = {4};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Surface(10) = {5};
