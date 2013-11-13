l=20;
p=1.0;
Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1}; 
}
Extrude {0, 1.0, 0} {
  Line{1};
}
Transfinite Line{1} = l+1;
Transfinite Line{2} = l+1;
Transfinite Line{-3} = l+1 Using Progression p;
Transfinite Line{-4} = l+1 Using Progression p;
Transfinite Surface{5} = {1,2,3,4} Alternate;
// Bottom
Physical Line(4) = {1};
Physical Line(2) = {4};
// Top
Physical Line(3) = {2};
Physical Line(1) = {3};
// Internal
Physical Surface(15) = {5};
