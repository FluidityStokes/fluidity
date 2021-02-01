l=20;
p=0.1;
Point(1) = {0, 0, 0};
Extrude {1, 0, 0} {
  Point{1}; 
}
Extrude {0, 0.5, 0} {
  Line{1};
}
Extrude {0, 0.5, 0} {
  Line{2};
}
Transfinite Line{1} = l+1 Using Bump p;
Transfinite Line{2} = l+1 Using Bump p;
Transfinite Line{6} = l+1 Using Bump p;
Transfinite Line{3} = 0.5*l+1 Using Bump p;
Transfinite Line{4} = 0.5*l+1 Using Bump p;
Transfinite Line{7} = 0.5*l+1 Using Bump p;
Transfinite Line{8} = 0.5*l+1 Using Bump p;
Transfinite Surface{5} = {1,2,3,4} Right;
Transfinite Surface{9} = {3,4,5,6} Right;
// Bottom
Physical Line(4) = {1};
// Right
Physical Line(2) = {4,8};
// Top
Physical Line(3) = {6};
// Left
Physical Line(1) = {3,7};
Physical Surface(0) = {5,9};
