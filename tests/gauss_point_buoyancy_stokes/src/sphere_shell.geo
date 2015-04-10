lr=4; // number of cells in the radial direction
Point(1) = {0, 0, 2.7, 1.0};
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{2};
}
Extrude {{-1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{1};
}
Extrude {{-1, 0, 0}, {0, 0, 0}, Pi/4} {
  Point{5};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{4, 3, 1, 2};
}
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{8, 5, 12, 16};
}
Extrude {{0, -1, 0}, {0, 0, 0}, Pi/4} {
  Line{4, 3, 1, 2};
}
Extrude {{0, -1, 0}, {0, 0, 0}, Pi/4} {
  Line{33, 36, 40, 44};
}
Transfinite Line {47} = lr Using Progression 1;
Transfinite Line {33} = lr Using Progression 1;
Transfinite Line {4}  = lr Using Progression 1;
Transfinite Line {5}  = lr Using Progression 1;
Transfinite Line {23} = lr Using Progression 1;
Transfinite Line {50} = lr Using Progression 1;
Transfinite Line {36} = lr Using Progression 1;
Transfinite Line {3}  = lr Using Progression 1;
Transfinite Line {8}  = lr Using Progression 1;
Transfinite Line {19} = lr Using Progression 1;
Transfinite Line {54} = lr Using Progression 1;
Transfinite Line {40} = lr Using Progression 1;
Transfinite Line {1}  = lr Using Progression 1;
Transfinite Line {12} = lr Using Progression 1;
Transfinite Line {26} = lr Using Progression 1;
Transfinite Line {58} = lr Using Progression 1;
Transfinite Line {44} = lr Using Progression 1;
Transfinite Line {2}  = lr Using Progression 1;
Transfinite Line {16} = lr Using Progression 1;
Transfinite Line {30} = lr Using Progression 1;
Transfinite Line {48} = lr Using Progression 1;
Transfinite Line {34} = lr Using Progression 1;
Transfinite Line {6}  = lr Using Progression 1;
Transfinite Line {21} = lr Using Progression 1;
Transfinite Line {51} = lr Using Progression 1;
Transfinite Line {37} = lr Using Progression 1;
Transfinite Line {9}  = lr Using Progression 1;
Transfinite Line {20} = lr Using Progression 1;
Transfinite Line {56} = lr Using Progression 1;
Transfinite Line {42} = lr Using Progression 1;
Transfinite Line {14} = lr Using Progression 1;
Transfinite Line {28} = lr Using Progression 1;


Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{23, 19, 26, 30};
}
Transfinite Line {61} = lr Using Progression 1;
Transfinite Line {64} = lr Using Progression 1;
Transfinite Line {68} = lr Using Progression 1;
Transfinite Line {72} = lr Using Progression 1;
Transfinite Line {62} = lr Using Progression 1;
Transfinite Line {65} = lr Using Progression 1;
Transfinite Line {70} = lr Using Progression 1;
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{61, 64, 68, 72};
}
Transfinite Line {75} = lr Using Progression 1;
Transfinite Line {78} = lr Using Progression 1;
Transfinite Line {82} = lr Using Progression 1;
Transfinite Line {86} = lr Using Progression 1;
Transfinite Line {76} = lr Using Progression 1;
Transfinite Line {79} = lr Using Progression 1;
Transfinite Line {84} = lr Using Progression 1;
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{75, 78, 82, 86};
}
Transfinite Line {89}  = lr Using Progression 1;
Transfinite Line {92}  = lr Using Progression 1;
Transfinite Line {96}  = lr Using Progression 1;
Transfinite Line {100} = lr Using Progression 1;
Transfinite Line {90}  = lr Using Progression 1;
Transfinite Line {93}  = lr Using Progression 1;
Transfinite Line {98}  = lr Using Progression 1;
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Line{89, 92, 96, 100};
}
Transfinite Line {104} = lr Using Progression 1;
Transfinite Line {107} = lr Using Progression 1;
Transfinite Line {112} = lr Using Progression 1;

// top
Transfinite Surface {63} Alternated;
Transfinite Surface {77} Alternated;
Transfinite Surface {91} Alternated;
Transfinite Surface {105} Alternated;
Transfinite Surface {49} Alternated;
Transfinite Surface {35} Alternated;
Transfinite Surface {7} Alternated;
Transfinite Surface {25} Alternated;
// upper middle
Transfinite Surface {67} Alternated;
Transfinite Surface {81} Alternated;
Transfinite Surface {95} Alternated;
Transfinite Surface {109} Alternated;
Transfinite Surface {53} Alternated;
Transfinite Surface {39} Alternated;
Transfinite Surface {11} Alternated;
Transfinite Surface {22} Alternated;
// lower middle
Transfinite Surface {71} Alternated;
Transfinite Surface {85} Alternated;
Transfinite Surface {99} Alternated;
Transfinite Surface {113} Alternated;
Transfinite Surface {57} Alternated;
Transfinite Surface {43} Alternated;
Transfinite Surface {15} Alternated;
Transfinite Surface {29} Alternated;
// bottom
Transfinite Surface {74} Alternated;
Transfinite Surface {88} Alternated;
Transfinite Surface {102} Alternated;
Transfinite Surface {116} Alternated;
Transfinite Surface {60} Alternated;
Transfinite Surface {46} Alternated;
Transfinite Surface {18} Alternated;
Transfinite Surface {32} Alternated;

// top
Physical Surface(33) = {49};
Physical Surface(34) = {35};
Physical Surface(35) = {7};
Physical Surface(36) = {25};
Physical Surface(133) = {63};
Physical Surface(134) = {77};
Physical Surface(135) = {91};
Physical Surface(136) = {105};
// upper middle
Physical Surface(25) = {53};
Physical Surface(26) = {39};
Physical Surface(27) = {11};
Physical Surface(28) = {22};
Physical Surface(125) = {67};
Physical Surface(126) = {81};
Physical Surface(127) = {95};
Physical Surface(128) = {109};
// lower middle
Physical Surface(29) = {57};
Physical Surface(30) = {43};
Physical Surface(31) = {15};
Physical Surface(32) = {29};
Physical Surface(129) = {71};
Physical Surface(130) = {85};
Physical Surface(131) = {99};
Physical Surface(132) = {113};
// bottom
Physical Surface(37) = {60};
Physical Surface(38) = {46};
Physical Surface(39) = {18};
Physical Surface(40) = {32};
Physical Surface(137) = {74};
Physical Surface(138) = {88};
Physical Surface(139) = {102};
Physical Surface(140) = {116};
