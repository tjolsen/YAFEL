lc = 0.02;

L = 1;
H = 0.2;

Point(1) = {0,0,0,lc};
Point(2) = {L,0,0,lc};
Point(3) = {L,H,0,lc};
Point(4) = {0,H,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Surface(1);
Recombine Surface(1);

Physical Surface(1) = {1};
Physical Line(100) = {1,3};
Physical Line(101) = {4};
Physical Line(102) = {2};