lc = .5;
L = 1;
H = .005;

Point(1) = {0,0,0,lc};
Point(2) = {L,0,0,lc};
Point(3) = {L,L,0,lc};
Point(4) = {0,L,0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};
Plane Surface(10) = {1};
Transfinite Surface(10);
Recombine Surface(10);
out[] = Extrude{0,0,H}{Surface{10}; Layers{1}; Recombine;};

Transfinite Surface(out[0]);
Transfinite Surface(out[2]);
Transfinite Surface(out[3]);
Transfinite Surface(out[4]);
Transfinite Surface(out[5]);
Recombine Surface(out[0]);
Recombine Surface(out[2]);
Recombine Surface(out[3]);
Recombine Surface(out[4]);
Recombine Surface(out[5]);

Transfinite Volume(out[1]);
Recombine Volume(out[1]);
Physical Volume(1) = {out[1]};