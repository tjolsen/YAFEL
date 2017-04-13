lc = 1;
L = 2;

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
//Transfinite Surface(10);
//Recombine Surface(10);
//Physical Surface(1) = {10};

out[] = Extrude{0,0,L}{Surface{10};};
//Transfinite Surface{out[0],out[2],out[3],out[4],out[5]};
//Recombine Surface(out[0]);
//Recombine Surface(out[2]);
//Recombine Surface(out[3]);
//Recombine Surface(out[4]);
//Recombine Surface(out[5]);

//Transfinite Volume(out[1]);
//Recombine Volume(out[1]);
Physical Volume(1) = {out[1]};
