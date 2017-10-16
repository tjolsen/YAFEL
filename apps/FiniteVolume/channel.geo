H = 0.1;
L = 1.0;
lc = 0.001;


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

//Physical Line(1) = {1,3}; // No-Slip walls
//Physical Line(2) = {2};   // Outlet
//Physical Line(3) = {4};   // Inlet

Physical Surface(10) = {1}; //Volume