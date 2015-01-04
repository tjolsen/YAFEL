lc = 100;

Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {2,0,0,lc};
Point(4) = {3,0,0,lc};
Point(5) = {1,1,0,lc};
Point(6) = {2,1,0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {1,5};
Line(5) = {4,6};
Line(6) = {2,5};
Line(7) = {3,6};
Line(8) = {2,6};
Line(9) = {3,5};
//Line(10) = {5,6};

Physical Point(1) = {1};
Physical Point(2) = {4};
Physical Line(3) = {1,2,3,4,5,6,7,8,9};//,10};