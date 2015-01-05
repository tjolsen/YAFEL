lc = 10;

Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {2,0,0,lc};
Line(1) = {1,2};
Line(2) = {2,3};

Physical Point(1) = {1};
Physical Point(2) = {3};
Physical Line(3) = {1,2};