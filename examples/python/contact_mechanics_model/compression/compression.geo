nb_points = 11;
nb_points_height = 11;
height = 2.0;
length = 1.0;

Point(1) = {-length/2,  -height/2, 0, 1};
Point(2) = {length/2, -height/2, 0, 1};
Point(3) = {length/2, 0, 0, 1};
Point(4) = {length/2, height/2, 0, 1};
Point(5) = {-length/2, height/2, 0, 1};
Point(6) = {-length/2, 0, 0, 1};
Point(7) = {length/2, 0, 0, 1};
Point(8) = {-length/2, 0, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {7, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {8, 1};
Line(8) = {3, 8};

Line Loop(1) = {1, 2, 8, 7};
Line Loop(2) = {3, 4, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line ("loading") = {4};
Physical Line ("fixed") = {1};

Physical Line ("sides") = {2, 7, 3, 5};

Physical Line ("contact") = {6, 8};

Physical Surface("lower") = {1};
Physical Surface("upper") = {2};

Transfinite Line {6,4} = nb_points Using Progression 1;
Transfinite Line {3,5} = nb_points_height Using Progression 1;

Transfinite Line {1, 8} = nb_points Using Progression 1;
Transfinite Line {7, 2} = nb_points_height Using Progression 1;

Transfinite Surface "*";
Recombine Surface "*";
