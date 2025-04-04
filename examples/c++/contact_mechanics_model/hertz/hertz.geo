cl1 = 0.0075;
cl2 = 0.005;
cl3 = 0.005;
Dy = 0.0;
radius = 0.1;
y = 0.1;
epsilon = 1e-3;

nb_points = 51;

Point(1) = {0, y, 0, 0.25*cl1};
Point(2) = {radius, radius + y, 0, cl1};

Point(4) = {0.2, y, 0, 0.5*cl1};
Point(5) = {0, y, 0, 0.5*cl1};
Point(6) = {0, 0.025, 0, 10*cl3};
Point(7) = {0.2, 0.025, 0, 10*cl3};
Point(8) = {0, radius + y, 0, cl1};

Circle(2) = {1, 8, 2};

Line(3) = {2, 8};
Line(13) = {8, 1};
Line(4) = {6, 7};
Line(5) = {7, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};

Line Loop(9) = {2, 3, 13};
Plane Surface(9) = {9};

Line Loop(11) = {6, 7, 4, 5};
Plane Surface(11) = {11};

Physical Line("contact_bottom") = {6};
Physical Line("contact_top") = {2};

Physical Line("loading") = {3};
Physical Line("fixed") = {4};

Physical Surface("upper") = {9};
Physical Surface("lower") = {11};

Physical Line("sides") = {5};

Physical Line("symmetry") = {13, 7};
Recombine Surface "*";
