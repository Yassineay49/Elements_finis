h = 1200.0e-8;
//+
Point(1) = {120.0e-8, -2760.0e-8, 0.0, h};
Point(2) = {120.0e-8, 72.0e-8, 0.0, h};
Point(3) = {2892.0e-8, 72.0e-8, 0.0, h};
Point(4) = {2892.0e-8, -2760.0e-8, 0.0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
x1 = 450.0e-8; y1=-240.0e-8; r1 = 300.0e-8;
Point(5) = {x1, y1, 0.0, h/2};
Point(6) = {x1-r1, y1, 0.0, h/4};
Circle(5) = {6,5,6};
//+
x2 = 750.0e-8; y2=-936.0e-8; r2 = 375.0e-8;
Point(7) = {x2, y2, 0.0, h/2};
Point(8) = {x2-r2, y2, 0.0, h/4};
Circle(6) = {8,7,8};
//+
x3 = 1425.0e-8; y3=-525.0e-8; r3 = 375.0e-8;
Point(9) = {x3, y3, 0.0, h/2};
Point(10) = {x3-r3, y3, 0.0, h/4};
Circle(7) = {10,9,10};
//+
x4 = 1281.0e-8; y4=-1503.0e-8; r4 = 360.0e-8;
Point(11) = {x4, y4, 0.0, h/2};
Point(12) = {x4-r4, y4, 0.0, h/4};
Circle(8) = {12,11,12};
//+
x5 = 597.0e-8; y5=-1893.0e-8; r5 = 390.0e-8;
Point(13) = {x5, y5, 0.0, h/2};
Point(14) = {x5-r5, y5, 0.0, h/4};
Circle(9) = {14,13,14};
//+
x6 = 1902.0e-8; y6=-1119.0e-8; r6 = 345.0e-8;
Point(15) = {x6, y6, 0.0, h/2};
Point(16) = {x6-r6, y6, 0.0, h/4};
Circle(10) = {16,15,16};
//+
x7 = 2556.0e-8; y7=-504.0e-8; r7 = 330.0e-8;
Point(17) = {x7, y7, 0.0, h/2};
Point(18) = {x7-r7, y7, 0.0, h/4};
Circle(11) = {18,17,18};
//+
x8 = 2406.0e-8; y8=-1692.0e-8; r8 = 315.0e-8;
Point(19) = {x8, y8, 0.0, h/2};
Point(20) = {x8-r8, y8, 0.0, h/4};
Circle(12) = {20,19,20};
//+
x9 = 1566.0e-8; y9=-2217.0e-8; r9 = 330.0e-8;
Point(21) = {x9, y9, 0.0, h/2};
Point(22) = {x9-r9, y9, 0.0, h/4};
Circle(13) = {22,21,22};
//+
x10 =2223.0e-8; y10=-2403.0e-8; r10 = 285.0e-8;
Point(23) = {x10, y10, 0.0, h/2};
Point(24) = {x10-r10, y10, 0.0, h/4};
Circle(14) = {24,23,24};
//+
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5};
Curve Loop(3) = {6};
Curve Loop(4) = {7};
Curve Loop(5) = {8};
Curve Loop(6) = {9};
Curve Loop(7) = {10};
Curve Loop(8) = {11};
Curve Loop(9) = {12};
Curve Loop(10) = {13};
Curve Loop(11) = {14};
//+
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; 
Plane Surface(2) = {2}; 
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5}; 
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8}; 
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
