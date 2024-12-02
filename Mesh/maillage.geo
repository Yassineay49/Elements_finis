h = 400.0;
//+
Point(1) = {40.0, -24.0, 0.0, h};
Point(2) = {964.0, -24.0, 0.0, h};
Point(3) = {964.0, 920.0, 0.0, h};
Point(4) = {40.0, 920.0, 0.0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
x1 = 150.0; y1=80.0; r1 = 100.0;
Point(5) = {x1, y1, 0.0, h/2};
Point(6) = {x1-r1, y1, 0.0, h/4};
Circle(5) = {6,5,6};
//+
x2 = 250.0; y2=312.0; r2 = 125.0;
Point(7) = {x2, y2, 0.0, h/2};
Point(8) = {x2-r2, y2, 0.0, h/4};
Circle(6) = {8,7,8};
//+
x3 = 475.0; y3=175.0; r3 = 125.0;
Point(9) = {x3, y3, 0.0, h/2};
Point(10) = {x3-r3, y3, 0.0, h/4};
Circle(7) = {10,9,10};
//+
x4 = 427.0; y4=501.0; r4 = 120.0;
Point(11) = {x4, y4, 0.0, h/2};
Point(12) = {x4-r4, y4, 0.0, h/4};
Circle(8) = {12,11,12};
//+
x5 = 199.0; y5=631.0; r5 = 130.0;
Point(13) = {x5, y5, 0.0, h/2};
Point(14) = {x5-r5, y5, 0.0, h/4};
Circle(9) = {14,13,14};
//+
x6 = 634.0; y6=373.0; r6 = 115.0;
Point(15) = {x6, y6, 0.0, h/2};
Point(16) = {x6-r6, y6, 0.0, h/4};
Circle(10) = {16,15,16};
//+
x7 = 852.0; y7=168.0; r7 = 110.0;
Point(17) = {x7, y7, 0.0, h/2};
Point(18) = {x7-r7, y7, 0.0, h/4};
Circle(11) = {18,17,18};
//+
x8 = 802.0; y8=564.0; r8 = 105.0;
Point(19) = {x8, y8, 0.0, h/2};
Point(20) = {x8-r8, y8, 0.0, h/4};
Circle(12) = {20,19,20};
//+
x9 = 522.0; y9=739.0; r9 = 110.0;
Point(21) = {x9, y9, 0.0, h/2};
Point(22) = {x9-r9, y9, 0.0, h/4};
Circle(13) = {22,21,22};
//+
x10 = 741.0; y10=801.0; r10 = 95.0;
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
//+
Physical Curve("ext", 1) = {1};
Physical Curve("int", 2) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
