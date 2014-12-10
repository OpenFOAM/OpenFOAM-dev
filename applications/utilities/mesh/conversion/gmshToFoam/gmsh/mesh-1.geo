// Parametres geometriques
r1 = 200*Cos(Pi/4)/1000;
r2 = 200*Cos(Pi/4)/1000;

h1 = 250/1000;
h2 = 360/1000;
h3 = 900/1000;
h4 = 1900/1000;

// Parametres de maillage
// selon le rayon
rCells = 10/2; 	rRatio = 0.85;
// selon S1
S1Cells = 30/2; 	S1ratio = 1;
// selon S2
S2Cells = 35/2;	S2ratio = 0.95;
// selon S3
S3Cells = 20/2;	S3ratio = 1;

Point(1) = {r1, r1, h4};
Point(2) = {r1, r1, h3};
Point(3) = {r2, r2, h2};
Point(4) = {r2, r2, h1};
Point(5) = {0, 0, h1};
Point(6) = {0, 0, h2};
Point(7) = {0, 0, h3};
Point(8) = {0, 0, h4};

Line(1) = {8, 1};
Line(2) = {1, 2};
Line(3) = {7, 2};
Line(4) = {8, 7};
Line(5) = {2, 3};
Line(6) = {6, 3};
Line(7) = {7, 6};
Line(8) = {3, 4};
Line(9) = {5, 4};
Line(10) = {6, 5};

Line Loop(11) = {1, 2, -3, -4};
Ruled Surface(12) = {11};
Line Loop(13) = {5, -6, -7, 3};
Ruled Surface(14) = {13};
Line Loop(15) = {8, -9, -10, 6};
Ruled Surface(16) = {15};
Transfinite Line {1, 3, 6, 9} = rCells Using Progression rRatio;
Transfinite Line {4, 2} = S1Cells Using Progression S1ratio;
Transfinite Line {7, 5} = S2Cells Using Progression S2ratio;
Transfinite Line {10, 8} = S3Cells Using Progression S3ratio;
Transfinite Surface {12} = {8, 1, 2, 7};
Transfinite Surface {14} = {7, 2, 3, 6};
Transfinite Surface {16} = {6, 3, 4, 5};
Recombine Surface {12, 14, 16};
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{12, 14, 16};
  Layers{25};
  Recombine;
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{33, 50, 67};
  Layers{25};
  Recombine;
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{84, 101, 118};
  Layers{25};
  Recombine;
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Surface{135, 152, 169};
  Layers{25};
  Recombine;
}
Physical Surface("entree") = {126, 75, 24, 177};
Physical Surface("S1") = {28, 181, 130, 79};
Physical Surface("S2") = {93, 42, 193, 144};
Physical Surface("S3") = {110, 59, 205, 161};
Physical Surface("fond") = {113, 62, 208, 164};
Physical Volume("fluide") = {4, 7, 10, 1, 5, 8, 11, 2, 9, 12, 3, 6};
