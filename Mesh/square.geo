// Adapté en fonction du raffinement de maillage souhaité
h = 0.04;

// Points du carré
Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {0.0, 1.0, 0.0, h};
Point(3) = {1.0, 0.0, 0.0, h};
Point(4) = {1.0, 1.0, 0.0, h};

// Lignes du carré
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};

// Boucle de lignes pour le carré
Line Loop(5) = {1, 2, 3, 4};

// Surface du carré
Plane Surface(1) = {5};

// Surface physique pour le carré
Physical Surface(1) = {1};

// Lignes physiques pour le carré
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

// Nouveaux points pour les cercles
Point(5) = {0.5, 0.5, 0.0, h};  // Centre du premier cercle
Point(6) = {0.75, 0.5, 0.0, h}; // Centre du deuxième cercle

// Cercles
Circle(1) = {5, 1, 6};  // Premier cercle
Circle(2) = {6, 2, 5};  // Deuxième cercle

// Boucle de lignes pour les cercles
Line Loop(6) = {7, 8};  // 7 pour le premier cercle, 8 pour le deuxième cercle

// Surface combinant le carré et les cercles
Surface Loop(7) = {1, 6};

// Surface physique pour les cercles
Physical Surface(2) = {7};

// Lignes physiques pour les cercles
Physical Line(5) = {7};  // Ligne associée au premier cercle
Physical Line(6) = {8};  // Ligne associée au deuxième cercle

