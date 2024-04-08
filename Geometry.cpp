#ifndef _GEOMETRY_CPP

#include<iostream>
#include "Geometry.h"

Geometry::Geometry(Mesh* mesh) : _vertices(mesh->getVertices()),
                                  _triangles(mesh->getTriangles()), _edgesNeumann(mesh->getNeumannEdges())
{
}

// Construction des 3 fonctions de base (phihat0, phihat1, phihat2)
double Geometry::phihat(int hati, Eigen::Vector2d hatX)
{
    double value;
    switch (hati)
    {
    case 0:
        value = 1 - hatX(0) - hatX(1);
        break;
    case 1:
        value = hatX(0);
        break;
    case 2:
        value = hatX(1);
        break;
    default:
        std::cout << "Veuillez choisir un numéro valide pour hati dans le triangle de référence (0, 1 ou 2)" << std::endl;
        abort();
    }
    return value;
}

// Construction des 3 gradients des fonctions de base (gradphihat0, gradphihat1, gradphihat2)
// (indépendant du vecteur hatX)
Eigen::Vector2d Geometry::gradphihat(int hati)
{
    Eigen::Vector2d gradphi;
    switch (hati)
    {
    case 0:
        gradphi(0) = -1;
        gradphi(1) = -1;
        break;
    case 1:
        gradphi(0) = 1;
        gradphi(1) = 0;
        break;
    case 2:
        gradphi(0) = 0;
        gradphi(1) = 1;
        break;
    default:
        std::cout << "Veuillez choisir un numéro valide pour hati dans le triangle de référence (0, 1 ou 2)" << std::endl;
        abort();
    }
    return gradphi;
}

// Construction de la fonction de transformation du triangle de référence alphahat en alpha
Eigen::Vector2d Geometry::Falpha(int refTrianglealpha, Eigen::Vector2d hatX)
{
    int n1(_triangles(refTrianglealpha, 0)), n2(_triangles(refTrianglealpha, 1)), n3(_triangles(refTrianglealpha, 2));
    double x1(_vertices(n1, 0)), y1(_vertices(n1, 1));
    double x2(_vertices(n2, 0)), y2(_vertices(n2, 1));
    double x3(_vertices(n3, 0)), y3(_vertices(n3, 1));

    Eigen::Vector2d X;
    X(0) = x1 * phihat(0, hatX) + x2 * phihat(1, hatX) + x3 * phihat(2, hatX);
    X(1) = y1 * phihat(0, hatX) + y2 * phihat(1, hatX) + y3 * phihat(2, hatX);

    return X;
}

// Construction de la Jacobienne de la fonction de transformation du triangle de référence alphahat en alpha
Eigen::Matrix2d Geometry::JFalpha(int refTrianglealpha)
{
    int n1(_triangles(refTrianglealpha, 0)), n2(_triangles(refTrianglealpha, 1)), n3(_triangles(refTrianglealpha, 2));
    double x1(_vertices(n1, 0)), y1(_vertices(n1, 1));
    double x2(_vertices(n2, 0)), y2(_vertices(n2, 1));
    double x3(_vertices(n3, 0)), y3(_vertices(n3, 1));

    Eigen::Matrix2d JF;
    JF(0, 0) = x2 - x1;
    JF(0, 1) = x3 - x1;
    JF(1, 0) = y2 - y1;
    JF(1, 1) = y3 - y1;

    return JF;
}

// Valeur absolue du déterminant de la Jacobienne
double Geometry::absoludetJ(int refTrianglealpha)
{
    return fabs((JFalpha(refTrianglealpha)).determinant());
}

// Construction de la matrice BoF
Eigen::MatrixXd Geometry::BoF_alpha(int refTrianglealpha)
{
    int n1(_triangles(refTrianglealpha, 0)), n2(_triangles(refTrianglealpha, 1)), n3(_triangles(refTrianglealpha, 2));
    double x1(_vertices(n1, 0)), y1(_vertices(n1, 1));
    double x2(_vertices(n2, 0)), y2(_vertices(n2, 1));
    double x3(_vertices(n3, 0)), y3(_vertices(n3, 1));

    Eigen::MatrixXd BoF(3,6);
    BoF(0, 0) = (y1 - y3) + (x3 - x1); BoF(0, 1) = 0; BoF(0, 2) = y3 - y1; BoF(0, 3) = 0; BoF(0, 4) = x1 - x3; BoF(0, 5) = 0;
    BoF(1, 0) = 0; BoF(1, 1) = (y2 - y1) + (x1 - x2); BoF(1, 2) = 0; BoF(1, 3) = y1 - y2;  BoF(1, 4) = 0; BoF(1, 5) = x2 - x1;
    BoF(2, 0) = BoF(1, 1); BoF(2, 1) = BoF(0, 0); BoF(2, 2) = BoF(1, 3); BoF(2, 3) = BoF(0, 2); BoF(2, 4) = BoF(1, 5); BoF(2, 5) = BoF(0, 4);

    return (1 / JFalpha(refTrianglealpha).determinant()) * BoF;
}

// Construction de la matrice des coefficients d'élasticité
Eigen::MatrixXd Geometry::coef_elastique(double E, double nu)
{
    Eigen::MatrixXd D(3,3);
    D(0, 0) = 1 - nu; D(0, 1) = nu; D(0,2) = 0;
    D(1, 0) = nu; D(1, 1) = 1 - nu; D(1,2) = 0;
    D(2,0)= 0; D(2,1)= 0; D(2,2)= 1-2*nu;

    return (E / ((1 + nu) * (1 - 2 * nu))) * D;
}

// Construction de la matrice élémentaire avant de passer à la quadrature pour le prochain cours
Eigen::MatrixXd Geometry::Elementaire(int refTrianglealpha, double mu, double E)
{
    Eigen::MatrixXd K_elem(3,6);
    K_elem = (BoF_alpha(refTrianglealpha).transpose()) * (coef_elastique(E, mu).transpose()) * (BoF_alpha(refTrianglealpha)) * (absoludetJ(refTrianglealpha));

    return K_elem;


}

// Points et poids de quadrature de la formule du milieu (intégration 2D)
void Geometry::quadraturePointsAndWeightsMidpointFormula(Eigen::VectorXd& weights, Eigen::Matrix<double, Eigen::Dynamic, 2>& points)
{
  weights.resize(3);
  weights(0) = 1./3.; weights(1) = 1./3.; weights(2) = 1./3.;
  points.resize(3,2);
  points(0,0) = 0.5; points(0,1) = 0;
  points(1,0) = 0.0; points(1,1) = 0.5;
  points(2,0) = 0.5; points(2,1) = 0.5;
}

// Points et poids de quadrature de la formule du milieu (intégration 1D)
void Geometry::quadraturePointsAndWeightsSimpsonFormula(Eigen::VectorXd& weights, Eigen::Matrix<double, Eigen::Dynamic, 2>& points, int numedge)
{
  weights.resize(3);
  weights(0) = 1./6.; weights(1) = 2./3.; weights(2) = 1./6.;
  points.resize(3,2);
  if (numedge == 1)
  {
    points(0,0) = 0.0; points(0,1) = 0.0;
    points(1,0) = 0.5; points(1,1) = 0.0;
    points(2,0) = 1.0; points(2,1) = 0.0;
  }
  else if (numedge == 2)
  {
    points(0,0) = 1.0; points(0,1) = 0.0;
    points(1,0) = 0.5; points(1,1) = 0.5;
    points(2,0) = 0.0; points(2,1) = 1.0;
  }
  else if (numedge == 3)
  {
    points(0,0) = 0.0; points(0,1) = 1.0;
    points(1,0) = 0.0; points(1,1) = 0.5;
    points(2,0) = 0.0; points(2,1) = 0.0;
  }
}


// traiter la conditon au bord de Neumann non homogene 
// Construction de la fonction de transformation
// de l'arête de référence hatE en E
Eigen::Vector2d Geometry::FE(int refEdgeE, Eigen::Vector2d hatX)
{
  int n1(_edgesNeumann(refEdgeE,0)), n2(_edgesNeumann(refEdgeE,1));
  double x1(_vertices(n1,0)), y1(_vertices(n1,1));
  double x2(_vertices(n2,0)), y2(_vertices(n2,1));

  Eigen::Vector2d X;
  X(0) = (x2-x1)*hatX(0)-(y2-y1)*hatX(1)+x1;
  X(1) = (y2-y1)*hatX(0)+(x2-x1)*hatX(1)+y1;

  return X;
}

// Construction de la jacobienne de la fonction de transformation
// de l'arête de référence hatE en E
// (indépendant du vecteur hatX)
Eigen::Matrix2d Geometry::JFE(int refEdgeE)
{
  int n1(_edgesNeumann(refEdgeE,0)), n2(_edgesNeumann(refEdgeE,1));
  double x1(_vertices(n1,0)), y1(_vertices(n1,1));
  double x2(_vertices(n2,0)), y2(_vertices(n2,1));

  Eigen::Matrix2d X;
  X(0,0) = (x2-x1); X(0,1) = -(y2-y1); X(1,0) = (y2-y1); X(1,1) = (x2-x1);

  return X;
}

// Construction de la valeur absolue du déterminant de la jacobienne
// de la fonction de transformation
// de l'arête de référence hatE en E
// (indépendant du vecteur hatX)
double Geometry::measE(int refEdgeE)
{
  return sqrt(fabs((JFE(refEdgeE)).determinant()));
}


#define _GEOMETRY_CPP
#endif

