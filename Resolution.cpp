#ifndef _RESOLUTION_H

#include "Resolution.h"
#include <fstream>
#include <iostream>

Resolution::Resolution(DataFile* dataFile, Assemblage* assemblage, Solver* solver, Mesh* mesh):
_solver(solver), _assemblage(assemblage), _vertices(mesh->getVertices()), _triangles(mesh->getTriangles())
{
}

Resolution:: ~Resolution() {}

void Resolution::saveSolution()
{
    saveSolution(_sol);
}

void Resolution::saveSolution(Eigen::SparseVector<double> sol)
{
    std::string name_file = "Resultats.vtk";
    int nb_vert = _vertices.rows();

  	std::ofstream solution;
  	solution.open(name_file, std::ios::out);
  	solution.precision(7);

    solution << "# vtk DataFile Version 3.0 " << std::endl;
    solution << "2D Unstructured Grid" << std::endl;
    solution << "ASCII" << std::endl;
    solution << "DATASET UNSTRUCTURED_GRID" << std::endl;

    solution << "POINTS " << nb_vert << " float " << std::endl;
    for (int i = 0 ; i < nb_vert ; i++)
    {
      solution << _vertices(i,0) << " " << _vertices(i,1) << " " << "0." << std::endl;
    }

    solution << "CELLS " << _triangles.rows() << " " << _triangles.rows()*4 << std::endl;
    for (int i = 0 ; i < _triangles.rows() ; i++)
    {
      solution << 3 << " " << _triangles(i,0) << " " << _triangles(i,1) << " " << _triangles(i,2) << std::endl;
    }

    solution << "CELL_TYPES " << _triangles.rows() << std::endl;
  	for (int i=0; i< _triangles.rows(); i++)
  	{
  		solution << 5 << std::endl;
  	}

    solution << "POINT_DATA " << nb_vert << std::endl; // Car solution aux sommets
    solution << "SCALARS sol float 1" << std::endl;
    solution << "LOOKUP_TABLE default" << std::endl;
    double eps = 1.0e-10;
  	for (int i = 0 ; i < nb_vert ; i++)
    {
  		solution << fmax(eps,sol.coeffRef(i)) << std::endl;
    }
    solution << std::endl;

  	solution.close();
}

void Resolution::systeme(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice)
{
    // Construction de la matrice de rigidité
    std::cout << "---- En cours d'assemblage -----"<<std::endl;
    _assemblage->matrice_globalK(E_fibre, E_matrice, nu_fibre, nu_matrice);
    Eigen::SparseMatrix<double, Eigen::RowMajor> systemMatrix(_assemblage->getStiffnessMatrix());

    // Applications des conditions aux bords sur la matrice
    _assemblage->applyBCToSystemMatrix(systemMatrix);

    // pour donner la matrice au solveur 
    std::cout<<"----- Application du solveur ----------"<<std::endl;
    _solver->setSystemMatrix(systemMatrix);

    // construction du second membre 
    _assemblage->assemblesourceNeumann();
    Eigen::SparseVector<double> RHS(_assemblage->getSourceAndNeumann());

    // mise a jour de la condition sur le terme de droite
    _assemblage->applyBCToRHS(RHS);
    // Resolution du système
    std::cout << "-------En cours de resolution--------"<<std::endl;
    _sol=_solver->solve(RHS);
    std::cout << "-------Termine--------"<<std::endl;
}


void Resolution::contrainte( double E_fibre, double E_matrice, double nu_fibre, double nu_matrice)
{
  //const double delta_x = 1e-6; // difference fini
   _contrainte.resize(_triangles.rows(), 3); // Redimensionner la matrice des contraintes
        
  // boucle sur les elements 
   for (int K = 0; K < _triangles.rows() ; K++){

    int n1(_triangles(K, 0)), n2(_triangles(K, 1)), n3(_triangles(K, 2));
    double x1(_vertices(n1, 0)), y1(_vertices(n1, 1));
    double x2(_vertices(n2, 0)), y2(_vertices(n2, 1));
    double x3(_vertices(n3, 0)), y3(_vertices(n3, 1));

    double u1 = _sol.coeff(n1);
    double u2 = _sol.coeff(n2);
    double u3 = _sol.coeff(n3);

    // calcul des gradients de deplacememnt (on le fait par approximation difference finie)
    // on note u la composante suivant x
    // et v la composante suivant y
    double du_dx = (u2-u1)/(x2-x1);
    double du_dy = (u2-u1)/(y2-y1);
    double dv_dx = (u3-u1)/(x3-x1);
    double dv_dy = (u3-u1)/(y3-y1);

    // calcul des composantes du tenseur des deformations 
    double epsilon_x = du_dx;
    double epsilon_y = du_dy;
    double gamma_xy = 0.5 * (du_dy + dv_dx);


    double sigma_x(0.0);
    double sigma_y(0.0);
    double tau_xy(0.0);
    // les contraintes 
    if(_triangles(K,3)=100){
        //K_elem=_geometry->Elementaire(K, nu_matrice, E_matrice);
        double a(E_matrice/((1-nu_matrice)*(1-2*nu_matrice)));
        sigma_x = (1-nu_matrice)*a*epsilon_x + nu_matrice*a*epsilon_y;
        sigma_x = a*nu_matrice*epsilon_x + (1-nu_matrice)*a*epsilon_y;
        tau_xy = 2*(E_matrice/(1+nu_matrice))*gamma_xy;
      }
      else{
        double a_fibre = E_fibre / ((1 - nu_fibre) * (1 - 2 * nu_fibre));
        double sigma_x = (1 - nu_fibre) * a_fibre * epsilon_x + nu_fibre * a_fibre * epsilon_y;
        double sigma_y = a_fibre * nu_fibre * epsilon_x + (1 - nu_fibre) * a_fibre * epsilon_y;
        double tau_xy = 2 * (E_fibre / (1 + nu_fibre)) * gamma_xy;
      }

    _contrainte.coeffRef(K, 0) = sigma_x; // Valeur de sigma_x
    _contrainte.coeffRef(K, 1) = sigma_y; // Valeur de sigma_y
    _contrainte.coeffRef(K, 2) = tau_xy;  // Valeur de tau_xy

  }
}

void Resolution::afficherContrainte()
{
  // Parcourir chaque élément du maillage
    for (int i = 0; i < _triangles.rows(); ++i) {
        // Afficher ou traiter les contraintes pour chaque élément
        double sigma_x = _contrainte.coeff(i, 0); // Valeur de sigma_x
        double sigma_y = _contrainte.coeff(i, 1); // Valeur de sigma_y
        double tau_xy = _contrainte.coeff(i, 2);  // Valeur de tau_xy
        
        // Afficher les valeurs de contrainte pour chaque élément
        std::cout << "Contraintes pour l'élément " << i << ":" << std::endl;
        std::cout << "Sigma_x = " << sigma_x << ", Sigma_y = " << sigma_y << ", Tau_xy = " << tau_xy << std::endl;
    }
}
#define _RESOLUTION_H
#endif