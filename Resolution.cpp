#ifndef _RESOLUTION_H

#include "Resolution.h"
#include <fstream>
#include <iostream>
#include <cmath>

Resolution::Resolution(DataFile* dataFile, Assemblage* assemblage, Solver* solver, Mesh* mesh):
_solver(solver), _assemblage(assemblage), _vertices(mesh->getVertices()), _triangles(mesh->getTriangles()), _datafile(dataFile)
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


// test sur tout le composé
void Resolution::contrainte(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice) {
    // const double delta_x = 1e-6; // différence finie
    _contrainte_Matrice.resize(_triangles.rows(), 3); // Redimensionner la matrice des contraintes
    _contrainte_Fibre.resize(_triangles.rows(), 3);

    // boucle sur les éléments
    for (int K = 0; K < _triangles.rows(); K++) {

        int n1 = _triangles(K, 0), n2 = _triangles(K, 1), n3 = _triangles(K, 2);
        // double x1(_vertices(n1, 0)), y1(_vertices(n1, 1));
        // double x2(_vertices(n2, 0)), y2(_vertices(n2, 1));
        // double x3(_vertices(n3, 0)), y3(_vertices(n3, 1));

        Eigen::MatrixXd B = _assemblage->_geometry->BoF_alpha(K);
        double u1x = _sol.coeff(n1);
        double u1y = _sol.coeff(n1+1);
        double u2x = _sol.coeff(n2);
        double u2y = _sol.coeff(n2+1);
        double u3x = _sol.coeff(n3);
        double u3y = 0.0;
        // if(n3 == _triangles.rows()-1){
        //   u3y = 100 ;
        // }
        // else{
        //   u3y = _sol.coeff(n3+1);
        // }
        

        Eigen::VectorXd U(6); 

        // Remplissage de U
        U << u1x, u1y, u2x, u2y, u3x, u3y;

        Eigen::VectorXd epsilon = B * U;

        double epsilon_x = epsilon[0];
        double epsilon_y = epsilon[1];
        double gamma_xy = epsilon[2];

        Eigen::VectorXd sigma_Matrice(3); 
        Eigen::VectorXd sigma_Fibre(3);   
        if (_triangles(K, 3) == 100) {   

            Eigen::MatrixXd D = _assemblage->_geometry->coef_elastique(E_matrice, nu_matrice);

            sigma_Matrice = D * epsilon;
            _contrainte_Matrice.row(K) = sigma_Matrice.transpose(); 
            _epsilon_xx_M += epsilon_x;
            _epsilon_yy_M += epsilon_y;
            _epsilon_xy_M += gamma_xy;
            _N_M++;
        }
        else {

            Eigen::MatrixXd D = _assemblage->_geometry->coef_elastique(E_fibre, nu_fibre);

            sigma_Fibre = D * epsilon;
            _contrainte_Fibre.row(K) = sigma_Fibre.transpose(); 
            _epsilon_xx_F += epsilon_x;
            _epsilon_yy_F += epsilon_y;
            _epsilon_xy_F += gamma_xy;
            _N_F++;
        }
    }
}



// // test uniquement sur le bord

// void Resolution::contrainte(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice) {
//     // const double delta_x = 1e-6; // différence finie
//     _contrainte_Matrice.resize(_triangles.rows(), 3); // Redimensionner la matrice des contraintes
//     _contrainte_Fibre.resize(_triangles.rows(), 3);

//     // boucle sur les éléments
//     for (int K = 0; K < _assemblage->getRefVerticesDirichletWithRef().size(); K++) {
//         int n1 = _assemblage->getRefVerticesDirichletWithRef()[K].first;
//         int n2 = _assemblage->getRefVerticesDirichletWithRef()[K].second;
//         int n3 = _triangles(K, 2);

//         Eigen::MatrixXd B = _assemblage->_geometry->BoF_alpha(K);
//         double u1x = _sol.coeff(n1);
//         double u1y = _sol.coeff(n1+1);
//         double u2x = _sol.coeff(n2);
//         double u2y = _sol.coeff(n2+1);
//         double u3x = _sol.coeff(n3);
//         double u3y = 0.0;

//         Eigen::VectorXd U(6); 

//         // Remplissage de U
//         U << u1x, u1y, u2x, u2y, u3x, u3y;

//         Eigen::VectorXd epsilon = B * U;

//         double epsilon_x = epsilon[0];
//         double epsilon_y = epsilon[1];
//         double gamma_xy = epsilon[2];

//         Eigen::VectorXd sigma_Matrice(3); 
//         Eigen::VectorXd sigma_Fibre(3);   
//         if (_triangles(K, 3) == 100) {   
//             Eigen::MatrixXd D = _assemblage->_geometry->coef_elastique(E_matrice, nu_matrice);

//             sigma_Matrice = D * epsilon;
//             _contrainte_Matrice.row(K) = sigma_Matrice.transpose(); 
//             _epsilon_xx_M += epsilon_x;
//             _epsilon_yy_M += epsilon_y;
//             _epsilon_xy_M += gamma_xy;
//             _N_M++;
//         }
//         else {
//             Eigen::MatrixXd D = _assemblage->_geometry->coef_elastique(E_fibre, nu_fibre);

//             sigma_Fibre = D * epsilon;
//             _contrainte_Fibre.row(K) = sigma_Fibre.transpose(); 
//             _epsilon_xx_F += epsilon_x;
//             _epsilon_yy_F += epsilon_y;
//             _epsilon_xy_F += gamma_xy;
//             _N_F++;
//         }
//     }
// }




void Resolution::afficherContrainte() {
    double moyenne_sigma_x_Matrice = 0.0;
    double moyenne_sigma_y_Matrice = 0.0;
    double moyenne_tau_xy_Matrice = 0.0;

    double moyenne_sigma_x_Fibre = 0.0;
    double moyenne_sigma_y_Fibre = 0.0;
    double moyenne_tau_xy_Fibre = 0.0;

    int nombre_elements = _assemblage->getRefVerticesDirichletWithRef().size();

    std::ofstream fichierMatrice("contrainte_matrice.txt");
    std::ofstream fichierFibre("contrainte_fibre.txt");

    if (fichierMatrice.is_open() && fichierFibre.is_open()) {
        for (int i = 0; i < nombre_elements; ++i) {
            // Accéder aux éléments de _contrainte_Matrice et _contrainte_Fibre
            double sigma_x_M = _contrainte_Matrice.coeff(i, 0);
            double sigma_y_M = _contrainte_Matrice.coeff(i, 1);
            double tau_xy_M = _contrainte_Matrice.coeff(i, 2);

            double sigma_x_F = _contrainte_Fibre.coeff(i, 0);
            double sigma_y_F = _contrainte_Fibre.coeff(i, 1);
            double tau_xy_F = _contrainte_Fibre.coeff(i, 2);

            // Ajouter les valeurs à la somme totale pour chaque matériau
            moyenne_sigma_x_Matrice += sigma_x_M;
            moyenne_sigma_y_Matrice += sigma_y_M;
            moyenne_tau_xy_Matrice += tau_xy_M;

            moyenne_sigma_x_Fibre += sigma_x_F;
            moyenne_sigma_y_Fibre += sigma_y_F;
            moyenne_tau_xy_Fibre += tau_xy_F;

            // Écrire les valeurs de contrainte pour chaque élément dans les fichiers
            fichierMatrice << "Contraintes pour l'élément " << i << ":" << std::endl;
            fichierMatrice << "Sigma_x = " << sigma_x_M << ", Sigma_y = " << sigma_y_M << ", Tau_xy = " << tau_xy_M << std::endl;

            fichierFibre << "Contraintes pour l'élément " << i << ":" << std::endl;
            fichierFibre << "Sigma_x = " << sigma_x_F << ", Sigma_y = " << sigma_y_F << ", Tau_xy = " << tau_xy_F << std::endl;
        }

        // Calculer les moyennes
        moyenne_sigma_x_Matrice /= _N_M;
        moyenne_sigma_y_Matrice /= _N_M;
        moyenne_tau_xy_Matrice /= _N_M;

        moyenne_sigma_x_Fibre /= _N_F;
        moyenne_sigma_y_Fibre /= _N_F;
        moyenne_tau_xy_Fibre /= _N_F;


        _contrainte_xx_M = moyenne_sigma_x_Matrice;
        _contrainte_yy_M = moyenne_sigma_y_Matrice;
        _contrainte_xy_M = moyenne_tau_xy_Matrice;

        _contrainte_xx_F = moyenne_sigma_x_Fibre;
        _contrainte_yy_F = moyenne_sigma_y_Fibre;
        _contrainte_xy_F = moyenne_tau_xy_Fibre;


        std::cout<< "sigma_x = "<< _contrainte_xx_F + _contrainte_xx_M<<std::endl;


        _epsilon_xx_M /= _N_M;
        _epsilon_yy_M /= _N_M;
        _epsilon_xy_M /= _N_M;

        _epsilon_xx_F /= _N_F;
        _epsilon_yy_F /= _N_F;
        _epsilon_xy_F /= _N_F;

        

        if(_datafile->getDirichletReferences()[0]==1 && _datafile->getDirichletReferences()[1]==3){

          std::cout<<"Traction dans la direction x"<<std::endl;
            // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la matrice
          _E_xx_M = _contrainte_xx_M / _epsilon_xx_M;
          _nu_xy_M = -(_epsilon_yy_M * _E_xx_M) / _contrainte_xx_M;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la fibre
          _E_xx_F = _contrainte_xx_F/_epsilon_xx_F ;
          _nu_xy_F = -(_epsilon_yy_F * _E_xx_F) / _contrainte_xx_F;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour le matériau composite
          _E_xx = abs(_E_xx_F)+abs(_E_xx_M);
          _nu_xy=abs(_nu_xy_F)+abs(_nu_xy_M);

           // Affichage des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la matrice
          std::cout << "Module d'élasticité (E_xx) pour la matrice : " << _E_xx_M << std::endl;
          std::cout << "Coefficient de Poisson (nu) pour la matrice : " << _nu_xy_M << std::endl;

          // Affichage des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la fibre
          std::cout << "Module d'élasticité (E_xx) pour la fibre : " << _E_xx_F << std::endl;
          std::cout << "Coefficient de Poisson (nu) pour la fibre : " << _nu_xy_F << std::endl;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour le matériau composite
          std::cout << "Module d'élasticité (E_xx) pour le matériau : " << _E_xx << std::endl;
          std::cout << "Coefficient de Poisson (nu) pour le matériau : " << _nu_xy << std::endl;

          // Calcul de la surface occupée par les fibres
          _Sf=static_cast<float>(_N_F)/(_N_F+_N_M) ;
          std::cout << "la surface occupée par les fibres : " << _Sf << std::endl;


        }
        else if(_datafile->getDirichletReferences()[0]==2 && _datafile->getDirichletReferences()[1]==4){
          std::cout<<"Traction dans la direction y"<<std::endl;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la matrice
          _E_yy_M = _contrainte_yy_M / _epsilon_yy_M;
          _nu_yx_M = -(_epsilon_xx_M * _E_yy_M) / _contrainte_yy_M;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la fibre
          _E_yy_F = _contrainte_yy_F / _epsilon_yy_F;
          _nu_yx_F = -(_epsilon_xx_F * _E_yy_F) / _contrainte_yy_F;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour le matériau composite
          _E_yy = abs(_E_yy_F)+abs(_E_yy_M)/2.0;
          _nu_yx=abs(_nu_yx_F)+abs(_nu_yx_M)/2.0;

           // Affichage des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la matrice
          std::cout << "Module d'élasticité (E_yy) pour la matrice : " << _E_yy_M << std::endl;
          std::cout << "Coefficient de Poisson (nu) pour la matrice : " << _nu_yx_M << std::endl;

          // Affichage des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour la fibre
          std::cout << "Module d'élasticité (E_yy) pour la fibre : " << _E_yy_F << std::endl;
          std::cout << "Coefficient de Poisson (nu) pour la fibre : " << _nu_yx_F << std::endl;

          // Calcul des modules d'élasticité (E_xx) et du coefficient de Poisson (nu) pour le matériau composite
          std::cout << "Module d'élasticité (E_yy) pour le matériau : " << _E_yy << std::endl;
          std::cout << "Coefficient de Poisson (nu_yx) pour le matériau : " << _nu_yx << std::endl;

          // Calcul de la surface occupée par les fibres
          _Sf=static_cast<float>(_N_F)/(_N_F+_N_M) ;
          std::cout << "la surface occupée par les fibres : " << _Sf << std::endl;

        }
        else{
          std::cout<<"Le cisaillement n'est pas possible "<<std::endl;
        }

        


       




        // Calcul analytique des modules d'élasticité (E_xx) pour le matériau composite 
          // double E_fibre=450e3;
          // double E_matrice=3.45e3;
          // double nu_fibre=0.32;
          // double nu_matrice=0.3;
          // _E_xx_th=1/((_Sf/E_fibre)+(1-_Sf)/E_matrice);
          // std::cout << "Module d'élasticité (E_xx) pour le matériau d'après la théorie des mélanges : " << _E_xx_th << std::endl;
        

       



        // Écrire les moyennes dans les fichiers
        fichierMatrice << "Moyenne des contraintes pour la matrice :" << std::endl;
        fichierMatrice << "Sigma_x = " << moyenne_sigma_x_Matrice << ", Sigma_y = " << moyenne_sigma_y_Matrice << ", Tau_xy = " << moyenne_tau_xy_Matrice << std::endl;

        fichierFibre << "Moyenne des contraintes pour la fibre :" << std::endl;
        fichierFibre << "Sigma_x = " << moyenne_sigma_x_Fibre << ", Sigma_y = " << moyenne_sigma_y_Fibre << ", Tau_xy = " << moyenne_tau_xy_Fibre << std::endl;

        // Fermer les fichiers
        fichierMatrice.close();
        fichierFibre.close();
        std::cout << "Les moyennes des contraintes ont été enregistrées dans les fichiers contrainte_matrice.txt et contrainte_fibre.txt." << std::endl;
        std::cout<<_N_F<<std::endl;
    } else {
        std::cerr << "Erreur : Impossible d'ouvrir les fichiers pour l'écriture." << std::endl;
    }
}


// void Resolution::afficherContrainte() {
//     double sigma_x_max = _contrainte.coeff(0, 0); // Initialisation avec la première valeur
//     double sigma_y_max = _contrainte.coeff(0, 1);
//     double tau_xy_max = _contrainte.coeff(0, 2);
    
//     // Calculer le nombre total d'éléments
//     int nombre_elements = _triangles.rows();

//     // Ouvrir le fichier en mode écriture
//     std::ofstream fichier("contrainte.txt");

//     // Vérifier si le fichier est correctement ouvert
//     if (fichier.is_open()) {
//         // Parcourir chaque élément du maillage
//         for (int i = 0; i < nombre_elements; ++i) {
//             // Écrire les valeurs de contrainte pour chaque élément dans le fichier
//             fichier << "Contraintes pour l'élément " << i << ":" << std::endl;
//             fichier << "Sigma_x = " << _contrainte.coeff(i, 0) << ", Sigma_y = " << _contrainte.coeff(i, 1) << ", Tau_xy = " << _contrainte.coeff(i, 2) << std::endl;

//             // Mettre à jour les valeurs maximales si nécessaire
//             sigma_x_max = std::max(sigma_x_max, _contrainte.coeff(i, 0));
//             sigma_y_max = std::max(sigma_y_max, _contrainte.coeff(i, 1));
//             tau_xy_max = std::max(tau_xy_max, _contrainte.coeff(i, 2));
//         }

//         // Écrire les valeurs maximales des contraintes dans le fichier
//         fichier << "Valeurs maximales des contraintes :" << std::endl;
//         fichier << "Sigma_x_max = " << sigma_x_max << ", Sigma_y_max = " << sigma_y_max << ", Tau_xy_max = " << tau_xy_max << std::endl;

//         // Fermer le fichier
//         fichier.close();
//         std::cout << "Les données des contraintes ont été enregistrées dans le fichier contrainte.txt." << std::endl;
//     } else {
//         // Afficher un message d'erreur si le fichier n'a pas pu être ouvert
//         std::cerr << "Erreur : Impossible d'ouvrir le fichier contrainte.txt pour l'écriture." << std::endl;
//     }
// }


// void Resolution::afficherContrainte() {
//     double sigma_x_max = _contrainte.coeff(0, 0); // Initialisation avec la première valeur
//     double sigma_y_max = _contrainte.coeff(0, 1);
//     double tau_xy_max = _contrainte.coeff(0, 2);
    
//     // Calculer le nombre total d'éléments
//     int nombre_elements = _triangles.rows();

//     // Parcourir chaque élément du maillage
//     for (int i = 0; i < nombre_elements; ++i) {
//         // Afficher les valeurs de contrainte pour chaque élément
//         std::cout << "Contraintes pour l'élément " << i << ":" << std::endl;
//         std::cout << "Sigma_x = " << _contrainte.coeff(i, 0) << ", Sigma_y = " << _contrainte.coeff(i, 1) << ", Tau_xy = " << _contrainte.coeff(i, 2) << std::endl;

//         // Mettre à jour les valeurs maximales si nécessaire
//         sigma_x_max = std::max(sigma_x_max, _contrainte.coeff(i, 0));
//         sigma_y_max = std::max(sigma_y_max, _contrainte.coeff(i, 1));
//         tau_xy_max = std::max(tau_xy_max, _contrainte.coeff(i, 2));
//     }

//     // Afficher les valeurs maximales des contraintes
//     std::cout << "Valeurs maximales des contraintes :" << std::endl;
//     std::cout << "Sigma_x_max = " << sigma_x_max << ", Sigma_y_max = " << sigma_y_max << ", Tau_xy_max = " << tau_xy_max << std::endl;
// }


#define _RESOLUTION_H
#endif
