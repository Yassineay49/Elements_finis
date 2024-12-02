#ifndef _ASSEMBLAGE_H

#include<vector>
#include<fstream>
#include "Dense"
#include "Sparse"
#include "mesh.h"
#include "Geometry.h"
#include "DataFile.h"

class Assemblage
{
    // coordonnées des sommets du maillge
    private:
        // Coordonnées des sommets du maillage
        Eigen::Matrix<double, Eigen::Dynamic, 2> _vertices;
        // Liste de tous les triangles et leur référence
        Eigen::Matrix<int, Eigen::Dynamic, 4> _triangles;
        // Liste des références de points concernés par Dirichlet
        std::vector<std::pair<int,int>> _refVerticesDirichletWithRef;
        // Liste de toutes les arêtes avec condition de Neumann, leur référence et le triangle associé
        Eigen::Matrix<int, Eigen::Dynamic, 5> _edgesNeumann;


        
        // Numerotation local
        Eigen::VectorXd _numlocal;
        // DataFile 
        DataFile* _dataFile;
        // Numerotation global 
        Eigen::VectorXd _numglobal;
        // Liste des points et des poids de quadrature pour la formule du Milieu (2D)
        Eigen::VectorXd _weightsMidPoints;
        Eigen::Matrix<double, Eigen::Dynamic, 2> _pointsMidPoints;
        // Liste des points et des poids de quadrature pour la formule du Simpson (1D)
        Eigen::VectorXd _weightsSimpson;
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> _pointsSimpson;



        // Matrice K global
        Eigen::SparseMatrix<double,Eigen::RowMajor> _Kglobal;
        // Second membre du systeme 
        Eigen::SparseVector<double> _sourceAndNeumann;

    public:
        // Geometrie
        Geometry * _geometry;
        Assemblage(Mesh*mesh, Geometry* geometry, DataFile*dataFile);
        void matrice_globalK(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice);
        void assemblesourceNeumann();
        void applyBCToSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> & systemMatrix);
        void applyBCToRHS(Eigen::SparseVector<double> & RHS);
        const Eigen::SparseMatrix<double,Eigen::RowMajor> & getStiffnessMatrix() const {return _Kglobal;};
        const Eigen::SparseVector<double> & getSourceAndNeumann() const {return _sourceAndNeumann;};
        const std::vector<std::pair<int,int>> & getRefVerticesDirichletWithRef() const {return _refVerticesDirichletWithRef;};



};

#define _ASSEMBLAGE_H
#endif