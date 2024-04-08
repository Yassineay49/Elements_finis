#ifndef _RESOLUTION_H

#include <vector>
#include <string>
#include "Dense"
#include "Sparse"
#include "DataFile.h"
#include "Assemblage.h"
#include "Solveur.h"

class Resolution
{
    protected:
        Solver * _solver;
        Assemblage* _assemblage;
        Eigen::Matrix<double, Eigen::Dynamic, 2> _vertices;
        Eigen::Matrix<int, Eigen::Dynamic, 4> _triangles;
        Eigen::SparseVector<double> _sol;
        
        //calcul de contrainte 

        Eigen::SparseMatrix<double> _contrainte;
    
    public:
        Resolution(DataFile* dataFile, Assemblage* assemblage, Solver* solver, Mesh* mesh);
        virtual ~Resolution();
        void saveSolution();
        void systeme(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice);
        void saveSolution(Eigen::SparseVector<double> sol);

        //calcul de contrainte
        void contrainte(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice);
        void afficherContrainte();
};

#define _RESOLUTION_H
#endif








    


