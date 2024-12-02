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
        DataFile* _datafile;
        Eigen::Matrix<double, Eigen::Dynamic, 2> _vertices;
        Eigen::Matrix<int, Eigen::Dynamic, 4> _triangles;
        Eigen::SparseVector<double> _sol;
        
        //calcul de contrainte 

        Eigen::MatrixXd _contrainte_Matrice;
        Eigen::MatrixXd _contrainte_Fibre;

        //recuperation des trois composante de deformation et de contrainte 
        double _epsilon_xx_F; 
        double _epsilon_yy_F;
        double _epsilon_xy_F;

        double _epsilon_xx_M; 
        double _epsilon_yy_M;
        double _epsilon_xy_M;



        double _contrainte_xx_M;
        double _contrainte_yy_M;
        double _contrainte_xy_M;

        double _contrainte_xx_F;
        double _contrainte_yy_F;
        double _contrainte_xy_F;

        double _E_xx_M;
        double _nu_xy_M;

        double _E_xx_F;
        double _nu_xy_F;

        double _E_xx;
        double _nu_xy;

        double _E_yy_M;
        double _nu_yx_M;

        double _E_yy_F;
        double _nu_yx_F;

        double _E_yy;
        double _nu_yx;


        int _N_F;
        int _N_M;

        double _Sf;
        double _E_xx_th;



        
    
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








    


