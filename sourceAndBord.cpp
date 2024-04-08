#ifndef _INIT_COND_H

#include "sourceAndBord.h"

// les forces volumique (considérer ici comme le terme de source)
double sourceTerm(Eigen::Vector2d X)
{
    // A modifier en fonction des parametre à considerer
    return 0;
}

// condition aux bords de Neumann
double neumannBC(Eigen::Vector2d X)
{
    double nx = 0;
    double ny = 0;
    if (abs(X(1))   < 1.0e-10) {ny = -1;}
    if (abs(X(0)-1) < 1.0e-10) {nx =  1;}
    if (abs(X(0))   < 1.0e-10) {nx = -1;}
    if (abs(X(1)-1) < 1.0e-10) {ny =  1;}
    nx = nx/sqrt(nx*nx+ny*ny);
    ny = ny/sqrt(nx*nx+ny*ny);

    // A definir en fonction des parametre du pb

    return 0;


}

// condition dirichlet
double DirichletBC(Eigen::Vector2d X, int ref, double traction)
{
    // dependra de la force de traction que nous allons appliquee
    if(ref>=5){
        return 100;
    }
    else
        return traction;
} 


#define _INIT_COND_H
#endif