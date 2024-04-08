#ifndef _INIT_COND_H

#include<cmath>
#include<iostream>
#include<vector>
#include "Dense"
#include "Sparse"

// le terme source et les cndition aux bords de Dirichlet et de Neumann
double sourceTerm(Eigen::Vector2d X);
double neumannBC(Eigen::Vector2d X);
double DirichletBC(Eigen::Vector2d X, int ref, double traction);

#define _INIT_COND_H
#endif