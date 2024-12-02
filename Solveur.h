#ifndef _SOLVER_H

#include "Sparse"
#include <SparseQR>
#include <IterativeLinearSolvers>



class Solver
{
public:
  Solver();
  virtual ~Solver();
  virtual void setSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> systemMatrix) = 0;
  virtual Eigen::SparseVector<double> solve(Eigen::SparseVector<double> RHS) = 0;
};

class EigenSolver : public Solver
{
private:
  Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> >  _solver;
public:
  void setSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> systemMatrix);
  Eigen::SparseVector<double> solve(Eigen::SparseVector<double> RHS);
};


class BiCGSTABSolver: public Solver // Define the BiCGSTAB solver class
{
private:
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> _solver; // Use BiCGSTAB solver
public:
    void setSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> systemMatrix);
    Eigen::SparseVector<double> solve(Eigen::SparseVector<double> RHS);
};


#endif