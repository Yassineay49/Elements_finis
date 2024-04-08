#ifndef _ASSEMBLAGE_H


#include "Assemblage.h"
#include "sourceAndBord.h"

Assemblage::Assemblage(Mesh* mesh, Geometry* geometry, DataFile* dataFile) :
_geometry(geometry), _triangles(mesh->getTriangles()), _vertices(mesh->getVertices()), _edgesNeumann(mesh->getNeumannEdges()),
_refVerticesDirichletWithRef(mesh->getRefVerticesDirichlet()), _Kglobal(_vertices.rows(), _vertices.rows()),
_sourceAndNeumann(_vertices.rows()), _dataFile(dataFile)
{
    // Weights and Points
  _geometry->quadraturePointsAndWeightsMidpointFormula(_weightsMidPoints, _pointsMidPoints);
  _pointsSimpson.resize(3);
  _geometry->quadraturePointsAndWeightsSimpsonFormula(_weightsSimpson, _pointsSimpson[0], 1);
  _geometry->quadraturePointsAndWeightsSimpsonFormula(_weightsSimpson, _pointsSimpson[1], 2);
  _geometry->quadraturePointsAndWeightsSimpsonFormula(_weightsSimpson, _pointsSimpson[2], 3);
}

void Assemblage::matrice_globalK(double E_fibre, double E_matrice, double nu_fibre, double nu_matrice)
{
    std::vector<Eigen::Triplet<double>> tripletsK;
    _Kglobal.setZero();

    Eigen::MatrixXd K_elem(3,6);

    // Boucle sur les elements K dans T_h
    for (int K = 0; K < _triangles.rows() ; K++)
    {
      //trouver  un moyen de traiter la matrice elemataire
      //sur la fibre et sur la matrice

      if(_triangles(K,3)=100){
        K_elem=_geometry->Elementaire(K, nu_matrice, E_matrice);
      }
      else 
        K_elem=_geometry->Elementaire(K,nu_fibre, E_fibre);

      // std::cout << "Progression de l'assemblage " << (K+1)/(1.0*_triangles.rows())*100 << "%" << std::endl;
      // Boucle sur les points de quadrature w_q, hatX_q
      for (int q = 0 ; q < _weightsMidPoints.rows() ; q++)
      {

        // Boucle sur les noeuds de l'élément de référence géométrique
        for (int hati = 0 ; hati < 3 ; hati++)
        {
          // Boucle sur les noeuds de l'élément de référence géométrique
          for (int hatj = 0 ; hatj < 3 ; hatj++)
          {
            // i = numéro global du noeud hati de K
            int i(_triangles(K,hati));
            // j = numéro global du noeud hatj de K
            int j(_triangles(K,hatj));

            // K(i,j) = K(i,j) + constribution de la matrice elementaire
            double Element(0.5 //aire du triangle de reference
                  *_weightsMidPoints(q)
                  *K_elem(hati,hatj)
            );
            tripletsK.push_back({i,j,Element});
          }
        }
      }
    }

    _Kglobal.setFromTriplets(tripletsK.begin(), tripletsK.end());

}


//second membre du systeme lineaire 

void Assemblage::assemblesourceNeumann()
{

  _sourceAndNeumann.setZero();
  // Boucle sur les elements K dans T_h
   for (int K = 0; K < _triangles.rows() ; K++)
  {
    // Boucle sur les points de quadrature w_q, hatX_q
    for (int q = 0 ; q < _weightsMidPoints.rows() ; q++)
    {
      // Boucle sur les noeuds de l'élément de référence géométrique
      for (int hati = 0 ; hati < 3 ; hati++)
      {
          // i = numéro global du noeud hati de K
          int i(_triangles(K,hati));

          //contribution du terme source 
          double termesource = 0;
          termesource = sourceTerm(_geometry->Falpha(K,_pointsMidPoints.row(q)));

          double contrib(0.5 // aire du triangle de reference
                *_weightsMidPoints(q)
                *_geometry->phihat(hati, _pointsMidPoints.row(q))
                *termesource
                *_geometry->absoludetJ(K));

          _sourceAndNeumann.coeffRef(i)=_sourceAndNeumann.coeffRef(i)+contrib;
      }
    }
  }

  //condition aux bord Neumann non homogene 

  // Boucle sur les arrets E de bord du domaine (ici on parle d'abord du grand domaine tout entier)

  Eigen::SparseVector<double> neumann(_sourceAndNeumann.size());
  neumann.setZero();

  for (int E = 0; E < _edgesNeumann.rows() ; E++)
  {
    //int ref = _edgesNeumann(E,2);
    int K = _edgesNeumann(E,3);
    int numedge = _edgesNeumann(E,4);

    // Boucle sur les points de quadrature w_q, hatX_q
    for (int q = 0 ; q < _weightsSimpson.rows() ; q++)
    {
      // Boucle sur les noeuds de l'élément de référence géométrique
      int hati;
      for (int temphati = 0 ; temphati < 2 ; temphati++)
      {
          // i = numéro global du noeud hati de K
          if (numedge==1) {hati = temphati;}
          else if (numedge==2) {if (temphati == 0) hati=1; if (temphati == 1) hati=2;}
          else if (numedge==3) {if (temphati == 0) hati=2; if (temphati == 1) hati=0;}

          int i(_triangles(K,hati));

          double neumannint = 0; 

          neumannint = neumannBC(_geometry->FE(E, _pointsSimpson[0].row(q)));

          double contrib(1.0 // longeur de l'element de reference
                *_weightsSimpson(q)
                *neumannint
                *_geometry->phihat(hati, _pointsSimpson[numedge-1].row(q))
                *_geometry->measE(E));

                neumann.coeffRef(i)=neumann.coeffRef(i)+contrib;
      }
    }
  }

  _sourceAndNeumann += neumann;

  
}

void Assemblage::applyBCToSystemMatrix(Eigen::SparseMatrix<double,Eigen::RowMajor> & systemMatrix)
{
  Eigen::SparseVector<double> zeroRow(systemMatrix.cols());
  for (int i = 0; i < _refVerticesDirichletWithRef.size() ; i++)
  {
    int refV = _refVerticesDirichletWithRef[i].first;
    systemMatrix.row(refV) = zeroRow;
    systemMatrix.coeffRef(refV,refV) = 1.;
  }
}

void Assemblage::applyBCToRHS(Eigen::SparseVector<double> & RHS)
{
  double traction = _dataFile->Get_traction();
  std::cout<<"Force de traction F= "<<traction<<std::endl;

  for (int i = 0; i < _refVerticesDirichletWithRef.size() ; i++)
  {
    // taiter separareent le cas des fibres et mettre des deplacements nuls ?
    
    int refV = _refVerticesDirichletWithRef[i].first;
    RHS.coeffRef(refV) = DirichletBC( _vertices.row(refV),_refVerticesDirichletWithRef[i].second, traction);
  }
}





#define _ASSEMBLAGE_H
#endif