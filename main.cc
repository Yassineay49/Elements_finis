#include <iostream>
#include <fstream>
#include <chrono>

#include "mesh.h"
#include "Geometry.h"
#include "DataFile.h"
#include "Assemblage.h"
#include "Resolution.h"



using namespace std;

int main(int argc, char** argv)
{ if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string dataFile_name = argv[1];

  std::cout << "-------------------------------------------------" << std::endl;
  // ----------------------- Fichier de données --------------------------------
  DataFile* dataFile = new DataFile(dataFile_name);
  dataFile->readDataFile();
  std::cout << "-------------------------------------------------" << std::endl;
  // ---------------------------------------------------------------------------

  // ----- Lecture du maillage et construction des entités géométriques --------
  Mesh* mesh = new Mesh(dataFile->getMeshName(),dataFile->getDirichletReferences(),
                                                 dataFile->getNeumannReferences());

  //std::cout<<"Programme"<<std::endl;

  Geometry* geometry = new Geometry(mesh);
  std::cout << "-------------------------------------------------" << std::endl;
  // ---------------------------------------------------------------------------//
  Assemblage* assemblage =  new Assemblage(mesh, geometry, dataFile);
  std::cout << "-------------------------------------------------" << std::endl;
  // ---------------------------------------------------------------------------//

  // --------------------------  Solveur ----------------------------//
  Solver* solver =new EigenSolver();
  // -------------------------------------------------------------------------//

  // --------------------------- Resolution -----------------------------
  Resolution* resolution = new Resolution(dataFile, assemblage, solver, mesh);
  // ---------------------------------------------------------------------------

  std::cout << "-------------------------------------------------" << std::endl;
  double E_fibre(311e3);
  double E_matrice(275e3);
  double nu_fibre(0.19);
  double nu_matrice(0.33);

  resolution->systeme(E_fibre, E_matrice, nu_fibre, nu_matrice);
  resolution->saveSolution();
  resolution->contrainte(E_fibre, E_matrice, nu_fibre, nu_matrice);
  resolution->afficherContrainte();






  //test de l'assemblage 
  // printf("\n");
  // 
  
  // assemblage->matrice_globalK(210000, 790000, 0.3, 0.4);
  // cout<<geometry->JFalpha(2)<<endl;
  // cout<<geometry->BoF_alpha(2)<<endl;
  // cout<<geometry->Elementaire(2,210000,0.3)<<endl;
  // cout<<assemblage->getStiffnessMatrix()<<endl;


  return 0;

}
