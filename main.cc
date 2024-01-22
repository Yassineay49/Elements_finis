#include <iostream>
#include <fstream>
#include <chrono>

#include "mesh.h"
#include "Geometry.h"
#include "DataFile.h"


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


  cout<<geometry->JFalpha(2)<<endl;
  cout<<geometry->BoF_alpha(2)<<endl;
  cout<<geometry->Elementaire(2,0.3,0.5)<<endl;


  return 0;

}
