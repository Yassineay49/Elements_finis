#ifndef _dataFile_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _fileName(file_name),  _ifMeshName(false), _ifResultsFolder(false),
_ifDirichlet(false), _ifNeumann(false),
_ifICAndSTFile(false), _itraction(false)
{

}

void DataFile::readDataFile()
{
  ifstream dataFile(_fileName.data());
  if (!dataFile.is_open())
  {
    cout << "Unable to open file " << _fileName << endl;
    abort();
  }
  else
  {
    cout << "Reading data file " << _fileName << endl;
  }

  string fileLine;

  while (!dataFile.eof())
  {
    getline(dataFile, fileLine);
    if (fileLine.find("mesh") != std::string::npos)
    {
      dataFile >> _meshName; _ifMeshName = true;
    }

    if (fileLine.find("traction") != std::string::npos)
    {
      dataFile >> _traction; _itraction = true;
    }

    // if (fileLine.find("initialCondition_and_sourceTerm_file") != std::string::npos)
    // {
    //   dataFile >> _ICAndSTFile; _ifICAndSTFile = true;
    //   //system(("cp -r ./" + _ICAndSTFile + " ./InitialConditionSourceTermFile.cpp").c_str());
    // }

    // if (fileLine.find("is_exact_sol") != std::string::npos)
    // {
    //   getline(dataFile, fileLine);
    //   if (fileLine == "true")
    //   {
    //     _isExactSol = true;
    //   }
    //   _ifIsExactSol = true;
    // }

    if (fileLine.find("dirichlet") != std::string::npos)
    {
      getline(dataFile, fileLine);
      std::string::size_type sz;   // alias of size_t
      for (int i = 0 ; i < fileLine.size() ; i++)
      {
        if (i%2==0)
        {
          int temp = std::stoi(fileLine.substr(i,1),&sz);
          _dirichlet.push_back(temp);
          _BCReferences.push_back(temp);
        }
      }
      _ifDirichlet = true;
    }

    if (fileLine.find("neumann") != std::string::npos)
    {
      getline(dataFile, fileLine);
      std::string::size_type sz;   // alias of size_t
      for (int i = 0 ; i < fileLine.size() ; i++)
      {
        if (i%2==0)
        {
          int temp = std::stoi(fileLine.substr(i,1),&sz);
          _neumann.push_back(temp);
          _BCReferences.push_back(temp);
        }
      }
      _ifNeumann = true;
    }
  }

  if (!_ifResultsFolder)
  {
    cout << "Beware - The default results folder name (results) is used." << endl;
    _resultsFolder = "results";
  }
  if (!_ifMeshName)
  {
    cout << "Do not forget to give the mesh name in the data file." << endl;
    abort();
  }
  // if (!_ifICAndSTFile)
  // {
  //   cout << "Do not forget to give the name of the file containing the initial condition and the source term." << endl;
  //   abort();
  // }
  
  if (!_ifDirichlet)
  {
    cout << "There is no Dirichlet boundary conditions." << endl;
  }
  else
  {
    cout << "Dirichlet BC(s) is (are) imposed on reference(s):";
    for (int i = 0 ; i < _dirichlet.size() ; i++)
    cout << " " << _dirichlet[i];
    cout << "." << endl;
  }
  if (!_ifNeumann)
  {
    cout << "There is no Neumann boundary conditions." << endl;
  }
  else
  {
    cout << "Neumann BC(s) is (are) imposed on reference(s):";
    for (int i = 0 ; i < _neumann.size() ; i++)
    cout << " " << _neumann[i];cout << "." << endl;
  }
}
#define _dataFile_CPP
#endif
