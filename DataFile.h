#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile {
private:
  std::string _fileName;
  std::string _meshName;
  std::string _resultsFolder;
  std::string _ICAndSTFile;
  std::vector<int> _dirichlet;
  std::vector<int> _neumann;
  std::vector<int> _BCReferences;

  bool  _ifMeshName, _ifResultsFolder, _ifDirichlet, _ifNeumann , _ifICAndSTFile;

public: // Méthodes et opérateurs de la classe
  DataFile(std::string file_name);
  void readDataFile();
  const std::string & getMeshName() const {return _meshName;};
  const std::string & getResultsFolder() const {return _resultsFolder;};
  const std::vector<int> & getDirichletReferences() const {return _dirichlet;};
  const std::vector<int> & getNeumannReferences() const {return _neumann;};
};

#define _DATA_FILE_H
#endif
