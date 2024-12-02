#include <iostream>
#include <fstream>
#include <chrono>

#include "mesh.h"
#include "Geometry.h"
#include "DataFile.h"
#include "Assemblage.h"
#include "Resolution.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv)
{
    if (argc < 2)
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

    Geometry* geometry = new Geometry(mesh);
    std::cout << "-------------------------------------------------" << std::endl;
    // ---------------------------------------------------------------------------//
    Assemblage* assemblage =  new Assemblage(mesh, geometry, dataFile);
    std::cout << "-------------------------------------------------" << std::endl;
    // ---------------------------------------------------------------------------//

    // --------------------------  Choix du solveur ----------------------------//
    cout << "Choose solver:" << endl;
    cout << "1. EigenSolver" << endl;
    cout << "2. BiCGSTABSolver" << endl;

    int solverChoice;
    cin >> solverChoice;

    Solver* solver;
    switch (solverChoice)
    {
        case 1:
            solver = new EigenSolver();
            break;
        case 2:
            solver = new BiCGSTABSolver();
            break;
        default:
            cout << "Invalid choice. Using default solver (EigenSolver)." << endl;
            solver = new EigenSolver();
    }
    // -------------------------------------------------------------------------//

    // --------------------------- Resolution -----------------------------
    Resolution* resolution = new Resolution(dataFile, assemblage, solver, mesh);
    // ---------------------------------------------------------------------------

    std::cout << "-------------------------------------------------" << std::endl;
    double E_fibre(311);
    double E_matrice(275);
    double nu_fibre(0.19);
    double nu_matrice(0.33);

    auto startResolution = high_resolution_clock::now(); // Enregistrer le temps de début de la résolution
    resolution->systeme(E_fibre, E_matrice, nu_fibre, nu_matrice);
    resolution->saveSolution();
    resolution->contrainte(E_fibre, E_matrice, nu_fibre, nu_matrice);
    resolution->afficherContrainte();
    auto endResolution = high_resolution_clock::now(); // Enregistrer le temps de fin de la résolution
    auto durationResolution = duration_cast<milliseconds>(endResolution - startResolution);
    cout << "Resolution time: " << durationResolution.count() << " milliseconds" << endl;

    delete solver; // Nettoyer la mémoire
    delete resolution;
    delete assemblage;
    delete geometry;
    delete mesh;
    delete dataFile;

    return 0;
}
