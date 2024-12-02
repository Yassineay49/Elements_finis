cmake_minimum_required(VERSION 3.0)
project(PG201_TP1 VERSION 2023)

set(CXX_STANDARD_REQUIRED On)
set(CXX_STANDARD C++11)

# Chemin vers la bibliothèque Eigen
set(LIBRARY_PATH_EIGEN "/home/ays/2A/PG201/pg201-students-m2-tp2a-yayoun/libraries/eigen/Eigen/")


# Options de compilation
set(CMAKE_CXX_FLAGS "-Wall")

# Options spécifiques à ce projet.
# Commentez ou décommentez la ligne selon le type de build désiré.
# set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_BUILD_TYPE "Release")

# Nom de l'exécutable
set(TARGET_NAME "run")

# Liste des fichiers à compiler
set(PROJECT_COMPILATION_FILE_LIST mesh.cpp DataFile.cpp Geometry.cpp sourceAndBord.cpp Assemblage.cpp Resolution.cpp Solveur.cpp main.cc)

# Assemblage des fichiers compilés avec l'exécutable
add_executable(${TARGET_NAME} ${PROJECT_COMPILATION_FILE_LIST})

# Inclusion du répertoire de la bibliothèque Eigen
target_include_directories(${TARGET_NAME} PUBLIC "${LIBRARY_PATH_EIGEN}")

# Fin des options spécifiques à ce projet.

