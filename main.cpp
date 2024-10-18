#include <iostream>
#include "MeshClass.h"

int main() {
    MeshClass mesh("/Users/lorenz_veithen/Desktop/Education/04-MPhil/01_Cambridge/01_MichaelmasTerm/Continuum/02_NumericalMethods_IncompressibleFluidDynamics/NMIF/inputFiles/points.txt",
                   "/Users/lorenz_veithen/Desktop/Education/04-MPhil/01_Cambridge/01_MichaelmasTerm/Continuum/02_NumericalMethods_IncompressibleFluidDynamics/NMIF/inputFiles/faces.txt",
                   "/Users/lorenz_veithen/Desktop/Education/04-MPhil/01_Cambridge/01_MichaelmasTerm/Continuum/02_NumericalMethods_IncompressibleFluidDynamics/NMIF/inputFiles/cells.txt",
                   "/Users/lorenz_veithen/Desktop/Education/04-MPhil/01_Cambridge/01_MichaelmasTerm/Continuum/02_NumericalMethods_IncompressibleFluidDynamics/NMIF/inputFiles/boundary.txt");

    unsigned int r1, r2, r3, r4;
    mesh.meshInfo(r1, r2, r3, r4);
    std::cout << r1 << std::endl;
    std::cout << r2 << std::endl;
    std::cout << r3 << std::endl;
    std::cout << r4 << std::endl;


    int id_test = 0;
    std::cout << "----- Face centroid -----" << std::endl;
    std::vector<double> faceCentroidTest = mesh.get_ithFaceCentroidCoordinates(id_test);
    for (double coord : faceCentroidTest){
        std::cout << coord << std::endl;
    }

    std::cout << "----- Face area -----" << std::endl;
    std::vector<double> faceAreaTest = mesh.get_ithFaceAreaVector(id_test);
    for (double coord : faceAreaTest){
        std::cout << coord << std::endl;
    }

    std::cout << "----- Cell centroid -----" << std::endl;
    std::vector<double> cellCentroidTest = mesh.get_ithCellCentroidCoordinates(id_test);
    for (double coord : cellCentroidTest){
        std::cout << coord << std::endl;
    }

    std::cout << "----- Cell volume -----" << std::endl;
    double cellVolumeTest = mesh.get_ithCellVolume(id_test);
    std::cout << cellVolumeTest << std::endl;

    std::cout << "----- Face owner -----" << std::endl;
    size_t ithFaceOwner = mesh.get_ithFaceOwner(id_test);
    std::cout << ithFaceOwner << std::endl;

    std::cout << "----- Face neighbours -----" << std::endl;
    std::vector<size_t> ithFaceNeighbours = mesh.get_ithFaceNeighbours(id_test);
    for (size_t nb : ithFaceNeighbours){
        std::cout << nb << std::endl;
    }

    std::cout << "----- Cell neighbours -----" << std::endl;
    std::vector<size_t> ithCellNeighbours = mesh.get_ithCellNeighbours(id_test);
    for (size_t nb : ithCellNeighbours){
        std::cout << nb << std::endl;
    }

    std::cout << "----- Patch Cell -----" << std::endl;
    std::vector<size_t> patchFaceID;
    std::vector<size_t> patchCellID;
    mesh.get_ithPatchCellIDs(id_test, patchFaceID, patchCellID);
    for (size_t i=0; i<patchFaceID.size(); i++){
        std::cout << "Face ID: " << patchFaceID[i] << ", Cell ID: " << patchCellID[i] << std::endl;
    }
}

