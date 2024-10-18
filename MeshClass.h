//
// Unstructured mesh class
//

#ifndef NMIF_MESHCLASS_H
#define NMIF_MESHCLASS_H
#include <iostream>
#include <tuple>
#include <string>
#include <fstream>
#include "linAlg.h"
#include <numeric>

class pointClass{
public:
    std::string label_;
    int pointDimension_;
    std::vector<double> coordinates_;
    size_t id_;

    pointClass(const std::string &label, const int &numberData, const std::vector<double> &v, const size_t &id);

};

class faceClass{
public:
    std::string label_;
    int numVertices_;
    std::vector<int> verticesID_;
    std::vector<pointClass> verticesObjects_;
    size_t id_;
    size_t cellOwnerID_=-1;
    std::vector<double> faceCentroidCoordinates_;
    std::vector<double> areaVector_;
    std::vector<size_t> faceCellNeighboursIDs_;

    faceClass(const std::string &label, const int &numberData, const std::vector<double> &v, const size_t &id);

    void assignPointsObjectsVector(const std::vector<pointClass> &v_point_objs);

    std::vector<double> faceCentroid();
    std::vector<double> areaVector();

private:
    void computeFaceCentroid();

    void computeAreaVector();
};

class cellClass{
public:
    std::string label_;
    int numFaces_;
    std::vector<int> facesID_;
    std::vector<faceClass> facesObjects_;
    std::vector<double> cellCentroidCoordinates_;
    double cellVolume_;
    size_t id_;

    cellClass(const std::string &label, const int &numberData, const std::vector<double> &v, const size_t &id);

    void assignFacesObjectsVector(const std::vector<faceClass> &v_faces_objs);

    std::vector<double> cellCentroid();

    double cellVolume();
private:
    void computeCellCentroid();
    void computeCellVolume();
};

class patchClass{
public:
    std::string label_;
    int numFaces_;
    std::vector<int> facesID_;
    std::vector<faceClass> facesObjects_;
    size_t id_;

    patchClass(const std::string &label, const int &numberData, const std::vector<double> &v, const size_t &id);

    void assignFacesObjectsVector(const std::vector<faceClass> &v_face_objs);



};


class MeshClass {

public:
    unsigned int numMeshPoints_, numMeshFaces_, numMeshCells_, numMeshBoundaryPatches_;
    std::vector<pointClass> listOfPoints_;
    std::vector<faceClass> listOfFaces_;
    std::vector<cellClass> listOfCells_;
    std::vector<patchClass> listOfPatches_;

    MeshClass(const std::string &pointsFileName,
              const std::string &facesFileName,
              const std::string &cellsFileName,
              const std::string &patchesFileName);

    void meshInfo(unsigned int &numMeshPoints,
                  unsigned int &numMeshFaces,
                  unsigned int &numMeshCells,
                  unsigned int &numMeshBoundaryPatches);

    std::vector<double> get_ithFaceCentroidCoordinates(const size_t faceID);
    std::vector<double> get_ithFaceAreaVector(const size_t faceID);

    std::vector<double> get_ithCellCentroidCoordinates(const size_t cellID);
    double get_ithCellVolume(const size_t cellID);

    size_t get_ithFaceOwner(const size_t faceID);

    std::vector<size_t> get_ithFaceNeighbours(const size_t faceID);

    std::vector<size_t> get_ithCellNeighbours(const size_t cellID);

    void get_ithPatchCellIDs(const size_t patchID, std::vector<size_t> &patchFaces, std::vector<size_t> &patchCells);

private:
    // Private methods
    template<class readObjectClass>
    int readInputFile(const std::string &filePath, std::vector<readObjectClass> &listOfObjects);

    static size_t split(const std::string &txt, std::vector<std::string> &strs, char ch);
};


#endif //NMIF_MESHCLASS_H
