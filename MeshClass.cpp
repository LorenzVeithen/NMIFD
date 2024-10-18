//
// Created by Lorenz Veithen on 13/10/2024.
//

#include "MeshClass.h"


/**
 * MeshClass
 */
MeshClass::MeshClass(const std::string &pointsFileName,
                     const std::string &facesFileName,
                     const std::string &cellsFileName,
                     const std::string &patchesFileName) {

    // Read all files and find number of elements for each
    readInputFile(pointsFileName, listOfPoints_);
    numMeshPoints_ = std::size(listOfPoints_);
    
    readInputFile(facesFileName, listOfFaces_);
    numMeshFaces_ = std::size(listOfFaces_);
    
    readInputFile(cellsFileName, listOfCells_);
    numMeshCells_ = std::size(listOfCells_);
    
    readInputFile(patchesFileName, listOfPatches_);
    numMeshBoundaryPatches_ = std::size(listOfPatches_);

    // Get owner and neighbour of each face
    for (cellClass &cellObj: listOfCells_) {
        for (const int &i: cellObj.facesID_) {
            // Make face aware that it has a cell owner and cell neighbours
            if (listOfFaces_[i].cellOwnerID_ == -1) {
                listOfFaces_[i].cellOwnerID_ = cellObj.id_;
            } else {
                listOfFaces_[i].faceCellNeighboursIDs_.push_back(cellObj.id_);
            }
        }
    }

    // Link objects to their parent classes
    for (faceClass &faceObj: listOfFaces_) {
        std::vector<pointClass> selectedObjects;
        int counter = 0;
        for (const int &i: faceObj.verticesID_) {
            pointClass selectedObject = listOfPoints_[i];
            selectedObjects.push_back(selectedObject);
            counter++;
            if (selectedObject.id_ != i) {
                std::cerr << "Error in assigning point objects, requested and received ID do not match: " << std::endl;
                std::cerr << selectedObject.id_ << "!=" << i << std::endl;
                exit(1);
            }
        }
        faceObj.assignPointsObjectsVector(selectedObjects);
    }

    for (cellClass &cellObj: listOfCells_) {
        std::vector<faceClass> selectedObjects;
        int counter = 0;
        for (const int &i: cellObj.facesID_) {
            faceClass selectedObject = listOfFaces_[i];
            selectedObjects.push_back(selectedObject);
            counter++;
            if (selectedObject.id_ != i) {
                std::cerr << "Error in assigning face objects, requested and received ID do not match: " << std::endl;
                std::cerr << selectedObject.id_ << "!=" << i << std::endl;
                exit(1);
            }
        }
        cellObj.assignFacesObjectsVector(selectedObjects);
    }
    
    for (patchClass &patchObj: listOfPatches_) {
        std::vector<faceClass> selectedObjects;
        int counter = 0;
        for (const int &i: patchObj.facesID_) {
            faceClass selectedObject = listOfFaces_[i];
            selectedObjects.push_back(selectedObject);
            counter++;
            if (selectedObject.id_ != i) {
                std::cerr << "Error in assigning face objects, requested and received ID do not match: " << std::endl;
                std::cerr << selectedObject.id_ << "!=" << i << std::endl;
                exit(1);
            }
        }
        patchObj.assignFacesObjectsVector(selectedObjects);
    }
}

template<class readObjectClass>
int MeshClass::readInputFile(const std::string &filePath, std::vector<readObjectClass> &listOfObjects) {
    std::ifstream inFile(filePath);
    if (!inFile) {
        std::cerr << "File could not be opened!" << std::endl;
        return 1;
    }
    int numItems;
    std::string line;

    std::string nextLabel = "";
    signed int nextNumberData = -1;
    std::vector<double> v_double;

    int blockLevel = 0;
    unsigned int itemCounter = 0;
    while (std::getline(inFile, line)) {
        // Strip down the line from the tab character
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());

        for (unsigned int i = 0; i < std::size(line); i++) {
            char c = line[i];
            if (isdigit(c) && blockLevel == 0) {
                // Get the number of elements in this section
                numItems = std::stoi(line);
            } else if (c == '(') {
                blockLevel++;
            } else if ((blockLevel == 2) && (c == ')')) {
                if ((nextNumberData != std::size(v_double)) && nextNumberData != -1) {
                    std::cerr << "Number of data points flag and number of received data points do not match"
                              << std::endl;
                }
                readObjectClass obj(nextLabel, nextNumberData, v_double, itemCounter);
                listOfObjects.push_back(obj);
                itemCounter++;

                // Re-initialise variables for next item
                blockLevel--;
                nextLabel = "";
                nextNumberData = -1;
                v_double.clear();
            } else if ((blockLevel == 1) && c == ' ') {
                continue;
            } else if ((blockLevel == 1) && isalpha(c)) {
                // Must be the label of the next section
                nextLabel += c;
            } else if ((blockLevel == 1) && isdigit(c)) {
                nextNumberData = std::stoi(std::string{c});
            } else if (blockLevel == 2) {
                // Data is here, find where the string ends
                size_t endData = line.find(')', i);
                size_t endInterval = (endData == std::string::npos) ? std::size(line) : endData;

                std::string dataOnly = line.substr(i, endInterval - i);
                std::vector<std::string> v_str;
                split(dataOnly, v_str, ' ');
                for (const std::string &v_str_j: v_str) {
                    if (!v_str_j.empty()) v_double.push_back(std::stod(v_str_j));
                }
                i = endInterval - 1;  // Skip the indices covered here
            }
        }
    }
    if (itemCounter != numItems) {
        std::cerr << "File level number of data points flag and number of received data points do not match"
                  << std::endl;
        std::cerr << itemCounter << "!=" << numItems << std::endl;
    }
    return 0;
}

size_t MeshClass::split(const std::string &txt, std::vector<std::string> &strs, char ch) {
    /**
     * From: https://stackoverflow.com/questions/5888022/split-string-by-single-spaces
     */
    size_t pos = txt.find(ch);
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while (pos != std::string::npos) {
        strs.push_back(txt.substr(initialPos, pos - initialPos));
        initialPos = pos + 1;

        pos = txt.find(ch, initialPos);
    }

    // Add the last one
    strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));

    return strs.size();
}

void MeshClass::meshInfo(unsigned int &numMeshPoints,
                         unsigned int &numMeshFaces,
                         unsigned int &numMeshCells,
                         unsigned int &numMeshBoundaryPatches) {
    /**
     * Return mesh information.
     *
     * Returns the number of mesh points, faces, cells and boundary patches.
     *
     * @return tuple containing the number of mesh points, faces, cells, and boundary patches.
     */
    numMeshPoints = numMeshPoints_;
    numMeshFaces = numMeshFaces_;
    numMeshCells = numMeshCells_;
    numMeshBoundaryPatches = numMeshBoundaryPatches_;
}

std::vector<double> MeshClass::get_ithFaceCentroidCoordinates(const size_t faceID){
    if (faceID > numMeshFaces_){
        std::cerr << "No face with id: " << faceID << std::endl;
    }
    return listOfFaces_[faceID].faceCentroidCoordinates_;
}
std::vector<double> MeshClass::get_ithFaceAreaVector(const size_t faceID){
    if (faceID > numMeshFaces_){
        std::cerr << "No face with id: " << faceID << std::endl;
    }
    return listOfFaces_[faceID].areaVector_;
}

std::vector<double> MeshClass::get_ithCellCentroidCoordinates(const size_t cellID){
    if (cellID > numMeshCells_){
        std::cerr << "No cell with id: " << cellID << std::endl;
    }
    return listOfCells_[cellID].cellCentroidCoordinates_;
}

double MeshClass::get_ithCellVolume(const size_t cellID){
    if (cellID > numMeshCells_){
        std::cerr << "No cell with id: " << cellID << std::endl;
    }
    return listOfCells_[cellID].cellVolume_;
}

size_t MeshClass::get_ithFaceOwner(const size_t faceID){
    if (faceID > numMeshFaces_){
        std::cerr << "No face with id: " << faceID << std::endl;
    }
    return listOfFaces_[faceID].cellOwnerID_;
}

std::vector<size_t> MeshClass::get_ithFaceNeighbours(const size_t faceID){
    if (faceID > numMeshFaces_){
        std::cerr << "No face with id: " << faceID << std::endl;
    }
    return listOfFaces_[faceID].faceCellNeighboursIDs_;
}

std::vector<size_t> MeshClass::get_ithCellNeighbours(const size_t cellID){
    std::vector<size_t> neighbouringCells;
    for (const faceClass &fObj : listOfCells_[cellID].facesObjects_){
        if (!fObj.faceCellNeighboursIDs_.empty()){
            // If the face has a neighbour, add it to the list
            for (const size_t neighbourCell : fObj.faceCellNeighboursIDs_){
                neighbouringCells.push_back(neighbourCell);
            }
        }
    }
    return neighbouringCells;
}

void MeshClass::get_ithPatchCellIDs(const size_t patchID, std::vector<size_t> &patchFaces, std::vector<size_t> &patchCells){
    if (!patchFaces.empty() || !patchCells.empty()){
        std::cerr << "Referenced vectors in get_ithPatchCellIDs should be empty" << std::endl;
    }

    for (const faceClass &fObj : listOfPatches_[patchID].facesObjects_){
        patchFaces.push_back(fObj.id_);
        patchCells.push_back(fObj.cellOwnerID_);
    }
}
/**
 * Point class
 */
pointClass::pointClass(const std::string &label, const int &numberData, const std::vector<double> &v,
                       const size_t &id) {
    label_ = label;
    pointDimension_ = numberData;
    coordinates_ = v;
    id_ = id;
}

/**
 * Face class
 */
faceClass::faceClass(const std::string &label, const int &numberData, const std::vector<double> &v,
                     const size_t &id) {
    label_ = label;
    numVertices_ = numberData;
    for (double num: v) {
        verticesID_.push_back(static_cast<int>(std::round(num)));
    }
    id_ = id;
}

void faceClass::assignPointsObjectsVector(const std::vector<pointClass> &v_point_objs) {
    for (const pointClass &pointObj: v_point_objs) {
        verticesObjects_.push_back(pointObj);
    }
    numVertices_ = verticesObjects_.size();

    // Compute the face centroid and store its value
    computeFaceCentroid();

    // Compute face area vector
    computeAreaVector();
}

void faceClass::computeFaceCentroid() {
    size_t coordinates_size = std::size(verticesObjects_[0].coordinates_);
    std::vector<double> centroid(coordinates_size);
    for (const pointClass &pointObj: verticesObjects_) {
        for (size_t i = 0; i < coordinates_size; i++) {
            centroid[i] += pointObj.coordinates_[i] / numVertices_;
        }
    }
    faceCentroidCoordinates_ = centroid;
}

std::vector<double> faceClass::faceCentroid() {
    return faceCentroidCoordinates_;
}

void faceClass::computeAreaVector() {
    std::vector<double> areaVector(verticesObjects_[0].coordinates_.size());
    // Sum the area vector of the different triangles making up the face
    for (int i = 0; i < numVertices_; i++) {
        pointClass curr_point = verticesObjects_[i];
        pointClass next_point = (i == (numVertices_-1)) ? verticesObjects_[0] : verticesObjects_[i + 1];     // If last point, comes back to the first one


        // Vectors relative to the face centroid
        std::vector<double> curr_face_centroid_to_point_vector(std::size(curr_point.coordinates_));
        std::vector<double> next_face_centroid_to_point_vector(std::size(curr_point.coordinates_));
        for (int j = 0; j < std::size(curr_point.coordinates_); j++) {
            curr_face_centroid_to_point_vector[j] = curr_point.coordinates_[j] - faceCentroidCoordinates_[j];
            next_face_centroid_to_point_vector[j] = next_point.coordinates_[j] - faceCentroidCoordinates_[j];
        }

        // Cross product of centroid to point vectors
        std::vector<double> currentPartialAreaVector(std::size(curr_point.coordinates_));
        crossProduct(curr_face_centroid_to_point_vector, next_face_centroid_to_point_vector, currentPartialAreaVector);

        // Sum the partial area vectors
        for (int j = 0; j < std::size(curr_point.coordinates_); j++) {
            areaVector[j] += currentPartialAreaVector[j]/2.;
        }
    }
    areaVector_ = areaVector;
}

std::vector<double> faceClass::areaVector() {
    return areaVector_;
}

/**
 * Cell class
 */
cellClass::cellClass(const std::string &label, const int &numberData, const std::vector<double> &v, const size_t &id) {
    label_ = label;
    numFaces_ = numberData;
    for (double num: v) {
        facesID_.push_back(static_cast<int>(std::round(num)));
    }
    id_ = id;
}

void cellClass::assignFacesObjectsVector(const std::vector<faceClass> &v_faces_objs) {
    for (const faceClass &faceObj: v_faces_objs) {
        facesObjects_.push_back(faceObj);
    }
    numFaces_ = std::size(facesObjects_);

    // Compute cell centroid
    computeCellCentroid();

    // Compute cell volume
    computeCellVolume();
}

void cellClass::computeCellCentroid() {
    int num_vertices = 0;
    std::vector<double> centroidCoordinates(facesObjects_[0].faceCentroidCoordinates_.size());
    std::vector<size_t> pointsUsed;
    for (const faceClass &fObj: facesObjects_) {
        for (const pointClass &pObj: fObj.verticesObjects_){
            if (std::find(pointsUsed.begin(), pointsUsed.end(), pObj.id_) != pointsUsed.end()) {
               continue;                   // The ID was already taken into account
            }
            num_vertices += 1;
            std::vector<double> currentPointCoordinates = pObj.coordinates_;
            for (int i=0; i<std::size(centroidCoordinates); i++){
                centroidCoordinates[i] += currentPointCoordinates[i];
            }
            pointsUsed.push_back(pObj.id_);
        }
    }
    for (double & centroidCoordinate : centroidCoordinates){
        centroidCoordinate /= num_vertices;
    }
    cellCentroidCoordinates_ = centroidCoordinates;
}

std::vector<double> cellClass::cellCentroid(){
    return cellCentroidCoordinates_;
}

void cellClass::computeCellVolume(){
    double cellVolume = 0.;
    for (faceClass &faceObj : facesObjects_){
        std::vector<double> currentFaceCentroid = faceObj.faceCentroidCoordinates_;

        // Face centroid position vector relative to cell centroid
        std::vector<double> cell_centroid_to_face_centroid_vector(std::size(currentFaceCentroid));
        for (int j = 0; j < std::size(currentFaceCentroid); j++) {
            cell_centroid_to_face_centroid_vector[j] = currentFaceCentroid[j] - cellCentroidCoordinates_[j];
        }
        for (int i=0; i<faceObj.numVertices_; i++){
            pointClass curr_point = faceObj.verticesObjects_[i];
            pointClass next_point = (i == (faceObj.numVertices_-1)) ? faceObj.verticesObjects_[0] : faceObj.verticesObjects_[i + 1];     // If last point, comes back to the first one

            // Vectors relative to cell centroid
            std::vector<double> curr_face_centroid_to_point_vector(std::size(currentFaceCentroid));
            std::vector<double> next_face_centroid_to_point_vector(std::size(currentFaceCentroid));
            for (int j = 0; j < std::size(curr_point.coordinates_); j++) {
                curr_face_centroid_to_point_vector[j] = curr_point.coordinates_[j] - cellCentroidCoordinates_[j];
                next_face_centroid_to_point_vector[j] = next_point.coordinates_[j] - cellCentroidCoordinates_[j];
            }

            // Compute tetrahedral volume
            std::vector<double> cross_product_point_vectors(std::size(currentFaceCentroid));
            crossProduct(curr_face_centroid_to_point_vector, next_face_centroid_to_point_vector, cross_product_point_vectors);
            double current_tetrahedral_volume = (1./6.) * std::inner_product(cross_product_point_vectors.begin(),
                                                                  cross_product_point_vectors.end(),
                                                                  cell_centroid_to_face_centroid_vector.begin(),
                                                                  0.0);
            cellVolume += current_tetrahedral_volume;
        }
    }
    cellVolume_ = cellVolume;
}

double cellClass::cellVolume(){
    return cellVolume_;
}

/**
 * Patch class
 */
patchClass::patchClass(const std::string &label, const int &numberData, const std::vector<double> &v,
                       const size_t &id) {
    label_ = label;
    numFaces_ = numberData;
    for (double num: v) {
        facesID_.push_back(static_cast<int>(std::round(num)));
    }
    id_ = id;
}

void patchClass::assignFacesObjectsVector(const std::vector<faceClass> &v_face_objs) {
    for (const faceClass &faceObj: v_face_objs) {
        facesObjects_.push_back(faceObj);
    }
}