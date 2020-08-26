/*
 * utils.cpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#include "couplingUtils.hpp"
#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

namespace coupling {
///*!
//    Build a scaled mesh from the input mesh, scaling the input mesh vertex coordinates by the radius factor
//
//    \param[in] unitRadiusSphereMesh the mesh to be scaled
//    \param[in] radius the scaling factor
//    \return the scaled mesh
//*/
//std::unique_ptr<SurfUnstructured> scale(const SurfUnstructured & unitRadiusSphereMesh, double radius) {
//
//    std::unique_ptr<SurfUnstructured> mesh = PatchKernel::clone(&unitRadiusSphereMesh);
//
//    PiercedVector<Vertex> & vertices = mesh->getVertices();
//
//    for(auto & v : vertices) {
//        v.setCoords(v.getCoords()*radius);
//    }
//
//    return mesh;
//}

/*!
    Interpolate data from the "fromMesh" to the "toMesh"

    \param[in] fromMesh the interpolation origin mesh
    \param[in] fromData the data of the origin to be interpolated
    \param[in] toMesh the interpolation destination mesh
    \param[in] toData the interpolated data on the destination mesh
    \return the scaled mesh
*/
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedVector<double> * fromData, SurfUnstructured * toMesh, PiercedVector<double> * toData){

    log::cout() << "Build Tree" << std::endl;

    SurfaceSkdTree fromTree(fromMesh);
    log::cout() << "Tree declared" << std::endl;
    fromTree.build(1);

    log::cout() << "Tree built" << std::endl;


    const PiercedVector<Vertex> & toVertices = toMesh->getVertices();
    const PiercedVector<Vertex> & fromVertices = fromMesh->getVertices();
    PiercedVector<double>::iterator toDataItEnd = toData->end();
    for (const Vertex & v : toVertices){

        // Find closest triangle of Mesh from for each vertex of Mesh to
        darray3 x = v.getCoords();
        long id;
        double dist;
        fromTree.findPointClosestCell(x, &id, &dist);

        // Recover vertices of triangle of Mesh 1
        ConstProxyVector<long> vIds = fromMesh->getCell(id).getVertexIds();

        // Compute barycentric coordinates of the projected vertex of Interpolation Mesh on the closest triangle of Reference Mesh
        darray3 lambda;
        darray3 xP = bitpit::CGElem::projectPointTriangle(x, fromVertices[vIds[0]].getCoords(), fromVertices[vIds[1]].getCoords(), fromVertices[vIds[2]].getCoords(), lambda);
        BITPIT_UNUSED(xP);
        // Interpolate data of Mesh 1 on vertex of Mesh 2
        double data = 0.;
        data += (*fromData)[vIds[0]]*lambda[0];
        data += (*fromData)[vIds[1]]*lambda[1];
        data += (*fromData)[vIds[2]]*lambda[2];

        // Update data of Interpolation Mesh
        PiercedVector<double>::iterator toDataIt = toData->find(v.getId());

        if( toDataIt == toDataItEnd) {
            toData->insert(v.getId(), data);
            throw std::runtime_error("Warning! New element inserted!");
        }
        else {
            *toDataIt = data;
        }
    }


};

/*!
    Interpolate cell collocated data from the "fromMesh" to the "toMesh"

    \param[in] fromMesh the interpolation origin mesh
    \param[in] fromData the data of the origin to be interpolated (PiercedStorage)
    \param[in] toMesh the interpolation destination mesh
    \param[in] toData the interpolated data on the destination mesh (PiercedStorage)
    \return the scaled mesh
*/
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedStorage<double,long> * fromData, SurfUnstructured * toMesh, PiercedStorage<double,long> * toData, int fieldIndex){

    assert(fieldIndex < fromData->getFieldCount() && fieldIndex < toData->getFieldCount() && fieldIndex >= 0);
    log::cout() << "Interpolating... " << std::endl;

    SurfaceSkdTree fromTree(fromMesh, false);//costruisce i ghost
    fromTree.build(1);

    //Collect cell center into array to be passed to parallel find PointClosest cell
    darray3 *toCellCenters = new darray3[toMesh->getInternalCount()];
    long *toCellIds = new long[toMesh->getInternalCount()];
    long toCellId = 0, cellCounter = 0;
    PiercedVector<Cell>::iterator beginInternalCells = toMesh->internalBegin();
    PiercedVector<Cell>::iterator endInternalCells = toMesh->internalEnd();
    for(PiercedVector<Cell>::iterator itCell = beginInternalCells; itCell != endInternalCells; ++itCell) {
        toCellId = itCell->getId();
        toCellCenters[cellCounter] = toMesh->evalCellCentroid(toCellId);
        toCellIds[cellCounter] = toCellId;
        ++cellCounter;
    }

    //Project cell centers to fromMesh getting the rank and the cell id owning the projection
    //Organize projection cell ids and centers by rank for communication
    long *fromClosestCellIds = new long[toMesh->getInternalCount()];
    int *fromClosestCellRanks = new int[toMesh->getInternalCount()];
    double *fromClosestCellDist = new double[toMesh->getInternalCount()];
    fromTree.findPointClosestGlobalCell(toMesh->getInternalCount(),toCellCenters,fromClosestCellIds,fromClosestCellRanks,fromClosestCellDist);

    //Counting points per rank
    std::vector<long> nofPointsPerRank(toMesh->getProcessorCount(),0);
    for(int i = 0; i < toMesh->getInternalCount(); ++i) {
        ++nofPointsPerRank[fromClosestCellRanks[i]];
    }
    //Map collecting (projectionId, originId, originCellCenter) per rank
    std::vector<std::vector<InterpolationInfo>> interpolationInfoPerRank;
    interpolationInfoPerRank.resize(toMesh->getProcessorCount());
    for(int r = 0; r < toMesh->getProcessorCount(); ++r) {
        interpolationInfoPerRank[r].reserve(nofPointsPerRank[r]);
    }
    for(int i = 0; i < toMesh->getInternalCount(); ++i) {
        InterpolationInfo buffer;
        buffer.originId = toCellIds[i];
        buffer.originCellCenter = toCellCenters[i];
        buffer.projectionId = fromClosestCellIds[i];
        interpolationInfoPerRank[fromClosestCellRanks[i]].push_back(buffer);
    }

    //Local interpolation points ready to be interpolated
    const std::vector<InterpolationInfo> & localInterpolationPoints = interpolationInfoPerRank[toMesh->getRank()];

    //Interpolation points from other ranks
    std::unordered_map<int,std::vector<InterpolationInfo>> otherRankIntepolationPoints;

    //Communicate demanding rank and points InterpolationInfo for non-local interpolations
    DataCommunicator interpolationComm(toMesh->getCommunicator());
    // demanding rank (int), nof points, InterpolationInfo (2 long and 3 doubles) per point
    for(int p = 0; p < toMesh->getProcessorCount(); ++p) {
        if(p==toMesh->getRank() || nofPointsPerRank[p] == 0) {
            continue;
        }
        std::size_t buffSize = sizeof(int) + sizeof(long) + (2*sizeof(long) * 3*sizeof(double))*nofPointsPerRank[p];
        interpolationComm.setSend(p,buffSize);
        SendBuffer &sendBuffer = interpolationComm.getSendBuffer(p);
        sendBuffer << toMesh->getRank();
        sendBuffer << nofPointsPerRank[p];
        for(long i = 0; i < nofPointsPerRank[p]; ++i){
            sendBuffer << interpolationInfoPerRank[p][i].originId;
            sendBuffer << interpolationInfoPerRank[p][i].projectionId;
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[0];
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[1];
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[2];
        }
    }

    interpolationComm.discoverRecvs();
    interpolationComm.startAllRecvs();
    interpolationComm.startAllSends();

    //Interpolate locals
    darray3 toCellCenter, neighCellCenter;
    long fromCellId;
    double weightSum, interpVal,cellCentersDist,weight;
    std::vector<long> neighs;
    for(const InterpolationInfo & info : localInterpolationPoints) {
        toCellId = info.originId;
        toCellCenter = info.originCellCenter;
        fromCellId = info.projectionId;
        neighs.clear();
        fromMesh->findCellNeighs(fromCellId,&neighs);
        //Interpolation
        weightSum = 0.0;
        interpVal = 0.0;
        //fromCell contribution
        neighCellCenter = fromMesh->evalCellCentroid(fromCellId);
        cellCentersDist = norm2(neighCellCenter-toCellCenter);
        weight = 1.0 / (cellCentersDist*cellCentersDist);
        weightSum += weight;
        interpVal += weight * fromData->at(fromCellId, fieldIndex);
        for(const long & neigh : neighs) {
            neighCellCenter = fromMesh->evalCellCentroid(neigh);
            cellCentersDist = norm2(neighCellCenter-toCellCenter);
            weight = 1.0 / (cellCentersDist*cellCentersDist);
            weightSum += weight;
            interpVal += weight * fromData->at(neigh,fieldIndex);
        }
        interpVal /= weightSum;

        //toData->set(toCellId,interpVal);
        toData->set(toCellId, fieldIndex, interpVal);
    }

    std::vector<int> recvRanks = interpolationComm.getRecvRanks();
    for(int rank : recvRanks) {
        interpolationComm.waitRecv(rank);
        RecvBuffer & recvBuffer = interpolationComm.getRecvBuffer(rank);
        int demandingRank;
        long nofPoints;
        InterpolationInfo pointInterpolationInfo;
        recvBuffer >> demandingRank;
        recvBuffer >> nofPoints;
        otherRankIntepolationPoints[demandingRank] = std::vector<InterpolationInfo>(nofPoints,pointInterpolationInfo);
        for(long p = 0; p < nofPoints; ++p) {
            recvBuffer >> pointInterpolationInfo.originId;
            recvBuffer >> pointInterpolationInfo.projectionId;
            recvBuffer >> pointInterpolationInfo.originCellCenter[0];
            recvBuffer >> pointInterpolationInfo.originCellCenter[1];
            recvBuffer >> pointInterpolationInfo.originCellCenter[2];
            otherRankIntepolationPoints[demandingRank][p] = pointInterpolationInfo;
        }
    }

    // Interpolate others and organize results by rank demanding interpolation
    std::unordered_map<int,std::vector<InterpolatedInfo>> otherInterpolatedValues;
    InterpolatedInfo emptyInterpolatedInfo;
    for(const auto & rankInfo : otherRankIntepolationPoints) {
        otherInterpolatedValues[rankInfo.first] = std::vector<InterpolatedInfo>(rankInfo.second.size(),emptyInterpolatedInfo);
    }

    for(const auto & rankInfo : otherRankIntepolationPoints) {
        long pointCounter = 0;
        for(const InterpolationInfo & info : rankInfo.second) {
            toCellId = info.originId;
            toCellCenter = info.originCellCenter;
            fromCellId = info.projectionId;
            neighs.clear();
            fromMesh->findCellNeighs(fromCellId,&neighs);
            //Interpolation
            weightSum = 0.0;
            interpVal = 0.0;
            //fromCell contribution
            neighCellCenter = fromMesh->evalCellCentroid(fromCellId);
            cellCentersDist = norm2(neighCellCenter-toCellCenter);
            weight = 1.0 / (cellCentersDist*cellCentersDist);
            weightSum += weight;
            interpVal += weight * fromData->at(fromCellId, fieldIndex);
            for(const long & neigh : neighs) {
                neighCellCenter = fromMesh->evalCellCentroid(neigh);
                cellCentersDist = norm2(neighCellCenter-toCellCenter);
                weight = 1.0 / (cellCentersDist*cellCentersDist);
                weightSum += weight;
                interpVal += weight * fromData->at(neigh,fieldIndex);
            }
            interpVal /= weightSum;
            InterpolatedInfo interpolatedInfo;
            interpolatedInfo.originId = toCellId;
            interpolatedInfo.value = interpVal;
            otherInterpolatedValues[rankInfo.first][pointCounter] = interpolatedInfo;
            ++pointCounter;
        }
    }

    //Communicate Interpolated values -
    //buffer = (originId(long) + interpolated value(double))*nofInterpolatedValues + nofInterpolatedValues(long)
    DataCommunicator interpolatedComm(toMesh->getCommunicator());
    for(const auto & rankInfo : otherInterpolatedValues) {
        std::size_t buffSize = sizeof(long) + (sizeof(long) + sizeof(double)) * rankInfo.second.size();
        interpolatedComm.setSend(rankInfo.first,buffSize);
        SendBuffer &sendBuffer = interpolatedComm.getSendBuffer(rankInfo.first);
        sendBuffer << rankInfo.second.size();
        for(size_t p = 0; p < rankInfo.second.size(); ++p) {
            sendBuffer << rankInfo.second[p].originId;
            sendBuffer << rankInfo.second[p].value;
        }
    }
    interpolatedComm.discoverRecvs();
    interpolatedComm.startAllRecvs();
    interpolatedComm.startAllSends();

//    std::stringstream ss;
//    ss << "interpolates_" << toMesh->getRank() << ".txt";
//    std::ofstream out(ss.str().c_str());
    recvRanks.clear();
    recvRanks = interpolatedComm.getRecvRanks();
    double value;
    for(int rank : recvRanks) {
        interpolatedComm.waitRecv(rank);
        RecvBuffer & recvBuffer = interpolatedComm.getRecvBuffer(rank);
        long nofPoints;
        recvBuffer >> nofPoints;
        for(long p = 0; p < nofPoints; ++p) {
            recvBuffer >> toCellId;
            recvBuffer >> value;
            toData->set(toCellId, fieldIndex, value);
//            out << toCellId << " " << value << " " << rank << std::endl;
        }
    }
//    out.close();

    delete [] toCellCenters;
    delete [] toCellIds;
    delete [] fromClosestCellIds;
    delete [] fromClosestCellRanks;
    delete [] fromClosestCellDist;
};


/*!
    Interpolate cell collocated data from the "fromMesh" to the "toMesh"

    \param[in] fromMesh the interpolation origin mesh
    \param[in] fromData the data of the origin to be interpolated (PiercedStorage)
    \param[in] toMesh the interpolation destination mesh
    \param[in] toData the interpolated data on the destination mesh (PiercedStorage)
    \return the scaled mesh
*/
void interpolateFromToMatrix(Mat * interpolationMatrix, SurfUnstructured * fromMesh, PatchNumberingInfo * fromNumberingInfo, SurfUnstructured * toMesh, PatchNumberingInfo * toNumberingInfo){

    log::cout() << "Compute Interpolation Matrix... " << std::endl;

    std::vector<MatrixRow> rowCollection(toMesh->getInternalCount());
    int rowCounter = 0;

    SurfaceSkdTree fromTree(fromMesh, false);//costruisce i ghost
    fromTree.build(1);

    //Collect cell center into array to be passed to parallel find PointClosest cell
    darray3 *toCellCenters = new darray3[toMesh->getInternalCount()];
    long *toCellIds = new long[toMesh->getInternalCount()];
    long toCellId = 0, cellCounter = 0;
    PiercedVector<Cell>::iterator beginInternalCells = toMesh->internalBegin();
    PiercedVector<Cell>::iterator endInternalCells = toMesh->internalEnd();
    for(PiercedVector<Cell>::iterator itCell = beginInternalCells; itCell != endInternalCells; ++itCell) {
        toCellId = itCell->getId();
        toCellCenters[cellCounter] = toMesh->evalCellCentroid(toCellId);
        toCellIds[cellCounter] = toCellId;
        ++cellCounter;
    }

    //Project cell centers to fromMesh getting the rank and the cell id owning the projection
    //Organize projection cell ids and centers by rank for communication
    long *fromClosestCellIds = new long[toMesh->getInternalCount()];
    int *fromClosestCellRanks = new int[toMesh->getInternalCount()];
    double *fromClosestCellDist = new double[toMesh->getInternalCount()];
    fromTree.findPointClosestGlobalCell(toMesh->getInternalCount(),toCellCenters,fromClosestCellIds,fromClosestCellRanks,fromClosestCellDist);

    //Counting points per rank
    std::vector<long> nofPointsPerRank(toMesh->getProcessorCount(),0);
    for(int i = 0; i < toMesh->getInternalCount(); ++i) {
        ++nofPointsPerRank[fromClosestCellRanks[i]];
    }
    //Map collecting (projectionId, originId, originCellCenter) per rank
    std::vector<std::vector<InterpolationInfo>> interpolationInfoPerRank;
    interpolationInfoPerRank.resize(toMesh->getProcessorCount());
    for(int r = 0; r < toMesh->getProcessorCount(); ++r) {
        interpolationInfoPerRank[r].reserve(nofPointsPerRank[r]);
    }
    for(int i = 0; i < toMesh->getInternalCount(); ++i) {
        InterpolationInfo buffer;
        buffer.originId = toCellIds[i];
        buffer.originCellCenter = toCellCenters[i];
        buffer.projectionId = fromClosestCellIds[i];
        interpolationInfoPerRank[fromClosestCellRanks[i]].push_back(buffer);
    }

    //Local interpolation points ready to be interpolated
    const std::vector<InterpolationInfo> & localInterpolationPoints = interpolationInfoPerRank[toMesh->getRank()];

    //Interpolation points from other ranks
    std::unordered_map<int,std::vector<InterpolationInfo>> otherRankIntepolationPoints;

    //Communicate demanding rank and points InterpolationInfo for non-local interpolations
    DataCommunicator interpolationComm(toMesh->getCommunicator());
    // demanding rank (int), nof points, InterpolationInfo (2 long and 3 doubles) per point
    for(int p = 0; p < toMesh->getProcessorCount(); ++p) {
        if(p==toMesh->getRank() || nofPointsPerRank[p] == 0) {
            continue;
        }
        std::size_t buffSize = sizeof(int) + sizeof(long) + (2*sizeof(long) * 3*sizeof(double))*nofPointsPerRank[p];
        interpolationComm.setSend(p,buffSize);
        SendBuffer &sendBuffer = interpolationComm.getSendBuffer(p);
        sendBuffer << toMesh->getRank();
        sendBuffer << nofPointsPerRank[p];
        for(long i = 0; i < nofPointsPerRank[p]; ++i){
            sendBuffer << interpolationInfoPerRank[p][i].originId;
            sendBuffer << interpolationInfoPerRank[p][i].projectionId;
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[0];
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[1];
            sendBuffer << interpolationInfoPerRank[p][i].originCellCenter[2];
        }
    }

    interpolationComm.discoverRecvs();
    interpolationComm.startAllRecvs();
    interpolationComm.startAllSends();

    //Interpolate locals
    darray3 toCellCenter, neighCellCenter;
    long fromCellId;
    double weightSum,cellCentersDist,weight;
    std::vector<long> neighs;
    std::vector<int> globalNeighs;
    std::vector<double> rowElements;
    for(const InterpolationInfo & info : localInterpolationPoints) {
        toCellId = info.originId;
        toCellCenter = info.originCellCenter;
        fromCellId = info.projectionId;
        neighs.clear();
        rowElements.clear();
        globalNeighs.clear();

        fromMesh->findCellNeighs(fromCellId,&neighs);
        globalNeighs.resize(neighs.size() + 1);
        rowElements.resize(neighs.size() + 1);

        //Interpolation
        weightSum = 0.0;
        //fromCell contribution
        long toCellGlobalId = toNumberingInfo->getCellConsecutiveId(toCellId);
        long fromCellGlobalId = fromNumberingInfo->getCellConsecutiveId(fromCellId);
        std::array<double,3> fromCellCenter = fromMesh->evalCellCentroid(fromCellId);
        cellCentersDist = norm2(fromCellCenter-toCellCenter);
        weight = 1.0 / (cellCentersDist*cellCentersDist);
        weightSum += weight;
        globalNeighs[0] = fromCellGlobalId;
        rowElements[0] = weight;

        int counter = 1;
        for(const long & neigh : neighs) {
            neighCellCenter = fromMesh->evalCellCentroid(neigh);
            cellCentersDist = norm2(neighCellCenter-toCellCenter);
            weight = 1.0 / (cellCentersDist*cellCentersDist);
            weightSum += weight;
            globalNeighs[counter] = fromNumberingInfo->getCellConsecutiveId(neigh);
            rowElements[counter] = weight;
            ++counter;
        }

        rowElements /= weightSum;

        MatrixRow matrixRow(toCellGlobalId,globalNeighs,rowElements);
        rowCollection[rowCounter] = matrixRow;
        ++rowCounter;
        //MatSetValues(*interpolationMatrix,1,row,1,globalNeighs.data(),rowElements.data(),INSERT_VALUES);
    }

    std::vector<int> recvRanks = interpolationComm.getRecvRanks();
    for(int rank : recvRanks) {
        interpolationComm.waitRecv(rank);
        RecvBuffer & recvBuffer = interpolationComm.getRecvBuffer(rank);
        int demandingRank;
        long nofPoints;
        InterpolationInfo pointInterpolationInfo;
        recvBuffer >> demandingRank;
        recvBuffer >> nofPoints;
        otherRankIntepolationPoints[demandingRank] = std::vector<InterpolationInfo>(nofPoints,pointInterpolationInfo);
        for(long p = 0; p < nofPoints; ++p) {
            recvBuffer >> pointInterpolationInfo.originId;
            recvBuffer >> pointInterpolationInfo.projectionId;
            recvBuffer >> pointInterpolationInfo.originCellCenter[0];
            recvBuffer >> pointInterpolationInfo.originCellCenter[1];
            recvBuffer >> pointInterpolationInfo.originCellCenter[2];
            otherRankIntepolationPoints[demandingRank][p] = pointInterpolationInfo;
        }
    }

    // Interpolate others and organize results by rank demanding interpolation
    std::unordered_map<int,std::vector<InterpolationMatrixInfo>> otherInterpolationElements;
    InterpolationMatrixInfo emptyInterpolatedInfo;
    for(const auto & rankInfo : otherRankIntepolationPoints) {
        otherInterpolationElements[rankInfo.first] = std::vector<InterpolationMatrixInfo>(rankInfo.second.size(),emptyInterpolatedInfo);
    }

    for(const auto & rankInfo : otherRankIntepolationPoints) {
        long pointCounter = 0;
        for(const InterpolationInfo & info : rankInfo.second) {
            toCellId = info.originId;
            toCellCenter = info.originCellCenter;
            fromCellId = info.projectionId;
            neighs.clear();
            rowElements.clear();
            globalNeighs.clear();

            fromMesh->findCellNeighs(fromCellId,&neighs);
            globalNeighs.resize(neighs.size() + 1);
            rowElements.resize(neighs.size() + 1);

            //Interpolation
            weightSum = 0.0;
            //fromCell contribution
            long toCellGlobalId = toNumberingInfo->getCellConsecutiveId(toCellId);
            long fromCellGlobalId = fromNumberingInfo->getCellConsecutiveId(fromCellId);
            std::array<double,3> fromCellCenter = fromMesh->evalCellCentroid(fromCellId);
            cellCentersDist = norm2(fromCellCenter-toCellCenter);
            weight = 1.0 / (cellCentersDist*cellCentersDist);
            weightSum += weight;
            globalNeighs[0] = fromCellGlobalId;
            rowElements[0] = weight;

            int counter = 1;
            for(const long & neigh : neighs) {
                neighCellCenter = fromMesh->evalCellCentroid(neigh);
                cellCentersDist = norm2(neighCellCenter-toCellCenter);
                weight = 1.0 / (cellCentersDist*cellCentersDist);
                weightSum += weight;
                globalNeighs[counter] = fromNumberingInfo->getCellConsecutiveId(neigh);
                rowElements[counter] = weight;
                ++counter;
            }
            rowElements /= weightSum;

            InterpolationMatrixInfo interpolationMatrixInfo;
            interpolationMatrixInfo.originConsecutiveId = toCellGlobalId;
            interpolationMatrixInfo.nofElements = globalNeighs.size();
            interpolationMatrixInfo.indices = globalNeighs;
            interpolationMatrixInfo.elements = rowElements;
            otherInterpolationElements[rankInfo.first][pointCounter] = interpolationMatrixInfo;
            ++pointCounter;
        }
    }

//    std::stringstream sss;
//    sss << "other_" << toMesh->getRank() << ".txt";
//    std::ofstream out1(sss.str().c_str());
//    for(const auto other : otherInterpolationElements) {
//        for(const auto point : other.second) {
//            out1 << point.originConsecutiveId << " - " << point.nofElements << " - " << point.indices << std::endl;
//        }
//    }
//    out1.close();

    //Communicate Interpolation elements -
    //buffer = (originId(long) + interpolated value(double))*nofInterpolatedValues + nofInterpolatedValues(long)
    DataCommunicator interpolationMatrixComm(toMesh->getCommunicator());
    for(const auto & rankInfo : otherInterpolationElements) {
        //std::size_t buffSize = sizeof(long) + (sizeof(long) + sizeof(double)) * rankInfo.second.size();
        std::size_t buffSize = sizeof(long); // number of points
        for(const auto point : rankInfo.second) {
            buffSize += sizeof(int); // nofElements
            buffSize += (sizeof(int) + sizeof(double)) * point.elements.size(); // indices and elements
            buffSize += sizeof(long); // orginConsecutiveId
        }
        interpolationMatrixComm.setSend(rankInfo.first,buffSize);
        SendBuffer &sendBuffer = interpolationMatrixComm.getSendBuffer(rankInfo.first);
        sendBuffer << rankInfo.second.size();
        for(size_t p = 0; p < rankInfo.second.size(); ++p) {
            sendBuffer << rankInfo.second[p].nofElements;
            sendBuffer << rankInfo.second[p].originConsecutiveId;
            for(size_t e = 0; e < rankInfo.second[p].indices.size(); ++e) {
                sendBuffer << rankInfo.second[p].indices[e];
                sendBuffer << rankInfo.second[p].elements[e];
            }
        }
    }
    interpolationMatrixComm.discoverRecvs();
    interpolationMatrixComm.startAllRecvs();
    interpolationMatrixComm.startAllSends();

    recvRanks.clear();
    recvRanks = interpolationMatrixComm.getRecvRanks();
    long row[1];
    int nofElements = 0;
    std::vector<int> indices;
    std::vector<double> elements;
    for(int rank : recvRanks) {
        interpolationMatrixComm.waitRecv(rank);
        RecvBuffer & recvBuffer = interpolationMatrixComm.getRecvBuffer(rank);
        long nofPoints;
        recvBuffer >> nofPoints;
        for(long p = 0; p < nofPoints; ++p) {
            recvBuffer >> nofElements;
            indices.resize(nofElements);
            elements.resize(nofElements);
            recvBuffer >> row[0];
            for(int e = 0; e < nofElements; ++e) {
                recvBuffer >> indices[e];
                recvBuffer >> elements[e];
            }

            MatrixRow matrixRow(int(row[0]),indices,elements);
            rowCollection[rowCounter] = matrixRow;
            ++rowCounter;
            //MatSetValues(*interpolationMatrix, 1, row, nofElements, indices.data(), elements.data(),INSERT_VALUES);
        }
    }

    //Pre-Allocate
    std::vector<int> d_nnz(toMesh->getInternalCount()), o_nnz(toMesh->getInternalCount());
    for( const auto r : rowCollection) {
        int rowLocalConsecutiveId = r.row - toNumberingInfo->getCellGlobalCountOffset();
        for(const auto i : r.indices) {
            if(fromNumberingInfo->getCellRankFromConsecutive(i) == fromMesh->getRank()) {
                ++d_nnz[rowLocalConsecutiveId];
            } else {
                ++o_nnz[rowLocalConsecutiveId];
            }
        }
    }
    MatCreateAIJ(toMesh->getCommunicator(), toMesh->getInternalCount(),fromMesh->getInternalCount(), PETSC_DETERMINE, PETSC_DETERMINE,
            0,d_nnz.data(),0,o_nnz.data(),interpolationMatrix);

    //Set elements
//    std::stringstream ss;
//    ss << "rows_" << toMesh->getRank() << ".txt";
//    std::ofstream out(ss.str().c_str());
    for(const auto r : rowCollection) {
        int rowArray[1];
        rowArray[0] = r.row;
//        out << rowArray[0] << " - " << r.indices << std::endl;
        MatSetValues(*interpolationMatrix,1,rowArray,r.indices.size(),r.indices.data(),r.elements.data(),INSERT_VALUES);
    }

    MatAssemblyBegin(*interpolationMatrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*interpolationMatrix, MAT_FINAL_ASSEMBLY);

    //DEBUG
//    PetscViewer matViewer;
//    PetscViewerCreate(toMesh->getCommunicator(), &matViewer);
//    PetscViewerSetType(matViewer, PETSCVIEWERASCII);
//    PetscViewerFileSetMode(matViewer, FILE_MODE_WRITE);
//    PetscViewerPushFormat(matViewer, PETSC_VIEWER_ASCII_MATLAB);
//
//    std::stringstream filePathStream;
//    filePathStream.str(std::string());
//    filePathStream << "./interpolationMatrix.txt";
//    PetscViewerFileSetName(matViewer, filePathStream.str().c_str());
//    MatView(*interpolationMatrix, matViewer);
//    PetscViewerDestroy(&matViewer);


//    out.close();

    delete [] toCellCenters;
    delete [] toCellIds;
    delete [] fromClosestCellIds;
    delete [] fromClosestCellRanks;
    delete [] fromClosestCellDist;
};



/*!
    Initialize a mesh-coherent PiercedVector Container with sin^2(sqrt(x_0^2+x_1^2+x_2^2))

    \param[in] mesh the mesh
    \return the initialized PiercedVector
*/
void initDoubleDataOnMesh(SurfUnstructured * mesh, PiercedVector<double>* data){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(const Vertex & v: vertices) {
        darray3 x = v.getCoords();
        double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        double datum = acos(x[2]/r) - M_PI/2;//sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        datum = sin(4*datum);
        datum *= datum;
        data->insert(v.getId(),datum);
    }
};

/*!
    Initialize a mesh-coherent PiercedVector Container with a C-array.
    It has to be noted that the array should contains data with the same order of the vertices of the associated mesh.

    \param[in] mesh the mesh
    \param[out] data PiercedVector container associated to the mesh that has to be filled
    \param[in] array a pointer to a C-array containing data to be inserted into the PiercedVector
    \param[in] arraySize the size of the C-array
*/
void initDataOnMeshFromArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    assert(vertices.size()==arraySize);
    size_t count = 0;
    std::cout << "Start inserting ...";
    for(const Vertex & v: vertices) {
        data->insert(v.getId(),array[count]);
        ++count;
    }
};

/*!
    Move data from a mesh-coherent PiercedVector Container to a C-array.
    It has to be noted that the array will contain data with the same order of the vertices of the associated mesh.

    \param[in] mesh the mesh
    \param[in] data PiercedVector container associated to the mesh containing data to be inserted into the C-array
    \param[out] array a pointer to a C-array to be filled with data coming from the PiercedVector
    \param[in] arraySize the size of the C-array
*/
void moveDataOnMeshToArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    assert(vertices.size()==arraySize);
    size_t count = 0;
    for(const Vertex & v: vertices) {
        array[count] = data->at(v.getId());
        ++count;
    }
};



/*!
    Write VTK file( (p)vtu ) containing the mesh

    \param[in] mesh the mesh
    \param[in] filename the name of the .(p)vtu file(s)
*/
void writeMesh(SurfUnstructured * mesh,std::string filename){
    mesh->getVTK().setName(filename);
    mesh->write();

};

/*!
    Write VTK file( (p)vtu ) containing the mesh and the passed data

    \param[in] mesh the mesh
    \param[in] filename the name of the .(p)vtu file(s)
    \param[in] data the PiercedVector containing the data to be plotted. NB only one scalar field is allowed
    \param[in] dataNames a list of the field data names. NB its size has to be 1.
*/
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> * data,const std::vector<std::string> & dataNames){

    mesh->getVTK().setName(filename);
    std::vector<double> vdata;
    vdata.reserve(data->size());
    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(const Vertex & v : vertices){
        //vdata.push_back(*(data.find(v.getId())));
        vdata.push_back(data->at(v.getId()));
    }
    mesh->getVTK().addData(dataNames[0], VTKFieldType::SCALAR, VTKLocation::POINT, vdata);
    mesh->write();

};

/*!
    Compute an ordered mesh partitioning, filling an empty long array (provided by the user) with the cell indices to be assigned to each rank.
    The cell numbering is given by the ordering in the mesh file. CAVEAT: The maximum global number of cells is bounded to 2 billions, due to the use of INT.
    bitpit can manage more cells by using special MPI datatype.

    \param[in] meshFile the file containing the mesh
    \param[in] comm the communicator used to partition the mesh
    \param[out] cellSizesPerRank a pointer to an empty long array provided by the user and filled with the cell indices for each rank.
*/
std::vector<int> computeMeshFilePartitioning(const std::string meshFile,MPI_Comm comm){

    std::vector<int> idRanks;
    //Ask communicator for its size
    int nofRanks;
    int rank;
    MPI_Comm_size(comm,&nofRanks);
    MPI_Comm_rank(comm,&rank);

    std::vector<int> sizes(nofRanks,0);
    int nofCells = 0;
    if(rank == 0){
        //Rank 0 reads mesh file and counts the cells
        SurfUnstructured mesh(2,3);
        mesh.importSTL(meshFile);
        mesh.deleteCoincidentVertices();
        nofCells = mesh.getCellCount();

        //Rank 0 compute the number of cells per rank
        int integerDivision = nofCells / nofRanks;
        int divisionReminder = nofCells % nofRanks;
        for(int r = 0; r < nofRanks; ++r) {
            sizes[r] = integerDivision;
        }
        for(int i = 0; i < divisionReminder; ++i) {
            ++sizes[i];
        }

        idRanks.resize(nofCells,0);
        int count = 0;
        for(int r = 0; r < nofRanks; ++r) {
            for(int i = 0; i < sizes[r]; ++i){
                idRanks[count] = r;
                ++count;
            }
        }
    }

    //Rank 0 broadcast nofCells
    MPI_Bcast(&nofCells,1,MPI_INT,0,comm);
    idRanks.resize(nofCells,0);
    MPI_Bcast(idRanks.data(),nofCells,MPI_INT,0,comm);

//    //DEBUG
//    for(int r = 0; r < nofRanks; ++r) {
//        if(rank == r) {
//            std::cout << "I'm rank " << rank << " and my idRanks is : " << std::endl;
//            for(int i = 0; i < nofCells; ++i) {
//                std::cout << "id " << i << " -> rank " << idRanks[i] << std::endl;
//            }
//        }
//        MPI_Barrier(comm);
//    }
//    //DEBUG

    return idRanks;
};


/*!
    Constructor

    \param[in] patch the patch (mesh) to be written in VTK format with attached data
    \param[in] scalarField
*/
FieldStreamer::FieldStreamer(const PatchKernel &patch, const PiercedStorage<double, long> &scalarField,std::vector<std::string> &fieldNames)
: m_patch(patch), m_scalarField(scalarField),m_fieldNames(fieldNames)
{
};


/*!
    Virtual method specialization to define how the data have to be written

    \param[in] stream the file stream used to write the VTK files
    \param[in] name the name of the field to be written
    \param[in] format the VTK format to write data
*/
void FieldStreamer::flushData(std::fstream &stream, const std::string & name, VTKFormat format)
{

    assert(format == VTKFormat::APPENDED);
    BITPIT_UNUSED(format);

    int fieldIndex;
    for(const std::string & fieldName : m_fieldNames) {
        if(name == fieldName) {
            if( name.find("temperature") != std::string::npos || name.find("Temperature") != std::string::npos ) {
                fieldIndex = 0;
            } else {
                fieldIndex = 1;
            }
            for (const Cell &cell : m_patch.getVTKCellWriteRange()) {
                long id = cell.getId();
                genericIO::flushBINARY(stream, m_scalarField.at(id,fieldIndex));
            }
        }
    }

};

const std::vector<std::string> & FieldStreamer::getFieldNames() const{

    return m_fieldNames;

}

}


