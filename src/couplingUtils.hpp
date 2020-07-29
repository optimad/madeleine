/*
 * couplingUtils.hpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#ifndef SRC_COUPLINGUTILS_HPP_
#define SRC_COUPLINGUTILS_HPP_

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

namespace coupling {
//SurfUnstructured scale(const SurfUnstructured & unitRadiusSphereMesh, double radius);
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedVector<double> * fromData, SurfUnstructured * toMesh, PiercedVector<double> * toData);
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedStorage<double,long> * fromData, SurfUnstructured * toMesh, PiercedStorage<double,long> * toData, int fieldIndex = -1);
void initDoubleDataOnMesh(SurfUnstructured * mesh, PiercedVector<double>* data);
void initDataOnMeshFromArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void moveDataOnMeshToArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void writeMesh(SurfUnstructured * mesh,std::string filename);
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> * data,const std::vector<std::string> & dataNames);
#if ENABLE_MPI==1
std::vector<int> computeMeshFilePartitioning(const std::string meshFile,MPI_Comm comm = MPI_COMM_WORLD);
#endif

// Auxiliary class to export fields in VTK format.
class FieldStreamer : public VTKBaseStreamer {

public:

    FieldStreamer(const PatchKernel &patch, const PiercedStorage<double, long> &scalarField,std::vector<std::string> &fieldNames);
    void flushData(std::fstream &stream, const std::string & name, VTKFormat format);
    const std::vector<std::string> & getFieldNames() const;

private:

    const PatchKernel &m_patch;
    const PiercedStorage<double, long> &m_scalarField;
    const std::vector<std::string> m_fieldNames;

};

struct InterpolationInfo {
    long projectionId;
    long originId;
    std::array<double,3> originCellCenter;

    InterpolationInfo(){projectionId = 0; originId = 0; originCellCenter = {{0.0,0.0,0.0}};};
    void display(std::ofstream & out) const {
        out << "Origin id - Origin center - projected id" << std::endl;
        out << originId << " - " << originCellCenter << " - " << projectionId << std::endl;
    }
    InterpolationInfo & operator=(InterpolationInfo & rhs) {
        this->originCellCenter = rhs.originCellCenter;
        this->projectionId = rhs.projectionId;
        this->originId = rhs.originId;
        return *this;
    }
};

struct InterpolatedInfo {
    long originId;
    double value;
    InterpolatedInfo(){originId = 0; value = 0.0;};
    InterpolatedInfo & operator=(InterpolatedInfo & rhs) {
        this->originId = rhs.originId;
        this->value = rhs.value;
        return *this;
    }
};

}
#endif /* SRC_COUPLINGUTILS_HPP_ */

