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
void initDoubleDataOnMesh(SurfUnstructured * mesh, PiercedVector<double>* data);
void initDataOnMeshFromArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void moveDataOnMeshToArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void writeMesh(SurfUnstructured * mesh,std::string filename);
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> * data,const std::vector<std::string> & dataNames);
#if ENABLE_MPI==1
void computeMeshFilePartitioning(const std::string meshFile,std::vector<int> & idRank,MPI_Comm comm = MPI_COMM_WORLD);
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


}
#endif /* SRC_COUPLINGUTILS_HPP_ */

