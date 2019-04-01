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
#include "coupling.hpp"

using namespace bitpit;

namespace coupling {
//SurfUnstructured scale(const SurfUnstructured & unitRadiusSphereMesh, double radius);
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedVector<double> * fromData, SurfUnstructured * toMesh, PiercedVector<double> * toData);
void initDoubleDataOnMesh(SurfUnstructured * mesh, PiercedVector<double>* data);
void initDataOnMeshFromArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void moveDataOnMeshToArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize);
void writeMesh(SurfUnstructured * mesh,std::string filename);
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> * data,const std::vector<std::string> & dataNames);
}
#endif /* SRC_COUPLINGUTILS_HPP_ */

