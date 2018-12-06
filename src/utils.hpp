/*
 * utils.hpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

//SurfUnstructured scale(const SurfUnstructured & unitRadiusSphereMesh, double radius);
void interpolateFromTo(SurfUnstructured * fromMesh, const PiercedVector<double> & fromData, SurfUnstructured * toMesh, PiercedVector<double> & toData);
void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector<double>* data);
void writeMesh(SurfUnstructured * mesh,std::string filename);
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> & data,const std::vector<std::string> & dataNames);
#endif /* SRC_UTILS_HPP_ */
