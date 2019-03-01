/*
 * coupling.hpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#ifndef _MADELEINE_COUPLING_HPP_
#define _MADELEINE_COUPLING_HPP_

// ========================================================================== //
// INCLUDE                                                                    //
// ========================================================================== //

#if ENABLE_MPI==1
#    include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

namespace coupling{

class MeshCoupling{

public:
    MeshCoupling();
    MeshCoupling(const vector<std::string> & inputNames, std::vector<std::string> & outputNames);
    void initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile, double sphereRadius);
    void compute(PiercedVector<double,long> * neutralData);
    const std::vector<std::string> & getInputDataNames();
    const std::vector<std::string> & getOutputDataNames();
    const SurfUnstructured * getDisciplineMesh();
    SurfUnstructured * getNeutralMesh();
    size_t getNeutralMeshSize();
    void close();

private:
    std::unique_ptr<SurfUnstructured> m_unitDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_unitNeutralMesh;
    double m_radius;
    std::unique_ptr<SurfUnstructured> m_scaledDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_scaledNeutralMesh;

    PiercedVector<double,long> m_disciplineData;

    std::vector<std::string> m_inputDataNames;
    std::vector<std::string> m_outputDataNames;
};

}
#endif /* SRC_COUPLING_HPP_ */
