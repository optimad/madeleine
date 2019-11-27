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

#include <bitpit_common.hpp>
#include <bitpit_IO.hpp>
#include <bitpit_CG.hpp>
#include <bitpit_surfunstructured.hpp>
#include "commons.hpp"

using namespace bitpit;

namespace coupling{

class MeshCoupling{

public:
#if ENABLE_MPI==1
    MeshCoupling(std::string disciplineName = "defaultDiscipline", MPI_Comm = MPI_COMM_WORLD);
    MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName, MPI_Comm = MPI_COMM_WORLD);
#else
    MeshCoupling(std::string disciplineName = "defaultDiscipline");
    MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName);
#endif
    void initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile, double sphereRadius, const std::vector<int> & globalNeutralId2MeshFileRank);
    void compute(PiercedVector<double,long> * neutralData);
    const std::vector<std::string> & getInputDataNames();
    const std::vector<std::string> & getOutputDataNames();
    const SurfUnstructured * getDisciplineMesh();
    SurfUnstructured * getNeutralMesh();
    size_t getNeutralMeshSize();
    void close();

private:
    void readUnitDisciplineMesh();
    void readUnitNeutralMesh();
    void staticPartitionNeutralMeshByMeshFileOrder(); //partitioning Nf

    void computeGlobalNeutralId2DisciplineRank();
    void computeGlobalDisciplineId2NeutralRank();

    std::unique_ptr<SurfUnstructured> m_unitDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_unitNeutralMesh;
    double m_radius;
    std::unique_ptr<SurfUnstructured> m_scaledDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_scaledNeutralMesh;

    PiercedVector<double,long> m_disciplineData;

    std::vector<std::string> m_inputDataNames;
    std::vector<std::string> m_outputDataNames;

    std::string m_name;

    std::string m_disciplineMeshFileName;
    std::string m_neutralMeshFileName;

    //oldies
    std::unordered_map<long,int> computeStaticPartitionByMetis(SurfUnstructured & mesh);
    void staticPartitionMeshByMetis(SurfUnstructured & mesh);
    void staticPartitionDisciplineMeshByMetis();
    void dynamicPartitionNeutralMeshByDiscipline();

#if ENABLE_MPI==1
    MPI_Comm m_comm;
    int m_nprocs;
    int m_rank;
    std::vector<int> m_globalNeutralId2DisciplineRank;
    std::vector<int> m_globalDisciplineId2NeutralRank;
    std::vector<int> m_globalNeutralId2MeshFileRank;
    std::unordered_map<long,int> m_neutralFile2DisciplineCellPerRanks;
    std::unordered_map<long,int> m_neutralDiscipline2FileCellPerRanks;
    std::unordered_map<long,int> m_discipline2FileNeutralCellPerRanks;
#endif

};

}
#endif /* SRC_COUPLING_HPP_ */
