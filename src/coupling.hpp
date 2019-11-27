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
    void staticPartitionDisciplineMeshByNeutralFile(); //partitioning D_Nf

    void computeGlobalDisciplineId2NeutralRank(); //compute global vector to get D_Nf
    void computeDiscipline2FileNeutralCellPerRanks(); // compute local vector to get D_Nf

    void buildScaledMeshes();

    void computeGlobalNeutralId2DisciplineRank(SurfUnstructured *serialNeutralMesh); //compute global vector to get N_{D_Nf}


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
    std::vector<int> m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank;             //global vector to build neutral cellPerRank to get partitioning overlapping the discipline neutral mesh file partitioning (N_{D_Nf})
    std::vector<int> m_globalDisciplineId2NeutralMeshFileRank;                               //global vector to build discipline cellPerRank to get partitioning overlapping neutral mesh file partitioning (D_Nf)
    std::vector<int> m_globalNeutralId2MeshFileRank;                                         //computed outside this class. Useful to build cellPerRanks to get neutral mesh file partitioning (Nf)
    std::unordered_map<long,int> m_discipline2FileNeutralCellPerRanks;                       //local map to partition the discipline mesh overlapping the neutral mesh file partitioning
    std::unordered_map<long,int> m_disciplineId2NeutralMeshFileCellPerRanks;                    //(D_Nf)     local map to partition the discipline mesh overlapping the neutral mesh file partitioning. It depends on the starting discipline mesh partitioning

    //oldies
    std::unordered_map<long,int> m_neutralFile2DisciplineCellPerRanks;

#endif

};

}
#endif /* SRC_COUPLING_HPP_ */
