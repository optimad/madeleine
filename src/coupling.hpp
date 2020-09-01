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
#include <bitpit_discretization.hpp>
#include <bitpit_LA.hpp>
#include <bitpit_surfunstructured.hpp>
#include "commons.hpp"
#include "couplingUtils.hpp"
#include "communications.hpp"
#include "matrixManager.hpp"

using namespace bitpit;

namespace coupling{

class MeshCoupling{

public:
#if ENABLE_MPI==1
    MeshCoupling(std::string disciplineName = "defaultDiscipline", MPI_Comm = MPI_COMM_WORLD);
    MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName,
            MPI_Comm = MPI_COMM_WORLD);
#else
    MeshCoupling(std::string disciplineName = "defaultDiscipline");
    MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName);
#endif
    void initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile,
            double sphereRadius, double sphereNeutralRadius,
            double sphereThickness, bool innerSphere, double sourceIntensity, std::vector<double> sourceDirection,
            double thermalDiffusivityCoefficient, double emissivity,
            const std::vector<int> & globalNeutralId2MeshFileRank, int kernel);
    void initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile,
            double sphereDisciplineRadius, double sphereNeutralRadius,
            double sphereThickness, bool innerSphere, double sourceIntensity, std::vector<double> sourceDirection,
            double thermalDiffusivityCoefficient, double emissivity,
            const std::vector<int> & globalNeutralId2MeshFileRank);
    void compute(double *neutralInputArray, std::size_t size, double newRadius, double otherRadius);
    const std::vector<std::string> & getInputDataNames();
    const std::vector<std::string> & getOutputDataNames();
    const SurfUnstructured * getDisciplineMesh();
    SurfUnstructured * getNeutralMesh();
    size_t getNeutralMeshSize();
    long getNeutralFirstCellId();
    long getNeutralGlobalConsecutiveId(long id);
    void close();
    void extractOutputInputJacobianRow(int index, int & cellGlobalId, std::vector<int> & columnIds, std::vector<double> & columnValues);
    void extractOutputControlJacobianRow(int index, int & cellGlobalId, std::vector<int> & columnIds, std::vector<double> & columnValues);
    JacobianMatricesManager* getJacobianManager();

private:
    void readUnitDisciplineMesh();
    void readUnitNeutralMesh();
    void staticPartitionNeutralMeshByMeshFileOrder(); //partitioning Nf
    void staticPartitionDisciplineMeshByNeutralFile(); //partitioning D_Nf

    void computeGlobalDisciplineId2NeutralRank(); //compute global vector to get D_Nf. It needs to be executed just once. The result is in m_globalDisciplineId2NeutralMeshFileRank.
    void computeDiscipline2FileNeutralCellPerRanks(); // compute local map to get D_Nf

    void buildScaledMeshes();

    void computeGlobalNeutralId2DisciplineRank(SurfUnstructured *serialNeutralMesh); //compute global vector to get N_{D_Nf}. It needs to be executed just once. The result is in m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank
    void computeNeutralId2DisciplineCellPerRanks(); //compute local map to get N_{D_Nf}

    void dynamicPartitionNeutralMeshByNeutralMeshFilePartitionedDiscipline();
    void dynamicPartitionNeutralMeshByNeutralMeshWithData();

    void synchronizeMeshData();
    void uniformlyInitAllData(double value);
    void prepareWritingData();

    void createGhostCommunicators(bool continuous);
    void initializeGhostCommunicators();
    void updateNeutralGhosts();
    void updateDisciplineGhosts();

    void disciplineKernel();
    void disciplineKernel1();
    void disciplineKernel2();
    const PatchNumberingInfo & getNeutraNumberingInfo();

    std::unique_ptr<SurfUnstructured> m_unitDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_unitNeutralMesh;
    double m_disciplineRadius,m_neutralRadius;
    double m_oldDisciplineRadius,m_oldNeutralRadius;
    std::unique_ptr<SurfUnstructured> m_scaledDisciplineMesh;
    std::unique_ptr<SurfUnstructured> m_scaledNeutralMesh;
    PatchNumberingInfo m_neutralNumberingInfo;
    PatchNumberingInfo m_disciplineNumberingInfo;

    PiercedStorage<double,long> m_disciplineData;
    PiercedStorage<double,long> m_neutralData;

    std::vector<std::string> m_inputDataNames;
    std::vector<std::string> m_outputDataNames;

    std::string m_name;

    int m_kernel;
    std::unique_ptr<StencilScalarSolverHandler> m_system;
    std::vector<StencilScalar> m_helmoltzStencils;
    std::unique_ptr<JacobianMatricesManager> m_jacobianManager;

    double m_thickness;
    bool m_innerSphere;
    double m_sourceMaxIntensity;
    std::array<double,3> m_sourceDirection;
    double m_thermalDiffusivityCoefficient; //linear radius dependency
    double m_emissivity;

    void scaleMeshToRadius(std::unique_ptr<SurfUnstructured> & mesh, double & oldRadius, const double & newRadius);

    double evalThermalDiffusivity();
    void assemblySimplifiedDiscreteHelmholtzSystem();
    void updateSimplifiedDiscreteHelmholtzSystem();
    void updateSystemRHS();
    void computeSimplifiedDiscreteLaplaceStencils(std::vector<StencilScalar> & laplaceStencils);
    void computeHelmholtzStencilsFromLaplaceStencils(std::vector<StencilScalar> & helmoltzStencils, const double & coefficient);
    double evalSourceIntensity(const std::array<double,3> & cellNormal);
    void solveSytem();
    void updateOutputField();
    void updateInputField(double *neutralInputArray, std::size_t size);

    std::string m_disciplineMeshFileName;
    std::string m_neutralMeshFileName;

    std::unique_ptr<coupling::FieldStreamer> m_disciplineVTKFieldStreamer;
    std::unique_ptr<coupling::FieldStreamer> m_neutralVTKFieldStreamer;

    std::unique_ptr<GhostCommunicator> m_neutralGhostCommunicator;
    std::unique_ptr<GhostCommunicator> m_disciplineGhostCommunicator;
    int m_neutralTag;
    int m_disciplineTag;
    std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>>> m_neutralGhostStreamer;
    std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>>> m_disciplineGhostStreamer;

    //oldies
    std::unordered_map<long,int> computeStaticPartitionByMetis(SurfUnstructured & mesh);
    void staticPartitionMeshByMetis(SurfUnstructured & mesh);
    void staticPartitionDisciplineMeshByMetis();
    void dynamicPartitionNeutralMeshByDiscipline();

#if ENABLE_MPI==1
    MPI_Comm m_comm;
    int m_nprocs;
    int m_rank;
    std::vector<int> m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank;                //(N_{D_Nf}) global vector to build neutral cellPerRank to get partitioning overlapping the discipline neutral mesh file partitioning
    std::vector<int> m_globalDisciplineId2NeutralMeshFileRank;                                  //(D_Nf)     global vector to build discipline cellPerRank to get partitioning overlapping neutral mesh file partitioning
    std::vector<int> m_globalNeutralId2MeshFileRank;                                            //(Nf)       computed outside this class. Useful to build cellPerRanks to get neutral mesh file partitioning
    std::unordered_map<long,int> m_disciplineId2NeutralMeshFileCellPerRanks;                    //(D_Nf)     local map to partition the discipline mesh overlapping the neutral mesh file partitioning. It depends on the starting discipline mesh partitioning
    std::unordered_map<long,int> m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks;  //(N_{D_Nf}) local map to partition the neutral mesh overlapping the discipline mesh file partitioning. It depends on the starting neutral mesh partitioning
    std::unordered_map<long,int> m_neutralMeshFilePartitionCellPerRank;                          //(Nf)       local map to partition the neutral mesh with mesh file partitioning. It depends on the starting neutral mesh partitioning

    std::unique_ptr<DataCommunicator> m_lbCommunicator;

    //oldies
    std::unordered_map<long,int> m_neutralFile2DisciplineCellPerRanks;
#endif

    //field access indices
    const int fid_temperature;
    const int fid_flux;
    int m_inputField;
    int m_outputField;
};

}
#endif /* SRC_COUPLING_HPP_ */
