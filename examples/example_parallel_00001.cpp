#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "bitpit_surfunstructured.hpp"
#include "commons.hpp"
#include "coupling.hpp"
#include "couplingUtils.hpp"
#include <exception>
#if ENABLE_MPI==1
#include "mpi.h"
#include <petscksp.h>
#endif

using namespace bitpit;

// =================================================================================== //

void test00001( int argc, char *argv[] ) {

    //
    // Initialization
    //

    // Initialize input files (default/passed by arguments)
    std::string disciplineFilemesh = "./data/unitSphere5.stl"; //NEUTRAL MESH SHOULD BE THE SAME MESH OF THE INNER DISCIPLINE TO AVOID RAY-TRACING FOR FLUX EXCHANGE
    std::string neutralFilemesh = "./data/unitSphere4.stl";
    if (argc > 2) {
        neutralFilemesh = argv[1];
        disciplineFilemesh = argv[2];
    }

    // Initialize parallel (each process runs independently from the others)
    int nProcessors;
    int rank, workRank;

#if ENABLE_MPI==1
    //MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcessors = 1;
    rank = 0;
#endif

    // Initialize logger
    log::manager().initialize(log::COMBINED, "madeleine", true, ".", nProcessors, rank);
#if ENABLE_DEBUG==1
    log::cout().setVisibility(log::GLOBAL);
#endif

    // Log file header
    log::cout() << "|=====================================================|" << std::endl;
    log::cout() << "|                                                     |" << std::endl;
    log::cout() << "|                 MADELEINE interpolation             |" << std::endl;
    log::cout() << "|                                                     |" << std::endl;
    log::cout() << "|=====================================================|" << std::endl;
    log::cout() << "| version:  v1.0                                      |" << std::endl;
    log::cout() << "|                                                     |" << std::endl;
    log::cout() << "| All rights reserved.                                |" << std::endl;
    log::cout() << "|=====================================================|" << std::endl;
    log::cout() << std::endl;

//    std::vector<std::string> inputs(1,"Flux"); //for inner sphere (outputs for outer one)
//    std::vector<std::string> outputs(1,"Temperature"); // for inner sphere (inputs for outer one)
    std::vector<std::string> inputs(1,"TemperatureIN"); //for outer sphere (outputs for outer one)
    std::vector<std::string> outputs(1,"TemperatureOUT"); // for outer sphere (inputs for outer one)

    //Create the world group
    MPI_Group worldGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

    //Select worker ranks in reduced communicator
    //int nofRankInWorkGroup = nProcessors != 1 ? nProcessors - 1 : 1;
    // PetscInitiliaze is called by bitpit mostly to pass custom arguments to petsc. PetscInitialize must be called by all the processes
    // If called from outside, bitpit cannot customize the linear solver using PETSc command-line arguments.
    // The default solver in bitpit is GMRES with ASM global preconditioner and sub-block ILU preconditioner
    int nofRankInWorkGroup = nProcessors;
    int *ranks = new int[nofRankInWorkGroup];
    for(int i = 0; i < nofRankInWorkGroup; ++i) {
        ranks[i] = i;
    }

    //DEBUG
    for(int i = 0; i < nofRankInWorkGroup; ++i) {
        std::cout << ranks[i] << std::endl;
    }
    std::cout << std::endl;
    //DEBUG

    //Create workers group
    MPI_Group workersGroup;
    MPI_Group_incl(worldGroup,nofRankInWorkGroup,ranks,&workersGroup);

    //Create workers communicator
    MPI_Comm workersComm;
    MPI_Comm_create_group(MPI_COMM_WORLD,workersGroup,0,&workersComm);

    //Make ONLY workers work
    if(workersComm != MPI_COMM_NULL) {
        MPI_Comm_rank(workersComm, &workRank);
        std::vector<int> cellIndicesPerRank = coupling::computeMeshFilePartitioning(neutralFilemesh,workersComm);
        //
        coupling::MeshCoupling parallelToyDiscipline1(inputs,outputs,"InnerSphere",workersComm);

        double radius = 2.0;
        double thickness = 1.0;
        bool isInnerSphere = false;
        double sourceIntensity = 1.0;
        //std::array<double,3> sourceDirection = {{1.0,0.0,0.0}};
        std::vector<double> sourceDirection(3,0.0);
        sourceDirection[0] = 1.0;
        double thermalDiffusivityCoefficient = 0.001;
        double emissivity = 0.001;
        double infinityTemperature = 274.0;

        parallelToyDiscipline1.initialize(disciplineFilemesh,neutralFilemesh,
                radius,radius,thickness,isInnerSphere,sourceIntensity,sourceDirection,
                thermalDiffusivityCoefficient, emissivity, infinityTemperature,
                cellIndicesPerRank);
        //
        //    PiercedVector<double> neutralData;
        //    coupling::initDoubleDataOnMesh(parallelToyDiscipline.getNeutralMesh(),&neutralData);

        //Homogeneous initialization of inputs from neutral mesh
        double initialTemperature = 274.0;
        double initialFlux = 0.274;

        std::vector<double> neutralData(parallelToyDiscipline1.getNeutralMesh()->getInternalCount(),initialTemperature); //Flux for inner, Temperature for outer DEBUG=rank
        parallelToyDiscipline1.compute(neutralData.data(),neutralData.size(),radius,radius);

        int globalRowId;
        std::vector<int> cols;
        std::vector<double> values;
        for(int i = 0; i < parallelToyDiscipline1.getNeutralMesh()->getInternalCount(); ++i) {
            //parallelToyDiscipline1.extractOutputInputJacobianRow(i, globalRowId, cols, values);
            //std::cout << "global row ID = " << globalRowId << " cols = " << cols << " values = " << values << std::endl;
            parallelToyDiscipline1.extractOutputControlJacobianRow(i, globalRowId, cols, values);
            std::cout << "global row ID = " << globalRowId << " cols = " << cols << " values = " << values << std::endl;
        }

        //
        //    parallelToyDiscipline.close();
    }

    //Free groups and workers communicator
    MPI_Group_free(&worldGroup);
    MPI_Group_free(&workersGroup);
    if(workersComm != MPI_COMM_NULL){
        MPI_Comm_free(&workersComm);
    }
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI_Init(&argc,&argv);
    const char help[] = "None";
    PetscInitialize(&argc, &argv, 0, help);

    {
#endif
        test00001(argc,argv) ;
#if ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;

}
