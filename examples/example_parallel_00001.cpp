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
#endif

using namespace bitpit;

// =================================================================================== //

void test00001( int argc, char *argv[] ) {

    //
    // Initialization
    //

    // Initialize input files (default/passed by arguments)
    std::string disciplineFilemesh = "./data/unitSphere5.stl";
    std::string neutralFilemesh = "./data/unitSphere4.stl";
    if (argc > 2) {
        neutralFilemesh = argv[1];
        disciplineFilemesh = argv[2];
    }

    // Initialize parallel (each process runs independently from the others)
    int nProcessors;
    int rank;

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

    std::vector<std::string> inputs(1,"Forces");

    std::vector<std::string> outputs(1,"Pressure");

    //Create the world group
    MPI_Group worldGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

    //Select worker ranks in reduced communicator
    int nofRankInWorkGroup = nProcessors != 1 ? nProcessors - 1 : 1;
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
        std::vector<int> cellIndicesPerRank;
        coupling::computeMeshFilePartitioning(neutralFilemesh,cellIndicesPerRank,workersComm);
        //
        coupling::MeshCoupling parallelToyDiscipline1(inputs,outputs,"toy1",workersComm);

        parallelToyDiscipline1.initialize(disciplineFilemesh,neutralFilemesh,2.0,cellIndicesPerRank);
        //
        //    PiercedVector<double> neutralData;
        //    coupling::initDoubleDataOnMesh(parallelToyDiscipline.getNeutralMesh(),&neutralData);
        //    parallelToyDiscipline.compute(&neutralData);
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

    {
#endif
        test00001(argc,argv) ;
#if ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;

}
