#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"
#include "coupling.hpp"
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
    std::string disciplineFilemesh = "./data/unitSphere1.stl";
    std::string neutralFilemesh = "./data/unitSphere2.stl";
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
    coupling::MeshCoupling neutralToDiscipline(inputs,outputs);

    neutralToDiscipline.initialize(disciplineFilemesh,neutralFilemesh,2.0);

    PiercedVector<double> neutralData;
    coupling::initDoubleDataOnMesh(neutralToDiscipline.getNeutralMesh(),&neutralData);
    neutralToDiscipline.compute(neutralData);

    neutralToDiscipline.close();

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI::Init(argc, argv);

    {
#endif
        //try{
            test00001(argc,argv) ;
//        }
//        catch(std::exception & e){
//            std::cout<<"example_00002 exited with an error of type : "<<e.what()<<std::endl;
//            return 1;
//        }
#if ENABLE_MPI==1
    }

    MPI_Finalize();
#endif

    return 0;

}
