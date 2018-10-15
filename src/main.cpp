#if ENABLE_MPI==1
#    include <mpi.h>
#endif

#include <bitpit_IO.hpp>

using namespace bitpit;

/*
 * Main
 */
int main(int argc, char *argv[])
{
    //
    // Initialization
    //

    // Initialize parallel
    int nProcessors;
    int rank;

#if ENABLE_MPI==1
    MPI_Init(&argc, &argv);

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

    //
    // Run
    //


    //
    // Finalization
    //

#if ENABLE_MPI==1
    // MPI finalization
    MPI_Finalize();
#endif

}
