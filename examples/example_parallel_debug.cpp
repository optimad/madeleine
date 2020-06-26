#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "bitpit_surfunstructured.hpp"
#include "commons.hpp"
#include <exception>
#if ENABLE_MPI==1
#include "mpi.h"
#endif

using namespace bitpit;

// =================================================================================== //

void test00002( int argc, char *argv[] ) {

    int nProcessors;
    int rank;

#if ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcessors = 1;
    rank = 0;
#endif

    SurfUnstructured mesh(2,3);
    mesh.setCommunicator(MPI_COMM_WORLD);

    std::stringstream meshDump;
    std::string folder = "./";
    meshDump << folder << "/neutralMeshFilePartitioned_" << rank << ".dat";
    std::ifstream meshStream;
    meshStream.open(meshDump.str().c_str());
    mesh.restore(meshStream);
    std::string meshIn("meshbefore");
    mesh.write(meshIn);

    std::stringstream mapFile;
    mapFile << folder << "/map_" << rank << ".dat";
    std::unordered_map<long,int> map;
    std::ifstream mapStream;
    mapStream.open(mapFile.str().c_str());
    long key;
    int value;
    if(mapStream.peek() != std::ifstream::traits_type::eof()) {
        while(mapStream.good()) {
            mapStream >> key;
            mapStream >> value;
            map[key] = value;
        }
    }

    PiercedStorage<double,long> data;

    data.setDynamicKernel(&(mesh.getCells()),PiercedVector<Cell>::SYNC_MODE_JOURNALED);

    data.fill(0.0);
    mesh.squeeze();
    mesh.sortCells();

//    for(auto & elem : map) {
//        std::cout << "rank " << rank << " " << elem.first << " " << elem.second << std::endl;
//    }

    std::vector<adaption::Info> partitionInfo = mesh.partition(map,true,true);
//    mesh.squeeze();
//    mesh.getCells().sync();
    std::string meshOut("meshafter");
    mesh.write(meshOut);
}


// =================================================================================== //

int main( int argc, char *argv[] ) {

    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
    MPI_Init(&argc,&argv);
    {
#endif
        test00002(argc,argv) ;
#if ENABLE_MPI==1
    }
    MPI_Finalize();
#endif

    return 0;

}
