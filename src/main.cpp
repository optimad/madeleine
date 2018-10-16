#if ENABLE_MPI==1
#    include <mpi.h>
#endif

#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

/*
 * Interpolation tool
 *
 * Usage:
 *
 * ./madeleine <path_mesh1>/<filename_mesh1.stl> <path_mesh2>/<filename_mesh2.stl>
 *
 * Default mesh files: ./data/sphere1.stl  ./data/sphere2.stl
 */
int main(int argc, char *argv[])
{
    //
    // Initialization
    //

    // Initialize input files (default/passed by arguments)
    std::string filemesh1 = "./data/sphere1.stl";
    std::string filemesh2 = "./data/sphere2.stl";
    if (argc > 2) {
        filemesh1 = argv[1];
        filemesh2 = argv[2];
    }

    // Initialize parallel (each process runs independently from the others)
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
    log::cout() << std::endl;


    // Initialize Mesh 1, vertices & data
    SurfUnstructured mesh1(2,3);
    mesh1.setExpert(true);
    mesh1.importSTL(filemesh1, false);

    PiercedVector<Vertex> vertices1;
    mesh1.deleteCoincidentVertices();
    vertices1 = mesh1.getVertices();

    PiercedVector<double> data1;
    for (auto v : vertices1){
        data1.insert(v.getId(), std::sin(v[0]));
    }


    // Initialize Mesh 2 & vertices
    SurfUnstructured mesh2(2,3);
    mesh2.setExpert(true);
    mesh2.importSTL(filemesh2, false);

    PiercedVector<Vertex> vertices2;
    mesh2.deleteCoincidentVertices();
    vertices2 = mesh2.getVertices();


    //
    // Run
    //

    // Output for Mesh 1 [Reference Mesh] & Mesh 2 [Interpolation Mesh]
    // Mesh 1 log info
    log::cout() << " First Mesh n vertices : " << mesh1.getVertexCount() << std::endl;
    log::cout() << " First Mesh n cells : " << mesh1.getInternalCount() << std::endl;

    // Write Mesh 1 with data (sphere1.vtu output file)
    {
        vector<double> vdata1;
        for (auto v : vertices1){
            vdata1.push_back(std::sin(v[0]));
        }
        mesh1.getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::POINT, vdata1) ;
        mesh1.getVTK().setName("mesh1");
        mesh1.write();
    }

    // Mesh 2 log info
    log::cout() << " Second Mesh n vertices : " << mesh2.getVertexCount() << std::endl;
    log::cout() << " Second Mesh n cells : " << mesh2.getInternalCount() << std::endl;



    // Build skd-tree of Mesh 1
    SurfaceSkdTree tree1(&mesh1);
    tree1.build();

    // Interpolate data from Reference Mesh [1] to Interpolation Mesh [2]
    // Barycentric coordinates of the projection of each node of Mesh 2 on the triangles of Mesh 1 are used
    PiercedVector<double> data2;
    for (auto v : vertices2){

        // Find closest triangle of Mesh 1 for each vertex of Mesh2
        darray3 x = v.getCoords();
        long id;
        double dist;
        tree1.findPointClosestCell(x, &id, &dist);

        // Recover vertices of triangle of Mesh 1
        ConstProxyVector<long> vIds = mesh1.getCell(id).getVertexIds();

        // Compute barycentric coordinates of the projected vertex of Interpolation Mesh on the closest triangle of Reference Mesh
        darray3 lambda;
        darray3 xP = CGElem::projectPointTriangle(x, vertices1[vIds[0]].getCoords(), vertices1[vIds[1]].getCoords(), vertices1[vIds[2]].getCoords(), lambda);

        // Interpolate data of Mesh 1 on vertex of Mesh 2
        double data = 0.;
        data += data1[vIds[0]]*lambda[0];
        data += data1[vIds[1]]*lambda[1];
        data += data1[vIds[2]]*lambda[2];

        // Update data of Interpolation Mesh
        data2.insert(v.getId(), data);

    }


    // Write Interpolation Mesh with interpolated data and exact data [analytically defined]
    {
        // Reference data on Mesh 2
        vector<double> dataRef;
        for (auto v : vertices2){
            dataRef.push_back(std::sin(v[0]));
        }

        // Add data to bitpit VTK object
        mesh2.getVTK().setName("mesh2");
        mesh2.getVTK().addData("dataRef", VTKFieldType::SCALAR, VTKLocation::POINT, dataRef) ;

        // Interpolated data [temporary copied in std::vector container to use base functions of bitpit VTK class]
        vector<double> vdata2;
        for (auto v : vertices2){
            vdata2.push_back(data2[v.getId()]);
        }

        // Add data to bitpit VTK object and write output file (sphere2.vtu)
        mesh2.getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::POINT, vdata2) ;
        mesh2.write();

    }


    //
    // Finalization
    //

#if ENABLE_MPI==1
    // MPI finalization
    MPI_Finalize();
#endif

}
