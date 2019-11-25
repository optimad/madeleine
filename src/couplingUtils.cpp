/*
 * utils.cpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#include "couplingUtils.hpp"
#include "bitpit_common.hpp"
#include "bitpit_IO.hpp"
#include "bitpit_CG.hpp"
#include "surfunstructured.hpp"
#include "commons.hpp"

using namespace bitpit;

namespace coupling {
///*!
//    Build a scaled mesh from the input mesh, scaling the input mesh vertex coordinates by the radius factor
//
//    \param[in] unitRadiusSphereMesh the mesh to be scaled
//    \param[in] radius the scaling factor
//    \return the scaled mesh
//*/
//std::unique_ptr<SurfUnstructured> scale(const SurfUnstructured & unitRadiusSphereMesh, double radius) {
//
//    std::unique_ptr<SurfUnstructured> mesh = PatchKernel::clone(&unitRadiusSphereMesh);
//
//    PiercedVector<Vertex> & vertices = mesh->getVertices();
//
//    for(auto & v : vertices) {
//        v.setCoords(v.getCoords()*radius);
//    }
//
//    return mesh;
//}

/*!
    Interpolate data from the "fromMesh" to the "toMesh"

    \param[in] fromMesh the interpolation origin mesh
    \param[in] fromData the data of the origin to be interpolated
    \param[in] toMesh the interpolation destination mesh
    \param[in] toData the interpolated data on the destination mesh
    \return the scaled mesh
*/
void interpolateFromTo(SurfUnstructured * fromMesh, PiercedVector<double> * fromData, SurfUnstructured * toMesh, PiercedVector<double> * toData){

    log::cout() << "Build Tree" << std::endl;

    SurfaceSkdTree fromTree(fromMesh);
    log::cout() << "Tree declared" << std::endl;
    fromTree.build(1);

    log::cout() << "Tree built" << std::endl;


    const PiercedVector<Vertex> & toVertices = toMesh->getVertices();
    const PiercedVector<Vertex> & fromVertices = fromMesh->getVertices();
    PiercedVector<double>::iterator toDataItEnd = toData->end();
    for (const Vertex & v : toVertices){

        // Find closest triangle of Mesh from for each vertex of Mesh to
        darray3 x = v.getCoords();
        long id;
        double dist;
        fromTree.findPointClosestCell(x, &id, &dist);

        // Recover vertices of triangle of Mesh 1
        ConstProxyVector<long> vIds = fromMesh->getCell(id).getVertexIds();

        // Compute barycentric coordinates of the projected vertex of Interpolation Mesh on the closest triangle of Reference Mesh
        darray3 lambda;
        darray3 xP = bitpit::CGElem::projectPointTriangle(x, fromVertices[vIds[0]].getCoords(), fromVertices[vIds[1]].getCoords(), fromVertices[vIds[2]].getCoords(), lambda);
        BITPIT_UNUSED(xP);
        // Interpolate data of Mesh 1 on vertex of Mesh 2
        double data = 0.;
        data += (*fromData)[vIds[0]]*lambda[0];
        data += (*fromData)[vIds[1]]*lambda[1];
        data += (*fromData)[vIds[2]]*lambda[2];

        // Update data of Interpolation Mesh
        PiercedVector<double>::iterator toDataIt = toData->find(v.getId());

        if( toDataIt == toDataItEnd) {
            toData->insert(v.getId(), data);
            throw std::runtime_error("Warning! New element inserted!");
        }
        else {
            *toDataIt = data;
        }
    }


};

/*!
    Initialize a mesh-coherent PiercedVector Container with sin^2(sqrt(x_0^2+x_1^2+x_2^2))

    \param[in] mesh the mesh
    \return the initialized PiercedVector
*/
void initDoubleDataOnMesh(SurfUnstructured * mesh, PiercedVector<double>* data){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(const Vertex & v: vertices) {
        darray3 x = v.getCoords();
        double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        double datum = acos(x[2]/r) - M_PI/2;//sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        datum = sin(4*datum);
        datum *= datum;
        data->insert(v.getId(),datum);
    }
};

/*!
    Initialize a mesh-coherent PiercedVector Container with a C-array.
    It has to be noted that the array should contains data with the same order of the vertices of the associated mesh.

    \param[in] mesh the mesh
    \param[out] data PiercedVector container associated to the mesh that has to be filled
    \param[in] array a pointer to a C-array containing data to be inserted into the PiercedVector
    \param[in] arraySize the size of the C-array
*/
void initDataOnMeshFromArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    assert(vertices.size()==arraySize);
    size_t count = 0;
    std::cout << "Start inserting ...";
    for(const Vertex & v: vertices) {
        data->insert(v.getId(),array[count]);
        ++count;
    }
};

/*!
    Move data from a mesh-coherent PiercedVector Container to a C-array.
    It has to be noted that the array will contain data with the same order of the vertices of the associated mesh.

    \param[in] mesh the mesh
    \param[in] data PiercedVector container associated to the mesh containing data to be inserted into the C-array
    \param[out] array a pointer to a C-array to be filled with data coming from the PiercedVector
    \param[in] arraySize the size of the C-array
*/
void moveDataOnMeshToArray(SurfUnstructured * mesh, PiercedVector<double>* data, double* array, size_t arraySize){

    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    assert(vertices.size()==arraySize);
    size_t count = 0;
    for(const Vertex & v: vertices) {
        array[count] = data->at(v.getId());
        ++count;
    }
};



/*!
    Write VTK file( (p)vtu ) containing the mesh

    \param[in] mesh the mesh
    \param[in] filename the name of the .(p)vtu file(s)
*/
void writeMesh(SurfUnstructured * mesh,std::string filename){
    mesh->getVTK().setName(filename);
    mesh->write();

};

/*!
    Write VTK file( (p)vtu ) containing the mesh and the passed data

    \param[in] mesh the mesh
    \param[in] filename the name of the .(p)vtu file(s)
    \param[in] data the PiercedVector containing the data to be plotted. NB only one scalar field is allowed
    \param[in] dataNames a list of the field data names. NB its size has to be 1.
*/
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> * data,const std::vector<std::string> & dataNames){

    mesh->getVTK().setName(filename);
    std::vector<double> vdata;
    vdata.reserve(data->size());
    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(const Vertex & v : vertices){
        //vdata.push_back(*(data.find(v.getId())));
        vdata.push_back(data->at(v.getId()));
    }
    mesh->getVTK().addData(dataNames[0], VTKFieldType::SCALAR, VTKLocation::POINT, vdata);
    mesh->write();

};

/*!
    Compute an ordered mesh partitioning, filling an empty long array (provided by the user) with the cell indices to be assigned to each rank.
    The cell numbering is given by the ordering in the mesh file. CAVEAT: The maximum global number of cells is bounded to 2 billions, due to the use of INT.
    bitpit can manage more cells by using special MPI datatype.

    \param[in] meshFile the file containing the mesh
    \param[in] comm the communicator used to partition the mesh
    \param[out] cellSizesPerRank a pointer to an empty long array provided by the user and filled with the cell indices for each rank.
*/
void computeMeshFilePartitioning(const std::string meshFile,std::vector<int> & idRanks,MPI_Comm comm){

    //Ask communicator for its size
    int nofRanks;
    int rank;
    MPI_Comm_size(comm,&nofRanks);
    MPI_Comm_rank(comm,&rank);

    std::vector<int> sizes(nofRanks,0);
    int nofCells = 0;
    if(rank == 0){
        //Rank 0 reads mesh file and counts the cells
        SurfUnstructured mesh(2,3);
        mesh.importSTL(meshFile);
        mesh.deleteCoincidentVertices();
        nofCells = mesh.getCellCount();

        //Rank 0 compute the number of cells per rank
        int integerDivision = nofCells / nofRanks;
        int divisionReminder = nofCells % nofRanks;
        for(int r = 0; r < nofRanks; ++r) {
            sizes[r] = integerDivision;
        }
        for(int i = 0; i < divisionReminder; ++i) {
            ++sizes[i];
        }

        idRanks.resize(nofCells,0);
        int count = 0;
        for(int r = 0; r < nofRanks; ++r) {
            for(int i = 0; i < sizes[r]; ++i){
                idRanks[count] = r;
                ++count;
            }
        }
    }

    //Rank 0 broadcast nofCells
    MPI_Bcast(&nofCells,1,MPI_INT,0,comm);
    idRanks.resize(nofCells,0);
    MPI_Bcast(idRanks.data(),nofCells,MPI_INT,0,comm);

//    //DEBUG
//    for(int r = 0; r < nofRanks; ++r) {
//        if(rank == r) {
//            std::cout << "I'm rank " << rank << " and my idRanks is : " << std::endl;
//            for(int i = 0; i < nofCells; ++i) {
//                std::cout << "id " << i << " -> rank " << idRanks[i] << std::endl;
//            }
//        }
//        MPI_Barrier(comm);
//    }
//    //DEBUG
};


}


