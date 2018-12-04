/*
 * utils.cpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#include "utils.hpp"
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
void interpolateFromTo(SurfUnstructured * fromMesh, const PiercedVector<double> & fromData, SurfUnstructured * toMesh, PiercedVector<double> & toData){

    SurfaceSkdTree fromTree(fromMesh);
    fromTree.build();

    const PiercedVector<Vertex> & toVertices = toMesh->getVertices();
    const PiercedVector<Vertex> & fromVertices = fromMesh->getVertices();
    PiercedVector<double>::iterator toDataItEnd = toData.end();
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

        // Interpolate data of Mesh 1 on vertex of Mesh 2
        double data = 0.;
        data += fromData[vIds[0]]*lambda[0];
        data += fromData[vIds[1]]*lambda[1];
        data += fromData[vIds[2]]*lambda[2];

        // Update data of Interpolation Mesh
        PiercedVector<double>::iterator toDataIt = toData.find(v.getId());

        if( toDataIt == toDataItEnd) {
            toData.insert(v.getId(), data);
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
void initDoubleDataOnMesh(const SurfUnstructured * mesh, PiercedVector<double>* data){

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
void writeData(SurfUnstructured * mesh,std::string filename,const PiercedVector<double> & data,const std::vector<std::string> & dataNames){

    mesh->getVTK().setName(filename);
    vector<double> vdata;
    vdata.reserve(data.size());
    const PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(const Vertex & v : vertices){
        //vdata.push_back(*(data.find(v.getId())));
        vdata.push_back(data.at(v.getId()));
    }
    mesh->getVTK().addData(dataNames[0], VTKFieldType::SCALAR, VTKLocation::POINT, vdata);
    mesh->write();

};


}



