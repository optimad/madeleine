/*
 * coupling.cpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#if ENABLE_MPI==1
#include <mpi.h>
#endif

#include "coupling.hpp"
#include "couplingUtils.hpp"

using namespace bitpit;

namespace coupling{

/*!
    Default Constructor
*/
MeshCoupling::MeshCoupling() : m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)), m_radius(1.0) {

};

/*!
    Constructor

    \param[in] inputNames list of the names of the input data
    \param[in] outputNames list of the names of the output data
*/
MeshCoupling::MeshCoupling(const vector<std::string> & inputNames, std::vector<std::string> & outputNames) : m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)), m_radius(1.0) {

    m_inputDataNames = inputNames;
    m_outputDataNames = outputNames;
};

/*!
    Initialize the unit sphere meshes, the radius and the radius scaled meshes. m_disciplineData is initialize to zero.

    \param[in] unitDisciplineMeshFile is the name of the file (.stl) containing the discipline mesh of the unit sphere
    \param[in] unitNeutralMeshFile is the name of the file (.stl) containing the neutral mesh of the unit sphere
    \param[in] radius is the value of the radius of the sphere discretized by the scaled meshes
*/
void MeshCoupling::initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile, double radius){

    //initialize radius
    m_radius = radius;

    //initialize discipline mesh
    m_unitDisciplineMesh->reset();
    m_disciplineData.clear(true);
    m_unitDisciplineMesh->importSTL(unitDisciplineMeshFile);
    m_unitDisciplineMesh->deleteCoincidentVertices();
    m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
    PiercedVector<Vertex> & scaledDisciplineVertices = m_scaledDisciplineMesh->getVertices();

    //initialize discipline data structure and resize the unit sphere, scaling its vertices
    for(Vertex & v : scaledDisciplineVertices) {
        v.setCoords(v.getCoords()*radius);
        m_disciplineData.insert(v.getId(),0.0);
    }

    //initialize neutral mesh and vertices container
    m_unitNeutralMesh->reset();
    m_unitNeutralMesh->importSTL(unitNeutralMeshFile);
    m_unitNeutralMesh->deleteCoincidentVertices();
    m_scaledNeutralMesh = PatchKernel::clone(m_unitNeutralMesh.get());
    PiercedVector<Vertex> & scaledNeutralVertices = m_scaledNeutralMesh->getVertices();

    //resize the neutral unit sphere mesh, scaling its vertices
    for(Vertex & v : scaledNeutralVertices) {
        v.setCoords(v.getCoords()*radius);
    }

    //print at log output vertices and cells of the discipline and neutral meshes
    log::cout() << " Discipline Mesh n vertices : " << m_unitDisciplineMesh->getVertexCount() << std::endl;
    log::cout() << " Discipline Mesh n cells : " << m_unitDisciplineMesh->getInternalCount() << std::endl;
    log::cout() << " Neutral Mesh n vertices : " << m_unitNeutralMesh->getVertexCount() << std::endl;
    log::cout() << " Neutral Mesh n cells : " << m_unitNeutralMesh->getInternalCount() << std::endl;

    //print on VTK files both discipline and neutral meshes
    std::string disciplineMeshVTKFileName = "disciplineMesh";
    std::string neutralMeshVTKFileName = "neutralMesh";
    writeMesh(m_scaledDisciplineMesh.get(),disciplineMeshVTKFileName);
    writeMesh(m_scaledNeutralMesh.get(),neutralMeshVTKFileName);

};

/*!
    Interpolate from neutral to discipline, compute smoothed step function, interpolate from discipline to neutral

    \param[in] neutralData the container with data from neutral mesh (data have to be coherent with the mesh)
*/
void MeshCoupling::compute(PiercedVector<double,long>* neutralData) {

    log::cout() << "First Interpolation" << std::endl;
    writeData(m_scaledNeutralMesh.get(),"neutralInitialization",neutralData,getInputDataNames());
    interpolateFromTo(m_scaledNeutralMesh.get(),neutralData,m_scaledDisciplineMesh.get(),&m_disciplineData);
    writeData(m_scaledDisciplineMesh.get(),"disciplineAfteInterpolation",&m_disciplineData,getOutputDataNames());

//    for(double & x : m_disciplineData) {
//        x = (tanh(10.0 * (x - 0.5)) + 1.0) / 2.0;
//    }
    log::cout() << "Computation" << std::endl;
    const PiercedVector<Vertex> & vertices = m_scaledDisciplineMesh->getVertices();
    for(const Vertex & v : vertices) {
        darray3 x = v.getCoords();
        double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        double datum = acos(x[0]/r) - M_PI/2;//sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        datum = sin(4*datum);
        datum *= datum;
        *(m_disciplineData.find(v.getId())) = datum;
    }

    writeData(m_scaledDisciplineMesh.get(),"disciplineAfterComputation",&m_disciplineData,getOutputDataNames());

    log::cout() << "Second Interpolation" << std::endl;
    interpolateFromTo(m_scaledDisciplineMesh.get(),&m_disciplineData,m_scaledNeutralMesh.get(),neutralData);
    writeData(m_scaledNeutralMesh.get(),"neutralAfterInterpolation",neutralData,getInputDataNames());

    log::cout() << "Exiting.." << std::endl;

};

/*!
    Initialize meshes and radius

    \return the vector containing the names of the input data
*/
const std::vector<std::string> & MeshCoupling::getInputDataNames() {

    return m_inputDataNames;

};

/*!
    Initialize meshes and radius

    \return the vector containing the names of the output data
*/
const std::vector<std::string> & MeshCoupling::getOutputDataNames() {

    return m_outputDataNames;

};

/*!
    Get scaled discipline mesh

    \return the scaled discipline mesh
*/
const SurfUnstructured* MeshCoupling::getDisciplineMesh(){

    return m_scaledDisciplineMesh.get();

};

/*!
    Get scaled neutral mesh

    \return the scaled neutral mesh
*/
SurfUnstructured* MeshCoupling::getNeutralMesh(){

    return m_scaledNeutralMesh.get();

};

/*!
    Possibly perform closing actions

*/
void MeshCoupling::close(){
    m_unitDisciplineMesh.reset(nullptr);
    m_unitNeutralMesh.reset(nullptr);
    m_scaledDisciplineMesh.reset(nullptr);
    m_scaledNeutralMesh.reset(nullptr);
};


}
