/*
 * coupling.cpp
 *
 *  Created on: 23 nov 2018
 *      Author: marco
 */

#if ENABLE_MPI==1
#include <mpi.h>
#include "metis.h"
#endif

#include "coupling.hpp"
#include "couplingUtils.hpp"
#include <fstream>
#include <bitpit_operators.hpp>

using namespace bitpit;

namespace coupling{

/*!
    Default Constructor
*/
#if ENABLE_MPI==1
MeshCoupling::MeshCoupling(std::string disciplineName, MPI_Comm comm) :
        m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)),
        m_disciplineData(2), m_neutralData(2),
        m_fid_temperatureIN(0), m_fid_temperatureOUT(1),
        m_disciplineRadius(1.0), m_oldDisciplineRadius(1.0), m_neutralRadius(1.0), m_oldNeutralRadius(1.0),
        m_name(disciplineName), m_system(nullptr), m_thickness(1.0), m_innerSphere(true), m_sourceMaxIntensity(0.0),
        m_sourceDirection({{1.0,0.0,0.0}})
#else
MeshCoupling::MeshCoupling(std::string disciplineName) :
        m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)),
        m_neutralData(2), m_disciplineData(2),
        m_disciplineRadius(1.0), m_oldDisciplineRadius(1.0), m_neutralRadius(1.0), m_oldNeutralRadius(1.0),
        m_name(disciplineName), m_system(nullptr), m_thickness(1.0), m_innerSphere(true), m_sourceMaxIntensity(0.0),
        m_sourceDirection({{1.0,0.0,0.0}})
#endif
{
#if ENABLE_MPI==1
    //duplicate the workers communicator
    MPI_Comm_dup(comm,&m_comm);

    //Set rank and number of processor
    MPI_Comm_size(m_comm,&m_nprocs);
    MPI_Comm_rank(m_comm,&m_rank);


    //Set workers communicator name
    char workers_comm_name_char[MPI_MAX_OBJECT_NAME];
    std::string workers_comm_name = m_name + "_COMM";
    strcpy(workers_comm_name_char,workers_comm_name.c_str());
    MPI_Comm_set_name(m_comm,workers_comm_name_char);

    /* char worldnameout[MPI_MAX_OBJECT_NAME]; */
    /* char workersnameout[MPI_MAX_OBJECT_NAME]; */
    /* int rlen; */
    /* MPI_Comm_get_name( MPI_COMM_WORLD, worldnameout, &rlen ); */
    /* MPI_Comm_get_name( m_comm, workersnameout, &rlen ); */
    /* std::string workersName(workersnameout); */
    /* std::string worldName(worldnameout); */
    /* int worldSize; */
    /* int worldRank; */
    /* MPI_Comm_size(MPI_COMM_WORLD,&worldSize); */
    /* MPI_Comm_rank(MPI_COMM_WORLD,&worldRank); */
    /* std::cout << "I'm " << m_rank << " of " << m_nprocs  << " on " << workersName  << " and " << worldRank << " of " << worldSize << " on " << worldName << std::endl; */

    m_neutralTag = 0;
    m_disciplineTag = 0;

    m_lbCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(m_comm));

#else
    m_rank = 0;
    m_nprocs = 1;
#endif
    m_kernel = -1;
    m_inputField = -1;
    m_outputField = -1;

    m_switchOnOutput = false;
};

/*!
    Constructor

    \param[in] inputNames list of the names of the input data
    \param[in] outputNames list of the names of the output data
*/
#if ENABLE_MPI==1
MeshCoupling::MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName, MPI_Comm comm) :
        MeshCoupling(disciplineName, comm)
#else
MeshCoupling::MeshCoupling(const std::vector<std::string> & inputNames, std::vector<std::string> & outputNames, std::string disciplineName) :
        MeshCoupling(disciplineName)
#endif
{
    m_inputDataNames = inputNames;
    m_outputDataNames = outputNames;
    m_system = std::unique_ptr<StencilScalarSolverHandler>(new StencilScalarSolverHandler(disciplineName,false));
    m_jacobianManager = std::unique_ptr<JacobianMatricesManager>(new JacobianMatricesManager(disciplineName,comm,*(m_system.get()) ));
};


void MeshCoupling::initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile,
        double disciplineRadius, double neutralRadius, double thickness, bool innerSphere, double sourceIntensity, std::vector<double> sourceDirection,
        double thermalDiffusivityCoefficient, double emissivity, double infinityTemperature,
        const std::vector<int> & globalNeutralId2MeshFileRank){
    int kernel = 1;
    initialize(unitDisciplineMeshFile, unitNeutralMeshFile, disciplineRadius, neutralRadius, thickness, innerSphere, sourceIntensity, sourceDirection,
            thermalDiffusivityCoefficient, emissivity, infinityTemperature,
            globalNeutralId2MeshFileRank, kernel);
}

/*!
    Initialize the unit sphere meshes, the radius and the radius scaled meshes. m_disciplineData is initialize to zero.

    \param[in] unitDisciplineMeshFile is the name of the file (.stl) containing the discipline mesh of the unit sphere
    \param[in] unitNeutralMeshFile is the name of the file (.stl) containing the neutral mesh of the unit sphere
    \param[in] radius is the value of the radius of the sphere discretized by the scaled meshes
    \param[in] thickness is the value of the radius of the sphere discretized by the scaled meshes
    \param[in] innerSphere is true if the discipline surface is the inner one
    \param[in] sourceIntensity is the value of the external constant source flux intensity (only for outer discipline)
    \param[in] sourceDirection is an array for the direction of the external constant source flux
*/
void MeshCoupling::initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile,
        double disciplineRadius, double neutralRadius, double thickness, bool innerSphere, double sourceIntensity, std::vector<double> sourceDirection,
        double thermalDiffusivityCoefficient, double emissivity, double infinityTemperature,
        const std::vector<int> & globalNeutralId2MeshFileRank, int kernel){

    //initialize radius
    m_disciplineRadius = disciplineRadius;
    m_neutralRadius = neutralRadius;
    m_thickness = thickness;
    m_innerSphere = innerSphere;
    m_sourceMaxIntensity = sourceIntensity;
    m_thermalDiffusivityCoefficient = thermalDiffusivityCoefficient;
    m_emissivity = emissivity;
    m_infinityTemperature = infinityTemperature;

    assert(sourceDirection.size() >= 3);
    for(int i = 0; i < 3; ++i) {
        m_sourceDirection[i] = sourceDirection[i];
    }
    m_sourceDirection /= norm2(sourceDirection);

    m_kernel = kernel;

    //initialize mesh file names
    m_disciplineMeshFileName = unitDisciplineMeshFile;
    m_neutralMeshFileName = unitNeutralMeshFile;

    //initialize neutral mesh map ids -> ranks for mesh file ordered partitioning
    m_globalNeutralId2MeshFileRank = globalNeutralId2MeshFileRank;

#if ENABLE_MPI==1
    //Read the unit meshes
    readUnitDisciplineMesh();
    readUnitNeutralMesh();

    std::unique_ptr<SurfUnstructured> serialNeutralMesh = PatchKernel::clone(m_unitNeutralMesh.get());

    //Partition the neutral mesh by mesh file order
    staticPartitionNeutralMeshByMeshFileOrder();

    //Compute the discipline ids/ranks map to get discipline rank sub-domains overlapping the neutral ones with the same rank
    computeGlobalDisciplineId2NeutralRank();

    computeDiscipline2FileNeutralCellPerRanks(); //compute static m_discipline2FileNeutralCellPerRanks: from D_serial to D_Nf

    staticPartitionDisciplineMeshByNeutralFile(); //D_Nf

    computeGlobalNeutralId2DisciplineRank(serialNeutralMesh.get()); // compute global vector to get N_{D_Nf}

    computeNeutralId2DisciplineCellPerRanks(); //compute dynamic m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks: from Nf to N_{D_Nf}

    buildScaledMeshes();
    //ATTENTION from here on unit meshes should not be used anymore.

    //Initialize data container
    synchronizeMeshData();
    uniformlyInitAllData(294.0);
    //Set mesh VTK writer
    //ATTENTION from here on scaled meshes are written with attached data, therefore data container have to be coherent with the relative mesh
    prepareWritingData();

    //Initialize ghost communicators
    initializeGhostCommunicators();

    updateDisciplineGhosts();

//    std::string name = "D_Nf";
//    m_scaledDisciplineMesh->write(name);


    //Prepare discipline
    //set Input/Output field access indices to be interpolated from/to neutral
    //Both discipline take temperature as input as give temperature as output
    m_inputField = m_fid_temperatureIN;
    m_outputField = m_fid_temperatureOUT;

    //Prepare Linear System
    //Matrix preallocation assembly
    m_helmoltzStencils.resize(getDisciplineMesh()->getInternalCount());
    assemblySimplifiedDiscreteHelmholtzSystem();
    updateSystemRHS(); //qui non serve solo testing

    getJacobianManager()->computeInterpolationMatrices(m_scaledDisciplineMesh.get(),&m_disciplineNumberingInfo,
            m_scaledNeutralMesh.get(),&m_neutralNumberingInfo);

#else
    m_unitDisciplineMesh->reset();
    m_disciplineData.clear(true);
    m_unitDisciplineMesh->importSTL(unitDisciplineMeshFile);
    m_unitDisciplineMesh->deleteCoincidentVertices();
    m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
#endif

};

/*!
    Interpolate from neutral(Nf) to discipline(D_Nf), compute smoothed step function on discipline(D_Nf), partition neutral mesh(N_{D_Nf}), interpolate from discipline(D_Nf) to neutral(N_{D_Nf})

    \param[in/out] neutralInputArray contiguous C-array from NUMPY array ordered like id-ordered cells in Nf mesh file partitioning
    \param[in] newRadius the radius of the discipline sphere coming from GEMS computations
    \param[in] ohertRadius the radius of the other discipline sphere coming form GEMS computations (currently not used)
*/


void MeshCoupling::compute(double *neutralInputArray, std::size_t size, double newRadius, double otherRadius) {

    //Set new radii and rescale meshes
    BITPIT_UNUSED(otherRadius);
    scaleMeshToRadius(m_scaledDisciplineMesh,m_disciplineRadius,newRadius);
    scaleMeshToRadius(m_scaledNeutralMesh,m_neutralRadius,newRadius);

    log::cout() << "First Interpolation" << std::endl;
    //sort cells by id - neutalInputArray should have values ordered like the neutral mesh file partitioning
    m_scaledNeutralMesh->sortCells();

    updateInputField(neutralInputArray,size);

    std::string name;
    if(m_switchOnOutput) {
        name = "initilizedNF";
        m_scaledNeutralMesh->write(name);
    }

    //Update neutral ghost cell values
    //It should not be the case but call m_neutralGhostCommunicator->resetExchangeLists() if neutral mesh has changed.

    updateNeutralGhosts();

    //Interpolate from N_f to D_{N_f]
    std::cout << "N_f to D_{N_f} interpolation." << std::endl;
    interpolateFromTo(m_scaledNeutralMesh.get(),&m_neutralData,m_scaledDisciplineMesh.get(),&m_disciplineData, m_inputField);

    if(m_switchOnOutput) {
        name = "intepolatedD_Nf";
        m_scaledDisciplineMesh->write(name);
    }

    //Update discipline ghost cell values
    std::cout << "Updatings discipline ghosts" << std::endl;
    updateDisciplineGhosts();

    //Solve Radiation problem
    std::cout << "Solve radiation problem" << std::endl;
    disciplineKernel();

    updateDisciplineGhosts();

    if(m_switchOnOutput) {
        name = "solutionD_Nf";
        m_scaledDisciplineMesh->write(name);
    }

//    name = "Nf";
//    m_scaledNeutralMesh->write(name);
    //dynamicPartitionNeutralMeshByNeutralMeshFilePartitionedDiscipline();

    //Interpolate from D_{N_f} to N_{D_{N_f}}
    std::cout << "D_{N_f} to N_{D_{N_f}} interpolation." << std::endl;
    interpolateFromTo(m_scaledDisciplineMesh.get(),&m_disciplineData,m_scaledNeutralMesh.get(),&m_neutralData,m_outputField);
    //This interpolation is for data check only. Input data will be overwritten at the beginning of the next compute call
    interpolateFromTo(m_scaledDisciplineMesh.get(),&m_disciplineData,m_scaledNeutralMesh.get(),&m_neutralData,m_inputField);

    if(m_switchOnOutput) {
        name = "interpolatedSolutionN_{D_Nf}";
        m_scaledNeutralMesh->write(name);
    }

    //dynamicPartitionNeutralMeshByNeutralMeshWithData();

    if(m_switchOnOutput) {
        name = "unsortedNf";
        m_scaledNeutralMesh->write(name);
    }

    m_scaledNeutralMesh->sortCells();

    if(m_switchOnOutput) {
        name = "sortedNf";
        m_scaledNeutralMesh->write(name);
    }

    //Update C-array to pass data to NUMPY array used by GEMS
    std::size_t counter = 0;
    for(const Cell & cell : m_scaledNeutralMesh->getCells()) {
        long id = cell.getId();
        if(cell.isInterior()) {
            neutralInputArray[counter] = m_neutralData.at(id, m_outputField);
            ++counter;
        }
    }

    //Neutral mesh numbering info is needed in computeJacobianRow. No need to update it after any partitioning. But DO NOT use it if old.
    m_neutralNumberingInfo.update();

    getJacobianManager()->computeEllipticOperatorInverse(m_scaledDisciplineMesh.get(),m_disciplineNumberingInfo);
    getJacobianManager()->computeOutputInputJacobian(m_comm,m_emissivity,m_innerSphere);
    getJacobianManager()->computeOutputControlJacobian(m_comm, m_emissivity, m_innerSphere, m_disciplineRadius,
            m_disciplineData, m_scaledDisciplineMesh.get(), m_disciplineNumberingInfo,
            m_scaledNeutralMesh.get(), m_neutralNumberingInfo);
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
    Get scaled neutral mesh number of internal cells

    \return the scaled neutral mesh size
*/
size_t MeshCoupling::getNeutralMeshSize(){

    return m_scaledNeutralMesh->getInternalCount();

};


/*!
    Get scaled neutral mesh first cell Id

    \return the scaled neutral mesh first id
*/
long MeshCoupling::getNeutralFirstCellId(){

    return m_scaledNeutralMesh->getCells().front().getId();

};



/*!
    Possibly perform closing actions

*/
void MeshCoupling::close(){
//    m_unitDisciplineMesh.reset(nullptr);
//    m_unitNeutralMesh.reset(nullptr);
//    m_scaledDisciplineMesh.reset(nullptr);
//    m_scaledNeutralMesh.reset(nullptr);
};

/*!
    Read unit discipline mesh from file
*/
void MeshCoupling::readUnitDisciplineMesh() {

    //Rank 0 initializes discipline mesh
    m_unitDisciplineMesh->reset();
    if(m_rank == 0) {
        //m_disciplineData.clear(true);
        m_unitDisciplineMesh->importSTL(m_disciplineMeshFileName);
        m_unitDisciplineMesh->deleteCoincidentVertices();
        std::cout << "Discipline number of cells = " << m_unitDisciplineMesh->getCellCount() << std::endl;
    }
    m_unitDisciplineMesh->buildAdjacencies();

}

/*!
    Read unit neutral mesh from file
*/
void MeshCoupling::readUnitNeutralMesh() {

    //Rank 0 reads unit neutral mesh
    m_unitNeutralMesh->reset();
    if(m_rank == 0) {
        m_unitNeutralMesh->importSTL(m_neutralMeshFileName);
        m_unitNeutralMesh->deleteCoincidentVertices();
        std::cout << "Neutral number of cells = " << m_unitNeutralMesh->getCellCount() << std::endl;
    }
    m_unitNeutralMesh->buildAdjacencies();
}


#if ENABLE_MPI==1
/*!
    Statically partition the mesh dividing cells by using metis
*/
void MeshCoupling::staticPartitionMeshByMetis(SurfUnstructured & mesh) {

    std::unordered_map<long,int> cellRanks = computeStaticPartitionByMetis(mesh);
    std::vector<adaption::Info> partitionInfo = mesh.partition(cellRanks,true,true);
}

/*!
    Compute the static partition of the mesh by using metis. Returns an unordered_map<id,rank> filled only by rank 0.
    Rank 0 can use this unordered_map to distribute cells among the other ranks in the discipline communicator
*/
std::unordered_map<long,int> MeshCoupling::computeStaticPartitionByMetis(SurfUnstructured & mesh) {


    bitpit::PatchNumberingInfo patchInfo = bitpit::PatchNumberingInfo(&mesh);

    std::unordered_map<long,int> cellRanks;
    mesh.buildAdjacencies();

    if(m_rank == 0) {
        //build an unordered map consecutive index to id
        std::unordered_map<long,long> mapcell;
        const PiercedVector<Cell> & cells = mesh.getCells();
        for(const Cell & cell : cells) {
            long id = cell.getId();
            long consecutive = patchInfo.getCellConsecutiveId(id);
            mapcell[consecutive] = id;
        }

        // Use cell centers as vertices of graph

        //
        //  The number of vertices.
        //
        idx_t nvtxs = idx_t(mesh.getCellCount());
        //DEBUG
        std::cout << "Cell count = " << nvtxs << std::endl;
        //DEBUG

        //
        // Number of balancing constraints, which must be at least 1.
        //
        idx_t ncon = 1;


        //Build Adjacencies
        std::vector<idx_t> xadj(nvtxs+1);

        idx_t duenedges = 0;
        for (int i=0; i < mesh.getInternalCount(); i++){
            duenedges += mesh.getCell(mapcell[i]).getAdjacencyCount();
        }

        std::vector<idx_t> adjncy(duenedges);

        idx_t i;
        idx_t j=0;
        std::vector<long> neighs;
        for (i=0; i < mesh.getInternalCount(); i++){
            xadj[i] = j;
            neighs.clear();
            mesh.findCellFaceNeighs(mapcell[i], &neighs);
            for (long id : neighs){
                if (id > 0){
                    adjncy[j] = idx_t(patchInfo.getCellConsecutiveId(id));
                    j++;
                }
            }
        }
        xadj[i] = j;

        idx_t nParts = m_nprocs;

        //
        //  On return, the edge cut volume of the partitioning solution.
        //
        idx_t objval;

        //
        //  On return, the partition vector for the graph.
        //
        std::vector<idx_t> part(nvtxs);

        int ret = METIS_PartGraphKway ( &nvtxs, &ncon, xadj.data(), adjncy.data(), NULL, NULL,
                NULL, &nParts, NULL, NULL, NULL, &objval, part.data() );
        BITPIT_UNUSED(ret);

//        //DEBUG
//        std::cout << "Metis partition method error flag = " << ret << std::endl;
//        //DEBUG

        cellRanks.reserve(nvtxs);
        for (long i=0; i<nvtxs; i++){
            cellRanks[mapcell[i]] = part[i];
        }
//        //DEBUG
//        int k = 0;
//        for(auto & r : cellRanks){
//            std::cout << "id = " << r.first << " - rank = " << r.second << std::endl;
//            ++k;
//        }
//        //DEBUG
    }

    return cellRanks;
}
#endif

/*!
    Statically partition of the discipline mesh by using metis.
*/
void MeshCoupling::staticPartitionDisciplineMeshByMetis() {

    //Partition the discipline mesh
    staticPartitionMeshByMetis(*(m_unitDisciplineMesh.get()));

}

/*!
    Statically partition of the neutral mesh by using mesh file order.
*/
void MeshCoupling::staticPartitionNeutralMeshByMeshFileOrder() {

    std::unordered_map<long,int> staticNeutralMeshFilePartitionCellPerRank;
    if(m_rank == 0) {
        for(Cell & cell : m_unitNeutralMesh->getCells()) {
            long id = cell.getId();
            staticNeutralMeshFilePartitionCellPerRank[id] = m_globalNeutralId2MeshFileRank[id];
        }
    }
    //Mesh file ordered unit neutral mesh partitioning
    m_unitNeutralMesh->setCommunicator(m_comm);
    std::vector<adaption::Info> partitionInfo = m_unitNeutralMesh->partition(staticNeutralMeshFilePartitionCellPerRank,true,true);

}

/*!
    Compute a global map ids->ranks for the neutral mesh in order to partition the neutral mesh
    with rank sub-domains overlapping the discipline sub-domains with the same rank.
*/
void MeshCoupling::computeGlobalNeutralId2DisciplineRank(SurfUnstructured *serialNeutralMesh) {
    //Map the discipline rank distribution on neutral mesh
    //Rank 0 communicate neutral mesh ids and cell centers to the other ranks
    long *idxNeutral = nullptr; //ordered by neutral mesh id
    double * cellCentersNeutral = nullptr; //ordered by neutral mesh id
    int nofCellsNeutral = 0;
    if(m_rank == 0) {
        nofCellsNeutral = serialNeutralMesh->getCellCount();
        idxNeutral = new long[nofCellsNeutral];
        cellCentersNeutral = new double[nofCellsNeutral * 3];
        for(auto cell : serialNeutralMesh->getCells()) {
            long id = cell.getId();
            const std::array<double,3> & cellCenter = serialNeutralMesh->evalCellCentroid(id);
            idxNeutral[id] = id;
            for(int i = 0; i < 3; ++i) {
                cellCentersNeutral[id * 3 + i] = cellCenter[i];
            }
        }
    }
    MPI_Bcast(&nofCellsNeutral,1,MPI_INT,0,m_comm);
    if(m_rank != 0) {
        idxNeutral = new long[nofCellsNeutral];
        cellCentersNeutral = new double [nofCellsNeutral * 3];
    }
    MPI_Bcast(idxNeutral,nofCellsNeutral,MPI_LONG,0,m_comm);
    MPI_Bcast(cellCentersNeutral,3*nofCellsNeutral,MPI_DOUBLE,0,m_comm);

    //Project rank from discipline partition to cellRanks for neutral mesh
    //Build discipline mesh surface skd tree
    SurfaceSkdTree disciplineTree(m_unitDisciplineMesh.get());
    disciplineTree.build(1);
    //Initialiaze rank maps - array ordered by id and containing ranks
    m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank.resize(nofCellsNeutral,-1);
    //Loop over neutral cell center and check the discipline cell the neutral cell center belong to. Assign rank according to
    std::array<double,3> tempCellCenter;
    darray3 lambda;
    darray3 xP;
    BITPIT_UNUSED(lambda);
    bitpit::utils::DoubleFloatingEqual checkFloating;
    for(int i = 0; i < nofCellsNeutral; ++i) {
        BITPIT_UNUSED(xP);
        for(int j = 0; j < 3; ++j) {
            tempCellCenter[j] = cellCentersNeutral[i * 3 + j];
        }
        //compute closest cell distance
        long id;
        double closestCellDist, closestCellPlaneDist;
        disciplineTree.findPointClosestCell(tempCellCenter, &id, &closestCellDist);
        //compute closest cell plane distance
        if(m_unitDisciplineMesh->getCell(id).isInterior()) {
            closestCellPlaneDist = bitpit::CGElem::distancePointPlane(tempCellCenter,
                    m_unitDisciplineMesh->getCellVertexCoordinates(id)[0],m_unitDisciplineMesh->evalFacetNormal(id),xP);

            if(checkFloating(closestCellDist,closestCellPlaneDist)) {
                m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank[i] = m_rank;
            }
        }
    }
    //Reduce disciplineMappedRankForNeutral on Rank 0 for static partitioning
    MPI_Allreduce(MPI_IN_PLACE,m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank.data(),nofCellsNeutral,MPI_INT,MPI_MAX,m_comm);

    delete [] idxNeutral;
    delete [] cellCentersNeutral;
}

/*!
    Compute a global map ids->ranks for the discipline mesh in order to partition the discipline mesh
    with rank sub-domains overlapping the neutral sub-domains with the same rank.
*/
void MeshCoupling::computeGlobalDisciplineId2NeutralRank() {
    //Map the neutral rank distribution on discipline mesh
    //Rank 0 communicate discipline mesh ids and cell centers to the other ranks
    long *idxDiscipline = nullptr; //ordered by discipline mesh id
    double * cellCentersDiscipline = nullptr; //ordered by discipline mesh id
    int nofCellsDiscipline = 0;
    if(m_rank == 0) {
        nofCellsDiscipline = m_unitDisciplineMesh->getCellCount();
        idxDiscipline = new long[nofCellsDiscipline];
        cellCentersDiscipline = new double[nofCellsDiscipline * 3];
        for(auto cell : m_unitDisciplineMesh->getCells()) {
            long id = cell.getId();
            const std::array<double,3> & cellCenter = m_unitDisciplineMesh->evalCellCentroid(id);
            idxDiscipline[id] = id;
            for(int i = 0; i < 3; ++i) {
                cellCentersDiscipline[id * 3 + i] = cellCenter[i];
            }
        }
    }
    MPI_Bcast(&nofCellsDiscipline,1,MPI_INT,0,m_comm);
    if(m_rank != 0) {
        idxDiscipline = new long[nofCellsDiscipline];
        cellCentersDiscipline = new double [nofCellsDiscipline * 3];
    }
    MPI_Bcast(idxDiscipline,nofCellsDiscipline,MPI_LONG,0,m_comm);
    MPI_Bcast(cellCentersDiscipline,3*nofCellsDiscipline,MPI_DOUBLE,0,m_comm);

    //Project rank from neutral partition to cellRanks for discipline mesh
    //Build neutral mesh surface skd tree
    SurfaceSkdTree neutralTree(m_unitNeutralMesh.get());
    neutralTree.build(1);
    //Initialiaze rank maps - array ordered by id and containing ranks
    m_globalDisciplineId2NeutralMeshFileRank.resize(nofCellsDiscipline,-1);
    std::vector<int> candidates(nofCellsDiscipline,-1);
    //Loop over discipline cell center and check the neutral cell the discipline cell center belong to. Assign rank according to
    std::array<double,3> tempCellCenter;
    darray3 lambda;
    darray3 xP;
    BITPIT_UNUSED(lambda);
    bitpit::utils::DoubleFloatingEqual checkFloating;
    for(int i = 0; i < nofCellsDiscipline; ++i) {
        BITPIT_UNUSED(xP);
        for(int j = 0; j < 3; ++j) {
            tempCellCenter[j] = cellCentersDiscipline[i * 3 + j];
        }
        //compute closest cell distance
        long id;
        double closestCellDist = 0.0, closestCellPlaneDist = 0.0;
        neutralTree.findPointClosestCell(tempCellCenter, &id, &closestCellDist);
        //compute closest cell plane distance
        if(m_unitNeutralMesh->getCell(id).isInterior()) {
            closestCellPlaneDist = bitpit::CGElem::distancePointPlane(tempCellCenter,
                    m_unitNeutralMesh->getCellVertexCoordinates(id)[0],m_unitNeutralMesh->evalFacetNormal(id),xP);
            if(checkFloating(closestCellDist,closestCellPlaneDist)) {
                m_globalDisciplineId2NeutralMeshFileRank[i] = m_rank;
            } else {
                candidates[i] = m_rank;
            }
        }
    }
    //Reduce disciplineMappedRankForNeutral on Rank 0 for static partitioning
    MPI_Allreduce(MPI_IN_PLACE,m_globalDisciplineId2NeutralMeshFileRank.data(),nofCellsDiscipline,MPI_INT,MPI_MAX,m_comm);
    MPI_Allreduce(MPI_IN_PLACE,candidates.data(),nofCellsDiscipline,MPI_INT,MPI_MAX,m_comm);
    //check cellcenters candidates
    for(int c = 0; c < nofCellsDiscipline; ++c) {
        if(m_globalDisciplineId2NeutralMeshFileRank[c] == -1) {
            m_globalDisciplineId2NeutralMeshFileRank[c] = candidates[c];
        }
    }

}

/*!
    Compute a local map ids->ranks for the discipline mesh in order to partition the discipline mesh
    with rank sub-domains overlapping the neutral sub-domains with the same rank.
*/
void MeshCoupling::computeDiscipline2FileNeutralCellPerRanks() {

    for(const Cell & cell : m_unitDisciplineMesh->getCells()) {
        if(cell.isInterior()) {
            long id =  cell.getId();
            int rank = m_globalDisciplineId2NeutralMeshFileRank[id];
            if(rank != m_rank) {
                m_disciplineId2NeutralMeshFileCellPerRanks[id] = rank;
            }
        }
    }
    std::cout << "Rank " << m_rank << " m_disciplineId2NeutralMeshFileCellPerRanks size " << m_disciplineId2NeutralMeshFileCellPerRanks.size() << std::endl;


}

/*!
    Compute a local map ids->ranks for the neutral mesh in order to partition the neutral mesh
    with rank sub-domains overlapping the discipline sub-domains with the same rank.
*/

void MeshCoupling::computeNeutralId2DisciplineCellPerRanks() {

    for(const Cell & cell : m_unitNeutralMesh->getCells()) {
        if(cell.isInterior()) {
            long id = cell.getId();
            int rank = m_globalNeutralId2NeutralMeshFilePartitionedDisciplineRank[id];
            if(rank != m_rank) {
                m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks[id] = rank;
            }
        }
    }
}

/*!
    Dynamically partition neutral mesh by discipline partitioning
*/
void MeshCoupling::dynamicPartitionNeutralMeshByDiscipline() {

//    std::stringstream dumpFileStringStream;
//    dumpFileStringStream << "neutralMeshFilePartitioned_" << m_rank << ".dat";
//    std::ofstream dumpStream;
//    dumpStream.open(dumpFileStringStream.str().c_str());
//    m_unitNeutralMesh->dump(dumpStream);
//    dumpStream.close();
//    std::stringstream mapStringStream;
//    mapStringStream << "map_" << m_rank << ".dat";
//    std::ofstream mapStream;
//    mapStream.open(mapStringStream.str().c_str());
//
//    for(auto & elem : m_neutralFile2DisciplineCellPerRanks) {
//        mapStream << elem.first << " " << elem.second << std::endl;
//    }
//
//    mapStream.close();

//    for(int r = 0; r < m_nprocs; ++r) {
//        if(m_rank == r) {
//            for(const Cell & cell : m_unitNeutralMesh->getCells()) {
//                if(cell.isInterior()) {
//                    std::cout << "before interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                } else {
//                    std::cout << "before ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_unitNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                }
//                std::cout << std::flush;
//            }
//
//        }
//        MPI_Barrier(m_comm);
//        std::cout << std::flush;
//    }


    std::vector<adaption::Info> partitionInfo = m_unitNeutralMesh->partition(m_neutralFile2DisciplineCellPerRanks,true,true);
    //DEBUG

//    for(int r = 0; r < m_nprocs; ++r) {
//        if(m_rank == r) {
//            for(const Cell & cell : m_unitNeutralMesh->getCells()) {
//                if(cell.isInterior()) {
//                    std::cout << "after interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                } else {
//                    std::cout << "after ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_unitNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                }
//                std::cout << std::flush;
//            }
//
//        }
//        MPI_Barrier(m_comm);
//        std::cout << std::flush;
//    }

}

/*!
    Dynamically partition discipline mesh by neutral mesh file partitioning
*/
void MeshCoupling::staticPartitionDisciplineMeshByNeutralFile() {

    m_unitDisciplineMesh->setCommunicator(m_comm);
    std::vector<adaption::Info> partitionInfo = m_unitDisciplineMesh->partition(m_disciplineId2NeutralMeshFileCellPerRanks,true,true);

}

/*!
    Set scaled meshes
*/
void MeshCoupling::buildScaledMeshes() {

    std::cout << "Discipline Radius: " << m_disciplineRadius << std::endl;
    std::cout << "Neutral Radius: " << m_neutralRadius << std::endl;
    //discipline
    m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
//    std::array<double,3> origin = {{0.0,0.0,0.0}};
//    std::array<double,3> scaling = {{m_radius,m_radius,m_radius}};
//    m_scaledDisciplineMesh->scale(scaling);
//    PiercedVector<Vertex> & disciplineVertices = m_scaledDisciplineMesh->getVertices();
//    for(Vertex &v : disciplineVertices) {
//        v.scale(scaling,origin);
//    }
    scaleMeshToRadius(m_scaledDisciplineMesh,m_oldDisciplineRadius,m_disciplineRadius);
    m_disciplineNumberingInfo = PatchNumberingInfo(m_scaledDisciplineMesh.get());

    //neutral
    m_scaledNeutralMesh = PatchKernel::clone(m_unitNeutralMesh.get());
    //m_scaledNeutralMesh->scale(m_radius);
//    PiercedVector<Vertex> & neutralVertices = m_scaledNeutralMesh->getVertices();
//    for(Vertex & v : neutralVertices) {
//        v.scale(scaling,origin);
//    }
    scaleMeshToRadius(m_scaledNeutralMesh, m_oldNeutralRadius, m_neutralRadius);
    m_neutralNumberingInfo = PatchNumberingInfo(m_scaledNeutralMesh.get());
}

/*!
    Dynamically partition neutral mesh by discipline partitioning
*/
void MeshCoupling::dynamicPartitionNeutralMeshByNeutralMeshFilePartitionedDiscipline() {

//    std::stringstream dumpFileStringStream;
//    dumpFileStringStream << "neutralMeshFilePartitioned_" << m_rank << ".dat";
//    std::ofstream dumpStream;
//    dumpStream.open(dumpFileStringStream.str().c_str());
//    m_scaledNeutralMesh->dump(dumpStream);
//    dumpStream.close();
//    std::stringstream mapStringStream;
//    mapStringStream << "map_" << m_rank << ".dat";
//    std::ofstream mapStream;
//    mapStream.open(mapStringStream.str().c_str());
//
//    for(auto & elem : m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks) {
//        mapStream << elem.first << " " << elem.second << std::endl;
//    }
//
//    mapStream.close();
//
//    std::string name = "scaledMeshBefore";
//    m_scaledNeutralMesh->write(name);

//    for(int r = 0; r < m_nprocs; ++r) {
//        if(m_rank == r) {
//            for(const Cell & cell : m_scaledNeutralMesh->getCells()) {
//                if(cell.isInterior()) {
//                    std::cout << "before interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_scaledNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                } else {
//                    std::cout << "before ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_scaledNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_scaledNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
//                }
//                std::cout << std::flush;
//            }
//
//        }
//        MPI_Barrier(m_comm);
//        std::cout << std::flush;
//    }
//
//    for(int r = 0; r < m_nprocs; ++r) {
//        if(m_rank == r) {
//            for(auto & elem : m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks) {
//                std::cout << "map rank " << m_rank << " " << elem.first << " " << elem.second << std::endl;
//            }
//        }
//        MPI_Barrier(m_comm);
//        std::cout << std::flush;
//    }
//    std::cout << std::flush;
//    MPI_Barrier(m_comm);

    std::cout << "Partition prepare.." << std::endl;
    std::vector<adaption::Info> partitionInfo = m_scaledNeutralMesh->partitioningPrepare(m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks,true);
    //std::vector<adaption::Info> partitionInfo = m_scaledNeutralMesh->partition(m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks,true,false);
    //m_scaledNeutralMesh->getCells().squeeze();
    //m_scaledNeutralMesh->squeeze();

    //Communicate exchanged cell values during partitioning
    size_t singleCellByteSize = sizeof(long) + m_neutralData.getFieldCount() * sizeof(double) * 2;
    for(auto info : partitionInfo) {

        int sendRank = info.rank;
        size_t bufferSize = info.previous.size() * singleCellByteSize + sizeof(size_t);
        m_lbCommunicator->setSend(sendRank,bufferSize);

        SendBuffer &sendBuffer = m_lbCommunicator->getSendBuffer(sendRank);
        sendBuffer << info.previous.size();
        for(const long & sendId : info.previous) {
            sendBuffer << sendId;
            double *data = m_neutralData.data(sendId,0);
            for(size_t k = 0; k < m_neutralData.getFieldCount(); ++k) {
                sendBuffer << data[k];
            }
        }

    }
    m_lbCommunicator->discoverRecvs();
    m_lbCommunicator->startAllRecvs();
    m_lbCommunicator->startAllSends();

    partitionInfo = m_scaledNeutralMesh->partitioningAlter(true,false);
    m_scaledNeutralMesh->partitioningCleanup();

    long recvId;
    size_t nofCells;
    const std::vector<int> & recvRanks = m_lbCommunicator->getRecvRanks();
    for(int rank : recvRanks) {
        m_lbCommunicator->waitRecv(rank);
        RecvBuffer & recvBuffer = m_lbCommunicator->getRecvBuffer(rank);
        recvBuffer >> nofCells;
        for(size_t cell = 0; cell < nofCells; ++cell) {
            recvBuffer >> recvId;
            double *data = m_neutralData.data(recvId,0);
            for(size_t k = 0; k < m_neutralData.getFieldCount(); ++k) {
                recvBuffer >> data[k];
            }
        }
    }

}


/*!
    Dynamically partition neutral mesh by discipline partitioning
*/
void MeshCoupling::dynamicPartitionNeutralMeshByNeutralMeshWithData() {

    //TODO continue from here
    for(Cell & cell : m_scaledNeutralMesh->getCells()) {
        if(cell.isInterior()) {
            long id = cell.getId();
            m_neutralMeshFilePartitionCellPerRank[id] = m_globalNeutralId2MeshFileRank[id];
        }
    }

    std::cout << "Partition prepare.." << std::endl;
    std::vector<adaption::Info> partitionInfo = m_scaledNeutralMesh->partitioningPrepare(m_neutralMeshFilePartitionCellPerRank,true);

    //Communicate exchanged cell values during partitioning
    size_t singleCellByteSize = sizeof(long) + m_neutralData.getFieldCount() * sizeof(double) * 2;
    for(auto info : partitionInfo) {

        int sendRank = info.rank;
        size_t bufferSize = info.previous.size() * singleCellByteSize + sizeof(size_t);
        m_lbCommunicator->setSend(sendRank,bufferSize);

        SendBuffer &sendBuffer = m_lbCommunicator->getSendBuffer(sendRank);
        sendBuffer << info.previous.size();
        for(const long & sendId : info.previous) {
            sendBuffer << sendId;
            double *data = m_neutralData.data(sendId,0);
            for(size_t k = 0; k < m_neutralData.getFieldCount(); ++k) {
                sendBuffer << data[k];
            }
        }

    }
    m_lbCommunicator->discoverRecvs();
    m_lbCommunicator->startAllRecvs();
    m_lbCommunicator->startAllSends();

    partitionInfo = m_scaledNeutralMesh->partitioningAlter(true,true);
    m_scaledNeutralMesh->partitioningCleanup();

    long recvId;
    size_t nofCells;
    const std::vector<int> & recvRanks = m_lbCommunicator->getRecvRanks();
    for(int rank : recvRanks) {
        m_lbCommunicator->waitRecv(rank);
        RecvBuffer & recvBuffer = m_lbCommunicator->getRecvBuffer(rank);
        recvBuffer >> nofCells;
        for(size_t cell = 0; cell < nofCells; ++cell) {
            recvBuffer >> recvId;
            double *data = m_neutralData.data(recvId,0);
            for(size_t k = 0; k < m_neutralData.getFieldCount(); ++k) {
                recvBuffer >> data[k];
            }
        }
    }

}

/*!
    Synchronize each data container with relative mesh
*/
void MeshCoupling::synchronizeMeshData() {

    m_disciplineData.setDynamicKernel(&(m_scaledDisciplineMesh->getCells()),PiercedVector<Cell>::SYNC_MODE_JOURNALED);
    m_neutralData.setDynamicKernel(&(m_scaledNeutralMesh->getCells()),PiercedVector<Cell>::SYNC_MODE_JOURNALED);

}

/*!
    Set all data container to the same uniform value
*/
void MeshCoupling::uniformlyInitAllData(double value) {

    m_neutralData.fill(value);
    m_disciplineData.fill(value);

//    //DEBUG
//    double neutralValues[2];
////    neutralValues[fid_temperature] = double(m_rank); neutralValues[fid_flux] = double(10 + m_rank);
//    neutralValues[m_fid_temperatureIN] = 274.0; neutralValues[m_fid_temperatureOUT] = 274.0;
//    double disciplineValues[2];
//    disciplineValues[m_fid_temperatureIN] = 274.0; disciplineValues[m_fid_temperatureOUT] = 274.0;
//
//    for(const auto & cell : m_scaledNeutralMesh->getCells()) {
//        long id = cell.getId();
//        m_neutralData.set(id,2,0,neutralValues);
//    }
//    for(const auto & cell : m_scaledDisciplineMesh->getCells()) {
//        long id = cell.getId();
//        m_disciplineData.set(id,2,0,disciplineValues);
////        m_disciplineData.set(id,fid_temperature,double(m_rank+40));
////        m_disciplineData.set(id,fid_flux,double(m_rank+30));
//    }
//
//    //DEBUG
}

/*!
    Set VTK writer (name, counter and field to be written)
*/
void MeshCoupling::prepareWritingData() {

    std::vector<std::string> fieldDataNames(m_inputDataNames);
    for(std::string & name : fieldDataNames) {
        name = name + "_input";
    }
    fieldDataNames.insert(fieldDataNames.end(),m_outputDataNames.begin(),m_outputDataNames.end());
    for(size_t i = m_inputDataNames.size(); i < fieldDataNames.size(); ++i) {
        fieldDataNames[i] = fieldDataNames[i] + "_output";
    }
//    for(std::string & name : fieldDataNames) {
//        std::cout << name << " ";
//    }
//    std::cout << std::endl;

    m_neutralVTKFieldStreamer = std::unique_ptr<coupling::FieldStreamer>(new coupling::FieldStreamer(*(m_scaledNeutralMesh.get()),m_neutralData,fieldDataNames));
    m_disciplineVTKFieldStreamer = std::unique_ptr<coupling::FieldStreamer>(new coupling::FieldStreamer(*(m_scaledDisciplineMesh.get()),m_disciplineData,fieldDataNames));
    m_scaledNeutralMesh->getVTK().setName(m_name + "_neutralMesh");
    m_scaledNeutralMesh->getVTK().setCounter();
    m_scaledNeutralMesh->setVTKWriteTarget(bitpit::PatchKernel::WRITE_TARGET_CELLS_INTERNAL);
    for(const std::string & fieldName : m_neutralVTKFieldStreamer->getFieldNames()) {
        m_scaledNeutralMesh->getVTK().addData<double>(fieldName, VTKFieldType::SCALAR, VTKLocation::CELL, m_neutralVTKFieldStreamer.get());
    }
    m_scaledDisciplineMesh->getVTK().setName(m_name + "_disciplineMesh");
    m_scaledDisciplineMesh->getVTK().setCounter();
    m_scaledNeutralMesh->setVTKWriteTarget(bitpit::PatchKernel::WRITE_TARGET_CELLS_INTERNAL);
    for(const std::string & fieldName : m_disciplineVTKFieldStreamer->getFieldNames()) {
        m_scaledDisciplineMesh->getVTK().addData<double>(fieldName, VTKFieldType::SCALAR, VTKLocation::CELL, m_disciplineVTKFieldStreamer.get());
    }

}

/*!
    Creates a new ghost communicator.

    \param continuous defines if the communicator will be set in continuous mode
*/
void MeshCoupling::createGhostCommunicators(bool continuous) {
    // Create communicators
    GhostCommunicator *neutralGhostCommunicator = new GhostCommunicator(m_scaledNeutralMesh.get());
    GhostCommunicator *disciplineGhostCommunicator = new GhostCommunicator(m_scaledDisciplineMesh.get());
    neutralGhostCommunicator->resetExchangeLists();
    disciplineGhostCommunicator->resetExchangeLists();
    neutralGhostCommunicator->setRecvsContinuous(continuous);
    disciplineGhostCommunicator->setRecvsContinuous(continuous);

    // Communicator tag
    m_neutralTag = neutralGhostCommunicator->getTag();
    m_disciplineTag = disciplineGhostCommunicator->getTag();

    // Add ghost communicator
    m_neutralGhostCommunicator = std::unique_ptr<GhostCommunicator>(neutralGhostCommunicator);
    m_disciplineGhostCommunicator = std::unique_ptr<GhostCommunicator>(disciplineGhostCommunicator);
}

/*!
    Initialize the data communicators that will be used for exchanging ghost data.
*/
void MeshCoupling::initializeGhostCommunicators() {

    m_neutralGhostStreamer = std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>,double>>(new ListBufferStreamer<bitpit::PiercedStorage<double, long>,double>(&m_neutralData,2*sizeof(double)));
    m_disciplineGhostStreamer = std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>,double>>(new ListBufferStreamer<bitpit::PiercedStorage<double, long>,double>(&m_disciplineData,2*sizeof(double)));

    createGhostCommunicators(true);

    m_neutralGhostCommunicator->addData(m_neutralGhostStreamer.get());
    m_disciplineGhostCommunicator->addData(m_disciplineGhostStreamer.get());
}

/*!
    Update neutral data ghost values
*/
void MeshCoupling::updateNeutralGhosts() {

    m_neutralGhostCommunicator->resetExchangeLists();

    if (m_scaledNeutralMesh->isPartitioned()) {
        m_neutralGhostCommunicator->startAllExchanges();
    }

    // Receive pressure
    if (m_scaledNeutralMesh->isPartitioned()) {
        m_neutralGhostCommunicator->completeAllExchanges();
    }

}

/*!
    Update discipline data ghost values
*/
void MeshCoupling::updateDisciplineGhosts() {

    if (m_scaledDisciplineMesh->isPartitioned()) {
        m_disciplineGhostCommunicator->startAllExchanges();
    }

    // Receive pressure
    if (m_scaledDisciplineMesh->isPartitioned()) {
        m_disciplineGhostCommunicator->completeAllExchanges();
    }

}

/*!
    Use inputs and radius to perform computation on discipline mesh
*/
void MeshCoupling::disciplineKernel() {

    //Update Matrix element values
    updateSimplifiedDiscreteHelmholtzSystem();
    //Update RHS
    updateSystemRHS();

    std::cout << "Solving ..." << std::endl;
    solveSytem();

    updateOutputField();
}


void MeshCoupling::disciplineKernel1() {
    float res = 0;
    std::size_t counter = 0;
    for(const Cell & cell : m_scaledDisciplineMesh->getCells()) {
        long id = cell.getId();
        if(cell.isInterior()) {
	    res = 1 - 0.4*m_disciplineData.at(id);
            m_disciplineData.set(id,res);
            ++counter;
        }
    }
    updateDisciplineGhosts();
}

void MeshCoupling::disciplineKernel2() {
    float res = 0;
    std::size_t counter = 0;
    for(const Cell & cell : m_scaledDisciplineMesh->getCells()) {
        long id = cell.getId();
        if(cell.isInterior()) {
	    res = m_disciplineData.at(id);
            m_disciplineData.set(id,res);
            ++counter;
        }
    }


    updateDisciplineGhosts();
}

/*!
    Extract row elements of the Output Input Jacobian Matrix to be passed to PETSc4Py
    \param[in] index it is the cell order number from GEMS. ATTENTION: IT IS MANDATORY THAT PARTITIONED MESHES ARE SQUEEZED!!!!
    \param[out] cellGlobalId a global consecutive id to be used in PETSc
    \param[out] columnIds a vector long containing the indices of non-zero elements in the matrix row relative to cellGlobalId
    \param[out] columnValues a vector double containing the values of non-zero elements in the matrix row relative to cellGlobalId
*/
void MeshCoupling::extractOutputInputJacobianRow(int index, int & cellGlobalId, std::vector<int> & columnIds, std::vector<double> & columnValues) {


    columnIds.clear();
    columnValues.clear();
    //Access cell id from raw iterator of PiercedVector using ordered local index from Python
    long cellId = m_scaledNeutralMesh->getCells().rawFind(index).getId();
    cellGlobalId = getNeutralGlobalConsecutiveId(cellId);

    int ncols;
    const int *cols;
    const double *values;

    MatGetRow(getJacobianManager()->getOuputInputJacobian(), int(cellGlobalId), &ncols, &cols, &values);

    columnIds.resize(ncols);
    columnValues.resize(ncols);
    std::copy(cols,cols + ncols,columnIds.begin());
    std::copy(values,values + ncols,columnValues.begin());

    PetscFree(cols);
    PetscFree(values);

}

/*!
    Extract row elements of the Output Control Jacobian Matrix to be passed to PETSc4Py
    \param[in] index it is the cell order number from GEMS. ATTENTION: IT IS MANDATORY THAT PARTITIONED MESHES ARE SQUEEZED!!!!
    \param[out] cellGlobalId a global consecutive id to be used in PETSc
    \param[out] columnIds a vector long containing the indices of non-zero elements in the matrix row relative to cellGlobalId
    \param[out] columnValues a vector double containing the values of non-zero elements in the matrix row relative to cellGlobalId
*/
void MeshCoupling::extractOutputControlJacobianRow(int index, int & cellGlobalId, std::vector<int> & columnIds, std::vector<double> & columnValues) {


    columnIds.clear();
    columnValues.clear();
    //Access cell id from raw iterator of PiercedVector using ordered local index from Python
    long cellId = m_scaledNeutralMesh->getCells().rawFind(index).getId();
    cellGlobalId = getNeutralGlobalConsecutiveId(cellId);

    int ncols;
    const int *cols;
    const double *values;

    MatGetRow(getJacobianManager()->getOuputControlJacobian(), int(cellGlobalId), &ncols, &cols, &values);

    columnIds.resize(ncols);
    columnValues.resize(ncols);
    std::copy(cols,cols + ncols,columnIds.begin());
    std::copy(values,values + ncols,columnValues.begin());

    PetscFree(cols);
    PetscFree(values);

}


/*!
    Get the information structure on neutral mesh elements numbering
    \return a PatchNumberingInfo for the scaled neutral mesh
*/
const PatchNumberingInfo & MeshCoupling::getNeutraNumberingInfo() {

    return m_neutralNumberingInfo;

}

/*!
    Get the global consecutive index of a neutral mesh cell from PatchNumberingInfo
    \param[in] id local index of the cell
    \return the global consecutive index for cell with local index id
*/
long MeshCoupling::getNeutralGlobalConsecutiveId(long id){
    return m_neutralNumberingInfo.getCellConsecutiveId(id);
}


/*!
    Assemble the linear system matrix for a simplified version of the discrete Helmoltz operator
*/
void MeshCoupling::assemblySimplifiedDiscreteHelmholtzSystem() {

    //Compute Stencils
    computeSimplifiedDiscreteLaplaceStencils(m_helmoltzStencils);
    computeHelmholtzStencilsFromLaplaceStencils(m_helmoltzStencils, m_emissivity);

    //Assembly matrix and linear system
    KSPOptions &solverOptions = m_system->getKSPOptions();
    solverOptions.rtol    = 1e-7;//set PETSc KSP global tolerance
    solverOptions.subrtol = 1e-7;//set PETSc blocks tolerance
    solverOptions.maxits  = 1000; //DEBUG TODO
    m_system->assembly(m_scaledDisciplineMesh->getCommunicator(),m_scaledDisciplineMesh->isPartitioned(),m_helmoltzStencils);

}

/*!
    Update the linear system matrix for a simplified version of the discrete Helmoltz operator
*/
void MeshCoupling::updateSimplifiedDiscreteHelmholtzSystem() {

    //Compute Stencils
    computeSimplifiedDiscreteLaplaceStencils(m_helmoltzStencils);
    computeHelmholtzStencilsFromLaplaceStencils(m_helmoltzStencils, m_emissivity);

    //update matrix
    //m_system->update(m_rowIds,m_helmoltzStencils);
    m_system->update(m_helmoltzStencils);
}

/*!
    Compute the linear system matrix for a simplified version of the discrete Laplace operator times thickness and thermal conductivity
    \param[out] laplaceStencils a vector collecting Laplace bitpit stencil to be passed to Sysyem Solver assembly method.
*/
void MeshCoupling::computeSimplifiedDiscreteLaplaceStencils(std::vector<StencilScalar> & laplaceStencils) {

    const std::unordered_map<long, long> &scaledDisciplineConsecutiveMapping = m_disciplineNumberingInfo.getCellConsecutiveMap();
    std::vector<long> neighs;
    std::vector<double> weights;
    neighs.reserve(20);
    weights.reserve(20);
    for (const Cell & cell : m_scaledDisciplineMesh->getCells()) {
        if(cell.isInterior()) {
            double weightsSum = 0.0;
            neighs.clear();
            weights.clear();

            long cellId = cell.getId();
            long cellConsecutiveId = m_disciplineNumberingInfo.getCellConsecutiveId(cellId);
            long cellLocalConsecutiveId = cellConsecutiveId - m_disciplineNumberingInfo.getCellGlobalCountOffset();
            std::array<double,3> cellCentroid = m_scaledDisciplineMesh->evalCellCentroid(cellId);

            StencilScalar & stencil = laplaceStencils[cellLocalConsecutiveId];
            stencil.initialize(1,0.0);
            stencil.zero();

            neighs = m_scaledDisciplineMesh->findCellNeighs(cellId);
            for(const long & neigh : neighs) {
                std::array<double,3> neighCentroid = m_scaledDisciplineMesh->evalCellCentroid(neigh);
                double weight = 1.0 / norm2(neighCentroid - cellCentroid);
                weightsSum += weight;
                weights.push_back(weight);
            }
            weights /= weightsSum;

            stencil.reserve(neighs.size()+1);
            for(size_t i = 0; i < neighs.size(); ++i) {
                stencil.sumItem(neighs[i],weights[i]);
            }
            //stencil.reserve(1);
            stencil.sumItem(cellId,-1.0);
//            if(m_scaledDisciplineMesh->getRank() == 0) {
//                std::cout << "local consecutive ID " << cellLocalConsecutiveId << " - " << cellConsecutiveId << " - " << cellId << std::endl;
//                stencil.display(std::cout);
//            }
            stencil.renumber(scaledDisciplineConsecutiveMapping);
//            if(m_scaledDisciplineMesh->getRank() == 0) {
//                std::cout << "local consecutive ID " << cellLocalConsecutiveId << " - " << cellConsecutiveId << " - " << cellId << std::endl;
//                stencil.display(std::cout);
//            }
            stencil = (-1.0) * stencil * m_thickness * evalThermalDiffusivity();
        }
    }

}

/*!
    Compute Helmholtz matrix from the Laplace one.
    \param[out] helmholtzStencils a vector collecting Helmholtz bitpit stencil to be passed to Sysyem Solver assembly method.
    \param[int] coefficient radiation flux coefficient
*/
void MeshCoupling::computeHelmholtzStencilsFromLaplaceStencils(std::vector<StencilScalar> & helmoltzStencils, const double & coefficient) {

    double emissivityMultiplier = 2.0;
    if(m_innerSphere) {
        emissivityMultiplier = 1.0;
    }
    for (const Cell & cell : m_scaledDisciplineMesh->getCells()) {
        if(cell.isInterior()) {
            long cellId = cell.getId();
            long cellConsecutiveId = m_disciplineNumberingInfo.getCellConsecutiveId(cellId);
            long cellLocalConsecutiveId = cellConsecutiveId - m_disciplineNumberingInfo.getCellGlobalCountOffset();
            StencilScalar & stencil = helmoltzStencils[cellLocalConsecutiveId];
            stencil.sumItem(cellConsecutiveId,emissivityMultiplier * coefficient);//global consecutive is mandatory because the stencil has been renumbered
        }
    }
}

/*!
    Update linear system rhs by inserting external constant source flux contribution
    \param[out] cellNormal the surface normal at cell location
*/
void MeshCoupling::updateSystemRHS() {

    long nLocalRow = m_system->getRowCount();
    double *rhs = m_system->getRHSRawPtr();
    if(m_innerSphere) {
        for(const Cell & cell : m_scaledDisciplineMesh->getCells()) {
            if(cell.isInterior()) {
                long cellId = cell.getId();
                long cellConsecutiveId = m_disciplineNumberingInfo.getCellConsecutiveId(cellId);
                long cellLocalConsecutiveId = cellConsecutiveId - m_disciplineNumberingInfo.getCellGlobalCountOffset();
                assert(cellLocalConsecutiveId < nLocalRow);

                rhs[cellLocalConsecutiveId] = 0.0;
                rhs[cellLocalConsecutiveId] = m_emissivity * m_disciplineData.at(cellId,m_inputField); // emissivity * Temp_from_outer_discipline
            }
        }
    } else {
        for(const Cell & cell : m_scaledDisciplineMesh->getCells()) {
            if(cell.isInterior()) {
                long cellId = cell.getId();
                long cellConsecutiveId = m_disciplineNumberingInfo.getCellConsecutiveId(cellId);
                long cellLocalConsecutiveId = cellConsecutiveId - m_disciplineNumberingInfo.getCellGlobalCountOffset();
                assert(cellLocalConsecutiveId < nLocalRow);

                std::array<double,3> cellNormal = m_scaledDisciplineMesh->evalFacetNormal(cellId);
                rhs[cellLocalConsecutiveId] = evalSourceIntensity(cellNormal) + m_emissivity * m_disciplineData.at(cellId,m_inputField)
                        + m_emissivity * m_infinityTemperature; // source + emissivity * Temp_from_inner_discipline + emissivity * Temp_from_infinity
            }
        }
    }
    m_system->restoreRHSRawPtr(rhs);
}

/*!
    Evaluate diffusivity as function of radius
*/
double MeshCoupling::evalThermalDiffusivity() {

    double thermalDiffusivity = m_thermalDiffusivityCoefficient * m_disciplineRadius;//coeff*m_radius for one discipline, coeff*m_radius^2 for the other one

    if(!m_innerSphere) {
        thermalDiffusivity *= m_disciplineRadius;
    }

    return thermalDiffusivity;
}

/*!
    Evaluate radiative external source intensity as function of the cell orientation
    \param[out] cellNormal local cell normal to the surface
*/
double MeshCoupling::evalSourceIntensity(const std::array<double,3> & cellNormal) {

    double intensity = 0.0;

    double visibility = dotProduct(cellNormal,m_sourceDirection);
    if(visibility >= 0.0) {
        return 0.0;
    }

    intensity = fabs(visibility) * m_sourceMaxIntensity;

    return intensity;
}

/*!
    Solve linear system
*/
void MeshCoupling::solveSytem() {

    m_system->solve();
    std::cout << "Iterations = " << m_system->getKSPStatus().its << std::endl;

}

/*!
    Update output field in discipline pierced storage with linear system solution
*/
void MeshCoupling::updateOutputField() {

    long nLocalRow = m_system->getRowCount();
    const double *solution = m_system->getSolutionRawReadPtr();
    for(const Cell & cell : m_scaledDisciplineMesh->getCells()) {
        if(cell.isInterior()) {
            long cellId = cell.getId();
            long cellConsecutiveId = m_disciplineNumberingInfo.getCellConsecutiveId(cellId);
            long cellLocalConsecutiveId = cellConsecutiveId - m_disciplineNumberingInfo.getCellGlobalCountOffset();
            assert(cellLocalConsecutiveId < nLocalRow);

            m_disciplineData.set(cellId,m_outputField,solution[cellLocalConsecutiveId]);
        }
    }
    m_system->restoreSolutionRawReadPtr(solution);

}

/*!
    Update input field neutral pierced storage with data from the other discipline
*/
void MeshCoupling::updateInputField(double *neutralInputArray, std::size_t size) {

    //put Input data into neutral PiercedStorage - neutralInputArray should have the same number of elements of neutral mesh rank sub-domain(internals)
    std::size_t counter = 0;
    assert(size == m_scaledNeutralMesh->getInternalCount());
    for(const Cell & cell : m_scaledNeutralMesh->getCells()) {
        long id = cell.getId();
        if(cell.isInterior()) {
            m_neutralData.set(id,m_inputField,neutralInputArray[counter]);
            ++counter;
        }
    }

}

/*!
    Scale mesh to current radius. It updates the bounding box
*/
void MeshCoupling::scaleMeshToRadius(std::unique_ptr<SurfUnstructured> & mesh, double & oldRadius, const double & newRadius) {

    std::array<double,3> origin = {{0.0,0.0,0.0}};
    double scaleFactor = newRadius / oldRadius;
    std::array<double,3> scaling = {{scaleFactor, scaleFactor, scaleFactor}};
    PiercedVector<Vertex> & vertices = mesh->getVertices();
    for(Vertex & v : vertices) {
        v.scale(scaling,origin);
    }
    mesh->updateBoundingBox(true);
    oldRadius = newRadius;
}

/*!
    Get a pointer to the Jacobian Manager
*/
JacobianMatricesManager* MeshCoupling::getJacobianManager() {
    return m_jacobianManager.get();
}

}

