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

using namespace bitpit;

namespace coupling{

/*!
    Default Constructor
*/
#if ENABLE_MPI==1
MeshCoupling::MeshCoupling(std::string disciplineName, MPI_Comm comm) :
        m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)), m_radius(1.0), m_name(disciplineName)
#else
MeshCoupling::MeshCoupling(std::string disciplineName, MPI_Comm comm) :
        m_unitDisciplineMesh(new SurfUnstructured(2,3)), m_unitNeutralMesh(new SurfUnstructured(2,3)), m_radius(1.0), m_name(disciplineName)
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

    char worldnameout[MPI_MAX_OBJECT_NAME];
    char workersnameout[MPI_MAX_OBJECT_NAME];
    int rlen;
    MPI_Comm_get_name( MPI_COMM_WORLD, worldnameout, &rlen );
    MPI_Comm_get_name( m_comm, workersnameout, &rlen );
    std::string workersName(workersnameout);
    std::string worldName(worldnameout);
    int worldSize;
    int worldRank;
    MPI_Comm_size(MPI_COMM_WORLD,&worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&worldRank);
    std::cout << "I'm " << m_rank << " of " << m_nprocs  << " on " << workersName  << " and " << worldRank << " of " << worldSize << " on " << worldName << std::endl;

    m_neutralTag = 0;
    m_disciplineTag = 0;

    m_lbCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(m_comm));

#else
    m_rank = 0;
    m_nprocs = 1;
#endif
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
};

/*!
    Initialize the unit sphere meshes, the radius and the radius scaled meshes. m_disciplineData is initialize to zero.

    \param[in] unitDisciplineMeshFile is the name of the file (.stl) containing the discipline mesh of the unit sphere
    \param[in] unitNeutralMeshFile is the name of the file (.stl) containing the neutral mesh of the unit sphere
    \param[in] radius is the value of the radius of the sphere discretized by the scaled meshes
*/
void MeshCoupling::initialize(const std::string & unitDisciplineMeshFile, const std::string & unitNeutralMeshFile, double radius, const std::vector<int> & globalNeutralId2MeshFileRank){

    //initialize radius
    m_radius = radius;

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
    uniformlyInitAllData(0.0);
    //Set mesh VTK writer
    //ATTENTION from here on scaled meshes are written with attached data, therefore data container have to be coherent with the relative mesh
    prepareWritingData();

    //Initialize ghost communicators
    initializeGhostCommunicators();

    std::string name = "D_Nf";
    m_scaledDisciplineMesh->write(name);

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

    \param[in] neutralInputArray contiguous C-array from NUMPY array ordered like id-ordered cells in Nf mesh file partitioning
*/
void MeshCoupling::compute(double *neutralInputArray, std::size_t size) {

    log::cout() << "First Interpolation" << std::endl;
    //sort cells by id - neutalInputArray should have values ordered like the neutral mesh file partitioning
    m_scaledNeutralMesh->sortCells();

    //put data into neutral PiercedStorage - neutralInputArray should have the same number of elements of neutral mesh rank sub-domain(internals)
    std::size_t counter = 0;
    for(const Cell & cell : m_scaledNeutralMesh->getCells()) {
        long id = cell.getId();
        if(cell.isInterior()) {
            m_neutralData.set(id,neutralInputArray[counter]);
            m_neutralData.set(id,m_scaledNeutralMesh->evalCellCentroid(id)[0]);//DEBUG
        } else {
            m_neutralData.set(id,-1);
        }
    }

    //Update neutral ghost cell values
    //It should not be the case but call m_neutralGhostCommunicator->resetExchangeLists() if neutral mesh has changed.

    updateNeutralGhosts();

    //Interpolate from N_f to D_{N_f]
    std::cout << "N_f to D_{N_f} interpolation." << std::endl;
    interpolateFromTo(m_scaledNeutralMesh.get(),&m_neutralData,m_scaledDisciplineMesh.get(),&m_disciplineData);

    //Update discipline ghost cell values
    updateDisciplineGhosts();

    disciplineKernel();

    std::string name = "Nf";
    m_scaledNeutralMesh->write(name);
    //With Data partitioning, but not useful, just for testing
    dynamicPartitionNeutralMeshByNeutralMeshFilePartitionedDiscipline();

    //Interpolate from N_f to D_{N_f]
    std::cout << "D_{N_f} to N_{D_{N_f}} interpolation." << std::endl;
    interpolateFromTo(m_scaledDisciplineMesh.get(),&m_disciplineData,m_scaledNeutralMesh.get(),&m_neutralData);

    dynamicPartitionNeutralMeshByNeutralMeshWithData();

    name = "unsortedNf";
    m_scaledNeutralMesh->write(name);

    m_scaledNeutralMesh->sortCells();

    name = "sortedNf";
    m_scaledNeutralMesh->write(name);

    for(int r = 0; r < m_nprocs; ++r) {
        if(m_rank == r) {
            for(const Cell & cell : m_scaledNeutralMesh->getCells()) {
                if(cell.isInterior()) {
                    std::cout << "before interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_scaledNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                } else {
                    std::cout << "before ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_scaledNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_scaledNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                }
                std::cout << std::flush;
            }

        }
        MPI_Barrier(m_comm);
        std::cout << std::flush;
    }


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
    Get scaled neutral mesh

    \return the scaled neutral mesh
*/
size_t MeshCoupling::getNeutralMeshSize(){

    return m_scaledNeutralMesh->getInternalCount();

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

    std::cout << "Radius: " << m_radius << std::endl;
    //discipline
    m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
    //m_scaledDisciplineMesh->scale(m_radius);
    PiercedVector<Vertex> & disciplineVertices = m_scaledDisciplineMesh->getVertices();
    std::array<double,3> origin = {{0.0,0.0,0.0}};
    std::array<double,3> scaling = {{m_radius,m_radius,m_radius}};
    for(Vertex &v : disciplineVertices) {
        v.scale(scaling,origin);
    }

    //neutral
    m_scaledNeutralMesh = PatchKernel::clone(m_unitNeutralMesh.get());
    PiercedVector<Vertex> & neutralVertices = m_scaledNeutralMesh->getVertices();
    for(Vertex & v : neutralVertices) {
        v.scale(scaling,origin);
    }

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
    for(int r = 0; r < m_nprocs; ++r) {
        if(m_rank == r) {
            for(auto & elem : m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks) {
                std::cout << "map rank " << m_rank << " " << elem.first << " " << elem.second << std::endl;
            }
        }
        MPI_Barrier(m_comm);
        std::cout << std::flush;
    }
    std::cout << std::flush;
    MPI_Barrier(m_comm);

    std::cout << "Partition prepare.." << std::endl;
    std::vector<adaption::Info> partitionInfo = m_scaledNeutralMesh->partitioningPrepare(m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks,true);
    //std::vector<adaption::Info> partitionInfo = m_scaledNeutralMesh->partition(m_neutralId2NeutralMeshFilePartitionedDisciplineCellPerRanks,true,false);
    //m_scaledNeutralMesh->getCells().squeeze();
    //m_scaledNeutralMesh->squeeze();

    //Communicate exchanged cell values during partitioning
    size_t singleCellByteSize = sizeof(long) + m_neutralData.getFieldCount() * sizeof(double);
    std::cout << "size single cell " << singleCellByteSize << std::endl;
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

    std::string name = "N_{D_Nf}";
    m_scaledNeutralMesh->write(name);


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
    size_t singleCellByteSize = sizeof(long) + m_neutralData.getFieldCount() * sizeof(double);
    std::cout << "size single cell " << singleCellByteSize << std::endl;
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

}

/*!
    Set VTK writer (name, counter and field to be written)
*/
void MeshCoupling::prepareWritingData() {

    m_neutralVTKFieldStreamer = std::unique_ptr<coupling::FieldStreamer>(new coupling::FieldStreamer(*(m_scaledNeutralMesh.get()),m_neutralData,m_inputDataNames));
    m_disciplineVTKFieldStreamer = std::unique_ptr<coupling::FieldStreamer>(new coupling::FieldStreamer(*(m_scaledDisciplineMesh.get()),m_disciplineData,m_outputDataNames));
    m_scaledNeutralMesh->getVTK().setName("neutralMesh");
    m_scaledNeutralMesh->getVTK().setCounter();
    m_scaledNeutralMesh->setVTKWriteTarget(bitpit::PatchKernel::WRITE_TARGET_CELLS_INTERNAL);
    for(const std::string & fieldName : m_neutralVTKFieldStreamer->getFieldNames()) {
        m_scaledNeutralMesh->getVTK().addData<double>(fieldName, VTKFieldType::SCALAR, VTKLocation::CELL, m_neutralVTKFieldStreamer.get());
    }
    m_scaledDisciplineMesh->getVTK().setName("disciplineMesh");
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

    m_neutralGhostStreamer = std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>>>(new ListBufferStreamer<bitpit::PiercedStorage<double, long>>(&m_neutralData));
    m_disciplineGhostStreamer = std::unique_ptr<ListBufferStreamer<bitpit::PiercedStorage<double, long>>>(new ListBufferStreamer<bitpit::PiercedStorage<double, long>>(&m_disciplineData));

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




    updateDisciplineGhosts();
}


}


