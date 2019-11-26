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

    //DEBUG
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
    //DEBUG
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

    //Partition the neutral mesh by mesh file order
    staticPartitionNeutralMeshByMeshFileOrder();

    //Compute the discipline ids/ranks map to get discipline rank sub-domains overlapping the neutral ones with the same rank
    computeGlobalDisciplineId2NeutralRank();




//    //Compute the neutral ids/ranks map to get neutral rank sub-domains overlapping the discipline ones with the same rank
//    computeGlobalNeutralId2DisciplineRank();
//
//
//    //TODO prepare cellRanks per discipline
//
    for(const Cell & cell : m_unitDisciplineMesh->getCells()) {
        if(cell.isInterior()) {
            long id =  cell.getId();
            int rank = m_globalDisciplineId2NeutralRank[id];
            if(rank != m_rank) {
                m_discipline2FileNeutralCellPerRanks[id] = rank;
            }
        }
    }
    std::cout << "Rank " << m_rank << " m_discipline2FileNeutralCellPerRanks size " << m_discipline2FileNeutralCellPerRanks.size() << std::endl;
//
    //DEBUG

    for(int r = 0; r < m_nprocs; ++r) {
        if(m_rank == r) {
            for(auto idRank : m_discipline2FileNeutralCellPerRanks) {
                long id = idRank.first;
                int rank = idRank.second;
                std::cout << "map Rank " << m_rank << " id " << id << " rank " << rank << std::endl;
                std::cout << std::flush;
            }
//            for(const Cell & cell : m_unitNeutralMesh->getCells()) {
//                if(cell.isInterior()) {
//                    std::cout << "interior Rank " << m_rank << " id " << cell.getId() << std::endl;
//                } else {
//                    std::cout << "ghost Rank " << m_rank << " id " << cell.getId() << std::endl;
//                }
//                std::cout << std::flush;
//            }

        }
        MPI_Barrier(m_comm);
        std::cout << std::flush;
    }
//
//    //DEBUG
//
//    //TODO partiziona neutrale come disciplina
//    dynamicPartitionNeutralMeshByDiscipline();
//
//
//    //TODO prepara cellRanks per neutrale da disciplina a file
//    //TODO partiziona neutrale come file
//
//    //TODO scala griglie
//
//    //Inizializza dati in parallelo


//    //DEBUG
//    for(int r = 0; r < m_nprocs; ++r) {
//        if(m_rank == r) {
//            for(int i = 0; i < nofCellsNeutral; ++i) {
//                std::cout << "My rank: " << m_rank << " Id: " << i << " cell center: "
//                        << cellCentersNeutral[i * 3] << " " << cellCentersNeutral[i * 3 + 1] << " " << cellCentersNeutral[i * 3 + 2] << " "
//                        << " on discpline rank partition: " << m_globalNeutralId2DisciplineRank[i] << std::endl;
//            }
//            std::cout << std::flush;
//        }
//        MPI_Barrier(m_comm);
//        std::cout << std::flush;
//    }
//    //DEBUG


    //m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
#else
    m_unitDisciplineMesh->reset();
    m_disciplineData.clear(true);
    m_unitDisciplineMesh->importSTL(unitDisciplineMeshFile);
    m_unitDisciplineMesh->deleteCoincidentVertices();
    m_scaledDisciplineMesh = PatchKernel::clone(m_unitDisciplineMesh.get());
#endif






//
//    PiercedVector<Vertex> & scaledDisciplineVertices = m_scaledDisciplineMesh->getVertices();
//
//
//    //initialize discipline data structure and resize the unit sphere, scaling its vertices
//    for(Vertex & v : scaledDisciplineVertices) {
//        v.setCoords(v.getCoords()*radius);
//        m_disciplineData.insert(v.getId(),0.0);
//    }
//
//    //initialize neutral mesh and vertices container
//    m_unitNeutralMesh->reset();
//    m_unitNeutralMesh->importSTL(unitNeutralMeshFile);
//    m_unitNeutralMesh->deleteCoincidentVertices();
//    m_scaledNeutralMesh = PatchKernel::clone(m_unitNeutralMesh.get());
//    PiercedVector<Vertex> & scaledNeutralVertices = m_scaledNeutralMesh->getVertices();
//
//    //resize the neutral unit sphere mesh, scaling its vertices
//    for(Vertex & v : scaledNeutralVertices) {
//        v.setCoords(v.getCoords()*radius);
//    }
//
//    //print at log output vertices and cells of the discipline and neutral meshes
//    log::cout() << " Discipline Mesh n vertices : " << m_unitDisciplineMesh->getVertexCount() << std::endl;
//    log::cout() << " Discipline Mesh n cells : " << m_unitDisciplineMesh->getInternalCount() << std::endl;
//    log::cout() << " Neutral Mesh n vertices : " << m_unitNeutralMesh->getVertexCount() << std::endl;
//    log::cout() << " Neutral Mesh n cells : " << m_unitNeutralMesh->getInternalCount() << std::endl;
//
//    //print on VTK files both discipline and neutral meshes
//    std::string disciplineMeshVTKFileName = "disciplineMesh";
//    std::string neutralMeshVTKFileName = "neutralMesh";
//    writeMesh(m_scaledDisciplineMesh.get(),disciplineMeshVTKFileName);
//    writeMesh(m_scaledNeutralMesh.get(),neutralMeshVTKFileName);

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
    Get scaled neutral mesh

    \return the scaled neutral mesh
*/
size_t MeshCoupling::getNeutralMeshSize(){

    return m_scaledNeutralMesh.get()->getVertices().size();

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
    std::vector<adaption::Info> partitionInfo = mesh.partition(cellRanks,true,false);
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
    std::string name("metisPatitionedDiscipline");
    m_unitDisciplineMesh->write(name);

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
    std::vector<adaption::Info> partitionInfo = m_unitNeutralMesh->partition(staticNeutralMeshFilePartitionCellPerRank,true,false);
    //DEBUG
    std::string name("filePartitionedNeutral");
    m_unitNeutralMesh->write(name);
    //DEBUG

}

/*!
    Compute a global map ids->ranks for the neutral mesh in order to partition the neutral mesh
    with rank sub-domains overlapping the discipline sub-domains with the same rank.
*/
void MeshCoupling::computeGlobalNeutralId2DisciplineRank() {
    //Map the discipline rank distribution on neutral mesh
    //Rank 0 communicate neutral mesh ids and cell centers to the other ranks
    long *idxNeutral = nullptr; //ordered by neutral mesh id
    double * cellCentersNeutral = nullptr; //ordered by neutral mesh id
    int nofCellsNeutral = 0;
    if(m_rank == 0) {
        nofCellsNeutral = m_unitNeutralMesh->getCellCount();
        idxNeutral = new long[nofCellsNeutral];
        cellCentersNeutral = new double[nofCellsNeutral * 3];
        for(auto cell : m_unitNeutralMesh->getCells()) {
            long id = cell.getId();
            const std::array<double,3> & cellCenter = m_unitNeutralMesh->evalCellCentroid(id);
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

    //    //DEBUG
    //    for(int r = 0; r < m_nprocs; ++r) {
    //        if(r==m_rank) {
    //            for(int i = 0; i < nofCells; ++i) {
    //                std::cout << "id: " << idx[i] << " - cell center: ";
    //                for(int j = 0; j < 3; ++j) {
    //                    std::cout << cellCenters[i * 3 + j] << " ";
    //                }
    //                std::cout << std::endl;
    //            }
    //        }
    //        MPI_Barrier(m_comm);
    //    }
    //    //DEBUG

    //Project rank from discipline partition to cellRanks for neutral mesh
    //Build discipline mesh surface skd tree
    SurfaceSkdTree disciplineTree(m_unitDisciplineMesh.get());
    log::cout() << "Tree declared" << std::endl;
    disciplineTree.build(1);
    //Initialiaze rank maps - array ordered by id and containing ranks
    m_globalNeutralId2DisciplineRank.resize(nofCellsNeutral,-1);
    //Loop over neutral cell center and check the discipline cell the neutral cell center belong to. Assign rank according to
    std::array<double,3> tempCellCenter;
    darray3 lambda;
    darray3 xP;
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
                m_globalNeutralId2DisciplineRank[i] = m_rank;
            }
        }
    }
    //Reduce disciplineMappedRankForNeutral on Rank 0 for static partitioning
    MPI_Allreduce(MPI_IN_PLACE,m_globalNeutralId2DisciplineRank.data(),nofCellsNeutral,MPI_INT,MPI_MAX,m_comm);

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

    //DEBUG
//    for(int r = 0; r < m_nprocs; ++r) {
//        if(r==m_rank) {
//            for(int i = 0; i < nofCellsDiscipline; ++i) {
//                std::cout << "id: " << idxDiscipline[i] << " - cell center: ";
//                for(int j = 0; j < 3; ++j) {
//                    std::cout << cellCentersDiscipline[i * 3 + j] << " ";
//                }
//                std::cout << std::endl;
//            }
//        }
//        MPI_Barrier(m_comm);
//    }
    //DEBUG

    //Project rank from neutral partition to cellRanks for discipline mesh
    //Build neutral mesh surface skd tree
    SurfaceSkdTree neutralTree(m_unitNeutralMesh.get());
    log::cout() << "Tree declared" << std::endl;
    neutralTree.build(1);
    //Initialiaze rank maps - array ordered by id and containing ranks
    m_globalDisciplineId2NeutralRank.resize(nofCellsDiscipline,-1);
    std::vector<int> candidates(nofCellsDiscipline,-1);
    //Loop over discipline cell center and check the neutral cell the discipline cell center belong to. Assign rank according to
    std::array<double,3> tempCellCenter;
    darray3 lambda;
    darray3 xP;
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
        if(i == 10) {
            std::cout << "Rank " << m_rank << " cc " << std::setprecision(20)<< tempCellCenter << " closest cell dist " << std::setprecision(20) << closestCellDist<< std::endl;
        }
        if(m_unitNeutralMesh->getCell(id).isInterior()) {
            closestCellPlaneDist = bitpit::CGElem::distancePointPlane(tempCellCenter,
                    m_unitNeutralMesh->getCellVertexCoordinates(id)[0],m_unitNeutralMesh->evalFacetNormal(id),xP);
            if(checkFloating(closestCellDist,closestCellPlaneDist)) {
                m_globalDisciplineId2NeutralRank[i] = m_rank;
            } else {
                candidates[i] = m_rank;
            }
        }
    }
    //Reduce disciplineMappedRankForNeutral on Rank 0 for static partitioning
    MPI_Allreduce(MPI_IN_PLACE,m_globalDisciplineId2NeutralRank.data(),nofCellsDiscipline,MPI_INT,MPI_MAX,m_comm);
    MPI_Allreduce(MPI_IN_PLACE,candidates.data(),nofCellsDiscipline,MPI_INT,MPI_MAX,m_comm);
    //check cellcenters candidates
    for(int c = 0; c < nofCellsDiscipline; ++c) {
        if(m_globalDisciplineId2NeutralRank[c] == -1) {
            m_globalDisciplineId2NeutralRank[c] = candidates[c];
        }
    }


    //DEBUG
    for(int r = 0; r < m_nprocs; ++r) {
        if(r==m_rank) {
            for(int i = 0; i < m_globalDisciplineId2NeutralRank.size(); ++i) {
                std::cout << "id: " << i << " - rank: " << m_globalDisciplineId2NeutralRank[i];
                std::cout << std::endl;
            }
        }
        MPI_Barrier(m_comm);
    }
    //DEBUG

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

    for(int r = 0; r < m_nprocs; ++r) {
        if(m_rank == r) {
            for(const Cell & cell : m_unitNeutralMesh->getCells()) {
                if(cell.isInterior()) {
                    std::cout << "before interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                } else {
                    std::cout << "before ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_unitNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                }
                std::cout << std::flush;
            }

        }
        MPI_Barrier(m_comm);
        std::cout << std::flush;
    }


    std::vector<adaption::Info> partitionInfo = m_unitNeutralMesh->partition(m_neutralFile2DisciplineCellPerRanks,true,false);
    //DEBUG

    std::string name("disciplinePartitionedNeutral");

    for(int r = 0; r < m_nprocs; ++r) {
        if(m_rank == r) {
            for(const Cell & cell : m_unitNeutralMesh->getCells()) {
                if(cell.isInterior()) {
                    std::cout << "after interior Rank " << m_rank << " id " << cell.getId() << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                } else {
                    std::cout << "after ghost Rank " << m_rank << " id " << cell.getId() << " owner " << m_unitNeutralMesh->getCellRank(cell.getId()) << " cc = " << m_unitNeutralMesh->evalCellCentroid(cell.getId()) << std::endl;
                }
                std::cout << std::flush;
            }

        }
        MPI_Barrier(m_comm);
        std::cout << std::flush;
    }

    //m_unitNeutralMesh->write(name);
    //DEBUG

}
}
