#ifndef __MATRIXMANAGER__
#define __MATRIXMANAGER__

#include <string>
#include <bitpit_discretization.hpp>
#include <bitpit_surfunstructured.hpp>
#include <bitpit_LA.hpp>
#include <petscmat.h>

using namespace bitpit;

class StencilScalarSolverHandler : public StencilScalarSolver{

public:
    StencilScalarSolverHandler(const std::string & prefix, bool debug = false);
    const Mat & getEllipticOperatorMatrix() const;
};

class JacobianMatricesManager {

    std::string m_name;
    Mat m_NeutralToDisciplineInterpolatorJacobian;
    Mat m_DisciplineToNeutralInterpolatorJacobian;
    Mat m_EllipticOperatorInverse;
    Mat m_OutputInputJacobian;
    Mat m_OutputControlJacobian;
    MPI_Comm m_comm;
    StencilScalarSolverHandler & m_ellipticLinearSystem;
    bool m_reuseJacobianMatrix;

public:
    JacobianMatricesManager(std::string name, MPI_Comm comm,
            StencilScalarSolverHandler & ellipticLinearSystem);
    ~JacobianMatricesManager();
    void computeEllipticOperatorInverse(SurfUnstructured* mesh, PatchNumberingInfo & numberingInfo);
    void computeInterpolationMatrices(SurfUnstructured * disciplineMesh, PatchNumberingInfo * disciplineNumberingInfo,
            SurfUnstructured * neutralMesh, PatchNumberingInfo * neutralNumberingInfo);
    void computeOutputInputJacobian(MPI_Comm comm, double emissivity, bool isInnerDiscipline);
    void computeOutputControlJacobian(MPI_Comm comm, double emissivity, bool isInnerDiscipline, double radius,
            const PiercedStorage<double,long> & disciplineData, const SurfUnstructured * disciplineMesh, const PatchNumberingInfo & disciplineNumberingInfo,
            const SurfUnstructured * neutralMesh, const PatchNumberingInfo & neutralNumberingInfo);
    Mat & getNeutralToDisciplineInterpolatorJacobian();
    Mat & getDisciplineToNeutralInterpolatorJacobian();
    Mat & getOuputInputJacobian();
    Mat & getOuputControlJacobian();
};

#endif
