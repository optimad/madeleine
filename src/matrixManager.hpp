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
    Mat m_EllipticOperatorJacobian;
    Mat m_EllipticOperatorInverse;
    Mat m_Jacobian;
    MPI_Comm m_comm;
    StencilScalarSolverHandler & m_ellipticLinearSystem;

public:
    JacobianMatricesManager(std::string name, MPI_Comm comm,
            StencilScalarSolverHandler & ellipticLinearSystem);
    ~JacobianMatricesManager();
    void computeEllipticOperatorInverse(SurfUnstructured* mesh, PatchNumberingInfo & numberingInfo);
    void computeInterpolationMatrices(SurfUnstructured * disciplineMesh, PatchNumberingInfo * disciplineNumberingInfo,
            SurfUnstructured * neutralMesh, PatchNumberingInfo * neutralNumberingInfo);
    Mat & getNeutralToDisciplineInterpolatorJacobian();
    Mat & getDisciplineToNeutralInterpolatorJacobian();
    Mat & getEllipticOperatorJacobian();
};

#endif
