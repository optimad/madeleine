/*
 * coupling.cpp
 *
 *  Created on: 14 jul 2020
 *      Author: marco
 */

#include "matrixManager.hpp"
#include "couplingUtils.hpp"

JacobianMatricesManager::JacobianMatricesManager(std::string name, MPI_Comm comm,
        StencilScalarSolverHandler & ellipticLinearSystem)
: m_name(name), m_comm(comm), m_ellipticLinearSystem(ellipticLinearSystem) {
}

JacobianMatricesManager::~JacobianMatricesManager() {
    MatDestroy(&m_EllipticOperatorJacobian);
    MatDestroy(&m_EllipticOperatorInverse);
    MatDestroy(&m_DisciplineToNeutralInterpolatorJacobian);
    MatDestroy(&m_NeutralToDisciplineInterpolatorJacobian);
}

void JacobianMatricesManager::computeEllipticOperatorInverse(SurfUnstructured* mesh, PatchNumberingInfo & numberingInfo) {
    Mat B;
    Mat A;

    PetscInt nofDisciplineLocalCells = mesh->getInternalCount();
    //Create Elliptic Operator inversion auxiliary dense B matrix
    MatCreateDense(m_comm,nofDisciplineLocalCells,nofDisciplineLocalCells,PETSC_DECIDE, PETSC_DECIDE,NULL, &B);
    //Fill B with the identity
    for (const Cell & cell : mesh->getCells()) {
        if(cell.isInterior()) {
            long cellId = cell.getId();
            long cellConsecutiveId = numberingInfo.getCellConsecutiveId(cellId);
            MatSetValue(B,cellConsecutiveId,cellConsecutiveId,1.0,INSERT_VALUES);
        }
    }
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    //Create the elliptic operator inverse matrix  by duplication of B dense matrix
    MatDuplicate(B,MAT_COPY_VALUES,&m_EllipticOperatorInverse);

    //Create LU factorization of elliptic problem matrix
    //MatDuplicate(m_ellipticLinearSystem.getEllipticOperatorMatrix(),MAT_COPY_VALUES,&A);
    MatGetFactor(m_ellipticLinearSystem.getEllipticOperatorMatrix(),MATSOLVERMUMPS,MAT_FACTOR_LU,&A);
    MatLUFactorSymbolic(A,m_ellipticLinearSystem.getEllipticOperatorMatrix(),NULL,NULL,NULL);
    MatLUFactorNumeric(A,m_ellipticLinearSystem.getEllipticOperatorMatrix(),NULL);

    //Compute inverse
    MatMatSolve(A,B,m_EllipticOperatorInverse);

//    //DEBUG
//    PetscViewer matViewer;
//    PetscViewerCreate(m_comm, &matViewer);
//    PetscViewerSetType(matViewer, PETSCVIEWERASCII);
//    PetscViewerFileSetMode(matViewer, FILE_MODE_WRITE);
//    PetscViewerPushFormat(matViewer, PETSC_VIEWER_ASCII_MATLAB);
//
//    std::stringstream filePathStream;
//    filePathStream.str(std::string());
//    filePathStream << "./inverse.txt";
//    PetscViewerFileSetName(matViewer, filePathStream.str().c_str());
//    MatView(m_EllipticOperatorInverse, matViewer);
//    PetscViewerDestroy(&matViewer);

    MatDestroy(&A);
    MatDestroy(&B);
}

void JacobianMatricesManager::computeInterpolationMatrices(SurfUnstructured * disciplineMesh, PatchNumberingInfo * disciplineNumberingInfo,
        SurfUnstructured * neutralMesh, PatchNumberingInfo * neutralNumberingInfo) {

    coupling::interpolateFromToMatrix(&m_DisciplineToNeutralInterpolatorJacobian,disciplineMesh,disciplineNumberingInfo,
            neutralMesh,neutralNumberingInfo);
    coupling::interpolateFromToMatrix(&m_NeutralToDisciplineInterpolatorJacobian,neutralMesh,neutralNumberingInfo,
            disciplineMesh,disciplineNumberingInfo);

}


StencilScalarSolverHandler::StencilScalarSolverHandler(const std::string & prefix, bool debug) : StencilScalarSolver(prefix,debug) {
}

const Mat & StencilScalarSolverHandler::getEllipticOperatorMatrix() const {
    return m_A;
}
