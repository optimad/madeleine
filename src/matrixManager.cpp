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
: m_name(name), m_comm(comm), m_ellipticLinearSystem(ellipticLinearSystem), m_reuseJacobianMatrix(false) {
}

JacobianMatricesManager::~JacobianMatricesManager() {
    MatDestroy(&m_OutputInputJacobian);
    MatDestroy(&m_OutputControlJacobian);
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
    //MatConvert(m_DisciplineToNeutralInterpolatorJacobian,MATDENSE,MAT_INPLACE_MATRIX,&m_DisciplineToNeutralInterpolatorJacobian);
    //MatConvert(m_NeutralToDisciplineInterpolatorJacobian,MATDENSE,MAT_INPLACE_MATRIX,&m_NeutralToDisciplineInterpolatorJacobian);

}


/*!
    Compute derivatives of output with respect to input.
    NB: PETSc is not able to mix sparse and dense matrices during multiplication for every multiplication ordering.
    Inner discipline follows
    We need to perform m_DisciplineToNeutralInterpolatorJacobian * m_EllipticOperatorInverse * m_NeutralToDisciplineInterpolatorJacobian.
    Let's call the matrices involved A,B,C and therefore we want A*B*C, where
    A is sparse, B is dense and C is sparse.
    PETSc does not allow for triple matrix multiplication using MatMatMatMult if B is dense.
    But it allows double matrix multiplication with dense matrix as the left side one, i.e. B*C is allowed.
    Considering that this product is dense, we cannot perform A*(B*C) for the same reason above.
    The only way is to use matrix transposition as trick to do what we need.
    The procedure to get A*B*C is:
    - compute B*C
    - transpose it, i.e. (B*C)^T
    - transpose A, i.e. A^T
    - compute (B*C)^T * A^T
    - transpose is, i.e. ((B*C)^T * A^T)^T = A*B*C

    Outer discipline follows the same reasoning.

    \param[in] emissivity value of the emissivity for outer discipline
    \param[in] isInnerDiscipline boolean value for different derivative calculation for inner and oute disciplines, respectively
*/
void JacobianMatricesManager::computeOutputInputJacobian(MPI_Comm comm, double emissivity, bool isInnerDiscipline) {

    log::cout() << "Compute Output/Input Jacobian." << std::endl;

    if(isInnerDiscipline) {
        Mat EmND_T;
        MatMatMult(m_EllipticOperatorInverse,m_NeutralToDisciplineInterpolatorJacobian,
                MAT_INITIAL_MATRIX,PETSC_DEFAULT,&EmND_T);
        MatTranspose(EmND_T,MAT_INPLACE_MATRIX,&EmND_T);

        //build emissivity identity
        int nofRows,nofCols;
        Mat emissivityIdentity;
        MatGetLocalSize(m_EllipticOperatorInverse,&nofRows,&nofCols);
        MatCreateAIJ(comm,nofRows,nofRows,PETSC_DETERMINE,PETSC_DETERMINE,1,PETSC_NULL,0,PETSC_NULL,&emissivityIdentity);
        MatZeroEntries(emissivityIdentity);
        MatAssemblyBegin(emissivityIdentity,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(emissivityIdentity,MAT_FINAL_ASSEMBLY);
        MatShift(emissivityIdentity,emissivity);

        Mat emiss_EmND;
        MatMatMult(EmND_T,emissivityIdentity,
                       MAT_INITIAL_MATRIX,PETSC_DEFAULT,&emiss_EmND);
        Mat DN_T;
        MatTranspose(m_DisciplineToNeutralInterpolatorJacobian,MAT_INITIAL_MATRIX,&DN_T);
        MatReuse matReuse = MAT_INITIAL_MATRIX;
        if(m_reuseJacobianMatrix) {
            matReuse = MAT_REUSE_MATRIX;
        }
        MatMatMult(emiss_EmND,DN_T,matReuse,PETSC_DEFAULT,&m_OutputInputJacobian);
        m_reuseJacobianMatrix = true;
        MatTranspose(m_OutputInputJacobian,MAT_INPLACE_MATRIX,&m_OutputInputJacobian);

        MatDestroy(&emissivityIdentity);
        MatDestroy(&emiss_EmND);
        MatDestroy(&EmND_T);
        MatDestroy(&DN_T);
    } else {
        Mat emissivityIdentity;
        Mat emissivityIdentityNeutral;
        Mat E_emissivityIdentity;
        Mat E_emissivityIdentity_ND_T;
        Mat DN_T;

        MatDuplicate(m_DisciplineToNeutralInterpolatorJacobian,MAT_COPY_VALUES,&DN_T);
        MatTranspose(DN_T,MAT_INPLACE_MATRIX,&DN_T);

        //build emissivity identity
        int nofRows,nofCols;
        MatGetLocalSize(m_EllipticOperatorInverse,&nofRows,&nofCols);
        MatCreateAIJ(comm,nofRows,nofRows,PETSC_DETERMINE,PETSC_DETERMINE,1,PETSC_NULL,0,PETSC_NULL,&emissivityIdentity);
        MatZeroEntries(emissivityIdentity);
        MatAssemblyBegin(emissivityIdentity,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(emissivityIdentity,MAT_FINAL_ASSEMBLY);
        MatShift(emissivityIdentity,emissivity);

        //
        MatMatMult(m_EllipticOperatorInverse,emissivityIdentity,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&E_emissivityIdentity);
        MatMatMult(E_emissivityIdentity,m_NeutralToDisciplineInterpolatorJacobian,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&E_emissivityIdentity_ND_T);
        MatTranspose(E_emissivityIdentity_ND_T,MAT_INPLACE_MATRIX,&E_emissivityIdentity_ND_T);

        MatReuse matReuse = MAT_INITIAL_MATRIX;
        if(m_reuseJacobianMatrix) {
            matReuse = MAT_REUSE_MATRIX;
        }
        MatMatMult(E_emissivityIdentity_ND_T,DN_T,matReuse,PETSC_DEFAULT,&m_OutputInputJacobian);
        m_reuseJacobianMatrix = true;
        MatTranspose(m_OutputInputJacobian,MAT_INPLACE_MATRIX,&m_OutputInputJacobian);

        MatDestroy(&DN_T);
        MatDestroy(&emissivityIdentity);
        MatDestroy(&E_emissivityIdentity);
        MatDestroy(&E_emissivityIdentity_ND_T);
    }
    log::cout() << "Output/Input Jacobian computed." << std::endl;


}

void JacobianMatricesManager::computeOutputControlJacobian(MPI_Comm comm, double emissivity, bool isInnerDiscipline, double radius,
        const PiercedStorage<double,long> & disciplineData, const SurfUnstructured * disciplineMesh, const PatchNumberingInfo & disciplineNumberingInfo,
        const SurfUnstructured * neutralMesh, const PatchNumberingInfo & neutralNumberingInfo) {

    //Set temperature in a Vec
    int tid = 0;
    if(isInnerDiscipline) {
        tid = 1;
    }
    Vec temperature;
    PetscScalar * t;
    int nofRows,nofCols;
    MatGetLocalSize(m_ellipticLinearSystem.getEllipticOperatorMatrix(),&nofRows,&nofCols);
    VecCreateMPI(disciplineMesh->getCommunicator(), nofRows, PETSC_DETERMINE, &temperature);
    VecGetArray(temperature,&t);
    for(const auto & cell : disciplineMesh->getCells()) {
        if(cell.isInterior()) {
            long cellId = cell.getId();
            long cellConsecutiveId = disciplineNumberingInfo.getCellConsecutiveId(cellId);
            long cellLocalConsecutiveId = cellConsecutiveId - disciplineNumberingInfo.getCellGlobalCountOffset();

            //t[cellLocalConsecutiveId] = disciplineData.at(cellId,tid);

        }

    }
    VecRestoreArray(temperature,&t);
    VecAssemblyBegin(temperature);
    VecAssemblyEnd(temperature);

    //Compute elliptic operator derivative wrt control
    Mat ellipticOperatorControlDerivative;
    MatDuplicate(m_ellipticLinearSystem.getEllipticOperatorMatrix(),MAT_COPY_VALUES,&ellipticOperatorControlDerivative);
    MatScale(ellipticOperatorControlDerivative,1.0 / radius);

    //Compute Jacobian as a Vec
    Vec dAdR_times_T;
    VecDuplicate(temperature,&dAdR_times_T);
    MatMult(ellipticOperatorControlDerivative,temperature,dAdR_times_T);
    Vec ellipticInverse_times_dAdR_times_T;
    VecDuplicate(temperature,&ellipticInverse_times_dAdR_times_T);
    MatMult(m_EllipticOperatorInverse,dAdR_times_T,ellipticInverse_times_dAdR_times_T);
    VecScale(ellipticInverse_times_dAdR_times_T,-1.0);

    Vec controlJacobian;
    MatGetLocalSize(m_DisciplineToNeutralInterpolatorJacobian,&nofRows,&nofCols);
    VecCreateMPI(disciplineMesh->getCommunicator(), nofRows, PETSC_DETERMINE, &controlJacobian);
    MatMult(m_DisciplineToNeutralInterpolatorJacobian,ellipticInverse_times_dAdR_times_T,controlJacobian);

    MatCreateDense(m_comm,nofRows,PETSC_DECIDE,PETSC_DECIDE,1,NULL,&m_OutputControlJacobian);
    MatGetLocalSize(m_OutputControlJacobian,&nofRows,&nofCols);

    const double *values;
    VecGetArrayRead(controlJacobian,&values);
    for(int i = 0; i < nofRows; ++i) {
        long cellId = neutralMesh->getCells().rawFind(i).getId();
        long cellGlobalId = neutralNumberingInfo.getCellConsecutiveId(cellId);
        MatSetValue(m_OutputControlJacobian,cellGlobalId,0,values[i],INSERT_VALUES);

    }
    VecRestoreArrayRead(controlJacobian,&values);
    MatAssemblyBegin(m_OutputControlJacobian,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_OutputControlJacobian,MAT_FINAL_ASSEMBLY);

    VecDestroy(&temperature);
    VecDestroy(&dAdR_times_T);
    VecDestroy(&ellipticInverse_times_dAdR_times_T);
    MatDestroy(&ellipticOperatorControlDerivative);
    VecDestroy(&controlJacobian);
}

Mat & JacobianMatricesManager::getOuputInputJacobian() {
    return m_OutputInputJacobian;
}

Mat & JacobianMatricesManager::getOuputControlJacobian() {
    return m_OutputControlJacobian;
}

StencilScalarSolverHandler::StencilScalarSolverHandler(const std::string & prefix, bool debug) : StencilScalarSolver(prefix,debug) {
}

const Mat & StencilScalarSolverHandler::getEllipticOperatorMatrix() const {
    return m_A;
}
