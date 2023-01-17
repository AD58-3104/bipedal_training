/**
 * @file MPCExample.cpp
 * @author Giulio Romualdi
 * @copyright Released under the terms of the BSD 3-Clause License
 * @date 2018
 */

// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"

// eigen
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include "gnuplot.h"
#include <cxxabi.h>

static constexpr bool enable_sparse_display = false;
using namespace Eigen;

void sparseDisplay(Eigen::SparseMatrix<double> matrix)
{
    if (enable_sparse_display)
    {
        GnuplotPipe gp;
        std::ofstream ofs;
        ofs.open("sparse_data.dat", std::ios::out);
        size_t row = matrix.rows();
        size_t col = matrix.cols();
        for (int k = 0; k < matrix.outerSize(); ++k)
            for (SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
            {
                ofs << it.col() << " " << -it.row() << "\n";
            }
        gp.sendLine("set terminal wxt size 1280,960");
        // gp.sendLine("set terminal wxt size 640,480");
        gp.sendLine("set xrange [-1:" + std::to_string(col) + "]");
        gp.sendLine("set yrange [-" + std::to_string(row + 1) + ":1]");
        gp.sendLine("plot 'sparse_data.dat' using 1:2");
    }
}

template <typename T>
void showTypeName(T &&tp)
{
    int tmp = 0;
    std::cout << abi::__cxa_demangle(typeid(decltype(std::forward<T>(tp))).name(), 0, 0, &tmp) << std::endl;
}

void showResult()
{
    GnuplotPipe gp;
    gp.sendLine("set terminal wxt size 1280,960");
    gp.sendLine("set xrange [0:3]");
    gp.sendLine("set yrange [-0.5:0.5]");
    gp.sendLine("plot 'x_data.dat' using 1:2");
}

Eigen::VectorXd generateRefTrajectory(const int32_t &step, const int32_t &horizon_length)
{
    static constexpr double T = 0.01;              // サンプリング周期 (s)
    static constexpr int_fast64_t start_step = 2;  // 10 * T = 0.1秒後に歩き出す
    static constexpr int_fast64_t cycle_step = 40; // 何サイクル毎に足踏みするか、 cycle_step * T = 周期(s)
    static constexpr double step_length = 0.3;     // 一歩の大きさ(m)
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(horizon_length);
    for (int32_t i = 0; i < horizon_length; i++)
    {
        if ((step + i) > start_step)
        {
            if (((step + i) / cycle_step) % 2 == 0)
            {
                ret(i, 0) = step_length;
            }
            else
            {
                ret(i, 0) = -step_length;
            }
        }
    }
    return ret;
}

template <size_t X_SIZE, size_t U_SIZE>
void castMPCToQPHessian(const Eigen::DiagonalMatrix<double, X_SIZE> &Q, const Eigen::DiagonalMatrix<double, U_SIZE> &R, const int &mpcWindow,
                        Eigen::SparseMatrix<double> &hessianMatrix)
{

    hessianMatrix.resize(X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow);

    // populate hessian matrix
    for (int i = 0; i < X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow; i++)
    {
        if (i < X_SIZE * (mpcWindow + 1))
        {
            int posQ = i % X_SIZE;
            float value = Q.diagonal()[posQ];
            if (value != 0)
                hessianMatrix.insert(i, i) = value;
        }
        else
        {
            int posR = i % U_SIZE;
            float value = R.diagonal()[posR];
            if (value != 0)
                hessianMatrix.insert(i, i) = value;
        }
    }
    sparseDisplay(hessianMatrix);
}

// todo xRefの形変える
template <size_t X_SIZE, size_t U_SIZE, size_t Z_SIZE, size_t mpcWindow>
void castMPCToQPGradient(const Eigen::DiagonalMatrix<double, Z_SIZE * mpcWindow + 1> &Q, const Eigen::Matrix<double, Z_SIZE * mpcWindow + 1, 1> &zRef,
                         const Eigen::Matrix<double, Z_SIZE, X_SIZE> &C, Eigen::VectorXd &gradient)
{

    Eigen::Matrix<double, Z_SIZE * mpcWindow + 1, 1> Qz_ref;
    Qz_ref = Q * (-zRef);
    // std::cout << "Qz_ref\n" <<Qz_ref << "\nvalue" <<std::endl;
    // std::cout << "X_SISE "<< X_SIZE << " mpcWindow " << mpcWindow << std::endl;
    // populate the gradient vector
    gradient = Eigen::VectorXd::Zero(X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, 1); // ここでのX_SIZEはC_SIZEを表す
    for (int i = 0; i * X_SIZE < X_SIZE * (mpcWindow + 1); i++)
    {
        int posQ = i;
        float value = Qz_ref(posQ, 0);
        // std::cout << posQ << " ←posQ " << value << std::endl;
        gradient.block(i * X_SIZE, 0, X_SIZE, 1) = value * C.transpose();
        // std::cout << C << std::endl;
        // std::cout << " i " << i << std::endl;
    }
    // std::cout << "gradient\n"<< gradient << std::endl;
}

// template <size_t X_SIZE, size_t U_SIZE, size_t Z_SIZE>
// void castMPCToQPConstraintMatrix(const Eigen::Matrix<double, X_SIZE, X_SIZE> &dynamicMatrix, const Eigen::Matrix<double, X_SIZE, U_SIZE> &controlMatrix,
//                                  int mpcWindow, Eigen::SparseMatrix<double> &constraintMatrix)
// {
//     constraintMatrix.resize(X_SIZE * (mpcWindow + 1) + X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow);

//     // populate linear constraint matrix
//     for (int i = 0; i < X_SIZE * (mpcWindow + 1); i++)
//     {
//         constraintMatrix.insert(i, i) = -1;
//     }

//     for (int i = 0; i < mpcWindow; i++)
//         for (int j = 0; j < X_SIZE; j++)
//             for (int k = 0; k < X_SIZE; k++)
//             {
//                 float value = dynamicMatrix(j, k);
//                 if (value != 0)
//                 {
//                     constraintMatrix.insert(X_SIZE * (i + 1) + j, X_SIZE * i + k) = value;
//                 }
//             }

//     for (int i = 0; i < mpcWindow; i++)
//         for (int j = 0; j < X_SIZE; j++)
//             for (int k = 0; k < U_SIZE; k++)
//             {
//                 float value = controlMatrix(j, k);
//                 if (value != 0)
//                 {
//                     constraintMatrix.insert(X_SIZE * (i + 1) + j, U_SIZE * i + k + X_SIZE * (mpcWindow + 1)) = value;
//                 }
//             }

//     for (int i = 0; i < X_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow; i++)
//     {
//         constraintMatrix.insert(i + (mpcWindow + 1) * X_SIZE, i) = 1;
//     }
// }

// template <size_t X_SIZE, size_t U_SIZE, size_t Z_SIZE>
// void castMPCToQPConstraintVectors(const Eigen::Matrix<double, Z_SIZE, 1> &zMax, const Eigen::Matrix<double, Z_SIZE, 1> &zMin,
//                                   const Eigen::Matrix<double, U_SIZE, 1> &uMax, const Eigen::Matrix<double, U_SIZE, 1> &uMin,
//                                   const Eigen::Matrix<double, X_SIZE, 1> &x0,
//                                   const int &mpcWindow, Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound)
// {
//     // evaluate the lower and the upper inequality vectors
//     Eigen::VectorXd lowerInequality = Eigen::MatrixXd::Zero(Z_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, 1);
//     Eigen::VectorXd upperInequality = Eigen::MatrixXd::Zero(Z_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, 1);
//     for (int i = 0; i < mpcWindow + 1; i++)
//     {
//         lowerInequality.block(Z_SIZE * i, 0, Z_SIZE, 1) = zMin;
//         upperInequality.block(Z_SIZE * i, 0, Z_SIZE, 1) = zMax;
//     }
//     for (int i = 0; i < mpcWindow; i++)
//     {
//         lowerInequality.block(U_SIZE * i + Z_SIZE * (mpcWindow + 1), 0, U_SIZE, 1) = uMin;
//         upperInequality.block(U_SIZE * i + Z_SIZE * (mpcWindow + 1), 0, U_SIZE, 1) = uMax;
//     }

//     // evaluate the lower and the upper equality vectors
//     Eigen::VectorXd lowerEquality = Eigen::MatrixXd::Zero(X_SIZE * (mpcWindow + 1), 1);
//     Eigen::VectorXd upperEquality;
//     lowerEquality.block(0, 0, X_SIZE, 1) = -x0;
//     upperEquality = lowerEquality;
//     lowerEquality = lowerEquality;

//     // merge inequality and equality vectors
//     lowerBound = Eigen::MatrixXd::Zero(X_SIZE * (mpcWindow + 1) + Z_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, 1);
//     lowerBound << lowerEquality,
//         lowerInequality;

//     upperBound = Eigen::MatrixXd::Zero(X_SIZE * (mpcWindow + 1) + Z_SIZE * (mpcWindow + 1) + U_SIZE * mpcWindow, 1);
//     upperBound << upperEquality,
//         upperInequality;
// }

// template<size_t X_SIZE>
// void updateConstraintVectors(const Eigen::Matrix<double, X_SIZE, 1> &x0,
//                              Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound)
// {
//     lowerBound.block(0, 0, X_SIZE, 1) = -x0;
//     upperBound.block(0, 0, X_SIZE, 1) = -x0;
// }

// template<size_t Z_SIZE>
// double getErrorNorm(const Eigen::Matrix<double, Z_SIZE, 1> &z,
//                     const Eigen::Matrix<double, Z_SIZE, 1> &zRef)
// {
//     // evaluate the error
//     Eigen::Matrix<double, Z_SIZE, 1> error = z - zRef;

//     // return the norm
//     return error.norm();
// }

int main()
{
    // set the preview window
    static constexpr double hCoM = 0.6;
    static constexpr double g = 9.81;
    static constexpr double T = 0.01;
    static constexpr int32_t mpcWindow = 10; // horizon length
    // number of iteration steps
    static constexpr int32_t numberOfSteps = 300;
    static constexpr int32_t Mu = 1;
    static constexpr int32_t Nx = 3;
    static constexpr int32_t Zx = 1;
    static constexpr int32_t num_of_variables = Nx * (numberOfSteps + 1) + Mu * numberOfSteps;

    // allocate the dynamics matrices
    Eigen::Matrix<double, Nx, Nx> A;
    A << 1.0f, T, T * T / 2.0f,
        0, 1.0f, T,
        0, 0, 1.0f;
    Eigen::Matrix<double, Nx, Mu> B;
    B << T * T * T / 6.0f,
        T * T / 2.0f,
        T;
    Eigen::Matrix<double, Zx, Nx> C;
    C << 1.0f, 0, -hCoM / g;

    // allocate the constraints vector
    Eigen::Matrix<double, Zx, 1> zMax;
    Eigen::Matrix<double, Zx, 1> zMin;
    Eigen::Matrix<double, Mu, 1> uMax;
    Eigen::Matrix<double, Mu, 1> uMin;
    zMax << 0.8;
    zMin << -0.8;
    uMax << 1.3;
    uMin << -1.3;

    // allocate the weight matrices
    // ホライゾン長に渡る、書く予測ステップ毎のZxに対するコスト。ここではZxが1次元なので、mpcWindowの数がQのサイズになる
    Eigen::DiagonalMatrix<double, Zx * mpcWindow + 1> Q;
    // uに掛けるコスト。こちらはホライゾン長に渡って共通のものを掛ける。サイズ = Mu
    Eigen::DiagonalMatrix<double, Mu> R;
    Q.setIdentity();
    R.setIdentity();
    Q = Q * 1000.0;
    R = R * 1.0;

    // allocate the initial and the reference state space
    Eigen::Matrix<double, Nx, 1> x0;
    Eigen::Matrix<double, Zx * mpcWindow + 1, 1> zRef;

    // allocate QP problem matrices and vectores
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;

    // set the initial and the desired states
    x0 << 0, 0, 0;
    zRef = generateRefTrajectory(0, mpcWindow + 1);

    // cast the MPC problem as QP problem
    // castMPCToQPHessian(Q, R, mpcWindow, hessian);
    std::cout << "previous malloc" << std::endl;
    castMPCToQPGradient<Nx, Mu, Zx, mpcWindow>(Q, zRef, C, gradient);
    // castMPCToQPConstraintMatrix(A, B, mpcWindow, linearMatrix);
    // castMPCToQPConstraintVectors(zMax, zMin, uMax, uMin, x0, mpcWindow, lowerBound, upperBound);

    // // instantiate the solver
    // OsqpEigen::Solver solver;

    // // settings
    // // solver.settings()->setVerbosity(false);
    // solver.settings()->setWarmStart(true);

    // // set the initial data of the QP solver
    // solver.data()->setNumberOfVariables(12 * (mpcWindow + 1) + 4 * mpcWindow);
    // solver.data()->setNumberOfConstraints(2 * 12 * (mpcWindow + 1) + 4 * mpcWindow);
    // if (!solver.data()->setHessianMatrix(hessian))
    //     return 1;
    // if (!solver.data()->setGradient(gradient))
    //     return 1;
    // if (!solver.data()->setLinearConstraintsMatrix(linearMatrix))
    //     return 1;
    // if (!solver.data()->setLowerBound(lowerBound))
    //     return 1;
    // if (!solver.data()->setUpperBound(upperBound))
    //     return 1;

    // // instantiate the solver
    // if (!solver.initSolver())
    //     return 1;

    // // controller input and QPSolution vector
    // Eigen::Vector4d ctr;
    // Eigen::VectorXd QPSolution;

    // for (int i = 0; i < numberOfSteps; i++)
    // {

    //     // solve the QP problem
    //     if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
    //         return 1;

    //     // get the controller input
    //     QPSolution = solver.getSolution();
    //     ctr = QPSolution.block(Nx * (mpcWindow + 1), 0, Mu, 1);

    //     // save data into file
    //     auto x0Data = x0.data();

    //     // propagate the model
    //     x0 = A * x0 + B * ctr;

    //     // update the constraint bound
    //     updateConstraintVectors(x0, lowerBound, upperBound);
    //     if (!solver.updateBounds(lowerBound, upperBound))
    //         return 1;
    // }
    return 0;
}