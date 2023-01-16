#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "OsqpEigen/OsqpEigen.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "gnuplot.h"
#include <cxxabi.h>

using namespace Eigen;

static constexpr bool enable_sparse_display = false;

void sparseDisplay(Eigen::SparseMatrix<double> matrix)
{
    if (enable_sparse_display)
    {
        GnuplotPipe gp;
        std::fstream ofs;
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

namespace trajectory_free_LMPC
{

    Eigen::VectorXd generateRefTrajectory(const int32_t &step, const int32_t &horizon_length)
    {
        static constexpr double T = 0.01;              // サンプリング周期 (s)
        static constexpr int_fast64_t start_step = 10; // 10 * T = 0.1秒後に歩き出す
        static constexpr int_fast64_t cycle_step = 40; // 何サイクル毎に足踏みするか、 cycle_step * T = 周期(s)
        static constexpr double step_length = 0.3;     // (m)
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

    template <size_t X_SIZE, size_t U_SIZE, size_t Hp>
    void getConstraintMatrix(const Eigen::Matrix<double, X_SIZE, X_SIZE> &A, const Eigen::Matrix<double, X_SIZE, U_SIZE> &B, Eigen::SparseMatrix<double> &constraintMatrix)
    {
        // QPのxの大きさ、システムの状態の方のxではない。
        static constexpr size_t size_of_vector_x_in_QP = X_SIZE * (Hp + 1) + U_SIZE * Hp;
        using constraintMat_t = Eigen::Matrix<double, X_SIZE *(Hp + 1) * 2 + U_SIZE * Hp, size_of_vector_x_in_QP>;
        // こんな事をしたらsparseを使う意味が無くなるが、面倒なので...。
        Eigen::MatrixXd tmp_constraintMat = Eigen::MatrixXd::Zero(X_SIZE * (Hp + 1) * 2 + U_SIZE * Hp, size_of_vector_x_in_QP);
        constraintMatrix.resize(X_SIZE * (Hp + 1) * 2 + U_SIZE * Hp, size_of_vector_x_in_QP);
        // 　↑の1つめの　X_SIZE * (Horizon_length + 1) = QPのxの中に状態を持たせる為に使う変数の分。 x0を不等式制約とし、
        //                                         そのx0を使ってそれ以降のダイナミクスを順々に計算し、それを = 0の制約とすると、
        //                                         なんとx(n) = Ax(n-1) + Bu(n-1)の式を制約で作る事ができ、QPのxの中に状態を変数として導入出来る
        // 　↑の2つめの　X_SIZE * (Horizon_length + 1) = 状態の不等式制約
        // 　↑の3つめの　U_SIZE * Horizon_length = 入力の不等式制約
        for (size_t i = 0; i < Hp; i++)
        {
            tmp_constraintMat.block(i * X_SIZE, i * X_SIZE, X_SIZE, X_SIZE) = -Eigen::MatrixXd::Identity(X_SIZE, X_SIZE);
            tmp_constraintMat.block((i + Hp) * X_SIZE, i * X_SIZE, X_SIZE, X_SIZE) = Eigen::MatrixXd::Identity(X_SIZE, X_SIZE);
            tmp_constraintMat.block((Hp)*X_SIZE * 2 + i * U_SIZE, Hp * X_SIZE + i * U_SIZE, U_SIZE, U_SIZE) = Eigen::MatrixXd::Identity(U_SIZE, U_SIZE);
            if (i > 0)
            {
                tmp_constraintMat.block((i)*X_SIZE, (i - 1) * X_SIZE, X_SIZE, X_SIZE) = A;
                tmp_constraintMat.block((i)*X_SIZE, X_SIZE * (Hp) + (i - 1) * U_SIZE, X_SIZE, U_SIZE) = B;
            }
        }
        std::cout << "---- constraint -----\n";
        // std::cout << tmp_constraintMat << std::endl;

        constraintMatrix = tmp_constraintMat.sparseView();
        sparseDisplay(constraintMatrix);
        return;
    }

    template <size_t X_SIZE, size_t U_SIZE, size_t Hp>
    void getConstraintVector(const Eigen::Matrix<double, X_SIZE, 1> &xMax, const Eigen::Matrix<double, X_SIZE, 1> &xMin,
                             const Eigen::Matrix<double, U_SIZE, 1> &uMax, const Eigen::Matrix<double, U_SIZE, 1> &uMin,
                             const Eigen::Matrix<double, X_SIZE, 1> &x0, Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound)
    {
        Eigen::VectorXd lowerInequality = Eigen::MatrixXd::Zero(X_SIZE * (Hp + 1) + U_SIZE * Hp, 1);
        Eigen::VectorXd upperInequality = Eigen::MatrixXd::Zero(X_SIZE * (Hp + 1) + U_SIZE * Hp, 1);
        for (int i = 0; i < Hp + 1; i++)
        {
            lowerInequality.block(X_SIZE * i, 0, X_SIZE, 1) = xMin;
            upperInequality.block(X_SIZE * i, 0, X_SIZE, 1) = xMax;
        }
        for (int i = 0; i < Hp; i++)
        {
            lowerInequality.block(U_SIZE * i + X_SIZE * (Hp + 1), 0, U_SIZE, 1) = uMin;
            upperInequality.block(U_SIZE * i + X_SIZE * (Hp + 1), 0, U_SIZE, 1) = uMax;
        }
        Eigen::VectorXd lowerEquality = Eigen::MatrixXd::Zero(X_SIZE * (Hp + 1), 1);
        lowerEquality.block(0, 0, X_SIZE, 1) = -x0;
        Eigen::VectorXd upperEquality = Eigen::MatrixXd::Zero(X_SIZE * (Hp + 1), 1);
        upperEquality = lowerEquality;
        lowerBound = upperBound = Eigen::MatrixXd::Zero(lowerEquality.rows() + lowerInequality.rows(), 1);
        lowerBound
            << lowerEquality,
            lowerInequality;
        upperBound << upperEquality,
            upperInequality;
    }

    template <size_t X_SIZE>
    void updateConstraintVectors(const Eigen::Matrix<double, X_SIZE, 1> &x0,
                                 Eigen::VectorXd &lowerBound, Eigen::VectorXd &upperBound)
    {
        lowerBound.block(0, 0, X_SIZE, 1) = -x0;
        upperBound.block(0, 0, X_SIZE, 1) = -x0;
    }

    // todo sparseを返す様にするべき
    template <size_t VAL_NUM, size_t TOTAL_U_SIZE>
    // Eigen::Matrix<double, VAL_NUM, VAL_NUM> expandHessianSize(const Eigen::Matrix<double, TOTAL_U_SIZE, TOTAL_U_SIZE> &hessian)
    Eigen::SparseMatrix<double> expandHessianSize(const Eigen::MatrixXd &hessian)
    {
        const size_t original_hessian_location = VAL_NUM - TOTAL_U_SIZE;
        Eigen::SparseMatrix<double> expand_hessian;
        expand_hessian.resize(VAL_NUM, VAL_NUM);
        for (size_t i = 0; i < hessian.cols(); ++i)
        {
            Eigen::Matrix<double, VAL_NUM, 1> tmp;
            tmp << Eigen::VectorXd::Zero(original_hessian_location), hessian.col(i);
            expand_hessian.col(original_hessian_location + i) = tmp.sparseView();
        }
        // Eigen::MatrixXd bighessian = Eigen::MatrixXd::Zero(VAL_NUM, VAL_NUM);
        // bighessian.block(original_hessian_location, original_hessian_location, hessian.cols(), hessian.rows()) = hessian;
        Eigen::MatrixXd expH(expand_hessian);
        std::cout << "!!! original hessian !!!" << std::endl;
        std::cout << hessian << std::endl;        
        std::cout << "!!! expand hessian !!!" << std::endl;
        std::cout << expH << std::endl;
        return expand_hessian;
    }

    template <size_t VAL_NUM, size_t TOTAL_U_SIZE>
    Eigen::Matrix<double, 1, VAL_NUM> expandGradientSize(const Eigen::Matrix<double, 1, TOTAL_U_SIZE> &gradient)
    {
        const size_t expand_col_size = VAL_NUM - TOTAL_U_SIZE;
        Eigen::Matrix<double, 1, VAL_NUM> ret;
        ret << Eigen::MatrixXd::Zero(1, expand_col_size), gradient;
        std::cout << "=== original gradient ===" << std::endl
                    << gradient << std::endl;
        std::cout << "=== expand gradient ===" << std::endl
                    << ret << std::endl;
        return ret;
    }

    /**
     * @brief implementation qp in the trajectory free linear MPC
     *  ref paper https://ieeexplore.ieee.org/document/4115592
     */
    void solve_qp()
    {
        static constexpr double hCoM = 0.6;
        static constexpr double g = 9.81;
        static constexpr double T = 0.01;
        static constexpr int32_t Horizon_length = 100; // horizon length
        static constexpr int32_t N = Horizon_length;   // horizon length
        const int numberOfSteps = 300;
        static constexpr int32_t Mu = 1;
        static constexpr int32_t Nx = 3;
        static constexpr int32_t Zx = 1;
        static constexpr int32_t num_of_variables = Nx * (Horizon_length + 1) + Mu * Horizon_length;

        // x = "A"x + Bu
        Eigen::Matrix<double, Nx, Nx> A;

        // x = Ax + "B"u
        Eigen::Matrix<double, Nx, Mu> B;
        // Z = "C"x
        Eigen::Matrix<double, 1, Nx> C;
        // 現在の状態
        Eigen::Matrix<double, Nx, 1> x_k;
        Eigen::Matrix<double, Horizon_length, Nx> Px;
        Eigen::Matrix<double, Horizon_length, Horizon_length * Mu> Pu;
        Eigen::Matrix<double, Horizon_length, Mu> u;    // 探したい制御入力
        Eigen::Matrix<double, Horizon_length, 1> Z_ref; // 参照軌道
        Eigen::Matrix<double, Horizon_length, 1> q_vec; // 重み asDiagonal()で対角になってくれるのでベクトル
        Eigen::Matrix<double, Horizon_length, 1> r_vec; // 重み asDiagonal()で対角になってくれるのでベクトル

        Eigen::Matrix<double, Horizon_length * Zx, Horizon_length * Zx> Q;
        Eigen::Matrix<double, Horizon_length * Mu, Horizon_length * Mu> R;
        Eigen::Matrix<double, Nx, Nx> q_1;
        Eigen::Matrix<double, Mu, Mu> r_1;
        std::cout << " Q, R " << std::endl;
        Q.setZero();
        R.setZero();
        // Q.diagonal() << 1, 1, 1, 0.3, 1, 0.5, 0.4, 0.3, 0.2, 0.1;
        // R.diagonal() << 1, 0.9, 0.8, 0.3, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;
        Q = decltype(Q)::Identity() * 100000.0f;
        R = decltype(R)::Identity() * 1.0f;

        Eigen::SparseMatrix<double> hessian;
        Eigen::Matrix<double, 1, num_of_variables> gradient;

        Eigen::SparseMatrix<double> ConstraintMatrix;
        Eigen::VectorXd lowerBound;
        Eigen::VectorXd upperBound;

        static constexpr double umax = 0.4 * 1000;
        static constexpr double umin = -0.4 * 1000;
        static constexpr double zmax = 100.0;
        static constexpr double zmin = -100.0;

        Eigen::Matrix<double, Nx, 1> xMax = Eigen::MatrixXd::Constant(Nx, 1, zmax);
        Eigen::Matrix<double, Nx, 1> xMin = Eigen::MatrixXd::Constant(Nx, 1, zmin);

        Eigen::Matrix<double, Mu, 1> uMax = Eigen::MatrixXd::Constant(Mu, 1, umax);
        Eigen::Matrix<double, Mu, 1> uMin = Eigen::MatrixXd::Constant(Mu, 1, umin);

        std::cout << "\n--- Q ---" << std::endl;
        std::cout << Q;
        std::cout << "\n--- R ---" << std::endl;
        std::cout << R;
        A << 1.0f, T, T * T / 2.0f,
            0, 1.0f, T,
            0, 0, 1.0f;
        B << T * T * T / 6.0f,
            T * T / 2.0f,
            T;
        C << 1.0f, 0, -hCoM / g;

        x_k << 1, 1, 1; // 一番最初の状態で初期化
        // construct Px
        for (size_t n = 0; n < Horizon_length; n++)
        {
            Px(n, 0) = 1.0;
            Px(n, 1) = T * (n + 1);
            Px(n, 2) = T * T * (n + 1) * (n + 1) / 2 - hCoM / g;
        }

        for (size_t n_row = 0; n_row < Horizon_length; n_row++)
        {
            // static Eigen::Matrix<double, Nx, 1> accumAB = B;
            // static Eigen::Matrix<double, Nx, Nx> powerA = Eigen::Matrix<double, Nx, Nx>::Identity();
            for (size_t n_col = 0; n_col < Horizon_length; n_col++)
            {
                if (n_col == 0)
                {
                    Pu(n_row, n_col) = (1.0 + 3.0 * n_row + 3.0 * n_row * n_row) * T * T * T / 6.0 - T * hCoM / g;
                    std::cout << "Pu(" << n_row << "," << n_col << ") = " << Pu(n_row, n_col) << std::endl;
                    // if (n_row > 0)
                    // {
                    //     powerA = powerA * A;
                    //     accumAB = accumAB + powerA * B;
                    // }
                    // std::cout << "calc " << C * accumAB << std::endl;
                }
                else if (n_col <= n_row)
                {
                    Pu(n_row, n_col) = Pu(n_row - 1, n_col - 1);
                }
                else
                {
                    Pu(n_row, n_col) = 0;
                }
            }
        }
        std::cout << Pu << std::endl;

        getConstraintMatrix<Nx, Mu, Horizon_length>(A, B, ConstraintMatrix);
        getConstraintVector<Nx, Mu, Horizon_length>(xMax, xMin, uMax, uMin, x_k, lowerBound, upperBound);

        // std::cout << "--- lowerBound ---" << std::endl;
        // std::cout << lowerBound << std::endl;
        // std::cout << "--- upperBound ---" << std::endl;
        // std::cout << upperBound << std::endl;

        // Eigen::Matrix<double, Horizon_length * Mu, Horizon_length *Mu> H = 1.0f / 2.0f * (Pu.transpose() * Q * Pu + R);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(Horizon_length * Mu, Horizon_length * Mu);
        H = 1.0f / 2.0f * (Pu.transpose() * Q * Pu + R);
        // Eigen::Matrix<double, num_of_variables, num_of_variables> expand_H = expandHessianSize<num_of_variables, Horizon_length * Mu>(H);
        hessian = 2 * expandHessianSize<num_of_variables, Horizon_length * Mu>(H);
        // hessian = (2.0f * expand_H).sparseView();
        sparseDisplay(hessian);
        Eigen::Matrix<double, 1, Horizon_length * Mu> original_G;
        original_G = ((x_k.transpose() * Px.transpose() - Z_ref.transpose()) * Q * Pu);
        gradient = expandGradientSize<num_of_variables, Horizon_length * Mu>(original_G);
        // std::cout << "--- gradient ---\n";
        // std::cout << gradient << std::endl;

        OsqpEigen::Solver solver;
        solver.settings()->setWarmStart(true);
        solver.data()->setNumberOfVariables(num_of_variables);
        solver.data()->setNumberOfConstraints(Nx * (Horizon_length + 1) * 2 + Mu * Horizon_length);
        // 　↑の1つめの　Nx * (Horizon_length + 1) = QPのxの中に状態を持たせる為に使う変数の分。 x0を不等式制約とし、
        //                                         そのx0を使ってそれ以降のダイナミクスを順々に計算し、それを = 0の制約とすると、
        //                                         なんとx(n) = Ax(n-1) + Bu(n-1)の式を制約で作る事ができ、QPのxの中に状態を変数として導入出来る
        // 　↑の2つめの　Nx * (Horizon_length + 1) = 状態の不等式制約
        // 　↑の3つめの　Mu * Horizon_length = 入力の不等式制約
        if (!solver.data()->setHessianMatrix(hessian))
            return;
        if (!solver.data()->setGradient(gradient))
            return;
        if (!solver.data()->setLinearConstraintsMatrix(ConstraintMatrix))
            return;
        if (!solver.data()->setLowerBound(lowerBound))
            return;
        if (!solver.data()->setUpperBound(upperBound))
            return;
        if (!solver.initSolver())
            return;

        Eigen::Matrix<double, Mu, 1> ctr;
        Eigen::VectorXd QPSolution;

        std::vector<double> calculated_x_list;
        std::ofstream ofs;
        ofs.open("x_data.dat");
        auto updateGradient = [&](int current_step) -> void
        {
            Z_ref = generateRefTrajectory(current_step, Horizon_length);
            original_G = ((x_k.transpose() * Px.transpose() - Z_ref.transpose()) * Q * Pu);
            gradient = expandGradientSize<num_of_variables, Horizon_length * Mu>(original_G);
            if (!solver.data()->setGradient(gradient))
            {
                throw std::runtime_error(" some false in setGradient");
                return;
            }
        };
        for (int i = 0; i < numberOfSteps; i++)
        {

            // solve the QP problem
            if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
            {
                std::cout << "solver error" << std::endl;
                return;
            }

            // get the controller input
            QPSolution = solver.getSolution();
            ctr = QPSolution.block(Nx * (Horizon_length + 1), 0, Mu, 1);

            // save data into file
            auto x0Data = x_k.data();
            // ofs << i * T << " " << x_k(0, 0) << std::endl;
            ofs << i * T << " " << x_k.transpose() << std::endl;

            // propagate the model
            x_k = A * x_k + B * ctr;
            updateGradient(i);
            // update the constraint bound
            updateConstraintVectors<Nx>(x_k, lowerBound, upperBound);
            if (!solver.updateBounds(lowerBound, upperBound))
                return;
            if (i == numberOfSteps - 1)
            {
                std::cout << "----answer----" << std::endl;
                std::cout << QPSolution << std::endl;
                std::cout << "answer cols " << QPSolution.cols() << " rows " << QPSolution.rows() << std::endl;
                std::cout << " control \n"
                          << ctr << std::endl;
            }
        }

        showResult();
        // std::cout << lowerBound.transpose() << std::endl;
        // std::cout << upperBound.transpose() << std::endl;
        // std::cout << "--- Pu =---" << std::endl;
        // std::cout << Pu << std::endl;
        // std::cout << "\n---  Px  ---" << std::endl
        //           << Px << std::endl;
        // std::cout << "--- Z_ref --- " << std::endl
        //           << Z_ref.transpose() << std::endl;
        // std::cout << " H " << std::endl << H << std::endl;
        return;
    }
}