#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "OsqpEigen/OsqpEigen.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace Eigen;

void solve()
{
    float T = 0.01;
    Eigen::Matrix<float, 3, 3> A, ansA;
    Eigen::Matrix<float, 3, 1> B, ansB;
    A << 1, T, T * T / 2,
        0, 1, T,
        0, 0, 1;
    B << T * T * T / 6,
        T * T / 2,
        T;
    ansA = A;
    ansB = B;
    for (size_t i = 0; i <= 10; i++)
    {
        ansA = ansA * A;
        ansB = A * ansB;
    }
    std::cout << ansA << std::endl;
    std::cout << ansB << std::endl;
}

/**
 * @brief
 *
 * @param step 現在のステップ
 * @param horizon_length ホライゾンの長さ
 * @return std::vector<double>
 */
std::vector<double> generateRefTrajectory(const uint32_t &step, const uint32_t &horizon_length)
{
    constexpr double step_T = 1.5f; // 秒に一歩
    static constexpr double T = 0.005;
    constexpr uint32_t loop_len = std::round(step_T / T);
    constexpr double rad = 2.0f * M_PI / (double)loop_len;
    std::vector<double> ret;
    ret.reserve(horizon_length);
    for (uint32_t i = 0; i < horizon_length; i++)
    {
        ret.push_back(rad * static_cast<double>(step + i));
    }
    return ret;
}

/**
 * @brief implementation qp in the trajectory free linear MPC
 *  ref paper https://ieeexplore.ieee.org/document/4115592
 */
void solve_qp()
{
    static constexpr double hCoM = 0.8;
    static constexpr double g = 9.81;
    static constexpr double T = 0.005;
    static constexpr int32_t Horizon_length = 10; // horizon length
    static constexpr int32_t N = Horizon_length;  // horizon length
    static constexpr int32_t Mu = 1;
    static constexpr int32_t Nx = 3;

    Eigen::Matrix<double, Nx, Nx> A;
    Eigen::Matrix<double, Nx, Mu> B;
    Eigen::Matrix<double, 1, Nx> C;
    Eigen::Matrix<double, Nx, 1> x_k;
    Eigen::Matrix<double, Horizon_length, Nx> Px;
    Eigen::Matrix<double, Horizon_length, Horizon_length> Pu;
    Eigen::Matrix<double, Horizon_length, Mu> u;     // 探したい制御入力
    Eigen::Matrix<double, Horizon_length, 1> Z_ref; // 参照軌道
    Eigen::Matrix<double, Horizon_length, 1> q_vec; // 重み asDiagonal()で対角になってくれるのでベクトル
    Eigen::Matrix<double, Horizon_length, 1> r_vec; // 重み asDiagonal()で対角になってくれるのでベクトル
    Eigen::SparseMatrix<double> hessian;
    Eigen::Matrix<double, 1, Horizon_length> gradient;
    Eigen::Matrix<double, Horizon_length, Horizon_length> R;
    Eigen::Matrix<double, Horizon_length, Horizon_length> Q;
    R.setZero();
    Q.setZero();
    Q.diagonal() << 1, 1, 1, 0.3, 1, 0.5, 0.4, 0.3, 0.2, 0.1;
    R.diagonal() << 1, 0.9, 0.8, 0.3, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;
    std::cout << "\n--- Q ---" << std::endl;
    std::cout << Q;
    std::cout << "\n--- R ---" << std::endl;
    std::cout << R;
    A << 1, T, T * T / 2,
        0, 1, T,
        0, 0, 1;
    B << T * T * T / 6,
        T * T / 2,
        T;
    C << 1, 0, hCoM / g;

    x_k << 0, 0, 0; // 一番最初の状態で初期化
    // construct Px
    for (size_t n = 0; n < N; n++)
    {
        Px(n, 0) = 1, Px(n, 1) = T * (n + 1), Px(n, 2) = T * T * (n + 1) * (n + 1) - hCoM / g;
    }
    std::cout << "\n---  Px  ---\n"
              << Px << std::endl;
    for (size_t n_row = 0; n_row < N; n_row++)
    {
        for (size_t n_col = 0; n_col < N; n_col++)
        {
            if (n_col == 0)
            {
                Pu(n_row, n_col) = (1 + 3 * n_row + 3 * n_row * n_row) * T * T * T / 6 - T * hCoM / g;
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
    std::cout << "---  Pu  ---\n"
              << Pu << std::endl;

    Eigen::Matrix<double, N, N> H = 1.0f / 2.0f * (Pu.transpose() * Q * Pu + R);
    hessian = (2.0f * H).sparseView();
    std::cout << hessian << std::endl;
    gradient = ((x_k.transpose() * Px.transpose() - Z_ref.transpose()) * Q * Pu);
    gradient = ((x_k.transpose() * Px.transpose() - Z_ref.transpose()) * Q * Pu);
    std::cout << std::endl;

    // osqp-eigenのサンプル
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    solver.data()->setNumberOfVariables(Nx * (Horizon_length + 1) + Mu * Horizon_length);
    solver.data()->setNumberOfConstraints(Nx * (Horizon_length + 1) * 2 + Mu * Horizon_length);
    //　↑の1つめの　Nx * (Horizon_length + 1) = QPのxの中に状態を持たせる為に使う変数の分。 x0を不等式制約とし、
    //                                        そのx0を使ってそれ以降のダイナミクスを順々に計算し、それを = 0の制約とすると、
    //                                        なんとx(n) = Ax(n-1) + Bu(n-1)の式を制約で作る事ができ、QPのxの中に状態を変数として導入出来る
    //　↑の2つめの　Nx * (Horizon_length + 1) = 状態の不等式制約
    //　↑の3つめの　Mu * Horizon_length = 入力の不等式制約
    // if (!solver.data()->setHessianMatrix(hessian))
    //     return;
    // if (!solver.data()->setGradient(gradient))
    //     return;
    // if (!solver.data()->setLinearConstraintsMatrix(linearMatrix))
    //     return;
    // if (!solver.data()->setLowerBound(lowerBound))
    //     return;
    // if (!solver.data()->setUpperBound(upperBound))
    //     return;
    return;
}