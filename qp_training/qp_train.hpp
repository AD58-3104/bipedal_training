#include <Eigen/Dense>
#include <Eigen/Geometry>

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
 * @brief implementation qp in the trajectory free linear MPC
 *  ref paper https://ieeexplore.ieee.org/document/4115592
 */
void solve_qp()
{
    static constexpr double hCoM = 0.8;
    static constexpr double g = 9.81;
    static constexpr double T = 0.005;
    static constexpr int32_t N = 10; // horizon length
    Eigen::Matrix<double, 3, 3> A;
    Eigen::Matrix<double, 3, 1> B;
    Eigen::Matrix<double, 1, 3> C;
    Eigen::Matrix<double, 3, 1> x_k;
    Eigen::Matrix<double, N, 3> Px;
    Eigen::Matrix<double, N, N> Pu;
    Eigen::Matrix<double, N, 1> u;     // 探したい制御入力
    Eigen::Matrix<double, N, 1> Z_ref; // 参照軌道
    Eigen::Matrix<double, N, 1> q_vec; // 重み asDiagonal()で対角になってくれるのでベクトル
    Eigen::Matrix<double, N, 1> r_vec; // 重み asDiagonal()で対角になってくれるのでベクトル
    Eigen::Matrix<double, N, N> hessian;
    Eigen::Matrix<double, N, N> gradient;
    Eigen::Matrix<double, N, N> Q;
    Eigen::Matrix<double, N, N> R;
    q_vec << 1, 1, 1, 0.3, 1, 0.5, 0.4, 0.3, 0.2, 0.1;
    r_vec << 1, 0.9, 0.8, 0.3, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1;
    Q = q_vec.asDiagonal();
    R = r_vec.asDiagonal();
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
    hessian = 2.0f * H;
    std::cout << hessian << std::endl;
    // gradient = H * u + ((x_k.transpose() * Px.transpose() - Z_ref.transpose())*Q*Pu);
    std::cout << std::endl;
    return;
}