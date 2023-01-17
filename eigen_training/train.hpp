#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <typeinfo>
#include <cxxabi.h>

using namespace Eigen;

template <typename T>
void showTypeName(T &&tp)
{
    int tmp = 0;
    std::cout << abi::__cxa_demangle(typeid(decltype(std::forward<T>(tp))).name(), 0, 0, &tmp) << std::endl;
}

void play()
{
    Eigen::Matrix<double, 12, 12> mat = MatrixXd::Zero(12, 12);
    mat.block(0, 0, 3, 3) = -MatrixXd::Identity(3, 3);
    Eigen::Matrix<double, 3, 3> tmp, a;
    tmp << 1, 1, 1,
        1, 1, 1,
        1, 1, 1;
    a = -MatrixXd::Identity(3, 3);
    std::cout << "---mat---" << std::endl;
    std::cout << tmp * a << std::endl;
    Eigen::Matrix3d tr = Eigen::Matrix3d::Identity();
    std::cout << tr << std::endl;
    std::cout << tr.array().pow(0) << std::endl;

    int32_t horizon_length = 300;
    static constexpr double T = 0.01;              // サンプリング周期 (s)
    static constexpr int_fast64_t start_step = 30; // 30 * T = 0.3秒後に歩き出す
    static constexpr int_fast64_t cycle_step = 40; // 何サイクル毎に足踏みするか、 cycle_step * T = 周期(s)
    static constexpr double step_length = 0.3;     // (m)
    Eigen::VectorXd ret = Eigen::VectorXd::Zero(horizon_length);
    for (int32_t i = 0; i < horizon_length; i++)
    {
        if ((i) > start_step)
        {
            if ( ((i) / cycle_step) % 2 == 0)
            {
                ret(i, 0) = step_length;
            }
            else
            {
                ret(i, 0) = -step_length;
            }
        }
        // ret(i, 0) = std::floor((double)(i + cycle_step) / (cycle_step + cycle_step)) - std::floor((double)(i) / (cycle_step + cycle_step)) * 2.0f;
    }
    std::cout << ret.transpose() << std::endl;
}
