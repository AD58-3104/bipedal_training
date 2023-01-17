#include <Eigen/Dense>
#include <Eigen/Sparse>
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

template <typename T>
void sparseBlockAssignation(Eigen::SparseMatrix<T>& sparse_mat, const size_t &row_location, const size_t &col_location,const Eigen::MatrixXd& assign_mat)
{
    const size_t assign_row = assign_mat.rows();
    const size_t assign_col = assign_mat.cols();
    const size_t sparse_max_col = sparse_mat.cols();
    const size_t sparse_max_row = sparse_mat.rows();
    if ((sparse_max_col < assign_col + col_location) || (sparse_max_row < assign_row + row_location))
    {
        throw std::range_error("write block will exceed sparse size!!!");
    }
    for (size_t row = row_location; row < row_location + assign_row; row++)
    {
        for (size_t col = col_location; col < col_location + assign_col; col++)
        {
            assert((sparse_max_col >= assign_col + col_location) || (sparse_max_row >= assign_row + row_location)); // over sparse size
            sparse_mat.insert(row, col) = assign_mat(row - row_location, col - col_location);
        }
    }
}

void play()
{
    Eigen::Matrix<double, 12, 12> mat = MatrixXd::Zero(12, 12);
    mat.block(0, 0, 3, 3) = -MatrixXd::Identity(3, 3);
    Eigen::Matrix<double, 3, 3> tmp, a;
    tmp << 1, 1, 1,
        1, 1, 1,
        1, 1, 1;
    a = 7 * MatrixXd::Constant(3, 3,1);
    std::cout << "---mat---" << std::endl;
    std::cout << tmp * a << std::endl;
    Eigen::Matrix3d tr = Eigen::Matrix3d::Identity();
    std::cout << tr << std::endl;
    std::cout << tr.array().pow(0) << std::endl;
    std::cout << " ------ sparse mat   ------" << std::endl;
    Eigen::SparseMatrix<double> sparse;
    Eigen::Vector3d vec;
    vec << 1,2,3;
    sparse.resize(40,40);
    sparseBlockAssignation(sparse,2,1,a);
    sparseBlockAssignation(sparse,36,39,vec);
    sparseBlockAssignation(sparse,39,37,vec.transpose());
    std::cout << sparse << std::endl;
}
