#include "train.hpp"

int main(int argc, char const *argv[])
{
    Eigen::MatrixXf m(4, 4);
    m << 1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16;
    std::cout << m << std::endl;
    auto zero = MatrixXf::Zero(3,3);
    auto cons = MatrixXf::Constant(3,3,5);
    m.block(0,0,3,3) = cons;
    std::cout << "--- change ---" << std::endl << m << std::endl;
    play();
    return 0;
}