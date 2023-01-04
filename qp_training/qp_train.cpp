#include "qp_train.hpp"
#include <unistd.h>

int main(int argc, char const *argv[])
{
    solve_qp();
    system("sleep 1");
    // system("gnuplot-x11 -persist show.plt");
    return 0;
}
