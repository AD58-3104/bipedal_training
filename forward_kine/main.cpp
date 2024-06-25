#include "plot.hpp"
#include "robot_model.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char const *argv[])
{
    std::system("killall gnuplot_qt");
    RobotModel robot;
    // y,  r,        p,p, p
    std::vector<double> degs = {15.0_deg, 0, 15.0_deg, 0, 15.0_deg, 0, 15.0_deg};
    // 膝上ピッチの回転が逆、足首の回転が逆
    robot.setAllJointAngle(degs);
    auto result = robot.calcForwardKinematics(true);
    plot3darm(result);
    plotXZ(result);
    plotYZ(result);
    std::cout << result.getEndPositionVec() << std::endl;
    return 0;
}