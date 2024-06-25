#include "plot.hpp"
#include "robot_model.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char const *argv[])
{
    std::system("killall gnuplot_qt");
    RobotModel robot;
    // y,  r,        p,p, p
    // std::vector<double> degs = {0.0_deg,6.3_deg,0.0_deg, 17.8_deg, -17.8_deg, 1.0_deg, -6.3_deg}; //right
    std::vector<double> degs = {0.0_deg,-6.3_deg,-23.0_deg, 17.8_deg, -17.8_deg, 1.0_deg, 6.3_deg}; //left
    // 膝上ピッチの回転が逆、足首の回転が逆
    robot.setAllJointAngle(degs);
    auto result = robot.calcForwardKinematics(true);
    plot3darm(result);
    plotXZ(result);
    plotYZ(result);
    robot.printAllJointAnglesDeg();
    return 0;
}