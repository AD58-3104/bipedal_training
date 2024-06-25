#ifndef  PLOT_HPP_
#define PLOT_HPP_
#include <sciplot/sciplot.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include "robot_model.hpp"

void plotArm(const ResultPosition data)
{
    std::ofstream ofs("pltarm.m");
    ofs << "x = [";
    for (auto &x : data.x)
    {
        ofs << x << " ";
    }
    ofs << "];" << std::endl;
    ofs << "y = [";
    for (auto &y : data.y)
    {
        ofs << y << " ";
    }
    ofs << "];" << std::endl;
    ofs << "z = [";
    for (auto &z : data.z)
    {
        ofs << z << " ";
    }
    ofs << "];" << std::endl;
    ofs
        << std::endl
        << "grid on;"
        << std::endl
        << "plot3(x,y,z,'-o');"
        << std::endl

        << "title('My Plot');"
        << std::endl
        << "xlabel('x');"
        << std::endl
        << "ylabel('y');"
        << std::endl
        << "zlabel('z');"
        << std::endl
        << "xlim([-100 100]);"
        << std::endl
        << "% ylim([-5 10]);"
        // << std::endl << "plot(y,z,'-o');"
        // << std::endl << "ylabel('y');zlabel('z');"
        ;
    ofs.close();
    std::system("octave --persist pltarm.m > /dev/null 2>&1");

    return;
}

void plot3darm(const ResultPosition& result)
{
    using namespace sciplot;
    static Plot3D plot3d;
    plot3d.autoclean(false);

    plot3d.legend().hide();
    plot3d.xlabel("x");
    plot3d.ylabel("y");
    plot3d.zlabel("z");
    plot3d.xrange(-300.f, 300.f);
    plot3d.yrange(-200.f, 200.f);
    plot3d.zrange(-300.f, 30.f);
    plot3d.drawWithVecs("linespoints", result.x, result.y, result.z).lineColor("red");
    Figure fig3d = {{plot3d}};
    Canvas canvas3d = {{fig3d}};
    canvas3d.size(500, 500);
    canvas3d.show();
    return;
}

/**
 * @brief ロボットを側面から描画する
 *
 * @param result
 */
void plotXZ(const ResultPosition& result)
{
    using namespace sciplot;
    static Plot2D plot2d;
    plot2d.legend().hide();
    plot2d.xlabel("x");
    plot2d.ylabel("z");
    plot2d.xrange(-500.f, 500.f);
    plot2d.yrange(-500.f, 30.f);
    plot2d.drawWithVecs("linespoints", result.x, result.z).lineColor("blue");
    Figure fig2d = {{plot2d}};
    fig2d.title("XZ Plane");
    Canvas canvas2d = {{fig2d}};
    canvas2d.show();
    return;
}

/**
 * @brief ロボットを正面から描画する
 *
 * @param result
 */
void plotYZ(const ResultPosition& result)
{
    using namespace sciplot;
    static Plot2D plot2d;
    plot2d.legend().hide();
    plot2d.xlabel("y");
    plot2d.ylabel("z");
    plot2d.xrange(-500.f, 500.f);
    plot2d.yrange(-500.f, 30.f);
    plot2d.drawWithVecs("linespoints", result.y, result.z).lineColor("green");
    Figure fig2d = {{plot2d}};
    fig2d.title("YZ Plane");
    Canvas canvas2d = {{fig2d}};
    canvas2d.show();
    return;
}


#endif // PLOT_HPP_