#include <iostream>
#include <cmath>
#include <vector>

#include "gnuplot-iostream.h"

int main(int argc, char** argv)
{

    Gnuplot gp;

    double x = 0.;
    double y = 0.8;
    double y_1 = 2;
    double h = 0.1;

    std::vector<std::pair<double, double> > coords;

    for (x = 0.; x <= 1.; x += h)
    {
        double ynext;
        double y_1next;
        y_1next = h*(cos(3*x) - 4*y) + y_1;
        ynext = y + h*y_1;
        coords.push_back(std::make_pair(x, y));
        y = ynext;
        y_1 = y_1next;
    }

    gp << "plot '-' with lines title 'numerical solution',"
        <<" (6*cos(2*x))/5 + sin(2*x) + sin(2*x)*(sin(5*x)/20 + sin(x)/4) - (2*cos(2*x)*(10*tan(x/2)**2 - 20*tan(x/2)**4 + 30*tan(x/2)**6 - 5*tan(x/2)**8 + 1))/(5*(tan(x/2)**2 + 1)**5) title 'analytical solution'"
        << std::endl;
    gp.send1d(coords);
    std::getchar();
    return 0;
}
