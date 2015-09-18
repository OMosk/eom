#include <iostream>
#include <getopt.h>
#include <cmath>

using namespace std;

double f(double x)
{
    return std::pow(2, 3*x);
}

double F(double x)
{
    return std::pow(2, 3*x) / (3*std::log(2));
}

void print_usage()
{
    cout << "Calculate integral of y = 2 ^ (3x) in interval [l; r]" << endl << endl
        << "Usage: lab2 -i input.txt -o output.txt" << endl
        << "Options:" << endl
        << " --integral, -i     perform integration" << endl
        << " --left, -l number  set left integration border" << endl
        << " --right, -r number set right integration border" << endl
        << " --step, -s number  set integration step" << endl
        << endl
        << " --derivative -d    perform differentiation" << endl
        << " --step, -s number  derivative scheme step " << endl
        << " --point -x number  set point" << endl
        << endl
        << " --help, -h         print this help" << endl
        << endl;
}

int main(int argc, char **argv)
{
    double left = 0.;
    double right = 0.;
    double step = 0.1;
    double x = 0.;

    bool find_integral = false;
    bool find_derivative = false;

    struct option options[] = {
        {"left", required_argument, NULL, 'l'},
        {"right", required_argument, NULL, 'r'},
        {"step", required_argument, NULL, 's'},
        {"integral", no_argument, NULL, 'i'},
        {"derivative", no_argument, NULL, 'd'},
        {"point", required_argument, NULL, 'x'},
        {"help", no_argument, NULL, 'h'},
        {0,0,0,0}
    };

    int opt;
    while ( (opt = getopt_long(argc, argv, "l:r:s:hidx:", options, NULL)) != -1 )
    {
        if (opt == 'l')
            left = std::atof(optarg);
        else if (opt == 'r')
            right = std::atof(optarg);
        else if (opt == 's')
            step = std::atof(optarg);
        else if (opt == 'i')
            find_integral = true;
        else if (opt == 'd')
            find_derivative = true;
        else if (opt == 'x')
            x = std::atof(optarg);
        else 
        {
            print_usage();
            return 0;
        }
    }

    if (find_integral)
    {
        double integral = 0.;
        for (double i = left + step; i <= right; i+=step)
        {
            integral += step * (f(i - step) + f(i)) / 2.;
        }
        std::cout << "Integral computation" << std::endl;
        std::cout << "Trapezoidal rule:       " << integral << std::endl;
        std::cout << "Newton-Leupzig formula: " << F(right) - F(left) << std::endl;
        std::cout << std::endl;
    }

    if (find_derivative)
    {
        double derivative = (F(x + step/2.) - F(x - step/2.))/step;
        std::cout << "Derivative computation" << std::endl;
        std::cout << "Central derivative scheme:   " << derivative << endl;
        std::cout << "Analytical derivative value: " << f(x) << endl;
        std::cout << std::endl;
    }
    return 0;
}
