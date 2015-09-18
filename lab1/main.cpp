#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cmath>

using namespace std;

class matrix
{
    double *array;
public:
    const int dim;
    matrix(int n): dim(n)
    {
        array = new double[n*n];
    }
    double& operator()(int i, int j)
    {
        return array[i*dim + j];
    }
    ~matrix()
    {
        delete[] array;
    }
};

void cholesky(matrix &in, matrix &out)
{
    for (int i = 0; i < in.dim; i++)
        for (int j = 0; j < in.dim; j++)
            out(i, j) = 0.;
    
    for (int j = 0; j < in.dim; j++)
    {
        out(j, j) = sqrt(in(j, j));
        for (int k = j + 1; k < in.dim; k++)
            out(k, j) = in(k, j) / out(j, j);

        for (int k = j + 1; k < in.dim; k++)
            for (int i = k; i < in.dim; i++)
                in(i, k) = in(i, k) - out(i, j)*out(k, j);
    }
}

std::vector<double> 
solve_lower_triangle_linear_system(matrix &A, const std::vector<double> &b)
{
    std::vector<double> x(A.dim, 0);
    for (int i = 0; i < A.dim; i++)
    {
        double b_i = b[i];
        for (int j = 0; j < i; j++)
            b_i -= A(i, j) * x[j];
        x[i] = b_i / A(i, i);
    }
    return x;
}

std::vector<double> 
solve_upper_triangle_linear_system(matrix &A, const std::vector<double> &b)
{
    std::vector<double> x(A.dim, 0);
    for (int i = A.dim - 1; i >= 0; i--)
    {
        double b_i = b[i];
        for (int j = i + 1; j < A.dim; j++)
            b_i -= A(i, j) * x[j];
        x[i] = b_i / A(i, i);
    }
    return x;
}

void print_usage()
{
    cout << "Purpose: solving system of linear equations using Cholesky decomposition" << endl
        << "Usage: lab1 -i input.txt -o output.txt" << endl
        << "Options:" << endl
        << " --input, -i input.txt    set input file name" << endl
        << " --output, -o output.txt  set output file name" << endl
        << " --help, -h               print this help" << endl
        << endl 
        << "Input file format:" << endl
        << "First line:   1 number N(dimension)" << endl
        << "Next N lines: N numbers(Matrix)" << endl
        << "Next N line:  1 number(coefficient vector)" << endl;
}

int main(int argc, char **argv)
{
    std::string input_file_name = "input.txt";
    std::string output_file_name = "output.txt";
    struct option options[] = {
        {"input", required_argument, NULL, 'i'},
        {"output", required_argument, NULL, 'o'},
        {"help", no_argument, NULL, 'h'},
        {0,0,0,0}
    };
    int opt;
    while ( (opt = getopt_long(argc, argv, "i:o:h", options, NULL)) != -1 )
    {
        if (opt == 'i')
            input_file_name = optarg;
        else if (opt == 'o')
            output_file_name = optarg;
        else 
        {
            print_usage();
            return 0;
        }
    }

    ifstream in(input_file_name.c_str());
    ofstream out(output_file_name.c_str());

    int dimension;
    in >> dimension;

    matrix A(dimension), U(dimension);
    std::vector<double> b(dimension);
    for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
            in >> A(i, j);

    for (int i = 0; i < dimension; i++)
        in >> b[i];

    for (int i = 0; i < dimension; i++)
        for (int j = i; j < dimension; j++)
            if (A(i, j) != A(j, i))
            {
                cout << "Input matrix is not symetric" << endl;
                return 0;
            }

    cholesky(A, U);
/*
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
            std::cout << U(i, j) << ' ';
        std::cout << std::endl;
    }
*/

    vector<double> b2 = solve_lower_triangle_linear_system(U, b);
/*
    std::cout << "------" << std::endl;
    for (int i = 0; i < dimension; i++)
        std::cout << b2[i] << std::endl;
*/
    for (int i = 0; i < dimension; i++)
        for (int j = 0; j <= i; j++)
            std::swap(U(i, j), U(j, i));
/*
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
            std::cout << U(i, j) << ' ';
        std::cout << std::endl;
    }
*/
    vector<double> x = solve_upper_triangle_linear_system(U, b2);

    //std::cout << "------" << std::endl;
    for (int i = 0; i < dimension; i++)
        out << x[i] << std::endl;

    return 0;
}
