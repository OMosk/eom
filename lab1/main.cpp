#include <iostream>
#include <fstream>
#include <string>
#include <array>
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

void choleckiy(matrix &in, matrix &out)
{
    for (int i = 0; i < in.dim; i++)
        for (int j = 0; j < in.dim; j++)
            out(i, j) = 0.;
    
    for (int j = 0; j < in.dim - 1; j++)
    {
        out(j, j) = sqrt(in(j, j));
        for (int k = j; k < in.dim; k++)
            out(k, j) = in(k, j) / out(j, j);

        for (int k = j + 1; k < in.dim; k++)
            for (int i = k; i < in.dim; i++)
                in(i, k) = in(i, k) - out(i, j)*out(k, j);
    }
}

void print_usage()
{
    cout << "Usage: lab1 -i input.txt -o output.txt" << endl
        << "Options:" << endl
        << " --input, -i input.txt    set input file name" << endl
        << " --output, -o output.txt  set output file name" << endl
        << " --help, -h               print this help" << endl;
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
    for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
            in >> A(i, j);

    for (int i = 0; i < dimension; i++)
        for (int j = i; j < dimension; j++)
            if (A(i, j) != A(j, i))
            {
                cout << "Input matrix is not symetric" << endl;
                return 0;
            }

    choleckiy(A, U);
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
            out << U(i, j) << ' ';
        out << endl;
    }


    return 0;
}
