#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include "partitioner.h"
using namespace std;

static void printRuntime(const string &msg, const chrono::high_resolution_clock::time_point &start)
{
    auto end = chrono::high_resolution_clock::now();
    cout << msg << ": "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
}

int main(int argc, char **argv)
{
    auto start = chrono::high_resolution_clock::now();
    fstream input, output;

    if (argc == 3) {
        input.open(argv[1], ios::in);
        output.open(argv[2], ios::out);
        if (!input) {
            cerr << "Cannot open the input file \"" << argv[1]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./fm <input file> <output file>" << endl;
        exit(1);
    }

    Partitioner *partitioner = new Partitioner(input);
    printRuntime("Parsing input", start);

    partitioner->initPartition();
    printRuntime("Initializing partition", start);

    partitioner->partition();
    printRuntime("Partitioning", start);

    partitioner->printSummary();
    printRuntime("Printing summary", start);

    partitioner->writeResult(output);
    printRuntime("Writing result", start);

    return 0;
}
