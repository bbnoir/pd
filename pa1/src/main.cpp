#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include "partitioner.h"
using namespace std;

int main(int argc, char** argv)
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

    auto end = chrono::high_resolution_clock::now();
    cout << "Initialization: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;

    Partitioner* partitioner = new Partitioner(input);
    end = chrono::high_resolution_clock::now();
    cout << "Parsing input: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    partitioner->partition();
    end = chrono::high_resolution_clock::now();
    cout << "Partitioning: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    partitioner->printSummary();
    end = chrono::high_resolution_clock::now();
    cout << "Printing summary: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    partitioner->writeResult(output);
    end = chrono::high_resolution_clock::now();
    cout << "Writing result: "
         << chrono::duration_cast<chrono::milliseconds>(end - start).count()
         << " ms" << endl;

    return 0;
}
