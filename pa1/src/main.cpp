#include <iostream>
#include <fstream>
#include <vector>
#include "partitioner.h"
#include <omp.h>

using namespace std;

int main(int argc, char **argv)
{
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

    vector<int> seedList = {8457, 2196, 1672, 42111, 49081, 99295};
    Partitioner *partitioner = new Partitioner(input);
    input.close();
    vector<Partitioner*> partitioners(6, nullptr);

    #pragma omp parallel num_threads(6)
    {
        int i = omp_get_thread_num();
        Partitioner *parti = new Partitioner(*partitioner);
        partitioners[i] = parti;
        parti->initPartition(seedList[i]);
        parti->partition();
    }

    int minCutSize = partitioners[0]->getCutSize();
    int minIndex = 0;
    for (int i = 1; i < 6; ++i) {
        if (partitioners[i]->getCutSize() < minCutSize) {
            minCutSize = partitioners[i]->getCutSize();
            minIndex = i;
        }
    }
    partitioners[minIndex]->writeResult(output);

    delete partitioner;
    for (Partitioner* parti : partitioners) {
        delete parti;
    }

    return 0;
}
