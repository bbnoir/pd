#include <iostream>
#include <fstream>
#include "fp.h"

#include <string>

using namespace std;

int main(int argv, char** argc) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    if (argv < 5) {
        cout << "[ERROR] Usage: " << argc[0] << " <alpha> <input.block> <input.nets> <output>" << endl;
        return EXIT_FAILURE;
    }

    double alpha = stod(argc[1]);
    if (alpha < 0 || alpha > 1) {
        cout << "[ERROR] alpha must be in [0, 1] but got " << alpha << endl;
        return EXIT_FAILURE;
    }
    ifstream inFileBlock(argc[2]);
    if (!inFileBlock.is_open()) {
        cout << "[ERROR] cannot open file " << argc[2] << endl;
        return EXIT_FAILURE;
    }
    ifstream inFileNet(argc[3]);
    if (!inFileNet.is_open()) {
        cout << "[ERROR] cannot open file " << argc[3] << endl;
        return EXIT_FAILURE;
    }
    ofstream outFile(argc[4]);
    if (!outFile.is_open()) {
        cout << "[ERROR] cannot open file " << argc[4] << endl;
        return EXIT_FAILURE;
    }

    FloorPlanner *fp = new FloorPlanner(alpha);
    fp->readInput(inFileBlock, inFileNet);
    inFileBlock.close();
    inFileNet.close();
    fp->floorplanParallel();
    fp->writeOutput(outFile);
    outFile.close();
    delete fp;

    return EXIT_SUCCESS;
}