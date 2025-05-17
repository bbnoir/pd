#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

using namespace std;

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {
    ////////////////////////////////////////////////////////////////////
    // Global placement algorithm
    ////////////////////////////////////////////////////////////////////

    // Randomly initialize the positions of modules
    const size_t num_modules = _placement.numModules();
    std::vector<Point2<double>> positions(num_modules);
    const double outline_width = _placement.boundryRight() - _placement.boundryLeft();
    const double outline_height = _placement.boundryTop() - _placement.boundryBottom();
    const double mid_x = _placement.boundryLeft() + outline_width / 2;
    const double mid_y = _placement.boundryBottom() + outline_height / 2;
    const double random_width = outline_width * 0.1;
    const double random_height = outline_height * 0.1;
    srand(42);
    for (size_t i = 0, end_i = num_modules; i < end_i; ++i) {
        Module &module = _placement.module(i);
        if (module.isFixed()) {
            positions[i].x = module.x();
            positions[i].y = module.y();
        } else {
            positions[i].x = mid_x + double(rand()) / RAND_MAX * random_width - random_width / 2;
            positions[i].y = mid_y + double(rand()) / RAND_MAX * random_height - random_height / 2;
        }
    }

    const double kAlpha = 0.5;
    ObjectiveFunction objFunc(_placement, positions);
    SimpleConjugateGradient optimizer(objFunc, positions, kAlpha);
    optimizer.Initialize();
    double prev_overflow_ratio = 1;
    int overflow_stop_cnt = 0;
    for (size_t i = 0, end_i = 1000; i < end_i; ++i) {
        if (i == 100) { objFunc.set_init_lambda(); }
        double f = objFunc(positions), wl = objFunc.wl_value(), d = objFunc.d_value(), overflow_ratio = objFunc.overflow_ratio();
        printf("iter %6ld: f = %10.2e, alpha = %10.2e, lambda = %10.2e, wl = %10.2e (%3.1f%%), d = %10.2e (%3.1f%%), overflow_ratio = %3.2f\n",
               i, f, optimizer.alpha(), objFunc.lambda(), wl, wl / f * 100, d, d * objFunc.lambda() / f * 100, overflow_ratio);
        optimizer.Step();
        for (size_t j = 0, end_j = num_modules; j < end_j; ++j) {
            Module &module = _placement.module(j);
            if (module.isFixed()) {
                // set back pos of fixed modules
                positions[j].x = module.x();
                positions[j].y = module.y();
            } else {
                // check if the module is out of the bounding box
                const double half_bin_width = objFunc.bin_width() * 0.5;
                if (positions[j].x < _placement.boundryLeft()) {
                    positions[j].x = _placement.boundryLeft() + half_bin_width;
                } else if (positions[j].x + module.width() > _placement.boundryRight()) {
                    positions[j].x = _placement.boundryRight() - module.width() - half_bin_width;
                }
                if (positions[j].y < _placement.boundryBottom()) {
                    positions[j].y = _placement.boundryBottom() + half_bin_width;
                } else if (positions[j].y + module.height() > _placement.boundryTop()) {
                    positions[j].y = _placement.boundryTop() - module.height() - half_bin_width;
                }
            }
        }
        if (prev_overflow_ratio - overflow_ratio < 0.01) {
            ++overflow_stop_cnt;
        } else {
            overflow_stop_cnt = 0;
        }
        if (d * objFunc.lambda() / f > 0.8) {
            objFunc.scale_lambda(0.1);
        }
        if (i > 30 && overflow_stop_cnt > 1) {
            objFunc.scale_lambda(2);
        }
        prev_overflow_ratio = overflow_ratio;
    }

    ////////////////////////////////////////////////////////////////////
    // Write the placement result into the database. (You may modify this part.)
    ////////////////////////////////////////////////////////////////////
    size_t fixed_cnt = 0;
    for (size_t i = 0; i < num_modules; i++) {
        // If the module is fixed, its position should not be changed.
        // In this programing assignment, a fixed module may be a terminal or a pre-placed module.
        if (_placement.module(i).isFixed()) {
            fixed_cnt++;
            continue;
        }
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    printf("INFO: %lu / %lu modules are fixed.\n", fixed_cnt, num_modules);
}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt) {
    // Generate data files for matplotlib
    string boundaryFileName = outfilename + "_boundary.dat";
    string moduleFileName = outfilename + "_modules.dat";
    string infoFileName = outfilename + "_info.dat";
    
    // Write boundary data
    ofstream boundaryFile(boundaryFileName.c_str(), ios::out);
    plotBoxDAT(boundaryFile, _placement.boundryLeft(), _placement.boundryBottom(), 
               _placement.boundryRight(), _placement.boundryTop());
    boundaryFile.close();
    
    // Write module data
    ofstream moduleFile(moduleFileName.c_str(), ios::out);
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxDAT(moduleFile, module.x(), module.y(), 
                  module.x() + module.width(), module.y() + module.height());
        // Add blank line between modules for matplotlib to recognize separate polygons
        moduleFile << endl;
    }
    moduleFile.close();
    
    // Write placement info
    ofstream infoFile(infoFileName.c_str(), ios::out);
    infoFile << _placement.computeHpwl() << endl;
    infoFile.close();
    
    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "python plot.py %s %s %s", 
                boundaryFileName.c_str(), moduleFileName.c_str(), infoFileName.c_str());
        if (system(cmd) != 0) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

// Helper function for writing box coordinates to data files
void GlobalPlacer::plotBoxDAT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << " " << y1 << endl
           << x2 << " " << y1 << endl
           << x2 << " " << y2 << endl
           << x1 << " " << y2 << endl
           << x1 << " " << y1 << endl;
}
