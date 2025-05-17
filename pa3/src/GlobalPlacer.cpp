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

static void revert_illegal(vector<Point2<double>> &positions, Placement &placement) {
    // clamp positions
    const double left = placement.boundryLeft();
    const double right = placement.boundryRight();
    const double bottom = placement.boundryBottom();
    const double top = placement.boundryTop();
    for (size_t i = 0, end_i = positions.size(); i < end_i; ++i) {
        Module &module = placement.module(i);
        if (module.isFixed()) {
            continue;
        }
        if (positions[i].x < left) {
            positions[i].x = left;
        } else if (positions[i].x + module.width() > right) {
            positions[i].x = right - module.width();
        }
        if (positions[i].y < bottom) {
            positions[i].y = bottom;
        } else if (positions[i].y + module.height() > top) {
            positions[i].y = top - module.height();
        }
    }
    // revert fixed modules
    for (size_t i = 0, end_i = positions.size(); i < end_i; ++i) {
        Module &module = placement.module(i);
        if (module.isFixed()) {
            positions[i].x = module.x();
            positions[i].y = module.y();
        }
    }
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
    const double random_width = outline_width * 0.005;
    const double random_height = outline_height * 0.005;
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

    const double kAlpha = 0.1;
    ObjectiveFunction objFunc(_placement, positions);
    SimpleConjugateGradient optimizer(objFunc, positions, kAlpha);
    optimizer.Initialize();

    const int stage_limit = 100;

    // optimize wirelength
    int wl_iter = 0;
    double best_wl = objFunc(positions); // initial lambda = 0
    int wl_stop_cnt = 0;
    const int max_wl_stop_cnt = 3;
    while (wl_iter < stage_limit && wl_stop_cnt < max_wl_stop_cnt) {
        ++wl_iter;
        optimizer.Step();
        revert_illegal(positions, _placement);
        double wl = objFunc(positions);
        if (wl < best_wl) {
            best_wl = wl;
            wl_stop_cnt = 0;
        } else {
            ++wl_stop_cnt;
        }
    }

    vector<int> bin_size_list = {8, 16, 32};
    for (int i = 0; i < 3; ++i) {
        // optimize density
        objFunc.resize_bin(bin_size_list[i]);
        objFunc.set_init_lambda();
        objFunc(positions);
        double prev_overflow_ratio = objFunc.overflow_ratio();
        int d_iter = 0;
        while (d_iter < stage_limit) {
            ++d_iter;
            optimizer.Step();
            revert_illegal(positions, _placement);
            const double overflow_ratio = objFunc.overflow_ratio();
            if (d_iter > 10 && overflow_ratio >= prev_overflow_ratio) {
                objFunc.scale_lambda(1.5);
            }
            if (overflow_ratio < 0.1) { break; }
            prev_overflow_ratio = overflow_ratio;
        }
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
