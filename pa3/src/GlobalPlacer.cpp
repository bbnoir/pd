#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

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

    const double kAlpha = 0.3;
    ObjectiveFunction objFunc(_placement, positions);
    SimpleConjugateGradient optimizer(objFunc, positions, kAlpha);
    optimizer.Initialize();
    double prev_overflow_ratio = 1;
    for (size_t i = 0, end_i = 1000; i < end_i; ++i) {
        if (i == 30) { objFunc.set_init_lambda(); }
        double f = objFunc(positions), wl = objFunc.wl_value(), d = objFunc.d_value(), overflow_ratio = objFunc.overflow_ratio();
        printf("iter %6ld: f = %10.2e, lambda = %10.2e, wl = %10.2e (%3.1f%%), d = %10.2e (%3.1f%%), overflow_ratio = %3.2f\n",
               i, f, objFunc.lambda(), wl, wl / f * 100, d, d * objFunc.lambda() / f * 100, overflow_ratio);
        optimizer.Step();
        for (size_t j = 0, end_j = num_modules; j < end_j; ++j) {
            Module &module = _placement.module(j);
            if (module.isFixed()) {
                // set back pos of fixed modules
                positions[j].x = module.x();
                positions[j].y = module.y();
            } else {
                // check if the module is out of the bounding box
                if (positions[j].x < _placement.boundryLeft()) {
                    positions[j].x = _placement.boundryLeft();
                } else if (positions[j].x + module.width() > _placement.boundryRight()) {
                    positions[j].x = _placement.boundryRight() - module.width();
                }
                if (positions[j].y < _placement.boundryBottom()) {
                    positions[j].y = _placement.boundryBottom();
                } else if (positions[j].y + module.height() > _placement.boundryTop()) {
                    positions[j].y = _placement.boundryTop() - module.height();
                }
            }
        }
        if (i > 30 && prev_overflow_ratio < overflow_ratio) {
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
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd)) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}
