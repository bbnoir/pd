#pragma once

// constants
constexpr int INF = __INT32_MAX__;
constexpr double initP = 0.8;
constexpr double initP2 = 0.4;
constexpr double R = 0.85; // for cooling
constexpr double TIME_LIMIT = 240.0; // seconds
constexpr double ALPHA_MIN = 0.5;

constexpr int RANDOM_START = 1;
constexpr int GREEDY_CONVERGE = 2;
constexpr int ANNEALING = 3;
constexpr int RE_ANNEALING = 4;

constexpr double RANDOM_START_TIMES = 1.0; // * numBlocks

// variables
extern int numThreads;