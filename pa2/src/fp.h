#pragma once
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <random>
#include "config.h"

class BlockLib;
class BlockInst;
class NetLib;
class NetInst;
class ContourNode;
class TreeNode;

struct Solution {
    int cost = INF;
    double normCost = INF;
    int wireLength = INF;
    int area = INF;
    int width = INF;
    int height = INF;
    int iter = 0;
    std::vector<std::pair<int, int>> blockLowerLeft;
    std::vector<std::pair<int, int>> blockUpperRight;

    Solution() {}
    Solution(int nBlocks) : blockLowerLeft(nBlocks), blockUpperRight(nBlocks) {}

    bool isInOutline(int outlineWidth, int outlineHeight) const {
        return width <= outlineWidth && height <= outlineHeight;
    }
};

struct SAParam {
    unsigned seed;
    int maxIter = 200000;
    bool shuffleInit = true;
    double convergePow = 1.2;

    SAParam() {}
};

class FloorPlanner {
public:
    FloorPlanner(double alpha);
    ~FloorPlanner();
    void readInput(std::ifstream& inFileBlock, std::ifstream& inFileNet);
    void floorplanParallel();
    void writeOutput(std::ofstream& outFile);

private:

    TreeNode* genDefaultTree(std::vector<TreeNode*>& treeNodes) const;
    Solution* sa(const SAParam& param) const;

    TreeNode* removeTreeNode(TreeNode*& root, TreeNode* node, std::vector<int>& touchedNodeId) const;
    void insertTreeNode(TreeNode* root, TreeNode* node, TreeNode* parent, std::vector<int>& touchedNodeId, bool left) const;

    void insertBlockOnContour(BlockInst* block, ContourNode* contourHead) const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> _startTime;
    double getElapsedTime() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - _startTime).count() / 1000.0;
    }

    // sa parameters
    double _alpha = 0;
    double _energyAlpha = 0;
    double _targetAR = 0;
    int _normArea = 0;
    int _normWL = 0;

    // circuit information
    int _outlineWidth = 0;
    int _outlineHeight = 0;
    int _nBlocks = 0;
    int _nNets = 0;

    // library and default data structures
    std::vector<BlockLib*> _blockLibs;
    std::vector<NetLib*> _netLibs;

    Solution* _bestSolution = nullptr;
};