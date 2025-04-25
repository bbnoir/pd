#include "fp.h"
#include "block.h"
#include "net.h"
#include "treenode.h"
#include "contour.h"
#include "config.h"
#include <cassert>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <omp.h>

using namespace std;

FloorPlanner::FloorPlanner(double alpha) : _startTime(chrono::high_resolution_clock::now()), _alpha(alpha)
{
}

FloorPlanner::~FloorPlanner()
{
    for (auto block : _blockLibs)
    {
        delete block;
    }
    for (auto net : _netLibs)
    {
        delete net;
    }
    if (_bestSolution != nullptr)
    {
        delete _bestSolution;
    }
}

void FloorPlanner::readInput(ifstream &inFileBlock, ifstream &inFileNet)
{
    // parse blocks
    string s;
    inFileBlock >> s;
    assert(s == "Outline:");
    inFileBlock >> _outlineWidth >> _outlineHeight;
    inFileBlock >> s;
    assert(s == "NumBlocks:");
    inFileBlock >> _nBlocks;
    inFileBlock >> s;
    assert(s == "NumTerminals:");
    int nTerminals;
    inFileBlock >> nTerminals;
    unordered_map<string, int> blockMap;
    unordered_map<string, int> terminalMap;
    vector<pair<int, int>> terminals;
    string name;
    int w, h;
    for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
    {
        inFileBlock >> name >> w >> h;
        blockMap[name] = i;
        _blockLibs.push_back(new BlockLib(name, w, h));
    }
    int x, y;
    for (int i = 0, end_i = nTerminals; i < end_i; ++i)
    {
        inFileBlock >> name >> s >> x >> y;
        assert(s == "terminal");
        terminalMap[name] = i;
        terminals.push_back(make_pair(x, y));
    }
    // parse nets
    inFileNet >> s;
    assert(s == "NumNets:");
    inFileNet >> _nNets;
    for (int i = 0, end_i = _nNets; i < end_i; ++i)
    {
        inFileNet >> s;
        assert(s == "NetDegree:");
        NetLib *net = new NetLib();
        _netLibs.push_back(net);
        int n;
        inFileNet >> n;
        for (int j = 0; j < n; ++j)
        {
            inFileNet >> s;
            assert(s != "NetDegree:");
            if (blockMap.find(s) != blockMap.end())
            {
                net->addBlockIdx(blockMap[s]);
            }
            else if (terminalMap.find(s) != terminalMap.end())
            {
                net->addTerminal(terminals[terminalMap[s]]);
            }
            #ifdef DEBUG
            else
            {
                cerr << "[ERROR] <FloorPlanner::readInput()> unknown block or terminal: " << s << endl;
                exit(EXIT_FAILURE);
            }
            #endif
        }
    }
}

void FloorPlanner::floorplanParallel()
{
    // initSA();
    // _energyAlpha = max(0.5, _alpha); // handle extereme case
    _energyAlpha = _alpha;
    _targetAR = double(_outlineWidth) / double(_outlineHeight);
    srand(42);
    mt19937 mt(time(0));
    SAParam param[numThreads];
    param[0].seed = rand();
    bool everInOutline = false;
    _bestSolution = sa(param[0]);
    #ifdef LOG
    int try_count = 0;
    int found_count = 0;
    int sum_iter = 0;
    #endif
    for (int i = 1; i < numThreads-1; i++)
    {
        param[i].seed = rand();
    }
    param[numThreads-1].seed = mt();
    vector<bool> foundList(numThreads, false);
    vector<bool> foundBetterList(numThreads, false);
    #pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < numThreads; i++)
    {
        Solution* sol = sa(param[i]);
        if (sol->isInOutline(_outlineWidth, _outlineHeight))
        {
            foundList[i] = true;
            everInOutline = true;
        }
        if (_bestSolution->cost > sol->cost)
        {
            foundBetterList[i] = true;
            #pragma omp critical
            {
                if (_bestSolution->cost > sol->cost)
                {
                    delete _bestSolution;
                    _bestSolution = sol;
                }
            }
        }
        else
        {
            delete sol;
        }
    }
    if (param[0].convergePow < 2)
    {
        for (int i = 0; i < numThreads; i++)
        {
            param[i].convergePow += 0.1;
        }
    }
}

void FloorPlanner::floorplanSeed(const int seed)
{
    _energyAlpha = _alpha;
    _targetAR = double(_outlineWidth) / double(_outlineHeight);
    _mt.seed(seed);
    SAParam param;
    param.seed = seed;
    int try_count = 0;
    Solution* sol = sa(param);
    // while (sol->isInOutline(_outlineWidth, _outlineHeight) == false)
    // {
    //     delete sol;
    //     sol = sa(param);
    //     try_count++;
    // }
    _bestSolution = sol;
}

TreeNode* FloorPlanner::genDefaultTree(vector<TreeNode*>& treeNodes)
{
    for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
    {
        treeNodes.push_back(new TreeNode(i));
    }
    for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
    {
        const int leftChildIndex = i * 2 + 1;
        if (leftChildIndex < _nBlocks)
        {
            treeNodes[i]->setLeftChild(treeNodes[leftChildIndex]);
            treeNodes[leftChildIndex]->setParent(treeNodes[i]);
        }
        const int rightChildIndex = i * 2 + 2;
        if (rightChildIndex < _nBlocks)
        {
            treeNodes[i]->setRightChild(treeNodes[rightChildIndex]);
            treeNodes[rightChildIndex]->setParent(treeNodes[i]);
        }
    }
    return treeNodes[0];
}

Solution *FloorPlanner::sa(const SAParam &param)
{
    // initialize solution
    Solution* bestSolution = new Solution(_nBlocks);
    vector<BlockInst*> blocks;
    for (auto block : _blockLibs)
    {
        blocks.push_back(new BlockInst(block));
    }
    vector<NetInst*> nets;
    for (auto net : _netLibs)
    {
        nets.push_back(new NetInst(net, blocks));
    }
    vector<TreeNode*> treeNodes, copyTreeNodes;
    TreeNode* root = genDefaultTree(treeNodes);
    vector<int> assignBlockIdx(_nBlocks);
    for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
    {
        assignBlockIdx[i] = i;
    }
    if (param.shuffleInit)
    {
        random_shuffle(assignBlockIdx.begin(), assignBlockIdx.end());
    }
    for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
    {
        treeNodes[i]->setBlock(blocks[assignBlockIdx[i]]);
        copyTreeNodes.push_back(new TreeNode(*treeNodes[i]));
    }

    // SA initialization
    double curEnergy = 0;
    double prevEnergy = 0;
    double T = 4000000;
    int iter = 0;
    const int maxIter = param.maxIter;

    enum Stage
    {
        RANDOM_START,
        GREEDY_CONVERGE,
        ANNEALING,
        RE_ANNEALING
    };
    Stage stage = RANDOM_START;
    bool converge = false;
    int failCount = 0;
    const int randomStartTimes = RANDOM_START_TIMES * _nBlocks;
    const int convergeCount = min(_nBlocks*5, 6000);
    double sumDeltaEnergy = 0;
    int upHillCount = 0;
    const int reduceCountThr = 1.5 * _nBlocks;
    const int reduceCountThr2 = 1.5 * _nBlocks;
    int reduceCount = 0;
    const double minT = 0.00001;
    bool hasBetter = false;
    bool loop_end = false;

    // for normalization
    double sumArea = 0, sumWL = 0, normArea = 0, normWL = 0;

    #ifdef LOG
    ofstream logFile("sa.log");
    #endif

    // SA loop
    while (iter < maxIter)
    {
        // stage
        switch (stage)
        {
        case RANDOM_START:
            if (iter > randomStartTimes)
            {
                stage = GREEDY_CONVERGE;
            }
            break;
        case GREEDY_CONVERGE:
            if (converge)
            {
                failCount = 0;
                stage = ANNEALING;
                T = (sumDeltaEnergy / double(upHillCount)) / -log(initP);
            }
            break;
        case ANNEALING:
            reduceCount++;
            if (reduceCount > reduceCountThr)
            {
                if (T > minT)
                {
                    T *= R;
                }
                reduceCount = 0;
            }
            if (converge)
            {
                if (hasBetter)
                {
                    hasBetter = false;
                    failCount = 0;
                    stage = RE_ANNEALING;
                    T = (sumDeltaEnergy / double(upHillCount)) / -log(initP2);
                }
                else
                {
                    loop_end = true;
                }
            }
            break;
        case RE_ANNEALING:
            reduceCount++;
            if (reduceCount > reduceCountThr2)
            {
                if (T > minT)
                {
                    T *= R;
                }
                reduceCount = 0;
            }
            if (converge)
            {
                // if (hasBetter)
                // {
                //     hasBetter = false;
                //     failCount = 0;
                //     T = (sumDeltaEnergy / double(upHillCount)) / -log(initP2);
                // }
                // else
                // {
                    loop_end = true;
                // }
            }
            break;
        default:
            break;
        }

        if (loop_end)
        {
            break;
        }

        // perturb
        const int perturbMode = genRandom(0, 2); // 0: rotate, 1: swap, 2: move
        const int idx1 = genRandom(0, _nBlocks - 1);
        int idx2;
        if (!(stage == RANDOM_START && iter != randomStartTimes))
        {
            for (int i = 0, end_i = _nBlocks; i < end_i; ++i)
            {
                copyTreeNodes[i]->copyInfo(treeNodes[i]);
            }
        }
        TreeNode* const originalRoot = root;
        vector<int> touchedNodeId;
        if (perturbMode == 0) // rotate
        {
            blocks[idx1]->rotate();
        }
        else if (perturbMode == 1) // swap
        {
            idx2 = (idx1 + 1 + genRandom(0, _nBlocks - 2)) % _nBlocks;
            treeNodes[idx1]->swapBlock(treeNodes[idx2]);
        }
        else if (perturbMode == 2) // move
        {
            TreeNode* moveNode = treeNodes[idx1];
            touchedNodeId.push_back(idx1);
            TreeNode* removedNode = removeTreeNode(root, moveNode, touchedNodeId);
            idx2 = (removedNode->getId() + 1 + genRandom(0, _nBlocks - 2)) % _nBlocks;
            insertTreeNode(root, removedNode, treeNodes[idx2], touchedNodeId);
        }

        // compact
        root->getBlock()->setX(0);
        ContourNode* contourHead = new ContourNode(0, 0);
        ContourNode* contourTail = new ContourNode(INF, 0);
        contourHead->setNext(contourTail);
        contourTail->setPrev(contourHead);
        vector<TreeNode*> nodeStack;
        nodeStack.push_back(root);
        while(!nodeStack.empty())
        {
            TreeNode* node = nodeStack.back();
            nodeStack.pop_back();
            BlockInst* const curBlock = node->getBlock();
            insertBlockOnContour(curBlock, contourHead);
            TreeNode* rightChild = node->getRightChild();
            if (rightChild != nullptr)
            {
                assert(rightChild->getBlock() != nullptr);
                rightChild->getBlock()->setX(curBlock->getX());
                nodeStack.push_back(rightChild);
            }
            TreeNode* leftChild = node->getLeftChild();
            if (leftChild != nullptr)
            {
                assert(leftChild->getBlock() != nullptr);
                leftChild->getBlock()->setX(curBlock->getX()+curBlock->getWidth());
                nodeStack.push_back(leftChild);
            }
        }
        int w = 0, h = 0;
        for (auto b : blocks)
        {
            w = max(w, b->getX() + b->getWidth());
            h = max(h, b->getY() + b->getHeight());
        }
        for (ContourNode* curNode = contourHead; curNode != contourTail; curNode = curNode->getNext())
        {
            delete curNode;
        }

        // calculate cost
        const double curArea = w * h;
        double curWireLength = 0;
        for (auto net : nets)
        {
            curWireLength += net->calcHPWL();
        }
        double curCost = double(curArea) * _alpha + double(curWireLength) * (1 - _alpha);
        #ifdef LOG
        logFile << iter << ", " << curCost << ", " << curWireLength << ", " << curArea << ", " << T << "\n";
        #endif
        // save best solution
        bool isInOutline = w <= _outlineWidth && h <= _outlineHeight;
        bool saveBest = false;
        double normCost = INF;
        if (stage == RANDOM_START)
        {
            saveBest = (isInOutline && (curCost < bestSolution->cost));
        }
        else
        {
            normCost = double(curArea) / normArea * _alpha + double(curWireLength) / normWL * (1 - _alpha);
            saveBest = (isInOutline && (normCost < bestSolution->normCost));
        }
        if (saveBest)
        {
            hasBetter = true;
            bestSolution->cost = curCost;
            bestSolution->normCost = normCost;
            bestSolution->wireLength = curWireLength;
            bestSolution->area = curArea;
            bestSolution->width = w;
            bestSolution->height = h;
            for (int i = 0; i < _nBlocks; i++)
            {
                bestSolution->blockLowerLeft[i] = make_pair(blocks[i]->getX(), blocks[i]->getY());
                bestSolution->blockUpperRight[i] = make_pair(blocks[i]->getX() + blocks[i]->getWidth(), blocks[i]->getY() + blocks[i]->getHeight());
            }
        }

        if (stage == RANDOM_START && iter != randomStartTimes)
        {
            sumArea += curArea;
            sumWL += curWireLength;
            iter++;
            continue;
        }

        if (iter == randomStartTimes)
        {
            normArea = sumArea / double(iter);
            normWL = sumWL / double(iter);
            bestSolution->cost = (bestSolution->area) / normArea * _alpha + (bestSolution->wireLength) / normWL * (1 - _alpha);
        }

        double curEnergy = (_energyAlpha*curArea/normArea+(1-_energyAlpha)*curWireLength/normWL);
        if (w > _outlineWidth || h > _outlineHeight) // only consider ar when out of outline
        {
            double energyAR = double(w) / double(h) / _targetAR;
            if (energyAR < 1)
            {
                energyAR = 1 / energyAR;
            }
            if (energyAR > 1.05)
            {
                curEnergy *= energyAR * energyAR;
            }
        }
        double deltaEnergy = curEnergy - prevEnergy;

        // accept or reject

        bool reject = true;
        if (deltaEnergy <= 0 || iter == randomStartTimes)
        {
            reject = false;
        }
        else if (stage == ANNEALING || stage == RE_ANNEALING)
        {
            const double prob = exp(-deltaEnergy / T);
            const double r = genRandomDouble();
            if (r < prob)
            {
                reject = false;
            }
        }

        if (stage == GREEDY_CONVERGE && deltaEnergy > 0)
        {
            sumDeltaEnergy += deltaEnergy;
            upHillCount++;
        }

        if (reject) 
        {
            failCount++;
            if (perturbMode == 0) // rotate
            {
                blocks[idx1]->rotate();
            }
            else if (perturbMode == 1) // swap
            {
                treeNodes[idx1]->swapBlock(treeNodes[idx2]);
            }
            else if (perturbMode == 2) // move
            {
                for (auto nodeId : touchedNodeId)
                {
                    treeNodes[nodeId]->copyInfo(copyTreeNodes[nodeId]);
                }
                root = originalRoot;
            }
        }
        else
        {
            if (deltaEnergy != 0)
            {
                failCount = 0;
            }
            prevEnergy = curEnergy;
        }
        iter++;
        converge = failCount > convergeCount;
    }


    // clean up
    for (auto block : blocks)
    {
        delete block;
    }
    for (auto net : nets)
    {
        delete net;
    }
    for (auto treeNode : treeNodes)
    {
        delete treeNode;
    }
    bestSolution->iter = iter;
    return bestSolution;
}

TreeNode* FloorPlanner::removeTreeNode(TreeNode* &root, TreeNode *node, vector<int>& touchedNodeId)
{
    TreeNode* parent = node->getParent();
    TreeNode* leftChild = node->getLeftChild();
    TreeNode* rightChild = node->getRightChild();
    if (leftChild == nullptr && rightChild == nullptr)
    {
        if (parent == nullptr)
        {
            root = nullptr;
        }
        else if (parent->getLeftChild() == node)
        {
            touchedNodeId.push_back(parent->getId());
            parent->setLeftChild(nullptr);
        }
        else
        {
            touchedNodeId.push_back(parent->getId());
            parent->setRightChild(nullptr);
        }
        return node;
    }
    else if (leftChild == nullptr)
    {
        if (parent == nullptr)
        {
            root = rightChild;
        }
        else if (parent->getLeftChild() == node)
        {
            touchedNodeId.push_back(parent->getId());
            parent->setLeftChild(rightChild);
        }
        else
        {
            touchedNodeId.push_back(parent->getId());
            parent->setRightChild(rightChild);
        }
        touchedNodeId.push_back(rightChild->getId());
        rightChild->setParent(parent);
        return node;
    }
    else if (rightChild == nullptr)
    {
        if (parent == nullptr)
        {
            root = leftChild;
        }
        else if (parent->getLeftChild() == node)
        {
            touchedNodeId.push_back(parent->getId());
            parent->setLeftChild(leftChild);
        }
        else
        {
            touchedNodeId.push_back(parent->getId());
            parent->setRightChild(leftChild);
        }
        touchedNodeId.push_back(leftChild->getId());
        leftChild->setParent(parent);
        return node;
    }
    else
    {
        TreeNode* replaceNode = (genRandom(0, 1) == 0) ? leftChild : rightChild;
        touchedNodeId.push_back(replaceNode->getId());
        node->swapBlock(replaceNode);
        return removeTreeNode(root, replaceNode, touchedNodeId);
    }
}

void FloorPlanner::insertTreeNode(TreeNode* root, TreeNode *node, TreeNode *parent, vector<int>& touchedNodeId)
{
    assert(parent != nullptr);
    touchedNodeId.push_back(parent->getId());
    node->setParent(parent);
    if (genRandom(0, 1) == 0)
    {
        TreeNode* leftChild = parent->getLeftChild();
        node->setLeftChild(leftChild);
        node->setRightChild(nullptr);
        if (leftChild != nullptr)
        {
            touchedNodeId.push_back(leftChild->getId());
            leftChild->setParent(node);
        }
        parent->setLeftChild(node);
    }
    else
    {
        TreeNode* rightChild = parent->getRightChild();
        node->setRightChild(rightChild);
        node->setLeftChild(nullptr);
        if (rightChild != nullptr)
        {
            touchedNodeId.push_back(rightChild->getId());
            rightChild->setParent(node);
        }
        parent->setRightChild(node);
    }
}

void FloorPlanner::insertBlockOnContour(BlockInst *block, ContourNode *contourHead)
{
    const int blockX = block->getX();
    const int blockX2 = blockX + block->getWidth();
    ContourNode* curNode = contourHead;
    while (curNode->getX() < blockX)
    {
        curNode = curNode->getNext();
    }
    int highestY = 0;
    ContourNode* startNode = curNode; // the first node with x >= blockX
    assert(startNode->getX() == blockX);
    curNode = curNode->getNext();
    while (curNode->getX() < blockX2)
    {
        highestY = max(highestY, curNode->getY());
        curNode = curNode->getNext();
        delete curNode->getPrev();
    }
    block->setY(highestY);
    ContourNode* endNode = curNode;
    while (endNode->getNext() && endNode->getNext()->getX() == blockX2)
    {
        endNode = endNode->getNext();
        delete endNode->getPrev();
    }
    // now endNode is the last node with x == blockX2 or the first node with x > blockX2
    ContourNode* newNodeLeft = new ContourNode(blockX, highestY + block->getHeight());
    startNode->setNext(newNodeLeft);
    newNodeLeft->setPrev(startNode);
    ContourNode* newNodeRightTop = new ContourNode(blockX2, highestY + block->getHeight());
    newNodeLeft->setNext(newNodeRightTop);
    newNodeRightTop->setPrev(newNodeLeft);
    if (endNode->getX() == blockX2)
    {
        newNodeRightTop->setNext(endNode);
        endNode->setPrev(newNodeRightTop);
    }
    else
    {
        ContourNode* newNodeRight = new ContourNode(blockX2, highestY);
        newNodeRightTop->setNext(newNodeRight);
        newNodeRight->setPrev(newNodeRightTop);
        newNodeRight->setNext(endNode);
        endNode->setPrev(newNodeRight);
    }
}

void FloorPlanner::writeOutput(ofstream &outFile)
{
    if (!(_bestSolution->isInOutline(_outlineWidth, _outlineHeight)))
    {
        cout << "-1\n";
        return;
    }
    else
    {
        cout << _bestSolution->cost << "\n";
    }
    assert(this->_bestSolution != nullptr);
    const Solution *const sol = this->_bestSolution;
    outFile << sol->cost << "\n";
    outFile << sol->wireLength << "\n";
    outFile << sol->area << "\n";
    outFile << sol->width << " " << sol->height << "\n";
    const auto runTime = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - _startTime).count() / 1000.0;
    outFile << runTime << "\n";
    for (int i = 0, end_i = this->_nBlocks; i < end_i; ++i)
    {
        outFile << this->_blockLibs[i]->getName() << " " << sol->blockLowerLeft[i].first << " " << sol->blockLowerLeft[i].second << " " << sol->blockUpperRight[i].first << " " << sol->blockUpperRight[i].second << "\n";
    }
}