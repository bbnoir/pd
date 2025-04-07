#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    unordered_map<string, int> cellName2Id;

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            Net* net = new Net(netName, _netNum);
            _netArray.push_back(net);
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        Cell* cell = new Cell(cellName, 0, cellId);
                        _cellArray.push_back(cell);
                        cellName2Id[cellName] = cellId;
                        cell->addNet(net);
                        cell->incPinNum();
                        net->addCell(cell);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(cellName2Id.count(cellName) == 1);
                            int cellId = cellName2Id[cellName];
                            Cell* cell = _cellArray[cellId];
                            cell->addNet(net);
                            cell->incPinNum();
                            net->addCell(cell);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::initPartition()
{
    // set first half of cells to part 0 and the second half to part 1
    for (size_t i = _cellNum / 2, end = _cellNum; i < end; ++i) {
        _cellArray[i]->setPart(1);
        _partSize[1]++;
    }
    _partSize[0] = _cellArray.size() - _partSize[1];
}

void Partitioner::initPartition(const int seed)
{
    mt19937 gen(seed);
    uniform_int_distribution<> dis(0, 1);
    for (size_t i = 0, end = _cellNum; i < end; ++i) {
        const int part = dis(gen);
        _cellArray[i]->setPart(part);
        ++_partSize[part];
    }
    _seed = seed;
}

void Partitioner::partition()
{
    // set up balance factor
    const int minPartSize = ceil((1 - _bFactor) / 2 * _cellNum);
    const int maxPartSize = floor((1 + _bFactor) / 2 * _cellNum);
    // set up Pmax
    for (auto cell : _cellArray) {
        _maxPinNum = max(_maxPinNum, cell->getPinNum());
    }
    const int bListSize = 2 * _maxPinNum + 1;
    _bList[0].resize(bListSize, nullptr);
    _bList[1].resize(bListSize, nullptr);
    for (int i = 0, end_i = bListSize; i < end_i; ++i) {
        _bList[0][i] = new Node(nullptr);
        _bList[1][i] = new Node(nullptr);
    }
    int initCutSize = 0;
    for (auto net : _netArray) {
        for (auto cell : net->getCellList()) {
            net->incPartCount(cell->getPart());
        }
        if (net->getPartCount(0) > 0 && net->getPartCount(1) > 0) {
            ++initCutSize;
        }
    }
    _cutSize = initCutSize;
    // iterate FM
    _iterNum = 0;
    mt19937 gen(_seed);
    do
    {
        // reset bucket lists
        for (int i = 0, end_i = bListSize; i < end_i; ++i) {
            _bList[0][i]->setNext(nullptr);
            _bList[1][i]->setNext(nullptr);
        }
        // reset cell gains
        _maxGain[0] = _maxGain[1] = -_maxPinNum-1;
        vector<int> cellIndices(_cellArray.size());
        iota(cellIndices.begin(), cellIndices.end(), 0);
        if (_seed != -1) {
            shuffle(cellIndices.begin(), cellIndices.end(), gen);
        }
        for (int idx : cellIndices) {
            auto cell = _cellArray[idx];
            const int part = cell->getPart();
            int gain = 0;
            for (auto net : cell->getNetList()) {
                if (net->getPartCount(part) == 1) {
                    ++gain;
                }
                if (net->getPartCount(1-part) == 0) {
                    --gain;
                }
            }
            cell->setGain(gain);
            cell->unlock();
            insertBList(part, gain, cell->getNode());
        }
        // move cells
        _accGain = 0;
        _maxAccGain = 0;
        _bestMoveNum = 0;
        _moveStack.clear();
        _moveStack.resize(_cellNum, nullptr);
        for (_moveNum = 0; _moveNum < _cellNum; ++_moveNum) {
            int chosenPart = 0;
            if (_partSize[0] <= minPartSize) {
                chosenPart = 1;
            }
            else if (_partSize[1] <= minPartSize) {
                chosenPart = 0;
            }
            else if (_maxGain[0] >= _maxGain[1]) {
                chosenPart = 0;
            }
            else {
                chosenPart = 1;
            }
            Node* const maxNode = _bList[chosenPart][gain2blistId(_maxGain[chosenPart])]->getNext();
            if (maxNode == nullptr) {
                break;
            }
            Cell* const cell = maxNode->getCell();
            _moveStack[_moveNum] = cell;
            const int curGain = cell->getGain();
            _accGain += curGain;
            if (_accGain > _maxAccGain) {
                _maxAccGain = _accGain;
                _bestMoveNum = _moveNum;
            }
            // update gain
            const int F = cell->getPart();
            const int T = 1-F;
            cell->move();
            cell->lock();
            --_partSize[F];
            ++_partSize[T];
            removeBList(F, curGain, maxNode);
            for (auto net : cell->getNetList()) {
                // check critical nets before moving
                const int partCountT = net->getPartCount(T);
                if (partCountT == 0) {
                    for (auto updateCell : net->getCellList()) {
                        if (!updateCell->getLock()) {
                            removeBList(F, updateCell->getGain(), updateCell->getNode());
                            updateCell->incGain();
                            insertBList(F, updateCell->getGain(), updateCell->getNode());
                        }
                    }
                }
                else if (partCountT == 1) {
                    for (auto updateCell : net->getCellList()) {
                        if (updateCell->getPart() == T && updateCell != cell) {
                            if (!updateCell->getLock()) {
                                removeBList(T, updateCell->getGain(), updateCell->getNode());
                                updateCell->decGain();
                                insertBList(T, updateCell->getGain(), updateCell->getNode());
                            }
                            break;
                        }
                    }
                }
                net->decPartCount(F);
                net->incPartCount(T);
                const int partCountF = net->getPartCount(F);
                if (partCountF == 0) {
                    for (auto updateCell : net->getCellList()) {
                        if (!updateCell->getLock()) {
                            removeBList(T, updateCell->getGain(), updateCell->getNode());
                            updateCell->decGain();
                            insertBList(T, updateCell->getGain(), updateCell->getNode());
                        }
                    }
                }
                else if (partCountF == 1) {
                    for (auto updateCell : net->getCellList()) {
                        if (updateCell->getPart() == F) {
                            if (!updateCell->getLock()) {
                                removeBList(F, updateCell->getGain(), updateCell->getNode());
                                updateCell->incGain();
                                insertBList(F, updateCell->getGain(), updateCell->getNode());
                            }
                            break;
                        }
                    }
                }
            }
        }
        // backtrack
        for (int i = _moveNum - 1; i > _bestMoveNum; --i) {
            Cell* const cell = _moveStack[i];
            const int F = cell->getPart();
            const int T = 1-F;
            cell->move();
            _partSize[F]--;
            _partSize[T]++;
            for (auto net : cell->getNetList()) {
                net->decPartCount(F);
                net->incPartCount(T);
            }
        }
        _cutSize -= _maxAccGain;
        _iterNum++;
    } while (_maxAccGain > 0);
}

void Partitioner::insertBList(int part, int gain, Node* node)
{
    Node* head = _bList[part][gain2blistId(gain)];
    Node* next = head->getNext();
    node->setNext(next);
    node->setPrev(head);
    head->setNext(node);
    if (next != nullptr) {
        next->setPrev(node);
    }
    _maxGain[part] = max(_maxGain[part], gain);
}

void Partitioner::removeBList(int part, int gain, Node* node)
{
    Node* prev = node->getPrev();
    Node* next = node->getNext();
    prev->setNext(next);
    if (next != nullptr) {
        next->setPrev(prev);
    }
    while (_maxGain[part] > -_maxPinNum && _bList[part][gain2blistId(_maxGain[part])]->getNext() == nullptr) {
        --_maxGain[part];
    }
    node->setPrev(nullptr);
    node->setNext(nullptr);
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (auto cell : _cellArray) {
        if (cell->getPart() == 0) {
            outFile << cell->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (auto cell : _cellArray) {
        if (cell->getPart() == 1) {
            outFile << cell->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (auto& cell : _cellArray) {
        delete cell->getNode();
        delete cell;
    }
    for (auto& net : _netArray) {
        delete net;
    }
    for (auto& node : _bList[0]) {
        delete node;
    }
    for (auto& node : _bList[1]) {
        delete node;
    }
    return;
}
