#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
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

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
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

void Partitioner::partition()
{
    // set up balance factor
    const int minPartSize = ceil((1 - _bFactor) / 2 * _cellNum);
    const int maxPartSize = floor((1 + _bFactor) / 2 * _cellNum);
    // set up Pmax
    for (auto cell : _cellArray) {
        _maxPinNum = max(_maxPinNum, cell->getPinNum());
    }
    for (int i = -_maxPinNum; i <= _maxPinNum; ++i) {
        _bList[0][i] = new Node(-1);
        _bList[1][i] = new Node(-1);
    }
    // init partition
    for (size_t i = _cellNum / 2, end = _cellNum; i < end; ++i) {
        _cellArray[i]->setPart(1);
        _partSize[1]++;
    }
    _partSize[0] = _cellArray.size() - _partSize[1];
    int initCutSize = 0;
    for (auto net : _netArray) {
        for (int cellId : net->getCellList()) {
            net->incPartCount(_cellArray[cellId]->getPart());
        }
        if (net->getPartCount(0) > 0 && net->getPartCount(1) > 0) {
            ++initCutSize;
        }
    }
    cout << "Initial cutsize = " << initCutSize << ", partSize[0] = " << _partSize[0] << ", partSize[1] = " << _partSize[1] << "\n";
    _cutSize = initCutSize;
    // iterate FM
    _iterNum = 0;
    do
    {
        // reset bucket lists
        for (int i = -_maxPinNum; i <= _maxPinNum; ++i) {
            _bList[0][i]->setNext(nullptr);
            _bList[1][i]->setNext(nullptr);
        }
        // reset cell gains
        _maxGain[0] = _maxGain[1] = -_maxPinNum-1;
        for (auto cell : _cellArray) {
            const int part = cell->getPart();
            const vector<int> netList = cell->getNetList();
            Net* net = nullptr;
            int gain = 0;
            for (int netId : netList) {
                net = _netArray[netId];
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
        _moveStack.resize(_cellNum, 0);
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
            Node* const maxNode = _bList[chosenPart][_maxGain[chosenPart]]->getNext();
            if (maxNode == nullptr) {
                break;
            }
            const int cellId = maxNode->getId();
            _moveStack[_moveNum] = cellId;
            Cell* cell = _cellArray[cellId];
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
            const vector<int> netList = cell->getNetList();
            for (int netId : netList) {
                Net* net = _netArray[netId];
                // check critical nets before moving
                const int partCountT = net->getPartCount(T);
                if (partCountT == 0) {
                    const vector<int> cellList = net->getCellList();
                    for (int cellId : cellList) {
                        Cell* updateCell = _cellArray[cellId];
                        if (!updateCell->getLock()) {
                            removeBList(F, updateCell->getGain(), updateCell->getNode());
                            updateCell->incGain();
                            insertBList(F, updateCell->getGain(), updateCell->getNode());
                        }
                    }
                }
                else if (partCountT == 1) {
                    const vector<int> cellList = net->getCellList();
                    for (int cellId : cellList) {
                        Cell* updateCell = _cellArray[cellId];
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
                    const vector<int> cellList = net->getCellList();
                    for (int cellId : cellList) {
                        Cell* updateCell = _cellArray[cellId];
                        if (!updateCell->getLock()) {
                            removeBList(T, updateCell->getGain(), updateCell->getNode());
                            updateCell->decGain();
                            insertBList(T, updateCell->getGain(), updateCell->getNode());
                        }
                    }
                }
                else if (partCountF == 1) {
                    const vector<int> cellList = net->getCellList();
                    for (int cellId : cellList) {
                        Cell* updateCell = _cellArray[cellId];
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
            Cell* const cell = _cellArray[_moveStack[i]];
            const int F = cell->getPart();
            const int T = 1-F;
            cell->move();
            _partSize[F]--;
            _partSize[T]++;
            const vector<int> netList = cell->getNetList();
            for (int netId : netList) {
                Net* net = _netArray[netId];
                net->decPartCount(F);
                net->incPartCount(T);
            }
        }
        _cutSize -= _maxAccGain;
        cout << "Iteration " << _iterNum << ": maxAccGain = " << _maxAccGain << ", cutSize = " << _cutSize << ", moveNum = " << _bestMoveNum << ", partSize[0] = " << _partSize[0] << ", partSize[1] = " << _partSize[1] << "\n";
        _iterNum++;
    } while (_maxAccGain > 0);
}

void Partitioner::insertBList(int part, int gain, Node* node)
{
    Node* head = _bList[part][gain];
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
    while (_maxGain[part] > -_maxPinNum && _bList[part][_maxGain[part]]->getNext() == nullptr) {
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

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
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
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
