#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <fstream>
#include <vector>
#include <unordered_map>
#include "cell.h"
#include "net.h"
using namespace std;

class Partitioner
{
public:
    // constructor and destructor
    Partitioner(fstream& inFile) :
        _cutSize(0), _netNum(0), _cellNum(0), _maxPinNum(0), _bFactor(0),
        _accGain(0), _maxAccGain(0), _iterNum(0) {
        parseInput(inFile);
        _partSize[0] = 0;
        _partSize[1] = 0;
    }
    ~Partitioner() {
        clear();
    }

    // copy constructor
    Partitioner(const Partitioner& part) :
        _cutSize(0), _netNum(part._netNum), _cellNum(part._cellNum),
        _maxPinNum(part._maxPinNum), _bFactor(part._bFactor),
        _accGain(part._accGain), _maxAccGain(part._maxAccGain),
        _iterNum(part._iterNum), _moveNum(0), _bestMoveNum(0) {
            // copy net array
            for (auto net : part._netArray) {
                Net* newNet = new Net(*net);
                _netArray.push_back(newNet);
            }
            // copy cell array
            for (auto cell : part._cellArray) {
                Cell* newCell = new Cell(*cell);
                _cellArray.push_back(newCell);
                // set up the net list of the new cell
                for (auto net : cell->getNetList()) {
                    newCell->addNet(_netArray[net->getId()]);
                }
            }
            // set up the cell list of the new net
            for (int i = 0; i < _netNum; ++i) {
                Net* net = _netArray[i];
                for (auto cell : part._netArray[i]->getCellList()) {
                    net->addCell(_cellArray[cell->getId()]);
                }
            }
            _partSize[0] = _partSize[1] = 0;
    }

    // basic access methods
    int getCutSize() const          { return _cutSize; }
    int getNetNum() const           { return _netNum; }
    int getCellNum() const          { return _cellNum; }
    double getBFactor() const       { return _bFactor; }
    int getPartSize(int part) const { return _partSize[part]; }

    // modify method
    void parseInput(fstream& inFile);
    void initPartition();
    void initPartition(const int seed);
    void partition();

    // member functions about reporting
    void printSummary() const;
    void writeResult(fstream& outFile);
    
    void insertBList(int part, int gain, Node* node);
    void removeBList(int part, int gain, Node* node);

private:
    int                 _cutSize;       // cut size
    int                 _partSize[2];   // size (cell number) of partition A(0) and B(1)
    int                 _netNum;        // number of nets
    int                 _cellNum;       // number of cells
    int                 _maxPinNum;     // Pmax for building bucket list
    double              _bFactor;       // the balance factor to be met
    // int                 _maxGainPart;   // partition of max gain cell
    int                 _maxGain[2];    // max gain
    // Node*               _maxGainNode[2];   // pointer to max gain cell
    vector<Net*>        _netArray;      // net array of the circuit
    vector<Cell*>       _cellArray;     // cell array of the circuit
    vector<Node*>       _bList[2];      // bucket list of partition A(0) and B(1)

    int                 _accGain;       // accumulative gain
    int                 _maxAccGain;    // maximum accumulative gain
    int                 _moveNum;       // number of cell movements
    int                 _iterNum;       // number of iterations
    int                 _bestMoveNum;   // store best number of movements
    int                 _unlockNum[2];  // number of unlocked cells
    vector<Cell*>       _moveStack;     // history of cell movement

    // Clean up partitioner
    void clear();
    inline int gain2blistId(int gain) const { return gain + _maxPinNum; }
};

#endif  // PARTITIONER_H
