#ifndef CELL_H
#define CELL_H

#include <vector>
using namespace std;

class Net;
class Cell;

class Node
{
    friend class Cell;

public:
    // Constructor and destructor
    Node(Cell* cell) : _cell(cell), _prev(nullptr), _next(nullptr) { }
    ~Node() { }

    // Basic access methods
    Cell* getCell() const   { return _cell; }
    Node* getPrev() const   { return _prev; }
    Node* getNext() const   { return _next; }

    // Set functions
    void setPrev(Node* prev)  { _prev = prev; }
    void setNext(Node* next)  { _next = next; }

private:
    Cell*       _cell;  // pointer to the cell
    Node*       _prev;  // pointer to the previous node
    Node*       _next;  // pointer to the next node
};

class Cell
{
public:
    // Constructor and destructor
    Cell(string& name, bool part, int id) :
        _id(id), _gain(0), _pinNum(0), _part(part), _lock(false), _name(name) {
        _node = new Node(this);
    }
    ~Cell() { }
    Cell(const Cell& cell) :
        _id(cell._id), _gain(cell._gain), _pinNum(cell._pinNum),
        _part(cell._part), _lock(cell._lock), _name(cell._name) {
        _node = new Node(this);
    }

    // Basic access methods
    int getId() const      { return _id; }
    int getGain() const     { return _gain; }
    int getPinNum() const   { return _pinNum; }
    bool getPart() const    { return _part; }
    bool getLock() const    { return _lock; }
    Node* getNode() const   { return _node; }
    string getName() const  { return _name; }
    // int getFirstNet() const { return _netList[0]; }
    const vector<Net*>& getNetList() const  { return _netList; }

    // Set functions
    void setNode(Node* node)        { _node = node; }
    void setGain(const int gain)    { _gain = gain; }
    void setPart(const bool part)   { _part = part; }
    void setName(const string name) { _name = name; }

    // Modify methods
    void move()         { _part = !_part; }
    void lock()         { _lock = true; }
    void unlock()       { _lock = false; }
    void incGain()      { ++_gain; }
    void decGain()      { --_gain; }
    void incPinNum()    { ++_pinNum; }
    void decPinNum()    { --_pinNum; }
    void addNet(Net* net) { _netList.push_back(net); }

private:
    int             _id;        // id of the cell
    int             _gain;      // gain of the cell
    int             _pinNum;    // number of pins the cell are connected to
    bool            _part;      // partition the cell belongs to (0-A, 1-B)
    bool            _lock;      // whether the cell is locked
    Node*           _node;      // node used to link the cells together
    string          _name;      // name of the cell
    vector<Net*>     _netList;   // list of nets the cell is connected to
};

#endif  // CELL_H
