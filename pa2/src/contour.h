#pragma once

class ContourNode
{
public:
    ContourNode() : _x(0), _y(0), _prev(nullptr), _next(nullptr) {}
    ContourNode(int x, int y) : _x(x), _y(y), _prev(nullptr), _next(nullptr) {}
    ~ContourNode() {}

    inline int getX() const { return _x; }
    inline int getY() const { return _y; }
    inline ContourNode* getPrev() const { return _prev; }
    inline ContourNode* getNext() const { return _next; }

    inline void setPrev(ContourNode* prev) { _prev = prev; }
    inline void setNext(ContourNode* next) { _next = next; }

private:
    int _x;
    int _y;
    ContourNode* _prev;
    ContourNode* _next;
};