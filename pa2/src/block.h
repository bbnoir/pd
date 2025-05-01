#pragma once
#include <string>
#include <utility>

class BlockLib {
public:
    BlockLib(std::string name, int width, int height) : _name(name), _width(width), _height(height) {}
    BlockLib(BlockLib *block) : _name(block->getName()), _width(block->getWidth()), _height(block->getHeight()) {}
    ~BlockLib() {}

    inline std::string getName() const { return _name; }
    inline int getWidth() const { return _width; }
    inline int getHeight() const { return _height; }
    inline int getArea() const { return _width * _height; }

protected:
    std::string _name;
    int _width;
    int _height;
};

class BlockInst : public BlockLib {
public:
    BlockInst(BlockLib *block) : BlockLib(block), _x(0), _y(0) {}

    inline int getX() const { return _x; }
    inline int getY() const { return _y; }

    inline void setX(int x) { _x = x; }
    inline void setY(int y) { _y = y; }
    inline void rotate() { std::swap(_width, _height); }

private:
    int _x;
    int _y;
};