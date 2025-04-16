#pragma once
#include <vector>
#include "config.h"

class BlockInst;

class NetLib
{
public:
    NetLib() {}
    NetLib(NetLib* net) : _blockIdxs(), _termMinX(net->_termMinX), _termMinY(net->_termMinY), _termMaxX(net->_termMaxX), _termMaxY(net->_termMaxY) {}
    ~NetLib() {}

    inline const std::vector<int>& getBlockIdxs() const { return _blockIdxs; }

    inline void addBlockIdx(int blockIdx) { _blockIdxs.push_back(blockIdx); }
    void addTerminal(std::pair<int, int> terminal);

protected:
    std::vector<int> _blockIdxs;
    int _termMinX = INF;
    int _termMinY = INF;
    int _termMaxX = 0;
    int _termMaxY = 0;
};

class NetInst : public NetLib
{
public:
    NetInst(NetLib *net, const std::vector<BlockInst *> &allBlocks);

    int calcHPWL() const;
private:
    std::vector<BlockInst *> _blocks;
};