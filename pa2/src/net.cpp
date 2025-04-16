#include "net.h"
#include "block.h"
#include <cassert>

using namespace std;

void NetLib::addTerminal(pair<int, int> terminal)
{
    _termMinX = min(_termMinX, terminal.first);
    _termMinY = min(_termMinY, terminal.second);
    _termMaxX = max(_termMaxX, terminal.first);
    _termMaxY = max(_termMaxY, terminal.second);
}

NetInst::NetInst(NetLib *net, const vector<BlockInst *> &allBlocks) : NetLib(net), _blocks()
{
    for (int blockIdx : net->getBlockIdxs())
    {
        _blocks.push_back(allBlocks[blockIdx]);
    }
}

int NetInst::calcHPWL() const
{
    int minX = _termMinX;
    int minY = _termMinY;
    int maxX = _termMaxX;
    int maxY = _termMaxY;
    assert(!_blocks.empty());
    for (BlockInst *block : _blocks)
    {
        const int x = block->getX() + block->getWidth() / 2;
        const int y = block->getY() + block->getHeight() / 2;
        minX = min(minX, x);
        minY = min(minY, y);
        maxX = max(maxX, x);
        maxY = max(maxY, y);
    }
    return (maxX - minX) + (maxY - minY);
}
