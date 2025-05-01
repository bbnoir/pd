#pragma once
#include <iostream>
#include <utility>

class BlockInst;

class TreeNode {
public:
    TreeNode(int id) : _id(id) {}
    TreeNode(int id, BlockInst *block) : _id(id), _block(block) {}
    TreeNode(const TreeNode &other) : _id(other._id), _block(other._block), _leftChild(other._leftChild), _rightChild(other._rightChild), _parent(other._parent) {}
    
    ~TreeNode() {} // Note: The _block is not deleted in the destructor

    inline int getId() const { return _id; }
    inline BlockInst *getBlock() const { return _block; }
    inline TreeNode *getLeftChild() const { return _leftChild; }
    inline TreeNode *getRightChild() const { return _rightChild; }
    inline TreeNode *getParent() const { return _parent; }

    inline void setBlock(BlockInst *block) { _block = block; }
    inline void setLeftChild(TreeNode *leftChild) { _leftChild = leftChild; }
    inline void setRightChild(TreeNode *rightChild) { _rightChild = rightChild; }
    inline void setParent(TreeNode *parent) { _parent = parent; }

    inline void swapBlock(TreeNode *other) { std::swap(_block, other->_block); }
    inline void copyInfo(const TreeNode *other) { _id = other->_id; _block = other->_block; _leftChild = other->_leftChild; _rightChild = other->_rightChild; _parent = other->_parent; }

private:
    int _id = -1;
    BlockInst *_block = nullptr;
    TreeNode *_leftChild = nullptr;
    TreeNode *_rightChild = nullptr;
    TreeNode *_parent = nullptr;
};