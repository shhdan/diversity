/*
* rtreeNode.h
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/

#ifndef RTREENODE_H_
#define RTREENODE_H_

#include "Rectangle.h"
#include "IndexablePoint.h"

class RtreeNode {
private:
	/*
	*  For an internal node, listPtr is the pointer of the child node list.
	*  For a leaf node, listPtr is the pointer of the element list.
	*/
	char** listPtr;
	int listSize;

	long long computeGroupPerimeter(bool _isLeaf, char** _list, int _start, int _length, int _dim);
public:
	Rectangle mbr;
	RtreeNode* parent;

	RtreeNode();
	virtual ~RtreeNode();

	int addElement(IndexablePoint* _elem, int _dim);
	int deleteElement(IndexablePoint* elem);

	int addChildNode(RtreeNode* _child, int _dim);
	int deleteChildNode(RtreeNode* _child);
	inline int getListSize() {
		return listSize;
	}
	inline char** getListPtr() {
		return listPtr;
	}

	void initialize(bool _isLeaf, int _dim, int _branchFactor, int _leafCapacity);

	RtreeNode* chooseBestChildNode(Rectangle& _rect, int _dim);
	RtreeNode* split(bool _isLeaf, char* _toInsPtr, int _dim, int _branchFactor, int _leafCapacity);

	void releaseSpace();
};

#endif /* RTREENODE_H_ */
