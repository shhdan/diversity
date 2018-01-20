/*
* rtreeNode.cpp
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/

//#include "../Point.h"
/* added by SHD*/
//#include "stdafx.h"

#include <stdlib.h>
#include <algorithm>
#include <limits.h>

#include "RtreeNode.h"
using namespace std;

/*
*  The largest integer.
*/
//const int MYINFTY = 2147483647;
const long long MYINFTY = LLONG_MAX;

/*
*  axisIndex is in the range of [0, DIM - 1].
*/
int axisIndex = 0;

/*
*  The dimensionality of rtree nodes.
*/
int rtreeNode_dim = 2;

/*
*  Compare two given tree node pointers by the specified axisIndex.
*/
bool CompRtreeNodeLess(char* _u1, char* _u2) {
	RtreeNode* u1 = (RtreeNode*)_u1;
	RtreeNode* u2 = (RtreeNode*)_u2;
	int c1 = (u1->mbr.valueList[axisIndex + rtreeNode_dim] + u1->mbr.valueList[axisIndex]);
	int c2 = (u2->mbr.valueList[axisIndex + rtreeNode_dim] + u2->mbr.valueList[axisIndex]);
	return c1 < c2;
}

/*
*  Compare two given element pointers by the specified axisIndex.
*/
bool CompElementLess(char* _u1, char* _u2) {
	IndexablePoint* p1 = (IndexablePoint*)_u1;
	IndexablePoint* p2 = (IndexablePoint*)_u2;
	return p1->coords[axisIndex] < p2->coords[axisIndex];
}

RtreeNode::RtreeNode() {
	// TODO Auto-generated constructor stub
	listPtr = NULL;
	listSize = 0;
	parent = NULL;
}

RtreeNode::~RtreeNode() {
	// TODO Auto-generated destructor stub
	this->releaseSpace();
}

int RtreeNode::addElement(IndexablePoint* _elem, int _dim) {
	IndexablePoint* p = _elem;
	Rectangle rect;
	rect.valueList = new int[2 * _dim];
	for (int i = 0; i < _dim; i++) {
		rect.valueList[i] = p->coords[i];
		rect.valueList[i + _dim] = p->coords[i];
	}
	// Update the MBR for the parent node.
	if (this->listSize == 0) {
		this->mbr.setByRectangle(rect, _dim);
	}
	else {
		this->mbr.enlarge(rect, _dim);
	}
	this->listPtr[listSize++] = (char*)_elem;

	// Added by jhgan on Dec 19.
	rect.releaseSpace();
	//
	return listSize;
}

int RtreeNode::deleteElement(IndexablePoint* _elem) {
	if (listSize == 0) {
		// TODO Error message
		return listSize;
	}
	for (int i = 0; i < listSize; i++) {
		if (listPtr[i] == (char*)_elem) {
			// Move the last child node to index i and remove the last one.
			listPtr[i] = listPtr[listSize - 1];
			listPtr[listSize - 1] = NULL;
			listSize--;

			// delete this element
			// TODO here need to release the space of the element.
			//			delete _elem;
			//			_elem = NULL;
			break;
		}
	}
	return listSize;
}

int RtreeNode::addChildNode(RtreeNode* _child, int _dim) {
	// Update the MBR for the parent node.
	if (this->listSize == 0) {
		this->mbr.setByRectangle(_child->mbr, _dim);
	}
	else {
		this->mbr.enlarge(_child->mbr, _dim);
	}
	this->listPtr[listSize++] = (char*)_child;
	_child->parent = this;
	return listSize;
}

int RtreeNode::deleteChildNode(RtreeNode* _child) {
	if (listSize == 0) {
		// TODO Error message
		return listSize;
	}
	for (int i = 0; i < listSize; i++) {
		if (listPtr[i] == (char*)_child) {
			// Move the last child node to index i and remove the last one.
			listPtr[i] = listPtr[listSize - 1];
			listPtr[listSize - 1] = NULL;
			listSize--;

			// delete this child node
			_child->releaseSpace();
			delete _child;
			_child = NULL;
			break;
		}
	}
	return listSize;
}

//int RtreeNode::getListSize() {
//	return listSize;
//}
//
//char** RtreeNode::getListPtr() {
//	return listPtr;
//}

void RtreeNode::initialize(bool _isLeaf, int _dim, int _branchFactor, int _leafCapacity) {
	if (listPtr != NULL) {
		// TODO Error detected.
	}
	if (_isLeaf) {
		listPtr = new char*[_leafCapacity];
	}
	else {
		listPtr = new char*[_branchFactor];
	}
	listSize = 0;
}

RtreeNode* RtreeNode::chooseBestChildNode(Rectangle& _rect, int _dim) {
	int bestIndex = -1;
	long long bestValue = MYINFTY;
	long long p = 0;

	int childNodeNum = listSize;
	RtreeNode** childList = (RtreeNode**)listPtr;
	for (int i = 0; i < childNodeNum; ++i) {
		p = childList[i]->mbr.perimeterIncrement(_rect, _dim);
		if (bestValue > p) {
			bestValue = p;
			bestIndex = i;
		}
		else if (bestValue == p && bestIndex >= 0) {
			long long p1 = childList[i]->mbr.computePerimeter(_dim);
			long long p2 = childList[bestIndex]->mbr.computePerimeter(_dim);
			if (p1 < p2) {
				// If the increment of perimeter are the same, then pick the one whose perimeter is smaller.
				bestValue = p;
				bestIndex = i;
			}
		}
	}
	return childList[bestIndex];
}

RtreeNode* RtreeNode::split(bool _isLeaf, char* _toInsPtr, int _dim, int _branchFactor, int _leafCapacity) {
	int childNodeNum = listSize + 1;
	int half = childNodeNum / 2;
	int rest = childNodeNum - half;
	long long p1 = 0;
	long long p2 = 0;
	long long bestValue = MYINFTY;
	int bestIndex = -1;

	char** tempList = new char*[childNodeNum];
	for (int i = 0; i < listSize; i++) {
		tempList[i] = listPtr[i];
	}
	tempList[listSize] = _toInsPtr;

	// Set the dimensionality of the rtree nodes.
	rtreeNode_dim = _dim;
	// Set the compare function pointer.
	bool(*comp)(char*, char*) = CompRtreeNodeLess;
	if (_isLeaf) {
		comp = CompElementLess;
	}

	// Find the best split axis.
	for (int i = 0; i < _dim; i++) {
		// sort the child nodes on each dimension coordinate.
		axisIndex = i;
		sort(tempList, tempList + childNodeNum, comp);
		p1 = this->computeGroupPerimeter(_isLeaf, tempList, 0, half, _dim);
		p2 = this->computeGroupPerimeter(_isLeaf, tempList, half, rest, _dim);
		if (p1 + p2 < bestValue) {
			// The smaller value the better.
			bestValue = p1 + p2;
			bestIndex = axisIndex;
		}
	}

	// Split this node by bestIndex.
	axisIndex = bestIndex;
	sort(tempList, tempList + childNodeNum, comp);

	// Create a new node u.
	RtreeNode* u = new RtreeNode();
	u->initialize(_isLeaf, _dim, _branchFactor, _leafCapacity);
	u->parent = this->parent;

	if (!_isLeaf) {
		// Reset this node.
		listSize = 0;
		for (int i = 0; i < half; ++i) {
			this->addChildNode((RtreeNode*)(tempList[i]), _dim);
		}
		for (int i = half; i < childNodeNum; ++i) {
			u->addChildNode((RtreeNode*)(tempList[i]), _dim);
		}
	}
	else {
		// Reset this node.
		listSize = 0;
		for (int i = 0; i < half; ++i) {
			this->addElement((IndexablePoint*)(tempList[i]), _dim);
		}
		for (int i = half; i < childNodeNum; ++i) {
			u->addElement((IndexablePoint*)(tempList[i]), _dim);
		}
	}
	return u;
}

void RtreeNode::releaseSpace() {
	mbr.releaseSpace();
	if (listPtr != NULL) {
		delete[] listPtr;
		listPtr = NULL;
	}
	listSize = 0;
}

long long RtreeNode::computeGroupPerimeter(bool _isLeaf, char** _list, int _start, int _length, int _dim) {
	if (_start < 0 || _length < 1) {
		// TODO Error message.
		return -1;
	}

	long long rtn = 0;
	int* valueList = new int[_dim * 2];
	int* minList = valueList;
	int* maxList = &(valueList[_dim]);

	if (!_isLeaf) {
		RtreeNode** childList = (RtreeNode**)_list;
		for (int k = 0; k < _dim * 2; k++) {
			valueList[k] = childList[_start]->mbr.valueList[k];
		}

		int end = _start + _length;
		for (int j = _start + 1; j < end; ++j) {
			for (int k = 0; k < _dim; ++k) {
				if (childList[j]->mbr.valueList[k] < minList[k])
					minList[k] = childList[j]->mbr.valueList[k];
				if (childList[j]->mbr.valueList[k + _dim] > maxList[k])
					maxList[k] = childList[j]->mbr.valueList[k + _dim];
			}
		}
	}
	else {
		// This is a leaf node.
		IndexablePoint** pointList = (IndexablePoint**)_list;
		for (int k = 0; k < _dim; k++) {
			minList[k] = pointList[_start]->coords[k];
			maxList[k] = pointList[_start]->coords[k];
		}

		int end = _start + _length;
		for (int j = _start + 1; j < end; ++j) {
			for (int k = 0; k < _dim; ++k) {
				if (pointList[j]->coords[k] < minList[k])
					minList[k] = pointList[j]->coords[k];
				if (pointList[j]->coords[k] > maxList[k])
					maxList[k] = pointList[j]->coords[k];
			}
		}
	}

	// Compute the perimeter of the mbr defined by minList and maxList.
	for (int i = 0; i < _dim; i++) {
		rtn += (maxList[i] - minList[i]);
	}
	rtn *= 2;

	delete[] valueList;

	return rtn;
}



