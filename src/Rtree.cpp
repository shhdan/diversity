/*
* rtree.cpp
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/
/* added by SHD*/
#define _CRT_SECURE_NO_DEPRECATE
//#include "stdafx.h"

#include "Rtree.h"

#include <stdlib.h>
#include <stdio.h>
#include <queue>
#include <algorithm>
#include <math.h>

#include "zorder.h"
using namespace std;

// === Codes for debug.
long long VisitedNodeNum = 0;
long long* VisitedNodeNumList = new long long[20];
// ===

// This value indicates the number of child nodes of an internal node
// is in the range of [fanout/FanoutRatio, fanout].
const int FanoutRatio = 4;

/*
*  The tool for sorting points by zorder.
*/
struct PointZorderSorter {
	int dim;
	PointZorderSorter(int _dim) :
		dim(_dim) {
	}
	bool operator()(const IndexablePoint* p1, const IndexablePoint* p2) {
		return compZorderByCoord(p1->coords, p2->coords, dim);
	}
};

/*
*  The tool for sorting points by the specified dimension coordinates.
*/
struct PointDimSorter {
	int dim;
	int dimIndex; // in [0, dim-1]
	PointDimSorter(int _dim, int _dimLev) :
		dim(_dim), dimIndex(_dimLev - 1) {
	}
	bool operator()(const IndexablePoint* p1, const IndexablePoint* p2) {
		return p1->coords[dimIndex] < p2->coords[dimIndex];
	}
};

struct NodeDimSorter {
	int dim;
	int dimIndex; // in [0, dim -1]
	NodeDimSorter(int _dim, int _dimLev) :
		dim(_dim), dimIndex(_dimLev - 1) {

	}
	bool operator()(const RtreeNode* node1, const RtreeNode* node2) {
		int coord1 = node1->mbr.valueList[dimIndex + dim] + node1->mbr.valueList[dimIndex];
		int coord2 = node2->mbr.valueList[dimIndex + dim] + node2->mbr.valueList[dimIndex];
		return coord1 < coord2;
		//		int coord1 = node1->mbr.valueList[dimIndex];
		//		int coord2 = node2->mbr.valueList[dimIndex];
		//		return coord1 < coord2;
	}
};

struct heapElement {
	int lev;
	char* ptr;
	double dist;
};

struct greaterThanByDist {
	bool operator()(const heapElement& _e1, const heapElement& _e2) const {
		return _e1.dist > _e2.dist;
	}
};

/*
*  The structure used in Rtree query algorithm.
*/
struct QueElement {
	RtreeNode* node;
	int lev;
};

/***************************************************************************/

Rtree::Rtree() {
	// TODO Auto-generated constructor stub
	root = NULL;
	rootLev = 0;
	dim = -1;
	branchFactor = 12;
	leafCapacity = 12;
}

Rtree::Rtree(int _dim, int _branchFactor, int _leafCapacity) {
	root = NULL;
	rootLev = 0;
	dim = _dim;
	branchFactor = _branchFactor;
	leafCapacity = _leafCapacity;
}

Rtree::~Rtree() {
	// TODO Auto-generated destructor stub
	this->releaseSpace();
}
void Rtree::initialize(int _dim, int _branchFactor, int _leafCapacity) {
	root = NULL;
	rootLev = 0;
	dim = _dim;
	branchFactor = _branchFactor;
	leafCapacity = _leafCapacity;
}

void Rtree::constructEmptyTree(int _dim, int _branchFactor, int _leafCapacity) {
	if (root != NULL) {
		// TODO error message
		// The tree is not empty!
		printf("Error in Rtree::constructRtree: The tree is not empty!\n");
	}
	dim = _dim;
	branchFactor = _branchFactor;
	leafCapacity = _leafCapacity;
	rootLev = 0;
	root = new RtreeNode();
	root->initialize(true, dim, branchFactor, leafCapacity);
}

RtreeNode* Rtree::constructRtree(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
	int _leafCapacity) {
	if (root != NULL) {
		printf("Error in Rtree::constructRtree: The tree is not empty!\n");
		return NULL;
	}

	dim = _dim;
	branchFactor = _branchFactor;
	leafCapacity = _leafCapacity;
	rootLev = 0;

	int num = _elemNum;
	if (num == 0) {
		// TODO error message.
		printf("Error in Rtree::ConstructRtree: the query list is empty!\n");
		return NULL;
	}
	else if (num <= leafCapacity) {
		root = new RtreeNode();
		root->initialize(true, dim, branchFactor, leafCapacity);
		for (int i = 0; i < num; i++)
			root->addElement(_elementList[i], dim);
		return root;
	}

	vector<RtreeNode*> treeNodeList;
	vector<RtreeNode*> parentList;
	//	Point** ptList = (Point**) _elementList;

	IndexablePoint** ptList = new IndexablePoint*[num];
	for (int i = 0; i < num; i++) {
		ptList[i] = _elementList[i];
	}

	/*************** Sort the Zorder directly by the coordinates of queries ****************/
	//	rtree_dim = dim;
	//	sort(ptList, ptList + num, compPtByZorder);
	sort(ptList, ptList + num, PointZorderSorter(dim));

	int leafNum = ceil(num * 1.0 / leafCapacity);
	treeNodeList.resize(leafNum);
	int cnt = 0;
	int tempNum = leafCapacity;
	// Construct the leaf level.
	for (int i = 0; i < leafNum; i++) {
		treeNodeList[i] = new RtreeNode();
		treeNodeList[i]->initialize(true, dim, branchFactor, leafCapacity);
		if (num < cnt + leafCapacity)
			tempNum = num - cnt;
		for (int j = 0; j < tempNum; j++) {
			treeNodeList[i]->addElement(ptList[cnt++], dim);
		}
	}

	int level = 0;
	int treeNodeNum = 0;
	int parentNum = 0;
	int rest = 0;
	int index = 0;

	while (treeNodeList.size() > 1) {

		treeNodeNum = treeNodeList.size();
		parentNum = ceil(treeNodeNum * 1.0 / branchFactor);
		parentList.resize(parentNum);

		// The size of the last group of tree nodes.
		rest = treeNodeNum - branchFactor * (treeNodeNum / branchFactor);

		for (int i = 0; i < treeNodeNum - rest; i += branchFactor) {
			index = i / branchFactor;
			parentList[index] = new RtreeNode();
			parentList[index]->initialize(false, dim, branchFactor, leafCapacity);

			for (int j = 0; j < branchFactor; ++j) {
				parentList[index]->addChildNode(treeNodeList[i + j], dim);
			}
		}

		if (rest > 0) {
			// Add the rest tree nodes to parent.
			parentList[parentNum - 1] = new RtreeNode();
			parentList[parentNum - 1]->initialize(false, dim, branchFactor, leafCapacity);
			for (int i = treeNodeNum - rest; i < treeNodeNum; ++i) {
				parentList[parentNum - 1]->addChildNode(treeNodeList[i], dim);
			}
		}
		treeNodeList.swap(parentList);
		parentList.clear();
		level++;
	}
	rootLev = level;
	root = treeNodeList[0];

	delete[] ptList;
	return root;
}

int Rtree::insertNode(RtreeNode* _toInsNode, int _targetLevel) {
	if (root == NULL) {
		if (dim == -1) {
			printf("Error in Rtree::insertNode: dim = -1! Please initialize the rtree first!\n");
			return 1;
		}
		root = _toInsNode;
		return 0;
	}
	if (rootLev == _targetLevel) {
		// In this case, we need to create a new root.
		RtreeNode* u = new RtreeNode();
		u->initialize(false, dim, branchFactor, leafCapacity);
		u->addChildNode(this->root, dim);
		u->addChildNode(_toInsNode, dim);
		this->root = u;
		rootLev = _targetLevel + 1;
	}
	else {
		this->insertNode_SubRoutine(this->root, rootLev, _toInsNode, _targetLevel);
	}
	return 0;
}

int Rtree::insertElement(IndexablePoint* _elem) {
	if (root == NULL) {
		if (dim == -1) {
			printf("Error in Rtree::insertElement: dim = -1! Please initialize the rtree first!\n");
			return 1;
		}
		root = new RtreeNode();
		root->initialize(true, dim, branchFactor, leafCapacity);
		root->addElement(_elem, dim);
		return 0;
	}
	if (rootLev == 0) {
		if (root->getListSize() + 1 > leafCapacity) {
			// the root will overflow after inserting this element.
			RtreeNode* u = root->split(true, (char*)_elem, dim, branchFactor, leafCapacity);
			RtreeNode* r = new RtreeNode();
			r->initialize(false, dim, branchFactor, leafCapacity);
			r->addChildNode(root, dim);
			r->addChildNode(u, dim);
			root = r;
			rootLev++;
		}
		else {
			root->addElement(_elem, dim);
		}
		return 0;
	}
	IndexablePoint* p = _elem;
	Rectangle rect;
	rect.valueList = new int[2 * dim];
	for (int i = 0; i < dim; i++) {
		rect.valueList[i] = p->coords[i];
		rect.valueList[i + dim] = p->coords[i];
	}
	int rtn = this->insertElement_SubRoutine(root, rootLev, _elem, rect);

	// release space
	rect.releaseSpace();
	return rtn;
}

int Rtree::deleteElement(IndexablePoint* _elem) {
	if (_elem == NULL || root == NULL) {
		printf("Error in Rtree::deleteElement: _elem or root is NULL!\n");
		return 1;
	}
	IndexablePoint* p = _elem;
	Rectangle rect;
	rect.setByPoint(p->coords, dim);
	RtreeNode* targetLeaf = this->findLeafForDel_SubRoutine(root, rootLev, _elem, rect);

	// === Codes for debug.
	if (targetLeaf == NULL) {
		printf("Error in Rtree::deleteElement: _elem deletion failed!\n");
	}
	// ===

	rect.releaseSpace();

	if (targetLeaf != NULL) {
		if (targetLeaf->deleteElement(_elem) < leafCapacity / FanoutRatio) {
			// targetLeaf underflow after deleting _elem.
			this->underflowHandler(targetLeaf, 0);
		}
		return 0;
	}
	return 1;
}

void Rtree::findKNN(IndexablePoint* _elem, int _topK, IndexablePoint** _targetPlace, int& _resultCnt) {
	if (root == NULL) {
		_resultCnt = 0;
		return;
	}
	priority_queue<heapElement, vector<heapElement>, greaterThanByDist> minHeap;
	_resultCnt = 0;
	IndexablePoint* pt = _elem;
	heapElement e;
	e.lev = rootLev;
	e.ptr = (char*)root;
	e.dist = root->mbr.computeMinDistToPt(pt->coords, dim);
	minHeap.push(e);

	int temp_num = 0;
	RtreeNode* temp_ptr = NULL;
	char** temp_list = NULL;
	heapElement min_e;

	while (minHeap.size() > 0 && _resultCnt < _topK) {
		min_e = minHeap.top();
		minHeap.pop();

		if (min_e.lev == -1) {
			// this is a leaf entry, i.e., a point.
			_targetPlace[_resultCnt++] = (IndexablePoint*)min_e.ptr;
		}
		else {
			// this is a tree node.

			// === Codes for debug.
			VisitedNodeNum++;
			VisitedNodeNumList[min_e.lev]++;
			// ===

			temp_ptr = (RtreeNode*)min_e.ptr;
			temp_num = temp_ptr->getListSize();
			temp_list = temp_ptr->getListPtr();
			e.lev = min_e.lev - 1;
			if (min_e.lev > 0) {
				for (int i = 0; i < temp_num; i++) {
					e.ptr = temp_list[i];
					e.dist = ((RtreeNode*)(e.ptr))->mbr.computeMinDistToPt(pt->coords, dim);
					minHeap.push(e);
				}
			}
			else {
				for (int i = 0; i < temp_num; i++) {
					e.ptr = temp_list[i];
					// Add points to heap according to selection criteria.
					e.dist = pt->computeDist(((IndexablePoint*)e.ptr)->coords, dim);
					//					e.dist = pt->getDist((IndexablePoint*) (e.ptr), dim);
					minHeap.push(e);
				}
			}
		}
	}
}

IndexablePoint* Rtree::findNN(IndexablePoint* _elem) {
	IndexablePoint* target = NULL;
	int cnt = 0;
	findKNN(_elem, 1, &target, cnt);
	return target;
}

/*
*  Return the number of points in ball B(_elem, _radius).
*/
int Rtree::rangeCntQuery(IndexablePoint* _elem, double _radius) {
	if (root == NULL) {
		printf("Warning in Rtree::rangeQuery: The root is NULL!\n");
		return 0;
	}
	RtreeNode** childList = NULL;
	int childNum = 0;
	IndexablePoint** ptList = NULL;
	IndexablePoint* q = _elem;
	int ptNum = 0;
	int state = 0;
	double sqrRadius = _radius * _radius;
	double sqrDist = 0;
	int cnt = 0;

	queue<QueElement> que;
	QueElement temp;
	QueElement curElement;
	curElement.node = root;
	curElement.lev = rootLev;
	que.push(curElement);
	while (!que.empty()) {
		curElement = que.front();
		que.pop();
		if (curElement.lev == 0) {
			// the current node is a leaf node.
			ptList = (IndexablePoint**)(curElement.node->getListPtr());
			ptNum = curElement.node->getListSize();
			for (int i = 0; i < ptNum; i++) {
				sqrDist = q->computeSqrDist(ptList[i]->coords, dim);
				if (sqrDist <= sqrRadius)
					cnt++;
			}
		}
		else {
			// the current node is an internal node.
			childList = (RtreeNode**)(curElement.node->getListPtr());
			childNum = curElement.node->getListSize();
			for (int i = 0; i < childNum; i++) {
				state = childList[i]->mbr.stateWithSphere(q->coords, _radius, dim);
				if (state == 0) {
					// intersect, then add this node to queue.
					temp.node = childList[i];
					temp.lev = curElement.lev - 1;
					que.push(temp);
				}
				else if (state == 1) {
					// fully inside, count the total number of points in the subtree of this node.
					cnt += countElementInSubtree(childList[i], curElement.lev - 1);
				}
			}
		}
	}
	return cnt;
}

/*
*  Added by jhgan on 2017-01-23.
*  The implementation of range reporting.
*/
int Rtree::rangeQuery(IndexablePoint* _elem, double _radius, vector<IndexablePoint*>& _targetPlace) {
	if (root == NULL) {
		printf("Warning in Rtree::rangeQuery: The root is NULL!\n");
		return 0;
	}
	// Clear the targetPlace.
	_targetPlace.clear();

	RtreeNode** childList = NULL;
	int childNum = 0;
	IndexablePoint** ptList = NULL;
	int ptNum = 0;
	int state = 0;
	double sqrRadius = _radius * _radius;
	double sqrDist = 0;

	queue<QueElement> que;
	QueElement temp;
	QueElement curElement;
	curElement.node = root;
	curElement.lev = rootLev;
	que.push(curElement);
	while (!que.empty()) {
		curElement = que.front();
		que.pop();

		if (curElement.lev == 0) {
			// the current node is a leaf node.
			ptList = (IndexablePoint**)(curElement.node->getListPtr());
			ptNum = curElement.node->getListSize();
			for (int i = 0; i < ptNum; i++) {
				sqrDist = _elem->computeSqrDist(ptList[i]->coords, dim);
				if (sqrDist <= sqrRadius)
					_targetPlace.push_back(ptList[i]);
			}
		}
		else {
			// the current node is an internal node.
			childList = (RtreeNode**)(curElement.node->getListPtr());
			childNum = curElement.node->getListSize();
			for (int i = 0; i < childNum; i++) {
				state = childList[i]->mbr.stateWithSphere(_elem->coords, _radius, dim);
				if (state != -1) {
					// intersect or fully inside, then add this node to queue.
					temp.node = childList[i];
					temp.lev = curElement.lev - 1;
					que.push(temp);
				}
			}
		}
	}
	return 0;
}

/*
*  Return the farthest point in B(_q, _radius).
*/
char* Rtree::rangeFarthestQuery(IndexablePoint* _elem, double _radius) {
	if (root == NULL) {
		printf("Warning in Rtree::rangeFarthestQuery: The root is NULL!\n");
		return NULL;
	}

	IndexablePoint* targetPoint = NULL;
	IndexablePoint* _q = _elem;

	RtreeNode** childList = NULL;
	int childNum = 0;
	IndexablePoint** ptList = NULL;
	int ptNum = 0;
	int state = 0;
	double sqrRadius = _radius * _radius;
	double sqrDist = 0;
	double curMax = 0;

	queue<QueElement> que;
	QueElement temp;
	QueElement curElement;
	curElement.node = root;
	curElement.lev = rootLev;
	que.push(curElement);
	while (!que.empty()) {
		curElement = que.front();
		que.pop();
		if (curElement.lev == 0) {
			// the current node is a leaf node.F
			ptList = (IndexablePoint**)(curElement.node->getListPtr());
			ptNum = curElement.node->getListSize();
			for (int i = 0; i < ptNum; i++) {
				sqrDist = _q->computeSqrDist(ptList[i]->coords, dim);
				if (sqrDist <= sqrRadius && curMax < sqrDist) {
					targetPoint = ptList[i];
					curMax = sqrDist;
				}
			}
		}
		else {
			// the current node is an internal node.
			childList = (RtreeNode**)(curElement.node->getListPtr());
			childNum = curElement.node->getListSize();
			for (int i = 0; i < childNum; i++) {
				state = childList[i]->mbr.stateWithSphere(_q->coords, _radius, dim);
				if (state != -1) {
					// intersect or fully inside, then add this node to queue.
					temp.node = childList[i];
					temp.lev = curElement.lev - 1;
					que.push(temp);
				}
			}
		}
	}
	return (char*)targetPoint;
}

int Rtree::getDim() {
	return dim;
}
bool Rtree::isValid() {
	return root != NULL;
}

void Rtree::releaseSpace() {
	if (root != NULL) {
		releaseSpace_SubRoutine(root, rootLev);
		delete root;
		root = NULL;
		dim = -1;
	}
}

/*
*  Private functions.
*/
int Rtree::insertElement_SubRoutine(RtreeNode* _subroot, int _subrootLevel, IndexablePoint* _elem,
	Rectangle& _elem_rect) {
	if (_subroot == NULL || _elem == NULL) {
		printf("Error in Rtree::insertElement_SubRoutine: _subroot or _elem is NULL!\n");
		return 1;
	}
	//	// This is a non-recursive version added by jhgan on 8 Jan, 2017.
	//	RtreeNode* subroot = _subroot;
	//	int subrootLevel = _subrootLevel;
	//	while (true) {
	//		if (subrootLevel == 0) {
	//			// subroot is a leaf node.
	//			if (subroot->getListSize() + 1 > leafCapacity) {
	//				// the _subroot will overflow after inserting this element.
	//				RtreeNode* u = subroot->split(true, (char*) _elem, dim, branchFactor, leafCapacity);
	//				this->insertNode_SubRoutine(subroot->parent, 1, u, 0);
	//			} else {
	//				subroot->addElement(_elem, dim);
	//			}
	//			break;
	//		} else {
	//			// subroot is an internal node.
	//			subroot->mbr.enlarge(_elem_rect, dim);
	//			subroot = subroot->chooseBestChildNode(_elem_rect, dim);
	//			subrootLevel--;
	//		}
	//	}
	//	return 0;

	// Version before Jan 8, 2017
	if (_subrootLevel == 0) {
		// subroot is a leaf node.
		if (_subroot->getListSize() + 1 > leafCapacity) {
			// the _subroot will overflow after inserting this element.
			RtreeNode* u = _subroot->split(true, (char*)_elem, dim, branchFactor, leafCapacity);
			this->insertNode_SubRoutine(_subroot->parent, 1, u, 0);
		}
		else {
			_subroot->addElement(_elem, dim);
		}
		return 0;
	}
	// subroot is an internal node.
	_subroot->mbr.enlarge(_elem_rect, dim);
	RtreeNode* bestChild = _subroot->chooseBestChildNode(_elem_rect, dim);
	int rtn = this->insertElement_SubRoutine(bestChild, _subrootLevel - 1, _elem, _elem_rect);
	return rtn;
}

void Rtree::insertNode_SubRoutine(RtreeNode* _subroot, int _subrootLevel, RtreeNode* _toInsNode, int _targetLevel) {
	if (_subroot == NULL) {
		printf("Error in Rtree::InsertNode_SubRoutine: subroot is NULL!\n");
		return;
	}
	if (_toInsNode == NULL) {
		printf("Error in Rtree::InsertNode_SubRoutine: toInsNode is NULL!\n");
		return;
	}

	if (_targetLevel >= _subrootLevel) {
		printf("Error in Rtree::InsertNode_SubRoutine: targetLevel >= subroot->level!\n");
		return;
	}

	//	// This is a non-recursive version added by jhgan on 8 Jan, 2017.
	//	RtreeNode* subroot = _subroot;
	//	int subrootLevel = _subrootLevel;
	//	while (true) {
	//		if (subrootLevel == _targetLevel + 1) {
	//			if (subroot->getListSize() + 1 > branchFactor) {
	//				// this subroot overflows after this insertion. It needs to split.
	//				this->overflowHandler(subroot, _toInsNode);
	//			} else {
	//				subroot->addChildNode(_toInsNode, dim);
	//			}
	//			break;
	//		} else {
	//			// Enlarge the mbr of subroot.
	//			subroot->mbr.enlarge(_toInsNode->mbr, dim);
	//			subroot = subroot->chooseBestChildNode(_toInsNode->mbr, dim);
	//			subrootLevel--;
	//		}
	//	}

	// Version before 8 Jan, 2017.
	if (_subrootLevel == _targetLevel + 1) {
		if (_subroot->getListSize() + 1 > branchFactor) {
			// this subroot overflows after this insertion. It needs to split.
			this->overflowHandler(_subroot, _toInsNode);
		}
		else {
			_subroot->addChildNode(_toInsNode, dim);
		}
	}
	else {
		// Enlarge the mbr of subroot.
		_subroot->mbr.enlarge(_toInsNode->mbr, dim);
		RtreeNode* bestChild = _subroot->chooseBestChildNode(_toInsNode->mbr, dim);
		this->insertNode_SubRoutine(bestChild, _subrootLevel - 1, _toInsNode, _targetLevel);
	}
}

void Rtree::overflowHandler(RtreeNode* _overflowNode, RtreeNode* _toInsNode) {
	RtreeNode* u = _overflowNode->split(false, (char*)_toInsNode, dim, branchFactor, leafCapacity);

	if (_overflowNode->parent != NULL) {
		if (_overflowNode->parent->getListSize() + 1 > branchFactor) {
			// the parent node overflows after inserting u.
			this->overflowHandler(_overflowNode->parent, u);
		}
		else {
			_overflowNode->parent->addChildNode(u, dim);
		}
	}
	else {
		// overflowNode is the root of this tree.
		// then create a new root for this tree.
		RtreeNode* r = new RtreeNode();
		r->initialize(false, dim, branchFactor, leafCapacity);
		r->addChildNode(_overflowNode, dim);
		r->addChildNode(u, dim);
		this->root = r;
		rootLev++;
	}
}

RtreeNode* Rtree::findLeafForDel_SubRoutine(RtreeNode* _subroot, int _subrootLevel, IndexablePoint* _elem,
	Rectangle& _elem_rect) {
	if (_subroot == NULL) {
		return NULL;
	}
	if (_subrootLevel == 0) {
		// subroot is a leaf node.
		int elemNum = _subroot->getListSize();
		IndexablePoint** elemList = (IndexablePoint**)_subroot->getListPtr();
		for (int i = 0; i < elemNum; i++) {
			if (elemList[i] == _elem)
				return _subroot;
		}
	}
	else {
		RtreeNode* rtn = NULL;
		int childNum = _subroot->getListSize();
		RtreeNode** childList = (RtreeNode**)(_subroot->getListPtr());
		for (int i = 0; i < childNum; ++i) {
			if (childList[i]->mbr.contain(_elem_rect, dim)) {
				// the mbr of this child node fully contains the mbr of _elem.
				rtn = this->findLeafForDel_SubRoutine(childList[i], _subrootLevel - 1, _elem, _elem_rect);
				if (rtn != NULL) {
					return rtn;
				}
			}
		}
	}
	return NULL;
}

void Rtree::underflowHandler(RtreeNode* _underflowNode, int _nodeLevel) {
	if (_underflowNode == NULL) {
		return;
	}
	RtreeNode* parent = _underflowNode->parent;
	if (parent == NULL) {
		// underflowNode is the root of this tree.
		int num = _underflowNode->getListSize();
		if (num > 1) {
			// underflowNode is a legal root.
			return;
		}
		else {
			if (num == 1) {
				if (rootLev > 0) {
					// underflowNode is an illegal root.
					// Let the only child node be the root.
					this->root = (RtreeNode*)(root->getListPtr()[0]);
					this->root->parent = NULL;
					rootLev = rootLev - 1;
					_underflowNode->releaseSpace();
					delete _underflowNode;
					_underflowNode = NULL;
				}
			}
			else {
				this->root = NULL;
				rootLev = 0;
				_underflowNode->releaseSpace();
				delete _underflowNode;
				_underflowNode = NULL;
			}
		}
	}
	else {
		// underflowNode is not the root of this tree.

		// Get all the child nodes of underflowNode.
		int reInsNum = _underflowNode->getListSize();
		char** listPtr = _underflowNode->getListPtr();
		char** reInsList = new char*[reInsNum];
		for (int i = 0; i < reInsNum; i++)
			reInsList[i] = listPtr[i];

		// Delete underflowNode from its parent.
		if (parent->deleteChildNode(_underflowNode) < branchFactor / FanoutRatio) {
			// parent underflows after deleting underflowNode.
			// Handle the underflowing event of parent.
			this->underflowHandler(parent, _nodeLevel + 1);
		}
		if (_nodeLevel > 0) {
			// Reinsert all the child nodes of underflowNode.
			for (int i = 0; i < reInsNum; i++) {
				this->insertNode((RtreeNode*)(reInsList[i]), _nodeLevel - 1);
			}
		}
		else {
			// Reinsert all the elements of underflowNode.
			for (int i = 0; i < reInsNum; i++) {
				//				// === Codes for debug.
				//				if(((Point*)(reInsList[i]))->id == -1){
				//					printf("time to break.\n");
				//				}
				//				// ===
				this->insertElement((IndexablePoint*)(reInsList[i]));
			}
		}
		delete[] reInsList;
	}
}

void Rtree::releaseSpace_SubRoutine(RtreeNode* _subroot, int _subrootLev) {
	char** childList = _subroot->getListPtr();
	int childNum = _subroot->getListSize();
	RtreeNode* curChildNode = NULL;
	if (_subrootLev > 0) {
		for (int i = 0; i < childNum; i++) {
			curChildNode = (RtreeNode*)(childList[i]);
			releaseSpace_SubRoutine(curChildNode, _subrootLev - 1);
			delete curChildNode;
		}
	}
	_subroot->releaseSpace();

	// version on Nov 28, 2016
	//	char** childList = _subroot->getListPtr();
	//	int childNum = _subroot->getListSize();
	//	if (_subrootLev > 0) {
	//		for (int i = 0; i < childNum; i++) {
	//			releaseSpace_SubRoutine((RtreeNode*) (childList[i]), _subrootLev - 1);
	//		}
	//	}
	//	_subroot->releaseSpace();
}

/*
*  Count the total number of elements in the specified subtree.
*/
int Rtree::countElementInSubtree(RtreeNode* _subroot, int _subrootLevel) {
	int cnt = 0;
	if (_subrootLevel == 0) {
		// This is a leaf node.
		cnt += _subroot->getListSize();
		return cnt;
	}
	// This is an internal node.
	int childNum = _subroot->getListSize();
	RtreeNode** childList = (RtreeNode**)(_subroot->getListPtr());
	for (int i = 0; i < childNum; ++i) {
		cnt += countElementInSubtree(childList[i], _subrootLevel - 1);
	}
	return cnt;
}

/***************************************************
*  Added by jhgan on 2017-01-19.
*  STR bulk loading related.
**************************************************/
/*
*  Construct leaf level by STR algorithm.
*  Parameter List:
*  	_leafLevel:			the target place to store the pointers of the leaf nodes.
*  	_dimLev:			which dimension we used for sorting in the current recursion.
*/
int Rtree::constructLeafLevel_STR_SubRoutine(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
	int _leafCapacity, vector<RtreeNode*>& _leafLevel, int _dimLev) {
	//	// === Codes for debug.
	//	const int print = 0;
	//	// ===

	if (_elementList == NULL) {
		printf("Error in Rtree::constructLeafLevel_STR_SubRoutine: _elementList == NULL!\n");
		return 1;
	}
	if (_elemNum <= 0) {
		printf("Error in Rtree::constructLeafLevel_STR_SubRoutine: _elemNum <= 0!\n");
		return 1;
	}

	double exp = 1.0 / (_dim - _dimLev + 1);
	int slabNum = ceil(pow((_elemNum * 1.0 / _leafCapacity), exp));
	int slabSize = ceil(_elemNum * 1.0 / slabNum); // The number of points in each slab.
												   //	if (print) {
												   //		printf("dimLev = %d: _nodeNum = %d:\t slabNum = %d\t slabSize = %d\n", _dimLev, _elemNum, slabNum, slabSize);
												   //	}

	if (_elemNum <= _leafCapacity) {
		// There is only one node.
		RtreeNode* curNode = new RtreeNode();
		curNode->initialize(true, _dim, _branchFactor, _leafCapacity);
		for (int i = 0; i < _elemNum; ++i) {
			curNode->addElement(_elementList[i], _dim);
		}
		_leafLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			curNode->mbr.showRectangle(_dim);
		////			printf("\n");
		//		}
		//		// ===
		return 0;
	}

	if (slabSize <= _leafCapacity / FanoutRatio || _dimLev == _dim) {
		// There are too few points to make slabs. Or this is the last level in the recursion.
		// Sort all the points by the _dimLev-dim coordinates.
		sort(_elementList, _elementList + _elemNum, PointDimSorter(_dim, _dimLev));

		// Scan the sorted list and construct tree nodes.
		int nodeNum = ceil(_elemNum * 1.0 / _leafCapacity); // Note that nodeNum > 1.
		RtreeNode* curNode = NULL;

		int curIndex = 0;
		int tempNum = _leafCapacity;
		// We construct the first (nodeNum - 2) nodes because we may need to deal with the last two nodes specially.
		for (int i = 0; i < nodeNum - 2; ++i) {
			curNode = new RtreeNode();
			curNode->initialize(true, _dim, _branchFactor, _leafCapacity);
			for (int j = 0; j < tempNum; ++j) {
				curNode->addElement(_elementList[curIndex], _dim);
				curIndex++;
			}
			_leafLevel.push_back(curNode);
			//			// === Codes for debug.
			//			if (print) {
			//				curNode->mbr.showRectangle(_dim);
			////			printf("\n");
			//			}
			//			// ===
		}
		// We deal with the last two nodes specially.
		tempNum = (_elemNum - curIndex) / 2;
		curNode = new RtreeNode();
		curNode->initialize(true, _dim, _branchFactor, _leafCapacity);
		for (int j = 0; j < tempNum; ++j) {
			curNode->addElement(_elementList[curIndex], _dim);
			curIndex++;
		}
		_leafLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			curNode->mbr.showRectangle(_dim);
		////		printf("\n");
		//		}
		//		// ===

		tempNum = _elemNum - curIndex;
		curNode = new RtreeNode();
		curNode->initialize(true, _dim, _branchFactor, _leafCapacity);
		for (int j = 0; j < tempNum; ++j) {
			curNode->addElement(_elementList[curIndex], _dim);
			curIndex++;
		}
		_leafLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			curNode->mbr.showRectangle(_dim);
		//			printf("\n");
		//		}
		//		// ===
		return 0;
	}
	else {
		// This is not the last level in the recursion and there are enough points in each slab.
		// Make slabs.
		// Note that: There must be at least two slabs because _leafCapacity < _elemNum.
		// Sort all the points by the _dimLev-dim coordinates.
		sort(_elementList, _elementList + _elemNum, PointDimSorter(_dim, _dimLev));
		int rtn = 0;
		IndexablePoint** curPtr = _elementList;
		for (int i = 0; i < slabNum - 2; ++i) {
			rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,
				_dimLev + 1);
			if (rtn == 1) {
				printf("Error in Rtree::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
				return 1;
			}
			curPtr += slabSize;
		}
		// Process the last two slabs specially.
		slabSize = (_elemNum - (curPtr - _elementList)) / 2;
		rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,
			_dimLev + 1);
		if (rtn == 1) {
			printf("Error in Rtree::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
			return 1;
		}
		curPtr += slabSize;

		// Put all the rest points to the last slab.
		slabSize = _elemNum - (curPtr - _elementList);
		rtn = constructLeafLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _leafLevel,
			_dimLev + 1);
		if (rtn == 1) {
			printf("Error in Rtree::constructLeafLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
			return 1;
		}
		curPtr += slabSize;
		return 0;
	}
	return 0;
}
int Rtree::constructNextLevel_STR_SubRoutine(RtreeNode** _nodeList, int _nodeNum, int _dim, int _branchFactor,
	int _leafCapacity, vector<RtreeNode*>& _nextLevel, int _dimLev) {
	//	// === Codes for debug.
	//	const int print = 0;
	//	// ===

	if (_nodeList == NULL) {
		printf("Error in Rtree::constructNextLevel_STR_SubRoutine: _nodeList == NULL!\n");
		return 1;
	}
	if (_nodeNum <= 0) {
		printf("Error in Rtree::constructNextLevel_STR_SubRoutine: _nodeNum <= 0!\n");
		return 1;
	}

	double exp = 1.0 / (_dim - _dimLev + 1);
	int slabNum = ceil(pow((_nodeNum * 1.0 / _branchFactor), exp));
	int slabSize = ceil(_nodeNum * 1.0 / slabNum); // The number of points in each slab.
												   //	// === Codes for debug.
												   //	if (print) {
												   //		printf("dimLev = %d: _nodeNum = %d:\t slabNum = %d\t slabSize = %d\n", _dimLev, _nodeNum, slabNum, slabSize);
												   //	}
												   //	// ===

	if (_nodeNum <= _branchFactor) {
		// There is only one node.
		RtreeNode* curNode = new RtreeNode();
		curNode->initialize(false, _dim, _branchFactor, _leafCapacity);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//		}
		//		// ===
		for (int i = 0; i < _nodeNum; ++i) {
			curNode->addChildNode(_nodeList[i], _dim);
			//			// === Codes for debug.
			//			if (print) {
			//				_nodeList[i]->mbr.showRectangle(_dim);
			//			}
			//			// ===
		}
		_nextLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//			curNode->mbr.showRectangle(_dim);
		////			printf("\n");
		//		}
		//		// ===
		return 0;
	}

	if (slabSize <= _branchFactor / FanoutRatio || _dimLev == _dim) {
		// There are too few points to make slabs. Or this is the last level in the recursion.
		// Sort all the points by the _dimLev-dim coordinates.
		sort(_nodeList, _nodeList + _nodeNum, NodeDimSorter(_dim, _dimLev));
		// Scan the sorted list and construct tree nodes.
		int nextNodeNum = ceil(_nodeNum * 1.0 / _branchFactor); // Note that nextNodeNum > 1.
		RtreeNode* curNode = NULL;

		int curIndex = 0;
		int tempNum = _branchFactor;
		// We construct the first (nextNodeNum - 2) nodes because we may need to deal with the last two nodes specially.
		for (int i = 0; i < nextNodeNum - 2; ++i) {
			curNode = new RtreeNode();
			curNode->initialize(false, _dim, _branchFactor, _leafCapacity);
			for (int j = 0; j < tempNum; ++j) {
				curNode->addChildNode(_nodeList[curIndex], dim);
				curIndex++;
			}
			_nextLevel.push_back(curNode);
			//			// === Codes for debug.
			//			if (print) {
			//				curNode->mbr.showRectangle(_dim);
			////			printf("\n");
			//			}
			//			// ===
		}
		// We deal with the last two nodes specially.
		tempNum = (_nodeNum - curIndex) / 2;
		curNode = new RtreeNode();
		curNode->initialize(false, _dim, _branchFactor, _leafCapacity);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//		}
		//		// ===

		for (int j = 0; j < tempNum; ++j) {
			//			// === Codes for debug.
			//			if (print) {
			//				_nodeList[curIndex]->mbr.showRectangle(_dim);
			//			}
			//			// ===
			curNode->addChildNode(_nodeList[curIndex], dim);
			curIndex++;
		}
		_nextLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//			curNode->mbr.showRectangle(_dim);
		////			printf("\n");
		//		}
		//		// ===

		tempNum = _nodeNum - curIndex;
		curNode = new RtreeNode();
		curNode->initialize(false, _dim, _branchFactor, _leafCapacity);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//		}
		//		// ===
		for (int j = 0; j < tempNum; ++j) {
			//			// === Codes for debug.
			//			if (print) {
			//				_nodeList[curIndex]->mbr.showRectangle(_dim);
			//			}
			//			// ===
			curNode->addChildNode(_nodeList[curIndex], _dim);
			curIndex++;
		}
		_nextLevel.push_back(curNode);
		//		// === Codes for debug.
		//		if (print) {
		//			printf("==========================\n");
		//			curNode->mbr.showRectangle(_dim);
		//			printf("\n");
		//		}
		//		// ===

		return 0;
	}
	else {
		// This is not the last level in the recursion and there are enough points in each slab.
		// Make slabs.
		// Note that: There must be at least two slabs because _branchFactor < _nodeNum.
		// Sort all the points by the _dimLev-dim coordinates.
		sort(_nodeList, _nodeList + _nodeNum, NodeDimSorter(_dim, _dimLev));
		int rtn = 0;
		RtreeNode** curPtr = _nodeList;
		for (int i = 0; i < slabNum - 2; ++i) {
			rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel,
				_dimLev + 1);
			if (rtn == 1) {
				printf("Error in Rtree::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
				return 1;
			}
			curPtr += slabSize;
		}
		// Process the last two slabs specially.
		slabSize = (_nodeNum - (curPtr - _nodeList)) / 2;
		rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel,
			_dimLev + 1);
		if (rtn == 1) {
			printf("Error in Rtree::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
			return 1;
		}
		curPtr += slabSize;

		// Put all the rest points to the last slab.
		slabSize = _nodeNum - (curPtr - _nodeList);
		rtn = constructNextLevel_STR_SubRoutine(curPtr, slabSize, _dim, _branchFactor, _leafCapacity, _nextLevel,
			_dimLev + 1);
		if (rtn == 1) {
			printf("Error in Rtree::constructNextLevel_STR_SubRoutine: Something must be wrong in recursion!\n");
			return 1;
		}
		curPtr += slabSize;
		return 0;
	}
	return 0;
}

/*
* 	The STR bulk loading algorithm.
*/
RtreeNode* Rtree::constructRtree_STR(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
	int _leafCapacity) {
	if (root != NULL) {
		printf("Error in Rtree::constructRtree_STR: The tree is not empty!\n");
		return NULL;
	}
	if (_branchFactor < FanoutRatio || _leafCapacity < FanoutRatio) {
		printf("Error in Rtree::constructRtree_STR: _branchFactor and _leafCapacity should be >= %d!\n", FanoutRatio);
		return NULL;
	}

	dim = _dim;
	branchFactor = _branchFactor;
	leafCapacity = _leafCapacity;
	rootLev = 0;

	if (_elementList == NULL) {
		printf("Error in Rtree::constructRtree_STR: _elementList == NULL!\n");
		return NULL;
	}
	if (_elemNum <= 0) {
		printf("Error in Rtree::constructRtree_STR: _elemNum <= 0!\n");
		return NULL;
	}

	IndexablePoint** ptList = new IndexablePoint*[_elemNum];
	for (int i = 0; i < _elemNum; i++) {
		ptList[i] = _elementList[i];
	}
	//	// === Codes for debug.
	//	double exp = 1.0 / (_dim);
	//	//	int slabSize = pow(_nodeNum * _branchFactor, exp); // The number of points in each slab.
	//	int slabNum = ceil(pow(_elemNum * 1.0 / _leafCapacity, exp));
	//	int slabSize = ceil(_elemNum * 1.0 / slabNum); // The number of points in each slab.
	//	printf("_nodeNum = %d:\t slabNum = %d\t slabSize = %d\n", _elemNum, slabNum, slabSize);
	//	// ===

	int rtn = 0;
	vector<RtreeNode*> curNodeList;
	rtn = constructLeafLevel_STR_SubRoutine(ptList, _elemNum, _dim, _branchFactor, _leafCapacity, curNodeList, 1);
	if (rtn == 1) {
		printf("Error in Rtree::constructRtree_STR: Leaf level construction failed!\n");
		curNodeList.clear();
		delete[] ptList;
		ptList = NULL;
		return NULL;
	}
	delete[] ptList;
	ptList = NULL;

	int level = 0;
	vector<RtreeNode*> nextNodeList;
	RtreeNode** curNodeListPtr = NULL;
	int curNodeNum = 0;
	// Start to construct internal levels.
	while (curNodeList.size() > 1) {
		nextNodeList.clear();
		curNodeListPtr = curNodeList.data();
		curNodeNum = curNodeList.size();
		//		// === Codes for debug.
		//		double exp = 1.0 / (_dim);
		//		//	int slabSize = pow(_nodeNum * _branchFactor, exp); // The number of points in each slab.
		//		int slabNum = ceil(pow(curNodeNum * 1.0 / _branchFactor, exp));
		//		int slabSize = ceil(curNodeNum * 1.0 / slabNum); // The number of points in each slab.
		//		printf("_nodeNum = %d:\t slabNum = %d\t slabSize = %d\n", curNodeNum, slabNum, slabSize);
		//		// ===
		rtn = constructNextLevel_STR_SubRoutine(curNodeListPtr, curNodeNum, _dim, _branchFactor, _leafCapacity,
			nextNodeList, 1);
		if (rtn == 1) {
			printf("Error in Rtree::constructRtee_STR: Next level construction failed!\n");
			curNodeList.clear();
			nextNodeList.clear();
			return NULL;
		}
		curNodeList.swap(nextNodeList);
		level++;
	}
	root = curNodeList[0];
	rootLev = level;

	curNodeList.clear();
	nextNodeList.clear();
	return root;
}

/*
*  Functions for debug.
*/
void Rtree::generateRtreePlotFile(char* _filename, int _targetLev) {
	RtreeNode* curNode = this->root;
	int curLev = rootLev;
	int preLevSize = 1;
	int nextLevSize = 0;
	if (curNode == NULL) {
		// TODO error message.
		printf("Error in rtree::generateRtreePlotFile: The root is NULL!\n");
		return;
	}

	RtreeNode** childList = NULL;
	int childNum = 0;
	int objID = 1;

	FILE* file = fopen(_filename, "wt");
	
	if (file == NULL) {
		printf("Error in rtree::generateRtreePlotFile:\nResult file open failed! Please check your file path!\n");
		return;
	}

	fprintf(file, "reset\n");
	fprintf(file, "set terminal postscript eps enhanced \"Times-Roman\" 35\n");
	fprintf(file, "set out \"%s.eps\"\n", _filename);

	queue<RtreeNode*> que;
	que.push(curNode);
	while (!que.empty()) {
		curNode = que.front();
		que.pop();
		preLevSize--;

		if (curLev == _targetLev) {
			fprintf(file, "set object %d rect from ", objID++);
			fprintf(file, "%d", curNode->mbr.valueList[0]);
			for (int i = 1; i < dim; ++i) {
				fprintf(file, ", %d", curNode->mbr.valueList[i]);
			}
			fprintf(file, " to %d", curNode->mbr.valueList[dim]);
			for (int i = 1; i < dim; ++i) {
				fprintf(file, ", %d", curNode->mbr.valueList[i + dim]);
			}
			fprintf(file, " fs empty border lc rgb \'red\'\n");
		}
		if (curLev > 0) {
			childNum = curNode->getListSize();
			childList = (RtreeNode**)(curNode->getListPtr());
			for (int i = 0; i < childNum; i++) {
				que.push(childList[i]);
				nextLevSize++;
			}
		}
		else {
			nextLevSize += curNode->getListSize();
		}
		if (preLevSize == 0) {
			curLev--;
			preLevSize = nextLevSize;
			printf("lev = %d has %d nodes.\n", curLev, preLevSize);
			fflush(stdout);
			nextLevSize = 0;
		}
	}

	fprintf(file, "set xrange [0:100000]\n");
	fprintf(file, "set yrange [0:100000]\n");
	fprintf(file, "plot x\n");
	//	printf("e\n");
	//	fprintf(file, "e\n");
	fclose(file);
}

int Rtree::getHeight() {
	return rootLev + 1;
}

