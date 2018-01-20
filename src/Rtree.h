#pragma once
/*
* rtree.h
*
*  Created on: 12/07/2016
*      Author: uqjgan
*/

#ifndef RTREE_H_
#define RTREE_H_

#include "RtreeNode.h"
//#include "../MyVector.h"
#include<vector>

// === Codes for debug.
extern long long VisitedNodeNum;
extern long long* VisitedNodeNumList;
// ===

class Rtree {
protected:
	RtreeNode* root; /* Leaf nodes are at level-0 and elements are at level-(-1). */
	int rootLev;

	int dim;
	int branchFactor;
	int leafCapacity;

private:
	void insertNode_SubRoutine(RtreeNode* _subroot, int _subrootLevel, RtreeNode* _toInsNode, int _targetLevel);
	int insertElement_SubRoutine(RtreeNode* _subroot, int _subrootLevel, IndexablePoint* _elem, Rectangle& _elem_rect);

	RtreeNode* findLeafForDel_SubRoutine(RtreeNode* _subroot, int _subrootLevel, IndexablePoint* _elem,
		Rectangle& _elem_rect);

	void overflowHandler(RtreeNode* _overflowNode, RtreeNode* _toInsNode);
	void underflowHandler(RtreeNode* _underflowNode, int _nodeLevel);

	void releaseSpace_SubRoutine(RtreeNode* _subroot, int _subrootLev);

	/*
	*  Count the total number of elements in the specified subtree.
	*/
	int countElementInSubtree(RtreeNode* _subroot, int _subrootLevel);

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
	int constructLeafLevel_STR_SubRoutine(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
		int _leafCapacity, vector<RtreeNode*>& _leafLevel, int _dimLev);
	int constructNextLevel_STR_SubRoutine(RtreeNode** _nodeList, int _nodeNum, int _dim, int _branchFactor,
		int _leafCapacity, vector<RtreeNode*>& _nextLevel, int _dimLev);

public:
	Rtree();
	Rtree(int _dim, int _branchFactor, int _leafCapacity);
	virtual ~Rtree();

	void initialize(int _dim, int _branchFactor, int _leafCapacity);
	void constructEmptyTree(int _dim, int _branchFactor, int _leafCapacity);

	RtreeNode* constructRtree(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
		int _leafCapacity);
	/*
	* 	Added by jhgan on 2017-01-19.
	* 	The STR bulk loading algorithm.
	*/
	RtreeNode* constructRtree_STR(IndexablePoint** _elementList, int _elemNum, int _dim, int _branchFactor,
		int _leafCapacity);

	int insertNode(RtreeNode* _toInsNode, int _targetLevel);
	int insertElement(IndexablePoint* _elem);
	int deleteElement(IndexablePoint* _elem);

	void findKNN(IndexablePoint* _elem, int _topK, IndexablePoint** _targetPlace, int& _resultCnt);
	IndexablePoint* findNN(IndexablePoint* _elem);

	/*
	*  Return the number of points in ball B(_elem, _radius).
	*/
	int rangeCntQuery(IndexablePoint* _elem, double _radius);

	/*
	*  Added by jhgan on 2017-01-23.
	*  The implementation of range reporting.
	*/
	int rangeQuery(IndexablePoint* _elem, double _radius, vector<IndexablePoint*>& _targetPlace);

	/*
	*  Return the farthest point in B(_q, _radius).
	*/
	char* rangeFarthestQuery(IndexablePoint* _q, double _radius);

	int getDim();
	bool isValid();
	void releaseSpace();

	// Functions for debug.
	void generateRtreePlotFile(char* _filename, int _targetLev);
	int getHeight();
};

#endif /* RTREE_H_ */

