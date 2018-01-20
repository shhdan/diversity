/*
* IndexablePoint.h
*
*  Created on: 14/09/2016
*      Author: uqjgan
*/

#ifndef INDEXABLEPOINT_H_
#define INDEXABLEPOINT_H_

/*
*  A base class for the indexable points such as in R-tree and other data structures.
*/
class IndexablePoint {
public:
	int id;
	/*
	*  The coordinates of this point.
	*/
	int* coords;
	int tra_id;
	
public:
	IndexablePoint();
	IndexablePoint(int _id, int* _coords);
	IndexablePoint(int _id, int * _coords, int _dim);
	~IndexablePoint();
	void releaseSpace();

	double computeSqrDist(int* _coords, int _dim);
	double computeDist(int* _coords, int _dim);
	long long computeDotProduct(int* _coords, int _dim);
	long long computeDotProduct(long long* _coords, int _dim);
};

#endif /* INDEXABLEPOINT_H_ */
