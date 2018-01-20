/*
* rectangle.cpp
*
*  Created on: 11/07/2016
*      Author: uqjgan
*/
/* added by SHD*/
//#include "stdafx.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Rectangle.h"
using namespace std;

Rectangle::Rectangle() {
	// TODO Auto-generated constructor stub
	//	this->valueList = new int[2 * dim];
	this->valueList = NULL;
}

Rectangle::Rectangle(int* _valueList) {
	this->valueList = _valueList;
}

Rectangle::~Rectangle() {
	// TODO Auto-generated destructor stub
}

/*
*  Detect whether this rectangle intersects with the given rectangle *_rec*.
*  Return:
*  	If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::intersect(const Rectangle& _rec, int _dim) const {
	for (int i = 0; i < _dim; i++) {
		if (this->valueList[i] > _rec.valueList[i + _dim]) {
			return false;
		}
		if (this->valueList[i + _dim] < _rec.valueList[i]) {
			return false;
		}
	}
	return true;
}

/*
*  Detect whether this rectangle is fully inside the given rectangle *_rec*.
*  Return:
*  	If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::inside(const Rectangle& _rec, int _dim) const {
	for (int i = 0; i < _dim; i++) {
		if (this->valueList[i] < _rec.valueList[i]) {
			return false;
		}
		if (this->valueList[i + _dim] > _rec.valueList[i + _dim]) {
			return false;
		}
	}
	return true;
}

/*
*  Detect whether the given rectangle *_rec* is fully contained in this rectangle.
*  Return:
*  	If YES, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::contain(const Rectangle& _rec, int _dim) const {
	for (int i = 0; i < _dim; ++i) {
		if (_rec.valueList[i] < this->valueList[i]) {
			return false;
		}
		if (_rec.valueList[i + _dim] > this->valueList[i + _dim]) {
			return false;
		}
	}
	return true;
}

/*
*  Compute the increment of the "perimeter" if the given rectangle is merged with this rectangle.
*  Parameter List:
*  	_rec:		the given rectangle.
*  Return:
*  	the sum of increments of the length in each dimension.
*/
long long Rectangle::perimeterIncrement(const Rectangle& _rec, int _dim) const {
	long long rtn = 0;
	for (int i = 0; i < _dim; ++i) {
		rtn += (max(this->valueList[i + _dim], _rec.valueList[i + _dim]) - min(this->valueList[i], _rec.valueList[i]));
	}
	//	rtn *= 2;
	return rtn - (this->computePerimeter(_dim));
}

/*
*  Compute the "perimeter" of this rectangle.
*  Return:
*  	the sum of the length in each dimension.
*/
long long Rectangle::computePerimeter(int _dim) const {
	long long rtn = 0;
	for (int i = 0; i < _dim; ++i) {
		rtn += (this->valueList[i + _dim] - this->valueList[i]);
	}
	//	rtn *= 2;
	return rtn;
}

/*
*  Find the minimum bounding rectangle of this and the given rectangles.
*  Parameter List:
*  	_rec:			the given rectangle
*  Return:
*  	the resulting MBR
*/
Rectangle& Rectangle::enlarge(const Rectangle& _rec, int _dim) {
	for (int i = 0; i < _dim; i++) {
		this->valueList[i] = min(this->valueList[i], _rec.valueList[i]);
		this->valueList[i + _dim] = max(this->valueList[i + _dim], _rec.valueList[i + _dim]);
	}
	return *(this);
}

void Rectangle::setByPoint(const int* _pointCoords, int _dim) {
	if (this->valueList == NULL) {
		this->valueList = new int[2 * _dim];
	}
	for (int i = 0; i < _dim; i++) {
		this->valueList[i] = _pointCoords[i];
		this->valueList[i + _dim] = _pointCoords[i];
	}
}

/*
*  Set this rectangle by a given rectangle.
*/
void Rectangle::setByRectangle(const Rectangle& _rec, int _dim) {
	this->setCoords(_rec.valueList, _dim);
}
/*
*  Set the list of coordinate values.
*/
void Rectangle::setCoords(const int* _valueList, int _dim) {
	if (this->valueList == NULL) {
		this->valueList = new int[2 * _dim];
	}
	for (int i = 0; i < 2 * _dim; i++) {
		this->valueList[i] = _valueList[i];
	}
}

/*
*  Check whether this rectangle is stabbed by the given point.
*  Parameter List:
*  	coords:		the coordinates of the given point.
*  Return:
*  	If Yes, return TRUE. Otherwise, return FALSE.
*/
bool Rectangle::checkStab(vector<int>& coords, int _dim) {
	for (int i = 0; i < _dim; ++i) {
		if (coords[i] < this->valueList[i] || coords[i] > this->valueList[i + _dim])
			return false;
	}
	return true;
}

/*
*  Compute the minimum L2 distant to the given point.
*/
double Rectangle::computeMinDistToPt(const int* _coords, int _dim) {
	double dist = 0;
	double temp = 0;
	for (int i = 0; i < _dim; i++) {
		temp = 0;
		if (_coords[i] < valueList[i]) {
			temp = valueList[i] - _coords[i];
		}
		else if (_coords[i] > valueList[i + _dim]) {
			temp = _coords[i] - valueList[i + _dim];
		}
		dist += temp * temp;
	}
	return sqrt(dist);
}


/*
*  Compute the state of this rectangle with the specified sphere.
*  Return:
*  	1:		this rectangle is fully inside the sphere.
*  	0:		this rectangle intersects with the sphere.
*  	-1:		this rectangle is fully outside the sphere.
*/
int Rectangle::stateWithSphere(const int* _centerCoords, double _radius, int _dim) {
	int dim = _dim;
	long long closestDist = 0;
	long long farthestDist = 0;
	long long temp = 0;
	long long temp2 = 0;
	long long sqr_temp = 0;
	long long sqr_temp2 = 0;
	double sqr_r = _radius * _radius;

	int* maxList = valueList + dim;

	// Find the distances from the closest and farthest points to q in this grid cell.
	for (int i = 0; i < dim; i++) {
		temp = valueList[i] - _centerCoords[i];
		//		temp2 = valueList[i + dim] - _centerCoords[i];
		temp2 = maxList[i] - _centerCoords[i];
		sqr_temp = temp * temp;
		sqr_temp2 = temp2 * temp2;
		if (temp > 0) {
			// q is to the left of this rectangle in this dimension
			closestDist += sqr_temp;
		}
		else if (temp2 < 0) {
			// q is to the right of this rectangle in this dimension
			closestDist += sqr_temp2;
		}
		farthestDist += (sqr_temp <= sqr_temp2 ? sqr_temp2 : sqr_temp);
	}
	if (closestDist <= sqr_r) {
		if (farthestDist <= sqr_r)
			return 1; // fully inside
		return 0; // intersect
	}
	return -1; // fully outside
}

/*
*  Show the coordinates of this rectangle.
*/
void Rectangle::showRectangle(int _dim) {
	// print the range
	printf("The rectangle is bounded by [%d, %d]", this->valueList[0], this->valueList[_dim]);
	for (int j = 1; j < _dim; ++j) {
		printf(" * [%d, %d]", this->valueList[j], this->valueList[j + _dim]);
	}
	printf("\n");
}

/*
*  Release the space of this rectangle.
*/
void Rectangle::releaseSpace() {
	if (this->valueList != NULL) {
		delete[](this->valueList);
		this->valueList = NULL;
	}
}
