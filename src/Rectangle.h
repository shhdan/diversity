/*
* rectangle.h
*
*  Created on: 11/07/2016
*      Author: uqjgan
*/

#ifndef RECTANGLE_H_
#define RECTANGLE_H_

#include <vector>
using namespace std;

class Rectangle {
public:
	/*
	*  The list of min values and then the max values.
	*  The length of this list is 2 * dim.
	*/
	int* valueList;

	Rectangle();
	Rectangle(int* _valueList);

	~Rectangle();

	bool intersect(const Rectangle& _rec, int _dim) const;
	bool inside(const Rectangle& _rec, int _dim) const;
	bool contain(const Rectangle& _rec, int _dim) const;
	long long perimeterIncrement(const Rectangle& _rec, int _dim) const;
	long long computePerimeter(int _dim) const;
	Rectangle& enlarge(const Rectangle& _rec, int _dim);
	void setByRectangle(const Rectangle& _rec, int _dim);
	void setByPoint(const int* _pointCoords, int _dim);
	void setCoords(const int* _valueList, int _dim);
	bool checkStab(vector<int>& coords, int _dim);

	double computeMinDistToPt(const int* _coords, int _dim);

	/*
	*  Compute the state of this rectangle with the specified sphere.
	*  Return:
	*  	1:		this rectangle is fully inside the sphere.
	*  	0:		this rectangle intersects with the sphere.
	*  	-1:		this rectangle is fully outside the sphere.
	*/
	int stateWithSphere(const int* _centerCoords, double _radius, int _dim);

	void showRectangle(int _dim);
	void releaseSpace();
};

#endif /* RECTANGLE_H_ */
