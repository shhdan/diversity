#pragma once
#include<vector>
#include "math.h"
#include "IndexablePoint.h"
using namespace std;

class Point: public IndexablePoint
{

private:
	int nDimensions;
	//vector<long> coordinate;
	//int* invertedlist;
public:
	Point();
	//Point(int Dimensions);
	long getM_coordinate(int nth);
	int getM_nDimensions();
	//void setM_coordinate(int nth, long value);
	//double PointDistance(Point point);

	~Point();





};

