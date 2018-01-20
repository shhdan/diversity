//#include "stdafx.h"
#include "Point.h"

Point::Point() {
	id = 0;
	coords = nullptr;
	nDimensions = 0;
	//invertedlist = nullptr;
}

Point::~Point() {
	if (coords != nullptr) {
		delete[] coords;
		coords = nullptr;
	}
	//if (invertedlist != nullptr) {
	//	delete[] invertedlist;
	//	invertedlist = nullptr;
	//}
}
//Point::Point() {
//	this->nDimensions = 2;
//	this->coordinate.reserve(nDimensions);
//	this->coordinate[0] = coordinate[1] = 0;
//}
//
//Point::Point(int Dimensions) {
//	this->nDimensions = Dimensions;
//	this->coordinate.reserve(nDimensions);
//	for (int i = 0; i < nDimensions; i++) {
//		this->coordinate[i] = 0;
//	}
//}
//
long Point::getM_coordinate(int _dim) {
	return this->coords[_dim];
}

int Point::getM_nDimensions() {
	return this->nDimensions;
}
//
//void Point::setM_coordinate(int nth, long value) {
//	this->coordinate[nth] = value;
//}
//
//double Point::PointDistance(Point point) {
//	double distance = 0.0;
//	long x1, x2, y1, y2;
//	x1 = this->getM_coordinate(1);
//	y1 = this->getM_coordinate(2);
//	x2 = point.getM_coordinate(1);
//	y2 = point.getM_coordinate(2);
//	distance = sqrt((double)((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
//	return distance;
//}
//
//Point::~Point()
//{
//	this->coordinate.clear();
//}
