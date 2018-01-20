/*
* IndexablePoint.cpp
*
*  Created on: 14/09/2016
*      Author: uqjgan
*/
/* added by SHD*/
//#include "stdafx.h"

#include "IndexablePoint.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

IndexablePoint::IndexablePoint() {
	id = 0;
	coords = NULL;
}

IndexablePoint::IndexablePoint(int _id, int* _coords) {
	id = _id;
	coords = _coords;
}

IndexablePoint::~IndexablePoint() {
	// TODO Auto-generated destructor stub
	//	releaseSpace();
}

void IndexablePoint::releaseSpace() {
	if (coords != NULL) {
		delete[] coords;
		coords = NULL;
	}
}

double IndexablePoint::computeSqrDist(int* _coords, int _dim) {
	double sqrDist = 0;
	long long temp = 0;
	for (int i = 0; i < _dim; ++i) {
		temp = coords[i] - _coords[i];
		sqrDist += temp * temp;
	}
	return sqrDist;
}

double IndexablePoint::computeDist(int* _coords, int _dim) {
	return sqrt(computeSqrDist(_coords, _dim));
}

long long IndexablePoint::computeDotProduct(int* _coords, int _dim) {
	long long ans = 0;
	long long temp = 0;
	for (int i = 0; i < _dim; ++i) {
		temp = _coords[i];
		temp *= coords[i];
		ans += temp;
	}
	return ans;
}

long long IndexablePoint::computeDotProduct(long long* _coords, int _dim) {
	long long ans = 0;
	for (int i = 0; i < _dim; ++i)
		ans += _coords[i] * coords[i];
	return ans;
}

