//#include "stdafx.h"
#include "InvDIX.h"


InvDIX::InvDIX()
{
	this->PointId = -1;
	vector<IndexablePoint*> trajectorylist = vector<IndexablePoint*>();
	this->trajectoryArray = trajectorylist;
}

void InvDIX::addTrajectoryToArray(IndexablePoint* point) {
	this->trajectoryArray.push_back(point);
}

void InvDIX::set_invertedlistId(int id) {
	this->PointId = id;
}

int InvDIX::get_invertedlistId() {
	return this->PointId;
}

vector<IndexablePoint*> InvDIX::get_trajectoryArray() {
	return this->trajectoryArray;
}

void InvDIX::set_trajectoryArray(vector<IndexablePoint*> trajectorylist) {
	this->trajectoryArray = trajectorylist;
}

InvDIX::~InvDIX()
{
	this->trajectoryArray.clear();
}
