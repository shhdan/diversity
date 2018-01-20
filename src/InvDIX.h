#pragma once
#include<vector>
#include "Point.h"
using namespace std;

class InvDIX
{
private:
	int PointId;
	vector<IndexablePoint*> trajectoryArray;
public:
	InvDIX();
	void addTrajectoryToArray(IndexablePoint* point);
	void set_invertedlistId(int id);
	int get_invertedlistId();
	vector<IndexablePoint*> get_trajectoryArray();
	void set_trajectoryArray(vector<IndexablePoint*> trajectorylist);
	~InvDIX();
};

