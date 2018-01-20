#pragma once
#include<vector>
#include "Point.h"
using namespace std;

class Trajectory
{

private:
	int trajectoryId;
	vector<IndexablePoint*> pointArray;
	int pointArray_size;

public:
	Trajectory();
	void addPointToArray(IndexablePoint* point);
	void set_trajectoryId(int id);
	int get_trajectoryId();
	vector<IndexablePoint*> get_pointArray();
	void set_pointArray(vector<IndexablePoint*> pointArray);
	double EuclideanDistance(Trajectory* Tra, int _dim);
	double EuclideanDistanceODPoint(IndexablePoint * opoint, IndexablePoint * dpoint, int _dim);
	~Trajectory();
};
