#pragma once
#include<vector>
#include<unordered_map>
#include"Diversification.h"
#include"Trajectory.h"
#include"IndexablePoint.h"
#include"InvDIX.h"
#include"Rtree.h"
#include"stdio.h"

using namespace std;

class Util
{
public:
	Util();

	void UpperBoundTopkDiversifiedSearch(vector<vector<IndexablePoint*>> odpairs, double _oradius, double _dradius, int k, long st, long et, Rtree* rtree, vector<InvDIX*> invertedindex, vector<Trajectory*> TraList, vector<Diversification*> resultpairs);

	void BruteForceTopkDiversifiedSearch(vector<vector<IndexablePoint*>> odpairs, double _oradius, double _dradius, int k, long st, long et, Rtree* rtree, vector<InvDIX*> invertedindex, vector<Trajectory*> TraList, vector<Diversification*> resultpairs);

	vector<IndexablePoint*> getPointList(IndexablePoint * _elem, double _radius, Rtree * _index);

	vector<IndexablePoint*> getTrajectory(IndexablePoint * _elem, double _radius, Rtree * _index);	
	~Util();
};

