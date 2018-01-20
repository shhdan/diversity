#pragma once
#include<vector>
#include<ctime>
#include "Point.h"
#include "Trajectory.h"
#include "InvDIX.h"
#include "Rtree.h"
using namespace std;

class Diversification
{
private:
	IndexablePoint* original;  // the origin region center
	double oradius;
	IndexablePoint* destination; // the destination region center
	double dradius;
	long starttime;  // the start time of the time window (optional)
	long endtime;  // the end time of the time window (optional)
	double diversity;
	//int nTop;
	double upperbound;
	//double lowerbound;
	//long ntotalPoint;
	vector<IndexablePoint*> intersectionTrajectory; // the trajectories traversing from the origin to the destination
	int nintersection; // the number of intersecting trajectories, initialze when getting the intersectionTrajectory
public:
	Diversification();
	void set_OD(IndexablePoint* _original, double _oradius, IndexablePoint* _destination, double _dradius, long _stime, long _etime);

	vector<IndexablePoint*> getPointList(IndexablePoint* _elem, double _radius, Rtree* _index);
	vector<IndexablePoint*> getInvertedList(vector<InvDIX*> _invindex, int _pointId);
	vector<IndexablePoint*> mergeTwoInvertedList(vector<IndexablePoint*> _invertedlist1, vector<IndexablePoint*> _invertedlist2);
	vector<IndexablePoint*> mergeInvertedList(vector<vector<IndexablePoint*>> _invertedlists);
	vector<IndexablePoint*> getTrajectory(IndexablePoint* _elem, double _radius, Rtree* _index, vector<InvDIX*> _invindex, long _timestart, long _timeend);
	void set_intersectionTrajectory(Rtree* _index, vector<InvDIX*> _invindex);
	vector<IndexablePoint*> get_intersectionTrajectoryArray();
	void set_diversity_with_distance(vector<Trajectory*> tralist_all);
	void set_diversity_upper_bound(vector<Trajectory*> tralist_all);
	int get_intersection();
	double get_diversity();
	double get_diversity_upper_bound();
	//void set_nTop(int numk);
	~Diversification();
};


