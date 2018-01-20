//#include "stdafx.h"
#include "stdio.h"
#include "Trajectory.h"

const long long MYINFTY = 1000000000;

Trajectory::Trajectory() {
	this->trajectoryId = -1;
	vector<IndexablePoint*> pointlist = vector<IndexablePoint*>();
	this->pointArray = pointlist;
}

void Trajectory::addPointToArray(IndexablePoint* point) {
	this->pointArray.push_back(point);
}

void Trajectory::set_trajectoryId(int id) {
	this->trajectoryId = id;
}

int Trajectory::get_trajectoryId() {
	return this->trajectoryId;
}


vector<IndexablePoint*> Trajectory::get_pointArray() {
	return this->pointArray;
}

void Trajectory::set_pointArray(vector<IndexablePoint*> pointlist) {
	this->pointArray = pointlist;
}

double Trajectory::EuclideanDistance(Trajectory* Tra, int _dim) {
	double euclideandist = 0.0;
	
	if (this->trajectoryId == Tra->get_trajectoryId()) {
		euclideandist = 0.0;
		return euclideandist;
	}

	for (int i = 0; i < this->pointArray.size(); i++) {
		double mindist = MYINFTY;
		long node1_id = this->pointArray[i]->id;
		
		//debug
		//printf("the number of points in the other trajectory = %d \n", Tra->get_pointArray().size());
		//printf("i = %d \n", i);
		//printf("trajectory id = %d\n", Tra->trajectoryId);
		//debug

		for (int j = 0; j < Tra->get_pointArray().size(); j++) {
			double dist = 0.0;
			
			//debug
		//	printf("****\n");
		//	printf("node2_id = %d \n", Tra->get_pointArray()[j]->coords[0]);
		//	printf("node2_id = %d \n", Tra->get_pointArray()[j]->coords[1]);
			//debug
			long node2_id = Tra->get_pointArray()[j]->id;

			////debug
		//	printf("node2_id = %d \n", Tra->get_pointArray()[j]->id);
			//debug
			if (node1_id == node2_id) {
				mindist = 0.0;
				break;
			}
			else {
				dist = this->pointArray[i]->computeDist(Tra->get_pointArray()[j]->coords, _dim);
			}
			if (dist < mindist)
				mindist = dist;
		}

		euclideandist = euclideandist + mindist;
	}
	euclideandist = euclideandist / this->pointArray.size();

	return euclideandist;
}

/*
 * calculate the Eucliean Distance between trajectory and the od points
 * */
double Trajectory::EuclideanDistanceODPoint(IndexablePoint* opoint, IndexablePoint* dpoint, int _dim) {
	double euclideandist = 0.0;

	for (int i = 0; i < this->pointArray.size(); i++) {
		double mindist = MYINFTY;
		long node_id = this->pointArray[i]->id;

		if (node_id == opoint->id || node_id == dpoint->id) {
			mindist = 0;
		}

		else {
			double odistance = this->pointArray[i]->computeDist(opoint->coords, _dim);
			if (odistance < mindist)
				mindist = odistance;

			double ddistance = this->pointArray[i]->computeDist(dpoint->coords, _dim);
			if (ddistance < mindist)
				mindist = ddistance;
		}

		euclideandist = euclideandist + mindist;
	}
	euclideandist = euclideandist / this->pointArray.size();

	return euclideandist;
}

Trajectory::~Trajectory()
{
	pointArray.clear();
}
