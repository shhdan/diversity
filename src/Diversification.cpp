//#include "stdafx.h"
#include "Diversification.h"
#include "stdio.h"

/*
	To calculate the diversity of trajectories with respect 
	to two OD regions.
*/

const long DATESECOND = 86400;
const long HOURSECOND = 3600;

Diversification::Diversification()
{
	this->original = nullptr;
	this->oradius = 0.0;
	this->destination = nullptr;
	this->dradius = 0.0;
	this->diversity = -1;
	this->starttime = 0;
	this->endtime = 0;
	//nTop = -1;
	this->upperbound = -1;
	//lowerbound = 10000000.00;
	//ntotalPoint = -1;
	//intersectionTrajectory = nullptr;
	this->intersectionTrajectory = vector<IndexablePoint*>();
	this->nintersection = -1;
}

//void Diversification::set_nTop(int numk) {
//	this->nTop = numk;
//}

void Diversification::set_OD(IndexablePoint* _original, double _oradius, IndexablePoint* _destination, double _dradius, long _stime, long _etime) {
	this->original = _original;
	this->oradius = _oradius;
	this->destination = _destination;
	this->dradius = _dradius;
	this->starttime = _stime;
	this->endtime = _etime;
}

/*
	For 2D Rtree, to get point list by range query.
	The result list contains only spatial points, each of the 
	point contains a inverted list referred to trajectories.
*/
vector<IndexablePoint*> Diversification::getPointList(IndexablePoint* _elem, double _radius, Rtree* _index) {
	if (! _index->isValid()) {
		printf("Warning in Diversification::findPointList: The Rtree is nullptr!\n");
	}
	if (_radius == 0.0) {
		printf("Warning in Diversification::findPointList: The _radius is zero!\n");
	}
	//clean the target list
	vector<IndexablePoint*> toPointlist;

	//rangeQuery to find candidate pointlist
	_index->rangeQuery(_elem, _radius, toPointlist);

	return toPointlist;
}

/*
	Get the inverted list with Hashing.
	The hashing map the point id to the index of corresponding inverted list.
	FIXME: Here point id = index, will be fixed later.
*/
vector<IndexablePoint*> Diversification::getInvertedList(vector<InvDIX*> _invindex, int _pointId) {
	int idx = _pointId;
	return _invindex[idx]->get_trajectoryArray();
}

/*
	Merge two inverted lists by trajectory id. The inverted lists are 
	sorted by trajectory id. When doing the merging, for two elements
	with the same trajectory id, we will keep the element whose time
	is larger, namely the subsequent point.
*/
vector<IndexablePoint*> Diversification::mergeTwoInvertedList(vector<IndexablePoint*> _invertedlist1, vector<IndexablePoint*> _invertedlist2) {
	
	vector<IndexablePoint*> result_invertedlist;
	if (_invertedlist1.size() == 0) {
		result_invertedlist = _invertedlist2;
	}
	else if (_invertedlist2.size() == 0) {
		result_invertedlist = _invertedlist1;
	}
	else {
		int idx1 = 0, idx2 = 0;
		while (idx1 < _invertedlist1.size() && idx2 < _invertedlist2.size()) {
			if (_invertedlist1[idx1]->id < _invertedlist2[idx2]->id) {
				result_invertedlist.push_back(_invertedlist1[idx1]);
				idx1++;
			}
			else if (_invertedlist1[idx1]->id > _invertedlist2[idx2]->id) {
				result_invertedlist.push_back(_invertedlist2[idx2]);
				idx2++;
			}
			else {

				if (_invertedlist1[idx1]->coords[0] > _invertedlist2[idx2]->coords[0])
					result_invertedlist.push_back(_invertedlist1[idx1]);
				else
					result_invertedlist.push_back(_invertedlist2[idx2]);
				idx1++;
				idx2++;
			}
		}
		if (idx1 == _invertedlist1.size()) {
			for(int i = idx2; i < _invertedlist2.size(); i++)
				result_invertedlist.push_back(_invertedlist2[i]);
		}
		if (idx2 == _invertedlist2.size()) {
			for (int i = idx1; i < _invertedlist1.size(); i++)
				result_invertedlist.push_back(_invertedlist1[i]);
		}
	}

	return result_invertedlist;
}

/*
	Merge multiples inverted lists into one list. Refer to the above.
*/
vector<IndexablePoint*> Diversification::mergeInvertedList(vector<vector<IndexablePoint*> > _invertedlists) {
	vector<IndexablePoint*> result_invertedlist;
	vector<vector<IndexablePoint*>> invertedlists = _invertedlists;

	for (int i = 0; i < invertedlists.size(); i++) {
		result_invertedlist = mergeTwoInvertedList(invertedlists[i], result_invertedlist);
	}
	
	//while (invertedlists.size() > 1) {
	//	int newsize = 0;
	//	if (invertedlists.size() % 2 == 0)
	//		newsize = invertedlists.size() / 2;
	//	else
	//		newsize = invertedlists.size() / 2 + 1;
	//	for (int i = 0; i < newsize ; i++) {
	//		vector<IndexablePoint*> temp = mergeTwoInvertedList(invertedlists[i], invertedlists[invertedlists.size() - i - 1]);
	//		invertedlists[invertedlists.size() - i - 1].clear();
	//		invertedlists[i].clear();
	//		invertedlists[i] = temp;
	//		temp.clear();
	//	}
	//}
	//result_invertedlist = invertedlists[0];

//	invertedlists.clear();
	return result_invertedlist;
}


/*
 * get trajectory with 2D rtree and inverted index
 * */
vector<IndexablePoint*> Diversification::getTrajectory(IndexablePoint* _elem, double _radius, Rtree* _index, vector<InvDIX*> _invindex, long _timestart, long _timeend) {
	//safety protection needed


	//clock_t start, end;	
	/*Perform range query to get spatial points*/
	vector<IndexablePoint*> reslut_invertedlist;
	//start = clock();
	vector<IndexablePoint*> pointlist = getPointList(_elem, _radius, _index);
	//end = clock();
	//float time = ((float)end - (float)start)/CLOCKS_PER_SEC;
	//printf("************* the time to perform range query = %f\n", time);
	//debug
	printf("the number of point = %d \n", pointlist.size());
	//debug

	vector<vector<IndexablePoint*> > invertedlists;
	for (int i = 0; i < pointlist.size(); i++) {
		invertedlists.push_back(getInvertedList(_invindex, pointlist[i]->id));
	}

	
	/*Check the time dimension here for each of the inverted list*/
	vector<vector<IndexablePoint*> > timefilter_invertedlists;
	for (int i = 0; i < invertedlists.size(); i++) {
		
		vector<IndexablePoint*> trajectorylist;
		for (int j = 0; j < invertedlists[i].size(); j++) {			
			//FIXME: THIS IS FOR RTREE EXPERIMENT
		//	if (invertedlists[i][j]->coords[1] < _timestart || invertedlists[i][j]->coords[1] > _timeend)
			if ((invertedlists[i][j]->coords[1] % DATESECOND)/HOURSECOND < _timestart || (invertedlists[i][j]->coords[1] % DATESECOND)/HOURSECOND > _timeend)
				continue;
			else
				trajectorylist.push_back(invertedlists[i][j]);
		}
		timefilter_invertedlists.push_back(trajectorylist);
		
		trajectorylist.clear();
	}
	reslut_invertedlist = mergeInvertedList(timefilter_invertedlists);
	pointlist.clear();
	invertedlists.clear();
	timefilter_invertedlists.clear();
	
	return reslut_invertedlist;
}

/*
	Given origin region and destination region, find the trajectories
	starting from the origin and ending at the destination by intersection.
*/
void Diversification::set_intersectionTrajectory(Rtree* _index, vector<InvDIX*> _invindex) {
	
	vector<IndexablePoint*> original_trajectorylist = getTrajectory(this->original, this->oradius, _index, _invindex, this->starttime, this->endtime);
	
	vector<IndexablePoint*> destination_trajectorylist = getTrajectory(this->destination, this->dradius, _index, _invindex, this->starttime, this->endtime);
	
	int Olength = original_trajectorylist.size();
	int Dlength = destination_trajectorylist.size();

	//debug
	printf("Olength = %d \n", Olength);
	printf("Dlength = %d \n", Dlength);
	//debug
	if (Olength == 0 || Dlength == 0) {
		this->nintersection = 0;
		return;
	}

	if (original_trajectorylist[Olength - 1]->id < destination_trajectorylist[0]->id
		|| destination_trajectorylist[Dlength - 1]->id < original_trajectorylist[0]->id) {
		this->nintersection = 0;
		return;
	}

	int Olist_index = 0;
	int Dlist_index = 0;

	long currentOID;
	long currentDID;

	while (Olist_index < Olength && Dlist_index < Dlength) {


		IndexablePoint* intersect = new IndexablePoint();

		currentOID = original_trajectorylist[Olist_index]->id;
		currentDID = destination_trajectorylist[Dlist_index]->id;

		while ((Dlist_index < Dlength - 1) && (currentOID > currentDID)) {
			Dlist_index++;
			currentDID = destination_trajectorylist[Dlist_index]->id;
		}

		if (currentOID > currentDID) {
			Dlist_index++;
			continue;
		}

		if (currentOID < currentDID) {
			Olist_index++;
			continue;
		}

		else if (currentOID == currentDID) {
			if (original_trajectorylist[Olist_index]->coords[0] < destination_trajectorylist[Dlist_index]->coords[0]) {
				int* temp = new int[2];
				intersect->id = currentOID;
				temp[0] = original_trajectorylist[Olist_index]->coords[0];
				temp[1] = destination_trajectorylist[Dlist_index]->coords[0];
				intersect->coords = temp;
				this->intersectionTrajectory.push_back(intersect);
			}
			Olist_index++;
			Dlist_index++;
		}
	}

	this->nintersection = this->intersectionTrajectory.size();
	printf("intersection= %d \n", this->nintersection);

}

vector<IndexablePoint*> Diversification::get_intersectionTrajectoryArray() {
	return this->intersectionTrajectory;
}

/*
	Given a set of trajectories with the same origin and destination, 
	calculate the diversity based on pairwise distances.
	
*/
void Diversification::set_diversity_with_distance(vector<Trajectory*> tralist_all) {

	double distance = 0.0;

	if (this->nintersection == 0) {
		this->diversity = -1;
		return;
	}

	else if (this->nintersection == 1) {
		this->diversity = 0;
		return;
	}

	else {

		for (int i = 0; i < this->nintersection - 1; i++) {
		
//		if(i % 10 == 0)
//			printf("i = %d \n", i);

			int traid_current, start_current, end_current;
			traid_current = (int)this->intersectionTrajectory[i]->id;
			start_current = (int) this->intersectionTrajectory[i]->coords[0];
			end_current = (int) this->intersectionTrajectory[i]->coords[1];

			Trajectory* tra_current =new Trajectory();
			tra_current->set_trajectoryId(traid_current);

			vector<IndexablePoint*> tra_current_pointlist;
			IndexablePoint* point_current = new IndexablePoint();

			for (int k = start_current; k < end_current; k++) {
				point_current = tralist_all[traid_current]->get_pointArray()[k];
				tra_current_pointlist.push_back(point_current);
			}

			tra_current->set_pointArray(tra_current_pointlist);

			//point_current->releaseSpace();
			tra_current_pointlist.clear();

			for (int j = i + 1; j < this->nintersection; j++) {

				int traid_subsequent, start_subsequent, end_subsequent;
				traid_subsequent = (int)this->intersectionTrajectory[j]->id;
				start_subsequent = (int) this->intersectionTrajectory[j]->coords[0];
				end_subsequent = (int) this->intersectionTrajectory[j]->coords[1];

				Trajectory* tra_subsequent = new Trajectory();
				tra_subsequent->set_trajectoryId(traid_subsequent);

				vector<IndexablePoint*> tra_subsequent_pointlist;
				IndexablePoint* point_subsequent = new IndexablePoint();

				for (int k = start_subsequent; k < end_subsequent; k++) {
					point_subsequent = tralist_all[traid_subsequent]->get_pointArray()[k];
					tra_subsequent_pointlist.push_back(point_subsequent);
				}

				tra_subsequent->set_pointArray(tra_subsequent_pointlist);

				//point_subsequent->releaseSpace();
				tra_subsequent_pointlist.clear();

				distance = distance + tra_current->EuclideanDistance(tra_subsequent, 2);

				tra_subsequent->~Trajectory();
			}

			tra_current->~Trajectory();
		}
	}
	long pairwise = this->nintersection * (this->nintersection - 1) / 2;

	this->diversity = distance / pairwise;

}

/*
 * Given a set of trajectories with the same origin and destination,
 * calculate the diversity upper bound by scanning once trajectory set.
 *
 * */
void Diversification::set_diversity_upper_bound(vector<Trajectory*> tralist_all) {

	double distance = 0.0;

	if (this->nintersection == 0) {
		this->upperbound = -1;
		return;
	}

	else if (this->nintersection == 1) {
		this->upperbound = 0;
		return;
	}
	else {

		for (int i = 0; i < this->nintersection - 1; i++) {

			int traid_current, start_current, end_current;
			traid_current = (int)this->intersectionTrajectory[i]->id;
			start_current = (int) this->intersectionTrajectory[i]->coords[0];
			end_current = (int) this->intersectionTrajectory[i]->coords[1];

			Trajectory* tra_current = new Trajectory();
			tra_current->set_trajectoryId(traid_current);

			vector<IndexablePoint*> tra_current_pointlist;
			IndexablePoint* point_current = new IndexablePoint();

			for (int k = start_current; k < end_current; k++) {
				point_current = tralist_all[traid_current]->get_pointArray()[k];
				tra_current_pointlist.push_back(point_current);
			}

			tra_current->set_pointArray(tra_current_pointlist);

			tra_current_pointlist.clear();

			double current_distance = tra_current->EuclideanDistanceODPoint(this->original, this->destination, 2);

			distance = distance + current_distance;

			tra_current->~Trajectory();
		}
	}

	this->upperbound = distance / (this->nintersection - 1);

}

int Diversification::get_intersection() {
	return this->nintersection;
}

double Diversification::get_diversity() {
	return this->diversity;
}

double Diversification::get_diversity_upper_bound() {
	return this->upperbound;
}

Diversification::~Diversification() {
	delete[] this->original;
	delete[] this->destination;
	this->intersectionTrajectory.clear();
}
