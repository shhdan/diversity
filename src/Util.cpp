#include"Util.h"

Util::Util()
{
}

void Util::UpperBoundTopkDiversifiedSearch(vector<vector<IndexablePoint*>> odpairs, double _oradius, 
	double _dradius, int k, long st, long et, Rtree* rtree, vector<InvDIX*> invertedindex, 
	vector<Trajectory*> TraList, vector<Diversification*> resultpairs) {
	if (k > odpairs.size()) {
		printf("Warning::Util::UpperBoundTopkDiversifiedSearch()::The number of od pairs is less than k!\n");
		return;
	}

	vector<Diversification*> diversitylist;

	for (int i = 0; i < k; i++) {
		Diversification* diversity = new Diversification();
		diversity->set_OD(odpairs[i][0], _oradius, odpairs[i][1], _dradius, st, et);
		diversity->set_intersectionTrajectory(rtree, invertedindex);
		diversity->set_diversity_with_distance(TraList);

		//debug
		printf("****************the diversity = %f \n", diversity->get_diversity());
		//debug

		resultpairs.push_back(diversity);

	//	diversity->~Diversification();
	}

	for (int i = k; i < odpairs.size(); i++) {
		Diversification* diversity = new Diversification();
		diversity->set_OD(odpairs[i][0], _oradius, odpairs[i][1], _dradius, st, et);
		diversity->set_intersectionTrajectory(rtree, invertedindex);
		diversity->set_diversity_upper_bound(TraList);

		double minidiversity = 1000000000;
		int miniidx = 0;
		for (int j = 0; j < resultpairs.size(); j++) {
			if (resultpairs[j]->get_diversity() < minidiversity) {
				minidiversity = resultpairs[j]->get_diversity();
				miniidx = j;
			}
		}

		if (diversity->get_diversity_upper_bound() < minidiversity)
			continue;

		else {
			diversity->set_diversity_with_distance(TraList);
		//debug
		printf("****************the diversity = %f \n", diversity->get_diversity());
		//debug
			if (diversity->get_diversity() > minidiversity) {
				resultpairs.at(miniidx) = diversity;
			}
		}

	//	diversity->~Diversification();
	}
}

void Util::BruteForceTopkDiversifiedSearch(vector<vector<IndexablePoint*>> odpairs, double _oradius, 
	double _dradius, int k, long st, long et, Rtree* rtree, vector<InvDIX*> invertedindex, 
	vector<Trajectory*> TraList, vector<Diversification*> resultpairs) {

	if (k > odpairs.size()) {
		printf("Warning::Util::BruteForceTopkDiversifiedSearch()::The number of od pairs is less than k!\n");
		return;
	}

	vector<Diversification*> diversitylist;

	for (int i = 0; i < odpairs.size(); i++) {
		Diversification* diversity = new Diversification();
		diversity->set_OD(odpairs[i][0], _oradius, odpairs[i][1], _dradius, st, et);
		diversity->set_intersectionTrajectory(rtree, invertedindex);
		diversity->set_diversity_with_distance(TraList);
		//debug
		printf("****************the diversity = %f \n", diversity->get_diversity());
		//debug
		diversitylist.push_back(diversity);

	//	diversity->~Diversification();
	}

	for (int i = 0; i < k; i++) {
		resultpairs.push_back(diversitylist[i]);
	}

	for (int i = k; i < diversitylist.size(); i++) {
		for (int j = 0; j < resultpairs.size(); j++) {
			if (diversitylist[i]->get_diversity() > resultpairs[j]->get_diversity()) {
				resultpairs.at(j) = diversitylist[i];
				break;
			}
		}
	}
}


/*
 * For 3D Rtree, to get point list by range query.
 * The result list contains only spatial-temporal points.
 * */
vector<IndexablePoint*> Util::getPointList(IndexablePoint* _elem, double _radius, Rtree* _index) {
	if (!_index->isValid()) {
		printf("Warning in Diversification::findPointList: The Rtree is NULL!\n");
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

/* get trajectory with 3D rtree where each point contains trajecotry id information*/
vector<IndexablePoint*> Util::getTrajectory(IndexablePoint* _elem, double _radius, Rtree* _index) {
	/*Perform range query to get spatial points*/
	vector<IndexablePoint*> reslut_invertedlist;
	vector<IndexablePoint*> pointlist = getPointList(_elem, _radius, _index);

	unordered_map<int, IndexablePoint*> trajectorylist;

	for (int i = 0; i < pointlist.size(); i++) {
		int tra_id = pointlist[i]->tra_id;
		unordered_map<int, IndexablePoint*>::iterator it;
		it = trajectorylist.find(tra_id);

		if (it != trajectorylist.end()) {
			long time_this = pointlist[i]->coords[2];
			long time_map = trajectorylist[tra_id]->coords[2];
			if (time_this > time_map) {
				trajectorylist[tra_id] = pointlist[i];
			}
		}

		else
		{
			trajectorylist[tra_id] = pointlist[i];
		}

	}

	for (auto elem : trajectorylist) {
		reslut_invertedlist.push_back(elem.second);
	}

	return reslut_invertedlist;

}

Util::~Util()
{
}
