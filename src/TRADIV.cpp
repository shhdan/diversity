// TRADIV.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Rtree.h"
#include "Trajectory.h"
#include "InvDIX.h"
#include "IndexablePoint.h"
#include "Diversification.h"
#include "Util.h"

#include <math.h>
#include <cmath> 
#define earthRadiusKm 6371.0

using namespace std;
vector<string> split(string s, char delim);
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d);

int main()
{	

	/**********************************Build The Rtree From Input Pointlist************************************/
//	printf("******************Build The Rtree From Input Pointlist*************** \n");
//
//	Rtree rtree = Rtree(2, 12, 12);
//
//	ifstream pointlistIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/pointlist");
//
//	if (!pointlistIn) {
//		printf("Warning: cannot open pointlist for building Rtree! \n");
//		return 1;
//	}
//
//	string pointline;
//	
//	vector<IndexablePoint> pointlist;
//	
//
//	while (pointlistIn && getline(pointlistIn, pointline)) {
//		IndexablePoint point = IndexablePoint();
//		
//		vector<string> rowdata = split(pointline, ',');
//
//		point.id = stoi(rowdata[0], nullptr, 0);
//
//		int* coordinates = new int[2];
//		coordinates[0] = stoi(rowdata[1], nullptr, 0);
//		coordinates[1] = stoi(rowdata[2], nullptr, 0);
//		point.coords = coordinates;
//
//		pointlist.push_back(point);
//
//		point.~IndexablePoint();
//		rowdata.clear();
//	}
//
//	
//	
//	IndexablePoint** ptList = new IndexablePoint*[pointlist.size()];
//	for (int i = 0; i < pointlist.size(); i++) {
//		ptList[i] = &pointlist[i];
//	}
//	
//	RtreeNode* root = rtree.constructRtree_STR(ptList, pointlist.size(), 2, 12, 12);
//	
//
//	/**************************************Build The Inverted Index With InFile*******************************************/
//	printf("******************Build The Inverted Index With InFile*************** \n");
//
//	vector<InvDIX*> invertedindex;
//	ifstream invertedlistIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/invertedindex-oneday");
//
//	if (!invertedlistIn) {
//		printf("Warning: cannot open invertedlist for building InvertedIndex! \n");
//		return 1;
//	}
//
//	string invertedline;
//	int check = 0;
//	while (invertedlistIn && getline(invertedlistIn, invertedline)) {
//		
//		//if (check > 100)
//		//	break;
//		//check++;
//
//		InvDIX* invertedlist = new InvDIX();
//
//		vector<string> rowdata = split(invertedline, ',');
//
//		invertedlist->set_invertedlistId(stoi(rowdata[0], nullptr, 0));
//
//		if (rowdata.size() == 1) {
//			invertedindex.push_back(invertedlist);
//			rowdata.clear();
//			//invertedlist.~InvDIX();
//			continue;
//		}
//
//		vector<IndexablePoint*> trajectorylist;
//
//		for (int i = 1; i < rowdata.size() - 2; i = i + 3) {
//
//			/* The attributes of elem is as follow:
//				elem.id = trajectory.Id
//				coords[0] = position
//				coords[1] = time
//			*/
//			IndexablePoint* elem = new IndexablePoint();
//
//			elem->id = stoi(rowdata[i], nullptr, 0);
//
//			int* coordinates = new int[2];
//			coordinates[0] = stoi(rowdata[i + 1], nullptr, 0);
//			coordinates[1] = stoi(rowdata[i + 2], nullptr, 0);
//			elem->coords = coordinates;
//
//			trajectorylist.push_back(elem);
//
//		}
//
//		invertedlist->set_trajectoryArray(trajectorylist);
//		invertedindex.push_back(invertedlist);
//
//		rowdata.clear();
//		trajectorylist.clear();
//		//invertedlist.~InvDIX();
//	
//	}
//
	//printf("check inverted index = %d, %d", invertedindex[90]->get_invertedlistId(), invertedindex[90]->get_trajectoryArray()[0]->coords[0]);

	/*******************************Build The Trajectory List*********************************/
	printf("******************Build The Trajectory List*************** \n");

	vector<Trajectory*> TraList;
	ifstream trajectorylistIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/trajectory-node-point");

	if (!trajectorylistIn) {
		printf("Warning: cannot open trajectorylist for building Trajectory! \n");
		return 1;
	}

	string trajectoryline;


	long pointnum = 0;
	
	while (trajectorylistIn && getline(trajectorylistIn, trajectoryline)) {

		Trajectory* trajectory = new Trajectory();

		vector<string> rowdata = split(trajectoryline, ',');

		trajectory->set_trajectoryId(stoi(rowdata[0], nullptr, 0));

		vector<IndexablePoint*> pointlist;

		for (int i = 1; i < rowdata.size() - 3; i = i + 4) {

			/* The attributes of elem is as follow:
				elem.id = point.Id
				coords[0] = longitude
				coords[1] = latitude
				coords[2] = time
			*/
			IndexablePoint* elem = new IndexablePoint();

			elem->id = stoi(rowdata[i], nullptr, 0);

			int* coordinates = new int[3];
			coordinates[0] = stoi(rowdata[i + 1], nullptr, 0);
			coordinates[1] = stoi(rowdata[i + 2], nullptr, 0);
			coordinates[2] = stoi(rowdata[i + 3], nullptr, 0);
			elem->coords = coordinates;

			pointlist.push_back(elem);
			
			pointnum++;

		}

		trajectory->set_pointArray(pointlist);
		TraList.push_back(trajectory);

		rowdata.clear();
		pointlist.clear();
      	//trajectory.~Trajectory();
      
	}

	printf("the number of point in this trajectory set = %d \n", pointnum);
	/*******************************Build The 3D Rtree*********************************/
//	printf("******************Build The 3D Rtree Index*************** \n");
//
//	Rtree rtree_3d = Rtree(3, 12, 12);
//	ifstream trajectorypointIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/trajectory-node-point-oneday");
//
//	if (!trajectorypointIn) {
//		printf("Warning: cannot open trajectorylist for building 3D Rtree! \n");
//		return 1;
//	}
//
//	string trajectorypointline;
//
//	vector<IndexablePoint*> pointlist;
//	while (trajectorypointIn && getline(trajectorypointIn, trajectorypointline)) {
//
//		vector<string> rowdata = split(trajectorypointline, ',');
//
//		int tra_id = stoi(rowdata[0], nullptr, 0);
//
//		for (int i = 1; i < rowdata.size() - 3; i = i + 4) {
//
//			/* The attributes of elem is as follow:
//			   elem.id = point.Id
//			   coords[0] = longitude
//			   coords[1] = latitude
//			   coords[2] = time
//			   */
//			IndexablePoint* elem = new IndexablePoint();
//
//			elem->id = stoi(rowdata[i], nullptr, 0);
//			elem->tra_id = tra_id;
//
//			int* coordinates = new int[3];
//			coordinates[0] = stoi(rowdata[i + 1], nullptr, 0);
//			coordinates[1] = stoi(rowdata[i + 2], nullptr, 0);
//			coordinates[2] = stoi(rowdata[i + 3], nullptr, 0);
//			elem->coords = coordinates;
//
//			pointlist.push_back(elem);
//			elem->~IndexablePoint();
//		}
//
//		rowdata.clear();
//	}
//
//			clock_t start, end;
//			start = clock();
//	IndexablePoint** ptList_3d = new IndexablePoint*[pointlist.size()];
//	for (int i = 0; i < pointlist.size(); i++) {
//		ptList_3d[i] = pointlist[i];
//	}
//
//	RtreeNode* root_3d = rtree_3d.constructRtree_STR(ptList_3d, pointlist.size(), 3, 12, 12);
//
//			end = clock();
//			float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//			printf("*************the running time = %f\n", time);
//
	/************************************Ready For Query******************************************/
//	printf("******************Ready For Query*************** \n");
//	
//	ifstream odorderIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/od.txt");
//
//	if (!odorderIn) {
//		printf("Warning: cannot open od order for query processing! \n");
//		return 1;
//	}
//
//	string odorderline;
//
//	vector<IndexablePoint*> odorder;
//	while (odorderIn && getline(odorderIn, odorderline)) {
//
//		vector<string> rowdata = split(odorderline, '\t');
//
//		IndexablePoint* elem = new IndexablePoint();
//		int* coordinates = new int[2];
//		coordinates[0] = stoi(rowdata[0], nullptr, 0);
//		coordinates[1] = stoi(rowdata[1], nullptr, 0);
//
//		elem->coords = coordinates;
//
//		odorder.push_back(elem);
//
//		rowdata.clear();
//	}
//
//	ofstream out;
//	out.open("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/output-1d-bf", ios_base::out);
//
//	ifstream odpairsIn("/media/dragon_data/uqdhe/BeijingFiveDays/diversity/odpairs");
//
//	if (!odpairsIn) {
//		printf("Warning: cannot open odpairs for query processing! \n");
//		return 1;
//	}
//
//	string odpairsline;
//
//	vector<IndexablePoint*> odpairs;
//	while (odpairsIn && getline(odpairsIn, odpairsline)) {
//
//		vector<string> rowdata = split(odpairsline, ',');
//
//		IndexablePoint* elem = new IndexablePoint();
//		int* coordinates = new int[2];
//		coordinates[0] = stoi(rowdata[0], nullptr, 0);
//		coordinates[1] = stoi(rowdata[1], nullptr, 0);
//		//coordinates[2] = 1427875200;
//
//		elem->coords = coordinates;
//
//		odpairs.push_back(elem);
//
//		rowdata.clear();
//	}
//
//	vector<vector<IndexablePoint*>> allpairs;
//
//	for(int i = 0; i < odorder.size(); i++){
//		vector<IndexablePoint*> onepair;
//		onepair.push_back(odpairs[odorder[i]->coords[0]]);
//		onepair.push_back(odpairs[odorder[i]->coords[1]]);
//		allpairs.push_back(onepair);
//
//		onepair.clear();
//	}
//
//	printf("the number of pairs = %d\n", allpairs.size());
//	int* k = new int[4];
//	k[0] = 1;
//	k[1] = 4;
//	k[2] = 16;
//	k[3] = 64;
//
//	for(int j = 0; j < 4; j++){
//
//		Util* util = new Util();
//		vector<Diversification*> resultpairs;
//
//		clock_t start, end;
//		start = clock();
//		util->BruteForceTopkDiversifiedSearch(allpairs, 100, 100, k[j], 7, 10, &rtree, invertedindex, TraList, resultpairs);
//		//util->UpperBoundTopkDiversifiedSearch(allpairs, 100, 100, k[j], 7, 10, &rtree, invertedindex, TraList, resultpairs);
//
//		end = clock();
//		float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//		printf("*************the running time = %f\n", time);
//
//		if(out.is_open()){
//			out << k[j] << ",";
//			//out << resultpairs[0]->get_diversity() << ",";
//			out << time << "\n";
//		}
//		//diversity->~Diversification();
//	}
//	out.close();

	//	int* temp1 = new int[2];
	//	int* temp2 = new int[2];
//	temp1[0] = 1163025;
//	temp1[1] = 399824;
//	temp2[0] = 1163588;
//	temp2[1] = 399430;
//	IndexablePoint opoint = IndexablePoint(1, temp1);
//	IndexablePoint dpoint = IndexablePoint(2, temp2);
//	for(int i = 0; i < 22; i=i+2){
//
//		Diversification* diversity = new Diversification();
//		printf("******************set_odpair*************** \n");
//		diversity->set_OD(&opoint, 200, &dpoint, 200, i, i+2);
//
//		diversity->set_intersectionTrajectory(&rtree, invertedindex);
//		clock_t start, end;
//		start = clock();
//		diversity->set_diversity_with_distance(TraList);
//		end = clock();
//		float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//		printf("*************the running time = %f\n", time);
//		printf("******************set_diversity_upper_bound*************** \n");
//		diversity->set_diversity_upper_bound(TraList);
//		printf("******************get_diversity*************** \n");
//		printf("check the diversity = %lf \n", diversity->get_diversity());
//		printf("******************get_diversity_upper_bound*************** \n");
//		printf("check the diversity upper bound = %lf \n", diversity->get_diversity_upper_bound());
//
//		if(out.is_open()){
//			out << i << ",";
//			out << diversity->get_intersection() << ",";
//			out << diversity->get_diversity() << ",";
//			out << diversity->get_diversity_upper_bound() << "\n";
//			//	out << distance << "\n";
//		}
//	}

//	for(int i = 0; i < odpairs.size(); i++){
//		for(int j = 0; j < odpairs.size(); j++){
//			double distance = 0;
//			double x1 = ((double)odpairs[i]->coords[1])/10000.00;
//			double y1 = ((double)odpairs[i]->coords[0])/10000.00;
//			double x2 = ((double)odpairs[j]->coords[1])/10000.00;
//			double y2 = ((double)odpairs[j]->coords[0])/10000.00;
//			distance = distanceEarth(x1,y1,x2,y2);
//			distance = odpairs[i]->computeDist(odpairs[j]->coords, 2);
//			Diversification* diversity = new Diversification();
//			printf("******************set_odpair*************** \n");
//			diversity->set_OD(odpairs[i], 100, odpairs[j], 100, 7, 10);
//
//			diversity->set_intersectionTrajectory(&rtree, invertedindex);
//			clock_t start, end;
//			start = clock();
//			diversity->set_diversity_with_distance(TraList);
//			end = clock();
//			float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//			printf("*************the running time = %f\n", time);
//			printf("******************set_diversity_upper_bound*************** \n");
//			diversity->set_diversity_upper_bound(TraList);
//			printf("******************get_diversity*************** \n");
//			printf("check the diversity = %lf \n", diversity->get_diversity());
//			printf("******************get_diversity_upper_bound*************** \n");
//			printf("check the diversity upper bound = %lf \n", diversity->get_diversity_upper_bound());
//			
//			if(out.is_open()){
//				out << i << "," << j << ",";
//				out << diversity->get_intersection() << ",";
//				out << diversity->get_diversity() << ",";
//				out << diversity->get_diversity_upper_bound() << ",";
//				out << distance << "\n";
//			}
//			//diversity->~Diversification();
//		}
//	}
//	out.close();
//
//	for(int i = 0; i < odpairs.size(); i++){
////		Util* util = new Util();
//		Diversification* diversity = new Diversification();
//		printf("******************set_odpair*************** \n");
//		diversity->set_OD(odpairs[i], 100, odpairs[i], 100, 1427875200, 1427875300);
//		//	Diversification* diversity = new Diversification();
//		//	int* temp1 = new int[2];
//		//	int* temp2 = new int[2];
//		//	temp1[0] = 1164650;
//		//	temp1[1] = 402060;
//		//	temp2[0] = 1163790;
//		//	temp2[1] = 401750;
//		//	IndexablePoint opoint = IndexablePoint(1, temp1);
//		//	IndexablePoint dpoint = IndexablePoint(2, temp2);
//		//	printf("******************set_OD*************** \n");
//		//	diversity->set_OD(&opoint, 1000, &dpoint, 1000, 1427875200, 1427876200);
//		printf("******************get_intersectionTrajectory*************** \n");
//		clock_t start, end;
//		start = clock();
//		
//	//	vector<IndexablePoint*> ODTrajectory = util->getTrajectory(odpairs[i], 100, &rtree_3d);
//		vector<IndexablePoint*> ODTrajectory = diversity->getTrajectory(odpairs[i], 100, &rtree, invertedindex, 1427875200, 1427875300);
//		end = clock();
//		float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//		printf("*************the running time = %f\n", time);
//		printf("the number of trajectory in OD set = %d\n", ODTrajectory.size());
//
//		if(out.is_open()){
//			out << ODTrajectory.size() << ",";
//			out << time << "\n";
//		}
//		
//	}
//	out.close();
//
//
//
//
//	printf("******************set_intersectionTrajectory*************** \n");
//	diversity.set_intersectionTrajectory(&rtree, invertedindex);
//	printf("******************set_diversity_with_distance*************** \n");
//	clock_t start, end;
//	start = clock();
//	diversity.set_diversity_with_distance(TraList);
//	end = clock();
//	float time = ((float)end - (float) start)/CLOCKS_PER_SEC;
//	printf("*************the running time = %f\n", time);
//	printf("******************set_diversity_upper_bound*************** \n");
//	diversity.set_diversity_upper_bound(TraList);
//	printf("******************get_diversity*************** \n");
//	printf("check the diversity = %lf \n", diversity.get_diversity());
//	printf("******************get_diversity_upper_bound*************** \n");
//	printf("check the diversity upper bound = %lf \n", diversity.get_diversity_upper_bound());
//	printf("*******this is the end of this program********* \n");
	return 0;
}



vector<string> split(string s, char delim) {
	stringstream ss(s);
	string item;
	vector<string> tokens;
	while (getline(ss, item, delim)) {
		tokens.push_back(item);
	}
	return tokens;
}


// This function converts decimal degrees to radians
double deg2rad(double deg) {
	return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
	return (rad * 180 / M_PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
	double lat1r, lon1r, lat2r, lon2r, u, v;
	lat1r = deg2rad(lat1d);
	lon1r = deg2rad(lon1d);
	lat2r = deg2rad(lat2d);
	lon2r = deg2rad(lon2d);
	u = sin((lat2r - lat1r)/2);
	v = sin((lon2r - lon1r)/2);
	return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}
