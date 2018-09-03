#pragma once
//thi is the test2 file, can not include test1 file, but has to define the class here.
//#include "Anchor_Unit.h"
#include <vector>
#include <list>
using namespace std;
class Anchor_Unit;
class Grid_Unit;
class Polymer_Cluster;
class Polymer_Unit
{
public:
	//Polymer_Unit(int col, int row);
	Polymer_Unit(list<Polymer_Cluster>::iterator polymer_cluster_iter, int direction);
	Polymer_Unit();
	Polymer_Unit(int i);//this is only generate a basically empty polymer_unit with only information of hydrolysis, for check the sum of hydrolysis information
	//~Polymer_Unit();
	//Anchor_Unit *anchor_pointer;
	list<Anchor_Unit>::iterator anchor_pointer;
	int direction;
	//vector<int> polymer_sequence;
	Grid_Unit *polymer_grid_pointer;
	/*
	there are two options to record the hydrolysis detail:
	vector: slow. Maybe scan the vector after each seperation?
	list: how the mantain the list. Hard to get the list information,
	maybe rebuild the list after each seperation? And how detail the list should be?

	Use the vector solution for now 
	*/
	int Hydrolysis_mark;
	double Consider_mark;
	double lifetime;
	int action_mark;
	int attatch_mark;
	int frap_mark;
	//int Is_Anchored;
	//int pointer_true;//if there is a polymer unit actually be there
	list<Polymer_Cluster>::iterator polymer_cluster_pointer;
	//vector<list<Polymer_Unit>::iterator> bundled_neighbor_polymer_unit_pointer;//this should be initialize to polymer_unit.end()
	//should i use this or just the pair relation? 
	
	vector<list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator> anealing_polymer_unit_pair_list_iter;
	//vector<list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator> TD_anealing_polymer_unit_pair_list_iter;
	//vector<list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator> TT_anealing_polymer_unit_pair_list_iter;
	vector<list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator> bundling_polymer_unit_pair_list_iter;
	vector<list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator> debundling_polymer_unit_pair_list_iter;

};