#pragma once
//thi is the test2 file, can not include test1 file, but has to define the class here.
//#include "Anchor_Unit.h"
#include <vector>
#include <list>
using namespace std;
class Anchor_Unit;
class Grid_Unit;
class Polymer_Unit;
class Polymer_Bundle;
class Polymer_Cluster
{
public:
	Polymer_Cluster();
	/*~Polymer_Cluster();*/
	//Anchor_Unit *anchor_pointer;
	list<list<Anchor_Unit>::iterator> anchor_pointer_list;//it can be multiple anchor
	int direction;//horizontal or vertical
	int polarize;//direction first or direction second
	list<Polymer_Unit> polymer_sequence;//left to right, down to up
	//Grid
	int TTT_sum;
	//list<Polymer_Unit>::iterator polymer_unit_pointer;
	list<Polymer_Bundle>::iterator cluster_bundle_pointer;
	//list<Polymer_Unit> bundling_list;
	//list<Polymer_Unit> debundling_list;
	/*
	there are two options to record the hydrolysis detail:
	vector: slow. Maybe scan the vector after each seperation?
	list: how the mantain the list. Hard to get the list information,
	maybe rebuild the list after each seperation? And how detail the list should be?

	Use the list solution for now

	*/
	
	int TTD_sum=0;//TTD and DTT are the same
	int DTD_sum=0;
	//this is for fragmentation
	int TT_sum=0;
	int TD_sum=0;
	int DD_sum=0;
	//use the same direction trick as in polymerize_end_sum
	int TT_depolymerize_plus_end_sum=0;
	int TD_depolymerize_plus_end_sum=0;//TD and DT are the same here
	int DD_depolymerize_plus_end_sum=0;
	int TT_depolymerize_minus_end_sum = 0;
	int TD_depolymerize_minus_end_sum = 0;//TD and DT are the same here
	int DD_depolymerize_minus_end_sum = 0;
	//int DD_fragmentation_sum = 0;
	double Consider_mark = 0;
	int ring_mark = 0;


	int TT_polymerize_plus_end_sum=0;//directions
	int TT_polymerize_minus_end_sum = 0;//directions
	int TD_polymerize_plus_end_sum = 0;//directions
	int TD_polymerize_minus_end_sum = 0;//directions
	vector<int> availabe_polymerize_direction;//I can even use it to detect polymerize direction
};