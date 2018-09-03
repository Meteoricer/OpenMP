#include "Rule_Anchoring_Deanchoring.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
//#include "Polymer_Unit.h"
#include "Grid_Unit.h"
#include <vector>
#include <cmath>
#include <random>
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
//extern vector<vector<Grid_Unit>> Polymer_Grid;
//extern vector<vector<Grid_Unit>> Polymer_Cluster_Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

Polymer_Cluster::Polymer_Cluster()
{
	this->availabe_polymerize_direction.resize(4);
	/*Polymer_Unit polymer_unit;
	this->polymer_sequence.push_back(polymer_unit);*/
	this->direction = 0;
	this->availabe_polymerize_direction[direction_up] = false;
	this->availabe_polymerize_direction[direction_down] = false;
	this->availabe_polymerize_direction[direction_left] = false;
	this->availabe_polymerize_direction[direction_right] = false;
	this->cluster_bundle_pointer = rule_structure.polymer_bundle_list.end();
	//Now i assume all is verticle
	
}

//Polymer_Cluster::~Polymer_Cluster()
//{
//	rule_parameters.DD_depolymerize_end_sum = rule_parameters.DD_depolymerize_end_sum - this->DD_depolymerize_end_sum;
//	rule_parameters.DD_fragmentation_sum -= this->DD_fragmentation_sum;
//	rule_parameters.TT_depolymerize_end_sum -= this->TT_depolymerize_end_sum;
//	rule_parameters.TTT_sum_total -= this->TTT_sum;
//	rule_parameters.TTD_sum_total -= this->TTD_sum;
//	
//}