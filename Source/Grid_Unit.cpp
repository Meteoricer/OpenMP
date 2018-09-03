#include "Grid_Unit.h"
#include <iostream>
//#include "parameters.h"
using namespace std;


extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

Grid_Unit::Grid_Unit()
{
	this->grid_anchor_pointer = rule_structure.anchor_list.end();
	this->grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
	//this->grid_polymer_unit_pointer= rule_structure.polymer_unit_list_end.begin();
	//this->Is_Anchored = no;
}
void Grid_Unit::output()
{
	//cout << Grid.size();
}