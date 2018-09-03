#include "Rule_Diffusion.h"
#include "Grid_Unit.h"
#include <list>
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include <cmath>
#include <random>
#include "parameters.h"
#include "Rule_Bundling_Debundling.h"
#include "General_Functions.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Parameters rule_parameters;

Rule_Structure::Rule_Structure()
{
	Polymer_Unit polymer_unit;
	polymer_unit.Hydrolysis_mark = -1;
	this->polymer_unit_list_end.push_front(polymer_unit);
	/*Anchor_Unit anchor_unit;
	this->anchor_list_end.push_front(anchor_unit);
	Polymer_Bundle polymer_bundle;
	this->polymer_bundle_list_end.push_front(polymer_bundle);
	Polymer_Cluster polymer_cluster;
	this->polymer_cluster_list_end.push_front(polymer_cluster);*/
	//Mark.resize(column_num, vector<int>(row_num));
	cout << "finished" << endl;
}