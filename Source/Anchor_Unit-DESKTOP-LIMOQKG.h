#pragma once
//#include "Grid_Unit.h"
//#include "Polymer_Unit.h"
#include <vector>
#include <list>
using namespace std;
class Grid_Unit;
class Polymer_Unit;
class Polymer_Cluster;
class Anchor_Unit
{
public:
	Grid_Unit *anchor_grid_pointer;
	vector<int> availabe_diffusion_direction;
	int available_polymer_diffusion_direction_sum;
	list<Polymer_Unit>::iterator polymer_unit_pointer;
	Anchor_Unit();
	int consider_mark=0;
	int attatch_mark = 0;
	//list<Polymer_Cluster>::iterator polymer_cluster_pointer;

	
	//int IsAnchored;
	//int anchor_no;

};