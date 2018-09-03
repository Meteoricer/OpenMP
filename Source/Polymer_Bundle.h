#pragma once
//#include "Grid_Unit.h"


class Polymer_Cluster;
class Anchor_Unit;
class Polymer_Bundle
{
public:
	//Grid_Unit *polymer_grid_pointer;
	list<list<Polymer_Cluster>::iterator> bundle_cluster_pointer_list;
	list<list<Anchor_Unit>::iterator> polymer_anchor_list;
	//list<int> diffusion_direction;
	vector<int> availabe_diffusion_direction;
	int available_polymer_diffusion_direction_sum;
	Polymer_Bundle();
	double Consider_mark;
	int bundle_direction;
	/*~Polymer_Bundle();*/

};