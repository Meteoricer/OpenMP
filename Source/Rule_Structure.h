#pragma once
#include "Grid_Unit.h"
//#include "Diffusion_Unit.h"
//#include "Polymer_Unit.h"
#include <vector>
#include <list>
//#include "parameters.h"
//#include "Anchor_Unit.h"
using namespace std;



class Rule_Structure
{
public:
	//list<Grid_Unit> anchor_empty_list;
	//list<Grid_Unit> anchor_diffusion_list;
	list<Anchor_Unit> anchor_list;
	list<Polymer_Cluster> polymer_cluster_list;
	list<Polymer_Bundle> polymer_bundle_list;
	//list<Anchor_Unit> anchor_list_end;
	//list<Polymer_Cluster> polymer_cluster_list_end;
	//list<Polymer_Bundle> polymer_bundle_list_end;
	list<Polymer_Unit> polymer_unit_list_end;//do not do anything meaningful, just to put null pointer at the begining.
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> bundle_candidate_list;// i still need to construct this list in diffusion, 
	
	
	// i should use a map for those pairs or a list?
	// i need to find a method to search element quickly, if i use map, what is the key value?
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> anealing_polymer_unit_pair_list;
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> bundling_polymer_unit_pair_list;
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> bundling_polymer_unit_pair_inverse_list;
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> debundling_polymer_unit_pair_list;
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> debundling_polymer_unit_pair_inverse_list;
	//then i need to fix map no splice problem
	//multimap<int,pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>> TT_anealing_polymer_unit_pair_list;
	//tried multimap, too complecated to deal with iterators, decided to go the brutal way

	//rebuild the boundary condition in fragmentation.
	Rule_Structure();
	//vector<vector<int>> Mark;





	
};