#pragma once
#include <vector>
#include <list>
#include "Anchor_Unit.h"
#include "Polymer_Unit.h"
#include "Polymer_Cluster.h"
#include "Rule_Parameters.h"
#include "Rule_Structure.h"
#include "Polymer_Bundle.h"
#include "test.h"
#include "Rule_Polymerize_Depolymerize.h"
#include "Rule_Diffusion.h"
#include "Rule_Fragmentation_Annealing.h"
#include "Rule_Hydrolysis.h"
#include "Rule_Bundling_Debundling.h"
#include "Rule_Anchoring_Deanchoring.h"
#include "Rule_Roration.h"
#include "General_Functions.h"
#include <iostream>
#include <fstream>
#include <map>
#include <chrono>

using namespace std;


class Grid_Unit
{
	public:
		Grid_Unit();
		void output();
		list<Polymer_Unit> polymer_unit_initial;
		int row_position;
		int col_position;
		//int Is_Anchored;
		/*Anchor_Unit a;
		Anchor_Unit b;*/
		list<Anchor_Unit>::iterator grid_anchor_pointer;
		
		//list<Polymer_Cluster>::iterator grid_polymer_cluster_pointer;
		
		list<Polymer_Unit>::iterator grid_polymer_unit_pointer;//=polymer_unit_initial.begin();//problem

		//list<Polymer_Bundle>::iterator grid_polymer_bundle_pointer;
private:
		
		
};