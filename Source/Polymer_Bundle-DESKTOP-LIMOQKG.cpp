#include "Rule_Anchoring_Deanchoring.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include "Grid_Unit.h"
#include <vector>
#include <cmath>
#include <random>
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
Polymer_Bundle::Polymer_Bundle()
{
	this->availabe_diffusion_direction.resize(4);
	this->availabe_diffusion_direction[direction_up] = false;
	this->availabe_diffusion_direction[direction_down] = false;
	this->availabe_diffusion_direction[direction_left] = false;
	this->availabe_diffusion_direction[direction_right] = false;
	this->Consider_mark = -1;
}
//Polymer_Bundle::~Polymer_Bundle()
//{
//	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_seq_iter;
//	list<Polymer_Cluster>::iterator polymer_cluster_iter;
//
//
//	if (this->bundle_cluster_pointer_list.size() > 0)
//	{
//		polymer_cluster_seq_iter = this->bundle_cluster_pointer_list.begin();
//		polymer_cluster_iter = (*polymer_cluster_seq_iter);
//		//rule_structure.polymer_cluster_list.erase(rule_structure.polymer_cluster_list.begin());
//		/*if (this->bundle_cluster_pointer_list.size() > 0)
//		{
//			rule_structure.polymer_cluster_list.erase(*this->bundle_cluster_pointer_list.begin(), *(this->bundle_cluster_pointer_list.end()));
//		}*/
//
//		//rule_structure.polymer_cluster_list.erase(*polymer_cluster_iter);
//		while (polymer_cluster_seq_iter != this->bundle_cluster_pointer_list.end())
//		{
//			polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
//			rule_structure.polymer_cluster_list.erase(polymer_cluster_iter);
//			this->bundle_cluster_pointer_list.erase(polymer_cluster_seq_iter);
//			polymer_cluster_seq_iter = this->bundle_cluster_pointer_list.begin();
//			if (polymer_cluster_seq_iter != this->bundle_cluster_pointer_list.end())
//				polymer_cluster_iter = (*polymer_cluster_seq_iter);
//			//i don't know how to make here work, there is a iterator incompatitable issue, now go with the stupid way
//		}
//		//rule_parameters.polymer_diffusion_direction_sum = rule_parameters.polymer_diffusion_direction_sum - this->available_polymer_diffusion_direction_sum;
//		rule_parameters.polymer_diffusion_direction_sum -= this->available_polymer_diffusion_direction_sum;
//	}
//}

