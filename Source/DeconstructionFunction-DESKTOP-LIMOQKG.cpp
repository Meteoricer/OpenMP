#include "General_Functions.h"
#include <vector>
#include <list>
#include <iterator>
#include "Grid_Unit.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include "parameters.h"
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
class Polymer_Unit;
class Anchor_Unit;

void delete_polymer_unit(list<Polymer_Cluster>::iterator polymer_cluster_iter, list<Polymer_Unit>::iterator polymer_unit_iter)
{



	
	rule_parameters.total_ftsZ_lifetime = rule_parameters.t - polymer_unit_iter->lifetime;
	//remove it's neighbor pointer
	int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
	int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
	int up_col_position = col_position + 1;
	if (up_col_position > rule_parameters.column_num - 1)
		up_col_position = 0;
	int down_col_position = col_position - 1;
	if (down_col_position < 0)
		down_col_position = rule_parameters.column_num - 1;
	int right_row_position = row_position + 1;
	if (right_row_position > rule_parameters.row_num - 1)
		right_row_position = 0;
	int left_row_position = row_position - 1;
	if (left_row_position < 0)
		left_row_position = rule_parameters.row_num - 1;

	//consider anchor and attatch information
	if (Grid[col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
	{
		if (Grid[col_position][row_position].grid_anchor_pointer->attatch_mark == 1)
			//it can do attatch
		{//since it is a deconstruction function of a polymer_unit, there must be a polymer_unit here
		 //a polymer_unit linked to a anchor can't be deleted, it must first be deanchored
			rule_parameters.emtpy_anchor_attatch_sum--;
			Grid[col_position][row_position].grid_anchor_pointer->attatch_mark = 0;
		}
	}
	


	
	{
		if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator anneal_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second];
			if (anneal_iter->first->Hydrolysis_mark == T&&anneal_iter->second->Hydrolysis_mark == T)
			{
				rule_parameters.TT_annealing_sum--;
			}
			else if (anneal_iter->first->Hydrolysis_mark == D&&anneal_iter->second->Hydrolysis_mark == D)
			{

			}
			else
			{
				rule_parameters.TD_annealing_sum--;
			}
			anneal_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
			anneal_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
			rule_structure.anealing_polymer_unit_pair_list.erase(anneal_iter);

		}
		if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator anneal_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
			if (anneal_iter->first->Hydrolysis_mark == T&&anneal_iter->second->Hydrolysis_mark == T)
			{
				rule_parameters.TT_annealing_sum--;
			}
			else if (anneal_iter->first->Hydrolysis_mark == D&&anneal_iter->second->Hydrolysis_mark == D)
			{

			}
			else
			{
				rule_parameters.TD_annealing_sum--;
			}
			anneal_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
			anneal_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
			rule_structure.anealing_polymer_unit_pair_list.erase(anneal_iter);
		}
		if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.bundling_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second];
			if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
			{
				rule_parameters.bundling_sum--;
			}
			else
			{
				rule_parameters.bundling_inverse_sum--;
			}
			bundle_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
			bundle_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
			rule_structure.bundling_polymer_unit_pair_list.erase(bundle_pair_iter);
		}
		if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.bundling_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
			if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
			{
				rule_parameters.bundling_sum--;
			}
			else
			{
				rule_parameters.bundling_inverse_sum--;
			}
			bundle_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
			bundle_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
			rule_structure.bundling_polymer_unit_pair_list.erase(bundle_pair_iter);
		}
		if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundle_pair_iter = polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first];
			if (debundle_pair_iter->first->polymer_cluster_pointer->polarize == debundle_pair_iter->second->polymer_cluster_pointer->polarize)
			{
				rule_parameters.debundling_sum--;
			}
			else
			{
				rule_parameters.debundling_inverse_sum--;
			}
			debundle_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
			debundle_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
			rule_structure.debundling_polymer_unit_pair_list.erase(debundle_pair_iter);
		}

		if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundle_pair_iter = polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second];
			if (debundle_pair_iter->first->polymer_cluster_pointer->polarize == debundle_pair_iter->second->polymer_cluster_pointer->polarize)
			{
				rule_parameters.debundling_sum--;
			}
			else
			{
				rule_parameters.debundling_inverse_sum--;
			}
			debundle_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
			debundle_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
			rule_structure.debundling_polymer_unit_pair_list.erase(debundle_pair_iter);
		}
		
		//if (polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter_1 != rule_structure.TT_anealing_polymer_unit_pair_list.end())
		//{
		//	rule_structure.TT_anealing_polymer_unit_pair_list.erase(polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter_1);
		//}
		//if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter_1 != rule_structure.bundling_polymer_unit_pair_list.end())
		//{
		//	rule_structure.bundling_polymer_unit_pair_list.erase(polymer_unit_iter->bundling_polymer_unit_pair_list_iter_1);
		//}
		///*if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter_1 != rule_structure.debundling_polymer_unit_pair_list.end())
		//{
		//rule_structure.debundling_polymer_unit_pair_list.erase(polymer_unit_iter->debundling_polymer_unit_pair_list_iter_1);
		//}*/
		//if (polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter_2 != rule_structure.TT_anealing_polymer_unit_pair_list.end())
		//{
		//	rule_structure.TT_anealing_polymer_unit_pair_list.erase(polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter_2);
		//}
		//if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter_2 != rule_structure.bundling_polymer_unit_pair_list.end())
		//{
		//	rule_structure.bundling_polymer_unit_pair_list.erase(polymer_unit_iter->bundling_polymer_unit_pair_list_iter_2);
		//}
		///*if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter_2 != rule_structure.debundling_polymer_unit_pair_list.end())
		//{
		//rule_structure.debundling_polymer_unit_pair_list.erase(polymer_unit_iter->debundling_polymer_unit_pair_list_iter_2);
		//}*/
	}


	polymer_unit_iter->polymer_grid_pointer->grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();

	//delete its relation with bundle//there is no bundle to unit or unit to bundle iterator, nothing to be done
	//delete its relation with anchor//a polymer with anchor will not be depolymerized,it will more likely to be framentation

	list<Polymer_Bundle>::iterator bundle_iter_up;
	if (&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		bundle_iter_up = Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
	else
		bundle_iter_up = rule_structure.polymer_bundle_list.end();


	list<Polymer_Bundle>::iterator bundle_iter_down;
	if (&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		bundle_iter_down = Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
	else
		bundle_iter_down = rule_structure.polymer_bundle_list.end();


	list<Polymer_Bundle>::iterator bundle_iter_left;
	if (&*Grid[col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		bundle_iter_left = Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
	else
		bundle_iter_left = rule_structure.polymer_bundle_list.end();


	list<Polymer_Bundle>::iterator bundle_iter_right;
	if (&*Grid[col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		bundle_iter_right = Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
	else
		bundle_iter_right = rule_structure.polymer_bundle_list.end();





	if (bundle_iter_up != rule_structure.polymer_bundle_list.end() &&
		bundle_iter_up != polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer)
	{//if they do not belong to the same bundle, need to check bundle diffusion relation
		check_polymer_bundle_diffusable_direction(bundle_iter_up,1);
	}
	if (bundle_iter_down != rule_structure.polymer_bundle_list.end() &&
		bundle_iter_down != polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer &&
		bundle_iter_up != bundle_iter_down)
	{
		check_polymer_bundle_diffusable_direction(bundle_iter_down,1);
	}
	if (bundle_iter_left != rule_structure.polymer_bundle_list.end() &&
		bundle_iter_left != polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer &&
		bundle_iter_left != bundle_iter_down && bundle_iter_left != bundle_iter_up)
	{
		check_polymer_bundle_diffusable_direction(bundle_iter_left,1);
	}
	if (bundle_iter_right != rule_structure.polymer_bundle_list.end() &&
		bundle_iter_right != polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer &&
		bundle_iter_right != bundle_iter_down && bundle_iter_right != bundle_iter_up&& bundle_iter_right != bundle_iter_left)
	{
		check_polymer_bundle_diffusable_direction(bundle_iter_right,1);
	}
	polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
};

void delete_polymer_unit(list<Polymer_Cluster>::iterator polymer_cluster_iter, list<Polymer_Unit>::iterator polymer_unit_iter_start, list<Polymer_Unit>::iterator polymer_unit_iter_end)
{
	list<Polymer_Unit>::iterator polymer_unit_iter=polymer_unit_iter_start;
	list<Polymer_Unit>::iterator polymer_unit_iter_temp = polymer_unit_iter_start;
	while (polymer_unit_iter != polymer_unit_iter_end)
	{
		polymer_unit_iter_temp++;
		delete_polymer_unit(polymer_cluster_iter, polymer_unit_iter);
		polymer_unit_iter=polymer_unit_iter_temp;
	}
};

void delete_polymer_cluster(list<Polymer_Cluster>::iterator polymer_cluster_iter)
{
	//check_polymer_hydrolysis_sum(polymer_cluster_iter->polymer_sequence.begin());
	check_polymer_hydrolysis_sum(polymer_cluster_iter);
	rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - polymer_cluster_iter->TT_depolymerize_plus_end_sum;
	rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - polymer_cluster_iter->TD_depolymerize_plus_end_sum;
	rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - polymer_cluster_iter->DD_depolymerize_plus_end_sum;
	rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - polymer_cluster_iter->TT_depolymerize_minus_end_sum;
	rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - polymer_cluster_iter->TD_depolymerize_minus_end_sum;
	rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - polymer_cluster_iter->DD_depolymerize_minus_end_sum;
	rule_parameters.DD_fragmentation_sum -= polymer_cluster_iter->DD_sum;
	rule_parameters.TD_fragmentation_sum -= polymer_cluster_iter->TD_sum;
	
	rule_parameters.TT_fragmentation_sum -= polymer_cluster_iter->TT_sum;
	rule_parameters.TTT_sum_total -= polymer_cluster_iter->TTT_sum;
	rule_parameters.TTD_sum_total -= polymer_cluster_iter->TTD_sum;
	rule_parameters.DTD_sum_total -= polymer_cluster_iter->DTD_sum;
	if (polymer_cluster_iter->direction == horizontal)
	{
		rule_parameters.TT_polymerize_end_horizontal_plus_sum -= polymer_cluster_iter->TT_polymerize_plus_end_sum;
		rule_parameters.TT_polymerize_end_horizontal_minus_sum -= polymer_cluster_iter->TT_polymerize_minus_end_sum;
		rule_parameters.TD_polymerize_end_horizontal_plus_sum -= polymer_cluster_iter->TD_polymerize_plus_end_sum;
		rule_parameters.TD_polymerize_end_horizontal_minus_sum -= polymer_cluster_iter->TD_polymerize_minus_end_sum;
	}
	if (polymer_cluster_iter->direction == vertical)
	{
		rule_parameters.TT_polymerize_end_vertical_plus_sum -= polymer_cluster_iter->TT_polymerize_plus_end_sum;
		rule_parameters.TT_polymerize_end_vertical_minus_sum -= polymer_cluster_iter->TT_polymerize_minus_end_sum;
		rule_parameters.TD_polymerize_end_vertical_plus_sum -= polymer_cluster_iter->TD_polymerize_plus_end_sum;
		rule_parameters.TD_polymerize_end_vertical_minus_sum -= polymer_cluster_iter->TD_polymerize_minus_end_sum;
	}
	
	delete_polymer_unit(polymer_cluster_iter, polymer_cluster_iter->polymer_sequence.begin(), polymer_cluster_iter->polymer_sequence.end());
	rule_structure.polymer_cluster_list.erase(polymer_cluster_iter);
};
void delete_polymer_bundle(list<Polymer_Bundle>::iterator polymer_bundle_iter)
{
	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_seq_iter;
	list<Polymer_Cluster>::iterator polymer_cluster_iter;
	polymer_cluster_seq_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
	if(polymer_cluster_seq_iter!= polymer_bundle_iter->bundle_cluster_pointer_list.end())
		polymer_cluster_iter = (*polymer_cluster_seq_iter);
		//rule_structure.polymer_cluster_list.erase(rule_structure.polymer_cluster_list.begin());
		/*if (polymer_bundle_iter->bundle_cluster_pointer_list.size() > 0)
		{
		rule_structure.polymer_cluster_list.erase(*polymer_bundle_iter->bundle_cluster_pointer_list.begin(), *(polymer_bundle_iter->bundle_cluster_pointer_list.end()));
		}*/

		//rule_structure.polymer_cluster_list.erase(*polymer_cluster_iter);
		while (polymer_cluster_seq_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
		{
			//polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			delete_polymer_cluster(polymer_cluster_iter);
			//rule_structure.polymer_cluster_list.erase(polymer_cluster_iter);
			polymer_bundle_iter->bundle_cluster_pointer_list.erase(polymer_cluster_seq_iter);
			polymer_cluster_seq_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
			if (polymer_cluster_seq_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
				polymer_cluster_iter = (*polymer_cluster_seq_iter);
			//i don't know how to make here work, there is a iterator incompatitable issue, now go with the stupid way
		}
		//rule_parameters.polymer_diffusion_direction_sum = rule_parameters.polymer_diffusion_direction_sum - polymer_bundle_iter->available_polymer_diffusion_direction_sum;
		rule_parameters.polymer_diffusion_direction_sum -= polymer_bundle_iter->available_polymer_diffusion_direction_sum;
		rule_structure.polymer_bundle_list.erase(polymer_bundle_iter);
};