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
#include "DeconstructionFunction.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

void bundling()
{
	///*//random_device rd;
	//mt19937_64 gen(rd());*/
	//default_random_engine gen;
	//uniform_int_distribution<> Bundling_Rand(0, rule_parameters.TT_bundling_sum);
	//int bundling_rand = Bundling_Rand(gen);
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_polymer_unit_pair_iter;
	//bundle_polymer_unit_pair_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
	//int bundle_temp_sum = 0;
	////anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	//while (bundle_temp_sum < bundling_rand)
	//{
	//	bundle_polymer_unit_pair_iter++;
	//}
	//if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer ==
	//	bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer)
	//{//two polymer_unit already belong to the same bundle, don't need to do anything

	//}
	//else
	//{//do not belong to the same bundle, need to reconstruct bundle
	//	//first link polymer_unit together
	//	int col_dist = bundle_polymer_unit_pair_iter->first->polymer_grid_pointer->col_position -
	//		bundle_polymer_unit_pair_iter->second->polymer_grid_pointer->col_position;
	//	int row_dist= bundle_polymer_unit_pair_iter->first->polymer_grid_pointer->row_position -
	//		bundle_polymer_unit_pair_iter->second->polymer_grid_pointer->row_position;
	//	if (col_dist == column_num - 1)
	//	{
	//		col_dist = -1;
	//	}
	//	if (col_dist == 1 - column_num)
	//	{
	//		col_dist = 1;
	//	}
	//	if (row_dist == row_num - 1)
	//	{
	//		row_dist = -1;
	//	}
	//	if (row_dist == 1 - row_num)
	//	{
	//		row_dist = 1;
	//	}
	//	//if (col_dist == 1)
	//	//{//the first one is up than second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_down] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_up] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//	
	//	//}
	//	//if (col_dist == -1)
	//	//{//the first one is downside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_up] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_down] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//if (row_dist == 1)
	//	//{//the first one is rightside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_left] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_right] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//if (row_dist == -1)
	//	//{//the first one is rightside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_right] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_left] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//reconstruct relation for first bundle		
	//	//check_polymer_bundle_diffusable_direction(bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer);
	//	
	//	list<Polymer_Bundle>::iterator to_be_deleted_bundle_iter;
	//	to_be_deleted_bundle_iter = bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer;

	//	check_polymer_bundle_diffusable_direction(bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer);
	//	
	//	//splice all the polymer_cluster from second bundle to first bundle
	//	bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.splice(
	//		bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin(),
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list,
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin(),
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.end());
	//	//now second bundle contains nothing
	//	//maybe i need to delete the second bundle here

	//	rule_structure.polymer_bundle_list.erase(to_be_deleted_bundle_iter);
	//	//move this pair from bundling list to debundling list
	//	rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
	//}

	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//cout << " enter bundling" << endl;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	//default_random_engine gen;
	uniform_int_distribution<> bundling_Rand(1, rule_parameters.bundling_sum);
	int bundle_rand = bundling_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_polymer_unit_pair_iter;
	bundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();


	int bundle_temp_sum = 0;
	if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		bundle_temp_sum = 1;
	}

	
	while (bundle_temp_sum < bundle_rand)
	{
		if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			bundle_temp_sum++;
		}
		if (bundle_temp_sum == bundle_rand)
		{
			break;
		}
		
		bundle_polymer_unit_pair_iter++;
		if (bundle_polymer_unit_pair_iter == rule_structure.bundling_polymer_unit_pair_list.end())
		{
			cout << "out of bound" << endl;
			system("pause");
		}
	}
	if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		cout << "bunlding information error";
		system("pause");
	}


	
	//anneal_polymer_unit_iter--;
	list<Polymer_Bundle>::iterator lower_bundle_iter;//it can be on the downside and left side
	list<Polymer_Bundle>::iterator higher_bundle_iter;//it can be on the upperside and right side
														//list<Polymer_Unit>::iterator large_polymer_unit_begin;
														//list<Polymer_Unit>::iterator large_polymer_unit_end;
														//i already define the direction of pair, so don't need to determine directino here

	lower_bundle_iter = bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer;
	higher_bundle_iter = bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer;
	/*if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->direction == vertical)
	{
		cout << "bundle verticle" << endl;
	}*/

	/*if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}
	if (bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}*/


	//now judge which part is larger, and insert the samller part into the larger part
	//if (lower_bundle_iter->polymer_sequence.size()>higher_bundle_iter->polymer_sequence.size())
	//the lower part is larger//i don't judge it any more, because there is a lot of bundling need to do here, only one polymer_sequence makes no difference

	
	if (lower_bundle_iter != higher_bundle_iter)
	{//if it is not the same bundle

		list<Polymer_Cluster>::iterator polymer_cluster_iter;
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_cluster_iter_iter = higher_bundle_iter->bundle_cluster_pointer_list.begin();
		polymer_cluster_iter = *polymer_cluster_iter_iter;
		do
		{
			//polymer_unit_iter->polymer_cluster_pointer = lower_bundle_iter;
			//polymer_unit_iter++;
			(*polymer_cluster_iter_iter)->cluster_bundle_pointer = lower_bundle_iter;
			polymer_cluster_iter_iter++;
		} while (polymer_cluster_iter_iter != higher_bundle_iter->bundle_cluster_pointer_list.end());


		lower_bundle_iter->bundle_cluster_pointer_list.splice(lower_bundle_iter->bundle_cluster_pointer_list.end(), higher_bundle_iter->bundle_cluster_pointer_list, higher_bundle_iter->bundle_cluster_pointer_list.begin(), higher_bundle_iter->bundle_cluster_pointer_list.end());
		lower_bundle_iter->polymer_anchor_list.splice(lower_bundle_iter->polymer_anchor_list.end(), higher_bundle_iter->polymer_anchor_list, higher_bundle_iter->polymer_anchor_list.begin(), higher_bundle_iter->polymer_anchor_list.end());
		
		delete_polymer_bundle(higher_bundle_iter);


		bundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
		bundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
		rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
		bundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
		bundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = bundle_polymer_unit_pair_iter;
		bundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = bundle_polymer_unit_pair_iter;
		rule_parameters.bundling_sum--;
		rule_parameters.debundling_sum++;
		//i don't need this part, because for the same bundle bundling, it will rebuild the bundle everytime, and check diffusion, but it is not necessary
		recursion_polymer_unit(bundle_polymer_unit_pair_iter->first, lower_bundle_iter);//here is the problem
		check_polymer_bundle_diffusable_direction(lower_bundle_iter, 1);//here is the problem
		
	}
	else
	{//if it is the same bundle
		bundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
		bundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
		rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
		bundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
		bundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = bundle_polymer_unit_pair_iter;
		bundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = bundle_polymer_unit_pair_iter;
		rule_parameters.bundling_sum--;
		rule_parameters.debundling_sum++;
	}

	//first write a consider mark function to mark one bundle
	//it they are the same, just break the debundle
	//if not, creat a new bundle, modify the recursion functino to creat bundle_pointer_iter

	/*if (lower_bundle_iter->polymer_anchor_list.size() > 1)
	{
		lower_bundle_iter->available_polymer_diffusion_direction_sum = 0;
	}*/
	//check_polymer_hydrolysis_sum(lower_bundle_iter-> polymer_sequence.begin());
	

	//rule_structure.TT_anealing_polymer_unit_pair_list.erase(bundle_polymer_unit_pair_iter);
	//cout << "bundling finished" << endl;
	

	

};


void bundling_inverse()
{
	///*//random_device rd;
	//mt19937_64 gen(rd());*/
	//default_random_engine gen;
	//uniform_int_distribution<> Bundling_Rand(0, rule_parameters.TT_bundling_sum);
	//int bundling_rand = Bundling_Rand(gen);
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_polymer_unit_pair_iter;
	//bundle_polymer_unit_pair_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
	//int bundle_temp_sum = 0;
	////anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	//while (bundle_temp_sum < bundling_rand)
	//{
	//	bundle_polymer_unit_pair_iter++;
	//}
	//if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer ==
	//	bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer)
	//{//two polymer_unit already belong to the same bundle, don't need to do anything

	//}
	//else
	//{//do not belong to the same bundle, need to reconstruct bundle
	//	//first link polymer_unit together
	//	int col_dist = bundle_polymer_unit_pair_iter->first->polymer_grid_pointer->col_position -
	//		bundle_polymer_unit_pair_iter->second->polymer_grid_pointer->col_position;
	//	int row_dist= bundle_polymer_unit_pair_iter->first->polymer_grid_pointer->row_position -
	//		bundle_polymer_unit_pair_iter->second->polymer_grid_pointer->row_position;
	//	if (col_dist == column_num - 1)
	//	{
	//		col_dist = -1;
	//	}
	//	if (col_dist == 1 - column_num)
	//	{
	//		col_dist = 1;
	//	}
	//	if (row_dist == row_num - 1)
	//	{
	//		row_dist = -1;
	//	}
	//	if (row_dist == 1 - row_num)
	//	{
	//		row_dist = 1;
	//	}
	//	//if (col_dist == 1)
	//	//{//the first one is up than second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_down] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_up] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//	
	//	//}
	//	//if (col_dist == -1)
	//	//{//the first one is downside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_up] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_down] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//if (row_dist == 1)
	//	//{//the first one is rightside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_left] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_right] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//if (row_dist == -1)
	//	//{//the first one is rightside of the second one
	//	//	bundle_polymer_unit_pair_iter->first->bundled_neighbor_polymer_unit_pointer[direction_right] =
	//	//		bundle_polymer_unit_pair_iter->second;
	//	//	bundle_polymer_unit_pair_iter->second->bundled_neighbor_polymer_unit_pointer[direction_left] =
	//	//		bundle_polymer_unit_pair_iter->first;
	//	//}
	//	//reconstruct relation for first bundle		
	//	//check_polymer_bundle_diffusable_direction(bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer);
	//	
	//	list<Polymer_Bundle>::iterator to_be_deleted_bundle_iter;
	//	to_be_deleted_bundle_iter = bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer;

	//	check_polymer_bundle_diffusable_direction(bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer);
	//	
	//	//splice all the polymer_cluster from second bundle to first bundle
	//	bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.splice(
	//		bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin(),
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list,
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin(),
	//		bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.end());
	//	//now second bundle contains nothing
	//	//maybe i need to delete the second bundle here

	//	rule_structure.polymer_bundle_list.erase(to_be_deleted_bundle_iter);
	//	//move this pair from bundling list to debundling list
	//	rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
	//}

	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//cout << " enter bundling" << endl;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	//default_random_engine gen;
	uniform_int_distribution<> bundling_Rand(1, rule_parameters.bundling_inverse_sum);
	int bundle_rand = bundling_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_polymer_unit_pair_iter;
	bundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
	int bundle_temp_sum = 0;
	if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		bundle_temp_sum = 1;
	}
	
	//anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	while (bundle_temp_sum < bundle_rand)
	{
		if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			bundle_temp_sum++;
		}		
		//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum;
		if (bundle_temp_sum == bundle_rand)
		{
			break;
		}
		bundle_polymer_unit_pair_iter++;
		if (bundle_polymer_unit_pair_iter == rule_structure.bundling_polymer_unit_pair_list.end())
		{
			cout << "out of bound" << endl;
			system("pause");
		}
	}
	if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		cout<<"bunlding information error";
		system("pause");
	}
	//anneal_polymer_unit_iter--;
	list<Polymer_Bundle>::iterator lower_bundle_iter;//it can be on the downside and left side
	list<Polymer_Bundle>::iterator higher_bundle_iter;//it can be on the upperside and right side
													  //list<Polymer_Unit>::iterator large_polymer_unit_begin;
													  //list<Polymer_Unit>::iterator large_polymer_unit_end;
													  //i already define the direction of pair, so don't need to determine directino here

	lower_bundle_iter = bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer;
	higher_bundle_iter = bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer;
	/*if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->direction == vertical)
	{
	cout << "bundle verticle" << endl;
	}*/

	/*if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
	if (bundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() == 1)
	{
	rule_parameters.rotation_sum--;
	}
	}
	if (bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
	if (bundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size() == 1)
	{
	rule_parameters.rotation_sum--;
	}
	}*/


	//now judge which part is larger, and insert the samller part into the larger part
	//if (lower_bundle_iter->polymer_sequence.size()>higher_bundle_iter->polymer_sequence.size())
	//the lower part is larger//i don't judge it any more, because there is a lot of bundling need to do here, only one polymer_sequence makes no difference


	if (lower_bundle_iter != higher_bundle_iter)
	{//if it is not the same bundle

		list<Polymer_Cluster>::iterator polymer_cluster_iter;
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_cluster_iter_iter = higher_bundle_iter->bundle_cluster_pointer_list.begin();
		polymer_cluster_iter = *polymer_cluster_iter_iter;
		do
		{
			//polymer_unit_iter->polymer_cluster_pointer = lower_bundle_iter;
			//polymer_unit_iter++;
			(*polymer_cluster_iter_iter)->cluster_bundle_pointer = lower_bundle_iter;
			polymer_cluster_iter_iter++;
		} while (polymer_cluster_iter_iter != higher_bundle_iter->bundle_cluster_pointer_list.end());


		lower_bundle_iter->bundle_cluster_pointer_list.splice(lower_bundle_iter->bundle_cluster_pointer_list.end(), higher_bundle_iter->bundle_cluster_pointer_list, higher_bundle_iter->bundle_cluster_pointer_list.begin(), higher_bundle_iter->bundle_cluster_pointer_list.end());
		lower_bundle_iter->polymer_anchor_list.splice(lower_bundle_iter->polymer_anchor_list.end(), higher_bundle_iter->polymer_anchor_list, higher_bundle_iter->polymer_anchor_list.begin(), higher_bundle_iter->polymer_anchor_list.end());

		delete_polymer_bundle(higher_bundle_iter);


		bundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
		bundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
		rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
		bundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
		bundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = bundle_polymer_unit_pair_iter;
		bundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = bundle_polymer_unit_pair_iter;
		rule_parameters.bundling_inverse_sum--;
		rule_parameters.debundling_inverse_sum++;
		//i don't need this part, because for the same bundle bundling, it will rebuild the bundle everytime, and check diffusion, but it is not necessary
		recursion_polymer_unit(bundle_polymer_unit_pair_iter->first, lower_bundle_iter);//here is the problem
		check_polymer_bundle_diffusable_direction(lower_bundle_iter, 1);//here is the problem

	}
	else
	{//if it is the same bundle
		bundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.end();
		bundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.end();
		rule_structure.debundling_polymer_unit_pair_list.splice(rule_structure.debundling_polymer_unit_pair_list.begin(), rule_structure.bundling_polymer_unit_pair_list, bundle_polymer_unit_pair_iter);
		bundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
		bundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = bundle_polymer_unit_pair_iter;
		bundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = bundle_polymer_unit_pair_iter;
		rule_parameters.bundling_inverse_sum--;
		rule_parameters.debundling_inverse_sum++;
	}

	//first write a consider mark function to mark one bundle
	//it they are the same, just break the debundle
	//if not, creat a new bundle, modify the recursion functino to creat bundle_pointer_iter

	/*if (lower_bundle_iter->polymer_anchor_list.size() > 1)
	{
	lower_bundle_iter->available_polymer_diffusion_direction_sum = 0;
	}*/
	//check_polymer_hydrolysis_sum(lower_bundle_iter-> polymer_sequence.begin());


	//rule_structure.TT_anealing_polymer_unit_pair_list.erase(bundle_polymer_unit_pair_iter);
	//cout << "bundling finished" << endl;




};








void debundling()
{
	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//cout << " enter debundling" << endl;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	//default_random_engine gen;
	uniform_int_distribution<> debundling_Rand(1, rule_parameters.debundling_sum);
	int debundle_rand = debundling_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundle_polymer_unit_pair_iter;
	debundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
	//int debundle_temp_sum = 1;
	////anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	//while (debundle_temp_sum < debundle_rand)
	//{
	//	debundle_temp_sum++;
	//	//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum;
	//	debundle_polymer_unit_pair_iter++;
	//}
	////anneal_polymer_unit_iter--;


	int debundle_temp_sum = 0;
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		debundle_temp_sum = 1;
	}

	//anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	while (debundle_temp_sum < debundle_rand)
	{
		if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			debundle_temp_sum++;
		}
		//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum;
		if (debundle_temp_sum == debundle_rand)
		{
			break;
		}
		debundle_polymer_unit_pair_iter++;
		if (debundle_polymer_unit_pair_iter == rule_structure.debundling_polymer_unit_pair_list.end())
		{
			cout << "out of bound" << endl;
			system("pause");
		}
	}
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		cout << "bunlding information error";
		system("pause");
	}








	list<Polymer_Unit>::iterator first_polymer_unit_iter;//it can be on the downside and left side
	list<Polymer_Unit>::iterator second_polymer_unit_iter;//it can be on the upperside and right side
													  //list<Polymer_Unit>::iterator large_polymer_unit_begin;
													  //list<Polymer_Unit>::iterator large_polymer_unit_end;
													  //i already define the direction of pair, so don't need to determine directino here

	first_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
	second_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;

	int mark = 0;
	auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
	auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
	//int direction_mark = 0;
	//direction_mark==0 means left is larger than right, otherwise means right is larger than left
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() >= debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size())
	{
		//auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
		//auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
		auto polymer_unit_iter = small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
		while (polymer_unit_iter != small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
		{
			if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
			{
				if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first]->first->polymer_cluster_pointer == large_polymer_unit_iter->polymer_cluster_pointer &&
					polymer_unit_iter != small_polymer_unit_iter)
				{
					mark = 1;//mark==1, it will still be the same bundle
					break;
				}
			}
				
			polymer_unit_iter++;
		}
		
	}
	else
	{
		auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
		auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
		//direction_mark = 1;


		auto polymer_unit_iter = small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
		while (polymer_unit_iter != small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
		{
			if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
			{
				if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second]->second->polymer_cluster_pointer == large_polymer_unit_iter->polymer_cluster_pointer &&
					polymer_unit_iter != small_polymer_unit_iter)
				{
					mark = 1;//mark==1, it will still be the same bundle
					break;
				}
			}
			
			polymer_unit_iter++;
		}

	}
	
	//auto polymer_unit_iter = first_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//while (polymer_unit_iter != polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
	//{
	//	if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second]->second->polymer_cluster_pointer == second_polymer_unit_iter->polymer_cluster_pointer &&
	//		polymer_unit_iter != first_polymer_unit_iter)
	//	{
	//		mark = 1;//mark==1, it will still be the same bundle
	//		break;
	//	}
	//}
	if (mark == 1)
	{//if mark==1, then it is still belong to the same bundle, don't have to go through this complex case (but if mark!=1, it could still be the same bundle)
		debundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
		debundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
		rule_structure.bundling_polymer_unit_pair_list.splice(rule_structure.bundling_polymer_unit_pair_list.begin(), rule_structure.debundling_polymer_unit_pair_list, debundle_polymer_unit_pair_iter);
		debundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
		debundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = debundle_polymer_unit_pair_iter;
		debundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = debundle_polymer_unit_pair_iter;
		rule_parameters.debundling_sum--;
		rule_parameters.bundling_sum++;
	}
	else
	{
		list<Polymer_Bundle>::iterator old_second_bundle_iter = second_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer;



		debundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
		debundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
		rule_structure.bundling_polymer_unit_pair_list.splice(rule_structure.bundling_polymer_unit_pair_list.begin(), rule_structure.debundling_polymer_unit_pair_list, debundle_polymer_unit_pair_iter);
		debundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
		debundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = debundle_polymer_unit_pair_iter;
		debundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = debundle_polymer_unit_pair_iter;
		rule_parameters.debundling_sum--;
		rule_parameters.bundling_sum++;

		Polymer_Bundle polymer_bundle;
		rule_structure.polymer_bundle_list.push_front(polymer_bundle);
		list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		recursion_polymer_unit(first_polymer_unit_iter, polymer_bundle_iter);
		if (first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer != second_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer)
		{//if they do not belong to the same bundle
			recursion_polymer_unit(second_polymer_unit_iter, old_second_bundle_iter);
			//cout << "1" << endl;;
			if (old_second_bundle_iter->polymer_anchor_list.size() == 0)
			{
				//cout <<"2" << endl;
				delete_polymer_bundle(old_second_bundle_iter);//why run this then go to the else section?

			}
			else
			{
				//cout <<"3"<< endl;
				//check_polymer_hydrolysis_sum((*old_second_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.begin());
				//check_polymer_bundle_diffusable_direction(second_polymer_unit_iter);
				//check_polymer_hydrolysis_sum(second_polymer_unit_iter);
				check_polymer_bundle_diffusable_direction(old_second_bundle_iter, 2);
				/*if (old_second_bundle_iter->bundle_cluster_pointer_list.size() == 1)
				{
					if ((*old_second_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.size() == 1)
					{
						rule_parameters.rotation_sum++;
					}
				}*/
			}
			if (first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
			{
				delete_polymer_bundle(first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			}
			else
			{
				check_polymer_hydrolysis_sum(first_polymer_unit_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
				/*if (polymer_bundle_iter->bundle_cluster_pointer_list.size() == 1)
				{
					if ((*polymer_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.size() == 1)
					{
						rule_parameters.rotation_sum++;
					}
				}*/
			}
		}
		else
		{//if they still belong to the same bundle
			old_second_bundle_iter->bundle_cluster_pointer_list.clear();
			old_second_bundle_iter->polymer_anchor_list.clear();

			for (int i = 0; i < 4; i++)
			{
				polymer_bundle_iter->availabe_diffusion_direction[i] = old_second_bundle_iter->availabe_diffusion_direction[i];

			}
			polymer_bundle_iter->available_polymer_diffusion_direction_sum = old_second_bundle_iter->available_polymer_diffusion_direction_sum;

			rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
			delete_polymer_bundle(old_second_bundle_iter);
		}





		//cout << "debundling finished" << endl;
		//rebuild_TT_Fragmentation_sum();
		//i used to have it because there is some unkonw error here, now i need to debug it
	}
	







	
	
};









void debundling_inverse()
{
	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//cout << " enter debundling" << endl;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	//default_random_engine gen;
	uniform_int_distribution<> debundling_Rand(1, rule_parameters.debundling_inverse_sum);
	int debundle_rand = debundling_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundle_polymer_unit_pair_iter;
	debundle_polymer_unit_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
	int debundle_temp_sum = 0;
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		debundle_temp_sum = 1;
	}

	//anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	while (debundle_temp_sum < debundle_rand)
	{
		if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize != debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			debundle_temp_sum++;
		}
		//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum
		if (debundle_temp_sum == debundle_rand)
		{
			break;
		}
		debundle_polymer_unit_pair_iter++;
		if (debundle_polymer_unit_pair_iter == rule_structure.debundling_polymer_unit_pair_list.end())
		{
			cout << "out of bound" << endl;
			system("pause");
		}
	}
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polarize == debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polarize)
	{
		cout << "bunlding information error";
		system("pause");
	}
	//anneal_polymer_unit_iter--;
	list<Polymer_Unit>::iterator first_polymer_unit_iter;//it can be on the downside and left side
	list<Polymer_Unit>::iterator second_polymer_unit_iter;//it can be on the upperside and right side
														  //list<Polymer_Unit>::iterator large_polymer_unit_begin;
														  //list<Polymer_Unit>::iterator large_polymer_unit_end;
														  //i already define the direction of pair, so don't need to determine directino here

	first_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
	second_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;

	int mark = 0;
	auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
	auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
	//int direction_mark = 0;
	//direction_mark==0 means left is larger than right, otherwise means right is larger than left
	if (debundle_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() >= debundle_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size())
	{
		//auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
		//auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
		auto polymer_unit_iter = small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
		while (polymer_unit_iter != small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
		{
			if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
			{
				if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first]->first->polymer_cluster_pointer == large_polymer_unit_iter->polymer_cluster_pointer &&
					polymer_unit_iter != small_polymer_unit_iter)
				{
					mark = 1;//mark==1, it will still be the same bundle
					break;
				}
			}

			polymer_unit_iter++;
		}

	}
	else
	{
		auto large_polymer_unit_iter = debundle_polymer_unit_pair_iter->second;
		auto small_polymer_unit_iter = debundle_polymer_unit_pair_iter->first;
		//direction_mark = 1;


		auto polymer_unit_iter = small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
		while (polymer_unit_iter != small_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
		{
			if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
			{
				if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second]->second->polymer_cluster_pointer == large_polymer_unit_iter->polymer_cluster_pointer &&
					polymer_unit_iter != small_polymer_unit_iter)
				{
					mark = 1;//mark==1, it will still be the same bundle
					break;
				}
			}

			polymer_unit_iter++;
		}

	}

	//auto polymer_unit_iter = first_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//while (polymer_unit_iter != polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
	//{
	//	if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second]->second->polymer_cluster_pointer == second_polymer_unit_iter->polymer_cluster_pointer &&
	//		polymer_unit_iter != first_polymer_unit_iter)
	//	{
	//		mark = 1;//mark==1, it will still be the same bundle
	//		break;
	//	}
	//}
	if (mark == 1)
	{//if mark==1, then it is still belong to the same bundle, don't have to go through this complex case (but if mark!=1, it could still be the same bundle)
		debundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
		debundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
		rule_structure.bundling_polymer_unit_pair_list.splice(rule_structure.bundling_polymer_unit_pair_list.begin(), rule_structure.debundling_polymer_unit_pair_list, debundle_polymer_unit_pair_iter);
		debundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
		debundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = debundle_polymer_unit_pair_iter;
		debundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = debundle_polymer_unit_pair_iter;
		rule_parameters.debundling_inverse_sum--;
		rule_parameters.bundling_inverse_sum++;
	}
	else
	{
		list<Polymer_Bundle>::iterator old_second_bundle_iter = second_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer;



		debundle_polymer_unit_pair_iter->first->debundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.debundling_polymer_unit_pair_list.end();
		debundle_polymer_unit_pair_iter->second->debundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.debundling_polymer_unit_pair_list.end();
		rule_structure.bundling_polymer_unit_pair_list.splice(rule_structure.bundling_polymer_unit_pair_list.begin(), rule_structure.debundling_polymer_unit_pair_list, debundle_polymer_unit_pair_iter);
		debundle_polymer_unit_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
		debundle_polymer_unit_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = debundle_polymer_unit_pair_iter;
		debundle_polymer_unit_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = debundle_polymer_unit_pair_iter;
		rule_parameters.debundling_inverse_sum--;
		rule_parameters.bundling_inverse_sum++;

		Polymer_Bundle polymer_bundle;
		rule_structure.polymer_bundle_list.push_front(polymer_bundle);
		list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		recursion_polymer_unit(first_polymer_unit_iter, polymer_bundle_iter);
		if (first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer != second_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer)
		{//if they do not belong to the same bundle
			recursion_polymer_unit(second_polymer_unit_iter, old_second_bundle_iter);
			//cout << "1" << endl;;
			if (old_second_bundle_iter->polymer_anchor_list.size() == 0)
			{
				//cout <<"2" << endl;
				delete_polymer_bundle(old_second_bundle_iter);//why run this then go to the else section?

			}
			else
			{
				//cout <<"3"<< endl;
				//check_polymer_hydrolysis_sum((*old_second_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.begin());
				//check_polymer_bundle_diffusable_direction(second_polymer_unit_iter);
				//check_polymer_hydrolysis_sum(second_polymer_unit_iter);
				check_polymer_bundle_diffusable_direction(old_second_bundle_iter, 2);
				/*if (old_second_bundle_iter->bundle_cluster_pointer_list.size() == 1)
				{
				if ((*old_second_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.size() == 1)
				{
				rule_parameters.rotation_sum++;
				}
				}*/
			}
			if (first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
			{
				delete_polymer_bundle(first_polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			}
			else
			{
				check_polymer_hydrolysis_sum(first_polymer_unit_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
				/*if (polymer_bundle_iter->bundle_cluster_pointer_list.size() == 1)
				{
				if ((*polymer_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.size() == 1)
				{
				rule_parameters.rotation_sum++;
				}
				}*/
			}
		}
		else
		{//if they still belong to the same bundle
			old_second_bundle_iter->bundle_cluster_pointer_list.clear();
			old_second_bundle_iter->polymer_anchor_list.clear();

			for (int i = 0; i < 4; i++)
			{
				polymer_bundle_iter->availabe_diffusion_direction[i] = old_second_bundle_iter->availabe_diffusion_direction[i];

			}
			polymer_bundle_iter->available_polymer_diffusion_direction_sum = old_second_bundle_iter->available_polymer_diffusion_direction_sum;

			rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
			delete_polymer_bundle(old_second_bundle_iter);
		}





		//cout << "debundling finished" << endl;
		//rebuild_TT_Fragmentation_sum();
		//i used to have it because there is some unkonw error here, now i need to debug it
	}










};