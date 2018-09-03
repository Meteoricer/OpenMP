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
#include <map>
#include <time.h>
//#include <string>
#include <sstream>
#include <omp.h>

using namespace std;

extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
class Polymer_Unit;
class Anchor_Unit;




list<Polymer_Unit>::iterator operator+(list<Polymer_Unit>::iterator input_iter, int input)
{
	list<Polymer_Unit>::iterator temp;
	temp = input_iter;
	for (int i = 0; i < input; i++)
	{
		temp++;
	}
	return temp;
}

list<Polymer_Unit>::iterator operator-(list<Polymer_Unit>::iterator input_iter, int input)
{
	list<Polymer_Unit>::iterator temp;
	temp = input_iter;
	for (int i = 0; i < input; i++)
	{
		temp--;
	}
	return temp;
}

//int::operator++(int a, int b)
//{
//
//}






void rebuild_polymer_bc(list<Polymer_Cluster>::iterator input_polymer_cluster_iter)//bc:boundary condition
{

	
};

void recursion_polymer_unit(list<Polymer_Unit>::iterator polymer_unit_iter, list<Polymer_Bundle>::iterator polymer_bundle_iter)
{
	//int directions;
	//first delete every anchor
	//i might also need to delete every polymer_cluster_iter with in it?
	if (polymer_bundle_iter->Consider_mark != -rule_parameters.t)
	{
		
		if (polymer_bundle_iter->polymer_anchor_list.size() != 0)
		{
			polymer_bundle_iter->polymer_anchor_list.erase(polymer_bundle_iter->polymer_anchor_list.begin(), polymer_bundle_iter->polymer_anchor_list.end());
		}
		if (polymer_bundle_iter->bundle_cluster_pointer_list.size() != 0)
		{
			polymer_bundle_iter->bundle_cluster_pointer_list.erase(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter->bundle_cluster_pointer_list.end());
		}
		
		polymer_bundle_iter->Consider_mark = -rule_parameters.t;
		polymer_bundle_iter->bundle_direction = polymer_unit_iter->polymer_cluster_pointer->direction;
	}
	
	if (polymer_unit_iter->Consider_mark != rule_parameters.t)
	{
		
		polymer_unit_iter->Consider_mark = rule_parameters.t;
		//I might need to delete this from the bundle_cluster pointer
		//directions = polymer_unit_iter->polymer_cluster_pointer->direction;
		//polymer_cluster_pointer to polymer bundle
		polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer = polymer_bundle_iter;
		
		//if (polymer_bundle_iter->bundle_cluster_pointer_list.size() == 0 || polymer_unit_iter->polymer_cluster_pointer != *(polymer_bundle_iter->bundle_cluster_pointer_list.begin()))
		if(polymer_unit_iter->polymer_cluster_pointer->Consider_mark!=rule_parameters.t)
		{//if there is the cluster has not been consider before, add it into the bundle list
			polymer_unit_iter->polymer_cluster_pointer->Consider_mark = rule_parameters.t;
			polymer_unit_iter->polymer_cluster_pointer->anchor_pointer_list.clear();
			polymer_bundle_iter->bundle_cluster_pointer_list.push_front(polymer_unit_iter->polymer_cluster_pointer);
			polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer = polymer_bundle_iter;
			polymer_unit_iter->polymer_cluster_pointer->direction = polymer_bundle_iter->bundle_direction;
		}
		//check anchor here;
		if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
		{
			//there are anchor;
			polymer_bundle_iter->polymer_anchor_list.push_front(polymer_unit_iter->anchor_pointer);
			polymer_unit_iter->polymer_cluster_pointer->anchor_pointer_list.push_front(polymer_unit_iter->anchor_pointer);
			polymer_unit_iter->anchor_pointer->polymer_unit_pointer = polymer_unit_iter;
		}
		
		
		//there are four iterators for different direction
		//don't consider this for now.
		/*for (int i = direction_up; i < direction_right + 1; i++)
		{
			if (polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[i] != rule_structure.polymer_unit_list_end.begin())
			{
				recursion_polymer_unit(polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[i], polymer_bundle_iter);
			}
		}*/
		/*if (polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[0] != rule_structure.polymer_unit_list_end.begin())
		{
			recursion_polymer_unit(polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[0],polymer_bundle_iter);
		}
		if (polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[1] != rule_structure.polymer_unit_list_end.begin())
		{
			recursion_polymer_unit(polymer_unit_iter->bundled_neighbor_polymer_unit_pointer[1],polymer_bundle_iter);
		}*/
		if (polymer_unit_iter != (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
		{//along the sequence, make the iterator advance
			polymer_unit_iter++;
			recursion_polymer_unit(polymer_unit_iter,polymer_bundle_iter);
			polymer_unit_iter--;
		}
		
		if (polymer_unit_iter != polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
		{//along the sequence, make the iterator backward
			polymer_unit_iter--;
			recursion_polymer_unit(polymer_unit_iter,polymer_bundle_iter);
			polymer_unit_iter++;
		}
		

		//first do recursion within the polymer_sequence

		if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
		{//if there is a bundle on its first direction, left or down
		 //list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
		 //int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
		 //int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
		 //if ((col_distance*col_distance + row_distance*row_distance) != 1)
		 //{	//it is no longer in the annealing position,need to delete this annealing relation
		 //	bundle_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
		 //	bundle_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
		 //	rule_structure.bundling_polymer_unit_pair_list.erase(bundle_pair_iter);
			recursion_polymer_unit(polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first]->first, polymer_bundle_iter);


			//}
		}

		if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
		{//if there is a bundle on its second direction, right or up
		 //list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second];
		 //int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
		 //int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
		 //if ((col_distance*col_distance + row_distance*row_distance) != 1)
		 //{	//it is no longer in the annealing position,need to delete this annealing relation
		 //	bundle_pair_iter->first->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
		 //	bundle_pair_iter->second->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
		 //	rule_structure.bundling_polymer_unit_pair_list.erase(bundle_pair_iter);
			recursion_polymer_unit(polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second]->second, polymer_bundle_iter);


			//}
		}
	}
};

//void check_annealing(list<Polymer_Cluster>::iterator polymer_cluster_iter)
//{
//	list<Polymer_Unit>::iterator polymer_unit_iter;
//	if (polymer_cluster_iter->direction == vertical)
//	{//now there is only this case
//		//begin with bottom
//		polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
//		int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
//		int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
//		int up_col_position = col_position + 1;
//		if (up_col_position > column_num - 1)
//			up_col_position = 0;
//		int down_col_position = col_position - 1;
//		if (down_col_position < 0)
//			up_col_position = column_num - 1;
//		int right_row_position = row_position + 1;
//		if (right_row_position > row_num - 1)
//			right_row_position = 0;
//		int left_row_position = row_position - 1;
//		if (left_row_position < 0)
//			left_row_position = row_num - 1;
//		//then check those anneling stuff
//	}
//	
//
//}

void check_polymer_hydrolysis_sum(list<Polymer_Unit>::iterator polymer_unit_iter)
{
	check_polymer_hydrolysis_sum(polymer_unit_iter->polymer_cluster_pointer);
	///*if (&*polymer_unit_iter == &*rule_structure.polymer_unit_list_end.begin())
	//{

	//}*/
	////i might also need to consider those depolymerize hydrolysis thing here.
	//rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - polymer_unit_iter->polymer_cluster_pointer->TTT_sum;
	//rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - polymer_unit_iter->polymer_cluster_pointer->TTD_sum;
	//rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - polymer_unit_iter->polymer_cluster_pointer->DTD_sum;

	//rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - polymer_unit_iter->polymer_cluster_pointer->TT_sum;
	//rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - polymer_unit_iter->polymer_cluster_pointer->TD_sum;
	//rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - polymer_unit_iter->polymer_cluster_pointer->DD_sum;


	//rule_parameters.TT_depolymerize_end_sum = rule_parameters.TT_depolymerize_end_sum - polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum;
	//rule_parameters.TD_depolymerize_end_sum = rule_parameters.TD_depolymerize_end_sum - polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum;
	//rule_parameters.DD_depolymerize_end_sum = rule_parameters.DD_depolymerize_end_sum - polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum;
	////rule_parameters.TT_annealing_sum: i can write a function: check_annealing();
	////rule_parameters.TT_fragmentation_sum
	////rule_parameters.TT_depolymerize_end_sum
	//switch (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size())
	//{
	//	case 0:
	//	{
	//		std::cout << "error: Polymer_sequence.size()=0" << endl;
	//		system("pause");
	//		break;
	//	}
	//	case 1:
	//	{
	//		polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//		if (polymer_unit_iter->Hydrolysis_mark == T)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 1;
	//		}
	//		polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//		polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//		check_annealing(polymer_unit_iter->polymer_cluster_pointer);
	//		break;
	//	}
	//	case 2:
	//	{
	//		list<Polymer_Unit>::iterator first_unit_iter;
	//		//list<Polymer_Unit>::iterator center_unit_iter;
	//		list<Polymer_Unit>::iterator last_unit_iter;
	//		first_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//		last_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin()++;
	//		if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)//TT case
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 2;
	//			polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 1;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;

	//		}
	//		else
	//		{
	//			if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)//DD case
	//			{
	//				polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 1;
	//			}
	//			else//TD case
	//			{
	//				polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//			}
	//		}
	//		check_annealing(polymer_unit_iter->polymer_cluster_pointer);
	//		break;
	//	}

	//	case 3:
	//	{
	//		list<Polymer_Unit>::iterator first_unit_iter;
	//		list<Polymer_Unit>::iterator center_unit_iter;
	//		list<Polymer_Unit>::iterator last_unit_iter;
	//		first_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//		center_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin()+1;
	//		last_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin()+2;
	//		if (first_unit_iter->Hydrolysis_mark == T && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
	//		{//TTT case
	//			polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 3;
	//			polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 2;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;

	//		}
	//		
	//		else
	//		{
	//			if ((first_unit_iter->Hydrolysis_mark == D && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
	//				&& (first_unit_iter->Hydrolysis_mark == T && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D))
	//			{//DTT case or TTD case
	//				polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//			}


	//			else
	//			{

	//				if (first_unit_iter->Hydrolysis_mark == D && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D)
	//				{//DTD case
	//					polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 1;
	//					polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 2;
	//					polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//				}

	//				else//DDD case
	//				{
	//					polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 1;
	//					polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//					polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 2;
	//				}
	//			}
	//			
	//		}


	//	
	//		
	//		check_annealing(polymer_unit_iter->polymer_cluster_pointer);
	//		break;
	//	}

	//	default:
	//	{
	//		list<Polymer_Unit>::iterator first_unit_iter;
	//		list<Polymer_Unit>::iterator center_unit_iter;
	//		list<Polymer_Unit>::iterator last_unit_iter;
	//		/*center_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//		last_unit_iter = center_unit_iter++;
	//		last_unit_iter++;
	//		first_unit_iter= polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();*/

	//		first_unit_iter = polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin();
	//		center_unit_iter = first_unit_iter + 1;
	//		last_unit_iter = first_unit_iter + 2;
	//		//need to first judge the first one and the second one.
	//		if (first_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
	//		{//TTX
	//			polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 1;
	//			polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 1;
	//			polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//			polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//		}
	//		else
	//		{
	//			if (first_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
	//			{//DDX
	//				polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 1;
	//			}
	//			else
	//			{//TDX or DTX
	//				polymer_unit_iter->polymer_cluster_pointer->TTT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TTD_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->DTD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum = 0;
	//				polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum = 1;
	//				polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum = 0;
	//				
	//			}
	//		}
	//		/*if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TT_sum = 1;
	//		}
	//		if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TD_sum = 1;
	//		}
	//		if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == T)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TD_sum = 1;
	//		}
	//		if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->DD_sum = 1;
	//		}*/
	//		do
	//		{
	//			if (center_unit_iter->Hydrolysis_mark == T)//only consider center T case for XXX test
	//			{
	//				if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
	//				{
	//					polymer_unit_iter->polymer_cluster_pointer->TTT_sum++;
	//					//polymer_unit_iter->polymer_cluster_pointer->TT_sum++;
	//				}
	//				else
	//				{
	//					if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
	//					{
	//						polymer_unit_iter->polymer_cluster_pointer->DTD_sum++;
	//						//polymer_unit_iter->polymer_cluster_pointer->TD_sum++;
	//					}
	//					else
	//					{
	//						polymer_unit_iter->polymer_cluster_pointer->TTD_sum++;
	//					}
	//				}	
	//			}

	//			if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
	//			{
	//				polymer_unit_iter->polymer_cluster_pointer->TT_sum++;
	//			}
	//			else
	//			{
	//				if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
	//				{
	//					polymer_unit_iter->polymer_cluster_pointer->DD_sum++;
	//				}
	//				else
	//				{
	//					polymer_unit_iter->polymer_cluster_pointer->TD_sum++;
	//				}
	//			}




	//			first_unit_iter++;
	//			center_unit_iter++;
	//			last_unit_iter++;
	//		} while (last_unit_iter != --polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end());
	//		//consider the last one

	//		if (last_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
	//		{
	//			polymer_unit_iter->polymer_cluster_pointer->TTT_sum ++;
	//			//polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum ++;
	//		
	//	
	//		}
	//		else
	//		{
	//			if (last_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
	//			{
	//			
	//				
	//				//polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum ++;
	//			}
	//			else
	//			{
	//		
	//				
	//				//polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum ++;
	//			
	//			
	//			
	//			}
	//		}
	//		check_annealing(polymer_unit_iter->polymer_cluster_pointer);
	//		break;

	//	}


	//}
	//	
	//rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum + polymer_unit_iter->polymer_cluster_pointer->TT_sum;
	//rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum + polymer_unit_iter->polymer_cluster_pointer->TD_sum;
	//rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum + polymer_unit_iter->polymer_cluster_pointer->DD_sum;


	//rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total + polymer_unit_iter->polymer_cluster_pointer->TTT_sum;
	//rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total + polymer_unit_iter->polymer_cluster_pointer->TTD_sum;
	//rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total + polymer_unit_iter->polymer_cluster_pointer->DTD_sum;
	//rule_parameters.TT_depolymerize_end_sum = rule_parameters.TT_depolymerize_end_sum + polymer_unit_iter->polymer_cluster_pointer->TT_depolymerize_end_sum;
	//rule_parameters.TD_depolymerize_end_sum = rule_parameters.TD_depolymerize_end_sum + polymer_unit_iter->polymer_cluster_pointer->TD_depolymerize_end_sum;
	//rule_parameters.DD_depolymerize_end_sum = rule_parameters.DD_depolymerize_end_sum + polymer_unit_iter->polymer_cluster_pointer->DD_depolymerize_end_sum;

};



void check_polymer_hydrolysis_sum(list<Polymer_Cluster>::iterator polymer_cluster_iter)
{
	if (rule_parameters.TTT_hydrolysys_rate == 0)
	{//no hydrolysis
		if (polymer_cluster_iter->ring_mark == 1)
		{
			cout << "no hydrolysis case, should not happend" << endl;
			system("pause");
			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			polymer_cluster_iter->TTD_sum = 0;
			polymer_cluster_iter->DTD_sum = 0;
			polymer_cluster_iter->TTT_sum = 0;
			polymer_cluster_iter->TT_sum = 0;
			polymer_cluster_iter->TD_sum = 0;
			polymer_cluster_iter->DD_sum = 0;
			polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;


			polymer_cluster_iter->TTT_sum = polymer_cluster_iter->polymer_sequence.size();
			polymer_cluster_iter->TT_sum = polymer_cluster_iter->polymer_sequence.size();
			//polymer_cluster_iter->TT_depolymerize_end_sum = 2;


			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum + polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum + polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum + polymer_cluster_iter->DD_sum;


			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total + polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total + polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total + polymer_cluster_iter->DTD_sum;
			//rule_parameters.TT_depolymerize_end_sum = rule_parameters.TT_depolymerize_end_sum + polymer_cluster_iter->TT_depolymerize_end_sum;
			//rule_parameters.TD_depolymerize_end_sum = rule_parameters.TD_depolymerize_end_sum + polymer_cluster_iter->TD_depolymerize_end_sum;
			//rule_parameters.DD_depolymerize_end_sum = rule_parameters.DD_depolymerize_end_sum + polymer_cluster_iter->DD_depolymerize_end_sum;
		}
		else
		{// not ring case
			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			polymer_cluster_iter->TTD_sum = 0;
			polymer_cluster_iter->DTD_sum = 0;
			polymer_cluster_iter->TTT_sum = 0;
			polymer_cluster_iter->TT_sum = 0;
			polymer_cluster_iter->TD_sum = 0;
			polymer_cluster_iter->DD_sum = 0;
			polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;
			int length = polymer_cluster_iter->polymer_sequence.size();
			switch (length)
			{
				case 0:
				{

					break;
				}
				case 1:
				{
					polymer_cluster_iter->TTT_sum = 1;
					break;
				}
				case 2:
				{
					polymer_cluster_iter->TTT_sum = 2;
					polymer_cluster_iter->TT_sum = 0;
					polymer_cluster_iter->TT_depolymerize_plus_end_sum = 1;
					polymer_cluster_iter->TT_depolymerize_minus_end_sum = 1;
					break;
				}
				case 3:
				{
					polymer_cluster_iter->TTT_sum = polymer_cluster_iter->polymer_sequence.size();
					polymer_cluster_iter->TT_sum = 0;
					//polymer_cluster_iter->TT_depolymerize_end_sum = 2;
					break;
				}
				default:
				{
					polymer_cluster_iter->TTT_sum = polymer_cluster_iter->polymer_sequence.size();
					polymer_cluster_iter->TT_sum = polymer_cluster_iter->polymer_sequence.size() - 3;
					//polymer_cluster_iter->TT_depolymerize_end_sum = 2;
					break;

				}
			}

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum + polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum + polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum + polymer_cluster_iter->DD_sum;


			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total + polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total + polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total + polymer_cluster_iter->DTD_sum;
			//rule_parameters.TT_depolymerize_end_sum = rule_parameters.TT_depolymerize_end_sum + polymer_cluster_iter->TT_depolymerize_end_sum;
			//rule_parameters.TD_depolymerize_end_sum = rule_parameters.TD_depolymerize_end_sum + polymer_cluster_iter->TD_depolymerize_end_sum;
			//rule_parameters.DD_depolymerize_end_sum = rule_parameters.DD_depolymerize_end_sum + polymer_cluster_iter->DD_depolymerize_end_sum;

			
		}

	}
	else
	{
		if (polymer_cluster_iter->ring_mark == 1)
		{
			/*cout << "ring case check polymer hydrolysis sum, it is not finished" << endl;
			system("pause");*/
			//has nothing to do with polymerize


			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			polymer_cluster_iter->TTD_sum = 0;
			polymer_cluster_iter->DTD_sum = 0;
			polymer_cluster_iter->TTT_sum = 0;
			polymer_cluster_iter->TT_sum = 0;
			polymer_cluster_iter->TD_sum = 0;
			polymer_cluster_iter->DD_sum = 0;
			polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;

			
				list<Polymer_Unit>::iterator first_unit_iter;
				list<Polymer_Unit>::iterator center_unit_iter;
				list<Polymer_Unit>::iterator last_unit_iter;
				/*center_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				last_unit_iter = center_unit_iter++;
				last_unit_iter++;
				first_unit_iter= polymer_cluster_iter->polymer_sequence.begin();*/

				first_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				center_unit_iter = first_unit_iter + 1;
				last_unit_iter = first_unit_iter + 2;
				int triple_sum = 0;
				triple_sum = (--polymer_cluster_iter->polymer_sequence.end()->Hydrolysis_mark) * 100 + first_unit_iter->Hydrolysis_mark * 10 + center_unit_iter->Hydrolysis_mark;
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum++;

					polymer_cluster_iter->TT_sum++;
				}
				break;
				case TTD:
				{
					polymer_cluster_iter->TTD_sum++;

					polymer_cluster_iter->TD_sum++;
				}
				break;
				case TDT:
				{
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DTT:
				{
					polymer_cluster_iter->TTD_sum++;

					polymer_cluster_iter->TT_sum++;
				}
				break;
				case TDD:
				{
					polymer_cluster_iter->DD_sum++;
				}
				break;
				case DTD:
				{
					polymer_cluster_iter->DTD_sum++;
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDT:
				{
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDD:
				{
					polymer_cluster_iter->DD_sum++;
				}
				break;
				}
				while (last_unit_iter != polymer_cluster_iter->polymer_sequence.end())
				{
					triple_sum = first_unit_iter->Hydrolysis_mark * 100 + center_unit_iter->Hydrolysis_mark * 10 + last_unit_iter->Hydrolysis_mark;
					switch (triple_sum)
					{
					case TTT:
					{
						polymer_cluster_iter->TTT_sum++;

						polymer_cluster_iter->TT_sum++;
					}
					break;
					case TTD:
					{
						polymer_cluster_iter->TTD_sum++;

						polymer_cluster_iter->TD_sum++;
					}
					break;
					case TDT:
					{
						polymer_cluster_iter->TD_sum++;
					}
					break;
					case DTT:
					{
						polymer_cluster_iter->TTD_sum++;

						polymer_cluster_iter->TT_sum++;
					}
					break;
					case TDD:
					{
						polymer_cluster_iter->DD_sum++;
					}
					break;
					case DTD:
					{
						polymer_cluster_iter->DTD_sum++;
						polymer_cluster_iter->TD_sum++;
					}
					break;
					case DDT:
					{
						polymer_cluster_iter->TD_sum++;
					}
					break;
					case DDD:
					{
						polymer_cluster_iter->DD_sum++;
					}
					break;
					}
					first_unit_iter++;
					center_unit_iter++;
					last_unit_iter++;
				}
				
				triple_sum = first_unit_iter->Hydrolysis_mark * 100 + center_unit_iter->Hydrolysis_mark * 10 + polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark;
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum++;

					polymer_cluster_iter->TT_sum++;
				}
				break;
				case TTD:
				{
					polymer_cluster_iter->TTD_sum++;

					polymer_cluster_iter->TD_sum++;
				}
				break;
				case TDT:
				{
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DTT:
				{
					polymer_cluster_iter->TTD_sum++;

					polymer_cluster_iter->TT_sum++;
				}
				break;
				case TDD:
				{
					polymer_cluster_iter->DD_sum++;
				}
				break;
				case DTD:
				{
					polymer_cluster_iter->DTD_sum++;
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDT:
				{
					polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDD:
				{
					polymer_cluster_iter->DD_sum++;
				}
				break;
				}
				////need to first judge the first one and the second one.
				//if (first_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
				//{//TTX
				//	polymer_cluster_iter->TTT_sum = 1;
				//	polymer_cluster_iter->TTD_sum = 0;
				//	polymer_cluster_iter->DTD_sum = 0;
				//	polymer_cluster_iter->TT_sum = 1;
				//	polymer_cluster_iter->TD_sum = 0;
				//	polymer_cluster_iter->DD_sum = 0;
				//	polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//	polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//	polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//	/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum--;
				//	}*/
				//}
				//else
				//{
				//	if (first_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
				//	{//DDX
				//		polymer_cluster_iter->TTT_sum = 0;
				//		polymer_cluster_iter->TTD_sum = 0;
				//		polymer_cluster_iter->DTD_sum = 0;
				//		polymer_cluster_iter->TT_sum = 0;
				//		polymer_cluster_iter->TD_sum = 0;
				//		polymer_cluster_iter->DD_sum = 1;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//		/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->DD_depolymerize_end_sum--;
				//		}*/
				//	}
				//	else
				//	{//TDX or DTX
				//		polymer_cluster_iter->TTT_sum = 0;
				//		polymer_cluster_iter->TTD_sum = 1;
				//		polymer_cluster_iter->DTD_sum = 0;
				//		polymer_cluster_iter->TT_sum = 0;
				//		polymer_cluster_iter->TD_sum = 1;
				//		polymer_cluster_iter->DD_sum = 0;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//		/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum--;
				//		}*/

				//	}
				//}
				///*if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//{
				//polymer_cluster_iter->TT_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D)
				//{
				//polymer_cluster_iter->TD_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == T)
				//{
				//polymer_cluster_iter->TD_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//{
				//polymer_cluster_iter->DD_sum = 1;
				//}*/
				//do
				//{
				//	if (center_unit_iter->Hydrolysis_mark == T)//only consider center T case for XXX test
				//	{
				//		if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//		{
				//			polymer_cluster_iter->TTT_sum++;
				//			//polymer_cluster_iter->TT_sum++;
				//		}
				//		else
				//		{
				//			if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//			{
				//				polymer_cluster_iter->DTD_sum++;
				//				//polymer_cluster_iter->TD_sum++;
				//			}
				//			else
				//			{
				//				polymer_cluster_iter->TTD_sum++;
				//			}
				//		}
				//	}

				//	if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//	{
				//		polymer_cluster_iter->TT_sum++;
				//	}
				//	else
				//	{
				//		if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//		{
				//			polymer_cluster_iter->DD_sum++;
				//		}
				//		else
				//		{
				//			polymer_cluster_iter->TD_sum++;
				//		}
				//	}




				//	first_unit_iter++;
				//	center_unit_iter++;
				//	last_unit_iter++;
				//} while (last_unit_iter != --polymer_cluster_iter->polymer_sequence.end());
				////consider the last one

				//last_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				//center_unit_iter++;
				//first_unit_iter++;



				//if (center_unit_iter->Hydrolysis_mark == T)//only consider center T case for XXX test
				//{
				//	if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//	{
				//		polymer_cluster_iter->TTT_sum++;
				//		//polymer_cluster_iter->TT_sum++;
				//	}
				//	else
				//	{
				//		if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//		{
				//			polymer_cluster_iter->DTD_sum++;
				//			//polymer_cluster_iter->TD_sum++;
				//		}
				//		else
				//		{
				//			polymer_cluster_iter->TTD_sum++;
				//		}
				//	}
				//}

				//if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//{
				//	polymer_cluster_iter->TT_sum++;
				//}
				//else
				//{
				//	if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//	{
				//		polymer_cluster_iter->DD_sum++;
				//	}
				//	else
				//	{
				//		polymer_cluster_iter->TD_sum++;
				//	}
				//}


				//if (last_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
				//{
				//	
				//	polymer_cluster_iter->TTT_sum++;
				//	//polymer_cluster_iter->TT_depolymerize_end_sum++;
				//	polymer_cluster_iter->TT_sum++;
				//	/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum--;
				//	}*/


				//}
				//else
				//{
				//	if (last_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
				//	{


				//		//polymer_cluster_iter->DD_depolymerize_end_sum++;
				//		polymer_cluster_iter->DD_sum++;
				//		/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->DD_depolymerize_end_sum--;
				//		}*/
				//	}
				//	else
				//	{


				//		//polymer_cluster_iter->TD_depolymerize_end_sum++;
				//		polymer_cluster_iter->TD_sum++;
				//		/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum--;
				//		}*/



				//	}
				//}
				//check_annealing(polymer_cluster_iter);


			

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum + polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum + polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum + polymer_cluster_iter->DD_sum;


			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total + polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total + polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total + polymer_cluster_iter->DTD_sum;
			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum + polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum + polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum + polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum + polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum + polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum + polymer_cluster_iter->DD_depolymerize_minus_end_sum;

		}
		else
		{
			//i might also need to consider those depolymerize hydrolysis thing here.
			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - polymer_cluster_iter->DD_sum;

			
			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			//cout << rule_parameters.TT_depolymerize_minus_end_sum << endl;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - polymer_cluster_iter->DD_depolymerize_minus_end_sum;
			//rule_parameters.TT_annealing_sum: i can write a function: check_annealing();
			//rule_parameters.TT_fragmentation_sum
			//rule_parameters.TT_depolymerize_end_sum

			polymer_cluster_iter->TTD_sum = 0;
			polymer_cluster_iter->DTD_sum = 0;
			polymer_cluster_iter->TTT_sum = 0;
			polymer_cluster_iter->TT_sum = 0;
			polymer_cluster_iter->TD_sum = 0;
			polymer_cluster_iter->DD_sum = 0;
			polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;
			if (polymer_cluster_iter->direction == vertical)
			{
				rule_parameters.TT_polymerize_end_vertical_plus_sum -= polymer_cluster_iter->TT_polymerize_plus_end_sum;
				rule_parameters.TT_polymerize_end_vertical_minus_sum -= polymer_cluster_iter->TT_polymerize_minus_end_sum;
				rule_parameters.TD_polymerize_end_vertical_plus_sum -= polymer_cluster_iter->TD_polymerize_plus_end_sum;
				rule_parameters.TD_polymerize_end_vertical_minus_sum -= polymer_cluster_iter->TD_polymerize_minus_end_sum;
				polymer_cluster_iter->TT_polymerize_plus_end_sum = 0;
				polymer_cluster_iter->TT_polymerize_minus_end_sum = 0;
				polymer_cluster_iter->TD_polymerize_plus_end_sum = 0;
				polymer_cluster_iter->TD_polymerize_minus_end_sum = 0;


				if (polymer_cluster_iter->polarize == direction_first)
				{//down or left//saves as left to right
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_down] == true)
					{
						if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum++;
						}
					}

					if (polymer_cluster_iter->availabe_polymerize_direction[direction_up] == true)
					{
						if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum++;
						}

					}
				}
				else
				{
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_up] == true)
					{
						if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum++;
						}
					}
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_down] == true)
					{
						if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum++;
						}

					}

				}
				rule_parameters.TT_polymerize_end_vertical_plus_sum += polymer_cluster_iter->TT_polymerize_plus_end_sum;
				rule_parameters.TT_polymerize_end_vertical_minus_sum += polymer_cluster_iter->TT_polymerize_minus_end_sum;
				rule_parameters.TD_polymerize_end_vertical_plus_sum += polymer_cluster_iter->TD_polymerize_plus_end_sum;
				rule_parameters.TD_polymerize_end_vertical_minus_sum += polymer_cluster_iter->TD_polymerize_minus_end_sum;
			}
			else
			{
				rule_parameters.TT_polymerize_end_horizontal_plus_sum -= polymer_cluster_iter->TT_polymerize_plus_end_sum;
				rule_parameters.TT_polymerize_end_horizontal_minus_sum -= polymer_cluster_iter->TT_polymerize_minus_end_sum;
				rule_parameters.TD_polymerize_end_horizontal_plus_sum -= polymer_cluster_iter->TD_polymerize_plus_end_sum;
				rule_parameters.TD_polymerize_end_horizontal_minus_sum -= polymer_cluster_iter->TD_polymerize_minus_end_sum;
				polymer_cluster_iter->TT_polymerize_plus_end_sum = 0;
				polymer_cluster_iter->TT_polymerize_minus_end_sum = 0;
				polymer_cluster_iter->TD_polymerize_plus_end_sum = 0;
				polymer_cluster_iter->TD_polymerize_minus_end_sum = 0;


				if (polymer_cluster_iter->polarize == direction_first)
				{//down or left//saves as left to right
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_left] == true)
					{
						if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum++;
						}
					}

					if (polymer_cluster_iter->availabe_polymerize_direction[direction_right] == true)
					{
						if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum++;
						}

					}
				}
				else
				{
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_right] == true)
					{
						if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum++;
						}
					}
					if (polymer_cluster_iter->availabe_polymerize_direction[direction_left] == true)
					{
						if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum++;
						}
						else
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum++;
						}

					}

				}
				rule_parameters.TT_polymerize_end_horizontal_plus_sum += polymer_cluster_iter->TT_polymerize_plus_end_sum;
				rule_parameters.TT_polymerize_end_horizontal_minus_sum += polymer_cluster_iter->TT_polymerize_minus_end_sum;
				rule_parameters.TD_polymerize_end_horizontal_plus_sum += polymer_cluster_iter->TD_polymerize_plus_end_sum;
				rule_parameters.TD_polymerize_end_horizontal_minus_sum += polymer_cluster_iter->TD_polymerize_minus_end_sum;
			}
			


			
			
				






			switch (polymer_cluster_iter->polymer_sequence.size())
			{
			case 0:
			{

				break;
			}
			
			case 1:
			{
				//polymer_cluster_iter->TTD_sum = 0;
				//polymer_cluster_iter->DTD_sum = 0;
				if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
				{
					polymer_cluster_iter->TTT_sum = 1;
				}

				//polymer_cluster_iter->TT_sum = 0;
				//polymer_cluster_iter->TD_sum = 0;
				//polymer_cluster_iter->DD_sum = 0;
				//polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//check_annealing(polymer_cluster_iter);
				break;
			}
			
			case 2:
			{
				list<Polymer_Unit>::iterator first_unit_iter;
				//list<Polymer_Unit>::iterator center_unit_iter;
				list<Polymer_Unit>::iterator last_unit_iter;
				first_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				last_unit_iter = (++polymer_cluster_iter->polymer_sequence.begin());
				if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)//TT case
				{
					polymer_cluster_iter->TTT_sum = 2;
					//polymer_cluster_iter->TTD_sum = 0;
					//polymer_cluster_iter->DTD_sum = 0;
					//polymer_cluster_iter->TT_sum = 0;
					//polymer_cluster_iter->TD_sum = 0;
					//polymer_cluster_iter->DD_sum = 0;
					polymer_cluster_iter->TT_depolymerize_plus_end_sum = 1;
					polymer_cluster_iter->TT_depolymerize_minus_end_sum = 1;
					//polymer_cluster_iter->TD_depolymerize_end_sum = 0;
					//polymer_cluster_iter->DD_depolymerize_end_sum = 0;
					//if (polymer_cluster_iter->anchor_pointer_list.size() < 2||
					//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first]!=rule_structure.debundling_polymer_unit_pair_list.end() ||
					//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end() ||
					//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
					//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())

					//{//if there is only one anchor, then depolymerize can happen, otherwise it should be deanchoring
					//	polymer_cluster_iter->TT_depolymerize_end_sum = 1;
					//}
					//else
					//{
					//	polymer_cluster_iter->TT_depolymerize_end_sum = 0;
					//}



				}
				else
				{
					if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)//DD case
					{
						//polymer_cluster_iter->TTT_sum = 0;
						//polymer_cluster_iter->TTD_sum = 0;
						//polymer_cluster_iter->DTD_sum = 0;
						//polymer_cluster_iter->TT_sum = 0;
						//polymer_cluster_iter->TD_sum = 0;
						//polymer_cluster_iter->DD_sum = 0;
						//polymer_cluster_iter->TT_depolymerize_end_sum = 0;
						//polymer_cluster_iter->TD_depolymerize_end_sum = 0;
						polymer_cluster_iter->DD_depolymerize_plus_end_sum = 1;
						polymer_cluster_iter->DD_depolymerize_minus_end_sum = 1;
						//if (polymer_cluster_iter->anchor_pointer_list.size() < 2 ||
						//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						//{//if there is only one anchor, then depolymerize can happen, otherwise it should be deanchoring
						//	polymer_cluster_iter->DD_depolymerize_end_sum = 1;
						//}
						//else
						//{
						//	polymer_cluster_iter->DD_depolymerize_end_sum = 0;
						//}
					}
					else//TD case
					{
						//polymer_cluster_iter->TTT_sum = 0;
						polymer_cluster_iter->TTD_sum = 1;
						//polymer_cluster_iter->DTD_sum = 0;
						//polymer_cluster_iter->TT_sum = 0;
						//polymer_cluster_iter->TD_sum = 0;
						//polymer_cluster_iter->DD_sum = 0;

						//polymer_cluster_iter->TT_depolymerize_end_sum = 0;
						polymer_cluster_iter->TD_depolymerize_plus_end_sum = 1;
						polymer_cluster_iter->TD_depolymerize_minus_end_sum = 1;
						//polymer_cluster_iter->DD_depolymerize_end_sum = 0;
						//if (polymer_cluster_iter->anchor_pointer_list.size() < 2 ||
						//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
						//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						//{//if there is only one anchor, then depolymerize can happen, otherwise it should be deanchoring
						//	polymer_cluster_iter->TT_depolymerize_end_sum = 1;
						//}
						//else
						//{
						//	polymer_cluster_iter->TT_depolymerize_end_sum = 0;
						//}
					}
				}
				//check_annealing(polymer_cluster_iter);
				break;
			}
			
			case 3:
			{
				list<Polymer_Unit>::iterator first_unit_iter;
				list<Polymer_Unit>::iterator center_unit_iter;
				list<Polymer_Unit>::iterator last_unit_iter;
				first_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				center_unit_iter = polymer_cluster_iter->polymer_sequence.begin() + 1;
				last_unit_iter = polymer_cluster_iter->polymer_sequence.begin() + 2;
				int triple_sum = first_unit_iter->Hydrolysis_mark * 100 + center_unit_iter->Hydrolysis_mark*10 + last_unit_iter->Hydrolysis_mark;
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum=3;

					polymer_cluster_iter->TT_depolymerize_plus_end_sum=1;
					polymer_cluster_iter->TT_depolymerize_minus_end_sum = 1;
				}
				break;
				case TTD:
				{
					polymer_cluster_iter->TTT_sum++;
					polymer_cluster_iter->TTD_sum++;
					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->TT_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TT_depolymerize_minus_end_sum++;
					}
				
				}
				break;
				case TDT:
				{
					polymer_cluster_iter->TD_depolymerize_plus_end_sum = 1;
					polymer_cluster_iter->TD_depolymerize_minus_end_sum = 1;
				}
				break;
				case DTT:
				{
					//polymer_cluster_iter->TTT_sum++;
					polymer_cluster_iter->TTD_sum++;

					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->TT_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TT_depolymerize_minus_end_sum++;
					}
				}
					break;
				case TDD:
				{
					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->DD_depolymerize_minus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->DD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
				break;
				case DTD:
				{
					polymer_cluster_iter->DTD_sum++;
					polymer_cluster_iter->TD_depolymerize_plus_end_sum = 1;
					polymer_cluster_iter->TD_depolymerize_minus_end_sum = 1;
				}
				
					break;
				case DDT:
				{
					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->DD_depolymerize_minus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->DD_depolymerize_plus_end_sum++;
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
					break;
				case DDD:
				{
					polymer_cluster_iter->DD_depolymerize_plus_end_sum = 1;
					polymer_cluster_iter->DD_depolymerize_minus_end_sum = 1;
				}
					break;
				}

				//if (first_unit_iter->Hydrolysis_mark == T && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//{//TTT case
				//	polymer_cluster_iter->TTT_sum = 3;
				//	//polymer_cluster_iter->TTD_sum = 0;
				//	//polymer_cluster_iter->DTD_sum = 0;
				//	//polymer_cluster_iter->TT_sum = 0;
				//	//polymer_cluster_iter->TD_sum = 0;
				//	//polymer_cluster_iter->DD_sum = 0;
				//	polymer_cluster_iter->TT_depolymerize_end_sum = 2;
				//	//polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//	//polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//	/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end()||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())&&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end()||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum --;
				//	}
				//	if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum--;
				//	}*/

				//}

				//else
				//{
				//	if ((first_unit_iter->Hydrolysis_mark == D && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//		|| (first_unit_iter->Hydrolysis_mark == T && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D))
				//	{//DTT case or TTD case
				//		polymer_cluster_iter->TTT_sum = 1;
				//		polymer_cluster_iter->TTD_sum = 1;
				//		//polymer_cluster_iter->DTD_sum = 0;
				//		//polymer_cluster_iter->TT_sum = 0;
				//		//polymer_cluster_iter->TD_sum = 0;
				//		//polymer_cluster_iter->DD_sum = 0;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 1;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 1;
				//		//polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//		/*if (first_unit_iter->Hydrolysis_mark == D &&
				//		(first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum --;
				//		}
				//		if (first_unit_iter->Hydrolysis_mark == T &&
				//		(first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TT_depolymerize_end_sum--;
				//		}
				//		if (last_unit_iter->Hydrolysis_mark == D &&
				//		(last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum--;
				//		}
				//		if (last_unit_iter->Hydrolysis_mark == T &&
				//		(last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TT_depolymerize_end_sum--;
				//		}*/
				//	}


				//	else
				//	{

				//		if (first_unit_iter->Hydrolysis_mark == D && center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D)
				//		{//DTD case
				//			polymer_cluster_iter->TTT_sum = 0;
				//			polymer_cluster_iter->TTD_sum = 0;
				//			polymer_cluster_iter->DTD_sum = 1;
				//			polymer_cluster_iter->TT_sum = 0;
				//			polymer_cluster_iter->TD_sum = 0;
				//			polymer_cluster_iter->DD_sum = 0;
				//			polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//			polymer_cluster_iter->TD_depolymerize_end_sum = 2;
				//			polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//			/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//			(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//			{
				//			polymer_cluster_iter->TD_depolymerize_end_sum--;
				//			}
				//			if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//			(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//			{
				//			polymer_cluster_iter->TD_depolymerize_end_sum--;
				//			}*/
				//		}

				//		else//DDD case
				//		{
				//			polymer_cluster_iter->TTT_sum = 0;
				//			polymer_cluster_iter->TTD_sum = 0;
				//			polymer_cluster_iter->DTD_sum = 0;
				//			polymer_cluster_iter->TT_sum = 0;
				//			polymer_cluster_iter->TD_sum = 0;
				//			polymer_cluster_iter->DD_sum = 0;
				//			polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//			polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//			polymer_cluster_iter->DD_depolymerize_end_sum = 2;
				//			/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//			(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//			{
				//			polymer_cluster_iter->DD_depolymerize_end_sum--;
				//			}
				//			if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//			(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//			center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//			{
				//			polymer_cluster_iter->DD_depolymerize_end_sum--;
				//			}*/
				//		}
				//	}

				//}




				////check_annealing(polymer_cluster_iter);
				break;
			}

			default:
			{
				/*Polymer_Unit begin_one(1);
				Polymer_Unit end_one(1);
				polymer_cluster_iter->polymer_sequence.push_front(begin_one);
				polymer_cluster_iter->polymer_sequence.push_back(end_one);*/
				list<Polymer_Unit>::iterator first_unit_iter;
				list<Polymer_Unit>::iterator center_unit_iter;
				list<Polymer_Unit>::iterator last_unit_iter;
				//list<Polymer_Unit>::iterator test_unit_iter;
				/*center_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				last_unit_iter = center_unit_iter++;
				last_unit_iter++;
				first_unit_iter= polymer_cluster_iter->polymer_sequence.begin();*/

				first_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				center_unit_iter = polymer_cluster_iter->polymer_sequence.begin() + 1;
				last_unit_iter = polymer_cluster_iter->polymer_sequence.begin() + 2;
				//test_unit_iter = polymer_cluster_iter->polymer_sequence.begin() + 3;
				int triple_sum=first_unit_iter->Hydrolysis_mark*10+center_unit_iter->Hydrolysis_mark +T*100;//for XXX like case
				/*int TTT = T * 10 + T + T;
				int TTD = T * 10 + T + D;
				int TDT = D * 10 + T + T;
				int DTT = TTD;
				int TDD = D * 10 + T + D;
				int DTD = T * 10 + D + D;
				int DDT = TDD;
				int DDD = D * 10 + D + D;*/
				//int double_sum;//for XX like case
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum++;
					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->TT_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TT_depolymerize_minus_end_sum++;
					}
						
				}
					break;
				case TTD:
				{
					polymer_cluster_iter->TTD_sum++;

					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
					break;
				case TDT:
				{
					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
					break;
				case DTT:
					//not possible
					break;
				case TDD:
				{
					if (polymer_cluster_iter->polarize == direction_first)
					{
						polymer_cluster_iter->DD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->DD_depolymerize_minus_end_sum++;
					}
				}
					break;
				case DTD:
					//not possible
					break;
				case DDT:
					//not possible
					break;
				case DDD:
					//not possible
					break;
				}
				auto end_iter = (--polymer_cluster_iter->polymer_sequence.end());
				while (last_unit_iter != end_iter)
				{
					triple_sum = first_unit_iter->Hydrolysis_mark * 100 + center_unit_iter->Hydrolysis_mark * 10 + last_unit_iter->Hydrolysis_mark;
					//now i need to consider XXX the center for hydrolysis
					//                        XX the last two for fragmentation
					switch (triple_sum)
					{
					case TTT:
					{
						polymer_cluster_iter->TTT_sum++;

						polymer_cluster_iter->TT_sum++;
					}
					break;
					case TTD:
					{
						polymer_cluster_iter->TTD_sum++;

						polymer_cluster_iter->TD_sum++;
					}
					break;
					case TDT:
					{
						polymer_cluster_iter->TD_sum++;
					}
					break;
					case DTT:
					{
						polymer_cluster_iter->TTD_sum++;

						polymer_cluster_iter->TT_sum++;
					}
						break;
					case TDD:
					{
						polymer_cluster_iter->DD_sum++;
					}
					break;
					case DTD:
					{
						polymer_cluster_iter->DTD_sum++;
						polymer_cluster_iter->TD_sum++;
					}
						break;
					case DDT:
					{
						polymer_cluster_iter->TD_sum++;
					}
						break;
					case DDD:
					{
						polymer_cluster_iter->DD_sum++;
					}
						break;
					}
					first_unit_iter++;
					center_unit_iter++;
					last_unit_iter++;
					//test_unit_iter++;
				}
				//now it is at the end XXX
				//                     XXX-END
				triple_sum = first_unit_iter->Hydrolysis_mark * 100 + center_unit_iter->Hydrolysis_mark * 10 + last_unit_iter->Hydrolysis_mark;
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum++;

					//polymer_cluster_iter->TT_sum++;
				}
				break;
				case TTD:
				{
					polymer_cluster_iter->TTD_sum++;

					//polymer_cluster_iter->TD_sum++;
				}
				break;
				case TDT:
				{
					//polymer_cluster_iter->TD_sum++;
				}
				break;
				case DTT:
				{
					polymer_cluster_iter->TTD_sum++;

					//polymer_cluster_iter->TD_sum++;
				}
				break;
				case TDD:
				{
					//polymer_cluster_iter->DD_sum++;
				}
				break;
				case DTD:
				{
					polymer_cluster_iter->DTD_sum++;
					//polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDT:
				{
					//polymer_cluster_iter->TD_sum++;
				}
				break;
				case DDD:
				{
					//polymer_cluster_iter->DD_sum++;
				}
				break;
				}


















				triple_sum = center_unit_iter->Hydrolysis_mark * 100 + last_unit_iter->Hydrolysis_mark * 10 + T;
				switch (triple_sum)
				{
				case TTT:
				{
					polymer_cluster_iter->TTT_sum++;

					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->TT_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TT_depolymerize_minus_end_sum++;
					}
				}
				break;
				case TTD:
				//not possible
				break;
				case TDT:
				{
					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
				break;
				case DTT:
				{
					polymer_cluster_iter->TTD_sum++;

					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->TD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->TD_depolymerize_minus_end_sum++;
					}
				}
				break;
				case TDD:
					//not possible
				break;
				case DTD:
					//not possible
				break;
				case DDT:
				{
					if (polymer_cluster_iter->polarize == direction_second)
					{
						polymer_cluster_iter->DD_depolymerize_plus_end_sum++;
					}
					else
					{
						polymer_cluster_iter->DD_depolymerize_minus_end_sum++;
					}
				}
				break;
				case DDD:
					//not possible
				break;
				}
			/*	if (triple_sum == TTT)
				{
					polymer_cluster_iter->TTT_sum++;
					
					polymer_cluster_iter->TT_depolymerize_end_sum ++;
					
				}
				else if (triple_sum == TTD)
				{
					
					polymer_cluster_iter->TTD_sum ++;
					
					polymer_cluster_iter->TD_depolymerize_end_sum++;
					
				}*/
				//else if (triple_sum == DTD)
				//{//not possible since the first one must be T

				//}





				//list<Polymer_Unit>::iterator first_unit_iter;
				//list<Polymer_Unit>::iterator center_unit_iter;
				//list<Polymer_Unit>::iterator last_unit_iter;
				///*center_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				//last_unit_iter = center_unit_iter++;
				//last_unit_iter++;
				//first_unit_iter= polymer_cluster_iter->polymer_sequence.begin();*/

				//first_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				//center_unit_iter = first_unit_iter + 1;
				//last_unit_iter = first_unit_iter + 2;
				////need to first judge the first one and the second one.
				//if (first_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
				//{//TTX
				//	polymer_cluster_iter->TTT_sum = 1;
				//	polymer_cluster_iter->TTD_sum = 0;
				//	polymer_cluster_iter->DTD_sum = 0;
				//	polymer_cluster_iter->TT_sum = 0;
				//	polymer_cluster_iter->TD_sum = 0;
				//	polymer_cluster_iter->DD_sum = 0;
				//	polymer_cluster_iter->TT_depolymerize_end_sum = 1;
				//	polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//	polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//	/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum--;
				//	}*/
				//}
				//else
				//{
				//	if (first_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
				//	{//DDX
				//		polymer_cluster_iter->TTT_sum = 0;
				//		polymer_cluster_iter->TTD_sum = 0;
				//		polymer_cluster_iter->DTD_sum = 0;
				//		polymer_cluster_iter->TT_sum = 0;
				//		polymer_cluster_iter->TD_sum = 0;
				//		polymer_cluster_iter->DD_sum = 0;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->DD_depolymerize_end_sum = 1;
				//		/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->DD_depolymerize_end_sum--;
				//		}*/
				//	}
				//	else if (first_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == D)//TDX
				//	{//TDX 
				//		polymer_cluster_iter->TTT_sum = 0;
				//		polymer_cluster_iter->TTD_sum = 1;
				//		polymer_cluster_iter->DTD_sum = 0;
				//		polymer_cluster_iter->TT_sum = 0;
				//		polymer_cluster_iter->TD_sum = 0;
				//		polymer_cluster_iter->DD_sum = 0;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 1;
				//		polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//		/*if ((first_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		first_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum--;
				//		}*/

				//	}
				//	else//or DTX
				//	{
				//		polymer_cluster_iter->TTT_sum = 0;
				//		polymer_cluster_iter->TTD_sum = 0;
				//		polymer_cluster_iter->DTD_sum = 0;
				//		polymer_cluster_iter->TT_sum = 0;
				//		polymer_cluster_iter->TD_sum = 0;
				//		polymer_cluster_iter->DD_sum = 0;
				//		polymer_cluster_iter->TT_depolymerize_end_sum = 0;
				//		polymer_cluster_iter->TD_depolymerize_end_sum = 1;
				//		polymer_cluster_iter->DD_depolymerize_end_sum = 0;
				//	}
				//}
				///*if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//{
				//polymer_cluster_iter->TT_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == D)
				//{
				//polymer_cluster_iter->TD_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == T)
				//{
				//polymer_cluster_iter->TD_sum = 1;
				//}
				//if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//{
				//polymer_cluster_iter->DD_sum = 1;
				//}*/
				//do
				//{
				//	if (center_unit_iter->Hydrolysis_mark == T)//only consider center T case for XXX test
				//	{
				//		if (first_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//		{
				//			polymer_cluster_iter->TTT_sum++;
				//			//polymer_cluster_iter->TT_sum++;
				//		}
				//		else
				//		{
				//			if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//			{
				//				polymer_cluster_iter->DTD_sum++;
				//				//polymer_cluster_iter->TD_sum++;
				//			}
				//			else //if (first_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//				 //TTD or DTT
				//			{
				//				polymer_cluster_iter->TTD_sum++;
				//			}
				//		}
				//	}

				//	if (center_unit_iter->Hydrolysis_mark == T && last_unit_iter->Hydrolysis_mark == T)
				//	{
				//		polymer_cluster_iter->TT_sum++;
				//	}
				//	else
				//	{
				//		if (center_unit_iter->Hydrolysis_mark == D && last_unit_iter->Hydrolysis_mark == D)
				//		{
				//			polymer_cluster_iter->DD_sum++;
				//		}
				//		else
				//		{
				//			polymer_cluster_iter->TD_sum++;
				//		}
				//	}




				//	first_unit_iter++;
				//	center_unit_iter++;
				//	last_unit_iter++;
				//} while (last_unit_iter != (--polymer_cluster_iter->polymer_sequence.end()));
				////consider the last one

				//if (last_unit_iter->Hydrolysis_mark == T &&center_unit_iter->Hydrolysis_mark == T)
				//{
				//	polymer_cluster_iter->TTT_sum++;
				//	polymer_cluster_iter->TT_depolymerize_end_sum++;
				//	/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//	(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//	center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//	{
				//	polymer_cluster_iter->TT_depolymerize_end_sum--;
				//	}*/


				//}
				//else
				//{
				//	if (last_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == D)
				//	{


				//		polymer_cluster_iter->DD_depolymerize_end_sum++;
				//		/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->DD_depolymerize_end_sum--;
				//		}*/
				//	}
				//	else if (last_unit_iter->Hydrolysis_mark == D &&center_unit_iter->Hydrolysis_mark == T)//DT_T
				//	{

				//		polymer_cluster_iter->TTD_sum++;
				//		polymer_cluster_iter->TD_depolymerize_end_sum++;
				//		/*if ((last_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		last_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()) &&
				//		(center_unit_iter->anchor_pointer != rule_structure.anchor_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end() ||
				//		center_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end()))
				//		{
				//		polymer_cluster_iter->TD_depolymerize_end_sum--;
				//		}*/



				//	}
				//	else//TD_T
				//	{
				//		polymer_cluster_iter->TD_depolymerize_end_sum++;
				//	}
				//}
				////check_annealing(polymer_cluster_iter);
				//break;

			}


			}

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum + polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum + polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum + polymer_cluster_iter->DD_sum;


			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total + polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total + polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total + polymer_cluster_iter->DTD_sum;
			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum + polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum + polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum + polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum + polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum + polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum + polymer_cluster_iter->DD_depolymerize_minus_end_sum;

		}
	}


	
	
};

//void check_polymer_bundle_diffusable_direction(list<Polymer_Unit>::iterator polymer_unit_iter)
//{
//
//	list<Polymer_Bundle>::iterator diffusion_bundle_iter;
//	diffusion_bundle_iter = polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer;
//	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
//	polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
//	//up case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position + 1;
//			if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//down case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position - 1;
//			if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//left case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position - 1;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
//			if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//right case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position + 1;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
//			if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//
//
//
//
//
//
//};


//void check_polymer_bundle_diffusable_direction(list<Polymer_Bundle>::iterator diffusion_bundle_iter)
//{
//
//	//list<Polymer_Bundle>::iterator diffusion_bundle_iter;
//	//diffusion_bundle_iter = polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer;
//	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
//	polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
//	//up case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position + 1;
//			if (Grid[col_position][row_position].grid_polymer_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_bundle_pointer != rule_structure.polymer_bundle_list_end.begin())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//down case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position - 1;
//			if (Grid[col_position][row_position].grid_polymer_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_bundle_pointer != rule_structure.polymer_bundle_list_end.begin())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//left case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position - 1;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
//			if (Grid[col_position][row_position].grid_polymer_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_bundle_pointer != rule_structure.polymer_bundle_list_end.begin())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//	//right case
//	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
//	{
//		list<Polymer_Unit>::iterator polymer_unit_iter;
//		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
//		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
//		{
//			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position + 1;
//			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
//			if (Grid[col_position][row_position].grid_polymer_bundle_pointer != diffusion_bundle_iter &&
//				Grid[col_position][row_position].grid_polymer_bundle_pointer != rule_structure.polymer_bundle_list_end.begin())//there is something there
//			{
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
//				Grid[col_position][row_position].grid_polymer_bundle_pointer->available_polymer_diffusion_direction_sum--;
//			}
//
//			polymer_unit_iter++;
//		}
//		polymer_cluster_iter++;
//	}
//
//
//
//
//
//
//};






void check_polymer_bundle_diffusable_direction(list<Polymer_Bundle>::iterator diffusion_bundle_iter,int check_mark)
{

	//list<Polymer_Bundle>::iterator diffusion_bundle_iter;
	//diffusion_bundle_iter = polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer;
	//diffusion_bundle_iter->availabe_diffusion_direction[direction_up] = true;
	//diffusion_bundle_iter->availabe_diffusion_direction[direction_down] = true;
	//diffusion_bundle_iter->availabe_diffusion_direction[direction_left] = true;
	//diffusion_bundle_iter->availabe_diffusion_direction[direction_right] = true;
	

	//rule_parameters.polymer_diffusion_propensity -= diffusion_bundle_iter->available_polymer_diffusion_direction_sum;

	if (diffusion_bundle_iter->Consider_mark != rule_parameters.t && check_mark!=0)
	{
		diffusion_bundle_iter->Consider_mark = rule_parameters.t;
		
		
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
			//output();
			//use step=2 to enable possible diffusion and polymerize

			polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();


			while (polymer_cluster_iter != diffusion_bundle_iter->bundle_cluster_pointer_list.end())
			{
				if ((*polymer_cluster_iter)->polymer_sequence.size() == 0)
				{
					rule_structure.polymer_cluster_list.erase(*polymer_cluster_iter);
					polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.erase(polymer_cluster_iter);
				}
				else
				{
					list<Polymer_Unit>::iterator polymer_unit_iter;
					polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
					for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)
					{
						if (1)
						{
							if (Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->Hydrolysis_mark != -1)
							{//if there is hydrolisys_mark==-1 problem

							 //first check if the stored annealing information is still valid
								if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
								{//if it is indeed in one of the annealing pair//first direction
									list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
									int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
									int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
									if ((col_distance*col_distance + row_distance*row_distance) != 1 ||
										annealing_pair_iter->first->polymer_cluster_pointer->direction != annealing_pair_iter->second->polymer_cluster_pointer->direction ||
										annealing_pair_iter->first->polymer_cluster_pointer->polarize != annealing_pair_iter->second->polymer_cluster_pointer->polarize)
									{	//it is no longer in the annealing position,need to delete this annealing relation
										if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
										{
											rule_parameters.TT_annealing_sum--;
										}
										else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
										{

										}
										else
										{
											rule_parameters.TD_annealing_sum--;
										}
										annealing_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
										annealing_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
										rule_structure.anealing_polymer_unit_pair_list.erase(annealing_pair_iter);



									}
								}
								if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
								{//if it is indeed in one of the annealing pair//second direction
									list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second];
									int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
									int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
									if ((col_distance*col_distance + row_distance*row_distance) != 1 ||
										annealing_pair_iter->first->polymer_cluster_pointer->direction != annealing_pair_iter->second->polymer_cluster_pointer->direction ||
										annealing_pair_iter->first->polymer_cluster_pointer->polarize != annealing_pair_iter->second->polymer_cluster_pointer->polarize)
									{	//it is no longer in the annealing position,need to delete this annealing relation
										if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
										{
											rule_parameters.TT_annealing_sum--;
										}
										else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
										{

										}
										else
										{
											rule_parameters.TD_annealing_sum--;
										}
										annealing_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
										annealing_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
										rule_structure.anealing_polymer_unit_pair_list.erase(annealing_pair_iter);



									}
								}

								if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.bundling_polymer_unit_pair_list.end())
								{//if there is a bundle on its first direction, left or down
									list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter;
									bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
									int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
									int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
									if ((col_distance*col_distance + row_distance*row_distance) != 1 ||
										bundle_pair_iter->first->polymer_cluster_pointer->direction != bundle_pair_iter->second->polymer_cluster_pointer->direction)
									{	//it is no longer in the bundling position,need to delete this bundling relation
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
								}

								if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.bundling_polymer_unit_pair_list.end())
								{//if there is a bundle on its second direction, right or up
									list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second];
									int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
									int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
									if ((col_distance*col_distance + row_distance*row_distance) != 1 ||
										bundle_pair_iter->first->polymer_cluster_pointer->direction != bundle_pair_iter->second->polymer_cluster_pointer->direction)
									{	//it is no longer in the bundling position,need to delete this bundling relation
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
								}



								//need to check attatch anchor information here
								//if (Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
								//{//if there is a anchor on the grid
								//	list<Anchor_Unit>::iterator anchor_iter = Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer;
								//	
								//	if (anchor_iter->attatch_mark == 0 && &*anchor_iter->polymer_unit_pointer==&*rule_structure.polymer_unit_list_end.begin())
								//	{//it is not able to attatch before and it does not link to any polymer_unit yet
								//		anchor_iter->attatch_mark = 1;
								//		rule_parameters.emtpy_anchor_attatch_sum++;
								//		//here i need to modify @ anchoring situation: there is anchor and polymer_unit on the same place on grid, but they are not linked together
								//		//in this case anchoring can't happened here, it need to do attatch
								//		//rule_parameters.empty_anchor_unit_num--;
								//	}

								//	if (anchor_iter->attatch_mark == 1 && &*anchor_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
								//	{//it is able to attatch before and it does link to any polymer_unit yet
								//		anchor_iter->attatch_mark = 0;
								//		rule_parameters.emtpy_anchor_attatch_sum--;
								//		//rule_parameters.empty_anchor_unit_num++;
								//	}

								//}



								//check_anchor_diffusion(polymer_unit_iter->polymer_grid_pointer->grid_anchor_pointer);
								//}
								if (polymer_unit_iter->attatch_mark == 1)
								{//i can do attatch before
									if (Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
										polymer_unit_iter->anchor_pointer == rule_structure.anchor_list.end())
									{//it still has a anchor on its bottom and the polymer_unit does not link to anything
									 //do nothing
									}
									else
									{
										list<Anchor_Unit>::iterator anchor_iter = Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer;
										anchor_iter->attatch_mark = 0;
										polymer_unit_iter->attatch_mark = 0;
										rule_parameters.emtpy_anchor_attatch_sum--;
									}
								}

								if (polymer_unit_iter->attatch_mark == 0)
								{//i can't do attatch before
									if (Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
										polymer_unit_iter->anchor_pointer == rule_structure.anchor_list.end())
									{//it now has a anchor on its bottom
										list<Anchor_Unit>::iterator anchor_iter = Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer;
										anchor_iter->attatch_mark = 1;
										polymer_unit_iter->attatch_mark = 1;
										rule_parameters.emtpy_anchor_attatch_sum++;
									}
									else
									{
										//do nothing
									}
								}




								//up case
								{
									int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
									int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
									int new_col_position = col_position + 1;
									if (new_col_position > rule_parameters.column_num - 1)
										new_col_position = 0;
									int next_new_col_position = new_col_position + 1;
									if (next_new_col_position > rule_parameters.column_num - 1)
										next_new_col_position = 0;
									//check_anchor_diffusion(Grid[new_col_position][row_position].grid_anchor_pointer);


									if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
									{//if up one step is empty
										if (Grid[col_position][row_position].grid_polymer_unit_pointer == (--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end()) &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if it is at the top end
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == false)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = true;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//up but polarize is direction down,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
														rule_parameters.TT_polymerize_end_vertical_minus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
														rule_parameters.TD_polymerize_end_vertical_minus_sum++;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
														rule_parameters.TT_polymerize_end_vertical_plus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
														rule_parameters.TD_polymerize_end_vertical_plus_sum++;
													}
												}
												
											}
										}

										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == false)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = true;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
											rule_parameters.polymer_diffusion_propensity+=calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,direction_up);
										}
										if (&*Grid[next_new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())//step two, there are something there
										{
											if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer == Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											{//if up step 2 is at is bottom and it is vertical
												if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == false)
												{
													Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = true;
													auto temp_polymer_unit_iter = Grid[next_new_col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//down but polarize is direction up,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
															rule_parameters.TT_polymerize_end_vertical_minus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
															rule_parameters.TD_polymerize_end_vertical_minus_sum++;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
															rule_parameters.TT_polymerize_end_vertical_plus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
															rule_parameters.TD_polymerize_end_vertical_plus_sum++;
														}
													}
												}
											}

											if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == false &&
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 1)
											{
												/*Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = true;
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
												rule_parameters.polymer_diffusion_propensity++;*/
												check_polymer_bundle_diffusable_direction(Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, check_mark - 1);
											}
										}
									}
								}


								//down case
								{
									int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
									int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
									int new_col_position = col_position - 1;
									if (new_col_position < 0)
										new_col_position = rule_parameters.column_num - 1;
									int next_new_col_position = new_col_position - 1;
									if (next_new_col_position < 0)
										next_new_col_position = rule_parameters.column_num - 1;
									//check_anchor_diffusion(Grid[new_col_position][row_position].grid_anchor_pointer);
									if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
									{
										if (Grid[col_position][row_position].grid_polymer_unit_pointer == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if it is at the bottom begin
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == false)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = true;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//down but polarize is direction up,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
														rule_parameters.TT_polymerize_end_vertical_minus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
														rule_parameters.TD_polymerize_end_vertical_minus_sum++;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
														rule_parameters.TT_polymerize_end_vertical_plus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
														rule_parameters.TD_polymerize_end_vertical_plus_sum++;
													}
												}
											}
										}

										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == false)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = true;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
											rule_parameters.polymer_diffusion_propensity +=calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,direction_down);
										}
										if (&*Grid[next_new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())//step two, there are something there
										{
											if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer == Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											{
												if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == false)
												{
													Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = true;
													auto temp_polymer_unit_iter = Grid[next_new_col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//up but polarize is direction down,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
															rule_parameters.TT_polymerize_end_vertical_minus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
															rule_parameters.TD_polymerize_end_vertical_minus_sum++;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
															rule_parameters.TT_polymerize_end_vertical_plus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
															rule_parameters.TD_polymerize_end_vertical_plus_sum++;
														}
													}
												}
											}

											if (Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == false &&
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 1)
											{
												/*Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = true;
												Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
												rule_parameters.polymer_diffusion_propensity++;*/
												check_polymer_bundle_diffusable_direction(Grid[next_new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, check_mark - 1);
											}
										}
									}
								}

								//left case
								{
									int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
									int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
									int new_row_position = row_position - 1;
									/*if (new_row_position < 0)
									new_row_position = row_num - 1;*/
									int next_new_row_position = new_row_position - 1;
									/*if (next_new_row_position < 0)
									next_new_row_position = row_num - 1;*/
									/*if(new_row_position>=0)
									check_anchor_diffusion(Grid[col_position][new_row_position].grid_anchor_pointer);*/
									if (new_row_position >= 0 && &*Grid[col_position][new_row_position].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
									{
										if (Grid[col_position][row_position].grid_polymer_unit_pointer == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
										{//if it is at the left begin
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == false)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = true;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//left but polarize is direction right,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
														rule_parameters.TT_polymerize_end_horizontal_minus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
														rule_parameters.TD_polymerize_end_horizontal_minus_sum++;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
														rule_parameters.TT_polymerize_end_horizontal_plus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
														rule_parameters.TD_polymerize_end_horizontal_plus_sum++;
													}
												}
											}
										}

										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == false)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = true;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
											rule_parameters.polymer_diffusion_propensity +=  calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,direction_left);
										}
										if (next_new_row_position >= 0 && &*Grid[col_position][next_new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
										{//step two, there are something there
											if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer == Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{
												if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == false)
												{
													Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = true;
													auto temp_polymer_unit_iter = Grid[col_position][next_new_row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//right but polarize is direction left,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum++;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum++;
														}
													}
												}
											}

											if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == false &&
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 1)
											{
												/*Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = true;
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
												rule_parameters.polymer_diffusion_propensity++;*/
												check_polymer_bundle_diffusable_direction(Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, check_mark - 1);
											}
										}
									}
								}
								//right case

								{
									int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
									int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
									int new_row_position = row_position + 1;
									if (new_row_position > rule_parameters.row_num - 1)
										new_row_position = 0;
									int next_new_row_position = new_row_position + 1;
									/*if (next_new_row_position > row_num - 1)
									next_new_row_position = 0;*/
									/*if(new_row_position<row_num)
									check_anchor_diffusion(Grid[col_position][new_row_position].grid_anchor_pointer);*/
									if (new_row_position<rule_parameters.row_num && &*Grid[col_position][new_row_position].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
									{
										if (Grid[col_position][row_position].grid_polymer_unit_pointer == --Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
										{//if it is at the right end
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == false)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = true;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//right but polarize is direction left,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
														rule_parameters.TT_polymerize_end_horizontal_minus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
														rule_parameters.TD_polymerize_end_horizontal_minus_sum++;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
														rule_parameters.TT_polymerize_end_horizontal_plus_sum++;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
														rule_parameters.TD_polymerize_end_horizontal_plus_sum++;
													}
												}
											}
										}

										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == false)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = true;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
											rule_parameters.polymer_diffusion_propensity +=  calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,direction_right);
										}
										if (next_new_row_position<rule_parameters.row_num && &*Grid[col_position][next_new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())//step two, there are something there
										{
											if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer == Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{
												if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == false)
												{
													Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = true;
													auto temp_polymer_unit_iter = Grid[col_position][next_new_row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//left but polarize is direction right,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum++;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum++;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum++;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum++;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum++;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum++;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum++;
														}
													}
												}
											}

											if (Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == false &&
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 1)
											{
												/*Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = true;
												Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum++;
												rule_parameters.polymer_diffusion_propensity++;*/
												check_polymer_bundle_diffusable_direction(Grid[col_position][next_new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, check_mark - 1);
											}
										}
									}
								}
							}
						}

						//polymer_unit_iter->Consider_mark = rule_parameters.t + 1;
						polymer_unit_iter++;
					}
					polymer_cluster_iter++;
				}


			}

			//use step=1 to disable possible diffusion and polymerize
			polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
			for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
			{
				list<Polymer_Unit>::iterator polymer_unit_iter;
				int action_mark = 0;//use this to mark bundling relation;
				polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
				for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
				{
					if (1)
					{
						//first consider anchor that might collide under diffusion
						if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
						{//if there is indeed a anchor
							int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
							int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;

							int up_col_position = col_position + 1;
							if (up_col_position > rule_parameters.column_num - 1)
								up_col_position = 0;
							int next_up_col_position = up_col_position + 1;
							if (next_up_col_position > rule_parameters.column_num - 1)
								next_up_col_position = 0;
							int down_col_position = col_position - 1;
							if (down_col_position < 0)
								down_col_position = rule_parameters.column_num - 1;
							int next_down_col_position = down_col_position - 1;
							if (next_down_col_position < 0)
								next_down_col_position = rule_parameters.column_num - 1;
							int right_row_position = row_position + 1;
							if (right_row_position > rule_parameters.row_num - 1)
								right_row_position = 0;
							int next_right_row_position = right_row_position + 1;
							if (next_right_row_position > rule_parameters.row_num - 1)
								next_right_row_position = 0;
							int left_row_position = row_position - 1;
							if (left_row_position < 0)
								left_row_position = rule_parameters.row_num - 1;
							int next_left_row_position = left_row_position - 1;
							if (next_left_row_position < 0)
								next_left_row_position = rule_parameters.row_num - 1;
							if (Grid[up_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
								&*Grid[up_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
							{//if there is a empty anchor on the up direction, disable diffusion up
								if (polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == true)
								{
									rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer, direction_up);
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = false;
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									
								}

							}
							if (Grid[down_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
								&*Grid[down_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
							{//if there is a empty anchor on the down direction, disable diffusion up
								if (polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == true)
								{
									rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer, direction_down);
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = false;
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									
								}

							}
							if (Grid[col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
								&*Grid[col_position][left_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
							{//if there is a empty anchor on the left direction, disable diffusion up
								if (polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == true)
								{
									rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer, direction_left);
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = false;
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									
								}

							}
							if (Grid[col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
								&*Grid[col_position][right_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
							{//if there is a empty anchor on the right direction, disable diffusion up
								if (polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == true)
								{
									rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer, direction_right);
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = false;
									polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									
								}

							}


						}



						if (Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->Hydrolysis_mark != -1)
						{
							//up case
							{



								int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
								int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
								int new_col_position = col_position + 1;
								if (new_col_position > rule_parameters.column_num - 1)
									new_col_position = 0;






								//belong to different bundle, it is for diffusion
								if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
								{//if there is a monomer in the upper side and they don't belong to the same bundle
									if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == yes)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_up);
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = no;
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										
									}
									//I need to consider the upperside bundle
									if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == yes)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_down);
										Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
										Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										

									}

								}
								//belong to either the same or different budnel, it is for polymerize
								if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
								{//if there is some monomer in the up directoin, not necessary to be different bundle
									if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
									{//if they don't belong to the same cluster
									 //i need to consider polymerize here
									 //consider top case
										action_mark = 0;
										if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if current polymer_unit is the top end and it is vertical
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == yes)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = no;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//up but polarize is direction down,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_minus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_minus_sum--;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_plus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_plus_sum--;
													}
												}

											}
											action_mark++;//mark this pair invalid;
										}
										if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
											Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if the upper one is on the bottom end and it is vertical
											if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == yes)
											{
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = no;
												auto temp_polymer_unit_iter = Grid[new_col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//down but polarize is direction up,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_minus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_minus_sum--;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_plus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_plus_sum--;
													}
												}

											}
											action_mark++;//mark it invalid
														  //rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[new_col_position][row_position].grid_polymer_unit_pointer));
														  //Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
														  //Grid[new_col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
										}
										if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
										{
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] == rule_structure.bundling_polymer_unit_pair_list.end()
												&& Grid[col_position][row_position].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_second] == rule_structure.debundling_polymer_unit_pair_list.end()
												&& Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal
												&&Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{
												rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[new_col_position][row_position].grid_polymer_unit_pointer));
												Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.begin();
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.begin();
												auto bundle_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
												if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
												{
													rule_parameters.bundling_sum++;
												}
												else
												{
													rule_parameters.bundling_inverse_sum++;
												}

											}


										}
										//add annealing relation here
										//if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
										//	--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
										//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										//{//if current polymer_unit is the top end and it is vertical

										//	if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
										//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
										//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										//	{//if the upper one is on the bottom end and it is vertical

										//		
										//	}
										//}
										if (action_mark == 2 &&
											Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize)
										{
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] == rule_structure.anealing_polymer_unit_pair_list.end())
											{
												rule_structure.anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[new_col_position][row_position].grid_polymer_unit_pointer));
												Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.begin();
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.begin();

												list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second];

												if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
												{
													rule_parameters.TT_annealing_sum++;
												}
												else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
												{

												}
												else
												{
													rule_parameters.TD_annealing_sum++;
												}
											}

										}
									}
									else
									{//if they belong to the same cluster, the up side is bottom end and down side is up end
										if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if current polymer_unit is the top end and it is vertical

											if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											{//if the upper one is on the bottom end and it is vertical
											 //this case can only happened under periodic conditions//ring case
												if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == true)
												{
													Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = false;
													auto temp_polymer_unit_iter = Grid[new_col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//down but polarize is direction up,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_plus_sum--;
														}
													}
												}
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == true)
												{
													Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = false;
													auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//up but polarize is direction down,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_plus_sum--;
														}
													}
												}

											}
										}

									}
								}






							}


							//down case
							{
								int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
								int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
								int new_col_position = col_position - 1;
								if (new_col_position <0)
									new_col_position = rule_parameters.column_num - 1;

								//belong to different bundle, it is for diffusion
								if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
								{//if there is a monomer in the down side and they don't belong to the same bundle
									if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == yes)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_down);
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										
									}
									//cout << new_col_position << endl;
									//cout << row_position << endl;
									//I need to consider the downside bundle
									if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == yes)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_up);
										Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = no;
										Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										

									}

								}
								//belong to either the same or different budnel, it is for polymerize
								if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
								{//if there is some monomer in the down directoin, not necessary to be different bundle
									if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
									{//if they don't belong to the same cluster
									 //i need to consider polymerize here
									 //consider top case
										action_mark = 0;
										if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if current polymer_unit is the bottom end and it is vertical
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == yes)
											{
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = no;
												auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//down but polarize is direction up,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_minus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_minus_sum--;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_plus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_plus_sum--;
													}
												}

											}
											action_mark++;//mark this pair invalid;
										}
										if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
											--Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if the down one is on the upper end and it is vertical
											if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == yes)
											{
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = no;
												auto temp_polymer_unit_iter = Grid[new_col_position][row_position].grid_polymer_unit_pointer;
												if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{//up but polarize is direction down,so it is the minus end
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_minus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_minus_sum--;
													}

												}
												else
												{
													if (temp_polymer_unit_iter->Hydrolysis_mark == T)
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
														rule_parameters.TT_polymerize_end_vertical_plus_sum--;
													}
													else
													{
														temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
														rule_parameters.TD_polymerize_end_vertical_plus_sum--;
													}
												}

											}
											action_mark++;//mark it invalid
														  //rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[new_col_position][row_position].grid_polymer_unit_pointer));
														  //Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
														  //Grid[new_col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
										}
										if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
										{
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] == rule_structure.bundling_polymer_unit_pair_list.end()
												&& Grid[col_position][row_position].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_first] == rule_structure.debundling_polymer_unit_pair_list.end()
												&& Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal
												&&Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{

												rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
												Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.begin();
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.begin();
												auto bundle_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
												if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
												{
													rule_parameters.bundling_sum++;
												}
												else
												{
													rule_parameters.bundling_inverse_sum++;
												}

											}

										}
										//add annealing relation here
										//if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
										//	--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
										//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										//{//if current polymer_unit is the top end and it is vertical

										//	if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
										//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
										//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										//	{//if the upper one is on the bottom end and it is vertical

										//		
										//	}
										//}
										if (action_mark == 2 &&
											Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize)
										{
											if (Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] == rule_structure.anealing_polymer_unit_pair_list.end())
											{
												rule_structure.anealing_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
												Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.begin();
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.begin();

												list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];

												if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
												{
													rule_parameters.TT_annealing_sum++;
												}
												else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
												{

												}
												else
												{
													rule_parameters.TD_annealing_sum++;
												}
											}

										}
									}
									else
									{//if the belong to the same cluster, the up side is bottom end and down side is up end
										if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
										{//if current polymer_unit is the bottom end and it is vertical

											if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
												--Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
												Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											{//if the downside one is on the top end and it is vertical
											 //this case can only happened under periodic conditions
												if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == true)
												{
													Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = false;
													auto temp_polymer_unit_iter = Grid[new_col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//up but polarize is direction down,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_plus_sum--;
														}
													}
												}
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == true)
												{
													Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = false;
													auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//down but polarize is direction up,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_vertical_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_vertical_plus_sum--;
														}
													}
												}

											}
										}

									}
								}


								//commented material
								{
									//if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									//	Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
									//{//not the same bundle
									//	if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == yes)
									//	{
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//	}
									//	//I need to consider the downside bundle
									//	if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == yes)
									//	{
									//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = no;
									//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//	}
									//	//i need to consider polymerize here
									//	//consider down case
									//	action_mark = 0;
									//	if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
									//	{//if the upper one is on its bottom end
									//		if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] == yes)
									//		{
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = no;
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymerize_end_sum--;
									//			rule_parameters.polymerize_end_sum--;
									//		}
									//		action_mark++;

									//	}
									//	if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
									//		--Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
									//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
									//	{//the downside one is on its upper end
									//		if (Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] == yes)
									//		{
									//			Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = no;
									//			Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymerize_end_sum--;
									//			rule_parameters.polymerize_end_sum--;
									//		}
									//		action_mark++;
									//		//rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//		//Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//		//Grid[new_col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//	}
									//	if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
									//	{
									//		rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_up] = rule_structure.bundling_polymer_unit_pair_list.begin();
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_down] = rule_structure.bundling_polymer_unit_pair_list.begin();
									//	}
									//	if (&*Grid[new_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									//		Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
									//	{//if there is a monomer on the top and it does not belong to the same cluster, add annealing information
									//		if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
									//		{//if current polymer_unit is the bottom end and it is vertical

									//			if (Grid[new_col_position][row_position].grid_polymer_unit_pointer ==
									//				--Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
									//				Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
									//			{//if the downside one is on the top end and it is vertical

									//				rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//				Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//				Grid[new_col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//			}
									//		}
									//	}
									//}
								}
							}

							//left case
							{
								int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
								int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;

								if (row_position == 0)
								{//it is at its left boundary, disable diffusion to left and polymerize to left at any case
									if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == true)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_left);
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = false;
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										
									}
									if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									{
										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == true)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = false;
											auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
											if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
											{
												if (temp_polymer_unit_iter->Hydrolysis_mark == T)
												{//left but polarize is direction right,so it is the minus end
													temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
													rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
												}
												else
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
													rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
												}

											}
											else
											{
												if (temp_polymer_unit_iter->Hydrolysis_mark == T)
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
													rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
												}
												else
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
													rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
												}
											}

										}
									}



								}
								else
								{
									int new_row_position = row_position - 1;
									/*if (new_row_position <0)
									new_row_position = row_num - 1;*/
									//belong to different bundle, it is for diffusion
									if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
										Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
									{//if there is a monomer in the left side and they don't belong to the same bundle
										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == yes)
										{
											rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_left);
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = no;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
											
										}
										//I need to consider the leftside bundle
										if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == yes)
										{
											rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_right);
											Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = no;
											Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
											

										}

									}
									//belong to either the same or different budnel, it is for polymerize
									if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
									{//if there is some monomer in the left direction, not necessary to be different bundle
										if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
										{//if they don't belong to the same cluster
										 //i need to consider polymerize here
										 //consider top case
											action_mark = 0;
											if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if current polymer_unit is the left end and it is horizontal
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == yes)
												{
													Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = no;
													auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//left but polarize is direction right,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
														}
													}

												}
												action_mark++;//mark this pair invalid;
											}
											if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
												--Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
												Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if the left one is on the right end and it is horizontal
												if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == yes)
												{
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = no;
													auto temp_polymer_unit_iter = Grid[col_position][new_row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//right but polarize is direction left,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
														}
													}

												}
												action_mark++;//mark it invalid
															  //rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][new_row_position].grid_polymer_unit_pointer));
															  //Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
															  //Grid[col_position][new_row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
											}
											if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
											{
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] == rule_structure.bundling_polymer_unit_pair_list.end()
													&& Grid[col_position][row_position].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_first] == rule_structure.debundling_polymer_unit_pair_list.end()
													&& Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical
													&&Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
												{
													rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][new_row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
													Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.begin();
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.begin();
													auto bundle_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
													if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
													{
														rule_parameters.bundling_sum++;
													}
													else
													{
														rule_parameters.bundling_inverse_sum++;
													}
												}


											}
											//add annealing relation here
											//if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											//	--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											//{//if current polymer_unit is the top end and it is vertical

											//	if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
											//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											//	{//if the upper one is on the bottom end and it is vertical

											//		
											//	}
											//}
											if (action_mark == 2 &&
												Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize)
											{
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] == rule_structure.anealing_polymer_unit_pair_list.end())
												{
													rule_structure.anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][new_row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
													Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.begin();
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.begin();
													list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];

													if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
													{
														rule_parameters.TT_annealing_sum++;
													}
													else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
													{

													}
													else
													{
														rule_parameters.TD_annealing_sum++;
													}
												}

											}
										}
										else
										{//if the belong to the same cluster, the right side is left end and left side is right end
											if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if current polymer_unit is the left end and it is horizontal

												if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
													--Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
												{//if the left one is on the right end and it is horizontal
												 //this case can only happened under periodic conditions
													if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == true)
													{
														Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = false;
														auto temp_polymer_unit_iter = Grid[col_position][new_row_position].grid_polymer_unit_pointer;
														if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{//right but polarize is direction left,so it is the minus end
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
															}

														}
														else
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
															}
														}
													}
													if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == true)
													{
														Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = false;
														auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
														if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{//left but polarize is direction right,so it is the minus end
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
															}

														}
														else
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
															}
														}
													}

												}
											}

										}
									}

								}
								//commented material
								{
									//if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									//	Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
									//{
									//	if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == yes)
									//	{
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = no;
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//	}
									//	//I need to consider the leftside bundle
									//	if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == yes)
									//	{
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = no;
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//	}


									//	//i need to consider polymerize here
									//	//consider down case
									//	action_mark = 0;
									//	if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									//	{//if the right one is on its left end
									//		if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == yes)
									//		{
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = no;
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymerize_end_sum--;
									//			rule_parameters.polymerize_end_sum--;
									//		}
									//		action_mark++;

									//	}
									//	if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
									//		--Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									//	{//the leftside one is on its right end
									//		if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == yes)
									//		{
									//			Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = no;
									//			Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymerize_end_sum--;
									//			rule_parameters.polymerize_end_sum--;
									//		}
									//		action_mark++;
									//		//rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[new_col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//		//Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//		//Grid[new_col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//	}
									//	if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
									//	{//now it is the left-right| cases
									//		rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][new_row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_right] = rule_structure.bundling_polymer_unit_pair_list.begin();
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_left] = rule_structure.bundling_polymer_unit_pair_list.begin();
									//	}
									//	if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
									//	{//if there is a monomer on the left and right does not belong to the same cluster, add annealing information
									//		if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
									//			Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									//		{//if current polymer_unit is the left end and it is horizontal

									//			if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
									//				--Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
									//				Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									//			{//if the leftside one is on the right end and it is horizontal

									//				rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][new_row_position].grid_polymer_unit_pointer, Grid[col_position][row_position].grid_polymer_unit_pointer));
									//				Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//				Grid[col_position][new_row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
									//			}
									//		}
									//	}

									//}
								}
							}


							//right case
							{

								int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
								int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;



								if (row_position == rule_parameters.row_num - 1)
								{//it is at its right boundary, disable diffusion to right and polymerize to right at any case
									if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == true)
									{
										rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_right);
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = false;
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
										
									}
									if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
										--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
										Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
									{
										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == true)
										{
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = false;
											auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
											if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
											{
												if (temp_polymer_unit_iter->Hydrolysis_mark == T)
												{//right but polarize is direction left,so it is the minus end
													temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
													rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
												}
												else
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
													rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
												}

											}
											else
											{
												if (temp_polymer_unit_iter->Hydrolysis_mark == T)
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
													rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
												}
												else
												{
													temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
													rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
												}
											}

										}
									}



								}
								else
								{
									/*int new_row_position = row_position - 1;*/
									int new_row_position = row_position + 1;
									/*	if (new_row_position > row_num - 1)
									new_row_position = 0;*/


									/*if (new_row_position <0)
									new_row_position = row_num - 1;*/
									//belong to different bundle, it is for diffusion
									if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
										Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
									{//if there is a monomer in the right side and they don't belong to the same bundle
										if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == yes)
										{
											rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_right);
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = no;
											Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
											
										}
										//I need to consider the rightside bundle
										if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == yes)
										{
											rule_parameters.polymer_diffusion_propensity = rule_parameters.polymer_diffusion_propensity - calculate_diffusion_rate(Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_left);
											Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = no;
											Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
											

										}

									}
									//belong to either the same or different budnel, it is for polymerize
									if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
									{//if there is some monomer in the right direction, not necessary to be different bundle
										if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer != Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer)
										{//if they don't belong to the same cluster
										 //i need to consider polymerize here
										 //consider right case
											action_mark = 0;
											if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
												--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if current polymer_unit is the right end and it is horizontal
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == yes)
												{
													Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = no;
													auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//right but polarize is direction left,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
														}
													}

												}
												action_mark++;//mark this pair invalid;
											}
											if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
												Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
												Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if the left one is on the right end and it is horizontal
												if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == yes)
												{
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = no;
													auto temp_polymer_unit_iter = Grid[col_position][new_row_position].grid_polymer_unit_pointer;
													if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{//left but polarize is direction right,so it is the minus end
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
														}

													}
													else
													{
														if (temp_polymer_unit_iter->Hydrolysis_mark == T)
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
															rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
														}
														else
														{
															temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
															rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
														}
													}

												}
												action_mark++;//mark it invalid
															  //rule_structure.TT_anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][new_row_position].grid_polymer_unit_pointer));
															  //Grid[col_position][row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
															  //Grid[col_position][new_row_position].grid_polymer_unit_pointer->TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
											}
											if (action_mark == 0)//if it is still 0 means it is none of the cases above: it is not on end both for up and down
											{
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] == rule_structure.bundling_polymer_unit_pair_list.end()
													&& Grid[col_position][row_position].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_second] == rule_structure.debundling_polymer_unit_pair_list.end()
													&& Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical
													&&Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
												{
													rule_structure.bundling_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][new_row_position].grid_polymer_unit_pointer));
													Grid[col_position][row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] = rule_structure.bundling_polymer_unit_pair_list.begin();
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] = rule_structure.bundling_polymer_unit_pair_list.begin();
													auto bundle_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
													if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
													{
														rule_parameters.bundling_sum++;
													}
													else
													{
														rule_parameters.bundling_inverse_sum++;
													}
												}

											}
											//add annealing relation here
											//if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
											//	--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
											//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											//{//if current polymer_unit is the top end and it is vertical

											//	if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
											//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
											//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
											//	{//if the upper one is on the bottom end and it is vertical

											//		
											//	}
											//}
											if (action_mark == 2 &&
												Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize == Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polarize)
											{
												if (Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] == rule_structure.anealing_polymer_unit_pair_list.end())
												{
													rule_structure.anealing_polymer_unit_pair_list.push_front(make_pair(Grid[col_position][row_position].grid_polymer_unit_pointer, Grid[col_position][new_row_position].grid_polymer_unit_pointer));
													Grid[col_position][row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.begin();
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.begin();
													list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_second];

													if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
													{
														rule_parameters.TT_annealing_sum++;
													}
													else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
													{

													}
													else
													{
														rule_parameters.TD_annealing_sum++;
													}
												}

											}
										}
										else
										{//if the belong to the same cluster, the right side is left end and left side is right end
											if (Grid[col_position][row_position].grid_polymer_unit_pointer ==
												--Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end() &&
												Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
											{//if current polymer_unit is the left end and it is horizontal

												if (Grid[col_position][new_row_position].grid_polymer_unit_pointer ==
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin() &&
													Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
												{//if the left one is on the right end and it is horizontal
												 //this case can only happened under periodic conditions
													if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] == true)
													{
														Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_left] = false;
														auto temp_polymer_unit_iter = Grid[col_position][new_row_position].grid_polymer_unit_pointer;
														if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{//left but polarize is direction right,so it is the minus end
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
															}

														}
														else
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
															}
														}
													}
													if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] == true)
													{
														Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->availabe_polymerize_direction[direction_right] = false;
														auto temp_polymer_unit_iter = Grid[col_position][row_position].grid_polymer_unit_pointer;
														if (temp_polymer_unit_iter->polymer_cluster_pointer->polarize == direction_first)
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{//right but polarize is direction left,so it is the minus end
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_minus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_minus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
															}

														}
														else
														{
															if (temp_polymer_unit_iter->Hydrolysis_mark == T)
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TT_polymerize_plus_end_sum--;
																rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
															}
															else
															{
																temp_polymer_unit_iter->polymer_cluster_pointer->TD_polymerize_plus_end_sum--;
																rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
															}
														}
													}

												}
											}

										}
									}

								}




								//commented material
								{
									//if (&*Grid[col_position][new_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
									//	Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter)//Grid[new_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
									//{
									//	if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == yes)
									//	{
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = no;
									//		Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//	}
									//	//I need to consider the upperside bundle
									//	if (Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == yes)
									//	{
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = no;
									//		Grid[col_position][new_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
									//		rule_parameters.polymer_diffusion_propensity--;
									//		//i need to consider polymerize here
									//	}

									//}
								}

								//int col_position = polymer_unit_iter->polymer_grid_pointer->col_position + 1;
								//int row_position = polymer_unit_iter->polymer_grid_pointer->row_position;
								//if (Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != diffusion_bundle_iter &&
								//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer != rule_structure.polymer_bundle_list.end())//there is something there
								//{
								//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
								//	Grid[col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
								//}
							}

						}
					}

					/*if (polymer_unit_iter->Consider_mark == rule_parameters.t + 1)
					{
					polymer_unit_iter->Consider_mark = rule_parameters.t;
					}*/

					polymer_unit_iter++;
				}
				polymer_cluster_iter++;
				/*int test_sum = 0;
				for (int i = 0; i < 4; i++)
				{
				test_sum += diffusion_bundle_iter->availabe_diffusion_direction[i];
				}
				if (test_sum != diffusion_bundle_iter->available_polymer_diffusion_direction_sum)
				{
				std::cout << "error: diffusion direction sum is not the same as sum up" << endl;
				}*/


			//now i enable all polymer no mater how many anchor it get
			/*if (diffusion_bundle_iter->polymer_anchor_list.size() > 1)
			{
				rule_parameters.polymer_diffusion_propensity -= diffusion_bundle_iter->available_polymer_diffusion_direction_sum;
				diffusion_bundle_iter->available_polymer_diffusion_direction_sum = 0;
				for (int i = 0; i < 4; i++)
				{
					diffusion_bundle_iter->availabe_diffusion_direction[i] = false;
				}

			}*/
		}
		
		
		
		
	}
	diffusion_bundle_iter->Consider_mark = 0;

};

void check_anchor_diffusion(list<Anchor_Unit>::iterator anchor_iter)
{
	if (anchor_iter != rule_structure.anchor_list.end() && anchor_iter->consider_mark != -rule_parameters.step_num)
	{//only when anchor_iter is exist, it has not be considered for diffusion

		if (&*anchor_iter->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
		{//if the anchor is not linked to the polymer_unit
			int col_position = anchor_iter->anchor_grid_pointer->col_position;
			int row_position = anchor_iter->anchor_grid_pointer->row_position;


			if (&*Grid[col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
			{//if there is indeed a polymer_unit there
				if (anchor_iter->attatch_mark == 0)
				{
					Grid[col_position][row_position].grid_polymer_unit_pointer->attatch_mark = 1;
					anchor_iter->attatch_mark = 1;
					rule_parameters.emtpy_anchor_attatch_sum++;
					//rule_parameters.empty_anchor_unit_num--;
				}

			}
			else
			{//if there is no polymer_unit above
				if (anchor_iter->attatch_mark == 1)
				{
					Grid[col_position][row_position].grid_polymer_unit_pointer->attatch_mark = 0;
					anchor_iter->attatch_mark = 0;
					rule_parameters.emtpy_anchor_attatch_sum--;
					//rule_parameters.empty_anchor_unit_num++;
				}

			}
			

			int up_col_position = col_position + 1;
			if (up_col_position > rule_parameters.column_num - 1)
				up_col_position = 0;
			int next_up_col_position = up_col_position + 1;
			if (next_up_col_position > rule_parameters.column_num - 1)
				next_up_col_position = 0;
			int down_col_position = col_position - 1;
			if (down_col_position < 0)
				down_col_position = rule_parameters.column_num - 1;
			int next_down_col_position = down_col_position - 1;
			if (next_down_col_position < 0)
				next_down_col_position = rule_parameters.column_num - 1;
			int right_row_position = row_position + 1;
			if (right_row_position > rule_parameters.row_num - 1)
				right_row_position = 0;
			int next_right_row_position = right_row_position + 1;
			if (next_right_row_position > rule_parameters.row_num - 1)
				next_right_row_position = 0;
			int left_row_position = row_position - 1;
			if (left_row_position < 0)
				left_row_position = rule_parameters.row_num - 1;
			int next_left_row_position = left_row_position - 1;
			if (next_left_row_position < 0)
				next_left_row_position = rule_parameters.row_num - 1;
			//free the diffusion
			rule_parameters.empty_anchor_diffusion_direction_sum -= anchor_iter->available_polymer_diffusion_direction_sum;
			anchor_iter->available_polymer_diffusion_direction_sum = 0;
			for (int i = 0; i < 4; i++)
			{
				anchor_iter->availabe_diffusion_direction[i] = false;
			}
			if (Grid[up_col_position][row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{//if up direction is empty, it can diffuse up//without anchor unit there
				anchor_iter->availabe_diffusion_direction[direction_up] = true;
				anchor_iter->available_polymer_diffusion_direction_sum++;
			}
			if (Grid[down_col_position][row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{//if up direction is empty, it can diffuse up
				anchor_iter->availabe_diffusion_direction[direction_down] = true;
				anchor_iter->available_polymer_diffusion_direction_sum++;
			}
			if (row_position != rule_parameters.row_num - 1 && Grid[col_position][right_row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{//if up direction is empty, it can diffuse right
				anchor_iter->availabe_diffusion_direction[direction_right] = true;
				anchor_iter->available_polymer_diffusion_direction_sum++;
			}
			if (row_position != 0 && Grid[col_position][left_row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{//if up direction is empty, it can diffuse right
				anchor_iter->availabe_diffusion_direction[direction_left] = true;
				anchor_iter->available_polymer_diffusion_direction_sum++;
			}
			rule_parameters.empty_anchor_diffusion_direction_sum += anchor_iter->available_polymer_diffusion_direction_sum;


			//now consider the neighboring anchor diffusion relation


			if (Grid[up_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[up_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] == true)
					{
						Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] = false;
						Grid[up_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					if (Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == true)
					{
						rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_down);
						Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = false;
						Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						
					}
				}
				
			}
			if (Grid[down_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[down_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] == true)
					{
						Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] = false;
						Grid[down_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					if (Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == true)
					{
						rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_up);
						Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = false;
						Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						
					}
				}

			}
			if (Grid[col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[col_position][left_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] == true)
					{
						Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] = false;
						Grid[col_position][left_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					if (Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == true)
					{
						rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_right);
						Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = false;
						Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						
					}
				}

			}
			if (Grid[col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[col_position][right_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] == true)
					{
						Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] = false;
						Grid[col_position][right_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					if (Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == true)
					{
						rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, direction_left);
						Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = false;
						Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						
					}
				}

			}
			
			//why do i need to do recursion here?
			Grid[col_position][row_position].grid_anchor_pointer->consider_mark = -rule_parameters.step_num;
			/*check_anchor_diffusion(Grid[up_col_position][left_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[up_col_position][right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[down_col_position][left_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[down_col_position][right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[next_up_col_position][row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[next_down_col_position][row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[col_position][next_right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[col_position][next_left_row_position].grid_anchor_pointer);
			Grid[col_position][row_position].grid_anchor_pointer->consider_mark = 0;*/
		}
		else
		{//if anchor already linked to a polymer_unit
			
			int col_position = anchor_iter->anchor_grid_pointer->col_position;
			int row_position = anchor_iter->anchor_grid_pointer->row_position;
			if (anchor_iter->attatch_mark == 1)
			{
				Grid[col_position][row_position].grid_polymer_unit_pointer->attatch_mark = 0;
				anchor_iter->attatch_mark = 0;
				rule_parameters.emtpy_anchor_attatch_sum--;
				//rule_parameters.empty_anchor_unit_num++;
			}
			//anchor_iter->polymer_unit_pointer->attatch_mark = 0;
			//anchor_iter->attatch_mark = 0;
			/*if (anchor_iter->available_polymer_diffusion_direction_sum != 0)
			{
				cout << "anchored polymer_unit has none zero empty anchor diffusion direction sum" << endl;
				system("pause");
			}*/
			rule_parameters.empty_anchor_diffusion_direction_sum -= anchor_iter->available_polymer_diffusion_direction_sum;
			for (int i = 0; i < 4; i++)
			{
				anchor_iter->availabe_diffusion_direction[i] = false;
			}
			anchor_iter->available_polymer_diffusion_direction_sum = 0;
			

			int up_col_position = col_position + 1;
			if (up_col_position > rule_parameters.column_num - 1)
				up_col_position = 0;
			int next_up_col_position = up_col_position + 1;
			if (next_up_col_position > rule_parameters.column_num - 1)
				next_up_col_position = 0;
			int down_col_position = col_position - 1;
			if (down_col_position < 0)
				down_col_position = rule_parameters.column_num - 1;
			int next_down_col_position = down_col_position - 1;
			if (next_down_col_position < 0)
				next_down_col_position = rule_parameters.column_num - 1;
			int right_row_position = row_position + 1;
			if (right_row_position > rule_parameters.row_num - 1)
				right_row_position = 0;
			int next_right_row_position = right_row_position + 1;
			if (next_right_row_position > rule_parameters.row_num - 1)
				next_right_row_position = 0;
			int left_row_position = row_position - 1;
			if (left_row_position < 0)
				left_row_position = rule_parameters.row_num - 1;
			int next_left_row_position = left_row_position - 1;
			if (next_left_row_position < 0)
				next_left_row_position = rule_parameters.row_num - 1;
			
			//if (Grid[up_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
			//	&*Grid[up_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
			//{//if there is a anchor on the up direction, disable the upper to diffuse down
			//	if (Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] == true)
			//	{
			//		Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] == false;
			//		Grid[up_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
			//		rule_parameters.empty_anchor_diffusion_direction_sum--;
			//	}
			//}
			//if (Grid[down_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
			//	&*Grid[down_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
			//{//if there is a anchor on the down direction, disable the upper to diffuse up
			//	if (Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] == true)
			//	{
			//		Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] == false;
			//		Grid[down_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
			//		rule_parameters.empty_anchor_diffusion_direction_sum--;
			//	}
			//}

			//if (Grid[col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
			//	&*Grid[col_position][left_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
			//{//if there is a anchor on the left direction, disable the upper to diffuse right
			//	if (Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] == true)
			//	{
			//		Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] == false;
			//		Grid[col_position][left_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
			//		rule_parameters.empty_anchor_diffusion_direction_sum--;
			//	}
			//}
			//if (Grid[col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
			//	&*Grid[col_position][right_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
			//{//if there is a anchor on the right direction, disable the upper to diffuse left
			//	if (Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] == true)
			//	{
			//		Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] == false;
			//		Grid[col_position][right_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
			//		rule_parameters.empty_anchor_diffusion_direction_sum--;
			//	}
			//}
			if (Grid[up_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[up_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] == true)
					{
						Grid[up_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_down] = false;
						Grid[up_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					/*if (Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] == true)
					{
						Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = false;
						Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.polymer_diffusion_propensity--;
					}*/
				}

			}
			if (Grid[down_col_position][row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[down_col_position][row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] == true)
					{
						Grid[down_col_position][row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_up] = false;
						Grid[down_col_position][row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					/*if (Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] == true)
					{
						Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = false;
						Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.polymer_diffusion_propensity--;
					}*/
				}

			}
			if (Grid[col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[col_position][left_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] == true)
					{
						Grid[col_position][left_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_right] = false;
						Grid[col_position][left_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					/*if (Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] == true)
					{
						Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_right] = false;
						Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.polymer_diffusion_propensity--;
					}*/
				}

			}
			if (Grid[col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
			{
				if (&*Grid[col_position][right_row_position].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//if there is a empty anchor on the up direction, disable the upper to diffuse down
					if (Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] == true)
					{
						Grid[col_position][right_row_position].grid_anchor_pointer->availabe_diffusion_direction[direction_left] = false;
						Grid[col_position][right_row_position].grid_anchor_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.empty_anchor_diffusion_direction_sum--;
					}
				}
				else
				{//if there is a not empty anchor on the up direction, disable the bundle to duffuse down
					/*if (Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] == true)
					{
						Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_left] = false;
						Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
						rule_parameters.polymer_diffusion_propensity--;
					}*/
				}

			}
			Grid[col_position][row_position].grid_anchor_pointer->consider_mark = -rule_parameters.step_num;
			/*check_anchor_diffusion(Grid[up_col_position][left_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[up_col_position][right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[down_col_position][left_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[down_col_position][right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[next_up_col_position][row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[next_down_col_position][row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[col_position][next_right_row_position].grid_anchor_pointer);
			check_anchor_diffusion(Grid[col_position][next_left_row_position].grid_anchor_pointer);
			Grid[col_position][row_position].grid_anchor_pointer->consider_mark = 0;*/
		}

	}
};

void set_action_mark_true(list < Polymer_Bundle>::iterator polymer_bundle_iter)
{
	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
	polymer_cluster_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
	for (int i = 0; i < polymer_bundle_iter->bundle_cluster_pointer_list.size(); i++)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
		{

			polymer_unit_iter->action_mark = true;




			polymer_unit_iter++;
		}
		polymer_cluster_iter++;
	}
}


void set_action_mark_false(list < Polymer_Bundle>::iterator polymer_bundle_iter)
{
	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
	polymer_cluster_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
	for (int i = 0; i < polymer_bundle_iter->bundle_cluster_pointer_list.size(); i++)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
		{

			polymer_unit_iter->action_mark = false;




			polymer_unit_iter++;
		}
		polymer_cluster_iter++;
	}
}









void check_bundling_pairs(list<Polymer_Bundle>::iterator bundle_iter)
{
	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
	polymer_cluster_iter = bundle_iter->bundle_cluster_pointer_list.begin();
	//up case
	for (int i = 0; i < bundle_iter->bundle_cluster_pointer_list.size(); i++)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
		for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
		{
			int col_position = polymer_unit_iter->polymer_grid_pointer->col_position;
			int row_position = polymer_unit_iter->polymer_grid_pointer->row_position + 1;
			//if (Grid[col_position][row_position].grid_polymer_bundle_pointer != diffusion_bundle_iter &&
			//	Grid[col_position][row_position].grid_polymer_bundle_pointer != rule_structure.polymer_bundle_list_end.begin())//there is something there
			//{
			//	Grid[col_position][row_position].grid_polymer_bundle_pointer->availabe_diffusion_direction[direction_down] = no;
			//	Grid[col_position][row_position].grid_polymer_bundle_pointer->available_polymer_diffusion_direction_sum--;
			//}

			polymer_unit_iter++;
		}
		polymer_cluster_iter++;
	}
};





void select_and_execute()
{

#ifdef DEBUG_RANDOM
	//random_device rd;
	default_random_engine gen(rule_parameters.step_num);
	//mt19937_64 gen(rd());
	
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	
	uniform_real_distribution <> rule_control(0, rule_parameters.sumps);
	//uniform_int_distribution<> test(0, 100);
	//uniform_real_distribution <> rule_control(0, 1000);
	
	//default_random_engine rule_random;
	rule_parameters.cumps =rule_parameters.ps[0];
	int irule = 0;
	//std::cout << "Check1" << endl;
	double r1;
	do
	{
		r1 = rule_control(gen);//rule_control is just to generate a random number;
	} while (r1 == 0);
	//int r2 = test(gen);
	//cout << r2 << endl;
	//system("pause");
	//r1 = r1*rule_parameters.sumps;
	while (r1 > rule_parameters.cumps)
	{
		irule++;
		rule_parameters.cumps = rule_parameters.cumps + rule_parameters.ps[irule];
	}
	rule_parameters.reaction_mark = irule;
	double start, end = 0;
	start = omp_get_wtime();
	/*if (irule == 17)
	{
		cout << "lalal";
	}*/
	if (rule_parameters.event_analysis==4 && irule != 0 && irule != 18 && irule != 10 && irule != 11)
	{
		rule_parameters.event_output <<setprecision(15)<< rule_parameters.t << " " << irule << " ";
		
	}
	/*if (irule == 20)
	{
		cout << "polymerize" << endl;
		system("pause");
	}*/
	double temp_diffusion_propensity = rule_parameters.polymer_diffusion_propensity;
	switch (irule)
	{
	//case rule_diffusion_anchor: diffuse_anchor();    //0 we have empty anchor diffuse;
	//	break;
	case rule_diffusion_polymer: 
	{
		polymer_diffusion();
	}
		break;
	case rule_TT_polymerize_horizontal_plus: 
	{
		TT_polymerize_horizontal_plus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
		break;
	case rule_TT_polymerize_horizontal_minus:
	{
		TT_polymerize_horizontal_minus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TT_polymerize_vertical_plus:
	{
		TT_polymerize_vertical_plus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TT_polymerize_vertical_minus:
	{
		TT_polymerize_vertical_minus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TD_polymerize_horizontal_plus:
	{
		TD_polymerize_horizontal_plus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TD_polymerize_horizontal_minus:
	{
		TD_polymerize_horizontal_minus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TD_polymerize_vertical_plus:
	{
		TD_polymerize_vertical_plus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_TD_polymerize_vertical_minus:
	{
		TD_polymerize_vertical_minus();
		rule_parameters.Z_ctp--;
		//calculate_Z_ctp();
	}
	break;
	case rule_depolymerize_plus_TT: 
	{
		depolymerizeTT_plus();
		calculate_Z_ctp();
	}
		break;
	case rule_depolymerize_plus_TD: 
	{
		depolymerizeTD_plus();
		calculate_Z_ctp();
	}
		break;
	case rule_depolymerize_plus_DD: 
	{
		depolymerizeDD_plus();
		calculate_Z_ctp();
	}
		break;
	case rule_depolymerize_minus_TT:
	{
		depolymerizeTT_minus();
		calculate_Z_ctp();
	}
		break;
	case rule_depolymerize_minus_TD:
	{
		depolymerizeTD_minus();
		calculate_Z_ctp();
	}
		break;
	case rule_depolymerize_minus_DD:
	{
		depolymerizeDD_minus();
		calculate_Z_ctp();
	}
		break;
	case rule_annealing_TT: 
	{
		annealing_TT();
	}
		break;
	case rule_annealing_TD:
	{
		annealing_TD();
	}
	
		break;
	case rule_fragmentate_TT: 
	{
		fragmentation_TT();
		calculate_Z_ctp();
	}
		break;
	case rule_fragmentate_TD:
	{
		fragmentation_TD();
		calculate_Z_ctp();
	}
		break;
	case rule_fragmentate_DD:
	{
		fragmentation_DD();
		calculate_Z_ctp();
	}
		break;
	case rule_bundling: 
	{
		bundling();
	}
		break;
	case rule_debundling: 
	{
		debundling();
		calculate_Z_ctp();
	}
		break;
	case rule_attatch:
	{
		anchor_attatch();
	}
		break;
	//case rule_attach: attatch();//replaced by anchor_attatch                     
		//break;
	case rule_anchoring: 
	{
		anchoring();
		calculate_Z_ctp();
	}
		break;
	case rule_deanchoring: 
	{
		deanchoring();
		calculate_Z_ctp();
	}
		break;
	case rule_hydrolysis_TTT:
	{
		Hydrolysis_TTT();
	}
		break;	
	case rule_hydrolysis_TTD:
	{
		Hydrolysis_TTD();
	}
		break;
	case rule_hydrolysis_DTD:
	{
		Hydrolysis_DTD();
	}
		break;
	case rule_empty_anchor_diffusion:
	{
		empty_anchor_diffusion();
	}
		break;
	case rule_bundling_inverse:
	{
		bundling_inverse();
	}
		break;
	case rule_debundling_inverse:
	{
		debundling_inverse();
		calculate_Z_ctp();
	}
		break;
	//case rule_rotation:rotate();
		//break;
	
	}
	//rotate();
	//calculate_Z_ctp();
	//check_polymerize_information();
	if (rule_parameters.step_num % 100 == 0)
	{
		//calculate_Z_ctp();
		rotate();
	}
	//check_polymerize_information();
	//clock_t start_time = clock();

	if (effect_rate == 1)
	{
		if (rule_parameters.TTT_sum_total + rule_parameters.TTD_sum_total + rule_parameters.DTD_sum_total != 0)
		{
			if ((rule_parameters.TTT_sum_total*rule_parameters.TTT_hydrolysys_rate + rule_parameters.TTD_sum_total*rule_parameters.TTD_hydrolysys_rate + rule_parameters.DTD_sum_total*rule_parameters.DTD_hydrolysys_rate) /
				(rule_parameters.TTT_sum_total + rule_parameters.TTD_sum_total + rule_parameters.DTD_sum_total) > rule_parameters.hydrolysis_rate)
			{
				rule_parameters.TTT_hydrolysys_rate = rule_parameters.TTT_hydrolysys_rate*0.99;
				rule_parameters.TTD_hydrolysys_rate = rule_parameters.TTT_hydrolysys_rate * rule_parameters.hydrolysis_multiplier;//1;//20;//20;// 100;//100;
				rule_parameters.DTD_hydrolysys_rate = rule_parameters.TTD_hydrolysys_rate * rule_parameters.hydrolysis_multiplier;//2;//40;//10;// 100;//100;
			}
			else
			{
				rule_parameters.TTT_hydrolysys_rate = rule_parameters.TTT_hydrolysys_rate*1.01;
				rule_parameters.TTD_hydrolysys_rate = rule_parameters.TTT_hydrolysys_rate * rule_parameters.hydrolysis_multiplier;//1;//20;//20;// 100;//100;
				rule_parameters.DTD_hydrolysys_rate = rule_parameters.TTD_hydrolysys_rate * rule_parameters.hydrolysis_multiplier;//2;//40;//10;// 100;//100;
			}
		}
		//if (rule_parameters.TT_depolymerize_minus_end_sum + rule_parameters.TD_depolymerize_minus_end_sum + rule_parameters.DD_depolymerize_minus_end_sum != 0)
		//{
		//	if ((rule_parameters.TT_depolymerize_minus_end_sum*rule_parameters.TT_depolymerize_minus_rate + rule_parameters.TD_depolymerize_minus_end_sum*rule_parameters.TD_depolymerize_minus_rate + rule_parameters.DD_depolymerize_minus_end_sum*rule_parameters.DD_depolymerize_minus_rate) /
		//		(rule_parameters.TT_depolymerize_minus_end_sum + rule_parameters.TD_depolymerize_minus_end_sum + rule_parameters.DD_depolymerize_minus_end_sum) > rule_parameters.depolymerize_rate)
		//	{
		//		rule_parameters.TT_depolymerize_minus_rate = rule_parameters.TT_depolymerize_minus_rate*0.99;
		//		rule_parameters.TT_depolymerize_plus_rate = rule_parameters.TT_depolymerize_minus_rate / 10;//1;//8;//80#1
		//		rule_parameters.TD_depolymerize_plus_rate = rule_parameters.TT_depolymerize_plus_rate * 10;//5;//20;//#5
		//		rule_parameters.DD_depolymerize_plus_rate = rule_parameters.TD_depolymerize_plus_rate * 10;//10;//40;//#5
		//		rule_parameters.TD_depolymerize_minus_rate = rule_parameters.TT_depolymerize_minus_rate * 10;//5;//20;//#5
		//		rule_parameters.DD_depolymerize_minus_rate = rule_parameters.TD_depolymerize_minus_rate * 10;//10;//40;//#5
		//	}
		//	else
		//	{
		//		rule_parameters.TT_depolymerize_minus_rate = rule_parameters.TT_depolymerize_minus_rate*1.01;
		//		rule_parameters.TT_depolymerize_plus_rate = rule_parameters.TT_depolymerize_minus_rate / 10;//1;//8;//80#1
		//		rule_parameters.TD_depolymerize_plus_rate = rule_parameters.TT_depolymerize_plus_rate * 10;//5;//20;//#5
		//		rule_parameters.DD_depolymerize_plus_rate = rule_parameters.TD_depolymerize_plus_rate * 10;//10;//40;//#5
		//		rule_parameters.TD_depolymerize_minus_rate = rule_parameters.TT_depolymerize_minus_rate * 10;//5;//20;//#5
		//		rule_parameters.DD_depolymerize_minus_rate = rule_parameters.TD_depolymerize_minus_rate * 10;//10;//40;//#5
		//	}

		//}
	}
	
	if (irule != 0 && irule != 18 && irule != 10 && irule != 11)
	{
		rule_parameters.event_output << rule_parameters.polymer_diffusion_propensity - temp_diffusion_propensity <<" "<<rule_parameters.polymer_diffusion_propensity<< endl;
	}






	rule_parameters.ps[rule_anchoring] = (rule_parameters.empty_anchor_unit_num - rule_parameters.emtpy_anchor_attatch_sum)*rule_parameters.Z_ctp*rule_parameters.anchoring_rate;
	rule_parameters.ps[rule_deanchoring] = (rule_parameters.Anchor_num - rule_parameters.empty_anchor_unit_num)*rule_parameters.deanchoring_rate;
	if (rule_parameters.MinC == 0)
	{
		rule_parameters.ps[rule_TT_polymerize_horizontal_plus] = rule_parameters.TT_polymerize_end_horizontal_plus_sum*rule_parameters.TT_polymerize_horizontal_plus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TT_polymerize_horizontal_minus] = rule_parameters.TT_polymerize_end_horizontal_minus_sum*rule_parameters.TT_polymerize_horizontal_minus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TT_polymerize_vertical_plus] = rule_parameters.TT_polymerize_end_vertical_plus_sum*rule_parameters.TT_polymerize_vertical_plus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TT_polymerize_vertical_minus] = rule_parameters.TT_polymerize_end_vertical_minus_sum*rule_parameters.TT_polymerize_vertical_minus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TD_polymerize_horizontal_plus] = rule_parameters.TD_polymerize_end_horizontal_plus_sum*rule_parameters.TD_polymerize_horizontal_plus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TD_polymerize_horizontal_minus] = rule_parameters.TD_polymerize_end_horizontal_minus_sum*rule_parameters.TD_polymerize_horizontal_minus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TD_polymerize_vertical_plus] = rule_parameters.TD_polymerize_end_vertical_plus_sum*rule_parameters.TD_polymerize_vertical_plus_rate*rule_parameters.Z_ctp;
		rule_parameters.ps[rule_TD_polymerize_vertical_minus] = rule_parameters.TD_polymerize_end_vertical_minus_sum*rule_parameters.TD_polymerize_vertical_minus_rate*rule_parameters.Z_ctp;
	}
	//rule_parameters.ps[rule_polymerize_vertical] = rule_parameters.polymerize_end_vertical_sum*polymerize_vertical_rate*rule_parameters.Z_ctp;
	rule_parameters.ps[rule_diffusion_polymer] = rule_parameters.polymer_diffusion_propensity*rule_parameters.diffusion_rate;
	rule_parameters.ps[rule_annealing_TT] = rule_parameters.TT_annealing_sum*rule_parameters.TT_annealing_rate;
	rule_parameters.ps[rule_annealing_TD] = rule_parameters.TD_annealing_sum*rule_parameters.TD_annealing_rate;
	rule_parameters.ps[rule_fragmentate_TT] = rule_parameters.TT_fragmentation_sum*rule_parameters.TT_fragmentation_rate;
	rule_parameters.ps[rule_fragmentate_TD] = rule_parameters.TD_fragmentation_sum*rule_parameters.TD_fragmentation_rate;
	rule_parameters.ps[rule_fragmentate_DD] = rule_parameters.DD_fragmentation_sum*rule_parameters.DD_fragmentation_rate;
	rule_parameters.ps[rule_deanchoring] = (rule_parameters.Anchor_num - rule_parameters.empty_anchor_unit_num)*rule_parameters.deanchoring_rate;
	rule_parameters.ps[rule_depolymerize_plus_TT] = rule_parameters.TT_depolymerize_plus_end_sum*rule_parameters.TT_depolymerize_plus_rate;
	rule_parameters.ps[rule_depolymerize_plus_TD] = rule_parameters.TD_depolymerize_plus_end_sum*rule_parameters.TD_depolymerize_plus_rate;
	rule_parameters.ps[rule_depolymerize_plus_DD] = rule_parameters.DD_depolymerize_plus_end_sum*rule_parameters.DD_depolymerize_plus_rate;
	rule_parameters.ps[rule_depolymerize_minus_TT] = rule_parameters.TT_depolymerize_minus_end_sum*rule_parameters.TT_depolymerize_minus_rate;
	rule_parameters.ps[rule_depolymerize_minus_TD] = rule_parameters.TD_depolymerize_minus_end_sum*rule_parameters.TD_depolymerize_minus_rate;
	rule_parameters.ps[rule_depolymerize_minus_DD] = rule_parameters.DD_depolymerize_minus_end_sum*rule_parameters.DD_depolymerize_minus_rate;
	rule_parameters.ps[rule_empty_anchor_diffusion]=rule_parameters.empty_anchor_diffusion_direction_sum*rule_parameters.emtpy_anchor_diffusion_rate;
	rule_parameters.ps[rule_bundling] = rule_parameters.bundling_sum*rule_parameters.bundling_rate;
	rule_parameters.ps[rule_debundling] =rule_parameters.debundling_sum*rule_parameters.debundling_rate;
	rule_parameters.ps[rule_attatch] = rule_parameters.emtpy_anchor_attatch_sum*rule_parameters.attatch_rate;
	rule_parameters.ps[rule_hydrolysis_TTT] = rule_parameters.TTT_sum_total*rule_parameters.TTT_hydrolysys_rate;
	rule_parameters.ps[rule_hydrolysis_TTD] = rule_parameters.TTD_sum_total*rule_parameters.TTD_hydrolysys_rate;
	rule_parameters.ps[rule_hydrolysis_DTD] = rule_parameters.DTD_sum_total*rule_parameters.DTD_hydrolysys_rate;
	rule_parameters.ps[rule_rotation] = 0;
	rule_parameters.ps[rule_bundling_inverse] = rule_parameters.bundling_inverse_sum*rule_parameters.bundling_inverse_rate;
	rule_parameters.ps[rule_debundling_inverse] = rule_parameters.debundling_inverse_sum*rule_parameters.debundling_inverse_rate;



	


	//rule_parameters.ps[rule_rotation] = rule_parameters.rotation_sum*rotation_rate;
	end = omp_get_wtime();
	//clock_t end_time = clock();
	rule_parameters.count[irule]++;
	rule_parameters.time[irule]= rule_parameters.time[irule]+(end-start);
	//cout << (double(end_time) - double(start_time)) << endl;
	check_sumps();
	


	rule_parameters.step_num++;
	
	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//default_random_engine gen;
	//uniform_real_distribution<> tempf(0, 1);//here is for the time part
	//double temp = 0;
	//do
	//{
	//	temp = tempf(gen);
	//} while (temp == 0);

	//double deltat = (1.0 / rule_parameters.sumps)*log(1.0 / temp);
	//rule_parameters.t = rule_parameters.t + deltat;


	//if (rule_parameters.step_num % 1000 == 0)
	//{
	//	//calculate_Z_ctp();
	//	rotate();
	//	//int double_sum;
	//	/*for (auto i = rule_structure.polymer_cluster_list.begin(); i != rule_structure.polymer_cluster_list.end(); i++)
	//	{
	//		if (i->polymer_sequence.size()>1)
	//		{
	//			double_sum = i->polymer_sequence.begin()->Hydrolysis_mark * 10 + (--i->polymer_sequence.end())->Hydrolysis_mark;
	//			switch (double_sum)
	//			{
	//				case TD:
	//				rule_parameters.TD_case++;
	//				break;
	//				case DT:
	//				rule_parameters.DT_case++;
	//				break;
	//				case TT:
	//				rule_parameters.TT_case++;
	//				break;
	//				case DD:
	//				rule_parameters.DD_case++;
	//				break;
	//			}
	//		}
	//	}*/
	//}
	




	
};

//void output(ofstream &fout)
//{
//	//std::cout << "lalala" << endl;
//	for (int i = 0; i<column_num; i++)
//	{
//		for (int j = 0; j < row_num; j++)
//		{
//			if (Grid[i][j].grid_anchor_pointer!=rule_structure.anchor_list.end())//&*Grid[i][j].grid_polymer_unit_pointer!=&*rule_structure.polymer_unit_list_end.begin())
//			{
//				fout << "1 ";
//			}
//			else
//			{
//				fout << "0 ";
//			}
//		}
//		fout << endl;
//	}
//	fout << endl;
//};

void output(ofstream &fout)
{
	//ofstream fout;
	//fout.open("testMatrix.txt");
	//fout << endl;
	//fout.width(2);
	//fout.setf(ios::left);
	//fout << "M";
	//fout << " ";
	//for (int j = 0; j < row_num; j++)
	//{
	//	fout.width(2);
	//	fout.setf(ios::left);
	//	fout <<j;
	//}
	//fout << endl;
	//fout << endl;
	//for (int i = 0; i<column_num; i++)
	//{
	//	fout.width(2);
	//	fout.setf(ios::left);
	//	fout << i;
	//	fout<<" ";
	//	for (int j = 0; j < row_num; j++)
	//	{
	//		if (&*Grid[i][j].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
	//		{//if there is a monomer there
	//			
	//			if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == T)
	//			{
	//				fout.width(2);
	//				fout.setf(ios::left);
	//				fout <<"T";
	//			}
	//			if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == D)
	//			{
	//				fout.width(2);
	//				fout.setf(ios::left);
	//				fout << "D";
	//			}

	//			
	//		}
	//		else
	//		{//if it is empty
	//			if (Grid[i][j].grid_anchor_pointer!=rule_structure.anchor_list.end())
	//			{//if there is a empty anchor
	//				fout.width(2);
	//				fout.setf(ios::left);
	//				fout << "O";
	//			}
	//			else
	//			{//if it is trully nothing
	//				fout.width(2);
	//				fout.setf(ios::left);
	//				fout << ".";
	//			}
	//			
	//		}
	//			
	//		
	//	}
	//	fout << endl;
	//}
	//fout << "Time: " << rule_parameters.t << endl;
	//system("pause");
	//fout.close();
	rule_parameters.T_sum = 0;
	rule_parameters.D_sum = 0;
	fout << rule_parameters.t << " ";
	int count = 0;
	int zero_count = 1;
	for (int i = 0; i < rule_parameters.column_num; i++)
	{
		for (int j = 0; j < rule_parameters.row_num; j++)
		{
			if (Grid[i][j].grid_anchor_pointer != rule_structure.anchor_list.end())
			{//if ther is a anchor
				//fout << A << " ";
				////rule_parameters.T_sum++;
				//count++;
				



				if (&*Grid[i][j].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{//no FtsZ above it, only anchor
					
					fout << zero_count<<" "<< AnchorBelowWithoutBound << Grid[i][j].grid_anchor_pointer->frap_mark <<NoFtsZ<< NoFtsZ << NoLateralBound<<NoDirections<<NotHeadOrTail<< " ";
					zero_count = 1;
					count++;

				}
				else
				{//there is FtsZ above it
					//first decide anchor bound
					list<Polymer_Unit>::iterator polymer_unit_iter = Grid[i][j].grid_polymer_unit_pointer;
					if (polymer_unit_iter->anchor_pointer == Grid[i][j].grid_anchor_pointer)
					{//FtsZ and anchor has a bound
						fout << zero_count << " " << AnchorBelowWithBound<< Grid[i][j].grid_anchor_pointer->frap_mark;
					}
					else
					{//FtsZ and anchor has no bound
						fout << zero_count << " " << AnchorBelowWithoutBound<< Grid[i][j].grid_anchor_pointer->frap_mark;
					}
					if (polymer_unit_iter->frap_mark == 1)
					{
						fout << 1;
					}
					if (polymer_unit_iter->frap_mark == 0)
					{
						fout << 0;
					}
					if (polymer_unit_iter->Hydrolysis_mark == T)
					{
						fout << T;
					}
					if (polymer_unit_iter->Hydrolysis_mark == D)
					{
						fout << D;
					}
					if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
					{
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							fout << LateralBoundDIrectionBoth;
						}
						else
						{
							fout << LateralBoundDirectionFirst;
						}
					}
					else
					{
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							fout << LateralBoundDirectionSecond;
						}
						else
						{
							fout << NoLateralBound;
						}
					}
					//directions

					if (polymer_unit_iter->polymer_cluster_pointer->direction == vertical)
					{
						if (polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
						{
							fout << direction_up;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
							

						}
						else
						{
							fout << direction_down;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
							
						}
					}
					else
					{//horizontal
						if (polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
						{
							fout << direction_right;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
						}
						else
						{
							fout << direction_left;
							
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
						}
					}

					
					
					


					
					


				}


			}
			else
			{//no anchor above
				if (&*Grid[i][j].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
				{


					//fout << 0 << " ";
					count++;
					zero_count++;

				}
				else
				{//there is FtsZ
					fout << zero_count << " " << NoAnchorBelow;
					list<Polymer_Unit>::iterator polymer_unit_iter = Grid[i][j].grid_polymer_unit_pointer;
					if (polymer_unit_iter->frap_mark == 1)
					{
						fout << 1;
					}
					if (polymer_unit_iter->frap_mark == 0)
					{
						fout << 0;
					}
					if (polymer_unit_iter->Hydrolysis_mark == T)
					{
						fout << T;
					}
					if (polymer_unit_iter->Hydrolysis_mark == D)
					{
						fout << D;
					}
					if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
					{
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							fout << LateralBoundDIrectionBoth;
						}
						else
						{
							fout << LateralBoundDirectionFirst;
						}
					}
					else
					{
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							fout << LateralBoundDirectionSecond;
						}
						else
						{
							fout << NoLateralBound;
						}
					}
					//directions

					if (polymer_unit_iter->polymer_cluster_pointer->direction == vertical)
					{
						if (polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
						{
							fout << direction_up;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}


						}
						else
						{
							fout << direction_down;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}

						}
					}
					else
					{//horizontal
						if (polymer_unit_iter->polymer_cluster_pointer->polarize == direction_second)
						{
							fout << direction_right;
							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
						}
						else
						{
							fout << direction_left;

							if (polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
							{
								fout << HeadAndTail << " ";
								zero_count = 1;
								count++;
							}
							else
							{
								if (polymer_unit_iter == polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.begin())
								{
									fout << Head << " ";
									zero_count = 1;
									count++;
								}
								else if (polymer_unit_iter == (--polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end()))
								{
									fout << Tail << " ";
									zero_count = 1;
									count++;
								}
								else
								{
									fout << NotHeadOrTail << " ";
									zero_count = 1;
									count++;
								}
							}
						}
					}
					
					
					


				}
			}
			

		}
	}
	fout << endl;
	if (count != rule_parameters.row_num*rule_parameters.column_num)
	{
		cout << "file writting error" << endl;
		system("pause");
	}
};



void output()
{
	//ofstream std::cout;
	//std::cout.open("testMatrix.txt");
	//std::cout << endl;
	std::cout.width(2);
	std::cout.setf(ios::left);
	std::cout << "M";
	std::cout << " ";
	for (int j = 0; j < rule_parameters.row_num; j++)
	{
		//std::cout << " ";
		std::cout.width(2);
		std::cout.setf(ios::left);
		std::cout << j;
		std::cout << "  ";
	}
	std::cout << endl;
	std::cout << endl;
	for (int i = rule_parameters.column_num-1; i>=0; i--)
	{
		//std::cout << " ";
		std::cout.width(2);
		std::cout.setf(ios::left);
		std::cout << i;
		std::cout << " ";
		for (int j = 0; j < rule_parameters.row_num; j++)
		{
			
			if (&*Grid[i][j].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
			{//if there is a monomer there

				if (Grid[i][j].grid_anchor_pointer != rule_structure.anchor_list.end())
				{
					if (&*Grid[i][j].grid_anchor_pointer->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
					{//if there is a anchor and they are not linked
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "@";
						
					}
					else
					{//if there is a anchor and they are linked
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "A";
						
					}


				}
				else
				{


					if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == T)
					{
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "T";
						
					}
					if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == D)
					{
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "D";
						
					}



				}

				//now std::cout the link symbol
				int tempmark = 1;
				if (Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
				{//verticle case
					if (Grid[i][j].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.bundling_polymer_unit_pair_list.end())
					{//bundleing
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << " ";
						tempmark = tempmark * 0;
					}
					if (Grid[i][j].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
					{//bundling
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "~";
						/*if (tempmark == 0)
						{
							std::cout << "error: polymer has bundling and debundling information on the same direction" << endl;
							system("pause");
						}*/
						tempmark = tempmark * 0;
					}
				}
				else
				{//horizontal//should not be activated for now
					if (Grid[i][j].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
					{//anealing
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << " ";
						tempmark = tempmark * 0;
					}
					if (j < rule_parameters.row_num - 1)
					{
						int tempj = j + 1;
						if (&*Grid[i][tempj].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
						{//if there is a monomer on the right side
							if (Grid[i][tempj].grid_polymer_unit_pointer->polymer_cluster_pointer == Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer&&
								Grid[i][tempj].grid_polymer_unit_pointer != (--Grid[i][tempj].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end()))
							{//if they belong to the same cluster
								std::cout.width(2);
								std::cout.setf(ios::left);
								std::cout << "-";
								tempmark = tempmark * 0;
							}
						}
					}
				}
				if (tempmark == 1)
				{
					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";
				}
					
				
				
				/*if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == -1)
				{
					std::cout << endl;
					std::cout << "error" << endl;
					std::cout << "col: " << i << " row: " << j << endl;
					system("pause");
				}*/
				/*if (Grid[i][j].grid_anchor_pointer != rule_structure.anchor_list.end() &&
					Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.size()==0)
				{
					std::cout << "error:anchor number does not meet with actual anchor number" << endl;
					system("pause");
				}*/
				/*if ((*Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin())->anchor_grid_pointer->col_position != (*Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.begin())->anchor_grid_pointer->col_position ||
					(*Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin())->anchor_grid_pointer->row_position != (*Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.begin())->anchor_grid_pointer->row_position)
				{
					std::cout << "error, bundle anchor information does not meet with grid information" << endl;
					system("pause");
				}*///multiple anchors for now
			/*	if (Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size()>10000)
				{
					std::cout << "error: invalid bundle with polymer_unit on grid" << endl;
					system("pause");
				}
				

				int col = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin()->polymer_grid_pointer->col_position;
				int row = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin()->polymer_grid_pointer->row_position;
				if (col != i && row != j)
				{
					std::cout << "error: this polymer_unit to a cluster that should not include this polymer_unit" << endl;
					system("pause");
				}*/

				/*int mark = 1;
				list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin();
				while (polymer_cluster_iter_iter != Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.end())
				{
					list<Polymer_Unit>::iterator polymer_unit_iter = (*polymer_cluster_iter_iter)->polymer_sequence.begin();
					while (polymer_unit_iter != (*polymer_cluster_iter_iter)->polymer_sequence.end())
					{
						col =polymer_unit_iter->polymer_grid_pointer->col_position;
						row = polymer_unit_iter->polymer_grid_pointer->row_position;
						if (col == i && row == j)
						{
							mark = mark * 0;
						}
						polymer_unit_iter++;
					}
					
					
					polymer_cluster_iter_iter++;
				}

				if (mark == 1)
				{
					std::cout << "error: this polymer_unit to a bundle that should not include this polymer_unit" << endl;
					system("pause");
				}*/
				
				//list<Polymer_Unit>::iterator polymer_unit_iter = Grid[i][j].grid_polymer_unit_pointer;
				//if (polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter != rule_structure.TT_anealing_polymer_unit_pair_list.end())
				//{//if it is indeed in one of the annealing pair
				//	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->TT_anealing_polymer_unit_pair_list_iter;
				//	int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
				//	int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
				//	if ((col_distance*col_distance + row_distance*row_distance) != 1)
				//	{	//it is no longer in the annealing position,need to delete this annealing relation
				//		std::cout << "Annealing relation error" << endl;
				//		system("pause");


				//	}
				//}
				//
				//if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.bundling_polymer_unit_pair_list.end())
				//{//if there is a bundle on its first direction, left or down
				//	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
				//	int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
				//	int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
				//	if ((col_distance*col_distance + row_distance*row_distance) != 1)
				//	{	//it is no longer in the annealing position,need to delete this annealing relation
				//		std::cout << "bundling relation error" << endl;
				//		system("pause");



				//	}
				//}

				//if (polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.bundling_polymer_unit_pair_list.end())
				//{//if there is a bundle on its second direction, right or up
				//	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_second];
				//	int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
				//	int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
				//	if ((col_distance*col_distance + row_distance*row_distance) != 1)
				//	{	//it is no longer in the annealing position,need to delete this annealing relation
				//		std::cout << "Bundling relation error" << endl;
				//		system("pause");



				//	}
				//}
				//

			}
			else
			{//if it is empty
				if (Grid[i][j].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a empty anchor
					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << "O";


					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";


				}
				else
				{//if it is trully nothing
					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";

					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";

				}

			}


		}
		std::cout << endl;
		std::cout.width(2);
		std::cout.setf(ios::left);
		std::cout << "";
		std::cout << " ";
		for (int j = 0; j < rule_parameters.row_num; j++)
		{
			if (&*Grid[i][j].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
			{
				int mark = 1;
				//now std::cout the link symbol
				if (Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == horizontal)
				{//horizontal case
					if (Grid[i][j].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.bundling_polymer_unit_pair_list.end())
					{//bundleing
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << " ";
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << ".";
						mark = mark * 0;
					}
					if (Grid[i][j].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
					{//bundleing
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << "S";
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << ".";
						mark = mark * 0;
					}
				}
				else
				{//verticle
					if (Grid[i][j].grid_polymer_unit_pointer->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
					{//anealing
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << " ";
						std::cout.width(2);
						std::cout.setf(ios::left);
						std::cout << ".";
						mark = mark * 0;
					}
					int tempi = i - 1;
					if (tempi < 0)
					{
						tempi = rule_parameters.column_num - 1;
					}

					if (&*Grid[tempi][j].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					{//if there is a monomer on the down side

						if (Grid[tempi][j].grid_polymer_unit_pointer->polymer_cluster_pointer->ring_mark == 1)
						{
							if (Grid[tempi][j].grid_polymer_unit_pointer->polymer_cluster_pointer == Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer)
							{//if they belong to the same cluster
								std::cout.width(2);
								std::cout.setf(ios::left);
								std::cout << "|";
								std::cout.width(2);
								std::cout.setf(ios::left);
								std::cout << ".";
								mark = mark * 0;
							}
						}
						else
						{
							if (Grid[tempi][j].grid_polymer_unit_pointer->polymer_cluster_pointer == Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer&&
								Grid[tempi][j].grid_polymer_unit_pointer != (--Grid[tempi][j].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.end()))
							{//if they belong to the same cluster
								std::cout.width(2);
								std::cout.setf(ios::left);
								std::cout << "|";
								std::cout.width(2);
								std::cout.setf(ios::left);
								std::cout << ".";
								mark = mark * 0;
							}
						}
						
					}
					
				}
				if (mark == 1)
				{
					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";
					std::cout.width(2);
					std::cout.setf(ios::left);
					std::cout << ".";
				}
			}
			else 
			{
				std::cout.width(2);
				std::cout.setf(ios::left);
				std::cout << ".";
				std::cout.width(2);
				std::cout.setf(ios::left);
				std::cout << ".";
			}
		}
		std::cout << endl;
	}
	std::cout << "Time: " << rule_parameters.t << endl;
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator TT_anealing_polymer_unit_pair_list_iter = rule_structure.TT_anealing_polymer_unit_pair_list.begin();
	//{
	//	map<tuple<int, int, int, int>,int> temp_map;
	//	/*vector<vector<int>> markmap_col;
	//	vector<vector<int>> markmap_row;
	//	markmap_col .resize(row_num);
	//	markmap_row.resize(row_num);
	//	for (int i = 0; i < row_num; i++)
	//	{
	//		markmap_col[i].resize(row_num);
	//		markmap_row[i].resize(row_num);
	//	}*/
	//	while (TT_anealing_polymer_unit_pair_list_iter != rule_structure.TT_anealing_polymer_unit_pair_list.end())
	//	{
	//		int col_distance = (TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position - TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
	//		int row_distance = (TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position - TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
	//		if ((col_distance*col_distance + row_distance*row_distance) != 1)
	//		{	//it is no longer in the annealing position,need to delete this annealing relation
	//			std::cout << "Annealing relation error" << endl;
	//			system("pause");
	//		}
	//		if (TT_anealing_polymer_unit_pair_list_iter->first->TT_anealing_polymer_unit_pair_list_iter[direction_second] != TT_anealing_polymer_unit_pair_list_iter)
	//		{//one of the pointer from polymer_unit to anealing pair is missing
	//			std::cout << "anealing_polymer_unit_pair_list_iter missing" << endl;
	//			system("pause");
	//		}
	//		if (TT_anealing_polymer_unit_pair_list_iter->second->TT_anealing_polymer_unit_pair_list_iter[direction_first] != TT_anealing_polymer_unit_pair_list_iter)
	//		{//one of the pointer from polymer_unit to anealing pair is missing
	//			std::cout << "anealing_polymer_unit_pair_list_iter missing" << endl;
	//			system("pause");
	//		}
	//		//if (markmap_col[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position] == 1 &&
	//		//	markmap_row[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position] == 1)
	//		//{//duplicate annealing information
	//		//	std::cout << "duplicate annealing information" << endl;
	//		//	system("pause");
	//		//}
	//		//markmap_col[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position] = 1;
	//		//markmap_row[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position] = 1;
	//		tuple<int, int, int, int> temp_tuple = make_tuple(TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position,
	//			TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position,
	//			TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position,
	//			TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position
	//			);
	//		if (temp_map.count(temp_tuple) > 0)
	//		{//duplicate annealing information
	//			std::cout << "duplicate annealing information" << endl;
	//			system("pause");
	//		}
	//		temp_map.insert(pair<tuple<int,int,int,int>,int>(temp_tuple,1));

	//		TT_anealing_polymer_unit_pair_list_iter++;
	//	}
	//}
	//
	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundling_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
	//while (bundling_pair_iter != rule_structure.bundling_polymer_unit_pair_list.end())
	//{
	//	int col_distance = (bundling_pair_iter->second->polymer_grid_pointer->col_position - bundling_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
	//	int row_distance = (bundling_pair_iter->second->polymer_grid_pointer->row_position - bundling_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
	//	if ((col_distance*col_distance + row_distance*row_distance) != 1)
	//	{	//it is no longer in the annealing position,need to delete this annealing relation
	//		std::cout << "Bundling relation error" << endl;
	//		system("pause");


	//	}
	//	bundling_pair_iter++;
	//}

	//list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundling_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
	//while (debundling_pair_iter != rule_structure.debundling_polymer_unit_pair_list.end())
	//{
	//	int col_distance = (debundling_pair_iter->second->polymer_grid_pointer->col_position - debundling_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
	//	int row_distance = (debundling_pair_iter->second->polymer_grid_pointer->row_position - debundling_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
	//	if ((col_distance*col_distance + row_distance*row_distance) != 1)
	//	{	//it is no longer in the annealing position,need to delete this annealing relation
	//		std::cout << "debundling relation error" << endl;
	//		system("pause");


	//	}
	//	debundling_pair_iter++;
	//}
	//vector<vector<int>> anchor_map(column_num, vector<int>(row_num,0));
	//list<Anchor_Unit>::iterator anchor_iter = rule_structure.anchor_list.begin();
	//int empty_anchor_sum = 0;
	//int attatch_anchor_sum = 0;
	//while (anchor_iter != rule_structure.anchor_list.end())
	//{
	//	if (Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
	//	{//there is anchor information on anchor list but not on the grid
	//		std::cout << "error: anchor information error" << endl;
	//		system("pause");
	//	}
	//	
	//	if (anchor_map[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position] == 0)
	//	{
	//		anchor_map[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position] = 1;
	//	}
	//	else
	//	{
	//		std::cout << "error: two anchor uses the same position" << endl;
	//		system("pause");
	//	}
	//	if (&*anchor_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
	//	{//if anchor is linked to a monomer
	//		if (anchor_iter->polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
	//		{//but this monomer is not linked to a ancher
	//			std::cout << "anchor attatch error" << endl;
	//			system("pause");
	//		}

	//	}
	//	if (&*anchor_iter->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
	//	{
	//		if (&*Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
	//		{
	//			attatch_anchor_sum++;
	//		}
	//		empty_anchor_sum++;
	//	}
	//		
	//	anchor_iter++;
	//}
	//if (attatch_anchor_sum != rule_parameters.emtpy_anchor_attatch_sum)
	//{
	//	std::cout << "attatch anchor sum doesn't match" << endl;
	//	system("pause");
	//}
	//if (empty_anchor_sum != rule_parameters.empty_anchor_unit_num)
	//{
	//	std::cout << "empty anchor sum doesn't match" << endl;
	//	system("pause");
	//}
	//
	//list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	//while(polymer_bundle_iter!=rule_structure.polymer_bundle_list.end())
	//{
	//	int temp_sum = 0;
	//	for (int i = 0; i < 4; i++)
	//	{
	//		temp_sum += polymer_bundle_iter->availabe_diffusion_direction[i];
	//	}
	//	if (temp_sum != polymer_bundle_iter->available_polymer_diffusion_direction_sum)
	//	{
	//		std::cout << "error: diffusion direction sum is not the same as sum up" << endl;
	//		system("pause");
	//	}
	//	if (polymer_bundle_iter->polymer_anchor_list.size()>1)
	//	{
	//		if (polymer_bundle_iter->available_polymer_diffusion_direction_sum != 0)
	//		{
	//			std::cout << "error: diffusion with more than 1 anchor" << endl;
	//			system("pause");
	//		}
	//	}
	//	polymer_bundle_iter++;
	//	
	//}

	//list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	//while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	//{
	//	
	//	if (polymer_cluster_iter->polymer_sequence.size()==0)
	//	{
	//		std::cout << "error: cluster size=0" << endl;
	//		system("pause");
	//	}

	//	if (polymer_cluster_iter->ring_mark == 1)
	//	{
	//		if (polymer_cluster_iter->TT_sum != column_num)
	//		{
	//			std::cout << "error: cluster in ring case, TT hydrolysis count error" << endl;
	//			system("pause");
	//		}
	//	}
	//	polymer_cluster_iter++;
	//}



	

	//system("pause");
	//std::cout.close();
};
void check()
{
	
	
	//for (int i = column_num - 1; i >= 0; i--)
	//{
	//	
	//	for (int j = 0; j < row_num; j++)
	//	{

	//		if (&*Grid[i][j].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
	//		{//if there is a monomer there
	//			if (Grid[i][j].grid_polymer_unit_pointer->attatch_mark == 1)
	//			{
	//				if (Grid[i][j].grid_anchor_pointer == rule_structure.anchor_list.end())
	//				{
	//					cout << "error: empty anchor with attatch_mark=1" << endl;
	//					system("pause");
	//				}
	//			}
	//			
	//			if (Grid[i][j].grid_anchor_pointer != rule_structure.anchor_list.end())
	//			{//if there is a anchor here
	//				if (&*Grid[i][j].grid_anchor_pointer->polymer_unit_pointer == &*Grid[i][j].grid_polymer_unit_pointer)
	//				{//anchor and polymer_unit does link to the other
	//					if (Grid[i][j].grid_anchor_pointer->attatch_mark == 1)
	//					{
	//						cout << "error: linked anchor with attatch_mark=1" << endl;
	//						system("pause");
	//					}
	//					
	//				}
	//				

	//			}
	//			//now std::cout the link symbol
	//			int tempmark = 1;
	//			if (Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->direction == vertical)
	//			{//verticle case
	//				if (Grid[i][j].grid_polymer_unit_pointer->bundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.bundling_polymer_unit_pair_list.end())
	//				{//bundleing
	//					
	//					tempmark = tempmark * 0;
	//				}
	//				if (Grid[i][j].grid_polymer_unit_pointer->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
	//				{//bundling
	//					 if (tempmark == 0)
	//					 {
	//						 std::cout << "error: polymer has bundling and debundling information on the same direction" << endl;
	//						 system("pause");
	//					 }
	//					tempmark = tempmark * 0;
	//				}
	//			}



	//			if (Grid[i][j].grid_polymer_unit_pointer->Hydrolysis_mark == -1)
	//			{
	//				std::cout << endl;
	//				std::cout << "error" << endl;
	//				std::cout << "col: " << i << " row: " << j << endl;
	//				system("pause");
	//			}
	//			
	//			if (Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size()>10000)
	//			{
	//				std::cout << "error: invalid bundle with polymer_unit on grid" << endl;
	//				system("pause");
	//			}


	//			int col = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin()->polymer_grid_pointer->col_position;
	//			int row = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.begin()->polymer_grid_pointer->row_position;
	//			if (col != i && row != j)
	//			{
	//				std::cout << "error: this polymer_unit to a cluster that should not include this polymer_unit" << endl;
	//				system("pause");
	//			}

	//			int mark = 1;
	//			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter = Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.begin();
	//			while (polymer_cluster_iter_iter != Grid[i][j].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.end())
	//			{
	//				list<Polymer_Unit>::iterator polymer_unit_iter = (*polymer_cluster_iter_iter)->polymer_sequence.begin();
	//				while (polymer_unit_iter != (*polymer_cluster_iter_iter)->polymer_sequence.end())
	//				{
	//					col = polymer_unit_iter->polymer_grid_pointer->col_position;
	//					row = polymer_unit_iter->polymer_grid_pointer->row_position;
	//					if (col == i && row == j)
	//					{
	//						mark = mark * 0;
	//					}
	//					polymer_unit_iter++;
	//				}


	//				polymer_cluster_iter_iter++;
	//			}

	//			if (mark == 1)
	//			{
	//				std::cout << "error: this polymer_unit to a bundle that should not include this polymer_unit" << endl;
	//				system("pause");
	//			}

	//			




	//		}
	//		else
	//		{//if there is no monomer
	//			
	//		}
	//		


	//	}
	//
	//	
	//}
	
	
	map<tuple<int, int, int, int>, int> temp_map;
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundling_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
	while (bundling_pair_iter != rule_structure.bundling_polymer_unit_pair_list.end())
	{
		if (&*bundling_pair_iter->first->polymer_grid_pointer->grid_polymer_unit_pointer != &*bundling_pair_iter->first)
		{
			cout << "bundle information error, bundle pair iter point to a monomer that does not point to this bundle pair iter" << endl;
			system("pause");
		}
		if (&*bundling_pair_iter->second->polymer_grid_pointer->grid_polymer_unit_pointer != &*bundling_pair_iter->second)
		{
			cout << "bundle information error, bundle pair iter point to a monomer that does not point to this bundle pair iter" << endl;
			system("pause");
		}
		int col_distance = (bundling_pair_iter->second->polymer_grid_pointer->col_position - bundling_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
		int row_distance = (bundling_pair_iter->second->polymer_grid_pointer->row_position - bundling_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
		if ((col_distance*col_distance + row_distance*row_distance) != 1)
		{	//it is no longer in the annealing position,need to delete this annealing relation
			std::cout << "Bundling relation error" << endl;
			system("pause");


		}
		tuple<int, int, int, int> temp_tuple = make_tuple(bundling_pair_iter->first->polymer_grid_pointer->col_position,
			bundling_pair_iter->first->polymer_grid_pointer->row_position,
			bundling_pair_iter->second->polymer_grid_pointer->col_position,
			bundling_pair_iter->second->polymer_grid_pointer->row_position
			);
		if (temp_map.count(temp_tuple) > 0)
		{//duplicate annealing information
			std::cout << "duplicate bundling information" << endl;
			system("pause");
		}
		temp_map.insert(pair<tuple<int, int, int, int>, int>(temp_tuple, 1));
		bundling_pair_iter++;
	}
	temp_map.clear();
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator debundling_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
	while (debundling_pair_iter != rule_structure.debundling_polymer_unit_pair_list.end())
	{
		if (debundling_pair_iter->first->polymer_grid_pointer->grid_polymer_unit_pointer != debundling_pair_iter->first)
		{
			cout << "debundle information error, debundle pair iter point to a monomer that does not point to this debundle pair iter" << endl;
			system("pause");
		}
		if (debundling_pair_iter->second->polymer_grid_pointer->grid_polymer_unit_pointer != debundling_pair_iter->second)
		{
			cout << "debundle information error, debundle pair iter point to a monomer that does not point to this debundle pair iter" << endl;
			system("pause");
		}
		int col_distance = (debundling_pair_iter->second->polymer_grid_pointer->col_position - debundling_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
		int row_distance = (debundling_pair_iter->second->polymer_grid_pointer->row_position - debundling_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
		if ((col_distance*col_distance + row_distance*row_distance) != 1)
		{	//it is no longer in the annealing position,need to delete this annealing relation
			std::cout << "debundling relation error" << endl;
			system("pause");


		}
		tuple<int, int, int, int> temp_tuple = make_tuple(debundling_pair_iter->first->polymer_grid_pointer->col_position,
			debundling_pair_iter->first->polymer_grid_pointer->row_position,
			debundling_pair_iter->second->polymer_grid_pointer->col_position,
			debundling_pair_iter->second->polymer_grid_pointer->row_position
			);
		if (temp_map.count(temp_tuple) > 0)
		{//duplicate annealing information
			std::cout << "duplicate debundling information" << endl;
			system("pause");
		}
		temp_map.insert(pair<tuple<int, int, int, int>, int>(temp_tuple, 1));
		debundling_pair_iter++;
	}
	vector<vector<int>> anchor_map(rule_parameters.column_num, vector<int>(rule_parameters.row_num, 0));
	list<Anchor_Unit>::iterator anchor_iter = rule_structure.anchor_list.begin();
	int empty_anchor_sum = 0;
	int attatch_anchor_sum = 0;
	int empty_anchor_direction_sum=0;
	while (anchor_iter != rule_structure.anchor_list.end())
	{
		if (Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
		{//there is anchor information on anchor list but not on the grid
			std::cout << "error: anchor information error" << endl;
			system("pause");
		}

		if (anchor_map[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position] == 0)
		{
			anchor_map[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position] = 1;
		}
		else
		{
			std::cout << "error: two anchor uses the same position" << endl;
			system("pause");
		}
		if (&*anchor_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{//if anchor is linked to a monomer
			if (anchor_iter->polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
			{//but this monomer is not linked to a ancher
				std::cout << "anchor attatch error" << endl;
				system("pause");
			}

		}
		if (&*anchor_iter->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
		{//if it is a empty anchor
			if (&*Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
			{
				attatch_anchor_sum++;
			}
			int col_position = anchor_iter->anchor_grid_pointer->col_position;
			int row_position = anchor_iter->anchor_grid_pointer->row_position;
			int up_col_position = col_position + 1;
			if (up_col_position > rule_parameters.column_num - 1)
				up_col_position = 0;
			int down_col_position = col_position - 1;
			if (down_col_position < 0)
				down_col_position = rule_parameters.column_num - 1;
			int right_row_position = row_position + 1;
			
			int left_row_position = row_position - 1;
			
			empty_anchor_sum++;

			if (Grid[up_col_position][row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{
				empty_anchor_direction_sum++;
			}
			else
			{
				if (anchor_iter->availabe_diffusion_direction[direction_up] == yes)
				{
					cout << "empty anchor diffusion direction error" << endl;
					system("pause");
				}
			}
			if (Grid[down_col_position][row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
			{
				empty_anchor_direction_sum++;
			}
			else
			{
				if (anchor_iter->availabe_diffusion_direction[direction_down] == yes)
				{
					cout << "empty anchor diffusion direction error" << endl;
					system("pause");
				}
			}
			if (left_row_position >= 0)
			{
				if (Grid[col_position][left_row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
				{
					empty_anchor_direction_sum++;
				}
			}
			else
			{
				if (anchor_iter->availabe_diffusion_direction[direction_left] == yes)
				{
					cout << "empty anchor diffusion direction error" << endl;
					system("pause");
				}
			}
			if (right_row_position < rule_parameters.row_num)
			{
				if (Grid[col_position][right_row_position].grid_anchor_pointer == rule_structure.anchor_list.end())
				{
					empty_anchor_direction_sum++;
				}
			}
			else
			{
				if (anchor_iter->availabe_diffusion_direction[direction_right] == yes)
				{
					cout << "empty anchor diffusion direction error" << endl;
					system("pause");
				}
			}
		
		}
		if (&*anchor_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{//if it is linked to a monomer
			if (anchor_iter->available_polymer_diffusion_direction_sum != 0)
			{
				cout << "non empty anchor with non-zero diffusion direction sum" << endl;
				system("pause");
			}

		}
		if (&*anchor_iter->anchor_grid_pointer->grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
		{
			if (anchor_iter->attatch_mark == 1)
			{
			
				cout << "error: empty anchor with attatch_mark=1" << endl;
				system("pause");
				
			}
		}



		anchor_iter++;
	}
	if (attatch_anchor_sum != rule_parameters.emtpy_anchor_attatch_sum)
	{
		std::cout << "attatch anchor sum doesn't match" << endl;
		system("pause");
	}
	if (empty_anchor_sum != rule_parameters.empty_anchor_unit_num)
	{
		std::cout << "empty anchor sum doesn't match" << endl;
		system("pause");
	}
	if (empty_anchor_direction_sum != rule_parameters.empty_anchor_diffusion_direction_sum)
	{
		std::cout << "empty anchor diffusion direction sum doesn't match :" << empty_anchor_direction_sum<<" "<< rule_parameters.empty_anchor_diffusion_direction_sum <<endl;
		system("pause");
	}
	double diffusion_direction_temp_sum = 0;
	int rotation_sum_temp = 0;
	list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	while (polymer_bundle_iter != rule_structure.polymer_bundle_list.end())
	{
		/*if (polymer_bundle_iter->polymer_anchor_list.size() > 1)
		{
			if (polymer_bundle_iter->available_polymer_diffusion_direction_sum != 0)
			{
				std::cout << "error: nonzero diffusion direction sum with more than one anchor" << endl;
				system("pause");
			}
		}*/
		
		/*if (polymer_bundle_iter->bundle_cluster_pointer_list.size() == 1)
		{
			if ((*polymer_bundle_iter->bundle_cluster_pointer_list.begin())->polymer_sequence.size() == 1)
			{
				rotation_sum_temp++;
			}
		}*/
		int temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			//(polymer_bundle_iter->availabe_diffusion_direction[i];
			diffusion_direction_temp_sum += calculate_diffusion_rate(polymer_bundle_iter,i);
		}
		
		/*if (temp_sum != polymer_bundle_iter->available_polymer_diffusion_direction_sum)
		{
			std::cout << "error: diffusion direction sum is not the same as sum up" << endl;
			system("pause");
		}*/
		/*check_polymer_bundle_diffusable_direction(polymer_bundle_iter,1);
		if (temp_sum != polymer_bundle_iter->available_polymer_diffusion_direction_sum)
		{
			std::cout << "error: diffusion direction sum on record is not the same as it actually is" << endl;
			system("pause");
		}*/
		/*if (polymer_bundle_iter->polymer_anchor_list.size()>1)
		{
			if (polymer_bundle_iter->available_polymer_diffusion_direction_sum != 0)
			{
				std::cout << "error: diffusion with more than 1 anchor" << endl;
				system("pause");
			}
		}*/
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_cluster_iter_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
		while (polymer_cluster_iter_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
		{
			if ((*polymer_cluster_iter_iter)->cluster_bundle_pointer != polymer_bundle_iter)
			{
				cout << "error: bundle_cluster_iter and cluster_bundle_iter error" << endl;
				system("pause");
			}
			polymer_cluster_iter_iter++;
		}

		polymer_bundle_iter++;

	}
	/*if (rotation_sum_temp != rule_parameters.rotation_sum)
	{
		cout << "rotation sum doesn't match" << endl;
		system("pause");
	}*/
	/*if (diffusion_direction_temp_sum != rule_parameters.polymer_diffusion_propensity)
	{
		cout << "polymer diffusion direction sum doesn't match: "<<diffusion_direction_temp_sum << endl;
		system("pause");
	}*/
	int endsum_horizontal = 0;
	int endsum_vertical = 0;
	int TT_endsum = 0;
	int TD_endsum = 0;
	int TT_fragmentation_sum = 0;
	int TD_fragmentation_sum = 0;
	int DD_fragmentation_sum = 0;
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{

		list<Polymer_Unit>::iterator polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
		while (polymer_unit_iter != polymer_cluster_iter->polymer_sequence.end())
		{
			if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
			{//if there is a bundle on its first direction, left or down
				list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first];
				int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
				int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
				if ((col_distance*col_distance + row_distance*row_distance) != 1)
				{
					cout << "bundling with distance larger than 1" << endl;
					system("pause");


					//}
				}
			}

				if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
				{//if there is a bundle on its second direction, right or up
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second];
					int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
					int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
					if ((col_distance*col_distance + row_distance*row_distance) != 1)
					{	//it is no longer in the annealing position,need to delete this annealing relation
						cout << "bundling with distance larger than 1" << endl;
						system("pause");
					}
				}
				polymer_unit_iter++;
		}

	






		if (polymer_cluster_iter->direction != vertical&&polymer_cluster_iter->direction != horizontal)
		{
			cout << "polymer_cluster_direction error" << endl;
			system("pause");
		}
		if (polymer_cluster_iter->ring_mark == true)
		{
			if (polymer_cluster_iter->polymer_sequence.size() != rule_parameters.column_num)
			{
				cout << "ring mark=1 with not ring case" << endl;
				system("pause");
			}
		}
		if (polymer_cluster_iter->direction == vertical)
		{
			if (polymer_cluster_iter->ring_mark == 1)
			{//if it is a ring
				/*if (polymer_cluster_iter->polymerize_end_sum != 0)
				{
					cout << "error: try to do polymerize in ring case" << endl;
					system("pause");
				}*/
			}
			else
			{
				{//consider down case
					int col_position = polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position;
					int row_position = polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position;;
					int up_col_position = col_position + 1;
					if (up_col_position > rule_parameters.column_num - 1)
						up_col_position = 0;
					int down_col_position = col_position - 1;
					if (down_col_position < 0)
						down_col_position = rule_parameters.column_num - 1;
					int right_row_position = row_position + 1;

					int left_row_position = row_position - 1;
					if (&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					{//if there is something in the lower end
						if (polymer_cluster_iter->availabe_polymerize_direction[direction_down] == 1)
						{
							cout << "error: try to do polymerize while there is already a polymer_unit there" << endl;
							system("pause");
						}
					}
				}
				{//consider up case
					int col_position = (--polymer_cluster_iter->polymer_sequence.end())->polymer_grid_pointer->col_position;
					int row_position = (--polymer_cluster_iter->polymer_sequence.end())->polymer_grid_pointer->row_position;;
					int up_col_position = col_position + 1;
					if (up_col_position > rule_parameters.column_num - 1)
						up_col_position = 0;
					int down_col_position = col_position - 1;
					if (down_col_position < 0)
						down_col_position = rule_parameters.column_num - 1;
					int right_row_position = row_position + 1;

					int left_row_position = row_position - 1;
					if (&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					{//if there is something in the lower end
						if (polymer_cluster_iter->availabe_polymerize_direction[direction_up] == 1)
						{
							cout << "error: try to do polymerize while there is already a polymer_unit there" << endl;
							system("pause");
						}
					}
				}
				
			}
		}
		if (polymer_cluster_iter->availabe_polymerize_direction.size() == 0)
		{
			std::cout << "error: polymerize information error" << endl;
			system("pause");
		}
		/*if (polymer_cluster_iter->polymerize_end_sum < 0)
		{
			std::cout << "error: negative polymerize end" << endl;
			system("pause");
		}*/
		if (polymer_cluster_iter->polymer_sequence.size() == 0)
		{
			std::cout << "error: cluster size=0" << endl;
			system("pause");
		}

		/*if (polymer_cluster_iter->ring_mark == 1)
		{
			if (polymer_cluster_iter->TT_sum != column_num)
			{
				std::cout << "error: cluster in ring case, TT hydrolysis count error" << endl;
				system("pause");
			}
		}*/
		int temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			temp_sum += polymer_cluster_iter->availabe_polymerize_direction[i];
		}
		/*if (temp_sum != polymer_cluster_iter->polymerize_end_sum)
		{
			std::cout << "error: polymer end sum mismatch" << endl;
			system("pause");
		}*/
		/*if (polymer_cluster_iter->TT_depolymerize_end_sum > 2)
		{
			std::cout << "error:depolymerize_end_sum error" << endl;
			system("pause");

		}
		if (polymer_cluster_iter->TD_depolymerize_end_sum > 2)
		{
			std::cout << "error:depolymerize_end_sum error" << endl;
			system("pause");

		}
		if (polymer_cluster_iter->DD_depolymerize_end_sum > 2)
		{
			std::cout << "error:depolymerize_end_sum error" << endl;
			system("pause");

		}*/
		if (polymer_cluster_iter->direction == horizontal)
		{
			endsum_horizontal += temp_sum;
		}
		if (polymer_cluster_iter->direction == vertical)
		{
			endsum_vertical += temp_sum;
		}
		
		//TT_endsum += polymer_cluster_iter->TT_depolymerize_end_sum;
		//TD_endsum += polymer_cluster_iter->TD_depolymerize_end_sum;
		TT_fragmentation_sum += polymer_cluster_iter->TT_sum;
		TD_fragmentation_sum += polymer_cluster_iter->TD_sum;
		DD_fragmentation_sum += polymer_cluster_iter->DD_sum;
		polymer_cluster_iter++;
	}
	/*if (endsum_vertical != rule_parameters.polymerize_end_vertical_sum)
	{
		cout << "polymer polymerize end sum doesn't match" << endl;
		system("pause");
	}
	if (endsum_horizontal != rule_parameters.polymerize_end_horizontal_sum)
	{
		cout << "polymer polymerize end sum doesn't match" << endl;
		system("pause");
	}*/
	/*if (TT_endsum != rule_parameters.TT_depolymerize_end_sum)
	{
		cout << "polymer TT_depolymerize end sum doesn't match" << endl;
		system("pause");
	}
	if (TD_endsum != rule_parameters.TD_depolymerize_end_sum)
	{
		cout << "polymer TD_depolymerize end sum doesn't match" << endl;
		system("pause");
	}*/
	if (rule_parameters.empty_anchor_diffusion_direction_sum < 0)
	{
		std::cout << "error: empty_anchor_diffusion_direction_sum can't be negative" << endl;
		system("pause");
	}
	if (rule_parameters.TT_fragmentation_sum !=TT_fragmentation_sum)
	{
		std::cout << "TT count error" << endl;
		system("pause");
	}
	if (rule_parameters.TD_fragmentation_sum != TD_fragmentation_sum)
	{
		std::cout << "TD count error" << endl;
		system("pause");
	}
	if (rule_parameters.DD_fragmentation_sum != DD_fragmentation_sum)
	{
		std::cout << "DD count error" << endl;
		system("pause");
	}


	//system("pause");
	//std::cout.close();
};

void check_direction()
{
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{
		if (polymer_cluster_iter->direction == vertical)
		{

		}
		else if (polymer_cluster_iter->direction == horizontal)
		{

		}
		else
		{
			std::cout << "polymer_cluster direction is not vertical or horizontal" << endl;
			system("pause");
		}
		polymer_cluster_iter++;
	}
}


void check_sumps()
{
	double temp=0;
	for (int i = 0; i < num_of_rules; i++)
	{
		temp = temp + rule_parameters.ps[i];
		if (rule_parameters.ps[i] < 0)
		{
			if (abs(rule_parameters.ps[i]) < 0.000000001)
			{
				rule_parameters.ps[i] = 0;
			}
				
			else
			{
				std::cout << "error: negative reaction parameter: " << rule_parameters.reaction_mark<<" "<<i<< endl;
				system("pause");
			}
			
		}
	}
	rule_parameters.sumps = temp;
};
void calculate_Z_ctp()
{
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	rule_parameters.cluster_bound = 0;
	int tempZ = rule_parameters.Z;
	for (int i = 0; i < rule_structure.polymer_cluster_list.size(); i++)
	{
		tempZ= tempZ - polymer_cluster_iter->polymer_sequence.size();
		rule_parameters.cluster_bound += polymer_cluster_iter->polymer_sequence.size() - 1;
		polymer_cluster_iter++;
	}
	rule_parameters.Z_ctp = tempZ;
};

//int hextodec(char *str)
//{
//	int i;
//	sscanf(str, "%x", &i);
//	return i;
//}
void test_annealing()
{
	int TT_sum_temp = 0;
	int TD_sum_temp = 0;
	
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator anealing_polymer_unit_pair_list_iter = rule_structure.anealing_polymer_unit_pair_list.begin();
	{
		map<tuple<int, int, int, int>, int> temp_map;
		/*vector<vector<int>> markmap_col;
		vector<vector<int>> markmap_row;
		markmap_col .resize(row_num);
		markmap_row.resize(row_num);
		for (int i = 0; i < row_num; i++)
		{
		markmap_col[i].resize(row_num);
		markmap_row[i].resize(row_num);
		}*/
		while (anealing_polymer_unit_pair_list_iter != rule_structure.anealing_polymer_unit_pair_list.end())
		{
			if (anealing_polymer_unit_pair_list_iter->first->Hydrolysis_mark == T && anealing_polymer_unit_pair_list_iter->second->Hydrolysis_mark == T)
			{
				TT_sum_temp++;
			}
			else if (anealing_polymer_unit_pair_list_iter->first->Hydrolysis_mark == D && anealing_polymer_unit_pair_list_iter->second->Hydrolysis_mark == D)
			{

			}
			else
			{
				TD_sum_temp++;
			}
			if (anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->grid_polymer_unit_pointer != anealing_polymer_unit_pair_list_iter->first)
			{
				cout << "TT_anealing_polymer_unit_pair_list_iter information error, TT_anealing_polymer_unit_pair_list_iter point to a monomer that does not point to this TT_anealing_polymer_unit_pair_list_iter" << endl;
				system("pause");
			}
			if (anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->grid_polymer_unit_pointer != anealing_polymer_unit_pair_list_iter->second)
			{
				cout << "TT_anealing_polymer_unit_pair_list_iter information error, TT_anealing_polymer_unit_pair_list_iter point to a monomer that does not point to this TT_anealing_polymer_unit_pair_list_iter" << endl;
				system("pause");
			}
			int col_distance = (anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position - anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position + rule_parameters.column_num) % rule_parameters.column_num;
			int row_distance = (anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position - anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position + rule_parameters.row_num) % rule_parameters.row_num;
			if ((col_distance*col_distance + row_distance*row_distance) != 1)
			{	//it is no longer in the annealing position,need to delete this annealing relation
				std::cout << "Annealing relation error" << endl;
				system("pause");
			}
			if (anealing_polymer_unit_pair_list_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] != anealing_polymer_unit_pair_list_iter)
			{//one of the pointer from polymer_unit to anealing pair is missing
				std::cout << "anealing_polymer_unit_pair_list_iter missing case 1" << endl;
				system("pause");
			}
			if (anealing_polymer_unit_pair_list_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] != anealing_polymer_unit_pair_list_iter)
			{//one of the pointer from polymer_unit to anealing pair is missing
				std::cout << "anealing_polymer_unit_pair_list_iter missing case 2" << endl;
				system("pause");
			}
			//if (markmap_col[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position] == 1 &&
			//	markmap_row[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position] == 1)
			//{//duplicate annealing information
			//	std::cout << "duplicate annealing information" << endl;
			//	system("pause");
			//}
			//markmap_col[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position] = 1;
			//markmap_row[TT_anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position][TT_anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position] = 1;
			tuple<int, int, int, int> temp_tuple = make_tuple(anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->col_position,
				anealing_polymer_unit_pair_list_iter->first->polymer_grid_pointer->row_position,
				anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->col_position,
				anealing_polymer_unit_pair_list_iter->second->polymer_grid_pointer->row_position
				);
			if (temp_map.count(temp_tuple) > 0)
			{//duplicate annealing information
				std::cout << "duplicate annealing information" << endl;
				system("pause");
			}
			temp_map.insert(pair<tuple<int, int, int, int>, int>(temp_tuple, 1));

			anealing_polymer_unit_pair_list_iter++;
		}
		int col = 0; int row = 0;
		for (auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin(); polymer_cluster_iter != rule_structure.polymer_cluster_list.end(); polymer_cluster_iter++)
		{
			auto polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
			col = polymer_unit_iter->polymer_grid_pointer->col_position;
			row = polymer_unit_iter->polymer_grid_pointer->row_position;
			int up_col_position = col + 1;
			if (up_col_position > rule_parameters.column_num - 1)
				up_col_position = 0;
			int down_col_position = col - 1;
			if (down_col_position < 0)
				down_col_position = rule_parameters.column_num - 1;
			int right_row_position = row + 1;
			if (right_row_position > rule_parameters.row_num)
			{
				cout << "out of right boundary" << endl;
				system("pause");
			}
			int left_row_position = row - 1;
			if (left_row_position < -1)
			{
				cout << "out of left boundary" << endl;
				system("pause");
			}
			if (Grid[down_col_position][row].grid_polymer_unit_pointer != rule_structure.polymer_unit_list_end.begin())
			{
				auto down_polymer_unit_iter = Grid[down_col_position][row].grid_polymer_unit_pointer;
				if (polymer_unit_iter->polymer_cluster_pointer->direction == vertical&&down_polymer_unit_iter->polymer_cluster_pointer->direction == vertical&&
					polymer_unit_iter->polymer_cluster_pointer->polarize == down_polymer_unit_iter->polymer_cluster_pointer->polarize)
				{
					if (down_polymer_unit_iter != --down_polymer_unit_iter->polymer_cluster_pointer->polymer_sequence.end())
					{
						cout << "polymer unit information error" << endl;
						system("pause");
					}
					//it is in the annealing case, make temp tuple
					if (polymer_unit_iter->polymer_cluster_pointer == down_polymer_unit_iter->polymer_cluster_pointer)
					{
						std::cout << "ring case" << endl;
						system("pause");
					}
					tuple<int, int, int, int> temp_tuple = make_tuple(down_polymer_unit_iter->polymer_grid_pointer->col_position,
						down_polymer_unit_iter->polymer_grid_pointer->row_position,
						polymer_unit_iter->polymer_grid_pointer->col_position,
						polymer_unit_iter->polymer_grid_pointer->row_position
						);
					if (temp_map.count(temp_tuple) != 1)
					{//duplicate annealing information
						std::cout << "annealing information missing" << endl;
						system("pause");
					}
				}
			}
		}
		if (TT_sum_temp != rule_parameters.TT_annealing_sum || TD_sum_temp != rule_parameters.TD_annealing_sum)
		{
			cout << "annealing information error" << endl;
			system("pause");
		}
	}
	
}

void check_bundling()
{
	int bundling_sum = 0;
	int debundling_sum = 0;
	int bundling_inverse_sum = 0;
	int debundling_inverse_sum = 0;
	for (auto bundle_pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin(); bundle_pair_iter != rule_structure.bundling_polymer_unit_pair_list.end(); bundle_pair_iter++)
	{
		if (bundle_pair_iter->first->polymer_cluster_pointer->polarize == bundle_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			bundling_sum ++;
		}
		else
		{
			bundling_inverse_sum++;
		}
	}
	for (auto debundle_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin(); debundle_pair_iter != rule_structure.debundling_polymer_unit_pair_list.end(); debundle_pair_iter++)
	{
		if (debundle_pair_iter->first->polymer_cluster_pointer->polarize == debundle_pair_iter->second->polymer_cluster_pointer->polarize)
		{
			debundling_sum++;
		}
		else
		{
			debundling_inverse_sum++;
		}
	}
	if (bundling_sum != rule_parameters.bundling_sum)
	{
		cout << "bundling infromation error" << endl;
		system("pause");
	}
	if (bundling_inverse_sum != rule_parameters.bundling_inverse_sum)
	{
		cout << "bundling inverse infromation error" << endl;
		system("pause");
	}
	if (debundling_sum != rule_parameters.debundling_sum)
	{
		cout << "debundling infromation error" << endl;
		system("pause");
	}
	if (debundling_inverse_sum != rule_parameters.debundling_inverse_sum)
	{
		cout << "debundling inverse infromation error" << endl;
		system("pause");
	}
}
void check_debundling_information()
{
	int double_sum=0;
	rule_parameters.debundle_pair_DD = 0;
	rule_parameters.debundle_pair_TD = 0;
	rule_parameters.debundle_pair_TT = 0;
	
	for (auto debundle_pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin(); debundle_pair_iter != rule_structure.debundling_polymer_unit_pair_list.end(); debundle_pair_iter++)
	{
		double_sum = debundle_pair_iter->first->Hydrolysis_mark * 10 + debundle_pair_iter->second->Hydrolysis_mark;
		if (double_sum == TT)
		{
			rule_parameters.debundle_pair_TT++;
		}
		else if (double_sum == TD||double_sum==DT)
		{
			rule_parameters.debundle_pair_TD++;
		}
		else if (double_sum == DD)
		{
			rule_parameters.debundle_pair_DD++;
		}
	}
	
}

void check_polymerize_information()
{
	int temp_TT_polymerize_end_vertical_plus_sum = 0;
	int temp_TT_polymerize_end_vertical_minus_sum = 0;
	int temp_TT_polymerize_end_horizontal_plus_sum = 0;
	int temp_TT_polymerize_end_horizontal_minus_sum = 0;
	int temp_TD_polymerize_end_vertical_plus_sum = 0;
	int temp_TD_polymerize_end_vertical_minus_sum = 0;
	int temp_TD_polymerize_end_horizontal_plus_sum = 0;
	int temp_TD_polymerize_end_horizontal_minus_sum = 0;
	int TT_polymerize_plus_end_sum = 0;//directions
	int TT_polymerize_minus_end_sum = 0;//directions
	int TD_polymerize_plus_end_sum = 0;//directions
	int TD_polymerize_minus_end_sum = 0;//directions
	int TT_depolymerize_plus_end_sum = 0;
	int TD_depolymerize_plus_end_sum = 0;//TD and DT are the same here
	int DD_depolymerize_plus_end_sum = 0;
	int TT_depolymerize_minus_end_sum = 0;
	int TD_depolymerize_minus_end_sum = 0;//TD and DT are the same here
	int DD_depolymerize_minus_end_sum = 0;
	int TT_depolymerize_plus_end_sum_total = 0;
	int TD_depolymerize_plus_end_sum_total = 0;//TD and DT are the same here
	int DD_depolymerize_plus_end_sum_total = 0;
	int TT_depolymerize_minus_end_sum_total = 0;
	int TD_depolymerize_minus_end_sum_total = 0;//TD and DT are the same here
	int DD_depolymerize_minus_end_sum_total = 0;


	for (auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin(); polymer_cluster_iter != rule_structure.polymer_cluster_list.end(); polymer_cluster_iter++)
	{	
		TT_polymerize_plus_end_sum = 0;//directions
		TT_polymerize_minus_end_sum = 0;//directions
		TD_polymerize_plus_end_sum = 0;//directions
		TD_polymerize_minus_end_sum = 0;//directions
		TT_depolymerize_plus_end_sum = 0;
		TD_depolymerize_plus_end_sum = 0;//TD and DT are the same here
		DD_depolymerize_plus_end_sum = 0;
		TT_depolymerize_minus_end_sum = 0;
		TD_depolymerize_minus_end_sum = 0;//TD and DT are the same here
		DD_depolymerize_minus_end_sum = 0;
		if (polymer_cluster_iter->polarize == direction_first)
		{//down or left//saves as left to right
			switch (polymer_cluster_iter->polymer_sequence.size())
			{
				case 0:
				{
					break;
				}
				case 1:
				{
					break;
				}
				case 2:
				{
					auto polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
					auto polymer_unit_iter_back = ++polymer_cluster_iter->polymer_sequence.begin();
					int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
					switch (doublesum)
					{
						case TT:
						{
							TT_depolymerize_plus_end_sum++;
							TT_depolymerize_minus_end_sum++;
						}
							break;
						case TD:
						{
							TD_depolymerize_plus_end_sum++;
							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DT:
						{
							TD_depolymerize_plus_end_sum++;
							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DD:
						{
							DD_depolymerize_plus_end_sum++;
							DD_depolymerize_minus_end_sum++;
						}
						break;
					}
					
					break;
				}
				default:
				{
					auto polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
					auto polymer_unit_iter_back = ++polymer_cluster_iter->polymer_sequence.begin();

					if (polymer_cluster_iter->polarize == direction_first)
					{
						int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
						switch (doublesum)
						{
						case TT:
						{
							TT_depolymerize_plus_end_sum++;
							
						}
						break;
						case TD:
						{
							TD_depolymerize_plus_end_sum++;
							
						}
						break;
						case DT:
						{
							TD_depolymerize_plus_end_sum++;
							
						}
						break;
						case DD:
						{
							DD_depolymerize_plus_end_sum++;
							
						}
						break;
						}
					}
					else
					{
						int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
						switch (doublesum)
						{
						case TT:
						{
							
							TT_depolymerize_minus_end_sum++;
						}
						break;
						case TD:
						{
							
							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DT:
						{
							
							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DD:
						{
							
							DD_depolymerize_minus_end_sum++;
						}
						break;
						}
					}
					polymer_unit_iter_back = --polymer_cluster_iter->polymer_sequence.end();
					polymer_unit_iter = --(--polymer_cluster_iter->polymer_sequence.end());

					if (polymer_cluster_iter->polarize == direction_second)
					{
						int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
						switch (doublesum)
						{
						case TT:
						{
							TT_depolymerize_plus_end_sum++;

						}
						break;
						case TD:
						{
							TD_depolymerize_plus_end_sum++;

						}
						break;
						case DT:
						{
							TD_depolymerize_plus_end_sum++;

						}
						break;
						case DD:
						{
							DD_depolymerize_plus_end_sum++;

						}
						break;
						}
					}
					else
					{
						int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
						switch (doublesum)
						{
						case TT:
						{

							TT_depolymerize_minus_end_sum++;
						}
						break;
						case TD:
						{

							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DT:
						{

							TD_depolymerize_minus_end_sum++;
						}
						break;
						case DD:
						{

							DD_depolymerize_minus_end_sum++;
						}
						break;
						}
					}

				}
			}
			if (TT_depolymerize_plus_end_sum != polymer_cluster_iter->TT_depolymerize_plus_end_sum)
			{
				cout << "TT_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TT_depolymerize_minus_end_sum != polymer_cluster_iter->TT_depolymerize_minus_end_sum)
			{
				cout << "TT_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			if (TD_depolymerize_plus_end_sum != polymer_cluster_iter->TD_depolymerize_plus_end_sum)
			{
				cout << "TD_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TD_depolymerize_minus_end_sum != polymer_cluster_iter->TD_depolymerize_minus_end_sum)
			{
				cout << "TD_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			if (DD_depolymerize_plus_end_sum != polymer_cluster_iter->DD_depolymerize_plus_end_sum)
			{
				cout << "DD_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (DD_depolymerize_minus_end_sum != polymer_cluster_iter->DD_depolymerize_minus_end_sum)
			{
				cout << "DD_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			TT_depolymerize_plus_end_sum_total += TT_depolymerize_plus_end_sum;
			TT_depolymerize_minus_end_sum_total += TT_depolymerize_minus_end_sum;
			TD_depolymerize_plus_end_sum_total += TD_depolymerize_plus_end_sum;
			TD_depolymerize_minus_end_sum_total += TD_depolymerize_minus_end_sum;
			DD_depolymerize_plus_end_sum_total += DD_depolymerize_plus_end_sum;
			DD_depolymerize_minus_end_sum_total += DD_depolymerize_minus_end_sum;
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_down] == true)
			{
				if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
				{
					TT_polymerize_plus_end_sum++;
					temp_TT_polymerize_end_vertical_plus_sum++;
				}
				else
				{
					TD_polymerize_plus_end_sum++;
					temp_TD_polymerize_end_vertical_plus_sum++;
				}
			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_left] == true)
			{
				if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
				{
					TT_polymerize_plus_end_sum++;
					temp_TT_polymerize_end_horizontal_plus_sum++;
				}
				else
				{
					TD_polymerize_plus_end_sum++;
					temp_TD_polymerize_end_horizontal_plus_sum++;
				}
			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_up] == true)
			{
				if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
				{
					TT_polymerize_minus_end_sum++;
					temp_TT_polymerize_end_vertical_minus_sum++;
				}
				else
				{
					TD_polymerize_minus_end_sum++;
					temp_TD_polymerize_end_vertical_minus_sum++;
				}

			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_right] == true)
			{
				if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
				{
					TT_polymerize_minus_end_sum++;
					temp_TT_polymerize_end_horizontal_minus_sum++;
				}
				else
				{
					TD_polymerize_minus_end_sum++;
					temp_TD_polymerize_end_horizontal_minus_sum++;
				}

			}
		}
		else
		{
			switch (polymer_cluster_iter->polymer_sequence.size())
			{
			case 0:
			{
				break;
			}
			case 1:
			{
				break;
			}
			case 2:
			{
				auto polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				auto polymer_unit_iter_back = ++polymer_cluster_iter->polymer_sequence.begin();
				int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
				switch (doublesum)
				{
				case TT:
				{
					TT_depolymerize_plus_end_sum++;
					TT_depolymerize_minus_end_sum++;
				}
				break;
				case TD:
				{
					TD_depolymerize_plus_end_sum++;
					TD_depolymerize_minus_end_sum++;
				}
				break;
				case DT:
				{
					TD_depolymerize_plus_end_sum++;
					TD_depolymerize_minus_end_sum++;
				}
				break;
				case DD:
				{
					DD_depolymerize_plus_end_sum++;
					DD_depolymerize_minus_end_sum++;
				}
				break;
				}

				break;
			}
			default:
			{
				auto polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
				auto polymer_unit_iter_back = ++polymer_cluster_iter->polymer_sequence.begin();

				if (polymer_cluster_iter->polarize == direction_first)
				{
					int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
					switch (doublesum)
					{
					case TT:
					{
						TT_depolymerize_plus_end_sum++;

					}
					break;
					case TD:
					{
						TD_depolymerize_plus_end_sum++;

					}
					break;
					case DT:
					{
						TD_depolymerize_plus_end_sum++;

					}
					break;
					case DD:
					{
						DD_depolymerize_plus_end_sum++;

					}
					break;
					}
				}
				else
				{
					int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
					switch (doublesum)
					{
					case TT:
					{

						TT_depolymerize_minus_end_sum++;
					}
					break;
					case TD:
					{

						TD_depolymerize_minus_end_sum++;
					}
					break;
					case DT:
					{

						TD_depolymerize_minus_end_sum++;
					}
					break;
					case DD:
					{

						DD_depolymerize_minus_end_sum++;
					}
					break;
					}
				}
				polymer_unit_iter_back = --polymer_cluster_iter->polymer_sequence.end();
				polymer_unit_iter = --(--polymer_cluster_iter->polymer_sequence.end());

				if (polymer_cluster_iter->polarize == direction_second)
				{
					int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
					switch (doublesum)
					{
					case TT:
					{
						TT_depolymerize_plus_end_sum++;

					}
					break;
					case TD:
					{
						TD_depolymerize_plus_end_sum++;

					}
					break;
					case DT:
					{
						TD_depolymerize_plus_end_sum++;

					}
					break;
					case DD:
					{
						DD_depolymerize_plus_end_sum++;

					}
					break;
					}
				}
				else
				{
					int doublesum = polymer_unit_iter->Hydrolysis_mark * 10 + polymer_unit_iter_back->Hydrolysis_mark;
					switch (doublesum)
					{
					case TT:
					{

						TT_depolymerize_minus_end_sum++;
					}
					break;
					case TD:
					{

						TD_depolymerize_minus_end_sum++;
					}
					break;
					case DT:
					{

						TD_depolymerize_minus_end_sum++;
					}
					break;
					case DD:
					{

						DD_depolymerize_minus_end_sum++;
					}
					break;
					}
				}

			}
			}
			if (TT_depolymerize_plus_end_sum != polymer_cluster_iter->TT_depolymerize_plus_end_sum)
			{
				cout << "TT_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TT_depolymerize_minus_end_sum != polymer_cluster_iter->TT_depolymerize_minus_end_sum)
			{
				cout << "TT_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			if (TD_depolymerize_plus_end_sum != polymer_cluster_iter->TD_depolymerize_plus_end_sum)
			{
				cout << "TD_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TD_depolymerize_minus_end_sum != polymer_cluster_iter->TD_depolymerize_minus_end_sum)
			{
				cout << "TD_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			if (DD_depolymerize_plus_end_sum != polymer_cluster_iter->DD_depolymerize_plus_end_sum)
			{
				cout << "DD_depolymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (DD_depolymerize_minus_end_sum != polymer_cluster_iter->DD_depolymerize_minus_end_sum)
			{
				cout << "DD_depolymerize_minus_end_sum" << endl;
				system("pause");
			}
			TT_depolymerize_plus_end_sum_total += TT_depolymerize_plus_end_sum;
			TT_depolymerize_minus_end_sum_total += TT_depolymerize_minus_end_sum;
			TD_depolymerize_plus_end_sum_total += TD_depolymerize_plus_end_sum;
			TD_depolymerize_minus_end_sum_total += TD_depolymerize_minus_end_sum;
			DD_depolymerize_plus_end_sum_total += DD_depolymerize_plus_end_sum;
			DD_depolymerize_minus_end_sum_total += DD_depolymerize_minus_end_sum;
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_up] == true)
			{
				if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
				{
					TT_polymerize_plus_end_sum++;
					temp_TT_polymerize_end_vertical_plus_sum++;
				}
				else
				{
					TD_polymerize_plus_end_sum++;
					temp_TD_polymerize_end_vertical_plus_sum++;
				}
			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_right] == true)
			{
				if ((--polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
				{
					TT_polymerize_plus_end_sum++;
					temp_TT_polymerize_end_horizontal_plus_sum++;
				}
				else
				{
					TD_polymerize_plus_end_sum++;
					temp_TD_polymerize_end_horizontal_plus_sum++;
				}
			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_down] == true)
			{
				if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
				{
					TT_polymerize_minus_end_sum++;
					temp_TT_polymerize_end_vertical_minus_sum++;
				}
				else
				{
					TD_polymerize_minus_end_sum++;
					temp_TD_polymerize_end_vertical_minus_sum++;
				}

			}
			if (polymer_cluster_iter->availabe_polymerize_direction[direction_left] == true)
			{
				if (polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T)
				{
					TT_polymerize_minus_end_sum++;
					temp_TT_polymerize_end_horizontal_minus_sum++;
				}
				else
				{
					TD_polymerize_minus_end_sum++;
					temp_TD_polymerize_end_horizontal_minus_sum++;
				}

			}
		}
		if (rule_parameters.MinC == 0)
		{
			if (TT_polymerize_plus_end_sum != polymer_cluster_iter->TT_polymerize_plus_end_sum)
			{
				cout << "TT_polymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TT_polymerize_minus_end_sum != polymer_cluster_iter->TT_polymerize_minus_end_sum)
			{
				cout << "TT_polymerize_minus_end_sum" << endl;
				system("pause");
			}
			if (TD_polymerize_plus_end_sum != polymer_cluster_iter->TD_polymerize_plus_end_sum)
			{
				cout << "TD_polymerize_plus_end_sum" << endl;
				system("pause");
			}
			if (TD_polymerize_minus_end_sum != polymer_cluster_iter->TD_polymerize_minus_end_sum)
			{
				cout << "TD_polymerize_minus_end_sum" << endl;
				system("pause");
			}
		}
		

	}
	if (TT_depolymerize_plus_end_sum_total != rule_parameters.TT_depolymerize_plus_end_sum)
	{
		cout << "TT_depolymerize_plus_end_sum_total" << TT_depolymerize_plus_end_sum_total << endl;
		system("pause");
	}
	if (TT_depolymerize_minus_end_sum_total != rule_parameters.TT_depolymerize_minus_end_sum)
	{
		cout << "TT_depolymerize_minus_end_sum_total" << TT_depolymerize_minus_end_sum_total << endl;
		system("pause");
	}
	if (TD_depolymerize_plus_end_sum_total != rule_parameters.TD_depolymerize_plus_end_sum)
	{
		cout << "TD_depolymerize_plus_end_sum_total" << TD_depolymerize_plus_end_sum_total << endl;
		system("pause");
	}
	if (TD_depolymerize_minus_end_sum_total != rule_parameters.TD_depolymerize_minus_end_sum)
	{
		cout << "TD_depolymerize_minus_end_sum_total" << TD_depolymerize_minus_end_sum_total << endl;
		system("pause");
	}
	if (DD_depolymerize_plus_end_sum_total != rule_parameters.DD_depolymerize_plus_end_sum)
	{
		cout << "DD_depolymerize_plus_end_sum_total" << DD_depolymerize_plus_end_sum_total << endl;
		system("pause");
	}
	if (DD_depolymerize_minus_end_sum_total != rule_parameters.DD_depolymerize_minus_end_sum)
	{
		cout << "DD_depolymerize_minus_end_sum_total" << DD_depolymerize_minus_end_sum_total << endl;
		system("pause");
	}




	if (rule_parameters.MinC == 0)
	{
		if (temp_TT_polymerize_end_vertical_plus_sum != rule_parameters.TT_polymerize_end_vertical_plus_sum)
		{
			cout << " temp_TT_polymerize_end_vertical_plus_sum" << temp_TT_polymerize_end_vertical_plus_sum << endl;
			system("pause");
		}
		if (temp_TT_polymerize_end_vertical_minus_sum != rule_parameters.TT_polymerize_end_vertical_minus_sum)
		{
			cout << " temp_TT_polymerize_end_vertical_minus_sum" << temp_TT_polymerize_end_vertical_minus_sum << endl;
			system("pause");
		}
		if (temp_TT_polymerize_end_horizontal_plus_sum != rule_parameters.TT_polymerize_end_horizontal_plus_sum)
		{
			cout << " temp_TT_polymerize_end_horizontal_plus_sum" << temp_TT_polymerize_end_horizontal_plus_sum << endl;
			system("pause");
		}
		if (temp_TT_polymerize_end_horizontal_minus_sum != rule_parameters.TT_polymerize_end_horizontal_minus_sum)
		{
			cout << " temp_TT_polymerize_end_horizontal_minus_sum" << temp_TT_polymerize_end_horizontal_minus_sum << endl;
			system("pause");
		}
		if (temp_TD_polymerize_end_vertical_plus_sum != rule_parameters.TD_polymerize_end_vertical_plus_sum)
		{
			cout << " temp_TD_polymerize_end_vertical_plus_sum" << temp_TD_polymerize_end_vertical_plus_sum << endl;
			system("pause");
		}
		if (temp_TD_polymerize_end_vertical_minus_sum != rule_parameters.TD_polymerize_end_vertical_minus_sum)
		{
			cout << " temp_TD_polymerize_end_vertical_minus_sum" << temp_TD_polymerize_end_vertical_minus_sum << endl;
			system("pause");
		}
		if (temp_TD_polymerize_end_horizontal_plus_sum != rule_parameters.TD_polymerize_end_horizontal_plus_sum)
		{
			cout << " temp_TD_polymerize_end_horizontal_plus_sum" << temp_TD_polymerize_end_horizontal_plus_sum << endl;
			system("pause");
		}
		if (temp_TD_polymerize_end_horizontal_minus_sum != rule_parameters.TD_polymerize_end_horizontal_minus_sum)
		{
			cout << " temp_TD_polymerize_end_horizontal_minus_sum" << temp_TD_polymerize_end_horizontal_minus_sum << endl;
			system("pause");
		}
	}
	



}

void check_T_D_distance()
{

	double distance = 0;
	double temp_distance = 0;
	auto cluster_iter = rule_structure.polymer_cluster_list.begin();
	rule_parameters.average_T_distance = 0;
	rule_parameters.average_D_distance = 0;
	for (int i = 0; i < rule_structure.polymer_cluster_list.size(); i++)
	{
		if (cluster_iter->polarize == direction_first)
		{
			temp_distance = 0;
		}
		else
		{
			temp_distance = cluster_iter->polymer_sequence.size() + 1;
		}
		auto unit_iter = cluster_iter->polymer_sequence.begin();
		for (int j = 0; j < cluster_iter->polymer_sequence.size(); j++)
		{
			if (cluster_iter->polarize == direction_first)
			{
				temp_distance ++;
			}
			else
			{
				temp_distance --;
			}
			if (unit_iter->Hydrolysis_mark == T)
			{
				rule_parameters.average_T_distance += (double)temp_distance;
			}
			else
			{
				rule_parameters.average_D_distance += (double)temp_distance;
			}
			unit_iter++;
		}
		cluster_iter++;
	}
	rule_parameters.average_T_distance = (double)rule_parameters.average_T_distance / rule_parameters.T_sum;
	rule_parameters.average_D_distance = (double)rule_parameters.average_D_distance / rule_parameters.D_sum;
}
void check_length_distribution()
{
	for (int i = 0; i < rule_parameters.column_num; i++)
	{
		rule_parameters.length_distribution[i] = 0;
	}
	auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	for (auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin(); polymer_cluster_iter != rule_structure.polymer_cluster_list.end(); polymer_cluster_iter++)
	{
		rule_parameters.length_distribution[polymer_cluster_iter->polymer_sequence.size()-1]++;
	}
}

void check_average_length_per_column()
{
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		rule_parameters.average_length_per_column[i] = 0;
	}
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		rule_parameters.polymer_count[i] = 0;
	}
	vector<int> sum_per_column(rule_parameters.row_num);
	//vector<int> count_per_column(rule_parameters.row_num);
	auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	for (auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin(); polymer_cluster_iter != rule_structure.polymer_cluster_list.end(); polymer_cluster_iter++)
	{
		sum_per_column[(polymer_cluster_iter->polymer_sequence.begin())->polymer_grid_pointer->row_position] += polymer_cluster_iter->polymer_sequence.size();
		rule_parameters.polymer_count[(polymer_cluster_iter->polymer_sequence.begin())->polymer_grid_pointer->row_position]++;
	}
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		if (rule_parameters.polymer_count[i] != 0)
		{
			rule_parameters.average_length_per_column[i] = (double)sum_per_column[i] / rule_parameters.polymer_count[i];
		}
		
	}

}
void frap(int left,int right,int up,int down)
{
	//left, right, up, down
	//for (int i = 0; i < rule_parameters.column_num; i++)
	for (int i = down; i < up; i++)
	{
		for (int j = left-1; j < right-1; j++)
		{
			
			if (&Grid[i][j].grid_polymer_unit_pointer != &rule_structure.polymer_unit_list_end.begin())
			{//if there is a ftsZ, need to be marked
				Grid[i][j].grid_polymer_unit_pointer->frap_mark = 1;
			}
			if (&Grid[i][j].grid_anchor_pointer != &rule_structure.anchor_list.end())
			{
				Grid[i][j].grid_anchor_pointer->frap_mark = 1;
			}
			//auto grid_polymer_unit_pointer = &Grid[i][j].grid_polymer_unit_pointer;
			//auto polymer_unit_null = &rule_structure.polymer_unit_list_end.begin();
			//auto grid_anchor_pointer = &Grid[i][j].grid_anchor_pointer;
			//auto anchor_unit_null = &rule_structure.anchor_list.end();

			//if (grid_polymer_unit_pointer != polymer_unit_null)
			//{//if there is a ftsZ, need to be marked
			//	Grid[i][j].grid_polymer_unit_pointer->frap_mark = 1;
			//}
			//if ( grid_anchor_pointer!=anchor_unit_null )
			//{
			//	Grid[i][j].grid_anchor_pointer->frap_mark = 1;
			//}
		}
	}
}
double calculate_diffusion_rate(list<Polymer_Bundle>::iterator polymer_bundle_iter, int direction)
{
	double rate = 0;
	
	if (polymer_bundle_iter->availabe_diffusion_direction[direction] == true)
	{
		if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
		{
			cout << "FtsA num=0, calculation error" << endl;
			system("pause");
		}
		
		//double rate = 1.0 / FtsA_num*exp(-FtsA_num*rule_parameters.diffusion_rate/5000);
		//double rate = 1.0 / FtsA_num*exp(-FtsA_num*rule_parameters.diffusion_rate / 10000.0);
		double FtsZ_sum = 0;
		double FtsA_sum = 0;
		auto polymer_cluster_FtsZ_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
		while (polymer_cluster_FtsZ_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
		{
			FtsZ_sum += (*polymer_cluster_FtsZ_iter)->polymer_sequence.size();
			polymer_cluster_FtsZ_iter++;
		}
		if (direction == direction_up || direction == direction_down)
		{
			//FtsZ_sum = sqrt(FtsZ_sum);
			//FtsA_sum = sqrt(polymer_bundle_iter->polymer_anchor_list.size());
			FtsA_sum = polymer_bundle_iter->polymer_anchor_list.size();
			//FtsZ_sum = log(FtsZ_sum)+1;
			//FtsA_sum = log(polymer_bundle_iter->polymer_anchor_list.size())+1;
		}
		if (direction == direction_left || direction == direction_right)
		{
			/*auto polymer_cluster_FtsZ_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
			while (polymer_cluster_FtsZ_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
			{
				FtsZ_sum += (*polymer_cluster_FtsZ_iter)->polymer_sequence.size();
				polymer_cluster_FtsZ_iter++;
			}*/
			FtsZ_sum = pow(FtsZ_sum,1.5);
			FtsA_sum = pow(polymer_bundle_iter->polymer_anchor_list.size(),1.5);
			
			//FtsZ_sum = log(FtsZ_sum)+1;
			//FtsZ_sum = sqrt(FtsZ_sum);
			//FtsA_sum = polymer_bundle_iter->polymer_anchor_list.size();
			
			

		}
		//rate = (1.0 / polymer_bundle_iter->polymer_anchor_list.size()) / FtsZ_sum;
		rate = (1.0 / (FtsA_sum+FtsZ_sum));
	}
	else 
	{
		rate = 0;
	}
	
	return rate;
};
void validate_diffusion_propensity()
{
	double old_diffusion_propensity = rule_parameters.polymer_diffusion_propensity;
	double diffusion_direction_temp_sum = 0;
	list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	while (polymer_bundle_iter != rule_structure.polymer_bundle_list.end())
	{
		check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
		int temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			temp_sum += polymer_bundle_iter->availabe_diffusion_direction[i];
			diffusion_direction_temp_sum += calculate_diffusion_rate(polymer_bundle_iter, i);
		}

		if (temp_sum != polymer_bundle_iter->available_polymer_diffusion_direction_sum)
		{
			std::cout << "error: diffusion direction sum is not the same as sum up" <<temp_sum<<" "<< polymer_bundle_iter->available_polymer_diffusion_direction_sum << endl;
			system("pause");
		}



		polymer_bundle_iter++;

	}

	if (fabs(diffusion_direction_temp_sum - old_diffusion_propensity) > 0.000001 || fabs(diffusion_direction_temp_sum - rule_parameters.polymer_diffusion_propensity) > 0.000001)
	{
		cout << "polymer diffusion propensity sum doesn't match: " << diffusion_direction_temp_sum <<" "<< old_diffusion_propensity << endl;
		system("pause");
		rule_parameters.polymer_diffusion_propensity = diffusion_direction_temp_sum;
	}
	else
	{
		rule_parameters.polymer_diffusion_propensity = diffusion_direction_temp_sum;
	}

};
void initialize_diffusion_propensity()
{
	double diffusion_direction_temp_sum = 0;
	list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	while (polymer_bundle_iter != rule_structure.polymer_bundle_list.end())
	{
		check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
		int temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			temp_sum += polymer_bundle_iter->availabe_diffusion_direction[i];
			diffusion_direction_temp_sum += calculate_diffusion_rate(polymer_bundle_iter, i);
		}

		if (temp_sum != polymer_bundle_iter->available_polymer_diffusion_direction_sum)
		{
			std::cout << "error: diffusion direction sum is not the same as sum up" << endl;
			system("pause");
		}



		polymer_bundle_iter++;

	}

	
	rule_parameters.polymer_diffusion_propensity = diffusion_direction_temp_sum;
	

};

void analyse_biased_diffusion_propensity()
{
	double biased_diffusion_up_sum = 0;
	double biased_diffusion_down_sum = 0;
	list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	while (polymer_bundle_iter != rule_structure.polymer_bundle_list.end())
	{
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_cluster_iter_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
		while (polymer_cluster_iter_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
		{
			auto polymer_cluster_iter = *polymer_cluster_iter_iter;
			if (polymer_cluster_iter->direction == vertical)
			{
				if (polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position>=rule_parameters.event_left && polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position<=rule_parameters.event_right)
				{
					biased_diffusion_up_sum += calculate_diffusion_rate(polymer_bundle_iter, direction_up);
					biased_diffusion_up_sum += calculate_diffusion_rate(polymer_bundle_iter, direction_down);
					break;
				}
				else
				{
					break;
				}
			}
			else
			{
				break;
			}
			polymer_cluster_iter_iter++;
		}





		polymer_bundle_iter++;

	}
	double biased_current = abs(biased_diffusion_up_sum - biased_diffusion_down_sum);
	if (biased_current > rule_parameters.biased_previous && rule_parameters.reaction_mark!=0 && rule_parameters.reaction_mark!=18)
	{
		rule_parameters.event_output << rule_parameters.reaction_mark << " " << biased_current << endl;
	}
	rule_parameters.biased_previous = biased_current;

	


};