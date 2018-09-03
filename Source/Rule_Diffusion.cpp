#include "Rule_Diffusion.h"
#include "Grid_Unit.h"
#include <list>
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include <cmath>
#include <random>
#include "parameters.h"
#include "General_Functions.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
//keypoint: copy the links and rebuild the boundary consitions.

void polymer_diffusion()
{
	//cout << " enter diffusion" << endl;
	//output();
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);
#endif // debug_random
	////random_device rd;
	//mt19937_64 gen(rd());
	//default_random_engine gen;
	uniform_real_distribution<> Diffusion_Rand(0, rule_parameters.polymer_diffusion_propensity);
	//a empty anchor could also diffuse, there is a problem//solved by empty anchor diffusion function
	//i change Diffusion_Rand(1, rule_parameters.polymer_diffusion_propensity) from 0 to 1, potentially a problem
	double diffusion_rand = Diffusion_Rand(gen);
	list<Polymer_Bundle>::iterator diffusion_bundle_iter;
	double diffusion_temp_sum = 0;
	diffusion_bundle_iter = rule_structure.polymer_bundle_list.begin();
	int store_col;
	int store_row;
	double temp_sum = 0;
	do
	{
		/*if(polymerize_polymer_iter->polymerizable_end_sum==both)
		polymerize_temp_sum = polymerize_temp_sum + 2;
		if (polymerize_polymer_iter->polymerizable_end_sum == up || polymerize_polymer_iter->polymerizable_end_sum == down)
		polymerize_temp_sum = polymerize_temp_sum + 1;		*/
		/*temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			temp_sum = temp_sum + diffusion_bundle_iter->availabe_diffusion_direction[i];
		}*/
		temp_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			//temp_sum = diffusion_bundle_iter->available_polymer_diffusion_direction_sum*calculate_diffusion_rate(diffusion_bundle_iter, i);
			temp_sum += calculate_diffusion_rate(diffusion_bundle_iter, i);
		}
		
		diffusion_temp_sum = diffusion_temp_sum + temp_sum;
		diffusion_bundle_iter++;
		/*if (diffusion_bundle_iter == rule_structure.polymer_bundle_list.end())
		{
			cout << "error:diffusion bundle iter out of range" << endl;
			system("pause");
		}*/
	}while (diffusion_temp_sum < diffusion_rand);
	diffusion_bundle_iter--;
	if (diffusion_bundle_iter == rule_structure.polymer_bundle_list.end())
	{
		cout << "error:diffusion bundle iter out of range" << endl;
		system("pause");
	}
	//uniform_int_distribution<> Direction_Rand(1, diffusion_bundle_iter->available_polymer_diffusion_direction_sum);
	////list<int>::iterator direction_iter;//this is the iterater to the current bundle available diffusion direction list
	////direction_iter = diffusion_bundle_iter->diffusion_direction.begin();
	//int direction_rand=Direction_Rand(gen);
	//
	///*for (int i = 0; i < direction_rand; i++)
	//{
	//	direction_iter++;
	//}*/
	//int diffusion_direction = 0;
	//int direction_temp = 0;
	//for (int i = 0; i < 4; i++)
	//{
	//	if (diffusion_bundle_iter->availabe_diffusion_direction[i] == yes)
	//	{
	//		direction_temp++;
	//	}
	//	if (direction_temp == direction_rand)
	//	{
	//		diffusion_direction = i;
	//		break;
	//	}
	//}
	uniform_real_distribution<> Direction_Rand(0, temp_sum);
	double direction_rand = Direction_Rand(gen);
	double direction_sum = 0;
	int diffusion_direction = 0;
	do
	{
		direction_sum += calculate_diffusion_rate(diffusion_bundle_iter, diffusion_direction);
		diffusion_direction++;
	} while (direction_sum < direction_rand);
	diffusion_direction--;


	/*if (diffusion_bundle_iter->polymer_anchor_list.size()>1)
	{
		cout << "error: diffusion with more than 1 anchor" << endl;
		system("pause");
	}*/
	//check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);
	// i still need to consider how to mark bundle sites.
	list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
	polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
	//check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);
	//unlink all the bundle realtion, it will be rebuild after diffusion
	for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		list<Anchor_Unit>::iterator anchor_iter;
		polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
		while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
		{
			
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
			
			polymer_unit_iter++;

		}
		

		polymer_cluster_iter++;

	}


	if (diffusion_direction == direction_up)
	{
		//rule_parameters.diffusion_count++;
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		//first copy each element in the list
		//list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		int col;
		int row;
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
			//for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
			{
				/*if (polymer_unit_iter == (--(*polymer_cluster_iter)->polymer_sequence.end()))
				{
					polymer_unit_iter++;
				}*/
				polymer_unit_iter--;




				if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{//if it is indeed in one of the annealing pair//first direction
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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






				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer=rule_structure.anchor_list.end();
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;
				//up case, need to consider two corners
				// T T							T T
				// TTT             ->           T T
				//	T							 T		
				//								 T
				
				rule_parameters.diffusion_statistics[row][direction_up]++;


				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;
				

				store_col = col;
				store_row = row;
				//anchor_iter = polymer_unit_iter->anchor_pointer;

				//Grid[col][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				if (polymer_unit_iter->attatch_mark == 1)
				{//it can do attatch before
					/*if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
						Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
					{*/
					list<Anchor_Unit>::iterator anchor_iter = Grid[col][row].grid_anchor_pointer;
					anchor_iter->attatch_mark = 0;
					polymer_unit_iter->attatch_mark = 0;
					rule_parameters.emtpy_anchor_attatch_sum--;
					//}

				}

				//Grid[col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				col = col + 1;
				if (col > rule_parameters.column_num - 1)
				{
					col = 0;
				}
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = rule_structure.polymer_bundle_list.end();
				//polymer_unit_iter->polymer_grid_pointer = &Grid[polymer_unit_iter->polymer_grid_pointer->col_position+1][polymer_unit_iter->polymer_grid_pointer->row_position];
				polymer_unit_iter->polymer_grid_pointer = &Grid[col][row];
				/*if(anchor_iter!=rule_structure.anchor_list.end())
					anchor_iter->anchor_grid_pointer = &Grid[col][row];*/
				//I don't know what is happening here for the anchor
				//list<Anchor_Unit>::iterator anchor_iter = *polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin();
				//anchor_iter->anchor_grid_pointer = polymer_unit_iter->polymer_grid_pointer;
				
				//I think below it is already +1, I don't need to do +1 again

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_anchor_pointer = anchor_iter;
				Grid[col][row].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[col][row].grid_anchor_pointer = anchor_iter;


				//now consider anchor stuff here
				//if (Grid[store_col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster
				//	if (&*Grid[store_col][row].grid_anchor_pointer->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
				//	{//if this anchor is indeed attatch to this cluster, need to move together with this cluster
				//		anchor_iter = Grid[store_col][row].grid_anchor_pointer;
				//		//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				//		Grid[col][row].grid_anchor_pointer = anchor_iter;
				//		anchor_iter->anchor_grid_pointer = &Grid[col][row];
				//		


				//	}
				//	else
				//	{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

				//	}
				//}



				if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster
					
					anchor_iter = polymer_unit_iter->anchor_pointer;
					//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					Grid[col][row].grid_anchor_pointer = anchor_iter;
					anchor_iter->anchor_grid_pointer = &Grid[col][row];				
				}

				
				
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer = *polymer_cluster_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = diffusion_bundle_iter;
				
				
				
				
			} 
			//check diffusion availability
			






			//then rebuild the boundary conditions.//Maybe only need to rebuild the first one and the last one?
			//maybe i have to write a general function to do it.

			//this should not a function of polymer_unit, but of polymer_cluster.

			polymer_cluster_iter++;
			
		}
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
			{
				
				polymer_unit_iter--;
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;

				if (&*Grid[down_col_position][row].grid_polymer_unit_pointer == &*polymer_unit_iter)
				{
					Grid[down_col_position][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				}
				
				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					if (Grid[down_col_position][row].grid_anchor_pointer!=rule_structure.anchor_list.end()&&
						&*Grid[down_col_position][row].grid_anchor_pointer== &*Grid[col][row].grid_anchor_pointer)
					{//if this down grid point to current anchor					
						Grid[down_col_position][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					}
					else
					{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

					}

					//if (Grid[down_col_position][row].grid_anchor_pointer->attatch_mark == 1)
					//{//i can do attatch before
					//	if (Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end() &&
					//		Grid[down_col_position][row].grid_polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
					//	{//it still has a anchor on its bottom and the polymer_unit does not link to anything
					//	 //do nothing
					//	}
					//	else
					//	{
					//		list<Anchor_Unit>::iterator anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
					//		anchor_iter->attatch_mark = 0;
					//		polymer_unit_iter->attatch_mark = 0;
					//		rule_parameters.emtpy_anchor_attatch_sum--;
					//	}
					//}
					
				}
				//if (polymer_unit_iter->attatch_mark == 1)
				//{//it can do attatch before
				//	if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				//		Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//	{
				//		list<Anchor_Unit>::iterator anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
				//		anchor_iter->attatch_mark = 0;
				//		polymer_unit_iter->attatch_mark = 0;
				//		rule_parameters.emtpy_anchor_attatch_sum--;
				//	}

				//}

				
		
			

			}
			

			polymer_cluster_iter++;

		}

		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
			{

				polymer_unit_iter--;
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int next_up_col_position = up_col_position + 1;
				if (next_up_col_position > rule_parameters.column_num - 1)
					next_up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int next_down_col_position = down_col_position - 1;
				if (next_down_col_position < 0)
					next_down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int next_right_row_position = right_row_position + 1;
				if (next_right_row_position > rule_parameters.row_num - 1)
					next_right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;
				int next_left_row_position = left_row_position - 1;
				if (next_left_row_position < 0)
					next_left_row_position = rule_parameters.row_num - 1;

				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}

				if (Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[next_down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[next_down_col_position][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[down_col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[down_col_position][left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[down_col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[down_col_position][right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				
				//check down_left and down_right bundle
				if (&*Grid[down_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[down_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[down_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[down_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[next_down_col_position][row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[next_down_col_position][row].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, 1);



			}


			polymer_cluster_iter++;

		}
		

		//here I assume there are only one anchor
		//list<Anchor_Unit>::iterator anchor_iter;
		//anchor_iter = (*diffusion_bundle_iter->polymer_anchor_list.begin());
		//anchor_iter->anchor_grid_pointer= &Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position + 1];
		//Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//diffusion finished, I have to rebuild the boundary conditions here

		//first to check up case
		//the basic thinking is that to move every polymer_unit up and check if there is already a defferent element there.

		check_polymer_bundle_diffusable_direction(diffusion_bundle_iter,2);//important!!
		//check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);


	}
	if (diffusion_direction == direction_down)
	{
		
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		//first copy each element in the list
		//list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		int col;
		int row;
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
				//for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
			{
				/*if (polymer_unit_iter == (--(*polymer_cluster_iter)->polymer_sequence.end()))
				{
				polymer_unit_iter++;
				}*/
				
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer=rule_structure.anchor_list.end();




				if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{//if it is indeed in one of the annealing pair//first direction
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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




				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;
				//up case, need to consider two corners
				// T T							T T
				// TTT             ->           T T
				//	T							 T		
				//								 T

				rule_parameters.diffusion_statistics[row][direction_down]++;


				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;


				store_col = col;
				store_row = row;
				//anchor_iter = polymer_unit_iter->anchor_pointer;

				//Grid[col][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();

				if (polymer_unit_iter->attatch_mark == 1)
				{//it can do attatch before
				 /*if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				 Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				 {*/
					list<Anchor_Unit>::iterator anchor_iter = Grid[col][row].grid_anchor_pointer;
					anchor_iter->attatch_mark = 0;
					polymer_unit_iter->attatch_mark = 0;
					rule_parameters.emtpy_anchor_attatch_sum--;
					//}

				}
				//Grid[col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				col = col - 1;
				if (col < 0)
				{
					col = rule_parameters.column_num-1;
				}
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = rule_structure.polymer_bundle_list.end();
				//polymer_unit_iter->polymer_grid_pointer = &Grid[polymer_unit_iter->polymer_grid_pointer->col_position+1][polymer_unit_iter->polymer_grid_pointer->row_position];
				polymer_unit_iter->polymer_grid_pointer = &Grid[col][row];
				/*if(anchor_iter!=rule_structure.anchor_list.end())
				anchor_iter->anchor_grid_pointer = &Grid[col][row];*/
				//I don't know what is happening here for the anchor
				//list<Anchor_Unit>::iterator anchor_iter = *polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin();
				//anchor_iter->anchor_grid_pointer = polymer_unit_iter->polymer_grid_pointer;

				//I think below it is already +1, I don't need to do +1 again

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_anchor_pointer = anchor_iter;
				Grid[col][row].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[col][row].grid_anchor_pointer = anchor_iter;


				//now consider anchor stuff here
				//if (Grid[store_col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster
				//	if (&*Grid[store_col][row].grid_anchor_pointer->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
				//	{//if this anchor is indeed attatch to this cluster, need to move together with this cluster
				//		anchor_iter = Grid[store_col][row].grid_anchor_pointer;
				//		//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				//		Grid[col][row].grid_anchor_pointer = anchor_iter;
				//		anchor_iter->anchor_grid_pointer = &Grid[col][row];



				//	}
				//	else
				//	{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

				//	}
				//}
				if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster

					anchor_iter = polymer_unit_iter->anchor_pointer;
					//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					Grid[col][row].grid_anchor_pointer = anchor_iter;
					anchor_iter->anchor_grid_pointer = &Grid[col][row];
				}

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer = *polymer_cluster_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = diffusion_bundle_iter;



				polymer_unit_iter++;
			}
			//check diffusion availability







			//then rebuild the boundary conditions.//Maybe only need to rebuild the first one and the last one?
			//maybe i have to write a general function to do it.

			//this should not a function of polymer_unit, but of polymer_cluster.

			polymer_cluster_iter++;

		}
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
			{

				
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;

				if (&*Grid[up_col_position][row].grid_polymer_unit_pointer == &*polymer_unit_iter)
				{
					Grid[up_col_position][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				}

				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					if (Grid[up_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end() &&
						&*Grid[up_col_position][row].grid_anchor_pointer == &*Grid[col][row].grid_anchor_pointer)
					{//if this down grid point to current anchor					
						Grid[up_col_position][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					}
					else
					{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

					}

					//if (Grid[down_col_position][row].grid_anchor_pointer->attatch_mark == 1)
					//{//i can do attatch before
					//	if (Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end() &&
					//		Grid[down_col_position][row].grid_polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
					//	{//it still has a anchor on its bottom and the polymer_unit does not link to anything
					//	 //do nothing
					//	}
					//	else
					//	{
					//		list<Anchor_Unit>::iterator anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
					//		anchor_iter->attatch_mark = 0;
					//		polymer_unit_iter->attatch_mark = 0;
					//		rule_parameters.emtpy_anchor_attatch_sum--;
					//	}
					//}

				}
				//if (polymer_unit_iter->attatch_mark == 1)
				//{//it can do attatch before
				//	if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				//		Grid[up_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//	{
				//		list<Anchor_Unit>::iterator anchor_iter = Grid[up_col_position][row].grid_anchor_pointer;
				//		anchor_iter->attatch_mark = 0;
				//		polymer_unit_iter->attatch_mark = 0;
				//		rule_parameters.emtpy_anchor_attatch_sum--;
				//	}

				//}


				polymer_unit_iter++;


			}


			polymer_cluster_iter++;

		}

		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
			{

				
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int next_up_col_position = up_col_position + 1;
				if (next_up_col_position > rule_parameters.column_num - 1)
					next_up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int next_down_col_position = down_col_position - 1;
				if (next_down_col_position < 0)
					next_down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int next_right_row_position = right_row_position + 1;
				if (next_right_row_position > rule_parameters.row_num - 1)
					next_right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;
				int next_left_row_position = left_row_position - 1;
				if (next_left_row_position < 0)
					next_left_row_position = rule_parameters.row_num - 1;

				if (Grid[next_up_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[next_up_col_position][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[up_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[up_col_position][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[up_col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[up_col_position][right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[up_col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[up_col_position][left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}

				//check down_left and down_right bundle
				if (&*Grid[up_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[up_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[up_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[up_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[next_up_col_position][row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[next_up_col_position][row].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, 1);


				polymer_unit_iter++;
			}


			polymer_cluster_iter++;

		}


		//here I assume there are only one anchor
		//list<Anchor_Unit>::iterator anchor_iter;
		//anchor_iter = (*diffusion_bundle_iter->polymer_anchor_list.begin());
		//anchor_iter->anchor_grid_pointer= &Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position + 1];
		//Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//diffusion finished, I have to rebuild the boundary conditions here

		//first to check up case
		//the basic thinking is that to move every polymer_unit up and check if there is already a defferent element there.

		check_polymer_bundle_diffusable_direction(diffusion_bundle_iter,2);//important!!
																		 //check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);
	}
	if (diffusion_direction == direction_right)
	{
		//right case should be similar to up case
		//first copy each element in the list
		//list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		int col;
		int row;
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
				//for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
			{
				/*if (polymer_unit_iter == (--(*polymer_cluster_iter)->polymer_sequence.end()))
				{
				polymer_unit_iter++;
				}*/
				polymer_unit_iter--;




				if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{//if it is indeed in one of the annealing pair//first direction
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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



				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer=rule_structure.anchor_list.end();
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;
				//up case, need to consider two corners
				// T T							T T
				// TTT             ->           T T
				//	T							 T		
				//								 T

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;


				store_col = col;
				store_row = row;
				//anchor_iter = polymer_unit_iter->anchor_pointer;

				//Grid[col][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				if (polymer_unit_iter->attatch_mark == 1)
				{//it can do attatch before
				 /*if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				 Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				 {*/
					list<Anchor_Unit>::iterator anchor_iter = Grid[col][row].grid_anchor_pointer;
					anchor_iter->attatch_mark = 0;
					polymer_unit_iter->attatch_mark = 0;
					rule_parameters.emtpy_anchor_attatch_sum--;
					//}

				}

				//Grid[col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				row = row + 1;
				if (row == rule_parameters.row_num)
				{
					cout<<"error: try to diffuse right at right boundary";
					system("pause");
				}
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = rule_structure.polymer_bundle_list.end();
				//polymer_unit_iter->polymer_grid_pointer = &Grid[polymer_unit_iter->polymer_grid_pointer->col_position+1][polymer_unit_iter->polymer_grid_pointer->row_position];
				polymer_unit_iter->polymer_grid_pointer = &Grid[col][row];
				/*if(anchor_iter!=rule_structure.anchor_list.end())
				anchor_iter->anchor_grid_pointer = &Grid[col][row];*/
				//I don't know what is happening here for the anchor
				//list<Anchor_Unit>::iterator anchor_iter = *polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin();
				//anchor_iter->anchor_grid_pointer = polymer_unit_iter->polymer_grid_pointer;

				//I think below it is already +1, I don't need to do +1 again

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_anchor_pointer = anchor_iter;
				Grid[col][row].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[col][row].grid_anchor_pointer = anchor_iter;


				//now consider anchor stuff here
				//if (Grid[col][store_row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster
				//	if (&*Grid[col][store_row].grid_anchor_pointer->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
				//	{//if this anchor is indeed attatch to this cluster, need to move together with this cluster
				//		anchor_iter = Grid[col][store_row].grid_anchor_pointer;
				//		//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				//		Grid[col][row].grid_anchor_pointer = anchor_iter;
				//		anchor_iter->anchor_grid_pointer = &Grid[col][row];



				//	}
				//	else
				//	{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

				//	}
				//}
				if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster

					anchor_iter = polymer_unit_iter->anchor_pointer;
					//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					Grid[col][row].grid_anchor_pointer = anchor_iter;
					anchor_iter->anchor_grid_pointer = &Grid[col][row];
				}

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer = *polymer_cluster_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = diffusion_bundle_iter;




			}
			//check diffusion availability







			//then rebuild the boundary conditions.//Maybe only need to rebuild the first one and the last one?
			//maybe i have to write a general function to do it.

			//this should not a function of polymer_unit, but of polymer_cluster.

			polymer_cluster_iter++;

		}
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
			{

				polymer_unit_iter--;
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;

				if (&*Grid[col][left_row_position].grid_polymer_unit_pointer == &*polymer_unit_iter)
				{
					Grid[col][left_row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				}

				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					if (Grid[col][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
						&*Grid[col][left_row_position].grid_anchor_pointer == &*Grid[col][row].grid_anchor_pointer)
					{//if this down grid point to current anchor					
						Grid[col][left_row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
					}
					else
					{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

					}

					//if (Grid[down_col_position][row].grid_anchor_pointer->attatch_mark == 1)
					//{//i can do attatch before
					//	if (Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end() &&
					//		Grid[down_col_position][row].grid_polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
					//	{//it still has a anchor on its bottom and the polymer_unit does not link to anything
					//	 //do nothing
					//	}
					//	else
					//	{
					//		list<Anchor_Unit>::iterator anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
					//		anchor_iter->attatch_mark = 0;
					//		polymer_unit_iter->attatch_mark = 0;
					//		rule_parameters.emtpy_anchor_attatch_sum--;
					//	}
					//}

				}
				//if (polymer_unit_iter->attatch_mark == 1)
				//{//it can do attatch before
				//	if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				//		Grid[col][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				//	{
				//		list<Anchor_Unit>::iterator anchor_iter = Grid[col][left_row_position].grid_anchor_pointer;
				//		anchor_iter->attatch_mark = 0;
				//		polymer_unit_iter->attatch_mark = 0;
				//		rule_parameters.emtpy_anchor_attatch_sum--;
				//	}

				//}





			}


			polymer_cluster_iter++;

		}

		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.end();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.begin())
			{

				polymer_unit_iter--;
				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int next_up_col_position = up_col_position + 1;
				if (next_up_col_position > rule_parameters.column_num - 1)
					next_up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int next_down_col_position = down_col_position - 1;
				if (next_down_col_position < 0)
					next_down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int next_right_row_position = right_row_position + 1;
				if (next_right_row_position > rule_parameters.row_num - 1)
					next_right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;
				int next_left_row_position = left_row_position - 1;
				if (next_left_row_position < 0)
					next_left_row_position = rule_parameters.row_num - 1;


				if (Grid[col][next_left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][next_left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[col][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[up_col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[up_col_position][left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[down_col_position][left_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[down_col_position][left_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}

				//check up_left and down_left bundle
				if (&*Grid[down_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[down_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[up_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[up_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[col][next_left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[col][next_left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, 1);



			}


			polymer_cluster_iter++;

		}


		//here I assume there are only one anchor
		//list<Anchor_Unit>::iterator anchor_iter;
		//anchor_iter = (*diffusion_bundle_iter->polymer_anchor_list.begin());
		//anchor_iter->anchor_grid_pointer= &Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position + 1];
		//Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//diffusion finished, I have to rebuild the boundary conditions here

		//first to check up case
		//the basic thinking is that to move every polymer_unit up and check if there is already a defferent element there.

		check_polymer_bundle_diffusable_direction(diffusion_bundle_iter,2);//important!!
																		 //check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);

	}
	if (diffusion_direction == direction_left)
	{
		//first copy each element in the list
		//list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		int col;
		int row;
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
				//for (int j = 0; j < (*polymer_cluster_iter)->polymer_sequence.size(); j++)//can I do it like this?
			{
				/*if (polymer_unit_iter == (--(*polymer_cluster_iter)->polymer_sequence.end()))
				{
				polymer_unit_iter++;
				}*/

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_anchor_pointer=rule_structure.anchor_list.end();





				if (polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{//if it is indeed in one of the annealing pair//first direction
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_unit_iter->anealing_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (annealing_pair_iter->second->polymer_grid_pointer->col_position - annealing_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (annealing_pair_iter->second->polymer_grid_pointer->row_position - annealing_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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
					//int col_distance = (bundle_pair_iter->second->polymer_grid_pointer->col_position - bundle_pair_iter->first->polymer_grid_pointer->col_position + column_num) % column_num;
					//int row_distance = (bundle_pair_iter->second->polymer_grid_pointer->row_position - bundle_pair_iter->first->polymer_grid_pointer->row_position + row_num) % row_num;
					//if ((col_distance*col_distance + row_distance*row_distance) != 1)
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






				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;
				//up case, need to consider two corners
				// T T							T T
				// TTT             ->           T T
				//	T							 T		
				//								 T

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;


				store_col = col;
				store_row = row;
				//anchor_iter = polymer_unit_iter->anchor_pointer;

				//Grid[col][row].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				if (polymer_unit_iter->attatch_mark == 1)
				{//it can do attatch before
				 /*if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				 Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				 {*/
					list<Anchor_Unit>::iterator anchor_iter = Grid[col][row].grid_anchor_pointer;
					anchor_iter->attatch_mark = 0;
					polymer_unit_iter->attatch_mark = 0;
					rule_parameters.emtpy_anchor_attatch_sum--;
					//}

				}

				//Grid[col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				row = row - 1;
				if (row < 0)
				{
					cout << "error: try to diffuse left at left boundary";
					system("pause");
				}
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = rule_structure.polymer_bundle_list.end();
				//polymer_unit_iter->polymer_grid_pointer = &Grid[polymer_unit_iter->polymer_grid_pointer->col_position+1][polymer_unit_iter->polymer_grid_pointer->row_position];
				polymer_unit_iter->polymer_grid_pointer = &Grid[col][row];
				/*if(anchor_iter!=rule_structure.anchor_list.end())
				anchor_iter->anchor_grid_pointer = &Grid[col][row];*/
				//I don't know what is happening here for the anchor
				//list<Anchor_Unit>::iterator anchor_iter = *polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.begin();
				//anchor_iter->anchor_grid_pointer = polymer_unit_iter->polymer_grid_pointer;

				//I think below it is already +1, I don't need to do +1 again

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_anchor_pointer = anchor_iter;
				Grid[col][row].grid_polymer_unit_pointer = polymer_unit_iter;
				//Grid[col][row].grid_anchor_pointer = anchor_iter;


				//now consider anchor stuff here
				//if (Grid[col][store_row].grid_anchor_pointer != rule_structure.anchor_list.end())
				//{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster
				//	if (&*Grid[col][store_row].grid_anchor_pointer->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
				//	{//if this anchor is indeed attatch to this cluster, need to move together with this cluster
				//		anchor_iter = Grid[col][store_row].grid_anchor_pointer;
				//		//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
				//		Grid[col][row].grid_anchor_pointer = anchor_iter;
				//		anchor_iter->anchor_grid_pointer = &Grid[col][row];



				//	}
				//	else
				//	{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

				//	}
				//}
				if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor, we need to judge if this anchor is attatched with this polymer_cluster

					anchor_iter = polymer_unit_iter->anchor_pointer;
					//Grid[store_col][row].grid_anchor_pointer = rule_structure.anchor_list.end();
					Grid[col][row].grid_anchor_pointer = anchor_iter;
					anchor_iter->anchor_grid_pointer = &Grid[col][row];
				}

				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer = *polymer_cluster_iter;
				//Grid[polymer_unit_iter->polymer_grid_pointer->col_position][polymer_unit_iter->polymer_grid_pointer->row_position /*+ 1*/].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = diffusion_bundle_iter;



				polymer_unit_iter++;
			}
			//check diffusion availability







			//then rebuild the boundary conditions.//Maybe only need to rebuild the first one and the last one?
			//maybe i have to write a general function to do it.

			//this should not a function of polymer_unit, but of polymer_cluster.

			polymer_cluster_iter++;

		}
		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
			{


				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;

				if (&*Grid[col][right_row_position].grid_polymer_unit_pointer == &*polymer_unit_iter)
				{
					Grid[col][right_row_position].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
				}

				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					if (Grid[col][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end() &&
						&*Grid[col][right_row_position].grid_anchor_pointer == &*Grid[col][row].grid_anchor_pointer)
					{//if this down grid point to current anchor					
						Grid[col][right_row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
					}
					else
					{//if this anchor is not attatch to this polymer_cluster, do not need to move,we need to do nothing.

					}

					//if (Grid[down_col_position][row].grid_anchor_pointer->attatch_mark == 1)
					//{//i can do attatch before
					//	if (Grid[down_col_position][row].grid_anchor_pointer != rule_structure.anchor_list.end() &&
					//		Grid[down_col_position][row].grid_polymer_unit_pointer->anchor_pointer == rule_structure.anchor_list.end())
					//	{//it still has a anchor on its bottom and the polymer_unit does not link to anything
					//	 //do nothing
					//	}
					//	else
					//	{
					//		list<Anchor_Unit>::iterator anchor_iter = Grid[down_col_position][row].grid_anchor_pointer;
					//		anchor_iter->attatch_mark = 0;
					//		polymer_unit_iter->attatch_mark = 0;
					//		rule_parameters.emtpy_anchor_attatch_sum--;
					//	}
					//}

				}
				//if (polymer_unit_iter->attatch_mark == 1)
				//{//it can do attatch before
				//	if (Grid[col][row].grid_anchor_pointer == rule_structure.anchor_list.end() &&
				//		Grid[col][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				//	{
				//		list<Anchor_Unit>::iterator anchor_iter = Grid[col][right_row_position].grid_anchor_pointer;
				//		anchor_iter->attatch_mark = 0;
				//		polymer_unit_iter->attatch_mark = 0;
				//		rule_parameters.emtpy_anchor_attatch_sum--;
				//	}

				//}


				polymer_unit_iter++;


			}


			polymer_cluster_iter++;

		}

		polymer_cluster_iter = diffusion_bundle_iter->bundle_cluster_pointer_list.begin();
		//iterate the cluster list
		//before this i also need to change the boundary condition
		for (int i = 0; i < diffusion_bundle_iter->bundle_cluster_pointer_list.size(); i++)
		{
			list<Polymer_Unit>::iterator polymer_unit_iter;
			list<Anchor_Unit>::iterator anchor_iter;
			polymer_unit_iter = (*polymer_cluster_iter)->polymer_sequence.begin();
			while (polymer_unit_iter != (*polymer_cluster_iter)->polymer_sequence.end())
			{


				col = polymer_unit_iter->polymer_grid_pointer->col_position;
				row = polymer_unit_iter->polymer_grid_pointer->row_position;

				/*int up_col_position = col + 1;
				if (up_col_position > column_num - 1)
					up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > row_num - 1)
					right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = row_num - 1;*/


				int up_col_position = col + 1;
				if (up_col_position > rule_parameters.column_num - 1)
					up_col_position = 0;
				int next_up_col_position = up_col_position + 1;
				if (next_up_col_position > rule_parameters.column_num - 1)
					next_up_col_position = 0;
				int down_col_position = col - 1;
				if (down_col_position < 0)
					down_col_position = rule_parameters.column_num - 1;
				int next_down_col_position = down_col_position - 1;
				if (next_down_col_position < 0)
					next_down_col_position = rule_parameters.column_num - 1;
				int right_row_position = row + 1;
				if (right_row_position > rule_parameters.row_num - 1)
					right_row_position = 0;
				int next_right_row_position = right_row_position + 1;
				if (next_right_row_position > rule_parameters.row_num - 1)
					next_right_row_position = 0;
				int left_row_position = row - 1;
				if (left_row_position < 0)
					left_row_position = rule_parameters.row_num - 1;
				int next_left_row_position = left_row_position - 1;
				if (next_left_row_position < 0)
					next_left_row_position = rule_parameters.row_num - 1;




				if (Grid[col][next_right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][next_right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[col][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[up_col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[up_col_position][right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[down_col_position][right_row_position].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[down_col_position][right_row_position].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}
				if (Grid[col][row].grid_anchor_pointer != rule_structure.anchor_list.end())
				{//if there is a anchor
					anchor_iter = Grid[col][row].grid_anchor_pointer;
					check_anchor_diffusion(anchor_iter);
				}

				//check down_right and up_right bundle
				if (&*Grid[up_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[up_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[down_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[down_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
				if (&*Grid[col][next_right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
					check_polymer_bundle_diffusable_direction(Grid[col][next_right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, 1);


				polymer_unit_iter++;
			}


			polymer_cluster_iter++;

		}


		//here I assume there are only one anchor
		//list<Anchor_Unit>::iterator anchor_iter;
		//anchor_iter = (*diffusion_bundle_iter->polymer_anchor_list.begin());
		//anchor_iter->anchor_grid_pointer= &Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position + 1];
		//Grid[anchor_iter->anchor_grid_pointer->col_position][anchor_iter->anchor_grid_pointer->row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//diffusion finished, I have to rebuild the boundary conditions here

		//first to check up case
		//the basic thinking is that to move every polymer_unit up and check if there is already a defferent element there.

		check_polymer_bundle_diffusable_direction(diffusion_bundle_iter,2);//important!!
																		 //check_polymer_bundle_diffusable_direction(diffusion_bundle_iter);
	}
	//rule_parameters.ps[rule_diffusion_polymer] = rule_parameters.polymer_diffusion_propensity*diffusion_rate;
	

	//check_sumps;
	//cout << "diffusion finished" << endl;
	//cout << store_col+1 << " " << store_row << endl;
	//output();
	//system("pause");
};
void empty_anchor_diffusion()
{
	//cout << " enter empty anchor diffusion" << endl;
	//output();
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);;
#endif // debug_random
	////random_device rd;
	//mt19937_64 gen(rd());
	//default_random_engine gen;
	uniform_int_distribution<> Diffusion_Rand(1, rule_parameters.empty_anchor_diffusion_direction_sum);
	//a empty anchor could also diffuse, there is a problem
	int diffusion_rand = Diffusion_Rand(gen);
	list<Anchor_Unit>::iterator diffusion_anchor_iter;
	int diffusion_temp_sum = 0;
	diffusion_anchor_iter = rule_structure.anchor_list.begin();
	//int store_col;
	//int store_row;
	int temp_sum = 0;
	do
	{
		if (&*diffusion_anchor_iter->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
		{//empty
			temp_sum = diffusion_anchor_iter->available_polymer_diffusion_direction_sum;
		}
		else
			temp_sum = 0;
		diffusion_temp_sum = diffusion_temp_sum + temp_sum;
			
		diffusion_anchor_iter++;

	} while (diffusion_temp_sum < diffusion_rand);
	diffusion_anchor_iter--;
	if (&*diffusion_anchor_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
	{
		cout << "empty anchor diffusion function tries to diffuse a anchor with polymer" << endl;
		system("pause");
	}
	uniform_int_distribution<> Direction_Rand(1, temp_sum);
	
	int direction_rand = Direction_Rand(gen);
	
	int diffusion_direction = 0;
	int direction_temp = 0;
	for (int i = 0; i < 4; i++)
	{
		if (diffusion_anchor_iter->availabe_diffusion_direction[i] == yes)
		{
			direction_temp++;
		}
		if (direction_temp == direction_rand)
		{
			diffusion_direction = i;
			break;
		}
	}
	int col_position = diffusion_anchor_iter->anchor_grid_pointer->col_position;
	int row_position = diffusion_anchor_iter->anchor_grid_pointer->row_position;
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
	if (diffusion_anchor_iter->attatch_mark == 1)
	{
		Grid[col_position][row_position].grid_polymer_unit_pointer->attatch_mark = 0;
		diffusion_anchor_iter->attatch_mark = 0;
		rule_parameters.emtpy_anchor_attatch_sum--;
	}
	if (diffusion_direction == direction_up)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		diffusion_anchor_iter->anchor_grid_pointer = &Grid[up_col_position][row_position];
		Grid[up_col_position][row_position].grid_anchor_pointer = diffusion_anchor_iter;
		Grid[col_position][row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//i need to check all the possible cases
		//next up, left, right, up_left, up_right,down
		check_anchor_diffusion(diffusion_anchor_iter);
		check_anchor_diffusion(Grid[next_up_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][row_position].grid_anchor_pointer);
		//what is it doing here? Here i have to consider empty anchor collision, so i need to do diffusion with a empty anchor, and because it is complecated, i'll use my master diffusion function.
		/*if (&*Grid[next_up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() &&
			Grid[next_up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size()<2)
		{
			check_polymer_bundle_diffusable_direction(Grid[next_up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer, 1);
		}*/
		if (&*Grid[next_up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			check_polymer_bundle_diffusable_direction(Grid[next_up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}

		
	}
	if (diffusion_direction == direction_down)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		diffusion_anchor_iter->anchor_grid_pointer = &Grid[down_col_position][row_position];
		Grid[down_col_position][row_position].grid_anchor_pointer = diffusion_anchor_iter;
		Grid[col_position][row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//i need to check all the possible cases
		//next up, left, right, up_left, up_right,down
		check_anchor_diffusion(diffusion_anchor_iter);
		check_anchor_diffusion(Grid[next_down_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][row_position].grid_anchor_pointer);
		if (&*Grid[next_down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[next_down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
	}
	if (diffusion_direction == direction_right)
	{
		diffusion_anchor_iter->anchor_grid_pointer = &Grid[col_position][right_row_position];
		Grid[col_position][right_row_position].grid_anchor_pointer = diffusion_anchor_iter;
		Grid[col_position][row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//i need to check all the possible cases
		//next up, left, right, up_left, up_right,down
		check_anchor_diffusion(diffusion_anchor_iter);
		check_anchor_diffusion(Grid[col_position][next_right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][right_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][left_row_position].grid_anchor_pointer);

		if (&*Grid[col_position][next_right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][next_right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		//this is for frap
		/*if (right_row_position >= rule_parameters.right_bound)
		{
			diffusion_anchor_iter->frap_mark = 0;
		}*/
		
	}
	if (diffusion_direction == direction_left)
	{
		diffusion_anchor_iter->anchor_grid_pointer = &Grid[col_position][left_row_position];
		Grid[col_position][left_row_position].grid_anchor_pointer = diffusion_anchor_iter;
		Grid[col_position][row_position].grid_anchor_pointer = rule_structure.anchor_list.end();
		//i need to check all the possible cases
		//next up, left, right, up_left, up_right,down
		check_anchor_diffusion(diffusion_anchor_iter);
		check_anchor_diffusion(Grid[col_position][next_left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[up_col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[down_col_position][left_row_position].grid_anchor_pointer);
		check_anchor_diffusion(Grid[col_position][right_row_position].grid_anchor_pointer);
		if (&*Grid[col_position][next_left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][next_left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			check_polymer_bundle_diffusable_direction(Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		if (&*Grid[col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin() )
		{
			check_polymer_bundle_diffusable_direction(Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,1);
		}
		//this is for frap
		/*if (right_row_position <= rule_parameters.left_bound)
		{
			diffusion_anchor_iter->frap_mark = 0;
		}*/
	}
	//rule_parameters.ps[rule_diffusion_polymer] = rule_parameters.polymer_diffusion_propensity*diffusion_rate;


	//check_sumps();
	//cout << "empty anchor diffusion finished" << endl;
	
};