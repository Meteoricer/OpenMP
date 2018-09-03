#include "Rule_Hydrolysis.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include "Grid_Unit.h"
#include "parameters.h"
#include <vector>
#include <cmath>
#include <random>
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

void rotate()
{
	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//cout << " enter annealing_TT" << endl;
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
	
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int rotation_sum_temp = 0;
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{
		if (polymer_cluster_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1 && polymer_cluster_iter->polymer_sequence.size()==1)
		{
			//uniform_int_distribution<> Rotation_Rand(0, 1);
			//int rotation_rand = Rotation_Rand(gen);
			//if (rotation_rand == 0)
			{//do rotation


				uniform_real_distribution<> Direction_Rand(0, 1);
				uniform_int_distribution<> Polarize_Rand(0, 1);
				int polarize_rand = Polarize_Rand(gen);
				double direction_rand = Direction_Rand(gen);
				if (direction_rand <= rule_parameters.posibility_vertical)
				{
					if (polymer_cluster_iter->direction == horizontal)
					{



						list<Polymer_Unit>::iterator polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
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
							list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
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
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							cout << "error, rotation unit should not bundled with anything" << endl;
							system("pause");
						}
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							cout << "error, rotation unit should not bundled with anything" << endl;
							system("pause");
						}
						for (int direction = 0; direction < 4; direction++)
						{
							if (polymer_cluster_iter->availabe_polymerize_direction[direction] == true)
							{
								polymer_cluster_iter->availabe_polymerize_direction[direction] = false;
								
								/*polymer_cluster_iter->polymerize_end_sum--;
								if (polymer_cluster_iter->direction == vertical)
								{
									rule_parameters.polymerize_end_vertical_sum--;
								}
								if (polymer_cluster_iter->direction == horizontal)
								{
									rule_parameters.polymerize_end_horizontal_sum--;
								}*/
							}
						}
						if (polymer_cluster_iter->TT_polymerize_minus_end_sum == 1)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum--;
							rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
						}
						if (polymer_cluster_iter->TT_polymerize_plus_end_sum == 1)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum--;
							rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
						}
						if (polymer_cluster_iter->TD_polymerize_minus_end_sum == 1)
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum--;
							rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
						}
						if (polymer_cluster_iter->TD_polymerize_plus_end_sum == 1)
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum--;
							rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
						}
						


						//uniform_real_distribution<> Direction_Rand(0, 1);
						//uniform_int_distribution<> Polarize_Rand(0, 1);
						//int polarize_rand = Polarize_Rand(gen);
						//double direction_rand = Direction_Rand(gen);
						//if (direction_rand <= posibility_vertical)
						//{
						//	if (polymer_cluster_iter->direction == horizontal)
						//	{
						//		polymer_cluster_iter->direction = vertical;
						//		polymer_cluster_iter->polarize = polarize_rand;
						//		check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						//	}
						//	else//direction is the same, but the polarize could be different
						//	{
						//		if (polymer_cluster_iter->polarize != polarize_rand)
						//		{
						//			polymer_cluster_iter->polarize = polarize_rand;
						//			check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						//		}
						//	}

						//}
						//else
						//{
						//	if (polymer_cluster_iter->direction == vertical)
						//	{
						//		polymer_cluster_iter->direction = horizontal;
						//		polymer_cluster_iter->polarize = polarize_rand;
						//		check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						//	}
						//	else//direction is the same, but the polarize could be different
						//	{
						//		if (polymer_cluster_iter->polarize != polarize_rand)
						//		{
						//			polymer_cluster_iter->polarize = polarize_rand;
						//			check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						//		}
						//	}
						//}






						polymer_cluster_iter->direction = vertical;
						polymer_cluster_iter->polarize = polarize_rand;
						check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
					}
					else//direction is the same, but the polarize could be different
					{
						if (polymer_cluster_iter->polarize != polarize_rand)
						{//if polarize has been changed



							list<Polymer_Unit>::iterator polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
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
								list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
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
							if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
							{
								cout << "error, rotation unit should not bundled with anything" << endl;
								system("pause");
							}
							if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
							{
								cout << "error, rotation unit should not bundled with anything" << endl;
								system("pause");
							}
							for (int direction = 0; direction < 4; direction++)
							{
								if (polymer_cluster_iter->availabe_polymerize_direction[direction] == true)
								{
									polymer_cluster_iter->availabe_polymerize_direction[direction] = false;
									/*polymer_cluster_iter->polymerize_end_sum--;
									if (polymer_cluster_iter->direction == vertical)
									{
										rule_parameters.polymerize_end_vertical_sum--;
									}
									if (polymer_cluster_iter->direction == horizontal)
									{
										rule_parameters.polymerize_end_horizontal_sum--;
									}*/
								}
							}

							if (polymer_cluster_iter->TT_polymerize_minus_end_sum == 1)
							{
								polymer_cluster_iter->TT_polymerize_minus_end_sum--;
								rule_parameters.TT_polymerize_end_vertical_minus_sum--;
							}
							if (polymer_cluster_iter->TT_polymerize_plus_end_sum == 1)
							{
								polymer_cluster_iter->TT_polymerize_plus_end_sum--;
								rule_parameters.TT_polymerize_end_vertical_plus_sum--;
							}
							if (polymer_cluster_iter->TD_polymerize_minus_end_sum == 1)
							{
								polymer_cluster_iter->TD_polymerize_minus_end_sum--;
								rule_parameters.TD_polymerize_end_vertical_minus_sum--;
							}
							if (polymer_cluster_iter->TD_polymerize_plus_end_sum == 1)
							{
								polymer_cluster_iter->TD_polymerize_plus_end_sum--;
								rule_parameters.TD_polymerize_end_vertical_plus_sum--;
							}


							



							polymer_cluster_iter->polarize = polarize_rand;
							check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						}
						else
						{//if nothing changed, don't need to do rotation

						}
					}

				}
				else
				{
					if (polymer_cluster_iter->direction == vertical)
					{


						list<Polymer_Unit>::iterator polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
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
							list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
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
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							cout << "error, rotation unit should not bundled with anything" << endl;
							system("pause");
						}
						if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
						{
							cout << "error, rotation unit should not bundled with anything" << endl;
							system("pause");
						}
						for (int direction = 0; direction < 4; direction++)
						{
							if (polymer_cluster_iter->availabe_polymerize_direction[direction] == true)
							{
								polymer_cluster_iter->availabe_polymerize_direction[direction] = false;
								/*polymer_cluster_iter->polymerize_end_sum--;
								if (polymer_cluster_iter->direction == vertical)
								{
									rule_parameters.polymerize_end_vertical_sum--;
								}
								if (polymer_cluster_iter->direction == horizontal)
								{
									rule_parameters.polymerize_end_horizontal_sum--;
								}*/
							}
						}

						if (polymer_cluster_iter->TT_polymerize_minus_end_sum == 1)
						{
							polymer_cluster_iter->TT_polymerize_minus_end_sum--;
							rule_parameters.TT_polymerize_end_vertical_minus_sum--;
						}
						if (polymer_cluster_iter->TT_polymerize_plus_end_sum == 1)
						{
							polymer_cluster_iter->TT_polymerize_plus_end_sum--;
							rule_parameters.TT_polymerize_end_vertical_plus_sum--;
						}
						if (polymer_cluster_iter->TD_polymerize_minus_end_sum == 1)
						{
							polymer_cluster_iter->TD_polymerize_minus_end_sum--;
							rule_parameters.TD_polymerize_end_vertical_minus_sum--;
						}
						if (polymer_cluster_iter->TD_polymerize_plus_end_sum == 1)
						{
							polymer_cluster_iter->TD_polymerize_plus_end_sum--;
							rule_parameters.TD_polymerize_end_vertical_plus_sum--;
						}
						







						polymer_cluster_iter->direction = horizontal;
						polymer_cluster_iter->polarize = polarize_rand;
						check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
					}
					else//direction is the same, but the polarize could be different
					{
						if (polymer_cluster_iter->polarize != polarize_rand)
						{






							list<Polymer_Unit>::iterator polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
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
								list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator bundle_pair_iter = polymer_unit_iter->bundling_polymer_unit_pair_list_iter[direction_first];
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
							if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_first] != rule_structure.debundling_polymer_unit_pair_list.end())
							{
								cout << "error, rotation unit should not bundled with anything" << endl;
								system("pause");
							}
							if (polymer_unit_iter->debundling_polymer_unit_pair_list_iter[direction_second] != rule_structure.debundling_polymer_unit_pair_list.end())
							{
								cout << "error, rotation unit should not bundled with anything" << endl;
								system("pause");
							}
							for (int direction = 0; direction < 4; direction++)
							{
								if (polymer_cluster_iter->availabe_polymerize_direction[direction] == true)
								{
									polymer_cluster_iter->availabe_polymerize_direction[direction] = false;
									/*polymer_cluster_iter->polymerize_end_sum--;
									if (polymer_cluster_iter->direction == vertical)
									{
										rule_parameters.polymerize_end_vertical_sum--;
									}
									if (polymer_cluster_iter->direction == horizontal)
									{
										rule_parameters.polymerize_end_horizontal_sum--;
									}*/
								}
							}


							if (polymer_cluster_iter->TT_polymerize_minus_end_sum == 1)
							{
								polymer_cluster_iter->TT_polymerize_minus_end_sum--;
								rule_parameters.TT_polymerize_end_horizontal_minus_sum--;
							}
							if (polymer_cluster_iter->TT_polymerize_plus_end_sum == 1)
							{
								polymer_cluster_iter->TT_polymerize_plus_end_sum--;
								rule_parameters.TT_polymerize_end_horizontal_plus_sum--;
							}
							if (polymer_cluster_iter->TD_polymerize_minus_end_sum == 1)
							{
								polymer_cluster_iter->TD_polymerize_minus_end_sum--;
								rule_parameters.TD_polymerize_end_horizontal_minus_sum--;
							}
							if (polymer_cluster_iter->TD_polymerize_plus_end_sum == 1)
							{
								polymer_cluster_iter->TD_polymerize_plus_end_sum--;
								rule_parameters.TD_polymerize_end_horizontal_plus_sum--;
							}




							polymer_cluster_iter->polarize = polarize_rand;
							check_polymer_bundle_diffusable_direction(polymer_cluster_iter->cluster_bundle_pointer, 1);
						}
						else
						{//nothing changed, don't need to do rotation

						}
					}
				}
















				


			/*	if (polymer_cluster_iter->direction == vertical)
					polymer_cluster_iter->direction = horizontal;
				else
					polymer_cluster_iter->direction = vertical;*/
				
			}
			//else
			//{//don't do rotation

			//}
		

			
		}
		
		polymer_cluster_iter++;
	}
	
	

};