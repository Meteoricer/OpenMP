#include "Rule_Polymerize_Depolymerize.h"
#include "Rule_Anchoring_Deanchoring.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include "Grid_Unit.h"
#include <vector>
#include <cmath>
#include <random>
#include "DeconstructionFunction.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

void TT_polymerize_vertical_plus()
{
	//cout << " enter polymerize function" << endl;
	//system("pause");
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
	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TT_polymerize_end_vertical_plus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == vertical)
		{
			temp_sum = polymerize_polymer_iter->TT_polymerize_plus_end_sum;
			/*for (int i = 0; i < 4; i++)
			{
			temp_sum = temp_sum + polymerize_polymer_iter->availabe_polymerize_direction[i];
			}*/

			/*if(polymerize_polymer_iter->polymerize_end_sum==both)
			polymerize_temp_sum = polymerize_temp_sum + 2;
			if (polymerize_polymer_iter->polymerize_end_sum == up || polymerize_polymer_iter->polymerize_end_sum == down)
			polymerize_temp_sum = polymerize_temp_sum + 1;		*/
			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}
		
		polymerize_polymer_iter++;

	}while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != vertical)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	}
	else
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	}
	//int polymerize_temp_sum = polymerize_temp_sum - int(polymerize_polymer_iter->polymerize_end_sum);
	//now we just choose one polymer_iter to polymer
	//if (temp_sum==1)//if there are only one end that can be polymerized
	//{
	//	/*int coltemp = polymerize_polymer_iter->polymer_grid_pointer->col_position;
	//	int rowtemp = polymerize_polymer_iter->polymer_grid_pointer->row_position;*/
	//	if (polymerize_polymer_iter->availabe_polymerize_direction[direction_up]==yes)
	//	{
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = --polymerize_polymer_iter->polymer_sequence.end();
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;

	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col+1][end_row];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		////polymerize_polymer_iter->polymer_sequence.pop_back();//maybe i need to do a research on how to insert.
	//		//polymerize_polymer_iter->polymer_sequence.push_back(polymer_unit);
	//		//polymer_iter = --polymerize_polymer_iter->polymer_sequence.end();
	//		//Grid[end_col + 1][end_row].grid_polymer_unit_pointer = polymer_iter;
	//		//Grid[end_col + 1][end_row].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		//Grid[end_col + 1][end_row].grid_polymer_cluster_pointer = polymerize_polymer_iter;
	//		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	//		
	//		//check if it is polymerizeable at this end;
	//		//if (Grid[end_col + 2][end_row].grid_polymer_unit_pointer==NULL)
	//		//{
	//		//	
	//		//}
	//		//if (Grid[end_col + 2][end_row].grid_polymer_unit_pointer != NULL)//have to make a change, leave blank for now
	//		//{
	//		//	
	//		//
	//		//}
	//		
	//	}
	//	if (polymerize_polymer_iter->availabe_polymerize_direction[direction_down] == yes)
	//	{
	//		//cout << "lalala" << endl;
	//		//system("pause");
	//		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		////polymer_iter->polymer_grid_pointer->col_position
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;

	//		///*int end_col = polymerize_polymer_iter->polymer_sequence[polymerize_polymer_iter->polymer_sequence.size() - 1].polymer_grid_pointer->col_position;
	//		//int end_row = polymerize_polymer_iter->polymer_sequence[polymerize_polymer_iter->polymer_sequence.size() - 1].polymer_grid_pointer->row_position;*/
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;
	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col - 1][end_row];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		//polymerize_polymer_iter->polymer_sequence.push_front(polymer_unit);
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();//down, first(front) element
	//		//Grid[end_col - 1][end_row].grid_polymer_unit_pointer = polymer_iter ;
	//		//Grid[end_col - 1][end_row].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		////Grid[end_col - 1][end_row].grid_polymer_cluster_pointer = polymerize_polymer_iter;
	//		
	//		
	//	}
	//	if (polymerize_polymer_iter->availabe_polymerize_direction[direction_left]==true)
	//	{
	//		/*list<Polymer_Unit>::iterator polymer_iter;
	//		polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		Polymer_Unit polymer_unit;
	//		polymer_unit.Hydrolysis_mark = T;
	//		polymer_unit.polymer_grid_pointer = &Grid[end_col][end_row - 1];
	//		polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		Grid[end_col][end_row - 1].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		polymerize_polymer_iter->polymer_sequence.push_front(polymer_unit);
	//		polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		Grid[end_col][end_row - 1].grid_polymer_unit_pointer = polymer_iter;*/
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	//	}
	//	if (polymerize_polymer_iter->availabe_polymerize_direction[direction_right]==true)
	//	{
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.end();
	//		//polymer_iter--;
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;
	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col][end_row+1];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		//polymerize_polymer_iter->polymer_sequence.push_back(polymer_unit);
	//		//polymer_iter++;//last element
	//		//Grid[end_col][end_row + 1].grid_polymer_unit_pointer = polymer_iter;
	//		//Grid[end_col][end_row + 1].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	//		
	//		
	//	}
	//}
	//not possible since maximun we have plus and minus end polymerize
	//if (temp_sum == 2)//if there are two
	//{
	//	int diraction_rand;
	//	//choose one direction
	//	if (polymerize_polymer_iter->direction == vertical)
	//	{
	//		uniform_int_distribution<> Diraction_Rand(direction_up, direction_down);//mark here
	//		diraction_rand = Diraction_Rand(gen);
	//	}
	//	if (polymerize_polymer_iter->direction == horizontal)
	//	{
	//		uniform_int_distribution<> Diraction_Rand(direction_left, direction_right);
	//		diraction_rand = Diraction_Rand(gen);
	//	}
	//	if (diraction_rand == direction_up)
	//	{
	//		//cout << " enter polymerize" << endl;
	//		//system("pause");
	//		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.end();
	//		//polymer_iter--;
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;

	//		//
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;
	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col + 1][end_row];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		//polymerize_polymer_iter->polymer_sequence.push_back(polymer_unit);
	//		//polymer_iter++;//the last element
	//		//Grid[end_col + 1][end_row].grid_polymer_unit_pointer = polymer_iter;
	//		//Grid[end_col + 1][end_row].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;


	//	}
	//	if (diraction_rand == direction_down)
	//	{
	//		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;
	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col - 1][end_row];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		//polymerize_polymer_iter->polymer_sequence.push_front(polymer_unit);
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();//down, first(front) element
	//		//Grid[end_col - 1][end_row].grid_polymer_unit_pointer = polymer_iter;
	//		//Grid[end_col - 1][end_row].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//	}
	//	if (diraction_rand == direction_left)
	//	{
	//		/*list<Polymer_Unit>::iterator polymer_iter;
	//		polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		Polymer_Unit polymer_unit;
	//		polymer_unit.Hydrolysis_mark = T;
	//		polymer_unit.polymer_grid_pointer = &Grid[end_col][end_row - 1];
	//		polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		Grid[end_col][end_row - 1].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		polymerize_polymer_iter->polymer_sequence.push_front(polymer_unit);
	//		polymer_iter = polymerize_polymer_iter->polymer_sequence.begin();
	//		Grid[end_col][end_row - 1].grid_polymer_unit_pointer = polymer_iter;*/
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	//	}
	//	if (diraction_rand == direction_right)
	//	{
	//		//list<Polymer_Unit>::iterator polymer_iter;
	//		//polymer_iter = polymerize_polymer_iter->polymer_sequence.end();
	//		//polymer_iter--;
	//		//int end_col = polymer_iter->polymer_grid_pointer->col_position;
	//		//int end_row = polymer_iter->polymer_grid_pointer->row_position;
	//		//Polymer_Unit polymer_unit;
	//		//polymer_unit.Hydrolysis_mark = T;
	//		//polymer_unit.polymer_grid_pointer = &Grid[end_col][end_row + 1];
	//		//polymer_unit.polymer_cluster_pointer = polymerize_polymer_iter;
	//		//polymerize_polymer_iter->polymer_sequence.push_back(polymer_unit);
	//		//polymer_iter++;//last element
	//		//Grid[end_col][end_row + 1].grid_polymer_unit_pointer = polymer_iter;
	//		//Grid[end_col][end_row + 1].grid_polymer_unit_pointer->polymer_cluster_pointer = polymerize_polymer_iter;
	//		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	//	}
	//	
	//}
	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer,2);

	
	//cout << "polymerize finished" << endl;
	
	//output();
};



void TT_polymerize_vertical_minus()
{
	
#ifdef DEBUG_RANDOM
	
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	
#endif 
	
	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TT_polymerize_end_vertical_minus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == vertical)
		{
			temp_sum = polymerize_polymer_iter->TT_polymerize_minus_end_sum;
			
			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != vertical)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	}
	else
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	}
	
	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};

void TD_polymerize_vertical_plus()
{

#ifdef DEBUG_RANDOM

	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());

#endif 

	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TD_polymerize_end_vertical_plus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == vertical)
		{
			temp_sum = polymerize_polymer_iter->TD_polymerize_plus_end_sum;

			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != vertical)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	}
	else
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	}

	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};

void TD_polymerize_vertical_minus()
{

#ifdef DEBUG_RANDOM

	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());

#endif 

	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TD_polymerize_end_vertical_minus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == vertical)
		{
			temp_sum = polymerize_polymer_iter->TD_polymerize_minus_end_sum;

			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != vertical)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_up);
	}
	else
	{
		rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_down);
	}

	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};


void TT_polymerize_horizontal_plus()
{
	
#ifdef DEBUG_RANDOM
	
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	
#endif 
	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TT_polymerize_end_horizontal_plus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == horizontal)
		{
			temp_sum = polymerize_polymer_iter->TT_polymerize_plus_end_sum;
			
			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != horizontal)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_left]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	}
	else
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	}
	
	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);


	
};



void TT_polymerize_horizontal_minus()
{

#ifdef DEBUG_RANDOM

	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());

#endif 

	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TT_polymerize_end_horizontal_minus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == horizontal)
		{
			temp_sum = polymerize_polymer_iter->TT_polymerize_minus_end_sum;

			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != horizontal)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	}
	else
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	}

	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};

void TD_polymerize_horizontal_plus()
{

#ifdef DEBUG_RANDOM

	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());

#endif 

	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TD_polymerize_end_horizontal_plus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == horizontal)
		{
			temp_sum = polymerize_polymer_iter->TD_polymerize_plus_end_sum;

			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != horizontal)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	}
	else
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	}

	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};

void TD_polymerize_horizontal_minus()
{

#ifdef DEBUG_RANDOM

	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());

#endif 

	uniform_int_distribution<> Polymerize_Rand(1, rule_parameters.TD_polymerize_end_horizontal_minus_sum);//change this into begin with 1
	int polymerize_rand = Polymerize_Rand(gen);
	list<Polymer_Cluster>::iterator polymerize_polymer_iter;
	int polymerize_temp_sum = 0;
	polymerize_polymer_iter = rule_structure.polymer_cluster_list.begin();
	int temp_sum = 0;
	do
	{
		if (polymerize_polymer_iter->direction == horizontal)
		{
			temp_sum = polymerize_polymer_iter->TD_polymerize_minus_end_sum;

			polymerize_temp_sum = polymerize_temp_sum + temp_sum;
		}

		polymerize_polymer_iter++;

	} while (polymerize_temp_sum < polymerize_rand);
	polymerize_polymer_iter--;
	if (polymerize_polymer_iter->direction != horizontal)
	{
		cout << "polymerize cluster direction error" << endl;
		system("pause");
	}
	if (polymerize_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (polymerize_polymer_iter->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (polymerize_polymer_iter->polarize == direction_first)
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_right);
	}
	else
	{
		//rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
		Polymer_Unit polymer_unit(polymerize_polymer_iter, direction_left);
	}

	check_polymer_hydrolysis_sum(polymerize_polymer_iter->polymer_sequence.begin());
	check_polymer_bundle_diffusable_direction(polymerize_polymer_iter->cluster_bundle_pointer, 2);

};





void depolymerizeTT_plus()
{
	
	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.TT_depolymerize_plus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{

		//temp_sum = 0;
		/*if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
				(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
			{
				depolymerize_temp_sum++;
				temp_sum++;
			}
			if (--depolymerize_polymer_cluster_iter->polymer_sequence.end()->Hydrolysis_mark == T &&
				(----depolymerize_polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
			{
				depolymerize_temp_sum++;
				temp_sum++;
			}


		}

		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
				(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
			{
				depolymerize_temp_sum++;
			}
		}*/
		
		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
		
		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		temp_sum++;
		depolymerize_polymer_cluster_iter++;
	} while (depolymerize_temp_sum <depolymerize_rand || depolymerize_end_sum==0);
	depolymerize_polymer_cluster_iter--;

	


	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{

			//list<Polymer_Unit>::iterator polymer_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
			//{//there is an anchor on the first unit, cut the second unit
			//	//polymer_unit_iter++;
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(++polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, ++polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//	//everything about diffusion should be already take care of in ~polymer_unit()
			//}
			//else
			//{//there is an anchor on the last unit, cut the first unit
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//}
			////check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//check_polymer_hydrolysis_sum(polymer_iter);
			
			
			
			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;
			
			
			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer,2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{
				
				
				{
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());
					if (depolymerize_polymer_cluster_iter->polarize == direction_first)
					
					{//if it is at the begin side
					 
						rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;


						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer,2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}

							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);
								/*if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer != polymer_bundle_iter)
								{
									cout << "polymer_iter_back do not point to the same bundle as polymer_bundle_iter should be" << endl;
									system("pause");
								}*/
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									//cout << polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() << endl;
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										
										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);


						}

					}
					else
					{
						
						{//if it is at the end side
						 
							rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer,2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);


							}
						}
					}
				}
				
				

			}
		}
	}
	
	
	
	
	/*cluster_size = cluster_size - depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (cluster_size == 0)
	{
		output();
		cout << "error: the size change for depolymerize is 0." << endl;
		cout << "cluster size change: " << cluster_size << endl;
		system("pause");
	}*/
	
 	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
	//check_polymer_bundle_diffusable_direction(depolymerize_polymer_cluster_iter->cluster_bundle_pointer);
	//rule_parameters.ps[rule_depolymerize_TT] = rule_parameters.TT_depolymerize_end_sum*TT_fragmentation_rate ;
	//cout << "TT depolymerize finished" << endl;
	//output();
	//system("pause");
	
};






void depolymerizeTT_minus()
{

	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.TT_depolymerize_minus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{

		//temp_sum = 0;
		/*if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}
		if (--depolymerize_polymer_cluster_iter->polymer_sequence.end()->Hydrolysis_mark == T &&
		(----depolymerize_polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}


		}

		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		}
		}*/

		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;

		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		temp_sum++;
		depolymerize_polymer_cluster_iter++;
	} while (depolymerize_temp_sum < depolymerize_rand || depolymerize_end_sum == 0);
	depolymerize_polymer_cluster_iter--;




	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{

			//list<Polymer_Unit>::iterator polymer_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
			//{//there is an anchor on the first unit, cut the second unit
			//	//polymer_unit_iter++;
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(++polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, ++polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//	//everything about diffusion should be already take care of in ~polymer_unit()
			//}
			//else
			//{//there is an anchor on the last unit, cut the first unit
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//}
			////check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//check_polymer_hydrolysis_sum(polymer_iter);



			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;


			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{


				{
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());
					if (depolymerize_polymer_cluster_iter->polarize == direction_second)

					{//if it is at the begin side

						rule_parameters.direction_record[rule_parameters.reaction_mark][direction_up]++;
						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;


						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}

							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
								/*if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer != polymer_bundle_iter)
								{
								cout << "polymer_iter_back do not point to the same bundle as polymer_bundle_iter should be" << endl;
								system("pause");
								}*/
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									//cout << polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() << endl;
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{

										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


						}

					}
					else
					{

						{//if it is at the end side

							rule_parameters.direction_record[rule_parameters.reaction_mark][direction_down]++;
							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


							}
						}
					}
				}



			}
		}
	}
};
void depolymerizeTD_plus()
{
	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.TD_depolymerize_plus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{
		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		//temp_sum = 0;
		/*if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}
		if (--depolymerize_polymer_cluster_iter->polymer_sequence.end()->Hydrolysis_mark == T &&
		(----depolymerize_polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}


		}

		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		}
		}*/

		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
		
		
		temp_sum++;
		depolymerize_polymer_cluster_iter++;
	} while (depolymerize_temp_sum <depolymerize_rand || depolymerize_end_sum == 0);
	depolymerize_polymer_cluster_iter--;
	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{

			
			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;
			
			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{
				
				{//there are only one possible depolymerize end, then it is either these two cases;
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());

				
					if (depolymerize_polymer_cluster_iter->polarize == direction_first)
					{//if it is at the begin side
					
						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


						}

					}
					else
					{
						
						{//if it is at the end side
						 
							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


							}
						}
					}
				}
			


			}
		}
	}

};
void depolymerizeTD_minus()
{
	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.TD_depolymerize_minus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{
		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		//temp_sum = 0;
		/*if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}
		if (--depolymerize_polymer_cluster_iter->polymer_sequence.end()->Hydrolysis_mark == T &&
		(----depolymerize_polymer_cluster_iter->polymer_sequence.end())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		temp_sum++;
		}


		}

		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.begin()->Hydrolysis_mark == T &&
		(++depolymerize_polymer_cluster_iter->polymer_sequence.begin())->Hydrolysis_mark == T)
		{
		depolymerize_temp_sum++;
		}
		}*/

		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;


		temp_sum++;
		depolymerize_polymer_cluster_iter++;
		} while (depolymerize_temp_sum <depolymerize_rand || depolymerize_end_sum == 0);
	depolymerize_polymer_cluster_iter--;
	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{


			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;

			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{

				{//there are only one possible depolymerize end, then it is either these two cases;
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());


					if (depolymerize_polymer_cluster_iter->polarize == direction_second)
					{//if it is at the begin side

						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();


						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;



						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


						}

					}
					else
					{
						
						{//if it is at the end side

							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


							}
						}
					}
				}



			}
		}
	}

};
void depolymerizeDD_plus()
{
	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.DD_depolymerize_plus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{

		

		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
		
		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		temp_sum++;
		depolymerize_polymer_cluster_iter++;
	} while (depolymerize_temp_sum <depolymerize_rand || depolymerize_end_sum == 0);
	depolymerize_polymer_cluster_iter--;
	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{

			//list<Polymer_Unit>::iterator polymer_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
			//{//there is an anchor on the first unit, cut the second unit
			//	//polymer_unit_iter++;
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(++polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, ++polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//	//everything about diffusion should be already take care of in ~polymer_unit()
			//}
			//else
			//{//there is an anchor on the last unit, cut the first unit
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//}
			////check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//check_polymer_hydrolysis_sum(polymer_iter);
			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{
				
				{//there are only one possible depolymerize end, then it is either these two cases;
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());

					if (depolymerize_polymer_cluster_iter->polarize==direction_first)	
					{//if it is at the begin side
					
						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


						}

					}
					else
					{
						
						{//if it is at the end side
						 //list<Polymer_Unit>::iterator polymer_unit_iter = --depolymerize_polymer_cluster_iter->polymer_sequence.end();
						 //if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
						 //{//there is an anchor on the last unit
						 //	if (depolymerize_polymer_cluster_iter->anchor_pointer_list.size() == 1)
						 //	{//if there is only one anchor, then cut the rest polymer_unit
						 //		//polymer_unit_iter++;
						 //		delete_polymer_unit(depolymerize_polymer_cluster_iter, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_unit_iter);
						 //		//depolymerize_polymer_cluster_iter->polymer_sequence.erase(depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_unit_iter);
						 //		//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
						 //	}
						 //	else
						 //	{
						 //		if ((depolymerize_polymer_cluster_iter->anchor_pointer_list.size() > 1))
						 //		{//if there are more than one anchor, then it behaves like fragmentation
						 //		 //not sure yet, put it empty
						 //		}
						 //	}
						 //	
						 //}
						 //else
						 //{//there is no anchor at the end, just delete this element
						 //	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
						 //	delete_polymer_unit(depolymerize_polymer_cluster_iter,polymer_unit_iter);
						 //	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
						 //}
							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;






							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


							}
						}
					}
				}
			


			}
		}
	}




	
};



void depolymerizeDD_minus()
{
	//cout << " enter depolymerize function" << endl;
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
	uniform_int_distribution<> Depolymerize_Rand(1, rule_parameters.DD_depolymerize_minus_end_sum);
	int depolymerize_rand = Depolymerize_Rand(gen);
	list<Polymer_Cluster>::iterator depolymerize_polymer_cluster_iter;
	int depolymerize_temp_sum = 0;
	depolymerize_polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int depolymerize_end_sum = 0;
	int temp_sum = 1;
	do
	{



		depolymerize_temp_sum = depolymerize_temp_sum + depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;
		depolymerize_end_sum = depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

		if (temp_sum > rule_structure.polymer_cluster_list.size())
		{
			cout << "loop out of range" << endl;
			system("pause");
		}
		temp_sum++;
		depolymerize_polymer_cluster_iter++;
	} while (depolymerize_temp_sum <depolymerize_rand || depolymerize_end_sum == 0);
	depolymerize_polymer_cluster_iter--;
	//start to do depolymerize here
	int cluster_size = depolymerize_polymer_cluster_iter->polymer_sequence.size();
	if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 1)
	{
		cout << "error: try to do depolymerize with polymer length 1" << endl;
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->col_position << " ";
		cout << depolymerize_polymer_cluster_iter->polymer_sequence.begin()->polymer_grid_pointer->row_position << endl;
		output();
		system("pause");
	}
	else
	{
		if (depolymerize_polymer_cluster_iter->polymer_sequence.size() == 2)
		{

			//list<Polymer_Unit>::iterator polymer_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
			//{//there is an anchor on the first unit, cut the second unit
			//	//polymer_unit_iter++;
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(++polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, ++polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//	//everything about diffusion should be already take care of in ~polymer_unit()
			//}
			//else
			//{//there is an anchor on the last unit, cut the first unit
			//	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
			//	delete_polymer_unit(depolymerize_polymer_cluster_iter, polymer_unit_iter);
			//	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			//}
			////check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
			list<Polymer_Unit>::iterator polymer_iter;
			//list<Polymer_Unit>::iterator polymer_iter_front;
			list<Polymer_Unit>::iterator polymer_iter_back;

			polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = polymer_iter;
			//polymer_iter_back;
			polymer_iter++;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//check_polymer_hydrolysis_sum(polymer_iter);
			Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_direction_first;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



			rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
			rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
			rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

			rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
			rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
			rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


			rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
			rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
			rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
			rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
			rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
			rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

			depolymerize_polymer_cluster_iter->TTT_sum = 0;
			depolymerize_polymer_cluster_iter->TTD_sum = 0;
			depolymerize_polymer_cluster_iter->DTD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_sum = 0;
			depolymerize_polymer_cluster_iter->TD_sum = 0;
			depolymerize_polymer_cluster_iter->DD_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
			depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




			while (polymer_unit_iter_temp != polymer_iter)
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
				polymer_unit_iter_temp++;
			}
			polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
			//use this to store this bundle_iter, will use it to delete
			list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
			list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
			polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
			polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
			polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
			recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
			if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
			{//if they do not belong to the same bundle
				recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
					delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
				}
				else
				{
					check_polymer_hydrolysis_sum(polymer_iter);
					check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
					if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
				if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
				{
					//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
					delete_polymer_bundle(polymer_bundle_iter);
				}
				else
				{

					check_polymer_hydrolysis_sum(polymer_iter_back);
					check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
					if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
						{
							rule_parameters.rotation_sum++;
						}
					}
				}
			}
			else
			{//if they still belong to the same bundle
			 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
				for (int i = 0; i < 4; i++)
				{
					polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

				}
				polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
				delete_polymer_bundle(polymer_bundle_iter_store);
				rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
				check_polymer_hydrolysis_sum(polymer_iter_back);
				check_polymer_hydrolysis_sum(polymer_iter);
				check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


			}
		}
		else
		{
			if (depolymerize_polymer_cluster_iter->polymer_sequence.size() >= 3)
			{

				{//there are only one possible depolymerize end, then it is either these two cases;
					list<Polymer_Unit>::iterator first_unit_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
					list<Polymer_Unit>::iterator center_unit_iter = (++depolymerize_polymer_cluster_iter->polymer_sequence.begin());

					if (depolymerize_polymer_cluster_iter->polarize == direction_second)
					{//if it is at the begin side

						list<Polymer_Unit>::iterator polymer_iter;
						//list<Polymer_Unit>::iterator polymer_iter_front;
						list<Polymer_Unit>::iterator polymer_iter_back;

						polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.begin();
						//polymer_iter++;
						polymer_iter_back = polymer_iter;
						//polymer_iter_back;
						polymer_iter++;
						//now polymer_iter is at TT position, 
						//						  ^
						//					      |
						//check_polymer_hydrolysis_sum(polymer_iter);
						Polymer_Bundle polymer_bundle;
						Polymer_Cluster polymer_cluster_direction_first;
						rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
						list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
						list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



						rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
						rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
						rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

						rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
						rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
						rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


						rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
						rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
						rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
						rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
						rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
						rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

						depolymerize_polymer_cluster_iter->TTT_sum = 0;
						depolymerize_polymer_cluster_iter->TTD_sum = 0;
						depolymerize_polymer_cluster_iter->DTD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_sum = 0;
						depolymerize_polymer_cluster_iter->TD_sum = 0;
						depolymerize_polymer_cluster_iter->DD_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
						depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;




						while (polymer_unit_iter_temp != polymer_iter)
						{
							polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							polymer_unit_iter_temp++;
						}
						polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
						//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
						//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
						rule_structure.polymer_bundle_list.push_front(polymer_bundle);
						list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
						list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
						//use this to store this bundle_iter, will use it to delete
						list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
						list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
						polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
						polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
						polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
						recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
						//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
						if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
						{//if they do not belong to the same bundle
							recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
								delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							}
							else
							{
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
							if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
							{
								//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
								delete_polymer_bundle(polymer_bundle_iter);
							}
							else
							{

								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
								if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
								{
									if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
									{
										rule_parameters.rotation_sum++;
									}
								}
							}
						}
						else
						{//if they still belong to the same bundle
						 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
						 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
							for (int i = 0; i < 4; i++)
							{
								polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

							}
							polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
							delete_polymer_bundle(polymer_bundle_iter_store);
							rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
							check_polymer_hydrolysis_sum(polymer_iter_back);
							check_polymer_hydrolysis_sum(polymer_iter);
							check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


						}

					}
					else
					{

						{//if it is at the end side
						 //list<Polymer_Unit>::iterator polymer_unit_iter = --depolymerize_polymer_cluster_iter->polymer_sequence.end();
						 //if (polymer_unit_iter->anchor_pointer != rule_structure.anchor_list.end())
						 //{//there is an anchor on the last unit
						 //	if (depolymerize_polymer_cluster_iter->anchor_pointer_list.size() == 1)
						 //	{//if there is only one anchor, then cut the rest polymer_unit
						 //		//polymer_unit_iter++;
						 //		delete_polymer_unit(depolymerize_polymer_cluster_iter, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_unit_iter);
						 //		//depolymerize_polymer_cluster_iter->polymer_sequence.erase(depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_unit_iter);
						 //		//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
						 //	}
						 //	else
						 //	{
						 //		if ((depolymerize_polymer_cluster_iter->anchor_pointer_list.size() > 1))
						 //		{//if there are more than one anchor, then it behaves like fragmentation
						 //		 //not sure yet, put it empty
						 //		}
						 //	}
						 //	
						 //}
						 //else
						 //{//there is no anchor at the end, just delete this element
						 //	//depolymerize_polymer_cluster_iter->polymer_sequence.erase(polymer_unit_iter);
						 //	delete_polymer_unit(depolymerize_polymer_cluster_iter,polymer_unit_iter);
						 //	//check_polymer_hydrolysis_sum(depolymerize_polymer_cluster_iter->polymer_sequence.begin());
						 //}
							list<Polymer_Unit>::iterator polymer_iter;
							//list<Polymer_Unit>::iterator polymer_iter_front;
							list<Polymer_Unit>::iterator polymer_iter_back;

							polymer_iter = depolymerize_polymer_cluster_iter->polymer_sequence.end();
							polymer_iter--;
							polymer_iter--;
							polymer_iter_back = polymer_iter;
							//polymer_iter_back;
							polymer_iter++;
							//now polymer_iter is at TT position, 
							//						  ^
							//					      |
							//check_polymer_hydrolysis_sum(polymer_iter);
							Polymer_Bundle polymer_bundle;
							Polymer_Cluster polymer_cluster_direction_first;
							rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
							list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
							list<Polymer_Unit>::iterator polymer_unit_iter_temp = depolymerize_polymer_cluster_iter->polymer_sequence.begin();



							rule_parameters.TTT_sum_total = rule_parameters.TTT_sum_total - depolymerize_polymer_cluster_iter->TTT_sum;
							rule_parameters.TTD_sum_total = rule_parameters.TTD_sum_total - depolymerize_polymer_cluster_iter->TTD_sum;
							rule_parameters.DTD_sum_total = rule_parameters.DTD_sum_total - depolymerize_polymer_cluster_iter->DTD_sum;

							rule_parameters.TT_fragmentation_sum = rule_parameters.TT_fragmentation_sum - depolymerize_polymer_cluster_iter->TT_sum;
							rule_parameters.TD_fragmentation_sum = rule_parameters.TD_fragmentation_sum - depolymerize_polymer_cluster_iter->TD_sum;
							rule_parameters.DD_fragmentation_sum = rule_parameters.DD_fragmentation_sum - depolymerize_polymer_cluster_iter->DD_sum;


							rule_parameters.TT_depolymerize_plus_end_sum = rule_parameters.TT_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum;
							rule_parameters.TD_depolymerize_plus_end_sum = rule_parameters.TD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum;
							rule_parameters.DD_depolymerize_plus_end_sum = rule_parameters.DD_depolymerize_plus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum;
							rule_parameters.TT_depolymerize_minus_end_sum = rule_parameters.TT_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum;
							rule_parameters.TD_depolymerize_minus_end_sum = rule_parameters.TD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum;
							rule_parameters.DD_depolymerize_minus_end_sum = rule_parameters.DD_depolymerize_minus_end_sum - depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum;

							depolymerize_polymer_cluster_iter->TTT_sum = 0;
							depolymerize_polymer_cluster_iter->TTD_sum = 0;
							depolymerize_polymer_cluster_iter->DTD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_sum = 0;
							depolymerize_polymer_cluster_iter->TD_sum = 0;
							depolymerize_polymer_cluster_iter->DD_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_plus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TT_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->TD_depolymerize_minus_end_sum = 0;
							depolymerize_polymer_cluster_iter->DD_depolymerize_minus_end_sum = 0;






							while (polymer_unit_iter_temp != polymer_iter)
							{
								polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
								polymer_unit_iter_temp++;
							}
							polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), depolymerize_polymer_cluster_iter->polymer_sequence, depolymerize_polymer_cluster_iter->polymer_sequence.begin(), polymer_iter);
							//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
							//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), depolymerize_polymer_cluster_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
							rule_structure.polymer_bundle_list.push_front(polymer_bundle);
							list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
							list<Polymer_Bundle>::iterator polymer_bundle_iter_store = depolymerize_polymer_cluster_iter->cluster_bundle_pointer;
							//use this to store this bundle_iter, will use it to delete
							list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
							list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
							polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
							polymer_cluster_iter_direction_first->direction = depolymerize_polymer_cluster_iter->direction;
							polymer_cluster_iter_direction_first->polarize = depolymerize_polymer_cluster_iter->polarize;
							recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
							//strange connect between Grid[14][6] and cluster of [0][5]... why they are in the same bundle that grid{14][6]is not in it?
							if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
							{//if they do not belong to the same bundle
								recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this sentense is to delete the empty cluster;
									delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
								}
								else
								{
									check_polymer_hydrolysis_sum(polymer_iter);
									check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer, 2);
									if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
								if (polymer_bundle_iter->polymer_anchor_list.size() == 0)
								{
									//check_polymer_bundle_diffusable_direction(polymer_bundle_iter);
									delete_polymer_bundle(polymer_bundle_iter);
								}
								else
								{

									check_polymer_hydrolysis_sum(polymer_iter_back);
									check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);
									if (polymer_iter_back->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
									{
										if (polymer_iter_back->polymer_cluster_pointer->polymer_sequence.size() == 1)
										{
											rule_parameters.rotation_sum++;
										}
									}
								}
							}
							else
							{//if they still belong to the same bundle
							 //delete_polymer_bundle(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
							 //the bundle should already be deleted in recursion function, maybe no need to delete it here.
								for (int i = 0; i < 4; i++)
								{
									polymer_bundle_iter->availabe_diffusion_direction[i] = polymer_bundle_iter_store->availabe_diffusion_direction[i];

								}
								polymer_bundle_iter->available_polymer_diffusion_direction_sum = polymer_bundle_iter_store->available_polymer_diffusion_direction_sum;
								delete_polymer_bundle(polymer_bundle_iter_store);
								rule_parameters.polymer_diffusion_direction_sum += polymer_bundle_iter->available_polymer_diffusion_direction_sum;
								check_polymer_hydrolysis_sum(polymer_iter_back);
								check_polymer_hydrolysis_sum(polymer_iter);
								check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


							}
						}
					}
				}



			}
		}
	}





};