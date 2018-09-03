#include "Polymer_Unit.h"
#include <list>
#include <vector>
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include <random>

using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;


Polymer_Unit::Polymer_Unit(int i)
{
	this->Hydrolysis_mark = T;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);
#endif // debug_random
	uniform_int_distribution<> Frap_Rand(1, rule_parameters.Z_ctp);
	int frap_rand = Frap_Rand(gen);
	if (frap_rand <= rule_parameters.bleached_Z&&this->frap_mark==0)
	{
		this->frap_mark = 1;
		rule_parameters.bleached_Z--;
	}
	else
	{
		frap_mark = 0;
	}
}



Polymer_Unit::Polymer_Unit()
{
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);
#endif // debug_random
	if (rule_parameters.Z_ctp > 0)
	{
		uniform_int_distribution<> Frap_Rand(1, rule_parameters.Z_ctp);
		int frap_rand = Frap_Rand(gen);
		if (frap_rand <= rule_parameters.bleached_Z&&this->frap_mark == 0)
		{
			this->frap_mark = 1;
			rule_parameters.bleached_Z--;
		}
		else
		{
			frap_mark = 0;
		}
	}
	
	this->lifetime = rule_parameters.t;
	rule_parameters.total_ftsZ_count++;
	this->attatch_mark = 0;
	/*this->bundled_neighbor_polymer_unit_pointer.resize(4);
	this->bundled_neighbor_polymer_unit_pointer[0] = rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[1]= rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[2] = rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[3] = rule_structure.polymer_unit_list_end.begin();*/
	this->Hydrolysis_mark = T;
	this->anealing_polymer_unit_pair_list_iter.resize(2);
	this->anealing_polymer_unit_pair_list_iter[0] = rule_structure.anealing_polymer_unit_pair_list.end();
	this->anealing_polymer_unit_pair_list_iter[1] = rule_structure.anealing_polymer_unit_pair_list.end();
	/*this->TT_anealing_polymer_unit_pair_list_iter.resize(4);
	this->TT_anealing_polymer_unit_pair_list_iter[0]= rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[1] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[2] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[3] = rule_structure.TT_anealing_polymer_unit_pair_list.end();*/
	this->debundling_polymer_unit_pair_list_iter.resize(2);
	this->debundling_polymer_unit_pair_list_iter[0]= rule_structure.debundling_polymer_unit_pair_list.end();
	this->debundling_polymer_unit_pair_list_iter[1] = rule_structure.debundling_polymer_unit_pair_list.end();
	//this->debundling_polymer_unit_pair_list_iter[2] = rule_structure.debundling_polymer_unit_pair_list.end();
	//this->debundling_polymer_unit_pair_list_iter[3] = rule_structure.debundling_polymer_unit_pair_list.end();
	this->bundling_polymer_unit_pair_list_iter.resize(2);
	this->bundling_polymer_unit_pair_list_iter[0]= rule_structure.bundling_polymer_unit_pair_list.end();
	this->bundling_polymer_unit_pair_list_iter[1] = rule_structure.bundling_polymer_unit_pair_list.end();
	//this->bundling_polymer_unit_pair_list_iter[2] = rule_structure.bundling_polymer_unit_pair_list.end();
	//this->bundling_polymer_unit_pair_list_iter[3] = rule_structure.bundling_polymer_unit_pair_list.end();

	this->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
	this->anchor_pointer = rule_structure.anchor_list.end();
	this->polymer_grid_pointer = NULL;
	this->Consider_mark = -1;
	//int new_col, new_row = 0;
	//new_col = col + 1;
	//new_row = row;
	//if (new_col > column_num - 1)
	//	new_col = 0;
	//if (&*Grid[new_col][new_row].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())//nothing there
	//{
	//	this->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = true;
	//	this->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = true;
	//	rule_parameters.polymerize_end_sum_sum = rule_parameters.polymerize_end_sum_sum++;
	//	rule_parameters.polymer_diffusion_propensity++;
	//}
	//new_col = col - 1;
	//new_row = row;
	//if (new_col <0)
	//	new_col = column_num-1;
	//if (&*Grid[new_col][new_row].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())//nothing there
	//{
	//	this->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = true;
	//	this->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = true;
	//	rule_parameters.polymerize_end_sum_sum = rule_parameters.polymerize_end_sum_sum++;
	//	rule_parameters.polymer_diffusion_propensity++;
	//}

	//rule_parameters.polymerize_end_sum_sum = rule_parameters.polymerize_end_sum_sum + 2;
	//this->anchor_pointer = rule_structure.anchor_list_end.begin();
	//this->polymer_cluster_pointer = rule_structure.polymer_cluster_list_end.begin();
};

Polymer_Unit::Polymer_Unit(list<Polymer_Cluster>::iterator polymer_cluster_iter, int direction)
{
	rule_parameters.total_ftsZ_count++;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);
#else
	random_device rd;
	mt19937_64 gen(rd());
	//default_random_engine gen(rule_parameters.step_num);
#endif // debug_random
	uniform_int_distribution<> Frap_Rand(1, rule_parameters.Z_ctp);
	int frap_rand = Frap_Rand(gen);
	if (frap_rand <= rule_parameters.bleached_Z&&this->frap_mark == 0)
	{
		this->frap_mark = 1;
		rule_parameters.bleached_Z--;
	}
	else
	{
		frap_mark = 0;
	}
	this->lifetime = rule_parameters.t;
	this->Consider_mark = -1;
	//this->attatch_mark = 0;
	this->Hydrolysis_mark = T;
	/*this->bundled_neighbor_polymer_unit_pointer.resize(4);
	this->bundled_neighbor_polymer_unit_pointer[0] = rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[1] = rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[2] = rule_structure.polymer_unit_list_end.begin();
	this->bundled_neighbor_polymer_unit_pointer[3] = rule_structure.polymer_unit_list_end.begin();*/
	this->anealing_polymer_unit_pair_list_iter.resize(2);
	this->anealing_polymer_unit_pair_list_iter[0] = rule_structure.anealing_polymer_unit_pair_list.end();
	this->anealing_polymer_unit_pair_list_iter[1] = rule_structure.anealing_polymer_unit_pair_list.end();
	/*this->TT_anealing_polymer_unit_pair_list_iter.resize(4);
	this->TT_anealing_polymer_unit_pair_list_iter[0] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[1] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[2] = rule_structure.TT_anealing_polymer_unit_pair_list.end();
	this->TT_anealing_polymer_unit_pair_list_iter[3] = rule_structure.TT_anealing_polymer_unit_pair_list.end();*/
	this->debundling_polymer_unit_pair_list_iter.resize(2);
	this->debundling_polymer_unit_pair_list_iter[0] = rule_structure.debundling_polymer_unit_pair_list.end();
	this->debundling_polymer_unit_pair_list_iter[1] = rule_structure.debundling_polymer_unit_pair_list.end();
	//this->debundling_polymer_unit_pair_list_iter[2] = rule_structure.debundling_polymer_unit_pair_list.end();
	//this->debundling_polymer_unit_pair_list_iter[3] = rule_structure.debundling_polymer_unit_pair_list.end();
	this->bundling_polymer_unit_pair_list_iter.resize(2);
	this->bundling_polymer_unit_pair_list_iter[0] = rule_structure.bundling_polymer_unit_pair_list.end();
	this->bundling_polymer_unit_pair_list_iter[1] = rule_structure.bundling_polymer_unit_pair_list.end();
	//this->bundling_polymer_unit_pair_list_iter[2] = rule_structure.bundling_polymer_unit_pair_list.end();
	//this->bundling_polymer_unit_pair_list_iter[3] = rule_structure.bundling_polymer_unit_pair_list.end();
	this->polymer_cluster_pointer = rule_structure.polymer_cluster_list.end();
	this->anchor_pointer = rule_structure.anchor_list.end();
	this->polymer_grid_pointer = NULL;
	if (direction == direction_up)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = (--polymer_cluster_iter->polymer_sequence.end());
		int end_col = polymer_unit_iter->polymer_grid_pointer->col_position;
		int end_row = polymer_unit_iter->polymer_grid_pointer->row_position;
		polymer_cluster_iter->polymer_sequence.push_back(*this);
		polymer_unit_iter =(--polymer_cluster_iter->polymer_sequence.end());
		//Polymer_Unit polymer_unit;
		polymer_unit_iter->attatch_mark = 0;
		//polymer_unit_iter->Hydrolysis_mark = T;
		this->Hydrolysis_mark = T;
		//when and where to check the boundary conditions?
		//now i need to consider periodic condition
		int new_col, new_row = 0;
		new_col = end_col + 1;
		new_row = end_row;
		if (new_col > rule_parameters.column_num - 1)
			new_col = 0;

		rule_parameters.polymerize_statistics[new_row][direction_up]++;
		

		if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			cout << "error: try to do polymerize when there is already a polymer_unit there" << endl;
			system("pause");
		}
		polymer_unit_iter->polymer_grid_pointer = &Grid[new_col][new_row];
		polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;
		
		polymer_unit_iter = --polymer_cluster_iter->polymer_sequence.end();
		Grid[new_col][new_row].grid_polymer_unit_pointer = polymer_unit_iter;
		check_anchor_diffusion(polymer_unit_iter->polymer_grid_pointer->grid_anchor_pointer);
		/*new_col = new_col +1;
		new_row = end_row;
		if (new_col > column_num - 1)
			new_col = 0;
		if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			polymer_unit_iter->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = no;
			polymer_unit_iter->polymer_cluster_pointer->polymerize_end_sum--;
			rule_parameters.polymerize_end_sum--;
			cout << "collide" << endl;
			system("pause");
		}*/
		//check polymerizeable end here
		
	}
	if (direction == direction_down)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
		int end_col = polymer_unit_iter->polymer_grid_pointer->col_position;
		int end_row = polymer_unit_iter->polymer_grid_pointer->row_position;
		polymer_cluster_iter->polymer_sequence.push_front(*this);
		polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
		//Polymer_Unit polymer_unit;
		polymer_unit_iter->attatch_mark = 0;
		this->Hydrolysis_mark = T;
		//when and where to check the boundary conditions?
		//now i need to consider periodic condition
		int new_col, new_row = 0;
		new_col = end_col - 1;
		new_row = end_row;
		if (new_col <0)
			new_col = rule_parameters.column_num-1;

		rule_parameters.polymerize_statistics[new_row][direction_down]++;

		if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			cout << "error: try to do polymerize when there is already a polymer_unit there" << endl;
			system("pause");
		}
		polymer_unit_iter->polymer_grid_pointer = &Grid[new_col][new_row];
		polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;
		
		Grid[new_col][new_row].grid_polymer_unit_pointer = polymer_unit_iter;
		check_anchor_diffusion(polymer_unit_iter->polymer_grid_pointer->grid_anchor_pointer);
		//new_col = new_col - 1;
		//new_row = end_row;
		//if (new_col <0)
		//	new_col = column_num-1;
		//if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		//{
		//	polymer_unit_iter->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = no;
		//	polymer_unit_iter->polymer_cluster_pointer->polymerize_end_sum--;
		//	rule_parameters.polymerize_end_sum--;
		//	cout << "collide" << endl;
		//	//system("pause");
		//}
		//check polymerizeable end here
	}
	if (direction == direction_left)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
		int end_col = polymer_unit_iter->polymer_grid_pointer->col_position;
		int end_row = polymer_unit_iter->polymer_grid_pointer->row_position;
		polymer_cluster_iter->polymer_sequence.push_front(*this);
		polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();
		//Polymer_Unit polymer_unit;
		polymer_unit_iter->attatch_mark = 0;
		this->Hydrolysis_mark = T;
		//when and where to check the boundary conditions?
		//now i need to consider periodic condition
		int new_col, new_row = 0;
		new_col = end_col;
		new_row = end_row-1;
		if (new_row <0)
			new_row = rule_parameters.row_num - 1;

		rule_parameters.polymerize_statistics[new_row][direction_left]++;

		if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			cout << "error: try to do polymerize when there is already a polymer_unit there" << endl;
			system("pause");
		}
		polymer_unit_iter->polymer_grid_pointer = &Grid[new_col][new_row];
		polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;

		Grid[new_col][new_row].grid_polymer_unit_pointer = polymer_unit_iter;
		check_anchor_diffusion(polymer_unit_iter->polymer_grid_pointer->grid_anchor_pointer);
	}
	if (direction == direction_right)
	{
		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = (--polymer_cluster_iter->polymer_sequence.end());
		int end_col = polymer_unit_iter->polymer_grid_pointer->col_position;
		int end_row = polymer_unit_iter->polymer_grid_pointer->row_position;
		polymer_cluster_iter->polymer_sequence.push_back(*this);
		polymer_unit_iter = (--polymer_cluster_iter->polymer_sequence.end());
		//Polymer_Unit polymer_unit;
		polymer_unit_iter->attatch_mark = 0;
		polymer_unit_iter->Hydrolysis_mark = T;
		//when and where to check the boundary conditions?
		//now i need to consider periodic condition
		int new_col, new_row = 0;
		new_col = end_col;
		new_row = end_row+1;
		if (new_row > rule_parameters.row_num - 1)
			new_row = 0;

		rule_parameters.polymerize_statistics[new_row][direction_right]++;

		if (&*Grid[new_col][new_row].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
		{
			cout << "error: try to do polymerize when there is already a polymer_unit there" << endl;
			system("pause");
		}
		polymer_unit_iter->polymer_grid_pointer = &Grid[new_col][new_row];
		polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;

		polymer_unit_iter = --polymer_cluster_iter->polymer_sequence.end();
		Grid[new_col][new_row].grid_polymer_unit_pointer = polymer_unit_iter;
		check_anchor_diffusion(polymer_unit_iter->polymer_grid_pointer->grid_anchor_pointer);
	}
	
};


//Polymer_Unit::~Polymer_Unit()
//{//this is to make polymer_unit automaticlly remove the links between neighbor
//	//how do i know if it is in a pair or not,maybe use a map?
//	//can't use map, have to do it the old way: for loop iterator
//
// 	//if (this->polymer_cluster_pointer != rule_structure.polymer_cluster_list.end())
//	if(this->polymer_grid_pointer!=NULL)
//	{
//
//		
//		//first remove it's neighbor pointer
//		int col_position = this->polymer_grid_pointer->col_position;
//		int row_position = this->polymer_grid_pointer->row_position;
//		int up_col_position = col_position + 1;
//		if (up_col_position > column_num - 1)
//			up_col_position = 0;
//		int down_col_position = col_position - 1;
//		if (down_col_position < 0)
//			down_col_position = column_num - 1;
//		int right_row_position = row_position + 1;
//		if (right_row_position > row_num - 1)
//			right_row_position = 0;
//		int left_row_position = row_position - 1;
//		if (left_row_position < 0)
//			left_row_position = row_num - 1;
//
//		////up case//i should check those in depolymerize function, not here
//		//if (Grid[up_col_position][row_position].grid_polymer_unit_pointer != rule_structure.polymer_unit_list_end.begin()&&
//		//	Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer!=this->polymer_cluster_pointer)
//		//{//there is something there and it does not belong to the same polymer_cluster
//		//	//abandom to do all the cases diffusion stuff here, decide to use check_diffusion stuff to decide, might take long
//		//	//only consider no bundleing for depolymerize case for now
//		//	if (this->bundled_neighbor_polymer_unit_pointer[direction_up] != rule_structure.polymer_unit_list_end.begin())
//		//	{//no bundle from up, can be depolymerized
//
//		//	}
//
//		//}
//
//
//		//delete all pair list except debundling list
//		{
//			if (this->TT_anealing_polymer_unit_pair_list_iter_1 != rule_structure.TT_anealing_polymer_unit_pair_list.end())
//			{
//				rule_structure.TT_anealing_polymer_unit_pair_list.erase(this->TT_anealing_polymer_unit_pair_list_iter_1);
//			}
//			if (this->bundling_polymer_unit_pair_list_iter_1 != rule_structure.bundling_polymer_unit_pair_list.end())
//			{
//				rule_structure.bundling_polymer_unit_pair_list.erase(this->bundling_polymer_unit_pair_list_iter_1);
//			}
//			/*if (this->debundling_polymer_unit_pair_list_iter_1 != rule_structure.debundling_polymer_unit_pair_list.end())
//			{
//				rule_structure.debundling_polymer_unit_pair_list.erase(this->debundling_polymer_unit_pair_list_iter_1);
//			}*/
//			if (this->TT_anealing_polymer_unit_pair_list_iter_2 != rule_structure.TT_anealing_polymer_unit_pair_list.end())
//			{
//				rule_structure.TT_anealing_polymer_unit_pair_list.erase(this->TT_anealing_polymer_unit_pair_list_iter_2);
//			}
//			if (this->bundling_polymer_unit_pair_list_iter_2 != rule_structure.bundling_polymer_unit_pair_list.end())
//			{
//				rule_structure.bundling_polymer_unit_pair_list.erase(this->bundling_polymer_unit_pair_list_iter_2);
//			}
//			/*if (this->debundling_polymer_unit_pair_list_iter_2 != rule_structure.debundling_polymer_unit_pair_list.end())
//			{
//				rule_structure.debundling_polymer_unit_pair_list.erase(this->debundling_polymer_unit_pair_list_iter_2);
//			}*/
//		}
//
//
//		this->polymer_grid_pointer->grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
//
//		//delete its relation with bundle//there is no bundle to unit or unit to bundle iterator, nothing to be done
//		//delete its relation with anchor//a polymer with anchor will not be depolymerized,it will more likely to be framentation
//		
//		list<Polymer_Bundle>::iterator bundle_iter_up;
//		if (Grid[up_col_position][row_position].grid_polymer_unit_pointer->action_mark!=true && 
//			&*Grid[up_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
//			bundle_iter_up = Grid[up_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
//		else
//			bundle_iter_up = rule_structure.polymer_bundle_list.end();
//
//
//		list<Polymer_Bundle>::iterator bundle_iter_down;
//		if (Grid[down_col_position][row_position].grid_polymer_unit_pointer->action_mark != true &&
//			&*Grid[down_col_position][row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
//			bundle_iter_down = Grid[down_col_position][row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
//		else
//			bundle_iter_down = rule_structure.polymer_bundle_list.end();
//
//
//		list<Polymer_Bundle>::iterator bundle_iter_left;
//		if (Grid[col_position][left_row_position].grid_polymer_unit_pointer->action_mark != true &&
//			&*Grid[col_position][left_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
//			bundle_iter_left = Grid[col_position][left_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
//		else
//			bundle_iter_left = rule_structure.polymer_bundle_list.end();
//
//
//		list<Polymer_Bundle>::iterator bundle_iter_right;
//		if (Grid[col_position][right_row_position].grid_polymer_unit_pointer->action_mark != true &&
//			&*Grid[col_position][right_row_position].grid_polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())
//			bundle_iter_right = Grid[col_position][right_row_position].grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
//		else
//			bundle_iter_right = rule_structure.polymer_bundle_list.end();
//
//		
//		
//
//		
//		if (bundle_iter_up != rule_structure.polymer_bundle_list.end()&&
//			bundle_iter_up != this->polymer_cluster_pointer->cluster_bundle_pointer)
//		{//if they do not belong to the same bundle, need to check bundle diffusion relation
//			check_polymer_bundle_diffusable_direction(bundle_iter_up);
//		}
//		if (bundle_iter_down != rule_structure.polymer_bundle_list.end() &&
//			bundle_iter_down != this->polymer_cluster_pointer->cluster_bundle_pointer &&
//			bundle_iter_up != bundle_iter_down)
//		{
//			check_polymer_bundle_diffusable_direction(bundle_iter_down);
//		}
//		if (bundle_iter_left != rule_structure.polymer_bundle_list.end() &&
//			bundle_iter_left != this->polymer_cluster_pointer->cluster_bundle_pointer &&
//			bundle_iter_left != bundle_iter_down && bundle_iter_left != bundle_iter_up)
//		{
//			check_polymer_bundle_diffusable_direction(bundle_iter_left);
//		}
//		if (bundle_iter_right != rule_structure.polymer_bundle_list.end() &&
//			bundle_iter_right != this->polymer_cluster_pointer->cluster_bundle_pointer &&
//			bundle_iter_right != bundle_iter_down && bundle_iter_right != bundle_iter_up&& bundle_iter_right != bundle_iter_left)
//		{
//			check_polymer_bundle_diffusable_direction(bundle_iter_right);
//		}
//		
//	}
//};