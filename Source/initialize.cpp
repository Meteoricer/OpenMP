 #include "initialize.h"
#include <cmath>
#include <random>
#include "Rule_Parameters.h"
#include "Rule_Structure.h"
#include "Grid_Unit.h"
//#include "Diffusion_Unit.h"
//#include "Anchor_Unit.h"
#include <list>

using namespace std;
//#include "parameters.h"
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;



void Initialize()
{
	Initialize_rule_base();
	Initialize_Grid();
	Initialize_Anchor();
	
	
}
void Initialize_Grid()
{
	Grid.resize(rule_parameters.column_num, vector<Grid_Unit>(rule_parameters.row_num));
	for (int i = 0; i < rule_parameters.column_num;i++)
		for (int j = 0; j < rule_parameters.row_num;j++)
		{
			Grid[i][j].col_position = i;
			Grid[i][j].row_position = j;
			Grid[i][j].grid_anchor_pointer = rule_structure.anchor_list.end();
			Grid[i][j].grid_polymer_unit_pointer= rule_structure.polymer_unit_list_end.begin();
			Grid[i][j].MinC = sin((double)j*3.141592653 / (double)rule_parameters.row_num);
		}
}

void Initialize_Anchor()
{
	/*//random_device rd;
	mt19937_64 gen(rd());*/
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
	uniform_int_distribution<> column(0, rule_parameters.column_num-1);
	uniform_int_distribution<> row(0, rule_parameters.row_num-1);
	int colt;
	int rowt;
	list<Anchor_Unit>::iterator anchor_initial_temp;
	anchor_initial_temp = rule_structure.anchor_list.begin();
	for (int i = 0; i < rule_parameters.Anchor_num; i++)//here I link the empty anchor list with Grid and initialized it.
	{
		/*
		colt = column(gen);
		rowt = row(gen);
		Grid[colt][rowt].Is_Anchored = yes;
		list<Grid_Unit>::iterator temp;
		temp = rule_structure.anchor_empty_list.begin();
		rule_structure.anchor_empty_list.insert(temp, Grid[colt][rowt]);
		temp--;
		Grid[colt][rowt].p_anchor_empty_list = temp;
		temp++;
		



		//check four directions
		//if it is not at boundary
		if (colt > 0 && rowt > 0)
		{
			if (Grid[colt - 1][rowt].Is_Anchored == yes)
			{
				//delete this point in the duffusion list.
			}
			if (Grid[colt - 1][rowt].Is_Anchored == no)
			{
				list<Grid_Unit>::iterator temp;
				temp = rule_structure.anchor_diffusion_list.begin();
				//i need to write a construct function for diffusion unit;
				//Diffusion_Unit diffusion_unit(Grid[colt][rowt], up);
				rule_structure.anchor_diffusion_list.insert(temp, diffusion_unit);
			}
		
		}
		*/	
		do
		{
			colt = column(gen);
			rowt = row(gen);
		} while (Grid[colt][rowt].grid_anchor_pointer != rule_structure.anchor_list.end());//&*Grid[colt][rowt].grid_polymer_unit_pointer!=&*rule_structure.polymer_unit_list_end.begin());//when there is something, roll again
		Anchor_Unit anchor_unit;
		anchor_unit.anchor_grid_pointer = &Grid[colt][rowt];
		anchor_unit.polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
		//anchor_unit.IsAnchored = no;
		//anchor_unit.anchor_no = i;
		rule_structure.anchor_list.insert(anchor_initial_temp, anchor_unit);
		anchor_initial_temp--;
		Grid[colt][rowt].grid_anchor_pointer = anchor_initial_temp;
		anchor_initial_temp->polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
		anchor_initial_temp->consider_mark = 100;
		//anchor_initial_temp->attatch_mark = 0;
		check_anchor_diffusion(anchor_initial_temp);
		
		//Grid[colt][rowt].Is_Anchored = yes;
		anchor_initial_temp++;
		//initialize diffusion anchor 
	}
	rule_parameters.empty_anchor_unit_num = rule_parameters.Anchor_num;
	//now need to consider the sums
	//empty_anchor_num = Anchor_num; don't need this, we can use list.size()
	

	
	
	//difusion_candidate_num = Anchor_num; I think we need to test if there are space between every anchor first;
	//cout << Grid.size();


	//now I need to modify the anchoring rate;
	rule_parameters.ps[rule_anchoring] = rule_parameters.empty_anchor_unit_num*rule_parameters.anchoring_rate;
	check_sumps();

	
}
void Initialize_rule_base()
{
	rule_parameters.ps.resize(num_of_rules);//17 is the current rules i have
	for (int i = 0; i < num_of_rules; i++)
	{
		rule_parameters.ps[i] = 0;
	}
	rule_parameters.count.resize(num_of_rules);//17 is the current rules i have
	for (int i = 0; i < num_of_rules; i++)
	{
		rule_parameters.count[i] = 0;
	}
	rule_parameters.time.resize(num_of_rules + 1);//21 is the current rules i have
	for (int i = 0; i < num_of_rules; i++)
	{
		rule_parameters.time[i] = 0;
	}
	rule_parameters.length_distribution.resize(rule_parameters.column_num);
	for (int i = 0; i < rule_parameters.column_num; i++)
	{
		rule_parameters.length_distribution[i] = 0;
	}
	rule_parameters.direction_record.resize(num_of_rules, vector<int>(4));

	rule_parameters.average_length_per_column.resize(rule_parameters.row_num);
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		rule_parameters.average_length_per_column[i] = 0;
	}
	rule_parameters.polymer_count.resize(rule_parameters.row_num);
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		rule_parameters.polymer_count[i] = 0;
	}
	rule_parameters.diffusion_statistics.resize(rule_parameters.row_num, vector<int>(4));
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			rule_parameters.diffusion_statistics[i][j] = 0;
		}

	}

	rule_parameters.polymerize_statistics.resize(rule_parameters.row_num,vector<int>(4));

	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			rule_parameters.polymerize_statistics[i][j] = 0;
		}		

	}

	rule_parameters.depolymerization_statistics.resize(rule_parameters.row_num);
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		
		rule_parameters.depolymerization_statistics[i] = 0;
		

	}

	rule_parameters.fragmentation_statistics.resize(rule_parameters.row_num);
	for (int i = 0; i < rule_parameters.row_num; i++)
	{

		rule_parameters.fragmentation_statistics[i] = 0;


	}

	

	rule_parameters.remove_statistics.resize(rule_parameters.row_num);
	for (int i = 0; i < rule_parameters.row_num; i++)
	{

		rule_parameters.remove_statistics[i] = 0;


	}
}