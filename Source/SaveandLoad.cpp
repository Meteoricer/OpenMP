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
#include <cmath>

using namespace std;

extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
class Polymer_Unit;
class Anchor_Unit;


template<typename Unit>
int get_position(list<Unit> *unit_list, typename list<Unit>::iterator unit_iter)
{
	//cout << "i'm here" << endl;
	//cout << "size:" << unit_list->size() << endl;
	//system("pause");
	int position = 0;
	if (unit_list->size() != 0 && unit_iter != unit_list->end())
		//if (false)
	{
		//cout << "lalala" << endl;
		//system("pause");
		auto iter = unit_list->begin();

		while (iter != unit_iter)
		{
			//cout << "position:" << position << endl;
			//system("pause");
			/*if (position>1000)
			{
				cout << "comparison error" << endl;
				system("pause");
			}*/
			position++;
			iter++;
		}
	}
	else
	{
		//cout << "lululu" << endl;
		//system("pause");
		position = -1;
	}
	//cout << "position:" << position << endl;
	return position;
};


//int get_position(list<Anchor_Unit> *unit_list, list<Anchor_Unit>::iterator unit_iter)
//{
//	cout << "i'm here" << endl;
//	cout << "size:" <<unit_list->size()<<"position:"<<unit_iter->anchor_grid_pointer->col_position<< endl;
//	system("pause");
//	int position = 0;
//	if (unit_list->size() != 0 || unit_iter != unit_list->end())
//	//if (false)
//	{
//		cout << "lalala" << endl;
//		system("pause");
//		auto iter = unit_list->begin();
//
//		while (iter != unit_iter)
//		{
//			cout << "position:"<<position << endl;
//			//system("pause");
//			if (iter->anchor_grid_pointer->col_position == unit_iter->anchor_grid_pointer->col_position&&
//				iter->anchor_grid_pointer->row_position == unit_iter->anchor_grid_pointer->row_position)
//			{
//				cout << "comparison error" << endl;
//				system("pause");
//			}
//			position++;
//			iter++;
//		}
//	}
//	else
//	{
//		cout << "lululu" << endl;
//		system("pause");
//		position = -1;
//	}
//	cout << "position:" << position << endl;
//	return position;
//};

void SaveAllParameters()
{
	ofstream savemark;
	savemark.open("savemark.txt");
	savemark << 1 << endl;
	savemark.close();
	ofstream parametersout;
	parametersout.open("parameters.txt");
	parametersout.unsetf(ios::scientific);
	//int empty_anchor_unit_num=0;
	parametersout << rule_parameters.empty_anchor_unit_num << endl;
	//for hydrolysis
	//int bundling_sum=0;
	parametersout << rule_parameters.bundling_sum << endl;
	//int bundling_inverse_sum=0;
	parametersout << rule_parameters.bundling_inverse_sum << endl;
	//int debundling_sum=0;
	parametersout << rule_parameters.debundling_sum << endl;
	//int debundling_inverse_sum=0;
	parametersout << rule_parameters.debundling_inverse_sum << endl;
	//int TTT_sum_total=0;
	parametersout << rule_parameters.TTT_sum_total << endl;
	//int TTD_sum_total;//TTD and DTT are the same here
	parametersout << rule_parameters.TTD_sum_total << endl;
	//int DTD_sum_total;
	parametersout << rule_parameters.DTD_sum_total << endl;
	//int TT_annealing_sum = 0;
	parametersout << rule_parameters.TT_annealing_sum << endl;
	//int TD_annealing_sum = 0;
	parametersout << rule_parameters.TD_annealing_sum << endl;
	//int TT_polymerize_end_vertical_plus_sum=0;
	parametersout << rule_parameters.TT_polymerize_end_vertical_plus_sum << endl;
	//int TT_polymerize_end_vertical_minus_sum=0;
	parametersout << rule_parameters.TT_polymerize_end_vertical_minus_sum << endl;
	//int TT_polymerize_end_horizontal_plus_sum = 0;
	parametersout << rule_parameters.TT_polymerize_end_horizontal_plus_sum << endl;
	//int TT_polymerize_end_horizontal_minus_sum = 0;
	parametersout << rule_parameters.TT_polymerize_end_horizontal_minus_sum << endl;
	//int TD_polymerize_end_vertical_plus_sum = 0;
	parametersout << rule_parameters.TD_polymerize_end_vertical_plus_sum << endl;
	//int TD_polymerize_end_vertical_minus_sum = 0;
	parametersout<<rule_parameters.TD_polymerize_end_vertical_minus_sum<<endl;
	//int TD_polymerize_end_horizontal_plus_sum = 0;
	parametersout << rule_parameters.TD_polymerize_end_horizontal_plus_sum <<endl;
	//int TD_polymerize_end_horizontal_minus_sum = 0;
	parametersout << rule_parameters.TD_polymerize_end_horizontal_minus_sum << endl;
	//int TT_depolymerize_plus_end_sum=0;
	parametersout << rule_parameters.TT_depolymerize_plus_end_sum << endl;
	//int TD_depolymerize_plus_end_sum=0;//TD and DT are the same here
	parametersout << rule_parameters.TD_depolymerize_plus_end_sum << endl;
	//int DD_depolymerize_plus_end_sum=0;
	parametersout << rule_parameters.DD_depolymerize_plus_end_sum << endl;
	//int TT_depolymerize_minus_end_sum=0;
	parametersout << rule_parameters.TT_depolymerize_minus_end_sum << endl;
	//int TD_depolymerize_minus_end_sum=0;//TD and DT are the same here
	parametersout << rule_parameters.TD_depolymerize_minus_end_sum << endl;
	//int DD_depolymerize_minus_end_sum=0;
	parametersout << rule_parameters.DD_depolymerize_minus_end_sum << endl;
	//int polymer_diffusion_propensity;
	parametersout << rule_parameters.polymer_diffusion_propensity << endl;
	//int empty_anchor_diffusion_direction_sum;
	parametersout << rule_parameters.empty_anchor_diffusion_direction_sum << endl;
	//int emtpy_anchor_attatch_sum;
	parametersout << rule_parameters.emtpy_anchor_attatch_sum << endl;
	//int TT_fragmentation_sum;
	parametersout << rule_parameters.TT_fragmentation_sum << endl;
	//int TD_fragmentation_sum;
	parametersout << rule_parameters.TD_fragmentation_sum << endl;
	//int DD_fragmentation_sum;
	parametersout << rule_parameters.DD_fragmentation_sum << endl;
	
	//int rotation_sum;
	parametersout << rule_parameters.rotation_sum << endl;
	//double TD_case = 0;//T------D
	parametersout << rule_parameters.TD_case << endl;
	//double DT_case = 0;//D------T
	parametersout << rule_parameters.DT_case << endl;
	//double TT_case = 0;//T------T
	parametersout << rule_parameters.TT_case << endl;
	//double DD_case = 0;//D------D
	parametersout << rule_parameters.DD_case << endl;
	//int T_sum;
	parametersout << rule_parameters.T_sum << endl;
	//int D_sum;
	parametersout << rule_parameters.D_sum << endl;
	//int TT_annealing_sum=0;
	//int step_num;
	parametersout << rule_parameters.step_num << endl;
	//int Z_ctp = 0;
	parametersout << rule_parameters.Z_ctp << endl;
	//int reaction_mark;
	parametersout << rule_parameters.reaction_mark << endl;
	//double cumps;
	parametersout << rule_parameters.cumps << endl;
	//double sumps;
	parametersout << rule_parameters.sumps << endl;
	//double t;
	parametersout << rule_parameters.t << endl;
	//double polymer_sequence_average_length = 0;
	parametersout << rule_parameters.polymer_sequence_average_length << endl;
	//int cluster_bound = 0;
	parametersout << rule_parameters.cluster_bound << endl;
	//int up_larger_than_down_annealing = 0;
	parametersout << rule_parameters.up_larger_than_down_annealing << endl;
	//int up_larger_than_down_fragmentation = 0;
	parametersout << rule_parameters.up_larger_than_down_fragmentation << endl;
	//int diffusion_count = 0;
	parametersout << rule_parameters.diffusion_count << endl;
	//int TT_annealing_count_up_larger_than_below = 0;
	parametersout << rule_parameters.TT_annealing_count_up_larger_than_below << endl;
	//int TT_annealing_count_up_smaller_than_below = 0;
	parametersout << rule_parameters.TT_annealing_count_up_smaller_than_below << endl;
	//int TD_annealing_count_up_larger_than_below = 0;
	parametersout << rule_parameters.TD_annealing_count_up_larger_than_below << endl;
	//int TD_annealing_count_up_smaller_than_below = 0;
	parametersout << rule_parameters.TD_annealing_count_up_smaller_than_below << endl;
	//int DT_annealing_count_up_larger_than_below = 0;
	parametersout << rule_parameters.DT_annealing_count_up_larger_than_below << endl;
	//int DT_annealing_count_up_smaller_than_below = 0;
	parametersout << rule_parameters.DT_annealing_count_up_smaller_than_below << endl;
	//int debundle_pair_TT;
	parametersout << rule_parameters.debundle_pair_TT << endl;
	//int debundle_pair_TD;
	parametersout << rule_parameters.debundle_pair_TD << endl;
	//int debundle_pair_DD;
	parametersout << rule_parameters.debundle_pair_DD << endl;
	//int total_ftsZ_count = 0;
	parametersout << rule_parameters.total_ftsZ_count << endl;
	//double total_ftsZ_lifetime = 0;
	parametersout << rule_parameters.total_ftsZ_lifetime << endl;
	
	//vector<vector<int>> direction_record;
	for (int i = 0; i < num_of_rules; i++)
	{
		
		for (int j = 0; j < 4; j++)
		{
			parametersout << rule_parameters.direction_record[i][j] << " ";
		}
		parametersout << endl;
	}
	//vector<double> ps;
	for (int i = 0; i < rule_parameters.ps.size(); i++)
	{
		parametersout << rule_parameters.ps[i] << endl;
	}
	//vector<double> time;
	for (int i = 0; i < rule_parameters.time.size(); i++)
	{
		parametersout << rule_parameters.time[i] << endl;
	}
	//vector<double> count;
	for (int i = 0; i < rule_parameters.count.size(); i++)
	{
		parametersout << rule_parameters.count[i] << endl;
	}

	//int typical_column = 600;
	parametersout << rule_parameters.typical_column << endl;
	//int typical_row = 800;
	parametersout << rule_parameters.typical_row << endl;


	//int column_num = 600;//600//150
	parametersout << rule_parameters.column_num << endl;
	//int row_num = 563;//800//200
	parametersout << rule_parameters.row_num << endl;
	//int Anchor_num = 1000;//1000
	parametersout << rule_parameters.Anchor_num << endl;
	//int Z = 8000;//8000//600//1000//4000
	parametersout << rule_parameters.Z << endl;

	//double bending_energy = 1;
	parametersout << rule_parameters.bending_energy << endl;
	//double scale_factor = typical_row / (double)row_num*typical_column*typical_column / column_num / column_num;
	parametersout << rule_parameters.scale_factor << endl;

	//double anchoring_rate = 0.005*scale_factor;//0.0005;
	parametersout << rule_parameters.anchoring_rate << endl;
	//double deanchoring_rate = 10;//8;//4;//0.2;
	parametersout << rule_parameters.deanchoring_rate << endl;
	//double TT_polymerize_vertical_plus_rate = 0.01*scale_factor;//0.03;//0.05
	parametersout << rule_parameters.TT_polymerize_horizontal_plus_rate << endl;
	//double TT_polymerize_vertical_minus_rate = 0 * scale_factor;//0.01;//0.05
	parametersout << rule_parameters.TT_polymerize_vertical_minus_rate << endl;
	//double TT_polymerize_horizontal_plus_rate = TT_polymerize_vertical_plus_rate*exp(-bending_energy);
	parametersout << rule_parameters.TT_polymerize_horizontal_plus_rate << endl;
	//double TT_polymerize_horizontal_minus_rate = TT_polymerize_vertical_minus_rate*exp(-bending_energy);
	parametersout << rule_parameters.TT_polymerize_horizontal_minus_rate << endl;
	//double TD_polymerize_vertical_plus_rate = 0.005*scale_factor;//0.05
	parametersout << rule_parameters.TD_polymerize_vertical_plus_rate << endl;
	//double TD_polymerize_vertical_minus_rate = 0 * scale_factor;//0.05
	parametersout << rule_parameters.TD_polymerize_vertical_minus_rate << endl;
	//double TD_polymerize_horizontal_plus_rate = TD_polymerize_vertical_plus_rate*exp(-bending_energy);
	parametersout << rule_parameters.TD_polymerize_horizontal_plus_rate << endl;
	//double TD_polymerize_horizontal_minus_rate = TD_polymerize_vertical_minus_rate*exp(-bending_energy);
	parametersout << rule_parameters.TD_polymerize_horizontal_minus_rate << endl;
	//double diffusion_rate = 4000;//1000;//4000
	parametersout << rule_parameters.diffusion_rate << endl;
	//double emtpy_anchor_diffusion_rate = 4000;//1000;//0.1;
	parametersout << rule_parameters.emtpy_anchor_diffusion_rate << endl;
	//double Ebdl = 0.05;//0.03//0.07
	parametersout << rule_parameters.Ebdl << endl;
	//double Ebdl_inverse = 0;
	parametersout << rule_parameters.Ebdl_inverse << endl;
	//double bundling_rate = 2.5;//50;//100;//500
	parametersout << rule_parameters.bundling_rate << endl;
	//double bundling_inverse_rate = 0;
	parametersout << rule_parameters.bundling_inverse_rate << endl;
	//double debundling_rate = bundling_rate*exp(-Ebdl);
	parametersout << rule_parameters.debundling_rate << endl;
	//double debundling_inverse_rate = bundling_inverse_rate*exp(-Ebdl_inverse);
	parametersout << rule_parameters.debundling_inverse_rate << endl;
	//double TT_annealing_rate = 100;//2500;
	parametersout << rule_parameters.TT_annealing_rate << endl;
	//double TD_annealing_rate = 10;
	parametersout << rule_parameters.TD_annealing_rate << endl;
	//double TT_fragmentation_rate = 0.001;//1;//8; #
	parametersout << rule_parameters.TT_fragmentation_rate << endl;
	//double TD_fragmentation_rate = 0.0001;//3;// 8;#
	parametersout << rule_parameters.TD_fragmentation_rate << endl;
	//double DD_fragmentation_rate = 0.00001;//5;//8;#5
	parametersout << rule_parameters.DD_fragmentation_rate << endl;
	//double TT_depolymerize_minus_rate = 0.1;//1;//8;//80#1
	parametersout << rule_parameters.TT_depolymerize_minus_rate << endl;
	//double TD_depolymerize_minus_rate = TT_depolymerize_minus_rate * 40;//5;//20;//#5
	parametersout << rule_parameters.TD_depolymerize_minus_rate << endl;
	//double DD_depolymerize_minus_rate = TD_depolymerize_minus_rate * 40;//10;//40;//#5
	parametersout << rule_parameters.DD_depolymerize_minus_rate << endl;
	//double TT_depolymerize_plus_rate = TT_depolymerize_minus_rate / 10;//1;//8;//80#1
	parametersout << rule_parameters.TT_depolymerize_plus_rate << endl;
	//double TD_depolymerize_plus_rate = TT_depolymerize_plus_rate * 40;//5;//20;//#5
	parametersout << rule_parameters.TD_depolymerize_plus_rate << endl;
	//double DD_depolymerize_plus_rate = TD_depolymerize_plus_rate * 40;//10;//40;//#5
	parametersout << rule_parameters.DD_depolymerize_plus_rate << endl;
	//double depolymerize_rate = 20;
	parametersout << rule_parameters.depolymerize_rate << endl;

	//double attatch_rate = 400;//100;//100//0.000001;
	parametersout << rule_parameters.attatch_rate << endl;
	//double both = 200;

	//double hydrolysis_rate = 0.0156;
	parametersout << rule_parameters.hydrolysis_rate << endl;
	//double TTT_hydrolysys_rate = 0.0985;//0.5;//2;//10;//100;//100;
	parametersout << rule_parameters.hydrolysis_multiplier << endl;
	parametersout << rule_parameters.TTT_hydrolysys_rate << endl;
	//double TTD_hydrolysys_rate = TTT_hydrolysys_rate * 5;//1;//20;//20;// 100;//100;
	parametersout << rule_parameters.TTD_hydrolysys_rate << endl;
	//double DTD_hydrolysys_rate = TTD_hydrolysys_rate * 5;//2;//40;//10;// 100;//100;
	parametersout << rule_parameters.DTD_hydrolysys_rate << endl;
	//double rotation_rate = 0;
	parametersout << rule_parameters.rotation_rate << endl;

	//double posibility_vertical = 1 / (1 + exp(-bending_energy));
	parametersout << rule_parameters.posibility_vertical << endl;
	//double posibility_horizontal = exp(-bending_energy) / (1 + exp(-bending_energy));
	parametersout << rule_parameters.posibility_horizontal << endl;
	parametersout.close();
};
void SaveAllObject()
{
	//save polymer_bundle object
	//bundle id
	//bundle cluster pointer list
	//bundle anchor pointer list
	//bundle direction
	ofstream polymer_bundle_out;
	polymer_bundle_out.open("polymer_bundle.txt");
	polymer_bundle_out << rule_structure.polymer_bundle_list.size() << endl;
	int bundle_id = 0;
	auto polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
	while(polymer_bundle_iter != rule_structure.polymer_bundle_list.end())
	{		
		polymer_bundle_out << bundle_id << endl;
		//this output line is for bundle_cluster_pointer_list
		if (polymer_bundle_iter->bundle_cluster_pointer_list.size() != 0)
		{
			polymer_bundle_out << polymer_bundle_iter->bundle_cluster_pointer_list.size() << " ";
			auto bundle_cluster_pointer_list_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
			int position = 0;
			while (bundle_cluster_pointer_list_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
			{
				position = 0;
				auto position_iter = rule_structure.polymer_cluster_list.begin();
				while (position_iter != *bundle_cluster_pointer_list_iter)
				{
					position++;
					position_iter++;
				}
				polymer_bundle_out << position << " ";

				bundle_cluster_pointer_list_iter++;
			}
			polymer_bundle_out << endl;
		}
		else
		{
			polymer_bundle_out << polymer_bundle_iter->bundle_cluster_pointer_list.size() << endl;
		}
		
		//this output line is for bundle_anchor_list
		if (polymer_bundle_iter->polymer_anchor_list.size() != 0)
		{
			polymer_bundle_out << polymer_bundle_iter->polymer_anchor_list.size() << " ";
			auto polymer_anchor_list_iter = polymer_bundle_iter->polymer_anchor_list.begin();
			int position = 0;
			while (polymer_anchor_list_iter != polymer_bundle_iter->polymer_anchor_list.end())
			{
				position = 0;
				auto position_iter = rule_structure.anchor_list.begin();
				while (position_iter != *polymer_anchor_list_iter)
				{
					position++;
					position_iter++;
				}
				polymer_bundle_out << position << " ";

				polymer_anchor_list_iter++;
			}
			polymer_bundle_out << endl;
		}
		else
		{
			polymer_bundle_out << polymer_bundle_iter->polymer_anchor_list.size() << endl;
		}
		//this output line is for availabe_diffusion_direction
		for (int i = direction_up; i <= direction_right; i++)
		{
			polymer_bundle_out << polymer_bundle_iter->availabe_diffusion_direction[i] << " ";
		}
		polymer_bundle_out << endl;
		polymer_bundle_out << polymer_bundle_iter->available_polymer_diffusion_direction_sum << endl;
		polymer_bundle_out << polymer_bundle_iter->Consider_mark << endl;
		polymer_bundle_out << polymer_bundle_iter->bundle_direction << endl;
		bundle_id++;
		polymer_bundle_iter++;
	}
	polymer_bundle_out.close();


	//save polymer_cluster object
	//cluster id
	//anchor_pointer_list
	//direction
	//polarize
	//polymer_sequence
	//TTT_sum
	//cluster bundle pointer
	//bundling_list
	//debundling_list
	//TTD_sum
	//DTD_sum
	//TT_sum
	//TD_sum
	//DD_sum
	//TT_depolymerize_plus_end_sum = 0;
	//TD_depolymerize_plus_end_sum = 0;//TD and DT are the same here
	//DD_depolymerize_plus_end_sum = 0;
	//TT_depolymerize_minus_end_sum = 0;
	//TD_depolymerize_minus_end_sum = 0;//TD and DT are the same here
	//DD_depolymerize_minus_end_sum = 0;
	//Consider_mark = 0;
	//ring_mark = 0;


	//TT_polymerize_plus_end_sum = 0;//directions
	//TT_polymerize_minus_end_sum = 0;//directions
	//TD_polymerize_plus_end_sum = 0;//directions
	//TD_polymerize_minus_end_sum = 0;//directions
	//vector<int> availabe_polymerize_direction;





	ofstream polymer_unit_out;
	polymer_unit_out.open("polymer_unit.txt");
	
	ofstream polymer_cluster_out;
	polymer_cluster_out.open("polymer_cluster.txt");
	polymer_cluster_out << rule_structure.polymer_cluster_list.size() << endl;
	int cluster_id = 0;
	auto polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{
		polymer_cluster_out << cluster_id << endl;
		if (polymer_cluster_iter->anchor_pointer_list.size() != 0)
		{
			polymer_cluster_out << polymer_cluster_iter->anchor_pointer_list.size() << " ";
			auto cluster_anchor_pointer_list_iter = polymer_cluster_iter->anchor_pointer_list.begin();
			int position = 0;
			while (cluster_anchor_pointer_list_iter != polymer_cluster_iter->anchor_pointer_list.end())
			{
				position = 0;
				auto position_iter = rule_structure.anchor_list.begin();
				while (position_iter != *cluster_anchor_pointer_list_iter)
				{
					position++;
					position_iter++;
				}
				polymer_cluster_out << position << " ";

				cluster_anchor_pointer_list_iter++;
			}
			polymer_cluster_out << endl;
		}
		else
		{
			polymer_cluster_out << polymer_cluster_iter->anchor_pointer_list.size()<< endl;
		}
		polymer_cluster_out << polymer_cluster_iter->direction << endl;
		polymer_cluster_out << polymer_cluster_iter->polarize << endl;

		//this part is for polymer_sequence
		
		//polymer_cluster_out << "start" << endl;//do i really need this one?
		polymer_cluster_out << polymer_cluster_iter->polymer_sequence.size() << endl;
		if (polymer_cluster_iter->polymer_sequence.size() != 0)
		{
			auto polymer_sequence_iter = polymer_cluster_iter->polymer_sequence.begin();
			while (polymer_sequence_iter != polymer_cluster_iter->polymer_sequence.end())
			{
				//this is where to output polymer_unit line by line
				
				polymer_unit_out << get_position(&(rule_structure.anchor_list), polymer_sequence_iter->anchor_pointer)<<" ";//anchor pointer;
				polymer_unit_out << polymer_sequence_iter->direction << " ";
				polymer_unit_out << polymer_sequence_iter->polymer_grid_pointer->col_position << " " << polymer_sequence_iter->polymer_grid_pointer->row_position << " ";
				polymer_unit_out << polymer_sequence_iter->Hydrolysis_mark << " ";
				polymer_unit_out << polymer_sequence_iter->Consider_mark << " ";
				polymer_unit_out << polymer_sequence_iter->lifetime << " ";
				polymer_unit_out << polymer_sequence_iter->action_mark << " ";
				polymer_unit_out << polymer_sequence_iter->attatch_mark << " ";
				polymer_unit_out << polymer_sequence_iter->frap_mark << " ";
				//polymer_unit_out << get_position(rule_structure.polymer_cluster_list,polymer_sequence_iter->polymer_cluster_pointer) << " ";
				for (int i = direction_first; i <= direction_second; i++)
				{
					polymer_unit_out << get_position(&rule_structure.anealing_polymer_unit_pair_list,polymer_sequence_iter->anealing_polymer_unit_pair_list_iter[i]) << " ";
				}
				for (int i = direction_first; i <= direction_second; i++)
				{
					polymer_unit_out <<  get_position(&rule_structure.bundling_polymer_unit_pair_list, polymer_sequence_iter->bundling_polymer_unit_pair_list_iter[i]) << " ";
				}
				for (int i = direction_first; i <= direction_second; i++)
				{
					polymer_unit_out <<  get_position(&rule_structure.debundling_polymer_unit_pair_list, polymer_sequence_iter->debundling_polymer_unit_pair_list_iter[i]) << " ";
				}
				polymer_unit_out << endl;
				polymer_sequence_iter++;
			}
			
		}
		else
		{
			cout << "cluster_sequence size 0, there is a error";
			system("pause");
		}
		//polymer_cluster_out << "end" << endl;
		
		polymer_cluster_out << get_position(&rule_structure.polymer_bundle_list, polymer_cluster_iter->cluster_bundle_pointer) << endl;;


		//bundling_list
		//debundling_list
		//i'm not sure what does those used for

		
		polymer_cluster_out << polymer_cluster_iter->TTT_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TTD_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->DTD_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TT_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TD_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->DD_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TT_depolymerize_plus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TD_depolymerize_plus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->DD_depolymerize_plus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TT_depolymerize_minus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TD_depolymerize_minus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->DD_depolymerize_minus_end_sum << endl;
		
		polymer_cluster_out << polymer_cluster_iter->Consider_mark << endl;
		polymer_cluster_out << polymer_cluster_iter->ring_mark << endl;
		polymer_cluster_out << polymer_cluster_iter->TT_polymerize_plus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter->TT_polymerize_minus_end_sum << endl;
		polymer_cluster_out << polymer_cluster_iter-> TD_polymerize_plus_end_sum<< endl;
		polymer_cluster_out << polymer_cluster_iter-> TD_polymerize_minus_end_sum<< endl;


		for (int i = direction_first; i <= direction_right; i++)
		{
			polymer_cluster_out << polymer_cluster_iter->availabe_polymerize_direction[i] << " ";
		}
		polymer_cluster_out << endl;
		cluster_id++;
		polymer_cluster_iter++;
		
		
	}
	polymer_cluster_out.close();
	polymer_unit_out.close();

	ofstream anchor_unit_out;
	anchor_unit_out.open("anchor_unit.txt");
	int anchor_id = 0;
	auto anchor_iter = rule_structure.anchor_list.begin();
	while (anchor_iter != rule_structure.anchor_list.end())
	{
		anchor_unit_out << anchor_id << endl;
		
		//anchor_grid_pointer
		anchor_unit_out << anchor_iter->anchor_grid_pointer->col_position << " " << anchor_iter->anchor_grid_pointer->row_position << endl;
		anchor_unit_out << anchor_iter->frap_mark << endl;
		//this output line is for availabe_diffusion_direction
		for (int i = direction_up; i <= direction_right; i++)
		{
			anchor_unit_out<<anchor_iter->availabe_diffusion_direction[i] << " ";
		}
		anchor_unit_out << endl;
		anchor_unit_out << anchor_iter->available_polymer_diffusion_direction_sum << endl;
		//save polymer_unit_pointer on polymer_unit side
		anchor_unit_out << anchor_iter->consider_mark << endl;
		anchor_unit_out << anchor_iter->attatch_mark << endl;
		
		anchor_id++;
		anchor_iter++;
	}
	anchor_unit_out.close();




	//how to save grid unit?
	ofstream grid_unit_out;
	grid_unit_out.open("grid_unit.txt");
	int grid_unit_id = 0;
	for (int i = 0; i < rule_parameters.column_num; i++)
	{
		for (int j = 0; j < rule_parameters.row_num; j++)
		{
			grid_unit_out << grid_unit_id << " ";
			grid_unit_out << Grid[i][j].col_position << " ";
			grid_unit_out << Grid[i][j].row_position << " ";
			grid_unit_out << get_position(&rule_structure.anchor_list, Grid[i][j].grid_anchor_pointer);
			//need to reconstruct grid_polymer_unit_pointer 
			
			grid_unit_out << endl;
			grid_unit_id++;
		}
		
	}
	grid_unit_out.close();

	//annealing_polymer_unit_pair_list
	//bundling_polymer_unit_pair_list_iter
	//debundling_polymer_unit_pair_list_iter
	//these three do not need to save and can be reconstruct by exist information
	ofstream pair_out;
	pair_out.open("unit_pair.txt");
	pair_out << rule_structure.anealing_polymer_unit_pair_list.size() << endl;
	pair_out << rule_structure.bundling_polymer_unit_pair_list.size()<< endl;
	pair_out << rule_structure.debundling_polymer_unit_pair_list.size()<< endl;
	pair_out.close();
	
};



void LoadRuleParameters()
{
	ifstream parametersin;
	parametersin.open("rule_parameters.txt");

	//int typical_column = 600;
	parametersin >> rule_parameters.typical_column;
	//int typical_row = 800;
	parametersin >> rule_parameters.typical_row;


	//int column_num = 600;//600//150
	parametersin >> rule_parameters.column_num;
	//int row_num = 563;//800//200
	parametersin >> rule_parameters.row_num;
	//int Anchor_num = 1000;//1000
	parametersin >> rule_parameters.Anchor_num;
	//int Z = 8000;//8000//600//1000//4000
	parametersin >> rule_parameters.Z;

	//double bending_energy = 1;
	parametersin >> rule_parameters.bending_energy;
	rule_parameters.scale_factor = rule_parameters.typical_row / (double)rule_parameters.row_num*rule_parameters.typical_column*rule_parameters.typical_column / rule_parameters.column_num / rule_parameters.column_num;
	//parametersin >> rule_parameters.scale_factor;

	//double anchoring_rate = 0.005*scale_factor;//0.0005;
	parametersin >> rule_parameters.anchoring_rate;
	rule_parameters.anchoring_rate = rule_parameters.anchoring_rate*rule_parameters.scale_factor;
	//double deanchoring_rate = 10;//8;//4;//0.2;
	parametersin >> rule_parameters.deanchoring_rate;
	//double TT_polymerize_vertical_plus_rate = 0.01*scale_factor;//0.03;//0.05
	parametersin >> rule_parameters.TT_polymerize_horizontal_plus_rate;
	rule_parameters.TT_polymerize_horizontal_plus_rate = rule_parameters.TT_polymerize_horizontal_plus_rate*rule_parameters.scale_factor;
	//double TT_polymerize_vertical_minus_rate = 0 * scale_factor;//0.01;//0.05
	parametersin >> rule_parameters.TT_polymerize_vertical_minus_rate;
	rule_parameters.TT_polymerize_vertical_minus_rate = rule_parameters.TT_polymerize_vertical_minus_rate*rule_parameters.scale_factor;
	//double TT_polymerize_horizontal_plus_rate = TT_polymerize_vertical_plus_rate*exp(-bending_energy);
	//parametersin >> rule_parameters.TT_polymerize_horizontal_plus_rate;
	rule_parameters.TT_polymerize_horizontal_plus_rate = rule_parameters.TT_polymerize_vertical_plus_rate*exp(-rule_parameters.bending_energy);
	//double TT_polymerize_horizontal_minus_rate = TT_polymerize_vertical_minus_rate*exp(-bending_energy);
	//parametersin >> rule_parameters.TT_polymerize_horizontal_minus_rate;
	rule_parameters.TT_polymerize_horizontal_minus_rate=rule_parameters.TT_polymerize_vertical_minus_rate*exp(-rule_parameters.bending_energy);
	//double TD_polymerize_vertical_plus_rate = 0.005*scale_factor;//0.05
	parametersin >> rule_parameters.TD_polymerize_vertical_plus_rate;
	rule_parameters.TD_polymerize_vertical_plus_rate = rule_parameters.TD_polymerize_vertical_plus_rate*rule_parameters.scale_factor;
	//double TD_polymerize_vertical_minus_rate = 0 * scale_factor;//0.05
	parametersin >> rule_parameters.TD_polymerize_vertical_minus_rate;
	rule_parameters.TD_polymerize_vertical_minus_rate= rule_parameters.TD_polymerize_vertical_minus_rate*rule_parameters.scale_factor;
	//double TD_polymerize_horizontal_plus_rate = TD_polymerize_vertical_plus_rate*exp(-bending_energy);
	//parametersin >> rule_parameters.TD_polymerize_horizontal_plus_rate;
	rule_parameters.TD_polymerize_horizontal_plus_rate=rule_parameters.TD_polymerize_vertical_plus_rate*exp(-rule_parameters.bending_energy);
	//double TD_polymerize_horizontal_minus_rate = TD_polymerize_vertical_minus_rate*exp(-bending_energy);
	//parametersin >> rule_parameters.TD_polymerize_horizontal_minus_rate;
	rule_parameters.TD_polymerize_horizontal_minus_rate=rule_parameters.TD_polymerize_vertical_minus_rate*exp(-rule_parameters.bending_energy);
	//double diffusion_rate = 4000;//1000;//4000
	parametersin >> rule_parameters.diffusion_rate;
	//double emtpy_anchor_diffusion_rate = 4000;//1000;//0.1;
	parametersin >> rule_parameters.emtpy_anchor_diffusion_rate;
	//double Ebdl = 0.05;//0.03//0.07
	parametersin >> rule_parameters.Ebdl;
	//double Ebdl_inverse = 0;
	parametersin >> rule_parameters.Ebdl_inverse;
	//double bundling_rate = 2.5;//50;//100;//500
	parametersin >> rule_parameters.bundling_rate;
	//double bundling_inverse_rate = 0;
	parametersin >> rule_parameters.bundling_inverse_rate;
	//double debundling_rate = bundling_rate*exp(-Ebdl);
	//parametersin >> rule_parameters.debundling_rate;
	rule_parameters.debundling_rate=rule_parameters.bundling_rate*exp(-rule_parameters.Ebdl);
	//double debundling_inverse_rate = bundling_inverse_rate*exp(-Ebdl_inverse);
	//parametersin >> rule_parameters.debundling_inverse_rate;
	rule_parameters.debundling_inverse_rate=rule_parameters.bundling_inverse_rate*exp(-rule_parameters.Ebdl_inverse);
	//double TT_annealing_rate = 100;//2500;
	parametersin >> rule_parameters.TT_annealing_rate;
	//double TD_annealing_rate = 10;
	parametersin >> rule_parameters.TD_annealing_rate;
	//double TT_fragmentation_rate = 0.001;//1;//8; #
	parametersin >> rule_parameters.fragmentation_multiplier;
	parametersin >> rule_parameters.TT_fragmentation_rate;
	//double TD_fragmentation_rate = 0.0001;//3;// 8;#
	rule_parameters.TD_fragmentation_rate=rule_parameters.TT_fragmentation_rate*rule_parameters.fragmentation_multiplier;
	//double DD_fragmentation_rate = 0.00001;//5;//8;#5
	rule_parameters.DD_fragmentation_rate=rule_parameters.TD_fragmentation_rate*rule_parameters.fragmentation_multiplier;
	//double TT_depolymerize_minus_rate = 0.1;//1;//8;//80#1
	rule_parameters.depolymerize_multiplier = rule_parameters.fragmentation_multiplier;
	parametersin >> rule_parameters.TT_depolymerize_minus_rate;
	//double TD_depolymerize_minus_rate = TT_depolymerize_minus_rate * 40;//5;//20;//#5
	//parametersin >> rule_parameters.TD_depolymerize_minus_rate;
	rule_parameters.TD_depolymerize_minus_rate = rule_parameters.TT_depolymerize_minus_rate*rule_parameters.fragmentation_multiplier;
	//double DD_depolymerize_minus_rate = TD_depolymerize_minus_rate * 40;//10;//40;//#5
	rule_parameters.DD_depolymerize_minus_rate=rule_parameters.TD_depolymerize_minus_rate*rule_parameters.depolymerize_multiplier;
	parametersin >> rule_parameters.TT_depolymerize_plus_rate;
	//double TD_depolymerize_minus_rate = TT_depolymerize_minus_rate * 40;//5;//20;//#5
	//parametersin >> rule_parameters.TD_depolymerize_minus_rate;
	rule_parameters.TD_depolymerize_plus_rate = rule_parameters.TT_depolymerize_plus_rate*rule_parameters.fragmentation_multiplier;
	//double DD_depolymerize_minus_rate = TD_depolymerize_minus_rate * 40;//10;//40;//#5
	rule_parameters.DD_depolymerize_plus_rate = rule_parameters.TD_depolymerize_plus_rate*rule_parameters.depolymerize_multiplier;
	//double depolymerize_rate = 20;
	parametersin >> rule_parameters.depolymerize_rate;

	//double attatch_rate = 400;//100;//100//0.000001;
	parametersin >> rule_parameters.attatch_rate;
	//double both = 200;

	//double hydrolysis_rate = 0.0156;
	parametersin >> rule_parameters.hydrolysis_rate;
	//double TTT_hydrolysys_rate = 0.0985;//0.5;//2;//10;//100;//100;
	parametersin >> rule_parameters.hydrolysis_multiplier;
	rule_parameters.TTT_hydrolysys_rate=rule_parameters.hydrolysis_rate;
	//double TTD_hydrolysys_rate = TTT_hydrolysys_rate * 5;//1;//20;//20;// 100;//100;
	rule_parameters.TTD_hydrolysys_rate=rule_parameters.TTT_hydrolysys_rate*rule_parameters.hydrolysis_multiplier;
	//double DTD_hydrolysys_rate = TTD_hydrolysys_rate * 5;//2;//40;//10;// 100;//100;
	rule_parameters.DTD_hydrolysys_rate=rule_parameters.TTD_hydrolysys_rate*rule_parameters.hydrolysis_multiplier;
	//double rotation_rate = 0;
	parametersin >> rule_parameters.rotation_rate;

	rule_parameters.posibility_vertical = 1 / (1 + exp(-rule_parameters.bending_energy));
	//parametersin >> rule_parameters.posibility_vertical;
	rule_parameters.posibility_horizontal = exp(-rule_parameters.bending_energy) / (1 + exp(-rule_parameters.bending_energy));
	//parametersin >> rule_parameters.posibility_horizontal;
	parametersin.close();
};



void LoadAllParameters()
{
	ifstream parametersin;
	parametersin.open("parameters.txt");
	
	//int empty_anchor_unit_num=0;
	parametersin >> rule_parameters.empty_anchor_unit_num ;
	//for hydrolysis
	//int bundling_sum=0;
	parametersin >> rule_parameters.bundling_sum ;
	//int bundling_inverse_sum=0;
	parametersin >> rule_parameters.bundling_inverse_sum ;
	//int debundling_sum=0;
	parametersin >> rule_parameters.debundling_sum ;
	//int debundling_inverse_sum=0;
	parametersin >> rule_parameters.debundling_inverse_sum ;
	//int TTT_sum_total=0;
	parametersin >> rule_parameters.TTT_sum_total ;
	//int TTD_sum_total;//TTD and DTT are the same here
	parametersin >> rule_parameters.TTD_sum_total ;
	//int DTD_sum_total;
	parametersin >> rule_parameters.DTD_sum_total ;
	//int TT_annealing_sum = 0;
	parametersin >> rule_parameters.TT_annealing_sum ;
	//int TD_annealing_sum = 0;
	parametersin >> rule_parameters.TD_annealing_sum ;
	//int TT_polymerize_end_vertical_plus_sum=0;
	parametersin >> rule_parameters.TT_polymerize_end_vertical_plus_sum ;
	//int TT_polymerize_end_vertical_minus_sum=0;
	parametersin >> rule_parameters.TT_polymerize_end_vertical_minus_sum ;
	//int TT_polymerize_end_horizontal_plus_sum = 0;
	parametersin >> rule_parameters.TT_polymerize_end_horizontal_plus_sum ;
	//int TT_polymerize_end_horizontal_minus_sum = 0;
	parametersin >> rule_parameters.TT_polymerize_end_horizontal_minus_sum ;
	//int TD_polymerize_end_vertical_plus_sum = 0;
	parametersin >> rule_parameters.TD_polymerize_end_vertical_plus_sum ;
	//int TD_polymerize_end_vertical_minus_sum = 0;
	parametersin >> rule_parameters.TD_polymerize_end_vertical_minus_sum ;
	//int TD_polymerize_end_horizontal_plus_sum = 0;
	parametersin >> rule_parameters.TD_polymerize_end_horizontal_plus_sum ;
	//int TD_polymerize_end_horizontal_minus_sum = 0;
	parametersin >> rule_parameters.TD_polymerize_end_horizontal_minus_sum ;
	//int TT_depolymerize_plus_end_sum=0;
	parametersin >> rule_parameters.TT_depolymerize_plus_end_sum ;
	//int TD_depolymerize_plus_end_sum=0;//TD and DT are the same here
	parametersin >> rule_parameters.TD_depolymerize_plus_end_sum ;
	//int DD_depolymerize_plus_end_sum=0;
	parametersin >> rule_parameters.DD_depolymerize_plus_end_sum ;
	//int TT_depolymerize_minus_end_sum=0;
	parametersin >> rule_parameters.TT_depolymerize_minus_end_sum ;
	//int TD_depolymerize_minus_end_sum=0;//TD and DT are the same here
	parametersin >> rule_parameters.TD_depolymerize_minus_end_sum ;
	//int DD_depolymerize_minus_end_sum=0;
	parametersin >> rule_parameters.DD_depolymerize_minus_end_sum ;
	//int polymer_diffusion_propensity;
	parametersin >> rule_parameters.polymer_diffusion_propensity ;
	//int empty_anchor_diffusion_direction_sum;
	parametersin >> rule_parameters.empty_anchor_diffusion_direction_sum ;
	//int emtpy_anchor_attatch_sum;
	parametersin >> rule_parameters.emtpy_anchor_attatch_sum ;
	//int TT_fragmentation_sum;
	parametersin >> rule_parameters.TT_fragmentation_sum ;
	//int TD_fragmentation_sum;
	parametersin >> rule_parameters.TD_fragmentation_sum ;
	//int DD_fragmentation_sum;
	parametersin >> rule_parameters.DD_fragmentation_sum ;

	//int rotation_sum;
	parametersin >> rule_parameters.rotation_sum ;
	//double TD_case = 0;//T------D
	parametersin >> rule_parameters.TD_case ;
	//double DT_case = 0;//D------T
	parametersin >> rule_parameters.DT_case ;
	//double TT_case = 0;//T------T
	parametersin >> rule_parameters.TT_case ;
	//double DD_case = 0;//D------D
	parametersin >> rule_parameters.DD_case ;
	//int T_sum;
	parametersin >> rule_parameters.T_sum ;
	//int D_sum;
	parametersin >> rule_parameters.D_sum ;
	//int TT_annealing_sum=0;
	//int step_num;
	parametersin >> rule_parameters.step_num ;
	//int Z_ctp = 0;
	parametersin >> rule_parameters.Z_ctp ;
	//int reaction_mark;
	parametersin >> rule_parameters.reaction_mark ;
	//double cumps;
	parametersin >> rule_parameters.cumps ;
	//double sumps;
	parametersin >> rule_parameters.sumps ;
	//double t;
	parametersin >> rule_parameters.t ;
	//double polymer_sequence_average_length = 0;
	parametersin >> rule_parameters.polymer_sequence_average_length ;
	//int cluster_bound = 0;
	parametersin >> rule_parameters.cluster_bound ;
	//int up_larger_than_down_annealing = 0;
	parametersin >> rule_parameters.up_larger_than_down_annealing ;
	//int up_larger_than_down_fragmentation = 0;
	parametersin >> rule_parameters.up_larger_than_down_fragmentation ;
	//int diffusion_count = 0;
	parametersin >> rule_parameters.diffusion_count ;
	//int TT_annealing_count_up_larger_than_below = 0;
	parametersin >> rule_parameters.TT_annealing_count_up_larger_than_below ;
	//int TT_annealing_count_up_smaller_than_below = 0;
	parametersin >> rule_parameters.TT_annealing_count_up_smaller_than_below ;
	//int TD_annealing_count_up_larger_than_below = 0;
	parametersin >> rule_parameters.TD_annealing_count_up_larger_than_below ;
	//int TD_annealing_count_up_smaller_than_below = 0;
	parametersin >> rule_parameters.TD_annealing_count_up_smaller_than_below ;
	//int DT_annealing_count_up_larger_than_below = 0;
	parametersin >> rule_parameters.DT_annealing_count_up_larger_than_below ;
	//int DT_annealing_count_up_smaller_than_below = 0;
	parametersin >> rule_parameters.DT_annealing_count_up_smaller_than_below ;
	//int debundle_pair_TT;
	parametersin >> rule_parameters.debundle_pair_TT ;
	//int debundle_pair_TD;
	parametersin >> rule_parameters.debundle_pair_TD ;
	//int debundle_pair_DD;
	parametersin >> rule_parameters.debundle_pair_DD ;
	//int total_ftsZ_count = 0;
	parametersin >> rule_parameters.total_ftsZ_count ;
	//double total_ftsZ_lifetime = 0;
	parametersin >> rule_parameters.total_ftsZ_lifetime ;



	rule_parameters.ps.resize(num_of_rules);//17 is the current rules i have
	
	rule_parameters.count.resize(num_of_rules);//17 is the current rules i have
	
	rule_parameters.time.resize(num_of_rules + 1);//21 is the current rules i have
	
	rule_parameters.direction_record.resize(num_of_rules, vector<int>(4));

	rule_parameters.diffusion_statistics.resize(rule_parameters.row_num, vector<int>(4));
	
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			rule_parameters.diffusion_statistics[i][j] = 0;
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

	//vector<vector<int>> direction_record;
	for (int i = 0; i < num_of_rules; i++)
	{

		for (int j = 0; j < 4; j++)
		{
			double lalala;
			//parametersin >> lalala;
			//rule_parameters.direction_record[i][j] = lalala;
			//cout << lalala << endl;
			parametersin >> rule_parameters.direction_record[i][j];
		}
		
	}
	//vector<double> ps;
	for (int i = 0; i < rule_parameters.ps.size(); i++)
	{
		double lalala;
		parametersin >> rule_parameters.ps[i];
		
	}
	//vector<double> time;
	for (int i = 0; i < rule_parameters.time.size(); i++)
	{
		parametersin >> rule_parameters.time[i];
	}
	//vector<double> count;
	for (int i = 0; i < rule_parameters.count.size(); i++)
	{
		parametersin >> rule_parameters.count[i];
	}

	rule_parameters.diffusion_statistics.resize(rule_parameters.row_num, vector<int>(4));
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			rule_parameters.diffusion_statistics[i][j] = 0;
		}

	}

	rule_parameters.polymerize_statistics.resize(rule_parameters.row_num, vector<int>(4));
	for (int i = 0; i < rule_parameters.row_num; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			rule_parameters.polymerize_statistics[i][j] = 0;
		}
		
	}

	//int typical_column = 600;
	parametersin >> rule_parameters.typical_column ;
	//int typical_row = 800;
	parametersin >> rule_parameters.typical_row ;


	//int column_num = 600;//600//150
	parametersin >> rule_parameters.column_num ;
	//int row_num = 563;//800//200
	parametersin >> rule_parameters.row_num ;
	//int Anchor_num = 1000;//1000
	parametersin >> rule_parameters.Anchor_num ;
	//int Z = 8000;//8000//600//1000//4000
	parametersin >> rule_parameters.Z ;

	//double bending_energy = 1;
	parametersin >> rule_parameters.bending_energy ;
	//double scale_factor = typical_row / (double)row_num*typical_column*typical_column / column_num / column_num;
	parametersin >> rule_parameters.scale_factor ;

	//double anchoring_rate = 0.005*scale_factor;//0.0005;
	parametersin >> rule_parameters.anchoring_rate ;
	//double deanchoring_rate = 10;//8;//4;//0.2;
	parametersin >> rule_parameters.deanchoring_rate ;
	//double TT_polymerize_vertical_plus_rate = 0.01*scale_factor;//0.03;//0.05
	parametersin >> rule_parameters.TT_polymerize_horizontal_plus_rate ;
	//double TT_polymerize_vertical_minus_rate = 0 * scale_factor;//0.01;//0.05
	parametersin >> rule_parameters.TT_polymerize_vertical_minus_rate ;
	//double TT_polymerize_horizontal_plus_rate = TT_polymerize_vertical_plus_rate*exp(-bending_energy);
	parametersin >> rule_parameters.TT_polymerize_horizontal_plus_rate ;
	//double TT_polymerize_horizontal_minus_rate = TT_polymerize_vertical_minus_rate*exp(-bending_energy);
	parametersin >> rule_parameters.TT_polymerize_horizontal_minus_rate ;
	//double TD_polymerize_vertical_plus_rate = 0.005*scale_factor;//0.05
	parametersin >> rule_parameters.TD_polymerize_vertical_plus_rate ;
	//double TD_polymerize_vertical_minus_rate = 0 * scale_factor;//0.05
	parametersin >> rule_parameters.TD_polymerize_vertical_minus_rate ;
	//double TD_polymerize_horizontal_plus_rate = TD_polymerize_vertical_plus_rate*exp(-bending_energy);
	parametersin >> rule_parameters.TD_polymerize_horizontal_plus_rate ;
	//double TD_polymerize_horizontal_minus_rate = TD_polymerize_vertical_minus_rate*exp(-bending_energy);
	parametersin >> rule_parameters.TD_polymerize_horizontal_minus_rate ;
	//double diffusion_rate = 4000;//1000;//4000
	parametersin >> rule_parameters.diffusion_rate ;
	//double emtpy_anchor_diffusion_rate = 4000;//1000;//0.1;
	parametersin >> rule_parameters.emtpy_anchor_diffusion_rate ;
	//double Ebdl = 0.05;//0.03//0.07
	parametersin >> rule_parameters.Ebdl ;
	//double Ebdl_inverse = 0;
	parametersin >> rule_parameters.Ebdl_inverse ;
	//double bundling_rate = 2.5;//50;//100;//500
	parametersin >> rule_parameters.bundling_rate ;
	//double bundling_inverse_rate = 0;
	parametersin >> rule_parameters.bundling_inverse_rate ;
	//double debundling_rate = bundling_rate*exp(-Ebdl);
	parametersin >> rule_parameters.debundling_rate ;
	//double debundling_inverse_rate = bundling_inverse_rate*exp(-Ebdl_inverse);
	parametersin >> rule_parameters.debundling_inverse_rate ;
	//double TT_annealing_rate = 100;//2500;
	parametersin >> rule_parameters.TT_annealing_rate ;
	//double TD_annealing_rate = 10;
	parametersin >> rule_parameters.TD_annealing_rate ;
	//double TT_fragmentation_rate = 0.001;//1;//8; #
	parametersin >> rule_parameters.TT_fragmentation_rate ;
	//double TD_fragmentation_rate = 0.0001;//3;// 8;#
	parametersin >> rule_parameters.TD_fragmentation_rate ;
	//double DD_fragmentation_rate = 0.00001;//5;//8;#5
	parametersin >> rule_parameters.DD_fragmentation_rate ;
	//double TT_depolymerize_minus_rate = 0.1;//1;//8;//80#1
	parametersin >> rule_parameters.TT_depolymerize_minus_rate ;
	//double TD_depolymerize_minus_rate = TT_depolymerize_minus_rate * 40;//5;//20;//#5
	parametersin >> rule_parameters.TD_depolymerize_minus_rate ;
	//double DD_depolymerize_minus_rate = TD_depolymerize_minus_rate * 40;//10;//40;//#5
	parametersin >> rule_parameters.DD_depolymerize_minus_rate ;
	//double TT_depolymerize_plus_rate = TT_depolymerize_minus_rate / 10;//1;//8;//80#1
	parametersin >> rule_parameters.TT_depolymerize_plus_rate ;
	//double TD_depolymerize_plus_rate = TT_depolymerize_plus_rate * 40;//5;//20;//#5
	parametersin >> rule_parameters.TD_depolymerize_plus_rate ;
	//double DD_depolymerize_plus_rate = TD_depolymerize_plus_rate * 40;//10;//40;//#5
	parametersin >> rule_parameters.DD_depolymerize_plus_rate ;
	//double depolymerize_rate = 20;
	parametersin >> rule_parameters.depolymerize_rate ;

	//double attatch_rate = 400;//100;//100//0.000001;
	parametersin >> rule_parameters.attatch_rate ;
	//double both = 200;

	//double hydrolysis_rate = 0.0156;
	parametersin >> rule_parameters.hydrolysis_rate ;
	//double TTT_hydrolysys_rate = 0.0985;//0.5;//2;//10;//100;//100;
	parametersin >> rule_parameters.hydrolysis_multiplier;
	parametersin >> rule_parameters.TTT_hydrolysys_rate ;
	//double TTD_hydrolysys_rate = TTT_hydrolysys_rate * 5;//1;//20;//20;// 100;//100;
	parametersin >> rule_parameters.TTD_hydrolysys_rate ;
	//double DTD_hydrolysys_rate = TTD_hydrolysys_rate * 5;//2;//40;//10;// 100;//100;
	parametersin >> rule_parameters.DTD_hydrolysys_rate ;
	//double rotation_rate = 0;
	parametersin >> rule_parameters.rotation_rate ;

	//double posibility_vertical = 1 / (1 + exp(-bending_energy));
	parametersin >> rule_parameters.posibility_vertical ;
	//double posibility_horizontal = exp(-bending_energy) / (1 + exp(-bending_energy));
	parametersin >> rule_parameters.posibility_horizontal ;
	parametersin.close();
};
void LoadAllObject()
{
	
	ifstream pair_in;
	pair_in.open("unit_pair.txt");
	int size = 0;
	pair_in >> size;
	for (int i = 0; i < size; i++)
	{
		rule_structure.anealing_polymer_unit_pair_list.push_back(make_pair(rule_structure.polymer_unit_list_end.begin(), rule_structure.polymer_unit_list_end.begin()));
	}
	pair_in >> size;
	for (int i = 0; i < size; i++)
	{
		rule_structure.bundling_polymer_unit_pair_list.push_back(make_pair(rule_structure.polymer_unit_list_end.begin(), rule_structure.polymer_unit_list_end.begin()));
	}
	pair_in >> size;
	for (int i = 0; i < size; i++)
	{
		rule_structure.debundling_polymer_unit_pair_list.push_back(make_pair(rule_structure.polymer_unit_list_end.begin(), rule_structure.polymer_unit_list_end.begin()));
	}
	pair_in.close();

	/*Anchor_Unit anchor_unit_temp;
	rule_structure.anchor_list.push_back(anchor_unit_temp);
	auto anchor_unit_iter = rule_structure.anchor_list.begin();*/

	//grid_unit
	ifstream grid_unit_in;
	grid_unit_in.open("grid_unit.txt");
	int grid_unit_id = 0;
	int anchor_temp_id = 0;
	Grid.resize(rule_parameters.column_num, vector<Grid_Unit>(rule_parameters.row_num));
	for (int i = 0; i < rule_parameters.column_num; i++)
		for (int j = 0; j < rule_parameters.row_num; j++)
		{
			grid_unit_in>>grid_unit_id;
			grid_unit_in >> Grid[i][j].col_position;
			grid_unit_in >> Grid[i][j].row_position;
			grid_unit_in >> anchor_temp_id;//Grid[i][j].grid_anchor_pointer
			Grid[i][j].grid_polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
			Grid[i][j].grid_anchor_pointer = rule_structure.anchor_list.end();
			//Grid[i][j].grid_anchor_pointer = anchor_unit_iter;
			
		}
	grid_unit_in.close();




	//anchor_list
	ifstream anchor_unit_in;
	anchor_unit_in.open("anchor_unit.txt");
	//auto anchor_unit_iter = rule_structure.anchor_list.begin();
	int anchor_id,col, row = 0;
	for (int i = 0; i < rule_parameters.Anchor_num; i++)
	{
		anchor_unit_in >> anchor_id;
		
		anchor_unit_in >> col;
		anchor_unit_in >> row;
		Anchor_Unit anchor_unit;
		//anchor_unit.frap_mark = 0;
		anchor_unit_in>>anchor_unit.frap_mark ;
		anchor_unit.anchor_grid_pointer = &Grid[col][row];
		anchor_unit.polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
		rule_structure.anchor_list.push_back(anchor_unit);
		auto anchor_unit_iter = (--rule_structure.anchor_list.end());
		Grid[col][row].grid_anchor_pointer = anchor_unit_iter;
		for (int i = direction_up; i <= direction_right; i++)
		{
			anchor_unit_in>>anchor_unit_iter->availabe_diffusion_direction[i];
		}
		anchor_unit_in >> anchor_unit_iter->available_polymer_diffusion_direction_sum;
		anchor_unit_in >> anchor_unit_iter->consider_mark;
		anchor_unit_in >> anchor_unit_iter->attatch_mark;
	
	}

	anchor_unit_in.close();



	//polymer_bundle_list
	ifstream polymer_bundle_in;
	polymer_bundle_in.open("polymer_bundle.txt");
	int bundle_id = 0;
	int bundle_list_size = 0;
	int position=0;
	int bundle_cluster_pointer_size = 0;
	int bundle_anchor_pointer_size = 0;
	polymer_bundle_in >> bundle_list_size;
	for (int i = 0; i < bundle_list_size; i++)
	{
		Polymer_Bundle polymer_bundle;
		rule_structure.polymer_bundle_list.push_back(polymer_bundle);
		auto polymer_bundle_iter = (--rule_structure.polymer_bundle_list.end());
		polymer_bundle_in >> bundle_id;
		
		polymer_bundle_in >> bundle_cluster_pointer_size;
		for (int j = 0; j < bundle_cluster_pointer_size; j++)
		{
			
			polymer_bundle_in >> position;
			//bundle_cluster_pointer is still missing//solved
		}
		polymer_bundle_in >> bundle_anchor_pointer_size;
		for (int j = 0; j < bundle_anchor_pointer_size; j++)
		{

			polymer_bundle_in >> position;
			auto anchor_unit_iter = rule_structure.anchor_list.begin();
			advance(anchor_unit_iter, position);
			polymer_bundle_iter->polymer_anchor_list.push_back(anchor_unit_iter);
			
		}
		for (int j = direction_up; j <= direction_right; j++)
		{
			polymer_bundle_in>>polymer_bundle_iter->availabe_diffusion_direction[j];
		}
		polymer_bundle_in >> polymer_bundle_iter->available_polymer_diffusion_direction_sum;
		polymer_bundle_in >> polymer_bundle_iter->Consider_mark;
		polymer_bundle_in >> polymer_bundle_iter->bundle_direction;
	}
	polymer_bundle_in.close();
	

	//polymer_cluster_list
	ifstream polymer_unit_in;
	ifstream polymer_cluster_in;
	polymer_unit_in.open("polymer_unit.txt");
	polymer_cluster_in.open("polymer_cluster.txt");
	int cluster_list_size = 0;
	int sequence_list_size = 0;
	int anchor_list_size = 0;
	polymer_cluster_in >> cluster_list_size;
	int cluster_id = 0;
	position=0;
	col, row = 0;
	
	for (int i = 0; i < cluster_list_size; i++)
	{
		polymer_cluster_in >> cluster_id;
		Polymer_Cluster polymer_cluster;
		rule_structure.polymer_cluster_list.push_back(polymer_cluster);
		auto polymer_cluster_iter = (--rule_structure.polymer_cluster_list.end());
		polymer_cluster_in >> anchor_list_size;
		for (int j = 0; j < anchor_list_size; j++)
		{
			auto anchor_unit_iter = rule_structure.anchor_list.begin();
			polymer_cluster_in >> position;
			advance(anchor_unit_iter, position);
			polymer_cluster_iter->anchor_pointer_list.push_back(anchor_unit_iter);
		}
		polymer_cluster_in >> polymer_cluster_iter->direction;
		polymer_cluster_in >> polymer_cluster_iter->polarize;
		polymer_cluster_in >> sequence_list_size;
		for (int j = 0; j < sequence_list_size; j++)
		{
			Polymer_Unit polymer_unit;
			polymer_cluster_iter->polymer_sequence.push_back(polymer_unit);
			auto polymer_unit_iter = (--polymer_cluster_iter->polymer_sequence.end());
			polymer_unit_in >> position;
			if (position != -1)
			{
				auto anchor_unit_iter = rule_structure.anchor_list.begin();
				advance(anchor_unit_iter, position);
				polymer_unit_iter->anchor_pointer = anchor_unit_iter;
				anchor_unit_iter->polymer_unit_pointer = polymer_unit_iter;
			}
			
			polymer_unit_in >> polymer_unit_iter->direction;
			//polymer_unit_grid_iter
			polymer_unit_in >> col;
			polymer_unit_in >> row;
			polymer_unit_iter->polymer_grid_pointer = &Grid[col][row];
			Grid[col][row].grid_polymer_unit_pointer = polymer_unit_iter;
			polymer_unit_in >> polymer_unit_iter->Hydrolysis_mark;
			polymer_unit_in >> polymer_unit_iter->Consider_mark;
			polymer_unit_in >> polymer_unit_iter->lifetime;
			polymer_unit_in >> polymer_unit_iter->action_mark;
			polymer_unit_in >> polymer_unit_iter->attatch_mark;
			polymer_unit_in >> polymer_unit_iter->frap_mark;
			//polymer_unit_iter->frap_mark = 0;
			polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;
			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.anealing_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->anealing_polymer_unit_pair_list_iter[0] = pair_iter;
				pair_iter->second = polymer_unit_iter;
			}
			
			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.anealing_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->anealing_polymer_unit_pair_list_iter[1] = pair_iter;
				pair_iter->first = polymer_unit_iter;
			}

			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->bundling_polymer_unit_pair_list_iter[0] = pair_iter;
				pair_iter->second = polymer_unit_iter;
			}
			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.bundling_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->bundling_polymer_unit_pair_list_iter[1] = pair_iter;
				pair_iter->first = polymer_unit_iter;
			}

			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->debundling_polymer_unit_pair_list_iter[0] = pair_iter;
				pair_iter->second = polymer_unit_iter;
			}
			polymer_unit_in >> position;
			if (position != -1)
			{
				auto pair_iter = rule_structure.debundling_polymer_unit_pair_list.begin();
				advance(pair_iter, position);
				polymer_unit_iter->debundling_polymer_unit_pair_list_iter[1] = pair_iter;
				pair_iter->first = polymer_unit_iter;
			}

		}

		auto polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		polymer_cluster_in >> position;
		advance(polymer_bundle_iter, position);
		polymer_cluster_iter->cluster_bundle_pointer = polymer_bundle_iter;
		polymer_bundle_iter->bundle_cluster_pointer_list.push_back(polymer_cluster_iter);

		polymer_cluster_in>> polymer_cluster_iter->TTT_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TTD_sum ;
		polymer_cluster_in >> polymer_cluster_iter->DTD_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TT_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TD_sum ;
		polymer_cluster_in >> polymer_cluster_iter->DD_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TT_depolymerize_plus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TD_depolymerize_plus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->DD_depolymerize_plus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TT_depolymerize_minus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TD_depolymerize_minus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->DD_depolymerize_minus_end_sum ;

		polymer_cluster_in >> polymer_cluster_iter->Consider_mark ;
		polymer_cluster_in >> polymer_cluster_iter->ring_mark ;
		polymer_cluster_in >> polymer_cluster_iter->TT_polymerize_plus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TT_polymerize_minus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TD_polymerize_plus_end_sum ;
		polymer_cluster_in >> polymer_cluster_iter->TD_polymerize_minus_end_sum ;
		for (int j = direction_up; j <= direction_right; j++)
		{
			polymer_cluster_in >> polymer_cluster_iter->availabe_polymerize_direction[j];
		}
	}

	polymer_unit_in.close();
	polymer_cluster_in.close();



	
	



};