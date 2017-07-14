#pragma once
//#define DEBUG_RANDOM 0;
class Rule_Parameters
{
public:
	int empty_anchor_unit_num=0;
	//for hydrolysis
	int bundling_sum=0;
	int bundling_inverse_sum=0;
	int debundling_sum=0;
	int debundling_inverse_sum=0;
	int TTT_sum_total=0;
	int TTD_sum_total;//TTD and DTT are the same here
	int DTD_sum_total;
	int TT_annealing_sum = 0;
	int TD_annealing_sum = 0;
	int TT_polymerize_end_vertical_plus_sum=0;
	int TT_polymerize_end_vertical_minus_sum=0;
	int TT_polymerize_end_horizontal_plus_sum = 0;
	int TT_polymerize_end_horizontal_minus_sum = 0;
	int TD_polymerize_end_vertical_plus_sum = 0;
	int TD_polymerize_end_vertical_minus_sum = 0;
	int TD_polymerize_end_horizontal_plus_sum = 0;
	int TD_polymerize_end_horizontal_minus_sum = 0;
	int TT_depolymerize_plus_end_sum=0;
	int TD_depolymerize_plus_end_sum=0;//TD and DT are the same here
	int DD_depolymerize_plus_end_sum=0;
	int TT_depolymerize_minus_end_sum=0;
	int TD_depolymerize_minus_end_sum=0;//TD and DT are the same here
	int DD_depolymerize_minus_end_sum=0;
	int polymer_diffusion_direction_sum;
	int empty_anchor_diffusion_direction_sum;
	int emtpy_anchor_attatch_sum;
	int TT_fragmentation_sum;
	int TD_fragmentation_sum;
	int DD_fragmentation_sum;
	/*int bundling_sum=0;
	int debundling_sum=0;*/
	int rotation_sum;
	double TD_case = 0;//T------D
	double DT_case = 0;//D------T
	double TT_case = 0;//T------T
	double DD_case = 0;//D------D
	int T_sum;
	int D_sum;
	//int TT_annealing_sum=0;
	int step_num;
	int Z_ctp=0;
	int reaction_mark;
	double cumps;
	double sumps;
	double t;
	double polymer_sequence_average_length=0;
	int cluster_bound=0;
	int up_larger_than_down_annealing=0;
	int up_larger_than_down_fragmentation = 0;
	int diffusion_count = 0;
	int TT_annealing_count_up_larger_than_below = 0;
	int TT_annealing_count_up_smaller_than_below = 0;
	int TD_annealing_count_up_larger_than_below = 0;
	int TD_annealing_count_up_smaller_than_below = 0;
	int DT_annealing_count_up_larger_than_below = 0;
	int DT_annealing_count_up_smaller_than_below = 0;
	int debundle_pair_TT;
	int debundle_pair_TD;
	int debundle_pair_DD;
	int total_ftsZ_count = 0;
	double total_ftsZ_lifetime = 0;
	vector<vector<int>> direction_record;
	vector<double> ps;
	vector<double> time;
	vector<double> count;



	int typical_column = 600;
	int typical_row = 800;


	int column_num = 600;//600//150
	int row_num = 563;//800//200
	int Anchor_num = 1000;//1000
	int Z = 8000;//8000//600//1000//4000

	double bending_energy = 1;
	double scale_factor = typical_row / (double)row_num*typical_column*typical_column / column_num / column_num;


	double anchoring_rate = 0.005*scale_factor;//0.0005;
	double deanchoring_rate = 10;//8;//4;//0.2;
	double TT_polymerize_vertical_plus_rate = 0.01*scale_factor;//0.03;//0.05
	double TT_polymerize_vertical_minus_rate = 0*scale_factor;//0.01;//0.05
	double TT_polymerize_horizontal_plus_rate = TT_polymerize_vertical_plus_rate*exp(-bending_energy);
	double TT_polymerize_horizontal_minus_rate = TT_polymerize_vertical_minus_rate*exp(-bending_energy);
	double TD_polymerize_vertical_plus_rate = 0.005*scale_factor;//0.05
	double TD_polymerize_vertical_minus_rate = 0*scale_factor;//0.05
	double TD_polymerize_horizontal_plus_rate = TD_polymerize_vertical_plus_rate*exp(-bending_energy);
	double TD_polymerize_horizontal_minus_rate = TD_polymerize_vertical_minus_rate*exp(-bending_energy);
	double diffusion_rate = 4000;//1000;//4000
	double emtpy_anchor_diffusion_rate = 4000;//1000;//0.1;
	double Ebdl = 0.05;//0.03//0.07
	double Ebdl_inverse = 0;
	double bundling_rate = 2.5;//50;//100;//500
	double bundling_inverse_rate = 0;
	double debundling_rate = bundling_rate*exp(-Ebdl);
	double debundling_inverse_rate = bundling_inverse_rate*exp(-Ebdl_inverse);
	double TT_annealing_rate = 100;//2500;
	double TD_annealing_rate = 10;
	double TT_fragmentation_rate = 0.001;//1;//8; #
	double TD_fragmentation_rate = 0.01;//3;// 8;#
	double DD_fragmentation_rate = 0.1;//5;//8;#5
	double TT_depolymerize_minus_rate = 2;//1;//8;//80#1
	double TD_depolymerize_minus_rate = TT_depolymerize_minus_rate * 40;//5;//20;//#5
	double DD_depolymerize_minus_rate = TD_depolymerize_minus_rate * 40;//10;//40;//#5
	double TT_depolymerize_plus_rate = TT_depolymerize_minus_rate / 10;//1;//8;//80#1
	double TD_depolymerize_plus_rate = TT_depolymerize_plus_rate * 40;//5;//20;//#5
	double DD_depolymerize_plus_rate = TD_depolymerize_plus_rate * 40;//10;//40;//#5
	double depolymerize_rate = 20;
	
	double attatch_rate = 400;//100;//100//0.000001;
	double both = 200;

	double hydrolysis_rate = 0.1150;
	double hydrolysis_multiplier = 2;
	double TTT_hydrolysys_rate = 0.0985;//0.5;//2;//10;//100;//100;
	double TTD_hydrolysys_rate = TTT_hydrolysys_rate * hydrolysis_multiplier;//1;//20;//20;// 100;//100;
	double DTD_hydrolysys_rate = TTD_hydrolysys_rate * hydrolysis_multiplier;//2;//40;//10;// 100;//100;
	
	double rotation_rate = 0;

	double posibility_vertical = 1 / (1 + exp(-bending_energy));
	double posibility_horizontal = exp(-bending_energy) / (1 + exp(-bending_energy));

};
