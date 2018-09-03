#include <vector>
#include <iterator>
#include <iostream>
#include "parameters.h"
#include "initialize.h"
#include "Grid_Unit.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include "test.h"
//#include <sys/types.h>
//#include <sys/stat.h>
//include <stdlib.h>
//#include <stdio.h>
#include "General_Functions.h"
#include <random>
#include <omp.h>
//#include <direct.h>//adding this for _mkdir

Rule_Structure rule_structure;
Rule_Parameters rule_parameters;
vector<vector<Grid_Unit>> Grid(rule_parameters.column_num, vector<Grid_Unit>(rule_parameters.row_num));//the grid we are using
//vector<vector<Grid_Unit>> Polymer_Grid(column_num, vector<Grid_Unit>(row_num));//the grid we are using
//vector<vector<Grid_Unit>> Polymer_Cluster_Grid(column_num, vector<Grid_Unit>(row_num));//the grid we are using




using namespace std;
int main()
{
	//vector<vector<int> > AnchorGrid(column_num, vector<int>(row_num));
	Initialize();
	cout << rule_parameters.sumps << endl;
	cout << rule_parameters.scale_factor << endl;
	cout << "Initialize complete" << endl;
	//system("pause");
	//cout << "empty anchor: " << rule_parameters.empty_anchor_unit_num << endl;
	//Grid grid;
	double fidx = 2400;//total time
	double tjudge = 0;
	double tbase = 300;
	double tdelta=0.5;//how much time it output another data
	double tdelta2 = 0.1;
	double ttemp = 0.000000001;
	rule_parameters.t = 0;
	char output_dir[120];
	//char output_file[120];
	char fname[120];
	cin;
	ofstream fout;
	//snprintf(output_dir, sizeof(output_dir), ".//%s%3.3lf/", "time", rule_parameters.t);
	/*_mkdir(output_dir);
	snprintf(fname, sizeof(fname), "%s%s", output_dir, "testmatrix.txt");*/
	fout.open("output.txt");

	ofstream statistics;
	statistics.open("statistics.txt");
	//output(fout);
	/*if (tjudge > tbase)
		tjudge = tbase;*/
	//tjudge = tbase + tjudge;
	while (rule_parameters.t <= fidx)
	{
		
		//cout << output_dir << endl;
	
		if (rule_parameters.t >= tjudge)
		{
			//snprintf(output_dir, sizeof(output_dir), ".//%s%3.3lf/", "time", rule_parameters.t);
			//_mkdir(output_dir);
			
			//snprintf(fname, sizeof(fname), "%s%s", output_dir, "testmatrix.txt");
			//cout << fname << endl;
			//system("pause");
			
			
			//cout << "output" << endl;


			output(fout);
		

			ofstream timeout;
			timeout.open("time.txt");
			for (int i = 0; i < num_of_rules+1; i++)
			{
				timeout << i<<" "<<rule_parameters.time[i] << endl;
			}
			timeout.close();

			ofstream countout;
			countout.open("count.txt");
			for (int i = 0; i < num_of_rules; i++)
			{
				countout << i<<" "<<rule_parameters.count[i] << endl;
			}
			countout.close();


			ofstream direction;
			direction.open("direction.txt");
			for (int i = 0; i < num_of_rules; i++)
			{
				direction << i << ": ";
				for (int j = 0; j < 4; j++)
				{
					direction << rule_parameters.direction_record[i][j] << " ";
				}
				direction << endl;
			}
			//int TD_case = 0;//T------D
			//int DT_case = 0;//D------T
			//int TT_case = 0;//T------T
			//int DD_case = 0;//D------D
			//int double_sum;
			
			direction << "TT: " << rule_parameters.TT_case << endl;
			direction << "TD: " << rule_parameters.TD_case << endl;
			direction << "DT: " << rule_parameters.DT_case << endl;
			direction << "DD: " << rule_parameters.DD_case << endl;
			direction << "TTT sum: " << rule_parameters.TTT_sum_total << endl;
			direction << "TTD sum: " << rule_parameters.TTD_sum_total << endl;
			direction << "DTD sum: " << rule_parameters.DTD_sum_total << endl;
			direction << "fragmentation TT: " << rule_parameters.TT_fragmentation_sum<<endl;
			direction << "fragmentation TD: " << rule_parameters.TD_fragmentation_sum << endl;;
			direction << "fragmentation DD: " << rule_parameters.DD_fragmentation_sum << endl;;
			direction << "depolymerization TT: " << rule_parameters.TT_depolymerize_minus_end_sum << endl;
			direction << "depolymerization TD: " << rule_parameters.TD_depolymerize_minus_end_sum << endl;
			direction << "depolymerization DD: " << rule_parameters.DD_depolymerize_minus_end_sum << endl;
			direction << "PS:" << endl;
			for (int i = 0; i < rule_parameters.ps.size(); i++)
			{
				direction << i << ": " << rule_parameters.ps[i] << endl;
			}
			direction << "up is larger than down for annealing: " << rule_parameters.up_larger_than_down_annealing << endl;
			direction << "up is larger than down for fragmentation: " << rule_parameters.up_larger_than_down_fragmentation << endl;
			direction << "TT up is larger than down for annealing: " << rule_parameters.TT_annealing_count_up_larger_than_below << endl;
			direction << "TT up is smaller than down for annealing: " << rule_parameters.TT_annealing_count_up_smaller_than_below << endl;
			direction << "TD up is larger than down for annealing: " << rule_parameters.TD_annealing_count_up_larger_than_below << endl;
			direction << "TD up is smaller than down for annealing: " << rule_parameters.TD_annealing_count_up_smaller_than_below << endl;
			direction << "DT up is larger than down for annealing: " << rule_parameters.DT_annealing_count_up_larger_than_below << endl;			
			direction << "DT up is smaller than down for annealing: " << rule_parameters.DT_annealing_count_up_smaller_than_below << endl;
			direction.close();

			check_debundling_information();
			rule_parameters.polymer_sequence_average_length = (double)(rule_parameters.Z-(double)rule_parameters.Z_ctp) / rule_structure.polymer_cluster_list.size();
			statistics << rule_parameters.t << " ";
			statistics << rule_parameters.Z_ctp << " ";
			statistics << rule_parameters.polymer_sequence_average_length << " ";
			statistics << rule_parameters.cluster_bound << " ";
			statistics << rule_parameters.TT_annealing_sum << " ";
			statistics << rule_structure.debundling_polymer_unit_pair_list.size() << " ";
			statistics << rule_parameters.debundle_pair_TT << " ";
			statistics << rule_parameters.debundle_pair_TD << " ";
			statistics << rule_parameters.debundle_pair_DD << " ";
			statistics << rule_parameters.T_sum << " ";
			statistics << rule_parameters.D_sum << " ";
			//statistics << rule_parameters.TT_fragmentation_sum+ rule_parameters.TT_depolymerize_end_sum << " ";
			//statistics << rule_parameters.TD_fragmentation_sum + rule_parameters.TD_depolymerize_end_sum<< " ";
			//statistics << rule_parameters.DD_fragmentation_sum + rule_parameters.DD_depolymerize_end_sum << " ";
			statistics << rule_parameters.total_ftsZ_lifetime / rule_parameters.total_ftsZ_count << " ";
			statistics << rule_parameters.bundling_sum <<" ";
			statistics << rule_parameters.bundling_inverse_sum << " ";
			statistics << rule_parameters.debundling_sum << " ";
			statistics << rule_parameters.debundling_inverse_sum << " ";
			statistics << rule_parameters.TTT_sum_total <<  " ";
			statistics << rule_parameters.TTD_sum_total << " ";
			statistics << rule_parameters.DTD_sum_total << " ";
			statistics << rule_parameters.TTT_hydrolysys_rate<< " ";
			statistics << rule_parameters.TT_depolymerize_minus_rate<< " ";
			statistics << rule_parameters.Z - rule_parameters.Z_ctp - rule_structure.polymer_cluster_list.size() << endl;
			

			//update in real time 
			/*double effective_hydrolysis_rate = ;
			double effective_depolymerize_rate = ;*/

			
			


			//cout << "time hit: " << rule_parameters.t;
			//output();
			//output(fout);
			//cout << "Time: " << rule_parameters.t << endl;
		/*	cout << "available anchoring anchor: " << rule_parameters.empty_anchor_unit_num << endl;
			cout << "available polymerize end: " << rule_parameters.polymerize_end_sum << endl;
			cout << "available TT depolymerize end: " << rule_parameters.TT_depolymerize_end_sum << endl;
			cout << "available TT fragmentation end: " << rule_parameters.TT_fragmentation_sum << endl;
			cout << "available diffusion direction sum: " << rule_parameters.polymer_diffusion_direction_sum << endl;
			cout << "available empty anchor diffusion direction sum: " << rule_parameters.empty_anchor_diffusion_direction_sum << endl;
			cout << "available attatch anchor sum: " << rule_parameters.emtpy_anchor_attatch_sum << endl;*/
			//system("pause");
			tjudge = rule_parameters.t + tdelta;
			//fout.close();
			
			//cout << rule_parameters.empty_anchor_unit_num << endl;
			//cout << tjudge << endl;
		}
		if (rule_parameters.sumps != 0)
		{
			
			
			//cout << "output:" << endl;
			double start = omp_get_wtime();
			//fout << "Time: " << rule_parameters.t << endl;
			if (rule_parameters.t > ttemp)
			{
				cout << "Time: " << rule_parameters.t << ";";
				cout << "effective hydrolysis rate:" << (rule_parameters.TTT_sum_total*rule_parameters.TTT_hydrolysys_rate + rule_parameters.TTD_sum_total*rule_parameters.TTD_hydrolysys_rate + rule_parameters.DTD_sum_total*rule_parameters.DTD_hydrolysys_rate) /
					(rule_parameters.TTT_sum_total + rule_parameters.TTD_sum_total + rule_parameters.DTD_sum_total) <<" hydrolysis_TTT rate:"<<rule_parameters.TTT_hydrolysys_rate<< ";";
				cout <<" effective depolymerization rate:"<< (rule_parameters.TT_depolymerize_minus_end_sum*rule_parameters.TT_depolymerize_minus_rate + rule_parameters.TD_depolymerize_minus_end_sum*rule_parameters.TD_depolymerize_minus_rate + rule_parameters.DD_depolymerize_minus_end_sum*rule_parameters.DD_depolymerize_minus_rate) /
					(rule_parameters.TT_depolymerize_minus_end_sum + rule_parameters.TD_depolymerize_minus_end_sum + rule_parameters.DD_depolymerize_minus_end_sum) <<" depolymerization TT minus rate:"<<rule_parameters.TT_depolymerize_minus_rate<< endl;
				ttemp += tdelta2;
			}
			
			select_and_execute();
			double end = omp_get_wtime();
			rule_parameters.time[num_of_rules] += (end - start);
			/*fout << "available anchoring anchor: " << rule_parameters.empty_anchor_unit_num << endl;
			fout << "available polymerize end: " << rule_parameters.polymerize_end_sum << endl;
			fout << "available TT depolymerize end: " << rule_parameters.TT_depolymerize_end_sum << endl;
			fout << "available diffusion direction sum: " << rule_parameters.polymer_diffusion_direction_sum << endl;*/
			//if (rule_parameters.t > tbase)
			//{
			//	break;
			//	cout << "available anchoring anchor: " << rule_parameters.empty_anchor_unit_num << endl;
			//	//cout << "available polymerize end: " << rule_parameters.polymerize_end_sum << endl;
			//	cout << "available TT depolymerize end: " << rule_parameters.TT_depolymerize_end_sum << endl;
			//	cout << "available TT fragmentation end: " << rule_parameters.TT_fragmentation_sum << endl;
			//	cout << "available diffusion direction sum: " << rule_parameters.polymer_diffusion_direction_sum << endl;
			//	cout << "available anchor diffusion direction sum: " << rule_parameters.empty_anchor_diffusion_direction_sum << endl;
			//	cout << "available attatch anchor sum: " << rule_parameters.emtpy_anchor_attatch_sum << endl;

			//	//output(fout);
			//	//output();
			//}
			//check_direction();
			//check();
			//check_polymerize_information();
			//test_annealing();
			//test_TT_Fragmentation_sum();
			//check_bundling();
			
			default_random_engine gen(rule_parameters.step_num);
			uniform_real_distribution<> tempf(0, 1);//here is for the time part
			double temp = 0;
			do
			{
				temp = tempf(gen);
			} while (temp == 0);

			double deltat = (1.0 / rule_parameters.sumps)*log(1.0 / temp);
			rule_parameters.t = rule_parameters.t + deltat;

			




			
			//system("pause");
			
		}
		else
		{
			////snprintf(output_dir, sizeof(output_dir), ".//%s%3.3lf/", "time", rule_parameters.t);
			//_mkdir(output_dir);
			//snprintf(fname, sizeof(fname), "%s%s", output_dir, "testmatrix.txt");
			//fout.open(fname);
			//output(fout);
			//fout.close();
			break;
		}
		
		
		

	}
	/*for (int i = 0; i < column_num; i++)
	{
		for (int j = 0; j < row_num; j++)
		{
			cout << rule_structure.Mark[i][j] <<" ";
		}
		cout << endl;
	}*/
	//system("pause");
	fout.close();
	statistics.close();
	


	//anchor_grid.output();
	//Initialize_Anchor();
	cin;
}