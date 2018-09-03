#include "Rule_Anchoring_Deanchoring.h"
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
//#include "Polymer_Unit.h"
#include "Grid_Unit.h"
#include <vector>
#include <cmath>
#include <random>
#include "DeconstructionFunction.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
//extern vector<vector<Grid_Unit>> Polymer_Grid;
//extern vector<vector<Grid_Unit>> Polymer_Cluster_Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
void anchoring()
{
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
	uniform_int_distribution<> Anchor_Rand(0, rule_parameters.empty_anchor_unit_num-rule_parameters.emtpy_anchor_attatch_sum-1);
	int anchor_rand = Anchor_Rand(gen);
	list<Anchor_Unit>::iterator anchoring_iter;
	anchoring_iter = rule_structure.anchor_list.begin();
	int anchoring_temp_sum=0;
	//cout << " enter anchoring function" << endl;
	//cout <<"anchor_rand: "<< anchor_rand << endl;
	int count = 0;
	//for (int i = 0; i < rule_structure.anchor_list.size(); i++)
	while(anchoring_iter != rule_structure.anchor_list.end())
	{
		//cout << std::addressof(anchoring_iter->polymer_unit_pointer) << " " << std::addressof(rule_structure.polymer_unit_list_end.begin()) << endl;
		//cout << *anchoring_iter->polymer_unit_pointer << endl;
		//system("Pause");
		
		if (&*anchoring_iter->polymer_unit_pointer==&*rule_structure.polymer_unit_list_end.begin()&&
			&*Grid[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())
			//when it is empty
			//when there is no monomer on this site. eg: no attatch
		{
			
			if (anchoring_temp_sum < anchor_rand)
			{
				anchoring_temp_sum++;
				anchoring_iter++;
				//cout << anchoring_temp_sum << endl;
			}
			else
			{//add a monomer to the anchor_unit
				//cout << "start to do anchoring" << endl;

				/*cout << "1: " << &*anchoring_iter->polymer_unit_pointer;
				cout << " 2: " << &*rule_structure.polymer_unit_list_end.begin() << endl;*/
				//anchoring_iter->polymer_unit_pointer != rule_structure.polymer_unit_list_end.begin();
				//cout << "before: " << rule_structure.Mark[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position] << endl;
				/*for (int i = 0; i < column_num; i++)
				{

				for (int j = 0; j < row_num; j++)
				{
				cout << rule_structure.Mark[i][j] << " ";
				}
				cout << endl;
				}*/



				Polymer_Unit polymer_unit;
				Polymer_Cluster polymer_cluster;
				Polymer_Bundle polymer_bundle;

				list<Polymer_Cluster>::iterator polymer_cluster_iter;
				list<Polymer_Bundle>::iterator polymer_bundle_iter;


				rule_structure.polymer_cluster_list.push_front(polymer_cluster);// polymer_unit polymer_cluster_iter
				polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
				/*polymer_cluster_iter->availabe_polymerize_direction[direction_up] = false;
				this->availabe_polymerize_direction[direction_down] = false;
				this->availabe_polymerize_direction[direction_left] = false;
				this->availabe_polymerize_direction[direction_right] = false;*/
				polymer_cluster_iter->anchor_pointer_list.push_front(anchoring_iter);

				//polymer_bundle.polymer_grid_pointer = anchoring_iter->anchor_grid_pointer;
				rule_structure.polymer_bundle_list.push_front(polymer_bundle);
				polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
				polymer_bundle_iter->bundle_cluster_pointer_list.push_front(polymer_cluster_iter);
				polymer_bundle_iter->polymer_anchor_list.push_front(anchoring_iter);




				polymer_cluster_iter->cluster_bundle_pointer = polymer_bundle_iter;//why

																					//anchoring_iter->polymer_cluster_pointer = polymer_cluster_iter;
																					//we still need to add the polymer_unit;
				polymer_cluster_iter->polymer_sequence.push_front(polymer_unit);
				list<Polymer_Unit>::iterator polymer_unit_iter;
				polymer_unit_iter = polymer_cluster_iter->polymer_sequence.begin();

				polymer_unit_iter->polymer_grid_pointer = anchoring_iter->anchor_grid_pointer;//&Grid[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position];//assign the polymer_grid
				polymer_unit_iter->anchor_pointer = anchoring_iter;
				polymer_unit_iter->polymer_cluster_pointer = polymer_cluster_iter;


				
				// i think i need to assign the directions here by some probability, don't know yet.


				anchoring_iter->polymer_unit_pointer = polymer_unit_iter;
				//anchoring_iter->IsAnchored = yes;
				//rule_structure.Mark[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position] = 1;
				Grid[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer = polymer_unit_iter;
				//int col = anchoring_iter->anchor_grid_pointer->col_position;
				//int row=anchoring_iter->anchor_grid_pointer->row_position;
				{//here i check the polymerize and diffusion boundary, but there is a problem, switch to use the general function

					
					uniform_real_distribution<> Direction_Rand(0, 1);
					uniform_int_distribution<> Polarize_Rand(0, 1);
					int polarize_rand = Polarize_Rand(gen);
					double direction_rand = Direction_Rand(gen);
					if (direction_rand <= rule_parameters.posibility_vertical)
					{
						polymer_cluster_iter->direction = vertical;
						polymer_cluster_iter->polarize = polarize_rand;
					}
					else
					{
						polymer_cluster_iter->direction = horizontal;
						polymer_cluster_iter->polarize = polarize_rand;
					}



					check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);
					check_polymer_hydrolysis_sum(polymer_unit_iter);
					//rule_parameters.empty_anchor_diffusion_direction_sum -= anchoring_iter->available_polymer_diffusion_direction_sum;
					rule_parameters.empty_anchor_diffusion_direction_sum -= anchoring_iter->available_polymer_diffusion_direction_sum;
					for (int i = 0; i < 4; i++)
					{
						anchoring_iter->availabe_diffusion_direction[i] = false;
					}
					anchoring_iter->available_polymer_diffusion_direction_sum = 0;
					//int new_col, new_row = 0;
					//new_col = col + 1;
					//new_row = row;
					//if (new_col > column_num - 1)
					//	new_col = 0;
					//if (&*Grid[new_col][new_row].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())//nothing there
					//{
					//	polymer_unit_iter->polymer_cluster_pointer->availabe_polymerize_direction[direction_up] = true;
					//	polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_up] = true;
					//	polymer_unit_iter->polymer_cluster_pointer->polymerize_end_sum++;
					//	rule_parameters.polymerize_end_sum = rule_parameters.polymerize_end_sum++;
					//	rule_parameters.polymer_diffusion_direction_sum++;
					//}
					//new_col = col - 1;
					//new_row = row;
					//if (new_col <0)
					//	new_col = column_num - 1;
					//if (&*Grid[new_col][new_row].grid_polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin())//nothing there
					//{
					//	polymer_unit_iter->polymer_cluster_pointer->availabe_polymerize_direction[direction_down] = true;
					//	polymer_unit_iter->polymer_cluster_pointer->cluster_bundle_pointer->availabe_diffusion_direction[direction_down] = true;
					//	polymer_unit_iter->polymer_cluster_pointer->polymerize_end_sum++;
					//	rule_parameters.polymerize_end_sum = rule_parameters.polymerize_end_sum++;
					//	rule_parameters.polymer_diffusion_direction_sum++;
					//}

					//int temp_sum = 0;
					//for (int i = 0; i < 4; i++)
					//{
					//temp_sum = temp_sum + polymer_bundle_iter->availabe_diffusion_direction[i];
					//}
					//polymer_bundle_iter->available_polymer_diffusion_direction_sum=temp_sum;


				}//check boundary considions

					//rule_parameters.polymerize_end_sum_sum = rule_parameters.polymerize_end_sum_sum + 2;
					//this->anchor_pointer = rule_structure.anchor_list_end.begin();
					//this->polymer_cluster_pointer = rule_structure.polymer_cluster_list_end.begin();}


					//anchoring_iter->anchor_grid_pointer->grid_polymer_unit_pointer = polymer_iter;
					//Grid[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer->polymer_cluster_pointer = polymer_cluster_iter;
					//Grid[anchoring_iter->anchor_grid_pointer->col_position][anchoring_iter->anchor_grid_pointer->row_position].grid_anchor_pointer = anchoring_iter;
					//cout << "jump" << endl;
					/*cout << "Col: " << anchoring_iter->anchor_grid_pointer->col_position << " Row: " << anchoring_iter->anchor_grid_pointer->row_position << endl;
					cout << "after anchoring" << endl;
					for (int i = 0; i < column_num; i++)
					{
					for (int j = 0; j < row_num; j++)
					{
					cout << rule_structure.Mark[i][j] << " ";
					}
					cout << endl;
					}



					cout << "1: " << &*anchoring_iter->polymer_unit_pointer << endl;*/
				break;

				
			}
			
			

		}
		else
		{
			anchoring_iter++;
			count++;
		}
	}
	//anchoring_iter = rule_structure.anchor_list.begin();

	//uniform_real_distribution<> tempf(0, 1);//here is for the time part
	//double temp = 0;
	//do
	//{
	//	temp = tempf(gen);
	//} while (temp == 0);

	//double deltat = (1.0 / rule_parameters.sumps)*log(1.0 / temp);
	//rule_parameters.t = rule_parameters.t + deltat;
	rule_parameters.rotation_sum++;
	rule_parameters.empty_anchor_unit_num--;
	//cout << rule_parameters.empty_anchor_unit_num << endl;
	//cout << rule_parameters.t << endl;
	//{
	//	int anchored_anchor_num = Anchor_num - rule_parameters.empty_anchor_unit_num;
	//	rule_parameters.ps[rule_deanchoring] = anchored_anchor_num*deanchoring_rate;
	//}//add deanchroing section
	/*rule_parameters.ps[rule_polymerize] = rule_parameters.polymerize_end_sum*polymerize_rate;
	rule_parameters.ps[rule_anchoring] = rule_parameters.empty_anchor_unit_num*anchoring_rate;
	rule_parameters.ps[rule_diffusion_polymer] = rule_parameters.polymerize_end_sum*diffusion_rate;
	rule_parameters.ps[rule_annealing_TT] = rule_structure.TT_anealing_polymer_unit_pair_list.size()*TT_annealing_rate;
	rule_parameters.ps[rule_fragmentate_TT] = rule_parameters.TT_fragmentation_sum*TT_fragmentation_rate;
	rule_parameters.ps[rule_deanchoring] = (Anchor_num - rule_parameters.empty_anchor_unit_num)*deanchoring_rate;*/

	//check_sumps();
	//cout << "anchoring finished" << endl;
	//output();
}

void deanchoring()
{
	//cout << " enter deanchoring function" << endl;
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
	uniform_int_distribution<> DeAnchor_Rand(0, rule_parameters.Anchor_num-rule_parameters.empty_anchor_unit_num);
	int deanchor_rand = DeAnchor_Rand(gen);
	list<Anchor_Unit>::iterator deanchoring_iter;
	deanchoring_iter = rule_structure.anchor_list.begin();
	int deanchoring_temp_sum = 0;
	
	int count = 0;
	//for (int i = 0; i < rule_structure.anchor_list.size(); i++)
	while (deanchoring_iter != rule_structure.anchor_list.end())
	{
		

		if (&*deanchoring_iter->polymer_unit_pointer != &*rule_structure.polymer_unit_list_end.begin())//when it is not empty
		{

			if (deanchoring_temp_sum < deanchor_rand-1)
			{
				deanchoring_temp_sum++;
				deanchoring_iter++;
				
			}
			else
			{//do the deanchoring stuff here
				//i have to judge whether there is only one anchor here
				if (deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 1)
				{//if the whole bundle has only one anchor, then the whole bundle has to be removed
					list<Polymer_Bundle>::iterator polymer_bundle_iter = deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
					/*list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter;
					list<Polymer_Unit>::iterator polymer_unit_iter;
					polymer_cluster_iter = polymer_bundle_iter->bundle_cluster_pointer_list.begin();
					while (polymer_cluster_iter != polymer_bundle_iter->bundle_cluster_pointer_list.end())
					{

					}*/

					//Maybe it is just:
					//set_action_mark_true(polymer_bundle_iter);
					if (deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->polymer_sequence.size() == 1 &&
						deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
					{
						rule_parameters.rotation_sum--;
					}
					deanchoring_iter->polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
					delete_polymer_bundle(polymer_bundle_iter);
					//rule_structure.polymer_bundle_list.erase(polymer_bundle_iter);
					
					check_anchor_diffusion(deanchoring_iter);
					//rule_parameters.empty_anchor_diffusion_direction_sum += deanchoring_iter->available_polymer_diffusion_direction_sum;
					rule_parameters.empty_anchor_unit_num++;
					
					//rule_parameters.emtpy_anchor_attatch_sum++;//add this because it will cause a error in 
					//deconstruction polymer_unit fuction. There i decrease every polymer_unit with anchor below it
					//but here is a special case

				}
				else
				{//i just need to remove this anchor
					list<Polymer_Bundle>::iterator polymer_bundle_iter = deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer;
					list<list<Anchor_Unit>::iterator>::iterator polymer_anchor_seq_iter;
					polymer_anchor_seq_iter = polymer_bundle_iter->polymer_anchor_list.begin();
					while (polymer_anchor_seq_iter != polymer_bundle_iter->polymer_anchor_list.end() && *polymer_anchor_seq_iter != deanchoring_iter)
					{//find this deanchoring iter from the bundle anchor list, erase it
						polymer_anchor_seq_iter++;
					}
					polymer_bundle_iter->polymer_anchor_list.erase(polymer_anchor_seq_iter);
					polymer_anchor_seq_iter = deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.begin();
					while (polymer_anchor_seq_iter != deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.end() && *polymer_anchor_seq_iter != deanchoring_iter)
					{//find this deanchoring iter from the cluster anchor list, erase it
						polymer_anchor_seq_iter++;
					}
					deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.erase(polymer_anchor_seq_iter);
					deanchoring_iter->polymer_unit_pointer->anchor_pointer = rule_structure.anchor_list.end();
					deanchoring_iter->polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
					rule_parameters.empty_anchor_unit_num++;
					//rule_parameters.emtpy_anchor_attatch_sum++;
					check_anchor_diffusion(deanchoring_iter);
					check_polymer_bundle_diffusable_direction(deanchoring_iter->anchor_grid_pointer->grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,2);
					//rule_parameters.empty_anchor_diffusion_direction_sum += deanchoring_iter->available_polymer_diffusion_direction_sum;

				}

				
				////if it is
				//if (true)//bundle has only this cluster
				//{
				//	if (true)//cluster has only this anchor
				//	{
				//		//delete bundle, cluster with every polymer_unit
				//		rule_structure.polymer_bundle_list.erase(deanchoring_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer);
				//		//I need to write a deconstruct function for all the classes that it will automaticly do every thing here.
				//	}
				//}

				//if the bundle has only this anchoring, then delete the whole bundle
				//if the bundle has multiple anchoring, then check for this cluster
				//if this cluster has only this anchoring, then delete this cluster, if not, just delete this cluster

				//I need to be carefull about delete a polymer_unit, i need to deal with boundary conditions:
				//delete all the neighbor_polymer_unit_links to it and also output ward;
				//free corresponding polymerize and diffusion
				break;

			}

		}
		else
		{
			deanchoring_iter++;
			count++;
		}
		
	}
	//deanchoring_iter = rule_structure.anchor_list.begin();
	
	//uniform_real_distribution<> tempf(0, 1);//here is for the time part
	//double temp = 0;
	//do
	//{
	//	temp = tempf(gen);
	//} while (temp == 0);

	//double deltat = (1.0 / rule_parameters.sumps)*log(1.0 / temp);
	//rule_parameters.t = rule_parameters.t + deltat;

	//rule_parameters.empty_anchor_unit_num--;
	//cout << rule_parameters.empty_anchor_unit_num << endl;
	//cout << rule_parameters.t << endl;
	//{
	//	int anchored_anchor_num = Anchor_num - rule_parameters.empty_anchor_unit_num;
	//	rule_parameters.ps[rule_deanchoring] = anchored_anchor_num*deanchoring_rate;
	//}//add deanchroing section
	//rule_parameters.ps[rule_anchoring] = rule_parameters.empty_anchor_unit_num*anchoring_rate;
	//rule_parameters.ps[rule_deanchoring] = (Anchor_num-rule_parameters.empty_anchor_unit_num)*deanchoring_rate;



	//rule_parameters.ps[rule_polymerize] = rule_parameters.polymerize_end_sum*polymerize_rate;
	////rule_parameters.ps[rule_anchoring] = rule_parameters.empty_anchor_unit_num*anchoring_rate;
	//rule_parameters.ps[rule_diffusion_polymer] = rule_parameters.polymerize_end_sum*diffusion_rate;
	//rule_parameters.ps[rule_annealing_TT] = rule_structure.TT_anealing_polymer_unit_pair_list.size()*TT_annealing_rate;
	//rule_parameters.ps[rule_fragmentate_TT] = rule_parameters.TT_fragmentation_sum*TT_fragmentation_rate;
	//rule_parameters.ps[rule_deanchoring] = (Anchor_num - rule_parameters.empty_anchor_unit_num)*deanchoring_rate;

	//check_sumps();

	//rule_parameters.empty_anchor_unit_num++;
	//cout << "deanchoring finished" << endl;
	//output();
	//system("pause");
};
void anchor_attatch()
{
	/*//random_device rd;
	mt19937_64 gen(rd());*/
	//default_random_engine gen;
#ifdef DEBUG_RANDOM
	//random_device rd;
	//mt19937_64 gen(rd());
	default_random_engine gen(rule_parameters.step_num);;
#else
	random_device rd;
	mt19937_64 gen(rd());
	//	gen(rule_parameters.step_num);;
#endif // debug_random
	uniform_int_distribution<> Attatch_Anchor_Rand(1, rule_parameters.emtpy_anchor_attatch_sum);
	int attatch_anchor_rand = Attatch_Anchor_Rand(gen);
	list<Anchor_Unit>::iterator attatch_anchor_iter;
	attatch_anchor_iter = rule_structure.anchor_list.begin();
	int attatch_temp_sum = 1;
	//cout << " enter attatch function" << endl;

	int count = 0;
	
	while (attatch_anchor_iter != rule_structure.anchor_list.end())
	{
		

		if (&*attatch_anchor_iter->polymer_unit_pointer == &*rule_structure.polymer_unit_list_end.begin()
			&& &*Grid[attatch_anchor_iter->anchor_grid_pointer->col_position][attatch_anchor_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer!=&*rule_structure.polymer_unit_list_end.begin())//when it is empty
		{

			if (attatch_temp_sum < attatch_anchor_rand)
			{
				attatch_temp_sum++;
				attatch_anchor_iter++;
				
			}
			else
			{//do attatchment here
				attatch_anchor_iter->polymer_unit_pointer = Grid[attatch_anchor_iter->anchor_grid_pointer->col_position][attatch_anchor_iter->anchor_grid_pointer->row_position].grid_polymer_unit_pointer;
				attatch_anchor_iter->polymer_unit_pointer->polymer_cluster_pointer->anchor_pointer_list.push_front(attatch_anchor_iter);
				attatch_anchor_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.push_front(attatch_anchor_iter);
				attatch_anchor_iter->polymer_unit_pointer->anchor_pointer = attatch_anchor_iter;
				//attatch_anchor_iter->attatch_mark = 0;
				//attatch_anchor_iter->polymer_unit_pointer->attatch_mark = 0;
				//rule_parameters.emtpy_anchor_attatch_sum--;
				//rule_parameters.empty_anchor_diffusion_direction_sum -= attatch_anchor_iter->available_polymer_diffusion_direction_sum;
				check_anchor_diffusion(attatch_anchor_iter);
				break;

			}

		}
		else
		{
			attatch_anchor_iter++;
			count++;
		}
	}
	

	rule_parameters.empty_anchor_unit_num--;
	//rule_parameters.emtpy_anchor_attatch_sum--;
	check_polymer_bundle_diffusable_direction(attatch_anchor_iter->polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer,2);
	
	//cout << "attatch finished" << endl;
	
};