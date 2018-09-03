#include "Rule_Diffusion.h"
#include "Grid_Unit.h"
#include <list>
#include "Rule_Structure.h"
#include "Rule_Parameters.h"
#include <cmath>
#include <random>
#include "parameters.h"
#include "Rule_Bundling_Debundling.h"
#include "General_Functions.h"
#include "DeconstructionFunction.h"
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;

void fragmentation_TT()
{
	//cout << " fragmentation_TT" << endl;
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
	uniform_int_distribution<> Fragmentation_Rand(1, rule_parameters.TT_fragmentation_sum);
	int fragmentation_rand = Fragmentation_Rand(gen);
	list<Polymer_Cluster>::iterator fragmentation_polymer_iter;
	int fragmentation_temp_sum = 0;
	fragmentation_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (fragmentation_temp_sum < fragmentation_rand)//use <= because there are cases that witnin_polymer_sum=0;//don't remember
	{
		if (fragmentation_polymer_iter == rule_structure.polymer_cluster_list.end())
		{
			cout << "error: iterator out of range" << endl;
			system("pause");
		}
		fragmentation_temp_sum = fragmentation_temp_sum + fragmentation_polymer_iter->TT_sum;
		fragmentation_polymer_iter++;
	}
	fragmentation_polymer_iter--;
	//int within_polymer_sum = fragmentation_temp_sum-fragmentation_rand;
	int within_polymer_sum = fragmentation_rand-(fragmentation_temp_sum - fragmentation_polymer_iter->TT_sum);
	
	//if ((fragmentation_polymer_iter->polymer_sequence.size() - 3) < within_polymer_sum)
	//{//
	//	cout << "error: within polymer sum larger than sequence size" << endl;
	//	system("pause");
	//}
	/*if (within_polymer_sum != fragmentation_polymer_iter->polymer_sequence.size() - 3)
	{
		cout << "error: fragmentation position number error" << endl;
		system("pause");
	}*/

	int within_polymer_sum_temp = 0;
	int double_sum;
	int step_count = 1;
	list<Polymer_Unit>::iterator polymer_iter;
	//list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;

	polymer_iter = fragmentation_polymer_iter->polymer_sequence.begin();
	polymer_iter++;
	polymer_iter++;
	polymer_iter_back = (++fragmentation_polymer_iter->polymer_sequence.begin());

	rule_parameters.fragmentation_statistics[polymer_iter->polymer_grid_pointer->row_position]++;

	//polymer_iter_back;
	//polymer_iter++;
	//polymer_iter_front = polymer_iter;
	//polymer_iter_front++;
	//back>>current>>front
	if (fragmentation_polymer_iter->ring_mark == 1)
	{
		polymer_iter--;
		polymer_iter_back--;//i decrease because it has to be different with depolymerize in normal case, here it doesn't hold any more
		if ((--fragmentation_polymer_iter->polymer_sequence.end())->Hydrolysis_mark == T&&polymer_iter_back->Hydrolysis_mark == T)
		{
			if (within_polymer_sum == 1)
			{//means the fragmentation point is just in the place
				fragmentation_polymer_iter->ring_mark = 0;
				check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
				return;
			}
			
			within_polymer_sum_temp++;
		}
		
		//if (within_polymer_sum == column_num) 
		//{//if the break need to happen in the ring break, don't need to fragment
		//	fragmentation_polymer_iter->ring_mark = 0;
		//	check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
		//}
		//else
		{//not in the already exist ring break, need to reconstruc this polymer_sequence
			while (within_polymer_sum_temp<within_polymer_sum)
			{
				if (polymer_iter->Hydrolysis_mark == T)
				{
					if (polymer_iter_back->Hydrolysis_mark == T)
					{
						within_polymer_sum_temp++;
					}
				}
				polymer_iter++;
				polymer_iter_back++;
				//polymer_iter_front++;
			}
			polymer_iter--;
			polymer_iter_back--;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_new;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_new);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_new = rule_structure.polymer_cluster_list.begin();

			polymer_cluster_iter_new->direction = fragmentation_polymer_iter->direction;
			polymer_cluster_iter_new->polarize = fragmentation_polymer_iter->polarize;
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
			do
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_new;
				polymer_unit_iter_temp++;
			} 
			while (polymer_unit_iter_temp != fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.end(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence.end());
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			//rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			//list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			//recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			polymer_cluster_iter_new->cluster_bundle_pointer = fragmentation_polymer_iter->cluster_bundle_pointer;
			fragmentation_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.push_front(polymer_cluster_iter_new);
			recursion_polymer_unit(polymer_iter, polymer_cluster_iter_new->cluster_bundle_pointer);
			fragmentation_polymer_iter->ring_mark = 0;
			delete_polymer_cluster(fragmentation_polymer_iter);
			polymer_cluster_iter_new->ring_mark = 0;
			check_polymer_hydrolysis_sum(polymer_cluster_iter_new);
			
			//polymer_cluster_iter_new->ring_mark = 0;
			
			
		}
	}
	else
	{//if it is not in the ring case, do usual fragmentation


	
		for (int i = 0; i <= 3; i++)
		{
			if (fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
			{
				rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(fragmentation_polymer_iter->cluster_bundle_pointer, i);
				fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
				fragmentation_polymer_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
				
			}
		}



		while (within_polymer_sum_temp<within_polymer_sum)
		{
			double_sum = polymer_iter_back->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
			/*if ((polymer_iter->Hydrolysis_mark == T && polymer_iter_back->Hydrolysis_mark == D) || (polymer_iter->Hydrolysis_mark == D && polymer_iter_back->Hydrolysis_mark == T))
			{

			within_polymer_sum_temp++;

			}*/

			if (double_sum == TT)
			{
				within_polymer_sum_temp++;
			}
			polymer_iter++;
			polymer_iter_back++;
			step_count++;
			//polymer_iter_front++;
		}
		polymer_iter--;
		polymer_iter_back--;

		if (step_count < fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation++;
		}
		else if (step_count > fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation--;
		}



		//now polymer_iter is at TT position, 
		//						  ^
		//					      |
		//check_polymer_hydrolysis_sum(polymer_iter);
		Polymer_Bundle polymer_bundle;
		Polymer_Cluster polymer_cluster_direction_first;
		rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
		list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
		list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
		while (polymer_unit_iter_temp != polymer_iter)
		{
			polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			polymer_unit_iter_temp++;
		}
		polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), polymer_iter);
		//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
		//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
		rule_structure.polymer_bundle_list.push_front(polymer_bundle);
		list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		list<Polymer_Bundle>::iterator polymer_bundle_iter_store = fragmentation_polymer_iter->cluster_bundle_pointer;
		//use this to store this bundle_iter, will use it to delete as it will be the old no use one
		//list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		//put every polymer_cluster_iter pointer from the old polymer_bundle to the new polymer_bundle
		polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(),polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
		polymer_cluster_iter_direction_first->direction = fragmentation_polymer_iter->direction;
		polymer_cluster_iter_direction_first->polarize = fragmentation_polymer_iter->polarize;
		//use polymer_iter_back to recursion the new polymer_bundle_iter, if it is the same bundle, the other side polymer_iter should also point to polymer_bundle_iter
		recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
		
		if (polymer_bundle_iter != polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer)
		{//if they do not belong to the same bundle
			recursion_polymer_unit(polymer_iter, polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);
			if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->polymer_anchor_list.size() == 0)
			{
				/*if (polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
				{
					if (polymer_iter->polymer_cluster_pointer->polymer_sequence.size() == 1)
					{
						rule_parameters.rotation_sum--;
					}
				}*/
				//check_polymer_bundle_diffusable_direction(polymer_iter->polymer_cluster_pointer->cluster_bundle_pointer);//actually this command is to delete the empty cluster;
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
			for (int i = 0; i < 4; i++)
			{
				rule_parameters.polymer_diffusion_propensity +=calculate_diffusion_rate(polymer_bundle_iter,i);
			}
			
			//now polymer_iter_back and polymer_iter belongs to two different polymer_cluster, but they belong to the same bundle, so only need to diffusion detection once
			check_polymer_hydrolysis_sum(polymer_iter_back);
			check_polymer_hydrolysis_sum(polymer_iter);
			check_polymer_bundle_diffusable_direction(polymer_bundle_iter,2);


		}
	}
	
	//rebuild_TT_Fragmentation_sum();
	//i used to have it because there is some unkonw error here, now i need to debug it
	//cout << "TT_fragmentatoin finished" << endl;



	//Polymer_Cluster sliced_polymer_cluster;
	//Polymer_Bundle polymer_bundle;
	//sliced_polymer_cluster.polymer_sequence.insert(sliced_polymer_cluster.polymer_sequence.begin(), polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());
	////need to move all the pointers to this new polymer_cluster
	////grid to polymer_unit
	//list<Polymer_Unit>::iterator unit_iter;
	//rule_structure.polymer_cluster_list.push_front(sliced_polymer_cluster);
	//rule_structure.polymer_bundle_list.push_front(polymer_bundle);
	//list<Polymer_Bundle>::iterator new_bundle_iter;
	//list<Polymer_Cluster>::iterator sliced_cluster_iter;
	//new_bundle_iter = rule_structure.polymer_bundle_list.begin();
	//sliced_cluster_iter = rule_structure.polymer_cluster_list.begin();
	////cluster_iter--;
	//
	//sliced_cluster_iter->polymer_sequence.splice(sliced_cluster_iter->polymer_sequence.begin(),fragmentation_polymer_iter->polymer_sequence,polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());


	//unit_iter = sliced_polymer_cluster.polymer_sequence.begin();//problems here, should not only this. Here i want to add a new list to save the latter list to a new polymer_cluster with new bundle
	//do {
	//	unit_iter->polymer_cluster_pointer = sliced_cluster_iter;
		//unit_iter->polymer_grid_pointer->grid_polymer_unit_pointer = unit_iter;
		////grid to polymer_cluster
		//unit_iter->polymer_grid_pointer->grid_polymer_unit_pointer->polymer_cluster_pointer = cluster_iter;
		////polymer_unit to polymer_cluster
		//unit_iter->polymer_cluster_pointer = cluster_iter;
		////grid to polymer_bundle//i should cancel this 
		////unit_iter->polymer_grid_pointer->grid_polymer_bundle_pointer = bundle_iter;
		//unit_iter->polymer_grid_pointer->grid_polymer_unit_pointer->polymer_cluster_pointer->cluster_bundle_pointer = bundle_iter;
		////bundle pointer to this iter
		//unit_iter->bundled_neighbor_polymer_unit_pointer[0]->bundled_neighbor_polymer_unit_pointer[1] = unit_iter;
		//unit_iter->bundled_neighbor_polymer_unit_pointer[1]->bundled_neighbor_polymer_unit_pointer[0] = unit_iter;
		//use different 0 and 1 to represent opposite direction.
	//	unit_iter++;

	//} while (unit_iter != sliced_polymer_cluster.polymer_sequence.end());
	////bundle pointer to unit
	
	//anchor to polymer_unit//this will be built by recusion_polymer_unit
	//anchor to polymer_cluster//this will be built by recusion_polymer_unit
	//still need to remove those at some stage.
	//fragmentation_polymer_iter->polymer_sequence.erase(polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());

	//now i am not sure which position polymer_iter at, but i don't need to worry about it.
	//polymer_iter--;
	//recursion_polymer_unit(polymer_iter,bundle_iter);
	//recursion_polymer_unit(fragmentation_polymer_iter->polymer_sequence.begin(), fragmentation_polymer_iter->cluster_bundle_pointer);
	//check_polymer_bundle_diffusable_direction(fragmentation_polymer_iter->cluster_bundle_pointer);

	//if (fragmentation_polymer_iter->cluster_bundle_pointer != sliced_cluster_iter->cluster_bundle_pointer)
	//{//the fragmented one is not point to the same bundle, i need to construct a new bundle to it
	//	recursion_polymer_unit(sliced_cluster_iter->polymer_sequence.begin(), new_bundle_iter);
	//	check_polymer_bundle_diffusable_direction(sliced_cluster_iter->cluster_bundle_pointer);
	//}





	
	// I need to read this 


	//now polymer_iter is at TT position, 
	//						 ^
	//		
	//check hydrolysis_sum
	//check_polymer_hydrolysis_sum(polymer_iter); 
	//check_polymer_hydrolysis_sum(fragmentation_polymer_iter->polymer_sequence.begin());
	//check_polymer_hydrolysis_sum(sliced_cluster_iter->polymer_sequence.begin());
	////check boundary conditions
	//check bundling condiiton
	//check diffusion condition
	//check_polymer_bundle_diffusable_direction(polymer_iter);
	

	//polymer_iter = sliced_polymer_cluster.polymer_sequence.begin();
	//now polymer_iter is at TT position, 
	//						  ^
	//			
	//check sum
	//check_polymer_hydrolysis_sum(polymer_iter);
	//check boundary conditions|
	//recursion_polymer_unit(polymer_iter,bundle_iter);
	//check diffusion directions
	//check_polymer_bundle_diffusable_direction(polymer_iter);



};

void fragmentation_TD()
{
	//cout << " fragmentation_TT" << endl;
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
	uniform_int_distribution<> Fragmentation_Rand(1, rule_parameters.TD_fragmentation_sum);
	int fragmentation_rand = Fragmentation_Rand(gen);
	list<Polymer_Cluster>::iterator fragmentation_polymer_iter;
	int fragmentation_temp_sum = 0;
	fragmentation_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (fragmentation_temp_sum < fragmentation_rand)//use <= because there are cases that witnin_polymer_sum=0;
	{
		if (fragmentation_polymer_iter == rule_structure.polymer_cluster_list.end())
		{
			cout << "error: iterator out of range" << endl;
			system("pause");
		}
		fragmentation_temp_sum = fragmentation_temp_sum + fragmentation_polymer_iter->TD_sum;
		fragmentation_polymer_iter++;
	}
	fragmentation_polymer_iter--;
	int within_polymer_sum = fragmentation_rand - (fragmentation_temp_sum - fragmentation_polymer_iter->TD_sum);
	//if ((fragmentation_polymer_iter->polymer_sequence.size() - 3) < within_polymer_sum)
	//{//
	//	cout << "error: within polymer sum larger than sequence size" << endl;
	//	system("pause");
	//}


	int within_polymer_sum_temp = 0;
	int double_sum;
	list<Polymer_Unit>::iterator polymer_iter;
	//list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;

	polymer_iter = fragmentation_polymer_iter->polymer_sequence.begin();
	polymer_iter++;
	polymer_iter++;
	polymer_iter_back = (++fragmentation_polymer_iter->polymer_sequence.begin());

	rule_parameters.fragmentation_statistics[polymer_iter->polymer_grid_pointer->row_position]++;

	int step_count = 1;
	//polymer_iter_back;
	//polymer_iter++;
	//polymer_iter_front = polymer_iter;
	//polymer_iter_front++;
	//back>>current>>front
	if (fragmentation_polymer_iter->ring_mark == 1)
	{
		polymer_iter--;
		polymer_iter_back--;


		if (((--fragmentation_polymer_iter->polymer_sequence.end())->Hydrolysis_mark == T&&polymer_iter_back->Hydrolysis_mark == D)||
			((--fragmentation_polymer_iter->polymer_sequence.end())->Hydrolysis_mark == D&&polymer_iter_back->Hydrolysis_mark == T))
		{
			if (within_polymer_sum == 1)
			{//means the fragmentation point is just in the place
				fragmentation_polymer_iter->ring_mark = 0;
				check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
				return;
			}

			within_polymer_sum_temp++;
		}


		//if (within_polymer_sum == column_num)
		//{//if the break need to happen in the ring break, don't need to fragment
		//	fragmentation_polymer_iter->ring_mark = 0;
		//	check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
		//}
		//else
		{//not in the already exist ring break, need to reconstruc this polymer_sequence
			while (within_polymer_sum_temp<within_polymer_sum)
			{
				if ((polymer_iter->Hydrolysis_mark == T&&polymer_iter_back->Hydrolysis_mark == D)|| (polymer_iter->Hydrolysis_mark == D&&polymer_iter_back->Hydrolysis_mark == T))
				{
					
					within_polymer_sum_temp++;
					
				}
				polymer_iter++;
				polymer_iter_back++;
				//polymer_iter_front++;
			}
			polymer_iter--;
			polymer_iter_back--;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_new;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_new);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_new = rule_structure.polymer_cluster_list.begin();
			polymer_cluster_iter_new->direction = fragmentation_polymer_iter->direction;
			polymer_cluster_iter_new->polarize = fragmentation_polymer_iter->polarize;
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
			do
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_new;
				polymer_unit_iter_temp++;
			} while (polymer_unit_iter_temp != fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.end(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence.end());
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			//rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			//list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			//recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			polymer_cluster_iter_new->cluster_bundle_pointer = fragmentation_polymer_iter->cluster_bundle_pointer;
			fragmentation_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.push_front(polymer_cluster_iter_new);
			recursion_polymer_unit(polymer_iter, polymer_cluster_iter_new->cluster_bundle_pointer);
			fragmentation_polymer_iter->ring_mark = 0;
			delete_polymer_cluster(fragmentation_polymer_iter);
			polymer_cluster_iter_new->ring_mark = 0;
			check_polymer_hydrolysis_sum(polymer_cluster_iter_new);

			//polymer_cluster_iter_new->ring_mark = 0;


		}
	}
	else
	{//if it is not in the ring case, do usual fragmentation

		for (int i = 0; i <= 3; i++)
		{
			if (fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
			{
				rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(fragmentation_polymer_iter->cluster_bundle_pointer, i);
				fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
				fragmentation_polymer_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
				
			}
		}

		while (within_polymer_sum_temp<within_polymer_sum)
		{
			double_sum = polymer_iter_back->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
			/*if ((polymer_iter->Hydrolysis_mark == T && polymer_iter_back->Hydrolysis_mark == D) || (polymer_iter->Hydrolysis_mark == D && polymer_iter_back->Hydrolysis_mark == T))
			{

				within_polymer_sum_temp++;

			}*/
			
			if (double_sum == TD || double_sum == DT)
			{
				within_polymer_sum_temp++;
			}
			
			polymer_iter++;
			polymer_iter_back++;
			step_count++;
			//polymer_iter_front++;
			if (polymer_iter == fragmentation_polymer_iter->polymer_sequence.end())
			{
				cout << "fragmentation error, do fragmentation on a depolymerize site" << endl;
				system("pause");
			}
		}
		
		polymer_iter--;
		polymer_iter_back--;
		if (step_count < fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation++;
		}
		else if (step_count > fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation--;
		}
		



		//now polymer_iter is at TT position, 
		//						  ^
		//					      |
		//check_polymer_hydrolysis_sum(polymer_iter);
		Polymer_Bundle polymer_bundle;
		Polymer_Cluster polymer_cluster_direction_first;
		rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
		list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
		list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
		while (polymer_unit_iter_temp != polymer_iter)
		{
			polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			polymer_unit_iter_temp++;
		}
		polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), polymer_iter);
		//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
		//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
		rule_structure.polymer_bundle_list.push_front(polymer_bundle);
		list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		list<Polymer_Bundle>::iterator polymer_bundle_iter_store = fragmentation_polymer_iter->cluster_bundle_pointer;
		//use this to store this bundle_iter, will use it to delete
		list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
		polymer_cluster_iter_direction_first->direction = fragmentation_polymer_iter->direction;
		polymer_cluster_iter_direction_first->polarize = fragmentation_polymer_iter->polarize;
		recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
		
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
			for (int i = 0; i < 4; i++)
			{
				rule_parameters.polymer_diffusion_propensity += calculate_diffusion_rate(polymer_bundle_iter,i);
			}
			
			check_polymer_hydrolysis_sum(polymer_iter_back);
			check_polymer_hydrolysis_sum(polymer_iter);
			check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


		}
	}
};

void fragmentation_DD()
{
	//cout << " fragmentation_TT" << endl;
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
	uniform_int_distribution<> Fragmentation_Rand(1, rule_parameters.DD_fragmentation_sum);
	int fragmentation_rand = Fragmentation_Rand(gen);
	list<Polymer_Cluster>::iterator fragmentation_polymer_iter;
	int fragmentation_temp_sum = 0;
	fragmentation_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (fragmentation_temp_sum < fragmentation_rand)//use <= because there are cases that witnin_polymer_sum=0;
	{
		if (fragmentation_polymer_iter == rule_structure.polymer_cluster_list.end())
		{
			cout << "error: iterator out of range" << endl;
			system("pause");
		}
		fragmentation_temp_sum = fragmentation_temp_sum + fragmentation_polymer_iter->DD_sum;
		fragmentation_polymer_iter++;
	}
	fragmentation_polymer_iter--;
	int within_polymer_sum = fragmentation_rand - (fragmentation_temp_sum - fragmentation_polymer_iter->DD_sum);
	//if ((fragmentation_polymer_iter->polymer_sequence.size() - 3) < within_polymer_sum)
	//{//
	//	cout << "error: within polymer sum larger than sequence size" << endl;
	//	system("pause");
	//}


	int within_polymer_sum_temp = 0;
	int double_sum;
	list<Polymer_Unit>::iterator polymer_iter;
	//list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;

	polymer_iter = fragmentation_polymer_iter->polymer_sequence.begin();
	polymer_iter++;
	polymer_iter++;
	polymer_iter_back =(++fragmentation_polymer_iter->polymer_sequence.begin());

	rule_parameters.fragmentation_statistics[polymer_iter->polymer_grid_pointer->row_position]++;

	int step_count=1;
	//polymer_iter_front = polymer_iter;
	//polymer_iter_front++;
	//back>>current>>front
	if (fragmentation_polymer_iter->ring_mark == 1)
	{
		polymer_iter--;
		polymer_iter_back--;


		if ((--fragmentation_polymer_iter->polymer_sequence.end())->Hydrolysis_mark == D&&polymer_iter_back->Hydrolysis_mark == D)
		{
			if (within_polymer_sum == 1)
			{//means the fragmentation point is just in the place
				fragmentation_polymer_iter->ring_mark = 0;
				check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
				return;
			}

			within_polymer_sum_temp++;
		}



		//if (within_polymer_sum == column_num)
		//{//if the break need to happen in the ring break, don't need to fragment
		//	fragmentation_polymer_iter->ring_mark = 0;
		//	check_polymer_hydrolysis_sum(fragmentation_polymer_iter);
		//}
		//else
		{//not in the already exist ring break, need to reconstruc this polymer_sequence
			while (within_polymer_sum_temp<within_polymer_sum)
			{
				if ((polymer_iter->Hydrolysis_mark == D&&polymer_iter_back->Hydrolysis_mark == D) )
				{

					within_polymer_sum_temp++;

				}
				polymer_iter++;
				polymer_iter_back++;
				//polymer_iter_front++;
			}
			polymer_iter--;
			polymer_iter_back--;
			//now polymer_iter is at TT position, 
			//						  ^
			//					      |
			//Polymer_Bundle polymer_bundle;
			Polymer_Cluster polymer_cluster_new;
			rule_structure.polymer_cluster_list.push_front(polymer_cluster_new);
			list<Polymer_Cluster>::iterator polymer_cluster_iter_new = rule_structure.polymer_cluster_list.begin();
			polymer_cluster_iter_new->direction = fragmentation_polymer_iter->direction;
			polymer_cluster_iter_new->polarize = fragmentation_polymer_iter->polarize;
			list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
			do
			{
				polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_new;
				polymer_unit_iter_temp++;
			} while (polymer_unit_iter_temp != fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, polymer_iter, fragmentation_polymer_iter->polymer_sequence.end());
			polymer_cluster_iter_new->polymer_sequence.splice(polymer_cluster_iter_new->polymer_sequence.end(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence.end());
			//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
			//rule_structure.polymer_bundle_list.push_front(polymer_bundle);
			//list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
			//recursion_polymer_unit(polymer_iter_back, polymer_bundle_iter);
			polymer_cluster_iter_new->cluster_bundle_pointer = fragmentation_polymer_iter->cluster_bundle_pointer;
			fragmentation_polymer_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.push_front(polymer_cluster_iter_new);
			recursion_polymer_unit(polymer_iter, polymer_cluster_iter_new->cluster_bundle_pointer);
			fragmentation_polymer_iter->ring_mark = 0;
			delete_polymer_cluster(fragmentation_polymer_iter);
			polymer_cluster_iter_new->ring_mark = 0;
			check_polymer_hydrolysis_sum(polymer_cluster_iter_new);

			//polymer_cluster_iter_new->ring_mark = 0;


		}
	}
	else
	{//if it is not in the ring case, do usual fragmentation



		for (int i = 0; i <= 3; i++)
		{
			if (fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
			{
				rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(fragmentation_polymer_iter->cluster_bundle_pointer, i);
				fragmentation_polymer_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
				fragmentation_polymer_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
				
			}
		}
		while (within_polymer_sum_temp<within_polymer_sum)
		{
			double_sum = polymer_iter_back->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
			/*if ((polymer_iter->Hydrolysis_mark == T && polymer_iter_back->Hydrolysis_mark == D) || (polymer_iter->Hydrolysis_mark == D && polymer_iter_back->Hydrolysis_mark == T))
			{

			within_polymer_sum_temp++;

			}*/

			if (double_sum == DD)
			{
				within_polymer_sum_temp++;
			}
			polymer_iter++;
			polymer_iter_back++;
			step_count++;
			//polymer_iter_front++;
		}
		polymer_iter--;
		polymer_iter_back--;
		if (step_count < fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation++;
		}
		else if (step_count > fragmentation_polymer_iter->polymer_sequence.size() / 2.0)
		{
			rule_parameters.up_larger_than_down_fragmentation--;
		}




		//now polymer_iter is at TT position, 
		//						  ^
		//					      |
		//check_polymer_hydrolysis_sum(polymer_iter);
		Polymer_Bundle polymer_bundle;
		Polymer_Cluster polymer_cluster_direction_first;
		rule_structure.polymer_cluster_list.push_front(polymer_cluster_direction_first);
		list<Polymer_Cluster>::iterator polymer_cluster_iter_direction_first = rule_structure.polymer_cluster_list.begin();
		list<Polymer_Unit>::iterator polymer_unit_iter_temp = fragmentation_polymer_iter->polymer_sequence.begin();
		while (polymer_unit_iter_temp != polymer_iter)
		{
			polymer_unit_iter_temp->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
			polymer_unit_iter_temp++;
		}
		polymer_cluster_iter_direction_first->polymer_sequence.splice(polymer_cluster_iter_direction_first->polymer_sequence.begin(), fragmentation_polymer_iter->polymer_sequence, fragmentation_polymer_iter->polymer_sequence.begin(), polymer_iter);
		//polymer_iter_back->polymer_cluster_pointer = polymer_cluster_iter_direction_first;
		//polymer_cluster_iter_direction_first->anchor_pointer_list.splice(polymer_cluster_iter_direction_first->anchor_pointer_list.begin(), fragmentation_polymer_iter->anchor_pointer_list, fragmentation_polymer_iter->anchor_pointer_list.begin(), polymer_iter_back);
		rule_structure.polymer_bundle_list.push_front(polymer_bundle);
		list<Polymer_Bundle>::iterator polymer_bundle_iter = rule_structure.polymer_bundle_list.begin();
		list<Polymer_Bundle>::iterator polymer_bundle_iter_store = fragmentation_polymer_iter->cluster_bundle_pointer;
		//use this to store this bundle_iter, will use it to delete
		list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
		list<list<Polymer_Cluster>::iterator>::iterator polymer_cluster_iter_iter;
		polymer_bundle_iter->bundle_cluster_pointer_list.splice(polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list, polymer_bundle_iter_store->bundle_cluster_pointer_list.begin(), polymer_bundle_iter_store->bundle_cluster_pointer_list.end());
		polymer_cluster_iter_direction_first->direction = fragmentation_polymer_iter->direction;
		polymer_cluster_iter_direction_first->polarize = fragmentation_polymer_iter->polarize;
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
			for (int i = 0; i < 4; i++)
			{
				rule_parameters.polymer_diffusion_propensity += calculate_diffusion_rate(polymer_bundle_iter,i);
			}
			
			check_polymer_hydrolysis_sum(polymer_iter_back);
			check_polymer_hydrolysis_sum(polymer_iter);
			check_polymer_bundle_diffusable_direction(polymer_bundle_iter, 2);


		}
	}
}
void annealing_TT()
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
	uniform_int_distribution<> Annealing_Rand(1, rule_parameters.TT_annealing_sum);
	int anneal_rand = Annealing_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator anneal_polymer_unit_pair_iter;
	anneal_polymer_unit_pair_iter = rule_structure.anealing_polymer_unit_pair_list.begin();
	int anneal_temp_sum = 0;
	//anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	//if (anneal_temp_sum < anneal_rand)
	//{
		while (anneal_temp_sum < anneal_rand)
		{
			if (anneal_polymer_unit_pair_iter->first->Hydrolysis_mark == T&&anneal_polymer_unit_pair_iter->second->Hydrolysis_mark == T)
			{
				anneal_temp_sum++;
			}
			//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum;
			anneal_polymer_unit_pair_iter++;
		}
		anneal_polymer_unit_pair_iter--;//why?
	//}
	
	//anneal_polymer_unit_iter--;
	list<Polymer_Cluster>::iterator lower_cluster_iter;//it can be on the downside and left side
	list<Polymer_Cluster>::iterator higher_cluster_iter;//it can be on the upperside and right side
	//list<Polymer_Unit>::iterator large_polymer_unit_begin;
	//list<Polymer_Unit>::iterator large_polymer_unit_end;
	//i already define the direction of pair, so don't need to determine directino here
	
	lower_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	higher_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;

	for (int i = 0; i <= 3; i++)
	{
		if (lower_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
		{
			rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(lower_cluster_iter->cluster_bundle_pointer, i);
			lower_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
			lower_cluster_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
			
		}
	}
	for (int i = 0; i <= 3; i++)
	{
		if (higher_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
		{
			rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(higher_cluster_iter->cluster_bundle_pointer, i);
			higher_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
			higher_cluster_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
			
		}
	}



	if (higher_cluster_iter->polymer_sequence.size() > lower_cluster_iter->polymer_sequence.size())
	{
		rule_parameters.up_larger_than_down_annealing++;
		rule_parameters.TT_annealing_count_up_larger_than_below++;
	}
	else if (higher_cluster_iter->polymer_sequence.size() < lower_cluster_iter->polymer_sequence.size())
	{
		rule_parameters.up_larger_than_down_annealing--;
		rule_parameters.TT_annealing_count_up_smaller_than_below++;
	}

	if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}
	if (anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (lower_cluster_iter == higher_cluster_iter)//to account for the ring case
	{
		//anneal_polymer_unit_pair_iter->second!=lower_cluster_iter->polymer_sequence.begin()
		/*cout << "error, try to anneal polymer_units belong to the same cluster" << endl;
		system("pause");*/
		lower_cluster_iter->ring_mark = 1;
		/*rule_parameters.polymerize_end_sum -= lower_cluster_iter->polymerize_end_sum;
		lower_cluster_iter->polymerize_end_sum = 0;
		for (int i = 0; i < 4; i++)
		{
			lower_cluster_iter->availabe_polymerize_direction[i] = 0;
		}*/
		//check_polymer_bundle_diffusable_direction(lower_cluster_iter->cluster_bundle_pointer);
		check_polymer_hydrolysis_sum(lower_cluster_iter->polymer_sequence.begin());
		//check_polymer_hydrolysis_sum(lower_cluster_iter);
		anneal_polymer_unit_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
		anneal_polymer_unit_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
		rule_structure.anealing_polymer_unit_pair_list.erase(anneal_polymer_unit_pair_iter);
		rule_parameters.TT_annealing_sum--;
	}
	//now judge which part is larger, and insert the samller part into the larger part
	//if (lower_cluster_iter->polymer_sequence.size()>higher_cluster_iter->polymer_sequence.size())
	else
	{//the lower part is larger//i don't judge it any more, because there is a lot of bundling need to do here, only one polymer_sequence makes no difference

		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = higher_cluster_iter->polymer_sequence.begin();
		do
		{
			polymer_unit_iter->polymer_cluster_pointer = lower_cluster_iter;
			polymer_unit_iter++;
		} while (polymer_unit_iter != higher_cluster_iter->polymer_sequence.end());
		lower_cluster_iter->polymer_sequence.splice(lower_cluster_iter->polymer_sequence.end(), higher_cluster_iter->polymer_sequence, higher_cluster_iter->polymer_sequence.begin(), higher_cluster_iter->polymer_sequence.end());
		lower_cluster_iter->anchor_pointer_list.splice(lower_cluster_iter->anchor_pointer_list.end(),higher_cluster_iter->anchor_pointer_list, higher_cluster_iter->anchor_pointer_list.begin(), higher_cluster_iter->anchor_pointer_list.end());

		//lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.splice(lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.end(), higher_cluster_iter->anchor_pointer_list, higher_cluster_iter->anchor_pointer_list.begin(), higher_cluster_iter->anchor_pointer_list.end());
		//i not only need to move polymer_unit insider this polymer_cluster, but i also need to consider rebuild the bundling relation
		//
		if (lower_cluster_iter->cluster_bundle_pointer != higher_cluster_iter->cluster_bundle_pointer)
		{//they do not belong to the same bundle,delete the upper bundle, link everything to the downside bundle
			list<Polymer_Bundle>::iterator polymer_bundle_iter = higher_cluster_iter->cluster_bundle_pointer;
			
			//delete_polymer_cluster(higher_cluster_iter);
			lower_cluster_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.splice(lower_cluster_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.end(), polymer_bundle_iter->bundle_cluster_pointer_list, polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter->bundle_cluster_pointer_list.end());
			delete_polymer_cluster(higher_cluster_iter);
			//lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.erase(lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.begin(), lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.end());
			recursion_polymer_unit(lower_cluster_iter->polymer_sequence.begin(), lower_cluster_iter->cluster_bundle_pointer);//this function is not tested yet;
			//delete_polymer_cluster(higher_cluster_iter);
			delete_polymer_bundle(polymer_bundle_iter);
		}
		else
		{//they belong to the same bundle, only need to delete this cluster
			//on this step there should already be nothing on high_cluster_iter
			delete_polymer_cluster(higher_cluster_iter);
			recursion_polymer_unit(lower_cluster_iter->polymer_sequence.begin(), lower_cluster_iter->cluster_bundle_pointer);
			//to rebuild this bundle
		}
		//delete_polymer_cluster(higher_cluster_iter);//this might cause problem
		check_polymer_bundle_diffusable_direction(lower_cluster_iter->cluster_bundle_pointer,2);
		check_polymer_hydrolysis_sum(lower_cluster_iter->polymer_sequence.begin());
		anneal_polymer_unit_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
		anneal_polymer_unit_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
		rule_structure.anealing_polymer_unit_pair_list.erase(anneal_polymer_unit_pair_iter);
		rule_parameters.TT_annealing_sum--;
		
	}

	//if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size()>=
	//	anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size())
	//{//the first polymer_cluster is larger than the second
	//	large_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	//	small_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;
	//	list<Polymer_Unit>::iterator polymer_unit_iter;
	//	polymer_unit_iter = small_cluster_iter->polymer_sequence.begin();
	//	do
	//	{
	//		polymer_unit_iter->polymer_cluster_pointer = large_cluster_iter;
	//		polymer_unit_iter++;
	//	} while (polymer_unit_iter != small_cluster_iter->polymer_sequence.end());
	//	//i already define the direction of pair, so don't need to determine directino here
	//	if (anneal_polymer_unit_pair_iter->first == large_cluster_iter->polymer_sequence.begin())
	//	{//small-large like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence,small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//	else
	//	{//large-small like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.end(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//		
	//	//now i need to consider the annealing direction
	//}
	//else
	//{//the first polymer_cluster is smaller than the second
	//	small_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	//	large_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;
	//	list<Polymer_Unit>::iterator polymer_unit_iter;
	//	polymer_unit_iter = small_cluster_iter->polymer_sequence.begin();
	//	do
	//	{
	//		polymer_unit_iter->polymer_cluster_pointer = large_cluster_iter;
	//		polymer_unit_iter++;
	//	} while (polymer_unit_iter != small_cluster_iter->polymer_sequence.end());
	//	//now i need to consider the annealing direction
	//	if (anneal_polymer_unit_pair_iter->second == large_cluster_iter->polymer_sequence.begin())
	//	{//small-large like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//	else
	//	{//large-small like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.end(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//}
	//check_polymer_hydrolysis_sum(large_cluster_iter->polymer_sequence.begin());
	//check_polymer_hydrolysis_sum(small_cluster_iter->polymer_sequence.begin());
	////remenber to delete those hydrolysis sum for the small one, the small one information will be deleted.
	////how can i delete the diffusion information for the small one?
	//check_polymer_bundle_diffusable_direction(large_cluster_iter->cluster_bundle_pointer);
	////check_polymer_bundle_diffusable_direction(small_cluster_iter->cluster_bundle_pointer);
	//cout << "annealing_TT finished" << endl;
};



void annealing_TD()
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
	uniform_int_distribution<> Annealing_Rand(1, rule_parameters.TD_annealing_sum);
	int anneal_rand = Annealing_Rand(gen);
	list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator anneal_polymer_unit_pair_iter;
	anneal_polymer_unit_pair_iter = rule_structure.anealing_polymer_unit_pair_list.begin();
	int anneal_temp_sum = 0;
	int double_sum;
	//anneal_polymer_unit_iter = rule_structure.polymer_cluster_list.begin();
	
	
		while (anneal_temp_sum < anneal_rand)
		{
			double_sum = anneal_polymer_unit_pair_iter->first->Hydrolysis_mark * 10 + anneal_polymer_unit_pair_iter->second->Hydrolysis_mark;
			if(double_sum==TD ||double_sum==DT)
			/*if (!(anneal_polymer_unit_pair_iter->first->Hydrolysis_mark == T && anneal_polymer_unit_pair_iter->second->Hydrolysis_mark == T) &&
				!(anneal_polymer_unit_pair_iter->first->Hydrolysis_mark == D && anneal_polymer_unit_pair_iter->second->Hydrolysis_mark == D))*/
			{
				anneal_temp_sum++;
			}
			//anneal_temp_sum = anneal_temp_sum + anneal_polymer_unit_iter->TTT_sum;
			anneal_polymer_unit_pair_iter++;
		}
		anneal_polymer_unit_pair_iter--;//why
	
	//anneal_polymer_unit_iter--;
	list<Polymer_Cluster>::iterator lower_cluster_iter;//it can be on the downside and left side
	list<Polymer_Cluster>::iterator higher_cluster_iter;//it can be on the upperside and right side
														//list<Polymer_Unit>::iterator large_polymer_unit_begin;
														//list<Polymer_Unit>::iterator large_polymer_unit_end;
														//i already define the direction of pair, so don't need to determine directino here

	lower_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	higher_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;

	for (int i = 0; i <= 3; i++)
	{
		if (lower_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
		{
			rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(lower_cluster_iter->cluster_bundle_pointer, i);
			lower_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
			lower_cluster_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
			
		}
	}
	for (int i = 0; i <= 3; i++)
	{
		if (higher_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] == true)
		{
			rule_parameters.polymer_diffusion_propensity -= calculate_diffusion_rate(higher_cluster_iter->cluster_bundle_pointer, i);
			higher_cluster_iter->cluster_bundle_pointer->availabe_diffusion_direction[i] = false;
			higher_cluster_iter->cluster_bundle_pointer->available_polymer_diffusion_direction_sum--;
			
		}
	}


	if (higher_cluster_iter->polymer_sequence.size() > lower_cluster_iter->polymer_sequence.size())
	{
		rule_parameters.up_larger_than_down_annealing++;
		if (double_sum == TD)
		{
			rule_parameters.TD_annealing_count_up_larger_than_below++;
		}
		else if (double_sum == DT)
		{
			rule_parameters.DT_annealing_count_up_larger_than_below++;
		}
	}
	else if (higher_cluster_iter->polymer_sequence.size() < lower_cluster_iter->polymer_sequence.size())
	{
		rule_parameters.up_larger_than_down_annealing--;
		if (double_sum == TD)
		{
			rule_parameters.TD_annealing_count_up_smaller_than_below++;
		}
		else if (double_sum == DT)
		{
			rule_parameters.DT_annealing_count_up_smaller_than_below++;
		}
	}
	

	if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}
	if (anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->cluster_bundle_pointer->bundle_cluster_pointer_list.size() == 1)
	{
		if (anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size() == 1)
		{
			rule_parameters.rotation_sum--;
		}
	}

	if (lower_cluster_iter == higher_cluster_iter)//to account for the ring case
	{
		//anneal_polymer_unit_pair_iter->second!=lower_cluster_iter->polymer_sequence.begin()
		/*cout << "error, try to anneal polymer_units belong to the same cluster" << endl;
		system("pause");*/
		lower_cluster_iter->ring_mark = 1;
		/*rule_parameters.polymerize_end_sum -= lower_cluster_iter->polymerize_end_sum;
		lower_cluster_iter->polymerize_end_sum = 0;
		for (int i = 0; i < 4; i++)
		{
		lower_cluster_iter->availabe_polymerize_direction[i] = 0;
		}*/
		//check_polymer_bundle_diffusable_direction(lower_cluster_iter->cluster_bundle_pointer);
		check_polymer_hydrolysis_sum(lower_cluster_iter->polymer_sequence.begin());
		//check_polymer_hydrolysis_sum(lower_cluster_iter);
		anneal_polymer_unit_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
		anneal_polymer_unit_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
		rule_structure.anealing_polymer_unit_pair_list.erase(anneal_polymer_unit_pair_iter);
		rule_parameters.TD_annealing_sum--;
	}
	//now judge which part is larger, and insert the samller part into the larger part
	//if (lower_cluster_iter->polymer_sequence.size()>higher_cluster_iter->polymer_sequence.size())
	else
	{//the lower part is larger//i don't judge it any more, because there is a lot of bundling need to do here, only one polymer_sequence makes no difference

		list<Polymer_Unit>::iterator polymer_unit_iter;
		polymer_unit_iter = higher_cluster_iter->polymer_sequence.begin();
		do
		{
			polymer_unit_iter->polymer_cluster_pointer = lower_cluster_iter;
			polymer_unit_iter++;
		} while (polymer_unit_iter != higher_cluster_iter->polymer_sequence.end());
		lower_cluster_iter->polymer_sequence.splice(lower_cluster_iter->polymer_sequence.end(), higher_cluster_iter->polymer_sequence, higher_cluster_iter->polymer_sequence.begin(), higher_cluster_iter->polymer_sequence.end());
		lower_cluster_iter->anchor_pointer_list.splice(lower_cluster_iter->anchor_pointer_list.end(), higher_cluster_iter->anchor_pointer_list, higher_cluster_iter->anchor_pointer_list.begin(), higher_cluster_iter->anchor_pointer_list.end());

		//lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.splice(lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.end(), higher_cluster_iter->anchor_pointer_list, higher_cluster_iter->anchor_pointer_list.begin(), higher_cluster_iter->anchor_pointer_list.end());
		//i not only need to move polymer_unit insider this polymer_cluster, but i also need to consider rebuild the bundling relation
		//
		if (lower_cluster_iter->cluster_bundle_pointer != higher_cluster_iter->cluster_bundle_pointer)
		{//they do not belong to the same bundle,delete the upper bundle, link everything to the downside bundle
			list<Polymer_Bundle>::iterator polymer_bundle_iter = higher_cluster_iter->cluster_bundle_pointer;

			//delete_polymer_cluster(higher_cluster_iter);
			lower_cluster_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.splice(lower_cluster_iter->cluster_bundle_pointer->bundle_cluster_pointer_list.end(), polymer_bundle_iter->bundle_cluster_pointer_list, polymer_bundle_iter->bundle_cluster_pointer_list.begin(), polymer_bundle_iter->bundle_cluster_pointer_list.end());
			delete_polymer_cluster(higher_cluster_iter);
			//lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.erase(lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.begin(), lower_cluster_iter->cluster_bundle_pointer->polymer_anchor_list.end());
			recursion_polymer_unit(lower_cluster_iter->polymer_sequence.begin(), lower_cluster_iter->cluster_bundle_pointer);//this function is not tested yet;
																															 //delete_polymer_cluster(higher_cluster_iter);
			delete_polymer_bundle(polymer_bundle_iter);
		}
		else
		{//they belong to the same bundle, only need to delete this cluster
		 //on this step there should already be nothing on high_cluster_iter
			delete_polymer_cluster(higher_cluster_iter);
			recursion_polymer_unit(lower_cluster_iter->polymer_sequence.begin(), lower_cluster_iter->cluster_bundle_pointer);
			//to rebuild this bundle
		}
		//delete_polymer_cluster(higher_cluster_iter);//this might cause problem
		check_polymer_bundle_diffusable_direction(lower_cluster_iter->cluster_bundle_pointer, 2);
		check_polymer_hydrolysis_sum(lower_cluster_iter->polymer_sequence.begin());
		anneal_polymer_unit_pair_iter->first->anealing_polymer_unit_pair_list_iter[direction_second] = rule_structure.anealing_polymer_unit_pair_list.end();
		anneal_polymer_unit_pair_iter->second->anealing_polymer_unit_pair_list_iter[direction_first] = rule_structure.anealing_polymer_unit_pair_list.end();
		rule_structure.anealing_polymer_unit_pair_list.erase(anneal_polymer_unit_pair_iter);
		rule_parameters.TD_annealing_sum--;

	}

	//if (anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer->polymer_sequence.size()>=
	//	anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer->polymer_sequence.size())
	//{//the first polymer_cluster is larger than the second
	//	large_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	//	small_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;
	//	list<Polymer_Unit>::iterator polymer_unit_iter;
	//	polymer_unit_iter = small_cluster_iter->polymer_sequence.begin();
	//	do
	//	{
	//		polymer_unit_iter->polymer_cluster_pointer = large_cluster_iter;
	//		polymer_unit_iter++;
	//	} while (polymer_unit_iter != small_cluster_iter->polymer_sequence.end());
	//	//i already define the direction of pair, so don't need to determine directino here
	//	if (anneal_polymer_unit_pair_iter->first == large_cluster_iter->polymer_sequence.begin())
	//	{//small-large like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence,small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//	else
	//	{//large-small like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.end(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//		
	//	//now i need to consider the annealing direction
	//}
	//else
	//{//the first polymer_cluster is smaller than the second
	//	small_cluster_iter = anneal_polymer_unit_pair_iter->first->polymer_cluster_pointer;
	//	large_cluster_iter = anneal_polymer_unit_pair_iter->second->polymer_cluster_pointer;
	//	list<Polymer_Unit>::iterator polymer_unit_iter;
	//	polymer_unit_iter = small_cluster_iter->polymer_sequence.begin();
	//	do
	//	{
	//		polymer_unit_iter->polymer_cluster_pointer = large_cluster_iter;
	//		polymer_unit_iter++;
	//	} while (polymer_unit_iter != small_cluster_iter->polymer_sequence.end());
	//	//now i need to consider the annealing direction
	//	if (anneal_polymer_unit_pair_iter->second == large_cluster_iter->polymer_sequence.begin())
	//	{//small-large like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//	else
	//	{//large-small like configuration
	//		large_cluster_iter->polymer_sequence.splice(large_cluster_iter->polymer_sequence.end(), small_cluster_iter->polymer_sequence, small_cluster_iter->polymer_sequence.begin(), small_cluster_iter->polymer_sequence.end());
	//	}
	//}
	//check_polymer_hydrolysis_sum(large_cluster_iter->polymer_sequence.begin());
	//check_polymer_hydrolysis_sum(small_cluster_iter->polymer_sequence.begin());
	////remenber to delete those hydrolysis sum for the small one, the small one information will be deleted.
	////how can i delete the diffusion information for the small one?
	//check_polymer_bundle_diffusable_direction(large_cluster_iter->cluster_bundle_pointer);
	////check_polymer_bundle_diffusable_direction(small_cluster_iter->cluster_bundle_pointer);
	//cout << "annealing_TT finished" << endl;
};

