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
// I can accelate this when the length of the polymer is higher than 5, use subvector of length 5 to check the changes;

void Hydrolysis_TTT()
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
	uniform_int_distribution<> Hydrolysis_Rand(1, rule_parameters.TTT_sum_total);
	int hydrolysis_rand = Hydrolysis_Rand(gen);
	//int hydrolysis_rand = Hydrolysis_Rand();
	list<Polymer_Cluster>::iterator hydrolysis_polymer_iter;
	int hydrolysis_temp_sum=0;
	hydrolysis_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (hydrolysis_temp_sum < hydrolysis_rand)
	{
		hydrolysis_temp_sum = hydrolysis_temp_sum + hydrolysis_polymer_iter->TTT_sum;
		hydrolysis_polymer_iter++;
		
	}
	hydrolysis_polymer_iter--;
	//int within_polymer_sum = hydrolysis_temp_sum - hydrolysis_rand;
	int within_polymer_sum = hydrolysis_rand - (hydrolysis_temp_sum - hydrolysis_polymer_iter->TTT_sum);
	//current_temp_sum is the unit to have reaction within the polymer;
	//now we need to iterate the polymer for the hydrolysis;
	int within_polymer_sum_temp=0;
	list<Polymer_Unit>::iterator polymer_iter;
	list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;
	//need to consider different cases as the length is 1 or 2 or more than 2



	if(hydrolysis_polymer_iter->ring_mark == 1)
	{//if it is the ring case//has nothing to do with polymerize
		
		polymer_iter_front = hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter = ++hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
		advance(polymer_iter_back, 2);

		int triple_sum = --hydrolysis_polymer_iter->polymer_sequence.end()->Hydrolysis_mark * 100 + polymer_iter_front->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
		if (triple_sum == TTT)
		{
			within_polymer_sum_temp++;
		}
		if (within_polymer_sum_temp == within_polymer_sum)
		{
			polymer_iter_front->Hydrolysis_mark = D;
		}
		else
		{
			while (polymer_iter_back != hydrolysis_polymer_iter->polymer_sequence.end())
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + polymer_iter_back->Hydrolysis_mark;
				if (triple_sum == TTT)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					break;
				}
				polymer_iter_front++;
				polymer_iter++;
				polymer_iter_back++;
			}
			if (within_polymer_sum_temp == within_polymer_sum)
			{
				polymer_iter->Hydrolysis_mark = D;
			}
			else
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark;
				if (triple_sum == TTT)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					polymer_iter->Hydrolysis_mark = D;
				}
				else
				{
					cout << "hydrolysis, ring case, loop to the end but still not meet the within_polymer_sum condition" << endl;
					system("pause");
				}
			}
		}
	}
	else if (hydrolysis_polymer_iter->polymer_sequence.size() == 1)
	{
		if (hydrolysis_polymer_iter->polymer_sequence.begin()->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = hydrolysis_polymer_iter->polymer_sequence.begin()->anealing_polymer_unit_pair_list_iter[direction_first];

			if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
			{
				rule_parameters.TT_annealing_sum--;
				rule_parameters.TD_annealing_sum++;
			}
			else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
			{
				cout << "error, should not be in this situation" << endl;
				system("pause");
			}
			else
			{
				rule_parameters.TD_annealing_sum--;
			}


		}
		if (hydrolysis_polymer_iter->polymer_sequence.begin()->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
		{
			list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = hydrolysis_polymer_iter->polymer_sequence.begin()->anealing_polymer_unit_pair_list_iter[direction_second];

			if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
			{
				rule_parameters.TT_annealing_sum--;
				rule_parameters.TD_annealing_sum++;
			}
			else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
			{
				cout << "error, should not be in this situation" << endl;
				system("pause");
			}
			else
			{
				rule_parameters.TD_annealing_sum--;
			}


		}
		if (hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark == D)
		{
			cout << "hydrolysis error, try to hydrolyse D :1 TTT" << endl;
			system("pause");
		}
		hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark = D;
		
	}
	else
	{
		if (hydrolysis_polymer_iter->polymer_sequence.size() == 2)
		{
			if (within_polymer_sum == 0)
			{
				cout << "strange case" << endl;
				system("pause");
			}
			polymer_iter = ++hydrolysis_polymer_iter->polymer_sequence.begin();
			//polymer_iter++;
			polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
			//polymer_iter_back--;
			
			if (within_polymer_sum == 1)
			{
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D: 2.1 TTT" << endl;
					system("pause");
				}
				polymer_iter_back->Hydrolysis_mark = D;
			}
			if (within_polymer_sum == 2)
			{
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D :2.2 TTT" << endl;
					system("pause");
				}
				polymer_iter->Hydrolysis_mark = D;
			}
		}
		else
		{//length is longer than 2
			//cout << rule_parameters.t << endl;
			polymer_iter = hydrolysis_polymer_iter->polymer_sequence.begin();
			polymer_iter++;
			polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();;
			//polymer_iter_back--;
			polymer_iter_front = polymer_iter;
			polymer_iter_front++;
			//int mark = 0;
			//cout << polymer_iter->Hydrolysis_mark << endl;
			//int triple_sum;

			/*while (within_polymer_sum_temp<within_polymer_sum)
			{
				if (polymer_iter->Hydrolysis_mark == T)
				{
					if (polymer_iter_back->Hydrolysis_mark == T && polymer_iter_front->Hydrolysis_mark == T)
					{
						within_polymer_sum_temp++;
					}
				}
				polymer_iter++;
				polymer_iter_back++;
				polymer_iter_front++;
			}*/
			if (within_polymer_sum == 0)
			{
				cout << "strange case" << endl;
				system("pause");
			}
			if (polymer_iter_back->Hydrolysis_mark == T)
			{
				if (polymer_iter->Hydrolysis_mark == T)
				{
					within_polymer_sum_temp++;
				}
			}
			if (within_polymer_sum_temp == within_polymer_sum)
			{
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}

				if (polymer_iter_back->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D " << rule_parameters.t <<":"<<hydrolysis_polymer_iter->polymer_sequence.size() <<"."<<within_polymer_sum<< endl;
					system("pause");
				}

				polymer_iter_back->Hydrolysis_mark = D;
			}
			else
			{
				//while (polymer_iter != hydrolysis_polymer_iter->polymer_sequence.end()&&within_polymer_sum_temp!=within_polymer_sum)
				while (polymer_iter != hydrolysis_polymer_iter->polymer_sequence.end())// && within_polymer_sum_temp<within_polymer_sum)
				{
					if (polymer_iter->Hydrolysis_mark == T)
					{
						if (polymer_iter_back->Hydrolysis_mark == T && (polymer_iter_front->Hydrolysis_mark == T || polymer_iter_front == hydrolysis_polymer_iter->polymer_sequence.end()))
						{
							within_polymer_sum_temp++;
							if (within_polymer_sum_temp == within_polymer_sum)
							{
								break;
							}
						}
					}
					polymer_iter++;
					polymer_iter_back++;
					polymer_iter_front++;
				}
				//if (mark == 1)
				//{
				//polymer_iter--;
				//polymer_iter_back--;
				//polymer_iter_front--;
				//}



				//polymer_iter--;


				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}

				if (polymer_iter->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D" << rule_parameters.t << endl;
					system("pause");
				}

				polymer_iter->Hydrolysis_mark = D;
			}
			
			
		}
	}

	
	
	check_polymer_hydrolysis_sum(hydrolysis_polymer_iter);
	//i need to change all related relation here
	//maybe i have to write a check function here; discard two sentences after
	/*hydrolysis_polymer_iter->TTT_sum--;
	rule_parameters.TTT_sum_total--;*/
};
void Hydrolysis_TTD()
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
	uniform_int_distribution<> Hydrolysis_Rand(1, rule_parameters.TTD_sum_total);
	int hydrolysis_rand = Hydrolysis_Rand(gen);
	list<Polymer_Cluster>::iterator hydrolysis_polymer_iter;
	int hydrolysis_temp_sum = 0;
	hydrolysis_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (hydrolysis_temp_sum < hydrolysis_rand)
	{
		hydrolysis_temp_sum = hydrolysis_temp_sum + hydrolysis_polymer_iter->TTD_sum;
		hydrolysis_polymer_iter++;

	}
	hydrolysis_polymer_iter--;
	int within_polymer_sum = hydrolysis_rand - (hydrolysis_temp_sum - hydrolysis_polymer_iter->TTD_sum);
	//current_temp_sum is the unit to have reaction within the polymer;
	//now we need to iterate the polymer for the hydrolysis;
	int within_polymer_sum_temp = 0;
	list<Polymer_Unit>::iterator polymer_iter;
	list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;
	//need to consider different cases as the length is 1 or 2 or more than 2



	if (hydrolysis_polymer_iter->ring_mark == 1)
	{//if it is the ring case

		polymer_iter_front = hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter = ++hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
		advance(polymer_iter_back, 2);

		int triple_sum = --hydrolysis_polymer_iter->polymer_sequence.end()->Hydrolysis_mark * 100 + polymer_iter_front->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
		if (triple_sum == TTD ||DTT)
		{
			within_polymer_sum_temp++;
		}
		if (within_polymer_sum_temp == within_polymer_sum)
		{
			polymer_iter_front->Hydrolysis_mark = D;
		}
		else
		{
			while (polymer_iter_back != hydrolysis_polymer_iter->polymer_sequence.end())
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + polymer_iter_back->Hydrolysis_mark;
				if (triple_sum == TTD||DTT)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					break;
				}
				polymer_iter_front++;
				polymer_iter++;
				polymer_iter_back++;
			}
			if (within_polymer_sum_temp == within_polymer_sum)
			{
				polymer_iter->Hydrolysis_mark = D;
			}
			else
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark;
				if (triple_sum == TTD||DDT)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					polymer_iter->Hydrolysis_mark = D;
				}
				else
				{
					cout << "hydrolysis, ring case, loop to the end but still not meet the within_polymer_sum condition" << endl;
					system("pause");
				}
			}
		}
	}




	else if (hydrolysis_polymer_iter->polymer_sequence.size() == 1)
	{
		//only 1 length can't be in TTD case
		cout << "should not be in this case" << endl;
		system("pause");
		
		//hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark = D;

	}
	else
	{
		if (hydrolysis_polymer_iter->polymer_sequence.size() == 2)
		{
			if (within_polymer_sum!=1)
			{
				cout << "strange case" << endl;
				system("pause");
			}
			
			polymer_iter = hydrolysis_polymer_iter->polymer_sequence.begin();
			polymer_iter++;
			polymer_iter_back = polymer_iter;
			polymer_iter_back--;
			if (polymer_iter_back->Hydrolysis_mark == D)//DT|T
			{
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D" << endl;
					system("pause");
				}
				polymer_iter->Hydrolysis_mark = D;
			}		
			else if (polymer_iter->Hydrolysis_mark == D)//T|TD
			{
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D" << endl;
					system("pause");
				}
				polymer_iter_back->Hydrolysis_mark = D;
			}		
		}
		else
		{//length is longer than 2
			polymer_iter = hydrolysis_polymer_iter->polymer_sequence.begin();
			polymer_iter++;
			polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
			//polymer_iter_back--;
			polymer_iter_front = polymer_iter;
			polymer_iter_front++;
			int triple_sum;

			/*while (within_polymer_sum_temp<within_polymer_sum)
			{
			if (polymer_iter->Hydrolysis_mark == T)
			{
			if (polymer_iter_back->Hydrolysis_mark == T && polymer_iter_front->Hydrolysis_mark == T)
			{
			within_polymer_sum_temp++;
			}
			}
			polymer_iter++;
			polymer_iter_back++;
			polymer_iter_front++;
			}*/
			if (polymer_iter_back->Hydrolysis_mark == T)
			{
				if (polymer_iter->Hydrolysis_mark == D)
				{
					within_polymer_sum_temp++;
				}
			}
			if (within_polymer_sum_temp == within_polymer_sum)
			{//starting case, it should consider polymer_iter_back
				if (within_polymer_sum != 1)
				{
					cout << "strange case" << endl;
					system("pause");
				}
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter_back->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}

				if (polymer_iter_back->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D" << endl;
					system("pause");
				}
				polymer_iter_back->Hydrolysis_mark = D;
			}
			else
			{
				while (polymer_iter != hydrolysis_polymer_iter->polymer_sequence.end())// && within_polymer_sum_temp<within_polymer_sum)
				{
					if (polymer_iter_front == hydrolysis_polymer_iter->polymer_sequence.end())
					{
						triple_sum = polymer_iter_back->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + T;
					}
					else
					{
						triple_sum = polymer_iter_back->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + polymer_iter_front->Hydrolysis_mark;
					}
					if (triple_sum == TTD || triple_sum == DTT)
					{
						within_polymer_sum_temp++;
						if (within_polymer_sum_temp == within_polymer_sum)
						{
							break;
						}

					}
					polymer_iter++;
					polymer_iter_back++;
					polymer_iter_front++;
				}




				//polymer_iter--;

				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_first];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}
				if (polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second] != rule_structure.anealing_polymer_unit_pair_list.end())
				{
					list<pair<list<Polymer_Unit>::iterator, list<Polymer_Unit>::iterator>>::iterator annealing_pair_iter = polymer_iter->anealing_polymer_unit_pair_list_iter[direction_second];

					if (annealing_pair_iter->first->Hydrolysis_mark == T&&annealing_pair_iter->second->Hydrolysis_mark == T)
					{
						rule_parameters.TT_annealing_sum--;
						rule_parameters.TD_annealing_sum++;
					}
					else if (annealing_pair_iter->first->Hydrolysis_mark == D&&annealing_pair_iter->second->Hydrolysis_mark == D)
					{
						cout << "error, should not be in this situation" << endl;
						system("pause");
					}
					else
					{
						rule_parameters.TD_annealing_sum--;
					}


				}

				if (polymer_iter->Hydrolysis_mark == D)
				{
					cout << "hydrolysis error, try to hydrolyse D" << endl;
					system("pause");
				}
				polymer_iter->Hydrolysis_mark = D;
			}
			
			
		}
	}



	check_polymer_hydrolysis_sum(hydrolysis_polymer_iter);
	//i need to change all related relation here
	//maybe i have to write a check function here; discard two sentences after
	/*hydrolysis_polymer_iter->TTT_sum--;
	rule_parameters.TTT_sum_total--;*/
	
};
void Hydrolysis_DTD()
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
	uniform_int_distribution<> Hydrolysis_Rand(1, rule_parameters.DTD_sum_total);
	int hydrolysis_rand = Hydrolysis_Rand(gen);
	list<Polymer_Cluster>::iterator hydrolysis_polymer_iter;
	int hydrolysis_temp_sum = 0;
	hydrolysis_polymer_iter = rule_structure.polymer_cluster_list.begin();
	while (hydrolysis_temp_sum < hydrolysis_rand)
	{
		hydrolysis_temp_sum = hydrolysis_temp_sum + hydrolysis_polymer_iter->DTD_sum;
		hydrolysis_polymer_iter++;

	}
	hydrolysis_polymer_iter--;
	int within_polymer_sum = hydrolysis_rand - (hydrolysis_temp_sum - hydrolysis_polymer_iter->DTD_sum);
	//current_temp_sum is the unit to have reaction within the polymer;
	//now we need to iterate the polymer for the hydrolysis;
	int within_polymer_sum_temp = 0;
	list<Polymer_Unit>::iterator polymer_iter;
	list<Polymer_Unit>::iterator polymer_iter_front;
	list<Polymer_Unit>::iterator polymer_iter_back;
	//need to consider different cases as the length is 1 or 2 or more than 2



	if (hydrolysis_polymer_iter->ring_mark == 1)
	{//if it is the ring case

		polymer_iter_front = hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter = ++hydrolysis_polymer_iter->polymer_sequence.begin();
		polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
		advance(polymer_iter_back, 2);

		int triple_sum = --hydrolysis_polymer_iter->polymer_sequence.end()->Hydrolysis_mark * 100 + polymer_iter_front->Hydrolysis_mark * 10 + polymer_iter->Hydrolysis_mark;
		if (triple_sum == DTD)
		{
			within_polymer_sum_temp++;
		}
		if (within_polymer_sum_temp == within_polymer_sum)
		{
			polymer_iter_front->Hydrolysis_mark = D;
		}
		else
		{
			while (polymer_iter_back != hydrolysis_polymer_iter->polymer_sequence.end())
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + polymer_iter_back->Hydrolysis_mark;
				if (triple_sum == DTD)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					break;
				}
				polymer_iter_front++;
				polymer_iter++;
				polymer_iter_back++;
			}
			if (within_polymer_sum_temp == within_polymer_sum)
			{
				polymer_iter->Hydrolysis_mark = D;
			}
			else
			{
				triple_sum = polymer_iter_front->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark;
				if (triple_sum == DTD)
				{
					within_polymer_sum_temp++;
				}
				if (within_polymer_sum_temp == within_polymer_sum)
				{
					polymer_iter->Hydrolysis_mark = D;
				}
				else
				{
					cout << "hydrolysis, ring case, loop to the end but still not meet the within_polymer_sum condition" << endl;
					system("pause");
				}
			}
		}
	}




	else if (hydrolysis_polymer_iter->polymer_sequence.size() == 1)
	{
		//hydrolysis_polymer_iter->polymer_sequence.begin()->Hydrolysis_mark = D;
		cout << "should not be in this case" << endl;
		system("pause");
	}
	else
	{
		if (hydrolysis_polymer_iter->polymer_sequence.size() == 2)
		{
			/*polymer_iter = hydrolysis_polymer_iter->polymer_sequence.begin();
			polymer_iter++;
			polymer_iter_back = polymer_iter;
			polymer_iter_back--;
			if (within_polymer_sum == 1)
			{
				polymer_iter_back->Hydrolysis_mark = D;
			}
			if (within_polymer_sum == 2)
			{
				polymer_iter->Hydrolysis_mark = D;
			}*/
			cout << "should not be in this case" << endl;
			system("pause");
		}
		else
		{//length is longer than 2
			polymer_iter = hydrolysis_polymer_iter->polymer_sequence.begin();
			polymer_iter++;
			polymer_iter_back = hydrolysis_polymer_iter->polymer_sequence.begin();
			//polymer_iter_back--;
			polymer_iter_front = polymer_iter;
			polymer_iter_front++;
			int triple_sum;

			/*while (within_polymer_sum_temp<within_polymer_sum)
			{
			if (polymer_iter->Hydrolysis_mark == T)
			{
			if (polymer_iter_back->Hydrolysis_mark == T && polymer_iter_front->Hydrolysis_mark == T)
			{
			within_polymer_sum_temp++;
			}
			}
			polymer_iter++;
			polymer_iter_back++;
			polymer_iter_front++;
			}*/
			while (polymer_iter != hydrolysis_polymer_iter->polymer_sequence.end())
			{
				if (polymer_iter_front == hydrolysis_polymer_iter->polymer_sequence.end())
				{
					triple_sum = polymer_iter_back->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + T;
				}
				else
				{
					triple_sum = polymer_iter_back->Hydrolysis_mark * 100 + polymer_iter->Hydrolysis_mark * 10 + polymer_iter_front->Hydrolysis_mark;
				}
				if (triple_sum == DTD)
				{
					within_polymer_sum_temp++;
					if (within_polymer_sum_temp == within_polymer_sum)
					{
						break;
					}

				}
				polymer_iter++;
				polymer_iter_back++;
				polymer_iter_front++;
			}




			//polymer_iter--;
			if (polymer_iter->Hydrolysis_mark == D)
			{
				cout << "hydrolysis error, try to hydrolyse D" << endl;
				system("pause");
			}
			polymer_iter->Hydrolysis_mark = D;
		}
	}



	check_polymer_hydrolysis_sum(hydrolysis_polymer_iter);
	//i need to change all related relation here
	//maybe i have to write a check function here; discard two sentences after
	/*hydrolysis_polymer_iter->TTT_sum--;
	rule_parameters.TTT_sum_total--;*/
};