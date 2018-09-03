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

using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;
class Polymer_Unit;
class Anchor_Unit;



void test_TT_Fragmentation_sum()
{
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int fragmentation_sum = 0;
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{

		int temp;
		if (polymer_cluster_iter->polymer_sequence.size() < 4)
		{
			temp = 0;
		}
		if (polymer_cluster_iter->polymer_sequence.size() >= 4)
		{
			temp = polymer_cluster_iter->polymer_sequence.size() - 3;
		}
		if (polymer_cluster_iter->ring_mark == 1)
		{
			temp = rule_parameters.column_num;
		}
		
		if (temp!= polymer_cluster_iter->TT_sum)
		{
			cout << "polymer_cluster fragmentation_TT sum doesn't match" << endl;
			system("pause");
		}


		fragmentation_sum += polymer_cluster_iter->TT_sum;
		polymer_cluster_iter++;
		
	}
	if (fragmentation_sum != rule_parameters.TT_fragmentation_sum)
	{
		cout << "fragmentation_TT sum doesn't match" << endl;
		system("pause");
	}
	

}


void rebuild_TT_Fragmentation_sum()
{
	list<Polymer_Cluster>::iterator polymer_cluster_iter = rule_structure.polymer_cluster_list.begin();
	int fragmentation_sum = 0;
	while (polymer_cluster_iter != rule_structure.polymer_cluster_list.end())
	{
		fragmentation_sum += polymer_cluster_iter->TT_sum;
		polymer_cluster_iter++;
	}
	rule_parameters.TT_fragmentation_sum = fragmentation_sum;

}