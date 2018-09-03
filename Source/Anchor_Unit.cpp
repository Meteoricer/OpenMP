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
using namespace std;
extern vector<vector<Grid_Unit>> Grid;
extern Rule_Structure rule_structure;
extern Rule_Parameters rule_parameters;


Anchor_Unit::Anchor_Unit()
{
	this->polymer_unit_pointer = rule_structure.polymer_unit_list_end.begin();
	this->attatch_mark = 0;
	this->availabe_diffusion_direction.resize(4);
	this->availabe_diffusion_direction[direction_up] = false;
	this->availabe_diffusion_direction[direction_down] = false;
	this->availabe_diffusion_direction[direction_left] = false;
	this->availabe_diffusion_direction[direction_right] = false;

};