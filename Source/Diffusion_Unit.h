#pragma once
#include "Grid_Unit.h"
#include <vector>
using namespace std;
class Diffusion_Unit
{
public:
	Diffusion_Unit();
	Diffusion_Unit(Grid_Unit input_grid_unit, int input_direction);
	Grid_Unit grid_unit;//deep copy?? Do i need a pointer here?
	int direction;

};