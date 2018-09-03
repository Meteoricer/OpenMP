#pragma once
#include <math.h>
using namespace std;
const int column_num = 600;//600//150
const int row_num = 800;//800//200
const int Anchor_num = 1000;//1000
const int yes = 1;
const int no = 0;
const int direction_up = 0;
const int direction_down = 1;
const int direction_left = 2;
const int direction_right = 3;
const int horizontal = 1;
const int vertical = 2;
const int direction_first=0;//it can be left or down
const int direction_second=1;//it can be up or right
const int T = 1;
const int D = 2;
const int Z = 8000;//8000//600//1000//4000
const double E_bdl = -0.07;
const int bending_energy = 1;
//const double up = 1.1;
//const double down = 1.2;
//const double left = 1.3;
//const double right = 1.4;

const double anchoring_rate = 0.000001;//100//200
const double deanchoring_rate = 0.2;
const double polymerize_vertical_rate=0.05;//0.05
const double polymerize_horizontal_rate = polymerize_vertical_rate*exp(-bending_energy);
const double diffusion_rate = 4000;//4000
const double emtpy_anchor_diffusion_rate = 0.1;
const double Ebdl = 0.07;
const double bundling_rate = 500;//500
const double debundling_rate = bundling_rate*exp(-Ebdl);
const double TT_annealing_rate = 2500;
const double TD_annealing_rate = 2500;
const double TT_fragmentation_rate = 8;
const double TD_fragmentation_rate = 8;
const double DD_fragmentation_rate = 8;
const double TT_depolymerize_rate = 80;//80
const double TD_depolymerize_rate = 80;
const double DD_depolymerize_rate = 80;
const double attatch_rate = 100;//0.000001;
const double both = 200;

const double TTT_hydrolysys_rate = 0;//100;//100;
const double TTD_hydrolysys_rate = 0;// 100;//100;
const double DTD_hydrolysys_rate = 0;// 100;//100;
const double rotation_rate = 0; 

const int num_of_rules=21;
const int rule_diffusion_polymer = 0;
const int rule_polymerize_horizontal = 1, rule_depolymerize_TT = 2, rule_depolymerize_TD = 3, rule_depolymerize_DD = 4;
const int rule_annealing_TT = 5, rule_annealing_TD = 6;
const int rule_fragmentate_TT = 7, rule_fragmentate_TD = 8, rule_fragmentate_DD = 9;
const int rule_bundling = 10, rule_debundling = 11, rule_attatch = 12, rule_anchoring = 13, rule_deanchoring = 14;
const int rule_hydrolysis_TTT = 15, rule_hydrolysis_TTD = 16, rule_hydrolysis_DTD = 17;
const int rule_empty_anchor_diffusion = 18;
const int rule_rotation = 19;
const int rule_polymerize_vertical = 20;
const int random_seed = 10;

const double posibility_vertical = 1 / (1 + exp(-bending_energy));
const double posibility_horizontal = exp(-bending_energy) / (1 + exp(-bending_energy));