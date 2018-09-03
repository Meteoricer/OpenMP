#pragma once
#include <math.h>
using namespace std;


const int yes = 1;
const int no = 0;
const int direction_up = 0;
const int direction_down = 1;
const int direction_left = 2;
const int direction_right = 3;
const int NoDirections=0;
const int horizontal = 1;
const int vertical = 2;
const int direction_first = 0;//it can be left or down
const int direction_second = 1;//it can be up or right
const int NoFtsZ = 0;
const int T = 1;
const int D = 2;
const int A = 3;
const int AT = 4;
const int AD = 5;
const int TN = 6;//N means down
const int DN = 7;
const int ATN = 8;
const int ADN = 9;
const int effect_rate = 0;
const int NoAnchorBelow = 0;
const int AnchorBelowWithoutBound = 1;
const int AnchorBelowWithBound = 2;
const int NoLateralBound = 0;
const int LateralBoundDirectionFirst = 1;
const int LateralBoundDirectionSecond = 2;
const int LateralBoundDIrectionBoth = 3;
const int Head = 1;
const int Tail = 2;
const int NotHeadOrTail = 0;
const int HeadAndTail = 3;


				   //const double E_bdl = 0.03;//-0.07


const int TTT = T * 100 + T * 10 + T;
const int TTD = T * 100 + T * 10 + D;
const int TDT = T * 100 + D * 10 + T;
const int DTT = D * 100 + T * 10 + T;
const int TDD = T * 100 + D * 10 + D;
const int DTD = D * 100 + T * 10 + D;
const int DDT = D * 100 + D * 10 + T;
const int DDD = D * 100 + D * 10 + D;
const int TT = T * 10 + T;
const int TD = T * 10 + D;
const int DT = D * 10 + T;
const int DD = D * 10 + D;
//const double up = 1.1;
//const double down = 1.2;
//const double left = 1.3;
//const double right = 1.4;


const int num_of_rules = 32;
const int rule_diffusion_polymer = 0;
const int rule_TT_polymerize_horizontal_plus = 1, rule_depolymerize_plus_TT = 2, rule_depolymerize_plus_TD = 3, rule_depolymerize_plus_DD = 4;
const int rule_annealing_TT = 5, rule_annealing_TD = 6;
const int rule_fragmentate_TT = 7, rule_fragmentate_TD = 8, rule_fragmentate_DD = 9;
const int rule_bundling = 10, rule_debundling = 11, rule_attatch = 12, rule_anchoring = 13, rule_deanchoring = 14;
const int rule_hydrolysis_TTT = 15, rule_hydrolysis_TTD = 16, rule_hydrolysis_DTD = 17;
const int rule_empty_anchor_diffusion = 18;
const int rule_rotation = 19;
const int rule_TT_polymerize_vertical_plus = 20;
const int rule_bundling_inverse = 21;
const int rule_debundling_inverse = 22;
const int rule_TT_polymerize_horizontal_minus = 23;
const int rule_TT_polymerize_vertical_minus = 24;
const int rule_TD_polymerize_vertical_plus = 25;
const int rule_TD_polymerize_vertical_minus = 26;
const int rule_TD_polymerize_horizontal_plus = 27;
const int rule_TD_polymerize_horizontal_minus = 28;
const int rule_depolymerize_minus_TT = 29, rule_depolymerize_minus_TD = 30, rule_depolymerize_minus_DD = 31;
const int random_seed = 10;

