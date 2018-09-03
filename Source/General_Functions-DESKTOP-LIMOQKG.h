#pragma once
#include "Grid_Unit.h"
#include <fstream>

//void check_polymer_hydrolysis_sum(list<Polymer_Unit>::iterator input_polymer_iterator);
void select_and_execute();
void rebuild_polymer_bc(list<Polymer_Cluster>::iterator input_polymer_cluster_iter);//bc:boundary condition
void recursion_polymer_unit(list<Polymer_Unit>::iterator polymer_unit_iter, list<Polymer_Bundle>::iterator polymer_bundle_iter);
void check_polymer_hydrolysis_sum(list<Polymer_Unit>::iterator polymer_unit_iter);
void check_polymer_hydrolysis_sum(list<Polymer_Cluster>::iterator polymer_cluster_iter);
//void check_polymer_bundle_diffusable_direction(list<Polymer_Unit>::iterator polymer_unit_iter);
void check_polymer_bundle_diffusable_direction(list<Polymer_Bundle>::iterator diffusion_bundle_iter,int check_count);
void check_bundling_pairs(list<Polymer_Bundle>::iterator bundle_iter);
void output(ofstream&);
void output();
void select_and_execute();
void check();
void check_sumps();
void check_annealing(list<Polymer_Cluster>::iterator polymer_cluster_iter);
void check_bundling();
void set_action_mark_true(list < Polymer_Bundle>::iterator polymer_bundle_iter);
void set_action_mark_false(list < Polymer_Bundle>::iterator polymer_bundle_iter);  
void check_anchor_diffusion(list<Anchor_Unit>::iterator anchor_iter);
void calculate_Z_ctp();
void test_annealing();
void check_direction();
void check_debundling_information();
void check_polymerize_information();
void check_T_D_distance();
void check_length_distribution();
void check_average_length_per_column();
void frap(int left,int right);


