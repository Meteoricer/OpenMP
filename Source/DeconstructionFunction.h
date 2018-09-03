#pragma once
#include "Grid_Unit.h"

void delete_polymer_unit(list<Polymer_Cluster>::iterator polymer_cluster_iter,list<Polymer_Unit>::iterator polymer_unit_iter);
void delete_polymer_unit(list<Polymer_Cluster>::iterator polymer_cluster_iter, list<Polymer_Unit>::iterator polymer_unit_iter_start, list<Polymer_Unit>::iterator polymer_unit_iter_end);
void delete_polymer_cluster(list<Polymer_Cluster>::iterator polymer_cluster_iter);
void delete_polymer_bundle(list<Polymer_Bundle>::iterator polymer_bundle_iter);
