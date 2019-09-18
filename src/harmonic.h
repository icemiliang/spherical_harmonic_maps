#include "Solid.h"
#include "iterators.h"
#include <math.h>
#include "Point.h"
#include "CTrait.h"

#define DELTA_T 1e-2
#define DELTA_E 1e-5

#define HARMONIC 1
#define TUETTE 0

using namespace MeshLib;


void star_map(Solid *cmesh);

double compute_energy(Solid *cmesh, int type);

void compute_gradient(Solid *cmesh, int type);

void conjugate_gradient(Solid *cmesh);

void update_mesh(Solid *cmesh);
void update_mass_center(Solid * cmesh);

void harmonic_map_conjugate_gd(Solid * pnmesh);

void harmonic_map_gd(Solid * pnmesh);
void tuette_map(Solid * pnmesh);
void set_up_normal(Solid *nmesh);

void set_up_kuv(Solid *nmesh);

void set_up(Solid * nmesh);
