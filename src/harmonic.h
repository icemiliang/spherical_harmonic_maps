#include "Solid.h"
#include "iterators.h"
#include <math.h>
#include "Point.h"
#include "CTrait.h"
#include <vector>

#define DELTA_T 1e-2
#define DELTA_E 1e-5

#define HARMONIC 1
#define TUETTE 0

using namespace MeshLib;

class Harmonic {

public:
	Harmonic(Solid *cmesh) {
		_nmesh = cmesh;
	};
	void star_map();
	void harmonic_map_conjugate_gd();
	void harmonic_map_gd();
	void tuette_map();
	void set_up();
	std::vector<double> compute_vertex_angles();
protected:
	void set_up_normal();
	void set_up_kuv();
	double compute_energy(int type);
	void compute_gradient(int type);
	void conjugate_gradient();
	void update_mesh();
	void update_mass_center();

	Solid * _nmesh;
};
