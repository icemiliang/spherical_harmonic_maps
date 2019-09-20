#include <fstream>
#include "OBJFileReader.h"
#include "Solid.h"
#include "harmonic.h"

using namespace MeshLib;
using namespace std;

int main(int argc, char *argv[]) {
	Solid nmesh;
	OBJFileReader of;
	ifstream in(argv[1]);
	of.readToSolid(&nmesh, in);

	Harmonic shm(&nmesh);
	cout << "--> Setting up" << endl;
	vector<double> angles_before = shm.compute_vertex_angles();
	shm.set_up();
	shm.star_map();
	string out = string(argv[2]) + "_star.obj";
	of.writeToObj(&nmesh, out);

	cout << "--> Computing Tuette map" << endl;
	shm.tuette_map();
	out = string(argv[2]) + "_tuette.obj";
	of.writeToObj(&nmesh, out);

	cout << "--> Computing harmonic map" << endl;
	// shm.harmonic_map_gd();
	shm.harmonic_map_conjugate_gd();
	out = string(argv[2]) + "_harmonic.obj";
	of.writeToObj(&nmesh, out);

	cout << "--> Checking angle distorsions" << endl;
	vector<double> angles_after = shm.compute_vertex_angles();
	double angle_before = accumulate(angles_before.begin(), angles_before.end(), 0);
	double angle_after = accumulate(angles_after.begin(), angles_after.end(), 0);
	cout << "--> Averaged angle difference: " 
		 << fabs(angle_after - angle_before) / double(angles_before.size()) 
		 << " degrees" << endl;

	return 0;
}
