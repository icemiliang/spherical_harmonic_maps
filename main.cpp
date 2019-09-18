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

	Solid * pnmesh = &nmesh;
	cout << "--> Setting up kuv" << endl;
	set_up(pnmesh);

	cout << "--> Gauss map" << endl;
	star_map(&nmesh);
	string out = string(argv[2]) + "_star.obj";
	of.writeToObj(&nmesh, out);

	cout << "--> Computing Tuette map" << endl;
	tuette_map(&nmesh);
	out = string(argv[2]) + "_tuette.obj";
	of.writeToObj(&nmesh, out);

	cout << "--> Computing harmonic map" << endl;
	// harmonic_map_gd(&nmesh);
	harmonic_map_conjugate_gd(&nmesh);

	out = string(argv[2]) + "_harmonic.obj";
	of.writeToObj(&nmesh, out);

	return 0;
}
