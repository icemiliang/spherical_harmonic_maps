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
	cout << "--> Setting up kuv" << endl;
	shm.set_up();

	cout << "--> Star map" << endl;
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
	// vector<double> errors = shm.check_conformality();

	

	return 0;
}
