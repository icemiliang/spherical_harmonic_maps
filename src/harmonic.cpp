#include "harmonic.h"



void Harmonic::star_map(Solid *cmesh) {
	Point center(0, 0, 0);
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		center += vertex->point();
	}
	center /= cmesh->numVertices();
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
		v_n(vertex) = vertex->point();
	}
}

double Harmonic::compute_energy(Solid *cmesh, int type) {
	double energy = 0;
	for (SolidEdgeIterator eiter(cmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Solid::tVertex v1 = cmesh->edgeVertex1(edge);
		Solid::tVertex v2 = cmesh->edgeVertex2(edge);
		Point uv = v1->point() - v2->point();
		if (type == HARMONIC) {
			energy += e_k(edge) * uv.norm2();
		}
		else energy += uv.norm2();
	}
	return energy;
}

void Harmonic::compute_gradient(Solid *cmesh, int type) {
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex sv = *viter;
		Point gradient = Point(0, 0, 0);
		for (VertexVertexIterator vviter(sv); !vviter.end(); ++vviter) {
			Vertex * tv = *vviter;
			if (type == HARMONIC)
				gradient += (sv->point() - tv->point()) * e_k(cmesh->vertexEdge(sv, tv));
			else gradient += (sv->point() - tv->point());
		}
		v_dv(sv) = gradient;
		v_abdv(sv) = v_dv(sv) - v_n(sv) * (v_dv(sv) * v_n(sv));
	}
}

void Harmonic::conjugate_gradient(Solid *cmesh) {
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		v_beta(vertex) = v_abdv(vertex)*v_abdv(vertex) / (v_fabdv(vertex)*v_fabdv(vertex));
		v_s(vertex) = v_s(vertex)*v_beta(vertex) - v_abdv(vertex);
	}
	double alpha = 1e-6;
	double e0 = compute_energy(cmesh, HARMONIC);
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v_mp(v) = v->point();
	}
	while (alpha<1e-2) {
		for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
			Solid::tVertex v = *viter;
			v->point() = v->point() + v_s(v)*alpha;
			v->point() /= v->point().norm();
		}
		double nenergy = compute_energy(cmesh, HARMONIC);
		if (nenergy < e0) {
			for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
				Solid::tVertex v = *viter;
				v_mp(v) = v->point();
			}
		}
		alpha *= 2;
	}
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v_mp(v);
	}
}

void Harmonic::update_mesh(Solid *cmesh) {
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Solid::tVertex v = *viter;
		v->point() = v->point()- v_abdv(v)*DELTA_T;
		v->point() /= v->point().norm();
	}
}

void Harmonic::update_mass_center(Solid * cmesh) {
	for (SolidFaceIterator fiter(cmesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;

		Point p[3];
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex * v = *fviter;
			p[i++] = v->point();
		}

		Point n = (p[1] - p[0]) ^ (p[2] - p[0]);
		f_a(face) = n.norm() / 2.0;
	}
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		double area = 0;
		for (VertexFaceIterator vfiter(vertex); !vfiter.end(); ++vfiter) {
			Face* vf = *vfiter;
			area += f_a(vf);
		}
		v_a(vertex) = area / 3.0;
	}

	Point center(0, 0, 0);
	double mass = 0;
	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		center += vertex->point() * v_a(vertex);
		mass += v_a(vertex);
	}
	center /= mass;

	for (SolidVertexIterator viter(cmesh); !viter.end(); ++viter) {
		Vertex * vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
	}
}

void Harmonic::harmonic_map_conjugate_gd(Solid * pnmesh) {
	double te = compute_energy(pnmesh, HARMONIC);
	double oe = 1000;
	compute_gradient(pnmesh, HARMONIC);
	update_mesh(pnmesh);
	for (SolidVertexIterator viter(pnmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		v_fabdv(vertex) = v_abdv(vertex);
		v_s(vertex) = v_abdv(vertex);
	}
	double he = compute_energy(pnmesh, HARMONIC);
	int i = 0;
	while (oe - he > DELTA_E) {
		compute_gradient(pnmesh, HARMONIC);
		conjugate_gradient(pnmesh);
		oe = he;
		he = compute_energy(pnmesh, HARMONIC);
		if (i++ < 10 or i % 100 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::harmonic_map_gd(Solid * pnmesh) {
	double he = compute_energy(pnmesh, HARMONIC);
	double oe = 1000;
	int i = 0;
	while (oe - he > DELTA_E) {
		compute_gradient(pnmesh, HARMONIC);
		update_mesh(pnmesh);
		update_mass_center(pnmesh);
		oe = he;
		he = compute_energy(pnmesh, HARMONIC);
		if (i++ < 10 or i % 100 == 0) {
			std::cout << i << ", energy: " << he << ", diff: "<< oe - he << std::endl;	
		}
	}
}

void Harmonic::tuette_map(Solid * pnmesh) {
	double te = compute_energy(pnmesh, TUETTE);
	double oe = 1000;
	int i = 0;
	while (oe-te > DELTA_E) {
		compute_gradient(pnmesh, TUETTE);
		update_mesh(pnmesh);
		oe = te;
		te = compute_energy(pnmesh, TUETTE);
		if (i++ < 10 or i % 100 == 0) {
			std::cout << i << ", energy: " << te << ", diff: "<< oe - te << std::endl;	
		}
	}
}

void Harmonic::set_up_normal(Solid *nmesh) {
	for (SolidVertexIterator viter(nmesh); !viter.end(); ++viter) {
		Solid::tVertex vertex = *viter;
		vertex->trait() = new CVertexTrait;
	}
	for (SolidFaceIterator fiter(nmesh); !fiter.end(); ++fiter) {
		Solid::tFace face = *fiter;
		face->trait() = new CFaceTrait;
		Point fp[3]; 
		int i = 0;
		for (FaceVertexIterator fviter(face); !fviter.end(); ++fviter) {
			Vertex* fv = *fviter;
			fp[i++] = fv->point();
		}
		Point normal = (fp[1] - fp[0]) ^ (fp[2] - fp[0]);
		f_n(face) = normal / normal.norm();
	}
}

void Harmonic::set_up_kuv(Solid *nmesh) {
	for (SolidEdgeIterator eiter(nmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		edge->trait() = new CEdgeTrait;
		e_k(edge) = 1.0;
	}

	double sum = 0;
	for (SolidEdgeIterator eiter(nmesh); !eiter.end(); ++eiter) {
		Solid::tEdge edge = *eiter;
		Point p1 = nmesh->edgeVertex1(edge)->point();
		Point p2 = nmesh->edgeVertex2(edge)->point();
		Point p3 = edge->halfedge(0)->he_next()->target()->point();

		double alpha, beta;
		alpha = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		p3 = edge->halfedge(0)->he_sym()->he_next()->target()->point();
		beta = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2;
		e_k(edge) = alpha + beta;
		sum += e_k(edge);
	}
}

void Harmonic::set_up(Solid * nmesh) {
	set_up_normal(nmesh);
	set_up_kuv(nmesh);
}
