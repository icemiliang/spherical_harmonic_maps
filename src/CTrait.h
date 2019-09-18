#ifndef _CTTAIT_H_
#define _CTTAIT_H_

#include "Trait.h"
#include <string>
#include <iterator>
#include "Point.h"

namespace MeshLib {
#define f_n(f) trait<CFaceTrait, Face>(f).fnormal()
#define f_a(f) trait<CFaceTrait, Face>(f).area()
#define e_k(e) trait<CEdgeTrait, Edge>(e).kuv()
#define e_l(e) trait<CEdgeTrait, Edge>(e).length()
#define v_n(v) trait<CVertexTrait, Vertex>(v).normal()
#define v_a(v) trait<CVertexTrait, Vertex>(v).area()
#define v_s(v) trait<CVertexTrait, Vertex>(v).s()
#define v_mp(v) trait<CVertexTrait, Vertex>(v).min_p()
#define v_dv(v) trait<CVertexTrait, Vertex>(v).dv()
#define v_abdv(v) trait<CVertexTrait, Vertex>(v).abdv()
#define v_beta(v) trait<CVertexTrait, Vertex>(v).beta()
#define v_fabdv(v) trait<CVertexTrait, Vertex>(v).fabdv()

class CEdgeTrait : public Trait {
public:
	CEdgeTrait(){};
	~CEdgeTrait() {};
	double & kuv() { return c_kuv; }
	double & length() { return c_length; }

private:
	double c_kuv;
	double c_length;
};

class CVertexTrait : public Trait {
public:
	CVertexTrait() {};
	~CVertexTrait() {};
	Point & normal() { return c_normal; }
	Point & dv() { return c_dv; }
	Point & abdv() { return c_abdv; }
	double & area() { return c_area; }
	double & beta() { return c_beta; }
	Point & fabdv() { return c_fabdv; }
	Point & s() { return c_s; }
	Point & min_p() { return c_minP; }

private:
	Point c_normal;
	Point c_dv;
	Point c_abdv;
	double c_area;
	double c_beta;
	Point c_fabdv;
	Point c_s;
	Point c_minP;
};

class CFaceTrait : public Trait {
public:
	CFaceTrait() {};
	~CFaceTrait() {};
	Point & fnormal() { return f_normal; }
	double & area() { return f_area; }
	
private:
	Point f_normal;
	double f_area;
};
}

#endif
