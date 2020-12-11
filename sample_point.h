#pragma once
#include <vector>
#include <KrisLibrary/math3d/primitives.h>
#include "helper.h"

class SamplePoint{
private:
public:
	Math3D::Vector3 normal;
	Math3D::Vector3 world_point;
	bool valid;
	SamplePoint():valid(false){}
	SamplePoint(Math3D::Vector3 n, Math3D::Vector3 w, bool v):normal(n),
		world_point(w), valid(v){}
};
