#pragma once
#include <vector>
#include <utility>
#include <KrisLibrary/math3d/primitives.h>
#include "sample_point.h"
#include "helper.h"

using Math3D::RigidTransform;
using Math3D::Vector3;

class Heightmap{
private:
public:
	std::vector<REAL> origin;
	std::vector<std::vector<REAL>> world_axes;
	REAL spacing;
	REAL border;
	// x, y; i.e., col, row
	std::pair<size_t, size_t> arr_origin;
	Math3D::RigidTransform trans;
	std::vector<std::vector<SamplePoint>> sample_points;
	Heightmap():spacing(0), border(0){}
	Heightmap(std::vector<REAL> o, std::vector<std::vector<REAL>> w,
		REAL s, REAL b, std::pair<size_t, size_t> ao):
		origin(o), world_axes(w), spacing(s), border(b), arr_origin(ao),
		trans(RigidTransform(Vector3(w[0].data()), Vector3(w[1].data()),
			Vector3(w[2].data()), Vector3(o.data()))){}
};
