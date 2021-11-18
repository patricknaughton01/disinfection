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
	REAL min_x, max_x, min_y, max_y;
	// x, y; i.e., col, row
	std::pair<size_t, size_t> arr_origin;
	Math3D::RigidTransform trans;
	std::vector<std::vector<SamplePoint>> sample_points;
	std::pair<std::vector<REAL>, std::vector<REAL>> interpolate(REAL x, REAL y);
	Heightmap():spacing(0), border(0), min_x(0), max_x(0), min_y(0), max_y(0){}
	Heightmap(std::vector<REAL> o, std::vector<std::vector<REAL>> w,
		REAL s, REAL b, REAL n_x, REAL x_x, REAL n_y, REAL x_y,
		std::pair<size_t, size_t> ao):
		origin(o), world_axes(w), spacing(s), border(b), min_x(n_x), max_x(x_x),
		min_y(n_y), max_y(x_y), arr_origin(ao),
		trans(RigidTransform(Vector3(w[0].data()), Vector3(w[1].data()),
			Vector3(w[2].data()), Vector3(o.data()))){}
};
