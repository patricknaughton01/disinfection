#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

std::vector<std::vector<std::vector<REAL>>> merge_triangle_mesh(
	const std::vector<std::vector<REAL>> &vertices,
	const std::vector<std::vector<size_t>> &inds,
	REAL threshold)
{
	PlaneFinder pf;
	pf.load_triangle_mesh(vertices, inds);
	plane_set out;
	pf.simplify_planes(out, threshold);
	// TODO: Can I return this by reference to python? Or via a shared_ptr?
	std::vector<std::vector<std::vector<REAL>>> ret;
	for(auto iter = out.begin(); iter != out.end(); iter++){
		std::vector<std::vector<REAL>> p;
		p.push_back((*iter)->get_norm());
		p.push_back((*iter)->get_centroid());
		ret.push_back(p);
	}
	return ret;
}
