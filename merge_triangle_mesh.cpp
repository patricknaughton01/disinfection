#include <iostream>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "merge_triangle_mesh.h"

std::vector<std::vector<std::vector<REAL>>> merge_triangle_mesh(
	const std::vector<std::vector<REAL>> &vertices,
	const std::vector<std::vector<size_t>> &inds,
	REAL threshold)
{
	std::cout << "Started" << std::endl;
	PlaneFinder pf;
	pf.load_triangle_mesh(vertices, inds);
	std::cout << "Loaded TM" << std::endl;
	plane_set out;
	pf.simplify_planes(out, threshold);
	std::cout << "Simplified planes" << std::endl;
	// TODO: Can I return this by reference to python? Or via a shared_ptr?
	std::vector<std::vector<std::vector<REAL>>> ret;
	for(auto iter = out.begin(); iter != out.end(); iter++){
		std::vector<std::vector<REAL>> p;
		p.push_back((*iter)->get_norm());
		p.push_back((*iter)->get_centroid());
		ret.push_back(p);
	}
	std::cout << "Reconstructed planes" << std::endl;
	return ret;
}
