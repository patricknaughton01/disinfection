#include <iostream>
#include <utility>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "merge_triangle_mesh.h"

std::pair<std::vector<std::vector<std::vector<REAL>>>,
	std::vector<std::vector<plane_id>>>
	merge_triangle_mesh(
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds,
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
	std::vector<std::vector<plane_id>> ret_tris;
	for(auto iter = out.begin(); iter != out.end(); iter++){
		std::vector<std::vector<REAL>> p;
		p.push_back((*iter)->get_norm());
		p.push_back((*iter)->get_centroid());
		ret.push_back(p);
		std::vector<plane_id> tmp;
		auto tris = (*iter)->get_triangles();
		for(auto t_it = tris.begin(); t_it != tris.end(); t_it++){
			tmp.push_back(*t_it);
		}
		ret_tris.push_back(tmp);
	}
	std::cout << "Reconstructed planes" << std::endl;
	return std::pair<std::vector<std::vector<std::vector<REAL>>>,
		std::vector<std::vector<plane_id>>>(ret, ret_tris);
}
