#include <iostream>
#include <utility>
#include <KrisLibrary/geometry/KDTree.h>
#include <KrisLibrary/math/VectorTemplate.h>
#include <KrisLibrary/math3d/primitives.h>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "heightmap.h"
#include "sample_point.h"
#include "merge_triangle_mesh.h"

PlaneFinder get_empty_plane_finder(){
	return PlaneFinder();
}

std::pair<std::vector<std::vector<std::vector<REAL>>>,
	std::vector<std::vector<plane_id>>>
	merge_triangle_mesh(
	PlaneFinder pf,
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds,
	REAL threshold)
{
	std::cout << "Started" << std::endl;
	pf.load_triangle_mesh(vertices, inds);
	std::cout << "Loaded TM" << std::endl;
	plane_set out;
	pf.simplify_planes(out, threshold);
	std::cout << "Simplified planes" << std::endl;
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

std::pair<std::vector<std::vector<std::vector<Math3D::Vector3>>>,
	std::vector<std::vector<std::vector<Math3D::Vector3>>>> get_heightmaps(
	PlaneFinder pf, REAL spacing, REAL border)
{
	// Return a pair containing the normal information for each plane and the
	// true world point for each plane
	std::vector<std::vector<std::vector<Math3D::Vector3>>> norms;
	std::vector<std::vector<std::vector<Math3D::Vector3>>> w_pts;
	pf.load_heightmaps(spacing, border);
	for(auto it = pf.planes.begin(); it != pf.planes.end(); it++){
		auto hm = pf.heightmaps[*it];
		std::vector<std::vector<Math3D::Vector3>> p_norms;
		std::vector<std::vector<Math3D::Vector3>> p_w_pts;
		for(size_t i = 0; i < hm->sample_points.size(); i++){
			std::vector<Math3D::Vector3> row_n;
			std::vector<Math3D::Vector3> row_w;
			for(size_t j = 0; j < hm->sample_points[i].size(); j++){
				SamplePoint sp = hm->sample_points[i][j];
				row_n.push_back(sp.normal);
				row_w.push_back(sp.world_point);
			}
			p_norms.push_back(row_n);
			p_w_pts.push_back(row_w);
		}
		norms.push_back(p_norms);
		w_pts.push_back(p_w_pts);
	}
	return std::pair<std::vector<std::vector<std::vector<Math3D::Vector3>>>,
		std::vector<std::vector<std::vector<Math3D::Vector3>>>>(norms, w_pts);
}

std::pair<std::vector<std::vector<REAL>>, std::vector<std::vector<size_t>>>
dedup_triangle_mesh(
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds)
{
	// Make a set of vertices (without duplicates)
	std::unordered_map<std::vector<REAL>, size_t, VectorHash> map_verts;
	std::unordered_map<size_t, size_t> ind_to_newind;
	// deduplicated vertices
	std::vector<std::vector<REAL>> dd_verts;
	Geometry::KDTree lookup_tree;
	for(auto it = vertices.begin(); it != vertices.end(); it++){
		size_t ind = it - vertices.begin();
		Math::VectorTemplate<REAL> test_vec(*it);
		REAL dist;
		bool repeat = false;
		size_t new_ind = dd_verts.size();
		if(lookup_tree.TreeSize() > 0){
			size_t lookup_id = lookup_tree.ClosestPoint(test_vec, dist);
			if(dist < 1e-5){
				repeat = true;
				new_ind = ind_to_newind[lookup_id];
			}
		}
		ind_to_newind[ind] = new_ind;
		lookup_tree.Insert(test_vec, ind);
		if(!repeat){
			dd_verts.push_back(*it);
		}
	}

	// Fixup the indices vector using the two maps:
	// old_index -> set_iterator -> new_index
	for(auto ind_it = inds.begin(); ind_it != inds.end(); ind_it++){
		for(auto it = ind_it->begin(); it != ind_it->end(); it++){
			*it = ind_to_newind[*it];
		}
	}
	std::cout << "Dedup Verts: " << std::endl;
	for(auto it = dd_verts.begin(); it != dd_verts.end(); it++){
		print(it->begin(), it->end());
	}
	std::cout << "Dedup Inds: " << std::endl;
	for(auto it = inds.begin(); it != inds.end(); it++){
		print(it->begin(), it->end());
	}
	return std::pair<std::vector<std::vector<REAL>>,
		std::vector<std::vector<size_t>>>(dd_verts, inds);
}
