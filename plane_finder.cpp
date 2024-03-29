#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <string>
#include <sstream>
#include <memory>
#include <stack>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <KrisLibrary/meshing/TriMesh.h>
#include <KrisLibrary/math3d/primitives.h>
#include <KrisLibrary/math3d/Ray3D.h>
#include <KrisLibrary/utils/IntTriple.h>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "sample_point.h"
#include "heightmap.h"

void PlaneFinder::load_heightmaps(REAL spacing,
	REAL border)
{
	typedef std::vector<REAL> vr;
	for(auto it = planes.begin(); it != planes.end(); it++){
		std::unordered_set<plane_id> triangles((*it)->get_triangles().begin(),
			(*it)->get_triangles().end());
		std::unordered_set<size_t> tested_v_inds;
		std::vector<vr> p_verts;
		vr x_axis;
		REAL max_len = 0;
		// Find the farthest point from the origin to make the x-axis
		for(auto t_iter = (*it)->get_triangles().begin();
			t_iter != (*it)->get_triangles().end(); t_iter++)
		{
			for(auto v_it = i_inds[*t_iter].begin();
				v_it != i_inds[*t_iter].end(); v_it++)
			{
				if(tested_v_inds.find(*v_it) == tested_v_inds.end()){
					vr target_v = i_vertices[*v_it];
					vr d = target_v - (*it)->get_centroid();
					vr n = (*it)->get_norm();
					vr perp_d = d - (d * n) * n;
					p_verts.push_back(perp_d);
					REAL dist = get_norm(perp_d);
					if(dist > max_len){
						max_len = dist;
						normalize_vector(perp_d);
						x_axis = perp_d;
					}
					tested_v_inds.insert(*v_it);
				}
			}
		}
		// Construct the right-handed y-axis
		vr y_axis = cross((*it)->get_norm(), x_axis);
		normalize_vector(y_axis);
		// Find bounding box of the points
		REAL max_x = 0, min_x = 0, max_y = 0, min_y = 0;
		for(size_t i = 0; i < p_verts.size(); i++){
			REAL x_val = p_verts[i] * x_axis;
			if(x_val > max_x){
				max_x = x_val;
			}
			if(x_val < min_x){
				min_x = x_val;
			}
			REAL y_val = p_verts[i] * y_axis;
			if(y_val > max_y){
				max_y = y_val;
			}
			if(y_val < min_y){
				min_y = y_val;
			}
		}
		int num_x = ((int)std::ceil((max_x - min_x + 2 * border)
			/ spacing));
		int num_y = ((int)std::ceil((max_y - min_y + 2 * border)
			/ spacing));
		int ao_x = ((int)((-min_x + border) / spacing));
		int ao_y = ((int)((max_y + border) / spacing));
		// std::cout << "ao x: " << ao_x << std::endl;
		// std::cout << "ao y: " << ao_y << std::endl;
		std::vector<std::vector<SamplePoint>> sample_points;
		std::vector<vr> axes{x_axis, y_axis, (*it)->get_norm()};
		std::shared_ptr<Heightmap> hm = std::make_shared<Heightmap>(
			(*it)->get_centroid(), axes, spacing, border, min_x, max_x, min_y,
			max_y, std::pair<size_t, size_t>(ao_x, ao_y));
		for(int i = 0; i < num_y; i++){
			std::vector<SamplePoint> row;
			for(int j = 0; j < num_x; j++){
				REAL x_val = (j - ao_x) * spacing;
				REAL y_val = (i - ao_y) * spacing;
				REAL z_val = 0;
				// std::cout << "Local point: " << x_val << " " << y_val << std::endl;
				Math3D::Vector3 plane_point(x_val, y_val, z_val);
				Math3D::Vector3 world_point;
				hm->trans.mulPoint(plane_point, world_point);
				Math3D::Vector3 direction(hm->world_axes[2].data());
				Math3D::Ray3D r1;
				Math3D::Ray3D r2;
				r1.source = world_point;
				r2.source = world_point;
				r1.direction = direction;
				r2.direction = -direction;
				Math3D::Vector3 hit1;
				Math3D::Vector3 hit2;
				int res1 = inflated_mesh.RayCast(r1, hit1);
				int res2 = inflated_mesh.RayCast(r2, hit2);
				bool valid_hit = false;
				REAL d1 = std::numeric_limits<double>::infinity();
				REAL d2 = d1;
				if(res1 > -1 && (i_normals[res1] * hm->world_axes[2]) > 0
					&& triangles.find(res1) != triangles.end())
				{
					valid_hit = true;
					d1 = (hit1 - world_point).norm();
				}
				if(res2 > -1 && (i_normals[res2] * hm->world_axes[2]) > 0
					&& triangles.find(res2) != triangles.end())
				{
					valid_hit = true;
					d2 = (hit2 - world_point).norm();
				}
				SamplePoint pt;
				if(valid_hit){
					if(d1 < d2){
						pt = SamplePoint(
							Math3D::Vector3(i_normals[res1].data()),
							hit1, true);
					}else{
						pt = SamplePoint(
							Math3D::Vector3(i_normals[res2].data()),
							hit2, true);
					}
					// std::cout << "Hitpoint: ";
					// std::cout << pt.world_point.x << " " << pt.world_point.y << " " << pt.world_point.z << std::endl;
				}
				row.push_back(pt);
			}
			sample_points.push_back(row);
		}
		hm->sample_points = sample_points;
		heightmaps[*it] = hm;
	}
}

void PlaneFinder::simplify_planes(plane_set &out, REAL thresh)
{
	// TODO: kind of sketchy that out gets shared_ptrs to planes inside the
	// internal (private) planes field, meaning someone else could modify the
	// planes in planes without permission. Perhaps could be fixed by making
	// out a set of pointers to const Planes so that it can't be used to modify
	// them? Then we have to definitely change the setting of planes to out.
	pq.clear();
	locs.clear();
	init_pq();
	if(!pq.size()){
		return;
	}
	auto min_iter = pq.begin();
	out = planes;

	// #ifdef DEBUG
	// 	std::cout << "Initial pq:" << std::endl;
	// 	for(auto iter = pq.begin(); iter != pq.end(); iter++){
	// 		std::cout << "(" << iter->first << "->(" << *(iter->second->first)
	// 			<< "," << *(iter->second->second) << ")) ";
	// 	}
	// 	std::cout << std::endl;
	// #endif // DEBUG

	// Merge until the best candidate is worse than the given threshold
	while(min_iter->first <= thresh){
		std::shared_ptr<Plane> p1 = min_iter->second->first;
		std::shared_ptr<Plane> p2 = min_iter->second->second;
		clean_neighbors(min_iter->second);
		std::shared_ptr<Plane> bigger_neighborhood(p1);
		std::shared_ptr<Plane> smaller_neighborhood(p2);
		if(p2->get_neighbors().size() > p1->get_neighbors().size()){
			bigger_neighborhood = p2;
			smaller_neighborhood = p1;
		}
		std::shared_ptr<Plane> nPlane(bigger_neighborhood);
		nPlane->merge(smaller_neighborhood);
		for(auto iter = nPlane->get_neighbors().begin();
			iter != nPlane->get_neighbors().end(); iter++)
		{
			if(auto spt = iter->lock()){
				std::shared_ptr<plane_pair> uPair = std::make_shared<plane_pair>(
					get_pair(nPlane, spt)
				);
				spt->add_neighbor(nPlane);
				pq_iter_t in_iter = pq.insert(
					std::make_pair(nPlane->score(spt), uPair));
				locs[uPair] = in_iter;
			}
		}
		out.erase(p1);
		out.erase(p2);
		out.insert(nPlane);
		locs.erase(min_iter->second);
		pq.erase(min_iter);
		if(!pq.size()){
			// We've merged all pairs, terminate
			break;
		}
		min_iter = pq.begin();
	}
	planes = out;
	#ifdef DEBUG
		std::cout << "Final pq:" << std::endl;
		for(auto iter = pq.begin(); iter != pq.end(); iter++){
			std::cout << "(" << iter->first << "->(" << *(iter->second->first)
				<< "," << *(iter->second->second) << ")) ";
		}
		std::cout << std::endl;
		std::cout << "Final planes:" << std::endl;
		for(auto iter = out.begin(); iter != out.end(); iter++){
			std::cout << **iter << std::endl;
		}
		std::cout << std::endl;
	#endif // DEBUG
}

void PlaneFinder::init_pq(){
	std::stack<std::shared_ptr<Plane>> open;
	plane_set visited;
	for(auto iter = planes.begin(); iter != planes.end(); iter++){
		if(visited.find(*iter) == visited.end()){
			// Haven't visited this plane yet, do a DFS from here to add
			// neighboring planes to the pq
			open.push(*iter);
			while(open.size()){
				std::shared_ptr<Plane> plane = open.top();
				visited.insert(plane);
				open.pop();
				for(auto piter = plane->get_neighbors().begin();
					piter != plane->get_neighbors().end(); piter++)
				{
					if(auto spt = piter->lock()){
						if (visited.find(spt) == visited.end()){
							open.push(spt);
							std::shared_ptr<plane_pair> p =
								std::make_shared<plane_pair>(
									get_pair(plane, spt));
							pq_iter_t in_iter = pq.insert(
								std::make_pair(plane->score(spt), p));
							locs[p] = in_iter;
						}
					}
				}
			}
		}
	}
}

void PlaneFinder::clean_neighbors(std::shared_ptr<plane_pair> ppair){
	clean_neighbors(ppair->first, ppair->second);
	clean_neighbors(ppair->second, ppair->first);
}

void PlaneFinder::clean_neighbors(std::shared_ptr<Plane> p1, std::shared_ptr<Plane> p2){
	for(auto iter = p1->get_neighbors().begin();
		iter != p1->get_neighbors().end(); iter++)
	{
		if(auto spt = iter->lock()){
			if(*spt != *p2){
				spt->remove_neighbor(p1);
				auto p1pair = std::make_shared<plane_pair>(get_pair(p1, spt));
				auto piter = locs[p1pair];
				locs.erase(p1pair);
				pq.erase(piter);
			}
		}
	}
}

/*
 * Each vector in vertices has three elements (spatial coordinates), each
 * vector in inds also has three elements, the three vertices of that triangle.
 * vertices and inds not const because they may be cleaned (remove duplicate
 * vertices). To make them constant, could make copies, but that would cause
 * some performance hit.
 */
void PlaneFinder::load_triangle_mesh(
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds)
{
	i_vertices = vertices;
	i_inds = inds;
	std::vector<size_t> tmp;
	// Map vertex i (row) to the triangle indices it borders
	std::vector<std::vector<size_t>> vmap(vertices.size(), tmp);
	// Map triangle indices to (pointers to ) planes
	std::unordered_map<size_t, std::shared_ptr<Plane>> t_to_plane;
	std::vector<Math3D::Vector3> krislibrary_verts;
	std::vector<IntTriple> krislibrary_inds;

	for(auto iter = inds.begin(); iter != inds.end(); iter++){
		size_t tri_ind = iter - inds.begin();
		// Compute geometric properties of planes
		size_t ind0 = (*iter)[0], ind1 = (*iter)[1], ind2 = (*iter)[2];
		krislibrary_inds.push_back(IntTriple(ind0, ind1, ind2));
		// IMPROVE: Could remove Eigen dependency by refactoring just this part.
		Eigen::Vector3d a(vertices[ind0].data());
		Eigen::Vector3d b(vertices[ind1].data());
		Eigen::Vector3d c(vertices[ind2].data());
		Eigen::Vector3d cen = (a + b + c) / 3.0;
		Eigen::Vector3d norm = (b - a).cross(c - a);
		REAL area = std::sqrt(norm.transpose() * norm) / 2;
		norm /= 2.0 * area;
		std::vector<REAL> vnorm(norm.data(), norm.data() + norm.size());
		std::vector<REAL> vcen(cen.data(), cen.data() + cen.size());
		std::list<plane_id> tri_list(1, tri_ind);
		std::shared_ptr<Plane> p = std::make_shared<Plane>(vnorm, vcen,
			tri_list, area);
		planes.insert(p);
		t_to_plane[tri_ind] = p;
		vmap[ind0].push_back(tri_ind);
		vmap[ind1].push_back(tri_ind);
		vmap[ind2].push_back(tri_ind);
		i_normals.push_back(vnorm);
		i_centroids.push_back(vcen);
	}

	for(size_t i = 0; i < inds.size(); i++){
		// Map triangle indices to the number of vertices of triangle i they
		// border
		std::unordered_map<size_t, int> count;
		for(size_t j = 0; j < inds[i].size(); j++){
			for(size_t k = 0; k < vmap[inds[i][j]].size(); k++){
				size_t tind = vmap[inds[i][j]][k];
				if(count.find(tind) == count.end()){
					count[tind] = 1;
				}else{
					count[tind]++;
					if(count[tind] == 2 && tind != i){
						// This triangle borders triangle i because tind != i
						// and it shares (at least) two vertices with triangle i
						t_to_plane[i]->add_neighbor(t_to_plane[tind]);
					}
				}
			}
		}
	}
	// Construct list of vertices for inflated_mesh
	for(auto it = vertices.begin(); it != vertices.end(); it++){
		krislibrary_verts.push_back(Math3D::Vector3(it->data()));
	}
	inflated_mesh.verts = krislibrary_verts;
	inflated_mesh.tris = krislibrary_inds;
}

PlaneFinder::PlaneFinder(){}

PlaneFinder::PlaneFinder(plane_set &p):planes(p){}
