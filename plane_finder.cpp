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
#include <Eigen/Dense>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

#define DEBUG

void PlaneFinder::simplify_planes(plane_set &out, REAL thresh)
{
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
	dedup_triangle_mesh(vertices, inds);
	std::vector<size_t> tmp;
	// Map vertex i (row) to the triangle indices it borders
	std::vector<std::vector<size_t>> vmap(vertices.size(), tmp);
	// Map triangle indices to (pointers to ) planes
	std::unordered_map<size_t, std::shared_ptr<Plane>> t_to_plane;

	for(auto iter = inds.begin(); iter != inds.end(); iter++){
		size_t tri_ind = iter - inds.begin();
		// Compute geometric properties of planes
		size_t ind0 = (*iter)[0], ind1 = (*iter)[1], ind2 = (*iter)[2];
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
		std::shared_ptr<Plane> p = std::make_shared<Plane>(vnorm, vcen, area);
		planes.insert(p);
		t_to_plane[tri_ind] = p;
		vmap[ind0].push_back(tri_ind);
		vmap[ind1].push_back(tri_ind);
		vmap[ind2].push_back(tri_ind);
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
	int defective = 0;
	for(auto it = planes.begin(); it != planes.end(); it++){
		if((*it)->get_neighbors().size() != 3){
			std::cout << "Found plane without three neighbors" << std::endl;
			std::cout << **it;
			std::cout << " has " << (*it)->get_neighbors().size() << " neighbors" << std::endl;
			defective++;
		}
	}
	std::cout << "Found " << defective << " defective triangles out of "
		<< planes.size() << std::endl;
}

void PlaneFinder::dedup_triangle_mesh(
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds)
{
	// Make a set of vertices (without duplicates)
	std::unordered_map<std::vector<REAL>, size_t, VectorHash> map_verts;
	std::unordered_map<size_t, size_t> ind_to_newind;
	// deduplicated vertices
	std::vector<std::vector<REAL>> dd_verts;
	for(auto it = vertices.begin(); it != vertices.end(); it++){
		size_t ind = it - vertices.begin();
		// N^2 algorithm for now, could speed up with kdtree
		bool repeat = false;
		size_t new_ind = 0;
		for(auto d_it = dd_verts.begin(); d_it != dd_verts.end(); d_it++){
			if(get_dist(*it, *d_it) < 1e-5){
				repeat = true;
				break;
			}
			new_ind++;
		}
		ind_to_newind[ind] = new_ind;
		if(!repeat){
			dd_verts.push_back(*it);
		}
		// if(map_verts.find(*it) == map_verts.end()){
		// 	map_verts[*it] = new_ind;
		// 	dd_verts.push_back(*it);
		// 	ind_to_newind[ind] = new_ind;
		// 	new_ind++;
		// }else{
		// 	std::cout << "------------------FOUND REPEAT VERTEX------------------" << std::endl;
		// 	ind_to_newind[ind] = map_verts[*it];
		// }
	}

	// Fixup the indices vector using the two maps:
	// old_index -> set_iterator -> new_index
	for(auto ind_it = inds.begin(); ind_it != inds.end(); ind_it++){
		for(auto it = ind_it->begin(); it != ind_it->end(); it++){
			*it = ind_to_newind[*it];
		}
	}
	// Update with the cleaned vertices
	vertices = dd_verts;
	std::cout << "Dedup Verts: " << std::endl;
	for(auto it = dd_verts.begin(); it != dd_verts.end(); it++){
		print(it->begin(), it->end());
	}
	std::cout << "Dedup Inds: " << std::endl;
	for(auto it = inds.begin(); it != inds.end(); it++){
		print(it->begin(), it->end());
	}
}

PlaneFinder::PlaneFinder(){}

PlaneFinder::PlaneFinder(plane_set &p):planes(p){}
