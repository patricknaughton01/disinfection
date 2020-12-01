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

int main(){
	plane_set planes;
	REAL area = 1.0;
	std::vector<REAL> n1{0.1, 0, 1.0};
	std::vector<REAL> n2{0, 0, 1.0};
	std::vector<REAL> n3{1.0, 0, 0};
	std::vector<REAL> c1{0, 0, 0};
	std::vector<REAL> c2{1.0, 0, 0};
	std::vector<REAL> c3{1.0, 0, -1.0};
	auto p1 = std::make_shared<Plane>(n1, c1, area);
	auto p2 = std::make_shared<Plane>(n2, c2, area);
	auto p3 = std::make_shared<Plane>(n3, c3, area);
	p1->add_neighbor(p2);
	p2->add_neighbor(p1);
	p2->add_neighbor(p3);
	p3->add_neighbor(p2);
	planes.insert(p1);
	planes.insert(p2);
	planes.insert(p3);
	plane_set out;
	PlaneFinder pf(planes);
	print_p(planes.begin(), planes.end());
	pf.simplify_planes(out, 0.5);
	print_p(planes.begin(), planes.end());
	print_p(out.begin(), planes.end());
	return 0;
}

void PlaneFinder::simplify_planes(plane_set &out, REAL thresh)
{
	pq.clear();
	locs.clear();
	init_pq();
	auto min_iter = pq.begin();
	out = planes;

	#ifdef DEBUG
		std::cout << "Initial pq:" << std::endl;
		for(auto iter = pq.begin(); iter != pq.end(); iter++){
			std::cout << "(" << iter->first << "->(" << *(iter->second->first)
				<< "," << *(iter->second->second) << ")) ";
		}
		std::cout << std::endl;
	#endif // DEBUG

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
			std::shared_ptr<plane_pair> uPair = std::make_shared<plane_pair>(
				get_pair(nPlane, *iter)
			);
			pq_iter_t in_iter = pq.insert(
				std::make_pair(nPlane->score(*iter), uPair));
			locs[uPair] = in_iter;
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
			visited.insert(*iter);
			while(open.size()){
				std::shared_ptr<Plane> plane = open.top();
				open.pop();
				for(auto piter = plane->get_neighbors().begin();
					piter != plane->get_neighbors().end(); piter++)
				{
					if (visited.find(*piter) == visited.end()){
						visited.insert(*piter);
						open.push(*piter);
						std::shared_ptr<plane_pair> p =
							std::make_shared<plane_pair>(
								get_pair(plane, *piter));
						pq_iter_t in_iter = pq.insert(
							std::make_pair(plane->score(*piter), p));
						locs[p] = in_iter;
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
		if(**iter != *p2){
			auto p1pair = std::make_shared<plane_pair>(get_pair(p1, *iter));
			auto piter = locs[p1pair];
			locs.erase(p1pair);
			pq.erase(piter);
		}
	}
}

/*
 * Each vector in vertices has three elements (spatial coordinates), each
 * vector in inds also has three elements, the three vertices of that triangle
 */
void PlaneFinder::load_triangle_mesh(
	const std::vector<std::vector<REAL>> &vertices,
	const std::vector<std::vector<size_t>> &inds)
{
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
		// Map triangle indices to the number of times they border the given
		// vertex
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
}

PlaneFinder::PlaneFinder(){}

PlaneFinder::PlaneFinder(plane_set &p):planes(p){}
