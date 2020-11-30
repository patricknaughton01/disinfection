#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <stack>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

int main(){
	plane_set planes;
	REAL area = 1.0;
	std::vector<REAL> n1{0, 0, 1.0};
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
	print_p(p1->get_neighbors().begin(), p1->get_neighbors().end());
	print_p(p2->get_neighbors().begin(), p2->get_neighbors().end());
	p1->merge(p2);
	print_p(p1->get_neighbors().begin(), p1->get_neighbors().end());
	plane_set out;
	PlaneFinder pf(planes);
	//pf.simplify_planes(out, 1.0);
	return 0;
}

void PlaneFinder::simplify_planes(plane_set &out, REAL thresh)
{
	pq.clear();
	locs.clear();
	init_pq();
	auto min_iter = pq.begin();
	while(min_iter->first <= thresh){
		std::shared_ptr<Plane> p1 = min_iter->second->first;
		std::shared_ptr<Plane> p2 = min_iter->second->second;
		for(auto iter = p1->get_neighbors().begin();
			iter != p1->get_neighbors().end(); iter++)
		{

		}
		pq.erase(min_iter);
		min_iter = pq.begin();
	}
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
				open.pop();
				for(auto piter = plane->get_neighbors().begin();
					piter != plane->get_neighbors().end(); piter++)
				{
					if (visited.find(*piter) == visited.end()){
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

PlaneFinder::PlaneFinder(){}

PlaneFinder::PlaneFinder(plane_set &p):planes(p){}
