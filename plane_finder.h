#pragma once

#include <unordered_set>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include "helper.h"
#include "plane.h"

class PlaneFinder{
private:
	// Use a multimap as a priority queue
	std::multimap<REAL, std::shared_ptr<plane_pair>> pq;
	typedef decltype(pq.begin()) pq_iter_t;
	std::unordered_map<std::shared_ptr<plane_pair>, pq_iter_t, PairPointerHash,
		PairDerefCompare> locs;
	void clean_neighbors(std::shared_ptr<plane_pair> ppair);
	void clean_neighbors(std::shared_ptr<Plane> p1, std::shared_ptr<Plane> p2);
	void dedup_triangle_mesh(std::vector<std::vector<REAL>> &vertices,
		std::vector<std::vector<size_t>> &inds);
public:
	plane_set planes;
	PlaneFinder();
	PlaneFinder(plane_set &p);
	void simplify_planes(plane_set &out, REAL thresh);
	void load_triangle_mesh(std::vector<std::vector<REAL>> &vertices,
		std::vector<std::vector<size_t>> &inds);
	void init_pq();
};
