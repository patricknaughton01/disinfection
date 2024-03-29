#pragma once

#include <unordered_set>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <KrisLibrary/meshing/TriMesh.h>
#include "helper.h"
#include "plane.h"
#include "heightmap.h"

class PlaneFinder{
private:
	// Use a multimap as a priority queue
	std::multimap<REAL, std::shared_ptr<plane_pair>> pq;
	typedef decltype(pq.begin()) pq_iter_t;
	std::unordered_map<std::shared_ptr<plane_pair>, pq_iter_t, PairPointerHash,
		PairDerefCompare> locs;
	plane_set planes;
	void clean_neighbors(std::shared_ptr<plane_pair> ppair);
	void clean_neighbors(std::shared_ptr<Plane> p1, std::shared_ptr<Plane> p2);
	std::vector<std::vector<REAL>> i_vertices;
	std::vector<std::vector<size_t>> i_inds;
	std::vector<std::vector<REAL>> i_normals;
	std::vector<std::vector<REAL>> i_centroids;
	Meshing::TriMesh inflated_mesh;
	std::unordered_map<std::shared_ptr<Plane>, std::shared_ptr<Heightmap>>
		heightmaps;
	friend std::pair<std::vector<std::vector<std::vector<std::vector<REAL>>>>,
		std::vector<std::vector<std::vector<std::vector<REAL>>>>> get_heightmaps(
		PlaneFinder &pf);
	friend std::vector<std::vector<REAL>> get_heightmap_metadata(
		PlaneFinder &pf);
public:
	PlaneFinder();
	PlaneFinder(plane_set &p);
	void simplify_planes(plane_set &out, REAL thresh);
	void load_heightmaps(REAL spacing, REAL border);
	void load_triangle_mesh(std::vector<std::vector<REAL>> &vertices,
		std::vector<std::vector<size_t>> &inds);
	void init_pq();
};
