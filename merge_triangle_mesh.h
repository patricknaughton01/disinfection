#include <vector>
#include <utility>
#include <KrisLibrary/math3d/primitives.h>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

PlaneFinder get_empty_plane_finder();

std::pair<std::vector<std::vector<std::vector<REAL>>>,
	std::vector<std::vector<plane_id>>> merge_triangle_mesh(
	PlaneFinder &pf,
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds,
	REAL threshold
);

std::pair<std::vector<std::vector<std::vector<std::vector<REAL>>>>,
	std::vector<std::vector<std::vector<std::vector<REAL>>>>> get_heightmaps(
	PlaneFinder &pf, REAL spacing, REAL border);

std::vector<REAL> to_std(Math3D::Vector3 v);

std::pair<std::vector<std::vector<REAL>>, std::vector<std::vector<size_t>>>
	dedup_triangle_mesh(
		std::vector<std::vector<REAL>> &vertices,
		std::vector<std::vector<size_t>> &inds);
