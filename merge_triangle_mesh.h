#include <vector>
#include <utility>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

std::pair<std::vector<std::vector<std::vector<REAL>>>,
	std::vector<std::vector<plane_id>>> merge_triangle_mesh(
	std::vector<std::vector<REAL>> &vertices,
	std::vector<std::vector<size_t>> &inds,
	REAL threshold
);
