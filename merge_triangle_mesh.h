#include <vector>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"

std::vector<std::vector<std::vector<REAL>>> merge_triangle_mesh(
	const std::vector<std::vector<REAL>> &vertices,
	const std::vector<std::vector<size_t>> &inds,
	REAL threshold
);
