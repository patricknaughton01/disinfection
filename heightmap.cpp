#include <vector>
#include <utility>
#include "heightmap.h"

std::pair<std::vector<REAL>, std::vector<REAL>>
	Heightmap::interpolate(REAL x, REAL y)
{
	int x_ind = (int)(x / spacing);
	int y_ind = -(int)(y / spacing);
	x_ind += arr_origin.first;
	y_ind += arr_origin.second;
}
