#pragma once
#include <vector>
#include "helper.h"

class SamplePoint{
private:
	std::vector<REAL> normal;
	std::vector<REAL> world_point;
public:
	SamplePoint(){};
	SamplePoint(std::vector<REAL> n, std::vector<REAL> w):normal(n),
		world_point(w){};
};
