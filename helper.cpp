#include <vector>
#include <cmath>
#include <sstream>
#include <exception>
#include "helper.h"

/*
 * Average vectors a and b into a, destroying the previous contents of a.
 * Throws length error if vectors are of different sizes. weight_a and
 * weight_b have default values of 1.0.
 */
void avg_vector(std::vector<REAL> &a, const std::vector<REAL> &b,
	REAL weight_a, REAL weight_b)
{
	assert_len_equal(a, b);
	REAL w_sum = weight_a + weight_b;
	for(size_t i = 0; i < a.size(); i++){
		a[i] = (weight_a * a[i] + weight_b * b[i]) / w_sum;
	}
}

REAL get_norm(const std::vector<REAL> &a){
	REAL sum = 0;
	for(auto iter = a.begin(); iter != a.end(); iter++){
		sum += *iter * *iter;
	}
	return std::sqrt(sum);
}

REAL get_dist(const std::vector<REAL> &a, const std::vector<REAL> &b){
	REAL sum = 0;
	if(a.size() != b.size()){ return -1; }
	for(size_t i = 0; i < a.size(); i++){
		REAL d = a[i] - b[i];
		sum += d*d;
	}
	return std::sqrt(sum);
}

void divide_vector(std::vector<REAL> &a, REAL d){
	if(d != 0){
		for(auto iter = a.begin(); iter != a.end(); iter++){
			*iter /= d;
		}
	}
}

void normalize_vector(std::vector<REAL> &a){
	REAL n = get_norm(a);
	divide_vector(a, n);
}

std::vector<REAL> cross(const std::vector<REAL> &a, const std::vector<REAL> &b){
	assert_len_equal(a, b);
	if(a.size() != 3){
		std::stringstream ss;
		ss << "Error in cross, length of a does "
			<< "not equal 3 (a.size = " << a.size() << ")";
		throw std::length_error(ss.str());
	}
	std::vector<REAL> r;
	r.push_back(a[1]*b[2] - a[2]*b[1]);
	r.push_back(a[2]*b[0] - a[0]*b[2]);
	r.push_back(a[0]*b[1] - a[1]*b[0]);
	return r;
}
