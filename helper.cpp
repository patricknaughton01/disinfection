
#include <vector>
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
