#pragma once

#include <unordered_set>
#include <vector>
#include <memory>

typedef long plane_id;

#define REAL double

void avg_vector(std::vector<REAL> &a, const std::vector<REAL> &b,
	REAL weight_a=1.0, REAL weight_b=1.0);

struct PointerHash{
	template <typename T>
    size_t operator()(const std::shared_ptr<T> &a) const;
};

struct DerefCompare {
	template <typename T>
	bool operator() (std::shared_ptr<T> const &a,
		std::shared_ptr<T> const &b) const;
};

class Plane{
private:
	std::vector<REAL> norm;
	std::vector<REAL> centroid;
	std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
		neighbors;
	REAL area;
	static plane_id gen_id;
	plane_id id;
	void init_id();
public:
	Plane();
	Plane(std::vector<REAL> &n, std::vector<REAL> &c, REAL a);
	Plane(std::vector<REAL> &n, std::vector<REAL> &c,
		std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
		&neigh, REAL a);
	void merge(const Plane &other);
	void neighbor_union(const Plane &other);
	void add_neighbor(std::shared_ptr<Plane> other);
	REAL score(const Plane &other) const;
	const std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>&
		get_neighbors() const;
	bool operator<(const Plane &other) const;
	bool operator==(const Plane &other) const;
	size_t hash() const;
};

typedef std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
	plane_set;

void simplify_planes(plane_set planes, plane_set out);
