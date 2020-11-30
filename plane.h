#pragma once
#include <unordered_set>
#include <vector>
#include <memory>
#include <iostream>
#include "helper.h"

class Plane;
struct PointerHash;
struct DerefCompare;

typedef long plane_id;
typedef std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
	plane_set;
typedef std::pair<std::shared_ptr<Plane>, std::shared_ptr<Plane>> plane_pair;

class Plane{
	friend std::ostream& operator<<(std::ostream &out, const Plane &p);
private:
	std::vector<REAL> norm;
	std::vector<REAL> centroid;
	std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
		neighbors;
	REAL area;
	static plane_id gen_id;
	void init_id();
public:
	plane_id id;
	Plane();
	Plane(std::vector<REAL> &n, std::vector<REAL> &c, REAL a);
	Plane(std::vector<REAL> &n, std::vector<REAL> &c,
		std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>
		&neigh, REAL a);
	void merge(std::shared_ptr<Plane> other);
	void neighbor_union(std::shared_ptr<Plane> other);
	void add_neighbor(std::shared_ptr<Plane> other);
	REAL score(std::shared_ptr<Plane> other) const;
	const std::unordered_set<std::shared_ptr<Plane>, PointerHash, DerefCompare>&
		get_neighbors() const;
	bool operator<(const Plane &other) const;
	bool operator==(const Plane &other) const;
	bool operator!=(const Plane &other) const;
	size_t hash() const;
};

plane_pair get_pair(std::shared_ptr<Plane> a,
	std::shared_ptr<Plane> b);

struct PairPointerHash{
    size_t operator()(const std::shared_ptr<plane_pair> &p) const;
};

struct PairDerefCompare{
	bool operator()(const std::shared_ptr<plane_pair> &p1,
		const std::shared_ptr<plane_pair> &p2) const;
};
