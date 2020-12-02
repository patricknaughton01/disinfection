#include <vector>
#include <iostream>
#include "plane.h"
#include "helper.h"

plane_id Plane::gen_id = 0;

Plane::Plane():area(0){
	init_id();
}

Plane::Plane(std::vector<REAL> &n, std::vector<REAL> &c, REAL a):norm(n),
	centroid(c), area(a){
	init_id();
	normalize_vector(norm);
}

Plane::Plane(std::vector<REAL> &n, std::vector<REAL> &c,
	plane_set &neigh, REAL a):norm(n), centroid(c), neighbors(neigh),
	area(a)
{
	init_id();
	normalize_vector(norm);
}

void Plane::init_id(){
	id = gen_id;
	gen_id++;
}

void Plane::merge(std::shared_ptr<Plane> other){
	avg_vector(norm, other->norm, area, other->area);
	normalize_vector(norm);
	avg_vector(centroid, other->centroid, area, other->area);
	area += other->area;
	neighbor_union(other);
	neighbors.erase(other);
}

void Plane::neighbor_union(std::shared_ptr<Plane> other){
	for(auto iter = other->neighbors.begin(); iter != other->neighbors.end();
		iter++)
	{
		if(*this != **iter){
			// Don't add self to list of neighbors
			neighbors.insert(*iter);
		}
	}
}

void Plane::add_neighbor(std::shared_ptr<Plane> other){
	neighbors.insert(other);
}

void Plane::remove_neighbor(std::shared_ptr<Plane> other){
	neighbors.erase(other);
}

REAL Plane::score(std::shared_ptr<Plane> other) const{
	assert_len_equal(norm, other->norm);
	REAL sum = 0;
	for(size_t i = 0; i < norm.size(); i++){
		sum += norm[i] * other->norm[i];
	}
	return 1.0 - sum;
}

const plane_set& Plane::get_neighbors() const{
	return neighbors;
}

plane_id Plane::get_id() const{
	return id;
}

bool Plane::operator<(const Plane &other) const{
	return id < other.id;
}

bool Plane::operator==(const Plane &other) const{
	return id == other.id;
}

bool Plane::operator!=(const Plane &other) const{
	return id != other.id;
}

std::ostream& operator<<(std::ostream &out, const Plane &p){
	out << "Plane(" << p.id << ", (";
	for(auto iter = p.norm.begin(); iter != p.norm.end(); iter++){
		out << *iter << ", ";
	}
	out << "), (";
	for(auto iter = p.centroid.begin(); iter != p.centroid.end(); iter++){
		out << *iter << ", ";
	}
	out << "))";
	return out;
}

size_t Plane::hash() const{
	return id;
}

plane_pair get_pair(std::shared_ptr<Plane> a,
	std::shared_ptr<Plane> b)
{
	plane_pair p;
	if (*a < *b){
		p.first = a; p.second = b;
	}else{
		p.first = b; p.second = a;
	}
	return p;
}

size_t PairPointerHash::operator()(const std::shared_ptr<plane_pair> &p) const{
	return p->first->hash() ^ p->second->hash();
}

bool PairDerefCompare::operator()(const std::shared_ptr<plane_pair> &p1,
	const std::shared_ptr<plane_pair> &p2) const
{
	return (*(p1->first) == *(p2->first))
		&& (*(p1->second) == *(p2->second));
}
