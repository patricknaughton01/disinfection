#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <exception>
#include <string>
#include <sstream>
#include <cmath>
#include <memory>
#include <stack>
#include "plane_finder.h"

#define MIN(x, y) ((x)<(y)?(x):(y))
#define MAX(x, y) ((x)>(y)?(x):(y))
#define ABS(x) ((x)>0?(x):(-(x)))

typedef std::pair<std::shared_ptr<Plane>, std::shared_ptr<Plane>> plane_pair;

plane_pair get_pair(const std::shared_ptr<Plane> a,
	const std::shared_ptr<Plane> b);
template <class T>
void print(T s, T e);
template <class T>
void print_d(T s, T e);
template <class T>
void assert_len_equal(T a, T b);

plane_id Plane::gen_id = 0;

int main(){
	plane_set planes;
	REAL area = 1.0;
	std::vector<REAL> n1{0, 0, 1.0};
	std::vector<REAL> n2{0, 0, 1.0};
	std::vector<REAL> n3{1.0, 0, 0};
	std::vector<REAL> c1{0, 0, 0};
	std::vector<REAL> c2{1.0, 0, 0};
	std::vector<REAL> c3{1.0, 0, -1.0};
	auto p1 = std::make_shared<Plane>(n1, c1, area);
	auto p2 = std::make_shared<Plane>(n2, c2, area);
	auto p3 = std::make_shared<Plane>(n3, c3, area);
	p1->add_neighbor(p2);
	p2->add_neighbor(p1);
	p2->add_neighbor(p3);
	p3->add_neighbor(p2);
	planes.insert(p1);
	planes.insert(p2);
	planes.insert(p3);
	plane_set out;
	simplify_planes(planes, out);
	return 0;
}


Plane::Plane():area(0){
	init_id();
}

Plane::Plane(std::vector<REAL> &n, std::vector<REAL> &c, REAL a):norm(n),
	centroid(c), area(a){
	init_id();
}

Plane::Plane(std::vector<REAL> &n, std::vector<REAL> &c,
	plane_set &neigh, REAL a):norm(n), centroid(c), neighbors(neigh),
	area(a)
{
	init_id();
}

void Plane::init_id(){
	id = gen_id;
	gen_id++;
}

struct PairPointerHash{
    size_t operator()(const std::shared_ptr<plane_pair> &p) const{
		return p->first->hash() ^ p->second->hash();
	}
};

struct PairDerefCompare{
	bool operator()(const std::shared_ptr<plane_pair> &p1,
		const std::shared_ptr<plane_pair> &p2) const
	{
		return (*(p1->first) == *(p2->first))
			&& (*(p1->second) == *(p2->second));
	}
};

void simplify_planes(plane_set planes, plane_set out)
{
	// Using a multimap as a priority queue
	std::multimap<REAL, std::shared_ptr<plane_pair>> pq;
	typedef decltype(pq.begin()) pq_iter_t;
	std::unordered_map<std::shared_ptr<plane_pair>, pq_iter_t, PairPointerHash,
		PairDerefCompare> locs;
}

void init_pq(std::multimap<REAL, std::shared_ptr<plane_pair>> &pq,
	plane_set planes)
{
	std::stack<std::shared_ptr<Plane>> open;
	plane_set visited;
	for(auto iter = planes.begin(); iter != planes.end(); iter++){
		if(visited.find(*iter) == visited.end()){
			// Haven't visited this plane yet, do a DFS from here to add
			// neighboring planes to the pq
			open.push(*iter);
			while(open.size()){
				std::shared_ptr<Plane> plane = open.top();
				open.pop();
				for(auto piter = plane->get_neighbors().begin();
					piter != plane->get_neighbors().end(); piter++)
				{
					if (visited.find(*piter) == visited.end()){
						open.push(*piter);
						std::shared_ptr<plane_pair> p =
							std::make_shared<plane_pair>(
								get_pair(plane, *piter));
						pq.insert(std::make_pair(plane->score(**piter), p));
					}
				}
			}
		}
	}
}

void Plane::merge(const Plane &other){
	avg_vector(norm, other.norm, area, other.area);
	avg_vector(centroid, other.centroid, area, other.area);
	area += other.area;
	neighbor_union(other);
}

void Plane::neighbor_union(const Plane &other){
	for(auto iter = other.neighbors.begin(); iter != other.neighbors.end();
		iter++)
	{
		neighbors.insert(*iter);
	}
}

void Plane::add_neighbor(std::shared_ptr<Plane> other){
	neighbors.insert(other);
}

REAL Plane::score(const Plane &other) const{
	assert_len_equal(norm, other.norm);
	REAL sum = 0;
	for(int i = 0; i < norm.size(); i++){
		sum += norm[i] * other.norm[i];
	}
	return sum;
}

const plane_set& Plane::get_neighbors() const{
	return neighbors;
}

plane_pair get_pair(const std::shared_ptr<Plane> a,
	const std::shared_ptr<Plane> b)
{
	plane_pair p;
	if (*a < *b){
		p.first = a; p.second = b;
	}else{
		p.first = b; p.second = a;
	}
	return p;
}

bool Plane::operator<(const Plane &other) const{
	return id < other.id;
}

bool Plane::operator==(const Plane &other) const{
	return id == other.id;
}

size_t Plane::hash() const{
	return id;
}

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
	for(int i = 0; i < a.size(); i++){
		a[i] = (weight_a * a[i] + weight_b * b[i]) / w_sum;
	}
}

template <class T>
void assert_len_equal(T a, T b){
	if(a.size() != b.size()){
		std::stringstream ss;
		ss << "Error in assert_len_equal, length of a does "
			<< "not match length of b ("
			<< a.size() << " != " << b.size() << ")";
		throw std::length_error(ss.str());
	}
}

template <class T>
void print(T s, T e){
	for(auto iter = s; iter != e; iter++){
		std::cout << *iter << " ";
	}
	std::cout << std::endl;
}

template <class T>
void print_d(T s, T e){
	for(auto iter = s; iter != e; iter++){
    	std::cout << "(" << iter->first << "->" << iter->second << ") ";
	}
	std::cout << std::endl;
}

template <typename T>
size_t PointerHash::operator()(const std::shared_ptr<T> &a) const{
	return a->hash();
}

template <typename T>
bool DerefCompare::operator() (std::shared_ptr<T> const &a,
	std::shared_ptr<T> const &b) const
{
	return *a == *b;
}
