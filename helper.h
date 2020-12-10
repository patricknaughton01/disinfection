#pragma once
#include <vector>
#include <memory>
#include <sstream>
#include <exception>
#include <iostream>

#define REAL double

#define MIN(x, y) ((x)<(y)?(x):(y))
#define MAX(x, y) ((x)>(y)?(x):(y))
#define ABS(x) ((x)>0?(x):(-(x)))

void avg_vector(std::vector<REAL> &a, const std::vector<REAL> &b,
	REAL weight_a=1.0, REAL weight_b=1.0);
REAL get_norm(const std::vector<REAL> &a);
REAL get_dist(const std::vector<REAL> &a, const std::vector<REAL> &b);
void divide_vector(std::vector<REAL> &a, REAL d);
void normalize_vector(std::vector<REAL> &a);
std::vector<REAL> cross(const std::vector<REAL> &a, const std::vector<REAL> &b);

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

template <class T>
void print_p(T s, T e){
	for(auto iter = s; iter != e; iter++){
		std::cout << **iter << " ";
	}
	std::cout << std::endl;
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

struct PointerHash{
	template <typename T>
    size_t operator()(const std::shared_ptr<T> &a) const{
		return a->hash();
	}
};

struct DerefCompare {
	template <typename T>
	bool operator() (std::shared_ptr<T> const &a,
		std::shared_ptr<T> const &b) const
	{
		return *a == *b;
	}
};

struct WeakPointerHash{
	template <typename T>
    size_t operator()(const std::weak_ptr<T> &a) const{
		if(auto spt = a.lock()){
			return spt->hash();
		}else{
			return 0;
		}
	}
};

struct WeakDerefCompare {
	template <typename T>
	bool operator() (std::weak_ptr<T> const &a,
		std::weak_ptr<T> const &b) const
	{
		auto spta = a.lock();
		auto sptb = b.lock();
		if(spta && sptb){
			return *spta == *sptb;
		}else{
			return false;
		}
	}
};

// Modified from https://stackoverflow.com/a/2595226/10193734
struct VectorHash {
	template <typename T>
	size_t operator()(const std::vector<T>& v) const {
		std::hash<T> hasher;
		size_t seed = 0;
		for (int i : v) {
			seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		}
		return seed;
	}
};

template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b){
	assert_len_equal(a, b);
	std::vector<T> r;
	for(size_t i = 0; i < a.size(); i++){
		r.push_back(a[i] + b[i]);
	}
	return r;
}

template<typename T>
std::vector<T> operator*(const T &a, const std::vector<T> &b){
	std::vector<T> r;
	for(size_t i = 0; i < b.size(); i++){
		r.push_back(a * b[i]);
	}
	return r;
}

template<typename T>
std::vector<T> operator*(const std::vector<T> &b, const T &a){
	return a*b;
}

template<typename T>
std::vector<T> operator/(const std::vector<T> &b, const T &a){
	return  (((T)1)/a) * b;
}

template<typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b){
	assert_len_equal(a, b);
	return a + (-((T)1) * b);
}

template<typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b){
	assert_len_equal(a, b);
	T r = 0;
	for(size_t i = 0; i < a.size(); i++){
		r += a[i] * b[i];
	}
	return r;
}
