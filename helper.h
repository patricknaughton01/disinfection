#pragma once
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
REAL normal(std::vector<REAL> &a);
void divide_vector(std::vector<REAL> &a, REAL d);
void normalize_vector(std::vector<REAL> &a);

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
