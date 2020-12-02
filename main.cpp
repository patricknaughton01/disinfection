#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "merge_planes.h"

void load_vector(std::vector<std::vector<REAL>> &v, int dim=0);

int main(){
	// plane_set planes;
	// REAL area = 1.0;
	// std::vector<REAL> n1{0.1, 0, 1.0};
	// std::vector<REAL> n2{0, 0, 1.0};
	// std::vector<REAL> n3{1.0, 0, 0};
	// std::vector<REAL> c1{0, 0, 0};
	// std::vector<REAL> c2{1.0, 0, 0};
	// std::vector<REAL> c3{1.0, 0, -1.0};
	// auto p1 = std::make_shared<Plane>(n1, c1, area);
	// auto p2 = std::make_shared<Plane>(n2, c2, area);
	// auto p3 = std::make_shared<Plane>(n3, c3, area);
	// p1->add_neighbor(p2);
	// p2->add_neighbor(p1);
	// p2->add_neighbor(p3);
	// p3->add_neighbor(p2);
	// planes.insert(p1);
	// planes.insert(p2);
	// planes.insert(p3);
	// plane_set out;
	// PlaneFinder pf(planes);
	// print_p(planes.begin(), planes.end());
	// pf.simplify_planes(out, 0.5);
	// print_p(planes.begin(), planes.end());
	// print_p(out.begin(), planes.end());


	std::vector<std::vector<REAL>> v;
	load_vector(v);
	std::vector<std::vector<size_t>> inds;
	for(int i = 0; i < 12; i++){
		std::vector<size_t> tmp;
		if(i < 6){
			tmp.push_back(0);
			tmp.push_back(i+1);
		}else{
			tmp.push_back(7);
			tmp.push_back(i-5);
		}
		inds.push_back(tmp);
	}
	inds[0].push_back(5);
	inds[1].push_back(3);
	inds[2].push_back(1);
	inds[3].push_back(6);
	inds[4].push_back(4);
	inds[5].push_back(2);

	inds[6].push_back(3);
	inds[7].push_back(6);
	inds[8].push_back(2);
	inds[9].push_back(5);
	inds[10].push_back(1);
	inds[11].push_back(4);
	// for(auto iter = inds.begin(); iter != inds.end(); iter++){
	// 	print(iter->begin(), iter->end());
	// }
	PlaneFinder pf;
	pf.load_triangle_mesh(v, inds);
	std::cout << "Initial Planes: " << std::endl;
	std::cout << pf.planes.size() << std::endl;
	for(auto iter = pf.planes.begin(); iter != pf.planes.end(); iter++){
		std::cout << **iter << std::endl;
	}
	plane_set out;
	pf.simplify_planes(out, 0.5);

	std::cout << "Final Planes: " << std::endl;
	std::cout << out.size() << std::endl;
	for(auto iter = out.begin(); iter != out.end(); iter++){
		std::cout << **iter << std::endl;
	}
	return 0;
}

void load_vector(std::vector<std::vector<REAL>> &v, int dim){
	if(dim == 2){
		for(int i = 0; i < 2; i++){
			std::vector<REAL> t(1, i);
			v.push_back(t);
		}
	}else{
		for(int i = 0; i < 2; i++){
			std::vector<std::vector<REAL>> tmp;
			load_vector(tmp, dim+1);
			for(size_t j = 0; j < tmp.size(); j++){
				tmp[j].push_back(i);
			}
			for(auto iter = tmp.begin(); iter != tmp.end(); iter++){
				v.push_back(*iter);
			}
		}
	}
}
