#include <iostream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include "helper.h"
#include "plane.h"
#include "plane_finder.h"
#include "merge_triangle_mesh.h"

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
	std::vector<REAL> tmp{0,0,0};
	v.push_back(tmp);
	std::cout << "Vertices: " << std::endl;
	for(auto it = v.begin(); it != v.end(); it++){
		print(it->begin(), it->end());
	}
	std::vector<std::vector<size_t>> inds;
	// std::vector<size_t> tmp{0, 1, 5};
	// inds.push_back(tmp);
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
	inds[0][0] = 8;
	dedup_triangle_mesh(v, inds);
	PlaneFinder pf;
	std::pair<std::vector<std::vector<std::vector<REAL>>>,
		std::vector<std::vector<plane_id>>> out_p = merge_triangle_mesh(
			pf, v, inds, 0.5
	);
	// std::cout << "Number of planes " << pf.planes.size() << std::endl;
	std::vector<std::vector<std::vector<REAL>>> out_planes = out_p.first;
	for(auto it = out_planes.begin(); it != out_planes.end(); it++){
		std::cout << "Plane " << (it - out_planes.begin()) << std::endl;
		print((*it)[0].begin(), (*it)[0].end());
		print((*it)[1].begin(), (*it)[1].end());
	}
	// std::cout << "Number of planes " << pf.planes.size() << std::endl;
	build_heightmaps(pf, 0.3, 0.1);
	std::pair<std::vector<std::vector<std::vector<std::vector<REAL>>>>,
		std::vector<std::vector<std::vector<std::vector<REAL>>>>> hm =
		get_heightmaps(pf);
	std::vector<std::vector<REAL>> data = get_heightmap_metadata(pf);
	// std::cout << "Number of heightmaps " << hm.first.size() << " "
	// 	<< hm.second.size() << std::endl;
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
