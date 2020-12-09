#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "merge_triangle_mesh.h"

PYBIND11_MODULE(merge_triangle_mesh, m){
	m.doc() = "Module for doing things with triangle meshes";
	m.def("merge_triangle_mesh", &merge_triangle_mesh,
		"Merge a triangle mesh");
	m.def("dedup_triangle_mesh", &dedup_triangle_mesh,
		"Deduplicate vertices in a triangle mesh");
}
