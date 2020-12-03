#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "merge_triangle_mesh.h"

PYBIND11_MODULE(merge_triangle_mesh, m){
	m.doc() = "Merge a triangle mesh";
	m.def("merge_triangle_mesh", &merge_triangle_mesh, "A function to merge a triangle mesh");
}
