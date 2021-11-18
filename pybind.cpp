#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "merge_triangle_mesh.h"
#include "plane_finder.h"
#include "heightmap.h"

PYBIND11_MODULE(merge_triangle_mesh, m){
	m.doc() = "Module for doing things with triangle meshes";
	m.def("get_empty_plane_finder", &get_empty_plane_finder,
		"Get an empty PlaneFinder");
	m.def("merge_triangle_mesh", &merge_triangle_mesh,
		"Merge a triangle mesh");
	m.def("build_heightmaps", &build_heightmaps,
		"Build and store heightmaps for a mesh");
	m.def("get_heightmaps", &get_heightmaps,
		"Get the heightmaps from a mesh");
	m.def("get_heightmap_metadata", &get_heightmap_metadata,
		"Get the metadata from each heightmap");
	m.def("heightmaps", &heightmaps,
		"Get the heightmap objects");
	m.def("dedup_triangle_mesh", &dedup_triangle_mesh,
		"Deduplicate vertices in a triangle mesh");
	pybind11::class_<PlaneFinder>(m, "PlaneFinder");
	pybind11::class_<Heightmap>(m, "Heightmap");
}
