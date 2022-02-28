#include"basis.h"
#include"curve.h"
#include <igl/opengl/glfw/Viewer.h>
#include<iostream>
#include<array>
#include"test.h"
#include<mesh_processing.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>
#include<exact_calculation.h>
#include<energy.h>
#include <fstream>
#include<integration.h>


std::string root_path(SI_MESH_DIR);
// run the example that insert lots of knot locally and refine it based on the energy.
void insert_redundant(const std::string path) {
	
	int model = 3;
	int nbr = 500;
	double par = 0.9;// ours
	double per = 0.5;
	ours_inserting_redundant(model, nbr, par, path, "", per, true);
}

// run the example which is spiky because of too few data points, and progressively smooth it by inserting knots.
void run_local_refine(const std::string path) {
	int pnbr = 30;
	double per = 0.5;
	double delta = 0.9;
	int model = 0;
	std::string tail = "";
	run_ours(model, pnbr, delta, path, tail, per,false);
	run_ours(model, pnbr, delta, path, tail, per, true);

}

// method = 0, peak; 
// method = 1, drop; 
// method = 2, hyperbolic; 
// method = 3, sinus;
// method = 4, bilinear
// method = 5, snail
// run the randomly sampled examples
void run_ours_over_all_models(double delta, double per, int pnbr, const std::string outpath) {
	std::vector<int> models = { {0,1,2,3,4,5} };
	
	for (int i = 0; i < models.size(); i++) {
		std::string tail = "";
		run_ours(models[i], pnbr, delta, outpath, tail, per, false);
	}
	std::cout << "files are written into " << outpath << std::endl;
}

void mesh_reconstruction(const std::string model, const std::string outpath) {
	

	double par = 0.9;// ours
	double per = 0.5;
	
	std::string tail = "";
	
	run_mesh_reconstruction(root_path,model, par, outpath, tail, per, false);
}


int main(int argc, char *argv[]) {
	//bool run_lofting_method = false;
	//if (run_lofting_method) {
	//	run_lofting(0, 500, 0.5, "./");
	//}
	//else {
	//	int model = 5;
	//	int nbr_pts = 100;
	//	int max_itr = 200;
	//	double threadshold = 1e-7;
	//	int nbr_sqrt = std::max(sqrt(nbr_pts), 4.0);// at least there should be 4 control points; 
	//	const std::string meshfile =
	//		//"D:\\vs\\sparse_data_interpolation\\sparse_data\\external\\meshes\\mask3kf.obj";
	//		"D:\\vs\\sparse_data_interpolation\\sparse_data\\external\\meshes\\tiger.obj";
	//	//run_pia_mesh_reconstruct( meshfile, max_itr, threadshold, 60, "./"); 
	//	run_pia(model, nbr_pts, max_itr, threadshold, nbr_sqrt, "./");
	//}
	mesh_reconstruction("mask3kf.obj", "./");
	

	/*double delta = 0.9;
	int algorithm = 0;
	int model = 1;
	run_ours(model, 100, delta, "./", "", 0.5, false);*/
	//run_other_sampling_strategies(model, delta, "./", "grid", 0.5, algorithm);
	
	std::cout << "done"<< std::endl;

	exit(0);
	std::string type = argv[1];

	// run the benchmark models, the arguments are:
	// ./Sparse_Interp_bin benchmarks <delta> <per> <number of points> <output path>
	if (type == "benchmarks") {
		double delta = std::stod(argv[2]);
		double per = std::stod(argv[3]);
		int nbr = std::stoi(argv[4]);
		std::string out = argv[5];
		run_ours_over_all_models(delta, per, nbr, out);
	}

	// run the mesh interpolation models, the arguments are:
	// ./Sparse_Interp_bin meshes <model> <output path>, where the <model> select from "mask3kf.obj" and "tiger.obj"
	if (type == "meshes") {
		std::string model = argv[2];
		std::string out = argv[3];
		mesh_reconstruction(model, out);
	}

	// improve the local non-smoothness caused by inserting too many knots in a small region. the arguments are:
	// ./Sparse_Interp_bin reproduce <output path>.
	if (type == "reproduce") {
		std::string out = argv[2];
		insert_redundant(out);
	}

	// improve the local non-smoothness caused by too few data points. the arguments are:
	// ./Sparse_Interp_bin local_energy <output path>.
	if (type == "local_energy") {
		std::string out = argv[2];
		run_local_refine(out);
	}
	std::cout << "done!" << std::endl;
	return 0;
}