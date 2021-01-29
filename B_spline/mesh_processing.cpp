#include <mesh_processing.h>


void test_read_mesh(const std::string &filename) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool read = igl::read_triangle_mesh(filename, V, F);
}