#include <mesh_processing.h>
#include<array>
#include<curve.h>

void test_read_mesh(const std::string &filename) {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool read = igl::read_triangle_mesh(filename, V, F);

}

void read_and_visual_mesh(const std::string &filename,Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	bool read = igl::read_triangle_mesh(filename, V, F);
}

void map_vertices_to_square(const Eigen::MatrixXd& V, const Eigen::VectorXi& bnd, Eigen::MatrixXd &bnd_uv){
	bnd_uv.resize(bnd.size(), 2);
	double total_dist = 0;
	std::vector<double> pdis(bnd.size());
	pdis[0] = (V.row(bnd[1]) - V.row(bnd[0])).norm();
	for (int i = 1; i < pdis.size() - 1; i++) {
		pdis[i] = (V.row(bnd[i + 1]) - V.row(bnd[i])).norm() + pdis[i - 1];
		
	}
	pdis[pdis.size() - 1] = pdis[pdis.size() - 2] + (V.row(bnd[0]) - V.row(bnd[bnd.size() - 1])).norm();
	total_dist = pdis[pdis.size() - 1];

	int k = 1;
	double quarter = total_dist / 4;
	std::array<double, 4> corners;
	corners[0] = 0;
	for (int i = 0; i < pdis.size()-1; i++) {
		if (pdis[i] <= quarter * k&&pdis[i + 1] > quarter*k) {
			if (i != corners[k - 1]) {
				corners[k] = i;
			}
			else {
				corners[k] = i + 1;
			}
			k += 1;
			if (k == 4) break;
		}
	}
	// u=0, v increase
	bnd_uv.row(0) << 0, 0;
	for (int i = 0; i < corners[1]; i++) {
		bnd_uv.row(i + 1) << 0, pdis[i] / pdis[corners[1]-1];
	}

	// v=1, u increase
	for (int i =corners[1]; i < corners[2]; i++) {
		bnd_uv.row(i + 1) << (pdis[i] - pdis[corners[1] - 1])/( pdis[corners[2] - 1]- pdis[corners[1] - 1]), 1;
	}

	// u=1, v decrease
	for (int i = corners[2]; i < corners[3]; i++) {
		bnd_uv.row(i + 1) << 1, 1 - (pdis[i] - pdis[corners[2] - 1]) / (pdis[corners[3] - 1] - pdis[corners[2] - 1]);
	}

	// v=0, u decrease
	for (int i = corners[3]; i < pdis.size()-1; i++) {
		bnd_uv.row(i + 1) << 1 - (pdis[i] - pdis[corners[3] - 1]) / (pdis[pdis.size()-1] - pdis[corners[3] - 1]), 0;
	}
	//std::cout << "the parameters\n"<<bnd_uv << std::endl;
	

}