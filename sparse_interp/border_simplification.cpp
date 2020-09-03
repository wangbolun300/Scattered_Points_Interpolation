#include"border_simplification.h"
#include<igl/triangle/cdt.h>

namespace workClean{
	namespace border {
		// get simplified points of the border. d1 is the expected distance of two points, d2 is the maximum search distance
		// constrain 1, the simplified points not too close to each other, constrain 2, the simplified shape is 
		// within an epsilon envelope.
		void get_simp_points(const std::vector<Vector3>& border, const double & d1, const double d2, std::vector<int>& pid) {
			std::vector<Scalar >cur;
			std::vector<Scalar> length; Scalar dis;
			get_local_curvature(border, cur, length, dis);
			
			// pick up the maximal curvature point
			Scalar cmax = 0;
			int cid = 0;
			for (int i = 0; i < border.size(); i++) {
				if (fabs(cur[i]) > cmax) {
					cmax = cur[i];
					cid = i;
				}
			}

			
			


		}

		// give pid, find another point pid1, accroding to the two conditions
		void find_next_point(const int &pid, const std::vector<Vector3>& border,
			const std::vector<Scalar>& cur,
			const std::vector<Scalar>& length, const double  d1, const double d2, int &pid1) {
			
			double segdis = 0;
			double maxcur = 0;
			int id = pid;
			int i1;// current id
			for (int i = 0; i < border.size(); i++) {
				int cid = pid + i;
				if (cid >= border.size()) i1 = cid - border.size();
				else i1 = cid;
				segdis += length[pid + i];
				if (segdis > d2) break;
				if (fabs(cur[i]) > maxcur) {
					maxcur = cur[i];
					id = i;
					//TODO check envelope here, check segdis>d1 here
				}
			}
		}

		// get length of each segment of border, border.size()=length.size()
		// the ith segment is defined by two points border[i] and border[i+1]
		void get_border_pp_length(const std::vector<Vector3>&border, std::vector<Scalar>& length, Scalar& dis) {
			dis = 0;
			length.resize(border.size());
			for (int i = 0; i < border.size()-1; i++) {
				length[i] = (border[i] - border[i + 1]).norm();
				dis += length[i];
			}
			length[border.size() - 1] = (border[border.size() - 1] - border[0]).norm();
			dis += length[border.size() - 1];
		}

		void get_local_curvature(const std::vector<Vector3>& border, std::vector<Scalar> &cur, 
			std::vector<Scalar>& length, Scalar &dis) {
			int size = border.size();
			cur.resize(size);
			Scalar angle;
			std::vector<Scalar> curvature;
			get_border_pp_length(border, length, dis);
			curvature.resize(border.size());
			angle = (border[0] - border[size - 1]).dot(border[1] - border[0]) / length[0] / length[size - 1];
			//get the angle from cosine
			angle = std::acos(angle);
			curvature[0] = angle * (length[0] + length[size - 1]) / 2;
			for (int i = 1; i < size; i++) {
				// get the cosine angle
				angle = (border[i] - border[i - 1]).dot(border[i + 1] - border[i]) / length[i] / length[i - 1];
				//get the angle from cosine
				angle = std::acos(angle);
				curvature[i] = angle * (length[i] + length[i - 1]) / 2;
			}
		}


	
	
	
	
	
	}
}



	


