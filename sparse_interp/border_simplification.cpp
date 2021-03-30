#include"border_simplification.h"
#include<igl/triangle/cdt.h>
#include<assert.h>


namespace border {
	
	std::vector<Vector3d> border_remove_redundant_points(const std::vector<Vector3d>&border) {
		std::vector<Vector3d> result;
		result.reserve(border.size());
		
		for (int i = 0; i < border.size() - 1; i++) {
			double dist = (border[i + 1] - border[i]).norm();
			if (dist > SCALAR_ZERO) {
				result.push_back(border[i + 1]);
			}
		}
		double dist = (border.back() - border.front()).norm();
		if (dist > SCALAR_ZERO) {
			result.push_back(border[0]);
		}
		return result;
	}

	// get length of each segment of border, border.size()=length.size()
	// the ith segment is defined by two points border[i] and border[i+1]
	void get_border_pp_length(const std::vector<Vector3d>&border, std::vector<double>& length, double& dis) {
		dis = 0;
		length.resize(border.size());
		for (int i = 0; i < border.size() - 1; i++) {
			length[i] = (border[i] - border[i + 1]).norm();
			dis += length[i];
		}
		length[border.size() - 1] = (border[border.size() - 1] - border[0]).norm();
		dis += length[border.size() - 1];
	}
	void get_local_curvature(const std::vector<Vector3d>& border, std::vector<double> &cur,
		std::vector<double>& length, double &dis) {
		int size = border.size();
		cur.resize(size);
		double angle;
		std::vector<double> curvature;
		get_border_pp_length(border, length, dis);
		curvature.resize(border.size());
		angle = (border[0] - border[size - 1]).dot(border[1] - border[0]) / length[0] / length[size - 1];
		//get the angle from cosine
		angle = std::acos(angle);
		curvature[0] = angle * (length[0] + length[size - 1]) / 2;
		for (int i = 1; i < size - 1; i++) {
			// get the cosine angle
			angle = (border[i] - border[i - 1]).dot(border[i + 1] - border[i]) / length[i] / length[i - 1];
			//get the angle from cosine
			angle = std::acos(angle);
			curvature[i] = angle * (length[i] + length[i - 1]) / 2;
		}

		angle = (border[size-1] - border[size-2]).dot(border[0] - border[size-1]) / length[size-1] / length[size-2];
		//get the angle from cosine
		angle = std::acos(angle);
		curvature[size - 1] = angle * (length[size - 1] + length[size - 2]) / 2;
		std::cout << "curvature get" << std::endl;
		cur = curvature;
	}
	
	


	bool point_in_seg_envelope(const Vector3d &p, const Vector3d &s1, const Vector3d &s2,const double envsize) {
		if ((p - s1).norm() < envsize || (p - s2).norm() < envsize) {
			return true;
		}

		Vector3d dir1 = s2 - s1;
		Vector3d dir2 = p - s1;
		Vector3d dir3 = p - s2;
		double area = dir1.cross(dir2).norm();
		double d1 = fabs(dir1.norm()), d2 = fabs(dir2.norm());
		double dot1 = dir2.dot(dir1), dot2 = -1 * dir1.dot(dir3), dot0 = dir2.dot(dir3);
		if (fabs(dot1) < SCALAR_ZERO) {// p is (almost) on the line. since s1 and s2 are not close
			if (dot0 < 0) {// p is in the middle of the seg
				return true;
			}
			else {
				return false;
			}
		}
		if (dot1 > 0 && dot2 > 0) {// p can be project to the seg
			double cos = dot1 / (d2*d1);
			double sin = sqrt(1 - cos * cos);
			if (d2*sin > envsize) {
				return false;
			}
			return true;
		}
		else {
			return false;
		}


	}
	// pid2>=border.size() is allowed
	bool envelope_check(const std::vector<Vector3d>& border, const int pid1, const int pid2, const double tolerance) {
		int howmany;// how many numbers are used in total
		if (pid2 > pid1) {
			howmany = pid2 - pid1 + 1;
		}
		else {
			howmany = border.size() - pid1 + pid2;
		}
		if (howmany == 2) return true;// if two neighbouring points are the two endpoints, then no need to check
		for (int i = 1; i < howmany - 1; i++) {
			int id2;
			if (pid1+i>border.size()-1) {
				id2 = pid1 + i - border.size();
			}
			else {
				id2 = pid1 + i;
			}
			if (!point_in_seg_envelope(border[id2], border[pid1], border[pid2], tolerance)) {
				return false;
			}
			
		}
		return true;
	}

	// give pid, find another point pid1, accroding to the two conditions
	void find_next_point_with_features(const int &pid, const std::vector<Vector3d>& border,
		const std::vector<double>& cur,
		const std::vector<double>& length, const double  d1, const double d2, int &pid1) {

		double segdis = 0;
		double maxcur = 0;
		int id = pid+1;
		int i1;// current id
		for (int i = 1; i < border.size(); i++) {
			int cid = pid + i;
			if (cid >= border.size()) i1 = cid - border.size();
			else i1 = cid;
			segdis += length[i1];
			if (segdis > d2) break;
			// select the point with maximal curvature in a certain region
			if (segdis >= d1 && fabs(cur[i1]) >= maxcur) {
				maxcur = cur[i1];
				id = i1;
			}
		}
		pid1 = id;

	}
	class loopinfo {
	private:
		bool is_id_setted = false;
		int id;
		int size;
	
	public:
		loopinfo(const int bordersize) {
			size = bordersize;
		}
		void set_id(int i) {
			id = i;
			is_id_setted = true;

		}
		int real_id(){
			if (!is_id_setted) {
				assert(false);
				std::cout << "ERROR USAGE: ID REQUEST BEFORE DEFINING" << std::endl;
			}
			if (id <= size - 1) {
				return id;
			}
			else {
				return id - size;
			}
		}
		int extended_id() {
			if (!is_id_setted) {
				assert(false);
				std::cout << "ERROR USAGE: ID REQUEST BEFORE DEFINING" << std::endl;
			}
			return id + size;
		}
	};
	bool id_in_between(const int id, const int start, const int end) {
		if (start < end) {
			if (id > start&&id < end) {
				return true;
			}
			return false;
		}
		else {
			if (id > start || id < end) {
				return true;
			}
			return false;
		}
	}
	std::vector<int> insert_a_point_to_id_list(const int border_size, const std::vector<int>&list, const int pid) {
		loopinfo lp(border_size);
		lp.set_id(pid);
		int realid = lp.real_id();

		std::vector<int> result;
		result.reserve(list.size());
		
		for (int i = 0; i < list.size() - 1; i++) {
			result.push_back(list[i]);
			if (id_in_between(realid,list[i],list[i+1])) {
				result.push_back(realid);
			}
		}
		result.push_back(list.back());
		if (id_in_between(realid, list.back(), list.front())) {
			result.push_back(realid);
		}
		assert(result.size() == list.size() + 1);
		return result;
	}
	// give an original border, a simplified feature list features, and a segment whose endpoints are border[segid0], border[segid1]. 
	// insert a new point border[p], s.t. segment border[segid0],border[segid0]
	// howmany_in_between is how many points in (segid0, segid1)
	std::vector<int> insert_a_point_to_refine(const std::vector<Vector3d>& border, const std::vector<int> &features,
		const int segid0, const int segid1, const double tolerance, const int howmany_in_between) {
		
		int pid;
		int end = howmany_in_between * 3 / 4;// to make the inserted point not too close to segid1.
		for (int i = 1;; i++) {
			bool check = envelope_check(border, segid0, segid0 + i, tolerance);
			if (check) {
				pid = segid0 + i;
			}
			if (i > end || !check) {// if too close to segid2 or not in envelope, break
				break;
			}
		}
		// now, pid is the point that we want to insert
		std::vector<int> result = insert_a_point_to_id_list(border.size(), features, pid);
		assert(result.size() == features.size() + 1);
		assert(pid < segid1);
		return result;
	}


	// refine the border until the simplified border is within envelope
	void refine_border_with_envelope(const std::vector<Vector3d>& border, const std::vector<int>& id_in,
		std::vector<int>& id_out, const double tolerance) {
		std::vector<int> tmpid = id_in;

		for (int i = 0; i < tmpid.size(); i++) {
			int howmany;// how many points in between two feature points
			int t = i + 1 > tmpid.size() - 1 ? 0 : i + 1;
			loopinfo endpoint(tmpid.size());
			int id1 = tmpid[i];// the segment ids are [id1, id2]
			int id2 = tmpid[t];
			bool check = envelope_check(border, id1, id2, tolerance);
			if (check) {// if this seg is within tolerance, then skip
				continue;
			}
			endpoint.set_id(id2);
			if (id2 > id1) {
				howmany = id2 - id1 - 1;
			}
			else {
				howmany = endpoint.extended_id() - id1 - 1;
			}
			print_vector(tmpid);
			if(howmany<=0){
				std::cout << "wrong howmany, " << howmany <<" id1, "<<id1<<" id2, "<<id2<< " i= "<<i<< std::endl;
				std::cout << "tmpid size " << tmpid.size() << std::endl;
			}

			assert(howmany > 0);
			// down here: if this segment is out of envelope, refine this segment
			tmpid = insert_a_point_to_refine(border, tmpid, id1, id2, tolerance, howmany);
		}
		id_out = tmpid;
	}


	// get simplified points of the border. d1 is the expected distance of two points, d2 is the maximum search distance
	// constrain 1, the simplified points not too close to each other, constrain 2, the simplified shape is 
	// within an epsilon envelope.
	void get_simp_points(const std::vector<Vector3d>& border, const double & d1, const double d2,
		const double tolerance, std::vector<int>& pid) {
		
		std::vector<double >cur;
		std::vector<double> length; double dis;
		get_local_curvature(border, cur, length, dis);

		std::vector<int> features;
		
		int old_id = 0, new_id;
		while(1){
			find_next_point_with_features(old_id, border, cur, length, d1, d2, new_id);
			features.push_back(new_id);
			if (old_id > new_id) {// if old_id > new_id, then the whole loop is searched
				break;
			}
			old_id = new_id;

		}
		if (features.back() > features.front()) {// avoid overlap
			features.pop_back();
		}

		refine_border_with_envelope(border, features, pid, tolerance);

	}

	// use expected point number to decide the sample length
	void get_simp_points(const std::vector<Vector3d>& border, const int expect_p_nbr,
		const double tolerance, std::vector<int>& pids) {
		std::vector<double >cur;
		std::vector<double> length; double dis;
		std::cout << "before curvature get" << std::endl;
		get_local_curvature(border, cur, length, dis);
		std::cout << "curvature get" << std::endl;
		const double  d1 = dis/expect_p_nbr;
		const double d2 = 2 * d1;
		std::vector<int> features;

		int old_id = 0, new_id;
		while (1) {
			find_next_point_with_features(old_id, border, cur, length, d1, d2, new_id);
			features.push_back(new_id);
			if (old_id > new_id) {// if old_id < new_id, then the whole circle is searched
				break;
			}
			old_id = new_id;

		}
		if (features.back() > features.front()) {// avoid overlap
			features.pop_back();
		}
		std::cout << "features get, size "<<features.size() << std::endl;
		print_vector(features);
		refine_border_with_envelope(border, features, pids, tolerance);
	}






}







