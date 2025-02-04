#pragma once

#include <cgv/math/functions.h>

#include "quadratic_bezier_tube.h"
#include <cgv/math/fvec.h>
#include <cgv/math/fmat.h>

class hermite_spline_tube {
public:
	using vec3 = cgv::vec3;
	using vec4 = cgv::vec4;
	using box3 = cgv::box3;
	using mat43 = cgv::math::fmat<float, 4, 3>;

	struct node {
		vec3 pos;
		float rad;
		vec3 pos_tan;
		float rad_tan;
	};

	node a;
	node b;

private:
	static const cgv::mat4 M;

	template<typename T>
	void split(unsigned segment_idx, T v0, T d0, T v1, T d1, T& v0_out, T& h_out, T& v1_out) const {

		v0_out = v0;
		h_out = v0 + d0 * 0.25f;
		T h1 =  v1 - d1 * 0.25f;
		v1_out = (h_out + h1) * 0.5f;

		if(segment_idx == 1) {
			v0_out = v1_out;
			h_out = h1;
			v1_out = v1;
		}
	}

public:
	hermite_spline_tube(const node& a, const node& b) : a(a), b(b) {};

	hermite_spline_tube(const vec3& pa, const vec3& pb, const float ra, const float rb, const vec3& pta, const vec3& ptb, const float& rta, const float& rtb) {
		a = { pa, ra, pta, rta };
		b = { pb, rb, ptb, rtb };
	};

	box3 bounding_box(bool exact) {
		box3 box;
		
		quadratic_bezier_tube q_tube0 = split_to_quadratic_bezier_tube(0);
		quadratic_bezier_tube q_tube1 = split_to_quadratic_bezier_tube(1);
		
		box.add_axis_aligned_box(q_tube0.bounding_box(exact));
		box.add_axis_aligned_box(q_tube1.bounding_box(exact));
		
		return box;
	}

	quadratic_bezier_tube split_to_quadratic_bezier_tube(unsigned segment_id) const {
		vec3 pa, pb, pc;
		float ra, rb, rc;

		split(segment_id,
			a.pos, a.pos_tan, b.pos, b.pos_tan,
			pa, pb, pc
		);
		split(segment_id,
			a.rad, a.rad_tan, b.rad, b.rad_tan,
			ra, rb, rc
		);

		return quadratic_bezier_tube(pa, pb, pc, ra, rb, rc);
	}

	vec3 interpolate(float t) const {
		vec4 T(1.0f, t, t * t, t * t * t);

		mat43 B;
		B.set_row(0, a.pos);
		B.set_row(1, a.pos_tan);
		B.set_row(2, b.pos);
		B.set_row(3, b.pos_tan);

		return T * M * B;
	}

	std::vector<vec3> sample(size_t num_segments) const {
		num_segments = std::max(num_segments, size_t(1));

		std::vector<vec3> points;
		float step = 1.0f / static_cast<float>(num_segments);
		for (size_t i = 0; i <= num_segments; ++i) {
			float t = step * static_cast<float>(i);
			points.push_back(interpolate(t));
		}

		return points;
	}

};
