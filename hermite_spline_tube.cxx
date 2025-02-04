#pragma once

#include "hermite_spline_tube.h"


const cgv::mat4 hermite_spline_tube::M = transpose(
	cgv::mat4{
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		-3.0f, -2.0f, 3.0f, -1.0f,
		2.0f, 1.0f, -2.0f, 1.0f
	});
