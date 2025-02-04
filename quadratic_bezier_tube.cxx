#pragma once

#include "quadratic_bezier_tube.h"

const cgv::mat3 quadratic_bezier_tube::M = cgv::mat3{
		1.0f, -2.0f, 1.0f,
		0.0f, 2.0f, -2.0f,
		0.0f, 0.0f, 1.0f
};
