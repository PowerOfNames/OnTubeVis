#pragma once

// C++ STL
#include <iostream>
#include <vector>

// CGV framework core
#include <cgv/math/fvec.h>
#include <cgv/math/fmat.h>
#include <cgv/media/color.h>

// local includes
#include "traj_loader.h"


/// provides read and write capabilites for Hermite splines in .bezdat format
template <class flt_type>
struct bezdat_handler : public traj_format_handler<flt_type>
{
	/// real number type
	typedef traj_format_handler::real real;

	/// 2D vector type
	typedef traj_format_handler::Vec2 Vec2;

	/// 3D vector type
	typedef traj_format_handler::Vec3 Vec3;

	/// 4D vector type
	typedef traj_format_handler::Vec4 Vec4;

	/// rgb color type
	typedef traj_format_handler::Color Color;

	/// test if the given data stream appears to be a .bezdat file
	virtual bool can_handle (std::istream &contents) const;

	/// parse the given stream containing the .bezdat file contents and report whether any data was loaded
	virtual traj_dataset<real> read (std::istream &contents, DatasetOrigin source, const std::string &path);
};
