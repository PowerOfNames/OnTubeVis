#pragma once

// C++ STL
#include <functional>

// CGV framework core
#include <cgv/math/fvec.h>


/// simple implementation of a 3D regular grid acceleration data structure.
template <class flt_type>
class grid3D
{

public:

	////
	// Exported types

	/// real number type
	typedef flt_type real;

	/// 3D vector type
	typedef cgv::math::fvec<real, 3> vec3;

	/// 4D vector type
	typedef cgv::math::fvec<real, 4> vec4;

	/// point accessor type
	typedef std::function<bool(vec3*, size_t)> point_accessor;

	/// grid cell functionality wrapper
	struct cell
	{
		/// STL-compliant hash functor for grid cells.
		struct hash
		{
			/// uses the 3d vector hash function described in "Optimized Spatial Hashing for Collision
			/// Detection of Deformable Objects" available here:
			/// http://www.beosil.com/download/CollisionDetectionHashing_VMV03.pdf
			size_t operator() (const cell &cell) const
			{
				return size_t(
					cell.x*73856093L ^ cell.y*19349663L ^ cell.z*83492791L
					// ^ mod N is ommitted since it would be (size_t)-1 here, which is
					//   basically a no-op
				);
			}
		};

		/// the @a x grid coordinate of the cell
		long long x;

		/// the @a y grid coordinate of the cell
		long long y;

		/// the @a z grid coordinate of the cell
		long long z;

		/// component-wise addition
		inline void operator+= (const cell &other)
		{
			x += other.x; y += other.y; z += other.z;
			// we intentionally don't support returning the result of applying the operator
		}

		/// equality comparison
		inline bool operator == (const cell &other) const
		{
			return x == other.x && y == other.y && z == other.z;
		}

		/// determines the cell that a given position lies in, given a certain grid cell width
		inline static cell get (const vec3 &position, real cellwidth)
		{
			const bool nx = position.x()<0, ny = position.y()<0, nz = position.z()<0;
			return {(long long)(position.x() / cellwidth) - long(nx),
			        (long long)(position.y() / cellwidth) - long(ny),
			        (long long)(position.z() / cellwidth) - long(nz)};
		}

		/// stores the cell resulting from applying an offset to a reference
		inline static void offset (cell &result, const cell &reference, const cell &offset)
		{
			result.x = reference.x + offset.x;
			result.y = reference.y + offset.y;
			result.z = reference.z + offset.z;
		}
	};


private:

	/// implementation forward
	struct Impl;

	/// implementation handle
	Impl *pimpl;


public:

	/// Default constructor. The grid will not be usable until \ref reset is called to obtain a meaningful state!
	grid3D();

	/// Constructs the grid with the specified cell width and point accessor (move semnantics on the latter).
	grid3D(real cellwidth, point_accessor &&pa);

	/// Constructs the grid with the specified cell width, optionally also specifying the point accessor.
	/// If no point accessor is specified, it needs to be set later on.
	grid3D(real cellwidth, const point_accessor &pa=point_accessor());

	/// the destructor
	~grid3D();

	/// clears all points from the grid
	void clear (void);

	/// rebuilds the grid with the given cell width. Returns true if any points got re-inserted, false otherwise (e.g. if the new
	/// cellwidth was identical to the old or if the grid was empty to begin with).
	bool rebuild (real cellwidth);

	/// clears all points from the grid and changes the cell width to the specified size, optinally assigning a new point accessor.
	void reset (real cellwidth, const point_accessor &pa=point_accessor());

	/// clears all points from the grid and changes the cell width to the specified size and assigns a new point accessor.
	void reset (real cellwidth, point_accessor &&pa);

	/// inserts the point at the given access index into the grid
	void insert (size_t index);

	/// Collects the indices of all currently inserted points within the 27-neighborhood of the cell that
	/// the query point falls in. By default, the points will be reported in ascending order of their
	/// distance from the query point.
	bool query (std::vector<size_t> *out, const vec3 &query_point, bool distance_sort=true) const;

	/// sets the point accessor used to retrieve actual point data
	void set_point_accessor (const point_accessor &pa);

	/// sets the point accessor (move semantics) used to retrieve actual point data
	void set_point_accessor (point_accessor &&pa);
};
