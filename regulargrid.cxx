
// C++ STL
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>

// implemented header
#include "regulargrid.h"


////
// Class implementation

template<class flt_type>
struct grid3D<flt_type>::Impl
{
	// types
	typedef flt_type real;
	typedef typename grid3D::vec3 vec3;

	// fields
	static const typename grid3D<flt_type>::cell n26 [26];
	typename grid3D::point_accessor pnt_access;
	real cellwidth;
	std::unordered_map<
		typename grid3D::cell, std::vector<size_t>,
		typename grid3D::cell::hash
	> grid;

	// helper methods
	Impl(real cellwidth, const grid3D::point_accessor &point_accessor)
		: cellwidth(cellwidth), pnt_access(point_accessor)
	{}
	Impl(real cellwidth, grid3D::point_accessor &&point_accessor)
		: cellwidth(cellwidth), pnt_access(std::move(point_accessor))
	{}
};
template <class flt_type>
const typename grid3D<flt_type>::cell grid3D<flt_type>::Impl::n26 [26] = {
	// 1-distance cells (basically, the 3d-cross centered around the middle cell)
	{-1,  0,  0}, { 1,  0,  0}, { 0,  0, -1}, { 0,  0,  1}, { 0,  1,  0}, { 0, -1,  0},
	// 2-distance cells
	{-1,  0, -1}, { 1,  0,  1}, {-1,  0,  1}, { 1,  0, -1}, // xz-slice
	{-1,  1,  0}, { 1, -1,  0}, { 1,  1,  0}, {-1, -1,  0}, // xy-slice
	{ 0,  1, -1}, { 0, -1,  1}, { 0,  1,  1}, { 0, -1, -1}, // yz-slice
	// 3-distance cells (basically, the outer corners)
	{-1,  1, -1}, { 1, -1,  1}, {-1,  1,  1}, { 1, -1, -1}, {-1, -1, -1}, { 1,  1,  1}, {-1, -1,  1}, { 1,  1, -1}
};

template <class flt_type>
grid3D<flt_type>::grid3D() : pimpl(nullptr)
{}

template <class flt_type>
grid3D<flt_type>::grid3D(real cellwidth, const point_accessor &pa) : pimpl(nullptr)
{
	pimpl = new Impl(cellwidth, pa);
}

template <class flt_type>
grid3D<flt_type>::grid3D(real cellwidth, point_accessor &&pa) : pimpl(nullptr)
{
	pimpl = new Impl(cellwidth, std::move(pa));
}

template <class flt_type>
grid3D<flt_type>::~grid3D()
{
	if (pimpl)
	{
		delete pimpl;
		pimpl = nullptr;
	}
}

template <class flt_type>
void grid3D<flt_type>::clear (void)
{
	pimpl->grid.clear();
}

template <class flt_type>
bool grid3D<flt_type>::rebuild (real cellwidth)
{
	// shortcut for saving one indirection
	auto &impl = *pimpl;

	// don't do anything if nothing changed
	if (impl.cellwidth == cellwidth)
		return false;

	// insert all points currently in the database into the resized grid
	auto grid_old = std::move(impl.grid);
	impl.cellwidth = cellwidth;
	for (const auto &oldcell : grid_old)
		for (const auto index : oldcell.second)
		{
			// retrieve the point
			vec3 point;
			impl.pnt_access(&point, index);

			// insert into appropriate cell
			const cell c = cell::get(point, cellwidth);
			impl.grid[c].push_back(index);
		}

	//done!
	return false;
}

template <class flt_type>
void grid3D<flt_type>::reset (real cellwidth, const point_accessor &pa)
{
	if (pimpl)
	{
		// shortcut for saving one indirection
		auto &impl = *pimpl;

		// reset state
		clear();
		impl.cellwidth = cellwidth;
		if (pa)
			impl.pnt_access = pa;
	}
	else
		pimpl = new Impl(cellwidth, pa);
}

template <class flt_type>
void grid3D<flt_type>::reset (real cellwidth, point_accessor &&pa)
{
	if (pimpl)
	{
		// shortcut for saving one indirection
		auto &impl = *pimpl;

		// reset state
		clear();
		impl.cellwidth = cellwidth;
		if (pa)
			impl.pnt_access = std::move(pa);
	}
	else
		pimpl = new Impl(cellwidth, std::move(pa));
}

template <class flt_type>
void grid3D<flt_type>::insert (size_t index)
{
	// shortcut for saving one indirection
	auto &impl = *pimpl;

	// retrieve the point
	vec3 point;
	impl.pnt_access(&point, index);

	// insert into appropriate cell
	const cell c = cell::get(point, impl.cellwidth);
	impl.grid[c].push_back(index);
}

template <class flt_type>
bool grid3D<flt_type>::query (std::vector<size_t> *out, const vec3 &query_point, bool distance_sort) const
{
	// shortcut for saving one indirection
	const auto &impl = *pimpl;

	// keep track of first new element position in output array
	const size_t sort_start = out->size();

	// determine query cell
	const cell qc = cell::get(query_point, impl.cellwidth);

	// look inside the query cell
	bool found = false;
	{
		auto c = impl.grid.find(qc);
		if (c != impl.grid.end())
		{
			found = true;
			out->insert(out->end(), c->second.begin(), c->second.end());
		}
	}

	// look in the 26-neighborhood of the query cell
	cell q_cur; for (unsigned i=0; i<26; i++)
	{
		cell::offset(q_cur, qc, Impl::n26[i]);
		auto c = impl.grid.find(q_cur);
		if (c != impl.grid.end())
		{
			found = true;
			out->insert(out->end(), c->second.begin(), c->second.end());
		}
	}

	// do exact distance sort if requested
	if (distance_sort)
		std::sort(out->begin()+sort_start, out->end(),
		[&impl, &query_point] (size_t l, size_t r) {
			vec3 lpoint, rpoint;
			impl.pnt_access(&lpoint, l);
			impl.pnt_access(&rpoint, r);
			return   (lpoint - query_point).sqr_length()
			       < (rpoint - query_point).sqr_length();
		});

	return found;
}

template <class flt_type>
void grid3D<flt_type>::set_point_accessor (const point_accessor &pa)
{
	pimpl->pnt_access = pa;
}

template <class flt_type>
void grid3D<flt_type>::set_point_accessor (point_accessor &&pa)
{
	pimpl->pnt_access = std::move(pa);
}



//////
//
// Explicit template instantiations
//

// Only float and double variants are intended
template class grid3D<float>;
template class grid3D<double>;
