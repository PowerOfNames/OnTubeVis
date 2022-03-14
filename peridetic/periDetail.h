//
//
// MIT License
//
// Copyright (c) 2020 Stellacore Corporation.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
// KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
// AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//


#ifndef periDetail_INCL_
#define periDetail_INCL_


#include <array>
#include <cmath>
#include <numeric>


//! \file periDetail.h Implementation details for Peridetic.


// utilities
namespace peri
{
	//! Use to indicate invalid data values
	constexpr double sNan{ std::numeric_limits<double>::quiet_NaN() };

	//! Invalid triplet
	constexpr std::array<double, 3u> sNull{ sNan, sNan, sNan };

	//! Classic square operation (value times itself)
	template <typename Type>
	inline
	Type
	sq  // peri::
		( Type const & value
		)
	{
		return (value * value);
	}

	//! "Vector addition" for two std::array data types
	inline
	std::array<double, 3u>
	operator+  // peri::
		( std::array<double, 3u> const & valuesA
		, std::array<double, 3u> const & valuesB
		)
	{
		return
			{ valuesA[0] + valuesB[0]
			, valuesA[1] + valuesB[1]
			, valuesA[2] + valuesB[2]
			};
	}

	//! "Vector 'subtraction'" for two std::array data types
	inline
	std::array<double, 3u>
	operator-  // peri::
		( std::array<double, 3u> const & valuesA
		, std::array<double, 3u> const & valuesB
		)
	{
		return
			{ valuesA[0] - valuesB[0]
			, valuesA[1] - valuesB[1]
			, valuesA[2] - valuesB[2]
			};
	}

	//! Unitary negation
	inline
	std::array<double, 3u>
	operator-  // peri::
		( std::array<double, 3u> const & values
		)
	{
		return
			{ -values[0]
			, -values[1]
			, -values[2]
			};
	}

	//! "scalar-Vector" multiplication for two std::array data types
	inline
	std::array<double, 3u>
	operator*  // peri::
		( double const & scale
		, std::array<double, 3u> const & values
		)
	{
		return
			{ scale * values[0]
			, scale * values[1]
			, scale * values[2]
			};
	}

	//! Vector dot product of two arrays
	inline
	double
	dot  // peri::
		( std::array<double, 3u> const & vecA
		, std::array<double, 3u> const & vecB
		)
	{
		return std::inner_product
			( vecA.begin(), vecA.end()
			, vecB.begin(), 0.
			);
		/*
		// C++17 syntax
		return std::inner_product
			( std::cbegin(vecA), std::cend(vecA)
			, std::cbegin(vecB), 0.
			);
		*/
	}

	//! Squared magnitude of vec Sum of squared components
	inline
	double
	magSq  // peri::
		( XYZ const & vec
		)
	{
		return dot(vec, vec);
	}


	//! Magnitude of vec (square root of sum of squared components)
	inline
	double
	magnitude  // peri::
		( XYZ const & vec
		)
	{
		return std::sqrt(magSq(vec));
	}

	/*! \brief Unitary direction associated with non-zero orig
	 *
	 * \note Specialzied for non-zero orig vector. There is no check
	 * for zero magnitude inputs.
	 *
	 */
	inline
	std::array<double, 3u>
	unit  // peri::
		( std::array<double, 3u> const & orig
		)
	{
		return { (1./magnitude(orig)) * orig };
	}


	//! Geodetic (Lon/Par) angles for local ellipsoid gradient (or up dir)
	inline
	std::pair<double, double>
	anglesLonParOf // Note: units are unimportant since angles are ratios
		( XYZ const & anyVec
			//!< Arbitrary: for geodetic lon/par use gradient at pVec
		)
	{
		// note computations are ratios and are independent of units
		double const & xx = anyVec[0];
		double const & yy = anyVec[1];
		double const & zz = anyVec[2];
		// radius of parallel circle
		double const hh{ std::sqrt(sq(xx) + sq(yy)) };
		// compute conventional lon/par angles
		double lon{ 0. };
		if (! (0. == hh)) // if small hh, somewhat random longitude
		{
			lon = std::atan2(yy, xx);
		}
		double const par{ std::atan2(zz, hh) };
		return { lon, par };
	}

} // [peri]

namespace peri
{
	/*! \brief Local vertical "up" unit direction at Lon/Par(Lat) location.
	 *
	 * Conventional definition with components interpreted as
	 * \arg [0] : component in equator positive toward lon=0, par=0
	 * \arg [1] : component dextrally orthogonal to [2],[0] components
	 * \arg [2] : component orthogonal to equator, positive to North pole
	 */
	inline
	XYZ
	upDirAtLpa  // peri::
		( LPA const & lpa //!< Only Lon,Par are used: Alt is ignored
		)
	{
		double const & lon = lpa[0];
		double const & par = lpa[1];
		return XYZ
			{ std::cos(par) * std::cos(lon)
			, std::cos(par) * std::sin(lon)
			, std::sin(par)
			};
	}


	/*! \brief Container for parameters describing oblate ellipsoidal shape.
	 *
	 * Represents an oblate spheroid of revolution generated by an ellipse
	 * rotated around its minor axis.
	 *
	 * Internal representation involves the semi-major axis, theRadA and
	 * semi-minor axis, theRadB, assuming (theRadB <= theRadA). Interpretively
	 * theRadA is the "equatorial radius" and theRadB is the "polar radius".
	 * (The term "radius" is used equivalent to "semi-axis-magnitude").
	 *
	 * For geodetic applications, typical usage is to create an instance
	 * from equatorial radius and inverse (first) flattening factor via:
	 * \arg fromMajorInvFlat() - Create shape from common geodetic values
	 *
	 * \note There is nothing to prevent constructing instances with polar
	 * radius larger than equatorial one (i.e. a prolate ellipsoid) or by
	 * providing an negative inverse flattening factor. However, doing so
	 * is inconsistent with basic geodesy use-cases.
	 *
	 * Several data members provide values for derived parameters including:
	 * 
	 * Characteristic size:
	 * \arg #theLambda - The geometric mean of the two radii: sqrt(radA*radB)
	 *
	 * Shape coefficients useful in math expressions (ref doc/PerideticMath):
	 * \arg #theMuSqs - The squared radii values ordered by index
	 *
	 * Data normalization:
	 * \arg normalizedShape() - Conforming shape with unit characteristic size
	 */
	struct Shape
	{
		//! Equatorial radius
		double const theRadA{ sNan };

		//! Polar radius
		double const theRadB{ sNan };

		//! Characteristic length (geometric mean: sqrt(theRadA*theRadB))
		double const theLambda{ sNan };

		//! Coefficients describing geometric shape (i.e. {a^2, a^2, b^2})
		std::array<double, 3u> const theMuSqs{ sNan, sNan, sNan };

	private:

		//! Value construction
		inline
		explicit
		Shape  // Shape::
			( double const & radA
				//!< Equatorial semi-axis magnitude
			, double const & radB
				//!< Polar semi-axis magnitude
			)
			: theRadA{ radA }
			, theRadB{ radB }
			, theLambda{ std::sqrt(theRadA * theRadB) }
			, theMuSqs{ sq(theRadA), sq(theRadA), sq(theRadB) }
		{ }

	public:

		//! Create an instance using the semi-major axis and inverse-Flattening
		static
		inline
		Shape
		fromMajorInvFlat  // Shape::
			( double const & equatorialRadius
				//!< Equatorial (semi-major) radius: for Earth~=6.378e6
			, double const & invFlatFactor
				//!< Inverse (first) flattening factor (aka 1/f): for Earth~=298
			)
		{
			double const & aa = equatorialRadius;
			double const ff{ 1. / invFlatFactor };
			double const bb{ (1.-ff) * aa };
			return Shape(aa, bb);
		}

		//! A null instance (nan data member values)
		Shape  // Shape::
			() = default;

		//! Algebraic (mis)closure relative to ellipsoid level surface
		inline
		XYZ
		gradientAt  // Shape::
			( XYZ const & pVecOnSurface
				//!< A point **ON** surface (i.e. assumes 0==funcValueAt(pVec))
			) const
		{
			return
				{ 2. * pVecOnSurface[0] / theMuSqs[0]
				, 2. * pVecOnSurface[1] / theMuSqs[1]
				, 2. * pVecOnSurface[2] / theMuSqs[2]
				};
		}

		//! A shape conformal to this one but with unit characteristic length.
		inline
		Shape
		normalizedShape  // Shape::
			() const
		{
			double const normPerOrig{ 1. / theLambda };
			return Shape(normPerOrig*theRadA, normPerOrig*theRadB);
		}

	}; // Shape


	/*! \brief Merit function to evaluate altitude scaling closure.
	 *
	 * Provides evaluation of ellipsoidal constraint function and
	 * its derivatives (with respect to 'sigma' scale parameter).
	 *
	 * For best numeric stability, construct and utilize with a shape
	 * that has a characteristic size near unity (e.g. radii near 1).
	 *
	 */
	struct ShapeClosure
	{
		//! Parameters describing the underlying shape.
		Shape theShape{};

		//! Value construction.
		inline
		explicit
		ShapeClosure // ShapeClosure::
			( Shape const & shape
			)
			: theShape{ shape }
		{
		}

		//! Default creates a null instance (member values are NaN)
		ShapeClosure // ShapeClosure::
			() = default;

		/*! \brief Ellipsoid constraint function and derivative values.
		 *
		 * Elements are:
		 * \arg [0]: Function value - ellipsoid "misclosure"
		 * \arg [1]: First derivative (with respect to sigma)
		 *
		 */
		 // * \arg [2]: Second derivative (with respect to sigma)
		inline
		std::array<double, 2u>
		funcDerivs // ShapeClosure::
			( double const & sigma
				//!< Free parameter at which to evaluate merit function
			, XYZ const & xVec
				//!< Point of interest location (in same units as theShape)
			) const
		{
			std::array<double, 2u> fdfs;
			XYZ const muPlusSigmas
				{ (theShape.theMuSqs[0] + sigma)
				, (theShape.theMuSqs[1] + sigma)
				, (theShape.theMuSqs[2] + sigma)
				};
			XYZ const muXSqs
				{ theShape.theMuSqs[0] * sq(xVec[0])
				, theShape.theMuSqs[1] * sq(xVec[1])
				, theShape.theMuSqs[2] * sq(xVec[2])
				};
			// function value - for ellipsoid condition equation
			XYZ terms
				{ muXSqs[0] / sq(muPlusSigmas[0])
				, muXSqs[1] / sq(muPlusSigmas[1])
				, muXSqs[2] / sq(muPlusSigmas[2])
				};
			fdfs[0] = (terms[0] + terms[1] + terms[2]) - 1.;
			// first derivative - for ellipsoid condition equation
			terms[0] /= muPlusSigmas[0];
			terms[1] /= muPlusSigmas[1];
			terms[2] /= muPlusSigmas[2];
			fdfs[1] = -2.*(terms[0] + terms[1] + terms[2]);
			/*
			// second derivative - for ellipsoid condition equation
			terms[0] /= muPlusSigmas[0];
			terms[1] /= muPlusSigmas[1];
			terms[2] /= muPlusSigmas[2];
			fdfs[2] = 6.*(terms[0] + terms[1] + terms[2]);
			*/
			return fdfs;
		}

	}; // ShapeClosure


	/*! \brief Math description of ellipsoid surface (as a scalar field)
	 *
	 * Has members:
	 * \arg theShapeOrig - Shape parameters in orig units
	 * \arg theShapeNorm - Shape parameters in normalized units
	 *
	 * Characteristic size (used for normalization and restoration):
	 * \arg lambdaOrig() - Characteristic size parameter is forward from
	 * 		theShapeOrig (is lambda==1 for theShapeNorm)
	 *
	 * Provides: Data normalization/restoration functions:
	 * \arg xyzNormFrom - normalized vector (components ranging [-1,1])
	 * \arg xyzOrigFrom - restored vector (components in physical units [m])
	 *
	 * Evaluations using normalized shape include:
	 * \arg gradientAt() - Vector gradient of shape field
	 *		(proportional to ellipsoid normal direction vector)
	 *
	 * Employed notation includes:
	 * \arg xVecOrig - an arbitrary point in space ('orig' physical units)
	 * \arg xVecNorm - normalized expression for xVec ('norm' units near 1)
	 */
	struct Ellipsoid
	{
		//! Original magnitude shape parameters
		Shape const theShapeOrig{};

		//! Normalized equivalent shape (1==theShapeNorm.theLambda())
		Shape const theShapeNorm{};

		//! A null instance
		Ellipsoid
			() = default;

		//! Value construction
		inline
		explicit
		Ellipsoid  // Ellipsoid::
			( Shape const & shapeOrig
				//!< Ellipsoid shape described by physical units (e.g. [m])
			)
			: theShapeOrig{ shapeOrig }
			, theShapeNorm{ theShapeOrig.normalizedShape() }
		{ }

		//! Characteristic size (geometric mean of original shape semi-axes)
		inline
		double
		lambdaOrig  // Ellipsoid::
			() const
		{
			return theShapeOrig.theLambda;
		}

		//! Cartesian vector normalized to working dimensions
		inline
		XYZ
		xyzNormFrom  // Ellipsoid::
			( XYZ const & xVecOrig
			) const
		{
			double const scl{ 1. / lambdaOrig() };
			return
				{ scl*xVecOrig[0]
				, scl*xVecOrig[1]
				, scl*xVecOrig[2]
				};
		}

		//! Cartesian vector restored to original units
		inline
		XYZ
		xyzOrigFrom  // Ellipsoid::
			( XYZ const & xVecNorm
			) const
		{
			double const scl{ lambdaOrig() };
			return
				{ scl*xVecNorm[0]
				, scl*xVecNorm[1]
				, scl*xVecNorm[2]
				};
		}

	}; // Ellipsoid

	/*! \brief Provide geodetic transforms at Earth scale (units of [m])
	 *
	 * Represents spatial configuration of Earth ellipsoidal shape
	 * and the relevant geometry in vicinity of its surface.
	 *
	 * All linear interface data are expressed in physical units (e.g. [m]).
	 * These data values are normalized before being used in computations
	 * (performed with a unit sized ellipsoid) after which linear values
	 * are restored to expression in physical units for return.
	 *
	 * Methods include:
	 * \arg lpaForXyz() - Geodetic coordinates from Cartesian
	 * \arg xyzForLpa() - Cartesian coordinates from Geodetic
	 * \arg nearEllipsoidPointFor() - Point on ellipsoid nearest point in space
	 */
	struct EarthModel
	{

		//! Geometric representation of surface
		Ellipsoid const theEllip{};

	private:

		//! Mathematical level condition associated with surface
		ShapeClosure const theMeritFuncNorm{}; // normalized units for stability

	public: // Note: public functions interface with physical units

		//! A null instance
		EarthModel
			() = default;

		//! Construct to match physical geometry description
		inline
		explicit
		EarthModel  // EarthModel::
			( Shape const & shape
				//!< Figure of Earth ellipsoid shape: physical units (e.g. [m])
			)
			: theEllip(shape)
			// construct with normalized values for stable numerics
			, theMeritFuncNorm(theEllip.theShapeNorm)
		{ }

		//! Geodetic coordinates associated with Cartesian coordinates xVec
		inline
		LPA
		lpaForXyz  // EarthModel::
			( XYZ const & xLocXyz
			) const
		{
			XYZ const & xVecOrig = xLocXyz;
			// normalize data values to facilitate stable computation
			XYZ const xVecNorm{ theEllip.xyzNormFrom(xVecOrig) };
			// find point, pVec, on ellipsoid closest to world point, xVec
			XYZ const pVecNorm{ poeNormFor(xVecNorm) };
			// compute local vertical direction from gradient
			XYZ const pGrad{ theEllip.theShapeNorm.gradientAt(pVecNorm) };
			XYZ const pUp{ unit(pGrad) };
			// compute altitude as directed distance from ellipsoid at pVec
			double const altNorm{ dot((xVecNorm - pVecNorm), pUp) };
			// rescale altitude to original units
			double const lambdaOrig{ theEllip.lambdaOrig() };
			double const altOrig{ lambdaOrig * altNorm };
			// find point, pVec, on ellipsoid closest to world point, xVec
			// extract LP (at A=0.) from vertical direction at pVec
			std::pair<double, double> const pairLonPar{ anglesLonParOf(pGrad) };
			// angles are invariant to scale (unaffected by normalization)
			double const & pLonOrig = pairLonPar.first;
			double const & pParOrig = pairLonPar.second;
			// return value as combo of LP and A computed results
			return LPA{ pLonOrig, pParOrig, altOrig };
		}

		//! Cartesian coordinates for geodetic location lpa
		inline
		XYZ
		xyzForLpa  // EarthModel::
			( LPA const & xLocLpa
			) const
		{
			double const & alt = xLocLpa[2];
			// determine vertical direction at LP location
			XYZ const up{ upDirAtLpa(xLocLpa) };
			// compute scaling coefficient
			std::array<double, 3u> const & muSqs
				= theEllip.theShapeNorm.theMuSqs;
			double const sumMuUpSq // positive since all mu values are positive
				{ muSqs[0]*sq(up[0])
				+ muSqs[1]*sq(up[1])
				+ muSqs[2]*sq(up[2])
				};
			double const scl{ theEllip.lambdaOrig() / std::sqrt(sumMuUpSq) };
			// compute Cartesian location as displacement along normal dir
			return
				{ (scl*muSqs[0] + alt) * up[0]
				, (scl*muSqs[1] + alt) * up[1]
				, (scl*muSqs[2] + alt) * up[2]
				};
		}

		//! Perpendicular projection (pVec) from xVec onto ellipsoid
		inline
		XYZ
		nearEllipsoidPointFor  // EarthModel::
			( XYZ const & xVecOrig
			) const
		{
			XYZ const xVecNorm{ theEllip.xyzNormFrom(xVecOrig) };
			XYZ const pVecNorm{ poeNormFor(xVecNorm) };
			XYZ const pVecOrig{ theEllip.xyzOrigFrom(pVecNorm) };
			return pVecOrig;
		}

	private: // Note: private functions operate with normalized data units

		//! Initial estimate for sigma factor (based on sphere approximation)
		inline
		double
		sigmaNormWrtSphere  // EarthModel::
			( XYZ const & xVecNorm
			) const
		{
			double const xMagNorm{ std::sqrt(dot(xVecNorm, xVecNorm)) };
			return (xMagNorm - 1.);
		}

		//! A linearly refined improvement to altitude scale factor currSigma
		inline
		double
		nextSigmaNormFor  // EarthModel::
			( double const & currSigmaNorm
			, XYZ const & xVecNorm
			) const
		{
			// evaluate function and derivatives at expansion point
			std::array<double, 2u> const fdfs
				{ theMeritFuncNorm.funcDerivs(currSigmaNorm, xVecNorm) };
			// linear upate
			double const num{ fdfs[0] };
			double const den{ fdfs[1] };
			double const nextSigma{ currSigmaNorm - num/den };
			return nextSigma;
		}

		//! Refined altitude scale factor at normalized point location xVecNorm
		inline
		double
		sigmaNormFor  // EarthModel::
			( XYZ const & xVecNorm
			) const
		{
			// linearized iteration
			double sigmaNorm{ sigmaNormWrtSphere(xVecNorm) };
			double currTestVal{ 1. + sigmaNorm };
			// Convergence is extremely quick within operational range
			// e.g. 3 or 2 iterations typically sufficient
			constexpr std::size_t nnMax{ 8u };
			for (std::size_t nn{0u} ; nn < nnMax ; ++nn)
			{
				sigmaNorm = nextSigmaNormFor(sigmaNorm, xVecNorm);
				double const nextTestVal{ 1. + sigmaNorm };
				// Tolerance suitable for 64-bit double type
				constexpr double tolDiff{ 1.e-15 };
				if (std::abs(currTestVal - nextTestVal) < tolDiff)
				{
					break;
				}
				currTestVal = nextTestVal;
			}
			return sigmaNorm;
		}

		//! Point-on-ellipsoid: pVec = perp projection onto ellipsoid from xVec
		inline
		XYZ
		poeNormFor  // EarthModel::
			( XYZ const & xVecNorm
				//!< Point of interest in normalized coordinates
			) const
		{
			double const sigmaNorm{ sigmaNormFor(xVecNorm) };
			std::array<double, 3u> const & muSqNorms
				= theEllip.theShapeNorm.theMuSqs;
			XYZ const pVecNorm
				{ muSqNorms[0] * xVecNorm[0] / (muSqNorms[0] + sigmaNorm)
				, muSqNorms[1] * xVecNorm[1] / (muSqNorms[1] + sigmaNorm)
				, muSqNorms[2] * xVecNorm[2] / (muSqNorms[2] + sigmaNorm)
				};
			return pVecNorm;
		}

	}; // EarthModel


} // [peri]



// Definitions for functions
// #include "periDetail.inl"

#endif // periDetail_INCL_

