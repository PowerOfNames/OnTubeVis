#pragma once

#include <cgv/math/fmat.h>
#include <cgv/math/quaternion.h>
#include <math.h>

namespace JacobiEigen {
	#define EPSILON 1e-10

	//THESIS: (from chatGPT)
	static std::pair<uint32_t, uint32_t> FindLargestOffDiagonal(const cgv::math::fmat<float, 3, 3>& matrix)
	{
		std::pair<uint32_t, uint32_t> result = { 0, 1 };
		for (uint32_t column = 0; column < 3; column++)
		{
			for (uint32_t row = column + 1; row < 3; row++)
			{
				if (fabs(matrix(row, column)) >= fabs(matrix(result.first, result.second)))
				{
					result.first = row;
					result.second = column;
				}
			}
		}
		//std::printf("Max: %i;%i \n", result.first, result.second);
		return result;
	}
	//THESIS: (from chatGPT, but with adjustment to the theta calculation (switch to matrix(q,q) - matrix(p,p) from matrix(p,p) - matrix(q,q) because of colum-major order) and equivalent
	//			for the eigenvector rotation
	static void JacobiRotation(cgv::math::fmat<float, 3, 3>& matrix, cgv::math::fmat<float, 3, 3>& eigenvectors, uint32_t p, uint32_t q)
	{
		if (fabs(matrix(p, q)) < EPSILON) return; // Already close to zero

		float theta = 0.5 * atan2(2 * matrix(p, q), matrix(q, q) - matrix(p, p));
		//std::printf("Current p: %i current q: %i CurrentRotationAngle: %f\n", p, q, theta);
		//std::printf("CurrentMaxOffDiagValue: %f ", matrix(p,q));


		float c = cos(theta);
		float s = sin(theta);
		// Apply the rotation to the matrix
		float app = c * c * matrix(p, p) - 2 * s * c * matrix(p, q) + s * s * matrix(q, q);
		float aqq = s * s * matrix(p, p) + 2 * s * c * matrix(p, q) + c * c * matrix(q, q);
		float apq = 0.0; // This will become zero

		float temp[3];
		for (int i = 0; i < 3; ++i) {
			if (i != p && i != q) {
				temp[i] = (c * matrix(i, p)) - (s * matrix(i, q));
				matrix(i, q) = (s * matrix(i, p)) + (c * matrix(i, q));
				matrix(i, p) = temp[i];
				matrix(p, i) = matrix(i, p); // Symmetric update
				matrix(q, i) = matrix(i, q); // Symmetric update
			}
		}

		matrix(p, p) = app;
		matrix(q, q) = aqq;
		matrix(p, q) = apq;
		matrix(q, p) = apq;

		// Rotate eigenvectors
		for (int i = 0; i < 3; ++i) {
			temp[i] = c * eigenvectors(p, i) - s * eigenvectors(q, i);
			eigenvectors(q, i) = s * eigenvectors(p, i) + c * eigenvectors(q, i);
			eigenvectors(p, i) = temp[i];
		}
	}

	static void PrintMatrix(const cgv::math::fmat<float, 3, 3>& matrix)
	{
		//std::printf("CurrentMatrixRotation: \n");
		for (uint32_t row = 0; row < 3; row++)
		{
			std::printf("{");
			for (uint32_t column = 0; column < 3; column++)
				std::printf(" %f ", matrix(row, column));
			std::printf("}\n");
		}
		std::printf("\n");
	}
	static void CheckEigen(const cgv::math::fmat<float, 3, 3>& matrix, const cgv::math::fmat<float, 3, 3>& eigenvectors, const cgv::math::fvec<float, 3>& eigenvalues)
	{
		for (int i = 0; i < 3; ++i) {
			double lambda = eigenvalues[i];
			double v[3] = { eigenvectors(i,0), eigenvectors(i,1), eigenvectors(i,2) };
			double Av[3] = { 0, 0, 0 };

			// Multiply A by v
			for (int j = 0; j < 3; ++j) {
				Av[j] = matrix(j, 0) * v[0] + matrix(j, 1) * v[1] + matrix(j, 2) * v[2];
			}

			// Check if Av is approximately equal to lambda * v
			for (int j = 0; j < 3; ++j) {
				if (fabs(Av[j] - lambda * v[j]) > EPSILON) {
					std::cout << "Eigenvector check failed for eigenvalue " << lambda << std::endl;
				}
				else {
					std::cout << "Eigenvalue " << lambda << " validated with eigenvector." << std::endl;
				}
			}
		}
	}
	//THESIS: (from chatGPT)	
	static void JacobiEigen(cgv::math::fmat<float, 3, 3>& matrix, cgv::math::fmat<float, 3, 3>& eigenvectors, cgv::math::fvec<float, 3>& eigenvalues)
	{
		//Eigenvector initialization to identity matrix
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				eigenvectors(i, j) = (i == j) ? 1.0 : 0.0;
			}
		}

		// Iterate until off-diagonal elements are sufficiently small
		for (int iter = 0; iter < 100; ++iter) {
			auto [p, q] = FindLargestOffDiagonal(matrix);
			if (fabs(matrix(p, q)) < EPSILON) break; // Converged
			JacobiRotation(matrix, eigenvectors, p, q);
			//PrintMatrix(matrix);
		}

		// Extract eigenvalues from the diagonal
		for (int i = 0; i < 3; ++i) {
			eigenvalues[i] = matrix(i, i);
		}
	}

	static void NormalizeEigenvectors(cgv::math::fmat<float, 3, 3>& eigenvectors)
	{
		for (int i = 0; i < 3; ++i) {
			double norm = sqrt(eigenvectors(i, 0) * eigenvectors(i, 0) +
				eigenvectors(i, 1) * eigenvectors(i, 1) +
				eigenvectors(i, 2) * eigenvectors(i, 2));
			for (int j = 0; j < 3; ++j) {
				eigenvectors(i, j) /= norm; // Normalize the eigenvector
			}
		}
	}

	static void CalculateAnglesFromEigenVectors(cgv::math::fmat<float, 3, 3>& eigenvectors, cgv::math::fvec<float, 3>& angles)
	{
		constexpr float PI = 3.1415926535;
		float rad1 = atan2(eigenvectors(1, 0), eigenvectors(0, 0));
		float rad2 = asin(-eigenvectors(2, 0));
		float rad3 = atan2(eigenvectors(2, 1), eigenvectors(2, 2));
		angles[0] = rad1 < 0.0 ? rad1 + 2 * PI : rad1; //yaw in x/y plane
		angles[0] /= 2 * PI;
		angles[1] = rad2 + PI; //pitch - tilt from z axis
		angles[1] /= PI;
		angles[2] = rad3 < 0.0 ? rad2 + 2 * PI : rad3;//roll
		angles[2] /= 2 * PI;

		//std::printf("Yaw: %f; Pitch: %f; Roll: %f \n", angles[0], angles[1], angles[2]);
		//std::printf("Yaw: %f; Pitch: %f; Roll: %f \n", (angles[0] * 180.0) / PI, (angles[1] * 180.0) / PI, (angles[2] * 180.0) / PI);
	}

	static cgv::math::quaternion<float> CalculateQuaternionFromEigenVectors(cgv::math::fmat<float, 3, 3>& eigenvectors)
	{
		return cgv::math::quaternion<float>(eigenvectors);
	}

	static void PrintQuaternion(cgv::math::quaternion<float> quat)
	{
		std::printf("Eigen-Quaternion: r:%f, i:%f, j:%f, k:%f \n", quat[0], quat[1], quat[2], quat[3]);
	}

}
