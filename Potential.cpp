#include "Potential.h"
#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;

double
Potential::Energy (
	const array3d& /*rij*/
) const noexcept
{
	return 0.;
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return 0.;
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return 0.;
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const std::vector<array3d>& /*rik*/
) const noexcept
{
	return 0.;
}

//---------------------------------------------------------------------------//

array3d
Potential::Force (
	const array3d& /*rij*/
) const noexcept
{
	return array3d{};
}


array3d
Potential::Force (
	const array3d& /*rij*/,
	FiniteDifference_t
) const noexcept
{
	return array3d{};
}


std::array<array3d, 2>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return {array3d{}, array3d{}};
}


std::array<array3d, 2>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rik*/,
	FiniteDifference_t
) const noexcept
{
	return {array3d{}, array3d{}};
}


array<array3d, 3>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return {array3d{}, array3d{}, array3d{}};
}


array<array3d, 3>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/,
	FiniteDifference_t
) const noexcept
{
	return {array3d{}, array3d{}, array3d{}};
}


vector<array3d>
Potential::Force (
	const array3d& /*rij*/,
	const vector<array3d>& rik
) const noexcept
{
	return vector<array3d>(rik.size()+1, array3d{});
}


vector<array3d>
Potential::Force (
	const array3d& /*rij*/,
	const vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	return vector<array3d>(rik.size()+1, array3d{});
}

//---------------------------------------------------------------------------//

matrix3d
Potential::Hessian (
	const array3d& /*rij*/
) const noexcept
{
	return matrix3d{};
}


matrix3d
Potential::Hessian (
	const array3d& /*rij*/,
	FiniteDifference_t
) const noexcept
{
	return matrix3d{};
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return {matrix3d{}, matrix3d{}, matrix3d{}};
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rik*/,
	FiniteDifference_t
) const noexcept
{
	return {matrix3d{}, matrix3d{}, matrix3d{}};
}


std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return {matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}};
}


std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/,
	FiniteDifference_t
) const noexcept
{
	return {matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}, matrix3d{}};
}


std::vector<matrix3d>
Potential::Hessian (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik
) const noexcept
{
	return vector<matrix3d>(rik.size()+1, matrix3d{});
}


std::vector<matrix3d>
Potential::Hessian (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	return vector<matrix3d>(rik.size()+1, matrix3d{});
}

