#include "Potential.h"
#include <utility>
#include <cmath>
#include <iostream>


constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;

double
Potential::Energy (
	const array3d& /*rij*/
) const noexcept
{
	return ZeroEnergy();
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return ZeroEnergy();
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return ZeroEnergy();
}


double
Potential::Energy (
	const array3d& /*rij*/,
	const std::vector<array3d>& /*rik*/
) const noexcept
{
	return ZeroEnergy();
}

//---------------------------------------------------------------------------//

array3d
Potential::Force (
	const array3d& /*rij*/
) const noexcept
{
	return ZeroForce( Pair_t() );
}


array3d
Potential::Force (
	const array3d& /*rij*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroForce( Pair_t() );
}


std::array<array3d, 2>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return ZeroForce( Angle_t() );
}


std::array<array3d, 2>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rik*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroForce( Angle_t() );
}


std::array<array3d, 3>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return ZeroForce( Dihedral_t() );
}


std::array<array3d, 3>
Potential::Force (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroForce( Dihedral_t() );
}


std::vector<array3d>
Potential::Force (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik
) const noexcept
{
	return ZeroForce( Manybody_t(), rik.size()+1 );
}


std::vector<array3d>
Potential::Force (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	return ZeroForce( Manybody_t(), rik.size()+1 );
}

//---------------------------------------------------------------------------//

matrix3d
Potential::Hessian (
	const array3d& /*rij*/
) const noexcept
{
	return ZeroHessian( Pair_t() );
}


matrix3d
Potential::Hessian (
	const array3d& /*rij*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroHessian( Pair_t() );
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rik*/
) const noexcept
{
	return ZeroHessian( Angle_t() );
}


std::array<matrix3d, 3>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rik*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroHessian( Angle_t() );
}


std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/
) const noexcept
{
	return ZeroHessian( Dihedral_t() );
}


std::array<matrix3d, 6>
Potential::Hessian (
	const array3d& /*rij*/,
	const array3d& /*rjk*/,
	const array3d& /*rkl*/,
	FiniteDifference_t
) const noexcept
{
	return ZeroHessian( Dihedral_t() );
}


std::vector<matrix3d>
Potential::Hessian (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik
) const noexcept
{
	return ZeroHessian( Manybody_t(), rik.size()+1 );
}


std::vector<matrix3d>
Potential::Hessian (
	const array3d& /*rij*/,
	const std::vector<array3d>& rik,
	FiniteDifference_t
) const noexcept
{
	return ZeroHessian( Manybody_t(), rik.size()+1 );
}

//---------------------------------------------//

constexpr double
Potential::ZeroEnergy () const noexcept
{
	return 0.;
}


constexpr array3d
Potential::ZeroForce () const noexcept
{
	return array3d{0., 0., 0.};
}


constexpr matrix3d
Potential::ZeroHessian () const noexcept
{
	return matrix3d {
		ZeroForce(), ZeroForce(), ZeroForce()
	};
}

//---------------------------------------------//

constexpr array3d
Potential::ZeroForce (
	Pair_t
) const noexcept
{
	return ZeroForce();
}


constexpr std::array<array3d, 2>
Potential::ZeroForce (
	Angle_t
) const noexcept
{
	return {ZeroForce(), ZeroForce()};
}

constexpr std::array<array3d, 3>
Potential::ZeroForce (
	Dihedral_t
) const noexcept
{
	return {ZeroForce(), ZeroForce(), ZeroForce()};
}


std::vector<array3d>
Potential::ZeroForce (
	Manybody_t,
	unsigned size
) const noexcept
{
	return std::vector<array3d>(size, ZeroForce());
}


constexpr matrix3d
Potential::ZeroHessian (
	Pair_t
) const noexcept
{
	return ZeroHessian();
}


constexpr std::array<matrix3d, 3>
Potential::ZeroHessian (
	Angle_t
) const noexcept
{
	return {ZeroHessian(), ZeroHessian(), ZeroHessian()};
}


constexpr std::array<matrix3d, 6>
Potential::ZeroHessian (
	Dihedral_t
) const noexcept
{
	return {ZeroHessian(), ZeroHessian(), ZeroHessian(),
		ZeroHessian(), ZeroHessian(), ZeroHessian()};
}


std::vector<matrix3d>
Potential::ZeroHessian (
	Manybody_t,
	unsigned size
) const noexcept
{
	return std::vector<matrix3d>(size, ZeroHessian());
}

