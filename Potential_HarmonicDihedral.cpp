#include "Potential_HarmonicDihedral.h"
#include <cmath>

using namespace Dihedral_Namespace;


Potential_HarmonicDihedral::Potential_HarmonicDihedral (
	double k,
	int d,
	unsigned n
) noexcept :
	m_k( k ),
	m_d( d ),
	m_n( n )
{}


double
Potential_HarmonicDihedral::ObjectiveFunction (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	return m_k*( 1. + m_d*CosNPhi(m_n, cosPhi) );
}


PotentialDihedral::_1stDerivative_t
Potential_HarmonicDihedral::_1stDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	_1stDerivative_t _1stDeri;
	_1stDeri[DCosPhi] = m_k*m_d*_1stDCosNPhi(m_n, cosPhi);
	return _1stDeri;
}


PotentialDihedral::_2ndDerivative_t
Potential_HarmonicDihedral::_2ndDerivative (
	double cosPhi,
	double /*rij_*/,
	double /*rjk_*/,
	double /*rkl_*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri[DCosPhi_DCosPhi] = m_k*m_d*_2ndDCosNPhi(m_n, cosPhi);
	return _2ndDeri;
}

