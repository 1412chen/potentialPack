#include "Potential_HarmonicBond.h"

using namespace std;


Potential_HarmonicBond::Potential_HarmonicBond (
	double k,
	double r0
) noexcept :
	m_k( k ),
	m_r0( r0 )
{}


double
Potential_HarmonicBond::ObjectiveFunction (
	double rij
) const noexcept
{
	auto dr = rij - m_r0;
	return 0.5 * m_k * dr * dr;
}


Potential::Pair_1stDerivative_t
Potential_HarmonicBond::_1stDerivative (
	double rij
) const noexcept
{
	return m_k * (rij - m_r0);
}


Potential::Pair_2ndDerivative_t
Potential_HarmonicBond::_2ndDerivative (
	double /*rij*/
) const noexcept
{
	return m_k;
}

