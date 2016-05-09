#include "Potential_CosineDihedral.h"


Potential_CosineDihedral::Potential_CosineDihedral (
	double k,
	int n,
	double d
) noexcept :
	m_k( k ),
	m_n( n ),
	m_d( d )
{}


double
Potential_CosineDihedral::ObjectiveFunction (
	double cosPhi,
	double /*rij*/,
	double /*rjk*/,
	double /*rkl*/,
	double /*cosThetaIJK*/,
	double /*cosThetaJKL*/
) const noexcept
{
	return m_k*( 1.+cosNPhi_d(cosPhi) );
}

