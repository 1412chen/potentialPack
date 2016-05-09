#include "Potential_HarmonicAngle.h"
#include <cmath>

using namespace Angle_Namespace;
constexpr double FiniteDistance = 1.e-13;
constexpr double InverseFiniteDistance = 1./FiniteDistance;


Potential_HarmonicAngle::Potential_HarmonicAngle (
	double k,
	double theta0
) noexcept :
	m_k( k ),
	m_theta0( theta0 )
{}


double
Potential_HarmonicAngle::ObjectiveFunction (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	return m_k*dTheta*dTheta;
}


PotentialAngle::_1stDerivative_t
Potential_HarmonicAngle::_1stDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	if ( 1.-fabs(cosTheta) < FiniteDistance )
		_1stDeri[DCosTheta] = 2.*m_k*m_theta0 * InverseFiniteDistance;
	else
	{
		auto theta = acos(cosTheta);
		auto dTheta = theta - m_theta0;

		_1stDeri[DCosTheta] = -2.*m_k*dTheta/sin(theta);
	}
	return _1stDeri;
}


PotentialAngle::_2ndDerivative_t
Potential_HarmonicAngle::_2ndDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	if ( 1.-fabs(cosTheta) < FiniteDistance )
		_2ndDeri[DCosTheta_DCosTheta] =  2.*m_k*(1.+m_theta0*cosTheta) * InverseFiniteDistance;
	else
	{
		auto sinThetaSq = 1. - cosTheta*cosTheta;
		auto theta = acos(cosTheta);
		auto dTheta = theta - m_theta0;
		auto cotTheta = cosTheta / sqrt(sinThetaSq);

		_2ndDeri[DCosTheta_DCosTheta] = 2.*m_k/sinThetaSq * (1.-dTheta*cotTheta);
	}
	return _2ndDeri;
}

