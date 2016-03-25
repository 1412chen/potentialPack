#include "Potential_HarmonicAngle.h"
#include <cmath>

using namespace Angle_Namespace;

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


_1stDerivative_t
Potential_HarmonicAngle::_1stDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;

	_1stDerivative_t _1stDeri;
	_1stDeri.fill(0.);
	_1stDeri[DCosTheta] = -2.*m_k*dTheta/sin(theta);
	return _1stDeri;
}


_2ndDerivative_t
Potential_HarmonicAngle::_2ndDerivative (
	double cosTheta,
	double /*rij*/,
	double /*rik*/
) const noexcept
{
	auto sinThetaSq = 1. - cosTheta*cosTheta;
	auto theta = acos(cosTheta);
	auto dTheta = theta - m_theta0;
	auto cotTheta = cosTheta / sqrt(sinThetaSq);

	_2ndDerivative_t _2ndDeri;
	_2ndDeri.fill(0.);
	_2ndDeri[DCosTheta_DCosTheta] = 2.*m_k/sinThetaSq * (1.-dTheta*cotTheta);
	return _2ndDeri;
}

