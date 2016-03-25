#ifndef POTENTIAL_HARMONICANGLE_H_INCLUDED
#define POTENTIAL_HARMONICANGLE_H_INCLUDED

#include "Potential.h"

class Potential_HarmonicAngle : public Potential {
public:
	Potential_HarmonicAngle (
		double k = 1.,
		double theta0 = 0.5
	) noexcept;

protected:
	double
	ObjectiveFunction (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;


	Angle_Namespace::_1stDerivative_t
	_1stDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;


	Angle_Namespace::_2ndDerivative_t
	_2ndDerivative (
		double cosTheta,
		double rij,
		double rik
	) const noexcept override;

private:
	double m_k;
	double m_theta0;
};

#endif // POTENTIAL_HARMONICANGLE_H_INCLUDED
