#ifndef POTENTIAL_HARMONICBOND_H_INCLUDED
#define POTENTIAL_HARMONICBOND_H_INCLUDED

#include "Potential.h"

class Potential_HarmonicBond : public Potential {
public:
	Potential_HarmonicBond (
		double k = 1.,
		double r0 = 1.
	) noexcept;

protected:
	double
	ObjectiveFunction (
		double rij
	) const noexcept override;


	Potential::Pair_1stDerivative_t
	_1stDerivative (
		double rij
	) const noexcept override;


	Potential::Pair_2ndDerivative_t
	_2ndDerivative (
		double rij
	) const noexcept override;

private:
	double m_k;
	double m_r0;
};

#endif // POTENTIAL_HARMONICBOND_H_INCLUDED
