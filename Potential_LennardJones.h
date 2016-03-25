#ifndef POTENTIAL_LENNARDJONES_H_INCLUDED
#define POTENTIAL_LENNARDJONES_H_INCLUDED

#include "Potential.h"

class Potential_LennardJones : public Potential {
public:
	Potential_LennardJones (
		double epsilon = 1.,
		double sigma = 1.
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

	//-------------------------------------//

	double
	Factor6 (
		double inverseRij
	) const noexcept;


	double
	Factor6 (
		double inverseRijSq,
		std::nullptr_t
	) const noexcept;

private:
	double m_epsilon;
	double m_sigma;
};

#endif // POTENTIAL_LENNARDJONES_H_INCLUDED
