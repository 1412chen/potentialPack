#ifndef POTENTIAL_BUCKINGHAM_H_INCLUDED
#define POTENTIAL_BUCKINGHAM_H_INCLUDED

#include "Potential.h"

class Potential_Buckingham : public Potential {
public:
	Potential_Buckingham (
		double A = 0.13169812355066,
		double rho = 2.993,
		double C = 12.
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
	double
	R6 (
		double rij
	) const noexcept;


	double
	R8 (
		double rij
	) const noexcept;

private:
	double m_A;
	double m_B;
	double m_C;
};

#endif // POTENTIAL_BUCKINGHAM_H_INCLUDED
