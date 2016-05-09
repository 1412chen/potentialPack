#ifndef POTENTIAL_COSINEDIHEDRAL_H_INCLUDED
#define POTENTIAL_COSINEDIHEDRAL_H_INCLUDED

#include "PotentialDihedral.h"


class Potential_CosineDihedral : public PotentialDihedral {
public:
	Potential_CosineDihedral (
		double k,
		int n,
		double d
	) noexcept;

private:
	double
	ObjectiveFunction (
		double cosPhi,
		double /*rij_*/,
		double /*rjk_*/,
		double /*rkl_*/,
		double /*cosThetaIJK*/,
		double /*cosThetaJKL*/
	) const noexcept override;


	PotentialDihedral::_1stDerivative_t
	_1stDerivative_ (
		double cosPhi,
		double /*rij_*/,
		double /*rjk_*/,
		double /*rkl_*/,
		double /*cosThetaIJK*/,
		double /*cosThetaJKL*/
	) const noexcept override;


	PotentialDihedral::_2ndDerivative_t
	_2ndDerivative (
		double cosPhi,
		double /*rij_*/,
		double /*rjk_*/,
		double /*rkl_*/,
		double /*cosThetaIJK*/,
		double /*cosThetaJKL*/
	) const noexcept override;

private:
	double m_k;
	int m_n;
	double m_d;
};

#endif // POTENTIAL_COSINEDIHEDRAL_H_INCLUDED
