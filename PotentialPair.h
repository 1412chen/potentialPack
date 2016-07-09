#ifndef POTENTIALPAIR_H_INCLUDED
#define POTENTIALPAIR_H_INCLUDED

#include "Potential.h"
#include <functional>


class PotentialPair : public Potential {
public:
	using _1stDerivative_t = Pair_Namespace::_1stDerivative_t;
	using _2ndDerivative_t = Pair_Namespace::_2ndDerivative_t;
	using _1stD_Func = std::function< _1stDerivative_t(double) >;
	using _2ndD_Func = std::function< _2ndDerivative_t(double) >;


	double
	Energy (
		const array3d& rij
	) const noexcept override;

	//-------------------------------------//

	array3d // fij
	Force (
		const array3d& rij
	) const noexcept override;


	array3d
	Force (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept override;

	//-------------------------------------//

	matrix3d // Hij
	Hessian (
		const array3d& rij
	) const noexcept override;


	matrix3d
	Hessian (
		const array3d& rij,
		FiniteDifference_t
	) const noexcept override;

private:
	array3d
	ForceImp (
		const array3d& rij,
		_1stD_Func _1stDeriFunc
	) const noexcept;


	matrix3d
	HessianImp (
		const array3d& rij,
		_1stD_Func _1stDeriFunc,
		_2ndD_Func _2ndDeriFunc
	) const noexcept;

	//-------------------------------------//
	//-------------------------------------//
	//-------------------------------------//

	virtual double
	ObjectiveFunction (
		double /*rij_*/
	) const noexcept = 0;


	virtual _1stDerivative_t // drij
	_1stDerivative (
		double rij_
	) const noexcept;


	virtual _2ndDerivative_t
	_2ndDerivative (
		double rij_
	) const noexcept;
};

#endif // POTENTIALPAIR_H_INCLUDED
